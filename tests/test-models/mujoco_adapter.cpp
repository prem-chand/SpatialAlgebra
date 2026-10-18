#include "mujoco_adapter.h"

#include <mujoco/mujoco.h>

#include <algorithm>
#include <fcntl.h>
#include <limits>
#include <unistd.h>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace test_models {

// -----------------------------------------------------------------------
// Internal implementation (PIMPL)
// -----------------------------------------------------------------------
struct MujocoAdapter::Impl {
    mjModel* model = nullptr;
    mjData* data = nullptr;
    std::vector<int> dof_body_id; // MuJoCo body ID per model DOF
    bool fkDirty = true;

    ~Impl() {
        if (data) mj_deleteData(data);
        if (model) mj_deleteModel(model);
    }
};

// -----------------------------------------------------------------------
// Helper: build MJCF XML string from RobotModel
// -----------------------------------------------------------------------
namespace {

void emitBodyXML(std::stringstream& ss, const RobotModel& model,
                 int idx, const std::vector<std::vector<int>>& children,
                 const std::string& indent) {
    const auto& js = model.joints[idx];

    // Translation and rotation from parentToJoint
    Eigen::Vector3d t = js.parentToJoint.block<3,1>(0, 3);
    Eigen::Matrix3d R = js.parentToJoint.block<3,3>(0, 0);
    Eigen::Quaterniond quat(R);

    ss << indent << "<body name=\"link_" << idx << "\""
       << " pos=\"" << t(0) << " " << t(1) << " " << t(2) << "\""
       << " quat=\"" << quat.w() << " " << quat.x() << " " << quat.y() << " " << quat.z() << "\">\n";

    if (js.type != JointType::FIXED) {
        std::string jointType = js.type == JointType::REVOLUTE ? "hinge" : "slide";
        ss << indent << "  <joint name=\"j" << idx << "\""
           << " type=\"" << jointType << "\""
           << " axis=\"" << js.jointAxis(0) << " " << js.jointAxis(1) << " " << js.jointAxis(2) << "\"/>\n";
    }

    // Inertial (mass, COM, rotational inertia about COM)
    const auto& I = js.inertia;
    if (js.mass > 1e-6) {
        ss << indent << "  <inertial"
           << " mass=\"" << js.mass << "\""
           << " pos=\"" << js.com(0) << " " << js.com(1) << " " << js.com(2) << "\""
           << " fullinertia=\""
           << I(0,0) << " " << I(1,1) << " " << I(2,2) << " "
           << I(0,1) << " " << I(0,2) << " " << I(1,2) << "\"/>\n";
    }

    for (int child : children[idx]) {
        emitBodyXML(ss, model, child, children, indent + "  ");
    }

    ss << indent << "</body>\n";
}

std::string buildMJCF(const RobotModel& model) {
    int n = static_cast<int>(model.joints.size());
    if (n == 0) {
        return "<mujoco model=\"empty\"><worldbody/></mujoco>";
    }

    // Build tree from parent indices
    std::vector<std::vector<int>> children(n);
    std::vector<int> roots;
    for (int i = 0; i < n; i++) {
        int p = model.joints[i].parent;
        if (p == -1) {
            roots.push_back(i);
        } else {
            if (p >= i)
                throw std::invalid_argument(
                    "RobotModel not in topological order: parent "
                    + std::to_string(p) + " >= child " + std::to_string(i));
            children[p].push_back(i);
        }
    }

    std::stringstream ss;
    ss.precision(std::numeric_limits<double>::max_digits10);
    ss << "<mujoco model=\"sa_robot\">\n";
    // MuJoCo rejects massless moving bodies, such as the virtual floating-base links.
    ss << "  <compiler boundmass=\"1e-12\" boundinertia=\"1e-12\"/>\n";
    ss << "  <worldbody>\n";
    for (int root : roots) {
        emitBodyXML(ss, model, root, children, "    ");
    }
    ss << "  </worldbody>\n";
    ss << "</mujoco>\n";
    return ss.str();
}

} // anonymous namespace

// -----------------------------------------------------------------------
// Factory
// -----------------------------------------------------------------------
std::unique_ptr<MujocoAdapter> MujocoAdapter::build(const RobotModel& model) {
    auto adapter = std::unique_ptr<MujocoAdapter>(new MujocoAdapter());
    auto impl = std::make_unique<Impl>();

    mjModel* m = nullptr;
    char error[2000] = {};

    // MuJoCo 3.4's mj_loadXML only loads from file (not XML string),
    // so we write the MJCF to a temp file first.
    char tmpname[] = "/tmp/sa_mujoco_XXXXXX";
    int fd = mkstemp(tmpname);
    if (fd == -1) {
        throw std::runtime_error("MujocoAdapter: failed to create temp file");
    }
    std::string xml_str = buildMJCF(model);
    ssize_t written = write(fd, xml_str.c_str(), xml_str.size());
    close(fd);
    if (written != static_cast<ssize_t>(xml_str.size())) {
        unlink(tmpname);
        throw std::runtime_error("MujocoAdapter: failed to write temp file");
    }

    m = mj_loadXML(tmpname, nullptr, error, 2000);
    unlink(tmpname);

    if (!m) {
        throw std::runtime_error(
            "MujocoAdapter: failed to build MuJoCo model:\n" +
            std::string(error) + "\nGenerated MJCF:\n" + xml_str);
    }
    impl->model = m;

    impl->data = mj_makeData(impl->model);

    // Build joint-to-DOF mapping
    int n = static_cast<int>(model.joints.size());
    impl->dof_body_id.reserve(static_cast<size_t>(model.getDOF()));
    for (int i = 0; i < n; ++i) {
        std::string name = "link_" + std::to_string(i);
        int bodyId = mj_name2id(impl->model, mjOBJ_BODY, name.c_str());
        if (bodyId < 0) {
            throw std::runtime_error("MujocoAdapter: generated body not found: " + name);
        }
        if (model.joints[static_cast<size_t>(i)].type != JointType::FIXED) {
            impl->dof_body_id.push_back(bodyId);
        }
    }
    if (impl->model->nv != model.getDOF()) {
        throw std::runtime_error(
            "MujocoAdapter: generated model DOF count does not match RobotModel");
    }

    adapter->impl_ = std::move(impl);
    return adapter;
}

MujocoAdapter::~MujocoAdapter() = default;

// -----------------------------------------------------------------------
// State
// -----------------------------------------------------------------------
void MujocoAdapter::setState(const Eigen::VectorXd& q,
                             const Eigen::VectorXd& qdot) {
    int n = getDOF();
    if (q.size() != n || qdot.size() != n) {
        throw std::invalid_argument("MujocoAdapter::setState: size mismatch");
    }
    // MuJoCo uses qpos for all position coordinates.
    // For a fixed-base model with 1-DOF joints, nq == nv == n.
    for (int i = 0; i < n; i++) {
        impl_->data->qpos[i] = q(i);
        impl_->data->qvel[i] = qdot(i);
    }
    impl_->fkDirty = true;
}

void MujocoAdapter::forwardKinematics() {
    mj_forward(impl_->model, impl_->data);
    impl_->fkDirty = false;
}

// -----------------------------------------------------------------------
// Inverse dynamics (RNEA)
// -----------------------------------------------------------------------
Eigen::VectorXd MujocoAdapter::computeTorques(
    const Eigen::VectorXd& qddot,
    const Eigen::Vector3d& gravity) {
    int n = getDOF();
    if (qddot.size() != n) {
        throw std::invalid_argument(
            "MujocoAdapter::computeTorques: size mismatch");
    }

    // Set gravity in model
    impl_->model->opt.gravity[0] = gravity(0);
    impl_->model->opt.gravity[1] = gravity(1);
    impl_->model->opt.gravity[2] = gravity(2);

    // Zero applied forces (so they don't pollute inverse dynamics)
    for (int i = 0; i < n; i++) {
        impl_->data->qfrc_applied[i] = 0.0;
    }

    // Update kinematics: mj_forward recomputes position-dependent quantities.
    // IMPORTANT: mj_forward OVERWRITES qacc with forward-dynamics result,
    // so we must set qacc AFTER this call.
    mj_forward(impl_->model, impl_->data);

    // Set desired accelerations for inverse dynamics
    for (int i = 0; i < n; i++) {
        impl_->data->qacc[i] = qddot(i);
    }

    // mj_rne with flg_acc=1 computes M(q)*qacc + C(q,qvel)
    std::vector<mjtNum> result_buf(n, 0.0);
    mj_rne(impl_->model, impl_->data, 1, result_buf.data());

    Eigen::VectorXd tau(n);
    for (int i = 0; i < n; i++) {
        tau(i) = result_buf[i];
    }
    return tau;
}

// -----------------------------------------------------------------------
// Forward dynamics (ABA)
// -----------------------------------------------------------------------
Eigen::VectorXd MujocoAdapter::computeAccelerations(
    const Eigen::VectorXd& tau,
    const Eigen::Vector3d& gravity) {
    int n = getDOF();
    if (tau.size() != n) {
        throw std::invalid_argument(
            "MujocoAdapter::computeAccelerations: size mismatch");
    }

    impl_->model->opt.gravity[0] = gravity(0);
    impl_->model->opt.gravity[1] = gravity(1);
    impl_->model->opt.gravity[2] = gravity(2);

    // Zero acceleration before forward dynamics
    for (int i = 0; i < n; i++) {
        impl_->data->qacc[i] = 0.0;
    }

    // Set applied torques
    for (int i = 0; i < n; i++) {
        impl_->data->qfrc_applied[i] = tau(i);
    }

    // mj_forward computes forward dynamics (qacc from applied forces)
    mj_forward(impl_->model, impl_->data);

    Eigen::VectorXd qddot(n);
    for (int i = 0; i < n; i++) {
        qddot(i) = impl_->data->qacc[i];
    }
    return qddot;
}

// -----------------------------------------------------------------------
// Mass matrix
// -----------------------------------------------------------------------
Eigen::MatrixXd MujocoAdapter::computeMassMatrix() {
    if (impl_->fkDirty) mj_forward(impl_->model, impl_->data);
    int n = getDOF();

    // MuJoCo stores the joint-space inertia matrix in d->qM (sparse format).
    // mj_fullM expands it to a full n×n matrix.
    std::vector<double> full(n * n, 0.0);
    mj_fullM(impl_->model, full.data(), impl_->data->qM);

    Eigen::MatrixXd H(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            H(i, j) = full[i * n + j];
        }
    }
    return H;
}

// -----------------------------------------------------------------------
// Gravity torques
// -----------------------------------------------------------------------
Eigen::VectorXd MujocoAdapter::computeGravityTorques(
    const Eigen::Vector3d& gravity) {
    int n = getDOF();
    impl_->model->opt.gravity[0] = gravity(0);
    impl_->model->opt.gravity[1] = gravity(1);
    impl_->model->opt.gravity[2] = gravity(2);

    // mj_rne reads the body velocities that mj_forward derives from qvel,
    // so qvel must be zero before mj_forward, not after it.
    std::vector<mjtNum> qvel(impl_->data->qvel, impl_->data->qvel + n);
    for (int i = 0; i < n; i++) {
        impl_->data->qfrc_applied[i] = 0.0;
        impl_->data->qvel[i] = 0.0;
    }
    mj_forward(impl_->model, impl_->data);

    std::vector<mjtNum> result_buf(n, 0.0);
    mj_rne(impl_->model, impl_->data, 0, result_buf.data());

    std::copy(qvel.begin(), qvel.end(), impl_->data->qvel);
    impl_->fkDirty = true;

    Eigen::VectorXd tau(n);
    for (int i = 0; i < n; i++) {
        tau(i) = result_buf[i];
    }
    return tau;
}

// -----------------------------------------------------------------------
// Jacobian
// -----------------------------------------------------------------------
Eigen::MatrixXd MujocoAdapter::computeJointSpaceJacobian(int idx) {
    int n = getDOF();
    if (idx < 0 || idx >= n) {
        throw std::out_of_range(
            "MujocoAdapter::computeJointSpaceJacobian: index out of range");
    }
    if (impl_->fkDirty) mj_forward(impl_->model, impl_->data);

    int body = impl_->dof_body_id[static_cast<size_t>(idx)];
    if (body < 0) {
        throw std::runtime_error("MujocoAdapter: body not found for Jacobian");
    }

    // Allocate temporaries (MuJoCo expects float* for jacp/jacr)
    std::vector<double> jacp(3 * n, 0.0);
    std::vector<double> jacr(3 * n, 0.0);

    // The world origin as the point gives the linear part of the spatial velocity.
    const mjtNum origin[3] = {0.0, 0.0, 0.0};
    mj_jac(impl_->model, impl_->data, jacp.data(), jacr.data(), origin, body);

    Eigen::MatrixXd J(6, n);
    J.setZero();
    // MuJoCo jacp = 3×n translational, jacr = 3×n rotational (column-major)
    for (int j = 0; j < n; j++) {
        J(0, j) = jacr[0 * n + j];
        J(1, j) = jacr[1 * n + j];
        J(2, j) = jacr[2 * n + j];
        J(3, j) = jacp[0 * n + j];
        J(4, j) = jacp[1 * n + j];
        J(5, j) = jacp[2 * n + j];
    }
    return J;
}

// -----------------------------------------------------------------------
// Joint transform
// -----------------------------------------------------------------------
Eigen::Matrix4d MujocoAdapter::getJointTransform(int idx) const {
    if (idx < 0 || idx >= getDOF()) {
        throw std::out_of_range(
            "MujocoAdapter::getJointTransform: index out of range");
    }
    int body = impl_->dof_body_id[static_cast<size_t>(idx)];
    if (body < 0) {
        throw std::runtime_error("MujocoAdapter: body not found");
    }

    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    for (int i = 0; i < 3; i++) {
        T(i, 3) = impl_->data->xpos[3 * body + i];
        for (int j = 0; j < 3; j++) {
            T(i, j) = impl_->data->xmat[9 * body + 3 * i + j];
        }
    }
    return T;
}

// -----------------------------------------------------------------------
// COM in world frame
// -----------------------------------------------------------------------
Eigen::Vector3d MujocoAdapter::getLinkCOM(int idx) const {
    if (idx < 0 || idx >= getDOF()) {
        throw std::out_of_range(
            "MujocoAdapter::getLinkCOM: index out of range");
    }
    int body = impl_->dof_body_id[static_cast<size_t>(idx)];
    if (body < 0) {
        throw std::runtime_error("MujocoAdapter: body not found");
    }

    // xipos = world-frame COM position of each body
    return Eigen::Vector3d(
        impl_->data->xipos[3 * body + 0],
        impl_->data->xipos[3 * body + 1],
        impl_->data->xipos[3 * body + 2]);
}

// -----------------------------------------------------------------------
// DOF
// -----------------------------------------------------------------------
int MujocoAdapter::getDOF() const {
    return impl_->model ? impl_->model->nv : 0;
}

} // namespace test_models
