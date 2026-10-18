#include "mjcf_parser.h"

#include <mujoco/mujoco.h>

#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace test_models {

RobotModel buildRobotModelFromMJCF(const std::string& xml_path) {
    char error[2000] = {};
    mjModel* m_raw = mj_loadXML(xml_path.c_str(), nullptr, error, 2000);
    if (!m_raw) {
        throw std::runtime_error(
            "buildRobotModelFromMJCF: failed to load '" + xml_path +
            "': " + error);
    }
    std::unique_ptr<mjModel, decltype(&mj_deleteModel)> m(m_raw, mj_deleteModel);

    int nb = m->nbody;
    int nj = m->njnt;

    // ----------------------------------------------------------------
    // Step 2: Detect floating base
    // ----------------------------------------------------------------
    // If body 1 has a free joint (mjJNT_FREE), the model is floating-base.
    // We create 6 virtual 1-DOF joints to represent the 6 base DOF:
    //   idx 0: prismatic X   idx 3: revolute X
    //   idx 1: prismatic Y   idx 4: revolute Y
    //   idx 2: prismatic Z   idx 5: revolute Z  (carries base body mass)
    // Actuated joints are shifted to idx 6+.
    int mj_floating_body = -1;  // MuJoCo body index of the floating base
    for (int j = 0; j < nj; j++) {
        if (m->jnt_type[j] == mjJNT_FREE) {
            mj_floating_body = m->jnt_bodyid[j];
            break;
        }
    }
    bool floating = (mj_floating_body > 0);

    // ----------------------------------------------------------------
    // Step 3: Build body → RobotModel index map
    // ----------------------------------------------------------------
    // Shift: 6 for floating base, 0 for fixed base
    std::vector<int> body_to_rmodel(nb, -1);
    int base_offset = floating ? 6 : 0;

    // Floating base body maps to index 5 (last virtual joint carries mass)
    if (floating) {
        body_to_rmodel[mj_floating_body] = 5;
    }

    // Map 1-DOF bodies
    struct ModelBody {
        int mj_body;
        int mj_jnt;
        JointType type;
    };
    std::vector<ModelBody> model_bodies;

    for (int b = 1; b < nb; b++) {
        if (b == mj_floating_body) continue;
        int jnt_start = m->body_jntadr[b];
        int jnt_count = m->body_jntnum[b];
        for (int k = 0; k < jnt_count; k++) {
            int j = jnt_start + k;
            if (m->jnt_type[j] == mjJNT_HINGE || m->jnt_type[j] == mjJNT_SLIDE) {
                if (body_to_rmodel[b] >= 0) {
                    throw std::runtime_error(
                        "buildRobotModelFromMJCF: multiple joints on one body are unsupported");
                }
                int idx = base_offset + static_cast<int>(model_bodies.size());
                body_to_rmodel[b] = idx;
                model_bodies.push_back({b, j,
                    m->jnt_type[j] == mjJNT_HINGE
                        ? JointType::REVOLUTE : JointType::PRISMATIC});
            } else {
                throw std::runtime_error(
                    "buildRobotModelFromMJCF: only a floating-base free joint, hinge, "
                    "slide, or fixed body is supported");
            }
        }
        if (body_to_rmodel[b] < 0) {
            int idx = base_offset + static_cast<int>(model_bodies.size());
            body_to_rmodel[b] = idx;
            model_bodies.push_back({b, -1, JointType::FIXED});
        }
    }

    // ----------------------------------------------------------------
    // Step 4: Helper functions
    // ----------------------------------------------------------------
    // Find the nearest RobotModel ancestor (walk up body tree)
    auto findParent = [&](int mj_body) -> int {
        int p = m->body_parentid[mj_body];
        while (p > 0 && body_to_rmodel[p] < 0) {
            p = m->body_parentid[p];
        }
        if (p <= 0) return -1;
        return body_to_rmodel[p];
    };

    // Compose transform across skipped bodies up to nearest RM ancestor
    auto composeTransform = [&](int mj_body) -> Eigen::Matrix4d {
        Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
        int b = mj_body;
        while (true) {
            Eigen::Vector3d pos(m->body_pos[3*b+0], m->body_pos[3*b+1], m->body_pos[3*b+2]);
            Eigen::Quaterniond quat(m->body_quat[4*b+0], m->body_quat[4*b+1],
                                    m->body_quat[4*b+2], m->body_quat[4*b+3]);
            Eigen::Matrix4d T_local = Eigen::Matrix4d::Identity();
            T_local.block<3,3>(0,0) = quat.toRotationMatrix();
            T_local.block<3,1>(0,3) = pos;
            T = T_local * T;
            int p = m->body_parentid[b];
            if (p <= 0 || body_to_rmodel[p] >= 0) break;
            b = p;
        }
        return T;
    };

    // Get 3×3 inertia matrix from MuJoCo body data
    auto getBodyInertia = [&](int b) -> Eigen::Matrix3d {
        double ixx = m->body_inertia[3*b+0], iyy = m->body_inertia[3*b+1], izz = m->body_inertia[3*b+2];
        Eigen::Quaterniond iq(m->body_iquat[4*b+0], m->body_iquat[4*b+1],
                              m->body_iquat[4*b+2], m->body_iquat[4*b+3]);
        Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
        I(0,0)=ixx; I(1,1)=iyy; I(2,2)=izz;
        if (!iq.isApprox(Eigen::Quaterniond::Identity(), 1e-6)) {
            Eigen::Matrix3d R = iq.toRotationMatrix();
            I = R * I * R.transpose();
        }
        return I;
    };

    // ----------------------------------------------------------------
    // Step 5: Build RobotModel entries
    // ----------------------------------------------------------------
    RobotModel robot;
    robot.floatingBase = floating;
    int total = base_offset + static_cast<int>(model_bodies.size());
    robot.joints.reserve(total);

    // Helper to add a joint spec
    auto addJoint = [&](JointSpec js) {
        robot.joints.push_back(std::move(js));
    };

    if (floating) {
        // ----- Virtual joints (0-5) -----
        // Joint 0: prismatic X
        { JointSpec js;
          js.parent = -1; js.parentToJoint = Eigen::Matrix4d::Identity();
          js.jointAxis = Eigen::Vector3d::UnitX(); js.type = JointType::PRISMATIC;
          js.mass = 0; js.com = Eigen::Vector3d::Zero(); js.inertia = Eigen::Matrix3d::Zero();
          js.name = "virtual_px"; addJoint(std::move(js)); }
        // Joint 1: prismatic Y
        { JointSpec js;
          js.parent = 0; js.parentToJoint = Eigen::Matrix4d::Identity();
          js.jointAxis = Eigen::Vector3d::UnitY(); js.type = JointType::PRISMATIC;
          js.mass = 0; js.com = Eigen::Vector3d::Zero(); js.inertia = Eigen::Matrix3d::Zero();
          js.name = "virtual_py"; addJoint(std::move(js)); }
        // Joint 2: prismatic Z
        { JointSpec js;
          js.parent = 1; js.parentToJoint = Eigen::Matrix4d::Identity();
          js.jointAxis = Eigen::Vector3d::UnitZ(); js.type = JointType::PRISMATIC;
          js.mass = 0; js.com = Eigen::Vector3d::Zero(); js.inertia = Eigen::Matrix3d::Zero();
          js.name = "virtual_pz"; addJoint(std::move(js)); }
        // Joint 3: revolute X
        { JointSpec js;
          js.parent = 2; js.parentToJoint = Eigen::Matrix4d::Identity();
          js.jointAxis = Eigen::Vector3d::UnitX(); js.type = JointType::REVOLUTE;
          js.mass = 0; js.com = Eigen::Vector3d::Zero(); js.inertia = Eigen::Matrix3d::Zero();
          js.name = "virtual_rx"; addJoint(std::move(js)); }
        // Joint 4: revolute Y
        { JointSpec js;
          js.parent = 3; js.parentToJoint = Eigen::Matrix4d::Identity();
          js.jointAxis = Eigen::Vector3d::UnitY(); js.type = JointType::REVOLUTE;
          js.mass = 0; js.com = Eigen::Vector3d::Zero(); js.inertia = Eigen::Matrix3d::Zero();
          js.name = "virtual_ry"; addJoint(std::move(js)); }
        // Joint 5: revolute Z (carries base body mass and placement)
        int fb = mj_floating_body;
        Eigen::Matrix4d T_base = composeTransform(fb);
        Eigen::Vector3d com_base_body(m->body_ipos[3*fb+0], m->body_ipos[3*fb+1], m->body_ipos[3*fb+2]);
        { JointSpec js;
          js.parent = 4; js.parentToJoint = T_base;
          js.jointAxis = Eigen::Vector3d::UnitZ(); js.type = JointType::REVOLUTE;
          js.mass = m->body_mass[fb]; js.com = com_base_body; js.inertia = getBodyInertia(fb);
          js.name = "virtual_rz"; addJoint(std::move(js)); }
    }

    // ----- Articulated and fixed bodies -----
    for (const auto& body : model_bodies) {
        int b = body.mj_body, j = body.mj_jnt;
        JointSpec js;
        js.type = body.type;
        if (j >= 0) {
            js.jointAxis = Eigen::Vector3d(
                m->jnt_axis[3*j+0], m->jnt_axis[3*j+1], m->jnt_axis[3*j+2]);
        }
        js.parent = findParent(b);
        js.parentToJoint = composeTransform(b);
        js.mass = m->body_mass[b];
        js.com = Eigen::Vector3d(m->body_ipos[3*b+0], m->body_ipos[3*b+1], m->body_ipos[3*b+2]);
        js.inertia = getBodyInertia(b);
        js.name = m->names + m->name_bodyadr[b];
        addJoint(std::move(js));
    }

    return robot;
}

} // namespace test_models
