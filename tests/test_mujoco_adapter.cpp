#include <gtest/gtest.h>
#include "mujoco_adapter.h"
#include "sa_adapter.h"
#include "mjcf_parser.h"
#include <mujoco/mujoco.h>
#include <Eigen/Dense>
#include <filesystem>
#include <memory>
#include <random>
#include <vector>

using namespace Eigen;
using namespace test_models;

namespace {

const std::string MODELS_DIR = std::string(MODELS_ROOT) + "/menagerie";
const Vector3d GRAVITY(0, 0, -9.81);

std::string modelPath(const std::string& relative) {
    return MODELS_DIR + "/" + relative;
}

Matrix4d pose(const Vector3d& t, const Matrix3d& R) {
    Matrix4d T = Matrix4d::Identity();
    T.topLeftCorner<3, 3>() = R;
    T.topRightCorner<3, 1>() = t;
    return T;
}

JointSpec link(int parent, const Matrix4d& T, const Vector3d& axis, JointType type) {
    JointSpec js;
    js.parent = parent;
    js.parentToJoint = T;
    js.jointAxis = axis;
    js.type = type;
    js.mass = 1.5;
    js.com = Vector3d(0.05, 0.02, -0.1);
    js.inertia = Vector3d(0.02, 0.03, 0.04).asDiagonal();
    return js;
}

// A branching fixed-base chain with rotated frames, a fixed body, and a
// prismatic joint, so that each joint convention affects the result.
RobotModel rotatedChain() {
    Matrix3d R = AngleAxisd(0.7, Vector3d(1, 2, 3).normalized()).toRotationMatrix();
    Vector3d t(0.1, -0.2, 0.3);
    RobotModel model;
    model.joints.push_back(link(-1, pose(t, R), Vector3d::UnitY(), JointType::REVOLUTE));
    model.joints.push_back(link(0, pose(t, R.transpose()), Vector3d::UnitX(), JointType::FIXED));
    model.joints.push_back(link(1, pose(-t, R), Vector3d::UnitZ(), JointType::REVOLUTE));
    model.joints.push_back(link(0, pose(t, R), Vector3d::UnitZ(), JointType::PRISMATIC));
    return model;
}

double relErr(const MatrixXd& a, const MatrixXd& b) {
    return (a - b).norm() / b.norm();
}

void expectSolversAgree(const RobotModel& model, double tol) {
    auto sa = SpatialAlgebraAdapter::build(model);
    auto mj = MujocoAdapter::build(model);
    const int n = model.getDOF();
    ASSERT_EQ(sa->getDOF(), n);
    ASSERT_EQ(mj->getDOF(), n);

    std::mt19937 rng(7);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (int trial = 0; trial < 10; ++trial) {
        VectorXd q(n), qdot(n), qddot(n);
        for (int i = 0; i < n; ++i) {
            q(i) = dist(rng);
            qdot(i) = dist(rng);
            qddot(i) = dist(rng);
        }

        sa->setState(q, qdot);
        mj->setState(q, qdot);
        EXPECT_LT(relErr(sa->computeMassMatrix(), mj->computeMassMatrix()), tol);
        EXPECT_LT(relErr(sa->computeGravityTorques(GRAVITY),
                         mj->computeGravityTorques(GRAVITY)), tol);

        sa->setState(q, qdot);
        mj->setState(q, qdot);
        VectorXd tau = mj->computeTorques(qddot, GRAVITY);
        EXPECT_LT(relErr(sa->computeTorques(qddot, GRAVITY), tau), tol);

        sa->setState(q, qdot);
        mj->setState(q, qdot);
        EXPECT_LT(relErr(sa->computeAccelerations(tau, GRAVITY), qddot), tol);
        EXPECT_LT(relErr(mj->computeAccelerations(tau, GRAVITY), qddot), tol);

        sa->setState(q, qdot);
        mj->setState(q, qdot);
        mj->forwardKinematics();
        for (int i = 0; i < n; ++i) {
            EXPECT_LT((sa->getJointTransform(i) - mj->getJointTransform(i)).norm(), tol);
            EXPECT_LT((sa->getLinkCOM(i) - mj->getLinkCOM(i)).norm(), tol);
            EXPECT_LT((sa->computeJointSpaceJacobian(i) -
                       mj->computeJointSpaceJacobian(i)).norm(), tol);
        }
    }
}

} // namespace

TEST(MujocoAdapterTest, RotatedChainMatchesSpatialAlgebra) {
    expectSolversAgree(rotatedChain(), 1e-10);
}

TEST(MujocoAdapterTest, Go1MatchesSpatialAlgebra) {
    std::string path = modelPath("unitree_go1/go1.xml");
    if (!std::filesystem::exists(path)) GTEST_SKIP() << "Run scripts/fetch_models.sh";
    expectSolversAgree(buildRobotModelFromMJCF(path), 1e-6);
}

// The parsed RobotModel must describe the same robot as the original MJCF.
// With the floating base at its reference pose and zero velocity, the
// actuated rows of the gravity torques and the mass matrix do not depend
// on how the six base coordinates are chosen.
TEST(MujocoAdapterTest, ParsedModelsMatchNativeMujoco) {
    for (const char* relative : {"unitree_go1/go1.xml", "unitree_go2/go2.xml",
                                 "unitree_h1/h1.xml", "unitree_g1/g1.xml"}) {
        SCOPED_TRACE(relative);
        std::string path = modelPath(relative);
        if (!std::filesystem::exists(path)) GTEST_SKIP() << "Run scripts/fetch_models.sh";

        char error[1000] = {};
        std::unique_ptr<mjModel, decltype(&mj_deleteModel)> m(
            mj_loadXML(path.c_str(), nullptr, error, sizeof(error)), mj_deleteModel);
        ASSERT_TRUE(m) << error;
        // RobotModel has no armature, so remove it from the native model.
        for (int i = 0; i < m->nv; ++i) m->dof_armature[i] = 0.0;
        std::unique_ptr<mjData, decltype(&mj_deleteData)> d(
            mj_makeData(m.get()), mj_deleteData);

        RobotModel model = buildRobotModelFromMJCF(path);
        ASSERT_TRUE(model.floatingBase);
        auto parsed = MujocoAdapter::build(model);
        const int n = model.getDOF();
        const int nact = m->nv - 6;
        ASSERT_EQ(n - 6, nact);

        std::mt19937 rng(3);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        VectorXd q = VectorXd::Zero(n);
        for (int i = 0; i < nact; ++i) {
            q(6 + i) = dist(rng);
            d->qpos[7 + i] = q(6 + i);
        }
        mj_forward(m.get(), d.get());

        std::vector<mjtNum> g(m->nv);
        mj_rne(m.get(), d.get(), 0, g.data());
        std::vector<mjtNum> M(m->nv * m->nv);
        mj_fullM(m.get(), M.data(), d->qM);
        VectorXd gNative(nact);
        MatrixXd hNative(nact, nact);
        for (int i = 0; i < nact; ++i) {
            gNative(i) = g[6 + i];
            for (int j = 0; j < nact; ++j) hNative(i, j) = M[(6 + i) * m->nv + 6 + j];
        }

        parsed->setState(q, VectorXd::Zero(n));
        EXPECT_LT(relErr(parsed->computeGravityTorques(GRAVITY).tail(nact), gNative), 1e-10);
        EXPECT_LT(relErr(parsed->computeMassMatrix().bottomRightCorner(nact, nact), hNative), 1e-6);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
