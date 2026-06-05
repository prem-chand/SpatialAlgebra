/**
 * @file model_factory.cpp
 * @brief Implementation of ModelFactory class for constructing solver objects
 * @details Implements serial chain, branching, and per-link-config overloads
 *          for both ForwardDynamics and InverseDynamics solver objects.
 *          Follows the link construction pattern from examples/dynamics.cpp.
 */

#include "model_factory.h"
#include "PluckerTransform.h"
#include "RigidBodyInertia.h"
#include "LowerTriangular.h"

namespace SpatialAlgebra::Bench {

using lt = LowerTriangular;
using Vector3d = Eigen::Vector3d;

ForwardDynamics ModelFactory::createFD(int nDOF, const JointConfig& cfg) {
    ForwardDynamics fd;
    fd.links.reserve(nDOF);

    for (int i = 0; i < nDOF; ++i) {
        Link link;
        link.parent = i - 1;  // -1 for base (index 0), i-1 for children

        // Transform: identity rotation + configurable translation (D-02)
        link.X = PluckerTransform(
            Rotation(Eigen::Matrix3d::Identity()),
            cfg.translation
        );

        // Inertia: configurable mass, COM, isotropic unit inertia (D-04)
        link.I = RigidBodyInertia(cfg.mass, cfg.com, lt::Identity(3));

        // Joint axis: configurable screw axis (D-02)
        link.S = cfg.axis;

        // State initialized to zero — benchmark sets q/qdot via RandomState
        link.q = 0.0;
        link.qdot = 0.0;

        // External forces initialized to zero
        link.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());

        fd.links.push_back(std::move(link));
    }

    return fd;
}

InverseDynamics ModelFactory::createID(int nDOF, const JointConfig& cfg) {
    InverseDynamics id;
    id.links.reserve(nDOF);

    for (int i = 0; i < nDOF; ++i) {
        InverseDynamicsLink link;
        link.parent = i - 1;  // -1 for base, i-1 for children

        // Transform: identity rotation + configurable translation (D-02)
        link.X = PluckerTransform(
            Rotation(Eigen::Matrix3d::Identity()),
            cfg.translation
        );

        // Inertia: configurable mass, COM, isotropic unit inertia (D-04)
        link.I = RigidBodyInertia(cfg.mass, cfg.com, lt::Identity(3));

        // Joint axis: configurable screw axis (D-02)
        link.S = cfg.axis;

        // State initialized to zero — benchmark sets q/qdot/qddot via RandomState
        link.q = 0.0;
        link.qdot = 0.0;
        link.qddot = 0.0;

        id.links.push_back(std::move(link));
    }

    return id;
}

ForwardDynamics ModelFactory::createFD(const std::vector<JointConfig>& configs) {
    ForwardDynamics fd;
    const int nDOF = static_cast<int>(configs.size());
    fd.links.reserve(nDOF);

    for (int i = 0; i < nDOF; ++i) {
        const JointConfig& cfg = configs[i];
        Link link;
        link.parent = i - 1;

        link.X = PluckerTransform(
            Rotation(Eigen::Matrix3d::Identity()),
            cfg.translation
        );

        link.I = RigidBodyInertia(cfg.mass, cfg.com, lt::Identity(3));
        link.S = cfg.axis;

        link.q = 0.0;
        link.qdot = 0.0;
        link.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());

        fd.links.push_back(std::move(link));
    }

    return fd;
}

InverseDynamics ModelFactory::createID(const std::vector<JointConfig>& configs) {
    InverseDynamics id;
    const int nDOF = static_cast<int>(configs.size());
    id.links.reserve(nDOF);

    for (int i = 0; i < nDOF; ++i) {
        const JointConfig& cfg = configs[i];
        InverseDynamicsLink link;
        link.parent = i - 1;

        link.X = PluckerTransform(
            Rotation(Eigen::Matrix3d::Identity()),
            cfg.translation
        );

        link.I = RigidBodyInertia(cfg.mass, cfg.com, lt::Identity(3));
        link.S = cfg.axis;

        link.q = 0.0;
        link.qdot = 0.0;
        link.qddot = 0.0;

        id.links.push_back(std::move(link));
    }

    return id;
}

ForwardDynamics ModelFactory::createFDBranching(int nDOF, int branchPoint,
                                                 const JointConfig& cfg) {
    ForwardDynamics fd;
    fd.links.reserve(nDOF);

    for (int i = 0; i < nDOF; ++i) {
        Link link;
        link.parent = i - 1;  // Default serial parent

        // Override for side branch: link at branchPoint+2 has parent = branchPoint
        // (same parent as the main branch child at branchPoint+1)
        if (i == branchPoint + 2) {
            link.parent = branchPoint;
        }

        link.X = PluckerTransform(
            Rotation(Eigen::Matrix3d::Identity()),
            cfg.translation
        );

        link.I = RigidBodyInertia(cfg.mass, cfg.com, lt::Identity(3));
        link.S = cfg.axis;

        link.q = 0.0;
        link.qdot = 0.0;
        link.f = ForceVector(Vector3d::Zero(), Vector3d::Zero());

        fd.links.push_back(std::move(link));
    }

    return fd;
}

InverseDynamics ModelFactory::createIDBranching(int nDOF, int branchPoint,
                                                  const JointConfig& cfg) {
    InverseDynamics id;
    id.links.reserve(nDOF);

    for (int i = 0; i < nDOF; ++i) {
        InverseDynamicsLink link;
        link.parent = i - 1;  // Default serial parent

        // Override for side branch: link at branchPoint+2 has parent = branchPoint
        // (same parent as the main branch child at branchPoint+1)
        if (i == branchPoint + 2) {
            link.parent = branchPoint;
        }

        link.X = PluckerTransform(
            Rotation(Eigen::Matrix3d::Identity()),
            cfg.translation
        );

        link.I = RigidBodyInertia(cfg.mass, cfg.com, lt::Identity(3));
        link.S = cfg.axis;

        link.q = 0.0;
        link.qdot = 0.0;
        link.qddot = 0.0;

        id.links.push_back(std::move(link));
    }

    return id;
}

}  // namespace SpatialAlgebra::Bench
