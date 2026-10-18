/**
 * @file ForwardDynamics.cpp
 * @brief Implementation of Articulated Body Algorithm (ABA) for forward dynamics
 * @details This file implements the core ABA recursive algorithms:
 *          - outwardPass(): velocity propagation and bias acceleration
 *          - inwardPass(): articulated inertia accumulation and force propagation
 *          - computeAccelerations(): main solver combining both passes
 *
 *          The implementation follows Featherstone (2008) Algorithm 7.3
 *          with O(n) computational complexity for serial chains.
 */

#include "ForwardDynamics.h"
#include "LowerTriangular.h"
#include "SpatialUtils.h"
#include <stdexcept>

namespace SpatialAlgebra
{
    namespace {

        ArticulatedBodyInertia transformInertiaToParent(const PluckerTransform &X,
                                                        const ArticulatedBodyInertia &Ia) {
            // I_parent = X^T · Ia · X  using block 3×3 formulas (avoids 6×6 temp)
            //
            // X = [R, 0; -R·skew(t), R]
            // X^T = [R^T, skew(t)·R^T; 0, R^T]
            // Ia = [I, H; H^T, M]
            //
            // Block result:
            //   I'_00 = R^T·I·R - R^T·H·R·r̂ + r̂·R^T·H^T·R - r̂·R^T·M·R·r̂
            //   I'_01 = R^T·H·R + r̂·R^T·M·R
            //   I'_11 = R^T·M·R
            // where r̂ = skew(t)

            const Eigen::Matrix3d& R = static_cast<const Eigen::Matrix3d &>(X.getRotation());
            Eigen::Matrix3d r_skew = skew(X.getTranslation());
            Eigen::Matrix3d Rt = R.transpose();

            Eigen::Matrix3d I_sym = Ia.getInertia().getSymmetricMatrix();
            Eigen::Matrix3d M_sym = Ia.getM().getSymmetricMatrix();
            Eigen::Matrix3d H = Ia.getH();

            // Precompute rotated blocks
            Eigen::Matrix3d R_I_Rt = Rt * I_sym * R;      // R^T·I·R
            Eigen::Matrix3d R_H_Rt = Rt * H * R;           // R^T·H·R
            Eigen::Matrix3d R_M_Rt = Rt * M_sym * R;       // R^T·M·R

            // I'_01 = R^T·H·R + r̂·R^T·M·R
            Eigen::Matrix3d H_new = R_H_Rt + r_skew * R_M_Rt;

            // I'_11 = R^T·M·R
            Eigen::Matrix3d M_new = R_M_Rt;

            // I'_00 = R^T·I·R - R^T·H·R·r̂ + r̂·R^T·H^T·R - r̂·R^T·M·R·r̂
            Eigen::Matrix3d R_H_Rt_skew = R_H_Rt * r_skew;          // R^T·H·R·r̂
            Eigen::Matrix3d skew_R_Ht_Rt = r_skew * R_H_Rt.transpose(); // r̂·R^T·H^T·R
            Eigen::Matrix3d skew_R_M_Rt = r_skew * R_M_Rt;           // r̂·R^T·M·R
            Eigen::Matrix3d skew_R_M_Rt_skew = skew_R_M_Rt * r_skew;  // r̂·R^T·M·R·r̂

            Eigen::Matrix3d I_new = R_I_Rt - R_H_Rt_skew + skew_R_Ht_Rt - skew_R_M_Rt_skew;

            // Symmetrize rotational and mass blocks
            I_new = 0.5 * (I_new + I_new.transpose());
            M_new = 0.5 * (M_new + M_new.transpose());

            return ArticulatedBodyInertia(
                lt::fromFullMatrix(I_new), H_new, lt::fromFullMatrix(M_new));
        }

    }

    void ForwardDynamics::updateKinematics()
    {
        for (auto& link : links)
        {
            // Skip FIXED joints — X is user-managed
            if (link.type == JointType::FIXED) continue;

            Rotation jointRot = Rotation(Eigen::Matrix3d::Identity());
            Vector3d jointTrans = Vector3d::Zero();

            switch (link.type)
            {
                case JointType::REVOLUTE:
                    jointRot = Rotation(Eigen::AngleAxisd(link.q, link.jointAxis));
                    break;
                case JointType::PRISMATIC:
                    jointTrans = link.jointAxis * link.q;
                    break;
                default:
                    break;
            }

            link.X = link.X0.multiply(PluckerTransform(jointRot, jointTrans));
        }
    }

    void ForwardDynamics::outwardPass()
    {
        for (int i = 0; i < static_cast<int>(links.size()); i++)
        {
            int parent = links[i].parent;
            MotionVector vJ = links[i].S * links[i].qdot;

            if (parent == -1)
            {
                links[i].v = vJ;
            }
            else
            {
                MotionVector vParent = links[parent].v;
                links[i].v = links[i].X.transformMotion(vParent) + vJ;
            }

            links[i].c = cross(links[i].v, vJ);
        }
    }

    void ForwardDynamics::inwardPass(const Eigen::VectorXd& tau)
    {
        // Phase 1: Initialize Ia and pa from rigid body inertia
        for (int i = 0; i < static_cast<int>(links.size()); i++)
        {
            double mass = links[i].I.getMass();
            Vector3d com = links[i].I.getCom();
            // Compute inertia about the joint origin (parallel axis theorem):
            // I_joint = I_cm + m*(|c|²*I₃ - c*c^T)
            Eigen::Matrix3d I_cm = links[i].I.getInertiaMatrixLT().getSymmetricMatrix();
            Eigen::Matrix3d I_joint = I_cm + mass * (
                com.dot(com) * Eigen::Matrix3d::Identity() - com * com.transpose()
            );
            links[i].Ia = ArticulatedBodyInertia(
                lt::fromFullMatrix(I_joint),
                skew(com) * mass,
                lt::Identity(3) * mass
            );

            ForceVector IaV = links[i].Ia.apply(links[i].v);
            links[i].pa = cross(links[i].v, IaV);

            links[i].pa = ForceVector(
                links[i].pa.getAngular() - links[i].f.getAngular(),
                links[i].pa.getLinear() - links[i].f.getLinear()
            );
        }

        // Phase 2: Tip-to-base sweep — condense Ia/pa and propagate to parent
        constexpr double EPSILON = 1e-10;
        struct JointSolveData {
            ForceVector U;
            double d;
            double u;
        };
        std::vector<JointSolveData> jointData(links.size());
        std::vector<int> tauIndex(links.size(), -1);
        int nextTau = 0;
        for (size_t i = 0; i < links.size(); ++i) {
            if (fixedLinksHaveDof || links[i].type != JointType::FIXED)
                tauIndex[i] = nextTau++;
        }

        for (int i = static_cast<int>(links.size()) - 1; i >= 0; i--)
        {
            int parent = links[i].parent;
            if (tauIndex[i] < 0) {
                links[i].qddot = 0.0;
                if (parent != -1) {
                    links[parent].Ia = links[parent].Ia
                        + transformInertiaToParent(links[i].X, links[i].Ia);
                    ForceVector paToParent = links[i].X.inverseTransformForce(
                        links[i].pa + links[i].Ia.apply(links[i].c));
                    links[parent].pa = links[parent].pa + paToParent;
                }
                continue;
            }

            ForceVector IaS = links[i].Ia.apply(links[i].S);
            double D = dot(links[i].S, IaS);

            if (std::abs(D) < EPSILON)
            {
                throw std::runtime_error(
                    "ForwardDynamics::inwardPass: Near-zero inertia at joint "
                    + std::to_string(i) + " (D=" + std::to_string(D) + ")"
                );
            }

            double u = tau[tauIndex[i]] - dot(links[i].S, links[i].pa);
            jointData[i] = {IaS, D, u};

            if (parent != -1)
            {
                double invD = 1.0 / D;

                Vector3d t = IaS.getAngular();
                Vector3d f = IaS.getLinear();

                lt inertiaCorr = lt::fromFullMatrix(t * t.transpose()) * invD;
                Eigen::Matrix3d HCorr = t * f.transpose() * invD;
                lt massCorr = lt::fromFullMatrix(f * f.transpose()) * invD;

                ArticulatedBodyInertia IaReduced = links[i].Ia
                    + (ArticulatedBodyInertia(inertiaCorr, HCorr, massCorr) * (-1.0));

                links[i].Ia = IaReduced;
                links[i].pa = links[i].pa + IaReduced.apply(links[i].c)
                            + (IaS * (u / D));

                ArticulatedBodyInertia IaTransformed =
                    transformInertiaToParent(links[i].X, links[i].Ia);
                links[parent].Ia = links[parent].Ia + IaTransformed;

                ForceVector paTransformed = links[i].X.inverseTransformForce(links[i].pa);
                links[parent].pa = ForceVector(
                    links[parent].pa.getAngular() + paTransformed.getAngular(),
                    links[parent].pa.getLinear() + paTransformed.getLinear()
                );
            }
        }

        // Phase 3: Forward pass — solve joint acceleration using the reduced
        // articulated quantities retained from the inward pass.
        std::vector<MotionVector> a(links.size());
        for (int i = 0; i < static_cast<int>(links.size()); i++)
        {
            int parent = links[i].parent;

            MotionVector aParentInChild = (parent == -1)
                ? links[i].X.transformMotion(
                    MotionVector(Vector3d::Zero(), -this->gravity))
                : links[i].X.transformMotion(a[parent]);

            if (tauIndex[i] < 0) {
                links[i].qddot = 0.0;
            } else {
                const JointSolveData& solve = jointData[i];
                links[i].qddot = (solve.u - dot(aParentInChild + links[i].c, solve.U)) / solve.d;
            }
            MotionVector a_prime = aParentInChild + links[i].c;
            a[i] = a_prime + links[i].S * links[i].qddot;
        }
    }

    void ForwardDynamics::computeAccelerations(const Eigen::VectorXd& tau, const Vector3d& gravity)
    {
        int dof = 0;
        for (const auto& link : links) {
            if (link.type != JointType::FIXED) ++dof;
        }
        if (tau.size() != dof && tau.size() != static_cast<int>(links.size()))
        {
            throw std::invalid_argument(
                "ForwardDynamics::computeAccelerations: tau size ("
                + std::to_string(tau.size()) + ") does not match DOF count ("
                + std::to_string(dof) + ")"
            );
        }
        fixedLinksHaveDof = tau.size() == static_cast<int>(links.size());

        for (int i = 0; i < tau.size(); i++)
        {
            if (std::isnan(tau[i]) || std::isinf(tau[i]))
            {
                throw std::invalid_argument(
                    "ForwardDynamics::computeAccelerations: Invalid torque at joint "
                    + std::to_string(i)
                );
            }
        }

        this->gravity = gravity;

        
        // Update kinematics: recompute X transforms from current q values
        updateKinematics();
        outwardPass();

        inwardPass(tau);
    }
}
