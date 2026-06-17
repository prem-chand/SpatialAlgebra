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
            const Eigen::Matrix3d& R = static_cast<const Eigen::Matrix3d &>(X.getRotation());
            Eigen::Matrix3d r_skew = skew(X.getTranslation());

            Eigen::Matrix3d I_full = Ia.getInertia().getSymmetricMatrix();
            Eigen::Matrix3d M_full = Ia.getM().getSymmetricMatrix();
            Eigen::Matrix3d H = Ia.getH();

            Eigen::Matrix<double, 6, 6> Ia_6x6;
            Ia_6x6.block<3,3>(0,0) = I_full;
            Ia_6x6.block<3,3>(0,3) = H;
            Ia_6x6.block<3,3>(3,0) = H.transpose();
            Ia_6x6.block<3,3>(3,3) = M_full;

            // Construct X directly based on transformMotion: v_c = X * v_p
            Eigen::Matrix<double, 6, 6> X_6x6;
            X_6x6.block<3,3>(0,0) = R;
            X_6x6.block<3,3>(0,3) = Eigen::Matrix3d::Zero();
            X_6x6.block<3,3>(3,0) = -R * r_skew;
            X_6x6.block<3,3>(3,3) = R;

            // Transform to parent frame: X^T * Ia * X
            Eigen::Matrix<double, 6, 6> result = X_6x6.transpose() * Ia_6x6 * X_6x6;

            Eigen::Matrix3d I_new = 0.5 * (result.block<3,3>(0,0) + result.block<3,3>(0,0).transpose());
            Eigen::Matrix3d M_new = 0.5 * (result.block<3,3>(3,3) + result.block<3,3>(3,3).transpose());
            Eigen::Matrix3d H_new = result.block<3,3>(0,3);

            return ArticulatedBodyInertia(
                lt::fromFullMatrix(I_new), H_new, lt::fromFullMatrix(M_new));
        }

    }

    void ForwardDynamics::outwardPass()
    {
        for (int i = 0; i < static_cast<int>(links.size()); i++)
        {
            int parent = links[i].parent;

            if (parent == -1)
            {
                links[i].v = links[i].S * links[i].qdot;

                links[i].c = MotionVector(Vector3d::Zero(), -this->gravity);
            }
            else
            {
                MotionVector vParent = links[parent].v;
                links[i].v = links[i].X.transformMotion(vParent) +
                             links[i].S * links[i].qdot;

                MotionVector cParent = links[parent].c;
                links[i].c = links[i].X.transformMotion(cParent) +
                             cross(links[i].v, links[i].S) * links[i].qdot;
            }
        }
    }

    void ForwardDynamics::inwardPass(const Eigen::VectorXd& tau)
    {
        // Phase 1: Initialize Ia and pa from rigid body inertia
        for (int i = 0; i < static_cast<int>(links.size()); i++)
        {
            double mass = links[i].I.getMass();
            Vector3d com = links[i].I.getCom();
            links[i].Ia = ArticulatedBodyInertia(
                links[i].I.getInertiaMatrixLT(),
                skew(com) * mass,
                lt::Identity(3) * mass
            );

            ForceVector IaV = links[i].Ia.apply(links[i].v);
            links[i].pa = links[i].Ia.apply(links[i].c) + cross(links[i].v, IaV);

            links[i].pa = ForceVector(
                links[i].pa.getAngular() + links[i].f.getAngular(),
                links[i].pa.getLinear() + links[i].f.getLinear()
            );
        }

        // Phase 2: Tip-to-base sweep — condense Ia/pa and propagate to parent
        constexpr double EPSILON = 1e-10;
        std::vector<double> D_store(links.size());
        std::vector<ArticulatedBodyInertia> Ia_unc(links.size());

        for (int i = static_cast<int>(links.size()) - 1; i >= 0; i--)
        {
            Ia_unc[i] = links[i].Ia;

            ForceVector IaS = links[i].Ia.apply(links[i].S);
            double D = dot(links[i].S, IaS);
            D_store[i] = D;

            if (std::abs(D) < EPSILON)
            {
                throw std::runtime_error(
                    "ForwardDynamics::inwardPass: Near-zero inertia at joint "
                    + std::to_string(i) + " (D=" + std::to_string(D) + ")"
                );
            }

            double u = tau[i] - dot(links[i].S, links[i].pa);
            links[i].qddot = u / D;

            int parent = links[i].parent;
            if (parent != -1)
            {
                double invD = 1.0 / D;

                Vector3d t = IaS.getAngular();
                Vector3d f = IaS.getLinear();

                lt inertiaCorr = lt::fromFullMatrix(t * t.transpose()) * invD;
                Eigen::Matrix3d HCorr = t * f.transpose() * invD;
                lt massCorr = lt::fromFullMatrix(f * f.transpose()) * invD;

                links[i].Ia = links[i].Ia
                    + (ArticulatedBodyInertia(inertiaCorr, HCorr, massCorr) * (-1.0));

                links[i].pa = links[i].pa + (IaS * links[i].qddot);

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

        // Phase 3: Forward pass — correct qddot for parent acceleration,
        //           using ENHANCED correction per Featherstone Algorithm 7.3
        std::vector<MotionVector> a(links.size());
        for (int i = 0; i < static_cast<int>(links.size()); i++)
        {
            int parent = links[i].parent;

            MotionVector aParentInChild = (parent == -1)
                ? MotionVector(Vector3d::Zero(), Vector3d::Zero())
                : links[i].X.transformMotion(a[parent]);

            // NIEMI: per Featherstone, the correction should account for
            // parent acceleration but NOT re-subtract Ia*c (which is
            // already included in pa via Phase 1 initialization).
            ForceVector accelBias = Ia_unc[i].apply(aParentInChild);
            double correction = dot(links[i].S, accelBias);
            links[i].qddot = links[i].qddot - (correction / D_store[i]);

            MotionVector a_prime = aParentInChild + links[i].c;
            a[i] = a_prime + links[i].S * links[i].qddot;
        }
    }

    void ForwardDynamics::computeAccelerations(const Eigen::VectorXd& tau, const Vector3d& gravity)
    {
        if (tau.size() != static_cast<int>(links.size()))
        {
            throw std::invalid_argument(
                "ForwardDynamics::computeAccelerations: tau size ("
                + std::to_string(tau.size()) + ") does not match link count ("
                + std::to_string(links.size()) + ")"
            );
        }

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

        outwardPass();

        inwardPass(tau);
    }
}
