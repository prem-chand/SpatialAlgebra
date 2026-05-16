/**
 * @file InverseDynamics.cpp
 * @brief Implementation of Recursive Newton-Euler Algorithm (RNEA) for inverse dynamics
 * @details This file implements the core RNEA recursive algorithms:
 *          - outwardPass(): velocity and acceleration propagation
 *          - inwardPass(): force propagation and torque computation
 *          - computeTorques(): main solver combining both passes
 * 
 *          The implementation follows Featherstone (2008) Algorithm 7.1
 *          with O(n) computational complexity for serial chains.
 */

#include "InverseDynamics.h"
#include <stdexcept>
#include <cmath>

namespace SpatialAlgebra
{
    void InverseDynamics::outwardPass()
    {
        // Iterate from base (0) to tip (n-1)
        for (int i = 0; i < static_cast<int>(links.size()); i++)
        {
            int parent = links[i].parent;
            
            if (parent == -1)
            {
                // Base link: velocity and acceleration from joint only
                // v₀ = S₀·q̇₀
                links[i].v = links[i].S * links[i].qdot;
                
                // a₀ = S₀·q̈₀ (no parent acceleration or velocity products)
                links[i].a = links[i].S * links[i].qddot;
            }
            else
            {
                // Propagate velocity from parent
                // vᵢ = Xᵢ⁻¹·v_parent + Sᵢ·q̇ᵢ
                MotionVector vParent = links[parent].v;
                links[i].v = links[i].X.inverseTransformMotion(vParent) + 
                             links[i].S * links[i].qdot;
                
                // Propagate acceleration from parent
                // aᵢ = Xᵢ⁻¹·a_parent + Sᵢ·q̈ᵢ + vᵢ × Sᵢ·q̇ᵢ
                // The cross product term is the Coriolis/centrifugal acceleration
                MotionVector aParent = links[parent].a;
                MotionVector coriolis = cross(links[i].v, links[i].S) * links[i].qdot;
                links[i].a = links[i].X.inverseTransformMotion(aParent) + 
                             links[i].S * links[i].qddot + 
                             coriolis;
            }
        }
    }

    Eigen::VectorXd InverseDynamics::inwardPass()
    {
        // Vector to store joint torques
        Eigen::VectorXd tau = Eigen::VectorXd::Zero(links.size());
        
        // Vector to store spatial forces (initialized to zero)
        std::vector<ForceVector> f(links.size(), 
            ForceVector(Vector3d::Zero(), Vector3d::Zero()));
        
        // Iterate from tip (n-1) to base (0)
        // This ensures children are processed before their parent
        for (int i = static_cast<int>(links.size()) - 1; i >= 0; i--)
        {
            // Compute spatial force at link i
            // fᵢ = Iᵢ·aᵢ + vᵢ × Iᵢ·vᵢ
            // First term: inertial force from acceleration
            // Second term: Coriolis/centrifugal force from velocity
            ForceVector inertialForce = links[i].I.apply(links[i].a);
            ForceVector coriolisForce = cross(links[i].v, links[i].I.apply(links[i].v));
            f[i] = ForceVector(
                inertialForce.getAngular() + coriolisForce.getAngular(),
                inertialForce.getLinear() + coriolisForce.getLinear()
            );
            
            // Add transformed forces from children
            // Find all children of this link and add their transformed forces
            for (int j = 0; j < static_cast<int>(links.size()); j++)
            {
                if (links[j].parent == i)
                {
                    // Child j: transform its force to parent frame and add
                    // fᵢ += Xⱼ·fⱼ (transform force from child to parent)
                    ForceVector fChildTransformed = links[j].X.transformForce(f[j]);
                    f[i] = ForceVector(
                        f[i].getAngular() + fChildTransformed.getAngular(),
                        f[i].getLinear() + fChildTransformed.getLinear()
                    );
                }
            }
            
            // Compute joint torque by projecting force onto joint axis
            // τᵢ = fᵢ·Sᵢ (dot product in 6D)
            tau[i] = f[i].getAngular().dot(links[i].S.getAngular()) + 
                     f[i].getLinear().dot(links[i].S.getLinear());
        }
        
        return tau;
    }

    Eigen::VectorXd InverseDynamics::computeTorques(const Eigen::VectorXd& qddot)
    {
        // Validate input dimensions
        if (qddot.size() != static_cast<int>(links.size()))
        {
            throw std::invalid_argument(
                "InverseDynamics::computeTorques: qddot size (" + 
                std::to_string(qddot.size()) + ") does not match link count (" + 
                std::to_string(links.size()) + ")"
            );
        }
        
        // Validate for NaN/Inf in input (threat mitigation T-09-01)
        for (int i = 0; i < qddot.size(); i++)
        {
            if (std::isnan(qddot[i]) || std::isinf(qddot[i]))
            {
                throw std::invalid_argument(
                    "InverseDynamics::computeTorques: Invalid acceleration at joint " + 
                    std::to_string(i)
                );
            }
        }
        
        // Set joint accelerations from input
        for (int i = 0; i < static_cast<int>(links.size()); i++)
        {
            links[i].qddot = qddot[i];
        }
        
        // Phase 1: Outward pass - propagate velocities and accelerations
        outwardPass();
        
        // Phase 2: Inward pass - propagate forces, compute joint torques
        Eigen::VectorXd tau = inwardPass();
        
        return tau;
    }
}
