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
#include <stdexcept>
#include <cmath>

namespace SpatialAlgebra
{
    void ForwardDynamics::outwardPass()
    {
        // Iterate from base (0) to tip (n-1)
        for (int i = 0; i < static_cast<int>(links.size()); i++)
        {
            int parent = links[i].parent;
            
            if (parent == -1)
            {
                // Base link: velocity from joint only
                // v₀ = S₀·q̇₀
                links[i].v = links[i].S * links[i].qdot;
                
                // Bias acceleration is zero for base
                // c₀ = 0
                links[i].c = MotionVector(Vector3d::Zero(), Vector3d::Zero());
            }
            else
            {
                // Propagate velocity from parent
                // vᵢ = Xᵢ⁻¹·v_parent + Sᵢ·q̇ᵢ
                MotionVector vParent = links[parent].v;
                links[i].v = links[i].X.inverseTransformMotion(vParent) + 
                             links[i].S * links[i].qdot;
                
                // Compute bias acceleration (Coriolis/centrifugal terms)
                // cᵢ = vᵢ × Sᵢ · q̇ᵢ
                // Using spatial cross product: cross(v, S)
                links[i].c = cross(links[i].v, links[i].S) * links[i].qdot;
            }
        }
    }

    void ForwardDynamics::inwardPass(const Eigen::VectorXd& tau)
    {
        // Iterate from tip (n-1) to base (0)
        // This ensures children are processed before their parent
        for (int i = static_cast<int>(links.size()) - 1; i >= 0; i--)
        {
            // Initialize articulated inertia with rigid body inertia
            // Iₐᵢ = Iᵢ (before adding child contributions)
            // RigidBodyInertia converts to ArticulatedBodyInertia as:
            // Ia.M = m*I (mass times identity), Ia.H = skew(com), Ia.Inertia = I_LT
            links[i].Ia = ArticulatedBodyInertia(
                links[i].I.getInertiaMatrixLT(),  // Rotational inertia
                skew(links[i].I.getCom()),         // Coupling matrix from COM
                lt::Identity(3) * links[i].I.getMass()  // Mass matrix
            );
            
            // Initialize bias force
            // pₐᵢ = Iₐᵢ·cᵢ (bias force from acceleration)
            links[i].pa = links[i].Ia.apply(links[i].c);
            
            // Add external forces if present
            // pₐᵢ = pₐᵢ + fᵢ
            links[i].pa = ForceVector(
                links[i].pa.getAngular() + links[i].f.getAngular(),
                links[i].pa.getLinear() + links[i].f.getLinear()
            );
            
            // Propagate to parent
            int parent = links[i].parent;
            if (parent != -1)
            {
                // Transform articulated inertia to parent frame
                // Iₐ_parent += Xᵢ·Iₐᵢ (where Xᵢ transforms from parent to child)
                // Note: tformABI transforms inertia from child to parent coordinates
                ArticulatedBodyInertia IaTransformed = links[i].X.tformABI(links[i].Ia);
                
                // Add to parent's articulated inertia
                links[parent].Ia = links[parent].Ia + IaTransformed;
                
                // Transform bias force to parent frame
                // pₐ_parent += Xᵢ·pₐᵢ
                ForceVector paTransformed = links[i].X.transformForce(links[i].pa);
                links[parent].pa = ForceVector(
                    links[parent].pa.getAngular() + paTransformed.getAngular(),
                    links[parent].pa.getLinear() + paTransformed.getLinear()
                );
            }
        }
        
        // Solve for joint accelerations
        // q̈ᵢ = (τᵢ - Sᵢᵀ·pₐᵢ) / (Sᵢᵀ·Iₐᵢ·Sᵢ)
        constexpr double EPSILON = 1e-10;
        
        for (int i = static_cast<int>(links.size()) - 1; i >= 0; i--)
        {
            // Compute Iₐᵢ·Sᵢ
            MotionVector IaS = links[i].Ia.apply(links[i].S);
            
            // Compute denominator: Sᵢᵀ·Iₐᵢ·Sᵢ
            // Dot product in 6D: angular·angular + linear·linear
            double denom = links[i].S.getAngular().dot(IaS.getAngular()) + 
                          links[i].S.getLinear().dot(IaS.getLinear());
            
            // Check for singular configuration
            if (std::abs(denom) < EPSILON)
            {
                throw std::runtime_error(
                    "ForwardDynamics::inwardPass: Near-zero inertia at joint " + 
                    std::to_string(i) + " (denom=" + std::to_string(denom) + ")"
                );
            }
            
            // Compute numerator: τᵢ - Sᵢᵀ·pₐᵢ
            double numer = tau[i] - (
                links[i].S.getAngular().dot(links[i].pa.getAngular()) + 
                links[i].S.getLinear().dot(links[i].pa.getLinear())
            );
            
            // Solve for acceleration
            links[i].qddot = numer / denom;
        }
    }

    void ForwardDynamics::computeAccelerations(const Eigen::VectorXd& tau)
    {
        // Validate input dimensions
        if (tau.size() != static_cast<int>(links.size()))
        {
            throw std::invalid_argument(
                "ForwardDynamics::computeAccelerations: tau size (" + 
                std::to_string(tau.size()) + ") does not match link count (" + 
                std::to_string(links.size()) + ")"
            );
        }
        
        // Validate for NaN/Inf in input
        for (int i = 0; i < tau.size(); i++)
        {
            if (std::isnan(tau[i]) || std::isinf(tau[i]))
            {
                throw std::invalid_argument(
                    "ForwardDynamics::computeAccelerations: Invalid torque at joint " + 
                    std::to_string(i)
                );
            }
        }
        
        // Phase 1: Outward pass - propagate velocities, compute bias accelerations
        outwardPass();
        
        // Phase 2: Inward pass - accumulate articulated inertias, solve for accelerations
        inwardPass(tau);
    }
}
