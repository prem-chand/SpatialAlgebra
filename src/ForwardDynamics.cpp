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
                
                // Bias acceleration includes gravity for base
                // c₀ = -g (Featherstone D-07: base link bias acceleration = -gravity)
                links[i].c = MotionVector(Vector3d::Zero(), -this->gravity);
            }
            else
            {
                // Propagate velocity from parent (X transforms parent→child)
                // vᵢ = Xᵢ·v_parent + Sᵢ·q̇ᵢ
                MotionVector vParent = links[parent].v;
                links[i].v = links[i].X.transformMotion(vParent) + 
                             links[i].S * links[i].qdot;
                
                // Compute bias acceleration (Coriolis/centrifugal terms)
                // cᵢ = Xᵢ·c_parent + vᵢ × Sᵢ · q̇ᵢ (Featherstone Algorithm 7.3)
                MotionVector cParent = links[parent].c;
                links[i].c = links[i].X.transformMotion(cParent) + 
                             cross(links[i].v, links[i].S) * links[i].qdot;
            }
        }
    }

void ForwardDynamics::inwardPass(const Eigen::VectorXd& tau)
{
    // Phase 1: Initialize Ia and pa for all links from rigid body inertia
    for (int i = 0; i < static_cast<int>(links.size()); i++)
    {
        double mass = links[i].I.getMass();
        Vector3d com = links[i].I.getCom();
        links[i].Ia = ArticulatedBodyInertia(
            links[i].I.getInertiaMatrixLT(),
            skew(com) * mass,
            lt::Identity(3) * mass
        );

        MotionVector IaV = links[i].Ia.apply(links[i].v);
        links[i].pa = links[i].Ia.apply(links[i].c) + cross(links[i].v, IaV);

        links[i].pa = ForceVector(
            links[i].pa.getAngular() + links[i].f.getAngular(),
            links[i].pa.getLinear() + links[i].f.getLinear()
        );
    }
    
    // Phase 2: Accumulate child contributions from tip to base
    // This must NOT re-initialize Ia/pa so child contributions are preserved
    for (int i = static_cast<int>(links.size()) - 1; i >= 0; i--)
    {
        int parent = links[i].parent;
        if (parent != -1)
        {
            // Transform child ABI from child to parent frame: X^{-1} * Ia * X^{-T}
            ArticulatedBodyInertia IaTransformed = links[i].X.invtformABI(links[i].Ia);
            links[parent].Ia = links[parent].Ia + IaTransformed;

            // Transform child bias force from child to parent: X^{-T} * pa
            ForceVector paTransformed = links[i].X.inverseTransformForce(links[i].pa);
            links[parent].pa = ForceVector(
                links[parent].pa.getAngular() + paTransformed.getAngular(),
                links[parent].pa.getLinear() + paTransformed.getLinear()
            );
        }
    }
    
    // Phase 3: Solve for joint accelerations
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

void ForwardDynamics::computeAccelerations(const Eigen::VectorXd& tau, const Vector3d& gravity)
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
    
    // Store gravity for use in outward pass
    this->gravity = gravity;
    
    // Phase 1: Outward pass - propagate velocities, compute bias accelerations
    outwardPass();
    
    // Phase 2: Inward pass - accumulate articulated inertias, solve for accelerations
    inwardPass(tau);
}
}
