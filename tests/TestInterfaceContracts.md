# Interface Contracts for Dynamics Tests

## InverseDynamics Interface

```cpp
namespace SpatialAlgebra {

struct Link {
  int parent;                     // Parent link index (-1 for base)
  PluckerTransform X;             // Transform from parent to this link
  RigidBodyInertia I;             // Rigid body inertia
  MotionVector S;                 // Joint motion axis (screw axis)
  
  double q;                       // Joint position
  double qdot;                    // Joint velocity
  double qddot;                   // Joint acceleration (input for RNEA)
  
  MotionVector v;                 // Spatial velocity (computed)
  MotionVector a;                 // Spatial acceleration (computed)
};

class InverseDynamics {
  std::vector<Link> links;        // Kinematic tree
  
  Eigen::VectorXd computeTorques(const Eigen::VectorXd& qddot);
  // Input: joint accelerations (n x 1)
  // Output: joint torques (n x 1)
  // Throws: std::invalid_argument if size mismatch or NaN/Inf
  
  void outwardPass();             // Propagate velocities/accelerations
  Eigen::VectorXd inwardPass();   // Propagate forces, compute torques
};

}
```

## ForwardDynamics Interface

```cpp
namespace SpatialAlgebra {

struct Link {
  int parent;                     // Parent link index (-1 for base)
  PluckerTransform X;             // Transform from parent to this link
  RigidBodyInertia I;             // Rigid body inertia
  MotionVector S;                 // Joint motion axis (screw axis)
  
  double q;                       // Joint position
  double qdot;                    // Joint velocity
  double qddot;                   // Joint acceleration (output for ABA)
  
  MotionVector v;                 // Spatial velocity (computed)
  MotionVector c;                 // Bias acceleration (computed)
  ForceVector f;                  // Spatial force (external forces)
  
  ArticulatedBodyInertia Ia;      // Articulated body inertia (computed)
  ForceVector pa;                 // Bias force (computed)
};

class ForwardDynamics {
  std::vector<Link> links;        // Kinematic tree
  
  void computeAccelerations(const Eigen::VectorXd& tau);
  // Input: joint torques (n x 1)
  // Output: sets links[i].qddot for all links
  // Throws: std::invalid_argument if size mismatch, std::runtime_error if singular
  
  void outwardPass();             // Propagate velocities, compute bias accelerations
  void inwardPass(const Eigen::VectorXd& tau);  // Accumulate inertias, solve
};

}
```

## Round-Trip Test Pattern

```cpp
// ABA → RNEA consistency
Eigen::VectorXd qddot_input = ...;  // Random accelerations
InverseDynamics id;
// Setup links with same q, qdot as ForwardDynamics
Eigen::VectorXd tau = id.computeTorques(qddot_input);

ForwardDynamics fd;
// Setup identical links
fd.computeAccelerations(tau);

// Verify: fd.links[i].qddot ≈ qddot_input[i]
```

```cpp
// RNEA → ABA consistency
Eigen::VectorXd tau_input = ...;  // Random torques
ForwardDynamics fd;
// Setup links
fd.computeAccelerations(tau_input);

InverseDynamics id;
// Setup identical links with qddot from fd
Eigen::VectorXd tau_output = id.computeTorques(qddot_from_fd);

// Verify: tau_output ≈ tau_input
```
