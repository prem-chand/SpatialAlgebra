# Phase 20: Test Model Library — Discussion Log

**Gathered:** 2026-06-17

## Discussion Summary

### Adapter Interface Design
- **Q:** Method scope? → **A:** Full solver API (9+ methods)
- **Q:** State management? → **A:** Hybrid — model holds reference, adapter caches internal
- **Q:** Construction pattern? → **A:** Builder pattern — adapter.build(Model) → ready solver
- **Q:** Concrete methods? → **A:** Core + FK + mass + Jacobians (setState, computeTorques, computeAccelerations, forwardKinematics, getJointTransform, computeMassMatrix, computeGravityTorques, computeJointSpaceJacobian, getLinkCOM, getDOF)
- **Q:** C++ interface pattern? → **A:** Abstract base class with virtual methods
- **Q:** Eigen types? → **A:** Dynamic Eigen (VectorXd, MatrixXd)

### Model Definition Format
- **Q:** Format? → **A:** Pure data POD structs (no methods, Eigen types only)
- **Q:** Directory? → **A:** tests/test-models/

### CMake Strategy
- **Q:** Integration? → **A:** Separate INTERFACE library target

### Naming
- **Q:** Class names? → **A:** RobotSolver (adapter ABC) + RobotModel (model POD)

### Key Decisions

| ID | Category | Decision |
|----|----------|----------|
| D-01 | Adapter | Abstract base class with virtual methods |
| D-02 | Adapter | Full API: setState, computeTorques, computeAccelerations, forwardKinematics, getJointTransform, computeMassMatrix, computeGravityTorques, computeJointSpaceJacobian, getLinkCOM, getDOF |
| D-03 | Adapter | Dynamic Eigen types (VectorXd, MatrixXd) |
| D-04 | Adapter | Hybrid state: model reference + adapter cache |
| D-05 | Adapter | Builder pattern: adapter.build(RobotModel) |
| D-06 | Naming | RobotSolver, RobotModel, JointSpec |
| D-07 | Model | Pure data POD structs, Eigen types only |
| D-08 | Model | JointSpec struct with parent, transform, axis, type, mass, com, inertia |
| D-09 | Model | JointType enum: REVOLUTE, PRISMATIC, FIXED |
| D-10 | Dir | tests/test-models/ with chains/ subdirectory |
| D-11 | Dir | One header per test domain in chains/ |
| D-12 | Dir | SA adapter is only file with SpatialAlgebra includes |
| D-13 | CMake | INTERFACE library test_models → Eigen3::Eigen |
| D-14 | CMake | SA adapter as separate sa_test_adapter library |
| D-15 | CMake | Existing tests NOT converted — smoke test only |
