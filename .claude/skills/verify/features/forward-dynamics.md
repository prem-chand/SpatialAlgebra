# Forward dynamics

The forward dynamics example shows a library user applying joint torques to a two-link model and reading the ABA joint accelerations.

## Sub-features

- `aba-torque-response` computes accelerations for the initial torque pair.
- `aba-opposing-torques` computes accelerations for opposing torques.
- `aba-zero-torque` confirms the stationary zero-torque case returns zero acceleration.

## How to get to it (user POV)

- Build and run the `example_dynamics` program from the repository root.

## Driving it with the shell

Preconditions:

- `build/examples/example_dynamics` exists and is executable.
- Configure and build it after source changes with `cmake -S . -B build && cmake --build build --target example_dynamics`.

- **Run the ABA example.** Execute `./build/examples/example_dynamics`. Require exit code `0`, the `SpatialAlgebra Forward Dynamics (ABA) Example` heading, and the `Example Complete` marker.
- **Check initial torque response.** Confirm the output reports `[1.0, 0.5] N·m` as the applied torques and prints both joint accelerations.
- **Check zero torque.** Confirm the third case reports `[0, 0] N·m` and both joint accelerations are `0`.
- **Save proof.** Capture the exact command, output, and exit code in `.verification/evidence/<run-id>/example_dynamics.txt`.

## Gotchas

- The example fixes both joint positions at zero. It does not prove q-dependent kinematics or every ABA topology.
- The example's printed coupling description is explanatory text. Validate numerical output for changed dynamics behavior against an independent reference when the feature map for that change requires it.
