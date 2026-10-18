# Rigid body inertia

The inertia example shows a library user creating mass properties, applying them to motion, and transforming the inertia to another frame.

## Sub-features

- `inertia-create` constructs centered and offset-center-of-mass inertias.
- `inertia-apply` maps rotational and linear motion to force.
- `inertia-transform` transforms inertia properties between frames.

## How to get to it (user POV)

- Build and run the `example_inertia` program from the repository root.

## Driving it with the shell

Preconditions:

- `build/examples/example_inertia` exists and is executable.
- Configure and build it after source changes with `cmake -S . -B build && cmake --build build --target example_inertia`.

- **Run inertia example.** Execute `./build/examples/example_inertia`. Require exit code `0`, the `SpatialAlgebra Inertia Example` heading, and the `Example Complete` marker.
- **Check inertia application.** Confirm the output shows force from rotational motion and force from linear acceleration.
- **Save proof.** Capture the exact command, output, and exit code in `.verification/evidence/<run-id>/example_inertia.txt`.

## Gotchas

- The example prints approximate decimal values. Compare mathematical results using the precision shown rather than exact full-line text.
- Do not treat the printed example as a substitute for testing changed code paths it does not exercise.
