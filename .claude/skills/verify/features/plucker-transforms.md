# Plucker transforms

The transform example shows a library user rotating and translating spatial motion and force values, then applying the inverse transform.

## Sub-features

- `transform-motion` transforms a motion vector between frames.
- `transform-force` transforms a force vector between frames.
- `transform-round-trip` applies a transform and its inverse.

## How to get to it (user POV)

- Build and run the `example_transforms` program from the repository root.

## Driving it with the shell

Preconditions:

- `build/examples/example_transforms` exists and is executable.
- Configure and build it after source changes with `cmake -S . -B build && cmake --build build --target example_transforms`.

- **Run transform example.** Execute `./build/examples/example_transforms`. Require exit code `0`, the `SpatialAlgebra Plücker Transforms Example` heading, and the `Example Complete` marker.
- **Check the round trip.** Confirm the transcript includes `Transform back to source frame` and the displayed value matches the original source motion within the example's printed precision.
- **Save proof.** Capture the exact command, output, and exit code in `.verification/evidence/<run-id>/example_transforms.txt`.

## Gotchas

- Check the inverse round-trip output as well as process exit status.
- Do not rely on exact spacing in Eigen's printed matrices.
