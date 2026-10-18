# Spatial vectors

The vector examples show a library user creating twists and wrenches, combining them, and evaluating spatial products.

## Sub-features

- `vector-motion` creates and combines motion vectors.
- `vector-force` creates and combines force vectors.
- `vector-products` evaluates motion-force and motion-motion products.

## How to get to it (user POV)

- Build and run the `example_vectors` program from the repository root.

## Driving it with the shell

Preconditions:

- `build/examples/example_vectors` exists and is executable.
- Configure and build it after source changes with `cmake -S . -B build && cmake --build build --target example_vectors`.

- **Run vector example.** Execute `./build/examples/example_vectors`. Require exit code `0`, the `SpatialAlgebra Basic Vectors Example` heading, and the `Example Complete` marker.
- **Check the shown operations.** Confirm the transcript includes motion-vector and force-vector output and a `Power =` result.
- **Save proof.** Capture the exact command, output, and exit code in `.verification/evidence/<run-id>/example_vectors.txt`.

## Gotchas

- Output is a demonstration, not a machine-readable API contract. Assert the headings and operation results, not every whitespace character.
- Do not replace this path with unit tests when proving the public example workflow.
