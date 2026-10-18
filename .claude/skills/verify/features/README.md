# SpatialAlgebra verification map

Read this index before driving the library. The C++ examples are short-lived programs under `examples/`. Build the relevant executable from the repository root, run it, and retain its transcript under `.verification/evidence/<run-id>/`.

## Baseline preconditions

- Run commands from the repository root.
- Configure and build with the Launch commands in `../SKILL.md` after source or CMake changes.
- Run the Doctor command in `../SKILL.md` before driving an existing build.
- Use one build at a time in the shared `build/` directory.
- Give every evidence run a unique `RUN_ID`.

## Driving conventions

- Read each feature file and drive every entry point it lists when the change affects that feature.
- Capture the literal command, combined stdout and stderr, and exit code.
- Check the feature's stated output. A successful build alone is not behavior evidence.
- Examples use fixed inputs and have no persistent data or external service dependencies.
- Do not use `example_robot_3link` output as verification. It contains stale limitation claims.

## Proof and cleanup

- Keep transcripts in `.verification/evidence/` after each run.
- Each transcript must show the action and resulting output.
- Do not delete `build/` or evidence during cleanup.
- Examples terminate by themselves. Stop only a process started by the current run if interruption leaves it running.

## Features

- [Spatial vectors](./spatial-vectors.md) covers motion and force vector examples.
- [Plucker transforms](./plucker-transforms.md) covers coordinate transforms and inverse round trips.
- [Rigid body inertia](./rigid-body-inertia.md) covers inertia application and transformed properties.
- [Forward dynamics](./forward-dynamics.md) covers ABA joint accelerations.
