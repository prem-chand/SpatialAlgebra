---
name: verify
description: Verify SpatialAlgebra's C++ library through its compiled examples. Use when checking spatial vectors, Plucker transforms, inertias, or forward dynamics after a code change.
---

# Verify SpatialAlgebra

Use the compiled examples as the user-facing surface for this C++ library. Read `features/README.md` first, then follow the recipe for every entry point relevant to the change.

## Launch

From the repository root, configure and build the examples:

```sh
cmake -S . -B build
cmake --build build --target example_vectors example_transforms example_inertia example_dynamics example_robot_2link
```

The build is ready when each requested target reports `Built target`. No server, port, environment variable, authentication, or seed data is required. Each example is a short-lived executable. Do not build concurrently in the shared `build/` directory.

## Doctor

Run this read-only check before driving an existing build:

```sh
test -x build/examples/example_dynamics && test -f build/libSpatialAlgebra.a && printf 'SpatialAlgebra example and library are present\n'
```

This only checks that the executable and library exist. After any source or CMake change, rerun Launch before relying on them.

## Drive

Run the mapped example from the repository root. These are the compiled examples documented in `README.md` and wired in `examples/CMakeLists.txt`. Do not replace the example path with direct calls into internal classes or test-only endpoints.

## Evidence

Create a unique evidence directory for each run and preserve it after cleanup:

```sh
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)-$$"
EVIDENCE_DIR=".verification/evidence/$RUN_ID"
mkdir -p "$EVIDENCE_DIR"
```

Capture the exact command, combined output, and exit code. For example:

```sh
EVIDENCE_FILE="$EVIDENCE_DIR/example_dynamics.txt"
printf '%s\n' '$ ./build/examples/example_dynamics' > "$EVIDENCE_FILE"
./build/examples/example_dynamics >> "$EVIDENCE_FILE" 2>&1
STATUS=$?
printf '\nexit code: %s\n' "$STATUS" >> "$EVIDENCE_FILE"
cat "$EVIDENCE_FILE"
test "$STATUS" -eq 0
```

Check both the exit code and the feature-specific output listed in the map. The examples do not persist application data or make network calls, so there are no external side effects to inspect. Do not claim a feature is verified from build success alone.

## Cleanup

Examples exit on their own. Do not kill processes by name. If a run is interrupted, stop only the exact process started by that run. Keep `build/` because it may contain the developer's shared build state. Remove only scratch files created outside the evidence directory. Never delete `.verification/evidence/`; proof transcripts must survive cleanup.

## Features

- [Spatial vectors](features/spatial-vectors.md) covers twists, wrenches, and spatial products.
- [Plucker transforms](features/plucker-transforms.md) covers frame transforms and inverse round trips.
- [Rigid body inertia](features/rigid-body-inertia.md) covers inertia application and transformed mass properties.
- [Forward dynamics](features/forward-dynamics.md) covers accelerations from joint torques with ABA.
