#!/usr/bin/env python3
"""
Phase 23 — Pinocchio Cross-Validation Python Harness

Invokes the C++ Pinocchio comparison benchmark and generates a structured
JSON report with per-model, per-check results (PCP-04).

No native Python pinocchio bindings required — the script delegates to the
already-working C++ benchmark executable and parses its standardized output.

Usage:
    python3 robot_dynamics/pinocchio_comparison.py            # runs benchmark
    python3 robot_dynamics/pinocchio_comparison.py --json     # writes results.json
    python3 robot_dynamics/pinocchio_comparison.py --help     # full options
"""

import argparse
import json
import os
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
BUILD_DIR = PROJECT_ROOT / "build"
DEFAULT_BENCH = BUILD_DIR / "tests" / "test-models" / "bench_pinocchio_comparison"
DEFAULT_OUTPUT = PROJECT_ROOT / "robot_dynamics" / "pinocchio_comparison_results.json"


# ── Test model definitions (PCP-02) ────────────────────────────────────────
# Python equivalents of the C++ test model factories. Each entry describes
# the chain layout, joint types, and expected mathematical behavior.
# These mirror chains defined in tests/test-models/chains/{inverse,forward}_dynamics.h

TEST_MODELS = {
    "SingleLinkChain": {
        "dof": 1,
        "layout": "1-DOF Z-revolute at origin",
        "links": [
            {"parent": -1, "type": "REVOLUTE", "axis": [0, 0, 1],
             "mass": 1.0, "com": [0, 0, 0], "inertia": "I_3x3",
             "parentToJoint": "identity"}
        ],
        "reference": "Featherstone (2008) §5.1 — simplest non-trivial RNEA",
    },
    "TwoLinkSerialChain": {
        "dof": 2,
        "layout": "2-DOF serial, Z-revolute, 1m X offset",
        "links": [
            {"parent": -1, "type": "REVOLUTE", "axis": [0, 0, 1],
             "mass": 1.0, "com": [0, 0, 0], "parentToJoint": "identity"},
            {"parent": 0, "type": "REVOLUTE", "axis": [0, 0, 1],
             "mass": 1.0, "com": [0, 0, 0], "parentToJoint": "translate(1,0,0)"},
        ],
        "reference": "Two-link serial; M[0,0]=3 at q=0, M[1,1]=1",
    },
    "ThreeLinkSerialChain": {
        "dof": 3,
        "layout": "3-DOF serial, Z-revolute, 1m X spacing",
        "links": [
            {"parent": -1, "type": "REVOLUTE", "axis": [0, 0, 1],
             "mass": 1.0, "com": [0, 0, 0], "parentToJoint": "identity"},
            {"parent": 0, "type": "REVOLUTE", "axis": [0, 0, 1],
             "mass": 1.0, "com": [0, 0, 0], "parentToJoint": "translate(1,0,0)"},
            {"parent": 1, "type": "REVOLUTE", "axis": [0, 0, 1],
             "mass": 1.0, "com": [0, 0, 0], "parentToJoint": "translate(1,0,0)"},
        ],
        "reference": "3-link; M[0,0]=8 at q=0 (parallel axis accumulation)",
    },
    "SingleLinkXAxisCOMOffset": {
        "dof": 1,
        "layout": "1-DOF X-revolute with COM offset at (0, 0.5, 0)",
        "links": [
            {"parent": -1, "type": "REVOLUTE", "axis": [1, 0, 0],
             "mass": 1.0, "com": [0, 0.5, 0], "inertia": "I_3x3",
             "parentToJoint": "identity"}
        ],
        "reference": "Parallel axis: I_joint(XX) = 1 + 1*0.25 = 1.25",
    },
    "BranchingYConfiguration": {
        "dof": 3,
        "layout": "3-DOF tree: base → left (Z-rev) + right (Z-rev), each offset",
        "links": [
            {"parent": -1, "type": "REVOLUTE", "axis": [0, 0, 1],
             "mass": 1.0, "com": [0, 0, 0], "parentToJoint": "identity"},
            {"parent": 0, "type": "REVOLUTE", "axis": [0, 0, 1],
             "mass": 1.0, "com": [0, 0, 0], "parentToJoint": "translate(1,0,0)"},
            {"parent": 0, "type": "REVOLUTE", "axis": [0, 0, 1],
             "mass": 1.0, "com": [0, 0, 0], "parentToJoint": "translate(0,1,0)"},
        ],
        "reference": "Branching topology test",
    },
    "TwoLinkWithCOMOffset": {
        "dof": 2,
        "layout": "2-DOF Z-revolute, each link COM=(0.1, 0, 0)",
        "links": [
            {"parent": -1, "type": "REVOLUTE", "axis": [0, 0, 1],
             "mass": 1.0, "com": [0.1, 0, 0], "parentToJoint": "identity"},
            {"parent": 0, "type": "REVOLUTE", "axis": [0, 0, 1],
             "mass": 1.0, "com": [0.1, 0, 0],
             "parentToJoint": "translate(1,0,0)"},
        ],
        "reference": "2-link with parallel axis from COM drift",
    },
    "ThreeLinkSerialNonZeroCOM": {
        "dof": 3,
        "layout": "3-DOF Z-revolute, each link COM=(0.1, 0, 0), 1m spacing",
        "links": [
            {"parent": -1, "type": "REVOLUTE", "axis": [0, 0, 1],
             "mass": 1.0, "com": [0.1, 0, 0], "parentToJoint": "identity"},
            {"parent": 0, "type": "REVOLUTE", "axis": [0, 0, 1],
             "mass": 1.0, "com": [0.1, 0, 0],
             "parentToJoint": "translate(1,0,0)"},
            {"parent": 1, "type": "REVOLUTE", "axis": [0, 0, 1],
             "mass": 1.0, "com": [0.1, 0, 0],
             "parentToJoint": "translate(1,0,0)"},
        ],
        "reference": "3-link with parallel axis corrections",
    },
    "TwoLinkWithGravity": {
        "dof": 2,
        "layout": "2-DOF Z-revolute, gravity compensation test",
        "links": [
            {"parent": -1, "type": "REVOLUTE", "axis": [0, 0, 1],
             "mass": 1.0, "com": [0, 0, 0], "parentToJoint": "identity"},
            {"parent": 0, "type": "REVOLUTE", "axis": [0, 0, 1],
             "mass": 1.0, "com": [0, 0, 0],
             "parentToJoint": "translate(1,0,0)"},
        ],
        "reference": "Gravity compensation: g=(0,0,-9.81), tau_g = J^T * F_g",
    },
}

# ── Parsing ────────────────────────────────────────────────────────────────

RE_CHECK_LINE = re.compile(
    r"\s*(?P<status>[✓✗])\s+(?P<name>\w+)\s+"
    r"max_err=(?P<max>[0-9eE.+-]+)\s+mean_err=(?P<mean>[0-9eE.+-]+)"
)

RE_MODEL_HEADER = re.compile(
    r"\s*(?P<status>PASS|FAIL)\s+\|\s+(?P<name>\w+)\s+\|\s+DOF=(?P<dof>\d+)"
)


def parse_benchmark_output(text: str):
    """Parse the C++ benchmark output into structured results."""
    results = {}
    current_model = None

    for line in text.splitlines():
        m = RE_MODEL_HEADER.search(line)
        if m:
            current_model = m.group("name")
            results[current_model] = {
                "pass": m.group("status") == "PASS",
                "checks": [],
            }
            continue

        if current_model is None:
            continue

        m = RE_CHECK_LINE.search(line)
        if m:
            results[current_model]["checks"].append({
                "name": m.group("name"),
                "pass": m.group("status") == "✓",
                "max_rel_err": float(m.group("max")),
                "mean_rel_err": float(m.group("mean")),
            })

    return results


# ── JSON report ────────────────────────────────────────────────────────────

def build_report(bench_results: dict) -> dict:
    """Build the structured JSON report from raw benchmark results."""
    report = {
        "schema_version": "1.0",
        "generated": datetime.utcnow().isoformat() + "Z",
        "toolchain": {
            "benchmark": str(DEFAULT_BENCH),
            "benchmark_type": "C++ GTest",
            "solver_a": "SpatialAlgebra (C++)",
            "solver_b": "Pinocchio (C++)",
            "tolerance": "1e-6",
        },
        "models": {},
        "summary": {"total": 0, "passed": 0, "failed": 0},
    }

    for model_name, model_data in TEST_MODELS.items():
        bench = bench_results.get(model_name, {})
        checks = bench.get("checks", [])
        model_pass = bench.get("pass", False)

        report["models"][model_name] = {
            "definition": model_data,
            "result": {
                "pass": model_pass,
                "checks": {
                    c["name"]: {
                        "pass": c["pass"],
                        "max_rel_err": c["max_rel_err"],
                        "mean_rel_err": c["mean_rel_err"],
                    }
                    for c in checks
                },
            },
        }
        report["summary"]["total"] += 1
        if model_pass:
            report["summary"]["passed"] += 1
        else:
            report["summary"]["failed"] += 1

    return report


# ── CLI ────────────────────────────────────────────────────────────────────

def run_benchmark(bench_path: Path) -> str:
    """Run the C++ comparison benchmark and return stdout."""
    if not bench_path.exists():
        raise FileNotFoundError(
            f"Benchmark executable not found: {bench_path}\n"
            f"Build with: cmake -B build -DSA_BUILD_PINOCCHIO_BENCHMARKS=ON "
            f"&& cmake --build build"
        )
    print(f"Running: {bench_path}", file=sys.stderr)
    result = subprocess.run(
        ["/usr/bin/env", "DYLD_LIBRARY_PATH=/usr/local/lib", str(bench_path)],
        capture_output=True, text=True, timeout=120, cwd=str(BUILD_DIR),
    )
    if result.returncode not in (0, 1):  # 0=pass, 1=test failures, >1=crash
        print(result.stderr, file=sys.stderr)
        print(result.stdout, file=sys.stderr)
        raise RuntimeError(f"Benchmark crashed (exit {result.returncode})")
    return result.stdout


def main():
    parser = argparse.ArgumentParser(
        description="Pinocchio cross-validation Python harness (Phase 23)",
    )
    parser.add_argument(
        "--bench", type=Path, default=DEFAULT_BENCH,
        help="Path to bench_pinocchio_comparison executable",
    )
    parser.add_argument(
        "--output", "-o", type=Path, default=DEFAULT_OUTPUT,
        help="JSON output path",
    )
    parser.add_argument(
        "--print", "-p", action="store_true",
        help="Print the report as formatted text instead of writing JSON",
    )
    args = parser.parse_args()

    # Run benchmark
    stdout = run_benchmark(args.bench)

    # Parse
    bench_results = parse_benchmark_output(stdout)

    # Build report
    report = build_report(bench_results)

    # Output
    if args.print:
        print(json.dumps(report, indent=2))
    else:
        with open(args.output, "w") as f:
            json.dump(report, f, indent=2)
        print(f"Wrote: {args.output}", file=sys.stderr)

    # Exit code matches the aggregate result
    sys.exit(0 if report["summary"]["failed"] == 0 else 1)


if __name__ == "__main__":
    main()
