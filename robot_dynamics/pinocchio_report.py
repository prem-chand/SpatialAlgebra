#!/usr/bin/env python3
"""
Phase 24 — Cross-Library Comparison Report Generator (RPT-01, RPT-02, RPT-03)

Reads the JSON results from Phase 23 and produces:
  - A formatted markdown report with per-model, per-check comparison tables
  - CI-integrated pass/fail assertions against configured tolerance thresholds
  - Summary statistics across all models and checks

Usage:
    python3 robot_dynamics/pinocchio_report.py                              # writes .md
    python3 robot_dynamics/pinocchio_report.py --ci                         # CI mode (exit 1 on failure)
    python3 robot_dynamics/pinocchio_report.py -o docs/comparison_report.md
"""

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_JSON = PROJECT_ROOT / "robot_dynamics" / "pinocchio_comparison_results.json"
DEFAULT_REPORT = PROJECT_ROOT / "robot_dynamics" / "pinocchio_comparison_report.md"

# ── Tolerance thresholds (RPT-03) ─────────────────────────────────────────
# Kinematics checks (joint transform, COM):   1e-8
# Dynamics checks (torques, accelerations):    1e-6
# Gravity checks:                              1e-6
# Mass matrix:                                 1e-6

TOLERANCES = {
    "getJointTransform":     1e-8,
    "getLinkCOM":            1e-8,
    "computeTorques":        1e-6,
    "computeAccelerations":  1e-6,
    "computeGravityTorques": 1e-6,
    "computeMassMatrix":     1e-6,
}

CHECK_ORDER = [
    "computeTorques",
    "computeAccelerations",
    "computeMassMatrix",
    "computeGravityTorques",
    "getJointTransform",
    "getLinkCOM",
]

CHECK_LABELS = {
    "computeTorques": "RNEA (τ = ID(qddot, g))",
    "computeAccelerations": "ABA (qddot = FD(τ, g))",
    "computeMassMatrix": "Mass matrix (CRBA / RNEA-column)",
    "computeGravityTorques": "Gravity compensation torques",
    "getJointTransform": "Forward kinematics (joint transforms)",
    "getLinkCOM": "Center of mass (world-frame COM)",
}

CHECK_CATEGORY = {
    "getJointTransform": "kinematics",
    "getLinkCOM": "kinematics",
    "computeTorques": "dynamics",
    "computeAccelerations": "dynamics",
    "computeGravityTorques": "dynamics",
    "computeMassMatrix": "dynamics",
}


def load_results(path: Path) -> dict:
    with open(path) as f:
        return json.load(f)


# ── Markdown report (RPT-01, RPT-02) ──────────────────────────────────────

def fmt_err(val: float) -> str:
    """Format error value: compact notation for scientific range."""
    if val < 1e-12:
        return "0.00e+00"
    return f"{val:.2e}"


def fmt_status(pass_: bool, max_err: float, tol: float) -> str:
    """Checkmark or cross with threshold context."""
    if pass_:
        return "✓"
    return f"✗ (max={max_err:.2e}, tol={tol:.0e})"


def generate_report(results: dict) -> str:
    """Build the full markdown report string."""
    lines = []
    ts = results.get("generated", "unknown")
    summary = results["summary"]

    # ── Header ────────────────────────────────────────────────────────
    lines.append("# SpatialAlgebra × Pinocchio — Cross-Validation Report")
    lines.append("")
    lines.append(f"**Generated**: {ts}")
    lines.append(f"**Solver A**: {results['toolchain']['solver_a']}")
    lines.append(f"**Solver B**: {results['toolchain']['solver_b']}")
    lines.append(f"**Benchmark**: `{results['toolchain']['benchmark']}`")
    lines.append(f"**Tolerance**: kinematics ≤ 1e-8, dynamics ≤ 1e-6")
    lines.append("")
    lines.append(f"**Summary**: {summary['passed']}/{summary['total']} models pass, "
                  f"{summary['failed']} fail")
    lines.append("")

    # ── Per-model tables ──────────────────────────────────────────────
    for model_name in sorted(results["models"].keys()):
        model = results["models"][model_name]
        defn = model["definition"]
        checks = model["result"]["checks"]
        model_pass = model["result"]["pass"]

        lines.append(f"---")
        lines.append(f"## {model_name}")
        lines.append("")
        lines.append(f"**DOF**: {defn['dof']}  |  "
                      f"**Layout**: {defn['layout']}  |  "
                      f"**Status**: {'✅ PASS' if model_pass else '❌ FAIL'}")
        lines.append("")
        if defn.get("reference"):
            lines.append(f"*{defn['reference']}*")
            lines.append("")

        # Detail table
        lines.append("| Check | Max Rel Err | Mean Rel Err | Threshold | Status |")
        lines.append("|-------|-------------|-------------|-----------|--------|")

        for name in CHECK_ORDER:
            if name not in checks:
                continue
            c = checks[name]
            tol = TOLERANCES.get(name, 1e-6)
            cat = CHECK_CATEGORY.get(name, "dynamics")
            label = CHECK_LABELS.get(name, name)
            max_s = fmt_err(c["max_rel_err"])
            mean_s = fmt_err(c["mean_rel_err"])
            tol_s = f"{tol:.0e}"
            status = "✓" if c["pass"] else "✗"
            lines.append(f"| {label} | {max_s} | {mean_s} | {tol_s} | {status} |")

        lines.append("")

    # ── Summary table ─────────────────────────────────────────────────
    lines.append("---")
    lines.append("## Cross-Model Summary")
    lines.append("")
    lines.append("| Model | Overall | RNEA | ABA | Mass Matrix | Gravity | FK | COM |")
    lines.append("|-------|---------|------|-----|-------------|---------|-----|-------|")

    for model_name in sorted(results["models"].keys()):
        checks = results["models"][model_name]["result"]["checks"]
        overall = "✅" if results["models"][model_name]["result"]["pass"] else "❌"
        cells = [overall]
        for name in CHECK_ORDER:
            if name in checks:
                c = checks[name]
                cells.append("✓" if c["pass"] else "✗")
            else:
                cells.append("—")
        lines.append(f"| {model_name} | " + " | ".join(cells) + " |")

    lines.append("")
    lines.append("---")
    lines.append(f"*Report generated by Phase 24 — Result Reporting*")
    lines.append("")

    return "\n".join(lines)


# ── CI assertions (RPT-03) ────────────────────────────────────────────────

def ci_check(results: dict, verbose: bool = True) -> int:
    """
    Validate all results against configured tolerance thresholds.
    Returns 0 if all checks pass, 1 otherwise (suitable for CI exit code).
    """
    failures = []
    total_checks = 0
    failed_checks = 0

    for model_name in sorted(results["models"].keys()):
        checks = results["models"][model_name]["result"]["checks"]
        for name, check in checks.items():
            total_checks += 1
            tol = TOLERANCES.get(name, 1e-6)
            if not check["pass"]:
                failed_checks += 1
                failures.append(
                    f"  ✗ {model_name}.{name}: "
                    f"max_rel_err={check['max_rel_err']:.2e} > {tol:.0e}"
                )

    if verbose:
        print(f"CI check: {total_checks} checks, {failed_checks} failures")
        for f in failures:
            print(f)

    if failed_checks == 0:
        if verbose:
            print("✅ All checks pass within tolerance thresholds.")
        return 0

    return 1


# ── CLI ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Phase 24 — Cross-library comparison report generator",
    )
    parser.add_argument(
        "--input", "-i", type=Path, default=DEFAULT_JSON,
        help="Input JSON results file",
    )
    parser.add_argument(
        "--output", "-o", type=Path, default=DEFAULT_REPORT,
        help="Output markdown report path",
    )
    parser.add_argument(
        "--ci", action="store_true",
        help="CI mode: exit with code 1 if any check exceeds tolerance",
    )
    parser.add_argument(
        "--print", "-p", action="store_true",
        help="Print report to stdout instead of writing to file",
    )
    args = parser.parse_args()

    results = load_results(args.input)

    if args.ci:
        code = ci_check(results)
        sys.exit(code)

    report = generate_report(results)

    if args.print:
        print(report)
    else:
        with open(args.output, "w") as f:
            f.write(report)
        print(f"Wrote: {args.output}", file=sys.stderr)

    # Also show CI status
    code = ci_check(results, verbose=True)
    if code != 0:
        print(f"\n⚠️  {code} models fail tolerance thresholds — see report.", file=sys.stderr)


if __name__ == "__main__":
    main()
