#!/usr/bin/env python3
"""Regression test suite for GFoil — forward and AD solvers."""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
GOLDEN_DIR = Path(__file__).parent / "golden"
TEST_INPUT = Path(__file__).parent / "input.json"

FWD_SCALAR_KEYS = ["CL", "CD", "CM", "OASPL"]
AD_SCALAR_KEYS  = ["d cl / d alpha", "d cd / d alpha", "d OASPL / d alpha"]
AD_ARRAY_KEYS   = ["d cl / d ycoords", "d cd / d ycoords", "d OASPL / d ycoords"]


# ---------------------------------------------------------------------------
# Build / run helpers
# ---------------------------------------------------------------------------

def do_build(build_dir: Path) -> None:
    print(f"Running cmake in {build_dir} ...")
    r = subprocess.run(["cmake", ".."], cwd=build_dir, capture_output=True, text=True)
    if r.returncode != 0:
        print("cmake failed:\n", r.stderr)
        sys.exit(1)

    nproc = subprocess.run(["nproc"], capture_output=True, text=True).stdout.strip() or "1"
    print(f"Running make -j{nproc} ...")
    r = subprocess.run(["make", f"-j{nproc}"], cwd=build_dir, capture_output=True, text=True)
    if r.returncode != 0:
        print("make failed:\n", r.stderr)
        sys.exit(1)
    print("Build complete.\n")


def run_solvers(build_dir: Path) -> None:
    """Run forward and AD solvers via the pybind11 module and write output files."""
    if not TEST_INPUT.exists():
        print(f"ERROR: test input not found: {TEST_INPUT}")
        sys.exit(1)

    # Make the GFoil package importable from the repo root.
    if str(REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(REPO_ROOT))

    try:
        import importlib
        import GFoil.gfoil_cpp as cpp
        importlib.reload(cpp)
    except ImportError as e:
        print(f"ERROR: could not import GFoil.gfoil_cpp — is the module built?\n{e}")
        sys.exit(1)

    inp = json.loads(TEST_INPUT.read_text())

    # ── forward solve ─────────────────────────────────────────────────────────
    print("  Running forward solver (pybind11) ...")
    fwd = cpp.run_forward(inp)
    if not fwd.get("conv", 0):
        print(f"ERROR: forward solver did not converge (failure_mode={fwd.get('failure_mode', '')})")
        sys.exit(1)

    out = {"CL": fwd["CL"], "CD": fwd["CD"], "CM": fwd["CM"], "OASPL": fwd["OASPL"]}
    (REPO_ROOT / "out.json").write_text(json.dumps(out, indent=2))

    # ── AD solve ──────────────────────────────────────────────────────────────
    print("  Running AD solver (pybind11) ...")
    g = cpp.run_AD(inp, fwd["jacobian"])

    ad = {
        "d cl / d alpha":    g["dCL_dalpha"],
        "d cd / d alpha":    g["dCD_dalpha"],
        "d OASPL / d alpha": g["dOASPL_dalpha"],
        "d cl / d ycoords":  g["dCL_dy"],
        "d cd / d ycoords":  g["dCD_dy"],
        "d OASPL / d ycoords": g["dOASPL_dy"],
    }
    (REPO_ROOT / "ad_gradients.json").write_text(json.dumps(ad, indent=2))


def load_outputs():
    out = json.loads((REPO_ROOT / "out.json").read_text())
    ad  = json.loads((REPO_ROOT / "ad_gradients.json").read_text())
    return out, ad


# ---------------------------------------------------------------------------
# Comparison helpers
# ---------------------------------------------------------------------------

def rel_err(current: float, golden: float) -> float:
    if abs(golden) > 1e-14:
        return abs(current - golden) / abs(golden)
    return abs(current - golden)


def compare_scalar(label: str, current: float, golden: float, tol: float) -> bool:
    err = rel_err(current, golden)
    passed = err < tol
    tag = "PASS" if passed else "FAIL"
    print(f"  [{tag}] {label:<34s}  current={current:.8g}  golden={golden:.8g}  rel_err={err:.3e}")
    return passed


def compare_array(label: str, current: list, golden: list, tol: float) -> bool:
    max_err = 0.0
    worst_idx = 0
    worst_curr = 0.0
    worst_gold = 0.0

    for i, (c, g) in enumerate(zip(current, golden)):
        e = rel_err(c, g)
        if e > max_err:
            max_err = e
            worst_idx = i
            worst_curr = c
            worst_gold = g

    passed = max_err < tol
    tag = "PASS" if passed else "FAIL"
    print(
        f"  [{tag}] {label:<34s}  max_rel_err={max_err:.3e}"
        f"  (worst: idx={worst_idx}, curr={worst_curr:.6g}, gold={worst_gold:.6g})"
    )
    return passed


# ---------------------------------------------------------------------------
# Main actions
# ---------------------------------------------------------------------------

def create_golden(build_dir: Path) -> None:
    print("Running solvers to generate golden reference ...\n")
    run_solvers(build_dir)
    out, ad = load_outputs()

    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)

    (GOLDEN_DIR / "fwd_scalars.json").write_text(
        json.dumps({k: out[k] for k in FWD_SCALAR_KEYS}, indent=2)
    )
    (GOLDEN_DIR / "ad_scalars.json").write_text(
        json.dumps({k: ad[k] for k in AD_SCALAR_KEYS}, indent=2)
    )
    (GOLDEN_DIR / "ad_gradients.json").write_text(
        json.dumps({k: ad[k] for k in AD_ARRAY_KEYS}, indent=2)
    )

    print(f"\nGolden files written to {GOLDEN_DIR}")
    print("Commit them to the repository before making any code changes.")


def run_tests(build_dir: Path, tol: float) -> None:
    for gf in [GOLDEN_DIR / "fwd_scalars.json",
               GOLDEN_DIR / "ad_scalars.json",
               GOLDEN_DIR / "ad_gradients.json"]:
        if not gf.exists():
            print(f"ERROR: golden file not found: {gf}")
            print("Run --create-golden first on known-good code.")
            sys.exit(1)

    print("Running solvers ...\n")
    run_solvers(build_dir)
    out, ad = load_outputs()

    golden_fwd        = json.loads((GOLDEN_DIR / "fwd_scalars.json").read_text())
    golden_ad_scalars = json.loads((GOLDEN_DIR / "ad_scalars.json").read_text())
    golden_ad_arrays  = json.loads((GOLDEN_DIR / "ad_gradients.json").read_text())

    results = []

    print("\nForward solver scalars:")
    for key in FWD_SCALAR_KEYS:
        results.append(compare_scalar(key, out[key], golden_fwd[key], tol))

    print("\nAD solver scalars:")
    for key in AD_SCALAR_KEYS:
        results.append(compare_scalar(key, ad[key], golden_ad_scalars[key], tol))

    print("\nAD solver arrays:")
    for key in AD_ARRAY_KEYS:
        results.append(compare_array(key, ad[key], golden_ad_arrays[key], tol))

    n_pass  = sum(results)
    n_total = len(results)
    print(f"\n{n_pass}/{n_total} checks passed.")
    sys.exit(0 if all(results) else 1)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="GFoil regression test suite",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python3 tests/regression_test.py --create-golden\n"
            "  python3 tests/regression_test.py --test\n"
            "  python3 tests/regression_test.py --build --test\n"
        ),
    )
    parser.add_argument("--create-golden", action="store_true",
                        help="Run solvers and save outputs as golden reference")
    parser.add_argument("--test", action="store_true",
                        help="Run solvers and compare against golden reference")
    parser.add_argument("--build", action="store_true",
                        help="Run cmake + make before testing")
    parser.add_argument("--build-dir", type=Path, default=None,
                        help="Path to build directory (default: ./build)")
    parser.add_argument("--tol", type=float, default=1e-8,
                        help="Relative tolerance for comparisons (default: 1e-8)")
    args = parser.parse_args()

    build_dir = args.build_dir if args.build_dir is not None else REPO_ROOT / "build"

    if args.build:
        do_build(build_dir)

    if args.create_golden:
        create_golden(build_dir)
    elif args.test:
        run_tests(build_dir, args.tol)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
