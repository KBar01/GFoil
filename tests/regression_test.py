#!/usr/bin/env python3
"""Regression test suite for GFoil — forward and AD solvers.

Covers three groups:
  1. Free-transition golden case (tests/input.json, NACA 0012 sharp-TE coords,
     alpha=2, Re=2e6, nCrit=5, model roz): forward scalars + AD scalars/arrays.
  2. Forced-transition golden case (same conditions, transition forced at
     x/c = 0.1 on both surfaces): exercises the taped-xift adjoint path
     (see CHANGELOG "Forced-transition adjoint fix").
  3. Windowed-Amiet anchors (config embedded below, ANALYTIC open-TE NACA 0012
     — NOT the input.json coords): scalar-0.98 and window-[0.95,0.99] OASPL
     reference values from CHANGELOG "Windowed Amiet TE sample / Reference
     values", tolerance 1e-6 relative.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
GOLDEN_DIR = Path(__file__).parent / "golden"
TEST_INPUT = Path(__file__).parent / "input.json"

FWD_SCALAR_KEYS = ["CL", "CD", "CM", "OASPL"]
AD_SCALAR_KEYS  = ["d cl / d alpha", "d cd / d alpha", "d OASPL / d alpha"]
AD_ARRAY_KEYS   = ["d cl / d ycoords", "d cd / d ycoords", "d OASPL / d ycoords"]

# Golden cases: (name, golden-file prefix, input-dict overrides)
CASES = [
    ("free transition",   "",        {}),
    ("forced transition", "forced_", {"toptrans": 0.1, "bottrans": 0.1,
                                      "forcetrans": 1}),
]

# ---------------------------------------------------------------------------
# Windowed-Amiet anchor checks (CHANGELOG "Reference values")
# ---------------------------------------------------------------------------
# Full config embedded here on purpose — these anchors do NOT use the
# input.json coordinates. The recorded values were produced on an analytic
# open-TE NACA 0012 (classic 4-digit thickness polynomial with the -0.1015
# trailing coefficient, cosine spacing, 301 points, TE-lower -> LE -> TE-upper
# ordering). Verified: the input.json sharp-TE coords give 63.05499/63.24251
# instead — the foil is part of the anchor config.

ANCHOR_TOL = 1e-6
ANCHORS = [
    ("anchor: scalar TEsample=0.98",    0.98, 0.98, 63.18598),
    ("anchor: window [0.95,0.99]",      0.95, 0.99, 63.36702),
]


def naca0012_analytic(n_half: int = 150):
    """Analytic open-TE NACA 0012, cosine-spaced, 2*n_half+1 points."""
    import numpy as np
    beta = np.linspace(0.0, np.pi, n_half + 1)
    xs = 0.5 * (1 + np.cos(beta))                  # 1 -> 0
    def yt(x):
        return 5 * 0.12 * (0.2969 * np.sqrt(x) - 0.1260 * x
                           - 0.3516 * x**2 + 0.2843 * x**3 - 0.1015 * x**4)
    x = np.concatenate([xs, xs[::-1][1:]])
    y = np.concatenate([-yt(xs), yt(xs[::-1][1:])])
    return x.tolist(), y.tolist()


def anchor_input(sample_lo: float, sample_hi: float) -> dict:
    x, y = naca0012_analytic()
    return {
        "xcoords": x, "ycoords": y,
        "alpha_degrees": 2.0,
        "Re": 2.0e6, "Ma": 0.0, "rho": 1.225, "nu": 1.5e-5,
        "restart": 0,
        "sampleTE": sample_lo, "sampleTE_hi": sample_hi,
        "X": [1.0], "Y": [0.0], "Z": [1.0],
        "S": 3.0,
        "ncrit": 5.0, "rtol": 1e-6,
        "Ufac": 1.0, "TEfac": 0.09,
        "toptrans": 1.0, "bottrans": 1.0, "forcetrans": 0,
        "model": "kam", "aWeighting": 0,
        "chord": 1.0, "f_min": 200.0, "f_max": 20000.0,
        "verbose": False,
    }


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


def import_module():
    if str(REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(REPO_ROOT))
    try:
        import GFoil.gfoil_cpp as cpp
    except ImportError as e:
        print(f"ERROR: could not import GFoil.gfoil_cpp — is the module built?\n{e}")
        sys.exit(1)
    return cpp


def run_case(cpp, overrides: dict, label: str):
    """Forward + AD solve for one golden case. Returns (out, ad) dicts."""
    if not TEST_INPUT.exists():
        print(f"ERROR: test input not found: {TEST_INPUT}")
        sys.exit(1)
    inp = json.loads(TEST_INPUT.read_text())
    inp.update(overrides)

    print(f"  [{label}] forward solve ...")
    fwd = cpp.run_forward(inp)
    if not fwd.get("conv", 0):
        print(f"ERROR: forward solver did not converge "
              f"(failure_mode={fwd.get('failure_mode', '')})")
        sys.exit(1)
    out = {k: fwd[k] for k in FWD_SCALAR_KEYS}

    print(f"  [{label}] AD solve ...")
    g = cpp.run_AD(inp, fwd["jacobian"])
    ad = {
        "d cl / d alpha":      g["dCL_dalpha"],
        "d cd / d alpha":      g["dCD_dalpha"],
        "d OASPL / d alpha":   g["dOASPL_dalpha"],
        "d cl / d ycoords":    g["dCL_dy"],
        "d cd / d ycoords":    g["dCD_dy"],
        "d OASPL / d ycoords": g["dOASPL_dy"],
    }
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
    print(f"  [{tag}] {label:<40s}  current={current:.8g}  golden={golden:.8g}  rel_err={err:.3e}")
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
        f"  [{tag}] {label:<40s}  max_rel_err={max_err:.3e}"
        f"  (worst: idx={worst_idx}, curr={worst_curr:.6g}, gold={worst_gold:.6g})"
    )
    return passed


# ---------------------------------------------------------------------------
# Main actions
# ---------------------------------------------------------------------------

def create_golden(build_dir: Path) -> None:
    cpp = import_module()
    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)

    for name, prefix, overrides in CASES:
        print(f"Generating golden reference: {name} ...")
        out, ad = run_case(cpp, overrides, name)
        (GOLDEN_DIR / f"{prefix}fwd_scalars.json").write_text(
            json.dumps(out, indent=2)
        )
        (GOLDEN_DIR / f"{prefix}ad_scalars.json").write_text(
            json.dumps({k: ad[k] for k in AD_SCALAR_KEYS}, indent=2)
        )
        (GOLDEN_DIR / f"{prefix}ad_gradients.json").write_text(
            json.dumps({k: ad[k] for k in AD_ARRAY_KEYS}, indent=2)
        )

    print(f"\nGolden files written to {GOLDEN_DIR}")
    print("Commit them to the repository before making any code changes.")


def run_anchor_checks(cpp) -> list:
    print("\nWindowed-Amiet anchors (analytic NACA 0012, kam, obs (1,0,1), span 3):")
    results = []
    for label, lo, hi, expected in ANCHORS:
        r = cpp.run_forward(anchor_input(lo, hi))
        if not r.get("conv", 0):
            print(f"  [FAIL] {label:<40s}  forward solve did not converge")
            results.append(False)
            continue
        results.append(compare_scalar(label, r["OASPL"], expected, ANCHOR_TOL))
    return results


def run_tests(build_dir: Path, tol: float) -> None:
    needed = []
    for _, prefix, _ in CASES:
        needed += [GOLDEN_DIR / f"{prefix}fwd_scalars.json",
                   GOLDEN_DIR / f"{prefix}ad_scalars.json",
                   GOLDEN_DIR / f"{prefix}ad_gradients.json"]
    for gf in needed:
        if not gf.exists():
            print(f"ERROR: golden file not found: {gf}")
            print("Run --create-golden first on known-good code.")
            sys.exit(1)

    cpp = import_module()
    results = []

    for name, prefix, overrides in CASES:
        print(f"\nRunning solvers: {name} ...")
        out, ad = run_case(cpp, overrides, name)

        golden_fwd        = json.loads((GOLDEN_DIR / f"{prefix}fwd_scalars.json").read_text())
        golden_ad_scalars = json.loads((GOLDEN_DIR / f"{prefix}ad_scalars.json").read_text())
        golden_ad_arrays  = json.loads((GOLDEN_DIR / f"{prefix}ad_gradients.json").read_text())

        print(f"\n{name} — forward solver scalars:")
        for key in FWD_SCALAR_KEYS:
            results.append(compare_scalar(key, out[key], golden_fwd[key], tol))

        print(f"\n{name} — AD solver scalars:")
        for key in AD_SCALAR_KEYS:
            results.append(compare_scalar(key, ad[key], golden_ad_scalars[key], tol))

        print(f"\n{name} — AD solver arrays:")
        for key in AD_ARRAY_KEYS:
            results.append(compare_array(key, ad[key], golden_ad_arrays[key], tol))

    results += run_anchor_checks(cpp)

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
                        help="Relative tolerance for golden comparisons (default: 1e-8)")
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
