#!/usr/bin/env python3
"""Mid-span reduction gate for the general oblique-gust Amiet kernel.

The phase-2 `_vec` kernel (src/include/newAmiet.hpp, "General oblique-gust
kernel") must reduce to the mid-span kernel it replaced, exactly, at x2 = 0.
This is the primary acceptance gate: it is an algebraic identity, not an
approximation, because kappa_bar(xi=0) = mu_bar, l_y(K2=0) = b_c*U_c/omega and
S0(y=0) is the mid-span S0.

Usage:
    python3 tests/midspan_reduction_check.py --dump ref.json   # on the OLD build
    python3 tests/midspan_reduction_check.py --check ref.json  # on the NEW build

The reference file is a throwaway (not committed): it is regenerated from the
pre-change build whenever the gate is re-run from scratch. The permanent
regression cover for this kernel is tests/golden/amiet_kernel_scalars.json.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from GFoil import noise_run  # noqa: E402

# Plausible turbulent TE BL states (same fixtures as tests/rotor_noise_test.py).
BL_SUCTION = np.array([1.80e-3, 2.60e-3, 5.2, 62.0, 1400.0, 3.1, 1.85e-2])
BL_PRESSURE = np.array([1.30e-3, 1.70e-3, 4.1, 58.0, -450.0, 3.6, 1.40e-2])

# Mid-span observers (y = 0 exactly), spanning fore/aft, near/far and overhead.
OBSERVERS = [
    (1.0, 0.0, 1.0),
    (0.25, 0.0, 3.0),
    (-2.0, 0.0, 5.0),
    (0.5, 0.0, 50.0),
    (-0.1, 0.0, 0.8),
    (3.0, 0.0, 0.4),
]

CASES = [
    # label,             chord, Re,     nu,      span, model
    ("kam Re=6e5", 0.15, 6.0e5, 1.48e-5, 2.0, "kam"),
    ("roz Re=2e6", 0.30, 2.0e6, 1.48e-5, 3.0, "roz"),
    ("goo Re=1e6", 0.10, 1.0e6, 1.48e-5, 1.0, "goo"),
    ("lee Re=3e6", 0.50, 3.0e6, 1.48e-5, 2.5, "lee"),
    ("tno Re=8e5", 0.20, 8.0e5, 1.48e-5, 1.5, "tno"),
]

FREQS = np.logspace(np.log10(200.0), np.log10(20000.0), 64)


def _values() -> dict:
    out = {}
    for label, chord, Re, nu, span, model in CASES:
        for iObs, obs in enumerate(OBSERVERS):
            # (a) BL-state path
            r = noise_run(BL_SUCTION, BL_PRESSURE, freqs_Hz=FREQS,
                          observerXYZ=np.array(obs), Re=Re, nu=nu, chord=chord,
                          span=span, alphaDeg=0.0, rho=1.225, model=model)
            out[f"BL|{label}|{iObs}"] = r.FF_spectra[0].tolist()

    # (b) custom-WPS path: an analytic spectrum, one surface then both.
    wps = 1e-3 * (FREQS / 1000.0) ** -1.7
    for iObs, obs in enumerate(OBSERVERS):
        for tag, cw in (("upper_only", np.column_stack([wps, np.zeros_like(wps)])),
                        ("both", np.column_stack([wps, 0.6 * wps]))):
            BL_t, BL_b = np.zeros(7), np.zeros(7)
            BL_t[3], BL_b[3] = 62.0, 58.0
            r = noise_run(BL_t, BL_b, freqs_Hz=FREQS, observerXYZ=np.array(obs),
                          Re=1.2e6, nu=1.48e-5, chord=0.18, span=2.0,
                          alphaDeg=0.0, rho=1.225, custom_WPS=cw)
            out[f"CW|{tag}|{iObs}"] = r.FF_spectra[0].tolist()
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--dump", metavar="FILE", help="write reference (run on OLD build)")
    g.add_argument("--check", metavar="FILE", help="compare (run on NEW build)")
    a = p.parse_args()

    vals = _values()
    if a.dump:
        Path(a.dump).write_text(json.dumps(vals))
        n = sum(len(v) for v in vals.values())
        print(f"wrote {len(vals)} spectra ({n} values) to {a.dump}")
        return 0

    ref = json.loads(Path(a.check).read_text())
    if set(ref) != set(vals):
        print("FAIL: case sets differ between reference and current")
        return 1

    worst, worst_key, n_bit, n_tot = 0.0, "", 0, 0
    for k in sorted(ref):
        r = np.array(ref[k])
        c = np.array(vals[k])
        if not np.all(np.isfinite(c)):
            print(f"FAIL: non-finite values in {k}")
            return 1
        nz = r != 0.0
        rel = np.zeros_like(r)
        rel[nz] = np.abs(c[nz] - r[nz]) / np.abs(r[nz])
        rel[~nz] = np.abs(c[~nz])
        n_bit += int(np.sum(c == r))
        n_tot += c.size
        m = float(rel.max())
        if m > worst:
            worst, worst_key = m, k

    print(f"  cases compared        {len(ref)}  ({n_tot} spectral values)")
    print(f"  bit-identical values  {n_bit}/{n_tot}  ({100.0 * n_bit / n_tot:.2f}%)")
    print(f"  worst relative error  {worst:.3e}   (case: {worst_key})")
    ok = worst <= 1e-13
    print(f"\n  [{'PASS' if ok else 'FAIL'}] mid-span reduction, rtol <= 1e-13")
    if n_bit == n_tot:
        print("  [INFO] reduction is BIT-IDENTICAL across every case.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
