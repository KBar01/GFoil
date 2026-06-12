#!/usr/bin/env python3
"""A/B benchmark for the run_AD tape reset strategy (resetHard vs reset).

Measures, in a single process, on the golden case (tests/input.json):
  1. AD wall time over N consecutive run_AD calls (report median).
  2. VmRSS (from /proc/self/status) after every call; plateau = median of
     the last 10 readings.
  3. Optionally (--alternating) the same with alternating fwd_run/run_AD
     pairs — the realistic optimisation pattern.

Stores the FULL gradient output of every call so variants can be compared
bit-for-bit afterwards (correctness gate: any difference on 2nd+ calls
between variants is a failure).

Usage:
  python3 bench/tape_reset_ab.py --label A --reps 30
  python3 bench/tape_reset_ab.py --label B --reps 30 --alternating
Writes bench/results/tape_reset_<label>.json (and _alt.json).
"""

import argparse
import json
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT))

import GFoil.gfoil_cpp as cpp  # noqa: E402

GRAD_KEYS = ["dCL_dy", "dCD_dy", "dOASPL_dy",
             "dCL_dalpha", "dCD_dalpha", "dOASPL_dalpha"]


def vm_rss_kb() -> int:
    for line in Path("/proc/self/status").read_text().splitlines():
        if line.startswith("VmRSS:"):
            return int(line.split()[1])
    raise RuntimeError("VmRSS not found")


def snapshot_grads(g) -> dict:
    out = {}
    for k in GRAD_KEYS:
        v = g[k]
        out[k] = list(v) if hasattr(v, "__len__") else v
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--reps", type=int, default=30)
    ap.add_argument("--alternating", action="store_true",
                    help="also run alternating fwd_run/run_AD measurement")
    args = ap.parse_args()

    inp = json.loads((REPO_ROOT / "tests" / "input.json").read_text())
    fwd = cpp.run_forward(inp)
    assert fwd.get("conv", 0), "forward solve did not converge"
    jac = fwd["jacobian"]

    # ── consecutive run_AD calls ────────────────────────────────────────────
    times_ms, rss_kb, grads = [], [], []
    for i in range(args.reps):
        t0 = time.perf_counter()
        g = cpp.run_AD(inp, jac)
        times_ms.append((time.perf_counter() - t0) * 1e3)
        rss_kb.append(vm_rss_kb())
        grads.append(snapshot_grads(g))

    result = {
        "label": args.label,
        "reps": args.reps,
        "ad_times_ms": times_ms,
        "ad_median_ms": sorted(times_ms)[len(times_ms) // 2],
        "rss_kb": rss_kb,
        "rss_plateau_kb": sorted(rss_kb[-10:])[5],
        "grads": grads,
    }
    out = REPO_ROOT / "bench" / "results" / f"tape_reset_{args.label}.json"
    out.write_text(json.dumps(result))
    print(f"[{args.label}] AD median {result['ad_median_ms']:.1f} ms, "
          f"RSS plateau {result['rss_plateau_kb'] / 1024:.1f} MB -> {out}")

    # ── alternating fwd/AD pattern ──────────────────────────────────────────
    if args.alternating:
        rss_alt = []
        for i in range(args.reps):
            f = cpp.run_forward(inp)
            assert f.get("conv", 0)
            cpp.run_AD(inp, f["jacobian"])
            rss_alt.append(vm_rss_kb())
        alt = {"label": args.label, "rss_kb": rss_alt,
               "rss_plateau_kb": sorted(rss_alt[-10:])[5]}
        out_alt = REPO_ROOT / "bench" / "results" / f"tape_reset_{args.label}_alt.json"
        out_alt.write_text(json.dumps(alt))
        print(f"[{args.label}] alternating fwd/AD RSS plateau "
              f"{alt['rss_plateau_kb'] / 1024:.1f} MB -> {out_alt}")


if __name__ == "__main__":
    main()
