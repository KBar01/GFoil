#!/usr/bin/env python3
"""B.3(iv) — 30-rep forward median wall-time on the golden case under a fix mode.

Run: GFOIL_BFIX=<off|A|B> python3 bench/b3_timing.py
"""
import json
import os
import statistics
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
import GFoil.gfoil_cpp as cpp

inp = json.loads((REPO / "tests" / "input.json").read_text())
# warmup
for _ in range(3):
    cpp.run_forward(inp)
ts = []
for _ in range(30):
    t0 = time.perf_counter()
    cpp.run_forward(inp)
    ts.append((time.perf_counter() - t0) * 1e3)
print(f"GFOIL_BFIX={os.environ.get('GFOIL_BFIX','B')}: "
      f"median {statistics.median(ts):.1f} ms  "
      f"min {min(ts):.1f}  p90 {sorted(ts)[int(0.9*29)]:.1f}")
