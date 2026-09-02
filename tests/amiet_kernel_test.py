#!/usr/bin/env python3
"""Test suite for the general oblique-gust Amiet kernel (phase 2).

Follows the conventions of tests/regression_test.py: a plain script (no pytest
dependency), PASS/FAIL lines, exit code 0 iff every check passes.

Covers the `_vec`-only general Roger & Moreau (2005) kernel restored in
CHANGELOG "General oblique-gust Amiet kernel". Deterministic and self-contained:
the acoustic kernel is driven directly, with no aero solve anywhere.

Groups:
  1. Algebra identities — xi = beta*|x2|/S0, kappa_bar = mu_bar*sqrt(1-xi^2),
     and the Corcos length l_y (Eq. 19) reduction / monotonicity / positivity.
  2. Cut behaviour — the unbridged narrow cut at xi = 1, continuity and endpoint
     slope matching of the bridged |I|, and subcritical decay for xi >> 1.
  3. Golden regression case (tests/golden/amiet_kernel_scalars.json).

The mid-span reduction gate (the general kernel must reproduce the mid-span
kernel it replaced, at x2 = 0, to rtol <= 1e-13) is a separate before/after
comparison against the pre-change build and lives in
tests/midspan_reduction_check.py.

Quick reference
---------------
    python3 tests/amiet_kernel_test.py --test
    python3 tests/amiet_kernel_test.py --create-golden

Regenerating the golden requires JUSTIFICATION and a CHANGELOG entry, the same
rule as the other goldens (see tests/README.md, "Updating the golden files").
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from GFoil import gfoil_cpp  # noqa: E402

GOLDEN_DIR = Path(__file__).parent / "golden"
GOLDEN_FILE = GOLDEN_DIR / "amiet_kernel_scalars.json"
GOLDEN_TOL = 1e-10

C0 = 340.0

# Paper Fig. 11 conditions (R&M 2005): c = 0.13 m, M = 0.05, f = 200 / 1000 Hz,
# mid-span observer swept over polar angle.
FIG11 = dict(chord=0.13, M=0.05, Ue=0.05 * C0)
FIG11_FREQS = [200.0, 1000.0]
FIG11_THETAS = [20.0, 60.0, 90.0, 120.0, 160.0]
R_OBS = 2.0


def _report(label: str, passed: bool, detail: str = "") -> bool:
    print(f"  [{'PASS' if passed else 'FAIL'}] {label:<50s}  {detail}")
    return passed


def _rel(a, b):
    return abs(a - b) / abs(b) if abs(b) > 1e-14 else abs(a - b)


def _obs(theta_deg, R=R_OBS):
    t = np.deg2rad(theta_deg)
    return R * np.cos(t), R * np.sin(t)


def kernel(chord, M, Ue, x, z, freqs, K2_bar, bridged=False):
    """Validation-only entry point: |I| with K2_bar a free parameter."""
    f = np.atleast_1d(np.asarray(freqs, dtype=float))
    k = np.atleast_1d(np.asarray(K2_bar, dtype=float))
    if k.size == 1 and f.size > 1:
        k = np.full(f.size, float(k[0]))
    if f.size == 1 and k.size > 1:
        f = np.full(k.size, float(f[0]))
    r = gfoil_cpp.amiet_kernel_I(dict(
        chord=chord, M=M, Ue=Ue, x=x, z=z, bridged=bool(bridged),
        freqs_Hz=list(f), K2_bar=list(k)))
    return {kk: np.asarray(v) for kk, v in r.items()}


# --------------------------------------------------------------------------- #
# 1. Algebra identities                                                       #
# --------------------------------------------------------------------------- #


def test_algebra() -> list:
    print("\n1. Algebra identities:")
    out = []

    chord, M, Ue = FIG11["chord"], FIG11["M"], FIG11["Ue"]
    beta = np.sqrt(1.0 - M * M)
    b = chord / 2.0

    # (a) xi = K2_bar/(beta*mu_bar) as the kernel reports it, against an
    #     independent computation from the wavenumber definitions.
    freqs = np.logspace(np.log10(200.0), np.log10(20000.0), 24)
    x, z = _obs(55.0)
    K2b = 0.4 * np.ones_like(freqs)
    r = kernel(chord, M, Ue, x, z, freqs, K2b)
    U = M * C0
    K_bar = (2.0 * np.pi * freqs / U) * b
    mu_bar = K_bar * M / (beta ** 2)
    out.append(_report("mu_bar = K_bar*M/beta^2 as kernel reports",
                       np.allclose(r["mu_bar"], mu_bar, rtol=1e-14, atol=0),
                       f"max_rel={np.max(np.abs(r['mu_bar']/mu_bar - 1)):.2e}"))
    xi_indep = np.abs(K2b) / (beta * mu_bar)
    out.append(_report("xi = |K2_bar|/(beta*mu_bar)",
                       np.allclose(r["xi"], xi_indep, rtol=1e-14, atol=0),
                       f"max_rel={np.max(np.abs(r['xi']/xi_indep - 1)):.2e}"))

    # (b) THE geometric identity: with K2_bar set by Eq. 18's selection
    #     (K2_bar = k_bar*x2/S0), xi must equal beta*|x2|/S0 exactly. This is
    #     what forces the production path to be supercritical everywhere.
    worst = 0.0
    for theta in (10.0, 45.0, 90.0, 135.0, 170.0):
        for y in (0.0, 0.3, 1.0, 5.0, 50.0):
            x, z = _obs(theta)
            S0 = np.sqrt(x * x + beta ** 2 * z * z + beta ** 2 * y * y)
            k_bar = (2.0 * np.pi * freqs / C0) * b
            K2b_sel = k_bar * y / S0
            rr = kernel(chord, M, Ue, x, z, freqs, K2b_sel)
            xi_geom = beta * abs(y) / S0
            worst = max(worst, float(np.max(np.abs(rr["xi"] - xi_geom))))
    out.append(_report("xi = beta*|x2|/S0 under the Eq. 18 selection",
                       worst < 1e-14, f"max_abs_err={worst:.2e} (25 geometries)"))

    # (c) ... and that identity's consequence: xi < 1 strictly, i.e. production
    #     never leaves the supercritical branch.
    xis = []
    for theta in (0.5, 30.0, 90.0, 150.0, 179.5):
        for y in (0.0, 0.1, 2.0, 1e3, 1e6):
            x, z = _obs(theta)
            S0 = np.sqrt(x * x + beta ** 2 * z * z + beta ** 2 * y * y)
            xis.append(beta * abs(y) / S0)
    out.append(_report("Eq.18 selection gives xi < 1 always (supercritical)",
                       max(xis) < 1.0, f"max xi = {max(xis):.6f}"))

    # (d) kappa_bar = mu_bar*sqrt(1 - xi^2) on the supercritical side, and
    #     kappa' = mu_bar*sqrt(xi^2 - 1) on the subcritical side. The kernel
    #     reports kappa signed: >0 supercritical, <0 for -kappa'.
    for xi_t, name in ((0.0, "xi=0"), (0.6, "xi=0.6"), (0.98, "xi=0.98")):
        x, z = _obs(70.0)
        rr = kernel(chord, M, Ue, x, z, freqs, xi_t * beta * mu_bar)
        want = mu_bar * np.sqrt(1.0 - xi_t ** 2)
        out.append(_report(f"kappa_bar = mu_bar*sqrt(1-xi^2)  [{name}]",
                           np.allclose(rr["kappa"], want, rtol=1e-13, atol=0),
                           f"max_rel={np.max(np.abs(rr['kappa']/want - 1)):.2e}"))
    for xi_t in (1.5, 4.0):
        x, z = _obs(70.0)
        rr = kernel(chord, M, Ue, x, z, freqs, xi_t * beta * mu_bar)
        want = -mu_bar * np.sqrt(xi_t ** 2 - 1.0)
        out.append(_report(f"kappa' = mu_bar*sqrt(xi^2-1)  [xi={xi_t}]",
                           np.allclose(rr["kappa"], want, rtol=1e-13, atol=0),
                           f"max_rel={np.max(np.abs(rr['kappa']/want - 1)):.2e}"))

    # (e) kappa_bar(xi=0) == mu_bar EXACTLY -- the bit-exactness this relies on
    #     is what makes the mid-span reduction bit-identical in I1.
    x, z = _obs(70.0)
    rr = kernel(chord, M, Ue, x, z, freqs, np.zeros_like(freqs))
    out.append(_report("kappa_bar(xi=0) == mu_bar bit-exactly",
                       bool(np.all(rr["kappa"] == rr["mu_bar"])), ""))

    # (e2) The G_d / G_e denominators D +/- 2*kappa_bar cannot vanish under the
    #      Eq. 18 selection. Radiation_integral2* guards them by clipping to
    #      1e-10 (a pre-existing TODO in the mid-span code), which near a zero
    #      would produce a huge spike rather than the analytic limit. At mid-span
    #      kappa_bar = mu_bar makes D + 2*kappa_bar = mu_bar*(3 - x1/S0) >= 2*mu_bar,
    #      so the mid-span kernel structurally never approached it. With general
    #      kappa_bar the clip IS reachable -- but only through the validation
    #      entry point, which drives K2_bar free of the geometry. Under Eq. 18,
    #      kappa_bar = mu_bar*sqrt(x1^2 + beta^2 x3^2)/S0, hence
    #        D + 2*kappa_bar = (mu_bar/S0)*(3*sqrt(x1^2+beta^2 x3^2) - x1) > 0
    #        D - 2*kappa_bar = -(mu_bar/S0)*(sqrt(x1^2+beta^2 x3^2) + x1) < 0
    #      strictly, since sqrt(x1^2+beta^2 x3^2) >= |x1|. This pins that.
    rng = np.random.default_rng(2005)
    n = 200000
    Mr = rng.uniform(0.01, 0.6, n)
    br = np.sqrt(1.0 - Mr * Mr)
    x1 = rng.normal(0.0, 50.0, n)
    x2 = rng.normal(0.0, 50.0, n)
    x3 = rng.normal(0.0, 50.0, n)
    S0r = np.sqrt(x1 ** 2 + br ** 2 * (x2 ** 2 + x3 ** 2))
    kom = np.sqrt(x1 ** 2 + br ** 2 * x3 ** 2) / S0r      # kappa_bar/mu_bar
    dplus = 3.0 * kom - x1 / S0r                          # (D+2kappa)*S0/mu_bar
    dminus = -kom - x1 / S0r                              # (D-2kappa)*S0/mu_bar
    dzero = kom - x1 / S0r                                # D*S0/mu_bar
    out.append(_report("Eq.18 selection keeps D + 2*kappa_bar > 0",
                       bool(np.all(dplus > 0.0)),
                       f"min={dplus.min():.3e} over {n} random geometries"))
    out.append(_report("Eq.18 selection keeps D - 2*kappa_bar < 0",
                       bool(np.all(dminus < 0.0)),
                       f"max={dminus.max():.3e} over {n} random geometries"))
    # D itself is the third guarded quantity (the G_e branch, which returns 0
    # for |D| < 1e-10 and takes a complex sqrt whose branch flips at D = 0).
    # D = (mu_bar/S0)*(sqrt(x1^2 + beta^2 x3^2) - x1) >= 0, vanishing only for
    # x3 = 0 with x1 > 0 -- an observer in the plate plane downstream, where
    # Eq. 18's prefactor (~x3^2) makes S_pp vanish anyway.
    out.append(_report("Eq.18 selection keeps D > 0",
                       bool(np.all(dzero > 0.0)),
                       f"min={dzero.min():.3e} over {n} random geometries"))

    # (f) Corcos length l_y (Eq. 19).
    U_c = 0.7 * Ue
    l_c = 1.47 * U_c / (2.0 * np.pi * freqs)
    rr = kernel(chord, M, Ue, x, z, freqs, np.zeros_like(freqs))
    out.append(_report("l_y(K2=0) == b_c*U_c/omega bit-identically",
                       bool(np.all(rr["l_y"] == l_c)),
                       "the mid-span Corcos length, exactly"))

    K2_grid = np.linspace(0.0, 400.0, 60)
    f0 = 2000.0
    rr = kernel(chord, M, Ue, x, z, np.full(K2_grid.size, f0), K2_grid * b)
    out.append(_report("l_y > 0 for all K2", bool(np.all(rr["l_y"] > 0.0)),
                       f"min={rr['l_y'].min():.3e}"))
    out.append(_report("l_y monotone decreasing in |K2|",
                       bool(np.all(np.diff(rr["l_y"]) < 0.0)),
                       f"l_y(0)={rr['l_y'][0]:.4e} -> l_y(400)={rr['l_y'][-1]:.4e}"))
    # even in K2: only K2^2 enters
    rr_neg = kernel(chord, M, Ue, x, z, np.full(K2_grid.size, f0), -K2_grid * b)
    out.append(_report("l_y even in K2 (sign of x2 irrelevant)",
                       bool(np.all(rr_neg["l_y"] == rr["l_y"])), ""))
    # closed form
    w = 2.0 * np.pi * f0 / (1.47 * U_c)
    want = w / (K2_grid ** 2 + w ** 2)
    out.append(_report("l_y == (w)/(K2^2 + w^2), w = omega/(b_c*U_c)",
                       np.allclose(rr["l_y"], want, rtol=1e-13, atol=0),
                       f"max_rel={np.max(np.abs(rr['l_y']/want - 1)):.2e}"))
    return out


# --------------------------------------------------------------------------- #
# 2. Cut behaviour                                                            #
# --------------------------------------------------------------------------- #


def test_cut() -> list:
    print("\n2. Cut behaviour at xi = 1 (paper Fig. 11 conditions):")
    out = []
    chord, M, Ue = FIG11["chord"], FIG11["M"], FIG11["Ue"]
    beta = np.sqrt(1.0 - M * M)

    # (a) The unbridged kernel must SHOW the cut: R&M sec.4.1 describe a deep,
    #     narrow dip at kappa_bar = 0 from non-convergence of the two-step
    #     Schwarzschild iteration. If it were absent the bridge would be
    #     pointless, so assert it is there.
    found = 0
    depths = []
    for f in FIG11_FREQS:
        for th in FIG11_THETAS:
            x, z = _obs(th)
            r0 = kernel(chord, M, Ue, x, z, [f], [0.0])
            mu = float(r0["mu_bar"][0])
            # sample tightly across xi = 1
            xi = np.concatenate([np.linspace(0.90, 0.9999, 60),
                                 np.linspace(1.0001, 1.10, 60)])
            r = kernel(chord, M, Ue, x, z, np.full(xi.size, f),
                       xi * beta * mu, bridged=False)
            I = r["I_abs"]
            # depth of the dip at the cut relative to the window edges
            edge = 0.5 * (I[0] + I[-1])
            dip = float(np.min(I)) / edge
            depths.append(dip)
            if dip < 0.9:
                found += 1
    out.append(_report("unbridged |I| shows a dip at xi = 1",
                       found >= 8,
                       f"{found}/10 cases dip >10%; min depth "
                       f"{min(depths):.3f}x edge"))

    # (b) The bridge must join the outer branches without a seam. A cubic
    #     Hermite has h00(0) = 1 and h01(1) = 1, so the joins should be exact to
    #     round-off; anything larger means the window bounds the kernel uses and
    #     the ones it reports have drifted apart.
    worst_seam, worst_lbl = 0.0, ""
    for f in FIG11_FREQS:
        for th in FIG11_THETAS:
            x, z = _obs(th)
            r0 = kernel(chord, M, Ue, x, z, [f], [0.0])
            mu = float(r0["mu_bar"][0])
            xa, xb = float(r0["xi_a"][0]), float(r0["xi_b"][0])
            d = 1e-9 * (xb - xa)
            for anchor, inward in ((xa, +1.0), (xb, -1.0)):
                pts = np.array([anchor, anchor + inward * d])
                rr = kernel(chord, M, Ue, x, z, np.full(2, f),
                            pts * beta * mu, bridged=True)
                # pts[0] takes the raw branch (xi <= xi_a / xi >= xi_b),
                # pts[1] is a hair inside and takes the bridge.
                e = abs(rr["I_abs"][1] - rr["I_abs"][0]) / abs(rr["I_abs"][0])
                if e > worst_seam:
                    worst_seam, worst_lbl = e, f"f={f:.0f} Th={th:.0f}"
    out.append(_report("bridge joins the outer branches with no seam",
                       worst_seam < 1e-7,
                       f"max relative seam jump {worst_seam:.2e} ({worst_lbl})"))

    # (c) Endpoint slopes: the bridge is a cubic Hermite whose end tangents ARE
    #     the outer branches' one-sided finite differences over its own step
    #     (0.01*(xi_b - xi_a), stepping away from the cut). Recompute those
    #     differences here from the raw branch and check the bridge honours them.
    worst_sl, worst_lbl = 0.0, ""
    for f in FIG11_FREQS:
        for th in FIG11_THETAS:
            x, z = _obs(th)
            r0 = kernel(chord, M, Ue, x, z, [f], [0.0])
            mu = float(r0["mu_bar"][0])
            xa, xb = float(r0["xi_a"][0]), float(r0["xi_b"][0])
            L = xb - xa
            h = 0.01 * L
            # raw one-sided slopes, stepping AWAY from the cut
            pts = np.array([xa, xa - h, xb, xb + h])
            rw = kernel(chord, M, Ue, x, z, np.full(4, f), pts * beta * mu,
                        bridged=False)["I_abs"]
            da_ref = (rw[0] - rw[1]) / h
            db_ref = (rw[3] - rw[2]) / h
            # bridge's own end tangents, by a tiny FD just inside each anchor
            d = 1e-7 * L
            pb = np.array([xa, xa + d, xb, xb - d])
            bg = kernel(chord, M, Ue, x, z, np.full(4, f), pb * beta * mu,
                        bridged=True)["I_abs"]
            da_got = (bg[1] - bg[0]) / d
            db_got = (bg[2] - bg[3]) / d
            for got, ref in ((da_got, da_ref), (db_got, db_ref)):
                e = abs(got - ref) / max(abs(ref), 1e-12)
                if e > worst_sl:
                    worst_sl, worst_lbl = e, f"f={f:.0f} Th={th:.0f}"
    # 1e-3: this compares a cubic's end tangent against the finite difference
    # that defined it, probed by a second finite difference. The residual is the
    # cubic's own curvature over the probe step, not an error in the bridge.
    out.append(_report("bridge end tangents == outer-branch FD slopes",
                       worst_sl < 1e-3,
                       f"max relative slope mismatch {worst_sl:.2e} ({worst_lbl})"))

    # (c2) No spike inside the window, in the non-degenerate case mu_bar >
    #      2*kappa_reg. There the window is narrow (its half-width in xi is
    #      ~kappa_reg/mu_bar), the anchor slopes are local, and the cubic cannot
    #      run away. See (c3) for the degenerate low-mu_bar case.
    worst_ov, worst_lbl, n_nd = 0.0, "", 0
    for f in FIG11_FREQS:
        for th in FIG11_THETAS:
            x, z = _obs(th)
            r0 = kernel(chord, M, Ue, x, z, [f], [0.0])
            mu = float(r0["mu_bar"][0])
            if mu <= 2.0 * 0.125:
                continue                      # degenerate; covered by (c3)
            n_nd += 1
            xa, xb = float(r0["xi_a"][0]), float(r0["xi_b"][0])
            xi = np.linspace(xa, xb, 200)
            I = kernel(chord, M, Ue, x, z, np.full(xi.size, f),
                       xi * beta * mu, bridged=True)["I_abs"]
            if not np.all(np.isfinite(I)):
                worst_ov, worst_lbl = np.inf, f"non-finite f={f} th={th}"
                break
            ov = float(np.max(I)) / max(I[0], I[-1])
            if ov > worst_ov:
                worst_ov, worst_lbl = ov, f"f={f:.0f} Th={th:.0f}"
    out.append(_report("bridged |I| does not spike (mu_bar > 2*kappa_reg)",
                       worst_ov < 1.5,
                       f"max |I|/max(endpoints) = {worst_ov:.3f} over {n_nd} "
                       f"cases ({worst_lbl})"))

    # (c3) KNOWN LIMITATION, pinned deliberately. Two things must coincide for
    #      the bridge to overshoot: (i) mu_bar <= 2*kappa_reg, where the
    #      degeneracy cap kappa_reg_eff = mu_bar/2 makes the window wide in xi
    #      (xi_a = sqrt(3)/2), so the anchor slope -- taken over 0.01*(xi_b-xi_a)
    #      per the regularisation design -- is levered across a 100x longer arm;
    #      and (ii) a D + 2*kappa_bar = 0 pole of Eq. 14 inside that wide window,
    #      which makes the anchor slope enormous (measured da ~ 250 at f=200 Hz,
    #      Theta=60 deg, against endpoint values ~1).
    #
    #      This is unreachable from any production path. The pole needs
    #      3*kappa_bar = mu_bar*x1/S0, which test_algebra (e2) shows the Eq. 18
    #      selection forbids; and the bridge's own anchors sit at xi <= xi_obs,
    #      where kappa_bar is LARGER than the observer's, so D + 2*kappa_bar is
    #      even more positive there. Only amiet_kernel_I, which drives K2_bar
    #      free of the geometry, can place a pole inside the window.
    #
    #      Asserted here so the limitation cannot be forgotten, and so that a
    #      future monotone-limited interpolant makes this test fail loudly
    #      rather than silently change the bridge.
    x, z = _obs(60.0)
    r0 = kernel(chord, M, Ue, x, z, [200.0], [0.0])
    mu = float(r0["mu_bar"][0])
    xa, xb = float(r0["xi_a"][0]), float(r0["xi_b"][0])
    xi = np.linspace(xa, xb, 200)
    I = kernel(chord, M, Ue, x, z, np.full(xi.size, 200.0),
               xi * beta * mu, bridged=True)["I_abs"]
    ov = float(np.max(I)) / max(I[0], I[-1])
    out.append(_report("degenerate wide window overshoots (known, validation-only)",
                       (mu <= 0.25) and np.all(np.isfinite(I)) and ov > 2.0,
                       f"mu_bar={mu:.4f} <= 2*kappa_reg, window [{xa:.3f},{xb:.3f}], "
                       f"overshoot {ov:.1f}x -- finite, and off-manifold only"))

    # (c4) The PRODUCTION route stays finite and positive all the way out to the
    #      spanwise axis, including deep inside the bridge window.
    #
    #      This must go through noise_run, not amiet_kernel_I: the validation
    #      hook deliberately keeps the mid-span S0 while K2_bar is driven freely,
    #      so the (K2_bar, S0) pair it forms is off-manifold by construction and
    #      the (e2) proof does not apply to it. noise_run builds both from the
    #      same general S0, which is what the proof needs.
    #
    #      Only finiteness and positivity are asserted -- deliberately. |I|^2
    #      rises one to two orders of magnitude towards the axis (xi -> 1 drives
    #      the selected gust critical, kappa_bar -> 0, where the radiation
    #      integral peaks) AND oscillates in x2: at mu_bar ~ 23 the phase 4*kappa
    #      sweeps tens of radians across this sweep. So neither a magnitude bound
    #      nor a neighbour-ratio bound is well posed here. The guarantee that no
    #      POLE is reached is algebraic, and is pinned by (e2); the physical
    #      sanity of the resulting directivity is pinned by the rotor suite's
    #      axis-vs-plane ordering test.
    from GFoil import noise_run  # noqa: E402  (local: keeps group 1 import-free)

    nu, span, Ue_p = 1.48e-5, 2.0, 62.0
    bad, n_pts, saw_bridge = 0, 0, False
    for chord_p, M_p in ((0.13, 0.05), (0.30, 0.20)):
        beta_p = np.sqrt(1.0 - M_p * M_p)
        U = M_p * C0
        Re = U * chord_p / nu
        b = chord_p / 2.0
        freqs = np.array([200.0, 1000.0, 8000.0])
        omega = 2.0 * np.pi * freqs
        K_bar = (omega / U) * b
        mu_p = K_bar * M_p / beta_p ** 2
        kr = np.where(mu_p > 2.0 * 0.125, 0.125, 0.5 * mu_p)
        xi_a_p = np.sqrt(1.0 - (kr / mu_p) ** 2)
        wps = np.full(freqs.size, 0.05)
        BL_t, BL_b = np.zeros(7), np.zeros(7)
        BL_t[3] = BL_b[3] = Ue_p
        for th in (25.0, 90.0, 165.0):
            x, z = _obs(th)
            yv = np.logspace(-1.0, 7.0, 120)
            for y in yv:
                r = noise_run(BL_t, BL_b, freqs_Hz=freqs,
                              observerXYZ=np.array([x + 0.75 * chord_p, y, z]),
                              Re=Re, nu=nu, chord=chord_p, span=span,
                              alphaDeg=0.0, rho=1.225,
                              custom_WPS=np.column_stack([wps,
                                                          np.zeros_like(wps)]))
                ff = r.FF_spectra[0]
                n_pts += ff.size
                if not np.all(np.isfinite(ff)) or np.any(ff <= 0.0):
                    bad += 1
            xi_sw = beta_p * yv / np.sqrt(x * x + beta_p ** 2 * z * z
                                          + beta_p ** 2 * yv ** 2)
            if np.any(xi_sw[:, None] > xi_a_p[None, :]):
                saw_bridge = True
    out.append(_report("production route finite & positive out to the axis",
                       bad == 0 and saw_bridge,
                       f"{n_pts} spectral values over a 8-decade x2 sweep, "
                       f"{bad} bad; bridge entered: {saw_bridge}"))

    # (d) Subcritical decay for xi >> 1, steeper at the higher frequency.
    #     Assert the sign of the trend, not digitised paper values.
    slopes = []
    for f in FIG11_FREQS:
        x, z = _obs(90.0)
        r0 = kernel(chord, M, Ue, x, z, [f], [0.0])
        mu = float(r0["mu_bar"][0])
        xi = np.logspace(np.log10(1.3), np.log10(20.0), 40)
        r = kernel(chord, M, Ue, x, z, np.full(xi.size, f),
                   xi * beta * mu, bridged=True)
        I = r["I_abs"]
        mono = bool(np.all(np.diff(I) < 0.0))
        sl = float(np.polyfit(np.log(xi[-15:]), np.log(I[-15:]), 1)[0])
        slopes.append(sl)
        out.append(_report(f"subcritical |I| monotone decreasing  [f={f:.0f}]",
                           mono, f"|I|: {I[0]:.3e} -> {I[-1]:.3e}, "
                                 f"log-log slope {sl:+.3f}"))
    out.append(_report("decay steeper at the higher frequency",
                       slopes[1] < slopes[0],
                       f"slope(200Hz)={slopes[0]:+.3f}  "
                       f"slope(1000Hz)={slopes[1]:+.3f}"))

    # (e) The subcritical branch must stay finite far past any physical need:
    #     the e^{-2iA1'} * F* product is O(1) but its factors are e^{-+2 kappa'}.
    x, z = _obs(90.0)
    r0 = kernel(2.0, 0.167, 56.7, -3.0, 900.0, [5000.0], [0.0])
    mu = float(r0["mu_bar"][0])
    beta2 = np.sqrt(1.0 - 0.167 ** 2)
    xi = np.array([12.0, 50.0, 200.0])
    r = kernel(2.0, 0.167, 56.7, -3.0, 900.0, np.full(3, 5000.0),
               xi * beta2 * mu, bridged=False)
    kp = -r["kappa"]
    out.append(_report("no overflow at extreme kappa' (erfcx path)",
                       bool(np.all(np.isfinite(r["I_abs"]))),
                       f"kappa' up to {kp.max():.0f}, |I| all finite"))
    return out


# --------------------------------------------------------------------------- #
# 3. Golden regression case                                                   #
# --------------------------------------------------------------------------- #


def _golden_case():
    """Fixed off-midspan kernel sweep exercising BOTH branches and the bridge.

    Do not edit without regenerating (and justifying) the golden.
    """
    chord, M, Ue = 0.13, 0.05, 17.0
    beta = np.sqrt(1.0 - M * M)
    freqs = np.logspace(np.log10(250.0), np.log10(16000.0), 12)
    cases = []
    for th in (20.0, 70.0, 115.0, 165.0):
        x, z = _obs(th)
        for xi_t in (0.0, 0.45, 0.93, 0.999, 1.02, 1.6, 5.0):
            cases.append((f"th{th:.0f}_xi{xi_t}", chord, M, Ue, x, z,
                          freqs, xi_t, beta))
    return cases


def _golden_values() -> dict:
    out = {}
    for lbl, chord, M, Ue, x, z, freqs, xi_t, beta in _golden_case():
        r0 = kernel(chord, M, Ue, x, z, freqs, np.zeros_like(freqs))
        K2b = xi_t * beta * r0["mu_bar"]
        for tag, br in (("raw", False), ("bridged", True)):
            r = kernel(chord, M, Ue, x, z, freqs, K2b, bridged=br)
            out[f"{lbl}|{tag}"] = r["I_abs"].tolist()
    return out


def create_golden() -> None:
    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
    GOLDEN_FILE.write_text(json.dumps(_golden_values(), indent=1))
    v = _golden_values()
    print(f"Golden file written to {GOLDEN_FILE} "
          f"({len(v)} sweeps, {sum(len(x) for x in v.values())} values)")
    print("Commit it, and record the justification in CHANGELOG.md.")


def test_golden() -> list:
    print(f"\n3. Golden regression case (rtol {GOLDEN_TOL:g}):")
    if not GOLDEN_FILE.exists():
        print(f"  [FAIL] golden file not found: {GOLDEN_FILE}")
        print("         run --create-golden first, on known-good code.")
        return [False]

    golden = json.loads(GOLDEN_FILE.read_text())
    current = _golden_values()
    if set(golden) != set(current):
        print("  [FAIL] golden case set differs from the current sweep set")
        return [False]

    worst, worst_key = 0.0, ""
    for k, g in golden.items():
        g = np.array(g)
        c = np.array(current[k])
        nz = np.abs(g) > 1e-14
        e = np.zeros_like(g)
        e[nz] = np.abs(c[nz] - g[nz]) / np.abs(g[nz])
        e[~nz] = np.abs(c[~nz])
        m = float(e.max())
        if m > worst:
            worst, worst_key = m, k
    n = sum(len(v) for v in golden.values())
    return [_report(f"|I| over {len(golden)} sweeps ({n} values)",
                    worst < GOLDEN_TOL,
                    f"max_rel_err={worst:.2e} (worst: {worst_key})")]


# --------------------------------------------------------------------------- #
# Entry point                                                                 #
# --------------------------------------------------------------------------- #


def run_tests() -> None:
    results = []
    results += test_algebra()
    results += test_cut()
    results += test_golden()
    print(f"\n{sum(results)}/{len(results)} checks passed.")
    sys.exit(0 if all(results) else 1)


def main() -> None:
    p = argparse.ArgumentParser(
        description="GFoil general Amiet kernel test suite",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=("Examples:\n"
                "  python3 tests/amiet_kernel_test.py --test\n"
                "  python3 tests/amiet_kernel_test.py --create-golden\n"),
    )
    p.add_argument("--test", action="store_true",
                   help="Run all checks against the golden file")
    p.add_argument("--create-golden", action="store_true",
                   help="Regenerate the golden file (requires justification)")
    a = p.parse_args()

    if a.create_golden:
        create_golden()
    elif a.test:
        run_tests()
    else:
        p.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
