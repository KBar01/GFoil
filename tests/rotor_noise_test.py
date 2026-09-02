#!/usr/bin/env python3
"""Test suite for GFoil.rotor_noise — Schlinker-Amiet rotor TE-noise wrapper.

Follows the conventions of tests/regression_test.py: a plain script (no pytest
dependency), PASS/FAIL lines, exit code 0 iff every check passes.

Deterministic and self-contained: no aero solve anywhere. Strips carry either
synthetic-but-plausible turbulent TE boundary-layer states or an analytic
custom wall-pressure spectrum, so nothing here depends on solver convergence.

Groups:
  1. Rotation matrices Rz/Ry against hand-computed values.
  2. Emission-time solve: zero-flow exactness, quadratic residual on random
     subsonic cases, and convergence to the paper's far-field Eq. 4.6.
  3. Doppler ratio: approach/recede/perpendicular limits and the positivity
     guards.
  4. Static anchor — the key regression gate. With Omega = 0 the wrapper must
     reduce to a bare noise_run call, bit-for-bit (rtol 1e-12).
  5. Azimuth convergence: 72 vs 144 samples agree to < 0.02 dB.
  6. Far-field spherical spreading: -6.02 +/- 0.1 dB per distance doubling.
  7. Blade count: B=4 vs B=1 differ by exactly 10*log10(4).
  8. kC = 1 crossover frequency.
  9. Custom-WPS path: source-frequency evaluation, the assembly (Doppler
     exponent +2, azimuth average, blade count), and the Ue trap -- noise_run's
     Amiet stage divides by U_c = 0.7*Ue even with a custom spectrum, so Ue = 0
     returns NaN rather than silence.
 10. Directivity ordering -- the phase-2 gate. The rotor axis must be LOUDER
     than the rotor plane on the paper's Table 2 wind-turbine element. This is
     the sign the phase-1 mid-span kernel got backwards; see the test docstring
     and CHANGELOG "General oblique-gust Amiet kernel".
 11. Golden regression case (tests/golden/rotor_noise_scalars.json).

Quick reference
---------------
    python3 tests/rotor_noise_test.py --test
    python3 tests/rotor_noise_test.py --create-golden   # see below

Regenerating the golden file requires JUSTIFICATION and a CHANGELOG entry, the
same rule as tests/golden/*_scalars.json (see tests/README.md, "Updating the
golden files"). The rotor wrapper is pure post-processing, so its golden values
can only move if the wrapper's own maths changes or the underlying Amiet/WPS
kernel changes. Either way that is a physics change, not a refactor: record why.
"""

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from GFoil import noise_run  # noqa: E402
from GFoil.rotor_noise import (  # noqa: E402
    RotorConfig, RotorStrip, rotor_noise_run, strip_relative_speed,
    _Rz, _Ry, _emission_time, _doppler_ratio, _kC1_frequency,
    _blade_frame_geometry, _TE_OFFSET_FRAC,
)

GOLDEN_DIR = Path(__file__).parent / "golden"
GOLDEN_FILE = GOLDEN_DIR / "rotor_noise_scalars.json"
GOLDEN_TOL = 1e-10


# --------------------------------------------------------------------------- #
# Shared synthetic fixtures                                                   #
# --------------------------------------------------------------------------- #

# Plausible turbulent trailing-edge BL state, ordering exactly as
# FwdResult.verbose_data.BL_top:
#   [theta, delta_star, tau_max, Ue, dpdx, tau_wall, delta99]
BL_SUCTION = np.array([1.80e-3, 2.60e-3, 5.2, 62.0, 1400.0, 3.1, 1.85e-2])
BL_PRESSURE = np.array([1.30e-3, 1.70e-3, 4.1, 58.0, -450.0, 3.6, 1.40e-2])


def _report(label: str, passed: bool, detail: str = "") -> bool:
    tag = "PASS" if passed else "FAIL"
    print(f"  [{tag}] {label:<44s}  {detail}")
    return passed


def _rel(a: float, b: float) -> float:
    return abs(a - b) / abs(b) if abs(b) > 1e-14 else abs(a - b)


# --------------------------------------------------------------------------- #
# 1. Rotation matrices                                                        #
# --------------------------------------------------------------------------- #


def test_rotations() -> list:
    print("\n1. Rotation matrices (paper Appendix B):")
    out = []

    out.append(_report("Rz(0) == I", np.allclose(_Rz(0.0), np.eye(3)), ""))
    out.append(_report("Ry(0) == I", np.allclose(_Ry(0.0), np.eye(3)), ""))

    # Rz(pi/2): x_hat -> y_hat, y_hat -> -x_hat, z_hat fixed.
    Rz90 = _Rz(0.5 * np.pi)
    ok = (np.allclose(Rz90 @ [1, 0, 0], [0, 1, 0], atol=1e-15)
          and np.allclose(Rz90 @ [0, 1, 0], [-1, 0, 0], atol=1e-15)
          and np.allclose(Rz90 @ [0, 0, 1], [0, 0, 1], atol=1e-15))
    out.append(_report("Rz(pi/2) maps x->y, y->-x, z->z", ok, ""))

    # Ry(pi/2) with this sign convention: x_hat -> z_hat, z_hat -> -x_hat.
    Ry90 = _Ry(0.5 * np.pi)
    ok = (np.allclose(Ry90 @ [1, 0, 0], [0, 0, 1], atol=1e-15)
          and np.allclose(Ry90 @ [0, 0, 1], [-1, 0, 0], atol=1e-15)
          and np.allclose(Ry90 @ [0, 1, 0], [0, 1, 0], atol=1e-15))
    out.append(_report("Ry(pi/2) maps x->z, z->-x, y->y", ok, ""))

    # Arbitrary angle, explicit hand-written entries.
    th = 0.37
    c, s = np.cos(th), np.sin(th)
    out.append(_report(
        "Rz(0.37) entries",
        np.allclose(_Rz(th), [[c, -s, 0], [s, c, 0], [0, 0, 1]], atol=1e-15), ""))
    out.append(_report(
        "Ry(0.37) entries",
        np.allclose(_Ry(th), [[c, 0, -s], [0, 1, 0], [s, 0, c]], atol=1e-15), ""))

    # Both must be proper rotations: orthogonal, det +1.
    for name, R in [("Rz", _Rz(1.1)), ("Ry", _Ry(-0.8))]:
        ok = (np.allclose(R @ R.T, np.eye(3), atol=1e-14)
              and abs(np.linalg.det(R) - 1.0) < 1e-14)
        out.append(_report(f"{name} orthogonal, det = +1", ok, ""))

    return out


# --------------------------------------------------------------------------- #
# 2. Emission-time solve                                                      #
# --------------------------------------------------------------------------- #


def test_emission_time() -> list:
    print("\n2. Emission-time solve:")
    out = []

    # (a) No flow -> pure geometric propagation, c0*Te = |xo - xe|.
    rng = np.random.default_rng(20130065)
    xo = np.array([12.0, -4.0, 7.5])
    xe = np.array([1.0, 0.5, 0.0])
    c0Te = _emission_time(xo, xe, np.zeros(3))
    err = _rel(float(c0Te), float(np.linalg.norm(xo - xe)))
    out.append(_report("M_FO = 0  =>  c0*Te = |xo - xe|", err < 1e-15,
                       f"rel_err={err:.3e}"))

    # (b) Quadratic residual |xo - xe - M_FO c0 Te| - c0 Te == 0 for random
    #     subsonic flow Mach vectors and random geometry.
    worst = 0.0
    for _ in range(200):
        xo = rng.uniform(-50, 50, 3)
        xe = rng.uniform(-3, 3, 3)
        # Random direction, random subsonic magnitude.
        m = rng.normal(size=3)
        m /= np.linalg.norm(m)
        M_FO = m * rng.uniform(0.0, 0.95)
        R = np.linalg.norm(xo - xe)
        c0Te = float(_emission_time(xo, xe, M_FO))
        resid = abs(np.linalg.norm(xo - xe - M_FO * c0Te) - c0Te)
        worst = max(worst, resid / R)
    out.append(_report("random subsonic: residual < 1e-12 * R", worst < 1e-12,
                       f"worst={worst:.3e} (200 cases)"))

    # (c) Far-field limit. With xe = 0 our solve must reproduce the paper's
    #     Eq. 4.6 exactly:
    #       Re = R(-M cosTheta + sqrt(1 - M^2 sin^2 Theta)) / (1 - M^2)
    #     with Theta the angle between M_FO and xo.
    def eq46(xo, M_FO):
        R = np.linalg.norm(xo)
        M = np.linalg.norm(M_FO)
        if M == 0.0:
            return R
        cosT = float(M_FO @ xo) / (M * R)
        sin2T = max(0.0, 1.0 - cosT * cosT)
        return R * (-M * cosT + np.sqrt(1.0 - M * M * sin2T)) / (1.0 - M * M)

    worst = 0.0
    for _ in range(100):
        xo = rng.uniform(-80, 80, 3)
        M_FO = np.array([0.0, 0.0, -rng.uniform(0.0, 0.9)])
        got = float(_emission_time(xo, np.zeros(3), M_FO))
        worst = max(worst, _rel(got, eq46(xo, M_FO)))
    out.append(_report("xe = 0 reproduces paper Eq. 4.6", worst < 1e-14,
                       f"worst_rel_err={worst:.3e}"))

    # ... and with xe finite it converges to Eq. 4.6 as |xo| >> |xe|.
    xe = np.array([1.0, 0.0, 0.0])
    M_FO = np.array([0.0, 0.0, -0.3])
    errs = []
    for scale in (1e2, 1e3, 1e4, 1e5):
        xo = np.array([0.0, 0.6, 0.8]) * scale
        errs.append(_rel(float(_emission_time(xo, xe, M_FO)), eq46(xo, M_FO)))
    # Monotone decay, ~1/|xo|: each decade of distance buys a decade of accuracy.
    ok = all(errs[i + 1] < 0.2 * errs[i] for i in range(len(errs) - 1)) \
        and errs[-1] < 1e-4
    out.append(_report("finite xe -> Eq. 4.6 as |xo|/r -> inf", ok,
                       f"rel_err={[f'{e:.2e}' for e in errs]}"))

    return out


# --------------------------------------------------------------------------- #
# 3. Doppler ratio                                                            #
# --------------------------------------------------------------------------- #


def test_doppler() -> list:
    print("\n3. Doppler ratio (paper Eq. 4.10):")
    out = []

    CO = np.array([0.0, 0.0, 1.0])
    M_FO0 = np.zeros(3)
    M = 0.25

    # NOTE ON THE EXPECTED VALUES.
    # Eq. 4.10 is  doppler = 1 + (M_BO.CO)/(1 + (M_FO - M_BO).CO), which for
    # M_FO = 0 and a head-on approach collapses to
    #     1 + M/(1 - M) = 1/(1 - M)
    # -- the textbook moving-source factor, NOT 1 + M. (1 + M is the first-order
    # expansion of it: 1/(1-M) = 1 + M + M^2 + ...) The implementation prompt's
    # sanity note claimed 1 +/- |M_BO|; that contradicts the very equation it
    # quotes, so we assert the equation. See the final report / CHANGELOG.

    got = float(_doppler_ratio(M * CO, M_FO0, CO))
    out.append(_report("M_FO=0, approaching  =>  1/(1 - M)",
                       _rel(got, 1.0 / (1.0 - M)) < 1e-15,
                       f"got={got:.8f} expected={1.0 / (1.0 - M):.8f}"))

    got = float(_doppler_ratio(-M * CO, M_FO0, CO))
    out.append(_report("M_FO=0, receding     =>  1/(1 + M)",
                       _rel(got, 1.0 / (1.0 + M)) < 1e-15,
                       f"got={got:.8f} expected={1.0 / (1.0 + M):.8f}"))

    M_perp = np.array([M, 0.0, 0.0])   # perpendicular to CO
    got = float(_doppler_ratio(M_perp, M_FO0, CO))
    out.append(_report("M_BO perpendicular to CO  =>  1", _rel(got, 1.0) < 1e-15,
                       f"got={got:.16f}"))

    # A static source in a moving medium must show no Doppler shift.
    got = float(_doppler_ratio(np.zeros(3), np.array([0.0, 0.0, -0.3]), CO))
    out.append(_report("M_BO = 0 (static source in flow)  =>  1",
                       _rel(got, 1.0) < 1e-15, f"got={got:.16f}"))

    # Equivalent closed form: (1 + M_FO.CO)/(1 + (M_FO - M_BO).CO).
    rng = np.random.default_rng(4711)
    worst = 0.0
    for _ in range(200):
        v = rng.normal(size=3)
        co = v / np.linalg.norm(v)
        M_BO = rng.normal(size=3) * 0.1
        M_FO = np.array([0.0, 0.0, -rng.uniform(0.0, 0.4)])
        got = float(_doppler_ratio(M_BO, M_FO, co))
        expect = (1.0 + M_FO @ co) / (1.0 + (M_FO - M_BO) @ co)
        worst = max(worst, _rel(got, expect))
    out.append(_report("matches (1+M_FO.CO)/(1+(M_FO-M_BO).CO)", worst < 1e-14,
                       f"worst_rel_err={worst:.3e}"))

    # Positivity guards must fire for contrived sonic / supersonic input.
    fired = 0
    for M_BO in (1.0 * CO, 1.2 * CO):     # denominator -> 0, then negative
        try:
            _doppler_ratio(M_BO, M_FO0, CO)
        except ValueError:
            fired += 1
    out.append(_report("guard fires at sonic and supersonic M_BO.CO",
                       fired == 2, f"{fired}/2 raised ValueError"))

    return out


# --------------------------------------------------------------------------- #
# 4. Static anchor — the key regression gate                                  #
# --------------------------------------------------------------------------- #


def test_static_anchor() -> list:
    r"""With Omega = 0 the wrapper must reduce to a bare noise_run call.

    Configuration: one strip, Omega = 0, Uz > 0, B = 1, n_azimuth = 1,
    pitch_rad = 0.

    Derivation of the observer placement.
      n_azimuth = 1  =>  the single azimuth is  gamma_0 = 2*pi*0/1 = 0.
      (The prompt said gamma = pi/2, but gamma_j = 2*pi*j/N_gamma gives 0 for
      N_gamma = 1; we use the grid the module actually generates. The chain is
      then a rotation composed with a translation rather than a pure
      translation, which is inverted analytically below -- exactly as sound.)

      Omega = 0  =>  M_BO = 0  =>  xp = xe  and  doppler = 1 + 0/denom = 1
      exactly, for any Uz. So the source frequencies equal the observer
      frequencies and the weight 1/doppler^2 is exactly 1.

      xe   = (r, 0, 0)                                       (gamma = 0)
      x1   = xo - xp = xo - xe
      Rz(pi/2 - 0) = Rz(pi/2) = [[0,-1,0],[1,0,0],[0,0,1]]
      x2   = Rz(pi/2) x1 = (-x1_y,  x1_x,  x1_z)
      X    = Ry(0) x2 = x2                                   (pitch = 0)

      So to land on a chosen blade-frame observer X = (X_c, Y_s, Z_n):
          -x1_y = X_c  and  x1_x = Y_s  and  x1_z = Z_n
      =>   x1 = (Y_s, -X_c, Z_n)
      =>   xo = xe + x1 = (r + Y_s, -X_c, Z_n)

      The wrapper then calls noise_run at observer
      (X_c + 0.75*chord, Y_s, Z_n) with alphaDeg = 0, which the kernel maps to
      TE-local (X_c, Y_s, Z_n) -- so a direct call at the same position must
      return the identical spectrum, and Spp = B * (1/1) * 1 * FF = FF.
    """
    print("\n4. Static anchor (Omega = 0: wrapper reduces to noise_run):")
    out = []

    radius, dr, chord = 2.5, 0.4, 0.18
    Uz, nu, rho = 45.0, 1.48e-5, 1.225
    model = "kam"

    X_c, Y_s, Z_n = 0.9, 0.35, 1.6      # target blade-frame observer
    xo = np.array([radius + Y_s, -X_c, Z_n])

    freqs = np.logspace(np.log10(400.0), np.log10(12000.0), 48)

    cfg = RotorConfig(Omega=0.0, Uz=Uz, B=1, rho=rho, nu=nu, n_azimuth=1)
    strip = RotorStrip(radius=radius, dr=dr, chord=chord, pitch_rad=0.0,
                       BL_top=BL_SUCTION, BL_bot=BL_PRESSURE, wps_model=model)
    res = rotor_noise_run(cfg, [strip], freqs, xo)

    # doppler must be exactly 1, so the source grid is the observer grid.
    ok = (res.diagnostics["doppler_min"] == 1.0
          and res.diagnostics["doppler_max"] == 1.0)
    out.append(_report("doppler identically 1", ok,
                       f"[{res.diagnostics['doppler_min']!r}, "
                       f"{res.diagnostics['doppler_max']!r}]"))

    U_rel = strip_relative_speed(cfg, radius)
    out.append(_report("U_rel == Uz when Omega = 0", U_rel == Uz,
                       f"U_rel={U_rel!r}"))

    direct = noise_run(
        BL_SUCTION, BL_PRESSURE,
        freqs_Hz=freqs,
        observerXYZ=np.array([X_c + _TE_OFFSET_FRAC * chord, Y_s, Z_n]),
        Re=U_rel * chord / nu, nu=nu, chord=chord, span=dr,
        alphaDeg=0.0, rho=rho, model=model,
    )

    # The kernel really did receive the paper-frame observer.
    ok = np.allclose(direct.obsXYZ_TElocal[0], [X_c, Y_s, Z_n], rtol=0, atol=1e-14)
    out.append(_report("kernel TE-local frame == (X_c, Y_s, Z_n)", ok,
                       f"got={direct.obsXYZ_TElocal[0]}"))

    got = res.Spp[0]
    want = direct.FF_spectra[0]
    ok = np.allclose(got, want, rtol=1e-12, atol=0.0)
    max_err = float(np.max(np.abs(got - want) / np.abs(want)))
    out.append(_report("Spp == direct noise_run FF (rtol 1e-12)", ok,
                       f"max_rel_err={max_err:.3e}"))

    return out


# --------------------------------------------------------------------------- #
# 5-7. Rotor-level behaviour                                                  #
# --------------------------------------------------------------------------- #


def _demo_case(n_azimuth: int = 72, B: int = 1, scale: float = 1.0):
    """A representative rotating case; observer distance scaled by `scale`."""
    cfg = RotorConfig(Omega=120.0, Uz=15.0, B=B, n_azimuth=n_azimuth)
    strips = [RotorStrip(radius=1.6, dr=0.5, chord=0.14, pitch_rad=0.12,
                         BL_top=BL_SUCTION, BL_bot=BL_PRESSURE, wps_model="kam")]
    freqs = np.logspace(np.log10(600.0), np.log10(12000.0), 32)
    xo = np.array([0.0, 60.0, 45.0]) * scale
    return cfg, strips, freqs, xo


def test_azimuth_convergence() -> list:
    print("\n5. Azimuth convergence (72 vs 144 samples):")
    o72 = rotor_noise_run(*_demo_case(n_azimuth=72)).OASPL_perObs[0]
    o144 = rotor_noise_run(*_demo_case(n_azimuth=144)).OASPL_perObs[0]
    d = abs(o72 - o144)
    return [_report("|OASPL(72) - OASPL(144)| < 0.02 dB", d < 0.02,
                    f"72={o72:.5f} dB  144={o144:.5f} dB  diff={d:.2e} dB")]


def test_spherical_spreading() -> list:
    print("\n6. Far-field spherical spreading (distance doubling):")
    near = rotor_noise_run(*_demo_case(scale=1.0)).OASPL_perObs[0]
    far = rotor_noise_run(*_demo_case(scale=2.0)).OASPL_perObs[0]
    d = far - near
    return [_report("dOASPL = -6.02 +/- 0.1 dB", abs(d - (-6.02)) < 0.1,
                    f"near={near:.4f} dB  far={far:.4f} dB  delta={d:.4f} dB")]


def test_blade_count() -> list:
    print("\n7. Blade count (incoherent power sum):")
    o1 = rotor_noise_run(*_demo_case(B=1)).OASPL_perObs[0]
    o4 = rotor_noise_run(*_demo_case(B=4)).OASPL_perObs[0]
    expect = 10.0 * np.log10(4.0)
    err = abs((o4 - o1) - expect)
    return [_report("OASPL(B=4) - OASPL(B=1) == 10*log10(4)", err < 1e-10,
                    f"delta={o4 - o1:.12f} dB  expected={expect:.12f} dB  "
                    f"abs_err={err:.2e}")]


# --------------------------------------------------------------------------- #
# 8. kC = 1 crossover                                                         #
# --------------------------------------------------------------------------- #


def test_kC1() -> list:
    print("\n8. kC = 1 crossover frequency:")
    out = []
    c0, chord = 340.0, 0.1

    f = _kC1_frequency(c0, chord)
    out.append(_report("f_kC1 = c0/(2*pi*C)", _rel(f, c0 / (2 * np.pi * chord)) < 1e-15,
                       f"f={f:.6f} Hz"))

    # The defining property: kC = (omega/c0)*C == 1 at that frequency.
    kC = (2.0 * np.pi * f / c0) * chord
    out.append(_report("kC(f_kC1) == 1", _rel(kC, 1.0) < 1e-15, f"kC={kC:.16f}"))

    # It is NOT c0/(pi*C) -- that expression is a factor 2 out (kC would be 2).
    out.append(_report("distinct from the c0/(pi*C) form",
                       _rel(f, c0 / (np.pi * chord)) > 0.49,
                       f"c0/(pi*C)={c0 / (np.pi * chord):.4f} Hz gives kC=2"))

    # Band fraction below kC = 1 is reported and warns above 20%.
    cfg = RotorConfig(Omega=100.0, Uz=10.0, B=1, n_azimuth=4)
    strips = [RotorStrip(radius=1.0, dr=0.2, chord=chord, pitch_rad=0.0,
                         BL_top=BL_SUCTION, BL_bot=BL_PRESSURE)]
    freqs = np.linspace(200.0, 10000.0, 50)
    d = rotor_noise_run(cfg, strips, freqs, np.array([0.0, 30.0, 20.0])) \
        .diagnostics["strips"][0]
    expect = float(np.mean(freqs < f))
    out.append(_report("frac_band_below_kC1 reported correctly",
                       _rel(d["frac_band_below_kC1"], expect) < 1e-15,
                       f"frac={d['frac_band_below_kC1']:.4f}  f_kC1={d['f_kC1_Hz']:.2f} Hz"))
    return out


# --------------------------------------------------------------------------- #
# 9. Custom-WPS path, source frequencies, and the Ue trap                     #
# --------------------------------------------------------------------------- #


def test_custom_wps_path() -> list:
    print("\n9. Custom-WPS path (source frequencies, Ue, assembly):")
    out = []

    cfg = RotorConfig(Omega=110.0, Uz=14.0, B=2, n_azimuth=4)
    freqs = np.logspace(np.log10(800.0), np.log10(6000.0), 8)
    omega = 2.0 * np.pi * freqs
    xo = np.array([0.0, 30.0, 22.0])

    # (a) custom_WPS_func must be evaluated at the DOPPLER-SHIFTED SOURCE
    #     frequencies, never the observer grid. The wall-pressure spectrum is a
    #     function of source frequency; getting this wrong evaluates the WPS at
    #     the wrong frequencies (prompt section 2.5).
    seen = []

    def spy(omega_src):
        seen.append(np.array(omega_src, copy=True))
        return np.full_like(omega_src, 0.05), np.full_like(omega_src, 0.02)

    strip = RotorStrip(radius=1.4, dr=0.3, chord=0.12, pitch_rad=0.0,
                       custom_WPS_func=spy, Ue_custom=(55.0, 52.0))
    res = rotor_noise_run(cfg, [strip], freqs, xo)

    _X, doppler = _blade_frame_geometry(cfg, strip, xo)
    out.append(_report("custom_WPS_func called once per azimuth",
                       len(seen) == cfg.n_azimuth,
                       f"{len(seen)} calls, expected {cfg.n_azimuth}"))

    worst = max(float(np.max(np.abs(seen[j] - omega / doppler[j])
                             / (omega / doppler[j])))
                for j in range(cfg.n_azimuth))
    out.append(_report("evaluated at omega/doppler (source freqs)",
                       worst < 1e-15, f"max_rel_err={worst:.2e}"))

    # And those grids really are shifted -- otherwise the check above is vacuous.
    shifted = max(float(np.max(np.abs(seen[j] / omega - 1.0)))
                  for j in range(cfg.n_azimuth))
    out.append(_report("source grids differ from the observer grid",
                       shifted > 1e-3,
                       f"max|omega_src/omega - 1| = {shifted:.4f}"))

    # (b) Assembly: independently rebuild Spp from per-azimuth noise_run calls.
    #     Pins the Doppler exponent (+2), the 1/N_gamma average, the factor B,
    #     the 0.75*chord offset and Re = U_rel*chord/nu all at once.
    U_rel = strip_relative_speed(cfg, strip.radius)
    BL_t, BL_b = np.zeros(7), np.zeros(7)
    BL_t[3], BL_b[3] = 55.0, 52.0
    expect = np.zeros(freqs.size)
    for j in range(cfg.n_azimuth):
        omega_src = omega / doppler[j]
        wu, wl = spy(omega_src)
        nr = noise_run(
            BL_t, BL_b, freqs_Hz=omega_src / (2.0 * np.pi),
            observerXYZ=np.array([_X[j, 0] + _TE_OFFSET_FRAC * strip.chord,
                                  _X[j, 1], _X[j, 2]]),
            Re=U_rel * strip.chord / cfg.nu, nu=cfg.nu, chord=strip.chord,
            span=strip.dr, alphaDeg=0.0, rho=cfg.rho,
            custom_WPS=np.column_stack([wu, wl]),
        )
        expect += nr.FF_spectra[0] / doppler[j] ** 2      # exponent +2
    expect *= cfg.B / cfg.n_azimuth

    err = float(np.max(np.abs(res.Spp[0] - expect) / np.abs(expect)))
    out.append(_report("Spp == B/Ng * sum FF/doppler^2 (rtol 1e-12)",
                       err < 1e-12, f"max_rel_err={err:.2e}"))

    # (c) The Ue trap. noise_run's Amiet stage divides by U_c = 0.7*Ue even when
    #     a custom WPS is supplied, so Ue = 0 yields a NaN far field -- which
    #     np.where(Spp > 0) would silently floor to -200 dB and read as "quiet".
    #     Ue_custom is validated up front; a zero-Ue BL state must raise.
    try:
        RotorStrip(radius=1.0, dr=0.2, chord=0.1, pitch_rad=0.0,
                   custom_WPS_func=spy, Ue_custom=(0.0, 50.0))
        ok = False
    except ValueError:
        ok = True
    out.append(_report("Ue_custom = 0 rejected at construction", ok, ""))

    BL_noue = BL_SUCTION.copy()
    BL_noue[3] = 0.0          # Ue = 0 -> kernel returns NaN
    bad = RotorStrip(radius=1.0, dr=0.2, chord=0.1, pitch_rad=0.0,
                     BL_top=BL_noue, BL_bot=BL_noue)
    try:
        rotor_noise_run(cfg, [bad], freqs, xo)
        ok, detail = False, "no error raised — NaN would reach the dB floor"
    except RuntimeError as e:
        ok = "non-finite" in str(e)
        detail = "RuntimeError raised"
    out.append(_report("NaN far field raises, never floors to -200 dB",
                       ok, detail))

    # (d) Ue_custom defaults to U_rel, matching the kernel's own fallback for a
    #     skipped surface (noise_run.hpp: Ue_top = upper_skip ? Uinf : Ue[0]).
    s_default = RotorStrip(radius=1.4, dr=0.3, chord=0.12, pitch_rad=0.0,
                           custom_WPS_func=spy)
    s_explicit = RotorStrip(radius=1.4, dr=0.3, chord=0.12, pitch_rad=0.0,
                            custom_WPS_func=spy,
                            Ue_custom=(U_rel, U_rel))
    a = rotor_noise_run(cfg, [s_default], freqs, xo).Spp
    b = rotor_noise_run(cfg, [s_explicit], freqs, xo).Spp
    out.append(_report("Ue_custom defaults to U_rel",
                       np.array_equal(a, b), f"U_rel={U_rel:.3f} m/s"))

    # (e) Mutual exclusion of the two source-spectrum paths.
    try:
        RotorStrip(radius=1.0, dr=0.2, chord=0.1, pitch_rad=0.0,
                   BL_top=BL_SUCTION, BL_bot=BL_PRESSURE, custom_WPS_func=spy)
        ok = False
    except ValueError:
        ok = True
    out.append(_report("BL states + custom_WPS_func rejected", ok, ""))

    return out


# --------------------------------------------------------------------------- #
# 10. Directivity ordering — the phase-2 headline gate                        #
# --------------------------------------------------------------------------- #


def _chou_george_delta_star(chord: float, chi_deg: float) -> float:
    """Displacement thickness [m] — Sinayoko et al. (2013) Eq. 3.10."""
    if chi_deg <= 4.0:
        return chord * (24.3 + 0.6625 * chi_deg) * 1e-4
    d = chi_deg - 4.0
    return chord * (26.95 + 0.6625 * d + 0.3044 * d ** 2 + 0.0104 * d ** 3) * 1e-4


def _chou_george_Sqq(omega, rho, U_X, delta_star):
    """Wall-pressure spectrum [Pa^2/(rad/s)] — Sinayoko et al. (2013) Eqs. 3.9/3.11."""
    w = np.asarray(omega * delta_star / U_X, dtype=float)
    lo = 1.732e-3 * w / (1.0 - 5.489 * w + 36.74 * w ** 2 + 0.1505 * w ** 5)
    hi = 1.4216e-3 * w / (0.3261 + 4.1837 * w + 22.818 * w ** 2
                          + 0.0013 * w ** 3 + 0.0028 * w ** 5)
    F = np.where(w < 0.06, lo, hi)
    return (0.5 * rho * U_X ** 2) ** 2 * (delta_star / U_X) * F


def test_directivity_ordering() -> list:
    r"""The rotor axis must be LOUDER than the rotor plane.

    This is the gate on the phase-2 headline fix, and the one thing phase 1 got
    backwards. The trailing-edge source is a compact edge dipole with its axis
    along the plate normal; at the Sinayoko et al. (2013) Table 2 wind-turbine
    element the chord line lies close to the rotor plane (chi = 10 deg), so the
    normal is close to the rotor AXIS. The directivity must therefore peak near
    Theta = 0 / 180 deg and null in the rotor plane (Theta = 90 deg).

    The phase-1 mid-span kernel inverted this: it formed
    S0 = sqrt(x_loc^2 + beta^2 z_loc^2), never reading y_loc, so at the azimuths
    where an observer in a plane containing the rotor axis crosses a blade's
    spanwise direction it discarded nearly the whole source-observer separation
    (a 1000 m observer treated as ~5 m). Those azimuths dominated the average
    and buried the null, leaving the rotor plane ~21 dB LOUDER than the axis.

    Coarse on purpose: a > 3 dB contrast, not a level match. The point is to pin
    the SIGN of the ordering against regression forever, and the sign is what
    changed. (Measured at the time of writing: +12.5 dB with the general kernel,
    against -21 dB with the mid-span one.)
    """
    print("\n10. Directivity ordering (axis louder than rotor plane):")
    out = []

    C0 = 340.0
    radius, chord = 21.75, 2.0
    M_BO, M_FO, chi_deg = 0.165, 0.029, 10.0
    span = radius / 3.0
    Omega = M_BO * C0 / radius
    Uz = M_FO * C0
    f_hz = (5.0 * C0 / chord) / (2.0 * np.pi)      # kC = 5

    cfg = RotorConfig(Omega=Omega, Uz=Uz, B=1, rho=1.225, nu=1.48e-5,
                      c0=C0, n_azimuth=72)
    U_X = strip_relative_speed(cfg, radius)
    dstar = _chou_george_delta_star(chord, chi_deg)

    def custom_WPS(omega_src):
        upper = _chou_george_Sqq(omega_src, cfg.rho, U_X, dstar)
        return upper, np.zeros_like(upper)

    strip = RotorStrip(radius=radius, dr=span, chord=chord,
                       pitch_rad=np.deg2rad(chi_deg),
                       custom_WPS_func=custom_WPS, Ue_custom=(U_X, U_X))

    freqs = np.array([0.99 * f_hz, f_hz, 1.01 * f_hz])
    R = 1000.0
    # Theta measured from the rotor axis (+z), observers in the y-z plane.
    thetas = [0.0, 90.0, 180.0]
    obs = np.array([[0.0, R * np.sin(np.deg2rad(t)), R * np.cos(np.deg2rad(t))]
                    for t in thetas])

    res = rotor_noise_run(cfg, [strip], freqs, obs)
    spl = 10.0 * np.log10(2.0 * np.pi * res.Spp[:, 1] * R ** 2 / (20e-6) ** 2)
    fore, plane, aft = spl[0], spl[1], spl[2]
    contrast = max(fore, aft) - plane

    out.append(_report("SPL finite at all three observers",
                       bool(np.all(np.isfinite(spl))),
                       f"axis0={fore:.2f}  plane={plane:.2f}  axis180={aft:.2f} dB/Hz"))
    out.append(_report("axis louder than rotor plane by > 3 dB",
                       contrast > 3.0,
                       f"contrast = {contrast:.2f} dB "
                       f"(phase-1 mid-span kernel gave about -21 dB)"))
    out.append(_report("both axial lobes exceed the rotor plane",
                       fore - plane > 3.0 and aft - plane > 3.0,
                       f"fore-plane={fore - plane:.2f} dB  "
                       f"aft-plane={aft - plane:.2f} dB"))
    return out


# --------------------------------------------------------------------------- #
# 11. Golden regression case                                                  #
# --------------------------------------------------------------------------- #

# Pinned (observer, frequency) indices sampled from Spp for the golden file.
GOLDEN_PINS = [(0, 0), (0, 7), (0, 15), (1, 3), (1, 11), (2, 9), (2, 19)]


def _golden_case():
    """Fixed multi-strip rotor case. Do not edit without regenerating."""
    cfg = RotorConfig(Omega=95.0, Uz=12.0, B=3, rho=1.225, nu=1.48e-5,
                      c0=340.0, n_azimuth=36)
    strips = [
        RotorStrip(radius=1.20, dr=0.40, chord=0.160, pitch_rad=0.20,
                   BL_top=BL_SUCTION, BL_bot=BL_PRESSURE, wps_model="kam"),
        RotorStrip(radius=1.60, dr=0.40, chord=0.140, pitch_rad=0.15,
                   BL_top=BL_SUCTION, BL_bot=BL_PRESSURE, wps_model="kam"),
        RotorStrip(radius=2.00, dr=0.40, chord=0.120, pitch_rad=0.10,
                   BL_top=1.15 * BL_SUCTION, BL_bot=0.9 * BL_PRESSURE,
                   wps_model="roz"),
    ]
    freqs = np.logspace(np.log10(500.0), np.log10(15000.0), 20)
    obs = np.array([[0.0, 40.0, 30.0],
                    [0.0, 0.0, 50.0],
                    [25.0, 25.0, 35.0]])
    return cfg, strips, freqs, obs


def _golden_values() -> dict:
    res = rotor_noise_run(*_golden_case())
    return {
        "OASPL_perObs": res.OASPL_perObs.tolist(),
        "Spp_pinned": {f"{i},{j}": float(res.Spp[i, j]) for i, j in GOLDEN_PINS},
    }


def create_golden() -> None:
    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
    GOLDEN_FILE.write_text(json.dumps(_golden_values(), indent=2))
    print(f"Golden file written to {GOLDEN_FILE}")
    print("Commit it, and record the justification in CHANGELOG.md.")


def test_golden() -> list:
    print(f"\n11. Golden regression case (rtol {GOLDEN_TOL:g}):")
    if not GOLDEN_FILE.exists():
        print(f"  [FAIL] golden file not found: {GOLDEN_FILE}")
        print("         run --create-golden first, on known-good code.")
        return [False]

    golden = json.loads(GOLDEN_FILE.read_text())
    current = _golden_values()
    out = []

    for i, (c, g) in enumerate(zip(current["OASPL_perObs"],
                                   golden["OASPL_perObs"])):
        e = _rel(c, g)
        out.append(_report(f"OASPL obs {i}", e < GOLDEN_TOL,
                           f"current={c:.10f}  golden={g:.10f}  rel_err={e:.2e}"))

    worst, worst_key = 0.0, ""
    for k, g in golden["Spp_pinned"].items():
        e = _rel(current["Spp_pinned"][k], g)
        if e > worst:
            worst, worst_key = e, k
    out.append(_report(f"Spp at {len(GOLDEN_PINS)} pinned (obs,freq)",
                       worst < GOLDEN_TOL,
                       f"max_rel_err={worst:.2e} (worst: [{worst_key}])"))
    return out


# --------------------------------------------------------------------------- #
# Entry point                                                                 #
# --------------------------------------------------------------------------- #


def run_tests() -> None:
    # Several cases deliberately sit in the mid-span kernel's advisory-warning
    # territory (spanwise neglect, kC, l_S). Those warnings are the module
    # working as intended and would drown the PASS/FAIL lines; every check here
    # asserts numbers, never warning text.
    warnings.simplefilter("ignore", UserWarning)

    results = []
    results += test_rotations()
    results += test_emission_time()
    results += test_doppler()
    results += test_static_anchor()
    results += test_azimuth_convergence()
    results += test_spherical_spreading()
    results += test_blade_count()
    results += test_kC1()
    results += test_custom_wps_path()
    results += test_directivity_ordering()
    results += test_golden()

    print(f"\n{sum(results)}/{len(results)} checks passed.")
    sys.exit(0 if all(results) else 1)


def main() -> None:
    p = argparse.ArgumentParser(
        description="GFoil rotor-noise test suite",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=("Examples:\n"
                "  python3 tests/rotor_noise_test.py --test\n"
                "  python3 tests/rotor_noise_test.py --create-golden\n"),
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
