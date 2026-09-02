#!/usr/bin/env python3
"""Qualitative paper check for GFoil.rotor_noise — NOT a gated test.

Reproduces the observer geometry of Fig. 9(b) of

    S. Sinayoko, M. Kingan & A. Agarwal (2013),
    "Trailing edge noise theory for rotating blades in uniform flow",
    Proc. R. Soc. A 469:20130065

for their Table 2 **wind turbine** blade element, driving GFoil's rotor wrapper
with the Chou & George (1984) wall-pressure model the paper itself uses
(their Eqs. 3.9-3.11), supplied through ``RotorStrip.custom_WPS_func``.


THIS IS A SHAPE CHECK ONLY — exact agreement with the paper is NOT expected
--------------------------------------------------------------------------
The two calculations still do not solve quite the same scattering problem.
GFoil's kernel now uses the general Roger & Moreau (2005) oblique-gust
formulation, but in the Eq. 18 large-aspect-ratio limit, where the spanwise
wavenumber integral collapses to a delta selecting the single gust
K2 = k*x2/S0 — rather than the paper's finite-span sinc over many gusts.
Absolute levels also depend on the assumed one-sidedness of Chou & George's
Sqq (see ``chou_george_Sqq``). The directivity SHAPE is what this checks.


HISTORY: phase 1 produced a quantified NEGATIVE result — an inverted pattern
----------------------------------------------------------------------------
Kept on the record, because it is why the phase-2 kernel exists.

Expected from the physics (what the paper shows): the TE source is a compact
edge dipole with its axis along the plate normal. At this element the chord line
lies close to the rotor plane (chi = 10 deg), so the plate normal is close to
the rotor axis. The directivity should therefore have **maxima near the rotor
axis** (Theta ~ 0, 180 deg) and a **null in the rotor plane** (Theta ~ 90, 270
deg), with a modest fore/aft asymmetry from the axial convection
(M_FO = -0.029 z_hat) and the blade-motion Doppler.

What phase 1 produced: the rotor plane ~21 dB **louder** than the axis — the
pattern inverted, not merely inaccurate. The mid-span kernel formed

    S0 = sqrt(x_loc^2 + beta^2 z_loc^2)          # y_loc never appeared

so the distance it used was sqrt(X_c^2 + Z_n^2), not |X|. An observer in any
plane containing the rotor axis crosses each blade's spanwise direction twice
per revolution, and at those azimuths almost the whole source-observer
separation sits in Y_s and was silently discarded: at its worst azimuth the
kernel placed this 1000 m observer at ~5 m, a 231x under-estimate, inflating
that azimuth's PSD (~1/S0^4) by ~95 dB. Those few azimuths dominated the
azimuthal average and buried the true dipole null.


RESOLVED IN PHASE 2 (July 2026) — see CHANGELOG, "General oblique-gust Amiet
kernel (phase 2)"
----------------------------------------------------------------------------
``TE_noise_outer_vec`` now forms the general S0 = sqrt(x1^2 + beta^2(x2^2+x3^2)),
selects the Eq. 18 gust K2_bar = k_bar*x2/S0, branches on the criticality
parameter xi = beta*|x2|/S0 (supercritical / subcritical, with a near-cutoff
regularisation), and uses the spanwise-wavenumber-corrected Corcos length.

Measured on this case, before -> after:

    Theta = 0   (upstream axis)   64.58 -> 64.58  dB/Hz    (unchanged)
    Theta = 180 (downstream axis) 65.18 -> 65.17  dB/Hz    (unchanged)
    Theta = 90  (rotor plane)     86.43 -> 52.62  dB/Hz    (-33.8 dB)
    axis-to-plane contrast       -21.25 -> +12.55 dB       (sign corrected)

The axial lobes barely move: near the rotor axis the observer never approaches a
blade's spanwise direction, so the mid-span kernel was already right there. The
whole 33.8 dB correction lands in the rotor plane, exactly where the phase-1
analysis said it would. The ordering is now gated as a test
(tests/rotor_noise_test.py, group 10).

``spanwise_neglect_ratio`` is still reported, but it no longer measures an
error — it is now an obliquity indicator and a proxy for xi. See
``GFoil/rotor_noise.py``.

Run:
    python3 bench/rotor_noise_paper_check.py
Writes:
    bench/results/rotor_noise_directivity.png   polar directivity
    bench/results/rotor_noise_cut.png           unbridged vs bridged |I(xi)|
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")          # headless
import matplotlib.pyplot as plt  # noqa: E402

REPO_ROOT = Path(__file__).parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from GFoil.rotor_noise import (  # noqa: E402
    RotorConfig, RotorStrip, rotor_noise_run, strip_relative_speed,
)

OUT_PNG = Path(__file__).parent / "results" / "rotor_noise_directivity.png"

P_REF = 20e-6

# --------------------------------------------------------------------------- #
# Paper Table 2 — wind turbine blade element                                  #
# --------------------------------------------------------------------------- #

C0 = 340.0            # speed of sound [m/s] (also what the kernel hard-codes)
RADIUS = 21.75        # blade element radius [m]
CHORD = 2.0           # chord [m]
M_BO = 0.165          # blade Mach number (rotational)
M_FO = 0.029          # flow (axial inflow) Mach number
CHI_DEG = 10.0        # section angle of attack [deg]
SPAN = RADIUS / 3.0   # spanwise extent of the element [m]

OMEGA_ROTOR = M_BO * C0 / RADIUS   # [rad/s]
UZ = M_FO * C0                     # axial inflow speed [m/s]

KC = 5.0                           # target acoustic chord parameter
# kC = (omega/c0)*C  =>  omega = kC*c0/C
OMEGA_ACOUSTIC = KC * C0 / CHORD
F_HZ = OMEGA_ACOUSTIC / (2.0 * np.pi)

R_OBS = 1000.0        # observer radius [m]; >> RADIUS and >> wavelength
N_THETA = 181


# --------------------------------------------------------------------------- #
# Chou & George (1984) wall-pressure spectrum — paper Eqs. 3.9-3.11           #
# --------------------------------------------------------------------------- #


def chou_george_delta_star(chord: float, chi_deg: float) -> float:
    """Displacement thickness [m] (paper Eq. 3.10).

    A pure empirical fit in the section angle of attack — no BL solve.
    """
    if chi_deg <= 4.0:
        return chord * (24.3 + 0.6625 * chi_deg) * 1e-4
    d = chi_deg - 4.0
    return chord * (26.95 + 0.6625 * d + 0.3044 * d ** 2
                    + 0.0104 * d ** 3) * 1e-4


def chou_george_F(omega_bar: np.ndarray) -> np.ndarray:
    """Non-dimensional spectral shape F(omega_bar) (paper Eq. 3.11).

    Piecewise at omega_bar = 0.06.
    """
    w = np.asarray(omega_bar, dtype=float)
    lo = 1.732e-3 * w / (1.0 - 5.489 * w + 36.74 * w ** 2 + 0.1505 * w ** 5)
    hi = 1.4216e-3 * w / (0.3261 + 4.1837 * w + 22.818 * w ** 2
                          + 0.0013 * w ** 3 + 0.0028 * w ** 5)
    return np.where(w < 0.06, lo, hi)


def chou_george_Sqq(omega: np.ndarray, rho: float, U_X: float,
                    delta_star: float) -> np.ndarray:
    """Wall-pressure spectrum Sqq(omega) (paper Eq. 3.9).

        Sqq = (0.5 rho U_X^2)^2 * (delta*/U_X) * F(omega_bar),
        omega_bar = omega delta*/U_X

    UNITS. The leading factor is a dynamic pressure squared [Pa^2] and
    delta*/U_X is a time [s], so Sqq is Pa^2 * s == Pa^2/(rad/s). That is
    already the convention ``noise_run``'s ``custom_WPS`` argument expects
    ("raw linear Pa^2/omega", NoiseResult docstring) — one-sided, per rad/s,
    NOT per Hz. So it passes through with NO conversion factor: no 2*pi (which
    would be needed for a per-Hz spectrum) and no factor 2 (which would be
    needed for a two-sided one). We assume the paper's Sqq is one-sided, as is
    standard for this family of empirical fits; that assumption sets the
    absolute level but not the directivity shape this script checks.
    """
    q = 0.5 * rho * U_X ** 2
    return (q ** 2) * (delta_star / U_X) * chou_george_F(omega * delta_star / U_X)


# --------------------------------------------------------------------------- #
# Directivity sweep                                                           #
# --------------------------------------------------------------------------- #


def main() -> int:
    cfg = RotorConfig(Omega=OMEGA_ROTOR, Uz=UZ, B=1, rho=1.225, nu=1.48e-5,
                      c0=C0, n_azimuth=72)

    U_X = strip_relative_speed(cfg, RADIUS)
    delta_star = chou_george_delta_star(CHORD, CHI_DEG)

    print("Sinayoko/Kingan/Agarwal (2013) Table 2 — wind turbine element")
    print(f"  radius       {RADIUS:8.3f} m      chord    {CHORD:6.3f} m")
    print(f"  M_BO         {M_BO:8.3f}        M_FO     {M_FO:6.3f}")
    print(f"  Omega        {OMEGA_ROTOR:8.4f} rad/s  Uz       {UZ:6.3f} m/s")
    print(f"  span         {SPAN:8.3f} m      chi      {CHI_DEG:6.1f} deg")
    print(f"  U_X (rel)    {U_X:8.3f} m/s")
    print(f"  delta* (C&G) {delta_star:8.5f} m")
    print(f"  kC = {KC:g}  ->  omega = {OMEGA_ACOUSTIC:.2f} rad/s"
          f"  ->  f = {F_HZ:.2f} Hz")
    print(f"  observer radius {R_OBS:g} m, {N_THETA} polar angles\n")

    # The Chou & George spectrum drives the SUCTION side only; the pressure
    # side is handed an all-zero column, which noise_run skips. One empirical
    # spectrum is one Amiet source — putting the same Sqq on both surfaces
    # would double the radiated power (+3 dB) with no physical basis.
    def custom_WPS(omega_src):
        upper = chou_george_Sqq(omega_src, cfg.rho, U_X, delta_star)
        return upper, np.zeros_like(upper)

    # Ue_custom: the Amiet stage still needs an edge velocity even with a
    # custom spectrum (U_c = 0.7*Ue). Chou & George is built on U_X, so use it
    # on both surfaces — which is also the module's default, stated explicitly.
    strip = RotorStrip(radius=RADIUS, dr=SPAN, chord=CHORD,
                       pitch_rad=np.deg2rad(CHI_DEG),
                       custom_WPS_func=custom_WPS,
                       Ue_custom=(U_X, U_X))

    # rotor_noise_run needs >= 2 frequencies (it integrates an OASPL); take a
    # tight bracket about the target and read the centre bin for the PSD.
    freqs = np.array([0.99 * F_HZ, F_HZ, 1.01 * F_HZ])
    i_mid = 1

    # Observer ring in the y-z plane (contains the rotor axis, +z).
    # Theta measured from the rotor axis.
    theta = np.linspace(0.0, 2.0 * np.pi, N_THETA)
    obs = R_OBS * np.stack([np.zeros_like(theta),
                            np.sin(theta),
                            np.cos(theta)], axis=-1)

    res = rotor_noise_run(cfg, [strip], freqs, obs, verbose=True)

    # Normalise to 1 m as the paper's Fig. 9(b) does: the far-field PSD falls
    # as 1/R^2, so scaling by R^2 removes the spreading and leaves directivity.
    Spp_1m = res.Spp[:, i_mid] * R_OBS ** 2
    with np.errstate(divide="ignore"):
        SPL_1m = np.where(Spp_1m > 0.0,
                          10.0 * np.log10(2.0 * np.pi * Spp_1m / P_REF ** 2),
                          -200.0)

    d = res.diagnostics
    print(f"\n  doppler range   [{d['doppler_min']:.4f}, {d['doppler_max']:.4f}]")
    print(f"  omega/Omega     {d['omega_over_Omega_min']:.1f}")
    print(f"  f_kC1           {d['strips'][0]['f_kC1_Hz']:.2f} Hz "
          f"(band fraction below: {d['strips'][0]['frac_band_below_kC1']:.0%})")
    print(f"  SPL_1m range    [{SPL_1m.min():.2f}, {SPL_1m.max():.2f}] dB/Hz")

    # Fore/aft asymmetry: compare the two axial lobes.
    i_fore = int(np.argmin(np.abs(theta - 0.0)))
    i_aft = int(np.argmin(np.abs(theta - np.pi)))
    i_plane = int(np.argmin(np.abs(theta - 0.5 * np.pi)))
    print(f"  SPL_1m at Theta=0   (upstream axis) {SPL_1m[i_fore]:8.3f} dB/Hz")
    print(f"  SPL_1m at Theta=180 (downstream)    {SPL_1m[i_aft]:8.3f} dB/Hz")
    print(f"  SPL_1m at Theta=90  (rotor plane)   {SPL_1m[i_plane]:8.3f} dB/Hz")
    print(f"  fore/aft asymmetry                  "
          f"{SPL_1m[i_fore] - SPL_1m[i_aft]:8.3f} dB")
    contrast = max(SPL_1m[i_fore], SPL_1m[i_aft]) - SPL_1m[i_plane]
    print(f"  axis-to-plane contrast              {contrast:8.3f} dB")

    ratio = d["spanwise_neglect_ratio"]
    print(f"\n  obliquity ratio |X|/sqrt(Xc^2+Zn^2)  {ratio:8.1f}x  "
          f"(worst azimuth; xi -> 1 there)")
    if contrast > 3.0:
        print("  -> axis louder than the rotor plane: the expected dipole "
              "pattern.\n"
              "     Phase 1 gave -21.25 dB here (inverted); the general "
              "oblique-gust kernel\n"
              "     restored the ordering. See this file's header.")
    else:
        print("  -> rotor plane louder than the axis: the dipole null is "
              "INVERTED.\n"
              "     Phase 2 was supposed to fix exactly this — treat as a "
              "REGRESSION and check\n"
              "     tests/rotor_noise_test.py group 10 and the CHANGELOG "
              "entry.")

    # ---- figure ---------------------------------------------------------- #
    fig, ax = plt.subplots(figsize=(7.6, 8.2),
                           subplot_kw={"projection": "polar"})
    fig.subplots_adjust(top=0.78, bottom=0.06)
    ax.plot(theta, SPL_1m, lw=1.8, color="#1f6feb", zorder=3)

    ax.set_theta_zero_location("N")   # rotor axis (+z) points up
    ax.set_theta_direction(-1)
    lo = np.floor((SPL_1m.max() - 40.0) / 10.0) * 10.0
    ax.set_ylim(lo, np.ceil(SPL_1m.max() / 5.0) * 5.0)

    fig.suptitle(
        "GFoil rotor TE noise — directivity in a plane containing the rotor axis\n"
        f"Sinayoko et al. (2013) Table 2 wind turbine element, kC = {KC:g} "
        f"(f = {F_HZ:.1f} Hz)\n"
        "SPL normalised to 1 m [dB/Hz re 20 uPa] — QUALITATIVE SHAPE CHECK ONLY",
        fontsize=10, y=0.985,
    )

    ax.annotate(
        f"rotor plane (90/270 deg):\ndipole NULL, as expected\n"
        f"phase 1 put a {-21.25:.0f} dB artifact peak here\n"
        f"(mid-span kernel dropped Y_s)",
        xy=(np.deg2rad(90), SPL_1m[i_plane]), xycoords="data",
        xytext=(0.985, 0.845), textcoords="figure fraction",
        fontsize=7.8, color="#1a7f37", ha="right", va="top",
        arrowprops=dict(arrowstyle="->", color="#1a7f37", lw=1.1,
                        connectionstyle="arc3,rad=-0.3"),
    )
    ax.annotate(
        f"rotor axis (0/180 deg):\ndipole MAXIMUM, as expected\n"
        f"axis-to-plane contrast +{contrast:.1f} dB\n"
        f"(phase 1: -21.3 dB, inverted)",
        xy=(np.deg2rad(0), SPL_1m[i_fore]), xycoords="data",
        xytext=(0.015, 0.845), textcoords="figure fraction",
        fontsize=7.8, color="#444", ha="left", va="top",
        arrowprops=dict(arrowstyle="->", color="#444", lw=1.1,
                        connectionstyle="arc3,rad=0.3"),
    )
    ax.grid(alpha=0.35)

    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=140, bbox_inches="tight")
    print(f"\nFigure written to {OUT_PNG}")

    plot_cut()
    return 0


# --------------------------------------------------------------------------- #
# Second figure: the near-cutoff regularisation                               #
# --------------------------------------------------------------------------- #


def plot_cut() -> None:
    """Unbridged vs bridged |I(xi)| at the paper's Fig. 11 conditions.

    Roger & Moreau sec. 4.1 report that the two-step Schwarzschild solution
    leaves deep, narrow cuts exactly at the cut-off (kappa_bar = 0, i.e. xi = 1),
    from non-convergence of the iteration there, and regularise by matching
    d|I|/dK2 from both sides. This shows GFoil's realisation of that: a cubic
    Hermite in xi on |I| between anchors at kappa_bar = +/- 0.125 (the paper's
    own stated accuracy threshold for the back-scatter approximation).
    """
    from GFoil import gfoil_cpp

    chord, M, Ue = 0.13, 0.05, 0.05 * C0
    beta = np.sqrt(1.0 - M * M)
    R = 2.0

    # Observer angles with x1 <= 0 only. Eq. 14's three guarded denominators are
    # D and D +/- 2*kappa_bar; sweeping K2_bar free of the geometry breaks the
    # Eq. 18 tie between xi and (x1, x3) that normally keeps them sign-definite,
    # and for x1 > 0 the sweep then drives them through zero, adding off-manifold
    # spikes that have nothing to do with the cut this figure is about. For
    # x1 <= 0 they stay sign-definite even off-manifold, so these curves show the
    # cut and nothing else. (See tests/amiet_kernel_test.py group 1, and the
    # CHANGELOG note on the guarded denominators.)
    ANGLES = ((90.0, "#b3261e"), (140.0, "#1f6feb"))

    def probe(x, z, f):
        r = gfoil_cpp.amiet_kernel_I(dict(
            chord=chord, M=M, Ue=Ue, x=x, z=z, bridged=False,
            freqs_Hz=[f], K2_bar=[0.0]))
        return (float(r["mu_bar"][0]), float(r["xi_a"][0]), float(r["xi_b"][0]))

    def kern(x, z, f, mu, xi, bridged):
        r = gfoil_cpp.amiet_kernel_I(dict(
            chord=chord, M=M, Ue=Ue, x=x, z=z, bridged=bridged,
            freqs_Hz=list(np.full(xi.size, f)),
            K2_bar=list(xi * beta * mu)))
        return np.asarray(r["I_abs"])

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.4))
    for col, f in enumerate((200.0, 1000.0)):
        top, bot = axes[0, col], axes[1, col]
        for th, colour in ANGLES:
            t = np.deg2rad(th)
            x, z = R * np.cos(t), R * np.sin(t)
            mu, xa, xb = probe(x, z, f)

            # top row: zoom on the cut (linear xi)
            pad = 2.5 * (xb - xa)
            xz = np.linspace(xa - pad, xb + pad, 1200)
            top.plot(xz, kern(x, z, f, mu, xz, False), lw=3.0, color=colour,
                     alpha=0.30, label=f"$\\Theta$={th:.0f}$\\degree$, unbridged")
            top.plot(xz, kern(x, z, f, mu, xz, True), lw=1.2, color=colour,
                     label=f"$\\Theta$={th:.0f}$\\degree$, bridged")
            top.axvspan(xa, xb, color="#999", alpha=0.18, lw=0)

            # bottom row: full sweep (log-log), showing the subcritical decay.
            # The bridge window is only ~0.01 wide in xi at 1 kHz, so a plain
            # log sweep aliases it into a spurious notch: sample it explicitly.
            xf = np.unique(np.concatenate([
                np.logspace(np.log10(0.05), np.log10(20.0), 900),
                np.linspace(xa - pad, xb + pad, 400)]))
            bot.loglog(xf, kern(x, z, f, mu, xf, True), lw=1.2, color=colour,
                       label=f"$\\Theta$={th:.0f}$\\degree$")

        top.axvline(1.0, color="k", ls=":", lw=0.9)
        top.set_title(f"f = {f:.0f} Hz — cut-off region", fontsize=10)
        top.set_xlabel(r"$\xi = \bar{K}_2/(\beta\bar{\mu})$")
        top.grid(alpha=0.3)
        bot.axvline(1.0, color="k", ls=":", lw=0.9)
        bot.set_title(f"f = {f:.0f} Hz — full sweep (bridged)", fontsize=10)
        bot.set_xlabel(r"$\xi = \bar{K}_2/(\beta\bar{\mu})$")
        bot.grid(alpha=0.3, which="both")
    axes[0, 0].set_ylabel(r"$|I|$")
    axes[1, 0].set_ylabel(r"$|I|$")
    axes[0, 0].legend(fontsize=7.5)
    axes[1, 0].legend(fontsize=7.5)
    fig.suptitle(
        "Near-cutoff regularisation of the Roger & Moreau radiation integral\n"
        "c = 0.13 m, M = 0.05, mid-span observer (their Fig. 11 conditions).  "
        r"Shaded = bridge window, $\bar{\kappa}_{reg}=0.125$;  dotted = cut-off "
        r"$\xi=1$." "\n"
        "Top: the two-step Schwarzschild solution's cut at the cut-off, and the "
        "cubic-Hermite bridge that removes it.\n"
        "Bottom: supercritical plateau, and the subcritical decay for "
        r"$\xi \gg 1$ (steeper at the higher frequency).",
        fontsize=9.5)
    fig.tight_layout(rect=(0, 0, 1, 0.89))
    out = OUT_PNG.parent / "rotor_noise_cut.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"Figure written to {out}")


if __name__ == "__main__":
    sys.exit(main())
