"""Rotor broadband trailing-edge noise via the corrected Schlinker-Amiet procedure.

Post-processing only. This module never touches C++, is never taped, and is
never on the AD/optimisation path — it wraps the acoustics-only ``noise_run``
entry point in the rotating-blade formulation of

    S. Sinayoko, M. Kingan & A. Agarwal (2013),
    "Trailing edge noise theory for rotating blades in uniform flow",
    Proc. R. Soc. A 469:20130065.

Their central result is that the *fixed-aerofoil* Amiet PSD — evaluated in the
blade frame at the Doppler-shifted **source** frequency and weighted by the
**squared** Doppler ratio, then averaged over azimuth — reproduces the exact
rotating-frame formulation to within 1 dB for kC > 1 and omega >> Omega, up to
chordwise Mach 0.95.


Equation map
------------
=========================  ====================================================
Paper                      Implementation
=========================  ====================================================
Eq. 4.6 (far-field Te)     ``_emission_time`` — GENERALISED. The paper takes
                           xe ~ 0 (source at the hub) to get a closed form. We
                           instead solve the exact quadratic with xe finite, so
                           the module stays valid for observers at moderate
                           distance. Reduces to Eq. 4.6 as |xo| >> r (checked in
                           the test suite).
Eq. 4.1 (present posn xp)  ``_source_positions``
Eq. 4.10 (Doppler ratio)   ``_doppler_ratio``
Eq. 4.9 + App. B (frame)   ``_blade_frame_observer``, ``_Rz``, ``_Ry``
Eq. 4.12 (time-averaged
  rotor PSD, exponent 2)   ``rotor_noise_run``
Eqs. 3.9-3.11 (Chou &
  George wall pressure)    NOT here — see ``bench/rotor_noise_paper_check.py``
=========================  ====================================================


The five-step procedure, per (strip, azimuth, observer)
-------------------------------------------------------
1. Solve the emission-time (propagation) equation for the retarded time Te.
2. Form the convected source position xc and the *present* source position xp.
3. Form the Doppler ratio from the convected source-to-observer direction, and
   with it the source frequencies omega' = omega / doppler.
4. Rotate the observer into the blade section frame about xp, and evaluate the
   fixed-aerofoil PSD there at omega' via ``noise_run``.
5. Weight by (omega'/omega)^2 = 1/doppler^2 and average over azimuth; sum
   strips incoherently and multiply by the blade count B.

Because the wall-pressure spectrum is a function of *source* frequency, each
``noise_run`` call is handed the Doppler-shifted grid directly — that evaluates
the WPS model and the radiation integral at the same, correct frequencies. A WPS
precomputed on the observer grid and reused would be evaluated at the wrong
frequencies, so this module never does that.


Doppler exponent
----------------
The weight is ``(omega'/omega)**2``, i.e. exponent **+2**. This is the paper's
headline correction to the original Schlinker & Amiet (1981) procedure and is
what their §5 validation against the exact rotating-frame solution supports. Do
not "fix" this to exponent 1 or -2.


Coordinate frames
-----------------
Hub-fixed, right-handed, origin at the hub centre; rotor turns about +z at
Omega. Uniform axial inflow of speed Uz flows in the **-z** direction, so the
flow Mach vector seen by a stationary observer is ``M_FO = -(Uz/c0) z_hat``. A
strip at radius r and azimuth gamma sits at ``xe = (r cos g, r sin g, 0)`` and
moves along ``gamma_hat = (-sin g, cos g, 0)`` at ``M_BO = (Omega r/c0) g_hat``.

The blade section frame produced by Eq. 4.9 is ``X = (X_c, Y_s, Z_n)``: X_c
chordwise and pointing downstream, Y_s spanwise, Z_n plate-normal. That is the
Roger-Moreau frame the Amiet stage of ``noise_run`` expects.

Section angle of attack chi does NOT enter the rotation. It enters only through
the BL states the caller supplies (from their own ``fwd_run`` at AoA chi). This
matches GFoil's existing convention — its Amiet frame is chord-aligned and its
Amiet Mach uses the full freestream speed — and is second-order-equivalent to
the paper's Appendix-D ``cos chi`` projection at small AoA.


Frame-mapping contract with ``noise_run`` (verified against the source)
-----------------------------------------------------------------------
All rotations are done here, so every ``noise_run`` call passes
``alphaDeg = 0``. At alpha = 0 the kernel's global -> TE-local transform
(``src/include/noise_run.hpp``, lines 142-153) reduces to::

    const Real te_offset = Real(0.75) * Real(chord);
    Real x_loc = obsX*cos_a - obsZ*sin_a - te_offset;   // cos_a=1, sin_a=0
    Real y_loc = obsY;
    Real z_loc = obsX*sin_a + obsZ*cos_a;

i.e. ``x_loc = obsX - 0.75*chord``, ``y_loc = obsY``, ``z_loc = obsZ``. So to
make the kernel's TE-local coordinates come out as exactly the paper-frame
``(X_c, Y_s, Z_n)`` we pass ``observerXYZ = (X_c + 0.75*chord, Y_s, Z_n)``.
(``noise_run``'s global frame has its origin at quarter-chord; the trailing edge
is 0.75 chord downstream of it.) Verified numerically: observer (1, 0, 1) with
chord 1 returns ``obsXYZ_TElocal = (0.25, 0, 1)``.

The kernel's Uinf is not a direct argument — it is rederived internally as
``Re*nu/chord``. That identity is therefore the ONLY handle this wrapper has on
the kernel's Uinf and Amiet Mach, so we set ``Re = U_rel*chord/nu`` and assert
the round-trip.

One further contract, which ``noise_run``'s own docstring does not spell out:
supplying ``custom_WPS`` skips the BL -> WPS model but **not** the use of
``Ue``. The Amiet stage still forms ``U_c = 0.7*Ue`` and from it the Corcos
correlation length and the convective wavenumber, so ``Ue = 0`` produces a
**NaN** far field, not silence. ``RotorStrip.Ue_custom`` therefore exists, and
defaults to ``U_rel``. Every ``noise_run`` result is checked for finiteness
before it reaches the dB conversion, because ``NaN > 0`` is False and a NaN
would otherwise be floored to -200 dB and read as a quiet observer — the silent
blank failure the forward path names ``acoustic_nan``.


Validity regime
---------------
* ``kC > 1``   — chord acoustically compact below this; the Amiet scattering
  solution is a high-frequency (Schwarzschild) approximation.
* ``omega >> Omega`` — the azimuthal average assumes the source spectrum is
  stationary over a revolution.
* ``l_S << r`` — strip theory: the spanwise correlation length must be small
  against the radius.

``rotor_noise_run`` reports all three in ``RotorNoiseResult.diagnostics`` and
warns when they look marginal; none of them ever fails a run.


Limitations
-----------
The phase-1 **mid-span kernel** limitation is **RESOLVED**. ``noise_run``'s
Amiet stage now implements the general three-dimensional oblique-gust
formulation of Roger & Moreau (2005): it forms
``S0 = sqrt(x1^2 + beta^2*(x2^2 + x3^2))``, selects the gust
``K2_bar = k_bar*x2/S0`` of their Eq. 18, branches on criticality, and applies
the spanwise-wavenumber-corrected Corcos length. The rotor directivity that
phase 1 inverted is correct: on the paper's Table 2 wind-turbine element the
rotor axis is now 12.5 dB louder than the rotor plane, against -21.3 dB before.
See CHANGELOG "General oblique-gust Amiet kernel (phase 2)". The remaining
limitations are:

* **Eq. 18 gust selection.** The kernel takes the paper's large-aspect-ratio
  limit, in which the spanwise wavenumber integral collapses to a delta at
  ``K2 = k*x2/S0`` — one oblique gust per (observer, frequency) rather than a
  finite-span sinc over many. That holds for ``L/(2b) >> 1``; strips of
  near-unit aspect ratio are outside it.
* **Axial inflow only.** No cross-flow and no shaft angle, so the section
  relative Mach is azimuth-independent and one BL state per strip is reused at
  every azimuth. The paper's Appendix-D azimuth-dependent M_X is deferred.
* **Kernel c0 fixed at 340 m/s.** ``noise_run.hpp`` hard-codes
  ``M_amiet = Uinf/340``. ``RotorConfig.c0`` drives all wrapper-side kinematics
  but cannot reach the kernel; a mismatch > 1 m/s warns.
* **No A-weighting.** ``Spp`` and ``OASPL_perObs`` are unweighted. (A future
  option — the forward path already has an A-weighting filter in ``sound.hpp``
  that could be mirrored here.)
"""

import warnings
from dataclasses import dataclass
from typing import Callable, Optional, Sequence

import numpy as np

from .gfoil import noise_run
from .inputs import _ResultMixin, _as_1d_float_array, _as_float_array

__all__ = [
    "RotorConfig",
    "RotorStrip",
    "RotorNoiseResult",
    "rotor_noise_run",
    "strip_relative_speed",
]


# Reference pressure for dB re 20 micro-Pa.
_P_REF = 20e-6
_P_REF2 = _P_REF * _P_REF

# The kernel places the global-frame origin at quarter-chord, so the trailing
# edge sits 0.75*chord downstream (noise_run.hpp: te_offset = 0.75*chord).
_TE_OFFSET_FRAC = 0.75

# Speed of sound hard-coded in the Amiet kernel (noise_run.hpp: Uinf/340.0).
_KERNEL_C0 = 340.0

# PSD dB floor for non-positive linear values (10*log10(0) = -inf otherwise).
_SPP_DB_FLOOR = -200.0

# OASPL floor for a silent observer. Matches sound.hpp's calc_OASPL, which
# returns 10*log10(1e-30) = -300 dB rather than -inf when there is no source.
_OASPL_FLOOR = 10.0 * np.log10(1e-30)


# --------------------------------------------------------------------------- #
# Inputs                                                                      #
# --------------------------------------------------------------------------- #


@dataclass
class RotorConfig:
    """Rotor operating point and hub-frame kinematics.

    Attributes
    ----------
    Omega : float
        Rotor angular velocity [rad/s], rotation about +z.
    Uz : float
        Axial inflow speed [m/s]; the flow runs in the -z direction, so the
        flow Mach vector is ``M_FO = -(Uz/c0) z_hat``. Must be >= 0.
    B : int
        Blade count, >= 1. Broadband noise from different blades is
        uncorrelated, so blades add in power: a linear factor B on the PSD.
    rho : float
        Density [kg/m^3].
    nu : float
        Kinematic viscosity [m^2/s]. Physical — it sets the per-strip Reynolds
        number handed to ``noise_run``, and with it the kernel's Uinf (see the
        frame-mapping contract in the module docstring).
    c0 : float
        Speed of sound [m/s]. Drives every wrapper-side kinematic quantity
        (Mach vectors, emission time, Doppler). It does NOT reach the Amiet
        kernel, which hard-codes 340 m/s; a deviation > 1 m/s warns.
    n_azimuth : int
        Azimuth samples per revolution. 72 is comfortably converged for the
        cases tested (72 vs 144 differ by < 0.02 dB OASPL).
    """

    Omega: float
    Uz: float
    B: int
    rho: float = 1.225
    nu: float = 1.48e-5
    c0: float = 340.0
    n_azimuth: int = 72

    def __post_init__(self):
        self.Omega = float(self.Omega)
        self.Uz = float(self.Uz)
        self.B = int(self.B)
        self.rho = float(self.rho)
        self.nu = float(self.nu)
        self.c0 = float(self.c0)
        self.n_azimuth = int(self.n_azimuth)

        if self.Uz < 0.0:
            raise ValueError(
                f"RotorConfig.Uz must be >= 0 (axial inflow speed; the -z "
                f"direction is already built in); got {self.Uz}"
            )
        if self.B < 1:
            raise ValueError(f"RotorConfig.B must be >= 1, got {self.B}")
        if self.n_azimuth < 1:
            raise ValueError(
                f"RotorConfig.n_azimuth must be >= 1, got {self.n_azimuth}"
            )
        if self.rho <= 0.0 or self.nu <= 0.0 or self.c0 <= 0.0:
            raise ValueError("RotorConfig.rho, .nu and .c0 must all be > 0")

        if self.Uz / self.c0 >= 1.0:
            raise ValueError(
                f"axial inflow is supersonic (Uz/c0 = {self.Uz / self.c0:.3f}); "
                "the emission-time solve assumes |M_FO| < 1"
            )
        if abs(self.c0 - _KERNEL_C0) > 1.0:
            warnings.warn(
                f"RotorConfig.c0 = {self.c0} m/s is used for all wrapper-side "
                f"kinematics (Mach vectors, emission time, Doppler), but the "
                f"underlying Amiet kernel hard-codes c0 = {_KERNEL_C0} m/s "
                f"(noise_run.hpp: M_amiet = Uinf/340). The two are "
                f"inconsistent; results are approximate.",
                stacklevel=2,
            )


@dataclass
class RotorStrip:
    """One radial blade element.

    The BL states are the caller's responsibility: run your own ``fwd_run`` at
    this strip's sectional conditions (relative speed ``strip_relative_speed``,
    section AoA chi, sectional Reynolds number) and pass
    ``result.verbose_data.BL_top`` / ``.BL_bot`` straight through. This module
    performs no aero solves.

    Exactly one of (``BL_top`` and ``BL_bot``) or ``custom_WPS_func`` must be
    supplied.

    Attributes
    ----------
    radius : float
        Strip centre radius [m], > 0.
    dr : float
        Strip radial extent [m], > 0. Passed to ``noise_run`` as ``span``; the
        far-field PSD is linear in it.
    chord : float
        Section chord [m], > 0.
    pitch_rad : float
        Geometric pitch (stagger) angle of the *chord line* relative to the
        rotor plane [rad]. Note this is not the section angle of attack — see
        the module docstring.
    BL_top, BL_bot : ndarray, shape (7,), optional
        Upper / lower surface trailing-edge BL states, ordered
        ``[theta, delta_star, tau_max, Ue, dpdx, tau_wall, delta99]`` — exactly
        ``FwdResult.verbose_data.BL_top`` / ``.BL_bot``. A surface with
        ``tau_max <= 0`` (fully-laminar TE) is skipped by the kernel.
    wps_model : str
        ``noise_run`` WPS model key: one of 'roz', 'goo', 'lee', 'kam', 'tno'.
        Used only when BL states are supplied.
    custom_WPS_func : callable, optional
        Alternative to BL states: ``f(omega_src) -> (WPS_upper, WPS_lower)``,
        where ``omega_src`` is the Doppler-shifted **source** angular-frequency
        array [rad/s] and both returned arrays are raw linear one-sided
        wall-pressure spectra in Pa^2/(rad/s) — the units ``noise_run``'s
        ``custom_WPS`` argument expects (see ``NoiseResult``: "Spectra are raw
        linear Pa^2/omega"). Called once per (strip, azimuth, observer), always
        at the source frequencies. An all-zero column skips that surface.
    Ue_custom : (float, float), optional
        Boundary-layer edge velocities ``(upper, lower)`` [m/s] for the
        custom-WPS path only. Defaults to ``U_rel`` on both surfaces.

        This is NOT redundant with ``custom_WPS_func``. Supplying ``custom_WPS``
        makes ``noise_run`` skip the BL -> WPS model, but its *Amiet* stage
        still reads ``Ue``: ``TE_noise_outer_vec`` forms the convection velocity
        ``U_c = 0.7*Ue`` and from it both the Corcos spanwise correlation length
        ``l_y = l_c/(1 + (K2*l_c)^2)``, ``l_c = 1.47*U_c/omega``, and the
        convective wavenumber ratio ``alpha = U/U_c``. So ``Ue`` is load-bearing
        even with a custom spectrum, and ``Ue = 0`` yields a NaN far field
        rather than silence.
        The ``U_rel`` default mirrors ``noise_run.hpp``'s own fallback for a
        skipped surface (``Ue_top = upper_skip ? Uinf : Ue[0]``).
    """

    radius: float
    dr: float
    chord: float
    pitch_rad: float
    BL_top: Optional[np.ndarray] = None
    BL_bot: Optional[np.ndarray] = None
    wps_model: str = "kam"
    custom_WPS_func: Optional[Callable] = None
    Ue_custom: Optional[Sequence[float]] = None

    def __post_init__(self):
        self.radius = float(self.radius)
        self.dr = float(self.dr)
        self.chord = float(self.chord)
        self.pitch_rad = float(self.pitch_rad)

        if self.radius <= 0.0:
            raise ValueError(f"RotorStrip.radius must be > 0, got {self.radius}")
        if self.dr <= 0.0:
            raise ValueError(f"RotorStrip.dr must be > 0, got {self.dr}")
        if self.chord <= 0.0:
            raise ValueError(f"RotorStrip.chord must be > 0, got {self.chord}")

        has_bl = self.BL_top is not None and self.BL_bot is not None
        has_custom = self.custom_WPS_func is not None

        if has_bl and has_custom:
            raise ValueError(
                "RotorStrip: supply EITHER BL_top/BL_bot OR custom_WPS_func, "
                "not both (the custom spectrum would silently win, since "
                "noise_run skips the BL->WPS path whenever custom_WPS is given)"
            )
        if not has_bl and not has_custom:
            raise ValueError(
                "RotorStrip: supply either BL_top and BL_bot (both, shape (7,)) "
                "or custom_WPS_func"
            )
        if has_custom and not callable(self.custom_WPS_func):
            raise ValueError("RotorStrip.custom_WPS_func must be callable")
        if (self.BL_top is None) != (self.BL_bot is None):
            raise ValueError(
                "RotorStrip: BL_top and BL_bot must be supplied together"
            )
        if self.Ue_custom is not None:
            if has_bl:
                raise ValueError(
                    "RotorStrip.Ue_custom applies to the custom_WPS_func path "
                    "only; with BL states the edge velocities come from "
                    "BL_top[3] / BL_bot[3]"
                )
            ue = _as_1d_float_array(self.Ue_custom, "RotorStrip.Ue_custom")
            if ue.size != 2:
                raise ValueError(
                    f"RotorStrip.Ue_custom must be (Ue_upper, Ue_lower); "
                    f"got {ue.size} element(s)"
                )
            if np.any(ue <= 0.0):
                raise ValueError(
                    f"RotorStrip.Ue_custom must be strictly positive "
                    f"(the Amiet stage divides by U_c = 0.7*Ue); got {ue}"
                )
            self.Ue_custom = ue

        if has_bl:
            self.BL_top = _as_1d_float_array(self.BL_top, "RotorStrip.BL_top")
            self.BL_bot = _as_1d_float_array(self.BL_bot, "RotorStrip.BL_bot")
            if self.BL_top.size != 7 or self.BL_bot.size != 7:
                raise ValueError(
                    f"RotorStrip.BL_top/BL_bot must have length 7 "
                    f"[theta, delta_star, tau_max, Ue, dpdx, tau_wall, delta99]; "
                    f"got {self.BL_top.size} and {self.BL_bot.size}"
                )


# --------------------------------------------------------------------------- #
# Output                                                                      #
# --------------------------------------------------------------------------- #


@dataclass(repr=False)
class RotorNoiseResult(_ResultMixin):
    """Returned by ``rotor_noise_run``.

    ``Spp`` is the time-averaged far-field broadband TE-noise PSD of the whole
    rotor, on the caller's **observer**-frame frequency grid.
    """

    _REPR_GROUPS = [
        (None, ["freqs_Hz", "Spp", "Spp_dB", "OASPL_perObs"]),
        ("diagnostics", ["diagnostics"]),
        ("verbose", ["Spp_perStrip"]),
    ]

    freqs_Hz: np.ndarray      # (N,)        echo of the input observer frequencies
    Spp: np.ndarray           # (nObs, N)   rotor PSD [Pa^2/(rad/s)], linear
    Spp_dB: np.ndarray        # (nObs, N)   10*log10(2*pi*Spp/pref^2)  [dB/Hz]
    OASPL_perObs: np.ndarray  # (nObs,)     [dB re 20 uPa]
    diagnostics: dict
    Spp_perStrip: Optional[np.ndarray] = None  # (nStrips, nObs, N), verbose only
                                               # pre-B-multiplication


# --------------------------------------------------------------------------- #
# Geometry / kinematics primitives                                            #
# --------------------------------------------------------------------------- #


def _Rz(theta: float) -> np.ndarray:
    """Rotation about z (paper Appendix B)."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s, 0.0],
                     [s,  c, 0.0],
                     [0.0, 0.0, 1.0]])


def _Ry(theta: float) -> np.ndarray:
    """Rotation about y (paper Appendix B)."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, 0.0, -s],
                     [0.0, 1.0, 0.0],
                     [s, 0.0,  c]])


def strip_relative_speed(cfg: RotorConfig, radius: float) -> float:
    """U_rel = sqrt((Omega*r)^2 + Uz^2) — sectional relative speed.

    Valid for axial inflow only, where the rotational and inflow components are
    orthogonal and the result is azimuth-independent. This is the speed to use
    for the strip's own ``fwd_run``, and it is what sets the kernel's Uinf here
    (via ``Re = U_rel*chord/nu``).
    """
    return float(np.hypot(cfg.Omega * radius, cfg.Uz))


def _emission_time(xo: np.ndarray, xe: np.ndarray, M_FO: np.ndarray) -> np.ndarray:
    """Solve c0*Te = |xo - xe - M_FO*c0*Te| for the propagation distance c0*Te.

    Generalises the paper's Eq. 4.6, which assumes xe ~ 0 (far field). Keeping
    xe finite makes the module valid at moderate observer distances; the two
    agree as |xo| >> |xe|.

    With d = xo - xe, R = |d| and M = M_FO, squaring gives

        (1 - |M|^2)(c0 Te)^2 + 2 (M.d)(c0 Te) - R^2 = 0

    whose positive root (the physical one, and the only one for |M| < 1) is

        c0 Te = [ -(M.d) + sqrt((M.d)^2 + (1 - |M|^2) R^2) ] / (1 - |M|^2)

    Parameters
    ----------
    xo : ndarray, shape (3,)
    xe : ndarray, shape (..., 3)
    M_FO : ndarray, shape (3,)

    Returns
    -------
    ndarray, shape (...,)
        c0*Te, the propagation distance [m].
    """
    d = xo - xe                       # (..., 3)
    R2 = np.sum(d * d, axis=-1)
    Md = d @ M_FO
    M2 = float(M_FO @ M_FO)
    one_minus_M2 = 1.0 - M2
    if one_minus_M2 <= 0.0:
        raise ValueError(
            f"emission-time solve needs |M_FO| < 1, got |M_FO| = {np.sqrt(M2):.4f}"
        )
    return (-Md + np.sqrt(Md * Md + one_minus_M2 * R2)) / one_minus_M2


def _doppler_ratio(M_BO: np.ndarray,
                   M_FO: np.ndarray,
                   CO: np.ndarray) -> np.ndarray:
    """Doppler ratio omega/omega' (paper Eq. 4.10).

        doppler = 1 + (M_BO.CO) / (1 + (M_FO - M_BO).CO)

    which rearranges to the equivalent, more obviously classical

        doppler = (1 + M_FO.CO) / (1 + (M_FO - M_BO).CO)

    Sanity limits: a static source in a moving medium (M_BO = 0) gives exactly
    1; a source in still air (M_FO = 0) approaching the observer head-on gives
    1/(1 - |M_BO|), the textbook moving-source factor.

    Parameters
    ----------
    M_BO : ndarray, shape (..., 3)   blade Mach vector
    M_FO : ndarray, shape (3,)       flow Mach vector
    CO : ndarray, shape (..., 3)     unit vector, convected source -> observer

    Returns
    -------
    ndarray, shape (...,)
    """
    a = np.sum(M_BO * CO, axis=-1)              # M_BO . CO
    denom = 1.0 + (CO @ M_FO) - a               # 1 + (M_FO - M_BO).CO
    if np.any(denom <= 0.0):
        raise ValueError(
            "Doppler denominator 1 + (M_FO - M_BO).CO <= 0: the source is "
            "sonic or supersonic relative to the observer direction, where "
            "Eq. 4.10 breaks down (multiple emission times). "
            f"min = {float(np.min(denom)):.6g}"
        )
    doppler = 1.0 + a / denom
    if np.any(doppler <= 0.0):
        raise ValueError(
            f"non-positive Doppler ratio (min = {float(np.min(doppler)):.6g}); "
            "source frequencies omega/doppler would be non-physical"
        )
    return doppler


def _blade_frame_geometry(cfg: RotorConfig,
                          strip: RotorStrip,
                          xo: np.ndarray):
    """Per-azimuth blade-frame observer positions and Doppler ratios.

    Vectorised over the azimuth samples gamma_j = 2*pi*j/n_azimuth,
    j = 0 .. n_azimuth-1.

    Returns
    -------
    X : ndarray, shape (n_azimuth, 3)
        Observer in the blade section frame (X_c, Y_s, Z_n), origin at the
        present source position xp.
    doppler : ndarray, shape (n_azimuth,)
        omega/omega' at each azimuth.
    """
    n = cfg.n_azimuth
    gamma = 2.0 * np.pi * np.arange(n) / n

    cg, sg = np.cos(gamma), np.sin(gamma)
    zero = np.zeros(n)

    xe = strip.radius * np.stack([cg, sg, zero], axis=-1)      # (n,3) emission posn
    ghat = np.stack([-sg, cg, zero], axis=-1)                  # (n,3) tangential
    M_BO = (cfg.Omega * strip.radius / cfg.c0) * ghat          # (n,3) blade Mach
    M_FO = np.array([0.0, 0.0, -cfg.Uz / cfg.c0])              # (3,)  flow Mach

    c0Te = _emission_time(xo, xe, M_FO)                        # (n,)

    # Convected and present source positions (paper Eq. 4.1).
    xc = xe + M_FO * c0Te[:, None]
    xp = xe + M_BO * c0Te[:, None]

    # Convected source -> observer. |xo - xc| == c0*Te by construction, but
    # normalise explicitly rather than lean on that identity.
    r_co = xo - xc
    CO = r_co / np.linalg.norm(r_co, axis=-1, keepdims=True)

    doppler = _doppler_ratio(M_BO, M_FO, CO)                   # (n,)

    # Observer into the blade section frame (paper Eq. 4.9).
    #   x1 = xo - xp;  x2 = Rz(pi/2 - gamma) x1;  X = Ry(pitch) x2
    x1 = xo - xp                                               # (n,3)
    Rz_stack = np.stack([_Rz(0.5 * np.pi - g) for g in gamma]) # (n,3,3)
    x2 = np.einsum("nij,nj->ni", Rz_stack, x1)
    X = x2 @ _Ry(strip.pitch_rad).T                            # (n,3)

    return X, doppler


# --------------------------------------------------------------------------- #
# Diagnostics                                                                 #
# --------------------------------------------------------------------------- #


def _kC1_frequency(c0: float, chord: float) -> float:
    """Frequency at which the acoustic chord parameter kC reaches 1.

    kC = (omega/c0)*C = (2*pi*f/c0)*C, so kC > 1  <=>  f > c0/(2*pi*C).
    """
    return c0 / (2.0 * np.pi * chord)


def _strip_diagnostics(cfg: RotorConfig,
                       strip: RotorStrip,
                       freqs_Hz: np.ndarray,
                       omega: np.ndarray,
                       U_rel: float,
                       Re: float) -> dict:
    """Advisory validity metrics for one strip. Warns; never raises."""
    f_kC1 = _kC1_frequency(cfg.c0, strip.chord)
    frac_below = float(np.mean(freqs_Hz < f_kC1))

    d = {
        "radius": strip.radius,
        "chord": strip.chord,
        "U_rel": U_rel,
        "Re": Re,
        "f_kC1_Hz": f_kC1,
        "frac_band_below_kC1": frac_below,
        "l_S": None,
        "l_S_over_radius": None,
    }

    if frac_below > 0.2:
        warnings.warn(
            f"strip r={strip.radius:g} m: {100.0 * frac_below:.0f}% of the "
            f"requested band lies below kC = 1 (f < {f_kC1:.1f} Hz), where the "
            f"Amiet high-frequency scattering solution is outside its validity "
            f"regime.",
            stacklevel=3,
        )

    # Corcos spanwise correlation length at the band minimum (the worst case:
    # l_S falls as 1/omega). Needs an edge velocity, so skip for custom WPS.
    #
    # This is the K2 = 0 length. The phase-2 kernel uses the general
    # l_y(K2) = l_c/(1 + (K2*l_c)^2) <= l_c, so l_c remains the correct
    # WORST CASE for the strip-theory check l_S << r: obliquity only shortens
    # the correlation length, never lengthens it.
    if strip.BL_top is not None:
        Ue_top = float(strip.BL_top[3])
        if Ue_top > 0.0 and omega[0] > 0.0:
            U_c = 0.7 * Ue_top
            l_S = 1.47 * U_c / omega[0]
            d["l_S"] = float(l_S)
            d["l_S_over_radius"] = float(l_S / strip.radius)
            if l_S / strip.radius > 0.2:
                warnings.warn(
                    f"strip r={strip.radius:g} m: Corcos spanwise correlation "
                    f"length l_S = {l_S:.3g} m is {l_S / strip.radius:.2f} of "
                    f"the radius at the lowest requested frequency; strip "
                    f"theory assumes l_S << r.",
                    stacklevel=3,
                )
    return d


# --------------------------------------------------------------------------- #
# Main entry point                                                            #
# --------------------------------------------------------------------------- #


def rotor_noise_run(cfg: RotorConfig,
                    strips: Sequence[RotorStrip],
                    freqs_Hz: np.ndarray,
                    observerXYZ: np.ndarray,
                    verbose: bool = False) -> RotorNoiseResult:
    """Time-averaged far-field broadband rotor TE-noise PSD (paper Eq. 4.12).

        S_pp(xo, omega) = B * sum_strips (1/Ng) sum_j (omega'_j/omega)^2
                                                * S'_pp(X_j, omega'_j)

    Blades add in power (broadband, mutually uncorrelated); strips add
    incoherently (strip theory); the azimuth average is a uniform Ng-point mean.
    ``S'_pp`` is the fixed-aerofoil PSD from ``noise_run``, evaluated in the
    blade frame at the Doppler-shifted source frequencies.

    One ``noise_run`` call is made per (strip, azimuth, observer). Observers
    cannot be batched into a single call even though ``noise_run`` accepts many:
    each observer sees its own Doppler ratio and therefore its own source
    frequency grid, and a call carries only one frequency array.

    Parameters
    ----------
    cfg : RotorConfig
    strips : sequence of RotorStrip
        At least one. Strips are independent; nothing checks that they tile the
        blade without gaps or overlap.
    freqs_Hz : array-like, shape (N,)
        **Observer**-frame frequencies [Hz], ascending, all > 0. Any spacing.
    observerXYZ : array-like, shape (3,) or (nObs, 3)
        Observer position(s) in the hub-fixed Cartesian frame [m].
    verbose : bool
        Populate ``Spp_perStrip`` and print per-strip progress.

    Returns
    -------
    RotorNoiseResult
    """
    if len(strips) == 0:
        raise ValueError("strips must contain at least one RotorStrip")
    for i, s in enumerate(strips):
        if not isinstance(s, RotorStrip):
            raise TypeError(f"strips[{i}] is {type(s).__name__}, not RotorStrip")

    freqs = _as_1d_float_array(freqs_Hz, "freqs_Hz")
    N = freqs.size
    if N < 2:
        raise ValueError(
            f"freqs_Hz needs at least 2 points to integrate an OASPL, got {N}"
        )
    if np.any(freqs <= 0.0):
        raise ValueError("freqs_Hz must be strictly positive")
    if np.any(np.diff(freqs) <= 0.0):
        raise ValueError("freqs_Hz must be strictly ascending")

    obs = _as_float_array(observerXYZ, "observerXYZ").astype(float)
    if obs.ndim == 1 and obs.size == 3:
        obs = obs.reshape(1, 3)
    elif not (obs.ndim == 2 and obs.shape[1] == 3):
        raise ValueError(
            f"observerXYZ must be shape (3,) or (nObs,3), got {obs.shape}"
        )
    nObs = obs.shape[0]

    omega = 2.0 * np.pi * freqs
    nStrips = len(strips)

    Spp_perStrip = np.zeros((nStrips, nObs, N))
    strip_diags = []
    dop_min, dop_max = np.inf, -np.inf
    span_neglect = 0.0

    for iStrip, strip in enumerate(strips):
        # Axial inflow => U_rel is azimuth-independent, so the sectional Re
        # (and with it the kernel's Uinf) is computed once per strip.
        U_rel = strip_relative_speed(cfg, strip.radius)
        if U_rel <= 0.0:
            raise ValueError(
                f"strips[{iStrip}]: relative speed is zero (Omega = {cfg.Omega}, "
                f"Uz = {cfg.Uz}) — the kernel's Uinf would be zero"
            )
        Re = U_rel * strip.chord / cfg.nu

        # The kernel takes no Uinf argument; it rederives Uinf = Re*nu/chord.
        # That identity is our only handle on it, so verify the round-trip.
        U_kernel = Re * cfg.nu / strip.chord
        if abs(U_kernel - U_rel) > 1e-12 * abs(U_rel):
            raise AssertionError(
                f"strips[{iStrip}]: Re round-trip failed — the kernel would use "
                f"Uinf = Re*nu/chord = {U_kernel!r} m/s instead of the intended "
                f"U_rel = {U_rel!r} m/s"
            )

        # BL states handed to noise_run. On the custom-WPS path the BL -> WPS
        # model is skipped, but the Amiet stage still reads Ue (index 3) for
        # U_c = 0.7*Ue, so a bare zeros(7) would give a NaN far field — hence
        # the U_rel default (matching the kernel's own skipped-surface
        # fallback). Every other slot is genuinely unused there.
        if strip.custom_WPS_func is not None:
            ue_u, ue_l = (strip.Ue_custom if strip.Ue_custom is not None
                          else (U_rel, U_rel))
            BL_top = np.zeros(7)
            BL_bot = np.zeros(7)
            BL_top[3] = ue_u
            BL_bot[3] = ue_l
        else:
            BL_top, BL_bot = strip.BL_top, strip.BL_bot

        strip_diags.append(
            _strip_diagnostics(cfg, strip, freqs, omega, U_rel, Re)
        )

        if verbose:
            print(f"[rotor_noise] strip {iStrip + 1}/{nStrips}: "
                  f"r={strip.radius:g} m  c={strip.chord:g} m  "
                  f"U_rel={U_rel:.2f} m/s  Re={Re:.3g}  "
                  f"({cfg.n_azimuth} azimuths x {nObs} obs)")

        for iObs in range(nObs):
            xo = obs[iObs]
            X, doppler = _blade_frame_geometry(cfg, strip, xo)

            dop_min = min(dop_min, float(doppler.min()))
            dop_max = max(dop_max, float(doppler.max()))

            # Obliquity metric, |X| / sqrt(X_c^2 + Z_n^2). This NO LONGER
            # measures a kernel error: the phase-2 kernel consumes the spanwise
            # coordinate, so its S0 is the true |X| (up to the beta weighting).
            # It is kept because it is a monotone proxy for the criticality
            # parameter the kernel actually branches on,
            #     xi = beta*|x2|/S0,
            # since |X|/sqrt(X_c^2 + Z_n^2) = 1/sqrt(1 - (|Y_s|/|X|)^2) and xi
            # rises with |Y_s|/|X|. 1.0 means the observer is in the blade's
            # chord-normal plane (xi = 0, mid-span gust); large values mean it
            # is near the blade's spanwise axis (xi -> 1), where the selected
            # gust approaches critical and the near-cutoff regularisation
            # engages. Reported, never warned on.
            seen = np.hypot(X[:, 0], X[:, 2])
            true = np.linalg.norm(X, axis=1)
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = np.where(seen > 0.0, true / seen, np.inf)
            span_neglect = max(span_neglect, float(np.max(ratio)))

            acc = np.zeros(N)
            for j in range(cfg.n_azimuth):
                # Source frequencies: the WPS and the radiation integral must
                # BOTH be evaluated here, which is exactly what handing this
                # grid to noise_run does.
                omega_src = omega / doppler[j]
                freqs_src = omega_src / (2.0 * np.pi)

                # Undo the kernel's TE offset so its TE-local coordinates come
                # out as the paper-frame (X_c, Y_s, Z_n). See module docstring.
                obs_kernel = np.array([X[j, 0] + _TE_OFFSET_FRAC * strip.chord,
                                       X[j, 1],
                                       X[j, 2]])

                custom = None
                if strip.custom_WPS_func is not None:
                    wps_u, wps_l = strip.custom_WPS_func(omega_src)
                    wps_u = _as_1d_float_array(wps_u, "custom_WPS_func upper")
                    wps_l = _as_1d_float_array(wps_l, "custom_WPS_func lower")
                    if wps_u.size != N or wps_l.size != N:
                        raise ValueError(
                            f"strips[{iStrip}].custom_WPS_func returned arrays of "
                            f"length {wps_u.size}/{wps_l.size}; expected {N} to "
                            f"match the frequency grid"
                        )
                    custom = np.column_stack([wps_u, wps_l])

                nr = noise_run(
                    BL_top, BL_bot,
                    freqs_Hz=freqs_src,
                    observerXYZ=obs_kernel,
                    Re=Re,
                    nu=cfg.nu,
                    chord=strip.chord,
                    span=strip.dr,
                    alphaDeg=0.0,      # every rotation is done above
                    rho=cfg.rho,
                    model=strip.wps_model,
                    custom_WPS=custom,
                )

                ff = nr.FF_spectra[0]
                if not np.all(np.isfinite(ff)):
                    # Never let a non-finite kernel result reach the dB
                    # conversion: NaN > 0 is False, so it would be silently
                    # floored to -200 dB and read as "quiet" (the same silent
                    # blank failure the forward path names acoustic_nan).
                    raise RuntimeError(
                        f"strips[{iStrip}]: noise_run returned a non-finite "
                        f"far-field spectrum at azimuth {j} "
                        f"(gamma = {2 * np.pi * j / cfg.n_azimuth:.4f} rad), "
                        f"observer {iObs}. Check the trailing-edge BL state — "
                        f"in particular Ue > 0 on the radiating surface, which "
                        f"the Amiet stage divides by."
                    )

                # Doppler exponent +2 — the paper's headline correction.
                acc += ff / (doppler[j] ** 2)

            Spp_perStrip[iStrip, iObs] = acc / cfg.n_azimuth

    # Blades are uncorrelated broadband sources -> linear power sum.
    Spp = cfg.B * Spp_perStrip.sum(axis=0)

    with np.errstate(divide="ignore", invalid="ignore"):
        Spp_dB = np.where(Spp > 0.0,
                          10.0 * np.log10(2.0 * np.pi * Spp / _P_REF2),
                          _SPP_DB_FLOOR)

    p2 = np.trapezoid(Spp, omega, axis=1) if hasattr(np, "trapezoid") \
        else np.trapz(Spp, omega, axis=1)
    OASPL = np.where(p2 / _P_REF2 > 1e-30,
                     10.0 * np.log10(np.maximum(p2, 1e-300) / _P_REF2),
                     _OASPL_FLOOR)

    omega_over_Omega_min = (float(omega[0] / abs(cfg.Omega))
                            if cfg.Omega != 0.0 else np.inf)
    if np.isfinite(omega_over_Omega_min) and omega_over_Omega_min < 10.0:
        warnings.warn(
            f"omega_min/Omega = {omega_over_Omega_min:.1f} < 10: the azimuthal "
            f"average assumes omega >> Omega, so the lowest requested "
            f"frequencies are outside the model's validity regime.",
            stacklevel=2,
        )

    # No warning on span_neglect: with the phase-2 general kernel this ratio is
    # a description of the geometry, not a defect. See the comment at its
    # computation above, and CHANGELOG "General oblique-gust Amiet kernel".

    diagnostics = {
        "n_azimuth": cfg.n_azimuth,
        "doppler_min": dop_min,
        "doppler_max": dop_max,
        "omega_over_Omega_min": omega_over_Omega_min,
        # max over (strip, azimuth, observer) of |X| / sqrt(X_c^2 + Z_n^2).
        # Obliquity indicator and a monotone proxy for the kernel's criticality
        # parameter xi = beta*|x2|/S0: 1.0 = observer in the blade's
        # chord-normal plane (xi = 0); large = near its spanwise axis (xi -> 1).
        # Historically this measured the phase-1 mid-span kernel's error
        # (~40*log10(ratio) dB of PSD inflation); it no longer does.
        "spanwise_neglect_ratio": span_neglect,
        "strips": strip_diags,
    }

    return RotorNoiseResult(
        freqs_Hz=freqs,
        Spp=Spp,
        Spp_dB=Spp_dB,
        OASPL_perObs=OASPL,
        diagnostics=diagnostics,
        Spp_perStrip=Spp_perStrip if verbose else None,
    )
