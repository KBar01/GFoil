import numpy as np
import os
from .inputs import (Aerofoil, Acoustics, OperatingConds, FwdResult, GradResult, VerboseResult, NoiseResult, _as_1d_float_array, _as_float_array)
from . import gfoil_cpp


def _build_input_dict(aerofoil: Aerofoil, operating: OperatingConds, acoustics: Acoustics, alphaDeg: float = None, fromRestart: int = 0, verbose: bool = False) -> dict:
    
    
    force = 1 if (operating.transition[0] != 1.0 or operating.transition[1] != 1.0) else 0
    alpha = alphaDeg if alphaDeg is not None else float(operating.alpha)
    # TESampleLoc is either a scalar x/c (single-point sample) or a length-2
    # [x_lo, x_hi] window (BL-averaged). Pass both endpoints; the C++ side treats
    # sampleTE_hi <= sampleTE as the scalar (degenerate-window) path.
    te = acoustics.TESampleLoc
    if isinstance(te, (list, tuple, np.ndarray)) and len(te) == 2:
        te_lo, te_hi = float(te[0]), float(te[1])
    else:
        te_lo = te_hi = float(te)
    return {
        "xcoords":       aerofoil.xcoords.tolist(),
        "ycoords":       aerofoil.ycoords.tolist(),
        "alpha_degrees": alpha,
        "Re":            float(operating.Re),
        "Ma":            float(operating.Ma),
        "rho":           float(operating.rho),
        "nu":            float(operating.nu),
        "restart":       fromRestart,
        "sampleTE":      te_lo,
        "sampleTE_hi":   te_hi,
        "X":             acoustics.observerXYZ[:, 0].tolist(),
        "Y":             acoustics.observerXYZ[:, 1].tolist(),
        "Z":             acoustics.observerXYZ[:, 2].tolist(),
        "S":             float(aerofoil.span),
        "ncrit":         float(operating.nCrit),
        "rtol":          float(operating.rtol),
        "Ufac":          float(aerofoil.panelUniformity),
        "TEfac":         float(aerofoil.panelTEspacing),
        "toptrans":      float(operating.transition[0]),
        "bottrans":      float(operating.transition[1]),
        "forcetrans":    force,
        "model":         acoustics.model,
        "aWeighting":    int(acoustics.aWeighting),
        "chord":         float(aerofoil.chord),
        "f_min":         float(acoustics.f_min),
        "f_max":         float(acoustics.f_max),
        "verbose":       verbose,
    }


def _call_forward(inp: dict, prev_result: "FwdResult" = None) -> "FwdResult":
    """Single forward solve via pybind11. Returns FwdResult."""
    if prev_result is not None and prev_result.converged:
        # Pass the donor's converged stagnation bracket too, so the C++ warm path
        # can seed stagpoint_move from it and reproduce the donor configuration
        # (avoids the one-node stag reindex that perturbs the warm entry state).
        jac_in = {"states": prev_result.states, "turb": prev_result.turb,
                  "stag": prev_result.stag}
        r = gfoil_cpp.run_forward(inp, jac_in)
    else:
        r = gfoil_cpp.run_forward(inp)

    if r["conv"] == 0:
        return FwdResult(
            converged=False,
            failure_mode=r.get("failure_mode", ""),
            newton_iterations=r.get("newton_iterations", 0),
        )

    jac = r["jacobian"]

    verb = None
    if r.get("innerFoilX") is not None:
        try:
            NS   = len(r["freq_Hz"])
            nObs = r["nObs"]
            verb = VerboseResult(
                x          = np.array(r["innerFoilX"]),
                y          = np.array(r["innerFoilY"]),
                Cp         = np.array(r["Cp_dist"]),
                delta_star = np.array(r["delta_star"]),
                theta      = np.array(r["theta"]),
                tau_wall   = np.array(r["tau_wall"]),
                tau_max    = np.array(r["tau_max"]),
                Ue         = np.array(r["Ue"]),
                dpdx       = np.array(r["dpdx"]),
                is_turb    = np.array(r["is_turb"], dtype=bool),
                topTransX  = float(r["topTransX"]),
                botTransX  = float(r["botTransX"]),
                BL_top     = np.array(r["BL_top"]),
                BL_bot     = np.array(r["BL_bot"]),
                freq_Hz    = np.array(r["freq_Hz"]),
                WPS_upper  = np.array(r["WPS_upper"]),
                WPS_lower  = np.array(r["WPS_lower"]),
                FF_spectra = np.array(r["FF_spectra"]).reshape(nObs, NS),
                OASPL_perObs   = np.array(r["OASPL_perObs"]),
                obsXYZ_TElocal = np.array(r["obsXYZ_TElocal"]).reshape(nObs, 3),
            )
        except (KeyError, TypeError, ValueError) as e:
            raise RuntimeError(
                f"verbose output from C++ solver is incomplete or malformed: {e}"
            ) from e

    return FwdResult(
        converged=True,
        CL=r["CL"], CD=r["CD"], CM=r["CM"], OASPL=r["OASPL"],
        states=jac["states"], turb=jac["turb"], stag=jac["stag"],
        RVvals=jac["RVvals"], RVrows=jac["RVrows"], RVcols=jac["RVcols"],
        RVnz=jac["RVnz"],
        ycoords=np.array(inp["ycoords"]),
        alpha=inp["alpha_degrees"],
        verbose_data=verb,
        failure_mode=r.get("failure_mode", ""),
        newton_iterations=r.get("newton_iterations", 0),
    )


def _is_stale_warm_accept(result: "FwdResult",
                          call_alpha: float,
                          donor_alpha: float) -> bool:
    """Detect a warm-restart stale iteration-0 accept.

    run_forward(inp, restart) can return converged with newton_iterations==0 and
    the DONOR's solution unchanged when warm-started across an alpha change:
    solve_coupled's entry convergence test only covers the BL-station rows, so an
    alpha shift (which enters through the ue-coupling rows) is invisible at
    iteration 0 and the donor state is accepted verbatim. Within standard_run
    geometry/Re/nCrit are invariant, so alpha is the only varying input — an
    it==0 accept at a CHANGED alpha is therefore always stale. Keys on converged
    AND it==0 AND changed alpha together, so a genuine 0-iteration accept at an
    identical alpha (true same-alpha restart) is NOT rejected.
    """
    return (result.converged and result.newton_iterations == 0
            and abs(call_alpha - donor_alpha) > 1e-12)


def standard_run(aerofoil: Aerofoil, operating: OperatingConds, acoustics: Acoustics, verbose: bool = False):
    """
    Forward solve with backstepping when it fails intial convergence.
    """
    inp = _build_input_dict(aerofoil, operating, acoustics, verbose=verbose)
    result = _call_forward(inp)
    if result.converged:
        return result

    initial_failure_mode = result.failure_mode

    print("Initial run failed. Starting backstepping ...")

    alphaDeg       = float(operating.alpha)
    step_direction = -1 if alphaDeg >= 0 else 1
    tempalf        = round(alphaDeg, 1) + step_direction * 1.0
    min_alpha      = alphaDeg + step_direction * 5.0
    small_step     = 0.5
    back_converged = False
    last_converged = None

    for _ in range(5):
        if abs(tempalf) < 2.0:
            small_step = 0.1
        if (step_direction < 0 and tempalf < min_alpha) or \
           (step_direction > 0 and tempalf > min_alpha):
            print("Minimum backstep AoA reached. Cannot continue.")
            break
        bs_inp = _build_input_dict(aerofoil, operating, acoustics,
                                   alphaDeg=tempalf, fromRestart=0)
        r = _call_forward(bs_inp)
        if r.converged:
            print(f"Backstep converged at {tempalf}")
            back_converged = True
            last_converged = r
            break
        tempalf += step_direction * small_step

    if not back_converged:
        print("Backstepping failed. No converged base solution.")
        return FwdResult(converged=False, failure_mode=initial_failure_mode)

    # Step forward toward original alphaDeg, warm-starting each step from the
    # previous converged solution.
    print("Starting forward stepping...")
    stepsize     = 0.5
    fwdalf       = tempalf - step_direction * stepsize
    attemptCount = 0
    overallCount = 0
    completed    = False

    while not completed and overallCount <= 10:
        print(f"Trying forward step to: {fwdalf:.3f}")
        is_final = abs(fwdalf - alphaDeg) < 1e-3
        fs_inp = _build_input_dict(aerofoil, operating, acoustics,
                                   alphaDeg=fwdalf,
                                   verbose=verbose if is_final else False)
        r = _call_forward(fs_inp, prev_result=last_converged)

        # Warm-restart stale-accept guard (see _is_stale_warm_accept): reject an
        # it==0 accept that merely echoes the donor at a changed alpha and retry
        # the SAME alpha cold. If the cold retry converges it flows through the
        # stepping logic normally; if it fails, the existing step-shrink/failure
        # logic below handles it unchanged.
        if _is_stale_warm_accept(r, fwdalf, float(last_converged.alpha)):
            print(f"[gfoil] warm-restart stale accept rejected at "
                  f"alpha={fwdalf:.4f} (donor {float(last_converged.alpha):.4f}); "
                  f"cold retry")
            r = _call_forward(fs_inp)   # cold solve at the same alpha

        if r.converged:
            last_converged = r
            if abs(fwdalf - alphaDeg) < 1e-3:
                completed = True
                break
            diff     = alphaDeg - fwdalf
            nextStep = fwdalf - step_direction * stepsize
            fwdalf   = alphaDeg if abs(diff) < abs(stepsize) else nextStep
            attemptCount = 0
        else:
            attemptCount += 1
            if attemptCount > 6:
                print("Forward stepping failed repeatedly.")
                break
            last_good = float(last_converged.alpha)
            fwdalf = last_good + (fwdalf - last_good) * 0.5
        overallCount += 1

    if completed:
        return last_converged
    return FwdResult(converged=False, failure_mode=initial_failure_mode)


def fwd_run(aerofoil: Aerofoil, operating: OperatingConds, acoustics: Acoustics, repanel: bool = False, verbose: bool = False):
    """
    Run forward solver. Returns FwdResult.
    result.converged is False on failure.
    result.CL, result.CD, result.CM, result.OASPL give the scalar outputs.
    Pass result to grad_run() to compute gradients.
    If verbose=True, result.verbose_data is populated with per-node aero
    and acoustic spectral data.
    """
    if repanel:
        result = standard_run(aerofoil, operating, acoustics, verbose=verbose)
        if result.converged:
            return result
        for count, (uf, tef) in enumerate(
                [(1.8, 0.1), (2.1, 0.09), (2.6, 0.09),
                 (1.0, 0.09), (1.0, 1.1), (1.5, 0.09)], 1):
            foil2 = Aerofoil(
                xcoords=aerofoil.xcoords.copy(),
                ycoords=aerofoil.ycoords.copy(),
                chord=aerofoil.chord,
                span=aerofoil.span,
                panelUniformity=uf,
                panelTEspacing=tef,
            )
            print(f"Trying different panel distribution ({count}/6)")
            result = standard_run(foil2, operating, acoustics, verbose=verbose)
            if result.converged:
                return result
        return FwdResult(converged=False)
    else:
        return standard_run(aerofoil, operating, acoustics, verbose=verbose)


def noise_run(BL_top,
              BL_bot,
              freqs_Hz,
              observerXYZ,
              Re: float,
              nu: float,
              chord: float,
              span: float,
              alphaDeg: float,
              rho: float = 1.225,
              Ma: float = 0.0,
              model: str = "kam",
              custom_WPS=None) -> NoiseResult:
    """Acoustics-only forward run: WPS models + Amiet, with NO aero solve.

    Supply trailing-edge boundary-layer states directly (or a custom
    wall-pressure spectrum) and get back the raw linear wall-pressure spectra
    and far-field PSD. This path is never differentiated.

    Parameters
    ----------
    BL_top, BL_bot : array-like, shape (7,)
        Upper / lower surface trailing-edge BL states, ordered
        [theta, delta_star, tau_max, Ue, dpdx, tau_wall, delta99] — exactly the
        ordering of FwdResult.verbose_data.BL_top / .BL_bot. A surface with
        tau_max <= 0 (fully-laminar TE) is skipped (its WPS column is zeros).
    freqs_Hz : array-like, shape (N,)
        Frequencies [Hz], ANY length, ANY spacing.
    observerXYZ : array-like, shape (N,3) or (3,)
        Observer position(s) in the global frame, origin at quarter-chord
        (same convention as Acoustics.observerXYZ).
    Re, nu, chord, span, alphaDeg : float
        Uinf is derived as Re*nu/chord; alphaDeg drives the global->TE-local
        observer rotation.
    rho, Ma : float
        Density and Mach (Ma is accepted for API symmetry; the Amiet Mach
        follows the Uinf/340 convention of the forward path).
    model : str
        WPS model key ('roz','goo','lee','kam','tno'). Ignored if custom_WPS
        is given.
    custom_WPS : array-like, shape (N,2), optional
        Columns [upper, lower] wall-pressure spectra [Pa^2/omega]. If given, the
        BL->WPS path is skipped for both surfaces and these are used directly.
        An all-zero column skips that surface's far-field contribution.

    Returns
    -------
    NoiseResult
        freqs_Hz (N,), WPS_upper (N,), WPS_lower (N,), FF_spectra (nObs,N),
        obsXYZ_TElocal (nObs,3). Spectra are raw linear Pa^2/omega.
    """
    BL_top = _as_1d_float_array(BL_top, "BL_top")
    BL_bot = _as_1d_float_array(BL_bot, "BL_bot")
    if BL_top.size != 7 or BL_bot.size != 7:
        raise ValueError(
            f"BL_top/BL_bot must have length 7 "
            f"[theta, delta_star, tau_max, Ue, dpdx, tau_wall, delta99]; "
            f"got {BL_top.size} and {BL_bot.size}"
        )

    freqs = np.asarray(freqs_Hz, dtype=float).ravel()
    N = freqs.size
    if N < 1:
        raise ValueError("freqs_Hz must be non-empty")

    obs = _as_float_array(observerXYZ, "observerXYZ").astype(float)
    if obs.ndim == 1 and obs.size == 3:
        obs = obs.reshape(1, 3)
    elif obs.ndim == 2 and obs.shape[1] == 3:
        pass
    else:
        raise ValueError(f"observerXYZ must be shape (3,) or (N,3), got {obs.shape}")
    nObs = obs.shape[0]

    inp = {
        "alphaDeg": float(alphaDeg),
        "Re":       float(Re),
        "rho":      float(rho),
        "nu":       float(nu),
        "Ma":       float(Ma),
        "chord":    float(chord),
        "span":     float(span),
        "X":        obs[:, 0].tolist(),
        "Y":        obs[:, 1].tolist(),
        "Z":        obs[:, 2].tolist(),
        # 7 BL quantities as [upper, lower] pairs
        "theta":     [float(BL_top[0]), float(BL_bot[0])],
        "deltaStar": [float(BL_top[1]), float(BL_bot[1])],
        "tauMax":    [float(BL_top[2]), float(BL_bot[2])],
        "Ue":        [float(BL_top[3]), float(BL_bot[3])],
        "dpdx":      [float(BL_top[4]), float(BL_bot[4])],
        "tauWall":   [float(BL_top[5]), float(BL_bot[5])],
        "delta99":   [float(BL_top[6]), float(BL_bot[6])],
        "freqs_Hz":  freqs.tolist(),
        "model":     str(model),
    }

    if custom_WPS is not None:
        cw = np.asarray(custom_WPS, dtype=float)
        if cw.shape != (N, 2):
            raise ValueError(
                f"custom_WPS must have shape (N,2)=({N},2) matching len(freqs_Hz); "
                f"got {cw.shape}"
            )
        inp["custom_WPS"] = cw.tolist()

    r = gfoil_cpp.noise_run(inp)

    return NoiseResult(
        freqs_Hz       = np.array(r["freqs_Hz"]),
        WPS_upper      = np.array(r["WPS_upper"]),
        WPS_lower      = np.array(r["WPS_lower"]),
        FF_spectra     = np.array(r["FF_spectra"]).reshape(nObs, N),
        obsXYZ_TElocal = np.array(r["obsXYZ_TElocal"]).reshape(nObs, 3),
    )


def grad_run(fwd_result: FwdResult,
             aerofoil: Aerofoil,
             operating: OperatingConds,
             acoustics: Acoustics) -> GradResult:
    """
    Run AD pass using the Jacobian state from fwd_run.
    Returns GradResult with dCL_dy, dCD_dy, dOASPL_dy arrays and alpha scalars.
    """
    if not fwd_result.converged:
        return GradResult(converged=False)

    inp = _build_input_dict(aerofoil, operating, acoustics)
    jacobian = {
        "states": fwd_result.states,
        "turb":   fwd_result.turb,
        "stag":   fwd_result.stag,
        "RVvals": fwd_result.RVvals,
        "RVrows": fwd_result.RVrows,
        "RVcols": fwd_result.RVcols,
        "RVnz":   fwd_result.RVnz,
    }

    g = gfoil_cpp.run_AD(inp, jacobian)
    return GradResult(
        converged=True,
        dCL_dy=np.array(g["dCL_dy"]),
        dCD_dy=np.array(g["dCD_dy"]),
        dOASPL_dy=np.array(g["dOASPL_dy"]),
        dCL_dalpha=g["dCL_dalpha"],
        dCD_dalpha=g["dCD_dalpha"],
        dOASPL_dalpha=g["dOASPL_dalpha"],
    )