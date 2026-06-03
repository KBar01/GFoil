import numpy as np
import os
from .inputs import Aerofoil, Acoustics, OperatingConds, FwdResult, GradResult, VerboseResult

from . import gfoil_cpp


def _build_input_dict(aerofoil: Aerofoil,
                      operating: OperatingConds,
                      acoustics: Acoustics,
                      alphaDeg: float = None,
                      fromRestart: int = 0,
                      verbose: bool = False) -> dict:
    force = 1 if (operating.transition[0] != 1.0 or operating.transition[1] != 1.0) else 0
    alpha = alphaDeg if alphaDeg is not None else float(operating.alpha)
    return {
        "xcoords":       aerofoil.xcoords.tolist(),
        "ycoords":       aerofoil.ycoords.tolist(),
        "alpha_degrees": alpha,
        "Re":            float(operating.Re),
        "Ma":            float(operating.Ma),
        "rho":           float(operating.rho),
        "nu":            float(operating.nu),
        "restart":       fromRestart,
        "sampleTE":      float(acoustics.TESampleLoc),
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
        jac_in = {"states": prev_result.states, "turb": prev_result.turb}
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


def standard_run(aerofoil: Aerofoil,
                 operating: OperatingConds,
                 acoustics: Acoustics,
                 verbose: bool = False) -> FwdResult:
    """
    Forward solve with backstepping/continuation on failure.
    Returns FwdResult; result.converged is False if all attempts fail.
    verbose=True populates result.verbose_data on the final converged solve.
    """
    inp = _build_input_dict(aerofoil, operating, acoustics, verbose=verbose)
    result = _call_forward(inp)
    if result.converged:
        return result

    # Preserve the failure_mode from the initial attempt so it can be returned
    # if all continuation attempts also fail.
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


def fwd_run(aerofoil: Aerofoil,
            operating: OperatingConds,
            acoustics: Acoustics,
            repanel: bool = False,
            verbose: bool = False) -> FwdResult:
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


