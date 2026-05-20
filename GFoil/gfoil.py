import json
import subprocess
import numpy as np
import os
from .inputs import Aerofoil, Acoustics, OperatingConds, WPSinfo, FwdResult, GradResult

# Try the compiled pybind11 module first (installed into this package dir).
# Falls back to the standalone binaries via subprocess if not available.
try:
    from . import gfoil_cpp
    _USE_BINDINGS = True
except ImportError:
    _USE_BINDINGS = False
    BIN_DIR = os.path.join(os.path.dirname(__file__), "bin")
    EXEC_FWD_codi = os.path.join(BIN_DIR, "GFoil_fwd_codi")
    EXEC_AD       = os.path.join(BIN_DIR, "GFoil_AD")


def _build_input_dict(aerofoil: Aerofoil,
                      operating: OperatingConds,
                      acoustics: Acoustics,
                      returnAllOutputs: bool,
                      alphaDeg: float = None,
                      fromRestart: int = 0) -> dict:
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
        "returnData":    int(returnAllOutputs),
        "ncrit":         float(operating.nCrit),
        "ncrithyst":     float(operating.ncrithyst),
        "Ufac":          float(aerofoil.panelUniformity),
        "TEfac":         float(aerofoil.panelTEspacing),
        "toptrans":      float(operating.transition[0]),
        "bottrans":      float(operating.transition[1]),
        "forcetrans":    force,
        "model":         acoustics.model,
        "aWeighting":    int(acoustics.aWeighting),
        "chord":         float(aerofoil.chord),
        "WPSonly":       0,
    }


def _call_forward(inp: dict, prev_result: "FwdResult" = None) -> "FwdResult":
    """Single forward solve — pybind11 or subprocess. Returns FwdResult."""
    if _USE_BINDINGS:
        if prev_result is not None and prev_result.converged:
            jac_in = {"states": prev_result.states, "turb": prev_result.turb}
            r = gfoil_cpp.run_forward(inp, jac_in)
        else:
            r = gfoil_cpp.run_forward(inp)
        if r["conv"] == 0:
            return FwdResult(converged=False)
        jac = r["jacobian"]
        return FwdResult(
            converged=True,
            CL=r["CL"], CD=r["CD"], CM=r["CM"], OASPL=r["OASPL"],
            states=jac["states"], turb=jac["turb"], stag=jac["stag"],
            RVvals=jac["RVvals"], RVrows=jac["RVrows"], RVcols=jac["RVcols"],
            RVnz=jac["RVnz"],
            ycoords=np.array(inp["ycoords"]),
            alpha=inp["alpha_degrees"],
        )
    else:
        cwd = os.getcwd()
        # subprocess warm-start: write restart.json from prev_result if available
        if prev_result is not None and prev_result.converged:
            rs_out = {
                "states": prev_result.states,
                "turb":   prev_result.turb,
                "stag":   prev_result.stag,
                "RVvals": prev_result.RVvals,
                "RVrows": prev_result.RVrows,
                "RVcols": prev_result.RVcols,
                "RVnz":   prev_result.RVnz,
            }
            with open(os.path.join(cwd, "restart.json"), "w") as f:
                json.dump(rs_out, f)
            inp = dict(inp)
            inp["restart"] = 1
        with open(os.path.join(cwd, "input.json"), "w") as f:
            json.dump(inp, f)
        result = subprocess.run([EXEC_FWD_codi], cwd=cwd, capture_output=True, text=True)
        if result.returncode != 1:
            return FwdResult(converged=False)
        with open(os.path.join(cwd, "out.json")) as f:
            out = json.load(f)
        if out.get("conv", 0) == 0:
            return FwdResult(converged=False)
        with open(os.path.join(cwd, "restart.json")) as f:
            rs = json.load(f)
        return FwdResult(
            converged=True,
            CL=out["CL"], CD=out["CD"], CM=out["CM"], OASPL=out["OASPL"],
            states=rs["states"], turb=rs["turb"], stag=rs["stag"],
            RVvals=rs["RVvals"], RVrows=rs["RVrows"], RVcols=rs["RVcols"],
            RVnz=rs["RVnz"],
            ycoords=np.array(inp["ycoords"]),
            alpha=inp["alpha_degrees"],
        )


def standard_run(aerofoil: Aerofoil,
                 operating: OperatingConds,
                 acoustics: Acoustics,
                 returnAllOutputs: bool = False) -> FwdResult:
    """
    Forward solve with backstepping/continuation on failure.
    Returns FwdResult; result.converged is False if all attempts fail.
    """
    inp = _build_input_dict(aerofoil, operating, acoustics, returnAllOutputs)
    result = _call_forward(inp)
    if result.converged:
        return result

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
        bs_inp = _build_input_dict(aerofoil, operating, acoustics, False,
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
        return FwdResult(converged=False)

    # Step forward toward original alphaDeg, warm-starting each step from the
    # previous converged solution (both pybind11 and subprocess paths).
    print("Starting forward stepping...")
    stepsize     = 0.5
    fwdalf       = tempalf - step_direction * stepsize
    attemptCount = 0
    overallCount = 0
    completed    = False

    while not completed and overallCount <= 6:
        print(f"Trying forward step to: {fwdalf:.2f}")
        fs_inp = _build_input_dict(aerofoil, operating, acoustics, returnAllOutputs,
                                   alphaDeg=fwdalf)
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
            fwdalf += step_direction * (stepsize / (2 ** attemptCount))
        overallCount += 1

    return last_converged if completed else FwdResult(converged=False)


def fwd_run(aerofoil: Aerofoil,
            operating: OperatingConds,
            acoustics: Acoustics,
            returnAllOutputs: bool = False,
            repanel: bool = False) -> FwdResult:
    """
    Run forward solver. Returns FwdResult.
    result.converged is False on failure.
    result.CL, result.CD, result.CM, result.OASPL give the scalar outputs.
    Pass result to grad_run() to compute gradients.
    """
    if repanel:
        result = standard_run(aerofoil, operating, acoustics, returnAllOutputs)
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
            result = standard_run(foil2, operating, acoustics, returnAllOutputs)
            if result.converged:
                return result
        return FwdResult(converged=False)
    else:
        return standard_run(aerofoil, operating, acoustics, returnAllOutputs)


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

    inp = _build_input_dict(aerofoil, operating, acoustics, False)
    jacobian = {
        "states": fwd_result.states,
        "turb":   fwd_result.turb,
        "stag":   fwd_result.stag,
        "RVvals": fwd_result.RVvals,
        "RVrows": fwd_result.RVrows,
        "RVcols": fwd_result.RVcols,
        "RVnz":   fwd_result.RVnz,
    }

    if _USE_BINDINGS:
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
    else:
        cwd = os.getcwd()
        with open(os.path.join(cwd, "input.json"), "w") as f:
            json.dump(inp, f)
        with open(os.path.join(cwd, "restart.json"), "w") as f:
            json.dump(jacobian, f)
        subprocess.run([EXEC_AD], cwd=cwd, capture_output=True, text=True)
        with open(os.path.join(cwd, "ad_gradients.json")) as f:
            g = json.load(f)
        return GradResult(
            converged=True,
            dCL_dy=np.array(g["d cl / d ycoords"]),
            dCD_dy=np.array(g["d cd / d ycoords"]),
            dOASPL_dy=np.array(g["d OASPL / d ycoords"]),
            dCL_dalpha=g["d cl / d alpha"],
            dCD_dalpha=g["d cd / d alpha"],
            dOASPL_dalpha=g["d OASPL / d alpha"],
        )


def WPS_run(data: WPSinfo) -> bool:
    cwd = os.getcwd()
    if _USE_BINDINGS:
        exec_path = os.path.join(os.path.dirname(__file__), "bin", "GFoil_fwd_codi")
    else:
        exec_path = EXEC_FWD_codi
    payload = {
        "Re":        data.Re,   "rho":  data.rho,   "nu":  data.nu,
        "X":         [data.observerXYZ[0]],
        "Y":         [data.observerXYZ[1]],
        "Z":         [data.observerXYZ[2]],
        "S":         data.span, "model": data.model, "chord": data.chord,
        "WPSonly":   1,
        "topdstar":  data.DispThick[0],  "toptheta":  data.MomThick[0],
        "topdelta":  data.BLHeight[0],   "toptauw":   data.wallShear[0],
        "toptaumax": data.maxShear[0],   "topue":     data.edgeVel[0],
        "topdpdx":   data.dpdx[0],
        "botdstar":  data.DispThick[1],  "bottheta":  data.MomThick[1],
        "botdelta":  data.BLHeight[1],   "bottauw":   data.wallShear[1],
        "bottaumax": data.maxShear[1],   "botue":     data.edgeVel[1],
        "botdpdx":   data.dpdx[1],
    }
    with open(os.path.join(cwd, "input.json"), "w") as f:
        json.dump(payload, f, indent=4)
    r = subprocess.run([exec_path], cwd=cwd, capture_output=True, text=True)
    return r.returncode == 1
