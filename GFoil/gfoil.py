import json
import subprocess
import numpy as np
import os
from dataclasses import dataclass, field
from .inputs import Aerofoil, Acoustics, OperatingConds, WPSinfo, customSpectrainfo


# Dynamically locate the installed executable path (in gradfoil/bin/)
BIN_DIR = os.path.join(os.path.dirname(__file__), "bin")
EXEC_FWD_codi = os.path.join(BIN_DIR, "GFoil_fwd_codi")
EXEC_AD = os.path.join(BIN_DIR, "GFoil_AD")




def standard_run(aerofoil: Aerofoil,
    operating: OperatingConds,
    acoustics: Acoustics,
    returnAllOutputs: bool = False):
    
    cwd = os.getcwd()
    in_json_path = os.path.join(cwd, "input.json")
    
    if (operating.transition[0] != 1.0) or (operating.transition[1] != 1.0):
        force = 1
    else:
        force = 0
    
    data = {
        "xcoords":       aerofoil.xcoords.tolist(),
        "ycoords":       aerofoil.ycoords.tolist(),
        "alpha_degrees": operating.alpha,
        "Re":            operating.Re,
        "Ma":            operating.Ma,
        "rho":           operating.rho,
        "nu":            operating.nu,
        "restart":        0,
        "sampleTE":      acoustics.TESampleLoc,
        "X":             acoustics.observerXYZ[0],
        "Y":             acoustics.observerXYZ[1],
        "Z":             acoustics.observerXYZ[2],
        "S":             aerofoil.span,
        "returnData":    returnAllOutputs,
        "ncrit":         operating.nCrit,
        "Ufac":          aerofoil.panelUniformity,
        "TEfac":         aerofoil.panelTEspacing,
        "toptrans":      operating.transition[0],
        "bottrans":      operating.transition[1],
        "forcetrans":    force,
        "model":  acoustics.model,
        "chord": aerofoil.chord,
        "WPSonly":0,
        "custSpectra":   0
    }

    # Write JSON input file
    with open(in_json_path, "w") as f:
        json.dump(data, f)

    initResult = subprocess.run([EXEC_FWD_codi],cwd=os.getcwd(), capture_output=True, text=True)
    initConvergence = initResult.returncode
    if initConvergence==1:
        return True

    print("Initial run failed. Starting backstepping ...")
    
    max_back_steps = 5
    stepsize = 1.0
    small_step = 0.5
    back_converged = False
    completed = False
    
    alphaDeg = float(operating.alpha)

    # Determine stepping direction based on sign of alphaDeg
    if alphaDeg >= 0:
        step_direction = -1
    else:
        step_direction = 1
    
    # take a step of 0.5degrees, set limit of backstep to start + 5.0
    tempalf = np.round(alphaDeg, decimals=1) + (step_direction * stepsize)
    min_alpha = alphaDeg + (step_direction * 5.0)
    
    for i in range(max_back_steps):
        
        if abs(tempalf) < 2.0:
            small_step = 0.1

        # Check if minimum/maximum alpha reached
        if (step_direction < 0 and tempalf < min_alpha) or (step_direction > 0 and tempalf > min_alpha):
            print("Minimum backstep AoA reached. Cannot continue.")
            break

        # Modify input JSON
        with open(in_json_path, "r") as f:
            data = json.load(f)

        data["alpha_degrees"] = tempalf
        data["restart"] = 0  # fresh run
        with open(in_json_path, "w") as f:
            json.dump(data, f)

        result = subprocess.run([EXEC_FWD_codi], cwd=os.getcwd(), capture_output=True, text=True)
        if result.returncode == 1:
            print(f"Backstep converged at {tempalf}")
            back_converged = True
            break

        tempalf += (step_direction * small_step)

    if not back_converged:
        print("Backstepping failed. No converged base solution.")
        return False
    
    
    # Step forward toward original alphaDeg using restart
    print("Starting forward stepping...")
    stepsize = 0.5
    fwdalf = tempalf - (step_direction * stepsize)
    attemptCount = 0
    overallCount = 0
    max_attempts = 6

    while (not completed) and (overallCount <= max_attempts):
        
        print(f"Trying forward step to: {fwdalf:.2f}")

        with open(in_json_path, "r") as f:
            data = json.load(f)

        data["restart"] = 1
        data["alpha_degrees"] = fwdalf
        with open(in_json_path, "w") as f:
            json.dump(data, f, indent=4)

        result = subprocess.run([EXEC_FWD_codi], cwd=os.getcwd(), capture_output=True, text=True)
        converged = result.returncode == 1

        if converged:
            if abs(fwdalf - alphaDeg) < 1e-3:
                completed = True
                break
            else:
                
                nextStep = fwdalf - (step_direction*stepsize)
                diff = alphaDeg - fwdalf

                if abs(diff) < abs(nextStep):
                    fwdalf = alphaDeg
                else:
                    fwdalf = nextStep

                attemptCount = 0
        else:
            attemptCount += 1
            if attemptCount > 6:
                print("Forward stepping failed repeatedly.")
                break

            fwdalf += step_direction * (stepsize / (2 ** attemptCount))

        overallCount += 1

    return completed
    

def fwd_run(
    aerofoil: Aerofoil,
    operating: OperatingConds,
    acoustics: Acoustics,
    returnAllOutputs: bool = False,
    repanel: bool = False
):

    if repanel:
        success = standard_run(aerofoil,operating,acoustics,returnAllOutputs)

        if success:
            return success
        else:
            count = 1 
            for uf, tef in [(1.8,0.1), (2.1,0.09), (2.6,0.09), (1.0,0.09), (1.0,1.1), (1.5,0.09)]:
                
                foil_panelChange = Aerofoil(
                    xcoords=aerofoil.xcoords.copy(),
                    ycoords=aerofoil.ycoords.copy(),
                    chord=aerofoil.chord,
                    span=aerofoil.span,
                    panelUniformity=uf,
                    panelTEspacing=tef,
                )
                print('trying different panel distribution ('+str(count)+'/6)')
                success = standard_run(foil_panelChange,operating,acoustics,returnAllOutputs)
                if success:
                    break
                count +=1
        
            return success
    
    else:
        success = standard_run(aerofoil,operating,acoustics,returnAllOutputs)
        return success


def grad_run():
    # Run the AD version of the code, using known solution from fwd run
    result = subprocess.run([EXEC_AD],cwd=os.getcwd(), capture_output=True, text=True)

def WPS_run(data: WPSinfo):
    cwd = os.getcwd()
    in_json_path = os.path.join(cwd, "input.json")
    data = {
        "Re":            data.Re,
        "rho":           data.rho,
        "nu":            data.nu,
        "X":             data.observerXYZ[0],
        "Y":             data.observerXYZ[1],
        "Z":             data.observerXYZ[2],
        "S":             data.span,
        "model":         data.model,
        "chord":         data.chord,
        "WPSonly":       1,
        "custSpectra":   0,
        
        # -------- TOP boundary layer states --------
        "topdstar":   data.DispThick[0],
        "toptheta":   data.MomThick[0],
        "topdelta":   data.BLHeight[0],
        "toptauw":    data.wallShear[0],
        "toptaumax":  data.maxShear[0],
        "topue":      data.edgeVel[0],
        "topdpdx":    data.dpdx[0],

        # -------- BOTTOM boundary layer states --------
        "botdstar":   data.DispThick[1],
        "bottheta":   data.MomThick[1],
        "botdelta":   data.BLHeight[1],
        "bottauw":    data.wallShear[1],
        "bottaumax":  data.maxShear[1],
        "botue":      data.edgeVel[1],
        "botdpdx":    data.dpdx[1]
    }

    # Write JSON input file
    with open(in_json_path, "w") as f:
        json.dump(data, f,indent=4)

    # Run the executable for first time, no restarting, use codi version to ensure output match to AD version of code
    initResult = subprocess.run([EXEC_FWD_codi],cwd=os.getcwd(), capture_output=True, text=True)
    initConvergence = initResult.returncode
    if initConvergence==1:
        return True
    

def custom_spectra_run(data: customSpectrainfo):
    cwd = os.getcwd()
    in_json_path = os.path.join(cwd, "input.json")
    data = {
        
        "custWPSLower":  data.WPSLower.tolist(),
        "custWPSUpper":  data.WPSUpper.tolist(),
        "custFreq":     data.Freq.tolist(),
        "Re":            data.Re,
        "rho":           data.rho,
        "nu":            data.nu,
        "X":             data.observerXYZ[0],
        "Y":             data.observerXYZ[1],
        "Z":             data.observerXYZ[2],
        "S":             data.span,
        "chord":         data.chord,
        "WPSonly":       0,
        "custSpectra":   1
    }

    # Write JSON input file
    with open(in_json_path, "w") as f:
        json.dump(data, f,indent=4)

    # Run the executable for first time, no restarting, use codi version to ensure output match to AD version of code
    initResult = subprocess.run([EXEC_FWD_codi],cwd=os.getcwd(), capture_output=True, text=True)
    initConvergence = initResult.returncode
    if initConvergence==1:
        return True