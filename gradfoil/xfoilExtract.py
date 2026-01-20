import json
import subprocess
import numpy as np
import os
#import tempfile
#from importlib.resources import files
import csv


def extract_foil_data(array):

    firstSection = []
    i = 0
    xcurr = array[i,0] 
    while xcurr <= 1.0:
        val = array[i,1]
        firstSection.append(val)
        i += 1 
        xcurr = array[i,0]
    
    #skip over wake section
    xprev = array[i-1,0]
    while xcurr > xprev:
        i += 1
        xcurr = array[i,0]
        xprev = array[i-1,0]

    #extract top theta
    secondSection = []
    while xcurr <= 1.0:
        val = array[i,1]
        secondSection.append(val)
        i += 1 
        if i >= array.shape[0]:
            break
        xcurr = array[i,0]

    return np.array(firstSection),np.array(secondSection)


def extract_ncrit_data(array):

    firstSection = []
    i = 0
    nprev = 0
    ncurr = array[i,1]
    while ncurr>=nprev:
        val = array[i,1]
        firstSection.append(val)
        i += 1 
        ncurr = array[i,1]
        nprev = array[i-1,1]


    #extract top theta
    secondSection = []
    while i < array.shape[0]:
        secondSection.append(array[i,1])
        i += 1

    return np.array(firstSection),np.array(secondSection)


def run_xfoil_get_BL_states(coords, alpha, Re, mach, xfoiPath, ncrit):
    
    """Run XFOIL on a 200-point airfoil coordinate array and extract boundary layer data.

    Args:
        coords (list of tuple): Airfoil surface coordinates.
        alpha (float): Angle of attack in degrees.
        Re (float): Reynolds number (excluding chord).
        mach (float): Mach number.
        xfoil_path (str): Full path to the XFOIL executable.
    """

    # Filenames
    dat_file = 'airfoil.dat'
    xfoil_in = 'xfoil_input.txt'
    files = {
        'n'  : 'bl_n.txt',
        'ct' : 'bl_ct.txt',
        'top': 'bl_top.txt',
        'bot': 'bl_bot.txt',
        'ue' : 'bl_ue.txt'
    }

    

    # Write airfoil.dat file
    with open(dat_file, 'w') as f:
        f.write('Foil\n')
        for x, y in coords:
            f.write(f"{x:.6f} {y:.6f}\n")

    # Write XFOIL input script
    with open(xfoil_in, 'w') as f:
        f.write(f"LOAD {dat_file}\n")
        f.write("OPER\n")
        f.write("VPAR\n")
        f.write("N\n")
        f.write(f"{ncrit}\n") 
        f.write("\n")          # Return from VPAR menu
        f.write(f"VISC {Re:.0f}\n")
        f.write(f"MACH {mach:.4f}\n")
        f.write("ITER 100\n")
        f.write(f"ALFA {alpha}\n")
        f.write("VPLO\n")
        f.write("N\nDUMP\n" + files['n'] + "\n")
        f.write("CT\nDUMP\n" + files['ct'] + "\n")
        f.write("DT\nDUMP\n" + files['top'] + "\n")
        f.write("DB\nDUMP\n" + files['bot'] + "\n")
        f.write("UE\nDUMP\n" + files['ue'] + "\n\nQUIT\n")

    # Run XFOIL
    #/mnt/c/Users/pa20830/XFOIL6.99/xfoil.exe
    result = subprocess.run([xfoiPath],input=open(xfoil_in).read(),text=True,stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)

    # Parse CL from output
    cl_value = None
    for line in result.stdout.splitlines():
        if "CL =" in line:
            try:
                cl_value = float(line.split("CL =")[1].split()[0])
                break
            except (IndexError, ValueError):
                continue
    
    if cl_value is None:
        print("Warning: CL value not found in XFOIL output.")


    def load_column_data(file):
        data = []
        with open(file, 'r') as f:
            for line in f:
                try:
                    a, b = map(float, line.strip().split())
                    data.append((a, b))
                except:
                    continue
        return np.array(data)

    # Load all BL data
    n_data   = load_column_data(files['n'])
    ct_data  = load_column_data(files['ct'])
    top_data = load_column_data(files['top'])
    bot_data = load_column_data(files['bot'])
    ue_data  = load_column_data(files['ue'])

    
    topDelta,topTheta = extract_foil_data(top_data)
    botDelta,botTheta = extract_foil_data(bot_data)
    topN,botN   = extract_ncrit_data(n_data)
    topCt,botCt = extract_foil_data(ct_data)
    topUe,botUe = extract_foil_data(ue_data)
   
    th = np.concatenate((np.flip(botTheta), topTheta))
    ds = np.concatenate((np.flip(botDelta),topDelta))
    sa = np.concatenate((np.flip(botCt), np.flip(botN), topN, topCt))
    ue = np.concatenate((np.flip(botUe), topUe))

    # uncorrect if mach used
    if mach > 0.0:
        ktb = np.sqrt(1-mach**2)
        l= (mach**2) / (1+(ktb))**2 
        uk = ue

        a = uk*l / 1 
        b = 1-l
        c = -uk

        ue_sols = np.array([(-b + np.sqrt((b**2)- 4*a*c))/(2*a),(-b - np.sqrt((b**2)- 4*a*c))/(2*a)])
        ue = ue_sols[0,:]
    
    out = np.vstack((th, ds, sa, ue))  # shape (4, N)

    # Create integer mask with 1s for Ct arrays and 0s for N arrays
    turb = np.concatenate((
        np.ones(len(botCt), dtype=int),   # 1s for flipped botCt
        np.zeros(len(botN), dtype=int),   # 0s for flipped botN
        np.zeros(len(topN), dtype=int),   # 0s for topN
        np.ones(len(topCt), dtype=int)    # 1s for topCt
    ))

    # Optionally clean up temporary files
    for f in [dat_file, xfoil_in] + list(files.values()):
        if os.path.exists(f):
            os.remove(f)

    return out,turb,cl_value  # shape (4, N)

def xfoil_start_run(EXEC_FWD,xcoords,ycoords,alphaDeg,Re,Ma,sampleTE,X,Y,Z,S,xfoilPath,Uinf,custUinf,trackCLs,returnFoilCps,ncrit,Ufac,TEfac,toptrans,bottrans,force,custNcrits,topncrit,botncrit,breakloopInt):
    

    cwd = os.getcwd()
    returnFoilCps = 0
    in_json_path = os.path.join(cwd, "input.json")
    data = {
        "xcoords":       xcoords,
        "ycoords":       ycoords,
        "alpha_degrees": alphaDeg,
        "Re":            Re,
        "Ma":            Ma,
        "restart":       0,
        "xfoilstart":    0,
        "xfoilgetpoints":1,
        "sampleTE":      sampleTE,
        "X":             X,
        "Y":             Y,
        "Z":             Z,
        "S":             S,
        "Uinf":          Uinf,
        "custUinf":      custUinf,
        "returnData":    returnFoilCps,
        "ncrit":         ncrit,
        "Ufac":          Ufac,
        "TEfac":         TEfac,
        "toptrans":      toptrans,
        "bottrans":      bottrans,
        "forcetrans":    force,
        "custncrits":    custNcrits,
        "topncrit":      topncrit,
        "botncrit":      botncrit,
        "breakloop":     breakloopInt
    }

    # Write JSON input
    with open(in_json_path, "w") as f:
        json.dump(data, f)

    # Run the executable to get inner foil coords
    initResult = subprocess.run([EXEC_FWD],cwd=os.getcwd(), capture_output=True, text=True)
    initConvergence = initResult.returncode

    foil_json_path = os.path.join(cwd, "innerfoilNodes.json")
    with open(foil_json_path, "r") as f:
        data = json.load(f)
    
    nodes = np.array(data["points"])
    nodes = nodes.reshape((2, 200), order='F')

    # Run xfoil with the foil for converged states
    innerFoil = nodes.T

    # convert chord Re to xfoil Re
    chord = np.max(xcoords) - np.min(xcoords)
    ReXfoil = Re/chord

    xfoilStates,xfoilturb,xfoilCL = run_xfoil_get_BL_states(innerFoil,alphaDeg,ReXfoil,Ma,xfoilPath,ncrit)
    xfoilCL /= chord
    # save them to json
    xfoil_json_path = os.path.join(cwd, "xfoilstates.json")
    data = {
        "states":(xfoilStates.flatten(order='F')).tolist(),
        "turb": xfoilturb.tolist()
    }
    with open(xfoil_json_path, "w") as f:
        json.dump(data, f)

    with open(in_json_path, "r") as f:
        data = json.load(f)
    data["xfoilgetpoints"] = 0 
    data["xfoilstart"] = 1      
    with open(in_json_path, "w") as f:
        json.dump(data, f, indent=4)
    
    Result = subprocess.run([EXEC_FWD],cwd=os.getcwd(), capture_output=True, text=True)
    converged = Result.returncode


    if trackCLs:
        # Get CL value from JSON output
        out_json_path = os.path.join(cwd, "out.json")
        with open(out_json_path, "r") as f:
            data = json.load(f)
        cppCL = data["CL"]

        # CSV file path
        cl_log_path = os.path.join(cwd, "trackCLs.csv")

        # Check if file exists
        file_exists = os.path.isfile(cl_log_path)

        # Open the CSV in append mode
        with open(cl_log_path, mode='a', newline='') as csvfile:
            writer = csv.writer(csvfile)

            # Write header if file doesn't exist
            if not file_exists:
                writer.writerow(["xfoilCL", "cppCL"])

            # Write the current CL values
            writer.writerow([xfoilCL, cppCL])



    return converged