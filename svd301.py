import os
import numpy as np
import matplotlib.pyplot as plt

# === CONFIGURATION ===


def doSVD():
    folder = "Smoothed_TEfixed_linear"
    baseline_file = "n0012_sharp.dat"
    #baseline_file = "n63-412blunt.dat"

    # === LOAD BASELINE AIRFOIL ===
    x_baseline, y_baseline = np.loadtxt(baseline_file, unpack=True)
    n_points = len(x_baseline)

    # === LOAD AIRFOIL LIBRARY ===
    Y = []  # list of all y-coordinate vectors

    for fname in os.listdir(folder):
        if not fname.endswith('.dat'):
            continue
        full_path = os.path.join(folder, fname)
        if fname == baseline_file:
            continue  # skip the baseline if it's in the same folder

        try:
            data = np.loadtxt(full_path, skiprows=3)
            x, y = data[:, 0], data[:, 1]

            # Sanity check
            if len(x) != n_points:
                print(f"Skipping {fname} due to point mismatch.")
                continue

            Y.append(y)
        except Exception as e:
            print(f"Error loading {fname}: {e}")

    Y = np.array(Y).T  # shape: (n_points, n_airfoils)

    # === COMPUTE DEFORMATIONS ===
    deltaY = Y - y_baseline[:, np.newaxis]  # shape: (n_points, n_airfoils)

    # === SVD ===
    U, S, Vt = np.linalg.svd(deltaY, full_matrices=False)

    print(S.shape)
    U_scaled = np.zeros_like(U)
    #scales = np.zeros(U.shape[1])

    energy = S / np.sum(S)
    for i in range(U.shape[1]):
        U_scaled[:, i] = U[:, i] * energy[i]

    
    return x_baseline,y_baseline,U_scaled

def doSVD_unWeighted():
    folder = "Smoothed_TEfixed"
    baseline_file = "n0012.dat"
    #baseline_file = "n63-412blunt.dat"

    # === LOAD BASELINE AIRFOIL ===
    x_baseline, y_baseline = np.loadtxt(baseline_file, unpack=True)
    n_points = len(x_baseline)

    # === LOAD AIRFOIL LIBRARY ===
    Y = []  # list of all y-coordinate vectors

    for fname in os.listdir(folder):
        if not fname.endswith('.dat'):
            continue
        full_path = os.path.join(folder, fname)
        if fname == baseline_file:
            continue  # skip the baseline if it's in the same folder

        try:
            data = np.loadtxt(full_path, skiprows=3)
            x, y = data[:, 0], data[:, 1]

            # Sanity check
            if len(x) != n_points:
                print(f"Skipping {fname} due to point mismatch.")
                continue

            Y.append(y)
        except Exception as e:
            print(f"Error loading {fname}: {e}")

    Y = np.array(Y).T  # shape: (n_points, n_airfoils)

    # === COMPUTE DEFORMATIONS ===
    deltaY = Y - y_baseline[:, np.newaxis]  # shape: (n_points, n_airfoils)

    # === SVD ===
    U, S, Vt = np.linalg.svd(deltaY, full_matrices=False)

    return x_baseline,y_baseline,U



