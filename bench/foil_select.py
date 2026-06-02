#!/usr/bin/env python3
"""
Geometry parsing, classification and stratified subset selection for the
GFoil cold-start convergence sweep (Phase 0 diagnostic harness).

Foil files are Tecplot-style: a 3-line header (TITLE / VARIABLES / ZONE) then
two whitespace-delimited columns (x, z), ordered TE -> LE -> TE.  We are lenient
about the header: any leading line that does not parse as two floats is skipped.

Classification is geometry-driven (max thickness, max camber, trailing-edge
reflex) with name hints used only to break ties / bias toward recognisable
representatives.  Five geometry classes are targeted:

    thin_low_re      thin glider / RC (AG, SD, SG ...), thickness < ~9%
    thick_symmetric  near-symmetric, thickness >= ~8%   (NACA 00xx ...)
    cambered_gp      moderate camber general purpose     (NACA 2412, Clark-Y ...)
    high_camber      high-lift / high-camber             (Eppler, Selig S ...)
    reflexed_6series 6-series or reflexed camber line
"""

import re
import numpy as np
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
FOIL_DIR = REPO_ROOT / "Smoothed_TEfixed_linear"

CLASSES = [
    "thin_low_re",
    "thick_symmetric",
    "cambered_gp",
    "high_camber",
    "reflexed_6series",
]


# --------------------------------------------------------------------------- #
# Loading                                                                      #
# --------------------------------------------------------------------------- #
def load_coords(path: Path):
    """Return (x, z) float arrays from a Tecplot-style .dat file.

    Header-tolerant: keeps only lines that parse as exactly two floats.
    """
    xs, zs = [], []
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) != 2:
                continue
            try:
                x = float(parts[0])
                z = float(parts[1])
            except ValueError:
                continue
            xs.append(x)
            zs.append(z)
    return np.asarray(xs), np.asarray(zs)


# --------------------------------------------------------------------------- #
# Geometry                                                                     #
# --------------------------------------------------------------------------- #
def _split_surfaces(x, z):
    """Split a TE->LE->TE contour into upper/lower surfaces interpolated onto a
    common cosine-spaced chordwise grid.  Returns (xq, z_upper, z_lower)."""
    le = int(np.argmin(x))
    x1, z1 = x[: le + 1], z[: le + 1]   # one surface (x decreasing 1 -> 0)
    x2, z2 = x[le:], z[le:]             # other surface (x increasing 0 -> 1)

    # Sort each surface by ascending x for interpolation.
    o1 = np.argsort(x1)
    o2 = np.argsort(x2)
    x1, z1 = x1[o1], z1[o1]
    x2, z2 = x2[o2], z2[o2]

    xq = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, 100)))  # 0..1 cosine
    z1q = np.interp(xq, x1, z1)
    z2q = np.interp(xq, x2, z2)

    # Designate the surface with the larger mean z as "upper".
    if np.mean(z1q) >= np.mean(z2q):
        return xq, z1q, z2q
    return xq, z2q, z1q


def geometry(x, z):
    """Return dict of geometric descriptors: t_max, camber_max (signed by the
    location of peak |camber|), reflex flag, te_camber."""
    xq, zu, zl = _split_surfaces(x, z)
    thickness = zu - zl
    camber = 0.5 * (zu + zl)

    t_max = float(np.max(thickness))
    i_cam = int(np.argmax(np.abs(camber)))
    camber_peak = float(camber[i_cam])

    # Reflex: meaningful camber whose line changes sign aft of the peak, or a
    # trailing-edge camber of opposite sign to the peak.
    te_camber = float(camber[-3])
    reflex = (abs(camber_peak) > 0.005 and
              te_camber * camber_peak < -1e-4 and
              xq[i_cam] < 0.6)

    return {
        "t_max": t_max,
        "camber_max": camber_peak,
        "abs_camber": abs(camber_peak),
        "te_camber": te_camber,
        "reflex": bool(reflex),
    }


# --------------------------------------------------------------------------- #
# Classification                                                               #
# --------------------------------------------------------------------------- #
_SIX_SERIES = re.compile(r"NACA\s*6[0-9]", re.I)
_THIN_NAMES = re.compile(r"\b(AG\d|SD\d|SG\d|MH\d|DAE\d)", re.I)
_HIGH_LIFT_NAMES = re.compile(r"\b(EPPLER|E\d{3}|^S\d{4}|SELIG|FX\b|WORTMANN|GO\b|GOE)", re.I)


def classify(name: str, g: dict) -> str:
    if _SIX_SERIES.search(name) or g["reflex"] or "SC(2)" in name.upper():
        return "reflexed_6series"
    if g["abs_camber"] < 0.006 and g["t_max"] >= 0.07:
        return "thick_symmetric"
    if (_THIN_NAMES.search(name) or g["t_max"] < 0.09) and g["abs_camber"] < 0.05:
        return "thin_low_re"
    if g["abs_camber"] >= 0.045 or _HIGH_LIFT_NAMES.search(name):
        return "high_camber"
    return "cambered_gp"


# Preferred recognisable representatives (regex) to seed each class so the
# printed subset is reviewable.
_PREFERRED = {
    "thin_low_re":      [r"^AG1\d", r"SD7037", r"SD7003", r"^AG0", r"MH\d"],
    "thick_symmetric":  [r"NACA 0012 AIR", r"NACA 0015", r"NACA 0018", r"NACA 0009", r"NACA 0008\b"],
    "cambered_gp":      [r"CLARK Y AIR", r"NACA 2412", r"NACA 4412", r"NACA 2415", r"USA-35"],
    "high_camber":      [r"EPPLER", r"^S1223", r"SELIG", r"^S\d{4}", r"FX "],
    "reflexed_6series": [r"NACA 63", r"NACA 64", r"NACA 65", r"SC\(2\)", r"E18\d"],
}


def build_catalog():
    """Parse every foil, return list of dicts {name, path, class, **geom}."""
    catalog = []
    for path in sorted(FOIL_DIR.glob("*.dat")):
        try:
            x, z = load_coords(path)
            if x.size < 20:
                continue
            g = geometry(x, z)
        except Exception:
            continue
        if not np.isfinite(g["t_max"]) or g["t_max"] <= 0.0 or g["t_max"] > 0.45:
            continue
        cls = classify(path.stem, g)
        catalog.append({"name": path.stem, "path": str(path), "class": cls, **g})
    return catalog


def select_subset(catalog, per_class=5, seed=12):
    """Stratified pick: prefer recognisable representatives, then fill randomly."""
    rng = np.random.default_rng(seed)
    by_class = {c: [] for c in CLASSES}
    for entry in catalog:
        by_class[entry["class"]].append(entry)

    chosen = []
    chosen_names = set()
    for cls in CLASSES:
        pool = by_class[cls]
        picks = []
        # 1) preferred recognisable names; cap each pattern so one foil family
        #    (e.g. the AG1x or Eppler 12xx block) cannot fill the whole class.
        for pat in _PREFERRED[cls]:
            rx = re.compile(pat, re.I)
            pat_count = 0
            for e in pool:
                if len(picks) >= per_class or pat_count >= 2:
                    break
                if e["name"] in chosen_names:
                    continue
                if rx.search(e["name"]):
                    picks.append(e)
                    chosen_names.add(e["name"])
                    pat_count += 1
        # 2) random fill from remaining
        remaining = [e for e in pool if e["name"] not in chosen_names]
        rng.shuffle(remaining)
        for e in remaining:
            if len(picks) >= per_class:
                break
            picks.append(e)
            chosen_names.add(e["name"])
        chosen.extend(picks)
    return chosen


if __name__ == "__main__":
    cat = build_catalog()
    print(f"Parsed {len(cat)} foils from {FOIL_DIR}\n")
    print("Class population:")
    for c in CLASSES:
        print(f"  {c:<18s} {sum(1 for e in cat if e['class'] == c)}")
    print("\nSelected stratified subset:")
    sub = select_subset(cat)
    for c in CLASSES:
        print(f"\n[{c}]")
        for e in sub:
            if e["class"] == c:
                print(f"  {e['name']:<42s} t={e['t_max']*100:5.1f}%  "
                      f"camber={e['camber_max']*100:+5.1f}%  reflex={e['reflex']}")
    print(f"\nTotal selected: {len(sub)}")
