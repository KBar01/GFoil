
import numpy as np
from dataclasses import dataclass, field
from typing import Optional



def _as_1d_float_array(a, name: str) -> np.ndarray:
    arr = np.asarray(a, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be 1D, got shape {arr.shape}")
    return arr

def _as_float_array(a, name: str) -> np.ndarray:
    try:
        return np.asarray(a, dtype=float)
    except Exception as e:
        raise ValueError(f"{name} could not be converted to float array") from e

def _require_length(arr: np.ndarray, n: int, name: str) -> np.ndarray:
    if arr.size != n:
        raise ValueError(f"{name} must have length {n}, got {arr.size}")
    return arr


@dataclass
class Aerofoil:
    xcoords: np.ndarray
    ycoords: np.ndarray
    chord: Optional[float] = 1.0
    span: Optional[float] = 2.0
    panelUniformity: Optional[float] = 1.0
    panelTEspacing: Optional[float] = 0.09

    def __post_init__(self):
        self.xcoords = _as_1d_float_array(self.xcoords, "Aerofoil.xcoords")
        self.ycoords = _as_1d_float_array(self.ycoords, "Aerofoil.ycoords")

        if self.xcoords.size != self.ycoords.size:
            raise ValueError(
                f"xcoords and ycoords must be the same length; "
                f"got {self.xcoords.size} and {self.ycoords.size}"
            )
        if self.xcoords.size < 2:
            raise ValueError("Need at least 2 points")

        x = self.xcoords
        y = self.ycoords

        te_x = 1.0
        tol = 1e-8

        te_mask = np.isclose(x, te_x, atol=tol, rtol=0.0)
        te_idx = np.flatnonzero(te_mask)

        if te_idx.size < 2:
            raise ValueError(
                "Could not find two TE points where x == 1.0. "
                f"Found {te_idx.size}. Ensure both lower and upper TE points have x=1.0."
            )

        first_te = int(te_idx.min())
        last_te = int(te_idx.max())

        if y[first_te] > y[last_te]:
            self.xcoords = self.xcoords[::-1].copy()
            self.ycoords = self.ycoords[::-1].copy()


@dataclass
class Acoustics:
    observerXYZ: np.ndarray
    TESampleLoc: Optional[float] = 0.97
    model: Optional[str] = "roz"

    def __post_init__(self):
        self.observerXYZ = _as_float_array(self.observerXYZ, "Acoustics.observerXYZ").astype(float)

        if self.observerXYZ.size != 3:
            raise ValueError(f"observerXYZ must have length 3, got shape {self.observerXYZ.shape}")
        self.observerXYZ = self.observerXYZ.reshape(3,)

        if self.TESampleLoc is not None:
            if not (0.0 <= self.TESampleLoc <= 1.0):
                raise ValueError("TESampleLoc must be between 0 and 1 (inclusive)")


@dataclass
class OperatingConds:
    alpha: Optional[float] = 0.0
    Re: Optional[float] = 2e6
    rho: Optional[float] = 1.225
    Ma: Optional[float] = 0.0
    nu: Optional[float] = 0.000015
    nCrit: Optional[float] = 9.0
    transition: np.ndarray = field(default_factory=lambda: np.array([1.0, 1.0], dtype=float))

    def __post_init__(self):
        self.transition = _as_float_array(self.transition, "OperatingConds.transition").astype(float)

        if self.transition.size != 2:
            raise ValueError(f"transition must have length 2, got shape {self.transition.shape}")
        self.transition = self.transition.reshape(2,)


@dataclass
class WPSinfo:
    Re: float
    observerXYZ: np.ndarray

    DispThick: np.ndarray
    MomThick: np.ndarray
    BLHeight: np.ndarray
    wallShear: np.ndarray
    maxShear: np.ndarray
    edgeVel: np.ndarray
    dpdx: np.ndarray

    chord: Optional[float] = 1.0
    span: Optional[float] = 2.0
    model: Optional[str] = "roz"
    rho: Optional[float] = 1.225
    nu: Optional[float] = 1.5e-5

    def __post_init__(self):
        self.observerXYZ = _as_1d_float_array(self.observerXYZ, "WPSinfo.observerXYZ")
        _require_length(self.observerXYZ, 3, "WPSinfo.observerXYZ")
        self.observerXYZ = self.observerXYZ.reshape(3,)

        for name in (
            "DispThick",
            "MomThick",
            "BLHeight",
            "wallShear",
            "maxShear",
            "edgeVel",
            "dpdx",
        ):
            arr = _as_1d_float_array(getattr(self, name), f"WPSinfo.{name}")
            arr = _require_length(arr, 2, f"WPSinfo.{name}")
            setattr(self, name, arr)


@dataclass
class customSpectrainfo:
    
    Freq: np.ndarray
    WPSUpper: np.ndarray
    WPSLower: np.ndarray
    Re: float
    observerXYZ: np.ndarray
    chord: Optional[float] = 1.0
    span: Optional[float] = 2.0
    rho: Optional[float] = 1.225
    nu: Optional[float] = 1.5e-5

    def __post_init__(self):
        self.observerXYZ = _as_1d_float_array(self.observerXYZ, "WPSinfo.observerXYZ")
        _require_length(self.observerXYZ, 3, "WPSinfo.observerXYZ")
        self.observerXYZ = self.observerXYZ.reshape(3,)

