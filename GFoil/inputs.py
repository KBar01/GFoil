import numpy as np
from dataclasses import dataclass, field, fields, is_dataclass
from typing import Optional, Sequence, Union


# --------------------------------------------------------------------------- #
# Result-object presentation / dict-access mixin                              #
# --------------------------------------------------------------------------- #



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

    """
    This is first key input, an Aerofoil dataclass that should be all info about aerofoil
    geometry, that includes the panelling coeffs as well 
    """
    xcoords: np.ndarray
    ycoords: np.ndarray
    chord: Optional[float] = 1.0
    span: Optional[float] = 3.0
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
    """
    This is acoustics input, so obervers, sampling location for WPS,
    what WPS model to use, and frequency range + A-weighting option 
    """
    observerXYZ: np.ndarray
    # TESampleLoc: scalar x/c float OR length-2 [x_lo, x_hi] window.
    #   scalar x   -> 0 < x < 1; stored as bare float; gfoil.py equal-packs to C++ single-point path.
    #   [x_lo,x_hi]-> 0 < x_lo < x_hi < 1; stored as (x_lo,x_hi) tuple; C++ BL-averaged window.
    #   x/c == 1.0 is rejected (TE node is degenerate for the interpolation stencil).
    TESampleLoc: Optional[Union[float, Sequence[float]]] = 0.98
    model: Optional[str] = "kam"
    aWeighting: bool = False
    f_min: float = 200.0   # lower acoustic frequency bound [Hz]
    f_max: float = 20000.0 # upper acoustic frequency bound [Hz]

    def __post_init__(self): # make sure its passed properly
        arr = _as_float_array(self.observerXYZ, "Acoustics.observerXYZ").astype(float)
        if arr.ndim == 1 and arr.size == 3:
            arr = arr.reshape(1, 3)
        elif arr.ndim == 2 and arr.shape[1] == 3:
            pass  # already (N, 3)
        else:
            raise ValueError(f"observerXYZ must be shape (3,) or (N,3), got {arr.shape}")
        self.observerXYZ = arr  # always (N, 3)

        if self.TESampleLoc is not None:
            te = np.asarray(self.TESampleLoc, dtype=float).ravel()
            if te.size == 1:
                x = float(te[0])
                if not (0.0 < x < 1.0):
                    raise ValueError(
                        f"TESampleLoc scalar must satisfy 0 < x < 1 (got {x}); "
                        "x/c == 1.0 is rejected (TE node is degenerate; "
                        "use a window for TE-region sampling)"
                    )
                self.TESampleLoc = x
            elif te.size == 2:
                x_lo, x_hi = float(te[0]), float(te[1])
                if not (0.0 < x_lo < x_hi < 1.0):
                    raise ValueError(
                        f"TESampleLoc window must satisfy 0 < x_lo < x_hi < 1 "
                        f"(got [{x_lo}, {x_hi}]); x/c == 1.0 is rejected (TE node is degenerate)"
                    )
                self.TESampleLoc = (x_lo, x_hi)
            else:
                raise ValueError(
                    f"TESampleLoc must be a scalar x/c or a length-2 [x_lo, x_hi] window; "
                    f"got {te.size} elements"
                )
           

        if self.f_min <= 0 or self.f_max <= self.f_min:
            raise ValueError("f_min must be > 0 and f_max > f_min")


@dataclass
class OperatingConds:

    """
    the conditions aerofoil is operating in, such as Re and angle of attack
    and other key params
    """

    alpha: Optional[float] = 0.0
    Re: Optional[float] = 2e6
    rho: Optional[float] = 1.225
    Ma: Optional[float] = 0.0
    nu: Optional[float] = 0.000015
    nCrit:     Optional[float] = 9.0
    transition: np.ndarray = field(default_factory=lambda: np.array([1.0, 1.0], dtype=float))
    rtol: Optional[float] = 1e-9   # solver RMS residual tolerance to be considered converged 

    def __post_init__(self):
        self.transition = _as_float_array(self.transition, "OperatingConds.transition").astype(float)

        if self.transition.size != 2:
            raise ValueError(f"transition must have length 2, got shape {self.transition.shape}")
        self.transition = self.transition.reshape(2,)

        self.rtol = float(self.rtol)
        if not (0.0 < self.rtol < 1.0):
            raise ValueError(f"rtol must be a positive float in (0, 1), got {self.rtol}")


##############################################################################################################

"""
Next below is output structures, where (through the use of some Claude), I have made them output how i want
them to, mainly to do with there structure in printing them etc
"""

def _short_descr(v):
    """One-line descriptor of a field value for the compact repr.

    Arrays/lists are summarised by shape/length (never dumped); scalars show
    their value; nested result dataclasses show "<ClassName> (set)".
    """
    if v is None:
        return "None"
    if isinstance(v, np.ndarray):
        return f"ndarray shape={tuple(v.shape)} dtype={v.dtype}"
    if isinstance(v, bool):               
        return str(v)
    if isinstance(v, list):
        return f"list len={len(v)}"
    if isinstance(v, tuple):
        return f"tuple len={len(v)}"
    if is_dataclass(v) and not isinstance(v, type):
        return f"{type(v).__name__} (set)"
    if isinstance(v, (float, np.floating)):
        return f"{float(v):.5g}"
    if isinstance(v, str):
        return f'"{v}"'
    return str(v)


class _ResultMixin:
    """Read-only dict-style access + compact repr for result dataclasses.

    Mixed into dataclasses declared with ``@dataclass(repr=False)`` so this
    ``__repr__`` is used instead of the verbose auto-generated one. Adds
    ``result["CL"]`` access alongside attribute access ``result.CL`` (read
    only — no ``__setitem__``), plus ``keys()``, ``in``, and ``to_dict()``.

    Field names and types are untouched, so direct attribute access (which the
    grad_run AD-upload path relies on) is unaffected.
    """

    # Subclasses may set _REPR_GROUPS = [(label_or_None, [field_name, ...]), ...]
    # to group the repr. Any field not listed is appended under an "other" group
    # so the repr stays complete if fields are added later.
    _REPR_GROUPS = None

    def _field_names(self):
        return [f.name for f in fields(self)]

    def _repr_layout(self):
        names = self._field_names()
        if self._REPR_GROUPS is None:
            return [(None, names)]
        listed = [n for _, grp in self._REPR_GROUPS for n in grp]
        missing = [n for n in names if n not in listed]
        layout = [(lbl, [n for n in grp if n in names])
                  for lbl, grp in self._REPR_GROUPS]
        if missing:
            layout.append(("other", missing))
        return layout

    def summary(self) -> str:
        """Return the compact multi-line summary string (same as repr).

        Handy for logging: ``logger.info(result.summary())``.
        """
        names = self._field_names()
        width = max((len(n) for n in names), default=0)
        lines = [type(self).__name__]
        for label, group in self._repr_layout():
            if label:
                lines.append(f"  -- {label} --")
            for n in group:
                lines.append(f"  {n:<{width}}  {_short_descr(getattr(self, n))}")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return self.summary()

    __str__ = __repr__

    # ---- read-only dict-style access --------------------------------------- #
    def __getitem__(self, key):
        if key in self._field_names():
            return getattr(self, key)
        raise KeyError(
            f"{key!r} is not a field of {type(self).__name__}. "
            f"Valid keys: {self._field_names()}"
        )

    def keys(self):
        """Field names, in declaration order (enables ``dict(result)``)."""
        return self._field_names()

    def __contains__(self, key) -> bool:
        return key in self._field_names()

    def to_dict(self) -> dict:
        """Plain ``{field_name: value}`` mapping.

        Shallow for arrays/lists (returned by reference, not copied). Nested
        result dataclasses (e.g. ``verbose_data``) are recursed into, so they
        become nested dicts rather than dataclass instances.
        """
        out = {}
        for f in fields(self):
            v = getattr(self, f.name)
            if is_dataclass(v) and not isinstance(v, type) and hasattr(v, "to_dict"):
                v = v.to_dict()
            out[f.name] = v
        return out


@dataclass(repr=False)
class VerboseResult(_ResultMixin):
    """
    This is verbose output showing entire flow solution of the solver,
    being Cp distribution, IBL values, transition location, as well as some
    of the acoustics stuff like WPS and FF-spectra
    """
    # Geometry
    x: np.ndarray        # panel node x-coords     shape (Ncoords,)
    y: np.ndarray        # panel node y-coords     shape (Ncoords,)

    # Per-node aero
    Cp: np.ndarray          # pressure coefficient          shape (Ncoords,)
    delta_star: np.ndarray  # displacement thickness [m]    shape (Ncoords,)
    theta: np.ndarray       # momentum thickness     [m]    shape (Ncoords,)
    tau_wall: np.ndarray    # wall shear stress      [Pa]   shape (Ncoords,)
    tau_max: np.ndarray     # max shear stress       [Pa]   shape (Ncoords,); 0 if laminar
    Ue: np.ndarray          # BL edge velocity       [m/s]  shape (Ncoords,)
    dpdx: np.ndarray        # pressure gradient      [Pa/m] shape (Ncoords,)
    is_turb: np.ndarray     # bool turbulence flag          shape (Ncoords,)

    # Transition
    topTransX: float     # upper surface transition x
    botTransX: float     # lower surface transition x

    # TE sampling: BL inputs to WPS / Amiet
    # Order: [theta, delta*, tau_max, Ue, dpdx, tau_wall, delta99]
    BL_top: np.ndarray   # shape (7,)  upper surface TE BL properties
    BL_bot: np.ndarray   # shape (7,)  lower surface TE BL properties

    # Acoustic spectra
    freq_Hz: np.ndarray    # frequency array [Hz]               shape (Nsound,)
    WPS_upper: np.ndarray  # upper surface WPS [dB/Hz re 20µPa] shape (Nsound,)
    WPS_lower: np.ndarray  # lower surface WPS [dB/Hz re 20µPa] shape (Nsound,)
    FF_spectra: np.ndarray # far-field spectra     [dB/Hz re 20µPa] shape (nObs, Nsound)

    # Per-observer integrated noise and observer geometry
    OASPL_perObs:   np.ndarray  # OASPL per observer [dB re 20µPa]            shape (nObs,)
    obsXYZ_TElocal: np.ndarray  # observer coords in TE-local Amiet frame [m] shape (nObs, 3) (chord-aligned, origin at the trailing edge)


@dataclass(repr=False)
class NoiseResult(_ResultMixin):
    """Returned by noise_run (acoustics-only, no aero solve).

    All spectra are RAW LINEAR wall-pressure / far-field PSD in Pa^2/omega
    (i.e. per rad/s — NOT per Hz, NOT in dB, NOT integrated to OASPL). A
    surface whose BL state has tau_max <= 0 (or an all-zero custom-WPS column)
    is skipped and its WPS column is zeros.
    """
    freqs_Hz:       np.ndarray  # frequency array [Hz]                 shape (N,)
    WPS_upper:      np.ndarray  # upper surface WPS [Pa^2/omega]       shape (N,)
    WPS_lower:      np.ndarray  # lower surface WPS [Pa^2/omega]       shape (N,)
    FF_spectra:     np.ndarray  # far-field PSD     [Pa^2/omega]       shape (nObs, N)
    obsXYZ_TElocal: np.ndarray  # observer coords in TE-local frame [m] shape (nObs, 3)


@dataclass(repr=False)
class FwdResult(_ResultMixin):
    """Returned by fwd_run. Pass to grad_run to get gradients."""

    _REPR_GROUPS = [
        (None, ["converged", "CL", "CD", "CM", "OASPL", "failure_mode", "newton_iterations"]),
        ("jacobian state", ["states", "turb", "stag", "RVvals", "RVrows", "RVcols", "RVnz"]),
        ("design point", ["ycoords", "alpha"]),
        (None, ["verbose_data"]),
    ]

    converged: bool
    CL:    float = 0.0
    CD:    float = 0.0
    CM:    float = 0.0
    OASPL: float = 0.0
    # Jacobian state for the AD pass
    states:  Optional[list] = None
    turb:    Optional[list] = None
    stag:    Optional[list] = None
    RVvals:  Optional[list] = None
    RVrows:  Optional[list] = None
    RVcols:  Optional[list] = None
    RVnz:    int = 0
    # Design point (for Jacobian-validity checks)
    ycoords: Optional[np.ndarray] = None
    alpha:   float = 0.0
    # Verbose output (None unless fwd_run called with verbose=True)
    verbose_data: Optional["VerboseResult"] = None
    newton_iterations: int = 0
    # Empty string when converged.  One of:
    #   "transition_front_oscillation" which is residual oscillating
    #   "diverged" which is residualNorm >= 1.0 at exit
    #   "no_convergence" is just anything else really, not very helpful
    failure_mode: str = ""


@dataclass(repr=False)
class GradResult(_ResultMixin):
    """Returned by grad_run."""
    converged: bool
    dCL_dy:        Optional[np.ndarray] = None
    dCD_dy:        Optional[np.ndarray] = None
    dOASPL_dy:     Optional[np.ndarray] = None
    dCL_dalpha:    float = 0.0
    dCD_dalpha:    float = 0.0
    dOASPL_dalpha: float = 0.0