
import numpy as np
from dataclasses import dataclass, field, fields, is_dataclass
from typing import Optional


# --------------------------------------------------------------------------- #
# Result-object presentation / dict-access mixin                              #
# --------------------------------------------------------------------------- #
def _short_descr(v) -> str:
    """One-line descriptor of a field value for the compact repr.

    Arrays/lists are summarised by shape/length (never dumped); scalars show
    their value; nested result dataclasses show "<ClassName> (set)".
    """
    if v is None:
        return "None"
    if isinstance(v, np.ndarray):
        return f"ndarray shape={tuple(v.shape)} dtype={v.dtype}"
    if isinstance(v, bool):                 # before int (bool is a subclass)
        return str(v)
    if isinstance(v, list):
        return f"list len={len(v)}"
    if isinstance(v, tuple):
        return f"tuple len={len(v)}"
    if is_dataclass(v) and not isinstance(v, type):
        return f"{type(v).__name__} (set)"
    if isinstance(v, (float, np.floating)):
        return f"{float(v):.5g}"            # repr only — stored value untouched
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
    """
    observerXYZ: Observer position(s) in the global freestream-aligned frame.
    Origin at the quarter-chord point. x=downstream, z=up, y=spanwise.
    Shape (N,3) or (3,) for a single observer. The code converts these to the
    TE-local chord-aligned Amiet frame accounting for angle of attack.
    """
    observerXYZ: np.ndarray
    TESampleLoc: Optional[float] = 0.97
    model: Optional[str] = "roz"
    aWeighting: bool = False
    f_min: float = 200.0   # lower acoustic frequency bound [Hz]
    f_max: float = 20000.0 # upper acoustic frequency bound [Hz]

    def __post_init__(self):
        arr = _as_float_array(self.observerXYZ, "Acoustics.observerXYZ").astype(float)
        if arr.ndim == 1 and arr.size == 3:
            arr = arr.reshape(1, 3)
        elif arr.ndim == 2 and arr.shape[1] == 3:
            pass  # already (N, 3)
        else:
            raise ValueError(f"observerXYZ must be shape (3,) or (N,3), got {arr.shape}")
        self.observerXYZ = arr  # always (N, 3)

        if self.TESampleLoc is not None:
            if not (0.0 <= self.TESampleLoc <= 1.0):
                raise ValueError("TESampleLoc must be between 0 and 1 (inclusive)")

        if self.f_min <= 0 or self.f_max <= self.f_min:
            raise ValueError("f_min must be > 0 and f_max > f_min")


@dataclass
class OperatingConds:
    alpha: Optional[float] = 0.0
    Re: Optional[float] = 2e6
    rho: Optional[float] = 1.225
    Ma: Optional[float] = 0.0
    nu: Optional[float] = 0.000015
    nCrit:     Optional[float] = 9.0
    transition: np.ndarray = field(default_factory=lambda: np.array([1.0, 1.0], dtype=float))

    def __post_init__(self):
        self.transition = _as_float_array(self.transition, "OperatingConds.transition").astype(float)

        if self.transition.size != 2:
            raise ValueError(f"transition must have length 2, got shape {self.transition.shape}")
        self.transition = self.transition.reshape(2,)


@dataclass(repr=False)
class VerboseResult(_ResultMixin):
    """
    Rich per-node and acoustic data returned when fwd_run(verbose=True).

    All per-node arrays are length Ncoords (200) and follow the internal
    panel ordering: lower surface TE -> LE, then upper surface LE -> TE.
    Physical units assume the chord length supplied in Aerofoil.chord.
    """
    # Geometry
    x: np.ndarray        # panel node x-coords    shape (Ncoords,)
    y: np.ndarray        # panel node y-coords     shape (Ncoords,)

    # Per-node aero
    Cp: np.ndarray          # pressure coefficient       shape (Ncoords,)
    delta_star: np.ndarray  # displacement thickness [m] shape (Ncoords,)
    theta: np.ndarray       # momentum thickness [m]     shape (Ncoords,)
    tau_wall: np.ndarray    # wall shear stress [Pa]     shape (Ncoords,)
    tau_max: np.ndarray     # max shear stress [Pa]      shape (Ncoords,); 0 if laminar
    Ue: np.ndarray          # BL edge velocity [m/s]     shape (Ncoords,)
    dpdx: np.ndarray        # pressure gradient [Pa/m]   shape (Ncoords,)
    is_turb: np.ndarray     # bool turbulence flag       shape (Ncoords,)

    # Transition
    topTransX: float     # upper surface transition x
    botTransX: float     # lower surface transition x

    # TE sampling: BL inputs to WPS / Amiet
    # Order: [theta, delta*, tau_max, Ue, dpdx, tau_wall, delta99]
    BL_top: np.ndarray   # shape (7,)  upper surface TE BL properties
    BL_bot: np.ndarray   # shape (7,)  lower surface TE BL properties

    # Acoustic spectra
    freq_Hz: np.ndarray    # frequency array [Hz]              shape (Nsound,)
    WPS_upper: np.ndarray  # upper surface WPS [dB/Hz re 20µPa] shape (Nsound,)
    WPS_lower: np.ndarray  # lower surface WPS [dB/Hz re 20µPa] shape (Nsound,)
    FF_spectra: np.ndarray # far-field PSD     [dB/Hz re 20µPa] shape (nObs, Nsound)

    # Per-observer integrated noise and observer geometry
    OASPL_perObs:   np.ndarray  # OASPL per observer [dB re 20µPa]            shape (nObs,)
    obsXYZ_TElocal: np.ndarray  # observer coords in TE-local Amiet frame [m] shape (nObs, 3)
                                # (chord-aligned, origin at the trailing edge)


@dataclass(repr=False)
class FwdResult(_ResultMixin):
    """Returned by fwd_run. Pass to grad_run to get gradients."""

    # repr grouping: headline scalars, then the bulky Jacobian-state arrays
    # (shape-only), the design point, and the optional verbose payload.
    _REPR_GROUPS = [
        (None, ["converged", "CL", "CD", "CM", "OASPL",
                "failure_mode", "newton_iterations"]),
        ("jacobian state", ["states", "turb", "stag",
                            "RVvals", "RVrows", "RVcols", "RVnz"]),
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
    #   "transition_front_oscillation" — ilam stable 15+ iters, residual oscillating
    #   "diverged"                     — residualNorm >= 1.0 at exit
    #   "no_convergence"               — catch-all
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


