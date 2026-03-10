from __future__ import annotations

from typing import TYPE_CHECKING, Optional

import numpy as np
import pandas as pd
from numpy import ndarray

if TYPE_CHECKING:
    from pynumad.objects.airfoil import Airfoil
    from pynumad.objects.station import Station
    from pynumad.objects.component import Component
    from pynumad.objects.material import Material


# Columns present in Definition.span_data (blade-stations grid).
_SPAN_DATA_COLS = [
    "chord",
    "degreestwist",
    "percentthick",
    "chordoffset",
    "aerocenter",
    "prebend",
    "sweep",
]

# Columns present in Definition.ispan_data (analysis-stations / ispan grid).
_ISPAN_DATA_COLS = [
    "leband",
    "teband",
    "sparcapwidth_hp",
    "sparcapwidth_lp",
    "sparcapoffset_hp",
    "sparcapoffset_lp",
]


class Definition:
    """Definition Class
    
    A Definition object is designed to contain the
    basic defining features of a blade, which can
    be populated manually or by reading a blade file.

    Attributes
    ----------
    components : dict[str, Component]
        Dictionary of components indexed by component name.
    shearweb : list
        List of shear webs. ``shearweb[i]`` gives an array defining the
        *i*-th shear web.
    materials : dict[str, Material]
        Dictionary of material objects indexed by material name.
    stacks : ndarray
    swstacks : ndarray
    span_data : pd.DataFrame
        Distributed aerodynamic properties at the blade-definition stations.
        Index: span [m].  Columns: ``chord``, ``degreestwist``,
        ``percentthick``, ``chordoffset``, ``aerocenter``, ``prebend``,
        ``sweep``.
    ispan_data : pd.DataFrame
        Distributed structural properties at the analysis (interpolated)
        stations.  Index: ispan [m].  Columns: ``leband``, ``teband``,
        ``sparcapwidth_hp``, ``sparcapwidth_lp``, ``sparcapoffset_hp``,
        ``sparcapoffset_lp``.
    stations : list[Station]
        List of :class:`~pynumad.objects.station.Station` objects.
    te_type : list[str]
    sparcapoffset : ndarray
        Single-array spar cap offset used by the Excel reader; the YAML
        reader populates ``ispan_data`` instead.
    sparcapwidth : ndarray
        Single-array spar cap width used by the Excel reader; the YAML
        reader populates ``ispan_data`` instead.
    """

    def __init__(self):
        self.components: dict[str, Component] = None
        self.shearweb: list = None
        self.materials: dict[str, Material] = None
        self.stations: list[Station] = None
        self.te_type: list[str] = None
        self.stacks: ndarray = None
        self.swstacks: ndarray = None

        # Blade-stations grid: aero distribution indexed by span [m]
        self.span_data: pd.DataFrame = pd.DataFrame(columns=_SPAN_DATA_COLS)

        # Analysis-stations grid: structural distribution indexed by ispan [m]
        self.ispan_data: pd.DataFrame = pd.DataFrame(columns=_ISPAN_DATA_COLS)

        # Kept for Excel-reader backward compat (single-value variants)
        self.sparcapoffset: Optional[ndarray] = None
        self.sparcapwidth: Optional[ndarray] = None

        # Configuration flags with validated setters
        self._natural_offset: int = 1
        self._rotorspin: int = 1
        self._swtwisted: int = 0

    # ------------------------------------------------------------------
    # Convenience properties: expose the DataFrame index/columns as
    # plain ndarrays so that existing downstream code that reads e.g.
    # definition.span, definition.chord continues to work without
    # modification.  Setters write back into the DataFrames.
    # ------------------------------------------------------------------

    @property
    def span(self) -> ndarray:
        """Spanwise locations of the blade-definition stations [m]."""
        return self.span_data.index.to_numpy()

    @span.setter
    def span(self, value: ndarray):
        value = np.asarray(value)
        self.span_data = pd.DataFrame(index=value, columns=_SPAN_DATA_COLS)

    @property
    def ispan(self) -> ndarray:
        """Spanwise locations of the analysis (interpolated) stations [m]."""
        return self.ispan_data.index.to_numpy()

    @ispan.setter
    def ispan(self, value: ndarray):
        value = np.asarray(value)
        self.ispan_data = pd.DataFrame(index=value, columns=_ISPAN_DATA_COLS)

    def _span_col(self, name: str) -> ndarray:
        """Return a span_data column as a float ndarray (None-safe)."""
        col = self.span_data[name]
        if col.isnull().all():
            return None
        return col.to_numpy(dtype=float, na_value=np.nan)

    def _set_span_col(self, name: str, value):
        arr = np.asarray(value, dtype=float)
        if arr.ndim == 0:
            # scalar: broadcast to current span length
            arr = np.full(len(self.span_data), float(arr))
        if self.span_data.empty or len(self.span_data) != len(arr):
            self.span_data = pd.DataFrame(
                index=pd.RangeIndex(len(arr)), columns=_SPAN_DATA_COLS
            )
        self.span_data[name] = arr

    def _ispan_col(self, name: str) -> ndarray:
        col = self.ispan_data[name]
        if col.isnull().all():
            return None
        return col.to_numpy(dtype=float, na_value=np.nan)

    def _set_ispan_col(self, name: str, value):
        arr = np.asarray(value, dtype=float)
        if arr.ndim == 0:
            # scalar: broadcast to current ispan length
            arr = np.full(len(self.ispan_data), float(arr))
        if self.ispan_data.empty or len(self.ispan_data) != len(arr):
            self.ispan_data = pd.DataFrame(
                index=pd.RangeIndex(len(arr)), columns=_ISPAN_DATA_COLS
            )
        self.ispan_data[name] = arr

    # Span-grid properties
    @property
    def chord(self) -> ndarray:
        return self._span_col("chord")

    @chord.setter
    def chord(self, v):
        self._set_span_col("chord", v)

    @property
    def degreestwist(self) -> ndarray:
        return self._span_col("degreestwist")

    @degreestwist.setter
    def degreestwist(self, v):
        self._set_span_col("degreestwist", v)

    @property
    def percentthick(self) -> ndarray:
        return self._span_col("percentthick")

    @percentthick.setter
    def percentthick(self, v):
        self._set_span_col("percentthick", v)

    @property
    def chordoffset(self) -> ndarray:
        return self._span_col("chordoffset")

    @chordoffset.setter
    def chordoffset(self, v):
        self._set_span_col("chordoffset", v)

    @property
    def aerocenter(self) -> ndarray:
        return self._span_col("aerocenter")

    @aerocenter.setter
    def aerocenter(self, v):
        self._set_span_col("aerocenter", v)

    @property
    def prebend(self) -> ndarray:
        return self._span_col("prebend")

    @prebend.setter
    def prebend(self, v):
        self._set_span_col("prebend", v)

    @property
    def sweep(self) -> ndarray:
        return self._span_col("sweep")

    @sweep.setter
    def sweep(self, v):
        self._set_span_col("sweep", v)

    # ispan-grid properties
    @property
    def leband(self) -> ndarray:
        return self._ispan_col("leband")

    @leband.setter
    def leband(self, v):
        self._set_ispan_col("leband", v)

    @property
    def teband(self) -> ndarray:
        return self._ispan_col("teband")

    @teband.setter
    def teband(self, v):
        self._set_ispan_col("teband", v)

    @property
    def sparcapwidth_hp(self) -> ndarray:
        return self._ispan_col("sparcapwidth_hp")

    @sparcapwidth_hp.setter
    def sparcapwidth_hp(self, v):
        self._set_ispan_col("sparcapwidth_hp", v)

    @property
    def sparcapwidth_lp(self) -> ndarray:
        return self._ispan_col("sparcapwidth_lp")

    @sparcapwidth_lp.setter
    def sparcapwidth_lp(self, v):
        self._set_ispan_col("sparcapwidth_lp", v)

    @property
    def sparcapoffset_hp(self) -> ndarray:
        return self._ispan_col("sparcapoffset_hp")

    @sparcapoffset_hp.setter
    def sparcapoffset_hp(self, v):
        self._set_ispan_col("sparcapoffset_hp", v)

    @property
    def sparcapoffset_lp(self) -> ndarray:
        return self._ispan_col("sparcapoffset_lp")

    @sparcapoffset_lp.setter
    def sparcapoffset_lp(self, v):
        self._set_ispan_col("sparcapoffset_lp", v)

    # ------------------------------------------------------------------

    def __eq__(self, other):
        if not isinstance(other, Definition):
            return NotImplemented
        attrs = vars(self).keys()
        for attr in attrs:
            self_attr = getattr(self, attr)
            other_attr = getattr(other, attr)
            if isinstance(self_attr, pd.DataFrame):
                if not self_attr.equals(other_attr):
                    return False
            elif isinstance(self_attr, (int, float, str, list, dict)):
                if self_attr != other_attr:
                    return False
            elif isinstance(self_attr, ndarray):
                if (self_attr != other_attr).any():
                    return False
        return True

    @property
    def natural_offset(self):
        """
        1 = offset by max thickness location,
        0 = do not offset to max thickness
        """
        return self._natural_offset

    @natural_offset.setter
    def natural_offset(self, new_natural_offset):
        if not (new_natural_offset == 0 or new_natural_offset == 1):
            raise Exception("natural_offset must be 0 or 1")
        else:
            self._natural_offset = new_natural_offset

    @property
    def rotorspin(self):
        """
        Rotor Spin,
        1 = CW rotation looking downwind,
        -1 = CCW rotation
        """
        return self._rotorspin

    @rotorspin.setter
    def rotorspin(self, new_rotorspin):
        if not (new_rotorspin == 1 or new_rotorspin == -1):
            raise Exception("rotorspin must be 1 (cw) or -1 (ccw)")
        else:
            self._rotorspin = new_rotorspin

    @property
    def swtwisted(self):
        """
        Shear Web,
        0 = planar shear webs,
        1 = shear webs twisted by blade twist
        """
        return self._swtwisted

    @swtwisted.setter
    def swtwisted(self, new_swtwisted):
        if not (new_swtwisted == 0 or new_swtwisted == 1):
            raise Exception("swtwisted must be 0 or 1")
        else:
            self._swtwisted = new_swtwisted

    def add_station(self, af: Airfoil, spanlocation: float):
        """Append a new station to this definition.

        Parameters
        ----------
        af : Airfoil or str
            Airfoil object or path to airfoil coordinate file.
        spanlocation : float
            Spanwise location of the station [m].
        """
        from pynumad.objects.station import Station  # deferred to avoid circular import
        new_station = Station(af)
        new_station.spanlocation = spanlocation
        self.stations.append(new_station)
        return self
