from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from pynumad.objects.airfoil import Airfoil


@dataclass(init=False, eq=False)
class Station:
    """A single spanwise blade station.

    Attributes
    ----------
    airfoil : Airfoil
        Airfoil profile at this station.
    spanlocation : float
        Spanwise location where the station is defined [m].
    """

    airfoil: Airfoil = field(default_factory=Airfoil)
    spanlocation: Optional[float] = None

    def __init__(self, af=None):
        """Create a Station, optionally initialising its airfoil.

        Parameters
        ----------
        af : str or Airfoil, optional
            If a ``str``, treat it as a path to an airfoil XML file.
            If an :class:`Airfoil`, use it directly.
            If ``None`` (default), create an empty :class:`Airfoil`.
        """
        self.spanlocation = None

        if isinstance(af, str):
            self.airfoil = Airfoil(filename=af)
        elif isinstance(af, Airfoil):
            self.airfoil = af
        else:
            self.airfoil = Airfoil()

    def __eq__(self, other):
        if not isinstance(other, Station):
            return NotImplemented
        return (
            self.airfoil == other.airfoil
            and self.spanlocation == other.spanlocation
        )
