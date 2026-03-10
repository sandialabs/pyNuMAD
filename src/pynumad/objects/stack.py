"""Composite stack and ply data classes.

These classes were previously defined in ``stackdb.py``.  They are
extracted here so that they can be imported independently without
pulling in the heavy :class:`~pynumad.objects.stackdb.StackDatabase`
logic.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from numpy import ndarray


@dataclass
class Ply:
    """A single ply (or consolidated ply group) within a :class:`Stack`.

    Attributes
    ----------
    component : str
        Name of the parent blade component.
    materialid : str
        Material identifier key into the blade material dictionary.
    thickness : float
        Nominal thickness of a single ply [mm].
    angle : float
        Fibre orientation angle [degrees].
    nPlies : int
        Number of plies in this group.
    """

    component: Optional[str] = None
    materialid: Optional[str] = None
    thickness: Optional[float] = None
    angle: Optional[float] = None
    nPlies: Optional[int] = None


class Stack:
    """A stack of composite plies for a single blade surface region.

    Parameters
    ----------
    name : str
        Human-readable label, e.g. ``'00_01_HP_SPAR'``.
    indices : list[int]
        ``[inboard station, outboard station, keypoint start, keypoint end]``
    plygroups : list[Ply]
        Ordered list of :class:`Ply` groups, inboard to outboard.

    Notes
    -----
    Consecutive plies from the *same* component at the *same* angle are
    automatically consolidated by :meth:`addply` (incrementing
    ``nPlies``).
    """

    def __init__(self):
        self.name: str = ""
        self.indices: list = []
        self.plygroups: list = []

    def __eq__(self, other):
        if not isinstance(other, Stack):
            return NotImplemented
        attrs = vars(self).keys()
        for attr in attrs:
            self_attr = getattr(self, attr)
            other_attr = getattr(other, attr)
            if isinstance(self_attr, (int, float, str, list, dict)):
                if self_attr != other_attr:
                    return False
            elif isinstance(self_attr, ndarray):
                if (self_attr != other_attr).any():
                    return False
        return True

    def addply(self, ply: Ply) -> "Stack":
        """Append *ply* to this stack, consolidating with the last group if compatible.

        Two plies are consolidated when they share the same component name
        and fabric angle; in that case ``nPlies`` on the last group is
        incremented rather than appending a new :class:`Ply`.

        Parameters
        ----------
        ply : Ply

        Returns
        -------
        self
        """
        plygroups = self.plygroups
        if (
            plygroups
            and ply.component == plygroups[-1].component
            and ply.angle == plygroups[-1].angle
        ):
            plygroups[-1].nPlies += 1
        else:
            plygroups.append(ply)
        return self

    def layer_thicknesses(self) -> ndarray:
        """Return total thickness contributed by each ply group.

        Returns
        -------
        ndarray
            1-D array of ``nPlies * thickness`` values, one per ply group.
        """
        return np.array(
            [pg.nPlies * pg.thickness for pg in self.plygroups]
        )
