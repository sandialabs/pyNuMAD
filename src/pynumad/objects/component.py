########################################################################
#                    Part of the SNL NuMAD Toolbox                     #
#  Developed by Sandia National Laboratories Wind Energy Technologies  #
#              See license.txt for disclaimer information              #
########################################################################

from matplotlib import pyplot as plt
import numpy as np
import matplotlib.pyplot as plt

from pynumad.utils.interpolation import interpolator_wrap


class Component:
    """ComponentDef:  A class definition for blade components.

    Parameters
    ----------
    None

    Attributes
    ----------
    group : int
        0 = blade, 1 = first shear web, 2 = second shear web, etc.
    name : str
        Name, such as 'spar'
    materialid : str
        Material id number from blade.materials
    fabricangle : float
        Fiber angle
    hpextents : list
        Array of keypoints such as {'b','c'}
    lpextents : list
        String Array: Array of keypoints such as {'b','c'}
    cp : np
        control points defining layer distribution
    imethod: str
        imethod
    pinnedends
    hCtrl
    hLine

    Examples:

        ``comp_self = ComponentDef();``

        ``comp_self = ComponentDef(comp_struct);``

    """

    def __init__(self):
        self.group: int = None
        self.name: str = None
        self.materialid: str = None
        self.fabricangle: float = None
        self.hpextents: list = None
        self.lpextents: list = None
        self.cp: np.ndarray = None
        self.imethod: str = "linear"
        self.pinnedends: bool = None
        self.hCtrl = None
        self.hLine = None
        self._keylabels = [
            "te",
            "e",
            "d",
            "c",
            "b",
            "a",
            "le",
            "a",
            "b",
            "c",
            "d",
            "e",
            "te",
        ]

    def getcp(self):
        if self.pinnedends:
            if np.any(self.cp[:, 0] < 0) or np.any(self.cp[:, 0] > 1):
                raise Exception(
                    'ComponentDef: first coordinate of control points must be in range [0,1] when using "pinned" ends'
                )
            cpx = np.concatenate(([-0.01], self.cp[:, 0], [1.01]))
            cpy = np.concatenate(([0], self.cp[:, 1], [0]))
        else:
            cpx = self.cp[:, 0]
            cpy = self.cp[:, 1]

        return cpx, cpy

    def get_num_layers(self, span):
        cpx, cpy = self.getcp()
        return interpolator_wrap(cpx, cpy, span, self.imethod, 0)
    
    def find_region_extents(self):
        """
        TODO docstring
        """
        le = self._keylabels.index("le")
        # "_keylabels" is expected to wrap from te on hp side around to te on lp side
        try:
            if len(self.hpextents) == 2:
                try:
                    hp1 = self._keylabels[:le + 1].index(self.hpextents[0])
                except KeyError:
                    print(f'HP extent label "{self.hpextents[0]}" not defined.')
                try:
                    hp2 = self._keylabels[:le + 1].index(self.hpextents[1])
                except KeyError:
                    print(f'HP extent label "{self.hpextents[1]}" not defined.')
                hpRegion = [hp1, hp2]
                hpRegion.sort()
            else:
                hpRegion = []
        except TypeError:
            hpRegion = []

        try:
            if len(self.lpextents) == 2:
                try:
                    lp1 = self._keylabels[le:].index(self.lpextents[0]) + le
                except KeyError:
                    print(f'HP extent label "{self.hpextents[0]}" not defined.')
                try:
                    lp2 = self._keylabels[le:].index(self.lpextents[1]) + le
                except KeyError:
                    print(f'HP extent label "{self.hpextents[1]}" not defined.')
                lpRegion = [lp1, lp2]
                lpRegion.sort()
            else:
                lpRegion = []
        except TypeError:
            lpRegion = []

        # if length(comp['hpextents'])==1 && length(comp['lpextents'])==1
        # sw1 = find(1==strcmpi(comp['hpextents']{1},_keylabels(1:le)));
        # assert(~isempty(sw1),'HP extent label "#s" not defined.',comp['hpextents']{1});
        #  w2 = find(1==strcmpi(comp['lpextents']{1},_keylabels(le:end))) + le-1;
        # assert(~isempty(sw2),'LP extent label "#s" not defined.',comp['lpextents']{1});
        # swRegion = [sw1 sw2];
        # else
        swRegion = []
        return hpRegion, lpRegion  # ,swRegion
