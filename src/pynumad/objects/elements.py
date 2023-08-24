import numpy as np

# typing
from numpy import ndarray

import numpy as np
import os
from pynumad.objects.airfoil import Airfoil
from pynumad.utils.interpolation import interpolator_wrap


class MatDBentry:
    """A simple class to organize the attributes of a material"""

    def __init__(self):
        self.type: str = None
        self.name: str = None
        self.reference: str = None
        self.dens: list = None
        self.nuxy: list = None
        self.ex: list = None
        self.ey: list = None
        self.ez: list = None
        self.gxy: list = None
        self.gyz: list = None
        self.gxz: list = None
        self.prxy: list = None
        self.pryz: list = None
        self.prxz: list = None
        self.xten: list = None
        self.xcmp: list = None
        self.yten: list = None
        self.ycmp: list = None
        self.zten: list = None
        self.zcmp: list = None
        self.xy: list = None
        self.yz: list = None
        self.xz: list = None
        self.xycp: list = None
        self.yzcp: list = None
        self.xzcp: list = None
        self.xzit: list = None
        self.xzic: list = None
        self.yzit: list = None
        self.yzic: list = None
        self.g1g2: list = None
        self.etal: list = None
        self.etat: list = None
        self.alp0: list = None
        self.thicknessType: list = None
        self.uniqueLayers: list = None
        self.symmetryType: list = None
        self.layer: list = None


class Layer:
    """A simple class to organize the attributes of a material layer.

    Attributes
    ----------

    self.layerName : str
    self.thicknessA : float
    self.thicknessB :float
    self.quantity : int
    self.theta : float
    """

    def __init__(self):
        self.layerName: str = None
        self.thicknessA: float = None
        self.thicknessB: float = None
        self.quantity: int = None
        self.theta: float = None


class Shearweb:
    """A simple class to organize the attributes of a Shearweb

    Attributes
    ----------

    Material : str
    BeginStation : int
    EndStation : int
    Corner : list
    """

    def __init__(self):
        self.Material: str = None
        self.BeginStation: int = None
        self.EndStation: int = None
        self.Corner: list = None


class BOM:
    """A simple class to organize the attributes of a Bill of Materials

    Attributes
    ----------

    layernum : int
        Layer
    materialid : int
        Material ID
    name : str
        Component or region name
    beginsta : float
        Begin station (m)
    endsta : float
        End station (m)
    maxwidth : float
        Max width (m)
    avgwidth : float
        Average width (m)
    area : float
        3D area (m^2)
    thickness : float
        Layer thickness (mm)
    weight : float
        Computed dry layer weight (g)
    """

    def __init__(self):
        self.layernum: int = None
        self.materialid: int = None
        self.name: str = None
        self.beginsta: float = None
        self.endsta: float = None
        self.maxwidth: float = None
        self.avgwidth: float = None
        self.area: float = None
        self.thickness: float = None
        self.weight: float = None


class Ply:
    """A simple class to organize the attributes of a ply

    Attributes
    ----------

    component : str
        parent component
    materialid : str
        Material id of ply
    thickness : float
        thickness of single ply (mm)
    angle : float
        ply angle
    nPlies : int
        number of plies
    """

    def __init__(self):
        self.component: str = None
        self.materialid: str = None
        self.thickness: float = None
        self.angle: float = None
        self.nPlies: int = None

    def _compare(self, other):
        pass


class SkinArea:
    def __init__(self):
        self.startIB = []
        self.endIB = []
        self.startOB = []
        self.endOB = []
        self.Material = []


class Stack:
    """A class definition for a stack of composite layers.

    Parameters
    ----------

    Attributes
    ----------
    name : string
        Name of the stack or composite material used by NuMAD, e.g. '000000_HP_LE_PANEL'
    indices : list
        Indices of stack, ``[in board station, out board station,
        1st kepoint, 2nd keypoint]``, e.g. ``[ibSta,obSta,keypt1,keypt2]``
    plygroups : list
        List of ``ply`` objects

    """

    def __init__(self):
        self.name: str = ""
        self.indices = []
        self.plygroups: list = []

    def addply(self, ply):
        """This method adds a Ply object to stack

        Parameters
        ----------
        ply : Ply object
            Ply object to be added
        Returns
        -------
        None

        """
        # modifies
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

    def layerThicknesses(self) -> ndarray:
        """Computes the thickness of each layer

        Returns:
            ndarray:
        """
        nLayers = len(self.plygroups)
        thicknesses = [
            self.plygroups[iLayer].nPlies * self.plygroups[iLayer].thickness
            for iLayer in range(nLayers)
        ]
        return np.array(thicknesses)
