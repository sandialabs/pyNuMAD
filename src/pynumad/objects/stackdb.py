from copy import deepcopy

import numpy as np
from numpy import ndarray

from pynumad.objects.keypoints import KeyPoints
from pynumad.objects.bom import BillOfMaterials

class StackDatabase:
    def __init__(self):
        self.stacks: ndarray = None
        self.swstacks: ndarray = None
        
    def __eq__(self, other):
        assert self.stacks.shape == other.stacks.shape
        
        assert self.swstacks.shape == other.swstacks.shape
        
        self_flat_stacks = self.stacks.flatten()
        other_flat_stacks = other.stacks.flatten()
        for i in range(self_flat_stacks.shape[0]):
            self_stack = self_flat_stacks[i]
            other_stack = other_flat_stacks[i]
            if self_stack != other_stack:
                return False
        self_flat_stacks = self.swstacks.flatten()
        other_flat_stacks = other.swstacks.flatten()
        for i in range(self_flat_stacks.shape[0]):
            self_stack = self_flat_stacks[i]
            other_stack = other_flat_stacks[i]
            if self_stack != other_stack:
                return False
        return True
    def generate(self, keypoints: KeyPoints, bom: BillOfMaterials):
        # build the material stack for each area
        n_segments = keypoints.key_areas.shape[0]
        n_stations = keypoints.key_areas.shape[1]
        n_webs = len(keypoints.web_points)
        segment_labels = [
            "HP_TE_FLAT",
            "HP_TE_REINF",
            "HP_TE_PANEL",
            "HP_SPAR",
            "HP_LE_PANEL",
            "HP_LE",
            "LP_LE",
            "LP_LE_PANEL",
            "LP_SPAR",
            "LP_TE_PANEL",
            "LP_TE_REINF",
            "LP_TE_FLAT",
        ]

        # initialize stack array
        self.stacks = np.empty(shape=(n_segments, n_stations), dtype=object)
        for k_seg in range(n_segments):
            for k_stat in range(n_stations):
                self.stacks[k_seg, k_stat] = Stack()
                self.stacks[k_seg, k_stat].name = "{:02d}_{:02d}_{}".format(
                    k_seg, k_stat, segment_labels[k_seg]
                )
                self.stacks[k_seg, k_stat].indices = [
                    k_stat,
                    k_stat + 1,
                    k_seg,
                    k_seg + 1,
                ]

        for k in range(len(bom["hp"])):
            # for each row in the BOM, get the ply definition ...
            cur_ply = Ply()
            cur_ply.component = bom["hp"][k].name
            cur_ply.materialid = bom["hp"][k].materialid
            cur_ply.thickness = bom["hp"][k].thickness
            cur_ply.angle = 0  # TODO, set to 0 for now, bom['lp'](k, );
            cur_ply.nPlies = 1  # default to 1, modified in addply() if necessary

            # ... and add the ply to every area that is part of the region
            ind = bom.indices["hp"][k]
            for k_seg in range(ind[2], ind[3]):
                for k_stat in range(ind[0], ind[1]):
                    # deepcopy is important to keep make ply object unique in each stack
                    self.stacks[k_seg, k_stat].addply(deepcopy(cur_ply))  

        for k in range(len(bom["lp"])):
            # for each row in the BOM, get the ply definition ...
            cur_ply = Ply()
            cur_ply.component = bom["lp"][k].name
            cur_ply.materialid = bom["lp"][k].materialid
            cur_ply.thickness = bom["lp"][k].thickness
            cur_ply.angle = 0  # TODO, set to 0 for now, bom['lp'](k, );
            cur_ply.nPlies = 1  # default to 1, modified in addply() if necessary

            # ... and add the ply to every area that is part of the region
            ind = bom.indices["lp"][k]
            for k_seg in range(ind[2], ind[3]):
                for k_stat in range(ind[0], ind[1]):
                    self.stacks[k_seg, k_stat].addply(deepcopy(cur_ply))

        self.swstacks = np.empty(shape=(n_webs, n_stations), dtype=object)
        for k_web in range(n_webs):
            for k_stat in range(n_stations):
                self.swstacks[k_web, k_stat] = Stack()
                # name the stacks <webnumber+1>_<stationnumber+1>_SW
                self.swstacks[k_web, k_stat].name = "{:02d}_{:02d}_SW".format(
                    k_web, k_stat
                )
                # currently, the shearweb indices do not change down the span
                ind = keypoints.web_indices[k_web]
                self.swstacks[k_web][k_stat].indices = [
                    k_stat,
                    k_stat + 1,
                    ind[0],
                    ind[1],
                ]
        for k_web in range(n_webs):
            for k in range(len(bom["sw"][k_web])):
                # for each row in the BOM, get the ply definition ...
                cur_ply = Ply()
                cur_ply.component = bom["sw"][k_web][k].name
                cur_ply.materialid = bom["sw"][k_web][k].materialid
                cur_ply.thickness = bom["sw"][k_web][k].thickness
                cur_ply.angle = 0  # TODO, set to 0 for now, bom['lp'](k, );
                cur_ply.nPlies = 1  # default to 1, modified in addply() if necessary
                
                # ... and add the ply to every area that is part of the region
                ind = bom.indices["sw"][k_web][k]
                for k_stat in range(ind[0], ind[1]):
                    self.swstacks[k_web, k_stat].addply(deepcopy(cur_ply))
        
    def edit_stacks_for_solid_mesh(self):
        """

        Returns
        -------
        Self
        """
        numSec, numStat = self.stacks.shape
        for i in range(numSec):
            for j in range(numStat):
                pg = self.stacks[i, j].plygroups
                if len(pg) == 4:
                    ply1 = deepcopy(pg[1])
                    ply2 = deepcopy(pg[2])
                    ply3 = deepcopy(pg[3])
                    newPg = np.array([ply1, ply2, ply3])
                else:
                    if len(pg) == 3:
                        # newPg = np.array([pg[1],pg[1],pg[2]])
                        ply1 = deepcopy(pg[1])
                        ply2 = deepcopy(pg[1])
                        ply3 = deepcopy(pg[2])
                        t2 = ply1.thickness
                        t3 = ply3.thickness
                        ply2.thickness = 0.3333333 * (t2 + t3)
                        ply1.thickness = 0.6666666 * t2
                        ply3.thickness = 0.6666666 * t3
                        newPg = np.array([ply1, ply2, ply3])
                    else:
                        if len(pg) == 2:
                            ply1 = deepcopy(pg[1])
                            ply2 = deepcopy(pg[1])
                            ply3 = deepcopy(pg[1])
                            # newPg = np.array([pg[0],pg[0],pg[1]])
                            t1 = ply1.thickness
                            t2 = ply3.thickness
                            ply2.thickness = 0.3333333 * (t1 + t2)
                            ply1.thickness = 0.6666666 * t1
                            ply3.thickness = 0.6666666 * t2
                            newPg = np.array([ply1, ply2, ply3])
                        else:
                            ply1 = deepcopy(pg[0])
                            ply2 = deepcopy(pg[0])
                            ply3 = deepcopy(pg[0])
                            # newPg = np.array([pg[0],pg[0],pg[0]])
                            t1 = ply1.thickness
                            ply2.thickness = 0.3333333 * t1
                            ply1.thickness = 0.3333333 * t1
                            ply3.thickness = 0.3333333 * t1
                            newPg = np.array([ply1, ply2, ply3])
                self.stacks[i, j].plygroups = newPg

        for i in range(2):
            stackLst = self.swstacks[i]
            for j in range(len(stackLst)):
                pg = stackLst[j].plygroups
                if len(pg) == 2:
                    ply1 = deepcopy(pg[0])
                    ply2 = deepcopy(pg[0])
                    ply3 = deepcopy(pg[1])
                    # newPg = np.array([pg[0],pg[0],pg[1]])
                    t1 = ply1.thickness
                    t2 = ply3.thickness
                    ply2.thickness = 0.3333333 * (t1 + t2)
                    ply1.thickness = 0.6666666 * t1
                    ply3.thickness = 0.6666666 * t2
                    newPg = np.array([ply1, ply2, ply3])
                    self.swstacks[i][j].plygroups = newPg
                elif len(pg) == 1:
                    ply1 = deepcopy(pg[0])
                    ply2 = deepcopy(pg[0])
                    ply3 = deepcopy(pg[0])
                    # newPg = np.array([pg[0],pg[0],pg[0]])
                    t1 = ply1.thickness
                    ply2.thickness = 0.3333333 * t1
                    ply1.thickness = 0.3333333 * t1
                    ply3.thickness = 0.3333333 * t1
                    newPg = np.array([ply1, ply2, ply3])
                    self.swstacks[i][j].plygroups = newPg
        return self
    
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
        
    def __eq__(self, other):
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

    def layer_thicknesses(self) -> ndarray:
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

    def __eq__(self, other):
        attrs = vars(self).keys()
        for attr in attrs:
            self_attr = getattr(self, attr)
            other_attr = getattr(other, attr)
            if isinstance(self_attr, (int, float, str)):
                if self_attr != other_attr:
                    return False
        return True