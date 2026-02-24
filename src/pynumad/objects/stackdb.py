from copy import deepcopy

import numpy as np
from numpy import ndarray

from pynumad.objects.keypoints import KeyPoints
from pynumad.objects.bom import BillOfMaterials
from pynumad.objects.stack import Stack, Ply

# Stack and Ply are defined in pynumad.objects.stack and re-exported here
# for backward compatibility with code that imports them from stackdb.
__all__ = ["StackDatabase", "Stack", "Ply"]


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
        """Build material stacks for every surface region and shear web.

        Parameters
        ----------
        keypoints : KeyPoints
        bom : BillOfMaterials
            Must have been generated prior to calling this method.
        """
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

        # HP layers — iterate over DataFrame rows
        for _, row in bom.hp.iterrows():
            cur_ply = Ply(
                component=row["name"],
                materialid=row["materialid"],
                thickness=row["thickness"],
                angle=row["angle"],
                nPlies=1,
            )
            for k_seg in range(int(row["seg_start"]), int(row["seg_end"])):
                for k_stat in range(int(row["sta_begin_idx"]), int(row["sta_end_idx"])):
                    self.stacks[k_seg, k_stat].addply(deepcopy(cur_ply))

        # LP layers
        for _, row in bom.lp.iterrows():
            cur_ply = Ply(
                component=row["name"],
                materialid=row["materialid"],
                thickness=row["thickness"],
                angle=row["angle"],
                nPlies=1,
            )
            for k_seg in range(int(row["seg_start"]), int(row["seg_end"])):
                for k_stat in range(int(row["sta_begin_idx"]), int(row["sta_end_idx"])):
                    self.stacks[k_seg, k_stat].addply(deepcopy(cur_ply))

        # Shear-web stacks
        self.swstacks = np.empty(shape=(n_webs, n_stations), dtype=object)
        for k_web in range(n_webs):
            for k_stat in range(n_stations):
                self.swstacks[k_web, k_stat] = Stack()
                self.swstacks[k_web, k_stat].name = "{:02d}_{:02d}_SW".format(
                    k_web, k_stat
                )
                ind = keypoints.web_indices[k_web]
                self.swstacks[k_web][k_stat].indices = [
                    k_stat,
                    k_stat + 1,
                    ind[0],
                    ind[1],
                ]

        # Shear-web layers
        for _, row in bom.sw.iterrows():
            k_web = int(row["web_id"])
            cur_ply = Ply(
                component=row["name"],
                materialid=row["materialid"],
                thickness=row["thickness"],
                angle=0,
                nPlies=1,
            )
            for k_stat in range(int(row["sta_begin_idx"]), int(row["sta_end_idx"])):
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
                    t1 = ply1.thickness
                    ply2.thickness = 0.3333333 * t1
                    ply1.thickness = 0.3333333 * t1
                    ply3.thickness = 0.3333333 * t1
                    newPg = np.array([ply1, ply2, ply3])
                    self.swstacks[i][j].plygroups = newPg
        return self
