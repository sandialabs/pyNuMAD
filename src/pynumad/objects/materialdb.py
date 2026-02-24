from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from pynumad.objects.stackdb import StackDatabase
from pynumad.objects.material import Layer  # re-exported from material.py


class MaterialDatabase:
    """Database of all blade materials in a solver-ready format.

    Entries are keyed by material/composite name and stored in
    ``entries``.  Use ``db[name]`` / ``db[name] = entry`` for convenient
    dict-style access.

    Attributes
    ----------
    entries : dict[str, MaterialDatabaseEntry]
        All material and composite stack entries indexed by name.
    """

    def __init__(self):
        self.entries: dict = {}

    # Dict-style access delegates to self.entries for convenience
    def __getitem__(self, key):
        return self.entries[key]

    def __setitem__(self, key, value):
        self.entries[key] = value

    def __contains__(self, key):
        return key in self.entries

    def generate(self, materials: dict, stackdb: StackDatabase):
        """Adds material and composites information to MatDB

        Parameters
        ----------
        materials : dict
        stackdb : StackDatabase

        Returns
        -------
        Self
        """
        mm_to_m = 0.001
        stacks = stackdb.stacks
        swstacks = stackdb.swstacks
        n_webs, n_stations = swstacks.shape

        # prepare material database ==========================================
        
        # add base materials
        for mat_name in materials:
            cur_entry = MaterialDatabaseEntry()
            cur_material = materials[mat_name]
            cur_entry.name = cur_material.name
            cur_entry.type = cur_material.type
            cur_entry.ex = cur_material.ex
            cur_entry.ey = cur_material.ey
            cur_entry.ez = cur_material.ez
            cur_entry.gxy = cur_material.gxy
            cur_entry.gyz = cur_material.gyz
            cur_entry.gxz = cur_material.gxz
            if cur_entry.type == "isotropic":
                cur_entry.nuxy = cur_material.prxy
            else:
                cur_entry.prxy = cur_material.prxy
                cur_entry.pryz = cur_material.pryz
                cur_entry.prxz = cur_material.prxz
            cur_entry.dens = cur_material.density
            cur_entry.reference = cur_material.reference
            self[mat_name] = cur_entry

        # add component stacks
        flat_stacks = stacks.flatten("F")
        for k_stack in range(flat_stacks.size):
            cur_stack = flat_stacks[k_stack]
            cur_entry = MaterialDatabaseEntry()
            cur_entry.name = cur_stack.name
            cur_entry.type = "composite"
            cur_entry.reference = "Reference text"
            cur_entry.thicknessType = "Constant"
            cur_entry.uniqueLayers = len(cur_stack.plygroups)
            cur_entry.symmetryType = "none"
            cur_entry.layer = [None] * cur_entry.uniqueLayers
            for j in range(cur_entry.uniqueLayers):
                cur_layer = Layer()
                plygroup = cur_stack.plygroups[j]
                matid = plygroup.materialid
                cur_layer.layer_name = self[matid].name
                cur_layer.thicknessA = mm_to_m * plygroup.thickness
                cur_layer.thicknessB = cur_layer.thicknessA
                cur_layer.quantity = plygroup.nPlies
                cur_layer.theta = plygroup.angle
                cur_entry.layer[j] = cur_layer
            self[cur_entry.name] = cur_entry

        # add shearweb stacks
        for k_web in range(n_webs):
            for k_stat in range(n_stations):
                cur_stack = swstacks[k_web, k_stat]
                cur_entry = MaterialDatabaseEntry()
                cur_entry.name = cur_stack.name
                cur_entry.type = "composite"
                cur_entry.reference = "Reference text"
                cur_entry.thicknessType = "Constant"
                try:
                    cur_entry.uniqueLayers = len(cur_stack.plygroups)
                except TypeError:
                    cur_entry.uniqueLayers = 0
                cur_entry.symmetryType = "none"
                cur_entry.layer = [None] * cur_entry.uniqueLayers
                for j in range(cur_entry.uniqueLayers):
                    cur_plygroup = cur_stack.plygroups[j]
                    cur_layer = Layer()
                    matid = cur_plygroup.materialid
                    cur_layer.layer_name = self[matid].name
                    cur_layer.thicknessA = (
                        mm_to_m * cur_plygroup.thickness
                    )
                    cur_layer.thicknessB = cur_layer.thicknessA
                    cur_layer.quantity = cur_plygroup.nPlies
                    cur_layer.theta = cur_plygroup.angle
                    cur_entry.layer[j] = cur_layer
                self[cur_entry.name] = cur_entry
        
        return self
    
    
@dataclass
class MaterialDatabaseEntry:
    """A simple class to organize the attributes of a material"""
    type: Optional[str] = None
    name: Optional[str] = None
    reference: Optional[str] = None
    dens: Optional[float] = None
    nuxy: Optional[float] = None
    ex: Optional[float] = None
    ey: Optional[float] = None
    ez: Optional[float] = None
    gxy: Optional[float] = None
    gyz: Optional[float] = None
    gxz: Optional[float] = None
    prxy: Optional[float] = None
    pryz: Optional[float] = None
    prxz: Optional[float] = None
    xten: Optional[float] = None
    xcmp: Optional[float] = None
    yten: Optional[float] = None
    ycmp: Optional[float] = None
    zten: Optional[float] = None
    zcmp: Optional[float] = None
    xy: Optional[float] = None
    yz: Optional[float] = None
    xz: Optional[float] = None
    xycp: Optional[float] = None
    yzcp: Optional[float] = None
    xzcp: Optional[float] = None
    xzit: Optional[float] = None
    xzic: Optional[float] = None
    yzit: Optional[float] = None
    yzic: Optional[float] = None
    g1g2: Optional[float] = None
    etal: Optional[float] = None
    etat: Optional[float] = None
    alp0: Optional[float] = None
    thicknessType: Optional[str] = None
    uniqueLayers: Optional[int] = None
    symmetryType: Optional[str] = None
    layer: Optional[list] = None
                 
# Layer is defined in pynumad.objects.material and re-exported here
# for backward compatibility.
__all__ = ["MaterialDatabase", "MaterialDatabaseEntry", "Layer"]


