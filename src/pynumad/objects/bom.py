
import numpy as np

from pynumad.objects.definition import Definition
from pynumad.objects.keypoints import KeyPoints

class BillOfMaterials(dict):
    """Bill of Materials class
    
    This class is designed replicate the bom
    blade attribute from NuMAD. It is currently
    not used in pyNuMAD analyses and isn't actively
    maintained.

    Attributes
    ----------
    indices : dict
    """
    def __init__(self):
        self.indices: dict = None
        pass
    
    def generate(self, definition: Definition, keypoints: KeyPoints):
        """This method generates the Bill-of-Materials
        See datatypes.BOM

        Returns
        -------

        None

        """
        # raise DeprecationWarning("update_bom currently deprecated. Please do not use.")
        # set conversion constants
        g_to_kg = 0.001
        m_to_mm = 1000.0
        mm_to_m = 0.001

        materials = definition.materials
        components = definition.components
        ispan = definition.ispan

        self["hp"] = []
        self["lp"] = []
        self["sw"] = []
        self["lebond"] = []
        self["tebond"] = []
        self["swbonds"] = []
        self["dryweight"] = []
        
        self.indices = {"hp": [], "lp": [], "sw": []}

        # calculate non-dimensional span
        ndspan = (ispan - ispan[0]) / (ispan[-1] - ispan[0])

        hprow = 0
        lprow = 0
        outer_shape_comps = [name for name in components if components[name].group == 0]
        for comp_name in outer_shape_comps:
            comp = components[comp_name]
            mat = materials[comp.materialid]
            hpRegion, lpRegion = comp.find_region_extents()
            num_layers = comp.get_num_layers(ndspan)
            num_layers = np.round(num_layers)

            for k_layer in range(1, int(np.max(num_layers)) + 1):
                begin_station, end_station = find_layer_extents(
                    num_layers, k_layer
                )
                ks_max = np.amin((len(begin_station), len(end_station)))
                # situation that beginSta/endSta is longer than 1
                for ks in range(ks_max):
                    ## END
                    if hpRegion:
                        areas = keypoints.key_areas[
                            hpRegion[0] : hpRegion[1],
                            begin_station[ks] : end_station[ks],
                        ]
                        regionarea = sum(areas.flatten())
                        arcs = (
                            keypoints.key_arcs[
                                hpRegion[1], begin_station[ks] : end_station[ks] + 1
                            ]
                            - keypoints.key_arcs[
                                hpRegion[0], begin_station[ks] : end_station[ks] + 1
                            ]
                        )
                        cur_bom = BillOfMaterialsEntry()
                        cur_bom.layernum = hprow
                        cur_bom.materialid = comp.materialid
                        cur_bom.name = comp.name
                        cur_bom.beginsta = ispan[begin_station[ks]]
                        cur_bom.endsta = ispan[end_station[ks]]
                        cur_bom.maxwidth = np.amax(arcs)
                        cur_bom.avgwidth = np.mean(arcs)
                        cur_bom.area = regionarea
                        cur_bom.thickness = mat.layerthickness
                        cur_bom.weight = mat.drydensity * regionarea
                        self.indices["hp"].append(
                            [begin_station[ks], end_station[ks], *hpRegion]
                        )
                        self["hp"].append(cur_bom)
                        hprow = hprow + 1

                    if lpRegion:
                        areas = keypoints.key_areas[
                            lpRegion[0] : lpRegion[1],
                            begin_station[ks] : end_station[ks],
                        ]
                        regionarea = sum(areas.flatten())
                        arcs = (
                            keypoints.key_arcs[
                                lpRegion[1], begin_station[ks] : end_station[ks] + 1
                            ]
                            - keypoints.key_arcs[
                                lpRegion[0], begin_station[ks] : end_station[ks] + 1
                            ]
                        )
                        cur_bom = BillOfMaterialsEntry()
                        cur_bom.layernum = lprow
                        cur_bom.materialid = comp.materialid
                        cur_bom.name = comp.name
                        cur_bom.beginsta = ispan[begin_station[ks]]
                        cur_bom.endsta = ispan[end_station[ks]]
                        cur_bom.maxwidth = np.amax(arcs)
                        cur_bom.avgwidth = np.mean(arcs)
                        cur_bom.area = regionarea
                        cur_bom.thickness = mat.layerthickness
                        cur_bom.weight = mat.drydensity * regionarea
                        self.indices["lp"].append(
                            [begin_station[ks], end_station[ks], *lpRegion]
                        )
                        self["lp"].append(cur_bom)
                        lprow = lprow + 1

        # shearwebs
        swnum = None
        swrow = 0
        sw_begin_station = []
        sw_end_station = []
        sw_comps = [comp for comp in components.values() if comp.group > 0]

        def sorter(e):
            return e.group

        # enforce an ordering on components based on group number
        # to ensure below loop functions correctly
        # probably should re-write the loop at some point...
        sw_comps.sort(key=sorter)
        for comp in sw_comps:
            mat = materials[comp.materialid]
            num_layers = comp.get_num_layers(ndspan)
            num_layers = np.round(num_layers)

            for k_layer in range(1, int(np.max(num_layers)) + 1):
                begin_station, end_station = find_layer_extents(
                    num_layers, k_layer
                )
                ks_max = np.amin((len(begin_station), len(end_station)))
                # situation that beginSta/endSta is longer than 1
                for ks in range(ks_max):
                    if swnum != comp.group - 1:
                        swnum = comp.group - 1
                        swrow = 0
                        sw_begin_station.append(begin_station[0])
                        sw_end_station.append(end_station[0])
                        self["sw"].append([])
                        self.indices["sw"].append([])
                    sw_begin_station[swnum] = np.amin(
                        [*begin_station, sw_begin_station[swnum]]
                    )
                    sw_end_station[swnum] = np.amax(
                        [*end_station, sw_end_station[swnum]]
                    )
                    areas = keypoints.web_areas[swnum][
                        begin_station[ks] : end_station[ks]
                    ]
                    regionarea = sum(areas.flatten())
                    cur_bom = BillOfMaterialsEntry()
                    cur_bom.layernum = swrow
                    cur_bom.materialid = comp.materialid
                    cur_bom.name = comp.name
                    cur_bom.beginsta = ispan[begin_station[ks]]
                    cur_bom.endsta = ispan[end_station[ks]]
                    cur_bom.maxwidth = np.amax(keypoints.web_width[swnum])
                    cur_bom.avgwidth = np.mean(keypoints.web_width[swnum])
                    cur_bom.area = regionarea
                    cur_bom.thickness = mat.layerthickness
                    cur_bom.weight = mat.drydensity * regionarea
                    self["sw"][swnum].append(cur_bom)
                    self.indices["sw"][swnum].append(
                        [begin_station[ks], end_station[ks]]
                    )
                    swrow = swrow + 1

        # compute lebond, tebond, and dryweight
        self["lebond"] = sum(keypoints.le_bond) * m_to_mm
        self["tebond"] = sum(keypoints.te_bond) * m_to_mm
        hp_dw = sum([L.weight for L in self["hp"]])
        lp_dw = sum([L.weight for L in self["lp"]])
        self["dryweight"] = g_to_kg * (hp_dw + lp_dw)

        nsw = len(self["sw"])
        self["swbonds"] = [None] * nsw
        for k in range(nsw):
            sw_dw = sum([L.weight for L in self["sw"][k]])
            self["dryweight"] = self["dryweight"] + sw_dw
            C = keypoints.web_bonds[k][:, sw_begin_station[k] : sw_end_station[k]]
            self["swbonds"][k] = m_to_mm * np.sum(C, 1)


        return self
    
    # Supporting function for update_bom
def find_layer_extents(layer_dist, layer_n):
    """
    TODO docstring
    """
    assert np.isscalar(layer_n), 'second argument "layer_n" must be a scalar'
    sta_logical = layer_dist >= layer_n
    prev = 0
    begin_station = []
    end_station = []
    for k in range(len(sta_logical)):
        if sta_logical[k] == 1 and prev == 0:
            begin_station.append(k)
        if sta_logical[k] == 0 and prev == 1:
            end_station.append(k)
        elif k == len(sta_logical) - 1 and prev == 1:
            end_station.append(k)
        prev = sta_logical[k]

    return begin_station, end_station

    
class BillOfMaterialsEntry:
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