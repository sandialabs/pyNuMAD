from __future__ import annotations

import numpy as np
import pandas as pd

from pynumad.objects.keypoints import KeyPoints


# Column names shared by hp, lp, and sw DataFrames.
_BOM_COLUMNS = [
    "layernum",
    "materialid",
    "name",
    "beginsta",
    "endsta",
    "maxwidth",
    "avgwidth",
    "area",
    "thickness",
    "weight",
    "angle",
    # Integer index bounds consumed by StackDatabase.generate()
    "sta_begin_idx",
    "sta_end_idx",
    "seg_start",
    "seg_end",
]

# The sw DataFrame adds one extra column.
_SW_COLUMNS = _BOM_COLUMNS + ["web_id"]


class BillOfMaterials:
    """Bill of Materials class.

    Contains per-layer material information for each component region of the
    blade (HP shell, LP shell, and each shear web), stored as
    :class:`pandas.DataFrame` objects.

    Attributes
    ----------
    hp : pd.DataFrame
        HP-side shell layers.  Columns: ``layernum``, ``materialid``,
        ``name``, ``beginsta``, ``endsta``, ``maxwidth``, ``avgwidth``,
        ``area``, ``thickness``, ``weight``, ``angle``,
        ``sta_begin_idx``, ``sta_end_idx``, ``seg_start``, ``seg_end``.
    lp : pd.DataFrame
        LP-side shell layers (same columns as ``hp``).
    sw : pd.DataFrame
        Shear-web layers.  Same columns as ``hp`` plus ``web_id`` (int).
    lebond : float
        Total leading-edge bond length [mm].
    tebond : float
        Total trailing-edge bond length [mm].
    swbonds : list[ndarray]
        Per-web bond-line length arrays [mm].
    dryweight : float
        Total dry material weight [kg].
    """

    def __init__(self):
        self.hp: pd.DataFrame = pd.DataFrame(columns=_BOM_COLUMNS)
        self.lp: pd.DataFrame = pd.DataFrame(columns=_BOM_COLUMNS)
        self.sw: pd.DataFrame = pd.DataFrame(columns=_SW_COLUMNS)
        self.lebond: float = None
        self.tebond: float = None
        self.swbonds: list = []
        self.dryweight: float = None

    def generate(self, ispan, components, materials, keypoints: KeyPoints):
        """Generate the Bill-of-Materials.

        Parameters
        ----------
        ispan : ndarray
            Interpolated span station locations [m].
        components : dict
            Component definitions (from ``Definition.components``).
        materials : dict
            Material definitions (from ``Definition.materials``).
        keypoints : KeyPoints

        Returns
        -------
        self
        """
        # set conversion constants
        g_to_kg = 0.001
        m_to_mm = 1000.0
        mm_to_m = 0.001

        # reset
        self.hp = pd.DataFrame(columns=_BOM_COLUMNS)
        self.lp = pd.DataFrame(columns=_BOM_COLUMNS)
        self.sw = pd.DataFrame(columns=_SW_COLUMNS)
        self.lebond = None
        self.tebond = None
        self.swbonds = []
        self.dryweight = None

        # calculate non-dimensional span
        ndspan = (ispan - ispan[0]) / (ispan[-1] - ispan[0])

        hp_rows = []
        lp_rows = []
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
                for ks in range(ks_max):
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
                        hp_rows.append({
                            "layernum": hprow,
                            "materialid": comp.materialid,
                            "name": comp.name,
                            "beginsta": ispan[begin_station[ks]],
                            "endsta": ispan[end_station[ks]],
                            "maxwidth": np.amax(arcs),
                            "avgwidth": np.mean(arcs),
                            "area": regionarea,
                            "thickness": mat.layerthickness,
                            "weight": mat.drydensity * regionarea,
                            "angle": comp.fabricangle * 180 / 3.141592653589793,
                            "sta_begin_idx": begin_station[ks],
                            "sta_end_idx": end_station[ks],
                            "seg_start": hpRegion[0],
                            "seg_end": hpRegion[1],
                        })
                        hprow += 1

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
                        lp_rows.append({
                            "layernum": lprow,
                            "materialid": comp.materialid,
                            "name": comp.name,
                            "beginsta": ispan[begin_station[ks]],
                            "endsta": ispan[end_station[ks]],
                            "maxwidth": np.amax(arcs),
                            "avgwidth": np.mean(arcs),
                            "area": regionarea,
                            "thickness": mat.layerthickness,
                            "weight": mat.drydensity * regionarea,
                            "angle": -1 * comp.fabricangle * 180 / 3.141592653589793,
                            "sta_begin_idx": begin_station[ks],
                            "sta_end_idx": end_station[ks],
                            "seg_start": lpRegion[0],
                            "seg_end": lpRegion[1],
                        })
                        lprow += 1

        self.hp = pd.DataFrame(hp_rows, columns=_BOM_COLUMNS)
        self.lp = pd.DataFrame(lp_rows, columns=_BOM_COLUMNS)

        # shearwebs
        swnum = None
        swrow = 0
        sw_begin_station = []
        sw_end_station = []
        sw_rows = []
        sw_comps = [comp for comp in components.values() if comp.group > 0]

        def sorter(e):
            return e.group

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
                for ks in range(ks_max):
                    if swnum != comp.group - 1:
                        swnum = comp.group - 1
                        swrow = 0
                        sw_begin_station.append(begin_station[0])
                        sw_end_station.append(end_station[0])
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
                    sw_rows.append({
                        "layernum": swrow,
                        "materialid": comp.materialid,
                        "name": comp.name,
                        "beginsta": ispan[begin_station[ks]],
                        "endsta": ispan[end_station[ks]],
                        "maxwidth": np.amax(keypoints.web_width[swnum]),
                        "avgwidth": np.mean(keypoints.web_width[swnum]),
                        "area": regionarea,
                        "thickness": mat.layerthickness,
                        "weight": mat.drydensity * regionarea,
                        "angle": comp.fabricangle * 180 / 3.141592653589793,
                        "sta_begin_idx": begin_station[ks],
                        "sta_end_idx": end_station[ks],
                        "seg_start": 0,  # not used for shear webs
                        "seg_end": 0,    # not used for shear webs
                        "web_id": swnum,
                    })
                    swrow += 1

        self.sw = pd.DataFrame(sw_rows, columns=_SW_COLUMNS)

        # compute lebond, tebond, and dryweight
        self.lebond = sum(keypoints.le_bond) * m_to_mm
        self.tebond = sum(keypoints.te_bond) * m_to_mm

        hp_dw = self.hp["weight"].sum() if not self.hp.empty else 0.0
        lp_dw = self.lp["weight"].sum() if not self.lp.empty else 0.0
        self.dryweight = g_to_kg * (hp_dw + lp_dw)

        nsw = len(sw_begin_station)
        self.swbonds = [None] * nsw
        for k in range(nsw):
            sw_web = self.sw[self.sw["web_id"] == k]
            sw_dw = sw_web["weight"].sum() if not sw_web.empty else 0.0
            self.dryweight = self.dryweight + sw_dw
            C = keypoints.web_bonds[k][:, sw_begin_station[k] : sw_end_station[k]]
            self.swbonds[k] = m_to_mm * np.sum(C, 1)

        return self


def find_layer_extents(layer_dist, layer_n):
    """Find the begin and end station indices where a layer count reaches *layer_n*.

    Parameters
    ----------
    layer_dist : array-like
        Per-station layer count distribution.
    layer_n : scalar
        Layer threshold to find extents for.

    Returns
    -------
    begin_station : list[int]
    end_station : list[int]
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
