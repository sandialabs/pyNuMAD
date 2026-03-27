from __future__ import annotations

import numpy as np
import pandas as pd

from pynumad.objects.keypoints import KeyPoints

_MOLD_COLS = ['Component', 'Area', 'Length', 'Start_station', 'End_station',
              'Max_width', 'Avg_width', 'Perimeter']

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
    "dryweight",
    "curedweight",
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
    mold : dict
        Dict of DataFrames with mold/prefab dimensional data (populated by
        generate_segmented). Keys: 'shell', 'spar', 'LEreinf', 'TEreinf',
        'root', 'web'.
    lebond : float
        Total leading-edge bond length [mm].
    tebond : float
        Total trailing-edge bond length [mm].
    swbonds : list[ndarray]
        Per-web bond-line length arrays [mm].
    total_dryweight : float
        Total dry (pre-cure fiber) weight [kg], computed from ``drydensity`` [g/m²].
    total_curedweight : float
        Total cured composite weight [kg], computed from ``density`` [kg/m³]
        times ``layerthickness`` [m] times area [m²].
    """

    def __init__(self):
        self.hp: pd.DataFrame = pd.DataFrame(columns=_BOM_COLUMNS)
        self.lp: pd.DataFrame = pd.DataFrame(columns=_BOM_COLUMNS)
        self.sw: pd.DataFrame = pd.DataFrame(columns=_SW_COLUMNS)
        self.mold: dict = {}
        self.lebond: float = None
        self.tebond: float = None
        self.swbonds: list = []
        self.total_dryweight: float = None
        self.total_curedweight: float = None

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
        self.total_dryweight = None
        self.total_curedweight = None

        # calculate non-dimensional span and interval midpoints
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
                            "dryweight": mat.drydensity * regionarea,
                            "curedweight": mat.density * mat.layerthickness * mm_to_m * regionarea,
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
                            "dryweight": mat.drydensity * regionarea,
                            "curedweight": mat.density * mat.layerthickness * mm_to_m * regionarea,
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
                        "dryweight": mat.drydensity * regionarea,
                        "curedweight": mat.density * mat.layerthickness * mm_to_m * regionarea,
                        "angle": comp.fabricangle * 180 / 3.141592653589793,
                        "sta_begin_idx": begin_station[ks],
                        "sta_end_idx": end_station[ks],
                        "seg_start": 0,  # not used for shear webs
                        "seg_end": 0,    # not used for shear webs
                        "web_id": swnum,
                    })
                    swrow += 1

        self.sw = pd.DataFrame(sw_rows, columns=_SW_COLUMNS)

        # compute lebond, tebond, total_dryweight, and total_curedweight
        self.lebond = sum(keypoints.le_bond) * m_to_mm
        self.tebond = sum(keypoints.te_bond) * m_to_mm

        hp_dw = self.hp["dryweight"].sum() if not self.hp.empty else 0.0
        lp_dw = self.lp["dryweight"].sum() if not self.lp.empty else 0.0
        self.total_dryweight = hp_dw + lp_dw

        hp_tw = self.hp["curedweight"].sum() if not self.hp.empty else 0.0
        lp_tw = self.lp["curedweight"].sum() if not self.lp.empty else 0.0
        self.total_curedweight = hp_tw + lp_tw

        nsw = len(sw_begin_station)
        self.swbonds = [None] * nsw
        for k in range(nsw):
            sw_web = self.sw[self.sw["web_id"] == k]
            sw_dw = sw_web["dryweight"].sum() if not sw_web.empty else 0.0
            sw_tw = sw_web["curedweight"].sum() if not sw_web.empty else 0.0
            self.total_dryweight = self.total_dryweight + sw_dw
            self.total_curedweight = self.total_curedweight + sw_tw
            C = keypoints.web_bonds[k][:, sw_begin_station[k] : sw_end_station[k]]
            self.swbonds[k] = m_to_mm * np.sum(C, 1)

        return self

    def generate_segmented(self, segments, ispan, components, materials, keypoints):
        """Generate a segmented Bill-of-Materials.

        Parameters
        ----------
        segments : list of dict
            Each entry describes one blade segment with keys:
              'span_nd'    : [start, end] normalized spanwise extents
              'hp_extents' : [key1, key2] HP chordwise bounds (e.g. ['le','te'])
              'lp_extents' : [key1, key2] LP chordwise bounds
        ispan : np.ndarray
            Interpolated spanwise stations (m)
        components : dict
            Component name -> Component objects
        materials : dict
            Material id -> Material objects
        keypoints : KeyPoints
            Computed keypoints (must already have key_bonds)
        """
        g_to_kg = 0.001
        m_to_mm = 1000.0
        mm_to_m = 0.001

        ndspan = (ispan - ispan[0]) / (ispan[-1] - ispan[0])
        n_segments = len(segments)
        key_labels = keypoints.key_labels

        hp_rows = []
        lp_rows = []
        sw_rows = []
        hprow = 0
        lprow = 0

        self.mold = {}
        self.lebond = []
        self.tebond = []

        n_webs = len(keypoints.web_areas)
        self.swbonds = [np.zeros((2, n_segments)) for _ in range(n_webs)]

        for ks, seg in enumerate(segments):
            span_nd = seg['span_nd']
            seg_hp_extents = seg['hp_extents']
            seg_lp_extents = seg['lp_extents']

            begin_sta = int(np.argmin(np.abs(ndspan - span_nd[0])))
            end_sta = int(np.argmin(np.abs(ndspan - span_nd[1])))
            station_ids = list(range(begin_sta, end_sta + 1))

            seg_data = {
                'begin_sta': begin_sta,
                'end_sta': end_sta,
                'hp_extents': seg_hp_extents,
                'lp_extents': seg_lp_extents,
                'sw_begin': {},
                'sw_end': {},
                'sw_extents': {},
                'sw_groups': {},
                'sw_labels': {},
            }

            for comp_name, comp in components.items():
                mat = materials[comp.materialid]
                hp_region, lp_region, sw_region = _find_region_extents_segmented(
                    key_labels, comp, seg_hp_extents, seg_lp_extents
                )
                num_layers = comp.get_num_layers(ndspan[begin_sta:end_sta + 1])
                num_layers = np.round(num_layers)

                seg_data = _calc_component_stations(comp, seg_data, num_layers, station_ids)

                klay_max = int(np.max(num_layers[:-1])) if len(num_layers) > 1 else 0

                for klay in range(1, klay_max + 1):
                    begin_stas, end_stas = _find_layer_extents_segmented(
                        num_layers, klay, station_ids
                    )
                    for ii in range(len(begin_stas)):
                        b_sta = begin_stas[ii]
                        e_sta = end_stas[ii]

                        if comp.group == 0 and hp_region:
                            areas = keypoints.key_areas[
                                hp_region[0]:hp_region[1], b_sta:e_sta
                            ]
                            arcs = (
                                keypoints.key_arcs[hp_region[1], b_sta:e_sta + 1]
                                - keypoints.key_arcs[hp_region[0], b_sta:e_sta + 1]
                            )
                            lengths = keypoints.key_bonds[
                                [hp_region[0], hp_region[1]], b_sta:e_sta
                            ]
                            arc_length = float(np.sum(lengths) / 2)
                            perimeter = float(np.sum(lengths) + arcs[0] + arcs[-1])
                            hp_rows.append({
                                'layernum': hprow,
                                'materialid': comp.materialid,
                                'name': comp.name,
                                'beginsta': ispan[b_sta],
                                'endsta': ispan[e_sta],
                                'maxwidth': float(np.amax(arcs)),
                                'avgwidth': float(np.mean(arcs)),
                                'area': float(np.sum(areas)),
                                'thickness': mat.layerthickness,
                                'dryweight': mat.drydensity * float(np.sum(areas)),
                                'curedweight': mat.density * mat.layerthickness * mm_to_m * float(np.sum(areas)),
                                'angle': comp.fabricangle * 180 / np.pi,
                                'sta_begin_idx': b_sta,
                                'sta_end_idx': e_sta,
                                'seg_start': hp_region[0],
                                'seg_end': hp_region[1],
                                'arc_length': arc_length,
                                'perimeter': perimeter,
                                'segment_id': f'HP{ks + 1}',
                            })
                            hprow += 1

                        if comp.group == 0 and lp_region:
                            areas = keypoints.key_areas[
                                lp_region[0]:lp_region[1], b_sta:e_sta
                            ]
                            arcs = (
                                keypoints.key_arcs[lp_region[1], b_sta:e_sta + 1]
                                - keypoints.key_arcs[lp_region[0], b_sta:e_sta + 1]
                            )
                            lengths = keypoints.key_bonds[
                                [lp_region[0], lp_region[1]], b_sta:e_sta
                            ]
                            arc_length = float(np.sum(lengths) / 2)
                            perimeter = float(np.sum(lengths) + arcs[0] + arcs[-1])
                            lp_rows.append({
                                'layernum': lprow,
                                'materialid': comp.materialid,
                                'name': comp.name,
                                'beginsta': ispan[b_sta],
                                'endsta': ispan[e_sta],
                                'maxwidth': float(np.amax(arcs)),
                                'avgwidth': float(np.mean(arcs)),
                                'area': float(np.sum(areas)),
                                'thickness': mat.layerthickness,
                                'dryweight': mat.drydensity * float(np.sum(areas)),
                                'curedweight': mat.density * mat.layerthickness * mm_to_m * float(np.sum(areas)),
                                'angle': -comp.fabricangle * 180 / np.pi,
                                'sta_begin_idx': b_sta,
                                'sta_end_idx': e_sta,
                                'seg_start': lp_region[0],
                                'seg_end': lp_region[1],
                                'arc_length': arc_length,
                                'perimeter': perimeter,
                                'segment_id': f'LP{ks + 1}',
                            })
                            lprow += 1

                        if comp.group > 0 and sw_region:
                            kw = comp.group - 1
                            areas = keypoints.web_areas[kw][b_sta:e_sta]
                            arcs = keypoints.web_width[kw][b_sta:e_sta + 1]
                            lengths = keypoints.web_bonds[kw][:, b_sta:e_sta]
                            arc_length = float(np.sum(lengths) / 2)
                            perimeter = float(np.sum(lengths) + arcs[0] + arcs[-1])

                            if kw not in seg_data['sw_labels']:
                                le = key_labels.index('le')
                                kl_surf = key_labels[:le + 1]
                                le_web_idxs = [
                                    i for i, k in enumerate(kl_surf) if k in ('a', 'b')
                                ]
                                web_kl_idx = kl_surf.index(comp.hpextents[0])
                                if web_kl_idx >= min(le_web_idxs):
                                    sw_label = f'LE{ks + 1}_SW{comp.group}'
                                else:
                                    sw_label = f'TE{ks + 1}_SW{comp.group}'
                                seg_data['sw_labels'][kw] = sw_label
                                seg_data['sw_extents'][kw] = [
                                    comp.hpextents[0], comp.lpextents[0]
                                ]
                                seg_data['sw_groups'][kw] = comp.group
                                seg_data['sw_begin'][kw] = b_sta
                                seg_data['sw_end'][kw] = e_sta
                            else:
                                seg_data['sw_begin'][kw] = min(
                                    seg_data['sw_begin'][kw], b_sta
                                )
                                seg_data['sw_end'][kw] = max(
                                    seg_data['sw_end'][kw], e_sta
                                )

                            sw_rows.append({
                                'layernum': len(sw_rows),
                                'materialid': comp.materialid,
                                'name': comp.name,
                                'beginsta': ispan[b_sta],
                                'endsta': ispan[e_sta],
                                'maxwidth': float(np.amax(arcs)),
                                'avgwidth': float(np.mean(arcs)),
                                'area': float(np.sum(areas)),
                                'thickness': mat.layerthickness,
                                'dryweight': mat.drydensity * float(np.sum(areas)),
                                'curedweight': mat.density * mat.layerthickness * mm_to_m * float(np.sum(areas)),
                                'angle': comp.fabricangle * 180 / np.pi,
                                'sta_begin_idx': b_sta,
                                'sta_end_idx': e_sta,
                                'web_id': kw,
                                'arc_length': arc_length,
                                'perimeter': perimeter,
                                'segment_id': seg_data['sw_labels'][kw],
                            })

            # LE/TE bond lengths (only for segments that include LE/TE boundaries)
            le_in_hp = 'le' in [x.lower() for x in seg_hp_extents]
            te_in_hp = 'te' in [x.lower() for x in seg_hp_extents]
            le_bond_ks = float(np.sum(keypoints.le_bond[begin_sta:end_sta])) * m_to_mm
            te_bond_ks = float(np.sum(keypoints.te_bond[begin_sta:end_sta])) * m_to_mm
            self.lebond.append(le_bond_ks if le_in_hp else 0.0)
            self.tebond.append(te_bond_ks if te_in_hp else 0.0)

            # SW bond lengths for this segment
            for kw in range(n_webs):
                if kw in seg_data['sw_begin']:
                    sw_b = seg_data['sw_begin'][kw]
                    sw_e = seg_data['sw_end'][kw]
                    C = keypoints.web_bonds[kw][:, sw_b:sw_e]
                else:
                    C = np.zeros((2, 1))
                self.swbonds[kw][:, ks] = m_to_mm * np.sum(C, 1)

            # Mold and prefab dimensions for this segment
            _calc_mold_segment(self, seg_data, ks, keypoints, ispan, key_labels)

        self.hp = pd.DataFrame(hp_rows)
        self.lp = pd.DataFrame(lp_rows)
        self.sw = pd.DataFrame(sw_rows)

        hp_dw = self.hp['dryweight'].sum() if not self.hp.empty else 0.0
        lp_dw = self.lp['dryweight'].sum() if not self.lp.empty else 0.0
        sw_dw = self.sw['dryweight'].sum() if not self.sw.empty else 0.0
        self.total_dryweight = hp_dw + lp_dw + sw_dw

        hp_tw = self.hp['curedweight'].sum() if not self.hp.empty else 0.0
        lp_tw = self.lp['curedweight'].sum() if not self.lp.empty else 0.0
        sw_tw = self.sw['curedweight'].sum() if not self.sw.empty else 0.0
        self.total_curedweight = hp_tw + lp_tw + sw_tw

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


def _find_layer_extents_segmented(layer_dist, layer_n, station_ids):
    """Find spanwise station extents within a segment where layer count >= layer_n.

    Returns begin_sta and end_sta as actual station indices from station_ids.
    """
    sta_logical = layer_dist >= layer_n
    begin_sta = []
    end_sta = []
    prev = 0
    for kk in range(len(sta_logical)):
        if sta_logical[kk] == 1 and prev == 0:
            begin_sta.append(station_ids[kk])
        if sta_logical[kk] == 0 and prev == 1:
            end_sta.append(station_ids[kk])
        elif kk == len(sta_logical) - 1 and prev == 1:
            end_sta.append(station_ids[kk])
        prev = sta_logical[kk]
    return begin_sta, end_sta


def _find_region_extents_segmented(key_labels, comp, seg_hp_extents, seg_lp_extents):
    """Find chordwise keypoint regions for a component, clipped to segment extents.

    Returns (hp_region, lp_region, sw_region) as sorted 2-element lists of
    0-indexed key_labels positions (or empty lists when not applicable).
    """
    le = key_labels.index('le')
    hp_region, lp_region, sw_region = [], [], []

    # Segment bounds on HP side (0-indexed positions in key_labels[:le+1])
    seg_hp_idxs = [key_labels[:le + 1].index(x) for x in seg_hp_extents]
    hp_max = max(seg_hp_idxs)  # closer to LE (higher index)
    hp_min = min(seg_hp_idxs)  # closer to TE (lower index)

    # Segment bounds on LP side (0-indexed positions in full key_labels)
    seg_lp_idxs = [key_labels[le:].index(x) + le for x in seg_lp_extents]
    lp_min = min(seg_lp_idxs)  # closer to LE (lower index)
    lp_max = max(seg_lp_idxs)  # closer to TE (higher index)

    if comp.hpextents and len(comp.hpextents) == 2:
        hp1, hp2 = [key_labels[:le + 1].index(x) for x in comp.hpextents]
        hp_start = min(max(hp1, hp2), hp_max)
        hp_end = max(min(hp1, hp2), hp_min)
        if max(hp1, hp2) > hp_min and min(hp1, hp2) < hp_max:
            hp_region = sorted([hp_start, hp_end])

    if comp.lpextents and len(comp.lpextents) == 2:
        lp1, lp2 = [key_labels[le:].index(x) + le for x in comp.lpextents]
        lp_start = max(min(lp1, lp2), lp_min)
        lp_end = min(max(lp1, lp2), lp_max)
        if min(lp1, lp2) < lp_max and max(lp1, lp2) > lp_min:
            lp_region = sorted([lp_start, lp_end])

    if (comp.hpextents and len(comp.hpextents) == 1
            and comp.lpextents and len(comp.lpextents) == 1):
        sw1 = key_labels[:le + 1].index(comp.hpextents[0])
        sw2 = key_labels[le:].index(comp.lpextents[0]) + le
        b_hp = key_labels[:le + 1].index('b')
        c_hp = key_labels[:le + 1].index('c')
        b_lp = key_labels[le:].index('b') + le
        c_lp = key_labels[le:].index('c') + le
        contains_hp_spar = hp_min <= min(b_hp, c_hp) and hp_max >= max(b_hp, c_hp)
        contains_lp_spar = lp_min <= min(b_lp, c_lp) and lp_max >= max(b_lp, c_lp)
        contains_hp_web = hp_min <= sw1 <= hp_max
        contains_lp_web = lp_min <= sw2 <= lp_max
        hp_spar_contains_web = min(b_hp, c_hp) <= sw1 <= max(b_hp, c_hp)
        lp_spar_contains_web = min(b_lp, c_lp) <= sw2 <= max(b_lp, c_lp)
        if contains_hp_web and contains_lp_web:
            if hp_spar_contains_web and lp_spar_contains_web:
                if contains_hp_spar and contains_lp_spar:
                    sw_region = [sw1, sw2]
            elif not hp_spar_contains_web and not lp_spar_contains_web:
                sw_region = [sw1, sw2]

    return hp_region, lp_region, sw_region


def _calc_component_stations(comp, seg_data, nlay, station_ids):
    """Track spanwise station bounds for mold-related component types."""
    valid = nlay > 0
    if not np.any(valid):
        return seg_data
    first_sta = station_ids[int(np.argmax(valid))]
    last_sta = station_ids[len(station_ids) - 1 - int(np.argmax(valid[::-1]))]
    name_lower = comp.name.lower()

    if 'spar' in name_lower:
        if (comp.hpextents and len(comp.hpextents) == 2
                and comp.lpextents and len(comp.lpextents) == 2):
            seg_data['hp_spar_begin'] = first_sta
            seg_data['hp_spar_end'] = last_sta
            seg_data['lp_spar_begin'] = first_sta
            seg_data['lp_spar_end'] = last_sta
        elif comp.lpextents and len(comp.lpextents) == 2:
            seg_data['lp_spar_begin'] = first_sta
            seg_data['lp_spar_end'] = last_sta
        elif comp.hpextents and len(comp.hpextents) == 2:
            seg_data['hp_spar_begin'] = first_sta
            seg_data['hp_spar_end'] = last_sta

    if 'reinf' in name_lower and 'le' in name_lower:
        seg_data['le_begin'] = min(seg_data.get('le_begin', first_sta), first_sta)
        seg_data['le_end'] = max(seg_data.get('le_end', last_sta), last_sta)

    if 'reinf' in name_lower and 'te' in name_lower:
        seg_data['te_begin'] = min(seg_data.get('te_begin', first_sta), first_sta)
        seg_data['te_end'] = max(seg_data.get('te_end', last_sta), last_sta)

    if 'root' in name_lower:
        seg_data['root_begin'] = min(seg_data.get('root_begin', first_sta), first_sta)
        seg_data['root_end'] = max(seg_data.get('root_end', last_sta), last_sta)

    return seg_data


def _calc_mold_segment(bom, seg_data, ks, keypoints, ispan, key_labels):
    """Compute mold and prefab dimensional characteristics for one blade segment."""
    le = key_labels.index('le')
    begin_sta = seg_data['begin_sta']
    end_sta = seg_data['end_sta']

    def _surf_kl(surf, labels):
        return sorted([i for i, k in enumerate(surf) if k.lower() in labels])

    def _add_mold_row(field, name, kl_idxs, sta_b, sta_e,
                      is_web=False, web_id=None):
        if sta_b >= sta_e:
            return
        if is_web:
            areas = keypoints.web_areas[web_id][sta_b:sta_e]
            arcs = keypoints.web_width[web_id][sta_b:sta_e + 1]
            lengths = keypoints.web_bonds[web_id][:, sta_b:sta_e]
        else:
            r0, r1 = min(kl_idxs), max(kl_idxs)
            areas = keypoints.key_areas[r0:r1, sta_b:sta_e]
            arcs = (keypoints.key_arcs[r1, sta_b:sta_e + 1]
                    - keypoints.key_arcs[r0, sta_b:sta_e + 1])
            lengths = keypoints.key_bonds[kl_idxs, sta_b:sta_e]
        row = {
            'Component': name,
            'Area': float(np.sum(areas)),
            'Length': float(np.sum(lengths) / 2),
            'Start_station': float(ispan[sta_b]),
            'End_station': float(ispan[sta_e]),
            'Max_width': float(np.amax(arcs)),
            'Avg_width': float(np.mean(arcs)),
            'Perimeter': float(np.sum(lengths) + arcs[0] + arcs[-1]),
        }
        if field not in bom.mold:
            bom.mold[field] = pd.DataFrame(columns=_MOLD_COLS)
        bom.mold[field] = pd.concat(
            [bom.mold[field], pd.DataFrame([row])], ignore_index=True
        )

    # HP surface molds
    if seg_data['hp_extents']:
        kl_surf = key_labels[:le + 1]
        spar_kl = _surf_kl(kl_surf, ['b', 'c'])
        le_reinf_kl = _surf_kl(kl_surf, ['le', 'a'])
        te_reinf_kl = _surf_kl(kl_surf, ['d', 'te'])
        le_web_kl = _surf_kl(kl_surf, ['a', 'b'])
        te_web_kl = _surf_kl(kl_surf, ['c', 'd'])
        root_kl = _surf_kl(kl_surf, ['le', 'te'])
        seg_kl = _surf_kl(kl_surf, [x.lower() for x in seg_data['hp_extents']])

        _add_mold_row('shell', f'HP{ks + 1}_shell_mold', seg_kl, begin_sta, end_sta)

        if min(seg_kl) <= min(spar_kl) and max(seg_kl) >= max(spar_kl):
            if 'hp_spar_begin' in seg_data:
                _add_mold_row('spar', f'HP{ks + 1}_spar_prefab', spar_kl,
                              seg_data['hp_spar_begin'], seg_data['hp_spar_end'])
            if begin_sta == 0 and 'root_begin' in seg_data:
                _add_mold_row('root', f'HP{ks + 1}_root_prefab', root_kl,
                              seg_data['root_begin'], seg_data['root_end'])

        if min(seg_kl) <= min(le_reinf_kl) and max(seg_kl) >= max(le_reinf_kl):
            if 'le_begin' in seg_data:
                _add_mold_row('LEreinf', f'HP{ks + 1}_LEreinf_preform', le_reinf_kl,
                              seg_data['le_begin'], seg_data['le_end'])

        if min(seg_kl) <= min(te_reinf_kl) and max(seg_kl) >= max(te_reinf_kl):
            if 'te_begin' in seg_data:
                _add_mold_row('TEreinf', f'HP{ks + 1}_TEreinf_prefab', te_reinf_kl,
                              seg_data['te_begin'], seg_data['te_end'])

        for kw, sw_b in seg_data['sw_begin'].items():
            sw_e = seg_data['sw_end'][kw]
            sw_label = seg_data['sw_labels'].get(kw, '')
            sw_hp_kl = kl_surf.index(seg_data['sw_extents'][kw][0])
            web_kl = le_web_kl if sw_hp_kl >= min(le_web_kl) else te_web_kl
            _add_mold_row('web', f'{sw_label}_mold', web_kl,
                          sw_b, sw_e, is_web=True, web_id=kw)

    # LP surface molds
    if seg_data['lp_extents']:
        kl_surf_lp = key_labels[le:]
        spar_kl_lp = [i + le for i in _surf_kl(kl_surf_lp, ['b', 'c'])]
        le_reinf_kl_lp = [i + le for i in _surf_kl(kl_surf_lp, ['le', 'a'])]
        te_reinf_kl_lp = [i + le for i in _surf_kl(kl_surf_lp, ['d', 'te'])]
        root_kl_lp = [i + le for i in _surf_kl(kl_surf_lp, ['le', 'te'])]
        seg_kl_lp = [
            i + le for i in _surf_kl(kl_surf_lp, [x.lower() for x in seg_data['lp_extents']])
        ]

        _add_mold_row('shell', f'LP{ks + 1}_shell_mold', seg_kl_lp, begin_sta, end_sta)

        if min(seg_kl_lp) <= min(spar_kl_lp) and max(seg_kl_lp) >= max(spar_kl_lp):
            if 'lp_spar_begin' in seg_data:
                _add_mold_row('spar', f'LP{ks + 1}_spar_prefab', spar_kl_lp,
                              seg_data['lp_spar_begin'], seg_data['lp_spar_end'])
            if begin_sta == 0 and 'root_begin' in seg_data:
                _add_mold_row('root', f'LP{ks + 1}_root_prefab', root_kl_lp,
                              seg_data['root_begin'], seg_data['root_end'])

        if min(seg_kl_lp) <= min(le_reinf_kl_lp) and max(seg_kl_lp) >= max(le_reinf_kl_lp):
            if 'le_begin' in seg_data:
                _add_mold_row('LEreinf', f'LP{ks + 1}_LEreinf_preform', le_reinf_kl_lp,
                              seg_data['le_begin'], seg_data['le_end'])

        if min(seg_kl_lp) <= min(te_reinf_kl_lp) and max(seg_kl_lp) >= max(te_reinf_kl_lp):
            if 'te_begin' in seg_data:
                _add_mold_row('TEreinf', f'LP{ks + 1}_TEreinf_prefab', te_reinf_kl_lp,
                              seg_data['te_begin'], seg_data['te_end'])


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
        self.angle: float = None
