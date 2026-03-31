import yaml
import numpy as np
import logging

from pynumad.utils.misc_utils import full_keys_from_substrings
from pynumad.io.airfoil_to_xml import airfoil_to_xml
from pynumad.io.yaml_to_blade import _add_materials
from pynumad.utils.interpolation import interpolator_wrap
from pynumad.objects.airfoil import Airfoil
from pynumad.objects.component import Component
from pynumad.objects.definition import Definition


def yaml_to_blade_v2(blade, filename: str, write_airfoils: bool = False):
    with open(filename) as f:
        data = yaml.load(f, Loader=yaml.FullLoader)

    definition = Definition()
    blade.definition = definition

    blade_data = data["components"]["blade"]
    blade_outer_shape = blade_data["outer_shape"]
    blade_reference_axis = blade_data["reference_axis"]
    blade_structure = blade_data["structure"]
    af_data = data["airfoils"]
    mat_data = data["materials"]
    
    ### Anchors ###
    # Resolve anchor references so all layers have inline {grid, values} dicts
    anchor_map = _build_anchor_map(blade_structure)
    _resolve_structure_anchors(blade_structure, anchor_map)

    ### Span ###
    L = blade_reference_axis["z"]["values"][-1]
    norm_span_grid = np.array(blade_reference_axis["z"]["grid"])
    full_span_grid = norm_span_grid * L
    
    ### Stations ###

    _add_stations_v2(
        definition, blade_outer_shape, blade_reference_axis, L, af_data, filename, write_airfoils
    )

    ### Spanwise props ###
    _add_spanwise_properties_v2(definition, blade_outer_shape, blade_reference_axis, L)
    blade.ispan = definition.ispan

    ### Materials ###
    _add_materials(definition, mat_data)

    ### Update blade structure ###
    # Ensure all components have full length grid/values
    blade_structure = update_internal_structure_v2(blade_structure, norm_span_grid)

    blade_structure_dict = {
        layer["name"].lower(): layer for layer in blade_structure["layers"]
    }
    
    ### Spar Caps ###

    _set_sparcap_arc_positions(definition, blade_structure_dict)
    
    ### TE/LE Bands ###
    definition.teband = _reinforcement_halfwidth_mm(["te", "reinf"], blade_structure_dict, anchor_map, full_span_grid)
    definition.leband = _reinforcement_halfwidth_mm(["le", "reinf"], blade_structure_dict, anchor_map, full_span_grid)
    
    ### Components ###
    definition.components = _build_component_dict(definition, blade_structure)
    _validate_components(definition.components)

    blade.update_blade()
    return blade


def _build_anchor_map(blade_structure):
    anchor_map = {}
    for anchor in blade_structure.get("anchors", []):
        anchor_map[anchor["name"]] = anchor
    for web in blade_structure.get("webs", []):
        for anchor in web.get("anchors", []):
            anchor_map[anchor["name"]] = anchor
    return anchor_map


def _resolve_nd_arc_field(field, anchor_map):
    if isinstance(field, dict) and "anchor" in field:
        anchor_name = field["anchor"]["name"]
        anchor_handle = field["anchor"]["handle"]
        anchor = anchor_map.get(anchor_name)
        if anchor is not None and anchor_handle in anchor:
            return anchor[anchor_handle]
        return {"grid": [0.0, 1.0], "values": [0.0, 0.0]}
    return field


def _resolve_structure_anchors(blade_structure, anchor_map):
    arc_keys = {"start_nd_arc", "end_nd_arc"}
    for layer in blade_structure.get("layers", []):
        for key in arc_keys:
            if key in layer:
                layer[key] = _resolve_nd_arc_field(layer[key], anchor_map)
    for web in blade_structure.get("webs", []):
        for key in arc_keys:
            if key in web:
                web[key] = _resolve_nd_arc_field(web[key], anchor_map)


def _add_stations_v2(
    definition, blade_outer_shape, blade_reference_axis, L, af_data, file, write_airfoils
):
    definition.span = np.multiply(np.transpose(blade_reference_axis["z"]["grid"]), L)
    definition.ispan = definition.span.copy()
    definition.stations = []

    af_dir_names = [af["name"] for af in af_data]
    tc_xL_list = []
    aero_cent = []

    for af_entry in blade_outer_shape["airfoils"]:
        af_name = af_entry["name"]
        tc_xL = af_entry["spanwise_position"]
        tc_xL_list.append(tc_xL)

        IAF = af_dir_names.index(af_name)
        aero_cent.append(af_data[IAF]["aerodynamic_center"])
        x = np.array(af_data[IAF]["coordinates"]["x"], dtype=float)
        y = np.array(af_data[IAF]["coordinates"]["y"], dtype=float)
        xf_coords = np.stack((x, y), 1)

        x, y = xf_coords[:, 0], xf_coords[:, 1]
        if np.sum(x[:-1] * y[1:] - x[1:] * y[:-1]) > 0:
            xf_coords = np.flipud(xf_coords)

        if write_airfoils:
            import os
            out_folder = "yaml2BladeDef_" + file.replace(".yaml", "")
            os.makedirs(out_folder + "/af_coords/", exist_ok=True)
            os.makedirs(out_folder + "/airfoil/", exist_ok=True)
            airfoil_to_xml(xf_coords, af_name, out_folder + "/af_coords/" + af_name + ".txt")

        af = Airfoil(coords=xf_coords, reference=af_name)
        af.resample(spacing="half-cosine")
        definition.add_station(af, tc_xL * L)

    definition.aerocenter = interpolator_wrap(np.multiply(tc_xL_list, L), aero_cent, definition.span)


def _add_spanwise_properties_v2(definition, blade_outer_shape, blade_reference_axis, L):
    span = definition.span

    definition.degreestwist = interpolator_wrap(
        np.multiply(np.transpose(blade_outer_shape["twist"]["grid"]), L),
        np.transpose(blade_outer_shape["twist"]["values"]),
        span,
    )
    definition.chord = interpolator_wrap(
        np.multiply(np.transpose(blade_outer_shape["chord"]["grid"]), L),
        np.transpose(blade_outer_shape["chord"]["values"]),
        span,
    )
    definition.percentthick = np.multiply(
        interpolator_wrap(
            np.multiply(blade_outer_shape["rthick"]["grid"], L),
            blade_outer_shape["rthick"]["values"],
            span,
        ),
        100,
    )
    # section_offset_y is absolute metres from LE; divide by chord for dimensionless fraction
    section_offset_y_at_span = interpolator_wrap(
        np.multiply(np.transpose(blade_outer_shape["section_offset_y"]["grid"]), L),
        np.transpose(blade_outer_shape["section_offset_y"]["values"]),
        span,
    )
    definition.chordoffset = section_offset_y_at_span / definition.chord

    definition.natural_offset = 0
    definition.prebend = interpolator_wrap(
        np.multiply(np.transpose(blade_reference_axis["x"]["grid"]), L),
        np.transpose(blade_reference_axis["x"]["values"]),
        span,
    )
    definition.sweep = interpolator_wrap(
        np.multiply(np.transpose(blade_reference_axis["y"]["grid"]), L),
        np.transpose(blade_reference_axis["y"]["values"]),
        span,
    )


def update_internal_structure_v2(blade_structure, norm_span_grid):
    """Interpolates all layer fields onto the full normalized span grid (0-1)."""
    nStations = len(norm_span_grid)
    keysToModify = {"thickness", "fiber_orientation", "start_nd_arc", "end_nd_arc"}

    for layer in blade_structure["layers"]:
        for currentKey in keysToModify.intersection(layer.keys()):
            field = layer[currentKey]
            if not isinstance(field, dict) or "grid" not in field:
                continue
            grid = field["grid"]
            values = field["values"]

            subIdx = np.where(
                (norm_span_grid >= grid[0]) & (norm_span_grid <= grid[-1])
            )[0]

            fullSpanValues = np.zeros(nStations)
            fullSpanValues[subIdx] = interpolator_wrap(grid, values, norm_span_grid[subIdx], "pchip")

            layer[currentKey]["grid"] = norm_span_grid
            layer[currentKey]["values"] = fullSpanValues

    return blade_structure


def _set_sparcap_arc_positions(definition, blade_structure_dict):
    sparCapKeys = full_keys_from_substrings(blade_structure_dict.keys(), ["spar"])
    if len(sparCapKeys) != 2:
        raise ValueError("Incorrect number of spar cap components")

    lpSideIndex = None
    hpSideIndex = None
    for i, name in enumerate(sparCapKeys):
        if "_ss" in name.lower():
            lpSideIndex = i
        elif "_ps" in name.lower():
            hpSideIndex = i

    if lpSideIndex is None or hpSideIndex is None:
        raise ValueError("Could not determine spar cap HP/LP sides from layer names (_SS/_PS suffix expected)")

    definition.sparcap_start_nd_arc_lp = blade_structure_dict[sparCapKeys[lpSideIndex]]["start_nd_arc"]["values"]
    definition.sparcap_end_nd_arc_lp   = blade_structure_dict[sparCapKeys[lpSideIndex]]["end_nd_arc"]["values"]
    definition.sparcap_start_nd_arc_hp = blade_structure_dict[sparCapKeys[hpSideIndex]]["start_nd_arc"]["values"]
    definition.sparcap_end_nd_arc_hp   = blade_structure_dict[sparCapKeys[hpSideIndex]]["end_nd_arc"]["values"]


def _reinforcement_halfwidth_mm(substrings, blade_structure_dict, anchor_map, full_span_grid):
    keys = full_keys_from_substrings(blade_structure_dict.keys(), substrings)
    L = full_span_grid[-1]
    total_width = np.zeros(len(full_span_grid))
    for key in keys:
        layer_name = blade_structure_dict[key]["name"]
        anchor = anchor_map.get(layer_name)
        if anchor and "width" in anchor:
            w_grid = np.multiply(anchor["width"]["grid"], L)
            total_width += interpolator_wrap(w_grid, anchor["width"]["values"], full_span_grid)
        else:
            logging.debug(f"No width found for band layer {layer_name}")
    return total_width * 1000 / 2



def _build_component_dict(definition, blade_structure):
    component_dict = {}
    for layer in blade_structure["layers"]:
        comp = Component()
        name = layer["name"].lower()
        comp.group = 0
        comp.name = layer["name"]
        comp.materialid = layer["material"]
        try:
            comp.fabricangle = np.mean(layer["fiber_orientation"]["values"])
        except (KeyError, TypeError):
            comp.fabricangle = 0
        comp.imethod = "pchip" if "spar" in name else "linear"
        grid = np.transpose(layer["thickness"]["grid"])
        n_plies = np.round(
            np.multiply(np.transpose(layer["thickness"]["values"]), 1000.0)
            / definition.materials[comp.materialid].layerthickness
        )
        comp.control_points = np.stack((grid, n_plies), axis=1)
        comp.pinnedends = 0

        _assign_extents_to_comp(comp, name)

        component_dict[comp.name] = comp
    return component_dict


def _assign_extents_to_comp(comp, name):
    if "spar" in name and "ps" in name:
        comp.hpextents = ["b", "c"]
    elif "spar" in name and "ss" in name:
        comp.lpextents = ["b", "c"]
    elif "shell" in name or "uv" in name or "gel" in name:
        comp.hpextents = ["le", "te"]
        comp.lpextents = ["le", "te"]
    elif "te_" in name and "reinf" in name:
        comp.hpextents = ["d", "te"] if "ps" in name else None
        comp.lpextents = ["d", "te"] if "ss" in name else None
        if "ps" not in name and "ss" not in name:
            comp.hpextents = ["d", "te"]
            comp.lpextents = ["d", "te"]
    elif "le_" in name and "reinf" in name:
        comp.hpextents = ["le", "a"] if "ps" in name else None
        comp.lpextents = ["le", "a"] if "ss" in name else None
        if "ps" not in name and "ss" not in name:
            comp.hpextents = ["le", "a"]
            comp.lpextents = ["le", "a"]
    elif "te_" in name and "ss" in name and "filler" in name:
        comp.lpextents = ["c", "d"]
    elif "le_" in name and "ss" in name and "filler" in name:
        comp.lpextents = ["a", "b"]
    elif "le_" in name and "ps" in name and "filler" in name:
        comp.hpextents = ["a", "b"]
    elif "te_" in name and "ps" in name and "filler" in name:
        comp.hpextents = ["c", "d"]
    elif "web" in name:
        if "fore" in name or (name.count("web") == 1 and "1" in name):
            comp.hpextents = ["b"]
            comp.lpextents = ["b"]
            comp.group = 1
        elif "aft" in name or "rear" in name or (name.count("web") == 1 and "0" in name):
            comp.hpextents = ["c"]
            comp.lpextents = ["c"]
            comp.group = 2


def _validate_components(component_dict):
    names = list(component_dict.keys())

    spar_ps = full_keys_from_substrings(names, ["spar", "ps"])
    spar_ss = full_keys_from_substrings(names, ["spar", "ss"])
    if len(spar_ps) != 1 or len(spar_ss) != 1:
        raise ValueError(f"Expected 1 PS and 1 SS spar cap, found {len(spar_ps)} PS and {len(spar_ss)} SS")

    coating = full_keys_from_substrings(names, ["uv"]) or full_keys_from_substrings(names, ["gel"])
    if len(coating) != 1:
        raise ValueError(f"Expected 1 UV/gelcoat component, found {len(coating)}")

    shell = full_keys_from_substrings(names, ["shell"])
    if len(shell) != 2:
        raise ValueError(f"Expected 2 shell components, found {len(shell)}")

    te_reinf = full_keys_from_substrings(names, ["te_", "reinf"])
    if len(te_reinf) not in (1, 2):
        raise ValueError(f"Expected 1 or 2 TE reinforcement layers, found {len(te_reinf)}")

    le_reinf = full_keys_from_substrings(names, ["le_", "reinf"])
    if len(le_reinf) not in (1, 2):
        raise ValueError(f"Expected 1 or 2 LE reinforcement layers, found {len(le_reinf)}")

    fore_web = (
        full_keys_from_substrings(names, ["web", "fore"])
        or full_keys_from_substrings(names, ["web", "1"])
    )
    if len(fore_web) == 0:
        raise ValueError("No fore web layers found")

    aft_web = (
        full_keys_from_substrings(names, ["web", "aft"])
        or full_keys_from_substrings(names, ["web", "0"])
        or full_keys_from_substrings(names, ["web", "rear"])
    )
    if len(aft_web) == 0:
        raise ValueError("No aft web layers found")
