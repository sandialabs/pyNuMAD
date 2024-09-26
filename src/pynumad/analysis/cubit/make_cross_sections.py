import warnings

try:
    from cubit import *
    from PyCubed_Main import *
except ModuleNotFoundError:
    warnings.warn("Cubit not configured, so cubit functionality will not be available.")

import numpy as np
import os

def debug():
    cubit.cmd(f"delete curve 1")
    cubit.cmd(f'save as "Debug.cub" overwrite')
def addColor(blade, volume_or_surface):
    # Adds color to volume or surfaces by material

    materials = blade.definition.materials
    color_dict = {}
    color_dict["adhesive"] = "yellow"
    color_dict["carbon"] = "grey"
    color_dict["uni"] = "seagreen"
    color_dict["triax"] = "lightgreen"
    color_dict["biax"] = "greenyellow"
    color_dict["foam"] = "khaki"

    for mat_name in materials:
        for color in color_dict.keys():
            if color in mat_name.lower():
                parse_string = f'with name "*{mat_name}*"'

                vol_ids = parse_cubit_list(volume_or_surface, parse_string)
                cubit.cmd(
                    f"color {volume_or_surface} {l2s(vol_ids)}  mesh {color_dict[color]}"
                )
                cubit.cmd(
                    f"color {volume_or_surface} {l2s(vol_ids)}  geometry {color_dict[color]}"
                )

                break


def surface_from_two_curves(top_curve, bottom_curve):
    v2Left, v2Right = selCurveVerts(top_curve)
    v1Left, v1Right = selCurveVerts(bottom_curve)
    cubit.cmd(f"create curve vertex {v1Left} {v2Left}")
    cubit.cmd(f"create curve vertex {v1Right} {v2Right}")
    cubit.cmd(
        f'create surface curve {l2s([get_last_id("curve")-1, bottom_curve,top_curve,get_last_id("curve")])}'
    )


def get_cs_normal_vector(xyz):
    npts, _ = xyz.shape

    # Create Referece line as a spline
    vertex_list = []
    for kcp in range(npts):
        vertex_list.append(create_vertex(xyz[kcp, 0], xyz[kcp, 1], xyz[kcp, 2]))

    c1 = create_curve(vertex_list[0], vertex_list[1])
    c2 = create_curve(vertex_list[0], vertex_list[2])

    midPoint = list(c1.position_from_fraction(0.5))
    tangent1 = c1.tangent(midPoint)
    midPoint = list(c2.position_from_fraction(0.5))
    tangent2 = c2.tangent(midPoint)
    cs_normal = vectNorm(crossProd(list(tangent1), list(tangent2)))
    cubit.cmd(f'delete curve {get_last_id("curve")}')
    cubit.cmd(f'delete curve {get_last_id("curve")-1}')
    return cs_normal


def get_blade_geometry_for_station(blade, i_station):
    geometry = blade.geometry
    return np.array(
        [
            geometry.coordinates[:, 0, i_station],
            geometry.coordinates[:, 1, i_station],
            geometry.coordinates[:, 2, i_station],
        ]
    ).transpose()
    # xcoords=blade.profiles[:,0,i_station]*blade.ichord[i_station]
    # ycoords=blade.profiles[:,1,i_station]*blade.ichord[i_station]
    # zcoord=blade.ispan[i_station]
    # zcoords=[zcoord]*len(xcoords)
    # return np.array([xcoords,ycoords,zcoords]).transpose()


def write_spline_from_coordinate_points(cubit, xyz):
    # xyz is npts by 3 array holding the coordinates of the points
    npts, _ = xyz.shape
    n_start = get_last_id("vertex") + 1
    for i_point in range(npts):
        coords = xyz[i_point, :]
        create_vertex(coords[0], coords[1], coords[2])
    n_end = get_last_id("vertex")
    vertex_list = list(range(n_start, n_end + 1))
    cubit.cmd(f"create curve spline vertex {l2s(vertex_list)}")


def extend_curve_at_vertex_to_length(curve_to_extend_id, extension_length, curve_start_or_end):
    # Extend all offset curves
    if curve_start_or_end.lower() == "start":
        temp_vert_id, _ = selCurveVerts(curve_to_extend_id)
        vectorDirectionSign = -1
    else:
        _, temp_vert_id = selCurveVerts(curve_to_extend_id)
        vectorDirectionSign = 1

    tangentLocationCoords = cubit.vertex(temp_vert_id).coordinates()

    x = cubit.curve(curve_to_extend_id).tangent(tangentLocationCoords)[0]
    y = cubit.curve(curve_to_extend_id).tangent(tangentLocationCoords)[1]
    z = cubit.curve(curve_to_extend_id).tangent(tangentLocationCoords)[2]
    tangent_direction = (
        vectorDirectionSign * extension_length * np.array(vectNorm([x, y, z]))
    )  # Unit vector of tangent *Scaled by offset curve length
    newVertexCoords = np.array(tangentLocationCoords) + tangent_direction

    v1 = cubit.create_vertex(newVertexCoords[0], newVertexCoords[1], newVertexCoords[2])

    # if statement needed to maintain original curve sense
    if curve_start_or_end.lower() == "start":
        cubit.create_curve(v1, cubit.vertex(temp_vert_id))
        first_curve_id=get_last_id("curve")
        second_curve_id=curve_to_extend_id
        keep_curve=2

    else:
        cubit.create_curve(cubit.vertex(temp_vert_id), v1)
        first_curve_id=curve_to_extend_id
        second_curve_id=get_last_id("curve")
        keep_curve=2


    cubit.cmd(f"split curve {first_curve_id} fraction 0.02 from end ")
    first_curve_id = get_last_id("curve") - 1
    cubit.cmd(f"split curve {second_curve_id} distance 0.02 from start ")
    second_curve_id = get_last_id("curve") 

    return splice_two_curves(first_curve_id, second_curve_id, keep_curve)


def removeBadTEgeometry(blade, i_station, curve_id, flatback_curve_id):
    definition = blade.definition
    cubit.cmd(
        f"split curve {curve_id} distance {definition.chord[i_station]*0.002} from start "
    )
    cubit.cmd(f'delete curve {get_last_id("curve")-1}')
    curve_id = get_last_id("curve")

    curve_start_or_end = "start"
    extension_length = 1 * cubit.curve(curve_id).length()
    curve_id = extend_curve_at_vertex_to_length(curve_id, extension_length, curve_start_or_end)

    _, v1 = selCurveVerts(curve_id)
    cubit.cmd(
        f"trim curve {curve_id} atintersection curve {flatback_curve_id} keepside vertex {v1}"
    )
    return get_last_id("curve")


def get_curve_offset_direction(curve_id, lp_hp_side, cs_normal):
    # This function is used to determine which way a curve offset will go. This is needed, for example,
    # to make sure the outer mold line curve is being offset towrads the interior of the blade.
    tol = 0.01
    cubit.cmd(f"create curve offset curve {curve_id} distance 0.0005")
    temp_id = get_last_id("curve")

    n_start = get_last_id("surface")   
    cubit.cmd(f"create surface skin curve {curve_id} {temp_id}")
    n_end = get_last_id("surface")
    
    if n_end - n_start==0:
      cubit.cmd(f'save as "Debug.cub" overwrite')
      raise Exception(
            f"Offset direction for curve {curve_id} was not able to be found due to failed skin curve operation.")
    
    cubit.cmd(f"delete curve {temp_id}")
    n = get_surface_normal(get_last_id("surface"))
    if (
        abs(n[0] - cs_normal[0]) < tol
        and abs(n[1] - cs_normal[1]) < tol
        and abs(n[2] - cs_normal[2]) < tol
    ):
        offset_sign = 1
    else:
        offset_sign = -1
    cubit.cmd(f'delete body {get_last_id("body")}')

    # Need to flip direction since LP curve is written to run clockwise
    if lp_hp_side.lower() == "lp":
        offset_sign = -1 * offset_sign
    return offset_sign


def offset_curve_and_combine_fragments_if_needed(curve_id, offset_distance):
    # Sometimes when offseting a curve, the offset is broken into multiple curves.
    # This happens when the curvature is to high for the offset distance. If this happens
    # this fuction combines the fragmented curves into one spline.
    n_start = get_last_id("curve")
    cubit.cmd(f"create curve offset curve {curve_id} distance {offset_distance} extended")
    n_end = get_last_id("curve")


    # cubit.cmd(f'save as "Debug.cub" overwrite')

    if n_end - n_start > 1:
        curveList = []
        curveList += list(range(n_start + 1, n_end + 1))
        vertex_list = []
        v1, _ = selCurveVerts(curveList[0])
        vertex_list.append(v1)


        for curve in curveList:   
            if cubit.curve(curve).length() > 1e-4: #Toss really small curves
                n_start = get_last_id("vertex") + 1
                cubit.cmd(f"create vertex on curve {curve} segment 200")
                n_end = get_last_id("vertex")
                vertex_list += list(range(n_start, n_end + 1))

        _, v1 = selCurveVerts(curveList[-1])
        cubit.cmd(f"create curve spline vertex {l2s(vertex_list)}")
        cubit.cmd(f"delete vertex {l2s(vertex_list[1:-1])}")
        cubit.cmd(f"delete curve {l2s(curveList)}")


def streamline_curve_intersections(first_curve_id, second_curve_id, keep_curve):
    # keep_curve is either 1 or 2
    n_start = get_last_id("vertex")
    cubit.cmd(f"create vertex atintersection curve {first_curve_id} {second_curve_id}")
    n_end = get_last_id("vertex")
    vertex_list = list(range(n_start + 1, n_end + 1))
    if len(vertex_list) > 0:
        cubit.cmd(f"split curve {first_curve_id} at vertex {vertex_list[-1]}")
        first_curve_id = get_last_id("curve") - 1
        cubit.cmd(f"split curve {second_curve_id} at vertex {vertex_list[-1]}")
        second_curve_id = get_last_id("curve")

        cubit.cmd(f"split curve {first_curve_id} distance 0.005 from end ")
        first_curve_id = get_last_id("curve") - 1
        cubit.cmd(f"split curve {second_curve_id} distance 0.005 from start ")
        second_curve_id = get_last_id("curve") 
        temp_id = splice_two_curves(first_curve_id, second_curve_id, keep_curve)
        return temp_id, vertex_list[-1]
    else:
        if keep_curve == 1:
            return first_curve_id, None
        else:
            return second_curve_id, None


def splice_two_curves(first_curve_id, second_curve_id, keep_curve):
    # Given two curves splice them into one curve. This should even work for curves that make corners.
    # Curve sence matters. The first curve's sense (tangent) should point towards the second curve.

    vertex_list = []
    v1, _ = selCurveVerts(first_curve_id)
    vertex_list.append(v1)
    n_start = get_last_id("vertex") + 1
    cubit.cmd(f"create vertex on curve {first_curve_id} segment 200")
    n_end = get_last_id("vertex")
    vertex_list += list(range(n_start, n_end + 1))
    _, v1 = selCurveVerts(first_curve_id)
    vertex_list.append(v1)
    v2, _ = selCurveVerts(second_curve_id)
    vertex_list.append(v2)
    n_start = get_last_id("vertex") + 1
    cubit.cmd(f"create vertex on curve {second_curve_id} segment 200")
    n_end = get_last_id("vertex")
    vertex_list += list(range(n_start, n_end + 1))
    _, v2 = selCurveVerts(second_curve_id)
    vertex_list.append(v2)
    cubit.cmd(f"create curve spline vertex {l2s(vertex_list)}")
    cubit.cmd(f"delete curve {first_curve_id} {second_curve_id}")
    second_curve_id = get_last_id("curve")
    cubit.cmd(f"delete vertex {l2s(vertex_list[1:-1])}")
    if keep_curve == 1:
        return first_curve_id
    elif keep_curve == 2:
        return second_curve_id


def extend_curve_past_curve_and_trim(
    curve_to_extend_id, curve_start_or_end, curve_idThatCutsExtendedCurve
):
    # Given two curves that are not necessarily intersecting extend {curve_to_extend_id} then trim it at
    # {curve_idThatCutsExtendedCurve}. {curve_start_or_end} defines which side of the curve to extend so
    # you need to know the curve sense

    extension_length = 2 * cubit.curve(curve_to_extend_id).length()
    curve_to_extend_id = extend_curve_at_vertex_to_length(
        curve_to_extend_id, extension_length, curve_start_or_end
    )

    n_start = get_last_id("vertex")
    cubit.cmd(
        f"create vertex AtIntersection curve {curve_idThatCutsExtendedCurve}  {curve_to_extend_id}"
    )
    splitVertexID = get_last_id("vertex")
    if n_start == splitVertexID:
        print(
            f"curve_to_extend_id {curve_to_extend_id} curve_idThatCutsExtendedCurve {curve_idThatCutsExtendedCurve}"
        )

        cubit.cmd(f'save as "Debug.cub" overwrite')
        raise Exception(
            f"Curve {curve_to_extend_id} was not able to be extended to curve {curve_idThatCutsExtendedCurve} because their intersection was not found."
        )

    cubit.cmd(f"split curve {curve_to_extend_id} at vertex {splitVertexID}")

    if curve_start_or_end.lower() == "start":
        curve_to_extend_id = get_last_id("curve")
        cubit.cmd(f'delete curve {get_last_id("curve")-1}')
    else:
        curve_to_extend_id = get_last_id("curve") - 1
        cubit.cmd(f'delete curve {get_last_id("curve")}')
    return curve_to_extend_id


# def rename_last_surface(part_name, i_station, i_modeled_layers, material_name, part_name_id):
#     # Every cross sectional surface that is created must be followed by a call to this function
#     part_name_id += 1


#     surface_name = (
#         part_name
#         + "Station"
#         + str(i_station).zfill(3)
#         + "_layer"
#         + str(i_modeled_layers)
#         + "_"
#         + material_name
#         + "_surface"
#         + str(part_name_id)
#     )
#     cubit.cmd(f'surface {get_last_id("surface")} rename "{surface_name}"')
#     return part_name_id


def add_surface_dict_entry(
    surface_dict, surface_object, my_curve_order, my_vert_order, material_name, ply_angle
):
    surface_dict[surface_object.id()] = {}

    # Curves:
    idList = []
    for curveObject in surface_object.curves():
        idList.append(curveObject.id())
    idList = [idList[i] for i in my_curve_order]

    surface_dict[surface_object.id()]["curves"] = idList

    # Verts:
    idList = []
    for vertObject in surface_object.vertices():
        idList.append(vertObject.id())
    idList = [idList[i] for i in my_vert_order]
    surface_dict[surface_object.id()]["verts"] = idList
    surface_dict[surface_object.id()]["material_name"] = material_name
    surface_dict[surface_object.id()]["ply_angle"] = ply_angle


def make_cross_section_surface(lp_hp_side,
    surface_dict,
    i_station,
    part_name,
    top_curve,
    bottom_curve,
    material_name,
    ply_angle,
    part_name_id,
    i_modeled_layers,
    materials_used,
    stack_id
):
    # Given two curves in a cross section, create a surface by connecting the end points then
    # rename the surface and add to the surface dictionary
    surface_from_two_curves(top_curve, bottom_curve)
    materials_used.add(material_name)

    curve_name = cubit.get_entity_name("curve", bottom_curve)

    if 'web_thickness' in curve_name:
        append_str='_web_thickness'
    else:
        append_str=''

    # Rename last surface
    part_name_id += 1
    if 'web_thickness' in curve_name:
        part_name+='_web_thickness'
        surface_name = ( 
            part_name+ lp_hp_side+"Station"+ str(i_station).zfill(3)+ "_layer"+ str(i_modeled_layers)+ "_"+ material_name+ "_surface"+ str(part_name_id))
    else:
        if stack_id>-1:
            surface_name = (
                part_name+ lp_hp_side+"Station"+ str(i_station).zfill(3)+ '_stack'+str(stack_id).zfill(3)+"_layer"+ str(i_modeled_layers)+ "_"+ 
                material_name+ "_surface"+ str(part_name_id))
        else:
            surface_name = (
                part_name+ lp_hp_side+"Station"+ str(i_station).zfill(3)+ "_layer"+ str(i_modeled_layers)+ "_"+ material_name+ "_surface"+ str(part_name_id))

    # surface_name = (
    #     part_name+ "Station"+ str(i_station).zfill(3)+ "_layer"+ str(i_modeled_layers)+ "_"+ material_name+ "_surface"+ str(part_name_id))
    cubit.cmd(f'surface {get_last_id("surface")} rename "{surface_name}"')

    # part_name_id = rename_last_surface(
    #     part_name+append_str, i_station, i_modeled_layers, material_name, part_name_id)
    add_surface_dict_entry(
        surface_dict,
        cubit.surface(get_last_id("surface")),
        [0, 1, 2, 3],
        [0, 1, 2, 3],
        material_name,
        ply_angle,
    )
    
    if i_modeled_layers ==1:
        layer_name_base = 'core_thickness'
    else:
        layer_name_base = 'face_thickness'


    cubit.cmd(f'curve {surface_dict[get_last_id("surface")]["curves"][1]} rename "{layer_name_base}{str(i_station).zfill(3)}_right"')
    cubit.cmd(f'curve {surface_dict[get_last_id("surface")]["curves"][-1]} rename "{layer_name_base}{str(i_station).zfill(3)}_left"')

    cubit.cmd(f'curve {surface_dict[get_last_id("surface")]["verts"][-1]} rename "layer_thickness"')
    
   # if i_modeled_layers == 0:
    curve_id = surface_dict[get_last_id("surface")]["curves"][0]
    curve_name = cubit.get_entity_name("curve", curve_id)

    if 'web_thickness' in curve_name:
        if 'face'in curve_name:
            append_str='_face_web_thickness'
        elif 'core' in curve_name:
            append_str='_core_web_thickness'

    else:
        if stack_id>-1:
            append_str = '_stack'+str(stack_id).zfill(3)
        else:
            append_str =''

    if i_modeled_layers==0:
        cubit.cmd(f'curve {surface_dict[get_last_id("surface")]["curves"][0]} rename "hoop_direction{str(i_station).zfill(3)+append_str}_oml"')
    else:
        cubit.cmd(f'curve {surface_dict[get_last_id("surface")]["curves"][0]} rename "hoop_direction{str(i_station).zfill(3)+append_str}"')
    cubit.cmd(f'curve {surface_dict[get_last_id("surface")]["curves"][2]} rename "hoop_direction{str(i_station).zfill(3)+append_str}"')

    return part_name_id, materials_used


def write_LE_adhesive_curves(
    HPstackThickness,
    LPstackThickness,
    adhesiveThickness,
    hp_key_curve,
    lp_key_curve,
    cs_normal,
):
    def unite_LE_geom_with_adhesive_curves(lp_hp_side, keyCurve):
        if lp_hp_side.lower() == "hp":
            # keyCurve=hp_key_curve
            adhesiveCurveID = le_adhesive_curve_ids[0][0]
            offset_curve = hpOffset
        else:
            # keyCurve=lp_key_curve
            adhesiveCurveID = le_adhesive_curve_ids[1][0]
            offset_curve = lpOffset


        v1, _ = selCurveVerts(keyCurve)

        cubit.cmd(
            f"trim curve {keyCurve} atintersection curve {adhesiveCurveID} keepside vertex {v1}"
        )
        keyCurve = get_last_id("curve")

        v1, _ = selCurveVerts(offset_curve)
        cubit.cmd(
            f"trim curve {offset_curve} atintersection curve {adhesiveCurveID} keepside vertex {v1}"
        )

        offset_curve = get_last_id("curve")
        _, v1 = selCurveVerts(keyCurve)
        _, v2 = selCurveVerts(offset_curve)

        cubit.cmd(f"create curve spline vertex {v1} {v2}")
        cubit.cmd(f"delete curve {adhesiveCurveID}")
        cubit.cmd(f"delete curve {offset_curve}")

        return get_last_id("curve"), keyCurve

    le_adhesive_curve_ids = [[], []]

    _, p1 = selCurveVerts(hp_key_curve)
    _, p2 = selCurveVerts(lp_key_curve)

    coords = []
    coords.append(list(cubit.vertex(p1).coordinates()))
    coords.append(list(cubit.vertex(p2).coordinates()))
    coords = np.array(coords)
    coords = np.mean(coords, 0)
    v1 = cubit.create_vertex(coords[0], coords[1], coords[2])

    cubit.cmd(f'vertex {get_last_id("vertex")} copy')
    midPointOML = get_last_id("vertex")

    # Offset OML curves to final layer offset
    cubit.cmd(f"curve {hp_key_curve} copy")
    cubit.cmd(f'split curve {get_last_id("curve")} fraction 0.05 from end')
    offset_curve_and_combine_fragments_if_needed(get_last_id("curve"), -1 * HPstackThickness)

    _, p1 = selCurveVerts(get_last_id("curve"))
    hpOffset = get_last_id("curve")

    cubit.cmd(f"curve {lp_key_curve} copy")
    cubit.cmd(f'split curve {get_last_id("curve")} fraction 0.05 from end')
    offset_curve_and_combine_fragments_if_needed(get_last_id("curve"), -1 * LPstackThickness)

    _, p2 = selCurveVerts(get_last_id("curve"))
    lpOffset = get_last_id("curve")

    coords = []
    coords.append(list(cubit.vertex(p1).coordinates()))
    coords.append(list(cubit.vertex(p2).coordinates()))
    coords = np.array(coords)
    coords = np.mean(coords, 0)
    v2 = cubit.create_vertex(coords[0], coords[1], coords[2])
    c1 = cubit.create_curve(v1, v2)
    adhesive_mid_line = get_last_id("curve")

    # Extend midline on both sides to make sure other curves eventually intersect with it
    curve_start_or_end = "start"
    extension_length = 2 * cubit.curve(adhesive_mid_line).length()
    adhesive_mid_line = extend_curve_at_vertex_to_length(
        adhesive_mid_line, extension_length, curve_start_or_end
    )

    curve_start_or_end = "end"
    extension_length = 1 * cubit.curve(adhesive_mid_line).length()
    adhesive_mid_line = extend_curve_at_vertex_to_length(
        adhesive_mid_line, extension_length, curve_start_or_end
    )

    # Copy and move since offset does not seem to work with strait lines

    # get offset vector

    axial_direction = cs_normal
    position = cubit.curve(adhesive_mid_line).position_from_fraction(1.0)
    tangent_direction = cubit.curve(adhesive_mid_line).tangent(position)

    normal_direction = crossProd(axial_direction, tangent_direction)
    normal_direction = (
        adhesiveThickness
        / 2
        * np.array(
            vectNorm([normal_direction[0], normal_direction[1], normal_direction[2]])
        )
    )
    cubit.cmd(
        f"curve {adhesive_mid_line} copy move x {normal_direction[0]} y {normal_direction[1]} z {normal_direction[2]} nomesh"
    )

    le_adhesive_curve_ids[0].append(get_last_id("curve"))
    normal_direction = -1 * normal_direction
    cubit.cmd(
        f"curve {adhesive_mid_line} copy move x {normal_direction[0]} y {normal_direction[1]} z {normal_direction[2]} nomesh"
    )
    le_adhesive_curve_ids[1].append(get_last_id("curve"))
    cubit.cmd(f"delete curve {adhesive_mid_line}")

    key_curves = [hp_key_curve, lp_key_curve]
    for iSide, lp_hp_side in enumerate(["HP", "LP"]):
        ###HP###
        le_adhesive_curve_ids[iSide][0], key_curves[iSide] = unite_LE_geom_with_adhesive_curves(
            lp_hp_side, key_curves[iSide]
        )

        # Make Copies
        cubit.cmd(f"curve {le_adhesive_curve_ids[iSide][0]} copy")
        le_adhesive_curve_ids[iSide].append(get_last_id("curve"))

        cubit.cmd(f"curve {le_adhesive_curve_ids[iSide][1]} copy")
        le_adhesive_curve_ids[iSide].append(get_last_id("curve"))
        # Extend
        curve_start_or_end = "end"

        extension_length = 1 * cubit.curve(le_adhesive_curve_ids[iSide][1]).length()
        le_adhesive_curve_ids[iSide][1] = extend_curve_at_vertex_to_length(
            le_adhesive_curve_ids[iSide][1], extension_length, curve_start_or_end
        )

    return key_curves[0], key_curves[1], le_adhesive_curve_ids


def split_curve_at_coordinte_points(coordinates_to_split_curve, curve_idToSplit):
    cubit.cmd(f"curve {curve_idToSplit} copy")
    tempCurveID = get_last_id("curve")

    nDPs, _ = coordinates_to_split_curve.shape
    id_start = get_last_id("vertex") + 1
    for kcp in range(nDPs):
        temp = cubit.curve(tempCurveID).closest_point(
            [
                coordinates_to_split_curve[kcp, 0],
                coordinates_to_split_curve[kcp, 1],
                coordinates_to_split_curve[kcp, 2],
            ]
        )
        create_vertex(temp[0], temp[1], temp[2])
    id_end = get_last_id("vertex")

    DPverticies = [i for i in range(id_start, id_end + 1)]

    id_start = get_last_id("curve") + 1
    cubit.cmd(f"split curve {tempCurveID} at vertex {l2s(DPverticies)}")
    id_end = get_last_id("curve")
    key_curves = [i for i in range(id_start, id_end + 1)]
    return key_curves


def split_key_curves(key_curves, aft_web_stack, fore_web_stack, web_adhesive_width):
    ###Do not split TE reinf
    temp_base_curve_ids = [key_curves[0]]

    ###split TE panel in half
    cubit.cmd(f"split curve {key_curves[1]} fraction 0.5")

    temp_base_curve_ids.append(get_last_id("curve") - 1)
    temp_base_curve_ids.append(get_last_id("curve"))

    ###Partition sparcap curve
    vertex_list = []
    web_layer_thickness = 0
    n_start = get_last_id("vertex") + 1
    for iLayer in reversed(range(len(aft_web_stack.plygroups))):
        web_layer_thickness += (
            aft_web_stack.plygroups[iLayer].thickness
            * aft_web_stack.plygroups[iLayer].nPlies
            / 1000
        )
        cubit.cmd(
            f"create vertex on curve {key_curves[2]}  distance {web_layer_thickness} from start"
        )
    cubit.cmd(
        f"create vertex on curve {key_curves[2]}  distance {web_layer_thickness+web_adhesive_width} from start"
    )

    # get total foreweb thickness
    web_layer_thickness = sum(fore_web_stack.layer_thicknesses()) / 1000
    cubit.cmd(
        f"create vertex on curve {key_curves[2]}  distance {web_layer_thickness+web_adhesive_width} from end"
    )
    for iLayer in reversed(range(len(fore_web_stack.plygroups))):
        cubit.cmd(
            f"create vertex on curve {key_curves[2]}  distance {web_layer_thickness} from end"
        )
        web_layer_thickness -= (
            fore_web_stack.plygroups[iLayer].thickness
            * fore_web_stack.plygroups[iLayer].nPlies
            / 1000
        )

    n_end = get_last_id("vertex")
    vertex_list += list(range(n_start, n_end + 1))

    n_start = get_last_id("curve") + 1
    cubit.cmd(f"split curve {key_curves[2]} at vertex {l2s(vertex_list)}")
    n_end = get_last_id("curve")

    cubit.cmd(f'curve {n_start} {n_start+2} rename "face_web_thickness"')
    cubit.cmd(f'curve {n_end-2} {n_end} rename "face_web_thickness"') 

    cubit.cmd(f'curve {n_start+1} rename "core_web_thickness"')
    cubit.cmd(f'curve {n_end-1} rename "core_web_thickness"') 

    temp_base_curve_ids.append(n_start)
    temp_base_curve_ids.append(n_end)
    spar_cap_base_curves = list(range(n_start + 1, n_end))

    ###split LE panel in half
    cubit.cmd(f"split curve {key_curves[3]} fraction 0.5")
    temp_base_curve_ids.append(get_last_id("curve") - 1)
    temp_base_curve_ids.append(get_last_id("curve"))

    ###Do not split LE reinf
    temp_base_curve_ids.append(key_curves[-1])

    return temp_base_curve_ids, spar_cap_base_curves


def get_mid_line(blade, iLE, i_station, geometry_scaling):
    geometry = blade.geometry
    X = geometry.coordinates[:, 0, i_station] * geometry_scaling
    Y = geometry.coordinates[:, 1, i_station] * geometry_scaling
    Z = geometry.coordinates[:, 2, i_station] * geometry_scaling

    ###### Get averge line
    xHP = X[1:iLE]
    xLP = np.flip(X[iLE - 1 : -1])
    yHP = Y[1:iLE]
    yLP = np.flip(Y[iLE - 1 : -1])
    zHP = Z[1:iLE]
    zLP = np.flip(Z[iLE - 1 : -1])
    midline = np.zeros((len(xHP), 3))
    for i_point in range(len(xHP)):
        midline[i_point, 0] = (xHP[i_point] + xLP[i_point]) / 2
        midline[i_point, 1] = (yHP[i_point] + yLP[i_point]) / 2
        midline[i_point, 2] = (zHP[i_point] + zLP[i_point]) / 2

    return midline


def get_adjustment_curve(curve_ids, layer_offset_dist, curve_start_or_end, end_layer_taper_curve):
    n_start = get_last_id("vertex") + 1
    curve_fraction = 1.0 / 3
    for i_curve, curve_id in enumerate(curve_ids):
        curve_length = cubit.curve(curve_id).length()
        if end_layer_taper_curve is not None and i_curve < end_layer_taper_curve - 1:
            if curve_length * curve_fraction < layer_offset_dist:
                cubit.cmd(
                    f"create vertex on curve {curve_id} fraction {curve_fraction} from {curve_start_or_end}"
                )
            else:
                cubit.cmd(
                    f"create vertex on curve {curve_id} distance {layer_offset_dist} from {curve_start_or_end}"
                )

        else:
            cubit.cmd(
                f"create vertex on curve {curve_id} distance {layer_offset_dist} from {curve_start_or_end}"
            )

    n_end = get_last_id("vertex")
    vertex_list = list(range(n_start, n_end + 1))
    cubit.cmd(f"create curve spline vertex {l2s(vertex_list)}")
    adjustment_curve = get_last_id("curve")
    cubit.cmd(f"delete vertex {l2s(vertex_list[1:-1])}")
    return adjustment_curve


def make_cs_perimeter_layer_areas(wt_name,
    surface_dict,
    i_station,
    station_stacks,
    cs_params,
    thickness_scaling,
    lp_hp_side,
    last_round_station,
    last_flat_station,
    part_name_id,
    n_modeled_layers,
    cs_normal,
    lp_hp_dict,
    materials_used,
    stack_ct
):
    # part_name = lp_hp_side + "shell"
    part_name = "shell"
    
    if i_station > last_round_station:
        is_flatback = True
    else:
        is_flatback = False

    # Assumes that #HP side is made first
    if lp_hp_side.lower() == "hp":
        lp_hp_side_index = 0
        camber_offsetSign = 1
    else:
        lp_hp_side_index = 1
        camber_offsetSign = -1

    offset_sign_camberID = get_curve_offset_direction(
        lp_hp_dict["camberID"], "LP", cs_normal
    )

    base_curve_index_ct = 0
    station_stacks.shape
    nStationLayups = len(station_stacks)
    station_stacks.shape

    last_perimeter = nStationLayups - 2

    for i_perimeter in range(nStationLayups - 1):  # Skip the last stack since the current and the next stack are generated at the same time.
        with open(f"{wt_name}.log", "a") as logFile:
            logFile.write(f"\tlp_hp_side {lp_hp_side}, i_perimeter={i_perimeter}\n")

        current_stack = station_stacks[i_perimeter]
        next_stack = station_stacks[i_perimeter + 1]

        current_stack_layer_thicknesses = np.array(current_stack.layer_thicknesses()) / 1000
        next_stack_layer_thicknesses = np.array(next_stack.layer_thicknesses()) / 1000

        cubit.cmd(
            f'curve {lp_hp_dict["base_curve_ids"][lp_hp_side_index][base_curve_index_ct]} copy'
        )
        current_base_curve_id = get_last_id("curve")
        base_curve_index_ct += 1
        cubit.cmd(
            f'curve {lp_hp_dict["base_curve_ids"][lp_hp_side_index][base_curve_index_ct]} copy'
        )
        next_base_curve_id = get_last_id("curve")
        base_curve_index_ct += 1

        current_stack_surface_list = []
        transition_stack_surface_list = []
        next_stack_surface_list = []

        current_stack_layer_offset = 0
        next_stack_layer_offset = 0
        layer_thickness_transition_lengths = []
        # Get all offsets and layer_thickness_transition_lengths
        thinest_layer_thickness_current_stack = 1e22  # initialize to a large value
        thinest_layer_thickness_next_stack = 1e22

        minimum_layer_transition_length=1.1*(sum(current_stack_layer_thicknesses)+sum(next_stack_layer_thicknesses))/2
        for i_modeled_layers in range(n_modeled_layers):
            current_stack_layer_offset += current_stack_layer_thicknesses[i_modeled_layers]
            next_stack_layer_offset += next_stack_layer_thicknesses[i_modeled_layers]

            adjacent_layer_missmatch = abs(current_stack_layer_offset - next_stack_layer_offset)

            if adjacent_layer_missmatch > minimum_layer_transition_length:
                layer_thickness_transition_lengths.append(
                    adjacent_layer_missmatch
                    / tan(math.radians(cs_params["layer_transition_angle"]))
                )
            else:
                layer_thickness_transition_lengths.append(minimum_layer_transition_length)

            # Also find the thinest layer in stack for meshing purposes
            if (
                current_stack_layer_thicknesses[i_modeled_layers] > 0
                and current_stack_layer_thicknesses[i_modeled_layers]
                < thinest_layer_thickness_current_stack
            ):
                thinest_layer_thickness_current_stack = current_stack_layer_thicknesses[
                    i_modeled_layers
                ]

            if (
                next_stack_layer_thicknesses[i_modeled_layers] > 0
                and next_stack_layer_thicknesses[i_modeled_layers]
                < thinest_layer_thickness_next_stack
            ):
                thinest_layer_thickness_next_stack = next_stack_layer_thicknesses[
                    i_modeled_layers
                ]

        max_layer_thickness_transition_length = max(layer_thickness_transition_lengths)

        if i_perimeter in [0, 2]:
            left_bottom_curve = current_base_curve_id
            cubit.cmd(
                f"split curve {next_base_curve_id} distance {max_layer_thickness_transition_length} from start "
            )
            transition_bottom_curve = get_last_id("curve") - 1
            right_bottom_curve = get_last_id("curve")
            transition_stack = next_stack
        elif i_perimeter in [1, 3]:
            right_bottom_curve = next_base_curve_id
            cubit.cmd(
                f"split curve {current_base_curve_id} distance {max_layer_thickness_transition_length} from end "
            )
            left_bottom_curve = get_last_id("curve") - 1
            transition_bottom_curve = get_last_id("curve")
            transition_stack = current_stack
        else:
            raise ValueError(f"i_perimeter {i_perimeter} not recognized")

        bottom_left_vertex_curve_left, bottom_right_vertex_curve_left = selCurveVerts(
            left_bottom_curve
        )
        bottom_left_vertex_curve_right, bottom_right_vertex_curve_right = selCurveVerts(
            right_bottom_curve
        )

        # This if statement prepares all layer curves such that they taper at the TE
        if i_perimeter == 0:
            current_stack_right_curves = []
            current_stack_left_curves = []

            # Base curve copy
            cubit.cmd(f"curve {left_bottom_curve} copy")
            base_curve_id_copy = get_last_id("curve")
            offset_sign_base_curve_id_copy = get_curve_offset_direction(
                base_curve_id_copy, lp_hp_side, cs_normal
            )

            # offset camber to make gap
            cubit.cmd(
                f'create curve offset curve {lp_hp_dict["camberID"]} distance {camber_offsetSign*offset_sign_camberID*cs_params["te_adhesive_thickness"][i_station]/2} extended'
            )
            camber_offset = get_last_id("curve")

            # Top Bounding Curve
            offset_distance = (
                1 * offset_sign_base_curve_id_copy * sum(current_stack_layer_thicknesses)
            )
            offset_curve_and_combine_fragments_if_needed(base_curve_id_copy, offset_distance)
            top_bounding_curve = get_last_id("curve")

            if is_flatback:
                curve_start_or_end = "start"
                extension_length = 1 * cubit.curve(top_bounding_curve).length()
                top_bounding_curve = extend_curve_at_vertex_to_length(
                    top_bounding_curve, extension_length, curve_start_or_end
                )
                keep_curve = 2

                (
                    top_bounding_curve,
                    begin_layer_taper_vertex_id,
                ) = streamline_curve_intersections(
                    camber_offset, top_bounding_curve, keep_curve
                )

            else:
                lp_hp_dict["flatback_curve_id"] = camber_offset

                curve_start_or_end = "start"
                extension_length = 1 * cubit.curve(base_curve_id_copy).length()
                base_curve_id_copy = extend_curve_at_vertex_to_length(
                    base_curve_id_copy, extension_length, curve_start_or_end
                )

                curve_start_or_end = "start"
                extension_length = 1 * cubit.curve(top_bounding_curve).length()
                top_bounding_curve = extend_curve_at_vertex_to_length(
                    top_bounding_curve, extension_length, curve_start_or_end
                )

                _, v1 = selCurveVerts(base_curve_id_copy)
                cubit.cmd(
                    f'trim curve {base_curve_id_copy} atintersection curve {lp_hp_dict["flatback_curve_id"]} keepside vertex {v1}'
                )
                base_curve_id_copy = get_last_id("curve")

            # Trim curve at TE.adhesive
            _, v1 = selCurveVerts(top_bounding_curve)
            cubit.cmd(
                f'trim curve {top_bounding_curve} atintersection curve {lp_hp_dict["flatback_curve_id"]} keepside vertex {v1}'
            )
            top_bounding_curve = get_last_id("curve")
            offset_sign_top_bounding_curve = get_curve_offset_direction(
                top_bounding_curve, lp_hp_side, cs_normal
            )

            if is_flatback:  # and begin_layer_taper_vertex_id is not None:
                # Make list of curves that will be used to taper each layer
                npts = 30
                n_start = get_last_id("curve") + 1

                if lp_hp_side.lower() == "hp":
                    sign_correction = 1
                else:
                    sign_correction = -1

                for i_point in range(npts):
                    if i_point == 0:
                        cubit.cmd(f"create vertex on curve {base_curve_id_copy} start")
                        cubit.cmd(f"create vertex on curve {top_bounding_curve} start")
                    else:
                        cubit.cmd(
                            f"create vertex on curve {base_curve_id_copy} fraction {(i_point)/(npts-1)} from start"
                        )
                        cubit.cmd(
                            f"create vertex on curve {top_bounding_curve} fraction {(i_point)/(npts-1)} from start"
                        )
                    cubit.cmd(
                        f'create curve vertex {get_last_id("vertex")-1} {get_last_id("vertex")}'
                    )

                n_end = get_last_id("curve")
                curve_ids = list(range(n_start, n_end + 1))

                # If layers are tapeded toward TE find the index in the list of curves (curve_ids) marks the end of the tapering
                if begin_layer_taper_vertex_id is not None:
                    found_flag = False
                    for curve_id in curve_ids:
                        temp = cubit.curve(curve_id).curve_center()
                        tangent_direction = cubit.curve(get_last_id("curve")).tangent(
                            temp
                        )

                        moment_arm = np.array(temp) - np.array(
                            cubit.vertex(begin_layer_taper_vertex_id).coordinates()
                        )

                        normal_direction = sign_correction * np.array(
                            crossProd(tangent_direction, moment_arm)
                        )

                        if (
                            not found_flag
                            and dotProd(normal_direction, cs_normal) < 0
                        ):
                            found_flag = True
                            end_layer_taper_curve = i_point
                else:
                    end_layer_taper_curve = None

            # TE Adhesive curve
            if not is_flatback:

                v1, _ = selCurveVerts(base_curve_id_copy)
                v2, _ = selCurveVerts(top_bounding_curve)
                cubit.cmd(
                    f'split curve {lp_hp_dict["flatback_curve_id"]} at vertex {v1} {v2}'
                )
                cubit.cmd(
                    f'delete curve {get_last_id("curve")} {get_last_id("curve")-2}'
                )
                lp_hp_dict["flatback_curve_id"] = get_last_id("curve") - 1
                cubit.cmd(f'curve {lp_hp_dict["flatback_curve_id"]} copy')

            iLayer = 0
            offset_distance = (
                1 * offset_sign_base_curve_id_copy * current_stack_layer_thicknesses[iLayer]
            )
            offset_curve_and_combine_fragments_if_needed(base_curve_id_copy, offset_distance)
            first_layer_offset = get_last_id("curve")

            curve_start_or_end = "start"
            first_layer_offset = extend_curve_past_curve_and_trim(
                first_layer_offset, curve_start_or_end, lp_hp_dict["flatback_curve_id"]
            )

            # Only do the following if all layer thicknesses are unequal

            if (is_flatback and abs(min(current_stack.layer_thicknesses())- max(current_stack.layer_thicknesses()))> 0.0001):
                layer_offset_dist = current_stack_layer_thicknesses[0]
                curve_start_or_end = "start"
                first_layer_offset = get_adjustment_curve(
                    curve_ids, layer_offset_dist, curve_start_or_end, end_layer_taper_curve
                )

            cubit.cmd(
                f'create curve offset curve {lp_hp_dict["camberID"]} distance {camber_offsetSign*offset_sign_camberID*cs_params["te_adhesive_thickness"][i_station]/2} extended'
            )
            camber_offset = get_last_id("curve")

            offset_sign_top_bounding_curve = get_curve_offset_direction(
                top_bounding_curve, lp_hp_side, cs_normal
            )

            offset_distance = (
                -1 * offset_sign_top_bounding_curve * current_stack_layer_thicknesses[-1]
            )
            offset_curve_and_combine_fragments_if_needed(top_bounding_curve, offset_distance)
            last_layer_offset = get_last_id("curve")

            curve_start_or_end = "start"
            last_layer_offset = extend_curve_past_curve_and_trim(
                last_layer_offset, curve_start_or_end, lp_hp_dict["flatback_curve_id"]
            )

            if (
                is_flatback
                and abs(
                    min(current_stack.layer_thicknesses())
                    - max(current_stack.layer_thicknesses())
                )
                > 0.0001
            ):
                layer_offset_dist = current_stack_layer_thicknesses[-1]
                curve_start_or_end = "end"
                last_layer_offset = get_adjustment_curve(
                    curve_ids, layer_offset_dist, curve_start_or_end, end_layer_taper_curve
                )


            # cubit.cmd(f'split curve {first_layer_offset} at vertex {lp_hp_dict["perimeter_verts_for_te_adhesive"][lp_hp_side]}')
            # current_stack_left_curves.append(get_last_id("curve") - 1)
            # current_stack_right_curves.append(get_last_id("curve"))
            # curve_len = cubit.curve(current_stack_left_curves[0]).length()

            # cubit.cmd(f'split curve {base_curve_id_copy} at vertex {lp_hp_dict["perimeter_verts_for_te_adhesive"][lp_hp_side]}')
            # current_stack_left_curves.insert(0, get_last_id("curve") - 1)
            # current_stack_right_curves.insert(0, get_last_id("curve"))

            # cubit.cmd(f'split curve {last_layer_offset} at vertex {lp_hp_dict["perimeter_verts_for_te_adhesive"][lp_hp_side]}')
            # current_stack_left_curves.append(get_last_id("curve") - 1)
            # current_stack_right_curves.append(get_last_id("curve"))
            
            # cubit.cmd(f'split curve {top_bounding_curve} at vertex {lp_hp_dict["perimeter_verts_for_te_adhesive"][lp_hp_side]}')
            # current_stack_left_curves.append(get_last_id("curve") - 1)
            # current_stack_right_curves.append(get_last_id("curve"))

            current_stack_right_curves=[base_curve_id_copy,first_layer_offset,last_layer_offset,top_bounding_curve]
            #Further subdived current_stack_left_curves for roundTEadhesive volumes

            
            current_stack_left_curves_splited=[]
            for vertex_id in lp_hp_dict["perimeter_verts_for_te_adhesive"][lp_hp_side]:
                temp_list_left=[]
                temp_list_right=[]
                for curve_id in current_stack_right_curves:
                    n_start = get_last_id("curve")
                    cubit.cmd(f'split curve {curve_id} at vertex {vertex_id}')
                    n_end = get_last_id("curve")
                    if n_end -  n_start< 2:
                        print('')
                    temp_list_left.append(get_last_id("curve") - 1)
                    temp_list_right.append(get_last_id("curve"))
                current_stack_right_curves=temp_list_right
                current_stack_left_curves_splited.append(temp_list_left)
            #current_stack_left_curves_splited.append(temp_list_right)
            
            #if is_flatback:
            if i_station in list(range(last_round_station+1,last_flat_station+1)):
                for i_split in range(len(current_stack_left_curves_splited)-1):
                    lp_hp_dict["round_te_adhesive_curve_list"][lp_hp_side_index].append(current_stack_left_curves_splited[i_split][-1])
            
            if i_station >= last_flat_station:
                lp_hp_dict["flat_te_adhesive_curve_list"][lp_hp_side_index].append(current_stack_left_curves_splited[-1][-1])
            
            #### Next Stack (the panel might intersect the camberline so the following is needed
            next_stack_curves = []
            cubit.cmd(f"curve {right_bottom_curve} copy")
            base_curve_id_copy = get_last_id("curve")
            offset_sign_base_curve_id_copy = get_curve_offset_direction(
                base_curve_id_copy, lp_hp_side, cs_normal
            )
            next_stack_curves.append(base_curve_id_copy)

            ### Offset camber to make gap
            # Offset is increased to create a larger clearance between HP LP shells so that the panels
            # to not self intersect during a simulation (this may not be needed)
            #cubit.cmd(f'create curve offset curve {lp_hp_dict["camberID"]} distance {camber_offsetSign*offset_sign_camberID*0.001*4} extended')
            cubit.cmd(f'create curve offset curve {lp_hp_dict["camberID"]} distance {camber_offsetSign*offset_sign_camberID*cs_params["te_adhesive_thickness"][i_station]/2} extended')

            camber_offset = get_last_id("curve")

            iLayer = 0
            offset_distance = (
                1 * offset_sign_base_curve_id_copy * next_stack_layer_thicknesses[0]
            )
            offset_curve_and_combine_fragments_if_needed(base_curve_id_copy, offset_distance)
            next_stack_curves.append(get_last_id("curve"))

            offset_distance = (
                1 * offset_sign_base_curve_id_copy * sum(next_stack_layer_thicknesses)
            )
            offset_curve_and_combine_fragments_if_needed(base_curve_id_copy, offset_distance)
            top_bounding_curve = get_last_id("curve")

            keep_curve = 2
            top_bounding_curve, intersectionVertex = streamline_curve_intersections(
                camber_offset, top_bounding_curve, keep_curve
            )

            v1, _ = selCurveVerts(base_curve_id_copy)
            cubit.cmd(f"split curve {top_bounding_curve} at vertex {v1}")
            top_bounding_curve = get_last_id("curve")
            cubit.cmd(f'delete curve {get_last_id("curve")-1}')
            next_stack_curves.append(top_bounding_curve)

            offset_sign_top_bounding_curve = get_curve_offset_direction(
                top_bounding_curve, lp_hp_side, cs_normal
            )

            offset_distance = (
                -1 * offset_sign_top_bounding_curve * next_stack_layer_thicknesses[-1]
            )
            offset_curve_and_combine_fragments_if_needed(top_bounding_curve, offset_distance)
            last_layer_offset = get_last_id("curve")
            next_stack_curves.insert(-1, last_layer_offset)

        for i_modeled_layers in range(n_modeled_layers):
            current_stackOffset = current_stack_layer_thicknesses[i_modeled_layers]
            next_stack_offset = next_stack_layer_thicknesses[i_modeled_layers]

            offset_sign_left_bottom_curve = get_curve_offset_direction(
                left_bottom_curve, lp_hp_side, cs_normal
            )

            offset_sign_right_bottom_curve = get_curve_offset_direction(
                right_bottom_curve, lp_hp_side, cs_normal
            )

            # Special Treatment for TE
            if i_perimeter == 0:
                # Create Left Areas Only

                material_name = current_stack.plygroups[i_modeled_layers].materialid
                ply_angle = current_stack.plygroups[i_modeled_layers].angle

                for i_split in range(len(current_stack_left_curves_splited)):
                    part_name_id, materials_used = make_cross_section_surface(lp_hp_side,
                        surface_dict,
                        i_station,
                        part_name,
                        current_stack_left_curves_splited[i_split][i_modeled_layers + 1],
                        current_stack_left_curves_splited[i_split][i_modeled_layers],
                        material_name,
                        ply_angle,
                        part_name_id,
                        i_modeled_layers,
                        materials_used,
                        stack_ct
                    )
                if i_station <= last_flat_station:
                #if not is_flatback:
                    lp_hp_dict["round_te_adhesive_curve_list"][lp_hp_side_index].append(
                        surface_dict[get_last_id("surface")]["curves"][-1]
                    )
                else:
                    surf_id=get_last_id("surface")
                    curve_id = cubit.surface(surf_id).curves()[3].id()
                    curve_name_split = cubit.get_entity_name("curve", curve_id).split('_')
                    curve_name=curve_name_split[0]+'_'+curve_name_split[1]+'_oml'
                    cubit.cmd(f'curve {curve_id} rename "{curve_name}"')
                if i_station == last_round_station - 1:
                    v1 = surface_dict[get_last_id("surface")]["verts"][0]
                    cubit.cmd(f'vertex {v1} rename "linear"')
                    v1 = surface_dict[get_last_id("surface")]["verts"][-1]
                    cubit.cmd(f'vertex {v1} rename "linear"')

                left_bottom_curve = current_stack_right_curves[i_modeled_layers]
                leftTopCurve = current_stack_right_curves[i_modeled_layers + 1]
                [bottom_left_vertex_curve_left, bottom_right_vertex_curve_left] = selCurveVerts(
                    current_stack_right_curves[i_modeled_layers]
                )
                [top_left_vertex_curve_left, top_right_vertex_curve_left] = selCurveVerts(
                    current_stack_right_curves[i_modeled_layers + 1]
                )

                # Create Left Areas Only
                right_bottom_curve = next_stack_curves[i_modeled_layers]
                right_top_curve = next_stack_curves[i_modeled_layers + 1]
                [
                    bottom_left_vertex_curve_right,
                    bottom_right_vertex_curve_right,
                ] = selCurveVerts(next_stack_curves[i_modeled_layers])
                [top_left_vertex_curve_right, top_right_vertex_curve_right] = selCurveVerts(
                    next_stack_curves[i_modeled_layers + 1]
                )
            else:


                cubit.cmd(
                    f"create curve offset curve {left_bottom_curve} distance {offset_sign_left_bottom_curve*current_stackOffset} extended"
                )
                leftTopCurve = get_last_id("curve")
                [top_left_vertex_curve_left, top_right_vertex_curve_left] = selCurveVerts(
                    get_last_id("curve")
                )
                
                offset_curve_and_combine_fragments_if_needed(
                    right_bottom_curve, offset_sign_right_bottom_curve * next_stack_offset
                )
                last_offset_curve = get_last_id("curve")
                right_top_curve = get_last_id("curve")
                [top_left_vertex_curve_right, top_right_vertex_curve_right] = selCurveVerts(
                    get_last_id("curve")
                )

            if i_perimeter == last_perimeter:
                curve_start_or_end = "end"
                
                last_offset_curve = extend_curve_past_curve_and_trim(
                    last_offset_curve,
                    curve_start_or_end,
                    lp_hp_dict["le_adhesive_curve_ids"][lp_hp_side_index][1],
                )
                right_top_curve = last_offset_curve
                [
                    bottom_left_vertex_curve_right,
                    bottom_right_vertex_curve_right,
                ] = selCurveVerts(right_bottom_curve)
                [top_left_vertex_curve_right, top_right_vertex_curve_right] = selCurveVerts(
                    right_top_curve
                )

            cubit.cmd(
                f"create curve vertex {top_right_vertex_curve_left} {top_left_vertex_curve_right}"
            )
            transition_top_curve = get_last_id("curve")
            cubit.cmd(
                f"create curve vertex {bottom_left_vertex_curve_left} {top_left_vertex_curve_left}"
            )
            cubit.cmd(
                f"create curve vertex {bottom_right_vertex_curve_left} {top_right_vertex_curve_left}"
            )
            n_curves_final = get_last_id("curve")
            left_side_curves = [i for i in range(n_curves_final - 1, n_curves_final + 1)]
            cubit.cmd(
                f"create curve vertex {bottom_left_vertex_curve_right} {top_left_vertex_curve_right}"
            )
            cubit.cmd(
                f"create curve vertex {bottom_right_vertex_curve_right} {top_right_vertex_curve_right}"
            )
            n_curves_final = get_last_id("curve")
            right_side_curves = [i for i in range(n_curves_final - 1, n_curves_final + 1)]

            if i_perimeter == last_perimeter:
                if is_flatback:
                    cubit.cmd(
                        f'split curve {lp_hp_dict["le_adhesive_curve_ids"][lp_hp_side_index][1]} at vertex {top_right_vertex_curve_right}'
                    )
                    right_side_curves[-1] = get_last_id("curve") - 1
                    lp_hp_dict["le_adhesive_curve_ids"][lp_hp_side_index][1] = get_last_id(
                        "curve"
                    )
                lp_hp_dict["le_adhesive_curve_list"][lp_hp_side_index].append(
                    right_side_curves[-1]
                )

            ### Create Surfaces ###
            # Sufaces for current_stack
            material_name = current_stack.plygroups[i_modeled_layers].materialid
            ply_angle = current_stack.plygroups[i_modeled_layers].angle
            part_name_id, materials_used = make_cross_section_surface(lp_hp_side,
                surface_dict,
                i_station,
                part_name,
                leftTopCurve,
                left_bottom_curve,
                material_name,
                ply_angle,
                part_name_id,
                i_modeled_layers,
                materials_used,
                stack_ct+1
            )
            current_stack_surface_list.append(get_last_id("surface"))

            # Surfaces for transition_stack
            material_name = transition_stack.plygroups[i_modeled_layers].materialid
            ply_angle = transition_stack.plygroups[i_modeled_layers].angle
            part_name_id, materials_used = make_cross_section_surface(lp_hp_side,
                surface_dict,
                i_station,
                part_name,
                transition_top_curve,
                transition_bottom_curve,
                material_name,
                ply_angle,
                part_name_id,
                i_modeled_layers,
                materials_used,
                stack_ct+2
            )
            transition_stack_surface_list.append(get_last_id("surface"))

            # Surfaces for next_stack
            material_name = next_stack.plygroups[i_modeled_layers].materialid
            ply_angle = next_stack.plygroups[i_modeled_layers].angle
            part_name_id, materials_used = make_cross_section_surface(lp_hp_side,
                surface_dict,
                i_station,
                part_name,
                right_top_curve,
                right_bottom_curve,
                material_name,
                ply_angle,
                part_name_id,
                i_modeled_layers,
                materials_used,
                stack_ct+3
            )
            next_stack_surface_list.append(get_last_id("surface"))

            ### Reset ###
            # Reset curves
            left_bottom_curve = leftTopCurve
            transition_bottom_curve = transition_top_curve
            right_bottom_curve = right_top_curve

            # Reset vertices
            bottom_left_vertex_curve_left = top_left_vertex_curve_left
            bottom_right_vertex_curve_left = top_right_vertex_curve_left
            bottom_left_vertex_curve_right = top_left_vertex_curve_right
            bottom_right_vertex_curve_right = top_right_vertex_curve_right

        # Build spar caps
        if i_perimeter == 1:
            lp_hp_dict["web_interface_curves"][lp_hp_side_index] = [right_top_curve]
            temp_ct=0
            for ic, current_curveID in enumerate(
                lp_hp_dict["spar_cap_base_curves"][lp_hp_side_index]
            ):
                bottom_curve = current_curveID
                offSetSign = get_curve_offset_direction(
                    bottom_curve, lp_hp_side, cs_normal
                )
                temp_ct+=1
                for it, thickness in enumerate(next_stack_layer_thicknesses):
                    cubit.cmd(
                        f"create curve offset curve {bottom_curve} distance {offSetSign*thickness} extended"
                    )
                    top_curve = get_last_id("curve")

                    material_name = next_stack.plygroups[it].materialid
                    ply_angle = next_stack.plygroups[it].angle
                    part_name_id, materials_used = make_cross_section_surface(lp_hp_side,
                        surface_dict,
                        i_station,
                        part_name,
                        top_curve,
                        bottom_curve,
                        material_name,
                        ply_angle,
                        part_name_id,
                        it,
                        materials_used,
                        stack_ct+temp_ct
                    )
                    next_stack_surface_list.append(get_last_id("surface"))

                    if it == 2 and ic != 3:
                        lp_hp_dict["web_interface_curves"][lp_hp_side_index].append(
                            top_curve
                        )
                    bottom_curve = top_curve
                    
            stack_ct+=1
        elif i_perimeter == 2:
            lp_hp_dict["web_interface_curves"][lp_hp_side_index].append(leftTopCurve)

        stack_ct+=3
    return part_name_id, lp_hp_dict, stack_ct


####################################################
####################################################
####################################################
####################################################
####################################################


def create_simplist_surface_for_TE_or_LE_adhesive(
    i_station,
    surface_dict,
    part_name,
    adhesive_curve_list,
    adhesiveMatID,
    part_name_id,
    n_modeled_layers,
    materials_used,
    stack_id
):

    for i_curve in range(len(adhesive_curve_list[0])):
        ply_angle = (
            0  # Ply angle is always zero since adhesive is always assumed as isotropic
        )
        if 'flat' in part_name:
            c_top  = adhesive_curve_list[1][i_curve]
            c_bot = adhesive_curve_list[0][i_curve]
        else:
            #Make input curves tangent to oml for material orientation puposes
            v_outer_top, v_inner_top = selCurveVerts(adhesive_curve_list[1][i_curve])
            v_outer_bot, v_inner_bot = selCurveVerts(adhesive_curve_list[0][i_curve])

            cubit.cmd(f'create curve vertex {v_inner_top} {v_inner_bot}')
            c_top = get_last_id("curve")

            cubit.cmd(f'create curve vertex {v_outer_top} {v_outer_bot}')
            c_bot = get_last_id("curve")   

        lp_hp_side=''
        part_name_id, materials_used = make_cross_section_surface(lp_hp_side,
            surface_dict,
            i_station,
            part_name,
            c_top,
            c_bot,
            adhesiveMatID,
            ply_angle,
            part_name_id,
            n_modeled_layers + 1,
            materials_used,
            stack_id
        )

    return part_name_id


def print_sine_curve_between_two_verts(vBot, vTop, amplitude, direction):
    n_sine_curve_sample_points = 7
    cubit.cmd(f"create curve vertex {vBot} {vTop}")

    idCurve = get_last_id("curve")
    vertex_list=[vBot,vTop]
    if round(amplitude, 3) > 0:
        n_start = get_last_id("vertex") + 1
        cubit.cmd(
            f'create vertex on curve {get_last_id("curve")} segment {n_sine_curve_sample_points-1}'
        )
        n_end = get_last_id("vertex")
        sine_curve_sample_points = np.linspace(0, pi, n_sine_curve_sample_points)
        vertex_offsets = amplitude * np.sin(sine_curve_sample_points)
        vertex_list = list(range(n_start, n_end + 1))
        for iVert, vertexOffset in enumerate(
            vertex_offsets[1:-1]
        ):  # skip first and last point since those are considered fixed and the offset is zero anyway
            cubit.cmd(f"move vertex {vertex_list[iVert]} {direction} {vertexOffset}")
        vertex_list.insert(0, vBot)
        vertex_list.append(vTop)
        cubit.cmd(f"create curve spline vertex {l2s(vertex_list)}")
        cubit.cmd(f"delete curve {idCurve}")
    return get_last_id("curve"),vertex_list


def make_cs_web_layer_areas(
    surface_dict,
    i_station,
    aft_web_stack,
    fore_web_stack,
    web_interface_curves,
    cs_params,
    part_name_id,
    cs_normal,
    n_modeled_layers,
    materials_used,
):
    
    aft_web_overwrap_thickness = (
        aft_web_stack.layer_thicknesses()[0] + aft_web_stack.layer_thicknesses()[-1]
    ) / 1000
    fore_web_overwrap_thickness = (
        fore_web_stack.layer_thicknesses()[0] + fore_web_stack.layer_thicknesses()[-1]
    ) / 1000
    part_name = "web"
    ### First create the first two layers. The first layer is the adhesive. The second layer is the web overwrap layer
    for i_curveList, curveList in enumerate(web_interface_curves):
        n_base_curves_web = len(curveList)
        if i_curveList == 0:
            lp_hp_side = "HP"
        else:
            lp_hp_side = "LP"
        for i_curve, bottom_curve in enumerate(curveList):
            offSetSign = get_curve_offset_direction(
                bottom_curve, lp_hp_side, cs_normal
            )

            if i_curve < n_base_curves_web / 2:
                layer_thicknesses = [
                    cs_params["web_aft_adhesive_thickness"][i_station],
                    aft_web_overwrap_thickness,
                ]
            else:
                layer_thicknesses = [
                    cs_params["web_fore_adhesive_thickness"][i_station],
                    fore_web_overwrap_thickness,
                ]

            for it, thickness in enumerate(layer_thicknesses):
                cubit.cmd(
                    f"create curve offset curve {bottom_curve} distance {offSetSign*thickness} extended"
                )
                top_curve = get_last_id("curve")

                if it == 0:
                    material_name = cs_params["adhesive_mat_name"]
                    ply_angle = 0

                else:
                    if i_curve < n_base_curves_web / 2:
                        material_name = aft_web_stack.plygroups[0].materialid
                        ply_angle = aft_web_stack.plygroups[0].angle
                    else:
                        material_name = fore_web_stack.plygroups[0].materialid
                        ply_angle = fore_web_stack.plygroups[0].angle
                
                part_name_id, materials_used = make_cross_section_surface(lp_hp_side,
                    surface_dict,
                    i_station,
                    part_name,
                    top_curve,
                    bottom_curve,
                    material_name,
                    ply_angle,
                    part_name_id,
                    n_modeled_layers + it,
                    materials_used,
                    -1
                )

                bottom_curve = top_curve

            # update web interface curves for vertical web retions
            web_interface_curves[i_curveList][i_curve] = top_curve

    ### Create vertical web regions
    lp_hp_side=''
    # remove curves that are not going to be part of the vertical web
    for i_curveList, curveList in enumerate(web_interface_curves):
        curveList.pop(3)
        curveList.pop(3)

    n_base_curves_web = len(web_interface_curves[0])
    for i_curve in range(n_base_curves_web):
        vHP, _ = selCurveVerts(web_interface_curves[0][i_curve])
        vLP, _ = selCurveVerts(web_interface_curves[1][i_curve])
        top_curve,_ = print_sine_curve_between_two_verts(
            vHP, vLP, cs_params["max_web_imperfection_distance"][i_station], "x"
        )
        _, vHP = selCurveVerts(web_interface_curves[0][i_curve])
        _, vLP = selCurveVerts(web_interface_curves[1][i_curve])
        bottom_curve,_ = print_sine_curve_between_two_verts(
            vHP, vLP, cs_params["max_web_imperfection_distance"][i_station], "x"
        )

        if i_curve < n_base_curves_web / 2:
            material_name = aft_web_stack.plygroups[i_curve].materialid
            ply_angle = aft_web_stack.plygroups[i_curve].angle
        else:
            material_name = fore_web_stack.plygroups[
                i_curve - int(n_base_curves_web / 2)
            ].materialid
            ply_angle = fore_web_stack.plygroups[i_curve - int(n_base_curves_web / 2)].angle
        part_name_id, materials_used = make_cross_section_surface(lp_hp_side,
            surface_dict,
            i_station,
            part_name,
            top_curve,
            bottom_curve,
            material_name,
            ply_angle,
            part_name_id,
            n_modeled_layers + it + 2 + i_curve,
            materials_used,
            -1
        )
        surf_id = get_last_id("surface")
        surf_name = cubit.get_entity_name("surface", surf_id).split('_')
        surf_name.insert(-1,'vertical')
        surf_name = '_'.join(surf_name)
        cubit.cmd(f'surface {surf_id} rename "{surf_name}"')
    return part_name_id, (vHP, vLP)


def make_a_cross_section(wt_name,
    surface_dict,
    i_station,
    i_station_geometry,
    blade,
    hasWebs,
    aft_web_stack,
    fore_web_stack,
    iLE,
    cs_params,
    geometry_scaling,
    thickness_scaling,
    last_round_station,
    last_flat_station,
    materials_used,
    cs_normal,
):
    geometry = blade.geometry
    stackdb = blade.stackdb
    keypoints = blade.keypoints
    
    if i_station > last_round_station:
        is_flatback=True
    else:
        is_flatback=False

    with open(f"{wt_name}.log", "a") as logFile:
        logFile.write(f"Working on Station: {i_station}\n")

    part_name_id = 0

    #### Step one create outer mold line
    xyz = get_blade_geometry_for_station(blade, i_station_geometry) * geometry_scaling

    # Start indexing from 1 (not 0) to ignore first point: because first point is not on the LP or HP surface but rather is the midpoint at the TE
    splinePoints = xyz[1:iLE, :]
    write_spline_from_coordinate_points(cubit, splinePoints)
    hp_key_curve = get_last_id("curve")

    xyz = np.flip(xyz, 0)
    splinePoints = xyz[1:iLE, :]
    write_spline_from_coordinate_points(cubit, splinePoints)
    lp_key_curve = get_last_id("curve")


    flatback_vBot, _ = selCurveVerts(hp_key_curve)
    flatback_vTop, _ = selCurveVerts(lp_key_curve)

    
    flatbackCurve = cubit.create_curve(
        cubit.vertex(flatback_vBot), cubit.vertex(flatback_vTop)
    )
    flatback_curve_id = flatbackCurve.id()

    #### Extend flatback ###
    curve_start_or_end = "start"
    extension_length = (
        100 * geometry.ichord[i_station_geometry] * cubit.curve(flatback_curve_id).length()
    )
    flatback_curve_id = extend_curve_at_vertex_to_length(
        flatback_curve_id, extension_length, curve_start_or_end
    )

    curve_start_or_end = "end"
    extension_length = 0.5 * cubit.curve(flatback_curve_id).length()
    flatback_curve_id = extend_curve_at_vertex_to_length(
        flatback_curve_id, extension_length, curve_start_or_end
    )

    if is_flatback:
        # Crate camber line
        offset_distance = 0
        npts = 100
        spacing = geometry.ichord[i_station_geometry] * 0.5 / npts

        cubit.cmd(f"curve {flatback_curve_id} copy")
        flatback_offset_curve_id = get_last_id("curve")

        vertex_list = []
        # Do first point manually outside of loop
        flatback_vBot, _ = selCurveVerts(hp_key_curve)
        flatback_vTop, _ = selCurveVerts(lp_key_curve)
        coords_hp = np.array(cubit.vertex(flatback_vBot).coordinates())
        coords_lp = np.array(cubit.vertex(flatback_vTop).coordinates())
        coords = np.mean(np.vstack((coords_hp, coords_lp)), 0)
        cubit.create_vertex(coords[0], coords[1], coords[2])
        vertex_list.append(get_last_id("vertex"))
        for i_point in range(npts):
            offset_distance += spacing
            axial_direction = cs_normal
            position = cubit.curve(flatback_offset_curve_id).position_from_fraction(1.0)
            tangent_direction = cubit.curve(flatback_offset_curve_id).tangent(position)

            normal_direction = crossProd(axial_direction, tangent_direction)
            normal_direction = (
                -1
                * spacing
                * np.array(
                    vectNorm(
                        [normal_direction[0], normal_direction[1], normal_direction[2]]
                    )
                )
            )
            cubit.cmd(
                f"curve {flatback_offset_curve_id} move x {normal_direction[0]} y {normal_direction[1]} z {normal_direction[2]} nomesh"
            )

            cubit.cmd(
                f"create vertex atintersection curve {flatback_offset_curve_id} {hp_key_curve}"
            )
            cubit.cmd(
                f"create vertex atintersection curve {flatback_offset_curve_id} {lp_key_curve}"
            )

            coords_hp = np.array(cubit.vertex(get_last_id("vertex") - 1).coordinates())
            coords_lp = np.array(cubit.vertex(get_last_id("vertex")).coordinates())

            coords = np.mean(np.vstack((coords_hp, coords_lp)), 0)
            cubit.create_vertex(coords[0], coords[1], coords[2])
            vertex_list.append(get_last_id("vertex"))

        cubit.cmd(f"create curve spline vertex {l2s(vertex_list)}")
        cubit.cmd(f"delete vertex {l2s(vertex_list[1:-1])}")
        cubit.cmd(f"delete curve {flatback_offset_curve_id}")

        camberID = get_last_id("curve")


    #     # cubit.cmd(f'save as "Debug.cub" overwrite')
    #     # foo
    else:
        xyz = get_mid_line(blade, iLE, i_station_geometry, geometry_scaling)
        npts, _ = xyz.shape
        npts = round(
            npts * 0.75
        )  # Only need write first 3/4 of camber line since LE is constructed another way

        splinePoints = xyz[0:npts, :]
        write_spline_from_coordinate_points(cubit, splinePoints)
        camberID = get_last_id("curve")

    #Special treatment for station that will connect round to flatback
    #if i_station == last_round_station+1:
    if i_station in list(range(last_round_station+1,last_flat_station+1)):
        cubit.create_curve(cubit.vertex(flatback_vBot), cubit.vertex(flatback_vTop))
        flatback_curve_id=get_last_id("curve")
        v1,v2 = selCurveVerts(flatback_curve_id)
        amplitude=cubit.curve(flatback_curve_id).length()/2
        flatback_curve_id,vertex_list = print_sine_curve_between_two_verts(v1, v2, amplitude, "x")
        
        #### Extend flatback ###
        curve_start_or_end = "start"
        extension_length = (
            100 * geometry.ichord[i_station_geometry] * cubit.curve(flatback_curve_id).length()
        )
        cubit.cmd(f"curve {flatback_curve_id} copy")
        curve_id = extend_curve_at_vertex_to_length(
            get_last_id("curve"), extension_length, curve_start_or_end
        )
        v1,_=selCurveVerts(curve_id)

        curve_start_or_end = "end"
        cubit.cmd(f"curve {flatback_curve_id} copy")
        curve_id = extend_curve_at_vertex_to_length(
            get_last_id("curve"), extension_length, curve_start_or_end
        )
        _,v2=selCurveVerts(curve_id)
        vertex_list=[v1]+vertex_list+[v2]
        cubit.cmd(f"create curve spline vertex {l2s(vertex_list)}")
        flatback_curve_id = get_last_id("curve")

        
    n_stacks = len(stackdb.stacks)

    le_hp_stack_thickness = (
        sum(stackdb.stacks[int(n_stacks / 2.0) - 1, i_station].layer_thicknesses()) / 1000
    )
    le_lp_stack_thickness = (
        sum(stackdb.stacks[int(n_stacks / 2.0), i_station].layer_thicknesses()) / 1000
    )

    # Define variables with HP side in index 0, LP side in index 1

    lp_hp_dict = {}
    lp_hp_dict["spar_cap_base_curves"] = [[], []]
    lp_hp_dict["web_interface_curves"] = [[], []]
    lp_hp_dict["base_curve_ids"] = [[], []]
    lp_hp_dict["round_te_adhesive_curve_list"] = [[], []]
    lp_hp_dict["flat_te_adhesive_curve_list"] = [[], []]
    lp_hp_dict["le_adhesive_curve_list"] = [[], []]
    lp_hp_dict["perimeter_verts_for_te_adhesive"] = {}

    hp_key_curve, lp_key_curve, lp_hp_dict["le_adhesive_curve_ids"] = write_LE_adhesive_curves(
        le_hp_stack_thickness,
        le_lp_stack_thickness,
        cs_params["le_adhesive_thickness"][i_station],
        hp_key_curve,
        lp_key_curve,
        cs_normal,
    )

    key_curves = split_curve_at_coordinte_points(
        keypoints.key_points[1:5, :, i_station_geometry], hp_key_curve
    )
    web_adhesive_width = cs_params["web_adhesive_width"][i_station]
    (
        lp_hp_dict["base_curve_ids"][0],
        lp_hp_dict["spar_cap_base_curves"][0],
    ) = split_key_curves(key_curves, aft_web_stack, fore_web_stack, web_adhesive_width)

    temp = np.flip(keypoints.key_points[:, :, i_station_geometry], 0)
    key_curves = split_curve_at_coordinte_points(temp[1:5, :], lp_key_curve)
    (
        lp_hp_dict["base_curve_ids"][1],
        lp_hp_dict["spar_cap_base_curves"][1],
    ) = split_key_curves(key_curves, aft_web_stack, fore_web_stack, web_adhesive_width)


    # Make sure that the adhesive width is the same on HP and LP sides. Also 
    # Ensure that the adhesive width does not exceed the minimum of the HP and
    # LP reinforcement widths.

    hpTEreinfCurveID=lp_hp_dict["base_curve_ids"][0][0]
    lpTEreinfCurveID=lp_hp_dict["base_curve_ids"][1][0]
    curve_lengths=[cubit.curve(hpTEreinfCurveID).length(),cubit.curve(lpTEreinfCurveID).length()]
    
    index=curve_lengths.index(min(curve_lengths))

    curve_len=min(curve_lengths)
    
    if index==0: #HP curve is the smallest
        stack_thicknesses=stackdb.stacks[1, i_station].layer_thicknesses()
    else:        #LP curve is the smallest
        stack_thicknesses=stackdb.stacks[-2, i_station].layer_thicknesses()

    minimumWidth=np.average(stack_thicknesses)/1000
    if curve_len -cs_params["te_adhesive_width"][i_station] > minimumWidth:
        split_length=cs_params["te_adhesive_width"][i_station]
    else:
        split_length=curve_len-minimumWidth
        with open(f"{wt_name}.log", "a") as logFile:
            logFile.write(f"    Warning: Adhesive width is wider than one of the TE reinforcements. Reducing adhesive width by about {minimumWidth}\n")
    
    if is_flatback:
        cubit.cmd(f"curve {camberID} copy")
        cubit.cmd(f"split curve {get_last_id('curve')}  distance {split_length} from start")
        curve_id=get_last_id("curve")-1

        #Further subdived current_stack_left_curves for roundTEadhesive volumes
        fraction_from_start=0
        n_start = get_last_id("vertex") + 1
        for fraction in []:
        #for fraction in [0.15,0.15,0.15]:
            fraction_from_start+=fraction
            cubit.cmd(f'create vertex on curve {curve_id} fraction {fraction_from_start} from start')
        n_end = get_last_id("vertex")
        vertex_list = list(range(n_start, n_end + 1))


        _,v1=selCurveVerts(curve_id)

        vertex_list.append(v1)
        lp_hp_dict["perimeter_verts_for_te_adhesive"]['HP']=vertex_list
        lp_hp_dict["perimeter_verts_for_te_adhesive"]['LP']=vertex_list

        #cs_params["round_te_adhesive_fractions_for_flatback"]=vertex_list
    else:
        cubit.cmd(f"create vertex on curve {hpTEreinfCurveID}  distance {split_length} from start")
        lp_hp_dict["perimeter_verts_for_te_adhesive"]['HP']=[get_last_id("vertex")]
        cubit.cmd(f"create vertex on curve {lpTEreinfCurveID}  distance {split_length} from start")
        lp_hp_dict["perimeter_verts_for_te_adhesive"]['LP']=[get_last_id("vertex")]

        #cs_params["round_te_adhesive_fractions_for_flatback"]=[]

    # Extend
    curve_start_or_end = "start"
    extension_length = 0.5 * cubit.curve(camberID).length()
    camberID = extend_curve_at_vertex_to_length(camberID, extension_length, curve_start_or_end)

    lp_hp_dict["camberID"] = camberID
    lp_hp_dict["flatback_curve_id"] = flatback_curve_id
    
    
    #For round sections make sure that base curves intersect the camber
    if not is_flatback:
        vHP, _ = selCurveVerts(lp_hp_dict["base_curve_ids"][0][0])
        vLP, _ = selCurveVerts(lp_hp_dict["base_curve_ids"][1][0])

        #Make two lines, each with the correct sense
        cubit.cmd(f'create curve vertex {vLP} {vHP}')
        hp_curve=get_last_id("curve")

        cubit.cmd(f'create curve vertex {vHP} {vLP}')
        lp_curve=get_last_id("curve")

        #HP
        cubit.cmd(f'split curve {hp_curve} fraction 0.4 from end')
        first_curve_id=get_last_id("curve")-1

        cubit.cmd(f'split curve {lp_hp_dict["base_curve_ids"][0][0]} fraction 0.05 from start')
        second_curve_id=get_last_id("curve")

        keep_curve=2
        lp_hp_dict["base_curve_ids"][0][0]=splice_two_curves(first_curve_id, second_curve_id, keep_curve)
        

        #LP
        cubit.cmd(f'split curve {lp_curve} fraction 0.4 from end')
        first_curve_id=get_last_id("curve")-1

        cubit.cmd(f'split curve {lp_hp_dict["base_curve_ids"][1][0]} fraction 0.025 from start')
        second_curve_id=get_last_id("curve")

        keep_curve=2
        lp_hp_dict["base_curve_ids"][1][0]=splice_two_curves(first_curve_id, second_curve_id, keep_curve)

  
    n_modeled_layers = 3

    lp_hp_side = "HP"
    
    stack_ct=0
    part_name_id, lp_hp_dict,stack_ct = make_cs_perimeter_layer_areas(wt_name,
        surface_dict,
        i_station,
        stackdb.stacks[1:6, i_station],
        cs_params,
        thickness_scaling,
        lp_hp_side,
        last_round_station,
        last_flat_station,
        part_name_id,
        n_modeled_layers,
        cs_normal,
        lp_hp_dict,
        materials_used,
        stack_ct
    )
    stack_ct+=1
    
    lp_hp_side = "LP"
    temp = stackdb.stacks[:, i_station]
    temp = np.flip(temp)
    part_name_id, lp_hp_dict,stack_ct = make_cs_perimeter_layer_areas(wt_name,
        surface_dict,
        i_station,
        temp[1:6],
        cs_params,
        thickness_scaling,
        lp_hp_side,
        last_round_station,
        last_flat_station,
        part_name_id,
        n_modeled_layers,
        cs_normal,
        lp_hp_dict,
        materials_used,
        stack_ct
    )

    part_name = "shell"
    part_name_id = create_simplist_surface_for_TE_or_LE_adhesive(
        i_station,
        surface_dict,
        part_name,
        lp_hp_dict["le_adhesive_curve_list"],
        stackdb.stacks[6, i_station].plygroups[0].materialid,
        part_name_id,
        n_modeled_layers,
        materials_used,
        -1
    )
    surf_id=get_last_id("surface")-2
    curve_id = cubit.surface(surf_id).curves()[0].id()
    curve_name_split = cubit.get_entity_name("curve", curve_id).split('_')
    curve_name=curve_name_split[0]+'_'+curve_name_split[1]+'_oml'
    cubit.cmd(f'curve {curve_id} rename "{curve_name}"')


    part_name_id = 0  # Reset since outer areoshell is complete (LE adhesive is accouted for as aeroshell)
    part_name = "roundTEadhesive"
    
    if lp_hp_dict["round_te_adhesive_curve_list"][0] and lp_hp_dict["round_te_adhesive_curve_list"][1]:
        round_te_present=True
        if i_station <= last_round_station:
            mat_name = cs_params["adhesive_mat_name"]
        else:
            mat_name = stackdb.stacks[1, i_station].plygroups[0].materialid
        part_name_id = create_simplist_surface_for_TE_or_LE_adhesive(
            i_station,
            surface_dict,
            part_name,
            lp_hp_dict["round_te_adhesive_curve_list"],
            mat_name,
            part_name_id,
            n_modeled_layers,
            materials_used,
            -1)
        
        surf_id=get_last_id("surface")-2
        curve_id = cubit.surface(surf_id).curves()[0].id()
        curve_name_split = cubit.get_entity_name("curve", curve_id).split('_')
        curve_name=curve_name_split[0]+'_'+curve_name_split[1]+'_oml'
        cubit.cmd(f'curve {curve_id} rename "{curve_name}"')
    else:
        round_te_present=False
        
    if is_flatback:
        part_name = "flatTEadhesive"
        part_name_id = 0  # Reset 
        part_name_id = create_simplist_surface_for_TE_or_LE_adhesive(
            i_station,
            surface_dict,
            part_name,
            lp_hp_dict["flat_te_adhesive_curve_list"],
            cs_params["adhesive_mat_name"],
            part_name_id,
            n_modeled_layers,
            materials_used,
            0
        )

        if not round_te_present:
            surf_id=get_last_id("surface")
            curve_id = cubit.surface(surf_id).curves()[3].id()
            curve_name_split = cubit.get_entity_name("curve", curve_id).split('_')
            curve_name=curve_name_split[0]+'_'+curve_name_split[1]+'_oml'
            cubit.cmd(f'curve {curve_id} rename "{curve_name}"')

    birds_mouth_verts = []
    if hasWebs:
        part_name_id = 0  # Reset since outer areoshell is complete (LE adhesive is accouted for as aeroshell)

        part_name_id, birds_mouth_verts = make_cs_web_layer_areas(
            surface_dict,
            i_station,
            aft_web_stack,
            fore_web_stack,
            lp_hp_dict["web_interface_curves"],
            cs_params,
            part_name_id,
            cs_normal,
            n_modeled_layers,
            materials_used,
        )

    parse_string = f'with name "*station{str(i_station).zfill(3)}*"'
    cs_surfaces = parse_cubit_list("surface", parse_string)
    for surface_id in cs_surfaces:
        n = get_surface_normal(surface_id)
        if n[2] < 0:
            cubit.cmd(f"reverse surface {surface_id}")

    cubit.cmd(f"delete vertex all with Is_Free")

    return birds_mouth_verts


def write_vabs_input(
    surface_dict,
    blade,
    cs_params,
    directory,
    file_name,
    surface_ids,
    materials_used,
    cs_normal,
):
    # Write VABS inputfile
    if cs_params["element_shape"].lower() == "quad":
        expandedConnectivityString = "face"
    elif cs_params["element_shape"].lower() == "tri":
        expandedConnectivityString = "tri"
    else:
        raise NameError(
            f'Element type: {cs_params["element_shape"]} not supported'
        )

    ######Write VABS input file
    nnodes = get_node_count()
    nelem = get_element_count()
    nmate = len(materials_used)
    nlayer = len(surface_ids)  # One VABS layer is defined for each surface

    path_name = directory + "/" + file_name

    if not os.path.exists(directory):
        os.makedirs(directory)

    with open(path_name, "w") as f:
        f.write(f"1 {nlayer}\n")  # New format. One layer definiton for each element.
        f.write("1 0 0\n")
        f.write("0 0 0 0\n")
        f.write(f"{nnodes} {nelem} {nmate}\n")
        # Write Nodes
        for iNd in range(nnodes):
            node_id = iNd + 1
            coords = list(get_nodal_coordinates(node_id))
            f.write(f"{node_id} {coords[0]} {coords[1]}\n")
        f.write("\n\n")
        # Write Elements
        maxNumberOfPossibleNodes = 9
        for iEl in range(nelem):
            element_id = iEl + 1
            nodesIDs = cubit.get_expanded_connectivity(
                cs_params["element_shape"], element_id
            )

            if nodesIDs[0] == 0 or nodesIDs[0] == 0.0:
                foo
            tempStr = str(nodesIDs)  # convert tuple to string
            tempStr = tempStr.replace("(", "")
            tempStr = tempStr.replace(")", "")
            tempStr = tempStr.replace(",", " ")
            nZeros = maxNumberOfPossibleNodes - len(nodesIDs)
            tempStr2 = str([0] * nZeros)
            tempStr2 = tempStr2.replace("[", "")
            tempStr2 = tempStr2.replace("]", "")
            tempStr2 = tempStr2.replace(",", " ")
            f.write(f"{element_id} {tempStr} {tempStr2}\n")
        # Write ply angle for all but the TE adhesive

        for i_surface, surface_id in enumerate(surface_ids):
            surface_name = cubit.get_entity_name("surface", surface_id)
            if 'LP' in surface_name:
                mat_ori_sign = -1.0
            else:
                mat_ori_sign = 1.0
            for iEl, element_id in enumerate(get_surface_quads(surface_id)):
                nodesIDs = cubit.get_expanded_connectivity("face", element_id)
                coords = []
                for iNd, node_id in enumerate(nodesIDs):
                    coords.append(list(get_nodal_coordinates(node_id)))
                coords = np.array(coords)
                # #######For Plotting - find the larges element side length #######
                # distances=[]
                # for iNd,node_id in enumerate(nodesIDs):
                #     for jNd,node_idj in enumerate(nodesIDs):
                #         distances.append(norm(vectSub(coords[iNd],coords[jNd])))
                # length=max(distances)
                # #######For Plotting - find the larges element side length #######
                coords = np.mean(coords, 0)
                # coords=cubit.get_center_point(cs_params['element_shape'], element_id)

                #                             minDist=inf #initialize
                #                             closestCurveID=nan #initialize
                #                             #Since there are possibly many curves due to the offset operation, see which curve is closeset to element center
                #                             for i_curve, curve_id in enumerate(curves):
                #                                 temp=cubit.curve(curve_id).closest_point(coords)
                #                                 distance=getDist(coords,temp)[0]
                #                                 if distance < minDist:
                #                                     minDist=distance
                #                                     closestCurveID=curve_id

                curve_id_for_mat_ori = cubit.surface(surface_id).curves()[0]
                curve_location_for_tangent = curve_id_for_mat_ori.closest_point(coords)
                x = mat_ori_sign*curve_id_for_mat_ori.tangent(curve_location_for_tangent)[0]
                y = mat_ori_sign*curve_id_for_mat_ori.tangent(curve_location_for_tangent)[1]
                z = mat_ori_sign*curve_id_for_mat_ori.tangent(curve_location_for_tangent)[2]
                tangent_direction = vectNorm([x, y, z])  # Unit vector of tangent

                # cross_prod = np.cross(coords,tangent_direction)

                # if np.dot([0,0,1],cross_prod) < 0: 
                #     tangent_direction = vectNorm([-x, -y, -z])  # Unit vector of tangent

                theta1 = math.atan2(tangent_direction[1], tangent_direction[0]) * 180 / pi

                f.write(f"{element_id} {i_surface+1} {theta1}\n")
                # # #######Only needed For Plotting Orientation Check#######
                # cubit.create_vertex(coords[0],coords[1],coords[2])
                # iVert1=get_last_id("vertex")
                # cubit.create_vertex(coords[0]+length*tangent_direction[0],coords[1]+length*tangent_direction[1],coords[2]+length*tangent_direction[2])
                # iVert2=get_last_id("vertex")
                # cubit.cmd(f'create curve vertex {iVert1} {iVert2}')
                # #######Only needed For Plotting Orientation Check#######
                # ##Normal to curve
                # #print(cs_normal)
                
                # # #######Only needed For Plotting Orientation Check#######
                # axial_direction = cs_normal  # There will be a slight error here for highly tapeded regions
                # normal_direction = crossProd(axial_direction, tangent_direction)
                # cubit.create_vertex(coords[0]+length*normal_direction[0],coords[1]+length*normal_direction[1],coords[2]+length*normal_direction[2])
                # cubit.cmd(f'create curve vertex {iVert1} {iVert2}')
                # #######Only needed For Plotting Orientation Check#######
        # Define Plies
        for i_surface, surface_id in enumerate(surface_ids):
            material_id = (
                list(materials_used).index(surface_dict[surface_id]["material_name"]) + 1
            )
            ply_angle = surface_dict[surface_id]["ply_angle"]
            f.write(f"{i_surface+1} {material_id} {ply_angle}\n")
        # Define Materials
        for imat, mat_name in enumerate(materials_used):
            material_id = imat + 1
            material = blade.definition.materials[mat_name]
            f.write(f"{material_id} {1} \n")
            f.write(f"{material.ex} {material.ey} {material.ez}\n")
            f.write(f"{material.gxy} {material.gxz} {material.gyz}\n")
            f.write(f"{material.prxy} {material.prxz} {material.pryz}\n")
            f.write(f"{material.density}\n")
    print("Done writing VABS input")
    return


# Main script fuctions


def get_te_angle(hp_key_curve, lp_key_curve, fraction):
    c1 = cubit.curve(hp_key_curve)
    c2 = cubit.curve(lp_key_curve)

    coords = list(c1.position_from_fraction(fraction))
    v1 = np.array(c1.tangent(coords))
    coords = list(c2.position_from_fraction(fraction))
    v2 = np.array(c2.tangent(coords))

    return math.degrees(
        math.acos(v1.dot(np.transpose(v2)) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    )


def get_mean_layer_thickness_at_station(i_stations):
        parse_string = f'with name "*_thickness{str(i_stations).zfill(3)}*"'
        thickness_curve_ids = parse_cubit_list("curve", parse_string)

        sum_length=0
        for curve_id in thickness_curve_ids:
            sum_length+=cubit.curve(curve_id).length()
        return sum_length/len(thickness_curve_ids)

def get_min_layer_thickness_at_station(i_stations):
        parse_string = f'with name "*layer_thickness{str(i_stations).zfill(3)}*"'
        thickness_curve_ids = parse_cubit_list("curve", parse_string)

        min_length=999999
        for curve_id in thickness_curve_ids:
            curve_len=cubit.curve(curve_id).length()
            if curve_len < min_length:
                min_length=curve_len
        return min_length

def get_locus_of_cross_sectional_centroids(station_list):
    
    # # Adding Nodesets
    x_bar_list=[]
    y_bar_list=[]
    z_bar_list=[]
    n_start=get_last_id("vertex")
    for iLoop, station_id in enumerate(station_list):
            
        parse_string = f'with name "*station{str(station_id).zfill(3)}*_surface*"'
        surface_ids = parse_cubit_list("surface", parse_string)

        #Find centroid
        x_bar=0
        y_bar=0
        z_bar=0
        total_area=0
        for surface_id in surface_ids:
            coords = get_surface_centroid(surface_id)
            area=cubit.surface(surface_id).area()

            total_area+=area
            x_bar+=coords[0]*area
            y_bar+=coords[1]*area
            z_bar+=coords[2]*area

        x_bar=x_bar/total_area
        y_bar=y_bar/total_area
        z_bar=z_bar/total_area

        x_bar_list.append(x_bar)
        y_bar_list.append(y_bar)
        z_bar_list.append(z_bar)
    n_end=get_last_id("vertex")


    centroidal_ref_line_coords = np.vstack([x_bar_list,y_bar_list,z_bar_list]).transpose()

    if len(station_list)>0:
        if len(station_list)==1:
            cubit.cmd(f"create vertex location {x_bar_list[0]} {y_bar_list[0]} {z_bar_list[0]}")
            return_index=get_last_id('vertex')
        else:
            write_spline_from_coordinate_points(cubit, centroidal_ref_line_coords)
            return_index=get_last_id('curve')


    return return_index,centroidal_ref_line_coords