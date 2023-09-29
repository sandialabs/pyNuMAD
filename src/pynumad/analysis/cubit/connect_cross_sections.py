from cubit import *
from PyCubed_Main import *
from pynumad.analysis.cubit.make_cross_sections import print_sine_curve_between_two_verts
import numpy as np
import re


def get_ordered_list(part_name):

    ordered_list = []
    surfaces_to_connect = [1]  # Initialize to enter loop
    i_surface = -1  # Initialize
    while surfaces_to_connect:
        i_surface += 1
        parse_string = f'with name "*{part_name}*surface{i_surface+1}"'
        surfaces_to_connect = parse_cubit_list("surface", parse_string)

        if surfaces_to_connect:
            ordered_list.append(surfaces_to_connect)

    return ordered_list


def make_spanwise_splines(surface_dict, ordered_list):
    spanwise_splines = []
    for aligned_surfaces in ordered_list:
        tempList = []
        for i_point in range(4):
            vertex_list = []
            for index, i_surface in enumerate(aligned_surfaces):
                vertex_id = surface_dict[i_surface]["verts"][i_point]
                vertex_list.append(vertex_id)
                vertex_name = cubit.get_entity_name("vertex", vertex_id)

            curve = cubit.cmd(f"create curve spline vertex {l2s(vertex_list)}")
            tempList.append(get_last_id("curve"))
        spanwise_splines.append(tempList)
    return spanwise_splines


def make_a_volume(
    current_surface_id, next_surface_id, spanwise_splines_for_a_volume, surface_dict, i_station_end
):
    cubit.cmd(f"surface {current_surface_id} copy")
    current_surface_id_copy = get_last_id("surface")

    cubit.cmd(f"surface {next_surface_id} copy")
    next_surface_id_copy = get_last_id("surface")

    current_surface = cubit.surface(current_surface_id)
    next_surface = cubit.surface(next_surface_id)

    current_surface_curves = surface_dict[current_surface_id]["curves"]
    next_surface_curves = surface_dict[next_surface_id]["curves"]

    current_surface_vertices = surface_dict[current_surface_id]["verts"]
    next_surface_vertices = surface_dict[next_surface_id]["verts"]

    spanwise_splines_for_a_volume.append(
        spanwise_splines_for_a_volume[0]
    )  # Make list circle back

    transverse_surface_ids = []
    for i_curve in range(len(current_surface_curves)):
        cubit.cmd(
            f"create surface curve {current_surface_curves[i_curve]} {spanwise_splines_for_a_volume[i_curve]} {next_surface_curves[i_curve]} {spanwise_splines_for_a_volume[i_curve+1]}"
        )
        transverse_surface_ids.append(get_last_id("surface"))

    surf_name = cubit.get_entity_name("surface", current_surface.id()).split("_")
    regex = re.compile("layer")
    layer_name = [string for string in surf_name if re.match(regex, string)][0]
    string_name = layer_name + "_bottomFace"
    cubit.cmd(f'surface {transverse_surface_ids[0]} rename "{string_name}"')
    string_name = layer_name + "_topFace"
    cubit.cmd(f'surface {transverse_surface_ids[2]} rename "{string_name}"')

    # cubit.cmd(f'save as "python1.cub" overwrite')
    # raise Exception(f'Volume "{volume_name}" creation failed')
    # Create Volume
    # n_start=get_last_id("volume")
    cubit.cmd(
        f"create volume surface {current_surface_id_copy} {next_surface_id_copy} {l2s(transverse_surface_ids)} noheal"
    )
    # n_end=get_last_id("volume")
    # print(f'n_start: {n_start}, n_end: {n_end}')

    if "Station" + str(i_station_end) in cubit.get_entity_name(
        "surface", next_surface_id
    ):  # This if statement is needed for componets that may have been droped between the last station and the second to last station
        volume_name = cubit.get_entity_name("surface", next_surface_id)
    else:
        volume_name = cubit.get_entity_name("surface", current_surface_id)
    if len(cubit.volume(get_last_id("volume")).surfaces()) < 6:
        print(
            f"\n\n ERROR with:\n\n create volume surface {current_surface_id_copy} {next_surface_id_copy} {l2s(transverse_surface_ids)} "
        )
        print(f"current_surface_id_copy: {current_surface_id_copy}")
        print(f"next_surface_id_copy: {next_surface_id_copy}")
        print(f"spanwise_splines_for_a_volume: {spanwise_splines_for_a_volume}")
        cubit.cmd(f'save as "python1.cub" overwrite')
        raise Exception(f'Volume "{volume_name}" creation failed')

    volume_name = volume_name.replace("surface", "volume")
    cubit.cmd(f'volume {get_last_id("volume")} rename "{volume_name}"')


def get_spanwise_splines_for_a_volume(
    i_span, n_cross_sections, spanwise_splines_for_a_surface, next_surface_vertices
):
    # Split off spanwise curves for a single volume and store them
    if i_span < n_cross_sections - 2:
        spanwise_splines_for_a_volume = []
        temp = []
        for i_curve, curve_id in enumerate(spanwise_splines_for_a_surface):
            cubit.cmd(f"split curve {curve_id} at vertex {next_surface_vertices[i_curve]}")
            temp.append(get_last_id("curve"))
            spanwise_splines_for_a_volume.append(get_last_id("curve") - 1)
        spanwise_splines_for_a_surface = temp
    else:
        spanwise_splines_for_a_volume = spanwise_splines_for_a_surface
    return spanwise_splines_for_a_volume, spanwise_splines_for_a_surface


# def assign_intervals(vol_id,nIntervals):
#     thickness_curve_id=cubit.volume(vol_id).curves()[1].id()
#     #cubit.cmd(f'locate curve {thickness_curve_id} ')
#     cubit.cmd(f'curve {thickness_curve_id} interval {nIntervals}')


def make_all_volumes_for_a_part(surface_dict, ordered_list, mesh_vol_list, i_station_end):
    # nIntervals=3
    spanwise_splines = make_spanwise_splines(surface_dict, ordered_list)
    n_cross_sections = len(ordered_list[0])
    nPartSurfaceIDs = len(ordered_list)
    if n_cross_sections > 1:
        for i_span in range(n_cross_sections - 1):
            for part_surface_ids in range(nPartSurfaceIDs):
                current_surface_id = ordered_list[part_surface_ids][i_span]
                next_surface_id = ordered_list[part_surface_ids][i_span + 1]
                (
                    spanwise_splines_for_a_volume,
                    spanwise_splines[part_surface_ids],
                ) = get_spanwise_splines_for_a_volume(
                    i_span,
                    n_cross_sections,
                    spanwise_splines[part_surface_ids],
                    surface_dict[next_surface_id]["verts"],
                )
                make_a_volume(
                    current_surface_id,
                    next_surface_id,
                    spanwise_splines_for_a_volume,
                    surface_dict,
                    i_station_end,
                )
                mesh_vol_list.append(get_last_id("volume"))
                # assign_intervals(get_last_id("volume"),nIntervals)
    else:
        raise ValueError("Can't make volumes with only one cross section.")

    return mesh_vol_list


def verify_web_cutting_amplitude(
    blade, amplitude, tolerance, i_station_first_web, i_station_last_web
):
    # Check to make sure that the amplitude does not result sharp volumes by cutting near a station location
    geometry = blade.geometry
    for i_station_check in range(i_station_first_web + 1, i_station_last_web + 1):
        blade_segment_length = (
            geometry.ispan[i_station_check] - geometry.ispan[i_station_first_web]
        )
        gap = blade_segment_length - amplitude
        # print(f'blade_segment_length: {blade_segment_length}\ngap {gap}')

        if abs(gap) > tolerance:
            break
        else:
            if gap > 0:
                amplitude = blade_segment_length - tolerance
            else:
                amplitude = blade_segment_length + tolerance
            break
    # print(f'new amplitude {amplitude} \nnew gap = {blade_segment_length-amplitude}')
    return amplitude


def make_birds_mouth(
    blade, birds_mouth_verts, birds_mouth_amplitude_fraction, i_station_first_web, i_station_last_web
):
    ### Make birds mouth volume that will cut the web volumes ###
    #############################################################

    # This function must be ran before merging the volumes since "birds_mouth_verts" will change during mergeing

    geometry = blade.geometry

    v1 = cubit.vertex(birds_mouth_verts[0])
    v2 = cubit.vertex(birds_mouth_verts[1])
    distance = getDist(v1.coordinates(), v2.coordinates())[0]
    create_curve(v1, v2)

    # Make the birds mouth cut-out start 5% from where the web meets the aeroshell
    cubit.cmd(
        f'create vertex on curve {get_last_id("curve")}  distance {0.05*distance} from start'
    )
    cubit.cmd(
        f'create vertex on curve {get_last_id("curve")}  distance {0.05*distance} from end'
    )
    v1 = cubit.vertex(get_last_id("vertex") - 1)
    v2 = cubit.vertex(get_last_id("vertex"))
    straight_line = create_curve(v1, v2)

    amplitude = birds_mouth_amplitude_fraction * distance
    tolerance = distance * 0.05

    amplitude = verify_web_cutting_amplitude(
        blade, amplitude, tolerance, i_station_first_web, i_station_last_web
    )

    curvedLine = cubit.curve(
        print_sine_curve_between_two_verts(v1.id(), v2.id(), amplitude, "z")
    )
    cubit.cmd(f"create surface skin curve {curvedLine.id()} {straight_line.id()}")
    baseSurface = get_last_id("surface")

    midPoint = list(curvedLine.position_from_fraction(0.5))
    tangent = straight_line.tangent(midPoint)

    # Get the cross-section normal
    parse_string = f'in surface with name "*webStation{i_station_first_web}*"'
    surface_id = parse_cubit_list("surface", parse_string)[
        0
    ]  # Pick the first surface in this list since all on same plane
    coords = cubit.get_center_point("surface", surface_id)
    surface_normal = cubit.surface(surface_id).normal_at(coords)
    cut_block_length = 5 * max(geometry.ichord)
    sweep_direction = np.array(vectNorm(crossProd(list(tangent), list(surface_normal))))

    cubit.cmd(
        f"sweep surface {baseSurface} direction {l2s(sweep_direction)} distance {cut_block_length}"
    )
    cubit.cmd(
        f'move volume {get_last_id("volume")} x {-cut_block_length/2*sweep_direction[0]} y {-cut_block_length/2*sweep_direction[1]} z {-cut_block_length/2*sweep_direction[2]}'
    )

    cutting_volume = get_last_id("volume")

    parse_string = f'with name "*webStation*"'
    web_volumes = parse_cubit_list("volume", parse_string)

    n_start = get_last_id("volume")
    cubit.cmd(f"subtract volume {cutting_volume} from volume {l2s(web_volumes)}")
    n_end = get_last_id("volume")
    
    #Add thickness_curves to ensure each volume has 4 of them
    vol_list = list(range(n_start + 1, n_end + 1))
    
    for volume_id in vol_list:
        
        #Get list of curves with name layer_thickness
        parse_string = f'with name "*layer_thickness*" in volume {volume_id}'
        thickness_curve_ids = parse_cubit_list("curve", parse_string)
        
        #All curves in volume
        parse_string = f'in volume {volume_id}'
        volume_curves=parse_cubit_list("curve", parse_string)

        if len(thickness_curve_ids) < 4:
            if len(thickness_curve_ids) > 0: #Make sure at least one thickness curve is defined

                #All curves in volume
                parse_string = f'in volume {volume_id}'
                volume_curves=set(parse_cubit_list("curve", parse_string))
                
                current_curve=cubit.curve(thickness_curve_ids[0])
                coords = current_curve.position_from_fraction(0.5)
                approximate_thickness_direction=current_curve.tangent(coords)

                for i_curve in set(volume_curves).difference(set(thickness_curve_ids)):
                    current_curve=cubit.curve(i_curve)
                    coords = current_curve.position_from_fraction(0.5)
                    curve_tangent=current_curve.tangent(coords)
                    if 1.0 - abs(np.dot(curve_tangent,approximate_thickness_direction)) < 0.15:
                        cubit.cmd(f'curve {i_curve} rename "layer_thickness"')
                
                #Get list of curves with name layer_thickness
                parse_string = f'with name "*layer_thickness*" in volume {volume_id}'
                thickness_curve_ids = parse_cubit_list("curve", parse_string)
                if len(thickness_curve_ids) !=4:
                    raise ValueError(
                        f"Something wrong with thickness curves in birds mounth volume"
                    )
            else:
                raise ValueError(
                    f"Zero thickness curves found in a web volume {volume_id} near birds mounth"
                )

        
    return


# cubit.cmd('open "/home/ecamare/myprojects/bar/cubitDev/python/python0.cub"')


def get_approximate_thickness_direction_for_volume(volume_id):
    # This function is used when assigning material orientations

    # Get thickness direction tangents
    approximate_thickness_direction = []
    for current_curve in cubit.volume(volume_id).curves():
        curve_name = cubit.get_entity_name("curve", current_curve.id())
        if "layer_thickness" in curve_name:
            coords = current_curve.position_from_fraction(0.5)
            approximate_thickness_direction.append(current_curve.tangent(coords))
    approximate_thickness_direction = np.array(approximate_thickness_direction)
    nThicknessCurves, _ = approximate_thickness_direction.shape

    if nThicknessCurves == 4:  # All other cases
        return np.mean(approximate_thickness_direction, 0)
    elif nThicknessCurves == 8:  # LE adhesive case and round TE adhesive
        return 0
    elif nThicknessCurves == 6:  # Web overwrap
        # Take the mean of all curves with name 'layer_thickness'
        mean = np.mean(approximate_thickness_direction, 0)

        errorList = []
        for i in range(nThicknessCurves):
            diff = approximate_thickness_direction[i] - mean

            errorList.append(sqrt(dotProd(diff, diff)))
        sortIndex = np.argsort(errorList)[
            :4
        ]  # Take the first four. This discards the two directions with the largest deviation from the average

        return np.mean(approximate_thickness_direction[sortIndex, :], 0)
    else:
        cubit.cmd(f'save as "Debug.cub" overwrite')
        raise ValueError(
            f"The number of thickness curves in volume is unexpected. Cannot assign material orientation. nThicknessCurves: {nThicknessCurves}"
        )

    return


def get_mat_ori_surface(volume_id, spanwise_mat_ori_curve):
    # This function is used when assigning material orientations
    # This gets returns the surface within a volume that will be used to get surface normals.
    # The sign +-1 is also returned since some of the surfaces are oriented the wrong way

    approximate_thickness_direction = get_approximate_thickness_direction_for_volume(volume_id)

    # Create a list of surface IDs in the given volume
    surface_ids = []
    volumeSurfaces = cubit.volume(volume_id).surfaces()
    for current_surface in volumeSurfaces:
        surface_ids.append(current_surface.id())

    # Eliminate surfaces that have two curves named thickness:
    surface_ct = 0
    for current_surface in volumeSurfaces:
        curve_ct = (
            0  # Counts the number of curves in the surface with name 'layer_thickness'
        )
        for current_curve in current_surface.curves():
            curve_name = cubit.get_entity_name("curve", current_curve.id())
            if "layer_thickness" in curve_name:
                curve_ct += 1

        if curve_ct >= 2:
            surface_ct += 1
            surface_ids.remove(current_surface.id())

    # surface_ids now has the list of surfaces w/o thickness curves
    if len(surface_ids) == 2 or len(surface_ids) == 1:
        if len(surface_ids) == 2:
            surface_name = cubit.get_entity_name("surface", surface_ids[0])
            if "topFace" in surface_name:
                surface_id = surface_ids[0]
            else:
                surface_id = surface_ids[-1]
        elif len(surface_ids) == 1:  # Web overwrap
            surface_id = surface_ids[0]

        coords = cubit.get_center_point("surface", surface_id)
        surface_normal = cubit.surface(surface_id).normal_at(coords)

        if dotProd(surface_normal, approximate_thickness_direction) > 0:
            sign = 1.0
        else:
            sign = -1.0
    elif (
        len(surface_ids) == 0
    ):  # LE adhesive and/or TE adhesive for round cross-sections
        # print(f'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~volume_id {volume_id}')
        surface_id = False
        sign = 1.0

    else:
        raise ValueError(
            "The number of thickness curves in volume is unexpected. Cannot assign material orientation"
        )

    return surface_id, sign
