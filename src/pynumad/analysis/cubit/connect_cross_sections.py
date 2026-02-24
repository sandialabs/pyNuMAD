import warnings

try:
    from cubit import *
    from PyCubed_Main import *
except ModuleNotFoundError:
    warnings.warn("Cubit not configured, so cubit functionality will not be available.")

from pynumad.analysis.cubit.make_cross_sections import print_sine_curve_between_two_verts
import numpy as np
import re

def debug():
    cubit.cmd(f"delete curve 1")
    cubit.cmd(f'save as "Debug.cub" overwrite')
def get_ordered_list(part_name):

    ordered_list = []
    surfaces_to_connect = [1]  # Initialize to enter loop
    i_surface = -1  # Initialize
    temp=[]
    while surfaces_to_connect:
        i_surface += 1
        parse_string = f'with name "{part_name}*Station*surface{i_surface+1}"'
        surfaces_to_connect = parse_cubit_list("surface", parse_string)
        if len(surfaces_to_connect)==4:
            temp.append(surfaces_to_connect)
        if surfaces_to_connect:
            ordered_list.append(surfaces_to_connect)

    return ordered_list


def make_spanwise_splines(surface_dict, ordered_list,spanwise_splines):
    flat_list = [x for xs in spanwise_splines for x in xs]

    for aligned_surfaces in ordered_list:

        tempList = []
        for i_point in range(4):
            #Check if first point lies on any prior splines
            vertex_id = surface_dict[aligned_surfaces[0]]["verts"][i_point]

            min_dist=1 #Initialize
            for spanwise_spline in flat_list:
                distance=get_distance_between_entities('curve', spanwise_spline, 'vertex',vertex_id)
                if distance < min_dist:
                    min_dist = distance
                    save_curve = spanwise_spline
                    #cubit.cmd(f"curve {spanwise_spline} copy")
                    #curve_id = get_last_id("curve")  
                    #break

            if min_dist < 0.0001:
                #cubit.cmd(f"curve {save_curve} copy")
                curve_id = save_curve
            else:
                vertex_list = []
                for index, i_surface in enumerate(aligned_surfaces):
                    vertex_id = surface_dict[i_surface]["verts"][i_point]
                    vertex_list.append(vertex_id)
                    vertex_name = cubit.get_entity_name("vertex", vertex_id)
                    if 'linear' in vertex_name:

                        try:
                            vertex_id2 = surface_dict[aligned_surfaces[index+1]]["verts"][i_point]
                            cubit.cmd(f"create curve vertex {vertex_id} {vertex_id2}")

                            n_start = get_last_id("vertex") + 1
                            cubit.cmd(f'create vertex on curve {get_last_id("curve")} segment 10')
                            n_end = get_last_id("vertex")
                            vertex_list += list(range(n_start, n_end + 1))
                        except:
                            pass

                cubit.cmd(f"create curve spline vertex {l2s(vertex_list)}")   
                curve_id = get_last_id("curve")
                #unique_spanwise_splines.append(curve_id)
                #cubit.cmd(f"curve {curve_id} copy")  
                #curve_id = get_last_id("curve")  
                cubit.cmd(f'curve {curve_id} rename "span_dir"')       

            tempList.append(curve_id)
        spanwise_splines.append(tempList)
    return spanwise_splines

def old_make_spanwise_splines(surface_dict, ordered_list):
    spanwise_splines = []
    for aligned_surfaces in ordered_list:

        tempList = []
        for i_point in range(4):
            vertex_list = []
            for index, i_surface in enumerate(aligned_surfaces):
                vertex_id = surface_dict[i_surface]["verts"][i_point]
                vertex_list.append(vertex_id)
                vertex_name = cubit.get_entity_name("vertex", vertex_id)
                if 'linear' in vertex_name:

                    try:
                        vertex_id2 = surface_dict[aligned_surfaces[index+1]]["verts"][i_point]
                        cubit.cmd(f"create curve vertex {vertex_id} {vertex_id2}")

                        n_start = get_last_id("vertex") + 1
                        cubit.cmd(f'create vertex on curve {get_last_id("curve")} segment 10')
                        n_end = get_last_id("vertex")
                        vertex_list += list(range(n_start, n_end + 1))
                    except:
                        pass

            cubit.cmd(f"create curve spline vertex {l2s(vertex_list)}") 
            cubit.cmd(f'curve {get_last_id("curve")} rename "span_dir"')       

            tempList.append(get_last_id("curve"))
        spanwise_splines.append(tempList)
    return spanwise_splines
def make_a_volume(
    i_span,current_surface_id, next_surface_id, spanwise_splines_for_a_volume, surface_dict, i_station_end,ply_angle
):
    cubit.cmd(f"surface {current_surface_id} copy")
    current_surface_id_copy = get_last_id("surface")

    cubit.cmd(f"surface {next_surface_id} copy")
    next_surface_id_copy = get_last_id("surface")

    current_surface = cubit.surface(current_surface_id)

    current_surface_curves = surface_dict[current_surface_id]["curves"]
    next_surface_curves = surface_dict[next_surface_id]["curves"]

    
    surf_name = cubit.get_entity_name("surface", current_surface_id)
    i_station = surf_name.split('_')[0][-3:]
    cubit.cmd(f'curve {l2s(spanwise_splines_for_a_volume)} rename "span_dir{str(i_station).zfill(3)}"')    

    spanwise_splines_for_a_volume.append(spanwise_splines_for_a_volume[0])  # Make list circle back

    transverse_surface_ids = []
    for i_curve in range(len(current_surface_curves)):
        cubit.cmd(
            f"create surface curve {current_surface_curves[i_curve]} {spanwise_splines_for_a_volume[i_curve]} {next_surface_curves[i_curve]} {spanwise_splines_for_a_volume[i_curve+1]}"
        )
        transverse_surface_ids.append(get_last_id("surface"))
    
    surf_name = cubit.get_entity_name("surface", current_surface.id())
    surf_name_split = surf_name.split("_")
    
    layer_name = surf_name_split[0]+'_'+surf_name_split[1]+'_'+surf_name_split[2]
    string_name = layer_name + "_topFace"
    cubit.cmd(f'surface {transverse_surface_ids[2]} rename "{string_name}"')

    if 'web_thickness' in surf_name:
        append_str='_web_thickness'
    else:
        append_str=''

    string_name = layer_name + "_bottomFace" + append_str
    cubit.cmd(f'surface {transverse_surface_ids[0]} rename "{string_name}"')

    cubit.cmd(f"create volume surface {current_surface_id_copy} {next_surface_id_copy} {l2s(transverse_surface_ids)} noheal")


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
        cubit.cmd(f'save as "Debug.cub" overwrite')
        raise Exception(f'Volume "{volume_name}" creation failed')

    volume_name = volume_name.replace("surface", "volume")
    cubit.cmd(f'volume {get_last_id("volume")} rename "{volume_name}"')

def exp_make_a_volume(
    i_span,current_surface_id, next_surface_id, spanwise_splines_for_a_volume, surface_dict, i_station_end
):
    cubit.cmd(f"surface {current_surface_id} copy")
    current_surface_id_copy = get_last_id("surface")

    cubit.cmd(f"surface {next_surface_id} copy")
    next_surface_id_copy = get_last_id("surface")

    current_surface = cubit.surface(current_surface_id)

    current_surface_curves = surface_dict[current_surface_id]["curves"]
    next_surface_curves = surface_dict[next_surface_id]["curves"]
    
    surf_name = cubit.get_entity_name("surface", current_surface_id)
    i_station = surf_name.split('_')[0][-3:]
    cubit.cmd(f'curve {l2s(spanwise_splines_for_a_volume)} rename "span_dir{str(i_station).zfill(3)}"')    

    transverse_surface_ids = []
    temp=[]

    set_overlap_max_gap(0.0005)
    for i_curve in range(len(current_surface_curves)):
        cubit.cmd(f"create surface curve {current_surface_curves[i_curve]} {spanwise_splines_for_a_volume[i_curve]} {next_surface_curves[i_curve]} {spanwise_splines_for_a_volume[(i_curve+1)%4]}")

        surf_id = get_last_id("surface")

        if False:
            overlapping_surfs = get_overlapping_surfaces_at_surface(surf_id,[])

            if current_surface_id_copy == 755:
                print()

            if len(overlapping_surfs) >0:

                #Vertecies of surf_id
                vertex_ids = selCurveVerts(cubit.curve(current_surface_curves[i_curve]).id())
                vertex_ids += selCurveVerts(cubit.curve(next_surface_curves[i_curve]).id())


                #Find out which surface shares all vertex coords with new surface: surf_id

                #First put all coords in all_coords
                all_coords =np.zeros((len(vertex_ids),3))

                for i_vert, vertex_id in enumerate(vertex_ids):
                    all_coords[i_vert,:]=cubit.vertex(vertex_id).coordinates()
                
                
                overlapping_verts_found = False
                for overlapping_surf in overlapping_surfs:
                    overlapping_verts=[]
                    for vertex in cubit.surface(overlapping_surf).vertices():
                        coords= np.array(vertex.coordinates())
                        diff = coords-all_coords
                        dist = []
                        for i in range(len(vertex_ids)):
                            dist.append(np.linalg.norm(diff[i,:]))
                        if len(np.where(np.array(dist)==0.0)[0]) == 1: 
                            overlapping_verts.append(vertex.id())
                        else:
                            break
                        print()
                    if len(overlapping_verts) ==4: 
                        overlapping_verts_found = True
                        cubit.cmd(f"surface {overlapping_surf} copy")
                        surf_id = get_last_id("surface")
                        break
                if overlapping_verts_found:
                    temp.append(surf_id)
                else:
                    temp.append(0.0)
                        
                    
            else:
                temp.append(0.0)
            transverse_surface_ids.append(surf_id)

        print()

        if sum(temp)>0:
            transverse_surface_ids=[]

            for i_surf,surface_id in enumerate(temp): 

                if surface_id:
                    if i_surf == 0:
                        index_2 = 0
                        index_1 = 2
                    elif i_surf == 1:
                        index_2 = 0
                        index_1 = 2
                    elif i_surf == 2:
                        index_1 = 0
                        index_2 = 2 
                    elif i_surf == 3:
                        index_1 = 0
                        index_2 = 2 
                    else:
                        foo


                    curve_id = cubit.surface(surface_id).curves()[index_1].id()
                    cubit.cmd(f'curve {curve_id} copy')
                    spanwise_splines_for_a_volume[(i_surf+1)%4] = get_last_id("curve")

                    curve_id = cubit.surface(surface_id).curves()[index_2].id()
                    cubit.cmd(f'curve {curve_id} copy')                
                    spanwise_splines_for_a_volume[i_surf] = get_last_id("curve")
    

            for i_curve in range(len(current_surface_curves)):
                cubit.cmd(f"create surface curve {current_surface_curves[i_curve]} {spanwise_splines_for_a_volume[i_curve]} {next_surface_curves[i_curve]} {spanwise_splines_for_a_volume[(i_curve+1)%4]}")
                transverse_surface_ids.append(get_last_id("surface"))

    surf_name = cubit.get_entity_name("surface", current_surface.id())
    surf_name_split = surf_name.split("_")
    


    layer_name = surf_name_split[0]+'_'+surf_name_split[1]+'_'+surf_name_split[2]
    string_name = layer_name + "_topFace"
    cubit.cmd(f'surface {transverse_surface_ids[2]} rename "{string_name}"')

    if 'web_thickness' in surf_name:
        append_str='_web_thickness'
    else:
        append_str=''

    string_name = layer_name + "_bottomFace" + append_str
    cubit.cmd(f'surface {transverse_surface_ids[0]} rename "{string_name}"')



    cubit.cmd(f"create volume surface {current_surface_id_copy} {next_surface_id_copy} {l2s(transverse_surface_ids)} noheal")
    vol_id=get_last_id("volume")

    if "Station" + str(i_station_end) in cubit.get_entity_name(
        "surface", next_surface_id
    ):  # This if statement is needed for componets that may have been droped between the last station and the second to last station
        volume_name = cubit.get_entity_name("surface", next_surface_id)
    else:
        volume_name = cubit.get_entity_name("surface", current_surface_id)
    if len(cubit.volume(vol_id).surfaces()) < 6:
        print(
            f"\n\n ERROR with:\n\n create volume surface {current_surface_id_copy} {next_surface_id_copy} {l2s(transverse_surface_ids)} "
        )
        print(f"current_surface_id_copy: {current_surface_id_copy}")
        print(f"next_surface_id_copy: {next_surface_id_copy}")
        print(f"spanwise_splines_for_a_volume: {spanwise_splines_for_a_volume}")
        cubit.cmd(f'save as "Debug.cub" overwrite')
        raise Exception(f'Volume "{volume_name}" creation failed')

    volume_name = volume_name.replace("surface", "volume")
    cubit.cmd(f'volume {get_last_id("volume")} rename "{volume_name}"')


def get_spanwise_splines_for_a_volume(current_surface_id,next_surface_id,spanwise_splines_for_a_surface,current_surface_vertices,next_surface_vertices):

    # Split off spanwise curves for a single volume and store them

    spanwise_splines_for_a_volume = []

    for i_curve, curve_id in enumerate(spanwise_splines_for_a_surface):
        
        cubit.cmd(f"curve {curve_id} copy")
        curve_id = get_last_id("curve")


        v_1,v_2 = selCurveVerts(curve_id)
        distance=get_distance_between_entities('surface', current_surface_id, 'vertex',v_1)
        if distance > 0.0001:
            cubit.cmd(f"split curve {curve_id} at vertex {current_surface_vertices[i_curve]}")
            curve_id = get_last_id("curve")
            v_1,v_2 = selCurveVerts(curve_id)

        distance=get_distance_between_entities('surface', next_surface_id, 'vertex',v_2)
        if distance > 0.0001:
            cubit.cmd(f"split curve {curve_id} at vertex {next_surface_vertices[i_curve]}")
            curve_id = get_last_id("curve") -1


        spanwise_splines_for_a_volume.append(curve_id)

    return spanwise_splines_for_a_volume



def make_all_volumes_for_a_part(surface_dict, volume_dict,ordered_list, i_station_end,spanwise_splines):
    # nIntervals=3
    vol_list=[]
    i_start = len(spanwise_splines)
    spanwise_splines = make_spanwise_splines(surface_dict, ordered_list,spanwise_splines)
    i_end = len(spanwise_splines)
    n_cross_sections = len(ordered_list[0])
    nPartSurfaceIDs = len(ordered_list)
    if n_cross_sections > 1:
        for i_span in range(n_cross_sections - 1):
            for i_surf, part_surface_ids in enumerate(range(i_start,i_end)):
                current_surface_id = ordered_list[i_surf][i_span]
                next_surface_id = ordered_list[i_surf][i_span + 1]

                spanwise_splines_for_a_volume = get_spanwise_splines_for_a_volume(current_surface_id,next_surface_id,spanwise_splines[part_surface_ids],
                surface_dict[current_surface_id]["verts"],surface_dict[next_surface_id]["verts"])
                make_a_volume(i_span,current_surface_id,next_surface_id,spanwise_splines_for_a_volume,surface_dict,i_station_end,surface_dict[current_surface_id]["ply_angle"])
                vol_id = get_last_id("volume")
                vol_list.append(vol_id)

                volume_dict[vol_id] = {}
                # volume_dict[vol_id]["material_name"] = material_name
                volume_dict[vol_id]["ply_angle"] = surface_dict[current_surface_id]["ply_angle"]
                del surface_dict[current_surface_id]
                
    else:
        raise ValueError("Can't make volumes with only one cross section.")

    return vol_list,spanwise_splines


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
        f'create vertex on curve {get_last_id("curve")}  distance {0.10*distance} from start'
    )
    cubit.cmd(
        f'create vertex on curve {get_last_id("curve")}  distance {0.10*distance} from end'
    )
    v1 = cubit.vertex(get_last_id("vertex") - 1)
    v2 = cubit.vertex(get_last_id("vertex"))
    straight_line = create_curve(v1, v2)

    amplitude = birds_mouth_amplitude_fraction * distance
    tolerance = distance * 0.05

    amplitude = verify_web_cutting_amplitude(
        blade, amplitude, tolerance, i_station_first_web, i_station_last_web
    )
    
    curve_id,_=print_sine_curve_between_two_verts(v1.id(), v2.id(), amplitude, "z")
    curvedLine = cubit.curve(curve_id)

    cubit.cmd(f"create surface skin curve {curvedLine.id()} {straight_line.id()}")
    baseSurface = get_last_id("surface")

    midPoint = list(curvedLine.position_from_fraction(0.5))
    tangent = straight_line.tangent(midPoint)

    # Get the cross-section normal
    parse_string = f'in surface with name "*web*Station{str(i_station_first_web).zfill(3)}*"'
    surface_id = parse_cubit_list("surface", parse_string)[0]  # Pick the first surface in this list since all on same plane
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

    parse_string = f'with name "*web*Station*"'
    web_volumes = parse_cubit_list("volume", parse_string)

    n_start = get_last_id("volume")
    cubit.cmd(f"subtract volume {cutting_volume} from volume {l2s(web_volumes)}")
    n_end = get_last_id("volume")
    
    #Add thickness_curves to ensure each volume has 4 of them
    vol_list = list(range(n_start + 1, n_end + 1))
    
    for volume_id in vol_list:
        
        #Get list of curves with name layer_thickness
        parse_string = f'with name "*thickness*" in volume {volume_id}'
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
                parse_string = f'with name "*thickness*" in volume {volume_id}'
                thickness_curve_ids = parse_cubit_list("curve", parse_string)
                if len(thickness_curve_ids) !=4:
                    raise ValueError(
                        f"Something wrong with thickness curves in birds mounth volume"
                    )
            else:
                raise ValueError(
                    f"Zero thickness curves found in a web volume {volume_id} near birds mounth"
                )
    parse_string = f'with name "*web*Station*"'
    
    return list(parse_cubit_list("volume", parse_string)) #update the web volumes




def get_mat_ori_surface(volume_id):


    parse_string = f'in vol {volume_id} with name "*topFace*"'
    surface_ids = list(parse_cubit_list("surface", parse_string))

    if len(surface_ids) > 0 and len(surface_ids)<3: 
        mat_ori_surface = surface_ids[0]
    elif len(surface_ids) >=3:           #Some web_thickness volumes have surfaces with topFace in the name. 
        for surface_id in surface_ids:   # These are filtered out in this loop
            surf_name = cubit.get_entity_name("surface", surface_id)
            if 'thickness' not in surf_name:
                mat_ori_surface = surface_id
                break
    else:
        cubit.cmd(f'save as "Debug.cub" overwrite')
        raise ValueError(f'Could not find mat_ori_surface. \n volume_id: {volume_id}')

    for curve in cubit.volume(volume_id).curves():
        curve_name = cubit.get_entity_name("curve", curve.id())
        if 'thickness' in curve_name:
            if 'web_thickness' not in curve_name:
                thickness_curve = curve
                break


    surface_normal=get_surface_normal(mat_ori_surface)
    midPoint = list(thickness_curve.position_from_fraction(0.5))

    approximate_thickness_direction = thickness_curve.tangent(midPoint)
    if dotProd(surface_normal, approximate_thickness_direction) > 0:
        sign = 1.0
    else:
        sign = -1.0             
    return mat_ori_surface,sign

