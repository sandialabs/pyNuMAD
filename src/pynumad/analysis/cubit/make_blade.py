from pynumad.analysis.cubit.make_cross_sections import *
from pynumad.analysis.cubit.connect_cross_sections import *
from pynumad.utils.orientations import *
import numpy as np
import os
import glob
import pickle
import time
import multiprocessing
def write_path_abscissas_to_file(set_verts,file_name,non_dim_span=[],directory='.'):
    if not os.path.exists(directory):
        os.makedirs(directory)
    file = open(f'{directory}/{file_name}', 'w')
    for set_name in set_verts.keys():
        all_set_coords=[]
        
        for vertex_id in set_verts[set_name]:
            node_id=get_vertex_node(vertex_id)
            coords=get_nodal_coordinates(node_id)
            all_set_coords.append(coords)

        non_dim_path_length=[]
        if 'span' not in set_name:
            all_set_coords=np.array(all_set_coords)
            segment_lengths=list((np.sqrt((all_set_coords[:-1]-all_set_coords[1:])**2)).sum(1))
            segment_lengths.insert(0,0.0)
            path_length=sum(segment_lengths)
            
            if path_length > 0:
                for i_seg in range(len(segment_lengths)):
                    non_dim_path_length.append(sum(segment_lengths[:i_seg+1])/path_length)
            else:
                non_dim_path_length = []
        else:
            non_dim_path_length=non_dim_span

        file.write(f'{set_name} {" ".join(map(str,non_dim_path_length))}\n')
    file.close()

def write_path_coords_to_file(set_verts,prepend,directory='.'):
    if not os.path.exists(directory):
        os.makedirs(directory)
    for set_name in set_verts.keys():
        file = open(f'{directory}/{prepend}{set_name}.coords', 'w')
        for vertex_id in set_verts[set_name]:
            node_id=get_vertex_node(vertex_id)
            coords=get_nodal_coordinates(node_id)
            file.write(f'{" ".join(map(str,coords))}\n')
        file.close()
def write_path_node_ids_to_file(set_verts,file_name,directory='.'):
    if not os.path.exists(directory):
        os.makedirs(directory)

    file = open(f'{directory}/{file_name}', 'w')
    for set_name in set_verts.keys():
        node_ids=[]
        for vertex_id in set_verts[set_name]:
            node_ids.append(get_vertex_node(vertex_id))
        file.write(f'{set_name} {" ".join(map(str,node_ids))}\n')
    file.close()
def write_path_node_angles_to_file(set_verts,prepend,directory='.'):
    if not os.path.exists(directory):
        os.makedirs(directory)

    for set_name in set_verts.keys():
        file = open(f'{directory}/{prepend}{set_name}.dcm', 'w')
        dcms=[ [] for _ in range(9) ]
        for vertex_id in set_verts[set_name]:
            node_id=get_vertex_node(vertex_id)
            parse_string = f'in node {node_id}'
            volume_id = parse_cubit_list("volume", parse_string)[0] #Just use first volume

            coords = get_nodal_coordinates(node_id)
            surf_id_for_mat_ori,sign = get_mat_ori_surface(volume_id)

            surface_normal = vectNorm(
                list(sign*np.array(get_surface_normal_at_coord(surf_id_for_mat_ori, coords))))

            ref_line_direction = [0,0,1]
            #https://www.maplesoft.com/support/help/maple/view.aspx?path=MathApps%2FProjectionOfVectorOntoPlane
            spanwise_direction = vectNorm(np.array(ref_line_direction)-np.dot(ref_line_direction,surface_normal)*np.array(surface_normal))

            perimeter_direction = vectNorm(np.cross(surface_normal, spanwise_direction))

            newCoordinateSystemVectors = [spanwise_direction,perimeter_direction,surface_normal]

            globalAxisBasisVectors = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            
            dcm = getDCM(globalAxisBasisVectors, newCoordinateSystemVectors)

            dcms[0].append(dcm[0,0])
            dcms[1].append(dcm[0,1])
            dcms[2].append(dcm[0,2])
            dcms[3].append(dcm[1,0])
            dcms[4].append(dcm[1,1])
            dcms[5].append(dcm[1,2])
            dcms[6].append(dcm[2,0])
            dcms[7].append(dcm[2,1])
            dcms[8].append(dcm[2,2])


        
        file.write(f'R_11 {" ".join(map(str,dcms[0]))}\n')
        file.write(f'R_12 {" ".join(map(str,dcms[1]))}\n')
        file.write(f'R_13 {" ".join(map(str,dcms[2]))}\n')
        file.write(f'R_21 {" ".join(map(str,dcms[3]))}\n')
        file.write(f'R_22 {" ".join(map(str,dcms[4]))}\n')
        file.write(f'R_23 {" ".join(map(str,dcms[5]))}\n')
        file.write(f'R_31 {" ".join(map(str,dcms[6]))}\n')
        file.write(f'R_32 {" ".join(map(str,dcms[7]))}\n')
        file.write(f'R_33 {" ".join(map(str,dcms[8]))}\n')
        file.close()
def get_hex_orientations_euler(volume_id):
    global_el_ids_in_vol=[]
    theta1s_in_vol=[]
    theta2s_in_vol=[]
    theta3s_in_vol=[]



    surf_id_for_mat_ori,sign = get_mat_ori_surface(volume_id)
    #volume_name = cubit.get_entity_name("volume", volume_id)
    #t0 = time.time()

    for el_id in get_volume_hexes(volume_id):
        coords = cubit.get_center_point("hex", el_id)
            
        surface_normal = vectNorm(
            list(sign*np.array(get_surface_normal_at_coord(surf_id_for_mat_ori, coords))))

        ref_line_direction = [0,0,1]
        #https://www.maplesoft.com/support/help/maple/view.aspx?path=MathApps%2FProjectionOfVectorOntoPlane
        spanwise_direction = vectNorm(np.array(ref_line_direction)-np.dot(ref_line_direction,surface_normal)*np.array(surface_normal))

        perimeter_direction = vectNorm(np.cross(surface_normal, spanwise_direction))

        newCoordinateSystemVectors = [spanwise_direction,perimeter_direction,surface_normal]

        globalAxisBasisVectors = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

        global_id=get_global_element_id('hex',el_id)
        
        dcm = getDCM(globalAxisBasisVectors, newCoordinateSystemVectors)

        temp1, temp2, temp3 = dcmToEulerAngles(dcm)

        global_el_ids_in_vol.append(global_id)
        theta1s_in_vol.append(-1*temp1)
        theta2s_in_vol.append(-1*temp2)
        theta3s_in_vol.append(-1*temp3)


    return global_el_ids_in_vol,theta1s_in_vol,theta2s_in_vol,theta3s_in_vol

def get_hex_orientations_two_points(volume_id):
    global_el_ids_in_vol=[]


    
    spanwise_directions_in_vol = []
    perimeter_directions_in_vol = []

    surf_id_for_mat_ori,sign = get_mat_ori_surface(volume_id)
    #volume_name = cubit.get_entity_name("volume", volume_id)
    #t0 = time.time()

    for el_id in get_volume_hexes(volume_id):
        coords = cubit.get_center_point("hex", el_id)
            
        surface_normal = vectNorm(
            list(sign*np.array(get_surface_normal_at_coord(surf_id_for_mat_ori, coords))))

        ref_line_direction = [0,0,1]
        #https://www.maplesoft.com/support/help/maple/view.aspx?path=MathApps%2FProjectionOfVectorOntoPlane
        spanwise_direction = vectNorm(np.array(ref_line_direction)-np.dot(ref_line_direction,surface_normal)*np.array(surface_normal))

        perimeter_direction = vectNorm(np.cross(surface_normal, spanwise_direction))

        global_id=get_global_element_id('hex',el_id)
        
        global_el_ids_in_vol.append(global_id)

        spanwise_directions_in_vol.append(spanwise_direction)
        perimeter_directions_in_vol.append(perimeter_direction)

    return global_el_ids_in_vol,spanwise_directions_in_vol,perimeter_directions_in_vol

def get_tet_orientations(volume_id):
    global_el_ids_in_vol=[]
    theta1s_in_vol=[]
    theta2s_in_vol=[]
    theta3s_in_vol=[]

    
    spanwise_directions_in_vol = []
    perimeter_directions_in_vol = []

    surf_id_for_mat_ori,sign = get_mat_ori_surface(volume_id)
    #volume_name = cubit.get_entity_name("volume", volume_id)
    #t0 = time.time()

    for el_id in get_volume_tets(volume_id):
        coords = cubit.get_center_point("tet", el_id)
            
        surface_normal = vectNorm(
            list(sign*np.array(get_surface_normal_at_coord(surf_id_for_mat_ori, coords))))

        ref_line_direction = [0,0,1]
        #https://www.maplesoft.com/support/help/maple/view.aspx?path=MathApps%2FProjectionOfVectorOntoPlane
        spanwise_direction = vectNorm(np.array(ref_line_direction)-np.dot(ref_line_direction,surface_normal)*np.array(surface_normal))

        perimeter_direction = vectNorm(np.cross(surface_normal, spanwise_direction))

        newCoordinateSystemVectors = [spanwise_direction,perimeter_direction,surface_normal]

        globalAxisBasisVectors = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

        global_id=get_global_element_id('tet',el_id)
        
        dcm = getDCM(globalAxisBasisVectors, newCoordinateSystemVectors)

        temp1, temp2, temp3 = dcmToEulerAngles(dcm)

        global_el_ids_in_vol.append(global_id)
        theta1s_in_vol.append(-1*temp1)
        theta2s_in_vol.append(-1*temp2)
        theta3s_in_vol.append(-1*temp3)

        spanwise_directions_in_vol.append(spanwise_direction)
        perimeter_directions_in_vol.append(perimeter_direction)

    return global_el_ids_in_vol,theta1s_in_vol,theta2s_in_vol,theta3s_in_vol,spanwise_directions_in_vol,perimeter_directions_in_vol
    
def assign_material_orientations(orientation_data):
    #Apply Material Orientation
    global_ids=orientation_data[0]
    n_el = len(global_ids)
    theta1s=orientation_data[1]
    theta2s=orientation_data[2]
    theta3s=orientation_data[3]

    cubit.set_element_variable(global_ids, 'rotation_angle_one', theta1s)
    cubit.set_element_variable(global_ids, 'rotation_angle_two', theta2s)
    cubit.set_element_variable(global_ids, 'rotation_angle_three', theta3s)

    cubit.set_element_variable(global_ids, 'rotation_axis_one', 1*np.ones(n_el))
    cubit.set_element_variable(global_ids, 'rotation_axis_two', 2*np.ones(n_el))
    cubit.set_element_variable(global_ids, 'rotation_axis_three', 3*np.ones(n_el))

    return

def compute_material_orientations(element_shape,output_format = 'euler',ncpus = 1):
    # # ####################################
    # # ### Assign material orientations ###
    # # ####################################
    
    parse_string = f'with name "*volume*"'
    all_volume_ids = parse_cubit_list("volume", parse_string)
    
    t0 = time.time()
    print(f'Calculating material orientations with {ncpus} CPU(s)...')

    if 'hex' in element_shape:
        if ncpus==1:
            ans = []
            if 'euler' in output_format:
                for vol_id in all_volume_ids:
                    ans.append(get_hex_orientations_euler(vol_id))
            elif 'two_points' in output_format:
                for vol_id in all_volume_ids:
                    ans.append(get_hex_orientations_two_points(vol_id))
            else:
                raise NameError(f'Material Orientation output format: {output_format} is not supported')
        else:
            pool_obj = multiprocessing.Pool(ncpus)
            if 'euler' in output_format:
                ans = pool_obj.map(get_hex_orientations_euler,all_volume_ids)
            elif 'two_points' in output_format:
                ans = pool_obj.map(get_hex_orientations_two_points,all_volume_ids)
            else:
                raise NameError(f'Material Orientation output format: {output_format} is not supported')
            
            pool_obj.close()
    else:
        raise NameError(f' element shape {element_shape} unsupported.')
    t1 = time.time()
    print(f'Total time for material orientations: {t1-t0}')

    ans=np.array(ans,dtype=object)
    global_ids=[]

    if 'euler' in output_format:
        theta1s=[]
        theta2s=[]
        theta3s=[]
        for i in range(len(all_volume_ids)):
            global_ids+=list(ans[i][0])
            theta1s+=list(ans[i][1])
            theta2s+=list(ans[i][2])
            theta3s+=list(ans[i][3])

        return [global_ids,theta1s,theta2s,theta3s]
    
    elif 'two_points' in output_format:
        spanwise_directions = []
        perimiter_directions = []

        for i in range(len(all_volume_ids)):
            global_ids+=list(ans[i][0])
            spanwise_directions+=list(ans[i][1])
            perimiter_directions+=list(ans[i][2])

        return [global_ids,spanwise_directions,perimiter_directions]



def order_path_points(points, ind):
    points_new = [ points.pop(ind) ]  # initialize a new list of points with the known first point
    pcurr      = points_new[-1]       # initialize the current point (as the known point)
    pointer=[ind]
    while len(points)>0:
        d      = np.linalg.norm(np.array(points) - np.array(pcurr), axis=1)  # distances between pcurr and all other remaining points
        ind    = d.argmin()                   # index of the closest point
        points_new.append( points.pop(ind) )  # append the closest point to points_new
        pointer.append(ind)
        pcurr  = points_new[-1]               # update the current point
    return pointer

def get_nodal_coordinates_from_set_of_nodes(nodeset_nodes):
    # nodeset_nodes is a list of ints
    coords=[]
    for node_id in nodeset_nodes:
        coords.append(get_nodal_coordinates(node_id))
    return coords

def get_nodeset_nodes_from_name(set_name):

        nodeset_id = parse_cubit_list('nodeset',set_name)
        if len(nodeset_id)==1:
            return get_nodeset_nodes_inclusive(nodeset_id[0])
        else:
            if len(nodeset_id)==0:
                raise RuntimeError(f'No nodeset named {set_name} found.')
            else:
                raise RuntimeError(f'More than one nodeset with name {set_name} found.')


def get_sideset_nodes_from_name(set_name):
        node_ids=[]
        sideset_id = parse_cubit_list('sideset',set_name)
        if len(sideset_id)==1:
            surface_ids=get_sideset_surfaces(sideset_id[0])
            for surface_id in surface_ids:
                node_ids+=list(get_surface_nodes(surface_id))
            return node_ids
        else:
            if len(sideset_id)==0:
                raise RuntimeError(f'No sideset named {set_name} found.')
            else:
                raise RuntimeError(f'More than one sideset with name {set_name} found.')
        

def sweep_volumes(vol_to_mesh):
    failed_volumes=[]
    for volume_id in vol_to_mesh:

        cubit.cmd(f"mesh vol {volume_id}")
        if not is_meshed('volume',volume_id): #Try and fix unmeshed volumes by explicitly sweeping
            volume_name=cubit.get_entity_name("volume", volume_id)
            source_target_side=[]
            source_target_side_is_meshed=[]

            for surface in cubit.volume(volume_id).surfaces():
                surface_name=cubit.get_entity_name("surface", surface.id())
                if volume_name.split('_')[0] == surface_name.split('_')[0]:
                    if 'Face' in surface_name:
                        source_target_side.append(surface.id())
                        if is_meshed('surface',surface.id()):
                            source_target_side_is_meshed.append(True)
                        else:
                            source_target_side_is_meshed.append(False)

            if len(source_target_side)==2:
                if source_target_side_is_meshed[0]:
                    source_side=source_target_side[0]
                    target_side=source_target_side[1]
                elif source_target_side_is_meshed[1]:
                    source_side=source_target_side[1]
                    target_side=source_target_side[0]
                else:
                    source_side=source_target_side[0]
                    target_side=source_target_side[1]

                cubit.cmd(f"volume {volume_id} scheme Sweep source surface {source_side} target surface {target_side} sweep transform least squares")
                cubit.cmd(f"mesh vol {volume_id}")
        if not is_meshed('volume',volume_id):
            failed_volumes.append(volume_id)
    return failed_volumes
def debug():
    #cubit.cmd(f"delete curve 1")
    cubit.cmd(f'save as "Debug.cub" overwrite')
def cubit_make_cross_sections(blade,wt_name,settings,cs_params,model2Dor3D,stationList=None,directory="."):
    
    """Directs Cubit to make cross sections from blade object data.

    Parameters
    ----------
    blade : blade object
        pyNuMAD blade object
    wt_name : str
        Used to name any files that are generated.
    settings : dict
        _description_
    cs_params : dict
        _description_
    model2Dor3D : str
        Users should set this '2d'. Functions such as cubit_make_solid_blade set this to '3d'.
    stationList : list, optional
        Integer list of stations user wants cross sections. By default None or empty list makes all the statations.
    directory : str
        Name of the directory to store all generated files.

    Returns
    -------
    cubit: cubit object
        cubit session data
    blade: blade object
        returns the modified blade object
    surface_dict: dict
        Keys are integers for the Cubit surface IDs for the cross sections. Each surface has
        it's own dictionary with the following keys: 'curves', 'verts', 'material_name', 'ply_angle'.
        e.g. 
        >> surface_dict[9]
        >>     {'curves': [131, 164, 129, 163], 'verts': [500, 501, 497, 496], 'material_name': 'glass_triax', 'ply_angle': 0}

    birds_mouth_verts: tuple
        Used internally.
    i_station_first_web: int
        Used internally.
    i_station_last_web: int
        Used internally. 
    materials_used: set
        Used for in FEA input file generation to define unique materials.
    spanwise_mat_ori_curve: int
        Cubit curve ID for the main spanwise spline corresponding to the curvilinear blade axis.


    Raises
    ------
    ValueError
       'Presweep is untested for cross-sectional meshing'
    ValueError
        'ANBA currently not supported'
    NameError
       f'Unknown beam cross-sectional solver: {settings["make_input_for"]}'
    NameError
        f'Unknown model export format: {settings["export"]}'
    """
    geometry = blade.geometry
    stackdb = blade.stackdb
    definition = blade.definition
    keypoints = blade.keypoints

    if stationList is None or len(stationList) == 0:
        stationList = list(range(len(geometry.ispan)))

    # Initialize variables
    surface_dict = {}
    # Uniquly track which materiall IDs are actuall used in blade model
    materials_used = set()
    iLE = geometry.LEindex + 1
    thickness_scaling = 0.001
    geometry_scaling = thickness_scaling * 1000

    # Set up Cubit
    cubit.init(["cubit", "-nojournal"])

    cubit.cmd("undo off")
    #cubit.cmd("set geometry accuracy 1e-6")
    # making numerus 3D volumes is very slow with autosize on
    cubit.cmd("set default autosize off")

    # Modify blade object to accomodate actual layer thicknesses

    expandTEthicknesses = list(
        cs_params["te_adhesive_thickness"]
        + 6 * cs_params["minimum_layer_thickness"]
    )
    blade.expand_blade_geometry_te(expandTEthicknesses)

    stackdb.edit_stacks_for_solid_mesh()


    all_layer_thicknesses = [[] for _ in range(3)] #three modeled layers
    _,n_stations = np.shape(stackdb.stacks)
    for i_station in range(n_stations):
        temp = stackdb.stacks[:, i_station]
        temp = np.flip(temp)

        stacks = list(stackdb.stacks[1:6, i_station]) + list(temp[1:6])

        for stack in stacks:
            for i_layer,layer_thicknesses in enumerate(all_layer_thicknesses):
                layer_thicknesses.append(stack.layer_thicknesses()[i_layer])

    if 'nel_per_core_layer' not in cs_params.keys() or cs_params['nel_per_core_layer'] < 1.0:
        #find number of elemements through core layer:
        nel_per_core_layer=[]
        for i in range(len(all_layer_thicknesses[0])):
            layer_t1=all_layer_thicknesses[0][i]
            layer_t2=all_layer_thicknesses[1][i]
            layer_t3=all_layer_thicknesses[2][i]

            element_t1= layer_t1/cs_params['nel_per_layer']
            element_t3= layer_t3/cs_params['nel_per_layer']

            nel_per_core_layer.append(mean([layer_t2/element_t1/cs_params['element_thickness_ar'],layer_t2/element_t3/cs_params['element_thickness_ar']]))
        
        cs_params['nel_per_core_layer']=int(max(round(mean(nel_per_core_layer)),1))

    

    hasWebs = []
    webNumber = 1
    for i_station in range(len(stackdb.swstacks[webNumber])):
        if not len(stackdb.swstacks[webNumber][i_station].plygroups) == 0:
            hasWebs.append(True)
        else:
            hasWebs.append(False)

    # WARNING - Last station never has webs. Fix later
    hasWebs.append(False)
    # WARNING - Last station never has webs. Fix later

    # Create Referece line as a spline

    ref_line_coords = np.vstack(
        ([definition.sweep, definition.prebend, definition.ispan])
    ).transpose()
    spanwise_mat_ori_curve = 1

    #if model2Dor3D.lower() == "3d":
        #write_spline_from_coordinate_points(cubit, ref_line_coords)
    #else: do it in cross section loop
    
    #Get last round station index
    # is_station_flatback = []
    # for i_station in range(len(blade.geometry.ispan)):
    #     if geometry.get_profile_te_type(i_station) == "flat":
    #         is_station_flatback.append(True)
    #     else:
    #         is_station_flatback.append(False)
    
    # is_station_flatback.append(True) #last station is never round
    # last_round_station=next((i-1 for i, x in enumerate(is_station_flatback) if x), None)





########################
    #### Step one create outer mold line
    excess_lengths=[]
    te_angles=[]
    for i_station_geometry in range(len(blade.geometry.ispan)-1): #-1 b/c fewer stacks than stations
        xyz = get_blade_geometry_for_station(blade, i_station_geometry) * geometry_scaling
        
        npts=5
        # Start indexing from 1 (not 0) to ignore first point: because first point is not on the LP or HP surface but rather is the midpoint at the TE
        splinePoints = xyz[1:npts, :]
        write_spline_from_coordinate_points(cubit, splinePoints)
        hp_key_curve = get_last_id("curve")

        xyz = np.flip(xyz, 0)
        splinePoints = xyz[1:npts, :]
        write_spline_from_coordinate_points(cubit, splinePoints)
        lp_key_curve = get_last_id("curve")


        first_point = xyz[-2, :]
        second_point = xyz[1, :]

        flatback_length=np.linalg.norm(second_point - first_point)

        athickness=cs_params["te_adhesive_thickness"][i_station_geometry]
        stack_thicknesses_hp=sum(stackdb.stacks[1, i_station_geometry].layer_thicknesses())/1000
        stack_thicknesses_lp=sum(stackdb.stacks[-2, i_station_geometry].layer_thicknesses())/1000

        excess_lengths.append(flatback_length-(stack_thicknesses_lp+stack_thicknesses_hp+athickness))



        curve_fraction = 0
        te_angles.append(get_te_angle(hp_key_curve, lp_key_curve, curve_fraction))
        # print(f"station {i_station}")
        # print(f"edgeLength={flatback_length*1000}")
        # print(cs_params)
        # print(f'athickness={cs_params["te_adhesive_thickness"][i_station]*1000}')
        # print(f'te_adhesive_width {cs_params["te_adhesive_width"][i_station]*1000}')
        # print(f"te_angle {te_angle}")

    last_round_station=next((i-1 for i, x in enumerate(te_angles) if x < 50.0), None)
    last_flat_station=next((i-1 for i, x in enumerate(te_angles) if x < 10.0), None)
    
    if last_round_station == None and last_round_station == None:
        last_round_station = len(te_angles)+1000 #Arbitrarily large
        last_flat_station = last_round_station+1 #Arbitrarily large
    else:

        last_10deg_station=last_flat_station
        for i_length, excess_length in enumerate(excess_lengths[last_10deg_station+1:]):
            athickness=cs_params["te_adhesive_thickness"][last_10deg_station+1+i_length]
            print(f'i {last_10deg_station+1+i_length},excess_length {excess_length*1000}, athickness{athickness*1000}')
            if (excess_length-athickness)/excess_length > 0.025:
                last_flat_station=last_flat_station+1+i_length
            else:
                break

########################





    with open(f"{wt_name}.log", "w") as logFile:
        logFile.write(f"Making cross sections for {wt_name}\n")

    path_name = directory + "/" + wt_name + "-crossSections"
    birds_mouth_verts = []



    for i_station in stationList:
        if model2Dor3D.lower() == "2d":
            cubit.cmd(
                "reset "
            )  # This is needed to restart node numbering for VABS. VABS neeeds every element and node starting from 1 to nelem/nnode should be present
        write_spline_from_coordinate_points(cubit, ref_line_coords)
        i_station_geometry = i_station
        if i_station == len(geometry.ispan) - 1:  # Only do this for the last station
            blade.add_interpolated_station(geometry.ispan[-1] * 0.999)
            stackdb.edit_stacks_for_solid_mesh()
            expandTEthicknesses.append(expandTEthicknesses[-1])
            blade.expand_blade_geometry_te(expandTEthicknesses)

            # adjustLastStackAfterNewTipStation(i_station)

            i_station_geometry = i_station + 1
        
        #is_flatback=is_station_flatback[i_station_geometry]


        i_station_first_web = np.argwhere(hasWebs)[0][0]
        i_station_last_web = np.argwhere(hasWebs)[-1][0]

        if hasWebs[i_station] == True:
            webNumber = 1
            aft_web_stack = stackdb.swstacks[webNumber][i_station]
            webNumber = 0
            fore_web_stack = stackdb.swstacks[webNumber][i_station]
        else:
            if i_station < i_station_first_web:
                iWebStation = i_station_first_web

            #         elif i_station_last_web == len(blade.ispan) - 1-1:
            else:
                iWebStation = i_station_last_web
            #         else:
            #             raise Exception('assuming web ends at last station for now. ')

            webNumber = 1
            aft_web_stack = stackdb.swstacks[webNumber][iWebStation]
            webNumber = 0
            fore_web_stack = stackdb.swstacks[webNumber][iWebStation]

        cs_normal = get_cs_normal_vector(
            np.array(
                [
                    keypoints.key_points[2, :, i_station_geometry],
                    keypoints.key_points[3, :, i_station_geometry],
                    keypoints.key_points[7, :, i_station_geometry],
                ]
            )
        )
        if model2Dor3D.lower() == "3d":
            make_webs = True
        else: 
            make_webs = hasWebs[i_station]

        # Only save birds_mouth_verts for the right cross-section
        if i_station == i_station_first_web:
            birds_mouth_verts = make_a_cross_section(wt_name,
                surface_dict,
                i_station,
                i_station_geometry,
                blade,
                make_webs,
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
            )
        else:
            make_a_cross_section(wt_name,
                surface_dict,
                i_station,
                i_station_geometry,
                blade,
                make_webs,
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
            )
            

        cubit.cmd(f"delete curve all with Is_Free except {spanwise_mat_ori_curve}")

        # Chord line for rotation of cross-section for homogenization
        if model2Dor3D.lower() == "2d":
            #         #Blocks
            if 'd_tube' in cs_params.keys() and cs_params['d_tube']:
                keep_list=[]

                cubit.cmd(f'delete surface with x_coord < 0"')
                cubit.cmd(f'delete surface with name "*layer9*"')
                cubit.cmd(f'delete surface with name "*layer10*"')
                cubit.cmd(f'delete surface with name "*layer11*"')

                delete_list=[]
                parse_string = f'with name "*layer3*"'
                delete_list += list(parse_cubit_list("surface", parse_string))
                parse_string = f'with name "*layer4*"'
                delete_list += list(parse_cubit_list("surface", parse_string))

                keep_list=[]
                #LE
                for i in [121,122,123]:
                    parse_string = f'with name "shell*Station*surface{i}"'
                    keep_list += list(parse_cubit_list("surface", parse_string))

                #Web
                for i in [1,2,3,4,5,6,17,18,19,20,21,22]:
                    parse_string = f'with name "web_web*surface{i}"'
                    keep_list += list(parse_cubit_list("surface", parse_string))

                vol_ids=set(delete_list).difference(set(keep_list))

                cubit.cmd(f'delete vol {l2s(vol_ids)}')
                cubit.cmd(f'delete vol with name "*Station005*"')

            for imat, material_name in enumerate(materials_used):
                cubit.cmd(f'block {imat+1} add surface with name "*{material_name}*"')
                cubit.cmd(f'block {imat+1} name "{material_name}"')

            addColor(blade, "surface")

            # create_vertex(blade.geometry[0,0,i_station]*geometry_scaling,blade.geometry[0,1,i_station]*geometry_scaling,blade.geometry[0,2,i_station]*geometry_scaling)
            # TEvert=get_last_id("vertex")
            # create_vertex(blade.geometry[iLE-1,0,i_station]*geometry_scaling,blade.geometry[iLE-1,1,i_station]*geometry_scaling,blade.geometry[iLE-1,2,i_station]*geometry_scaling)
            # LEvert=get_last_id("vertex")

            # cubit.cmd(f'create curve vertex {TEvert} {LEvert}')
            # coords=cubit.vertex(TEvert).coordinates()
            # tangent=cubit.curve(get_last_id("curve")).tangent(coords)
            # tangent_direction=vectNorm(list(tangent))  #Unit vector of tangent.
            # crossSectionRotationAngle=math.atan2(tangent_direction[1],tangent_direction[0])*180/pi

            parse_string = f'with name "*Station{str(i_station).zfill(3)}*"'
            volume_ids = parse_cubit_list("surface", parse_string)

            # Undo initial twist
            cubit.cmd(
                f"rotate Surface {l2s(volume_ids)} angle {definition.degreestwist[i_station]} about Z include_merged "
            )

            # Undo prebend
            if definition.prebend[i_station] != 0:
                cubit.cmd(f"move surface {l2s(volume_ids)} y {-1*definition.prebend[i_station]} include_merged")

            # Undo sweep
            if definition.sweep[i_station] != 0:
                raise ValueError("Presweep is untested for cross-sectional meshing")

            if 'ref_line_type' in cs_params and 'centroid' in cs_params['ref_line_type'].lower():
                centroidal_vert_id, centroidal_ref_line_coords=get_locus_of_cross_sectional_centroids([i_station])
                cubit.cmd(f"move surface {l2s(volume_ids)} x {-1*centroidal_ref_line_coords[0][0]} include_merged")
                cubit.cmd(f"move surface {l2s(volume_ids)} y {-1*centroidal_ref_line_coords[0][1]} include_merged")
                #cubit.cmd(f"move surface {l2s(volume_ids)} z {-1*centroidal_ref_line_coords[0][2]} include_merged")
                
                centroidal_vert_id, centroidal_ref_line_coords=get_locus_of_cross_sectional_centroids([i_station])


            # Mesh the cross-section
            cubit.cmd(f'curve with name "face_thickness*" interval {cs_params["nel_per_layer"]}')
            cubit.cmd(f'curve with name "*face_web_thickness*" interval {cs_params["nel_per_layer"]}')

            cubit.cmd(f'curve with name "core_thickness*" interval {cs_params["nel_per_core_layer"]}')
            cubit.cmd(f'curve with name "*core_web_thickness*" interval {cs_params["nel_per_core_layer"]}')


            cubit.cmd(f'curve with name "*hoop*" in surface with name "roundTEadhesive*" interval {cs_params["nel_per_layer"]}')

            # cubit.cmd(f'imprint volume {l2s(surface_ids)}')
            cubit.cmd(f"merge volume {l2s(volume_ids)}")
            cubit.cmd(f"set default autosize on")

            if cs_params["element_shape"].lower() == "tri":
                cubit.cmd(f"surface {l2s(volume_ids)} scheme tri")
            else:
                cubit.cmd(f"surface {l2s(volume_ids)} scheme map")


            t_1=get_mean_layer_thickness_at_station(i_station) #Find mean layer thickness for first station
            e_size_1=t_1/cs_params['nel_per_layer']*cs_params['element_ar'] #Get spanwise element size at station_id cross section
            
            cubit.cmd(f"surface all size {e_size_1}")
            cubit.cmd(f"mesh surface {l2s(volume_ids)}")

            file_name = wt_name + "-" + str(i_station) + "-t-0.in"

            if not os.path.exists(directory):
                os.makedirs(directory)

            if get_mesh_error_count() ==0:
                if settings["make_input_for"] is not None:
                    if "vabs" in settings["make_input_for"].lower():
                        write_vabs_input(
                            surface_dict,
                            blade,
                            cs_params,
                            directory,
                            file_name,
                            volume_ids,
                            materials_used,
                            cs_normal,
                        )

                elif "anba" in settings["make_input_for"].lower():
                    raise ValueError("ANBA currently not supported")
                else:
                    raise NameError(
                        f'Unknown beam cross-sectional solver: {settings["make_input_for"]}'
                    )
            else: 
                with open(f"{wt_name}.log", "a") as logFile:
                    logFile.write(f"    Warning: {get_mesh_error_count()} cross section mesh errors exist in station {i_station}\n")
    
            if 'd_tube' in cs_params.keys() and cs_params['d_tube']:

                set_verts={}
                set_verts[f'thickness_{str(i_station).zfill(3)}_s1']=[2804, 9817, 9823, 9831]
                set_verts[f'thickness_{str(i_station).zfill(3)}_s2']=[9985, 10001, 10039, 10083]
                set_verts[f'spanwise_s2']=[9985]
                set_verts[f'circumferential_{str(i_station).zfill(3)}']=[9650,2801,2802,2804,2806,2808,2810,2812,2814,9979,9983,9985,10117,10114,10115,6396,6395,6398,6266,6264,6260,2764,2762,2760,2758,2756,2754,2752,2751,5931,5951,5989,6033,11530,11538,11714,11706,9752,9708,9670,9650]
                set_verts[f'p1']=[9650]
                set_verts[f'p2']=[5931]

                file_name=f'beam_{str(i_station).zfill(3)}.nodes'
                write_path_node_ids_to_file(set_verts,file_name,directory)
                #write_path_coords_to_file(set_verts,prepend,dir_name)
                
                file_name=f'{directory}/beam_{str(i_station).zfill(3)}.abscissa'
                write_path_abscissas_to_file(set_verts,file_name)   
                

                # nodeset_id= cubit.get_next_nodeset_id()
                # cubit.cmd(f'nodeset {nodeset_id} add curve 925 928 932')
                # cubit.cmd(f'nodeset {nodeset_id} name "{node_set_name}"')

                # path_type='thickness'
                #node_order=get_path_node_order(node_set_name,path_type)
                #def get_path_node_order(node_set_name,path_type):


                # nodeset_nodes = get_nodeset_nodes_from_name(node_set_name)
                # coords=get_nodal_coordinates_from_set_of_nodes(nodeset_nodes)

                # pointer=order_path_points(coords, ind)



                # file = open(directory +'/'+ axisFileName, 'w')
                # file.write('--------- BEAMDYN with OpenFAST INPUT FILE -------------------------------------------\n')
                
                # nodeset_id= cubit.get_next_nodeset_id()
                # cubit.cmd(f'nodeset {nodeset_id} add curve 1040 1065 1091')
                # cubit.cmd(f'nodeset {nodeset_id} name "s2_thickness_path"')


                # nodeset_id= cubit.get_next_nodeset_id()
                # cubit.cmd(f'nodeset {nodeset_id} add vertex 9985')
                # cubit.cmd(f'nodeset {nodeset_id} name "spanwise_path"')


                # nodeset_id= cubit.get_next_nodeset_id()
                # cubit.cmd(f'nodeset {nodeset_id} add curve 73 74 75 76 77 78 79 92 93 94 95 96 97 98 335 361 476 478 479 556 557 558 824 848 873 899 1014 1016 1017 1094 1095 1096 1196 1220 1224 1324 1328 1422')
                # cubit.cmd(f'nodeset {nodeset_id} name "circumferential_path"')
                


                # nodeset_id= cubit.get_next_nodeset_id()
                # cubit.cmd(f'nodeset {nodeset_id} add curve with name "*oml*"')
                # cubit.cmd(f'nodeset {nodeset_id} name "{node_set_name}"')
                # oml_nodes = get_nodeset_nodes_from_name(node_set_name)


            if settings["export"] is not None:
                if (
                    "g" in settings["export"].lower()
                    or "cub" in settings["export"].lower()
                ):
                    if "g" in settings["export"].lower():
                        cubit.cmd(f'export mesh "{path_name}-{str(i_station)}.g" overwrite')
                    if "cub" in settings["export"].lower():
                        cubit.cmd(f"delete curve {spanwise_mat_ori_curve}")
                        cubit.cmd(f'save as "{path_name}-{str(i_station)}.cub" overwrite')
                        print('')
                elif len(settings["export"]) == 0:
                    pass
                else:
                    raise NameError(
                        f'Unknown model export format: {settings["export"]}'
                    )



    # Import all cross-sections into one cub file
    if model2Dor3D.lower() == "2d" and settings["export"] is not None and "cub" in settings["export"].lower():
        cubit.cmd("reset ")
        
        #Since cross sections were translated for cross sectional codes, remove prebend and sweep from ref axis.
        ref_line_coords[:,0]=np.zeros(len(ref_line_coords[:,0]))
        ref_line_coords[:,1]=np.zeros(len(ref_line_coords[:,0]))
        write_spline_from_coordinate_points(cubit, ref_line_coords)

        for i_station in stationList:
            cubit.cmd(f'import cubit "{path_name}-{str(i_station)}.cub"')
        addColor(blade, "surface")
        cubit.cmd(f"delete vertex with Is_Free")
        cubit.cmd(f'save as "{path_name}.cub" overwrite')

        # Remove unnecessary files to save space
        # for filePath in glob.glob(f"{path_name}-*.cub"):
        #     os.remove(filePath)
    return (cubit,blade,surface_dict,birds_mouth_verts,i_station_first_web,i_station_last_web,materials_used,spanwise_mat_ori_curve,hasWebs)


def cubit_make_solid_blade(
    blade, wt_name, settings, cs_params, stationList=None
):
    """_summary_

    Parameters
    ----------
    blade : blade object
        pyNuMAD blade object
    wt_name : str
        Used to name any files that are generated.
    settings : dict
        _description_
    cs_params : dict
        _description_
    stationList : list, optional
        Integer list of stations user wants cross sections. By default None or empty list makes all the statations.

    Returns
    -------
    materials_used: set
        Used for in FEA input file generation to define unique materials.

    Raises
    ------
    ValueError
        "Need more than one cross section to make a solid model"
    """
    if stationList is None or len(stationList) == 0:
        stationList = list(range(len(blade.ispan)))
    elif len(stationList) == 1:
        raise ValueError("Need more than one cross section to make a solid model")

    (cubit,blade,surface_dict,birds_mouth_verts,i_station_first_web,
     i_station_last_web,materials_used,spanwise_mat_ori_curve,hasWebs) = cubit_make_cross_sections(
        blade, wt_name, settings, cs_params, "3D", stationList)


    i_station_start = stationList[0]
    i_station_end = stationList[-1]
    ### ### ### ###
    # Make volumes along the span.
    ### ### ### ###
    mesh_vol_list = []

    part_name = "shell"
    ordered_list = get_ordered_list(part_name)
    spanwise_splines=[]
    if len(ordered_list) > 0:
        shell_vol_list,spanwise_splines = make_all_volumes_for_a_part(surface_dict, ordered_list, i_station_end,spanwise_splines)
    else:
        shell_vol_list=[]

    part_name = "web"
    ordered_list = get_ordered_list(part_name)
    ordered_list_web = ordered_list.copy()
    if ordered_list and len(ordered_list[0]) > 1:
        web_vol_list,spanwise_splines = make_all_volumes_for_a_part(surface_dict, ordered_list, i_station_end,spanwise_splines)
    else:
        web_vol_list=[]

    part_name = "roundTEadhesive"
    ordered_list = get_ordered_list(part_name)
    if ordered_list and len(ordered_list[0]) > 1:
        roundTEadhesive_vol_list,spanwise_splines = make_all_volumes_for_a_part(surface_dict, ordered_list, i_station_end,spanwise_splines)
    else:
        roundTEadhesive_vol_list=[]


    part_name = "flatTEadhesive"
    ordered_list = get_ordered_list(part_name)

    if ordered_list and len(ordered_list[0]) > 1:
        flatTEadhesive_vol_list,spanwise_splines = make_all_volumes_for_a_part(surface_dict, ordered_list, i_station_end,spanwise_splines)
    else:
        flatTEadhesive_vol_list=[]

    if (
        ordered_list_web
        and len(ordered_list_web[0]) > 1
        and cs_params["birds_mouth_amplitude_fraction"]
        and birds_mouth_verts
    ):
        web_vol_list=make_birds_mouth(
            blade,
            birds_mouth_verts,
            cs_params["birds_mouth_amplitude_fraction"],
            i_station_first_web,
            i_station_last_web,
        )
    mesh_vol_list=shell_vol_list+web_vol_list+roundTEadhesive_vol_list+flatTEadhesive_vol_list
    
    # cubit.cmd(f"merge tol 1e-3")
    cubit.cmd(f"delete surface with Is_Free")

    for i_station in stationList[0:-1]: 
    #for i_keep,keep_web in enumerate(hasWebs):
        if not hasWebs[i_station]:
            cubit.cmd(f"delete volume with name 'web*tation{str(i_station).zfill(3)}*'")



    cubit.cmd(f"merge volume all")
    

    addColor(blade, "volume")


# Mesh sizing
    if 'tet' in cs_params['element_shape']:
        cubit.cmd("set default autosize on")
        cubit.cmd("volume all scheme tetmesh")
        cubit.cmd("set trimesher geometry sizing off")
        cubit.cmd("volume all size auto factor 10")
        cubit.cmd("mesh vol all")
    elif float(cs_params['element_ar']) != 0.0:
        cubit.cmd("set default autosize on")
        omit_surf_mesh=[]
        hplp_max_stack_ct = 14 #The number of stacks in HP surface. Hard code for now.

        #cubit.cmd(f'surface with name "*shell*Station*layer0_bottomFace*" scheme map')
        cubit.cmd(f"surf with name '*Station*surface*' scheme map") #All cross sections are map

        ### Hoop direcrtion mesh spacing for every cross section
        

        #Find transition from round to flatback adhesive
        flatback_adhesive_station_list=[]
        round_adhesive_station_list=[]
        for i_station, station_id in enumerate(stationList):
            i_stack=0
            parse_string = f'with name "*flatTEadhesiveStation{str(station_id).zfill(3)}*surface*"'
            surface_ids = parse_cubit_list("surface", parse_string)
            if len(surface_ids)>0:
                flatback_adhesive_station_list.append(station_id)
            else:
                round_adhesive_station_list.append(station_id)

        e_size=[]
        remove_thickness_curves=[]
        #hoop spacing for flatback_adhesive_stations, if any
        for i_station, station_id in enumerate(flatback_adhesive_station_list):
            t_1=get_mean_layer_thickness_at_station(station_id) #Find mean layer thickness for first station
            e_size_1=t_1/cs_params['nel_per_layer']*cs_params['element_ar'] #Get spanwise element size at station_id cross section
            
            e_size.append(e_size_1)

            for i_stack in range(2*hplp_max_stack_ct):

                if i_stack == 0 or i_stack == hplp_max_stack_ct:

                    parse_string = f'with name "*shell*{str(station_id).zfill(3)}_stack{str(i_stack).zfill(3)}_*surface*"' 
                    surface_ids = parse_cubit_list("surface", parse_string)

                    parse_string = f'with name "*hoop_direction*" in surface {l2s(surface_ids)}'
                    curve_ids = parse_cubit_list("curve", parse_string)


                    tangents = []
                    for curve_id in curve_ids:
                        c1=cubit.curve(curve_id)
                        #list(c1.position_from_fraction(0.5))
                        tangents.append(c1.tangent(list(c1.position_from_fraction(0.5))))

                    mean_tangent = np.mean(tangents, 0)
                    curve_ids = set()
                    for surface_id in surface_ids:
                        for curve in cubit.surface(surface_id).curves():
                            print(np.dot(mean_tangent,curve.tangent(list(curve.position_from_fraction(0.5)))))
                            if round(np.dot(mean_tangent,curve.tangent(list(curve.position_from_fraction(0.5)))),1)>0.925:
                                curve_ids.add(curve.id())
                    curve_ids=list(curve_ids)
                    
                    parse_string = f'with name "*thickness*" in curve {l2s(curve_ids)}'
                    temp_ids = parse_cubit_list("curve", parse_string)
                    remove_thickness_curves+=temp_ids
                else:
                    parse_string = f'with name "hoop_direction{str(station_id).zfill(3)}_stack{str(i_stack).zfill(3)}*"' 
                    curve_ids = parse_cubit_list("curve", parse_string)


                if len(curve_ids) != 4:
                    raise ValueError('Expecting 4 curves for interval assignment.')

                if i_stack == 12 or i_stack == 26:  #Transition from thick core to LE reinf. for HP (i_stack=12) and LP (i_stack=26)
                    cubit.cmd(f"curve {l2s(curve_ids)} interval {cs_params['nel_per_layer']}") 
                else:
                    cubit.cmd(f"curve {curve_ids[0]} scheme bias fine size {e_size_1} coarse size {e_size_1} ")  ### SAME FOR NOW
                current_hoop_interval=get_mesh_intervals("curve" , curve_ids[0])

                #Adjust intervals if needed
                if i_station > 0:
                    #get prior stack interval count
                    parse_string = f'with name "hoop_direction{str(station_id-1).zfill(3)}_stack{str(i_stack).zfill(3)}*"' 
                    temp_curve_ids = parse_cubit_list("curve", parse_string)
                    prior_hoop_interval=get_mesh_intervals("curve" , temp_curve_ids[0])

                    if abs(current_hoop_interval-prior_hoop_interval)<2:
                        current_hoop_interval=prior_hoop_interval
                    
                    #Make sure that interval count is even in surface loop. Otherwise meshing is impossible. 
                    elif (prior_hoop_interval-current_hoop_interval) % 2:
                        current_hoop_interval+=1

                if i_stack == hplp_max_stack_ct and station_id != flatback_adhesive_station_list[-1]:
                    parse_string = f'with name "shell*Station{str(station_id).zfill(3)}_stack{str(i_stack).zfill(3)}_layer0_bottomFace"' 
                    omit_surf_mesh.append(parse_cubit_list("surface", parse_string)[0])

                cubit.cmd(f"curve {l2s(curve_ids)} interval {current_hoop_interval}")  ### SAME FOR NOW

                if i_stack == 0:
                    save_interval =  current_hoop_interval #Make sure TE intervals match
                elif i_stack == hplp_max_stack_ct:

                    cubit.cmd(f"curve {l2s(curve_ids)} interval {save_interval}")  


        #hoop spacing for round_adhesive_station_list, if any
        if len(round_adhesive_station_list)>1:
            round_adhesive_station_list.reverse()
        for i_station, station_id in enumerate(round_adhesive_station_list):

            
            t_1=get_mean_layer_thickness_at_station(station_id) #Find mean layer thickness for first station
            e_size_1=t_1/cs_params['nel_per_layer']*cs_params['element_ar'] #Get spanwise element size at station_id cross section
            
            e_size.insert(0,e_size_1)

            for i_stack in range(2*hplp_max_stack_ct):

                parse_string = f'with name "hoop_direction{str(station_id).zfill(3)}_stack{str(i_stack).zfill(3)}*"' 
                curve_ids = parse_cubit_list("curve", parse_string)
                if i_stack == 12 or i_stack == 26:
                    cubit.cmd(f"curve {l2s(curve_ids)} interval {cs_params['nel_per_layer']}") 
                else:
                    cubit.cmd(f"curve {curve_ids[0]} scheme bias fine size {e_size_1} coarse size {e_size_1} ")  ### SAME FOR NOW
                current_hoop_interval=get_mesh_intervals("curve" , curve_ids[0])

                #Adjust intervals if needed
                if i_station > 0 or len(flatback_adhesive_station_list)>0:
                    #get prior stack interval count
                    parse_string = f'with name "hoop_direction{str(station_id+1).zfill(3)}_stack{str(i_stack).zfill(3)}*"' 
                    temp_curve_ids = parse_cubit_list("curve", parse_string)
                    prior_hoop_interval=get_mesh_intervals("curve" , temp_curve_ids[0])

                    if abs(current_hoop_interval-prior_hoop_interval)<2:
                        current_hoop_interval=prior_hoop_interval
                    
                    #Make sure that interval count is even in surface loop. Otherwise meshing is impossible. 
                    elif (prior_hoop_interval-current_hoop_interval) % 2:
                        current_hoop_interval+=1

                cubit.cmd(f"curve {l2s(curve_ids)} interval {current_hoop_interval}")  ### SAME FOR NOW




        #Set each surface mesh scheme
        parse_string = f"shell*Station*_stack*_layer0_bottomFace"
        oml_to_mesh = parse_cubit_list("surf", parse_string)
        for surface_id in oml_to_mesh:
            curves=cubit.surface(surface_id).curves()
            hoop_intervals=[]
            if len(curves) == 4:
                for curve in curves:
                    curve_id=curve.id()
                    curve_name = cubit.get_entity_name("curve", curve_id)
                    if 'hoop_dir' in curve_name:
                        hoop_intervals.append(get_mesh_intervals("curve" , curve_id))


                if len(hoop_intervals)==2:
                    if hoop_intervals[0] == hoop_intervals[1]:  #Only pave the surfaces that actually need it. Otherwise, the mesh will be unnecessarily goofy. 
                        cubit.cmd(f'surface {surface_id} scheme map')
                    else:
                        cubit.cmd(f'surface {surface_id} scheme pave')
                else:
                    raise ValueError(f'Found {len(hoop_intervals)} hoop curves in surface {surface_id}. Expecting 2')
            else:
                raise ValueError(f'Found {len(curves)} curves in surface {surface_id}. Expecting 4')
        
        #Spanwise spacing for every segment then mesh OML surface
        for iLoop, station_id in enumerate(stationList[:-1]): #Skip the last station since there are n_stations-1 spanwise volumes
            
            e_size_1=e_size[iLoop]
            e_size_2=e_size[iLoop+1]
            # t_2=get_mean_layer_thickness_at_station(stationList[iLoop+1]) #Find mean layer thickness for station station_id+1
            # e_size_2=t_2/cs_params['nel_per_layer']*cs_params['element_ar'] #Get spanwise element size at station_id+1 cross section

            parse_string = f'with name "shell*Station{str(station_id).zfill(3)}*layer0*bottomFace*"' #WOrks for TE only (use later)
            oml_to_mesh = parse_cubit_list("surf", parse_string)
            
            surface_id=oml_to_mesh[0]
            curves=cubit.surface(surface_id).curves()
            if len(curves) == 4:
                for curve in curves:
                    curve_id=curve.id()
                    curve_name = cubit.get_entity_name("curve", curve_id)

                    if 'span_dir' in curve_name:
                        cubit.cmd(f"curve {curve_id} scheme bias fine size {e_size_1} coarse size {e_size_2} ") 
                        span_interval=get_mesh_intervals("curve" , curve_id) 
                        break
            else:
                raise ValueError(f'Found {len(curves)} curves in {volume_id}. Expecting 4')
            
            cubit.cmd(f"curve with name 'span_dir{str(station_id).zfill(3)}*' interval {span_interval}") 
            

            cubit.cmd(f"mesh surface {l2s(oml_to_mesh)} except surf {l2s(omit_surf_mesh)}")


        parse_string = f'with name "face_thickness*"'
        curve_ids = set(parse_cubit_list("curve", parse_string))
        
        cubit.cmd(f'curve {l2s(curve_ids.difference(set(remove_thickness_curves)))} interval {cs_params["nel_per_layer"]}')
        cubit.cmd(f'curve with name "*face_web_thickness*" interval {cs_params["nel_per_layer"]}')  #none

        cubit.cmd(f'curve with name "core_thickness*" interval {cs_params["nel_per_core_layer"]}')
        cubit.cmd(f'curve with name "*core_web_thickness*" interval {cs_params["nel_per_core_layer"]}') #seems good

        parse_string = f"with name 'shell*Station*layer0*' except vol with name '*stack{str(hplp_max_stack_ct).zfill(3)}*'"
        vol_to_mesh = parse_cubit_list("volume", parse_string)
        failed_volumes=sweep_volumes(vol_to_mesh)   

        parse_string = f"with name 'shell*Station*layer1*' except vol with name '*stack{str(hplp_max_stack_ct).zfill(3)}*'"
        vol_to_mesh = parse_cubit_list("volume", parse_string)
        failed_volumes+=sweep_volumes(vol_to_mesh)   

        parse_string = f"with name 'shell*Station*layer2*' except vol with name '*stack{str(hplp_max_stack_ct).zfill(3)}*'"
        vol_to_mesh = parse_cubit_list("volume", parse_string)
        failed_volumes+=sweep_volumes(vol_to_mesh)   

        cubit.cmd(f"mesh vol with name 'flatTEad*'")
        cubit.cmd(f"mesh vol with name '*stack{str(hplp_max_stack_ct).zfill(3)}*'")
        cubit.cmd(f"mesh vol with name 'shell_web_thicknessStation*'")
        cubit.cmd(f"mesh vol with name '*layer3*'")
        cubit.cmd(f"mesh vol with name '*layer4*'")

        #Spanwise spacing for every segment then mesh OML surface
        for iLoop, station_id in enumerate(stationList[:-1]): #Skip the last station since there are n_stations-1 spanwise volumes
            e_size_1=e_size[iLoop]
            e_size_2=e_size[iLoop+1]

            parse_string = f'with name "web*Station{str(station_id).zfill(3)}_layer6*_topFace"' #
            surf_to_mesh = parse_cubit_list("surf", parse_string)

            cubit.cmd(f'surface {l2s(surf_to_mesh)} scheme pave')
            cubit.cmd(f'surface {l2s(surf_to_mesh)} size {(e_size_1+e_size_2)/2}')
            cubit.cmd(f'mesh surface {l2s(surf_to_mesh)}')

            parse_string = f'with name "web*Station{str(station_id).zfill(3)}_layer9*_topFace"' #
            surf_to_mesh = parse_cubit_list("surf", parse_string)

            cubit.cmd(f'surface {l2s(surf_to_mesh)} scheme pave')
            cubit.cmd(f'surface {l2s(surf_to_mesh)} size {(e_size_1+e_size_2)/2}')
            cubit.cmd(f'mesh surface {l2s(surf_to_mesh)}')

        cubit.cmd(f"mesh vol with name 'web*Station*'")



    else:
        cubit.cmd(f"reset volume all")
        cubit.cmd(f"mesh volume all")
        

    if get_mesh_error_count():
        with open(f"{wt_name}.log", "a") as logFile:
            logFile.write(f"    Warning: There are {get_mesh_error_count()} mesh errors\n")

    parse_string = f'with not is_meshed'
    un_meshed_vol_ids=parse_cubit_list("volume", parse_string)
    if un_meshed_vol_ids:
        un_meshed_vol_names=[]
        for volume_id in un_meshed_vol_ids:
            un_meshed_vol_names.append(cubit.get_entity_name("volume", volume_id))
        raise ValueError(f'There are {len(un_meshed_vol_ids)} unmeshed volumes. The following volumes failed to mesh: {un_meshed_vol_names}\n\nVolume IDs: {un_meshed_vol_ids}')

    ## Remove materials that were used in one cross section but not the other
    material_not_used = []
    for imat, material_name in enumerate(materials_used):
        parse_string = f'with name "*_{material_name}_*'
        vol_ids=parse_cubit_list("volume", parse_string)

        if not vol_ids:
            material_not_used.append(material_name)

    materials_used = list(materials_used)
    for material_name in material_not_used:
        materials_used.remove(material_name)



    # Blocks
    for imat, material_name in enumerate(materials_used):
        cubit.cmd(f'block {imat+1} add volume with name "*_{material_name}_*"')
        cubit.cmd(f'block {imat+1} name "{material_name}"')

    

    # # Adding Nodesets
    for iLoop, station_id in enumerate(stationList):


        node_set_name=f'station{str(station_id).zfill(3)}'
            
            
        parse_string = f'with name "*station{str(station_id).zfill(3)}*_surface*"'
        surface_ids = parse_cubit_list("surface", parse_string)

        nodeset_id=cubit.get_next_nodeset_id()
        cubit.cmd(f"nodeset {nodeset_id} add surface {l2s(surface_ids)} ")
        cubit.cmd(f'nodeset {nodeset_id} name "{node_set_name}_ns"')

        sideset_id=cubit.get_next_sideset_id()
        cubit.cmd(f"sideset {sideset_id} add surface {l2s(surface_ids)} ")
        cubit.cmd(f'sideset {sideset_id} name "{node_set_name}_ss"')


    # Outer mold-line sideset

    parse_string = f'in curve with name "*oml*"'
    surface_ids = []

    for surf_id in parse_cubit_list("surface", parse_string):
        if 'surface' not in cubit.get_entity_name("surface", surf_id):
            surface_ids.append(surf_id)

    sideset_id=cubit.get_next_sideset_id()
    cubit.cmd(f"sideset {sideset_id} add surface {l2s(surface_ids)} ")
    cubit.cmd(f'sideset {sideset_id} name "oml_ss"')


    cubit.cmd(f"delete curve all with Is_Free except {spanwise_mat_ori_curve}")
    cubit.cmd(f"delete vertex all with Is_Free except {spanwise_mat_ori_curve}")
    
    # if settings["export"] is not None:
    #     if "g" in settings["export"].lower():
    #         cubit.set_element_variable(global_ids, 'rotation_angle_one', theta1s)
    #         cubit.set_element_variable(global_ids, 'rotation_angle_two', theta2s)
    #         cubit.set_element_variable(global_ids, 'rotation_angle_three', theta3s)

    #         cubit.set_element_variable(global_ids, 'rotation_axis_one', 1*np.ones(n_el))
    #         cubit.set_element_variable(global_ids, 'rotation_axis_two', 2*np.ones(n_el))
    #         cubit.set_element_variable(global_ids, 'rotation_axis_three', 3*np.ones(n_el))
    #         cubit.cmd(f'export mesh "{wt_name}.g" overwrite')
    #     if "cub" in settings["export"].lower():
    #         cubit.cmd(f'save as "{wt_name}.cub" overwrite')


    return materials_used
