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

            hoop_direction = vectNorm(np.cross(surface_normal, spanwise_direction))

            newCoordinateSystemVectors = [spanwise_direction,hoop_direction,surface_normal]

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
def get_orientations_euler(volume_id,element_shape_string):
    global_el_ids_in_vol=[]
    theta1s_in_vol=[]
    theta2s_in_vol=[]
    theta3s_in_vol=[]



    surf_id_for_mat_ori,sign = get_mat_ori_surface(volume_id)
    #volume_name = cubit.get_entity_name("volume", volume_id)
    #t0 = time.time()

    if 'hex' in element_shape_string:
        element_ids = get_volume_hexes(volume_id)
    elif 'tet' in element_shape_string:
        element_ids = get_volume_tets(volume_id)
        

    for el_id in element_ids:
        coords = cubit.get_center_point(element_shape_string, el_id)
            
        surface_normal = vectNorm(
            list(sign*np.array(get_surface_normal_at_coord(surf_id_for_mat_ori, coords))))

        ref_line_direction = [0,0,1]
        #https://www.maplesoft.com/support/help/maple/view.aspx?path=MathApps%2FProjectionOfVectorOntoPlane
        spanwise_direction = vectNorm(np.array(ref_line_direction)-np.dot(ref_line_direction,surface_normal)*np.array(surface_normal))

        hoop_direction = vectNorm(np.cross(surface_normal, spanwise_direction))

        newCoordinateSystemVectors = [spanwise_direction,hoop_direction,surface_normal]

        globalAxisBasisVectors = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

        global_id=get_global_element_id(element_shape_string,el_id)
        
        dcm = getDCM(globalAxisBasisVectors, newCoordinateSystemVectors)

        temp1, temp2, temp3 = dcmToEulerAngles(dcm)

        global_el_ids_in_vol.append(global_id)
        theta1s_in_vol.append(-1*temp1)
        theta2s_in_vol.append(-1*temp2)
        theta3s_in_vol.append(-1*temp3)


    return global_el_ids_in_vol,theta1s_in_vol,theta2s_in_vol,theta3s_in_vol

def get_element_orientations_vectors(element_ids,mat_ori_surfs,signs,thetas):

    spanwise_directions = []
    hoop_directions = []
    surface_normal_directions = []

    for el_id in element_ids:

        surf_id_for_mat_ori = int(mat_ori_surfs[el_id-1])
        sign = signs[el_id-1]

        node_ids = parse_cubit_list('node',f'in element {el_id}')
        coords = [list(get_nodal_coordinates(node_id)) for node_id in node_ids]
        coords = np.array(coords)
        coords = np.mean(coords, 0)

        surface_normal = vectNorm(
            list(sign*np.array(get_surface_normal_at_coord(surf_id_for_mat_ori, coords))))

        ref_line_direction = [0,0,1]
        #https://www.maplesoft.com/support/help/maple/view.aspx?path=MathApps%2FProjectionOfVectorOntoPlane
        spanwise_direction = vectNorm(np.array(ref_line_direction)-np.dot(ref_line_direction,surface_normal)*np.array(surface_normal))

        hoop_direction = vectNorm(np.cross(surface_normal, spanwise_direction))
       
        #Additional rotation about surface normal
        theta = thetas[el_id-1]
        c = math.cos(theta)
        s = math.sin(theta)
        beta = np.array([[c,s, 0],[-s, c, 0],[0,0,1]])

        new_directions = beta @ [spanwise_direction,hoop_direction,surface_normal]

        spanwise_directions.append(list(new_directions[0]))
        hoop_directions.append(list(new_directions[1]))
        surface_normal_directions.append(list(new_directions[2]))
    


    return spanwise_directions,hoop_directions,surface_normal_directions
def assign_material_orientation_angles(orientation_data):
    #Assigne Euler angles to exodus mesh variables
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
def assign_material_orientation_vectors(orientation_data):
    #Assign material orientation vectors to exodus mesh variables
    global_ids=orientation_data[0]
    n_el = len(global_ids)

    

    one_axis=np.array(orientation_data[1])
    two_axis=np.array(orientation_data[2])
    three_axis=np.array(orientation_data[3])

    cubit.set_element_variable(global_ids, 'matCoord_1_x', one_axis[:,0])
    cubit.set_element_variable(global_ids, 'matCoord_1_y', one_axis[:,1])
    cubit.set_element_variable(global_ids, 'matCoord_1_z', one_axis[:,2])

    cubit.set_element_variable(global_ids, 'matCoord_2_x', two_axis[:,0])
    cubit.set_element_variable(global_ids, 'matCoord_2_y', two_axis[:,1])
    cubit.set_element_variable(global_ids, 'matCoord_2_z', two_axis[:,2])

    cubit.set_element_variable(global_ids, 'matCoord_3_x', three_axis[:,0])
    cubit.set_element_variable(global_ids, 'matCoord_3_y', three_axis[:,1])
    cubit.set_element_variable(global_ids, 'matCoord_3_z', three_axis[:,2])

    return

def get_material_orientation_angles(orientation_vectors):
    # orientation_vectors is the output of assign_material_orientation_vectors()

    t0 = time.time()
    print(f'Converting orientation vectors to Euler Angles...')

    globalAxisBasisVectors = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    R_1_angles = []
    R_2_angles = []
    R_3_angles = []
    for i_el, element_id in enumerate(orientation_vectors[0]):
        spanwise_direction = orientation_vectors[1][i_el]
        hoop_direction = orientation_vectors[2][i_el]
        surface_normal = orientation_vectors[3][i_el]

        newCoordinateSystemVectors = [spanwise_direction,hoop_direction,surface_normal]

        dcm = getDCM(globalAxisBasisVectors, newCoordinateSystemVectors)

        temp1, temp2, temp3 = dcmToEulerAngles(dcm)

        R_1_angles.append(-1*temp1)
        R_2_angles.append(-1*temp2)
        R_3_angles.append(-1*temp3)

    t1 = time.time()
    print(f'Total time for Euler angles: {t1-t0}')
    return orientation_vectors[0],R_1_angles,R_2_angles,R_3_angles


def get_material_orientation_vectors(ncpus = 1):
    # # ####################################
    # # ### Get material orientations ###
    # # ####################################
    cubit.cmd("renumber element all start_id 1 uniqueids")

    parse_string = f'in volume with name "*volume*"'
    global_element_ids = parse_cubit_list("element", parse_string)

    t0 = time.time()
    mat_ori_surfs = np.zeros(len(global_element_ids))
    signs = np.zeros(len(global_element_ids))
    thetas = np.zeros(len(global_element_ids))
    
    parse_string = f'in volume with name "*volume*"'
    volume_ids = parse_cubit_list("volume", parse_string)

    for volume_id in volume_ids:
        surf_id_for_mat_ori,sign = get_mat_ori_surface(volume_id)

        parse_string = f'in volume {volume_id}'
        this_volume_element_ids = parse_cubit_list("element", parse_string)

        for el_id in this_volume_element_ids:
            mat_ori_surfs[el_id-1] = surf_id_for_mat_ori
            signs[el_id-1] = sign
            ply_angle = float(get_entity_name("volume", volume_id).split('_')[-2])
            thetas[el_id-1] = math.radians(ply_angle)
    t1 = time.time()
    print(f'Total time for material orientation arrays: {t1-t0}')

    t0 = time.time()
    print(f'Calculating material orientations with {ncpus} CPU(s)...')
    spanwise_directions,hoop_directions,surface_normal_directions = get_element_orientations_vectors(global_element_ids,mat_ori_surfs,signs,thetas)
    t1 = time.time()
    print(f'Total time for material orientations: {t1-t0}')

    return [global_element_ids,spanwise_directions,hoop_directions,surface_normal_directions]


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
        stationList = list(range(len(definition.ispan)))

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
    for i_station_geometry in range(len(blade.definition.ispan)-1): #-1 b/c fewer stacks than stations
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
                last_flat_station+=1
            else:
                break
        last_flat_station+=1

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
        if i_station == len(definition.ispan) - 1:  # Only do this for the last station
            blade.add_interpolated_station(definition.ispan[-1] * 0.999)
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

        if 'precomp' in settings['make_input_for']:
            make_a_precomp_cross_section(wt_name,
                                surface_dict,i_station,i_station_geometry,blade,make_webs,aft_web_stack,fore_web_stack,iLE,cs_params,
                                geometry_scaling,thickness_scaling,last_round_station,last_flat_station,materials_used,cs_normal)
        else:
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
        cubit.cmd(f"delete block all")
        #Since cross sections were translated for cross sectional codes, remove prebend and sweep from ref axis.
        ref_line_coords[:,0]=np.zeros(len(ref_line_coords[:,0]))
        ref_line_coords[:,1]=np.zeros(len(ref_line_coords[:,0]))
        write_spline_from_coordinate_points(cubit, ref_line_coords)

        for i_station in stationList:
            cubit.cmd(f'import cubit "{path_name}-{str(i_station)}.cub"')
        addColor(blade, "surface")

        for imat, material_name in enumerate(materials_used):
            cubit.cmd(f'block {imat+1} add surface with name "*{material_name}*"')
            cubit.cmd(f'block {imat+1} name "{material_name}"')
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
    volume_dict ={}

    part_name = "shell"
    ordered_list = get_ordered_list(part_name)
    spanwise_splines=[]
    if len(ordered_list) > 0:
        shell_vol_list,spanwise_splines = make_all_volumes_for_a_part(surface_dict, volume_dict, ordered_list, i_station_end,spanwise_splines)
    else:
        shell_vol_list=[]

    part_name = "web"
    ordered_list = get_ordered_list(part_name)
    ordered_list_web = ordered_list.copy()
    if ordered_list and len(ordered_list[0]) > 1:
        web_vol_list,spanwise_splines = make_all_volumes_for_a_part(surface_dict, volume_dict, ordered_list, i_station_end,spanwise_splines)
    else:
        web_vol_list=[]

    part_name = "roundTEadhesive"
    ordered_list = get_ordered_list(part_name)
    if ordered_list and len(ordered_list[0]) > 1:
        roundTEadhesive_vol_list,spanwise_splines = make_all_volumes_for_a_part(surface_dict, volume_dict, ordered_list, i_station_end,spanwise_splines)
    else:
        roundTEadhesive_vol_list=[]


    part_name = "flatTEadhesive"
    ordered_list = get_ordered_list(part_name)

    if ordered_list and len(ordered_list[0]) > 1:
        flatTEadhesive_vol_list,spanwise_splines = make_all_volumes_for_a_part(surface_dict, volume_dict, ordered_list, i_station_end,spanwise_splines)
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
        cubit.cmd(f'save as "debug.cub" overwrite')
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

        # sideset_id=cubit.get_next_sideset_id()
        # cubit.cmd(f"sideset {sideset_id} add surface {l2s(surface_ids)} ")
        # cubit.cmd(f'sideset {sideset_id} name "{node_set_name}_ss"')


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
    


    return materials_used

def yaml_mesh_to_cubit(dir_name,yaml_file_base,plot_mat_ori = True):
    import yaml
    
    print(f'Importing {dir_name}/{yaml_file_base} to cubit ...')
    with open(f'{dir_name}/{yaml_file_base}.yaml', 'r') as file:
        mesh_data = yaml.load(file, Loader=yaml.CLoader)

    print('    Making nodes ...')
    for node_coords in mesh_data['nodes']:
        cubit.silent_cmd(f'create node location {node_coords[0]}')

    print('    Making elements ...')
    n_node_per_el = len(mesh_data['elements'][0][0].split())

    if n_node_per_el == 3:
        element_type = 'tri'
    elif n_node_per_el == 4:
        element_type = 'face'
    elif n_node_per_el == 6:
        element_type = 'tet'
    elif n_node_per_el == 8:
        element_type = 'hex'

    for element_conn in mesh_data['elements']:

        cubit.silent_cmd(f'create {element_type} node {element_conn[0]}')

    print('    Making blocks ...')
    for i_mat, mat_data in enumerate(mesh_data['materials']):
        mat_name = mat_data['name']
        for el_set in mesh_data['sets']['element']:
            if mat_name.lower() == el_set['name'].lower():
                cubit.silent_cmd(f'block {i_mat+1} add {element_type} {str(el_set["labels"]).replace("[", "").replace("]", "")}')
                cubit.silent_cmd(f'block {i_mat+1} name "{mat_name}"')

    #Material Orientations
    spanwise_directions = []
    hoop_directions = []
    surface_normal_directions = []

    print('    Assigning Material orientations ...')
    # global_element_ids = tuple(np.arange(len(mesh_data['elements']))+1)
    global_element_ids = []
    for i_el, element_ori in enumerate(mesh_data['elementOrientations']):
        hex_id = i_el+1   
        global_element_ids.append(get_global_element_id('hex',hex_id))             

        spanwise_directions.append(element_ori[:3])
        hoop_directions.append(element_ori[3:6])
        surface_normal_directions.append(element_ori[6:])

    
    orientation_vectors=[global_element_ids,spanwise_directions,hoop_directions,surface_normal_directions]       
    
    # assign_material_orientation_vectors(orientation_vectors)

    if plot_mat_ori:
        print('    Plotting material orientation lines ...')
        dir_strings = ['xdir','ydir','zdir',]
        for i_el, element_ori in enumerate(mesh_data['elementOrientations']):
            hex_id = i_el+1                

            node_ids = cubit.get_expanded_connectivity(element_type, hex_id)
            coords = []
            for iNd, node_id in enumerate(node_ids):
                coords.append(list(get_nodal_coordinates(node_id)))
            coords = np.array(coords)

            # #######For Plotting - find the largest element side length #######
            distances=[]
            for iNd,node_id in enumerate(node_ids):
                for jNd,node_idj in enumerate(node_ids):
                    distances.append(norm(vectSub(coords[iNd],coords[jNd])))
            length=max(distances)

            coords = np.mean(coords, 0)

            cubit.create_vertex(coords[0],coords[1],coords[2])
            iVert1=get_last_id("vertex")
            for i_dir in range(3):
                index = 3*i_dir
                cubit.create_vertex(coords[0]+length*element_ori[index],coords[1]+length*element_ori[index+1],coords[2]+length*element_ori[index+2])
                iVert2=get_last_id("vertex")
                cubit.silent_cmd(f'create curve vertex {iVert1} {iVert2}')
                cubit.silent_cmd(f'curve {get_last_id("curve")} name "{dir_strings[i_dir]}"')


        cubit.silent_cmd(f'color curve with name "xdir*"  geometry seagreen')
        cubit.silent_cmd(f'color curve with name "ydir*"  geometry blue')
        cubit.silent_cmd(f'color curve with name "zdir*"  geometry red')

    print(f'Done importing! \n')

def get_all_cross_section_surface_ids(station_list):


    all_cross_section_surface_ids=[]

    for i_station in station_list:
        parse_string = f'with name "*Station{str(i_station).zfill(3)}*surface*"'
        all_cross_section_surface_ids.append(parse_cubit_list("surface", parse_string))
    return all_cross_section_surface_ids
def get_all_segments_volume_ids(station_list,grouped_segments):
    from copy import deepcopy
    station_list.pop(-1)
    all_segments_volume_ids=[[]]

    #Make sure grouped stations are in ascending order
    for group in grouped_segments:
        group.sort()

    #Make sure grouped stations are in ascending order
    for igroup, group in enumerate(grouped_segments):
        if len(group)>0:
            i_start = group[0]
            i_end =group[-1]
            grouped_segments[igroup]=list(range(i_start,i_end+1))

    ## Make sure grouped segments are in station list: 
    for i_staion in [ x for xs in grouped_segments for x in xs]:
        if i_staion not in station_list:
            raise IOError(f'Grouped station {i_staion} not in station_list')

    segment_ct=0
    group_ct=0
    for i_station in station_list:

        parse_string = f'with name "*Station{str(i_station).zfill(3)}*"'
        segment_volume_ids = parse_cubit_list("volume", parse_string)
        # my_list = [1, 2, 3, 4, 5]
        # segment_volume_ids = [i * (i_station+1) for i in my_list]
        if i_station not in grouped_segments[group_ct]:
            all_segments_volume_ids[segment_ct] = deepcopy(segment_volume_ids)
            all_segments_volume_ids.append([])
            segment_ct+=1
        else:
            if i_station == grouped_segments[group_ct][0]:
                all_segments_volume_ids[segment_ct] = deepcopy(segment_volume_ids)
            else:
                all_segments_volume_ids[segment_ct]+=segment_volume_ids
            if i_station == grouped_segments[group_ct][-1]:
                segment_ct+=1
                all_segments_volume_ids.append([])
                if group_ct < len(grouped_segments)-1:
                    group_ct+=1
    all_segments_volume_ids.pop(-1)

    return all_segments_volume_ids,grouped_segments

#################################
def export_segments_to_yaml(case_type,blade,station_list,grouped_segments,mesh_yaml_file_name_base,directory):

    materials_used  = ['medium_density_foam','Adhesive','glass_biax','glass_uni','glass_triax','carbon_uni_industry_baseline']

    block_ids = parse_cubit_list('block',f'in vol all')  
    materials_used =[cubit.get_exodus_entity_name('block', block_id) for block_id in block_ids]
    if station_list is None or len(station_list) == 0:
        station_list=[]
        for i_station in range(100): #Assuming there are 100 stations or less. 
            if len(parse_cubit_list('surface', f'with name "*tation{str(i_station).zfill(3)}*surface*')):
                station_list.append(i_station)

    orientation_vectors=get_material_orientation_vectors(ncpus = 1)
    if 'solid_cross_section' in case_type.lower(): 
        all_segments_volume_ids = get_all_cross_section_surface_ids(station_list)
        unit_type = 'surf'
        element_type = 'face'

        # flattened_list = [item for sublist in all_segments_volume_ids for item in sublist]
        # all_face_ids = parse_cubit_list(element_type, f'in {unit_type} {l2s(flattened_list)}')
        
        # for face_id in all_face_ids:
        #     hex_ids = parse_cubit_list('hex',f'in face {face_id}')
        #     if len(hex_ids) > 0:
        #         z_coord = []
        #         for hex_id in hex_id:
        #             z_coord.append(get_center_point("hex", hex_id)[-1])
        #         hex_id = hex_ids[np.argmax(z_coord)]
        #         element_id = get_hex_global_element_id(hex_id)

            # else:
            #     raise RuntimeError(f'No hex with face {face_id} found')

            # get_hex_global_element_id(1537)


    elif 'solid_segment' in case_type.lower():
        
        print('Partitioning blade volumes into segments ...')
        all_segments_volume_ids,grouped_segments = get_all_segments_volume_ids(station_list,grouped_segments)
        unit_type = 'vol'
        element_type = 'element'
        mesh_yaml_file_name_base+='_segment'

    for iseg,segment_volume_ids in enumerate(all_segments_volume_ids):

        print(f'Writing segment {iseg} ...')

        file_name=f'{mesh_yaml_file_name_base}-{station_list[iseg]}.yaml'
        path_name = directory + "/" + file_name
        
        segment_element_ids = parse_cubit_list(element_type, f'in {unit_type} {l2s(segment_volume_ids)}')
        segment_node_ids = parse_cubit_list('node', f'in {unit_type} {l2s(segment_volume_ids)}')

        # Write Nodes
        print(f'    Writing {len(segment_node_ids)} nodes...')
        data_string =''
        data_string+="nodes:\n"
        # data_string+=' '
        temp_str=[' - '+str(list(get_nodal_coordinates(node_id))) + '\n' for node_id in segment_node_ids]
        separator = ""
        temp_str = separator.join(temp_str).replace(',','')

        data_string+=temp_str
        with open(path_name, "w") as f:
            f.write(data_string)

        # Write Elements
        data_string =''
        print(f'    Writing {len(segment_element_ids)} elements...')
        data_string =''
        data_string+="elements:\n"
        # data_string+=' '

        node_index_dict = dict((value, idx) for idx,value in enumerate(segment_node_ids))

        temp_str=[' - '+str([node_index_dict[node_id]+1 for node_id in cubit.get_expanded_connectivity(element_type, segment_element_id)]) +'\n'for segment_element_id in segment_element_ids]
        separator = ""
        temp_str = separator.join(temp_str).replace(',','')

        data_string+=temp_str
        with open(path_name, "a") as f:
            f.write(data_string)

        # Write Element Sets
        print(f'    Writing element sets...')
        data_string =''
        data_string+="sets:\n"

        data_string+=f"  element:\n"
        index_dict = dict((value, idx) for idx,value in enumerate(segment_element_ids))
        for material_name in materials_used:
            
            if 'solid_cross_section' in case_type.lower(): 
                element_list = list(parse_cubit_list('face', f'in surf with name "*tation{str(station_list[iseg]).zfill(3)}_*{material_name}*surface*"'))

            elif 'solid_segment' in case_type.lower():
                block_elements = set(parse_cubit_list('element', f'in block with name "{material_name}"'))
                segment_elements = set(parse_cubit_list('element', f'in {unit_type} {l2s(segment_volume_ids)}'))
                element_list = block_elements.intersection(segment_elements)

            data_string+=f"  - name: {material_name}\n"
            data_string+=f"    labels:\n"


            temp_str = str([index_dict[iEl]+1 for iEl in element_list])
            temp_str = temp_str.replace("[", "")
            temp_str = temp_str.replace("]", "")

            temp_str = temp_str.replace(",", "\n    -")
            data_string+=f"    - {temp_str}\n"

        with open(path_name, "a") as f:
            f.write(data_string)

        # Write Element Orientations
        print(f'    Writing material orientations...')
        data_string =''
        data_string+="elementOrientations:\n"
        # data_string+=' '

        index_dict = dict((value, idx) for idx,value in enumerate(orientation_vectors[0]))

        # '- '+separator.join([str(orientation_vectors[i_dir][22]) for i_dir in [1,2,3]]).replace("], [", ", ")+'\n'
        # data_string+=separator.join(['- '+separator.join([str(orientation_vectors[i_dir][index_dict[segment_element_id]]) for i_dir in [1,2,3]]).replace("], [", ", ")+'\n' for segment_element_id in segment_element_ids]).replace("\n, - ", "\n - ")
        if 'solid_cross_section' in case_type.lower():
            # parse_string = f'in node in surf with name "*tation{str(iseg).zfill(3)}*surface*" except hex in node with z_coord < {blade.ispan[station_list[iseg]]}'
            # parse_string = f'in face {l2s(segment_element_ids)} and hex in node with z_coord < {blade.ispan[station_list[iseg]]}'
            if station_list[iseg] == station_list[-1]:
                parse_string = f'in face {l2s(segment_element_ids)} in hex with z_coord < {blade.ispan[station_list[iseg]]}'
            else:
                parse_string = f'in face {l2s(segment_element_ids)} in hex with z_coord > {blade.ispan[station_list[iseg]]}'

            #Overwrite segment_element_ids since not used downstream
            temp_element_ids =  list(parse_cubit_list('hex', parse_string))
            
            segment_element_ids=[get_hex_global_element_id(hex_id) for hex_id in temp_element_ids]


        temp_str=[f" - [{', '.join(str(e) for e in orientation_vectors[1][index_dict[segment_element_id]])}, {', '.join(str(e) for e in orientation_vectors[2][index_dict[segment_element_id]])}, {', '.join(str(e) for e in orientation_vectors[3][index_dict[segment_element_id]])}]\n" for segment_element_id in segment_element_ids]
        separator = ""
        temp_str = separator.join(temp_str)
        data_string+=temp_str
        with open(path_name, "a") as f:
            f.write(data_string)

        #Write Materials
        print(f'    Writing materials...')
        data_string =''
        data_string+="materials:\n"

        materials = blade.definition.materials
        for material_name in materials_used:
            material=materials[material_name]

            data_string+=f'   -  name: {material.name}\n'
            data_string+=f'      E: [{material.ex}, {material.ey}, {material.ez}]\n'
            data_string+=f'      G: [{material.gxy}, {material.gxz}, {material.gyz}]\n'
            data_string+=f'      nu: [{material.prxy}, {material.prxz}, {material.pryz}]\n'
            data_string+=f'      rho: {material.density}\n'
        with open(path_name, "a") as f:
            f.write(data_string)


