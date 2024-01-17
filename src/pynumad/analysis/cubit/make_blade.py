from pynumad.analysis.cubit.make_cross_sections import *
from pynumad.analysis.cubit.connect_cross_sections import *
from pynumad.utils.orientations import *
import numpy as np
import os
import glob
import pickle


def cubit_make_cross_sections(
    blade,
    wt_name,
    settings,
    cs_params,
    model2Dor3D,
    stationList=None,
    directory=".",
):
    
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
    cubit.cmd("set geometry accuracy 1e-6")
    # making numerus 3D volumes is very slow with autosize on
    cubit.cmd("set default autosize off")

    # Modify blade object to accomodate actual layer thicknesses

    expandTEthicknesses = list(
        cs_params["te_adhesive_thickness"]
        + 6 * cs_params["minimum_layer_thickness"]
    )
    blade.expand_blade_geometry_te(expandTEthicknesses)

    stackdb.edit_stacks_for_solid_mesh()

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

    refLineCoords = np.vstack(
        ([definition.sweep, definition.prebend, definition.ispan])
    ).transpose()
    spanwise_mat_ori_curve = 1

    if model2Dor3D.lower() == "3d":
        write_spline_from_coordinate_points(cubit, refLineCoords)
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
        write_spline_from_coordinate_points(cubit, refLineCoords)
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
        
        # Only save birds_mouth_verts for the right cross-section
        if i_station == i_station_first_web:
            birds_mouth_verts = make_a_cross_section(wt_name,
                surface_dict,
                i_station,
                i_station_geometry,
                blade,
                hasWebs[i_station],
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
                hasWebs[i_station],
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

            for imat, material_name in enumerate(materials_used):
                cubit.cmd(f'block {imat+1} add surface with name "*{material_name}*"')

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
                cubit.cmd(
                    f"move surface {l2s(volume_ids)} y {-1*definition.prebend[i_station]} include_merged"
                )

            # Undo sweep
            if definition.sweep[i_station] != 0:
                raise ValueError("Presweep is untested for cross-sectional meshing")

            # Mesh the cross-section
            cubit.cmd(
                f'curve with name "layer_thickness*" interval {cs_params["nel_per_layer"]}'
            )
            # cubit.cmd(f'imprint volume {l2s(surface_ids)}')
            cubit.cmd(f"merge volume {l2s(volume_ids)}")
            cubit.cmd(f"set default autosize on")

            if cs_params["element_shape"].lower() == "tri":
                cubit.cmd(f"surface {l2s(volume_ids)} scheme tri")
            else:
                cubit.cmd(f"surface {l2s(volume_ids)} scheme map")

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
        refLineCoords[:,0]=np.zeros(len(refLineCoords[:,0]))
        refLineCoords[:,1]=np.zeros(len(refLineCoords[:,0]))
        write_spline_from_coordinate_points(cubit, refLineCoords)

        for i_station in stationList:
            cubit.cmd(f'import cubit "{path_name}-{str(i_station)}.cub"')
        addColor(blade, "surface")
        cubit.cmd(f"delete vertex with Is_Free")
        cubit.cmd(f'save as "{path_name}.cub" overwrite')

        # Remove unnecessary files to save space
        for filePath in glob.glob(f"{path_name}-*.cub"):
            os.remove(filePath)
    return (
        cubit,
        blade,
        surface_dict,
        birds_mouth_verts,
        i_station_first_web,
        i_station_last_web,
        materials_used,
        spanwise_mat_ori_curve,
    )


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

    (
        cubit,
        blade,
        surface_dict,
        birds_mouth_verts,
        i_station_first_web,
        i_station_last_web,
        materials_used,
        spanwise_mat_ori_curve,
    ) = cubit_make_cross_sections(
        blade, wt_name, settings, cs_params, "3D", stationList
    )

    i_station_start = stationList[0]
    i_station_end = stationList[-1]
    ### ### ### ###
    # Make volumes along the span.
    ### ### ### ###
    mesh_vol_list = []

    part_name = "shell"
    ordered_list = get_ordered_list(part_name)
    if len(ordered_list) > 0:
        shell_vol_list = make_all_volumes_for_a_part(surface_dict, ordered_list, i_station_end)
    else:
        shell_vol_list=[]

    part_name = "web"
    ordered_list = get_ordered_list(part_name)
    ordered_list_web = ordered_list.copy()
    if ordered_list and len(ordered_list[0]) > 1:
        web_vol_list = make_all_volumes_for_a_part(surface_dict, ordered_list, i_station_end)
    else:
        web_vol_list=[]

    part_name = "roundTEadhesive"
    ordered_list = get_ordered_list(part_name)
    if ordered_list and len(ordered_list[0]) > 1:
        roundTEadhesive_vol_list = make_all_volumes_for_a_part(surface_dict, ordered_list, i_station_end)
    else:
        roundTEadhesive_vol_list=[]


    part_name = "flatTEadhesive"
    ordered_list = get_ordered_list(part_name)

    if ordered_list and len(ordered_list[0]) > 1:
        flatTEadhesive_vol_list = make_all_volumes_for_a_part(surface_dict, ordered_list, i_station_end)
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

    cubit.cmd(f"merge volume {l2s(mesh_vol_list)}")
    cubit.cmd(f"reset volume all")

    cubit.cmd(f"delete surface with Is_Free")
    #cubit.cmd("vol all size 0.2")
    # cubit.cmd(f'curve with name "layer_thickness*" interval {cs_params["nel_per_layer"]}')
    #cubit.cmd("set default autosize on")
    cubit.cmd(f"mesh volume {l2s(mesh_vol_list)}")
    cubit.cmd(f"draw volume {l2s(mesh_vol_list)}")

    if get_mesh_error_count() !=0:
        with open(f"{wt_name}.log", "a") as logFile:
            logFile.write(f"    Warning: {get_mesh_error_count()} cross section mesh errors exist in station {i_station}\n")
            
    # Blocks
    for imat, material_name in enumerate(materials_used):
        cubit.cmd(f'block {imat+1} add volume with name "*{material_name}*"')
        cubit.cmd(f'block {imat+1} name "{material_name}"')

    addColor(blade, "volume")


    # # Adding Nodesets
    for iLoop, i_station in enumerate(stationList):

        if iLoop ==0:
            node_set_name='root'
        elif iLoop==len(stationList)-1:
            node_set_name='tip'
        else:
            node_set_name=f'station{str(i_station).zfill(3)}'
            
            
        parse_string = f'with name "*station{str(i_station).zfill(3)}*"'
        surface_ids = parse_cubit_list("surface", parse_string)

        cubit.cmd(f"nodeset {iLoop+1} add surface {l2s(surface_ids)} ")
        cubit.cmd(f'nodeset {iLoop+1} name "{node_set_name}"')



    # Outer mold-line nodeset
    cubit.cmd('draw surf with name "*layer0_bottomFace*"')
    parse_string = f'with name "*layer0_bottomFace*"'
    surface_ids = parse_cubit_list("surface", parse_string)
    cubit.cmd(f"nodeset {iLoop+1} add surface {l2s(surface_ids)} ")
    cubit.cmd(f'nodeset {iLoop+1} name "oml"')

    # ####################################
    # ### Assign material orientations ###
    # ####################################

    parse_string = f'in volume with name "*volume*"'
    all_volume_ids = parse_cubit_list("volume", parse_string)

    global_ids=[]
    theta1s=[]
    theta2s=[]
    theta3s=[]

    for volume_id in all_volume_ids:
        surf_id_for_mat_ori, sign = get_mat_ori_surface(volume_id, spanwise_mat_ori_curve)
        volume_name = cubit.get_entity_name("volume", volume_id)
        for hex_id in get_volume_hexes(volume_id):
            coords = cubit.get_center_point("hex", hex_id)

            if surf_id_for_mat_ori:
                surface_normal = vectNorm(
                    list(
                        sign
                        * np.array(get_surface_normal_at_coord(surf_id_for_mat_ori, coords))
                    )
                )
                if 'web' in volume_name.lower():
                    spanwise_direction = [0,0,1]
                else:
                    curve_location_for_tangent = cubit.curve(
                        spanwise_mat_ori_curve
                    ).closest_point(coords)
                    x = cubit.curve(spanwise_mat_ori_curve).tangent(curve_location_for_tangent)[0]
                    y = cubit.curve(spanwise_mat_ori_curve).tangent(curve_location_for_tangent)[1]
                    z = cubit.curve(spanwise_mat_ori_curve).tangent(curve_location_for_tangent)[2]
                    spanwise_direction = vectNorm([x, y, z])

                perimeter_direction = vectNorm(
                    crossProd(surface_normal, spanwise_direction)
                )

                # Recalculate to garantee orthogonal system
                surface_normal = crossProd(spanwise_direction, perimeter_direction)
            else:
                perimeter_direction = [1, 0, 0]
                surface_normal = [0, 1, 0]
                spanwise_direction = [0, 0, 1]

            newCoordinateSystemVectors = [
                spanwise_direction,
                perimeter_direction,
                surface_normal,
            ]
            globalAxisBasisVectors = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]


            global_id=get_global_element_id('hex',hex_id)
            
            dcm = getDCM(globalAxisBasisVectors, newCoordinateSystemVectors)

            temp1, temp2, temp3 = dcmToEulerAngles(dcm)

            global_ids.append(global_id)
            theta1s.append(-1*temp1)
            theta2s.append(-1*temp2)
            theta3s.append(-1*temp3)

    n_el=len(global_ids)
    cubit.cmd(f"delete curve all with Is_Free except {spanwise_mat_ori_curve}")
    cubit.cmd(f"delete vertex all with Is_Free except {spanwise_mat_ori_curve}")
    if settings["export"] is not None:
        if "g" in settings["export"].lower():
            cubit.set_element_variable(global_ids, 'rotation_angle_one', theta1s)
            cubit.set_element_variable(global_ids, 'rotation_angle_two', theta2s)
            cubit.set_element_variable(global_ids, 'rotation_angle_three', theta3s)

            cubit.set_element_variable(global_ids, 'rotation_axis_one', 1*np.ones(n_el))
            cubit.set_element_variable(global_ids, 'rotation_axis_two', 2*np.ones(n_el))
            cubit.set_element_variable(global_ids, 'rotation_axis_three', 3*np.ones(n_el))
            cubit.cmd(f'export mesh "{wt_name}.g" overwrite')
        if "cub" in settings["export"].lower():
            cubit.cmd(f'save as "{wt_name}.cub" overwrite')


    return materials_used
