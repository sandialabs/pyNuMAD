from pynumad.analysis.cubit.utils import *
from pynumad.analysis.cubit.solid_model_utils import *
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
            {'curves': [131, 164, 129, 163], 'verts': [500, 501, 497, 496], 'material_name': 'glass_triax', 'ply_angle': 0}

    birds_mouth_verts: tuple
        Used internally.
    i_stationFirstWeb: int
        Used internally.
    i_stationLastWeb: int
        Used internally. 
    materials_used: set
        Used for in FEA input file generation to define unique materials.
    spanwiseMatOriCurve: int
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
    spanwiseMatOriCurve = 1

    te_types = [station.airfoil.te_type for station in definition.stations]
    roundStations = np.argwhere(np.array(te_types) == "round")
    roundStations = list(roundStations[:, 0])
    lastRoundStation = roundStations[-1]

    with open(f"{wt_name}.log", "w") as logFile:
        logFile.write(f"Making cross sections for {wt_name}\n")

    path_name = directory + "/" + wt_name + "-crossSections"

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

        if geometry.get_profile_te_type(i_station_geometry) == "flat":
            is_flatback = True
        else:
            is_flatback = False

        i_stationFirstWeb = np.argwhere(hasWebs)[0][0]
        i_stationLastWeb = np.argwhere(hasWebs)[-1][0]

        if hasWebs[i_station] == True:
            webNumber = 1
            aft_web_stack = stackdb.swstacks[webNumber][i_station]
            webNumber = 0
            fore_web_stack = stackdb.swstacks[webNumber][i_station]
        else:
            if i_station < i_stationFirstWeb:
                iWebStation = i_stationFirstWeb

            #         elif i_stationLastWeb == len(blade.ispan) - 1-1:
            else:
                iWebStation = i_stationLastWeb
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
        if i_station == i_stationFirstWeb:
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
                is_flatback,
                lastRoundStation,
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
                is_flatback,
                lastRoundStation,
                materials_used,
                cs_normal,
            )
            birds_mouth_verts = []

        cubit.cmd(f"delete curve all with Is_Free except {spanwiseMatOriCurve}")

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

            parse_string = f'with name "*Station{str(i_station)}*"'
            volumeIDs = parse_cubit_list("surface", parse_string)

            # Undo initial twist
            cubit.cmd(
                f"rotate Surface {l2s(volumeIDs)} angle {definition.degreestwist[i_station]} about Z include_merged "
            )

            # Undo prebend
            if definition.prebend[i_station] != 0:
                cubit.cmd(
                    f"move surface {l2s(volumeIDs)} y {-1*definition.prebend[i_station]} include_merged"
                )

            # Undo sweep
            if definition.sweep[i_station] != 0:
                raise ValueError("Presweep is untested for cross-sectional meshing")

            # Mesh the cross-section
            cubit.cmd(
                f'curve with name "layer_thickness*" interval {cs_params["nel_per_layer"]}'
            )
            # cubit.cmd(f'imprint volume {l2s(surface_ids)}')
            cubit.cmd(f"merge volume {l2s(volumeIDs)}")
            cubit.cmd(f"set default autosize on")

            if cs_params["element_shape"].lower() == "tri":
                cubit.cmd(f"surface {l2s(volumeIDs)} scheme tri")
            else:
                cubit.cmd(f"surface {l2s(volumeIDs)} scheme map")

            cubit.cmd(f"mesh surface {l2s(volumeIDs)}")

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
                            volumeIDs,
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
                        cubit.cmd(f"delete curve {spanwiseMatOriCurve}")
                        cubit.cmd(f'save as "{path_name}-{str(i_station)}.cub" overwrite')
                elif len(settings["export"]) == 0:
                    pass
                else:
                    raise NameError(
                        f'Unknown model export format: {settings["export"]}'
                    )
        # elif model2Dor3D.lower() == "3d":
        #     cubit.cmd(f"delete curve {spanwiseMatOriCurve}")
        #     cubit.cmd(f'save as "{path_name}-{str(i_station)}.cub" overwrite')

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
        i_stationFirstWeb,
        i_stationLastWeb,
        materials_used,
        spanwiseMatOriCurve,
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
        "Need more that one cross section to make a solid model"
    """
    if stationList is None or len(stationList) == 0:
        stationList = list(range(len(blade.ispan)))
    elif len(stationList) == 1:
        raise ValueError("Need more that one cross section to make a solid model")

    (
        cubit,
        blade,
        surface_dict,
        birds_mouth_verts,
        i_stationFirstWeb,
        i_stationLastWeb,
        materials_used,
        spanwiseMatOriCurve,
    ) = cubit_make_cross_sections(
        blade, wt_name, settings, cs_params, "3D", stationList
    )

    i_station_start = stationList[0]
    i_station_end = stationList[-1]
    ### ### ### ###
    # Make volumes along the span.
    ### ### ### ###
    meshVolList = []

    part_name = "shell"
    orderedList = get_ordered_list(part_name)
    if len(orderedList) > 0:
        meshVolList = make_all_volumes_for_a_part(surface_dict, orderedList, meshVolList, i_station_end)
    #     cubit.cmd(f'save as "python2.cub" overwrite')
    #     foo

    part_name = "web"
    orderedList = get_ordered_list(part_name)
    orderedListWeb = orderedList.copy()
    if orderedList and len(orderedList[0]) > 1:
        meshVolList = make_all_volumes_for_a_part(surface_dict, orderedList, meshVolList, i_station_end)

    part_name = "roundTEadhesive"
    orderedList = get_ordered_list(part_name)
    if orderedList and len(orderedList[0]) > 1:
        meshVolList = make_all_volumes_for_a_part(surface_dict, orderedList, meshVolList, i_station_end)

    part_name = "flatTEadhesive"
    orderedList = get_ordered_list(part_name)

    if orderedList and len(orderedList[0]) > 1:
        meshVolList = make_all_volumes_for_a_part(surface_dict, orderedList, meshVolList, i_station_end)

    if (
        orderedListWeb
        and len(orderedListWeb[0]) > 1
        and cs_params["birds_mouth_amplitude_fraction"]
        and birds_mouth_verts
    ):
        make_birds_mouth(
            blade,
            birds_mouth_verts,
            cs_params["birds_mouth_amplitude_fraction"],
            i_stationFirstWeb,
            i_stationLastWeb,
        )

    cubit.cmd(f"merge volume {l2s(meshVolList)}")
    cubit.cmd(f"reset volume all")

    cubit.cmd(f"delete surface with Is_Free")
    #cubit.cmd("vol all size 0.2")
    # cubit.cmd(f'curve with name "layer_thickness*" interval {cs_params["nel_per_layer"]}')
    #cubit.cmd("set default autosize on")
    cubit.cmd(f"mesh volume {l2s(meshVolList)}")
    cubit.cmd(f"draw volume {l2s(meshVolList)}")

    # Blocks
    for imat, material_name in enumerate(materials_used):
        cubit.cmd(f'block {imat+1} add volume with name "*{material_name}*"')
        cubit.cmd(f'block {imat+1} name "{material_name}"')

    addColor(blade, "volume")

    # Adding Nodesets
    # Root Nodeset
    parse_string = f'with name "*station{i_station_start}*"'
    print(f"parse_string{parse_string}")
    surface_ids = parse_cubit_list("surface", parse_string)
    cubit.cmd(f"nodeset 1 add surface {l2s(surface_ids)} ")
    cubit.cmd(f'nodeset 1 name "root"')

    for iLoop, i_station in enumerate(stationList[1:-1]):
        parse_string = f'with name "*station{i_station}*"'
        print(f"parse_string{parse_string}")
        surface_ids = parse_cubit_list("surface", parse_string)
        cubit.cmd(f"nodeset {iLoop+2} add surface {l2s(surface_ids)} ")
        cubit.cmd(f'nodeset {iLoop+2} name "station{i_station}"')
    if not stationList[1:-1]:
        iLoop = -1
    # Tip Nodeset
    parse_string = f'with name "*station{i_station_end}*"'
    surface_ids = parse_cubit_list("surface", parse_string)
    cubit.cmd(f"nodeset {iLoop+3} add surface {l2s(surface_ids)} ")
    cubit.cmd(f'nodeset {iLoop+3} name "tip"')

    # Outer mold-line nodeset
    cubit.cmd('draw surf with name "*layer0_bottomFace*"')
    parse_string = f'with name "*layer0_bottomFace*"'
    surface_ids = parse_cubit_list("surface", parse_string)
    cubit.cmd(f"nodeset {iLoop+4} add surface {l2s(surface_ids)} ")
    cubit.cmd(f'nodeset {iLoop+4} name "oml"')

    # ####################################
    # ### Assign material orientations ###
    # ####################################

    parse_string = f'in volume with name "*volume*"'
    allVolumeIDs = parse_cubit_list("volume", parse_string)
    theta2 = {}
    theta1 = {}
    theta3 = {}
    directions = {}

    for volumeID in allVolumeIDs:
        surfIDforMatOri, sign = get_mat_ori_surface(volumeID, spanwiseMatOriCurve)

        for hex_id in get_volume_hexes(volumeID):
            coords = cubit.get_center_point("hex", hex_id)

            if surfIDforMatOri:
                surfaceNormal = vectNorm(
                    list(
                        sign
                        * np.array(get_surface_normal_at_coord(surfIDforMatOri, coords))
                    )
                )

                curveLocationForTangent = cubit.curve(
                    spanwiseMatOriCurve
                ).closest_point(coords)
                x = cubit.curve(spanwiseMatOriCurve).tangent(curveLocationForTangent)[0]
                y = cubit.curve(spanwiseMatOriCurve).tangent(curveLocationForTangent)[1]
                z = cubit.curve(spanwiseMatOriCurve).tangent(curveLocationForTangent)[2]
                spanwiseDirection = vectNorm([x, y, z])

                perimeterDirection = vectNorm(
                    crossProd(surfaceNormal, spanwiseDirection)
                )

                # Recalculate to garantee orthogonal system
                surfaceNormal = crossProd(spanwiseDirection, perimeterDirection)
            else:
                perimeterDirection = [1, 0, 0]
                surfaceNormal = [0, 1, 0]
                spanwiseDirection = [0, 0, 1]

            newCoordinateSystemVectors = [
                spanwiseDirection,
                perimeterDirection,
                surfaceNormal,
            ]
            globalAxisBasisVectors = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]


            global_id=get_global_element_id('hex',hex_id)
            dcm = getDCM(globalAxisBasisVectors, newCoordinateSystemVectors)
            directions[global_id] = newCoordinateSystemVectors

            theta1[global_id], theta2[global_id], theta3[global_id] = dcmToEulerAngles(
                dcm
            )


    if settings["export"] is not None:
        if "g" in settings["export"].lower():
            cubit.cmd(f'export mesh "{wt_name}.g" overwrite')
        if "cub" in settings["export"].lower():
            cubit.cmd(f'save as "{wt_name}.cub" overwrite')
        with open("euler", "wb") as file:
            # A new file will be created
            pickle.dump((theta1, theta2, theta3), file)
        with open("directions", "wb") as file:
            # A new file will be created
            pickle.dump(directions, file)
    return materials_used
