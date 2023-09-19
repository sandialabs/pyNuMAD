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
    surfaceDict: dict
        Keys are integers for the Cubit surface IDs for the cross sections. Each surface has
        it's own dictionary with the following keys: 'curves', 'verts', 'materialName', 'plyAngle'.
        
        e.g. 
        >> surfaceDict[9]
            {'curves': [131, 164, 129, 163], 'verts': [500, 501, 497, 496], 'materialName': 'glass_triax', 'plyAngle': 0}

    birdsMouthVerts: tuple
        Used internally.
    iStationFirstWeb: int
        Used internally.
    iStationLastWeb: int
        Used internally. 
    materialsUsed: set
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
    surfaceDict = {}
    # Uniquly track which materiall IDs are actuall used in blade model
    materialsUsed = set()
    iLE = geometry.LEindex + 1
    thicknessScaling = 0.001
    geometryScaling = thicknessScaling * 1000

    # Set up Cubit
    cubit.init(["cubit", "-nojournal"])

    cubit.cmd("undo off")
    cubit.cmd("set geometry accuracy 1e-6")
    # making numerus 3D volumes is very slow with autosize on
    cubit.cmd("set default autosize off")

    # Modify blade object to accomodate actual layer thicknesses

    expandTEthicknesses = list(
        cs_params["TE_adhesive_thickness"]
        + 6 * cs_params["minimum_layer_thickness"]
    )
    blade.expand_blade_geometry_te(expandTEthicknesses)

    stackdb.edit_stacks_for_solid_mesh()

    hasWebs = []
    webNumber = 1
    for iStation in range(len(stackdb.swstacks[webNumber])):
        if not len(stackdb.swstacks[webNumber][iStation].plygroups) == 0:
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

    with open("make_blade.log", "w") as logFile:
        logFile.write(f"Making cross sections for {wt_name}\n")

    pathName = directory + "/" + wt_name + "-crossSections"

    for iStation in stationList:
        if model2Dor3D.lower() == "2d":
            cubit.cmd(
                "reset "
            )  # This is needed to restart node numbering for VABS. VABS neeeds every element and node starting from 1 to nelem/nnode should be present
        writeSplineFromCoordinatePoints(cubit, refLineCoords)
        iStationGeometry = iStation
        if iStation == len(geometry.ispan) - 1:  # Only do this for the last station
            blade.add_interpolated_station(geometry.ispan[-1] * 0.999)
            stackdb.edit_stacks_for_solid_mesh()
            expandTEthicknesses.append(expandTEthicknesses[-1])
            blade.expand_blade_geometry_te(expandTEthicknesses)

            # adjustLastStackAfterNewTipStation(iStation)

            iStationGeometry = iStation + 1

        if geometry.get_profile_te_type(iStationGeometry) == "flat":
            isFlatback = True
        else:
            isFlatback = False

        iStationFirstWeb = np.argwhere(hasWebs)[0][0]
        iStationLastWeb = np.argwhere(hasWebs)[-1][0]

        if hasWebs[iStation] == True:
            webNumber = 1
            aftWebStack = stackdb.swstacks[webNumber][iStation]
            webNumber = 0
            foreWebStack = stackdb.swstacks[webNumber][iStation]
        else:
            if iStation < iStationFirstWeb:
                iWebStation = iStationFirstWeb

            #         elif iStationLastWeb == len(blade.ispan) - 1-1:
            else:
                iWebStation = iStationLastWeb
            #         else:
            #             raise Exception('assuming web ends at last station for now. ')

            webNumber = 1
            aftWebStack = stackdb.swstacks[webNumber][iWebStation]
            webNumber = 0
            foreWebStack = stackdb.swstacks[webNumber][iWebStation]

        crossSectionNormal = getCrossSectionNormalVector(
            np.array(
                [
                    keypoints.key_points[2, :, iStationGeometry],
                    keypoints.key_points[3, :, iStationGeometry],
                    keypoints.key_points[7, :, iStationGeometry],
                ]
            )
        )

        # Only save birdsMouthVerts for the right cross-section
        if iStation == iStationFirstWeb:
            birdsMouthVerts = writeCubitCrossSection(
                surfaceDict,
                iStation,
                iStationGeometry,
                blade,
                hasWebs[iStation],
                aftWebStack,
                foreWebStack,
                iLE,
                cs_params,
                geometryScaling,
                thicknessScaling,
                isFlatback,
                lastRoundStation,
                materialsUsed,
                crossSectionNormal,
            )
        else:
            writeCubitCrossSection(
                surfaceDict,
                iStation,
                iStationGeometry,
                blade,
                hasWebs[iStation],
                aftWebStack,
                foreWebStack,
                iLE,
                cs_params,
                geometryScaling,
                thicknessScaling,
                isFlatback,
                lastRoundStation,
                materialsUsed,
                crossSectionNormal,
            )
            birdsMouthVerts = []

        cubit.cmd(f"delete curve all with Is_Free except {spanwiseMatOriCurve}")

        # Chord line for rotation of cross-section for homogenization
        if model2Dor3D.lower() == "2d":
            #         #Blocks

            for imat, materialName in enumerate(materialsUsed):
                cubit.cmd(f'block {imat+1} add surface with name "*{materialName}*"')

            addColor(blade, "surface")

            # create_vertex(blade.geometry[0,0,iStation]*geometryScaling,blade.geometry[0,1,iStation]*geometryScaling,blade.geometry[0,2,iStation]*geometryScaling)
            # TEvert=get_last_id("vertex")
            # create_vertex(blade.geometry[iLE-1,0,iStation]*geometryScaling,blade.geometry[iLE-1,1,iStation]*geometryScaling,blade.geometry[iLE-1,2,iStation]*geometryScaling)
            # LEvert=get_last_id("vertex")

            # cubit.cmd(f'create curve vertex {TEvert} {LEvert}')
            # coords=cubit.vertex(TEvert).coordinates()
            # tangent=cubit.curve(get_last_id("curve")).tangent(coords)
            # tangentDirection=vectNorm(list(tangent))  #Unit vector of tangent.
            # crossSectionRotationAngle=math.atan2(tangentDirection[1],tangentDirection[0])*180/pi

            parseString = f'with name "*Station{str(iStation)}*"'
            volumeIDs = parse_cubit_list("surface", parseString)

            # Undo initial twist
            cubit.cmd(
                f"rotate Surface {l2s(volumeIDs)} angle {definition.degreestwist[iStation]} about Z include_merged "
            )

            # Undo prebend
            if definition.prebend[iStation] != 0:
                cubit.cmd(
                    f"move surface {l2s(volumeIDs)} y {-1*definition.prebend[iStation]} include_merged"
                )

            # Undo sweep
            if definition.sweep[iStation] != 0:
                raise ValueError("Presweep is untested for cross-sectional meshing")

            # Mesh the cross-section
            cubit.cmd(
                f'curve with name "layerThickness*" interval {cs_params["nel_per_layer"]}'
            )
            # cubit.cmd(f'imprint volume {l2s(surfaceIDs)}')
            cubit.cmd(f"merge volume {l2s(volumeIDs)}")
            cubit.cmd(f"set default autosize on")

            if cs_params["element_shape"].lower() == "tri":
                cubit.cmd(f"surface {l2s(volumeIDs)} scheme tri")
            else:
                cubit.cmd(f"surface {l2s(volumeIDs)} scheme map")

            cubit.cmd(f"mesh surface {l2s(volumeIDs)}")

            fileName = wt_name + "-" + str(iStation) + "-t-0.in"

            if not os.path.exists(directory):
                os.makedirs(directory)

            if settings["make_input_for"] is not None:
                if "vabs" in settings["make_input_for"].lower():
                    writeVABSinput(
                        surfaceDict,
                        blade,
                        cs_params,
                        directory,
                        fileName,
                        volumeIDs,
                        materialsUsed,
                        crossSectionNormal,
                    )

            elif "anba" in settings["make_input_for"].lower():
                raise ValueError("ANBA currently not supported")
            else:
                raise NameError(
                    f'Unknown beam cross-sectional solver: {settings["make_input_for"]}'
                )

            if settings["export"] is not None:
                if (
                    "g" in settings["export"].lower()
                    or "cub" in settings["export"].lower()
                ):
                    if "g" in settings["export"].lower():
                        cubit.cmd(f'export mesh "{pathName}.g" overwrite')
                    if "cub" in settings["export"].lower():
                        cubit.cmd(f"delete curve {spanwiseMatOriCurve}")
                        cubit.cmd(f'save as "{pathName}-{str(iStation)}.cub" overwrite')
                elif len(settings["export"]) == 0:
                    pass
                else:
                    raise NameError(
                        f'Unknown model export format: {settings["export"]}'
                    )

            # Import all cross-sections into one cub file
            if settings["export"] is not None and "cub" in settings["export"].lower():
                cubit.cmd("reset ")
                writeSplineFromCoordinatePoints(cubit, refLineCoords)

                for iStation in stationList:
                    cubit.cmd(f'import cubit "{pathName}-{str(iStation)}.cub"')
                cubit.cmd(f'save as "{pathName}.cub" overwrite')

                # Remove unnecessary files to save space
                for filePath in glob.glob(f"{pathName}-*.cub"):
                    os.remove(filePath)
    return (
        cubit,
        blade,
        surfaceDict,
        birdsMouthVerts,
        iStationFirstWeb,
        iStationLastWeb,
        materialsUsed,
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
    materialsUsed: set
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
        surfaceDict,
        birdsMouthVerts,
        iStationFirstWeb,
        iStationLastWeb,
        materialsUsed,
        spanwiseMatOriCurve,
    ) = cubit_make_cross_sections(
        blade, wt_name, settings, cs_params, "3D", stationList
    )

    iStationStart = stationList[0]
    iStationEnd = stationList[-1]
    ### ### ### ###
    # Make volumes along the span.
    ### ### ### ###
    meshVolList = []

    partName = "shell"
    orderedList = getOrderedList(partName)
    if len(orderedList) > 0:
        meshVolList = makeAeroshell(surfaceDict, orderedList, meshVolList, iStationEnd)
    #     cubit.cmd(f'save as "python2.cub" overwrite')
    #     foo

    partName = "web"
    orderedList = getOrderedList(partName)
    orderedListWeb = orderedList.copy()
    if orderedList and len(orderedList[0]) > 1:
        meshVolList = makeAeroshell(surfaceDict, orderedList, meshVolList, iStationEnd)

    partName = "roundTEadhesive"
    orderedList = getOrderedList(partName)
    if orderedList and len(orderedList[0]) > 1:
        meshVolList = makeAeroshell(surfaceDict, orderedList, meshVolList, iStationEnd)

    partName = "flatTEadhesive"
    orderedList = getOrderedList(partName)

    if orderedList and len(orderedList[0]) > 1:
        meshVolList = makeAeroshell(surfaceDict, orderedList, meshVolList, iStationEnd)

    if (
        orderedListWeb
        and len(orderedListWeb[0]) > 1
        and cs_params["birds_mouth_amplitude_fraction"]
        and birdsMouthVerts
    ):
        makeBirdsMouth(
            blade,
            birdsMouthVerts,
            cs_params["birds_mouth_amplitude_fraction"],
            iStationFirstWeb,
            iStationLastWeb,
        )

    cubit.cmd(f"merge volume {l2s(meshVolList)}")
    cubit.cmd(f"reset volume all")

    cubit.cmd(f"delete surface with Is_Free")
    cubit.cmd("vol all size 0.2")
    # cubit.cmd(f'curve with name "layerThickness*" interval {cs_params["nel_per_layer"]}')
    cubit.cmd("set default autosize on")
    cubit.cmd(f"mesh volume {l2s(meshVolList)}")
    cubit.cmd(f"draw volume {l2s(meshVolList)}")

    # Blocks
    # for imat,material in enumerate(blade.materials):
    for imat, materialName in enumerate(materialsUsed):
        cubit.cmd(f'block {imat+1} add volume with name "*{materialName}*"')
        cubit.cmd(f'block {imat+1} name "{materialName}"')

    addColor(blade, "volume")

    # Adding Nodesets
    # Root Nodeset
    parseString = f'with name "*station{iStationStart}*"'
    print(f"parseString{parseString}")
    surfaceIDs = parse_cubit_list("surface", parseString)
    cubit.cmd(f"nodeset 1 add surface {l2s(surfaceIDs)} ")
    cubit.cmd(f'nodeset 1 name "root"')

    for iLoop, iStation in enumerate(stationList[1:-1]):
        parseString = f'with name "*station{iStation}*"'
        print(f"parseString{parseString}")
        surfaceIDs = parse_cubit_list("surface", parseString)
        cubit.cmd(f"nodeset {iLoop+2} add surface {l2s(surfaceIDs)} ")
        cubit.cmd(f'nodeset {iLoop+2} name "station{iStation}"')
    if not stationList[1:-1]:
        iLoop = -1
    # Tip Nodeset
    parseString = f'with name "*station{iStationEnd}*"'
    surfaceIDs = parse_cubit_list("surface", parseString)
    cubit.cmd(f"nodeset {iLoop+3} add surface {l2s(surfaceIDs)} ")
    cubit.cmd(f'nodeset {iLoop+3} name "tip"')

    # Outer mold-line nodeset
    cubit.cmd('draw surf with name "*layer0_bottomFace*"')
    parseString = f'with name "*layer0_bottomFace*"'
    surfaceIDs = parse_cubit_list("surface", parseString)
    cubit.cmd(f"nodeset {iLoop+4} add surface {l2s(surfaceIDs)} ")
    cubit.cmd(f'nodeset {iLoop+4} name "oml"')

    # ####################################
    # ### Assign material orientations ###
    # ####################################

    parseString = f'in volume with name "*volume*"'
    allVolumeIDs = parse_cubit_list("volume", parseString)
    theta2 = {}
    theta1 = {}
    theta3 = {}
    directions = {}

    for iVol, volumeID in enumerate(allVolumeIDs):
        surfIDforMatOri, sign = getMatOriSurface(volumeID, spanwiseMatOriCurve)

        for iEl, elementID in enumerate(get_volume_hexes(volumeID)):
            coords = cubit.get_center_point("hex", elementID)
            if elementID == 33:
                print()
            cubit.create_vertex(coords[0], coords[1], coords[2])
            iVert1 = get_last_id("vertex")
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

            # with open('spanwiseDirection.txt', 'a') as f:
            #     f.write(f' {spanwiseDirection[0]}, {spanwiseDirection[1]}, {spanwiseDirection[2]};')

            # with open('perimeterDirection.txt', 'a') as f:
            #     f.write(f' {perimeterDirection[0]}, {perimeterDirection[1]}, {perimeterDirection[2]};')
            # with open('surfaceNormal.txt', 'a') as f:
            #     f.write(f' {surfaceNormal[0]}, {surfaceNormal[1]}, {surfaceNormal[2]};')

            # print(f'iEl {iEl} elementID {elementID}')
            # print(f'spanwiseDirection: {spanwiseDirection}')
            # print(f'perimeterDirection: {perimeterDirection}')
            # print(f'surfaceNormal: {surfaceNormal}')

            dcm = getDCM(globalAxisBasisVectors, newCoordinateSystemVectors)
            directions[elementID] = newCoordinateSystemVectors
            theta1[elementID], theta2[elementID], theta3[elementID] = dcmToEulerAngles(
                dcm
            )

            # with open('theta.txt', 'a') as f:
            #     f.write(f' {theta2[elementID]}, {theta1[elementID]}, {theta3[elementID]};')

            # print(f'theta2: {theta2[elementID]}, theta1: {theta1[elementID]}, theta3: {theta3[elementID]}')

            length = 0.1
            cubit.create_vertex(
                coords[0] + length * perimeterDirection[0],
                coords[1] + length * perimeterDirection[1],
                coords[2] + length * perimeterDirection[2],
            )
            iVert2 = get_last_id("vertex")
            cubit.cmd(f"create curve vertex {iVert1} {iVert2}")
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
    return materialsUsed
