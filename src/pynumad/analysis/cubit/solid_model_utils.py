from cubit import *
from PyCubed_Main import *
from pynumad.analysis.cubit.utils import print_sine_curve_between_two_verts
import numpy as np
import re


def get_ordered_list(part_name):

    orderedList = []
    surfacesToConnect = [1]  # Initialize to enter loop
    iSurface = -1  # Initialize
    while surfacesToConnect:
        iSurface += 1
        parse_string = f'with name "*{part_name}*surface{iSurface+1}"'
        surfacesToConnect = parse_cubit_list("surface", parse_string)

        if surfacesToConnect:
            orderedList.append(surfacesToConnect)

    return orderedList


def make_spanwise_splines(surface_dict, orderedList):
    spanwiseSplines = []
    for alignedSurfaces in orderedList:
        tempList = []
        for i_point in range(4):
            vertex_list = []
            for index, iSurface in enumerate(alignedSurfaces):
                vertexID = surface_dict[iSurface]["verts"][i_point]
                vertex_list.append(vertexID)
                vertexName = cubit.get_entity_name("vertex", vertexID)

            curve = cubit.cmd(f"create curve spline vertex {l2s(vertex_list)}")
            tempList.append(get_last_id("curve"))
        spanwiseSplines.append(tempList)
    return spanwiseSplines


def make_a_volume(
    currentSurfaceID, nextSurfaceID, spanwiseSplinesForAvolume, surface_dict, i_station_end
):
    cubit.cmd(f"surface {currentSurfaceID} copy")
    currentSurfaceIDcopy = get_last_id("surface")

    cubit.cmd(f"surface {nextSurfaceID} copy")
    nextSurfaceIDcopy = get_last_id("surface")

    currentSurface = cubit.surface(currentSurfaceID)
    nextSurface = cubit.surface(nextSurfaceID)

    currentSurfaceCurves = surface_dict[currentSurfaceID]["curves"]
    nextSurfaceCurves = surface_dict[nextSurfaceID]["curves"]

    currentSurfaceVerteces = surface_dict[currentSurfaceID]["verts"]
    nextSurfaceVerteces = surface_dict[nextSurfaceID]["verts"]

    spanwiseSplinesForAvolume.append(
        spanwiseSplinesForAvolume[0]
    )  # Make list circle back

    transverseSurfaceIDs = []
    for iCurve in range(len(currentSurfaceCurves)):
        cubit.cmd(
            f"create surface curve {currentSurfaceCurves[iCurve]} {spanwiseSplinesForAvolume[iCurve]} {nextSurfaceCurves[iCurve]} {spanwiseSplinesForAvolume[iCurve+1]}"
        )
        transverseSurfaceIDs.append(get_last_id("surface"))

    surfName = cubit.get_entity_name("surface", currentSurface.id()).split("_")
    regex = re.compile("layer")
    layerName = [string for string in surfName if re.match(regex, string)][0]
    stringName = layerName + "_bottomFace"
    cubit.cmd(f'surface {transverseSurfaceIDs[0]} rename "{stringName}"')
    stringName = layerName + "_topFace"
    cubit.cmd(f'surface {transverseSurfaceIDs[2]} rename "{stringName}"')

    # cubit.cmd(f'save as "python1.cub" overwrite')
    # raise Exception(f'Volume "{volumeName}" creation failed')
    # Create Volume
    # n_start=get_last_id("volume")
    cubit.cmd(
        f"create volume surface {currentSurfaceIDcopy} {nextSurfaceIDcopy} {l2s(transverseSurfaceIDs)} noheal"
    )
    # n_end=get_last_id("volume")
    # print(f'n_start: {n_start}, n_end: {n_end}')

    if "Station" + str(i_station_end) in cubit.get_entity_name(
        "surface", nextSurfaceID
    ):  # This if statement is needed for componets that may have been droped between the last station and the second to last station
        volumeName = cubit.get_entity_name("surface", nextSurfaceID)
    else:
        volumeName = cubit.get_entity_name("surface", currentSurfaceID)
    if len(cubit.volume(get_last_id("volume")).surfaces()) < 6:
        print(
            f"\n\n ERROR with:\n\n create volume surface {currentSurfaceIDcopy} {nextSurfaceIDcopy} {l2s(transverseSurfaceIDs)} "
        )
        print(f"currentSurfaceIDcopy: {currentSurfaceIDcopy}")
        print(f"nextSurfaceIDcopy: {nextSurfaceIDcopy}")
        print(f"spanwiseSplinesForAvolume: {spanwiseSplinesForAvolume}")
        cubit.cmd(f'save as "python1.cub" overwrite')
        raise Exception(f'Volume "{volumeName}" creation failed')

    volumeName = volumeName.replace("surface", "volume")
    cubit.cmd(f'volume {get_last_id("volume")} rename "{volumeName}"')


def get_spanwise_splines_for_a_volume(
    iSpan, nCrossSections, spanwiseSplinesForOneSurface, nextSurfaceVerteces
):
    # Split off spanwise curves for a single volume and store them
    if iSpan < nCrossSections - 2:
        spanwiseSplinesForAvolume = []
        temp = []
        for iCurve, curve_id in enumerate(spanwiseSplinesForOneSurface):
            cubit.cmd(f"split curve {curve_id} at vertex {nextSurfaceVerteces[iCurve]}")
            temp.append(get_last_id("curve"))
            spanwiseSplinesForAvolume.append(get_last_id("curve") - 1)
        spanwiseSplinesForOneSurface = temp
    else:
        spanwiseSplinesForAvolume = spanwiseSplinesForOneSurface
    return spanwiseSplinesForAvolume, spanwiseSplinesForOneSurface


# def assignIntervals(volID,nIntervals):
#     thicknessCurveID=cubit.volume(volID).curves()[1].id()
#     #cubit.cmd(f'locate curve {thicknessCurveID} ')
#     cubit.cmd(f'curve {thicknessCurveID} interval {nIntervals}')


def make_all_volumes_for_a_part(surface_dict, orderedList, meshVolList, i_station_end):
    # nIntervals=3
    spanwiseSplines = make_spanwise_splines(surface_dict, orderedList)
    nCrossSections = len(orderedList[0])
    nPartSurfaceIDs = len(orderedList)
    if nCrossSections > 1:
        for iSpan in range(nCrossSections - 1):
            for partSurfaceIDs in range(nPartSurfaceIDs):
                currentSurfaceID = orderedList[partSurfaceIDs][iSpan]
                nextSurfaceID = orderedList[partSurfaceIDs][iSpan + 1]
                (
                    spanwiseSplinesForAvolume,
                    spanwiseSplines[partSurfaceIDs],
                ) = get_spanwise_splines_for_a_volume(
                    iSpan,
                    nCrossSections,
                    spanwiseSplines[partSurfaceIDs],
                    surface_dict[nextSurfaceID]["verts"],
                )
                make_a_volume(
                    currentSurfaceID,
                    nextSurfaceID,
                    spanwiseSplinesForAvolume,
                    surface_dict,
                    i_station_end,
                )
                meshVolList.append(get_last_id("volume"))
                # assignIntervals(get_last_id("volume"),nIntervals)
    else:
        raise ValueError("Can't make volumes with only one cross section.")

    return meshVolList


def verify_web_cutting_amplitude(
    blade, amplitude, tolerance, i_stationFirstWeb, i_stationLastWeb
):
    # Check to make sure that the amplitude does not result sharp volumes by cutting near a station location
    geometry = blade.geometry
    for i_stationCheck in range(i_stationFirstWeb + 1, i_stationLastWeb + 1):
        bladeSegmentLength = (
            geometry.ispan[i_stationCheck] - geometry.ispan[i_stationFirstWeb]
        )
        gap = bladeSegmentLength - amplitude
        # print(f'bladeSegmentLength: {bladeSegmentLength}\ngap {gap}')

        if abs(gap) > tolerance:
            break
        else:
            if gap > 0:
                amplitude = bladeSegmentLength - tolerance
            else:
                amplitude = bladeSegmentLength + tolerance
            break
    # print(f'new amplitude {amplitude} \nnew gap = {bladeSegmentLength-amplitude}')
    return amplitude


def make_birds_mouth(
    blade, birds_mouth_verts, birds_mouth_amplitude_fraction, i_stationFirstWeb, i_stationLastWeb
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
    straightLine = create_curve(v1, v2)

    amplitude = birds_mouth_amplitude_fraction * distance
    tolerance = distance * 0.05

    amplitude = verify_web_cutting_amplitude(
        blade, amplitude, tolerance, i_stationFirstWeb, i_stationLastWeb
    )

    curvedLine = cubit.curve(
        print_sine_curve_between_two_verts(v1.id(), v2.id(), amplitude, "z")
    )
    cubit.cmd(f"create surface skin curve {curvedLine.id()} {straightLine.id()}")
    baseSurface = get_last_id("surface")

    midPoint = list(curvedLine.position_from_fraction(0.5))
    tangent = straightLine.tangent(midPoint)

    # Get the cross-section normal
    parse_string = f'in surface with name "*webStation{i_stationFirstWeb}*"'
    surface_id = parse_cubit_list("surface", parse_string)[
        0
    ]  # Pick the first surface in this list since all on same plane
    coords = cubit.get_center_point("surface", surface_id)
    surfaceNormal = cubit.surface(surface_id).normal_at(coords)
    cutBlockLength = 5 * max(geometry.ichord)
    sweepDirection = np.array(vectNorm(crossProd(list(tangent), list(surfaceNormal))))

    cubit.cmd(
        f"sweep surface {baseSurface} direction {l2s(sweepDirection)} distance {cutBlockLength}"
    )
    cubit.cmd(
        f'move volume {get_last_id("volume")} x {-cutBlockLength/2*sweepDirection[0]} y {-cutBlockLength/2*sweepDirection[1]} z {-cutBlockLength/2*sweepDirection[2]}'
    )

    cuttingVolume = get_last_id("volume")

    parse_string = f'with name "*webStation*"'
    webVolumes = parse_cubit_list("volume", parse_string)

    cubit.cmd(f"subtract volume {cuttingVolume} from volume {l2s(webVolumes)}")

    return


# cubit.cmd('open "/home/ecamare/myprojects/bar/cubitDev/python/python0.cub"')


def get_approximate_thickness_direction_for_volume(volumeID):
    # This function is used when assigning material orientations

    # Get thickness direction tangents
    approximateThicknessDirection = []
    for currentCurve in cubit.volume(volumeID).curves():
        curveName = cubit.get_entity_name("curve", currentCurve.id())
        if "layer_thickness" in curveName:
            coords = currentCurve.position_from_fraction(0.5)
            approximateThicknessDirection.append(currentCurve.tangent(coords))
    approximateThicknessDirection = np.array(approximateThicknessDirection)
    nThicknessCurves, _ = approximateThicknessDirection.shape

    if nThicknessCurves == 4:  # All other cases
        return np.mean(approximateThicknessDirection, 0)
    elif nThicknessCurves == 8:  # LE adhesive case and round TE adhesive
        return 0
    elif nThicknessCurves == 6:  # Web overwrap
        # Take the mean of all curves with name 'layer_thickness'
        mean = np.mean(approximateThicknessDirection, 0)

        errorList = []
        for i in range(nThicknessCurves):
            diff = approximateThicknessDirection[i] - mean

            errorList.append(sqrt(dotProd(diff, diff)))
        sortIndex = np.argsort(errorList)[
            :4
        ]  # Take the first four. This discards the two directions with the largest deviation from the average

        return np.mean(approximateThicknessDirection[sortIndex, :], 0)
    else:
        cubit.cmd(f'save as "Debug.cub" overwrite')
        raise ValueError(
            f"The number of thickness curves in volume is unexpected. Cannot assign material orientation. nThicknessCurves: {nThicknessCurves}"
        )

    return


def get_mat_ori_surface(volumeID, spanwiseMatOriCurve):
    # This function is used when assigning material orientations
    # This gets returns the surface within a volume that will be used to get surface normals.
    # The sign +-1 is also returned since some of the surfaces are oriented the wrong way

    approximateThicknessDirection = get_approximate_thickness_direction_for_volume(volumeID)

    # Create a list of surface IDs in the given volume
    surface_ids = []
    volumeSurfaces = cubit.volume(volumeID).surfaces()
    for currentSurface in volumeSurfaces:
        surface_ids.append(currentSurface.id())

    # Eliminate surfaces that have two curves named thickness:
    surfaceCT = 0
    for currentSurface in volumeSurfaces:
        curveCT = (
            0  # Counts the number of curves in the surface with name 'layer_thickness'
        )
        for currentCurve in currentSurface.curves():
            curveName = cubit.get_entity_name("curve", currentCurve.id())
            if "layer_thickness" in curveName:
                curveCT += 1

        if curveCT >= 2:
            surfaceCT += 1
            surface_ids.remove(currentSurface.id())

    # surface_ids now has the list of surfaces w/o thickness curves
    if len(surface_ids) == 2 or len(surface_ids) == 1:
        if len(surface_ids) == 2:
            surfaceName = cubit.get_entity_name("surface", surface_ids[0])
            if "topFace" in surfaceName:
                surface_id = surface_ids[0]
            else:
                surface_id = surface_ids[-1]
        elif len(surface_ids) == 1:  # Web overwrap
            surface_id = surface_ids[0]

        coords = cubit.get_center_point("surface", surface_id)
        surfaceNormal = cubit.surface(surface_id).normal_at(coords)

        if dotProd(surfaceNormal, approximateThicknessDirection) > 0:
            sign = 1.0
        else:
            sign = -1.0
    elif (
        len(surface_ids) == 0
    ):  # LE adhesive and/or TE adhesive for round cross-sections
        # print(f'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~volumeID {volumeID}')
        surface_id = False
        sign = 1.0

    else:
        raise ValueError(
            "The number of thickness curves in volume is unexpected. Cannot assign material orientation"
        )

    return surface_id, sign
