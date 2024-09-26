import numpy as np
import warnings
from os import getcwd
from os.path import join
import subprocess
import copy as cp
from scipy import interpolate

import pynumad
from pynumad.utils.interpolation import interpolator_wrap

##from pynumad.mesh_gen.shellClasses import shellRegion, elementSet, NuMesh3D, spatialGridList2D, spatialGridList3D
from pynumad.mesh_gen.boundary2d import *
from pynumad.mesh_gen.surface import Surface
from pynumad.mesh_gen.mesh2d import *
from pynumad.mesh_gen.mesh3d import Mesh3D
from pynumad.mesh_gen.shell_region import ShellRegion
from pynumad.mesh_gen.mesh_tools import *
from pynumad.mesh_gen.element_utils import *
#from pynumad.analysis.ansys.write import writeAnsysShellModel


def shell_mesh_general(blade, forSolid, includeAdhesive, elementSize):
    """
    This method generates a finite element shell mesh for the blade, based on what is
    stored in blade.geometry.coordinates, blade.keypoints.key_points, 
    and blade.geometry.profiles.  Output is given as a python dictionary.

    Parameters
    -----------
    blade: Blade
    forSolid: bool
    includeAdhesive: bool
    elementSize: float

    Returns
    -------
    meshData:
    Nodes and elements for outer shell and shear webs:
    nodes:
        - [x, y, z]
        - [x, y, z]
        ...
    elements:
        - [n1,n2,n3,n4]
        - [n1,n2,n3,n4]
        ...
    Set list and section list for the outer shell and shear webs.
    These are companion lists with the same length and order,
    so meshData['sets']['element'][i] corresponds to meshData['sections'][i]
    sets:
        element:
            - name: set1Name
              labels: [e1, e2, e3 ...]
            - name: set2Name
              labels: [e1, e2, e3 ...]
            ...
    sections:
        - type: 'shell'
          layup:
             - [materialid,thickness,angle] ## layer 1
             - [materialid,thickness,angle] ## layer 2
             ...
          elementSet: set1Name
        - type: 'shell'
          layup:
             - [materialid,thickness,angle] ## layer 1
             - [materialid,thickness,angle] ## layer 2
             ...
          elementSet: set2Name
    Nodes, elements and set for adhesive elements
    adhesiveNds:
        - [x, y, z]
        ...
    adhesiveEls:
        - [n1,n2,n3,n4,n5,n6,n7,n8]
        ...
    adhesiveElSet:
        name: 'adhesiveElements'
        labels: [0,1,2,3....(number of adhesive elements)]
    """
    geometry = blade.geometry
    coordinates = geometry.coordinates
    profiles = geometry.profiles
    key_points = blade.keypoints.key_points
    stacks = blade.stackdb.stacks
    swstacks = blade.stackdb.swstacks
    
    geomSz = coordinates.shape
    lenGeom = geomSz[0]
    numXsec = geomSz[2]
    XSCurvePts = np.array([], dtype=int)

    ## Determine the key curve points along the OML at each cross section
    for i in range(numXsec):
        keyPts = np.array([0])
        minDist = 1
        lePt = 0
        for j in range(lenGeom):
            prof = profiles[j, :, i]
            mag = np.linalg.norm(prof)
            if mag < minDist:
                minDist = mag
                lePt = j

        for j in range(5):
            kpCrd = key_points[j, :, i]
            minDist = geometry.ichord[i]
            pti = 1
            for k in range(lePt):
                ptCrd = coordinates[k, :, i]
                vec = ptCrd - kpCrd
                mag = np.linalg.norm(vec)
                if mag < minDist:
                    minDist = mag
                    pti = k
            keyPts = np.concatenate((keyPts, [pti]))
            coordinates[pti, :, i] = np.array(kpCrd)

        keyPts = np.concatenate((keyPts, [lePt]))
        for j in range(5, 10):
            kpCrd = key_points[j, :, i]
            minDist = geometry.ichord[i]
            pti = 1
            for k in range(lePt, lenGeom):
                ptCrd = coordinates[k, :, i]
                vec = ptCrd - kpCrd
                mag = np.linalg.norm(vec)
                if mag < minDist:
                    minDist = mag
                    pti = k
            keyPts = np.concatenate((keyPts, [pti]))
            coordinates[pti, :, i] = np.array(kpCrd)

        keyPts = np.concatenate((keyPts, [lenGeom - 1]))
        allPts = np.array([keyPts[0]])
        for j in range(0, len(keyPts) - 1):
            secPts = np.linspace(keyPts[j], keyPts[j + 1], 4)
            secPts = np.round(secPts).astype(int)
            allPts = np.concatenate((allPts, secPts[1:4]))

        XSCurvePts = np.vstack((XSCurvePts, allPts)) if XSCurvePts.size else allPts
    rws, cls = XSCurvePts.shape

    ## Create longitudinal splines down the blade through each of the key X-section points

    splineX = coordinates[XSCurvePts[0, :], 0, 0]
    splineY = coordinates[XSCurvePts[0, :], 1, 0]
    splineZ = coordinates[XSCurvePts[0, :], 2, 0]
    for i in range(1, rws):
        Xrow = coordinates[XSCurvePts[i, :], 0, i]
        splineX = np.vstack((splineX, Xrow.T))
        Yrow = coordinates[XSCurvePts[i, :], 1, i]
        splineY = np.vstack((splineY, Yrow.T))
        Zrow = coordinates[XSCurvePts[i, :], 2, i]
        splineZ = np.vstack((splineZ, Zrow.T))

    spParam = np.transpose(np.linspace(0, 1, rws))
    nSpi = rws + 2 * (rws - 1)
    spParami = np.transpose(np.linspace(0, 1, nSpi))
    splineXi = interpolator_wrap(spParam, splineX[:, 0], spParami, "pchip")
    splineYi = interpolator_wrap(spParam, splineY[:, 0], spParami, "pchip")
    splineZi = interpolator_wrap(spParam, splineZ[:, 0], spParami, "pchip")
    for i in range(1, cls):
        splineXi = np.vstack(
            [splineXi, interpolator_wrap(spParam, splineX[:, i], spParami, "pchip")]
        )
        splineYi = np.vstack(
            [splineYi, interpolator_wrap(spParam, splineY[:, i], spParami, "pchip")]
        )
        splineZi = np.vstack(
            [splineZi, interpolator_wrap(spParam, splineZ[:, i], spParami, "pchip")]
        )
    splineXi = splineXi.T
    splineYi = splineYi.T
    splineZi = splineZi.T
    ## Determine the first spanwise section that needs adhesive
    if includeAdhesive == 1:
        stPt = 0
        frstXS = 0
        while frstXS == 0 and stPt < splineXi.shape[0]:
            v1x = splineXi[stPt, 6] - splineXi[stPt, 4]
            v1y = splineYi[stPt, 6] - splineYi[stPt, 4]
            v1z = splineZi[stPt, 6] - splineZi[stPt, 4]
            v2x = splineXi[stPt, 30] - splineXi[stPt, 32]
            v2y = splineYi[stPt, 30] - splineYi[stPt, 32]
            v2z = splineZi[stPt, 30] - splineZi[stPt, 32]
            mag1 = np.sqrt(v1x * v1x + v1y * v1y + v1z * v1z)
            mag2 = np.sqrt(v2x * v2x + v2y * v2y + v2z * v2z)
            dp = (1 / (mag1 * mag2)) * (v1x * v2x + v1y * v2y + v1z * v2z)
            if dp > 0.7071:
                frstXS = stPt
            stPt = stPt + 3

        if frstXS == 0:
            frstXS = splineXi.shape[0]
    else:
        frstXS = splineXi.shape[0]

    ## Generate the mesh using the splines as surface guides
    
    ## --------------------------TEMP
    ## secNel = [1,3,15,5,8,1,1,13,5,12,3,1]
    ## secNel = [1,4,10,10,8,2,2,8,10,10,4,1]
    ## --------------------------END TEMP
    
    bladeSurf = Surface()
    ## Outer shell sections
    outShES = set()
    secList = list()
    stPt = 0
    for i in range(rws - 1):
        # if stPt < frstXS:
        #     stSec = 0
        #     endSec = 11
        #     stSp = 0
        # else:
        #     stSec = 1
        #     endSec = 10
        #     stSp = 3
        stSec = 0
        endSec = 11
        stSp = 0
        for j in range(stSec, endSec + 1):
            shellKp = np.array(
                [
                    [splineXi[stPt, stSp], splineYi[stPt, stSp], splineZi[stPt, stSp]],
                    [
                        splineXi[stPt, stSp + 3],
                        splineYi[stPt, stSp + 3],
                        splineZi[stPt, stSp + 3],
                    ],
                    [
                        splineXi[stPt + 3, stSp + 3],
                        splineYi[stPt + 3, stSp + 3],
                        splineZi[stPt + 3, stSp + 3],
                    ],
                    [
                        splineXi[stPt + 3, stSp],
                        splineYi[stPt + 3, stSp],
                        splineZi[stPt + 3, stSp],
                    ],
                    [
                        splineXi[stPt, stSp + 1],
                        splineYi[stPt, stSp + 1],
                        splineZi[stPt, stSp + 1],
                    ],
                    [
                        splineXi[stPt, stSp + 2],
                        splineYi[stPt, stSp + 2],
                        splineZi[stPt, stSp + 2],
                    ],
                    [
                        splineXi[stPt + 1, stSp + 3],
                        splineYi[stPt + 1, stSp + 3],
                        splineZi[stPt + 1, stSp + 3],
                    ],
                    [
                        splineXi[stPt + 2, stSp + 3],
                        splineYi[stPt + 2, stSp + 3],
                        splineZi[stPt + 2, stSp + 3],
                    ],
                    [
                        splineXi[stPt + 3, stSp + 2],
                        splineYi[stPt + 3, stSp + 2],
                        splineZi[stPt + 3, stSp + 2],
                    ],
                    [
                        splineXi[stPt + 3, stSp + 1],
                        splineYi[stPt + 3, stSp + 1],
                        splineZi[stPt + 3, stSp + 1],
                    ],
                    [
                        splineXi[stPt + 2, stSp],
                        splineYi[stPt + 2, stSp],
                        splineZi[stPt + 2, stSp],
                    ],
                    [
                        splineXi[stPt + 1, stSp],
                        splineYi[stPt + 1, stSp],
                        splineZi[stPt + 1, stSp],
                    ],
                    [
                        splineXi[stPt + 1, stSp + 1],
                        splineYi[stPt + 1, stSp + 1],
                        splineZi[stPt + 1, stSp + 1],
                    ],
                    [
                        splineXi[stPt + 1, stSp + 2],
                        splineYi[stPt + 1, stSp + 2],
                        splineZi[stPt + 1, stSp + 2],
                    ],
                    [
                        splineXi[stPt + 2, stSp + 2],
                        splineYi[stPt + 2, stSp + 2],
                        splineZi[stPt + 2, stSp + 2],
                    ],
                    [
                        splineXi[stPt + 2, stSp + 1],
                        splineYi[stPt + 2, stSp + 1],
                        splineZi[stPt + 2, stSp + 1],
                    ],
                ]
            )
            # vec = shellKp[1, :] - shellKp[0, :]
            # mag = np.linalg.norm(vec)
            # nEl = np.array([], dtype=int)
            # nEl = np.concatenate([nEl, [np.ceil(mag / elementSize).astype(int)]])
            # vec = shellKp[2, :] - shellKp[1, :]
            # mag = np.linalg.norm(vec)
            # nEl = np.concatenate([nEl, [np.ceil(mag / elementSize).astype(int)]])
            # vec = shellKp[3, :] - shellKp[2, :]
            # mag = np.linalg.norm(vec)
            # nEl = np.concatenate([nEl, [np.ceil(mag / elementSize).astype(int)]])
            # vec = shellKp[0, :] - shellKp[3, :]
            # mag = np.linalg.norm(vec)
            # nEl = np.concatenate([nEl, [np.ceil(mag / elementSize).astype(int)]])
            
            vec = shellKp[1, :] - shellKp[0, :]
            mag = np.linalg.norm(vec)
            
            nEl1 = int(np.ceil(mag/elementSize))
            ##
            ## nEl1 = secNel[j]
            ##
            
            vec = shellKp[2, :] - shellKp[1, :]
            mag = np.linalg.norm(vec)
            nEl2 = int(np.ceil(mag/elementSize))
            vec = shellKp[3, :] - shellKp[2, :]
            mag = np.linalg.norm(vec)
            
            nEl3 = int(np.ceil(mag/elementSize))
            ##
            ## nEl3 = secNel[j]
            ##
            
            vec = shellKp[0, :] - shellKp[3, :]
            mag = np.linalg.norm(vec)
            nEl4 = int(np.ceil(mag/elementSize))
            nEl = np.array([nEl1,nEl2,nEl3,nEl4])
            
            bladeSurf.addShellRegion(
                "quad3",
                shellKp,
                nEl,
                name=stacks[j, i].name,
                elType="quad",
                meshMethod="structured",
            )
            outShES.add(stacks[j,i].name)
            newSec = dict()
            newSec["type"] = "shell"
            layup = list()
            for pg in stacks[j, i].plygroups:
                totThick = 0.001*pg.thickness * pg.nPlies
                ply = [pg.materialid, totThick, pg.angle]
                layup.append(ply)
            newSec["layup"] = layup
            newSec["elementSet"] = stacks[j, i].name
            newSec["xDir"] = (shellKp[3,:] - shellKp[0,:]) + (shellKp[2,:] - shellKp[1,:])
            newSec["xyDir"] = (shellKp[1,:] - shellKp[0,:]) + (shellKp[2,:] - shellKp[3,:])
            secList.append(newSec)
            stSp = stSp + 3
        stPt = stPt + 3

    ## Shift the appropriate splines if the mesh is for a solid model seed
    if forSolid == 1:
        caseIndex = np.array([[9, 27, 3], [12, 24, 3], [24, 12, 8], [27, 9, 8]])
        for i in range(caseIndex.shape[0]):
            spl = caseIndex[i, 0]
            tgtSp = caseIndex[i, 1]
            sec = caseIndex[i, 2]
            stPt = 0
            for j in range(rws - 1):
                totalThick = 0
                for k in range(3):
                    tpp = 0.001 * stacks[sec, j].plygroups[k].thickness
                    npls = stacks[sec, j].plygroups[k].nPlies
                    totalThick = totalThick + tpp * npls
                for k in range(3):
                    vx = splineXi[stPt, tgtSp] - splineXi[stPt, spl]
                    vy = splineYi[stPt, tgtSp] - splineYi[stPt, spl]
                    vz = splineZi[stPt, tgtSp] - splineZi[stPt, spl]
                    magInv = 1 / np.sqrt(vx * vx + vy * vy + vz * vz)
                    ux = magInv * vx
                    uy = magInv * vy
                    uz = magInv * vz
                    splineXi[stPt, spl] = splineXi[stPt, spl] + 0.1*totalThick * ux
                    splineYi[stPt, spl] = splineYi[stPt, spl] + 0.1*totalThick * uy
                    splineZi[stPt, spl] = splineZi[stPt, spl] + 0.1*totalThick * uz
                    stPt = stPt + 1

    ## Shear web sections
    swES = set()
    stPt = 0
    web1Sets = np.array([])
    web2Sets = np.array([])
    for i in range(rws - 1):
        if swstacks[0][i].plygroups:
            shellKp = np.zeros((16, 3))
            shellKp[0, :] = np.array(
                [splineXi[stPt, 12], splineYi[stPt, 12], splineZi[stPt, 12]]
            )
            shellKp[1, :] = np.array(
                [splineXi[stPt, 24], splineYi[stPt, 24], splineZi[stPt, 24]]
            )
            shellKp[2, :] = np.array(
                [splineXi[stPt + 3, 24], splineYi[stPt + 3, 24], splineZi[stPt + 3, 24]]
            )
            shellKp[3, :] = np.array(
                [splineXi[stPt + 3, 12], splineYi[stPt + 3, 12], splineZi[stPt + 3, 12]]
            )
            shellKp[6, :] = np.array(
                [splineXi[stPt + 1, 24], splineYi[stPt + 1, 24], splineZi[stPt + 1, 24]]
            )
            shellKp[7, :] = np.array(
                [splineXi[stPt + 2, 24], splineYi[stPt + 2, 24], splineZi[stPt + 2, 24]]
            )
            shellKp[10, :] = np.array(
                [splineXi[stPt + 2, 12], splineYi[stPt + 2, 12], splineZi[stPt + 2, 12]]
            )
            shellKp[11, :] = np.array(
                [splineXi[stPt + 1, 12], splineYi[stPt + 1, 12], splineZi[stPt + 1, 12]]
            )
            shellKp[4, :] = 0.6666 * shellKp[0, :] + 0.3333 * shellKp[1, :]
            shellKp[5, :] = 0.3333 * shellKp[0, :] + 0.6666 * shellKp[1, :]
            shellKp[8, :] = 0.6666 * shellKp[2, :] + 0.3333 * shellKp[3, :]
            shellKp[9, :] = 0.3333 * shellKp[2, :] + 0.6666 * shellKp[3, :]
            shellKp[12, :] = 0.6666 * shellKp[11, :] + 0.3333 * shellKp[6, :]
            shellKp[13, :] = 0.3333 * shellKp[11, :] + 0.6666 * shellKp[6, :]
            shellKp[14, :] = 0.6666 * shellKp[7, :] + 0.3333 * shellKp[10, :]
            shellKp[15, :] = 0.3333 * shellKp[7, :] + 0.6666 * shellKp[10, :]

            # vec = shellKp[1, :] - shellKp[0, :]
            # mag = np.linalg.norm(vec)

            # nEl = np.array([], dtype=int)
            # nEl = np.concatenate([nEl, [np.ceil(mag / elementSize).astype(int)]])
            # vec = shellKp[2, :] - shellKp[1, :]
            # mag = np.linalg.norm(vec)
            # nEl = np.concatenate([nEl, [np.ceil(mag / elementSize).astype(int)]])
            # vec = shellKp[3, :] - shellKp[2, :]
            # mag = np.linalg.norm(vec)
            # nEl = np.concatenate([nEl, [np.ceil(mag / elementSize).astype(int)]])
            # vec = shellKp[0, :] - shellKp[3, :]
            # mag = np.linalg.norm(vec)
            # nEl = np.concatenate([nEl, [np.ceil(mag / elementSize).astype(int)]])
            
            vec = shellKp[1, :] - shellKp[0, :]
            mag = np.linalg.norm(vec)
            
            nEl1 = int(np.ceil(mag/elementSize))
            ##
            ## nEl1 = 8
            ##
            
            vec = shellKp[2, :] - shellKp[1, :]
            mag = np.linalg.norm(vec)
            nEl2 = int(np.ceil(mag/elementSize))
            vec = shellKp[3, :] - shellKp[2, :]
            mag = np.linalg.norm(vec)
            
            nEl3 = int(np.ceil(mag/elementSize))
            ##
            ## nEl3 = 8
            ##
            
            vec = shellKp[0, :] - shellKp[3, :]
            mag = np.linalg.norm(vec)
            nEl4 = int(np.ceil(mag/elementSize))
            nEl = np.array([nEl1,nEl2,nEl3,nEl4])

            bladeSurf.addShellRegion(
                "quad3",
                shellKp,
                nEl,
                name=swstacks[0][i].name,
                elType="quad",
                meshMethod="structured",
            )
            swES.add(swstacks[0][i].name)
            newSec = dict()
            newSec["type"] = "shell"
            layup = list()
            for pg in swstacks[0][i].plygroups:
                totThick = 0.001*pg.thickness * pg.nPlies
                ply = [pg.materialid, totThick, pg.angle]
                layup.append(ply)
            newSec["layup"] = layup
            newSec["elementSet"] = swstacks[0][i].name
            newSec["xDir"] = np.array([0.0,0.0,1.0])
            newSec["xyDir"] = (shellKp[1,:] - shellKp[0,:]) + (shellKp[2,:] - shellKp[3,:])
            secList.append(newSec)
        if swstacks[1][i].plygroups:
            shellKp = np.zeros((16, 3))
            shellKp[0, :] = np.array(
                [splineXi[stPt, 27], splineYi[stPt, 27], splineZi[stPt, 27]]
            )
            shellKp[1, :] = np.array(
                [splineXi[stPt, 9], splineYi[stPt, 9], splineZi[stPt, 9]]
            )
            shellKp[2, :] = np.array(
                [splineXi[stPt + 3, 9], splineYi[stPt + 3, 9], splineZi[stPt + 3, 9]]
            )
            shellKp[3, :] = np.array(
                [splineXi[stPt + 3, 27], splineYi[stPt + 3, 27], splineZi[stPt + 3, 27]]
            )
            shellKp[6, :] = np.array(
                [splineXi[stPt + 1, 9], splineYi[stPt + 1, 9], splineZi[stPt + 1, 9]]
            )
            shellKp[7, :] = np.array(
                [splineXi[stPt + 2, 9], splineYi[stPt + 2, 9], splineZi[stPt + 2, 9]]
            )
            shellKp[10, :] = np.array(
                [splineXi[stPt + 2, 27], splineYi[stPt + 2, 27], splineZi[stPt + 2, 27]]
            )
            shellKp[11, :] = np.array(
                [splineXi[stPt + 1, 27], splineYi[stPt + 1, 27], splineZi[stPt + 1, 27]]
            )
            shellKp[4, :] = 0.6666 * shellKp[0, :] + 0.3333 * shellKp[1, :]
            shellKp[5, :] = 0.3333 * shellKp[0, :] + 0.6666 * shellKp[1, :]
            shellKp[8, :] = 0.6666 * shellKp[2, :] + 0.3333 * shellKp[3, :]
            shellKp[9, :] = 0.3333 * shellKp[2, :] + 0.6666 * shellKp[3, :]
            shellKp[12, :] = 0.6666 * shellKp[11, :] + 0.3333 * shellKp[6, :]
            shellKp[13, :] = 0.3333 * shellKp[11, :] + 0.6666 * shellKp[6, :]
            shellKp[14, :] = 0.6666 * shellKp[7, :] + 0.3333 * shellKp[10, :]
            shellKp[15, :] = 0.3333 * shellKp[7, :] + 0.6666 * shellKp[10, :]

            # vec = shellKp[1, :] - shellKp[0, :]
            # mag = np.linalg.norm(vec)

            # nEl = np.array([], dtype=int)
            # nEl = np.concatenate([nEl, [np.ceil(mag / elementSize).astype(int)]])
            # vec = shellKp[2, :] - shellKp[1, :]
            # mag = np.linalg.norm(vec)
            # nEl = np.concatenate([nEl, [np.ceil(mag / elementSize).astype(int)]])
            # vec = shellKp[3, :] - shellKp[2, :]
            # mag = np.linalg.norm(vec)
            # nEl = np.concatenate([nEl, [np.ceil(mag / elementSize).astype(int)]])
            # vec = shellKp[0, :] - shellKp[3, :]
            # mag = np.linalg.norm(vec)
            # nEl = np.concatenate([nEl, [np.ceil(mag / elementSize).astype(int)]])
            
            vec = shellKp[1, :] - shellKp[0, :]
            mag = np.linalg.norm(vec)
            
            nEl1 = int(np.ceil(mag/elementSize))
            ##
            ## nEl1 = 8
            ##
            
            vec = shellKp[2, :] - shellKp[1, :]
            mag = np.linalg.norm(vec)
            nEl2 = int(np.ceil(mag/elementSize))
            vec = shellKp[3, :] - shellKp[2, :]
            mag = np.linalg.norm(vec)
            
            nEl3 = int(np.ceil(mag/elementSize))
            ##
            ## nEl3 = 8
            ##
            
            vec = shellKp[0, :] - shellKp[3, :]
            mag = np.linalg.norm(vec)
            nEl4 = int(np.ceil(mag/elementSize))
            nEl = np.array([nEl1,nEl2,nEl3,nEl4])

            bladeSurf.addShellRegion(
                "quad3",
                shellKp,
                nEl,
                name=swstacks[1][i].name,
                elType="quad",
                meshMethod="structured",
            )
            swES.add(swstacks[1][i].name)
            newSec = dict()
            newSec["type"] = "shell"
            layup = list()
            for pg in swstacks[1][i].plygroups:
                totThick = 0.001*pg.thickness * pg.nPlies
                ply = [pg.materialid, totThick, pg.angle]
                layup.append(ply)
            newSec["layup"] = layup
            newSec["elementSet"] = swstacks[1][i].name
            newSec["xDir"] = np.array([0.0,0.0,1.0])
            newSec["xyDir"] = (shellKp[1,:] - shellKp[0,:]) + (shellKp[2,:] - shellKp[3,:])
            secList.append(newSec)
        stPt = stPt + 3

    ## Generate Shell mesh

    print('getting blade mesh')
    shellData = bladeSurf.getSurfaceMesh()
    shellData["sections"] = secList
    
    ## Get local direction cosine orientations for individual elements
    print('getting element orientations')
    nodes = shellData["nodes"]
    elements = shellData["elements"]
    numEls = len(shellData["elements"])
    elOri = np.zeros((numEls,9),dtype=float)
    
    esi = 0
    for sec in secList:
        dirCos = get_direction_cosines(sec["xDir"],sec["xyDir"])
        es = shellData["sets"]["element"][esi]
        for ei in es["labels"]:
            eX = list()
            eY = list()
            eZ = list()
            for ndi in elements[ei]:
                if(ndi > -1):
                    eX.append(nodes[ndi,0])
                    eY.append(nodes[ndi,1])
                    eZ.append(nodes[ndi,2])
            if(len(eX) == 3):
                elType = "shell3"
            else:
                elType = "shell4"
            elCrd = np.array([eX,eY,eZ])
            elDirCos = correct_orient(dirCos,elCrd,elType)
            elOri[ei,0:3] = elDirCos[0]
            elOri[ei,3:6] = elDirCos[1]
            elOri[ei,6:9] = elDirCos[2]
        esi = esi + 1
        
    shellData["elementOrientations"] = elOri
    
    ## Get all outer shell and all shear web element sets
    
    outerLab = list()
    swLab = list()
    for es in shellData["sets"]["element"]:
        nm = es["name"]
        if(nm in outShES):
            outerLab.extend(es["labels"])
        elif(nm in swES):
            swLab.extend(es["labels"])
    outerSet = dict()
    outerSet["name"] = "allOuterShellEls"
    outerSet["labels"] = outerLab
    shellData["sets"]["element"].append(outerSet)
    swSet = dict()
    swSet["name"] = "allShearWebEls"
    swSet["labels"] = swLab
    shellData["sets"]["element"].append(swSet)
    
    ## Get root (Zmin) node set
    minZ = np.min(splineZi)
    rootLabs = list()
    lab = 0
    for nd in nodes:
        if(np.abs(nd[2] - minZ) < 0.25*elementSize):
            rootLabs.append(lab)
        lab = lab + 1
    newSet = dict()
    newSet["name"] = "RootNodes"
    newSet["labels"] = rootLabs
    try:
        shellData["sets"]["node"].append(newSet)
    except:
        nodeSets = list()
        nodeSets.append(newSet)
        shellData["sets"]["node"] = nodeSets

    ## Generate mesh for trailing edge adhesive if requested
    print('getting adhesive mesh')
    if includeAdhesive == 1:
        stPt = frstXS
        v1x = splineXi[stPt, 6] - splineXi[stPt, 4]
        v1y = splineYi[stPt, 6] - splineYi[stPt, 4]
        v1z = splineZi[stPt, 6] - splineZi[stPt, 4]
        mag1 = np.sqrt(v1x * v1x + v1y * v1y + v1z * v1z)
        v2x = splineXi[stPt, 30] - splineXi[stPt, 32]
        v2y = splineYi[stPt, 30] - splineYi[stPt, 32]
        v2z = splineZi[stPt, 30] - splineZi[stPt, 32]
        mag2 = np.sqrt(v2x * v2x + v2y * v2y + v2z * v2z)
        v3x = splineXi[stPt, 6] - splineXi[stPt, 30]
        v3y = splineYi[stPt, 6] - splineYi[stPt, 30]
        v3z = splineZi[stPt, 6] - splineZi[stPt, 30]
        mag3 = np.sqrt(v3x * v3x + v3y * v3y + v3z * v3z)
        v4x = splineXi[stPt, 4] - splineXi[stPt, 32]
        v4y = splineYi[stPt, 4] - splineYi[stPt, 32]
        v4z = splineZi[stPt, 4] - splineZi[stPt, 32]
        mag4 = np.sqrt(v4x * v4x + v4y * v4y + v4z * v4z)
        nE1 = np.ceil(mag1 / elementSize).astype(int)
        nE2 = np.ceil(mag3 / elementSize).astype(int)
        nE3 = np.ceil(mag2 / elementSize).astype(int)
        nE4 = np.ceil(mag4 / elementSize).astype(int)
        nEl = np.array([nE1, nE2, nE3, nE4])
        gdLayer = 0
        sweepElements = []
        guideNds = []
        while stPt < splineXi.shape[0]:
            shellKp = np.zeros((9, 3))
            shellKp[0, :] = np.array(
                [splineXi[stPt, 4], splineYi[stPt, 4], splineZi[stPt, 4]]
            )
            shellKp[1, :] = np.array(
                [splineXi[stPt, 6], splineYi[stPt, 6], splineZi[stPt, 6]]
            )
            shellKp[2, :] = np.array(
                [splineXi[stPt, 30], splineYi[stPt, 30], splineZi[stPt, 30]]
            )
            shellKp[3, :] = np.array(
                [splineXi[stPt, 32], splineYi[stPt, 32], splineZi[stPt, 32]]
            )
            shellKp[4, :] = np.array(
                [splineXi[stPt, 5], splineYi[stPt, 5], splineZi[stPt, 5]]
            )
            shellKp[5, :] = 0.5 * shellKp[1, :] + 0.5 * shellKp[2, :]
            shellKp[6, :] = np.array(
                [splineXi[stPt, 31], splineYi[stPt, 31], splineZi[stPt, 31]]
            )
            shellKp[7, :] = 0.5 * shellKp[0, :] + 0.5 * shellKp[3, :]
            shellKp[8, :] = 0.5 * shellKp[4, :] + 0.5 * shellKp[6, :]
            sReg = ShellRegion("quad2", shellKp, nEl, elType="quad", meshMethod="free")
            regMesh = sReg.createShellMesh()

            if stPt == frstXS:
                adhesMesh = Mesh3D(regMesh["nodes"], regMesh["elements"])
            else:
                guideNds.append(regMesh["nodes"])
                layerSwEl = np.ceil(
                    (splineZi[stPt, 4] - splineZi[(stPt - 3), 4]) / elementSize
                ).astype(int)
                sweepElements.append(layerSwEl)
            stPt = stPt + 3

        adMeshData = adhesMesh.createSweptMesh(
            "toDestNodes", sweepElements, destNodes=guideNds, interpMethod="smooth"
        )
        shellData["adhesiveEls"] = adMeshData["elements"]
        shellData["adhesiveNds"] = adMeshData["nodes"]
        adEls = len(adMeshData["elements"])
        adhesSet = dict()
        adhesSet["name"] = "adhesiveElements"
        labList = list(range(0, adEls))
        adhesSet["labels"] = labList
        adMeshData["sets"] = {"element": [adhesSet], "node": []}
        shellData["adhesiveElSet"] = adhesSet
        
        print('getting constraints')
        nDir = np.array([0.0,1.0,0.0])
        adMeshData = get_surface_nodes(adMeshData,"adhesiveElements","LP_AdNodes",nDir,normTol=45.0)
        nDir = np.array([0.0,-1.0,0.0])
        adMeshData = get_surface_nodes(adMeshData,"adhesiveElements","HP_AdNodes",nDir,normTol=45.0)
        lpEls = list()
        hpEls = list()
        for es in shellData["sets"]["element"]:
            if("LP_TE_REINF" in es["name"]):
                lpEls.extend(es["labels"])
            elif("HP_TE_REINF" in es["name"]):
                hpEls.extend(es["labels"])
        lpSet = {"name": "LP_TE_REINF", "labels": lpEls}
        shellData = add_element_set(shellData,lpSet)
        hpSet = {"name": "HP_TE_REINF", "labels": hpEls}
        shellData = add_element_set(shellData,hpSet)
        constraints = tie_2_meshes_constraints(adMeshData,"LP_AdNodes",shellData,"LP_TE_REINF",0.5*elementSize)
        hpConst = tie_2_meshes_constraints(adMeshData, "HP_AdNodes", shellData, "HP_TE_REINF", 0.5*elementSize)
        constraints.extend(hpConst)
        shellData["constraints"] = constraints
        
    matList = list()
    for mn in blade.definition.materials:
        newMat = dict()
        mat = blade.definition.materials[mn]
        newMat['name'] = mat.name
        newMat['density'] = mat.density
        newMat['elastic'] = dict()
        newMat['elastic']['E'] = [mat.ex,mat.ey,mat.ez]
        if(mat.type == "isotropic"):
            nu = mat.prxy
            newMat['elastic']['nu'] = [nu,nu,nu]
        else:
            newMat['elastic']['nu'] = [mat.prxy,mat.prxz,mat.pryz]
        newMat['elastic']['G'] = [mat.gxy,mat.gxz,mat.gyz]
        matList.append(newMat)

    shellData['materials'] = matList

    return shellData



def solidMeshFromShell(blade, shellMesh, layerNumEls, elementSize):
    shNodes = shellMesh["nodes"]
    shElements = shellMesh["elements"]
    elSets = shellMesh["sets"]["element"]
    sectns = shellMesh["sections"]

    print('building solid mesh')
    ## Initialize 3D solid mesh from the shell mesh
    bladeMesh = Mesh3D(shNodes, shElements)
    ## Calculate unit normal vectors for all nodes
    numNds = len(shNodes)
    numShEls = len(shElements)
    nodeNorms = np.zeros((numNds, 3))
    for i in range(0, len(shElements)):
        n1 = shElements[i, 0]
        n2 = shElements[i, 1]
        n3 = shElements[i, 2]
        n4 = shElements[i, 3]
        if n4 == -1:
            v1 = shNodes[n3, :] - shNodes[n1, :]
            v2 = shNodes[n2, :] - shNodes[n1, :]
        else:
            v1 = shNodes[n4, :] - shNodes[n2, :]
            v2 = shNodes[n3, :] - shNodes[n1, :]
        v3x = v1[1] * v2[2] - v1[2] * v2[1]
        v3y = v1[2] * v2[0] - v1[0] * v2[2]
        v3z = v1[0] * v2[1] - v1[1] * v2[0]
        v3 = np.array([v3x, v3y, v3z])
        mag = np.linalg.norm(v3)
        uNorm = (1.0 / mag) * v3
        for j in range(4):
            nj = shElements[i, j]
            if nj != -1:
                nodeNorms[nj, :] = nodeNorms[nj, :] + uNorm

    for i in range(numNds):
        mag = np.linalg.norm(nodeNorms[i])
        nodeNorms[i] = (1.0 / mag) * nodeNorms[i]

    ## Extrude shell mesh into solid mesh
    if len(layerNumEls) == 0:
        layerNumEls = np.array([1, 1, 1])

    prevLayer = shNodes.copy()
    guideNds = list()
    for i in range(0, len(layerNumEls)):
        nodeDist = np.zeros(numNds)
        nodeHitCt = np.zeros(numNds, dtype=int)
        numSec, numStat = blade.stackdb.stacks.shape
        j = 0
        for sec in sectns:
            layerThick = sec["layup"][i][1]
            #layerThick = sectns[j]["layup"][i][1]
            for el in elSets[j]["labels"]:
                for nd in shElements[el]:
                    if nd != -1:
                        nodeDist[nd] = nodeDist[nd] + layerThick
                        nodeHitCt[nd] = nodeHitCt[nd] + 1
            j = j + 1
        newLayer = np.zeros((numNds, 3))
        for j in range(0, numNds):
            if nodeHitCt[j] != 0:
                nodeDist[j] = nodeDist[j] / nodeHitCt[j]
                newLayer[j] = prevLayer[j] + nodeDist[j] * nodeNorms[j]

        guideNds.append(newLayer)
        prevLayer = newLayer.copy()

    solidMesh = bladeMesh.createSweptMesh(
        sweepMethod="toDestNodes",
        sweepElements=layerNumEls,
        destNodes=guideNds,
        interpMethod="linear",
    )

    print('getting element sets')
    ## Construct the element set list, extrapolated from the shell model
    newSetList = list()
    newSectList = list()
    esi = 0
    for sec in sectns:
        es = elSets[esi]
        elArray = np.array(es["labels"])
        elLayer = 0
        li = 1
        for lne in layerNumEls:
            newSet = dict()
            newSet["name"] = es["name"] + "layer_" + str(li)
            newLabels = list()
            for i in range(0, lne):
                newLabels.extend(elArray + numShEls * elLayer)
                elLayer = elLayer + 1
            newSet["labels"] = newLabels
            newSetList.append(newSet)
            newSec = dict()
            newSec["type"] = "solid"
            newSec["elementSet"] = newSet["name"]
            newSec["material"] = sec["layup"][li - 1][0]
            newSec["xDir"] = sec["xDir"]
            newSec["xyDir"] = sec["xyDir"]
            newSectList.append(newSec)
            li = li + 1
        esi = esi + 1

    solidMesh["sets"] = dict()
    solidMesh["sets"]["element"] = newSetList
    solidMesh["sections"] = newSectList

    solidMesh["adhesiveNds"] = shellMesh["adhesiveNds"]
    solidMesh["adhesiveEls"] = shellMesh["adhesiveEls"]
    solidMesh["adhesiveElSet"] = shellMesh["adhesiveElSet"]
    
    numBdEls = len(solidMesh["elements"])
    totElLrs = sum(layerNumEls)
    
    ## Get orientations for all elements
    elOri = np.zeros((numBdEls,9),dtype=float)
    for li in range(0,totElLrs):
        shft = li*numShEls
        oi = 0
        for ori in shellMesh["elementOrientations"]:
            elOri[oi+shft] = ori
            oi = oi + 1
    solidMesh["elementOrientations"] = elOri
    
    ## Get extruded node and element sets
    extSets = get_extruded_sets(shellMesh, totElLrs)
    solidMesh["sets"]["node"] = extSets["node"]
    solidMesh["sets"]["element"].extend(extSets["element"])
    
    ## Get section node sets
    solidMesh = get_matching_node_sets(solidMesh)
    
    ## Get shear web edge set
    
    numBdNds = len(solidMesh["nodes"])
    ndElCt = np.zeros(numBdNds,dtype=int)
    for el in solidMesh["elements"]:
        for nd in el:
            if(nd > -1):
                ndElCt[nd] = ndElCt[nd] + 1
    
    edgeLab = list()
    for ns in solidMesh["sets"]["node"]:
        if(ns["name"] == "allShearWebEls"):
            for nd in ns["labels"]:
                if(ndElCt[nd] <= 2):
                    edgeLab.append(nd)
    
    newSet = dict()
    newSet["name"] = "shearWebEdges"
    newSet["labels"] = edgeLab
    
    solidMesh["sets"]["node"].append(newSet)
    
    print('getting shear web constraints.  This can take a few minutes...')
    ## Get tie constraints for shear webs
    
    constraints = shellMesh["constraints"]
    
    webConst = tie_2_sets_constraints(solidMesh, "shearWebEdges", "allOuterShellEls", 0.05*elementSize)
    
    constraints.extend(webConst)
    
    #print('getting adhesive constraints')
    ## Get tie constraints for adhesive
    
    # adMesh = dict()
    # adMesh["nodes"] = solidMesh["adhesiveNds"]
    # adMesh["elements"] = solidMesh["adhesiveEls"]
    # adConst = tie_2_meshes_constraints(adMesh,solidMesh,0.015*ndSp)
    
    #constraints.extend(adConst)
    
    solidMesh["constraints"] = constraints

    return solidMesh


def get_solid_mesh(blade, layerNumEls, elementSize):
    ## Edit stacks to be usable for 3D solid mesh
    blade.stackdb.edit_stacks_for_solid_mesh()
    ## Create shell mesh as seed
    ## Note the new output structure of shellMeshGeneral, as a single python dictionary  -E Anderson
    shellMesh = shell_mesh_general(blade, 1, 1, elementSize)
    print("finished shell mesh")
    solidMesh = solidMeshFromShell(blade, shellMesh, layerNumEls, elementSize)
    return solidMesh


def get_shell_mesh(blade, includeAdhesive, elementSize):
    meshData = shell_mesh_general(blade, 0, includeAdhesive, elementSize)
    return meshData

def get_root_mesh(axisPts, radius, thickness, elementSize, elLayers, config='inserts', boltRad=None, insertThk=None, adhesiveThk=None, numIns=None, tubeThk=None, tubeExtend=None, extNumEls=None):
    if(config == 'inserts'):
        ## Create insert cross section mesh
        thInc = 2.0*np.pi/numIns
        minR = radius[0] - thickness[0] ## inner radius of root
        midR = 0.5*(2.0*radius[0] - thickness[0]) ## mid-thickness radius of root
        
        bd = Boundary2D()
        x1 = (midR - boltRad)*np.cos(0.5*thInc)
        y1 = (midR - boltRad)*np.sin(0.5*thInc)
        x2 = (midR + boltRad)*np.cos(0.5*thInc)
        y2 = (midR + boltRad)*np.sin(0.5*thInc)
        kp = [[x1,y1],
              [x2,y2],
              [x1,y1]]
        nEls = int(np.ceil(2.0*np.pi*boltRad/elementSize))
        bd.addSegment('arc',kp,nEls)
        bdData = bd.getBoundaryMesh()
        
        mesher = Mesh2D(bdData['nodes'],bdData['elements'])
        boltMesh = mesher.createUnstructuredMesh('quad')
        boltMesh = make3D(boltMesh)
        # plotShellMesh(insMesh)
        
        ## Create insert cross section mesh
        
        insertRad = boltRad + insertThk
        x1 = (midR - insertRad)*np.cos(0.5*thInc)
        y1 = (midR - insertRad)*np.sin(0.5*thInc)
        x2 = (midR + insertRad)*np.cos(0.5*thInc)
        y2 = (midR + insertRad)*np.sin(0.5*thInc)
        kp = [[x1,y1],
              [x2,y2],
              [x1,y1]]
        bd.addSegment('arc',kp,nEls)
        bdData = bd.getBoundaryMesh()
        
        mesher = Mesh2D(bdData['nodes'],bdData['elements'])
        insMesh = mesher.createUnstructuredMesh('quad')
        insMesh = make3D(insMesh)
        # plotShellMesh(adhMesh)
        
        ## Create adhesive cross section mesh
        
        
        bd = Boundary2D()
        x1 = (midR - insertRad)*np.cos(0.5*thInc)
        y1 = (midR - insertRad)*np.sin(0.5*thInc)
        x2 = (midR + insertRad)*np.cos(0.5*thInc)
        y2 = (midR + insertRad)*np.sin(0.5*thInc)
        kp = [[x1,y1],
              [x2,y2],
              [x1,y1]]
        bd.addSegment('arc',kp,nEls)
        
        adhesiveRad = boltRad + insertThk + adhesiveThk
        x1 = (midR - adhesiveRad)*np.cos(0.5*thInc)
        y1 = (midR - adhesiveRad)*np.sin(0.5*thInc)
        x2 = (midR + adhesiveRad)*np.cos(0.5*thInc)
        y2 = (midR + adhesiveRad)*np.sin(0.5*thInc)
        kp = [[x1,y1],
              [x2,y2],
              [x1,y1]]
        bd.addSegment('arc',kp,nEls)
        bdData = bd.getBoundaryMesh()
        
        mesher = Mesh2D(bdData['nodes'],bdData['elements'])
        adhMesh = mesher.createUnstructuredMesh('quad')
        adhMesh = make3D(adhMesh)

        ## Create root fill material unit cell mesh
        
        bd = Boundary2D()
        x1 = (midR - adhesiveRad)*np.cos(0.5*thInc)
        y1 = (midR - adhesiveRad)*np.sin(0.5*thInc)
        x2 = (midR + adhesiveRad)*np.cos(0.5*thInc)
        y2 = (midR + adhesiveRad)*np.sin(0.5*thInc)
        kp = [[x1,y1],
              [x2,y2],
              [x1,y1]]
        bd.addSegment('arc',kp,nEls)
        
        kp = [[minR,0.],
              [radius[0],0.]]
        nEls = int(np.ceil(thickness[0]/elementSize))
        bd.addSegment('line',kp,nEls)
        
        x1 = radius[0]
        y1 = 0.
        x2 = radius[0]*np.cos(0.5*thInc)
        y2 = radius[0]*np.sin(0.5*thInc)
        x3 = radius[0]*np.cos(thInc)
        y3 = radius[0]*np.sin(thInc)
        kp = [[x1,y1],
              [x2,y2],
              [x3,y3]]
        nEls = int(np.ceil(radius[0]*thInc/elementSize))
        bd.addSegment('arc',kp,nEls)
        
        x1 = radius[0]*np.cos(thInc)
        y1 = radius[0]*np.sin(thInc)
        x2 = minR*np.cos(thInc)
        y2 = minR*np.sin(thInc)
        kp = [[x1,y1],
              [x2,y2]]
        nEls = int(np.ceil(thickness[0]/elementSize))
        bd.addSegment('line',kp,nEls)
        
        x1 = minR
        y1 = 0.
        x2 = minR*np.cos(0.5*thInc)
        y2 = minR*np.sin(0.5*thInc)
        x3 = minR*np.cos(thInc)
        y3 = minR*np.sin(thInc)
        kp = [[x1,y1],
              [x2,y2],
              [x3,y3]]
        nEls = int(np.ceil(radius[0]*thInc/elementSize))
        bd.addSegment('arc',kp,nEls)
        
        bdData = bd.getBoundaryMesh()
        
        mesher = Mesh2D(bdData['nodes'],bdData['elements'])
        fillMesh = mesher.createUnstructuredMesh('quad')
        fillMesh = make3D(fillMesh)
        # plotShellMesh(fillMesh)
        
        ## Create mesh of whole root cross section
        faceSurf = Surface()
        
        boltSets = list()
        insSets = list()
        adhSets = list()
        fillSets = list()
        for i in range(0,numIns):
            rotAng = thInc*i*180.0/np.pi
            newBolt = cp.deepcopy(boltMesh)
            newBolt = rotate_mesh(newBolt,[0.,0.,0.],[0.,0.,1.],rotAng)
            nm = 'bolt_' + str(i)
            faceSurf.addMesh(newBolt,name=nm)
            boltSets.append(nm)
            newIns = cp.deepcopy(insMesh)
            newIns = rotate_mesh(newIns,[0.,0.,0.],[0.,0.,1.],rotAng)
            nm = 'insert_' + str(i)
            faceSurf.addMesh(newIns,name=nm)
            insSets.append(nm)
            newAd = cp.deepcopy(adhMesh)
            newAd = rotate_mesh(newAd,[0.,0.,0.],[0.,0.,1.],rotAng)
            nm = 'adhesive_' + str(i)
            faceSurf.addMesh(newAd,name=nm)
            adhSets.append(nm)
            newFill = cp.deepcopy(fillMesh)
            newFill = rotate_mesh(newFill,[0.,0.,0.],[0.,0.,1.],rotAng)
            nm = 'fill_' + str(i)
            faceSurf.addMesh(newFill,name=nm)
            fillSets.append(nm)
            
        faceMesh = faceSurf.getSurfaceMesh()
        # plotShellMesh(faceMesh)
        
        tVec = [0.,0.,axisPts[0]]
        faceMesh = translate_mesh(faceMesh,tVec)
        
        nPts = len(axisPts)
        allZpts = np.linspace(axisPts[0],axisPts[nPts-1],elLayers+1)
        radFun = interpolate.interp1d(axisPts,radius,kind='linear')
        thkFun = interpolate.interp1d(axisPts,thickness,kind='linear')
        ptRad = radFun(allZpts)
        ptThk = thkFun(allZpts)
        
        rootMesher = Mesh3D(faceMesh['nodes'],faceMesh['elements'])
        
        gdNds = list()
        nEls = list()
        for i in range(1,len(allZpts)):
            newNds = faceMesh['nodes'].copy()
            thkFact = ptThk[i]/ptThk[0]
            shftDist = ptRad[i] - ptRad[0]*thkFact
            for j, nd in enumerate(newNds):
                mag = np.linalg.norm(nd[0:2])
                unitXY = (1.0/mag)*nd[0:2]
                newXY = thkFact*nd[0:2] + shftDist*unitXY
                newNds[j,0:2] = newXY
                newNds[j,2] = allZpts[i]
            gdNds.append(newNds)
            nEls.append(1)
            
        rootMesh = rootMesher.createSweptMesh('toDestNodes',nEls,destNodes=gdNds,interpMethod='linear')
        
        # rootMesher = Mesh3D(faceMesh['nodes'],faceMesh['elements'])
        # sD = axisPts[nPts-1] - axisPts[0]
        # rootMesh = rootMesher.createSweptMesh('inDirection',elLayers,sweepDistance=sD,axis=[0.,0.,1.])
        
        # plotSolidMesh(rootMesh)
        
        extSets = get_extruded_sets(faceMesh,elLayers)
        
        rootMesh['sets'] = extSets
        
        rootMesh = get_element_set_union(rootMesh,boltSets,'allBolts')
        rootMesh = get_element_set_union(rootMesh,insSets,'allInsert')
        rootMesh = get_element_set_union(rootMesh,adhSets,'allAdhesive')
        rootMesh = get_element_set_union(rootMesh,fillSets,'allFill')
        
        rootMesh = get_matching_node_sets(rootMesh)
        
        sections = list()
        
        newSec = dict()
        newSec['type'] = 'solid'
        newSec['elementSet'] = 'allBolts'
        newSec['material'] = 'boltMat'
        newSec['xDir'] = [0.,0.,1.]
        newSec['xyDir'] = [1.,0.,0.]
        sections.append(newSec)
        
        newSec = dict()
        newSec['type'] = 'solid'
        newSec['elementSet'] = 'allInsert'
        newSec['material'] = 'insertMat'
        newSec['xDir'] = [0.,0.,1.]
        newSec['xyDir'] = [1.,0.,0.]
        sections.append(newSec)
        
        newSec = dict()
        newSec['type'] = 'solid'
        newSec['elementSet'] = 'allAdhesive'
        newSec['material'] = 'adhesiveMat'
        newSec['xDir'] = [0.,0.,1.]
        newSec['xyDir'] = [1.,0.,0.]
        sections.append(newSec)
        
        newSec = dict()
        newSec['type'] = 'solid'
        newSec['elementSet'] = 'allFill'
        newSec['material'] = 'fillMat'
        newSec['xDir'] = [0.,0.,1.]
        newSec['xyDir'] = [1.,0.,0.]
        sections.append(newSec)
        
        rootMesh['sections'] = sections
        
        return rootMesh
        
    else:
        ## Create main root
        circEls = int(2.0*np.pi*radius[0]/elementSize)
        minR = radius[0] - thickness[0]
        hFillThk = 0.5*(thickness[0] - 2.0*adhesiveThk - tubeThk)
        faceSurf = Surface()
        
        ## Inner fill
        print('Inner fill')
        
        bd = Boundary2D()
        xi = 0.0
        
        y1 = -minR
        y2 = minR
        kp = [[xi,y1],
              [xi,y2],
              [xi,y1]]
        bd.addSegment('arc',kp,circEls)
        
        # y1 = -minR - hFillThk
        # y2 = minR + hFillThk
        # kp = [[xi,y1],
              # [xi,y2],
              # [xi,y1]]
        # bd.addSegment('arc',kp,circEls)
        
        bdData = bd.getBoundaryMesh()
        mesher = Mesh2D(bdData['nodes'],bdData['elements'])
        #layerMesh = mesher.createUnstructuredMesh('quad')
        sEls = int(np.ceil(hFillThk/elementSize))
        layerMesh = mesher.createSweptMesh('fromPoint',sEls,sweepDistance=hFillThk,point=[0.,0.])
        layerMesh = make3D(layerMesh)
        faceSurf.addMesh(layerMesh,name='innerFill')
        
        ## Inner adhesive
        print('Inner adhesive')
        
        bd = Boundary2D()
        xi = 0.0
        
        y1 = -minR - hFillThk
        y2 = minR + hFillThk
        kp = [[xi,y1],
              [xi,y2],
              [xi,y1]]
        bd.addSegment('arc',kp,circEls)
        
        # y1 = -minR - hFillThk - adhesiveThk
        # y2 = minR + hFillThk + adhesiveThk
        # kp = [[xi,y1],
              # [xi,y2],
              # [xi,y1]]
        # bd.addSegment('arc',kp,circEls)
        
        bdData = bd.getBoundaryMesh()
        mesher = Mesh2D(bdData['nodes'],bdData['elements'])
        #layerMesh = mesher.createUnstructuredMesh('quad')
        sEls = int(np.ceil(adhesiveThk/elementSize))
        layerMesh = mesher.createSweptMesh('fromPoint',sEls,sweepDistance=adhesiveThk,point=[0.,0.])
        layerMesh = make3D(layerMesh)
        faceSurf.addMesh(layerMesh,name='innerAdhesive')
        
        ## Tube
        print('Tube')
        
        bd = Boundary2D()
        xi = 0.0
        
        y1 = -minR - hFillThk - adhesiveThk
        y2 = minR + hFillThk + adhesiveThk
        kp = [[xi,y1],
              [xi,y2],
              [xi,y1]]
        bd.addSegment('arc',kp,circEls)
        
        # y1 = -minR - hFillThk - adhesiveThk - tubeThk
        # y2 = minR + hFillThk + adhesiveThk + tubeThk
        # kp = [[xi,y1],
              # [xi,y2],
              # [xi,y1]]
        # bd.addSegment('arc',kp,circEls)
        
        bdData = bd.getBoundaryMesh()
        mesher = Mesh2D(bdData['nodes'],bdData['elements'])
        #layerMesh = mesher.createUnstructuredMesh('quad')
        sEls = int(np.ceil(tubeThk/elementSize))
        layerMesh = mesher.createSweptMesh('fromPoint',sEls,sweepDistance=tubeThk,point=[0.,0.])
        layerMesh = make3D(layerMesh)
        faceSurf.addMesh(layerMesh,name='tube')
        
        ## Outer adhesive
        print('Outer adhesive')
        
        bd = Boundary2D()
        xi = 0.0
        
        y1 = -radius[0] + hFillThk + adhesiveThk
        y2 = -y1
        kp = [[xi,y1],
              [xi,y2],
              [xi,y1]]
        bd.addSegment('arc',kp,circEls)
        
        # y1 = -radius[0] + hFillThk
        # y2 = -y1
        # kp = [[xi,y1],
              # [xi,y2],
              # [xi,y1]]
        # bd.addSegment('arc',kp,circEls)
        
        bdData = bd.getBoundaryMesh()
        mesher = Mesh2D(bdData['nodes'],bdData['elements'])
        #layerMesh = mesher.createUnstructuredMesh('quad')
        sEls = int(np.ceil(adhesiveThk/elementSize))
        layerMesh = mesher.createSweptMesh('fromPoint',sEls,sweepDistance=adhesiveThk,point=[0.,0.])
        layerMesh = make3D(layerMesh)
        faceSurf.addMesh(layerMesh,name='outerAdhesive')
        
        ## Outer fill
        print('Outer fill')
        
        bd = Boundary2D()
        xi = 0.0
        
        y1 = -radius[0] + hFillThk
        y2 = -y1
        kp = [[xi,y1],
              [xi,y2],
              [xi,y1]]
        bd.addSegment('arc',kp,circEls)
        
        # y1 = -radius[0]
        # y2 = -y1
        # kp = [[xi,y1],
              # [xi,y2],
              # [xi,y1]]
        # bd.addSegment('arc',kp,circEls)
        
        bdData = bd.getBoundaryMesh()
        mesher = Mesh2D(bdData['nodes'],bdData['elements'])
        #layerMesh = mesher.createUnstructuredMesh('quad')
        sEls = int(np.ceil(hFillThk/elementSize))
        layerMesh = mesher.createSweptMesh('fromPoint',sEls,sweepDistance=hFillThk,point=[0.,0.])
        layerMesh = make3D(layerMesh)
        faceSurf.addMesh(layerMesh,name='outerFill')
        
        ## Build 3D mesh
        faceMesh = faceSurf.getSurfaceMesh()
        
        tVec = [0.,0.,axisPts[0]]
        faceMesh = translate_mesh(faceMesh,tVec)
        
        nPts = len(axisPts)
        allZpts = np.linspace(axisPts[0],axisPts[nPts-1],elLayers+1)
        radFun = interpolate.interp1d(axisPts,radius,kind='linear')
        thkFun = interpolate.interp1d(axisPts,thickness,kind='linear')
        ptRad = radFun(allZpts)
        ptThk = thkFun(allZpts)
        
        rootMesher = Mesh3D(faceMesh['nodes'],faceMesh['elements'])
        
        gdNds = list()
        nEls = list()
        for i in range(1,len(allZpts)):
            newNds = faceMesh['nodes'].copy()
            thkFact = ptThk[i]/ptThk[0]
            shftDist = ptRad[i] - ptRad[0]*thkFact
            for j, nd in enumerate(newNds):
                mag = np.linalg.norm(nd[0:2])
                unitXY = (1.0/mag)*nd[0:2]
                newXY = thkFact*nd[0:2] + shftDist*unitXY
                newNds[j,0:2] = newXY
                newNds[j,2] = allZpts[i]
            gdNds.append(newNds)
            nEls.append(1)
            
        rootMesh = rootMesher.createSweptMesh('toDestNodes',nEls,destNodes=gdNds,interpMethod='linear')
        
        extSets = get_extruded_sets(faceMesh,len(nEls))
        
        rootMesh['sets'] = extSets
        
        ## Build tube extension
        print('extension')
        
        faceSurf = Surface()
        bd = Boundary2D()
        xi = 0.0
        
        y1 = -minR - hFillThk - adhesiveThk
        y2 = minR + hFillThk + adhesiveThk
        kp = [[xi,y1],
              [xi,y2],
              [xi,y1]]
        bd.addSegment('arc',kp,circEls)
        
        # y1 = -minR - hFillThk - adhesiveThk - tubeThk
        # y2 = minR + hFillThk + adhesiveThk + tubeThk
        # kp = [[xi,y1],
              # [xi,y2],
              # [xi,y1]]
        # bd.addSegment('arc',kp,circEls)
        
        bdData = bd.getBoundaryMesh()
        mesher = Mesh2D(bdData['nodes'],bdData['elements'])
        #layerMesh = mesher.createUnstructuredMesh('quad')
        sEls = int(np.ceil(tubeThk/elementSize))
        layerMesh = mesher.createSweptMesh('fromPoint',sEls,sweepDistance=tubeThk,point=[0.,0.])
        layerMesh = make3D(layerMesh)
        faceSurf.addMesh(layerMesh,name='tubeExtension')
        
        faceMesh = faceSurf.getSurfaceMesh()
        
        tVec = [0.,0.,axisPts[0]]
        faceMesh = translate_mesh(faceMesh,tVec)
        
        extMesher = Mesh3D(layerMesh['nodes'],layerMesh['elements'])
        extMesh = extMesher.createSweptMesh('inDirection',extNumEls,sweepDistance=tubeExtend,axis=[0.,0.,-1.])
        
        extSets = get_extruded_sets(faceMesh,extNumEls)
        
        extMesh['sets'] = extSets
        
        fullMesh = merge_meshes(rootMesh,extMesh)
        
        sections = list()
        
        newSec = dict()
        newSec['type'] = 'solid'
        newSec['elementSet'] = 'innerFill'
        newSec['material'] = 'fillMat'
        newSec['xDir'] = [0.,0.,1.]
        newSec['xyDir'] = [1.,0.,0.]
        sections.append(newSec)
        
        newSec = dict()
        newSec['type'] = 'solid'
        newSec['elementSet'] = 'innerAdhesive'
        newSec['material'] = 'adhesiveMat'
        newSec['xDir'] = [0.,0.,1.]
        newSec['xyDir'] = [1.,0.,0.]
        sections.append(newSec)
        
        newSec = dict()
        newSec['type'] = 'solid'
        newSec['elementSet'] = 'tube'
        newSec['material'] = 'tubeMat'
        newSec['xDir'] = [0.,0.,1.]
        newSec['xyDir'] = [1.,0.,0.]
        sections.append(newSec)
        
        newSec = dict()
        newSec['type'] = 'solid'
        newSec['elementSet'] = 'outerAdhesive'
        newSec['material'] = 'adhesiveMat'
        newSec['xDir'] = [0.,0.,1.]
        newSec['xyDir'] = [1.,0.,0.]
        sections.append(newSec)
        
        newSec = dict()
        newSec['type'] = 'solid'
        newSec['elementSet'] = 'outerFill'
        newSec['material'] = 'fillMat'
        newSec['xDir'] = [0.,0.,1.]
        newSec['xyDir'] = [1.,0.,0.]
        sections.append(newSec)
        
        newSec = dict()
        newSec['type'] = 'solid'
        newSec['elementSet'] = 'tubeExtension'
        newSec['material'] = 'tubeMat'
        newSec['xDir'] = [0.,0.,1.]
        newSec['xyDir'] = [1.,0.,0.]
        sections.append(newSec)
        
        fullMesh['sections'] = sections
        
        return fullMesh