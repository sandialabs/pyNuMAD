import numpy as np
import plotly.graph_objects as go
from pynumad.mesh_gen.spatial_grid_list2d import *
from pynumad.mesh_gen.spatial_grid_list3d import *
from pynumad.mesh_gen.element_utils import *

def rotate_vector(vec,axis,angle):
    if(angle < 0.0000000001):
        return vec.copy()
    else:
        axAr = np.array(axis)
        mag = np.linalg.norm(axis)
        unitAxis = (1.0/mag)*axAr
        alp1 = np.zeros((3,3),dtype=float)
        alp1[0] = unitAxis
        i1 = 0
        if(abs(unitAxis[1]) < abs(unitAxis[0])):
            i1 = 1
        if(abs(unitAxis[2]) < abs(unitAxis[i1])):
            i1 = 2
        alp1[1,i1] = np.sqrt(1.0 - alp1[0,i1]*alp1[0,i1])
        for i2 in range(0,3):
            if(i2 != i1):
                alp1[1,i2] = -alp1[0,i1]*alp1[0,i2]/alp1[1,i1]
        alp1[2] = cross_prod(alp1[0], alp1[1])
        theta = angle*np.pi/180.0
        cs = np.cos(theta)
        sn = np.sin(theta)
        alp2 = np.array([[1.0,0.0,0.0],
                         [0.0,cs,-sn],
                         [0.0,sn,cs]])
        rV = np.matmul(alp1,vec)
        rV = np.matmul(alp2,rV)
        rV = np.matmul(rV,alp1)
        return rV

def translate_mesh(meshData,tVec):
    tAr = np.array(tVec)
    nLen = len(meshData['nodes'])
    newNds = np.zeros((nLen,3),dtype=float)
    for i, nd in enumerate(meshData['nodes']):
        newNds[i] = nd + tAr
    meshData['nodes'] = newNds
    return meshData

def rotate_mesh(meshData,pt,axis,angle):
    ptAr = np.array(pt)
    nLen = len(meshData['nodes'])
    newNds = np.zeros((nLen,3),dtype=float)
    for i, nd in enumerate(meshData['nodes']):
        tCrd = nd - ptAr
        rCrd = rotate_vector(tCrd,axis,angle)
        newNds[i] = ptAr + rCrd
    meshData['nodes'] = newNds
    return meshData

def get_direction_cosines(xDir,xyDir):
    mag = np.linalg.norm(xDir)
    a1 = (1.0/mag)*xDir
    zDir = cross_prod(xDir,xyDir)
    mag = np.linalg.norm(zDir)
    a3 = (1.0/mag)*zDir
    a2 = cross_prod(a3,a1)
    dirCos = np.array([a1,a2,a3])
    return dirCos

def get_average_node_spacing(nodes,elements):
    totDist = 0.0
    ct = 0
    for el in elements:
        for ndi in el:
            if(ndi > -1):
                nd1 = nodes[ndi]
                for ndi2 in el:
                    if(ndi2 > -1 and ndi2 != ndi):
                        nd2 = nodes[ndi2]
                        vec = nd1 - nd2
                        dist = np.linalg.norm(vec)
                        totDist = totDist + dist
                        ct = ct + 1
    return totDist/ct

def check_all_jacobians(nodes,elements):
    failedEls = set()
    ei = 0
    for el in elements:
        xC = []
        yC = []
        zC = []
        for nd in el:
            if(nd > -1):
                xC.append(nodes[nd,0])
                yC.append(nodes[nd,1])
                zC.append(nodes[nd,2])
        elCrd = np.array([xC,yC,zC])
        nn = len(xC)
        if(nn == 8):
            elType = 'brick8'
        elif(nn == 6):
            elType = 'wedge6'
        else:
            elType = ''
        passed = check_jacobian(elCrd,elType)
        if(not passed):
            failedEls.add(ei)
        ei = ei + 1
    return failedEls

def get_element_volumes(meshData):
    elVols = dict()
    elMats = dict()
    nodes = meshData['nodes']
    elements = meshData['elements']
    elSets = meshData['sets']['element']
    if(len(elements[0]) > 4):
        for i, sec in enumerate(meshData['sections']):
            for ei in elSets[i]['labels']:
                elCrd = get_el_coord(elements[ei],nodes)
                if(elements[ei][6] == -1):
                    elType = 'wedge6'
                else:
                    elType = 'brick8'
                eVol = get_volume(elCrd,elType)
                stei = str(ei+1)
                elVols[stei] = eVol
                elMats[stei] = sec['material']
    else:
        for i, sec in enumerate(meshData['sections']):
            for ei in elSets[i]['labels']:
                elCrd = get_el_coord(elements[ei],nodes)
                if(elements[ei][3] == -1):
                    elType = 'shell3'
                else:
                    elType = 'shell4'
                eVol = get_volume(elCrd,elType)
                layVol = list()
                layMat = list()
                for lay in sec['layup']:
                    layVol.append(eVol*lay[1])
                    layMat.append(lay[0])
                stei = str(ei+1)
                elVols[stei] = layVol
                elMats[stei] = layMat
    
    output = dict()
    output['elVols'] = elVols
    output['elMats'] = elMats
    return output

def get_mesh_spatial_list(nodes,xSpacing=0,ySpacing=0,zSpacing=0):
    totNds = len(nodes)
    spaceDim = len(nodes[0])

    maxX = np.amax(nodes[:,0])
    minX = np.amin(nodes[:,0])
    maxY = np.amax(nodes[:,1])
    minY = np.amin(nodes[:,1])
    nto1_2 = np.power(totNds,0.5)
    nto1_3 = np.power(totNds,0.333333333333)
    if(spaceDim == 3):
        maxZ = np.amax(nodes[:,2])
        minZ = np.amin(nodes[:,2])
        dimVec = np.array([(maxX-minX),(maxY-minY),(maxZ-minZ)])
        meshDim = np.linalg.norm(dimVec)
        maxX = maxX + 0.01*meshDim
        minX = minX - 0.01*meshDim
        maxY = maxY + 0.01*meshDim
        minY = minY - 0.01*meshDim
        maxZ = maxZ + 0.01*meshDim
        minZ = minZ - 0.01*meshDim
        if(xSpacing == 0):
            xS = 0.5*(maxX - minX)/nto1_3
        else:
            xS = xSpacing
        if(ySpacing == 0):
            yS = 0.5*(maxY - minY)/nto1_3
        else:
            yS = ySpacing
        if(zSpacing == 0):
            zS = 0.5*(maxZ - minZ)/nto1_3
        else:
            zS = zSpacing
        meshGL = spatial_grid_list3d(minX,maxX,minY,maxY,minZ,maxZ,xS,yS,zS)
        #tol = 1.0e-6*meshDim/nto1_3
    else:
        dimVec = np.array([(maxX-minX),(maxY-minY)])
        meshDim = np.linalg.norm(dimVec)
        maxX = maxX + 0.01*meshDim
        minX = minX - 0.01*meshDim
        maxY = maxY + 0.01*meshDim
        minY = minY - 0.01*meshDim
        if(xSpacing == 0):
            xS = 0.5*(maxX - minX)/nto1_2
        else:
            xS = xSpacing
        if(ySpacing == 0):
            yS = 0.5*(maxY - minY)/nto1_2
        else:
            yS = ySpacing
        meshGL = spatial_grid_list2d(minX,maxX,minY,maxY,xS,yS)
        #tol = 1.0e-6*meshDim/nto1_2
    return meshGL

## - Convert list of mesh objects into a single merged mesh, returning sets representing the elements/nodes from the original meshes
def mergeDuplicateNodes(meshData,tolerance=None):
    allNds = meshData['nodes']
    allEls = meshData['elements']
    totNds = len(allNds)
    totEls = len(allEls)
    elDim = len(allEls[0])

    avgSp = get_average_node_spacing(meshData['nodes'],meshData['elements'])
    sp = 2*avgSp
    nodeGL = get_mesh_spatial_list(allNds,xSpacing=sp,ySpacing=sp,zSpacing=sp)
    if(tolerance == None):
        tol = 1.0e-4*avgSp
    else:
        tol = tolerance
    
    # glDim = nodeGL.getDim()
    # mag = np.linalg.norm(glDim)
    # nto1_3 = np.power(len(allNds),0.3333333333)
    # tol = 1.0e-6*mag/nto1_3

    i = 0
    for nd in allNds:
        nodeGL.addEntry(i, nd)
        i = i + 1

    ndElim = -np.ones(totNds, dtype=int)
    ndNewInd = -np.ones(totNds, dtype=int)
    for n1i in range(0, totNds):
        if ndElim[n1i] == -1:
            nearNds = nodeGL.findInRadius(allNds[n1i], tol)
            for n2i in nearNds:
                if n2i > n1i and ndElim[n2i] == -1:
                    proj = allNds[n2i] - allNds[n1i]
                    dist = np.linalg.norm(proj)
                    if dist < tol:
                        ndElim[n2i] = n1i
    ndi = 0
    nodesFinal = list()
    for n1i in range(0, totNds):
        if ndElim[n1i] == -1:
            nodesFinal.append(allNds[n1i])
            ndNewInd[n1i] = ndi
            ndi = ndi + 1
    nodesFinal = np.array(nodesFinal)
    for eli in range(0, totEls):
        for j in range(0, elDim):
            nd = allEls[eli, j]
            if nd != -1:
                if ndElim[nd] == -1:
                    allEls[eli, j] = ndNewInd[nd]
                else:
                    allEls[eli, j] = ndNewInd[ndElim[nd]]

    meshData["nodes"] = nodesFinal
    meshData["elements"] = allEls

    return meshData
    
def merge_meshes(mData1,mData2,tolerance=None):
    mergedData = dict()
    nds1 = mData1['nodes']
    nLen1 = len(nds1)
    nds2 = mData2['nodes']
    nLen2 = len(nds2)
    totNds = nLen1 + nLen2
    mrgNds = np.zeros((totNds,3),dtype=float)
    mrgNds[0:nLen1] = nds1
    mrgNds[nLen1:totNds] = nds2
    els1 = mData1['elements']
    eLen1 = len(els1)
    els2 = mData2['elements']
    eLen2 = len(els2)
    totEls = eLen1 + eLen2
    eCols = len(els1[0])
    mrgEls = -1*np.ones((totEls,eCols),dtype=int)
    mrgEls[0:eLen1] = els1
    for i, el in enumerate(els2,start=eLen1):
        addVec = np.zeros(eCols,dtype=int)
        for j, nd in enumerate(el):
            if(nd > -1):
                addVec[j] = nLen1
        mrgEls[i] = el + addVec
    mergedData['nodes'] = mrgNds
    mergedData['elements'] = mrgEls
    mergedData['sets'] = dict()
    mergedData['sets']['node'] = list()
    mergedData['sets']['element'] = list()
    try:
        mergedData['sets']['node'].extend(mData1['sets']['node'])
    except:
        pass
    try:
        mergedData['sets']['element'].extend(mData1['sets']['element'])
    except:
        pass
    try:
        for ns in mData2['sets']['node']:
            newSet = dict()
            newSet['name'] = ns['name']
            labs = list()
            for nd in ns['labels']:
                labs.append(nd+nLen1)
            newSet['labels'] = labs
            mergedData['sets']['node'].append(newSet)
    except:
        pass
    try:
        for es in mData2['sets']['element']:
            newSet = dict()
            newSet['name'] = es['name']
            labs = list()
            for el in es['labels']:
                labs.append(el+eLen1)
            newSet['labels'] = labs
            mergedData['sets']['element'].append(newSet)
    except:
        pass
    return mergeDuplicateNodes(mergedData,tolerance)

def add_node_set(meshData,newSet):
    try:
        meshData['sets']['node'].append(newSet)
    except:
        nSets = list()
        nSets.append(newSet)
        try:
            meshData['sets']['node'] = nSets
        except:
            sets = dict()
            sets['node'] = nSets
            meshData['sets'] = sets
    return meshData

def add_element_set(meshData,newSet):
    try:
        meshData['sets']['element'].append(newSet)
    except:
        elSets = list()
        elSets.append(newSet)
        try:
            meshData['sets']['element'] = elSets
        except:
            sets = dict()
            sets['element'] = elSets
            meshData['sets'] = sets
    return meshData

def get_matching_node_sets(meshData):
    elements = meshData['elements']
    elSets = meshData['sets']['element']
    nodeSets = list()
    for es in elSets:
        ns = set()
        for ei in es['labels']:
            for elnd in elements[ei]:
                if(elnd > -1):
                    ns.add(elnd)
        newSet = dict()
        newSet['name'] = es['name']
        newSet['labels'] = list(ns)
        nodeSets.append(newSet)
    try:
        meshData['sets']['node'].extend(nodeSets)
    except:
        meshData['sets']['node'] = nodeSets
        
    return meshData

def get_extruded_sets(meshData,numLayers):
    numEls = len(meshData['elements'])
    numNds = len(meshData['nodes'])
    extSets = dict()
    try:
        elSets = meshData['sets']['element']
        extES = list()
        for es in elSets:
            labels = list()
            for lay in range(0,numLayers):
                for ei in es['labels']:
                    newLab = ei + numEls*lay
                    labels.append(newLab)
            newSet = dict()
            newSet['name'] = es['name']
            newSet['labels'] = labels
            extES.append(newSet)
        extSets['element'] = extES
    except:
        pass
    
    try:
        ndSets = meshData['sets']['node']
        extNS = list()
        for ns in ndSets:
            labels = list()
            for lay in range(0,(numLayers + 1)):
                for ni in ns['labels']:
                    newLab = ni + numNds*lay
                    labels.append(newLab)
            newSet = dict()
            newSet['name'] = ns['name']
            newSet['labels'] = labels
            extNS.append(newSet)
        extSets['node'] = extNS
    except:
        pass        
        
    return extSets

def get_element_set_union(meshData,setList,newSetName):
    un = set()
    for es in meshData['sets']['element']:
        if(es['name'] in setList):
            thisSet = set(es['labels'])
            un = un.union(thisSet)
    newSet = dict()
    newSet['name'] = newSetName
    newSet['labels'] = list(un)
    meshData['sets']['element'].append(newSet)
    return meshData

def make3D(meshData):
    numNodes = len(meshData["nodes"])
    nodes3D = np.zeros((numNodes, 3))
    nodes3D[:, 0:2] = meshData["nodes"]
    dataOut = dict()
    dataOut["nodes"] = nodes3D
    dataOut["elements"] = meshData["elements"]
    return dataOut

def get_surface_nodes(meshData,elSet,newSetName,normDir,normTol=5.0):
    nds = meshData['nodes']
    els = meshData['elements']
    mag = np.linalg.norm(normDir)
    unitNorm = (1.0/mag)*normDir
    cosTol = np.cos(normTol*np.pi/180.0)
    faceDic = dict()
    for es in meshData['sets']['element']:
        if(es['name'] == elSet):
            for ei in es['labels']:
                fcStr, globFc = get_sorted_face_strings(els[ei])
                for fi, fk in enumerate(fcStr):
                    try:
                        curr = faceDic[fk]
                        faceDic[fk] = None
                    except:
                        faceDic[fk] = globFc[fi]
    surfSet = set()
    for fk in faceDic:
        glob = faceDic[fk]
        if(glob is not None):
            gLen = len(glob)
            if(gLen == 3):
                v1 = nds[glob[1]] - nds[glob[0]]
                v2 = nds[glob[2]] - nds[glob[1]]
            elif(gLen == 4):
                v1 = nds[glob[2]] - nds[glob[0]]
                v2 = nds[glob[3]] - nds[glob[1]]
            cp = cross_prod(v1,v2)
            mag = np.linalg.norm(cp)
            fcNrm = (1.0/mag)*cp
            dp = np.dot(fcNrm,unitNorm)
            if(dp >= cosTol):
                for nd in glob:
                    surfSet.add(nd)
    newSet = dict()
    newSet['name'] = newSetName
    newSet['labels'] = list(surfSet)
    return add_node_set(meshData,newSet)

def tie_2_meshes_constraints(tiedMesh,tiedSetName,tgtMesh,tgtSetName,maxDist):
    tiedNds = tiedMesh['nodes']
    tgtNds = tgtMesh['nodes']
    tgtEls = tgtMesh['elements']
    for ns in tiedMesh['sets']['node']:
        if(ns['name'] == tiedSetName):
            tiedSet = ns['labels']
    for es in tgtMesh['sets']['element']:
        if(es['name'] == tgtSetName):
            tgtSet = es['labels']
    radius = get_average_node_spacing(tgtNds,tgtEls)
    elGL = get_mesh_spatial_list(tgtNds,radius,radius,radius)
    if(radius < maxDist):
        radius = maxDist
    for ei in tgtSet:
        el = tgtEls[ei]
        fstNd = tgtNds[el[0]]
        elGL.addEntry(ei,fstNd)
    solidStr = 'tet4 wedge6 brick8'
    constraints = list()
    for ni in tiedSet:
        nd = tiedNds[ni]
        nearEls = elGL.findInRadius(nd,radius)
        minDist = 1.0e+100
        minPO = dict()
        minEi = -1
        for ei in nearEls:
            if(len(tgtEls[ei]) <= 4):
                if(tgtEls[ei,3] == -1):
                    elType = 'shell3'
                else:
                    elType = 'shell4'
            elif(len(tgtEls[ei]) <= 8):
                if(tgtEls[ei,4] == -1):
                    elType = 'tet4'
                elif(tgtEls[ei,6] == -1):
                    elType = 'wedge6'
                else:
                    elType = 'brick8'
            else:
                pstr = 'Warning: encountered unsupported element type in tie2MeshesConstraints'
                print(pstr)
            xC = []
            yC = []
            zC = []
            for en in tgtEls[ei]:
                if(en > -1):
                    xC.append(tgtNds[en,0])
                    yC.append(tgtNds[en,1])
                    zC.append(tgtNds[en,2])
            elCrd = np.array([xC,yC,zC])
            pO = get_proj_dist(elCrd,elType,nd)
            if(elType in solidStr):
                if(pO['distance'] > 0.0):
                    solidPO = get_solid_surf_proj(elCrd,elType,nd)
                    if(solidPO['distance'] < minDist):
                        minDist = solidPO['distance']
                        minPO = solidPO
                        minEi = ei
                else:
                    minDist = 0.0
                    minPO = pO
                    minEi = ei
            else:
                if(pO['distance'] < minDist):
                    minDist = pO['distance']
                    minPO = pO
                    minEi = ei
        if(minDist < maxDist):
            newConst = dict()
            terms = list()
            newTerm = dict()
            newTerm['nodeSet'] = 'tiedMesh'
            newTerm['node'] = ni
            newTerm['coef'] = -1.0
            terms.append(newTerm)
            nVec = minPO['nVec']
            nVi = 0
            for en in tgtEls[minEi]:
                if(en > -1):
                    newTerm = dict()
                    newTerm['nodeSet'] = 'targetMesh'
                    newTerm['node'] = en
                    newTerm['coef'] = nVec[nVi]
                    terms.append(newTerm)
                    nVi = nVi + 1
            newConst['terms'] = terms
            newConst['rhs'] = 0.0
            constraints.append(newConst)
    return constraints

def tie_2_sets_constraints(mesh,tiedSetName,tgtSetName,maxDist):
    try:
        elements = mesh['elements']
        nodes = mesh['nodes']
        elSets = mesh['sets']['element']
        ndSets = mesh['sets']['node']
        for es in elSets:
            if(es['name'] == tgtSetName):
                tgtSet = es['labels']
        for ns in ndSets:
            if(ns['name'] == tiedSetName):
                tiedSet = ns['labels']
            if(ns['name'] == tgtSetName):
                tgtNdSet = ns['labels']
        fstEl = tgtSet[0]
        fstNd = tgtSet[0]
        fstTgtNd = tgtNdSet[0]
    except:
        raise Exception('There was a problem accessing the mesh data in tie2SetsConstraints().  Check the set names and make sure nodes, elements and sets exist in the input mesh')

    tgtNdCrd = []
    for ni in tgtNdSet:
        tgtNdCrd.append(nodes[ni])
    tgtNdCrd = np.array(tgtNdCrd)
    
    radius = get_average_node_spacing(nodes, elements)
    elGL = get_mesh_spatial_list(tgtNdCrd,radius,radius,radius)
    if(radius < maxDist):
        radius = maxDist
    for ei in tgtSet:
        fstNd = tgtNdCrd[elements[ei,0]]
        elGL.addEntry(ei,fstNd)
        ei = ei + 1    
    
    solidStr = 'tet4 wedge6 brick8'
    constraints = list()
    for ni in tiedSet:
        nd = nodes[ni]
        nearEls = elGL.findInRadius(nd,radius)
        minDist = 1.0e+100
        minPO = dict()
        minEi = -1
        for ei in nearEls:
            if(len(elements[ei]) <= 4):
                if(elements[ei,3] == -1):
                    elType = 'shell3'
                else:
                    elType = 'shell4'
            elif(len(elements[ei]) <= 8):
                if(elements[ei,4] == -1):
                    elType = 'tet4'
                elif(elements[ei,6] == -1):
                    elType = 'wedge6'
                else:
                    elType = 'brick8'
            else:
                pstr = 'Warning: encountered unsupported element type in tie2SetsConstraints'
                print(pstr)
            xC = []
            yC = []
            zC = []
            for en in elements[ei]:
                if(en > -1):
                    xC.append(nodes[en,0])
                    yC.append(nodes[en,1])
                    zC.append(nodes[en,2])
            elCrd = np.array([xC,yC,zC])
            pO = get_proj_dist(elCrd,elType,nd)
            if(elType in solidStr):
                if(pO['distance'] > 0.0):
                    solidPO = get_solid_surf_proj(elCrd,elType,nd)
                    if(solidPO['distance'] < minDist):
                        minDist = solidPO['distance']
                        minPO = solidPO
                        minEi = ei
                else:
                    minDist = 0.0
                    minPO = pO
                    minEi = ei
            else:
                if(pO['distance'] < minDist):
                    minDist = pO['distance']
                    minPO = pO
                    minEi = ei
        if(minDist < maxDist and (ni not in elements[minEi])):
            newConst = dict()
            terms = list()
            newTerm = dict()
            newTerm['nodeSet'] = tiedSetName
            newTerm['node'] = ni
            newTerm['coef'] = -1.0
            terms.append(newTerm)
            nVec = minPO['nVec']
            nVi = 0
            for en in elements[minEi]:
                if(en > -1):
                    newTerm = dict()
                    newTerm['nodeSet'] = tgtSetName
                    newTerm['node'] = en
                    newTerm['coef'] = nVec[nVi]
                    terms.append(newTerm)
                    nVi = nVi + 1
            newConst['terms'] = terms
            newConst['rhs'] = 0.0
            constraints.append(newConst)
    
    return constraints

def plotNodes(meshData):
    xLst = meshData['nodes'][:,0]
    yLst = meshData['nodes'][:,1]
    try:
        zLst = meshData['nodes'][:,2]
    except:
        zLst = np.zeros(len(xLst),dtype=float)
    
    fig = go.Figure(data=[go.Scatter3d(x=xLst,y=yLst,z=zLst,mode='markers')])
    
    fig.show()

def plotShellMesh(meshData):
    xLst = meshData["nodes"][:, 0]
    yLst = meshData["nodes"][:, 1]
    try:
        zLst = meshData["nodes"][:, 2]
    except:
        zLst = np.zeros(len(xLst))
    value = list()
    v1 = list()
    v2 = list()
    v3 = list()
    i = 0
    for el in meshData["elements"]:
        v1.append(el[0])
        v2.append(el[1])
        v3.append(el[2])
        value.append(np.sin(i))
        if el[3] != -1:
            v1.append(el[0])
            v2.append(el[2])
            v3.append(el[3])
            value.append(np.sin(i))
        i = i + 1
    fig = go.Figure(
        data=[
            go.Mesh3d(
                x=xLst,
                y=yLst,
                z=zLst,
                colorbar_title="",
                colorscale=[[0.0, "white"], [0.5, "gray"], [1.0, "black"]],
                intensity=value,
                intensitymode="cell",
                i=v1,
                j=v2,
                k=v3,
                name="",
                showscale=True,
            )
        ]
    )

    fig.show()

def plotSolidMesh(meshData):
    xLst = meshData["nodes"][:, 0]
    yLst = meshData["nodes"][:, 1]
    zLst = meshData["nodes"][:, 2]
    value = list()
    v1 = list()
    v2 = list()
    v3 = list()
    i = 0
    for el in meshData["elements"]:
        si = np.sin(i)
        if el[4] == -1:
            v1.append(el[0])
            v2.append(el[1])
            v3.append(el[2])
            value.append(si)
            v1.append(el[0])
            v2.append(el[1])
            v3.append(el[3])
            value.append(si)
            v1.append(el[0])
            v2.append(el[2])
            v3.append(el[3])
            value.append(si)
            v1.append(el[1])
            v2.append(el[2])
            v3.append(el[3])
            value.append(si)
        elif el[6] == -1:
            v1.append(el[0])
            v2.append(el[1])
            v3.append(el[2])
            value.append(si)
            v1.append(el[3])
            v2.append(el[4])
            v3.append(el[5])
            value.append(si)
            v1.append(el[0])
            v2.append(el[1])
            v3.append(el[3])
            value.append(si)
            v1.append(el[1])
            v2.append(el[3])
            v3.append(el[4])
            value.append(si)

            v1.append(el[0])
            v2.append(el[2])
            v3.append(el[3])
            value.append(si)
            v1.append(el[2])
            v2.append(el[3])
            v3.append(el[5])
            value.append(si)
            v1.append(el[1])
            v2.append(el[2])
            v3.append(el[4])
            value.append(si)
            v1.append(el[2])
            v2.append(el[4])
            v3.append(el[5])
            value.append(si)
        else:
            v1.append(el[0])
            v2.append(el[3])
            v3.append(el[4])
            value.append(si)
            v1.append(el[3])
            v2.append(el[4])
            v3.append(el[7])
            value.append(si)
            v1.append(el[1])
            v2.append(el[2])
            v3.append(el[5])
            value.append(si)
            v1.append(el[2])
            v2.append(el[5])
            v3.append(el[6])
            value.append(si)

            v1.append(el[0])
            v2.append(el[1])
            v3.append(el[4])
            value.append(si)
            v1.append(el[1])
            v2.append(el[4])
            v3.append(el[5])
            value.append(si)
            v1.append(el[2])
            v2.append(el[3])
            v3.append(el[6])
            value.append(si)
            v1.append(el[3])
            v2.append(el[6])
            v3.append(el[7])
            value.append(si)

            v1.append(el[0])
            v2.append(el[1])
            v3.append(el[2])
            value.append(si)
            v1.append(el[0])
            v2.append(el[2])
            v3.append(el[3])
            value.append(si)
            v1.append(el[4])
            v2.append(el[5])
            v3.append(el[6])
            value.append(si)
            v1.append(el[4])
            v2.append(el[6])
            v3.append(el[7])
            value.append(si)
        i = i + 1
    fig = go.Figure(
        data=[
            go.Mesh3d(
                x=xLst,
                y=yLst,
                z=zLst,
                colorbar_title="",
                colorscale=[[0.0, "white"], [0.5, "gray"], [1.0, "black"]],
                intensity=value,
                intensitymode="cell",
                i=v1,
                j=v2,
                k=v3,
                name="",
                showscale=True,
            )
        ]
    )

    fig.show()


## -Create node/element set within a spatial range or radius
