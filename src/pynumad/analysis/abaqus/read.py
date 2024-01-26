# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 15:00:28 2023

@author: evaande
"""
import csv
import yaml
from yaml import CLoader as Loader
from pynumad.mesh_gen.mesh_tools import get_element_volumes

## read_csv_model_output will likely be phased out, but currently remains for reference
def read_csv_modal_output(fileName):
    outData = dict()
    outData['modes'] = dict()
    inFile = open(fileName,'r')
    reader = csv.reader(inFile)
    colHd = ['Frame', 'Instance Name', 'Element Label', 'IntPt', 'Section Point',
             'E-E11', 'E-E22', 'E-E33', 'E-E12', 'E-E13', 'E-E23',
             'S-S11', 'S-S22', 'S-S33', 'S-S12', 'S-S13', 'S-S23']
    colInd = dict()
    for hd in colHd:
        colInd[hd] = -1
    Est = 5
    Sst = 11
    for row in reader:
        if(colInd['Element Label'] == -1):
            for j, hd in enumerate(row):
                sthd = hd.strip()
                if(sthd in colInd):
                    colInd[sthd] = j
        else:
            newEnt = dict()
            newEnt['label'] = int(row[colInd['Element Label']])
            j = colInd['Instance Name']
            if(j > -1):
                newEnt['part'] = row[j]
            j = colInd['IntPt']
            if(j > -1):
                newEnt['intPt'] = int(row[j])
            spStr = row[colInd['Section Point']]
            if('Layer =' in spStr):
                sLst = spStr.split('Layer =')
                newEnt['layer'] = int(sLst[1])
            stress = [0.,0.,0.,0.,0.,0.]
            strain = [0.,0.,0.,0.,0.,0.]
            for i in range(0,6):
                j = colInd[colHd[Est+i]]
                if(j > -1):
                    strain[i] = float(row[j])
                j = colInd[colHd[Sst+i]]
                if(j > -1):
                    strain[i] = float(row[j])
            newEnt['stress'] = stress
            newEnt['strain'] = strain
            mdStr = 'mode_unkn'
            j = colInd['Frame']
            if(j > -1):
                frmStr = row[j]
                if('Mode' in frmStr):
                    lst1 = frmStr.split('Mode')
                    lst2 = lst1[1].split(':')
                    mdStr = 'mode_' + lst2[0].strip()
            try:
                outData['modes'][mdStr].append(newEnt)
            except:
                outData['modes'][mdStr] = [newEnt]
    inFile.close()
    return outData
    
def build_damping_input(abq_results,material_loss_factors,meshData):
    inFile = open(abq_results,'r')
    dampInp = yaml.load(inFile, Loader=Loader)
    inFile.close()
    modes = dampInp['modes']
    bladeEV = get_element_volumes(meshData)
    try:
        adMeshData = dict()
        adMeshData['nodes'] = meshData['adhesiveNds']
        adMeshData['elements'] = meshData['adhesiveEls']
        es = dict()
        es['name'] = 'adhesiveEls'
        numAdEl = len(meshData['adhesiveEls'])
        es['labels'] = list(range(0,numAdEl))
        elSets = [es]
        sets = dict()
        sets['element'] = elSets
        adMeshData['sets'] = sets
        adMeshData['sections'] = [meshData['adhesiveSection']]
        adEV = get_element_volumes(adMeshData)
    except:
        pass
    for md in modes:
        thisMd = modes[md]
        volume = list()
        material = list()
        for i, lab in enumerate(thisMd['label']):
            labSt = str(lab)
            part = thisMd['part'][i]
            if('adhesive' in part.lower()):
                volume.append(adEV['elVols'][labSt])
                material.append(adEV['elMats'][labSt])
            else:
                lay = thisMd['layer'][i] - 1
                if(lay >= 0):
                    volume.append(bladeEV['elVols'][labSt][lay])
                    material.append(bladeEV['elMats'][labSt][lay])
                else:
                    volume.append(bladeEV['elVols'][labSt])
                    material.append(bladeEV['elMats'][labSt])
        modes[md]['volume'] = volume
        modes[md]['material'] = material
    dampInp['mat_loss_factors'] = material_loss_factors
    return dampInp