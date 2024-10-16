import yaml
from yaml import CLoader as Loader
import numpy as np
import os


def mesh_to_yaml(meshData, file_name):
    """
    TODO docstring
    """
    mDataOut = dict()
    nodes = list()
    for nd in meshData["nodes"]:
        ndstr = str(list(nd))
        nodes.append(ndstr)
    elements = list()
    for el in meshData["elements"]:
        elstr = str(list(el))
        elements.append(elstr)
    esList = list()
    for es in meshData["sets"]["element"]:
        newSet = dict()
        newSet["name"] = es["name"]
        labels = list()
        for el in es["labels"]:
            labels.append(int(el))
        newSet["labels"] = labels
        esList.append(newSet)
    nsList = list()
    try:
        for ns in meshData["sets"]["node"]:
            newSet = dict()
            newSet["name"] = ns["name"]
            labels = list()
            for nd in ns["labels"]:
                labels.append(int(nd))
            newSet["labels"] = labels
            nsList.append(newSet)
    except:
        pass
    sections = list()
    for sec in meshData["sections"]:
        newSec = dict()
        newSec["type"] = sec["type"]
        newSec["elementSet"] = sec["elementSet"]
        if sec["type"] == "shell":
            newLayup = list()
            for lay in sec["layup"]:
                laystr = str(lay)
                newLayup.append(laystr)
            newSec["layup"] = newLayup
            newSec["xDir"] = str(list(sec["xDir"]))
            newSec["xyDir"] = str(list(sec["xyDir"]))
        else:
            newSec["material"] = sec["material"]
        sections.append(newSec)
    elOri = list()
    for ori in meshData["elementOrientations"]:
        elOri.append(str(list(ori)))
    
    mDataOut["nodes"] = nodes
    mDataOut["elements"] = elements
    mDataOut["sets"] = dict()
    mDataOut["sets"]["element"] = esList
    mDataOut["sections"] = sections
    mDataOut["elementOrientations"] = elOri
    try:
        adNds = list()
        for nd in meshData["adhesiveNds"]:
            ndstr = str(list(nd))
            adNds.append(ndstr)
        adEls = list()
        for el in meshData["adhesiveEls"]:
            elstr = str(list(el))
            adEls.append(elstr)
        mDataOut["adhesiveNds"] = adNds
        mDataOut["adhesiveEls"] = adEls
        mDataOut["adhesiveElSet"] = meshData["adhesiveElSet"]
        mDataOut["adhesiveSection"] = meshData["adhesiveSection"]
        constraints = list()
        for c in meshData["constraints"]:
            terms = list()
            for t in c["terms"]:
                newTerm = dict()
                newTerm["nodeSet"] = t["nodeSet"]
                newTerm["node"] = int(t["node"])
                newTerm["coef"] = float(t["coef"])
                terms.append(newTerm)
            newConst = dict()
            newConst["terms"] = terms
            newConst["rhs"] = c["rhs"]
            constraints.append(newConst)
        mDataOut["constraints"] = constraints
    except:
        pass
    
    try:
        mDataOut["materials"] = meshData["materials"]
    except:
        pass
    
    fileStr = yaml.dump(mDataOut,sort_keys=False)
    
    fileStr = fileStr.replace("'","")
    fileStr = fileStr.replace('"','')
    
    outStream = open(file_name,'w')
    outStream.write(fileStr)
    outStream.close()
    
def yaml_to_mesh(fileName):
    inFile = open(fileName,'r')
    meshData = yaml.load(inFile,Loader=Loader)
    inFile.close()
    meshData['nodes'] = np.array(meshData['nodes'])
    meshData['elements'] = np.array(meshData['elements'])
    try:
        meshData['adhesiveNds'] = np.array(meshData['adhesiveNds'])
        meshData['adhesiveEls'] = np.array(meshData['adhesiveEls'])
    except:
        pass
    return meshData