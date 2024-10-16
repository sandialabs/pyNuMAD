import pynumad as pynu
import numpy as np
import os
from os.path import join

from pynumad.mesh_gen.mesh_gen import get_root_mesh
from pynumad.io.mesh_to_yaml import *
from pynumad.analysis.abaqus.write import *

##  This example builds a solid element model of the root of the IEA-15MW reference offshore wind turbine
##  using an imbedded solid tube insert.

## Define input yaml file and name of Abaqus input file to be generated
bladeYaml = 'example_data/IEA-15-240-RWT.yaml'
abaqusFile = 'IEA_15_tube_root.inp'

## Read the model into a blade object
blade = pynu.Blade()
blade.read_yaml(bladeYaml)

## Define blade root specifications
rootDiameter = blade.definition.chord[0]
rtRad = 0.5*rootDiameter

tubeThk = 0.048  ## Thickness of tube insert
adhThk = 0.01  ## Thickness of adhesive
axisPts = [0.,1.]  ## Range of points along the spanwise axis in meters
radius = [rtRad,rtRad]  ##  Radius corresponding to spanwise axis points
rtThk = 1.5*(tubeThk + 2.0*adhThk) ## Thickness of root
thickness = [rtThk,rtThk]  ## Axis point thickness
elementSize = adhThk  ## Nominal element size for the model
elLayers = 10  ##  Number of element layers in the spanwise direction
tubeExt = 1.0  ## Length of tube extension
extNE = 10  ## Number of element layer in tube extension

## Generate the root mesh
tubeMesh = get_root_mesh(axisPts,radius,thickness,elementSize,elLayers,config='tube',adhesiveThk=adhThk,tubeThk=tubeThk,tubeExtend=tubeExt,extNumEls=extNE)

## Define material properties specific to root design
materials = list()

mat = dict()
mat['name'] = 'adhesiveMat'
mat['density'] = 1100.0
mat['elastic'] = {'E': [4.56e+6,4.56e+6,4.56e+6], 
                  'nu': [0.49,0.49,0.49], 
                  'G': [1.52e+6,1.52e+6,1.52e+6]}
materials.append(mat)

mat = dict()
mat['name'] = 'fillMat'
mat['density'] = 1940.0
mat['elastic'] = {'E': [28.7e+9,16.6e+9,16.7e+9], 
                  'nu': [0.5,0.0,0.17], 
                  'G': [8.4e+9,3.49e+9,3.49e+9]}
materials.append(mat)

mat = dict()
mat['name'] = 'tubeMat'
mat['density'] = 7800.0
mat['elastic'] = {'E': [200.0e+9,200.0e+9,200.0e+9], 
                  'nu': [0.3,0.3,0.3], 
                  'G': [79.3e+9,79.3e+9,79.3e+9]}
materials.append(mat)

tubeMesh['materials'] = materials

##  Write the Abaqus input file

write_solid_general(abaqusFile,tubeMesh)