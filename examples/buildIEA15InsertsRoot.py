import pynumad as pynu
import numpy as np
import os
from os.path import join

from pynumad.mesh_gen.mesh_gen import get_root_mesh
from pynumad.io.mesh_to_yaml import *
from pynumad.analysis.abaqus.write import *

##  This example builds a solid element model of the root of the IEA-15MW reference offshore wind turbine
##  using a conventional bolt insert design.

## Define input yaml file and name of abaqus input file to be generated
bladeYaml = 'example_data/IEA-15-240-RWT.yaml'
abaqusFile = 'IEA_15_insert_root.inp'

## Read the model into a blade object
blade = pynu.Blade()
blade.read_yaml(bladeYaml)

## Define blade root specifications
rootDiameter = blade.definition.chord[0]
rtRad = 0.5*rootDiameter

axisPts = [0.,1.0]  ## Range of points along the spanwise axis in meters
radius = [rtRad,rtRad]  ## Radius corresponding to spanwise axis points
boltRad = 0.024  ## Radius of a single bolt
insThk = 0.01    ## Thickness of bolt insert layer
adhThk = 0.01    ## Thickness of bolt adhesive layer
rtThk = 3.0*(boltRad + adhThk + insThk)  ## Total thickness of wall at root
thickness = [rtThk,rtThk]  ## Axis point thickness
nIns = int(np.pi*rootDiameter/rtThk)  ## Number of inserts
elementSize = adhThk  ## Nominal element size for the model
elLayers = 10  ##  Number of element layers in the spanwise direction

## Generate the root mesh
insertMesh = get_root_mesh(axisPts,radius,thickness,elementSize,elLayers,config='inserts',boltRad=boltRad,insertThk=insThk,adhesiveThk=adhThk,numIns=nIns)

## Define material properties specific to root design
materials = list()

mat = dict()
mat['name'] = 'boltMat'
mat['density'] = 7800.0
mat['elastic'] = {'E': [200.0e+9,200.0e+9,200.0e+9], 
                  'nu': [0.3,0.3,0.3], 
                  'G': [79.3e+9,79.3e+9,79.3e+9]}
materials.append(mat)

mat = dict()
mat['name'] = 'insertMat'
mat['density'] = 7800.0
mat['elastic'] = {'E': [200.0e+9,200.0e+9,200.0e+9], 
                  'nu': [0.3,0.3,0.3], 
                  'G': [79.3e+9,79.3e+9,79.3e+9]}
materials.append(mat)

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

insertMesh['materials'] = materials

##  Write the Abaqus input file

write_solid_general(abaqusFile,insertMesh)