import pynumad as pynu
import numpy as np
import os
from os.path import join

from pynumad.mesh_gen.mesh_gen import get_shell_mesh
from pynumad.io.mesh_to_yaml import *
from pynumad.analysis.abaqus.write import *

## Define inputs
bladeYaml = join("example_data","blade.yaml")
meshYaml = "BAR0.yaml"
abqFileName = "BAR0.inp"
abqScriptName = "runModalAnalysis.py"
dampFileName = "BAR0DampingInp.yaml"
adhesiveMat = "Adhesive"

## Read blade data from yaml file
blade = pynu.Blade()
blade.read_yaml(bladeYaml)

## Set the airfoil point resolution
for stat in blade.definition.stations:
    stat.airfoil.resample(n_samples=300)
    
#blade.generate_geometry()
blade.update_blade()
nStations = blade.geometry.coordinates.shape[2]
minTELengths = 0.001*np.ones(nStations)
blade.expand_blade_geometry_te(minTELengths)

## Set the target element size for the mesh
elementSize = 0.2

## Generate mesh
adhes = 1
bladeMesh = get_shell_mesh(blade, adhes, elementSize)

## Add a section to define the adhesive line
sec = dict()
sec['type'] = 'solid'
sec['elementSet'] = 'adhesiveElSet'
sec['material'] = adhesiveMat
bladeMesh['adhesiveSection'] = sec

mesh_to_yaml(bladeMesh,meshYaml)

## Write Abaqus input

numModes = 3
write_shell_modal_input(abqFileName, blade, bladeMesh, adhesiveMat, numModes)
