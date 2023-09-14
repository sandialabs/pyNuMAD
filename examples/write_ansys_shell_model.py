import numpy as np
from os.path import join
from pynumad.analysis.ansys.write import writeAnsysShellModel
from pynumad.shell.shell import get_shell_mesh
from pynumad.objects.blade import Blade

"""
This example loads in a blade, gets the shell mesh, and writes the shell
model to an ansys .src file
"""

yamlpath = join("..","example_data","blade.yaml")
blade = Blade(yamlpath)
elementSize = 0.2
adhes = 1

meshData = get_shell_mesh(blade, includeAdhesive=adhes, elementSize=elementSize)
config = dict()
config["BoundaryCondition"] = 'cantilevered'
config["elementType"] = '181'
config["MultipleLayerBehavior"] = 'multiply'
config["dbgen"] = 1
config["dbname"] = 'master'

filename = 'buildAnsysShell.src'
includeAdhesive = 1

writeAnsysShellModel(blade,filename,meshData,config,includeAdhesive)