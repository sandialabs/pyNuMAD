import pynumad as pynu
import numpy as np
import pickle
from pprint import pprint
from os.path import join

from pynumad.shell.shell import get_shell_mesh

blade = pynu.Blade()
fileName = join("example_data","blade.yaml")
blade.read_yaml(fileName)

elementSize = 0.2
adhes = 1

meshData = get_shell_mesh(blade, adhes, elementSize)

# Print first 10 nodes coordinates
pprint(meshData['nodes'][:10,:])

# Print first 10 element connectivities
pprint(meshData['elements'][:10,:])
