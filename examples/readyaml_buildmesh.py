import pynumad as pynu
import numpy as np
import pickle
from pprint import pprint

from pynumad.shell.shell import get_shell_mesh

blade = pynu.Blade()
fileName = 'example_data/myBlade.yaml'
blade.read_yaml(fileName)

elementSize = 0.2
adhes = 1

meshData = get_shell_mesh(blade, adhes, elementSize)

# Print first 10 nodes coordinates
pprint(meshData['nodes'][:10,:])

# Print first 10 element connectivities
pprint(meshData['elements'][:10,:])
