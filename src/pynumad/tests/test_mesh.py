import unittest
from os.path import abspath, dirname, join
from pynumad.analysis.ansys.write import writeAnsysShellModel
from pynumad.shell.shell import getShellMesh
from pynumad.objects.Blade import Blade
from pynumad.paths import DATA_PATH

test_data_dir = DATA_PATH


class TestMesh(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.yamlfile = join(test_data_dir, "blade_yamls", "myBlade_modified.yaml")

    def test_mesh(self):
        blade = Blade(self.yamlfile)
        elementSize = 0.2
        adhes = 1
        meshData = getShellMesh(blade, includeAdhesive=adhes, elementSize=elementSize)
