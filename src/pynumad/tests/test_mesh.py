import unittest
import os
from pynumad.mesh_gen.mesh_gen import get_shell_mesh
from pynumad.objects.blade import Blade

test_data_dir = os.path.join(os.path.dirname(__file__), "test_data")


class TestMesh(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.yamlfile = os.path.join(test_data_dir, "blades", "blade.yaml")

    def test_mesh(self):
        blade = Blade(self.yamlfile)
        elementSize = 0.2
        adhes = 1
        meshData = get_shell_mesh(blade, includeAdhesive=adhes, elementSize=elementSize)
