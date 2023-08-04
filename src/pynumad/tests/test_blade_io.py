import unittest
import os

from pynumad.objects.blade import Blade

test_data_dir = os.path.join(os.path.dirname(__file__), 'test_data')


class TestBladeIO(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.xlsxfile = os.path.join(test_data_dir, "blades", "blade.xlsx")
        cls.yamlfile = os.path.join(test_data_dir, "blades", "blade.yaml")

    def test_xlsx_blade(self):
        xlsxblade = Blade(self.xlsxfile)

    def test_yaml_blade(self):
        yamlblade = Blade(self.yamlfile)


if __name__ == "__main__":
    unittest.main()
