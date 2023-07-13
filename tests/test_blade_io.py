import unittest
from os.path import abspath, dirname, join

from pynumad.objects.Blade import Blade
from pynumad.paths import DATA_PATH

test_data_dir = DATA_PATH

class TestBladeIO(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.xlsxfile = join(test_data_dir,"blade_excels","Excel2ObjectExample.xlsx")
        self.yamlfile = join(test_data_dir,"blade_yamls","BAR1_SNL_1_18_2021.yaml")
        
    def test_xlsx_blade(self):
        xlsxblade = Blade(self.xlsxfile)

    def test_yaml_blade(self):
        yamlblade = Blade(self.yamlfile)

if __name__ == "__main__":
    unittest.main()