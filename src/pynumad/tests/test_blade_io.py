import unittest
import os
import pickle

from pynumad.objects.blade import Blade

test_data_dir = os.path.join(os.path.dirname(__file__), "test_data")


class TestBladeIO(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.xlsxfile = os.path.join(test_data_dir, "blades", "blade.xlsx")
        cls.yamlfile = os.path.join(test_data_dir, "blades", "blade.yaml")
        
        with open(os.path.join(test_data_dir,"blades","yaml_blade.pkl"), "rb") as f:
            cls.yaml_blade_pkl = pickle.load(f)
            
        with open(os.path.join(test_data_dir,"blades","excel_blade.pkl"), "rb") as f:
            cls.excel_blade_pkl = pickle.load(f)

    def test_xlsx_blade(self):
        xlsxblade = Blade(self.xlsxfile)
        # assert xlsxblade.definition == self.excel_blade_pkl.definition
        

    def test_yaml_blade(self):
        yamlblade = Blade(self.yamlfile)
        # assert yamlblade.definition == self.yaml_blade_pkl.definition


if __name__ == "__main__":
    unittest.main()
