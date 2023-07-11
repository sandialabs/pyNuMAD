import unittest
from os.path import abspath, dirname, join

from pynumad.objects.Airfoil import Airfoil
from pynumad.paths import DATA_PATH

testdir = dirname(abspath(str(__file__)))
test_data_dir = DATA_PATH

class TestAirfoil(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.xmlfile = join(test_data_dir,"airfoils", "DU91-W-250.txt")

    def test_load_xml(self):
        x = Airfoil(filename = self.xmlfile)

        #check reference
        self.assertEqual(x.reference, "Points generated by BRR for NuMAD, 6/2/2011")

        #TODO: check coords

if __name__ == "__main__":
    unittest.main()