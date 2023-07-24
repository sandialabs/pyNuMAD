from os.path import abspath, dirname, join

ANSYS_PATH = ""
CUBIT_PATH = ""

package_dir = dirname(abspath(str(__file__)))
DATA_PATH = join(package_dir, "..", "data")
