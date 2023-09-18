from os.path import abspath, dirname, join
import json
# DATA_PATH is intended to be used by pynumad internals to find airfoil files
#   as needed by the excel_to_blade function
package_dir = dirname(abspath(str(__file__)))
DATA_PATH = join(package_dir, "..", "data")

with open(join(package_dir, "software_paths.json"), "rt") as f:
    SOFTWARE_PATHS = json.load(f)