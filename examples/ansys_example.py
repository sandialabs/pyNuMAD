import pynumad
import numpy as np
import matplotlib.pyplot as plt

import json
from pynumad.utils.misc_utils import setup_logging
from pynumad.shell.shell import shell_mesh_general
from pynumad.analysis.ansys.write import write_ansys_shell_model
from pynumad.analysis.ansys.run import call_ansys
from pynumad.utils.distributed_loading import loads_table_coordinate_trans
from pynumad.analysis.ansys.main_ansys_analysis import main_ansys_analysis






blade=pynumad.Blade()
yamlName='myBlade_Modified'
blade.read_yaml('example_data/'+yamlName+'.yaml') 

log=setup_logging(yamlName+'_ansys')


# create a shell model in ansys w/o adhesive
elementSize=0.45
mesh_data=shell_mesh_general(blade, forSolid=False, includeAdhesive=False, elementSize=elementSize)

config = dict()
config["elementType"] = '181'
config['blade_name'] = yamlName
filename=write_ansys_shell_model(blade, mesh_data, config)

call_ansys(filename,log)


# Load a previously built loads_table
with open('myBlade_loadsTable.json','r') as fid:
    loads_table = [json.load(fid)]

loads_table=loads_table_coordinate_trans(loads_table)

# Set up configuration for deflection run

analysis_config = dict()
analysis_config['meshFile'] = 'master.db'
analysis_config['analysisFileName'] = 'bladeAnalysis'
analysis_config['np'] = 2
analysis_config['analysisFlags'] = dict()
analysis_config['analysisFlags']['mass'] = True
analysis_config['analysisFlags']['deflection'] = True
analysis_config['analysisFlags']['failure'] = 'smax'
analysis_config['analysisFlags']['fatigue'] = 'all'
analysis_config['analysisFlags']['resultants'] = True
analysis_config['analysisFlags']['local_fields'] = 'all'
analysis_config['analysisFlags']['globalBuckling'] = 3

ansys_result = main_ansys_analysis(blade,mesh_data,loads_table,analysis_config,elementSize,log)

# Plot the flapwise deflection
y=ansys_result['deflection'][0][1]

plt.figure()
plt.plot(blade.ispan,y)


print('')