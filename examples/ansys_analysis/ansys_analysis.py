import pynumad
import numpy as np
import logging
import matplotlib.pyplot as plt
import pickle

import json
from pynumad.analysis.ansys.main_ansys_analysis import main_ansys_analysis
from pynumad.analysis.ansys.write import write_ansys_shell_model
from pynumad.shell.shell import shellMeshGeneral
from pynumad.analysis.ansys.run import call_ansys
from pynumad.utils.misc_utils import setup_logging
from pynumad.utils.distributed_loading import *
import matplotlib.pyplot as plt


blade=pynumad.Blade()
#yamlName='BAR4_SNL_1_18_2021'
#yamlName='BAR0_NREL_1_4_2021'
#yamlName='barURC'
#yamlName='myBlade'
#yamlName='IEA-15-240-RWT'
#yamlName='IEA-22-280-RWT'

yamlName='myBlade_Modified'

log=setup_logging(yamlName+'_ansys')




if False:
    blade.read_yaml(f'{yamlName}.yaml')
    file = open(yamlName+'.pkl', 'wb')
    pickle.dump(blade, file)
    file.close()

file = open(yamlName+'.pkl', 'rb')
blade = pickle.load(file)
file.close()

# create a shell model in ansys w/o adhesive
mesh_data=shellMeshGeneral(blade, forSolid=False, includeAdhesive=False, elementSize=0.45)

config = dict()
config["elementType"] = '181'
config["MultipleLayerBehavior"] = 'multiply'
config["dbgen"] = 1
config['blade_name'] = yamlName
filename=write_ansys_shell_model(blade, mesh_data, config)

call_ansys(filename,log)


# Load a previously built loads_table
with open('/home/ecamare/repos/fork_pynumad/examples/ansys_analysis/myBlade_loadsTable.json','r') as fid:
    loads_table1 = json.load(fid)
with open('/home/ecamare/repos/fork_pynumad/examples/ansys_analysis/myBlade_loadsTable.json','r') as fid:
    loads_table2 = json.load(fid)

loads_table2['Fxb']=list(np.array(loads_table2['Fxb'])*3)
loads_table=[loads_table1]
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

ansys_result = main_ansys_analysis(blade,mesh_data,loads_table,analysis_config,log)

# Plot the flapwise deflection
y=ansys_result['deflection'][0][1]

plt.figure()
plt.plot(blade.ispan,y)

print('')