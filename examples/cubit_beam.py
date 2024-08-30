import sys 
import pynumad

sys.path.append(pynumad.SOFTWARE_PATHS['cubit'])
sys.path.append(pynumad.SOFTWARE_PATHS['cubit_enhancements'])

import cubit
from pynumad.analysis.cubit.make_blade import *
import numpy as np

from pynumad.analysis.make_models import write_beam_model
import logging


def get_cs_params():
    cs_params = {}
    cs_params['nel_per_layer'] = 3 
    cs_params['element_ar'] = 5
    cs_params['element_shape'] = 'quad'
    cs_params['element_thickness_ar'] =5

    cs_params['layer_transition_angle'] = 30
    cs_params['birds_mouth_amplitude_fraction']=1.0
    cs_params['minimum_layer_thickness'] = 0.001



    cs_params['adhesive_mat_name'] = 'Adhesive'


    totalStations = np.asarray(blade.ispan).size


    root_adhesive_thickness=0.001
    tip_adhesive_thickness=0.001
    te_adhesive_width_percent_chord=5
    adhesive_array=np.linspace(root_adhesive_thickness, tip_adhesive_thickness, num=totalStations)

    cs_params['web_fore_adhesive_thickness'] = adhesive_array
    cs_params['web_aft_adhesive_thickness'] = adhesive_array  
    cs_params['le_adhesive_thickness'] = adhesive_array
    cs_params['te_adhesive_thickness'] = adhesive_array
    cs_params['web_adhesive_width']=np.zeros((totalStations,))
    cs_params['te_adhesive_width']=np.zeros((totalStations,))
    cs_params['max_web_imperfection_distance']=np.zeros((totalStations,))

    #Spanwise varying parameters
    thickness_scaling=0.001
    geometry_scaling=thickness_scaling*1000
    for i_station in range(totalStations):
    
        if blade.definition.leband[i_station]!=0.0 and cs_params['le_adhesive_thickness'][i_station] > blade.definition.leband[i_station]*0.85/1000: #Adhesive thickness can't be greater than 85% of the LE band
            raise ValueError(f'LE adhesive thickness of {cs_params["le_adhesive_thickness"][i_station]} for station {i_station} is too large for the specified LE band width of {blade.definition.leband[i_station]/1000}.')
    

        #Adhesive width parameters
        cs_params['web_adhesive_width'][i_station] = 0.03 * blade.geometry.ichord[i_station]
        cs_params['te_adhesive_width'][i_station] = te_adhesive_width_percent_chord/100.0 * blade.geometry.ichord[i_station] * geometry_scaling
        
        cs_params['max_web_imperfection_distance'][i_station] = 0.0001 * blade.geometry.ichord[i_station] * geometry_scaling
        
    return cs_params 

blade=pynumad.Blade()

yamlName='myBlade_Modified'
blade.read_yaml('example_data/'+yamlName+'.yaml') 




wt_name=yamlName
dir_name='.'


cs_params=get_cs_params()
settings={}
settings['make_input_for']='VABSbeamdyn'  #SM, VABS, ANBA, or None
settings['export']='cubg' #cub, g, or None

station_list = []
cubit_make_cross_sections(blade,wt_name,settings,cs_params,'2D',stationList=station_list,directory=dir_name) #Note an empty list for stationList will make all cross sections.
    
#Proportional damping values for BeamDyn file.
mu=[0.00257593, 0.0017469,  0.0017469,  0.0017469,  0.00257593, 0.0017469]


#Set up logging
log =logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
fh=logging.FileHandler(wt_name+'driver.log',mode='w')
log.addHandler(fh)

#Read in a fresh new blade
blade=pynumad.Blade()
blade.read_yaml('example_data/'+yamlName+'.yaml') 
file_names=write_beam_model(wt_name,station_list,settings,blade,mu,log,directory=dir_name)