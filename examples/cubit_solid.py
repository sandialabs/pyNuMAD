import sys 
import pynumad

sys.path.append(pynumad.SOFTWARE_PATHS['cubit'])
sys.path.append(pynumad.SOFTWARE_PATHS['cubitEnhancements'])

import cubit
from pynumad.analysis.cubit.make_blade import *
import numpy as np

from pynumad.analysis.make_models import write_sierra_model



def get_cs_params():
    cs_params = {}
    cs_params['nel_per_layer'] = 3 
    cs_params['element_ar'] = 5
    cs_params['element_shape'] = 'quad'


    cs_params['layer_transition_angle'] = 30
    cs_params['birds_mouth_amplitude_fraction']=1.0
    cs_params['minimum_layer_thickness'] = 0.001



    cs_params['adhesive_mat_name'] = 'Adhesive'


    totalStations = np.asarray(blade.ispan).size


    rootAdhesiveThickness=0.001
    tipAdhesiveThickness=0.001
    teAdhesiveWidthPercentChord=5
    adhesiveArray=np.linspace(rootAdhesiveThickness, tipAdhesiveThickness, num=totalStations)

    cs_params['web_fore_adhesive_thickness'] = adhesiveArray
    cs_params['web_aft_adhesive_thickness'] = adhesiveArray  
    cs_params['LE_adhesive_thickness'] = adhesiveArray
    cs_params['TE_adhesive_thickness'] = adhesiveArray
    cs_params['web_adhesive_width']=np.zeros((totalStations,))
    cs_params['TE_adhesive_width']=np.zeros((totalStations,))
    cs_params['max_web_imperfection_distance']=np.zeros((totalStations,))
    cs_params['minimum_layer_transition_length']=np.zeros((totalStations,))

    #Spanwise varying parameters
    thicknessScaling=0.001
    geometryScaling=thicknessScaling*1000
    for iStation in range(totalStations):
    
        if blade.definition.leband[iStation]!=0.0 and cs_params['LE_adhesive_thickness'][iStation] > blade.definition.leband[iStation]*0.85/1000: #Adhesive thickness can't be greater than 85% of the LE band
            raise ValueError(f'LE adhesive thickness of {cs_params["LE_adhesive_thickness"][iStation]} for station {iStation} is too large for the specified LE band width of {blade.definition.leband[iStation]/1000}.')
    

        #Adhesive width parameters
        cs_params['web_adhesive_width'][iStation] = 0.03 * blade.geometry.ichord[iStation]
        cs_params['TE_adhesive_width'][iStation] = teAdhesiveWidthPercentChord/100.0 * blade.geometry.ichord[iStation] * geometryScaling
        
        cs_params['max_web_imperfection_distance'][iStation] = 0.0001 * blade.geometry.ichord[iStation] * geometryScaling
        
        #Parameters that vary with span but are not part of the study 
        cs_params['minimum_layer_transition_length'][iStation]=blade.geometry.ichord[iStation]*0.002*geometryScaling;
    return cs_params 

blade=pynumad.Blade()

yamlName='myBlade_Modified'
blade.read_yaml('example_data/'+yamlName+'.yaml') 




wt_name=yamlName
dirName='.'


cs_params=get_cs_params()
settings={}
settings['make_input_for']='Sd'  #SM, VABS, ANBA, or None
settings['export']='cubg' #cub, g, or None

materialsUsed=cubit_make_solid_blade(blade, wt_name, settings, cs_params, stationList=[2,3])
    

write_sierra_model(wt_name,settings,blade,materialsUsed,'.') 