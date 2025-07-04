import sys 
import pynumad

sys.path.append(pynumad.SOFTWARE_PATHS['cubit'])
sys.path.append(pynumad.SOFTWARE_PATHS['cubit_enhancements'])

import cubit
from pynumad.analysis.cubit.make_blade import *
import numpy as np

from pynumad.analysis.make_models import write_sierra_sm_model
from pynumad.analysis.make_models import write_sierra_sd_model


def get_cs_params():
    cs_params = {}
    cs_params['nel_per_layer'] = 1 
    cs_params['element_ar'] = 5
    cs_params['element_shape'] = 'hex'
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
dirName='.'


cs_params=get_cs_params()
settings={}
settings['make_input_for']='SmSd'  #SM, VABS, ANBA, or None
settings['export']='cubg' #cub, g, or None


#Make Cubit Geometry
station_list = [2,3]
materials_used, volume_dict=cubit_make_solid_blade(blade, wt_name, settings, cs_params, stationList=station_list)

#Compute material orientation
orientation_vectors=get_material_orientation_vectors(volume_dict,ncpus = 1)
orientation_angles=get_material_orientation_angles(orientation_vectors)

#assign material orientation in Cubit
assign_material_orientation_vectors(orientation_vectors)
assign_material_orientation_angles(orientation_angles)

#Export mesh in Genisis format 
cubit.cmd(f'export mesh "{wt_name}.g" overwrite')


#Write Sierra input file
from pynumad.paths import SOFTWARE_PATHS
template_path=SOFTWARE_PATHS['pynumad']+'src/pynumad/data/templates/'

write_sierra_sm_model(template_path+'sm.i.template',wt_name,station_list,blade,materials_used,'.') 

write_sierra_sd_model(template_path+'sd.i.template',wt_name,station_list,blade,materials_used,'.') 