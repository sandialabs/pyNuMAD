#import logging
import subprocess
#import os
import glob
import numpy as np
from pynumad.utils.misc_utils import copy_and_replace
from pynumad.paths import SOFTWARE_PATHS
import os


def write_beam_model(wt_name,station_list,settings,blade,mu,log,directory='.'):
    import pynumad.analysis.beam_utils as beam_utils

    #     #Runs VABS or OpenSG to homogenize
    #     #Makes beamDyn or GEBT files
    geometry = blade.geometry


    # #Run input files

    if 'vabs' in settings['make_input_for'].lower():

        log.info(f'\n\n\nRunning VABS homogenization.')
        
        fileCount=0
        #First remove any lck files
        pattern=directory+'/'+wt_name+'*.in'
        if len(glob.glob(pattern))==0:
            raise RuntimeError(f'Could not find files with pattern: {pattern}. Beam model generation failed')
        MAXnLicenceTries=100
        for filePath in glob.glob(directory+'/'+wt_name+'*.in'):
            fileCount+=1
            try:
                this_cmd = SOFTWARE_PATHS['vabs']+' ' +filePath
                log.info(f' running: {this_cmd}')

                licenseAvailable=False
                nLicenceTries=0
                while not licenseAvailable and nLicenceTries <=MAXnLicenceTries-1:
                    subprocess.run(this_cmd, shell=True, check=True, capture_output=True)

                    with open(filePath+'.ech', 'r') as f:
                        lines = f.readlines()
                    #log the last line of .ech file:
                    
                    if 'Congratulations! No errors' in lines[-1]:
                        log.info(f'****************************\n{lines[-1]}\n******************************')
                        licenseAvailable=True
                        nLicenceTries=0
                    elif 'license' in lines[-1].lower():
                        nLicenceTries+=1
                        log.info(f'****************************\nnLicenceTries: {nLicenceTries}, {lines[-1]}\n******************************')

                    else:
                        log.error(f'****************************\n{lines[-1]}\n******************************')
                        raise Exception(f'Cross-sectional homogenization for file {filePath} failed due to: \n {lines[-1]} \n Beam model creation failed.') 
                if nLicenceTries ==MAXnLicenceTries:
                        string=f'License failed to be obtained after {MAXnLicenceTries} tries. Beam model creation failed.'
                        log.error(string)
                        raise Exception(string) 

            except subprocess.CalledProcessError as e:
                log.error(f'Error running {this_cmd}: {e}')
        
        # if fileCount != nStations:
        #     raise Exception('Error. Not enough VABS input files:')

    elif 'anba' in settings['make_input_for'].lower():
        raise ValueError('ANBA currently not supported')



### Read inputs
    extension='K'
    if station_list is None or len(station_list) == 0:
        station_list = list(range(len(geometry.ispan)))

    radial_stations=geometry.ispan/geometry.ispan[-1]
    radial_stations=radial_stations[station_list]

    if round(radial_stations[-1],2) ==1.0:
        radial_stations[-1]=1.0
    else:
        raise ValueError('The last radial station should be 1.0')
    beam_stiff = np.zeros([len(radial_stations), 6, 6])
    beam_inertia = np.zeros([len(radial_stations), 6, 6])

    for i_station, station in enumerate(station_list):
        file_name=f'{directory}/{wt_name}-{station}-t-0.in.K'
        beam_stiff[i_station,:,:],beam_inertia[i_station,:,:]=beam_utils.readVABShomogenization(file_name)
    # if len(glob.glob(directory+'/'+wt_name+"*." +extension)) >0:
    #     for file_name in glob.glob(directory+'/'+wt_name+"*." +extension):
    #         i_station=int(file_name.split('-')[-3].split('.')[0])
    #         print(f'file_name {file_name} i_station {i_station}')
    #         if i_station in range(len(station_list)):
    #             beam_stiff[i_station,:,:],beam_inertia[i_station,:,:]=beam_utils.readVABShomogenization(file_name)
    # else:
    #     raise RuntimeError(f'No VABS homogenization files found in directory: {directory}')
    

    if 'beamdyn' in settings['make_input_for'].lower():
        beam_stiff,beam_inertia=beam_utils.transformMatrixToBeamDyn(beam_stiff,beam_inertia)
        axisFileName=beam_utils.write_beamdyn_axis(directory, wt_name, blade)
        propFileName=beam_utils.write_beamdyn_prop(directory, wt_name, radial_stations, beam_stiff, beam_inertia, mu)
    return [axisFileName,propFileName]


def write_sierra_sm_model(template_file,wt_name,station_list,blade,materials_used,directory='.'):
    if not os.path.exists(directory):
        os.makedirs(directory)

    if station_list is None or len(station_list) == 0:
        station_list = list(range(len(blade.ispan)))
    elif len(station_list) == 1:
        raise ValueError("Need more than one cross section to make a solid model")
    # radial_stations=blade.ispan/blade.ispan[-1]
    # nStations=len(radial_stations)
    # #Run input files
    
    materials = blade.definition.materials

    sierra_file_name=f'{directory}/sm_{wt_name}.i'

    materialLines=f''
    blockLines=f''
    user_output_lines=f''

    for material_name in materials_used:
        material=materials[material_name]
        print(material.name)
        materialLines+=f'begin property specification for material {material.name}\n'
        materialLines+=f'   DENSITY      = {material.density}\n'
        materialLines+='    begin parameters for model elastic_orthotropic\n'
        materialLines+='        youngs modulus = 70e9\n'
        materialLines+='        poissons ratio  = 0.33\n'
        materialLines+=f'        E11          = {material.ex}\n'
        materialLines+=f'        E22          = {material.ey}\n'
        materialLines+=f'        E33          = {material.ez}\n'
        materialLines+=f'        NU12         = {material.prxy}\n'
        materialLines+=f'        NU13         = {material.prxz}\n'
        materialLines+=f'        NU23         = {material.pryz}\n'
        materialLines+=f'        G12          = {material.gxy}\n'
        materialLines+=f'        G13          = {material.gxz}\n'
        materialLines+=f'        G23          = {material.gyz}\n'

        materialLines+='\n' 
        materialLines+='        # Coordinate system\n'
        materialLines+='        COORDINATE SYSTEM = sysR\n' 
        materialLines+='    end parameters for model elastic_orthotropic\n' 
        materialLines+=f'end property specification for material {material.name}\n' 
        materialLines+='\n\n' 
        

        blockLines+=f'begin parameters for block {material.name}\n'
        blockLines+=f'    material {material.name}\n'
        blockLines+=f'    solid mechanics use model elastic_orthotropic\n'
        blockLines+=f'    section = hex_section\n'
        blockLines+=f'end parameters for block {material.name}\n'
        blockLines+='\n\n' 

        user_output_lines+=f'begin user output\n'
        user_output_lines+=f'    block = {material.name}\n'
        user_output_lines+=f'    extrapolate element variable mat_cauchy_stress_xx to nodal variable s_11_{material.name}\n'
        user_output_lines+=f'    extrapolate element variable mat_cauchy_stress_yy to nodal variable s_22_{material.name}\n'
        user_output_lines+=f'    extrapolate element variable mat_cauchy_stress_zz to nodal variable s_33_{material.name}\n'
        user_output_lines+=f'    extrapolate element variable mat_cauchy_stress_yz to nodal variable s_23_{material.name}\n'
        user_output_lines+=f'    extrapolate element variable mat_cauchy_stress_zx to nodal variable s_13_{material.name}\n'
        user_output_lines+=f'    extrapolate element variable mat_cauchy_stress_xy to nodal variable s_12_{material.name}\n\n'

        user_output_lines+=f'    extrapolate element variable mat_log_strain_xx to nodal variable e_11_{material.name}\n'
        user_output_lines+=f'    extrapolate element variable mat_log_strain_yy to nodal variable e_22_{material.name}\n'
        user_output_lines+=f'    extrapolate element variable mat_log_strain_zz to nodal variable e_33_{material.name}\n'
        user_output_lines+=f'    extrapolate element variable mat_log_strain_yz to nodal variable e_23_{material.name}\n'
        user_output_lines+=f'    extrapolate element variable mat_log_strain_zx to nodal variable e_13_{material.name}\n'
        user_output_lines+=f'    extrapolate element variable mat_log_strain_xy to nodal variable e_12_{material.name}\n'
        user_output_lines+=f'end user output\n'

    copy_and_replace(template_file, sierra_file_name,
        {
            'ROOT_STATION': 'station'+str(station_list[0]).zfill(3),
            'BLADE_MATERIALS': materialLines,
            'IN_MESH':wt_name+'.g',
            'OUT_MESH':wt_name+'.e',
            'BLADE_BLOCKS': blockLines,
            'BEGIN_USER_OUTPUT': user_output_lines
        })
    
def write_sierra_sd_model(template_file,wt_name,station_list,blade,materials_used,directory='.'):
    if station_list is None or len(station_list) == 0:
        station_list = list(range(len(blade.ispan)))
    elif len(station_list) == 1:
        raise ValueError("Need more than one cross section to make a solid model")
    
    if not os.path.exists(directory):
        os.makedirs(directory)
        
    materials = blade.definition.materials

    sierra_file_name=f'{directory}/sd_{wt_name}.i'

    materialLines=f''
    blockLines=f''
    for material_name in materials_used:
        material=materials[material_name]
        print(material.name)
        materialLines+=f'material {material.name}\n'
        materialLines+=f'orthotropic_prop\n'
        materialLines+=f'    E1          = {material.ex}\n'
        materialLines+=f'    E2          = {material.ey}\n'
        materialLines+=f'    E3          = {material.ez}\n'
        materialLines+=f'    nu12         = {material.prxy}\n'
        materialLines+=f'    nu23         = {material.pryz}\n'
        materialLines+=f'    nu13         = {material.prxz}\n'
        materialLines+=f'    G12          = {material.gxy}\n'
        materialLines+=f'    G23          = {material.gyz}\n'
        materialLines+=f'    G13          = {material.gxz}\n'
        materialLines+=f'    density       = {material.density}\n'
        materialLines+='\n' 

        materialLines+=f'end \n' 
        materialLines+='\n\n' 
        

        blockLines+=f'block {material.name}\n'
        blockLines+=f'    material {material.name}\n'
        blockLines+=f'end \n'
        blockLines+='\n\n' 

    copy_and_replace(template_file, sierra_file_name,
        {
            'BLADE_MATERIALS': materialLines,
            'IN_MESH':wt_name+'.g',
            'BLADE_BLOCKS': blockLines,
        })


