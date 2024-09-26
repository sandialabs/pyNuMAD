import numpy as np
from pynumad.utils.interpolation import interpolator_wrap
import os
def readVABShomogenization(file_name):


    #Stiffness
    beam_stiff=np.zeros((6,6))
    with open(file_name) as f:
        lines=f.readlines()
        
    for lineNumber,line in enumerate(lines):
        if 'Timoshenko Stiffness Matrix' in line:
            break
    lineStart=lineNumber+3
    lineEnd=lineStart+6
    ct=0
    for iLine in range(lineStart,lineEnd):
        dataList=[float(i) for i in lines[iLine].split() if i.strip()]
        beam_stiff[ct,:]=dataList
        ct+=1

    
    #mass
    beam_inertia=np.zeros((6,6))
    for lineNumber,line in enumerate(lines):
        if 'Mass Matrix' in line:
            break
    lineStart=lineNumber+3
    lineEnd=lineStart+6
    ct=0
    for iLine in range(lineStart,lineEnd):
        dataList=[float(i) for i in lines[iLine].split() if i.strip()]
        beam_inertia[ct,:]=dataList
        ct+=1
  
    return beam_stiff, beam_inertia

def transform_bd_vectors_to_VABS(V):
    B = np.transpose(np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]]))  # transformation matrix
    return np.matmul(B,V)

def transformMatrixToBeamDyn(beam_stiff,beam_inertia):
    beamDynData={}

    B = np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]])  # NEW transformation matrix
    T = np.dot(np.identity(3), np.linalg.inv(B))
    
    nStations, _,_=np.shape(beam_stiff)

    for i_station in range(nStations):
        beam_stiff[i_station,:,:]=trsf_sixbysix(beam_stiff[i_station,:,:], T)
        beam_inertia[i_station,:,:]=trsf_sixbysix(beam_inertia[i_station,:,:], T)
   
    return(beam_stiff,beam_inertia)

def trsf_sixbysix(M, T):
    """
    Transform six-by-six compliance/stiffness matrix. 
    change of reference frame in engineering (or Voigt) notation.
    
    Parameters
    ----------
    M : np.ndarray
        6x6 Siffness or Mass Matrix
    T : np.ndarray
        Transformation Matrix
        
    Returns
    ----------
    res : np.ndarray
        Transformed 6x6 matrix
    """

    TS_1 = np.dot(np.dot(T.T, M[0:3, 0:3]), T)
    TS_2 = np.dot(np.dot(T.T, M[3:6, 0:3]), T)
    TS_3 = np.dot(np.dot(T.T, M[0:3, 3:6]), T)
    TS_4 = np.dot(np.dot(T.T, M[3:6, 3:6]), T)

    tmp_1 = np.vstack((TS_1, TS_2))
    tmp_2 = np.vstack((TS_3, TS_4))
    res = np.hstack((tmp_1, tmp_2))
    return res

# --- Write BeamDyn file with blade reference line locations ---#
def write_beamdyn_axis(directory, wt_name, blade):
    definition = blade.definition

    geometry = blade.geometry
    n_pts = 50
    input_span = blade.geometry.ispan/blade.geometry.ispan[-1]
    target_span = np.linspace(0, 1, n_pts)

    kp_xr=interpolator_wrap(input_span,definition.prebend,target_span,'pchip', axis=1)
    kp_yr=interpolator_wrap(input_span,definition.sweep,target_span,'pchip', axis=1)
    kp_zr=interpolator_wrap(input_span,geometry.ispan,target_span,'pchip', axis=1)
    twist_interp=interpolator_wrap(input_span,geometry.idegreestwist,target_span,'pchip', axis=1)


    data = np.vstack((kp_xr, kp_yr, kp_zr, twist_interp)).T

    if not os.path.exists(directory):
        os.makedirs(directory)

    axisFileName='bd_primary_'+wt_name + '.inp'
    
    file = open(directory +'/'+ axisFileName, 'w')
    file.write('--------- BEAMDYN with OpenFAST INPUT FILE -------------------------------------------\n')
    file.write('%s blade\n' % (wt_name))
    file.write('---------------------- SIMULATION CONTROL --------------------------------------\n')
    file.write('True          Echo            - Echo input data to "<RootName>.ech" (flag)\n')
    file.write('True          QuasiStaticInit - Use quasistatic pre-conditioning with centripetal accelerations in initialization (flag) [dynamic solve only]\n')
    file.write(' 0            rhoinf          - Numerical damping parameter for generalized-alpha integrator\n')
    file.write(' 2            quadrature      - Quadrature method: 1=Gaussian; 2=Trapezoidal (switch)\n')
    file.write(' 2            refine          - Refinement factor for trapezoidal quadrature (-). DEFAULT = 1 [used only when quadrature=2]\n')
    file.write('"DEFAULT"     n_fact          - Factorization frequency (-). DEFAULT = 5\n')
    file.write('"DEFAULT"     DTBeam          - Time step size (s).\n')
    file.write(' 50           load_retries    - Number of factored load retries before quitting the aimulation\n')
    file.write('"DEFAULT"     NRMax           - Max number of iterations in Newton-Ralphson algorithm (-). DEFAULT = 10\n')
    file.write('"DEFAULT"     stop_tol        - Tolerance for stopping criterion (-)\n')
    file.write('"DEFAULT"     tngt_stf_fd     - Flag to use finite differenced tangent stiffness matrix (-)\n')
    file.write('"DEFAULT"     tngt_stf_comp   - Flag to compare analytical finite differenced tangent stiffness matrix  (-)\n')
    file.write('"DEFAULT"     tngt_stf_pert   - perturbation size for finite differencing (-)\n')
    file.write('"DEFAULT"     tngt_stf_difftol- Maximum allowable relative difference between analytical and fd tangent stiffness (-)\n')
    file.write('True          RotStates       - Orient states in the rotating frame during linearization? (flag) [used only when linearizing]\n')
    file.write('---------------------- GEOMETRY PARAMETER --------------------------------------\n')
    file.write('          1   member_total    - Total number of members (-)\n')
    file.write('         %u   kp_total        - Total number of key points (-) [must be at least 3]\n' % (n_pts))
    file.write('     1     %u                 - Member number; Number of key points in this member\n' % (n_pts))
    file.write('\t\t kp_xr \t\t\t kp_yr \t\t\t kp_zr \t\t initial_twist\n')
    file.write('\t\t  (m)  \t\t\t  (m)  \t\t\t  (m)  \t\t   (deg)\n')


    for i in range(n_pts):
        file.write('\t %.5e \t %.5e \t %.5e \t %.5e \n' % (data[i, 0], data[i, 1], data[i, 2], data[i, 3]))

    file.write('---------------------- MESH PARAMETER ------------------------------------------\n')
    file.write('          10   order_elem     - Order of interpolation (basis) function (-)\n')
    file.write('---------------------- MATERIAL PARAMETER --------------------------------------\n')
    file.write('"%s"    BldFile - Name of file containing properties for blade (quoted string)\n' % (wt_name + '_BeamDyn_Blade.dat'))
    file.write('---------------------- PITCH ACTUATOR PARAMETERS -------------------------------\n')
    file.write('False         UsePitchAct - Whether a pitch actuator should be used (flag)\n')
    file.write('        200   PitchJ      - Pitch actuator inertia (kg-m^2) [used only when UsePitchAct is true]\n')
    file.write('      2E+07   PitchK      - Pitch actuator stiffness (kg-m^2/s^2) [used only when UsePitchAct is true]\n')
    file.write('     500000   PitchC      - Pitch actuator damping (kg-m^2/s) [used only when UsePitchAct is true]\n')
    file.write('---------------------- OUTPUTS -------------------------------------------------\n')
    file.write('True          SumPrint       - Print summary data to "<RootName>.sum" (flag)\n')
    file.write('"ES10.3E2"    OutFmt         - Format used for text tabular output, excluding the time channel.\n')
    file.write('          9   NNodeOuts      - Number of nodes to output to file [0 - 9] (-)\n')
    file.write('2,3, 6, 9, 12, 15, 18, 21, 50, 59    OutNd          - Nodes whose values will be output  (-)\n')
    file.write('          OutList            - The next line(s) contains a list of output parameters. See OutListParameters.xlsx for a listing of available output channels, (-)\n')

    # coordinate={}
    # coordinate['F']='l'
    # coordinate['M']='l'
    # coordinate['RD']='r'
    # coordinate['TD']='r'
    
    # channelList=['F','M','RD','TD']
    # for iNode in range(9):
    #     for load in channelList:
    #         for dir in ['x','y','z']:
    #             file.write(f'"N{iNode+1}{load}{dir}{coordinate[load]}"\n')

    
    # #Root
    # coordinate={}
    # coordinate['F']='r'
    # coordinate['M']='r'
    
    # channelList=['F','M']
    # for iNode in ['Root']:
    #     for load in channelList:
    #         for dir in ['x','y','z']:
    #             file.write(f'"{iNode}{load}{dir}{coordinate[load]}"\n')

    # #Tip
    # coordinate={}
    # coordinate['RD']='r'
    # coordinate['TD']='r'
    
    # channelList=['RD','TD']
    # for iNode in ['Tip']:
    #     for load in channelList:
    #         for dir in ['x','y','z']:
    #             file.write(f'"{iNode}{load}{dir}{coordinate[load]}"\n')
                

    file.write('END of input file (the word "END" must appear in the first 3 columns of this last OutList line)\n')
    file.write('---------------------- NODE OUTPUTS --------------------------------------------\n')
    file.write('         99   BldNd_BlOutNd   - Blade nodes on each blade (currently unused)\n')
    file.write('              OutList     - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx, BeamDyn_Nodes tab for a listing of available output channels, (-)\n')
    file.write('"FxL"\n')
    file.write('"FyL"\n')
    file.write('"FzL"\n')
    file.write('"MxL"\n')
    file.write('"MyL"\n')
    file.write('"MzL"\n')
    file.write('"TDxr"\n')
    file.write('"TDyr"\n')
    file.write('"TDzr"\n')
    file.write('"RDxr"\n')
    file.write('"RDyr"\n')
    file.write('"RDzr"\n')
    file.write('"AbsXr"\n')
    file.write('"AbsYr"\n')
    file.write('"AbsZr"\n')
    file.write('END of input file (the word "END" must appear in the first 3 columns of this last OutList line)\n')
    file.write('---------------------------------------------------------------------------------------\n')
    file.write('PFxL\n')
    file.write('PFyL\n')
    file.write('PFzL\n')
    file.write('PMxL\n')
    file.write('PMyL\n')
    file.write('PMzL\n')
    file.write('DFxL\n')
    file.write('DFyL\n')
    file.write('DFzL\n')
    file.write('DMxL\n')
    file.write('DMyL\n')
    file.write('DMzL\n')
    file.write('DFxR\n')
    file.write('DFyR\n')
    file.write('DFzR\n')
    file.write('DMxR\n')
    file.write('DMyR\n')
    file.write('DMzR\n')

    file.close()

    print('Finished writing BeamDyn File')

    return axisFileName

# --- Write BeamDyn_Blade file with blade properties ---#
def write_beamdyn_prop(folder, wt_name, radial_stations, beam_stiff, beam_inertia, mu):
    n_pts = len(radial_stations)

    if not os.path.exists(folder):
        os.makedirs(folder)
        
    propFileName= 'bd_props_'+wt_name + '.inp'
    
    
    file = open(folder +'/'+propFileName, 'w')
    file.write(' ------- BEAMDYN V1.00.* INDIVIDUAL BLADE INPUT FILE --------------------------\n')
    file.write(' Test Format 1\n')
    file.write(' ---------------------- BLADE PARAMETERS --------------------------------------\n')
    file.write('%u   station_total    - Number of blade input stations (-)\n' % (n_pts))
    file.write(' 1   damp_type        - Damping type: 0: no damping; 1: damped\n')
    file.write('  ---------------------- DAMPING COEFFICIENT------------------------------------\n')
    file.write('   mu1        mu2        mu3        mu4        mu5        mu6\n')
    file.write('   (-)        (-)        (-)        (-)        (-)        (-)\n')
    file.write('\t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e\n' % (mu[0], mu[1], mu[2], mu[3], mu[4], mu[5])) 
    file.write(' ---------------------- DISTRIBUTED PROPERTIES---------------------------------\n')
    
    for i in range(n_pts):
        file.write('\t %.6f \n' % (radial_stations[i]))
        # write stiffness matrices
        for j in range(6):
            file.write('\t %.16e \t %.16e \t %.16e \t %.16e \t %.16e \t %.16e\n' % (
            beam_stiff[i, j, 0], beam_stiff[i, j, 1], beam_stiff[i, j, 2], beam_stiff[i, j, 3], beam_stiff[i, j, 4],
            beam_stiff[i, j, 5]))
        file.write('\n')

        # write inertia properties
        for j in range(6):
            file.write('\t %.16e \t %.16e \t %.16e \t %.16e \t %.16e \t %.16e\n' % (
            beam_inertia[i, j, 0], beam_inertia[i, j, 1], beam_inertia[i, j, 2], beam_inertia[i, j, 3],
            beam_inertia[i, j, 4], beam_inertia[i, j, 5]))
        file.write('\n')
        # ToDO: check correct translation of stiffness and mass matrices from VABS and anbax !!!
    file.close()

    print('Finished writing BeamDyn_Blade File')

    return propFileName


def writeBeamDynStandAlone(file_names,disrLoads,tipLoads,directory='.'):

    if not os.path.exists(directory):
        os.makedirs(directory)

    from pynumad.utils.misc_utils import copy_and_replace

    templateFileName='beamDynStandAlone.template'
    
    analysisFileName=file_names[0]+'_driver.inp'

    path_name=directory+'/'+analysisFileName



    
    copy_and_replace(templateFileName, path_name,
        {
            'DISTRLOAD1' : str(disrLoads[0]),
            'DISTRLOAD2' : str(disrLoads[1]),
            'DISTRLOAD3' : str(disrLoads[2]),
            'DISTRLOAD4' : str(disrLoads[3]),
            'DISTRLOAD5' : str(disrLoads[4]),
            'DISTRLOAD6' : str(disrLoads[5]),
            'TIPLOAD1' : str(tipLoads[0]),
            'TIPLOAD2' : str(tipLoads[1]),
            'TIPLOAD3' : str(tipLoads[2]),
            'TIPLOAD4' : str(tipLoads[3]),
            'TIPLOAD5' : str(tipLoads[4]),
            'TIPLOAD6' : str(tipLoads[5]),
            'AXIS FILE NAME': file_names[0],
        })
    return analysisFileName

def runBeamDynStandAlone(beamDynDriverFileName,log,directory='.'):
    from pynumad import path_data
    import subprocess
    try:
        this_cmd = path_data['openFast']+'beamdyn_driver '+directory+'/'+beamDynDriverFileName
        log.info(f' running: {this_cmd}')
        subprocess.run(this_cmd, shell=True, check=True, capture_output=True)

        # with open(filePath+'.ech', 'r') as f:
        #     lines = f.readlines()
        # #log the last line of .ech file:
        # log.error(f'****************************\n{lines[-1]}\n******************************')

    except subprocess.CalledProcessError as e:
        log.error(f'Error running {this_cmd}: {e}')

def read_fast_out(filename):

    if 'outb' in filename:
        raise ValueError('FAST binary output not yet supported')
    
    with open(filename) as f:
        lines=f.readlines()

    for iline,line in enumerate(lines):  
        if 'Time' in line:
            break
        
    #Initialize empty dict
    data={val.strip():[] for val in lines[iline].split('\t')}

    for line in lines[iline+2:]:
        vals = line.split('\t')
        for ikey, key in enumerate(data.keys()):
            data[key].append(float(vals[ikey]))
    return data
            
def get_bd_source_radial_stations(bd_sum_file_name):
    '''
    Returns 
        nnodes, number of quadrature points (qp)
        node_location, nnodes x 3 array of global coordinates
    '''
    import yaml

    with open(bd_sum_file_name) as blade_yaml:
        data = yaml.load(blade_yaml,Loader=yaml.FullLoader)
    node_location=list(np.array(data['Init_Nodes_E1'])[:,2]/np.array(data['Init_Nodes_E1'])[-1,2])
    return node_location

# def get_bd_source_radial_stations(qp_location,blade_length,nqp): #Adds root and top to qp data

#     source_radial_stations=[0]
#     for i in range(nqp):
#         source_radial_stations.append(qp_location[i,2]/blade_length)
#     source_radial_stations.append(1.0)
#     return source_radial_stations

def get_bd_target_radial_stations(bd_props_file_name):
    with open(bd_props_file_name) as f:
        lines=f.readlines()

    for iline,line in enumerate(lines):  
        if 'properties' in line.lower():
            line_start=iline+1
            break

    target_radial_station = []
    for iline,line in enumerate(lines[line_start:]): 
        if len(line.split('\t')) < 5:
            try:
                target_radial_station.append(float(line.strip()))
            except :
                pass

    return target_radial_station

# def get_bd_spanwise_source_resultants_at_time(out_data,nqp,time_index):
#     from pynumad.utils.misc_utils import full_keys_from_substrings

#     resultant_labels=['Fx','Fy','Fz','Mx','My','Mz']


#     all_resultants=[]

#     bd_frame='r' 
#     key_name_prepend='Root'
#     resultant=[]
#     for ir,resultant_label in enumerate(resultant_labels):
#         sub_string = [key_name_prepend,resultant_label,bd_frame]

#         key_name=full_keys_from_substrings(out_data.keys(),sub_string,ignore_case=False)[0]
#         resultant.append(out_data[key_name][time_index])
#     all_resultants.append(resultant)

#     bd_frame='L' #case matters
#     key_name_prepend='N'
#     for qp in range(1,nqp+1):
#         resultant=[]
#         for ir,resultant_label in enumerate(resultant_labels):
#             sub_string = [key_name_prepend+str(qp).zfill(3),resultant_label,bd_frame]

#             key_name=full_keys_from_substrings(out_data.keys(),sub_string,ignore_case=False)[0]
#             resultant.append(out_data[key_name][time_index])
#         all_resultants.append(resultant)


#     #Append zeros to tip
#     resultant=[]
#     for ir,resultant_label in enumerate(resultant_labels):
#         resultant.append(0.0)
#     all_resultants.append(resultant)
#     return all_resultants
def get_bd_spanwise_source_qois_at_time(time_index,out_data,qoi,nnodes,bd_frame):
    from pynumad.utils.misc_utils import full_keys_from_substrings

    if 'l' in bd_frame.lower():
        bd_frame=bd_frame.upper() #case matters
    elif 'r' in bd_frame.lower():
        bd_frame=bd_frame.lower() #case matters
    else:
        raise ValueError(f'unknown BeamDyn coordinate frame: {bd_frame}')
    
    if 'forces' in qoi.lower():
        data_labels=['F']
    elif 'moments' in qoi.lower():
        data_labels=['M']
    elif 'disp' in qoi.lower():
        data_labels=['TD']
    elif 'rot' in qoi.lower():
        data_labels=['RD']
    else:
        raise ValueError('Unknown QOI: {qoi}')

    spanwise_source_qois=[]

    key_name_prepend='N'
    for inode in range(1,nnodes+1):
        source_qois=[]
        for ir,source_qois_label in enumerate(data_labels):
            for direction in ['x','y','z']:
                sub_string = [key_name_prepend+str(inode).zfill(3),source_qois_label+direction+bd_frame]

                key_name=full_keys_from_substrings(out_data.keys(),sub_string,ignore_case=False)
                
                if len(key_name)==1:
                    key_name=key_name[0]
                elif len(key_name)==0:
                    raise ValueError(f'The strings: {sub_string} could not be found in the BeamDyn/OpenFAST results file.')
                else:
                    raise ValueError(f'Too many variables were found when searching BeamDyn/OpenFAST results file for {sub_string} ')
                
                

                source_qois.append(out_data[key_name][time_index])
        spanwise_source_qois.append(source_qois)
    return spanwise_source_qois
def get_bd_blade_length(bd_sum_file_name):
    '''
    Returns 
        blade length
    '''
    import yaml

    with open(bd_sum_file_name) as blade_yaml:
        data = yaml.load(blade_yaml,Loader=yaml.FullLoader)
    return data['Length']

def get_spanwise_target_qois_at_one_time(source_radial_stations,target_radial_stations,spanwise_source_qois):
    nqois=len(spanwise_source_qois[0])
    target_qois=[]
    for target_radial_station in np.array(target_radial_stations):
        target_qoi=[]
        for i in range(nqois): #three forces and three moments or three displacements or three rotations
            
            target_qoi.append(float(interpolator_wrap(source_radial_stations,np.transpose(np.array(spanwise_source_qois))[i],target_radial_station,'pchip', axis=1)))
        target_qois.append(target_qoi)
    return target_qois

def get_dcm_from_WM_params(wm):

    '''
    Gets the direction cosine matrix (dcm) from Wiener-Milenkovic (WM) parameters
     
    This dcm transforms tensors from the b frame to the B frame. These 
    frames are defined in: 

    Hodges, D. H. (2006). Nonlinear composite beam theory. American Institute of Aeronautics and Astronautics.

    
    Formulas are from https://github.com/OpenFAST/openfast/issues/10
    '''
    dcms=[]
    for i in range(len(wm)):
        c=np.array(wm[i])
        c0 = 2.0 - 1.0/8.0*np.matmul(np.transpose(c),c)
        rx=c[0]; ry=c[1]; rz=c[2]

        
        if c0 < 2:
            R11 = (c0**2 + rx**2 - ry**2 - rz**2)
            R12 = (2.0*(rx*ry - c0*rz))
            R13 = (2.0*(rx*rz + c0*ry))

            R21 = (2.0*(rx*ry + c0*rz))
            R22 = (c0**2 - rx**2 + ry**2 - rz**2)
            R23 = (2.0*(ry*rz - c0*rx))

            R31 = (2.0*(rx*rz - c0*ry))
            R32 = (2.0*(ry*rz + c0*rx))
            R33 = (c0**2 - rx**2 - ry**2 + rz**2)

            R = np.zeros((3, 3))

            R[0,0] = R11
            R[0,1] = R12
            R[0,2] = R13

            R[1,0] = R21
            R[1,1] = R22
            R[1,2] = R23

            R[2,0] = R31
            R[2,1] = R32
            R[2,2] = R33

            scale_factor = 1.0/(4.0 - c0)**2
            R = scale_factor * R
            dcm=np.transpose(R)
        else:
            dcm = np.identity(3) #No rotations (mainly for root)
        dcms.append(dcm)
    return dcms

def write_vabs_glb_input(folder,file_name_base,u,dcm,F,M,f,m,fp,mp,fpp,mpp,fppp,mppp):
    
    file = open(f'{folder}/{file_name_base}.glb', 'w')
    file.write(f'{u[0]} {u[1]} {u[2]}\n')
    file.write(f'{dcm[0,0]} {dcm[0,1]} {dcm[0,2]}\n')
    file.write(f'{dcm[1,0]} {dcm[1,1]} {dcm[1,2]}\n')
    file.write(f'{dcm[2,0]} {dcm[2,1]} {dcm[2,2]}\n')
    file.write(f'{F[0]} {M[0]} {M[1]} {M[2]}\n')
    file.write(f'{F[1]} {F[2]}\n')
    file.write(f'{f[0]} {f[1]} {f[2]} {m[0]} {m[1]} {m[2]}\n')
    file.write(f'{fp[0]} {fp[1]} {fp[2]} {mp[0]} {mp[1]} {mp[2]}\n')
    file.write(f'{fpp[0]} {fpp[1]} {fpp[2]} {mpp[0]} {mpp[1]} {mpp[2]}\n')
    file.write(f'{fppp[0]} {fppp[1]} {fppp[2]} {mppp[0]} {mppp[1]} {mppp[2]}\n')
    file.close()
    return

def run_vabs_recovery(dir_name,file_name):
    import subprocess
    from pynumad.paths import SOFTWARE_PATHS

    this_cmd = f'{SOFTWARE_PATHS["vabs"]} {dir_name}/{file_name} 1'
    try:
        print(f'Running: {this_cmd}')
        subprocess.run(this_cmd, shell=True, check=True, capture_output=True)
    except:

        raise RuntimeError(f'The command: {this_cmd} failed.')
    return