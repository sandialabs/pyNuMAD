import numpy as np
import os
import warnings
import matplotlib.pyplot as plt
from datetime import datetime
import subprocess
import time
import os,glob

from pynumad.analysis.ansys.utility import *
from pynumad.analysis.ansys.read import *
from pynumad.analysis.ansys.write import *
from pynumad.analysis.ansys.run import call_ansys

#from pynumad.analysis.utils.distributed_loading import 


    
def main_ansys_analysis(
        blade,
        mesh_data,
        loads_table,
        analysis_config,
        elementSize,
        log=False
    ):
    """
    Parameters
    ----------
    blade: Blade
        pynumad blade object
    mesh_data: dict
        pynumad mesh output
    loads_table: list
        list of loads
    analysis_config: dict
        #TODO define parameters
    varargin = None

    Returns
    -------
    """
    an_flag_names = analysis_config["analysisFlags"].keys()
    joined_flag_names=''.join(an_flag_names)
    if ('imperfection' in analysis_config["analysisFlags"]) and \
        analysis_config["analysisFlags"]["imperfection"] and \
        not analysis_config["analysisFlags"]["globalBuckling"]:
        raise Exception('Specify number of buckling modes when performing nonlinear buckling')
    
    # Original mesh file to analize
    if ('meshFile' not in analysis_config):
        analysis_config["meshFile"] = 'master'
    
    # File name base name for ansys analysis files
    if ('analysisFileName' in analysis_config):
        ansysFilename = analysis_config["analysisFileName"]
    else:
        ansysFilename = 'FEmodel'
    
    # Number of CPUs to use
    if 'np' in analysis_config and analysis_config["np"] > 0:
        ncpus = analysis_config["np"]
    else:
        ncpus = 1
    
    #Initialize
    results = dict()
    for key in an_flag_names:
        if key in ['globalBuckling', 'resultants', 'deflection', 'mass']:
            if analysis_config["analysisFlags"][key] != 0:
                results[key] = [None]*len(loads_table)
        elif key in [ 'failure', 'fatigue', 'imperfection', 'mass']:
            if analysis_config["analysisFlags"][key]:
                results[key] = [None]*len(loads_table)
    
    if not results:
        raise Exception('no analyses are configured in configuration st.')
    
    flags={}
    for flag_name in ['globalBuckling','imperfection','resultants','deflection']+ ['localBuckling','failure','fatigue','local_fields']:
        if flag_name in joined_flag_names and analysis_config['analysisFlags'][flag_name]:
            flags[flag_name]=analysis_config['analysisFlags'][flag_name]
        else:
            flags[flag_name]=False #Fillin missing flags with False

    # for flag_name in ['localBuckling','failure','fatigue','local_fields']:
    #     if flag_name in joined_flag_names and len(analysis_config['analysisFlags'][flag_name])>:
    #         flags[flag_name]=analysis_config['analysisFlags'][flag_name]
    #     else:
    #         flags[flag_name]=False #Fillin missing flags with False



    for iLoad, loads in enumerate(loads_table):
        if iLoad == 0 and 'mass' in joined_flag_names and analysis_config['analysisFlags']['mass']:
            mass_flag=True
        else:
            mass_flag=False
        ## ************************************************************************
        # ================= APPLY LOADS TO FEA MESH ================= #NOTE: Priority
        forcefilename = 'forces'
        # only want outershell VVV
        nodeData = np.concatenate([np.arange(mesh_data['nodes'].shape[0]).reshape((-1,1)),mesh_data['nodes']], axis=1)
        
        write_ansys_loads(nodeData, loads, forcefilename, analysis_config)

        ## ************************************************************************
        # ================= WRITE LINEAR STATIC ANALYSIS COMMANDS ================= #NOTE: Priority
        # run buckling computations in ansys
        print(' ')
        print('Running ANSYS analysis...')
        script_name = 'ansysAnalysis.mac'
        script_out = 'ansysAnalysisEcho.out'
        fid = open(script_name,'w+')
        fid.write('/NERR,,99999999\n' % ())
        fid.write('/CWD, %s\n' % (os.getcwd()))
        fid.write('resume,master,db\n' % ())
        #         fprintf(fid,'/FILNAME,''#s'',1\n',ansysFilename);   #From master, change the jobname
        fid.write('/FILNAME,%s,1\n' % (ansysFilename+'-Load'+str(iLoad)))
        #fprintf(fid,'resume\n');
        fid.write('! BEGIN LINEAR STATIC SCRIPT\n' % ())
        fid.write('esel,all\n' % ())
        fid.write('/prep7\n' % ())
        fid.write('fdel,all\n' % ())
        fid.write('/input,%s,src\n' % (forcefilename))
        #Linear Static Analysis
        fid.write('/solu\n' % ())
        fid.write('antype,static\n' % ())
        if 'StaticNonlinear' in analysis_config["analysisFlags"] and not \
            len(analysis_config["analysisFlags"].StaticNonlinear)==0  and \
            analysis_config["analysisFlags"].StaticNonlinear != 0:
            fid.write('nlgeom,1\n' % ())
            fid.write('OUTRES,all,ALL\n' % ())
        #         else
        #             fprintf(fid,'pstres,on\n');
        fid.write('irlf,-1\n' % ())
        fid.write('bcsoption,,incore\n' % ())
        fid.write('solve\n' % ())
        fid.write('finish\n' % ())
        #Only compute mass on the first load case
        if mass_flag:
            #Get Mass Here
            fid.write('*GET, Z_mass, ELEM, 0, MTOT, X\n' % ())
            fid.write('/output, results_mass,txt\n' % ())
            fid.write('*status,Z_mass\n' % ())
            fid.write('/output\n' % ())
            fid.write('finish\n' % ())
        ## ************************************************************************
        #================= WRITE Deflection ANALYSIS COMMANDS =================  
        if flags['deflection']:
            deflectionFilename = 'results_deflection'
            writeAnsysDeflections(blade,analysis_config,iLoad,fid,deflectionFilename)
        # # calculate face stresses for wrinkling NOTE: Skip
        # if flags['localBuckling']:
        #     #Check for wrinkling here in a linear analysis
        #     app,SkinAreas,compsInModel = writeAnsysGetFaceStresses(blade,fid,analysis_config["analysisFlags"].localBuckling)
        ### Output resultant force and moments to file NOTE: Skip
        if flags['resultants']:
            writeAnsysResultants(blade,analysis_config,iLoad,fid,elementSize)
        ## ************************************************************************
        # ================= WRITE FATIGUE ANALYSIS COMMANDS================= 
        if flags['fatigue']:
            writeAnsysFatigue(fid,iLoad)
        ## ************************************************************************
        # ================= CREATE LOCAL FIELD RESULTS FOR  ================= 
        if flags['local_fields']:
            writeAnsysLocalFields(blade,analysis_config,iLoad,fid)
        ## ************************************************************************
        # ================= WRITE FAILURE ANALYSIS COMMANDS ================= 
        # Initialize GUI commands from batch operation to identify maxima
        if flags['failure']:
            failureFilename = 'results_failure'
            writeAnsysRupture(analysis_config["analysisFlags"]['failure'],iLoad,fid,failureFilename)
        ## ************************************************************************
        # ================= WRITE BUCKLING ANALYSIS COMMANDS ================= 
        #Linear Buckling Analysis
        if flags['globalBuckling']:
            if analysis_config["analysisFlags"]['globalBuckling'] > 0:
                bucklingFilename = 'results_buckling'
                writeAnsysLinearBuckling(blade,analysis_config["analysisFlags"]['globalBuckling'],fid,bucklingFilename)
            else:
                raise Exception('analysis_config["analysisFlags"].globalBuckling must be greater than or equal to zero')
        fid.close()
        
        
        
        call_ansys(script_name,log,ncpus=analysis_config['np'])

    

        #   POST PROCESS ##########################################
        ## ************************************************************************
        # ================= READ MASS RESULTS  =================
        if mass_flag:
            results['mass'] = read_1_ANSYSoutput('results_mass.txt')
        ## ************************************************************************
        # ================= READ DEFLECTION RESULTS   =================
        if flags['deflection']:
            results['deflection'][iLoad] = read_ansys_deflections(blade,analysis_config,iLoad,deflectionFilename)
        ## ************************************************************************
        # ================= READ STRESS RESULTANTS   =================
        if flags['resultants']:
            file_name = 'results_resultants.txt'
            results['resultants'][iLoad] = txt2mat(file_name)

        ## ************************************************************************
        # ================= READ LINEAR BUCKLING RESULTS =================
        # read buckling results
        if flags['globalBuckling']:
            linearLoadFactors = readAnsysLinearBuckling(analysis_config["analysisFlags"]["globalBuckling"],bucklingFilename)
        ## ************************************************************************
        # ================= WRITE NON-LINEAR BUCKLING/WRINKLING ANALYSIS COMMANDS =================
        # Perform nonlinear buckling here if required (and writeANSYSgetFaceStresses
        # at the end of the nonlinear analysis for wrikling check
        if flags['imperfection']:
            warnings.warn('output results. Currently does not work for nonlinear cases')
            imperfection = analysis_config["analysisFlags"].imperfection / 1000
            nonlinearLoadFactors = np.zeros((len(linearLoadFactors),len(imperfection)))
            critDesignvar = np.zeros((len(imperfection),1))
            wrinklingLimitingElementData = np.zeros((len(linearLoadFactors),4,len(imperfection)))
            marker = np.array(['-ok','-sk','-dk','-*k','-^k','-<k','->k','-pk','-hk'])
            #SF=max(LLF); #Use one loads file for all buckling modes
            for jj in range(len(imperfection)):
                for ii in range(len(linearLoadFactors)):
                    # For each load factor, create a new jobname and database and run a nonlinear static analysis
                    nonlinearLoadFactors[ii,jj] = writeAnsysNonLinearBuckling(ansysFilename,ansysPath,ansys_product,analysis_config,ii,jj,ncpus,iLoad)
                    wrinklingLimitingElementData[ii,:,jj] = wrinklingForNonlinearBuckling(blade,analysis_config["analysisFlags"].localBuckling,settings,ncpus,ansysFilename,ii,jj)
                minnLLF,minnLLFMode = np.amin(nonlinearLoadFactors[:,jj])
                minWLF,minWLFMode = np.amin(wrinklingLimitingElementData[:,2,jj])
                critDesignvar[jj] = np.amin(minnLLF,minWLF)
            plt.figure(5)
            for k in range(len(linearLoadFactors)):
                #disp(strcat('-',marker(j),'k'))
                plt.plot(imperfection * 1000,nonlinearLoadFactors[k,:],marker[k])
            plt.legend('Mode-' + str(np.arange(len(linearLoadFactors))))
            plt.title('Imperfection Study (Linear Elements) SNL3p0-148-mk0p2-s1-fiberglass')
            plt.xlabel('Max Imperfection Size [mm]')
            plt.ylabel('Buckling Load Factors [ ]')
            #wrinklingLimitingElementData - [ansysSecNumber elno lf phicr]
            results.globalBuckling[iLoad] = np.amin(np.amin(critDesignvar))
        else:
            if flags['globalBuckling']:
                results['globalBuckling'][iLoad] = linearLoadFactors[0]
        ## ************************************************************************
        # ================= POST-PROCESS PANEL WRINKLING FACTORS =================
        # if flags['localBuckling']:
        #         #UNSUPPORTED AT THIS TIME
        #         writeAnsysNonLinearLocalBuckling(blade,analysis_config,iLoad,fid,ansysFilename,ii,jj)
        #     # perform wrinkling check
        #     wrinklingLimitingElementData = writeAnsysFagerberWrinkling(app,SkinAreas,compsInModel,analysis_config["analysisFlags"].localBuckling
        #     results.localBuckling[iLoad] = wrinklingLimitingElementData[2]
        #     os.path.delete('*faceAvgStresses.txt') # NOTE: I believe * is supposed to glob here, but I am not sure it is doing that -kb
        ## ************************************************************************
        # ================= READ FAILURE RESULTS   =================
        if flags['failure']:
            file_name = failureFilename+'.txt'
            results['failure'][iLoad] = read_1_ANSYSoutput(file_name)
            print('')
    
    ## ************************************************************************
    
    # # ================= RUN FATIGUE POST PROCESSOR =================
    # #After all load directions are solved compute fatige damage if needed
    # if flags['fatigue']:
    #     if not len(varargin)==0  and class_(varargin[0])=='IECDef':
    #         # cd ..
    #         IEC = varargin[0]
    #         wt,rccdata = getWindSpeedDistribution(IEC.avgws)
    #         # cd 'NuMAD'
    #         results.fatigue = postprocessANSYSfatigue(blade,mesh_data,wt,rccdata,IEC,loads_table,analysis_config)
    #     else:
    #         raise Exception('IECDef required to run fatigue analysis in main_ansys_analysis')
    for f in glob.glob("results*.txt"):
        os.remove(f)
    return results
    
    
def saveData(results = None,iLoad = None,airfoilSegmentName = None,iSpan = None,nodes = None,midNodei = None): 
    getattr[results.local_fields[iLoad],[airfoilSegmentName]].x[iSpan] = nodes[midNodei,1]
    getattr[results.local_fields[iLoad],[airfoilSegmentName]].y[iSpan] = nodes[midNodei,2]
    getattr[results.local_fields[iLoad],[airfoilSegmentName]].z[iSpan] = nodes[midNodei,3]
    getattr[results.local_fields[iLoad],[airfoilSegmentName]].data[iSpan] = nodes[midNodei,4]
    return results