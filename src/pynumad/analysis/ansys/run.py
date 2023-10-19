
import subprocess
from pynumad.paths import SOFTWARE_PATHS
import glob,os
def call_ansys(script_name,log=False,script_out='ansys.out',ncpus=1):
    for f in glob.glob("*.lock"):
        os.remove(f)
    ansys_path=SOFTWARE_PATHS['ansys']
    MAXnLicenceTries=100
    try:

        this_cmd = f'export KMP_STACKSIZE=2048k & {ansys_path} -b -p ANSYS -I {script_name} -o {script_out} -np {str(ncpus)}' 
        if log:
            log.info(f' running: {this_cmd}')

        licenseAvailable=False
        nLicenceTries=0
        while not licenseAvailable and nLicenceTries <=MAXnLicenceTries-1:
            subprocess.run(this_cmd, shell=True, check=True, capture_output=True)

            with open(script_out, 'r') as f:
                lines = f.readlines()
            if 'lapsed time' in ''.join(lines):
                licenseAvailable=True
                if log:
                    log.info(f' Complete: {this_cmd}')
            # #log the last line of .ech file:
            
            # if 'Congratulations! No errors' in lines[-1]:
            #     log.info(f'****************************\n{lines[-1]}\n******************************')
            #     licenseAvailable=True
            #     nLicenceTries=0
            # elif 'license' in lines[-1].lower():
            #     nLicenceTries+=1
            #     log.info(f'****************************\nnLicenceTries: {nLicenceTries}, {lines[-1]}\n******************************')

            # else:
            #     log.error(f'****************************\n{lines[-1]}\n******************************')
            #     raise Exception(f'Cross-sectional homogenization for file {filePath} failed due to: \n {lines[-1]} \n Beam model creation failed.') 
            if nLicenceTries ==MAXnLicenceTries:
                    string=f'License failed to be obtained after {MAXnLicenceTries} tries. ANSYS model creation failed.'
                    if log:
                        log.error(string)
                    raise Exception(string) 

    except subprocess.CalledProcessError as e:
        with open(script_out, 'r') as f:
                lines = f.readlines()
        if log:
            log.error(''.join(lines))
            log.error(f'Error running {this_cmd}: {e}')

