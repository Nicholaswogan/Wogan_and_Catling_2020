import numpy as np
import sys
sys.path.append('./python_scripts/')
from gen_inputs import *
import subprocess
import os
FNULL = open(os.devnull, 'w')
import time

# cycle through all of the photochemial outputs
# and calculate the gibbs free energy of each of
# each of them



# Loop through all files
# there are Q=1 and Q=25
QQ = range(1,26,2)
for Q in QQ:
    print('Q =',Q)
    start = time.time()
    direc = 'Volc_iter/Q'+str(int(Q))+'/'
    # #first do prebiotic case
    gen_inputs_atmos_only(direc+'PTZ_mixingratios_pre.dist','inputs_Volc_iter_atmos_only.txt')
    print('Calculating Delta_G for prebiotic model...')
    subprocess.call(['matlab', '-nodesktop', '-nosplash', '-r', "loadDatabase; loadDatabaseB; loadDatabaseC; loadDatabaseD; load_Pitzer; Main_script_iterate; exit"])#,stdout=FNULL, stderr=subprocess.STDOUT)
    end = time.time()
    print('Done. Iteration time:',end-start,'seconds')
    #save output
    subprocess.call(['mv','Gibbs_min_result.txt',direc+'Gibbs_result_pre_atmos_only.txt'])

    #first do prebiotic case
    start = time.time()
    gen_inputs_atmos_only(direc+'PTZ_mixingratios_eco.dist','inputs_Volc_iter_atmos_only.txt')
    print('Calculating Delta_G for ecosystem model...')
    subprocess.call(['matlab', '-nodesktop', '-nosplash', '-r', "loadDatabase; loadDatabaseB; loadDatabaseC; loadDatabaseD; load_Pitzer; Main_script_iterate; exit"])#,stdout=FNULL, stderr=subprocess.STDOUT)
    end = time.time()
    print('Done. Iteration time:',end-start,'seconds')
    #save output
    subprocess.call(['mv','Gibbs_min_result.txt',direc+'Gibbs_result_eco_atmos_only.txt'])
    print('\n')
