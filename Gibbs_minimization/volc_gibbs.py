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

QQ = range(1,26,2)
for Q in QQ:
    start = time.time()
    print('Q =',Q)
    direc = 'Volc_iter/Q'+str(int(Q))+'/'
    #first do prebiotic case
    gen_inputs(direc+'PTZ_mixingratios_pre.dist','inputs_Volc_iter.txt')
    print('Calculating Delta_G for prebiotic model...')
    subprocess.call(['matlab', '-nodesktop', '-nosplash', '-r', "loadDatabase; loadDatabaseB; loadDatabaseC; loadDatabaseD; load_Pitzer; Main_script_iterate; exit"])#,stdout=FNULL, stderr=subprocess.STDOUT)
    end = time.time()
    print('Done. Iteration time:',end-start)
    # save output
    subprocess.call(['mv','Gibbs_min_result.txt',direc+'Gibbs_result_pre.txt'])

    #first do prebiotic case
    start = time.time()
    gen_inputs(direc+'PTZ_mixingratios_eco.dist','inputs_Volc_iter.txt')
    print('Calculating Delta_G for ecosystem model...')
    subprocess.call(['matlab', '-nodesktop', '-nosplash', '-r', "loadDatabase; loadDatabaseB; loadDatabaseC; loadDatabaseD; load_Pitzer; Main_script_iterate; exit"])#,stdout=FNULL, stderr=subprocess.STDOUT)
    end = time.time()
    print('Done. Iteration time:',end-start)
    #save output
    subprocess.call(['mv','Gibbs_min_result.txt',direc+'Gibbs_result_eco.txt'])
    print('\n')
