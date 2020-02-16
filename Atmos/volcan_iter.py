import numpy as np
from matplotlib import pyplot as plt
import subprocess
import os
FNULL = open(os.devnull, 'w')

import sys
sys.path.append('./')
from volc_func import *

# The photochemical model needs a starting point called in.dist. Here we make
# sure the in.dist will allow the first model to converge.
subprocess.call(['rm','PHOTOCHEM/in.dist'])
subprocess.call(['cp','PHOTOCHEM/INPUTFILES/Disequilibrium/Start/in.dist','PHOTOCHEM/'])

# import fluxes!
(Q,FH2,FCO,FSO2,FH2S) = np.loadtxt('Data/volc_photochem_inputs.txt',skiprows=1).T

in_dir = 'PHOTOCHEM/INPUTFILES/'
out_dir = 'PHOTOCHEM/OUTPUT/'

for j in range(0,len(Q)):
    # First run prebiotic case:
    print('Q =',int(Q[j]))
    print('Solving the Prebiotic model...')


    change_volc(FH2[j],FCO[j],FSO2[j],FH2S[j])
    change_LBOUND_H(3,ecosystem='n')

    # Run photochem
    subprocess.call(['./Photo.run'],stdout=FNULL, stderr=subprocess.STDOUT)
    print('Completed')
    # Save the profile
    subprocess.call(['mkdir','Output/disequilibrium/Volc_iter/Q'+'%i' % Q[j]]) #open directory
    access = 'Output/disequilibrium/Volc_iter/Q'+'%i' % Q[j]+'/'
    subprocess.call(['cp',out_dir+'PTZ_mixingratios_out.dist',access])
    subprocess.call(['mv',access+'PTZ_mixingratios_out.dist',access+'PTZ_mixingratios_pre.dist'])
    subprocess.call(['cp',out_dir+'out.dist',access])
    subprocess.call(['mv',access+'out.dist',access+'out_pre.dist'])
    f_out = profile2dic(access+'PTZ_mixingratios_pre.dist')

    # Now run eco case:
    change_LBOUND_H(1,ecosystem='y')
    ft = (f_out['H2'][0]+2*f_out['CH4'][0])*1e6

    #change initial conditions to speed things up.
    if j>0:
        access1 = 'Output/disequilibrium/Volc_iter/Q'+'%i' % Q[j-1]+'/'
        subprocess.call(['cp',access1+'out_eco.dist','PHOTOCHEM/'])
        subprocess.call(['mv','PHOTOCHEM/out_eco.dist','PHOTOCHEM/in.dist'])
    solve_eco(ft)
    # save the profile
    subprocess.call(['cp',out_dir+'PTZ_mixingratios_out.dist',access])
    subprocess.call(['mv',access+'PTZ_mixingratios_out.dist',access+'PTZ_mixingratios_eco.dist'])
    subprocess.call(['cp',out_dir+'out.dist',access])
    subprocess.call(['mv',access+'out.dist',access+'out_eco.dist'])

    # prep for next iteration
    subprocess.call(['cp',access+'out_pre.dist','PHOTOCHEM/'])
    subprocess.call(['mv','PHOTOCHEM/out_pre.dist','PHOTOCHEM/in.dist'])
    print()
