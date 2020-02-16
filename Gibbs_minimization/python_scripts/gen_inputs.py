import numpy as np
import sys
sys.path.append('./')
from thermodynamics import *
from volc_func import *

def gen_inputs(input_file,output_file,CO2 = 0.2,N2=0.75):
    # This script generates a input file for the
    # Disequilibrium code from the output of a
    # photochemical model.

    # constants
    R = 8.314

    #Inputs
    # name of python dictionary
    # name of output file
    # Temperature
    pH = 6.6
    f = profile2dic(input_file)
    f['CO2'] = CO2*np.ones(len(f['H2']))
    f['N2'] = N2*np.ones(len(f['H2']))
    keys = []
    for key, value in f.items() :
        keys.append(key)
    #T = f1['TEMP'][0]
    T = 298
    P = 1
    fNH3 = 1e-10
    mSO4 = 1e-15
    mH2S = 1e-15
    prebiotic = "n"

    # Calculate Carbon chemistry equilibrium
    # Calculate equilibirum constants
    pK1 = 17.788-0.073104*T-0.0051087*35+1.1463e-4*T**2
    K1 = 10**-pK1
    pK2 = 20.919-0.064209*T-0.011887*35+8.7313e-5*T**2
    K2 = 10**-pK2
    mCO2 = f['CO2'][0]*henrys_coef('CO2',T,P)
    mH = 10**-pH
    mHCO3 = K1*mCO2/mH
    mCO3 = K2*mHCO3/mH

    # OH
    mOH = 10**-14/mH

    # NH3 NH4 chemsitry
    mNH3 = fNH3*henrys_coef('NH3',T,P)
    DGr = gibbsAQ('NH4(+)',T,P)+gibbsAQ('OH(-)',T,P)-(gibbsAQ('NH3',T,P)+(-237130))
    K = np.exp(-DGr/(R*T))
    mNH4 = K*mNH3/mOH



    file = open('inputs_Archean_max.txt','r')

    lines = file.readlines()
    file.close()
    ident = lines[1].strip('\n').split()
    molc = lines[2].strip('\n').split(',')
    concen = lines[3].strip('\n').split()


    for i in range(0,len(molc)):
        for key in keys:
            # gases
            if molc[i]=="'"+key+"  '" and ident[i]=='1':
                concen[i] = '%0.3e' % (f[key][0])
            # henry law
            if molc[i]=="'"+key+"'" and ident[i]=='4':
                concen[i] = '%0.3e' % (f[key][0]*henrys_coef(key,T,P))

            if molc[i]=="'"+key+"(0)'" and ident[i]=='4':
                concen[i] = '%0.3e' % (f[key][0]*henrys_coef(key,T,P))

        # set pH and carbon chemsitries
        if molc[i]=="'H(+)'":
            concen[i] = '%0.3e' % (mH)
        if molc[i]=="'OH(-)'":
            concen[i] = '%0.3e' % (mOH)
        if molc[i]=="'HCO3(-)'":
            concen[i] = '%0.3e' % (mHCO3)
        if molc[i]=="'CO3(-2)'":
            concen[i] = '%0.3e' % (mCO3)
        # Nitrogen stuff
        if molc[i]=="'NH3  '" and ident[i]=='1':
            concen[i] = '%0.3e' % (fNH3)
        if molc[i]=="'NH3'" and ident[i]=='4':
            concen[i] = '%0.3e' % (mNH3)
        if molc[i]=="'NH4(+)'":
            concen[i] = '%0.3e' % (mNH4)
        # SO4
        if molc[i]=="'SO4(-2)'" and ident[i]=='4':
            concen[i] = '%0.3e' % (mSO4)
        if molc[i]=="'H2S'" and ident[i]=='4':
            concen[i] = '%0.3e' % (mH2S)
        if molc[i]=="'H2O  '" and ident[i]=='1':
            concen[i] = '%0.3e' % (2.5e-2)

    #Write output file
    file = open(output_file,'w')
    file.write(lines[0])
    file.write(lines[1])
    file.write(lines[2])
    file.write(' '.join(concen)+'\n')
    file.write(lines[4])
    file.close()

def gen_inputs_atmos_only(input_file,output_file,CO2 = 0.2,N2=0.75):
    # This script generates a input file for the
    # Disequilibrium code from the output of a
    # photochemical model.

    # constants
    R = 8.314

    #Inputs
    # name of python dictionary
    # name of output file
    # Temperature
    pH = 6.6
    f = profile2dic(input_file)
    f['CO2'] = CO2*np.ones(len(f['H2']))
    f['N2'] = N2*np.ones(len(f['H2']))
    fNH3 = 1e-10
    f['NH3'] = fNH3*np.ones(len(f['H2']))
    keys = []
    for key, value in f.items() :
        keys.append(key)
    #T = f1['TEMP'][0]
    T = 298
    P = 1


    file = open('atmosphere_only.txt','r')

    lines = file.readlines()
    file.close()
    ident = lines[1].strip('\n').split()
    molc = lines[2].strip('\n').split(',')
    concen = lines[3].strip('\n').split()

    for i in range(0,len(molc)):
        for key in keys:
            # gases
            if molc[i]=="'"+key+"  '" and ident[i]=='1':
                concen[i] = '%0.3e' % (f[key][0])
        if molc[i]=="'H2O  '" and ident[i]=='1':
            concen[i] = '%0.3e' % (2.5e-2)

    #Write output file
    file = open(output_file,'w')
    file.write(lines[0])
    file.write(lines[1])
    file.write(lines[2])
    file.write(' '.join(concen)+'\n')
    file.write(lines[4])
    file.close()
