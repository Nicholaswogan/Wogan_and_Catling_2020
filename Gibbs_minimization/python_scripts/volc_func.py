import numpy as np
import subprocess
import os


#retrieve Ve and Vp
def find_Vs(out_dir='PHOTOCHEM/OUTPUT/'):
    f = open(out_dir+'out.out','r')
    lines = f.readlines()
    for i in range(0,len(lines)):
        if lines[i]==' DISSOLVED CH4 AND H2 (ECO MODEL)\n':
            line = lines[i+3]

    VP = float(line[68:80])
    VE = float(line[92:104])
    RAT = float(line[128:139])
    f.close()
    return VP,VE

# function for changing H2 and CH4
def change_H(fH2,fCH4,in_dir='PHOTOCHEM/INPUTFILES/'):

    f = open(in_dir+'species.dat','r')
    lines = f.readlines()
    f.close()

    lines_new = []
    for line in lines:
        if line[0:5].replace(' ','')=='H2':
            TMP = list(('%.3E' % (fH2*1e-6)).replace('E+0','E+').replace('E-0','E-'))
            line = list(line)
            line[45:53] = TMP
            line = ''.join(line)
        if line[0:5].replace(' ','')=='CH4':
            TMP = list(('%.3E' % (fCH4*1e-6)).replace('E+0','E+').replace('E-0','E-'))
            line = list(line)
            line[45:53] = TMP
            line = ''.join(line)
        lines_new.append(line)

    f = open(in_dir+'species.dat','w')
    for line in lines_new:
        f.write(line)
    f.close()

def change_volc(FH2,FCO,FSO2,FH2S,in_dir='PHOTOCHEM/INPUTFILES/'):
    #import fluxes into species.dat
    f = open(in_dir+'species.dat','r')
    lines = f.readlines()
    f.close()

    lines_new = []
    for line in lines:
        if line[0:5].replace(' ','')=='H2':
            TMP = list(('%.3E' % FH2))
            line = list(line)
            line[54:63] = TMP
            line = ''.join(line)

        if line[0:5].replace(' ','')=='CO':
            TMP = list(('%.3E' % FCO))
            line = list(line)
            line[54:63] = TMP
            line = ''.join(line)

        if line[0:5].replace(' ','')=='SO2':
            TMP = list(('%.3E' % FSO2))
            line = list(line)
            line[54:63] = TMP
            line = ''.join(line)

        if line[0:5].replace(' ','')=='H2S':
            TMP = list(('%.3E' % FH2S))
            line = list(line)
            line[54:63] = TMP
            line = ''.join(line)
        lines_new.append(line)


    f = open(in_dir+'species.dat','w')
    for line in lines_new:
        f.write(line)
    f.close()

def change_LBOUND_H(bound,in_dir='PHOTOCHEM/INPUTFILES/'):
    f = open(in_dir+'species.dat','r')
    lines = f.readlines()
    f.close()

    lines_new = []
    for line in lines:
        if line[0:5].replace(' ','')=='H2':
            line = list(line)
            line[30] = str(int(bound))
            line = ''.join(line)
        if line[0:5].replace(' ','')=='CH4':
            line = list(line)
            line[30] = str(int(bound))
            line = ''.join(line)

        lines_new.append(line)

    f = open(in_dir+'species.dat','w')
    for line in lines_new:
        f.write(line)
    f.close()


def profile2dic(filename):
    file = open(filename)
    lines = file.readlines()

    key = lines[0].split()

    #build dictionary of output
    out = []
    for i in range(0,len(lines)):
        if i>0:
            tmp = []
            for j in lines[i].split():
                try:
                    tmp.append(float(j))
                except ValueError:
                    tmp.append(0)
            out.append(tmp)
    out = np.array(out)

    f_out = {}
    for i in range(0,len(key)):
        f_out[key[i]]=out[:,i]
    file.close()

    return f_out

def solve_eco(ft,out_dir = 'PHOTOCHEM/OUTPUT/'):
    FNULL = open(os.devnull, 'w')
    # Initial conditions
    fH2_1 = ft/5
    fCH4_1 = find_CH4(fH2_1,ft)
    print('Solving the Eco Model...')
    # Run photochem once to get new file
    change_H(fH2_1,fCH4_1)
    subprocess.call(['./Photo.run'],stdout=FNULL, stderr=subprocess.STDOUT)
    subprocess.call(['cp',out_dir+'out.dist','PHOTOCHEM/'],)
    subprocess.call(['mv','PHOTOCHEM/out.dist','PHOTOCHEM/in.dist'])

    # tolerance and iterations
    iter_count = 0
    tol = 0.05
    alpha = .5
    iters = 20

    ve = []
    vp = []
    fCH4 = []
    fH2 = []
    Error = []
    M = []
    B = []
    while iter_count<iters:
        iter_count+=1
        # change input file
        if iter_count>1:
            change_H(fH2_1,fCH4_1)
            subprocess.call(['./Photo.run'],stdout=FNULL, stderr=subprocess.STDOUT)
        # pull out results
        VP1,VE1 = find_Vs()
        y1 = VE1/VP1
        # calculate error
        err = abs(1-y1)

        # Save things
        Error.append(err)
        ve.append(VE1)
        vp.append(VP1)
        fCH4.append(fCH4_1)
        fH2.append(fH2_1)

        # break the loop if convergence is reached
        if err<tol:
            print('Convergence reached on iteration',iter_count-1)
            print('Ve/Vp =','%0.2e' % (VE1/VP1),'Error =','%0.2e' % err)
            print('CH4 =','%.2f' % fCH4_1,3*' ','H2 =','%.2f' % fH2_1)
            break

        # else calculate a new point
        else:
            print('Iteration',str(iter_count-1)+':','Ve/Vp =','%0.2e' % (VE1/VP1),'Error =','%0.2e' % err)
            print(12*' ','CH4 =','%.2f' % fCH4_1,3*' ','H2 =','%.2f' % fH2_1)
            # calculate another point
            fCH4_2 = fCH4_1+alpha
            fH2_2 = find_H2(fCH4_2,ft)
            change_H(fH2_2,fCH4_2)
            subprocess.call(['./Photo.run'],stdout=FNULL, stderr=subprocess.STDOUT)
            # pull out results
            VP2,VE2 = find_Vs()
            y2 = VE2/VP2

            #Calculate slope and intercept and then solve for new CH4
            m = (y2-y1)/(fCH4_2-fCH4_1)
            b = y2-m*fCH4_2
            fCH4_new = (1-b)/m
            fH2_new = find_H2(fCH4_new,ft)
            M.append(m)
            B.append(b)

            # update fCH4_1 and fH2_1
            fCH4_1 = fCH4_new
            fH2_1 = fH2_new

    if iter_count==iters:
        print('No Convergence!')
        print('Ve/Vp =','%0.2e' % (VE1/VP1),'Err =','%0.2e' % err)

def find_CH4(f_H2,ft):
    return (ft-f_H2)/2

def find_H2(f_CH4,ft):
    return ft-2*f_CH4
