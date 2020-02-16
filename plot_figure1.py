import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib import cm

import sys
sys.path.append('./Atmos/')
from volc_func import *

QQ = np.arange(1,26,2)
folder = 'Gibbs_minimization/Volc_iter/'

DG_pre = []
DG_eco = []
DG_pre_a = []
DG_eco_a = []


H2_pre = []
CH4_pre = []
CO_pre = []
O2_pre = []
H2_eco = []
CH4_eco = []
CO_eco = []
O2_eco = []
for Q in QQ:
    file = open(folder+'Q'+str(Q)+'/'+'Gibbs_result_pre.txt')
    lines = file.readlines()
    DG_pre.append(float(lines[1].strip('\n')))
    H2_pre.append(float(lines[11].strip('\n').split()[0]))
    CH4_pre.append(float(lines[9].strip('\n').split()[0]))
    CO_pre.append(float(lines[10].strip('\n').split()[0]))
    O2_pre.append(float(lines[4].strip('\n').split()[0]))


    file.close()

    file = open(folder+'Q'+str(Q)+'/'+'Gibbs_result_pre_atmos_only.txt')
    lines = file.readlines()
    DG_pre_a.append(float(lines[1].strip('\n')))
    file.close()

    file = open(folder+'Q'+str(Q)+'/'+'Gibbs_result_eco.txt')
    lines = file.readlines()
    DG_eco.append(float(lines[1].strip('\n')))
    H2_eco.append(float(lines[11].strip('\n').split()[0]))
    CH4_eco.append(float(lines[9].strip('\n').split()[0]))
    CO_eco.append(float(lines[10].strip('\n').split()[0]))
    O2_eco.append(float(lines[4].strip('\n').split()[0]))
    file.close()

    file = open(folder+'Q'+str(Q)+'/'+'Gibbs_result_eco_atmos_only.txt')
    lines = file.readlines()
    DG_eco_a.append(float(lines[1].strip('\n')))
    file.close()

DG_pre = np.array(DG_pre)
DG_eco = np.array(DG_eco)
DG_pre_a = np.array(DG_pre_a)
DG_eco_a = np.array(DG_eco_a)

scale = 3
fnt = 25*scale
fnt1 = 20*scale

start = 0.0
stop = 0.55
number_of_lines= 2
cm_subsection = np.linspace(start, stop, number_of_lines)
colors = [ cm.inferno(x) for x in cm_subsection ]
mpl.rcParams['axes.linewidth'] = 1*scale

fig,ax = plt.subplots(1,1,figsize=[16*scale,10*scale])

color1 = colors[0]
color2 = colors[1]
color3 = 'k'

ax.plot(QQ,H2_pre,'-',color=color1,linewidth=2.5*scale,label = r'H$_2$ prebiotic')
ax.plot(QQ,CH4_pre,'--',color=color1,linewidth=2.5*scale,label = r'CH$_4$ prebiotic')
ax.plot(QQ,CO_pre,'-.',color=color1,linewidth=2.5*scale,label = r'CO prebiotic')
# ax.plot(QQ,O2_pre,color=colors[2],linewidth=2.5,label = r'O$_2$ prebiotic')

ax.plot(QQ,H2_eco,'-',color=color2,linewidth=2.5*scale,label = r'H$_2$ chemotrophic')
ax.plot(QQ,CH4_eco,'--',color=color2,linewidth=2.5*scale,label = r'CH$_4$ chemotrophic')
ax.plot(QQ,CO_eco,'-.',color=color2,linewidth=2.5*scale,label = r'CO chemotrophic')
# ax.plot(QQ,O2_eco,'--',color=colors[2],linewidth=2.5,label = r'O$_2$ chemosynthetic')

# AA = [5,5]
# AA1 = [1e-3,1e-3]
# ax.plot(AA,AA1,'k',label=r'Prebiotic')
# ax.plot(AA,AA1,'k--',label=r'Chemosynthetic')

ax.set_ylabel('Mixing ratio',fontsize=fnt)
ax.set_xlabel('Volcanic outgassing multiplier, $C$',fontsize=fnt)
ax.tick_params(axis='both', which='major', labelsize=fnt)
ax.tick_params(axis='both', which='minor', labelsize=fnt)

ax.set_xticks(np.arange(1,26,4))

ax.legend(fontsize=fnt1,loc=4,ncol=2)
# ax.set_ylabel('Avaliable Gibbs Energy (J/mol)')
# ax.set_xlabel('Heat flux relative to modern')
ax.set_yscale('log')
ax.set_yticks(10**np.arange(-9,-.9,1))
ax.xaxis.set_tick_params(which='major',width=1*scale,length=4*scale)
ax.xaxis.set_tick_params(which='minor',width=.8*scale,length=2*scale)
ax.yaxis.set_tick_params(which='major',width=1*scale,length=4*scale)
ax.yaxis.set_tick_params(which='minor',width=.8*scale,length=2*scale)

ax.tick_params(axis='x', which='minor')

locmin = mpl.ticker.LogLocator(base=10.0,subs=np.arange(.1,1,.1),numticks=12)
ax.yaxis.set_minor_locator(locmin)
ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())

#plt.savefig("surface_mixing_rat.jpg",bbox_inches='tight')
plt.savefig("surface_mixing_rat.pdf",bbox_inches='tight')
