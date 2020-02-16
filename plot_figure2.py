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
#DG_pre_a = []
#DG_eco_a = []


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

    # file = open(folder+'Q'+str(Q)+'/'+'Gibbs_result_pre_atmos_only.txt')
    # lines = file.readlines()
    # DG_pre_a.append(float(lines[1].strip('\n')))
    # file.close()

    file = open(folder+'Q'+str(Q)+'/'+'Gibbs_result_eco.txt')
    lines = file.readlines()
    DG_eco.append(float(lines[1].strip('\n')))
    H2_eco.append(float(lines[11].strip('\n').split()[0]))
    CH4_eco.append(float(lines[9].strip('\n').split()[0]))
    CO_eco.append(float(lines[10].strip('\n').split()[0]))
    O2_eco.append(float(lines[4].strip('\n').split()[0]))
    file.close()

    # file = open(folder+'Q'+str(Q)+'/'+'Gibbs_result_eco_atmos_only.txt')
    # lines = file.readlines()
    # DG_eco_a.append(float(lines[1].strip('\n')))
    # file.close()

DG_pre = np.array(DG_pre)
DG_eco = np.array(DG_eco)
# DG_pre_a = np.array(DG_pre_a)
# DG_eco_a = np.array(DG_eco_a)

scale = 1
mpl.rcParams['axes.linewidth'] = 1*scale
fig,ax = plt.subplots(1,1,figsize=[16*scale,10*scale])




fnt = 25*scale
fnt1 = 20*scale

start = 0.0
stop = .55
number_of_lines= 2
cm_subsection = np.linspace(start, stop, number_of_lines)
colors = [ cm.inferno(x) for x in cm_subsection ]

color1 = colors[0]
color2 = colors[1]

# plot the raw results
#ax.scatter(QQ,-DG_pre,c='b',marker='s',label=r'prebiotic')
#ax.scatter(QQ,-DG_eco,c='r',marker='s',label=r'chemosynthetic')

ax.plot(QQ,-DG_pre,linewidth=2.5*scale,c=color1,marker='o',markersize=7*scale)
ax.plot(QQ,-DG_eco,linewidth=2.5*scale,c=color2,marker='o',markersize=7*scale)

# ax.plot(QQ,-DG_pre_a,'--',linewidth=2.5*scale,c=color1,marker='o',markersize=7*scale)
# ax.plot(QQ,-DG_eco_a,'--',linewidth=2.5*scale,c=color2,marker='o',markersize=7*scale)

AA = [10,10]
AA1 = [.01,.01]

#ax.plot(AA,AA1,'-',color = color1,marker='o',label=r'Atmosphere-Ocean')
#ax.plot(AA,AA1,'--',color = color2,label=r'Atmosphere-only')

ax.text(7,3000,'Prebiotic',fontsize=fnt,color=color1)
ax.text(10,340,'Chemotrophic',fontsize=fnt,color=color2)

ax.text(20,350,'Atmosphere-Ocean',fontsize=fnt1,color=color2)
# ax.text(20,19,'Atmosphere-only',fontsize=fnt1,color=color2)
ax.text(12,3900,'Atmosphere-Ocean',fontsize=fnt1,color=color1)
# ax.text(12,1300,'Atmosphere-only',fontsize=fnt1,color=color1)

# plot the robust fit instead
#ax.plot(QQ,DG_pre_fit,'b',label=r'prebiotic')
#ax.plot(QQ,DG_eco_fit,'r',label=r'chemosynthetic')
#ax.plot(QQ,-DG_mod,label=r'$\Delta G_{modern}$')

#leg = ax.legend(loc=4,fontsize=fnt)
ax.set_ylabel('Avaliable Gibbs Energy (J/mol)',fontsize=fnt)
ax.set_xlabel('Volcanic outgassing multiplier, $C$',fontsize=fnt)
ax.set_yscale('log')
ax.tick_params(axis='both', which='major', labelsize=fnt)
ax.tick_params(axis='both', which='minor', labelsize=fnt)
ax.set_xticks(np.arange(1,26,4))
#ax.grid()


ax.xaxis.set_tick_params(width=1*scale,length=4*scale)
ax.xaxis.set_tick_params(which='minor',width=.8*scale,length=2*scale)
ax.yaxis.set_tick_params(width=1*scale,length=4*scale)
ax.yaxis.set_tick_params(which='minor',width=.8*scale,length=2*scale)

ax.set_ylim(1,8000)
#ax.set_title(r'Modern $\Delta$G = 2326 J/mol')

#plt.savefig("volc_iter.jpg",bbox_inches='tight')
plt.savefig("volc_iter.pdf",bbox_inches='tight',format='pdf')
