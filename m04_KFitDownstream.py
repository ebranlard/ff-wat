import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import welib.weio as weio
from welib.essentials import *
import sys
sys.path.append('amr')
from helper_functions import *
from plothelper import *
from t12_PlotK import *

setFigureFont(14)


export=False
export=True


outDir = 'amr/_out_all'
figDir = '_out_all'
Meander=True
syms = [True, False]
rms =['0WT', 'Outer']
#     rms =['0WT']

caseNames = [] 
caseNames += ['neutral2WT']
caseNames += ['stable2WT'] 
caseNames += ['unstable2WT'] 
#     caseNames += ['neutral1WT']
#     caseNames += ['stable1WT'] 
#     caseNames += ['unstable1WT'] 

Cases = dict((k,v) for k,v in AllCases.items() if k in caseNames)
# --- Load all
D={}
# for smooth in [0,2,4]:
for smooth in [2,4]:
    for nWT in [2]:
        D0 = loadallKFits(Meander, syms, Cases, rms, smooth=smooth, nWT=nWT, outDir=outDir)
        D.update(D0)
# --- Average parametric studies
D_avg={}
for stab in stabilities:
    keys = [k for k in D.keys() if k.find(stab)==5]
    n=len(keys)
    D_avg[stab] =np.zeros_like(D[keys[0]])
    for k in keys:
        D_avg[stab] +=D[k]/n


# --- Super KD average
xPlanes=np.arange(7) # TODO for 1WT..
KD_avg=np.zeros_like(D_avg[stabilities[0]][0,:,0]   )
for stab in stabilities:
    KD_avg += D_avg[stab][0,:,0]/len(stabilities)
print('KD_avg:',KD_avg)
print('KD_avg:',np.mean(KD_avg))
xPlanes2 = np.linspace(0,25,100)
kd_model = 0.74 + xPlanes2*0
kd_model= eddyFilter(xPlanes2, Dmin=0, Dmax=2, Fmin=0.00, Exp=1.00)* 0.74

CKF=2.5
kf_model1= eddyFilter(xPlanes2, **eddyFilter_Shr)* CKF
kf_model2= eddyFilter(xPlanes2, Dmin=0, Dmax=12, Fmin=0.00, Exp=0.65)*CKF





# --- Plot
fig,axes = plt.subplots(1, 2, sharex=True, figsize=(12.8,3.5)) # (6.4,4.8)
fig.subplots_adjust(left=0.07, right=0.985, top=0.97, bottom=0.140, hspace=0.20, wspace=0.265)
for iStab, stab in enumerate(stabilities):
    KS = D_avg[stab]
    iRow = 0
    kd = KS[iRow, :, 0] 
    kg = KS[iRow, :, 1] 

    keys = [k for k in D.keys() if k.find(stab)==5]
    for k in keys:
        axes[0].plot(xPlanes, D[k][iRow,:,0], c=colrsStab[iStab], alpha=0.2, lw=0.8)
        axes[1].plot(xPlanes, D[k][iRow,:,1], c=colrsStab[iStab], alpha=0.2, lw=0.8)


    axes[0].plot(xPlanes, kd, c=colrsStab[iStab], label=stab.capitalize(), lw=1.8)
    axes[1].plot(xPlanes, kg, c=colrsStab[iStab], label=stab.capitalize(), lw=1.8)


axes[0].plot(xPlanes2, kd_model, ls='--', c='k', label='Model', lw=1.8)
axes[1].plot(xPlanes2, kf_model1, ls=':', c='k', label='Model (shr)', lw=1.8)
axes[1].plot(xPlanes2, kf_model2, ls='--', c='k', label='Model2', lw=1.8)


axes[0].set_ylabel(r'$k_{def}$ [-]')
axes[1].set_ylabel(r'$k_{grad}$ [-]')

for ax in np.array(axes).flatten():
    ax.set_ylim([0,4.0])
    #ax.set_xlim([np.min(xPlanes),np.max(xPlanes)])
    ax.set_xlim([0,14])
    ax.set_xlabel('x/D [-]')
    ax.tick_params(direction='in', top=True, right=True, labelright=False, labeltop=False, which='both')

axes[0].legend(fontsize=13, ncol=1, loc='upper right')
fig.suptitle('KFitDownstream')

if export:
    export2pdf()

plt.show()
