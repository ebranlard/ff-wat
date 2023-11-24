import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatch
# Local 
import welib.weio as weio
from welib.essentials import *
import sys
sys.path.append('amr')
from helper_functions import *
from plothelper import *
# from t90_PostProMeanderingFrame import *

setFigureFont(14)
setFigurePath('../article_torque/figs/')

export=False
export=True


# --- Parameters
outDir='amr/_out_all/'
caseNames =[]
caseNames += ['stable2WT']     # NOT
caseNames += ['neutral2WT']    # NOTE: respect order!
# caseNames += ['unstable2WT']

symmetric=False
removeBG='0WT'
smooth=4
Meander=True
iTimeMin =500
oneRow=True

# Derived parametesr
D=240
for iCase, caseName in enumerate(caseNames):
    Case=AllCases[caseName]


    label = getLabel(caseName, Meander, symmetric, removeBG, smooth)
    outDirFits = os.path.join(outDir, 'kfit')
    prefix = '_Meander_' if Meander else '_Inertial_'
    # --- Derived parameters
    U0, dt, D, xyWT, HubHeight, xPlanes = getSimParamsAMR(Case['stability'])
    outDirHori = os.path.join(outDir, 'lines')

    # --- Load horizontal planes
    ds21 = readDataSet(os.path.join(outDirHori, caseName + prefix + 'WT1.nc_small')) # Wake behind turbine 1
    if oneRow:
        ds22 = readDataSet(os.path.join(outDirHori, caseName + prefix + 'WT2.nc_small')) # Wake behind turbine 2

    if removeBG=='0WT':
        caseName0 = caseName.replace('2WT','0WT').replace('1WT','0WT')
        if Meander:
            superPrefix = caseName0+ '_IN'+ caseName
        else:
            superPrefix = caseName0
        ds01 = readDataSet(os.path.join(outDirHori, superPrefix + prefix + 'WT1.nc_small')) # Wake behind turbine 1 - No turbine
        if oneRow:
            ds02 = readDataSet(os.path.join(outDirHori, superPrefix + prefix + 'WT2.nc_small')) # Wake behind turbine 2 - No turbine
    else:
        ds01=None
        ds02=None
    ds = ds21

    # --- Intersection of time
    # NOTE: not fully fair, but 0WT has shorter sim time for neutral
    if removeBG=='0WT':
        ITime = common_itime(ds01, ds21, ds02, ds22)
    else:
        ITime = common_itime(ds21, ds22)
    ITime = ITime[ITime>iTimeMin]
    print('ITimeSel', ITime[0], ITime[-1])

    # Storage
    KS = np.zeros((2, len(xPlanes), 2)) # iWT, xPlanes, (Kd, Kg)


    if oneRow:
        #fig,axes = plt.subplots(1, len(xPlanes), sharex=True, sharey=True, figsize=(12.8,3.8)) # (6.4,4.8)
        #fig.subplots_adjust(left=0.06, right=0.98, top=0.88, bottom=0.11, hspace=0.20, wspace=0.20)
        fig,axes = plt.subplots(1, len(xPlanes), sharex=True, sharey=True, figsize=(12.8,1.80)) # (6.4,4.8)
        fig.subplots_adjust(left=0.05, right=0.995, top=0.88, bottom=0.11, hspace=0.20, wspace=0.20)
        axes = np.atleast_2d(axes)
    else:
        fig,axes = plt.subplots(2, len(xPlanes), sharex=True, sharey=True, figsize=(12.8,3.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.06, right=0.98, top=0.881, bottom=0.11, hspace=0.20, wspace=0.20)

    for iP, x in enumerate(xPlanes):
        if iCase==0:
            axes[0,iP].set_title(r'$x={}D$'.format(int(x)), fontsize=12)
        # --- KFit 
        y, var, u, k_fit, ks_fit, fitter = ds_fitK(ds21, iP, ITime, D, U0, symmetric=symmetric, ds0=ds01, removeBG=removeBG, smooth=smooth)
        KS[0,iP, :] = ks_fit
        # plot
        axes[0,iP].plot(np.sqrt(var), y/D, c=colrsStab[iCase], label=r'$\sqrt{\sigma_u^2-\sigma_U^2 }$')
        axes[0,iP].plot(k_fit       , y/D,  'k--' ,label='k - Fit')
        if fitter is not None:
            #s = fitter.coeffsToString(sep='\n')
            s =  r'$k_{deg}=$' + r'${:.2f}$'.format(fitter.model['coeffs']['kd']) + '\n'
            s += r'$k_{grad}=$' +r'${:.2f}$'.format(fitter.model['coeffs']['kg'])
            axes[0,iP].text(0.32,-1.4, s, ha='left', va='center', fontsize=11)

    axes[0,0].set_ylabel('y/D [-]')
    if oneRow:
        axes[0,0].legend(fontsize=10, ncol=1, loc='upper right', borderpad=0)

    if not oneRow:
        for iP, x in enumerate(xPlanes):
            # --- KFit 
            y, var, u, k_fit, ks_fit, fitter = ds_fitK(ds22, iP, ITime, D, U0, symmetric=symmetric, ds0=ds02, removeBG=removeBG, smooth=smooth)
            KS[1,iP, :] = ks_fit
            # plot
            axes[1,iP].plot(np.sqrt(var), y/D,        label=r'$\sqrt{\sigma^2-o}$')
            axes[1,iP].plot(k_fit       , y/D, 'k--' ,label='fit')
            if fitter is not None:
                axes[1,iP].text(1,1, fitter.coeffsToString(sep='\n'), ha='center', va='center')

            # Compute some kind of K factors based on first row
            #print('---- iP', iP)
            kd =         KS[1,iP,0]    
            kg =         KS[1,iP,1]    
            axes[0,iP].set_xlim([0.0,1.9])
            axes[1,iP].set_xlim([0.0,1.9])
            axes[1,iP].set_xlabel('k [-]')
        axes[1,0].legend(fontsize=8, ncol=2, loc='upper center', borderpad=0)
        axes[1,0].set_ylabel('y/D [-]')

    for ax in axes.flatten():
        ax.tick_params(direction='in', top=True, right=True, labelright=False, labeltop=False, which='both')
        ax.set_ylim([-2,2])
        ax.set_xlim([ 0.0 ,1.2])


    fig.suptitle('WAT-FitExamples-{}'.format(caseName))

if export:
    export2pdf()

plt.show()
