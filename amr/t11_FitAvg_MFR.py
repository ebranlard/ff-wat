""" 
Fit wake average profiles in the non meandering frame of reference.
Just to get familiar with nonlinear fitting and gradients
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray
# Local 
from welib.weio.pickle_file import PickleFile
from welib.essentials import *
from welib.tools.curve_fitting import model_fit
from helper_functions import *
import warnings
warnings.simplefilter("ignore", category=RuntimeWarning)


def fitAvgProfiles(caseName, Case, Meander=True, symmetric=False, removeBG='Outer', smooth=0, outDir='_out', figDir='_figs', iTimeMin=500):
    label = getLabel(caseName, Meander, symmetric, removeBG, smooth)
    print(f'-----------------------------------------------------------------------')
    print(f'--- fitAvgProfiles - Case: {label}')
    print(f'-----------------------------------------------------------------------')
    if Case['nWT']==0:
        WARN('Skipping case (script only makes sense for one or two WT)')
        return

    outDirFits = os.path.join(outDir, 'kfit')
    figDir = os.path.join(figDir, '_figs_kfit')
    if not os.path.exists(figDir):
        os.makedirs(figDir)
    if not os.path.exists(outDirFits):
        os.makedirs(outDirFits)

    prefix = '_Meander_' if Meander else '_Inertial_'

    # --- Derived parameters
    U0, dt, D, xyWT, HubHeight, xPlanes = getSimParamsAMR(Case['stability'])
    outDirHori = os.path.join(outDir, 'lines')

    # --- Load horizontal planes
    ds21 = readDataSet(os.path.join(outDirHori, caseName + prefix + 'WT1.nc_small')) # Wake behind turbine 1
    ds22 = readDataSet(os.path.join(outDirHori, caseName + prefix + 'WT2.nc_small')) # Wake behind turbine 2

    if removeBG=='0WT':
        caseName0 = caseName.replace('2WT','0WT').replace('1WT','0WT')
        if Meander:
            superPrefix = caseName0+ '_IN'+ caseName
        else:
            superPrefix = caseName0
        ds01 = readDataSet(os.path.join(outDirHori, superPrefix + prefix + 'WT1.nc_small')) # Wake behind turbine 1 - No turbine
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

    fig,axes = plt.subplots(2, len(xPlanes), sharex=True, sharey=True, figsize=(12.8,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.06, right=0.98, top=0.88, bottom=0.11, hspace=0.20, wspace=0.20)
    for iP, x in enumerate(xPlanes):
        axes[0,iP].set_title('x= {}D'.format(int(x)))
        # --- KFit 
        y, var, u, k_fit, ks_fit, fitter = ds_fitK(ds21, iP, ITime, D, U0, symmetric=symmetric, ds0=ds01, removeBG=removeBG, smooth=smooth)
        KS[0,iP, :] = ks_fit
        # plot
        axes[0,iP].plot(np.sqrt(var), y/D,        label='sqrt(Var)')
        axes[0,iP].plot(k_fit       , y/D,  'k--' ,label='fit')
        if fitter is not None:
            axes[0,iP].text(1,1, fitter.coeffsToString(sep='\n'), ha='center', va='center')

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
        #print('kd={:.3f} kg={:.3f}'.format(kd, kg))
#         kd = np.sqrt(KS[0,iP,0]**2 + KS[0,-1,0]**2)
#         kg = np.sqrt(KS[0,iP,1]**2 + KS[0,-1,1]**2)
#         print('kd={:.3f} kg={:.3f}'.format(kd, kg))
#         K, KD, KG, du, gu = kWAT1D(y, U0, u, D, kd=kd, kg=kg)
#         axes[1,iP].plot(K           , y/D,  ':',     label='superp')
        axes[0,iP].set_xlim([0.0,1.9])
        axes[1,iP].set_xlim([0.0,1.9])
        axes[1,iP].set_xlabel('k [-]')
    axes[1,0].legend(fontsize=8, ncol=2, loc='upper center', borderpad=0)
    axes[0,0].set_ylabel('y/D [-]')
    axes[1,0].set_ylabel('y/D [-]')
    figname = label
    fig.suptitle(figname)
    fig.savefig(os.path.join(figDir,figname+'.png'))


    # --- Export
    outFile = os.path.join(outDirFits, label+'.pkl')
    print('Writting: ',outFile)
    pkl = PickleFile(data=KS)
    pkl.write(outFile)


if __name__ == '__main__':
    caseNames = [] 
    caseNames += ['neutral2WT']
#     caseNames += ['stable2WT'] 
#     caseNames += ['unstable2WT'] 
#     caseNames += ['neutral1WT']
#     caseNames += ['stable1WT'] 
#     caseNames += ['unstable1WT'] 

    Cases = dict((k,v) for k,v in AllCases.items() if k in caseNames)

    BG=['Outer','0WT']

    for caseName, Case in Cases.items():
        for smooth in [0,2,4]:
            for Meander in [True, False]:
                for symmetric in [True, False]:
                    for bg in BG:
                        label = getLabel(caseName, Meander, symmetric, bg, smooth)
                        with FailSafe(label, False):
                            fitAvgProfiles(caseName, Case, Meander=Meander, symmetric=symmetric, removeBG=bg, smooth=smooth)
