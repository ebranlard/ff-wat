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
import weio
from weio.pickle_file import PickleFile
from welib.essentials import *
from welib.tools.curve_fitting import model_fit
from helper_functions import *
import warnings
warnings.simplefilter("ignore", category=RuntimeWarning)


def fitAvgProfiles(case, Meander=True, symmetric=False, removeBG='Outer', smooth=0):
    figsPath='_figs_KFit/'
    outPath='_out_KFit/'

    # --- Derived parameters
    if Meander: 
        label='KFit_Meander_'+case
    else:
        label='KFit_Inertial_'+case
    if symmetric:
        label+='_sym'
    label+='_rm'+removeBG+'_smooth{:d}'.format(smooth)
    print('>>>>>>>>>>>>> ',label)

    U0, dt, D, xyWT1, xyWT2, xPlanes = getSimParamsAMR(case)
    datapath0 = os.path.join('02-iea15-0WT/', case, 'processedData')
    datapath2 = os.path.join('02-iea15-2WT/', case, 'processedData')

    # --- Load horizontal planes
    if Meander:
        ds21= xarray.open_dataset(os.path.join(datapath2,'MeanderWT1.nc_small')) # Wake behind turbine 1
        ds22= xarray.open_dataset(os.path.join(datapath2,'MeanderWT2.nc_small')) # Wake behind turbine 2
        try:
            ds01= xarray.open_dataset(os.path.join(datapath0,'MeanderWT1.nc_small')) # Wake behind turbine 1 - No turbine
            ds02= xarray.open_dataset(os.path.join(datapath0,'MeanderWT2.nc_small')) # Wake behind turbine 2 - No turbine
        except:
            ds01=ds21.copy()
            ds02=ds22.copy()
    else:
        ds21= xarray.open_dataset(os.path.join(datapath2,'HubHeightWT1.nc_small')) # Wake behind turbine 1
        ds22= xarray.open_dataset(os.path.join(datapath2,'HubHeightWT2.nc_small')) # Wake behind turbine 2
        try:
            ds01= xarray.open_dataset(os.path.join(datapath0,'HubHeightWT1.nc_small')) # Wake behind turbine 1 - No turbine
            ds02= xarray.open_dataset(os.path.join(datapath0,'HubHeightWT2.nc_small')) # Wake behind turbine 2 - No turbine
        except:
            ds01=ds21.copy()
            ds02=ds22.copy()
    ds = ds21

    # --- Intersection of time
    # NOTE: not fully fair, but 0WT has shorter sim time for neutral
    if removeBG=='0WT':
        ITime = common_itime(ds01, ds21, ds02, ds22)
    else:
        ITime = common_itime(ds21, ds22)
    ITime = np.arange(501, ITime[-1])
    print('ITime', ITime[0], ITime[-1])

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
        print('---- iP', iP)
        kd =         KS[1,iP,0]    
        kg =         KS[1,iP,1]    
        print('kd={:.3f} kg={:.3f}'.format(kd, kg))
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
    fig.savefig(figsPath+figname+'.png')


    # --- Export
    pkl = PickleFile(data=KS)
    pkl.write(outPath + label+'.pkl')


if __name__ == '__main__':
    # --- Parameters 
    cases=[]
    cases+=['neutral']
    cases+=['stable']
    cases+=['unstable']


    for smooth in [0,2,4]:
    #     for Meander in [False]:
        for Meander in [True, False]:
            for symmetric in [True, False]:
    #         for symmetric in [False]:
            #for symmetric in [True]:
                for case in cases:
                    if case!='unstable':
                        BG=['Outer','0WT']
                    else:
                        BG=['Outer']
    #                 BG=['Outer']
                    for bg in BG:
                        fitAvgProfiles(case, Meander=Meander, symmetric=symmetric, removeBG=bg, smooth=smooth)


#     plt.show()
