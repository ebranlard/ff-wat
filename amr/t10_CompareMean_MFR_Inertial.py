""" 
Plot mean horizontal velocity profile at each downstream plane between WT1&2
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray
# Local 
import weio
from welib.essentials import *
from helper_functions import *
import warnings
warnings.simplefilter("ignore", category=RuntimeWarning)

# --- Script parameters
figsPath='_figs_MFR_vs_Interial/'

# for nWT in [0,2]:
# for case in ['neutral']:
# for case in ['neutral']:
# for case in ['stable']:
# for case in ['unstable']:
for case in ['neutral','stable','unstable']:

    U0, dt, D, xyWT, HubHeight, xPlanes = getSimParamsAMR(case)
    # --- Derived parameters
    datapath = os.path.join('02-iea15-{}WT/'.format(2), case, 'processedData')
    figbase = 'MFR_vs_Inertial_Hori_{}'.format(case)

    ds1MFR = xarray.open_dataset(os.path.join(datapath,'MeanderWT1.nc_small'))
    ds2MFR = xarray.open_dataset(os.path.join(datapath,'MeanderWT2.nc_small'))
    ds1Ine = xarray.open_dataset(os.path.join(datapath,'HubHeightWT1.nc_small'))
    ds2Ine = xarray.open_dataset(os.path.join(datapath,'HubHeightWT2.nc_small'))
    datapath = os.path.join('02-iea15-{}WT/'.format(0), case, 'processedData')
    if case !='unstable':
        ds3MFR = xarray.open_dataset(os.path.join(datapath,'MeanderWT1.nc_small'))
        ds4MFR = xarray.open_dataset(os.path.join(datapath,'MeanderWT2.nc_small'))
        ds3Ine = xarray.open_dataset(os.path.join(datapath,'HubHeightWT1.nc_small'))
        ds4Ine = xarray.open_dataset(os.path.join(datapath,'HubHeightWT2.nc_small'))
    else:
        ds3MFR = ds1MFR.copy()
        ds4MFR = ds1MFR.copy()
        ds3Ine = ds1MFR.copy()
        ds4Ine = ds1MFR.copy()
        ds3MFR['u'] *= 0
        ds4MFR['u'] *= 0
        ds3Ine['u'] *= 0
        ds4Ine['u'] *= 0
    print('it',ds1MFR.it.values[0],ds1MFR.it.values[-1]  )

    y = ds1MFR.y.values
    t = ds1MFR.it.values * dt
    ITime = common_itime(ds1MFR,)
    print('ITime 2WT',ITime[0], ITime[-1])
    if case !='unstable':
        ITime = common_itime(ds3MFR,)
        print('ITime 0WT',ITime[0], ITime[-1])
        ITime = common_itime(ds1MFR, ds3MFR)
        print('ITime com',ITime[0], ITime[-1])
    ITime = np.arange(501, ITime[-1])
    print('ITime com',ITime[0], ITime[-1])


    Colrs =[fColrs(1), fColrs(2), fColrs(1), fColrs(2)]

    # --- Plot mean at hub height as function of time
    ysym=False
    fig,axes = plt.subplots(2, len(xPlanes), sharey=True, figsize=(12.8,8.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.91, bottom=0.07, hspace=0.15, wspace=0.20)
    for iP, x in enumerate(xPlanes):
        u1MFR = compute_u(ds1MFR, iP, ITime, ysym=ysym)
        u1Ine = compute_u(ds1Ine, iP, ITime, ysym=ysym)
        u3MFR = compute_u(ds3MFR, iP, ITime, ysym=ysym)
        u3Ine = compute_u(ds3Ine, iP, ITime, ysym=ysym)
        u2MFR = compute_u(ds2MFR, iP, ITime, ysym=ysym)
        u2Ine = compute_u(ds2Ine, iP, ITime, ysym=ysym)
        u4MFR = compute_u(ds4MFR, iP, ITime, ysym=ysym)
        u4Ine = compute_u(ds4Ine, iP, ITime, ysym=ysym)

        axes[0,iP].plot(u1MFR, y/D, c=Colrs[0], lw=2.0,           label='MFR'        )
        axes[0,iP].plot(u1Ine, y/D, c=Colrs[1], lw=2.0,           label='Inertial'   )
        axes[0,iP].plot(u3MFR, y/D, c=Colrs[2], lw=1.0, ls='--',  label='BG MFR'     )
        axes[0,iP].plot(u3Ine, y/D, c=Colrs[3], lw=1.0, ls='--',  label='BG Inertial')
        axes[1,iP].plot(u2MFR, y/D, c=Colrs[0], lw=2.0,            label='MFR')
        axes[1,iP].plot(u2Ine, y/D, c=Colrs[1], lw=2.0,            label='Inertial')
        axes[1,iP].plot(u4MFR, y/D, c=Colrs[2], lw=1.0, ls='--',   label='BG MFR')
        axes[1,iP].plot(u4Ine, y/D, c=Colrs[3], lw=1.0, ls='--',   label='BG Inertial')

        u_fitMFR,  sigMFR, fitterMFR =  fitGaussian01(y, u1MFR.values, U0, D)
        u_fitIne,  sigIne, fitterIne =  fitGaussian01(y, u1Ine.values, U0, D)
        axes[0,iP].text(5,1.2, r'$\sigma={}$'.format(np.around(sigMFR,1)),  c=Colrs[0], ha='center')
        axes[0,iP].text(5,1  , r'$\sigma={}$'.format(np.around(sigIne,1)),  c=Colrs[1], ha='center')
        #axes[1,iP].plot(u_fitMFR, y/D, c=Colrs[0], lw=0.5, marker='.',          label='MFR'        )
        #axes[1,iP].plot(u_fitIne, y/D, c=Colrs[1], lw=0.5, marker='.',          label='Inertial'   )

        u_fitMFR,  sigMFR, fitterMFR =  fitGaussian01(y, u2MFR.values, U0, D)
        u_fitIne,  sigIne, fitterIne =  fitGaussian01(y, u2Ine.values, U0, D)
        axes[1,iP].text(5,1.2, r'$\sigma={}$'.format(np.around(sigMFR,1)),  c=Colrs[0], ha='center')
        axes[1,iP].text(5,1  , r'$\sigma={}$'.format(np.around(sigIne,1)),  c=Colrs[1], ha='center')
        #axes[1,iP].plot(u_fitMFR, y/D, c=Colrs[0], lw=0.5, marker='.',          label='MFR'        )
        #axes[1,iP].plot(u_fitIne, y/D, c=Colrs[1], lw=0.5, marker='.',          label='Inertial'   )

        axes[0,iP].set_title('x= {}D'.format(int(x)))
        axes[1,iP].set_xlabel('y/D [-]')
        if iP==0:
            axes[0,0].legend(fontsize=8, ncol=2, loc='upper center')

    axes[0,0].set_ylabel('Time-averaged velocity\n horizontal plane [m/s]')
    axes[1,0].set_ylabel('Time-averaged velocity\n horizontal plane [m/s]')
    for ax in np.array(axes).flatten():
        ax.tick_params(direction='in', top=True, right=True, labelright=False, labeltop=False, which='both')
        ax.set_xlim([2,10])
    figname = figbase+ '__U'
    fig.suptitle(figname)
    fig.savefig(figsPath+figname+'.png')


    # --- Plot var at hub height as function of time
    ysym=False
    removeBG='Zero'
#     removeBG='Outer'
    fig,axes = plt.subplots(2, len(xPlanes), sharey=True, figsize=(12.8,8.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.91, bottom=0.07, hspace=0.15, wspace=0.20)
    for iP, x in enumerate(xPlanes):

        axes[0,iP].set_title('x= {}D'.format(int(x)))
        axes[0,iP].plot(np.sqrt(compute_var(ds1MFR, iP, ITime, D=D, removeBG=removeBG, ysym=ysym)), y/D, c=Colrs[0], lw=2.0,            label='MFR')
        axes[0,iP].plot(np.sqrt(compute_var(ds1Ine, iP, ITime, D=D, removeBG=removeBG, ysym=ysym)), y/D, c=Colrs[1], lw=2.0,            label='Inertial')
        axes[0,iP].plot(np.sqrt(compute_var(ds3MFR, iP, ITime, D=D, removeBG=removeBG, ysym=ysym)), y/D, c=Colrs[2], lw=1.0, ls='--',   label='BG MFR')
        axes[0,iP].plot(np.sqrt(compute_var(ds3Ine, iP, ITime, D=D, removeBG=removeBG, ysym=ysym)), y/D, c=Colrs[3], lw=1.0, ls='--',   label='BG Inertial')
        axes[1,iP].plot(np.sqrt(compute_var(ds2MFR, iP, ITime, D=D, removeBG=removeBG, ysym=ysym)), y/D, c=Colrs[0], lw=2.0,            label='MFR')
        axes[1,iP].plot(np.sqrt(compute_var(ds2Ine, iP, ITime, D=D, removeBG=removeBG, ysym=ysym)), y/D, c=Colrs[1], lw=2.0,            label='Inertial')
        axes[1,iP].plot(np.sqrt(compute_var(ds4MFR, iP, ITime, D=D, removeBG=removeBG, ysym=ysym)), y/D, c=Colrs[2], lw=1.0, ls='--',   label='BG MFR')
        axes[1,iP].plot(np.sqrt(compute_var(ds4Ine, iP, ITime, D=D, removeBG=removeBG, ysym=ysym)), y/D, c=Colrs[3], lw=1.0, ls='--',   label='BG Inertial')
        axes[1,iP].set_xlabel('y/D [-]')
        if iP==0:
            axes[0,0].legend(fontsize=8, ncol=2, loc='upper center')

    for ax in np.array(axes).flatten():
        ax.tick_params(direction='in', top=True, right=True, labelright=False, labeltop=False, which='both')
        ax.set_xlim([0,2.0])

    axes[0,0].set_ylabel(r'$\sigma$'+'\n horizontal plane [m/s]')
    axes[1,0].set_ylabel(r'$\sigma$'+'\n horizontal plane [m/s]')
    figname = figbase+ '__Sigma_'+removeBG
    fig.suptitle(figname)
    fig.savefig(figsPath+figname+'.png')

plt.show()
