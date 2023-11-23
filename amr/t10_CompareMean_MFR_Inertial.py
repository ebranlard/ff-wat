""" 
Plot mean horizontal velocity profile at each downstream plane between WT1&2
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from helper_functions import *
import warnings
warnings.simplefilter("ignore", category=RuntimeWarning)


def plotMFRInertialComp(caseName, Case, outDir='_out', figDir='_fig', iTimeMin=500):
    print(f'-----------------------------------------------------------------------')
    print(f'--- plotMFRInertialComp - Case: {caseName}')
    print(f'-----------------------------------------------------------------------')
    if Case['nWT']==0:
        WARN('Skipping case (script only makes sense for one or two WT)')
        return

    outDirHori = os.path.join(outDir  , 'lines')
    figDirComp  = os.path.join(figDir, '_figs_compMFR')
    if not os.path.exists(outDirHori):
        os.makedirs(outDirHori)
    if not os.path.exists(figDirComp):
        os.makedirs(figDirComp)

    U0, dt, D, xyWT, HubHeight, xPlanes = getSimParamsAMR(Case['stability'])
    # --- Derived parameters
    figbase = '{}_MFR_vs_Inertial_Hori_'.format(caseName)

    ds1MFR = readDataSet(os.path.join(outDirHori,caseName+'_Meander_WT1.nc_small'))
    ds2MFR = readDataSet(os.path.join(outDirHori,caseName+'_Meander_WT2.nc_small'))
    ds1Ine = readDataSet(os.path.join(outDirHori,caseName+'_Inertial_WT1.nc_small'))
    ds2Ine = readDataSet(os.path.join(outDirHori,caseName+'_Inertial_WT2.nc_small'))
    #  OWT case
    caseName0 = caseName.replace('2WT','0WT').replace('1WT','0WT')
    ds3MFR = readDataSet(os.path.join(outDirHori, caseName0 + '_IN'+ caseName + '_Meander_WT1.nc_small'))
    ds4MFR = readDataSet(os.path.join(outDirHori, caseName0 + '_IN'+ caseName + '_Meander_WT2.nc_small'))
    ds3Ine = readDataSet(os.path.join(outDirHori, caseName0 + '_Inertial_WT1.nc_small'))
    ds4Ine = readDataSet(os.path.join(outDirHori, caseName0 + '_Inertial_WT2.nc_small'))
    print('it',ds1MFR.it.values[0],ds1MFR.it.values[-1]  )

    y = ds1MFR.y.values
    t = ds1MFR.it.values * dt
    ITime = common_itime(ds1MFR,)
    print('ITime 2WT',ITime[0], ITime[-1])
    ITime = common_itime(ds3MFR,)
    print('ITime 0WT',ITime[0], ITime[-1])
    ITime = common_itime(ds1MFR, ds3MFR)
    print('ITime com',ITime[0], ITime[-1])
    ITime = np.arange(iTimeMin, ITime[-1])
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
    fig.savefig(os.path.join(figDirComp,figname+'.png'))


    # --- Plot var at hub height as function of time
    ysym=False
    for removeBG in ['Zero','Outer']:
#     removeBG='Zero'
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
        fig.savefig(os.path.join(figDirComp,figname+'.png'))


if __name__ == '__main__':
    caseNames = [] 
    caseNames += ['neutral2WT']
    caseNames += ['stable2WT'] 
    caseNames += ['unstable2WT'] 
#     caseNames += ['neutral1WT']
#     caseNames += ['stable1WT'] 
#     caseNames += ['unstable1WT'] 

    Cases = dict((k,v) for k,v in AllCases.items() if k in caseNames)

    for caseName, Case in Cases.items():
        with FailSafe(caseName, False):
            plotMFRInertialComp(caseName, Case, outDir='_out')

    plt.show()
