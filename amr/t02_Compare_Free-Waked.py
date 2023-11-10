""" 
Compare 
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray
# Local 
import weio
from welib.essentials import *
from welib.tools.colors import *
from helper_functions import *

IPlot = [  ]
# IPlot += [ 0] # Velocity at Point
# IPlot += [ 1] # Mean velocity
# IPlot += [ 2] # TKE
IPlot += [ 3] # Var
# IPlot += [ 4] # k^2 / Var

# --- Parameters 
figsPath='_figs/'

case = 'neutral'
case = 'stable' 

for case in ['neutral', 'stable']:
# for case in ['neutral']:

    # --- Derived parameters
    U0, dt, D, xyWT1, xyWT2, xPlanes = getSimParamsAMR(case)
    figbase = figsPath+'{}_02WT_'.format(case)
    datapath0 = os.path.join('03-iea15-0WT/', case, 'processedData')
    datapath2 = os.path.join('02-iea15-2WT/', case, 'processedData')

    ds01= xarray.open_dataset(os.path.join(datapath0,'HubHeightWT1.nc'))
    ds21= xarray.open_dataset(os.path.join(datapath2,'HubHeightWT1.nc'))
    ds02= xarray.open_dataset(os.path.join(datapath0,'HubHeightWT2.nc'))
    ds22= xarray.open_dataset(os.path.join(datapath2,'HubHeightWT2.nc'))

    # --- Intersection of time
    ITime = common_itime(ds01, ds21)

    y = ds01.y.values
    iy0 = np.argmin(np.abs(y-0    ))

    # --------------------------------------------------------------------------------}
    # --- Plot velocity at a given point
    # --------------------------------------------------------------------------------{
    if 0 in IPlot:
        # --- Plot mean at hub height as function of time
        iy=0
        for iy in [0, iy0, -1]:
            fig,axes = plt.subplots(2, len(xPlanes), sharex=True, sharey=True, figsize=(12.8,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.06, right=0.99, top=0.88, bottom=0.11, hspace=0.20, wspace=0.20)
            for iP, x in enumerate(xPlanes):
                tArrival  = (x*D)    /(U0*(1-0.3))
                tArrival2 = ((x+7)*D)/(U0*(1-0.3))
                axes[0,iP].set_title('x= {}D'.format(int(x)))
                axes[0,iP].plot(ds01.it*dt, ds01.isel(x=iP, y=iy).u, '-', label='No turbine')
                axes[0,iP].plot(ds21.it*dt, ds21.isel(x=iP, y=iy).u, '--', label='Turbine sim')
                axes[0,iP].axvline(tArrival, c='k', ls=':')
                axes[0,iP].set_xlabel('Time [s]')

                axes[1,iP].set_title('x= {}D'.format(int(x)))
                axes[1,iP].plot(ds02.it*dt, ds02.isel(x=iP, y=iy).u, '-', label='No turbine')
                axes[1,iP].plot(ds22.it*dt, ds22.isel(x=iP, y=iy).u, '--', label='Turbine sim')
                axes[1,iP].axvline(tArrival,  c='k', ls=':')
                axes[1,iP].axvline(tArrival2, c='k', ls='--')
                axes[1,iP].set_xlabel('Time [s]')
            axes[0,0].set_ylabel('Velocity [m/s]')
            axes[1,0].set_ylabel('Velocity [m/s]')
        #     axes[1,0].set_xlim([0,300])
            axes[0,0].legend()
            y_D=float((ds01.y[iy])/D)
            fig.suptitle('Velocity at point y={:.1f}D'.format(y_D))
            fig.savefig(figbase+ 'PointVelocity_{:.1f}D.png'.format(y_D))

    # --------------------------------------------------------------------------------}
    # --- Plot mean velocity at downstream planes
    # --------------------------------------------------------------------------------{

    if 1 in IPlot:
        # --- Plot mean at hub height as function of time
        fig,axes = plt.subplots(2, len(xPlanes), sharex=True, sharey=True, figsize=(12.8,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.06, right=0.98, top=0.88, bottom=0.11, hspace=0.20, wspace=0.20)
        for iP, x in enumerate(xPlanes):
            axes[0,iP].set_title('x= {}D'.format(int(x)))
            axes[0,iP].plot(ds01.isel(x=iP, it=ITime).mean(dim='it').u, y/D, label='No turbine')
            axes[0,iP].plot(ds21.isel(x=iP, it=ITime).mean(dim='it').u, y/D, label='Turbine')
            axes[1,iP].plot(ds02.isel(x=iP, it=ITime).mean(dim='it').u, y/D, label='No turbine')
            axes[1,iP].plot(ds22.isel(x=iP, it=ITime).mean(dim='it').u, y/D, label='Turbine')
            #axes[0,iP].set_xlabel('u [m/s]')
            axes[1,iP].set_xlabel('u [m/s]')
        axes[0,0].set_ylabel('y/D [-]')
        axes[1,0].set_ylabel('y/D [-]')
        figname = figbase+ 'MeanHoriProfilesHH'
        fig.suptitle(figname + ' n:{}/{}'.format(it_min,it_max))
        fig.savefig(figname+'.png')


    # --------------------------------------------------------------------------------}
    # --- Plot TKE
    # --------------------------------------------------------------------------------{
    if 2 in IPlot:
        # --- Plot mean at hub height as function of time
        fig,axes = plt.subplots(2, len(xPlanes), sharex=True, sharey=True, figsize=(12.8,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.06, right=0.98, top=0.88, bottom=0.11, hspace=0.20, wspace=0.20)

        y=ds21.y.values

        for iP, x in enumerate(xPlanes):
            tke0 = tke_ds(ds01.isel(x=iP, it=ITime),   axis='it').mean()
            axes[0,iP].set_title('x= {}D'.format(int(x)))
            axes[0,iP].plot(2*tke_ds(ds01.isel(x=iP, it=ITime),   axis='it'),        y/D, label='No turbine')
            axes[0,iP].plot(2*tke_ds(ds21.isel(x=iP, it=ITime),   axis='it'),        y/D, label='Turbine', c='k')
            axes[0,iP].plot(         ds21.isel(x=iP, it=ITime).var(dim='it').u+tke0, y/D, label='Turbine - Var u/2')
            axes[1,iP].plot(2*tke_ds(ds02.isel(x=iP, it=ITime),   axis='it'),        y/D, label='No turbine')
            axes[1,iP].plot(2*tke_ds(ds22.isel(x=iP, it=ITime),   axis='it'),        y/D, label='Turbine', c='k')
            axes[1,iP].plot(         ds22.isel(x=iP, it=ITime).var(dim='it').u+tke0, y/D, label='Turbine - Var u/2')

            #axes[0,iP].set_xlabel('u [m/s]')
            axes[0,0].legend(fontsize=9)

            axes[1,iP].set_xlabel('2 TKE [m^2/s^2]')
            axes[0,iP].set_xlim([-0.5,5])
        axes[0,0].set_ylabel('y/D [-]')
        axes[1,0].set_ylabel('y/D [-]')
        figname = figbase+ 'MeanTKEHori'
        fig.suptitle(figname + ' n:{}/{}'.format(it_min,it_max))
        fig.savefig(figname+'.png')

    # --------------------------------------------------------------------------------}
    # --- Plot Var
    # --------------------------------------------------------------------------------{
    if 3 in IPlot:
        # --- Plot mean at hub height as function of time
        fig,axes = plt.subplots(2, len(xPlanes), sharex=True, sharey=True, figsize=(12.8,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.06, right=0.98, top=0.88, bottom=0.11, hspace=0.20, wspace=0.20)

        y=ds21.y.values

        for iP, x in enumerate(xPlanes):
            tke0 = tke_ds(ds01.isel(x=iP, it=ITime),   axis='it').mean()
            axes[0,iP].set_title('x= {}D'.format(int(x)))
            axes[0,iP].plot(ds01.isel(x=iP).var(dim='it').u,   y/D, label='No turbine')
            axes[0,iP].plot(ds21.isel(x=iP).var(dim='it').u,   y/D, label='Turbine', c='k')
            axes[1,iP].plot(ds02.isel(x=iP).var(dim='it').u,  y/D, label='No turbine')
            axes[1,iP].plot(ds22.isel(x=iP).var(dim='it').u,  y/D, label='Turbine', c='k')

            #axes[0,iP].set_xlabel('u [m/s]')
            axes[0,0].legend(fontsize=9)

            axes[1,iP].set_xlabel('Var [m^2/s^2]')
            axes[0,iP].set_xlim([-0.5,5])
        axes[0,0].set_ylabel('y/D [-]')
        axes[1,0].set_ylabel('y/D [-]')
        figname = figbase+ 'MeanVAR'
        fig.suptitle(figname + ' n:{}/{}'.format(it_min,it_max))
        fig.savefig(figname+'.png')


    # --------------------------------------------------------------------------------}
    # --- Plot k**2 / var
    # --------------------------------------------------------------------------------{
    kDef  = 3.5
    kGrad = 3.5

    if 4 in IPlot:
        # --- Plot mean at hub height as function of time
        fig,axes = plt.subplots(2, len(xPlanes), sharex=True, sharey=True, figsize=(12.8,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.06, right=0.98, top=0.88, bottom=0.11, hspace=0.20, wspace=0.20)

        y=ds21.y.values

        for iP, x in enumerate(xPlanes):
            """ 

             k = k_def / U0  | ~u |  + k_grad * D /(2U0) |  du~/dr + 1/r du~/dth |

  - Extract LES results in meandering frame of reference (MFOR) using wake tracker
  - Extract mean deficit u_ and gradient du_ in MFOR
  - Subtract mean deficit u_ from wake time series in MFOR to get u_perturb(t)
  - Compute var of u_perturb(t,r)
  - Calibrate   var(u_perturb) = (k1 u + k2 du)^2-var(U0)  

            """
            # Instantaneous velocity for a given plane
            i_du1= 1-ds21.isel(x=iP, it=ITime).u/U0
            i_du2= 1-ds22.isel(x=iP, it=ITime).u/U0
            # Instantaneous gradients
            #     gU_r  = np.gradient(Uh_r, y)*(D/2)/U0

            # Time averaged velocity profiles
            du1= 1-ds21.isel(x=iP, it=ITime).mean(dim='it').u.values/U0
            du2= 1-ds22.isel(x=iP, it=ITime).mean(dim='it').u.values/U0
            # Time averaged gradients
            gu_r1 = np.gradient(du1, y)*(D/2)/U0
            gu_r2 = np.gradient(du1, y)*(D/2)/U0
            # Factor k
            k1 = kDef * np.abs(du1) + kGrad * np.abs(gu_r1)
            k2 = kDef * np.abs(du2) + kGrad * np.abs(gu_r2)
# 

            axes[0,iP].set_title('x= {}D'.format(int(x)))
#             axes[0,iP].plot(tke_ds(ds01.isel(x=iP, it=ITime), axis='it'), y/D, label='No turbine')
#             axes[0,iP].plot(tke_ds(ds21.isel(x=iP, it=ITime), axis='it'), y/D, label='Turbine')
#             axes[1,iP].plot(tke_ds(ds02.isel(x=iP, it=ITime), axis='it'), y/D, label='No turbine')
#             axes[1,iP].plot(tke_ds(ds22.isel(x=iP, it=ITime), axis='it'), y/D, label='Turbine')
            varbg1 =  ds01.isel(x=iP, it=ITime).var(dim='it').u
            varbg2 =  ds02.isel(x=iP, it=ITime).var(dim='it').u

            axes[0,iP].plot(ds01.isel(x=iP, it=ITime).var(dim='it').u - varbg1, y/D, label=None)
            axes[0,iP].plot(ds21.isel(x=iP, it=ITime).var(dim='it').u - varbg1, y/D, label='Turbine')
            axes[1,iP].plot(ds02.isel(x=iP, it=ITime).var(dim='it').u - varbg2, y/D, label=None)
            axes[1,iP].plot(ds22.isel(x=iP, it=ITime).var(dim='it').u - varbg2, y/D, label='Turbine')


            axes[0,iP].plot(du1*5, y/D, '--', c=python_colors(2) , label='Deficit')
            axes[1,iP].plot(du1*5, y/D, '--', c=python_colors(2) , label='Deficit')

            axes[0,iP].plot(gu_r1*5, y/D, ':' , c=(0.5,0.5,0.5),  label='Gradient')
            axes[1,iP].plot(gu_r1*5, y/D, ':' , c=(0.5,0.5,0.5),  label='Gradient')

            axes[0,iP].plot(k1**2, y/D, 'k--', label='k')
            axes[1,iP].plot(k2**2, y/D, 'k--', label='k')



            #axes[0,iP].set_xlabel('u [m/s]')
            axes[0,0].legend(ncol=2, fontsize=8)
            axes[0,iP].set_xlim([-0.5,5])

            axes[1,iP].set_xlabel('Var [m^2/s^2]')
            #axes[1,iP].set_xlabel('TKE [m^2/s^2]')
        axes[0,0].set_ylabel('y/D [-]')
        axes[1,0].set_ylabel('y/D [-]')
        figname = figbase+ 'MeanK_Hori'
        fig.suptitle(figname + ' n:{}/{}'.format(it_min,it_max))
        fig.savefig(figname+'.png')



# import pdb; pdb.set_trace()

plt.show()
