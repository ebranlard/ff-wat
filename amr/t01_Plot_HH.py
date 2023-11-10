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

# --- Script parameters
figsPath='_figs/'

Meander=True

# for nWT in [0,2]:
for nWT in [2]:
    #for case in ['neutral','stable']:
    for case in ['neutral']:

        # --- Derived parameters
        datapath = os.path.join('0{}-iea15-{}WT/'.format({2:2, 0:3}[nWT], nWT), case, 'processedData')
        U0, dt, D, xyWT1, xyWT2, xPlanes = getSimParamsAMR(case)
        if Meander:
            ds1 = xarray.open_dataset(os.path.join(datapath,'MeanderWT1.nc_small'))
            ds2 = xarray.open_dataset(os.path.join(datapath,'MeanderWT2.nc_small'))
            figbase = figsPath+'Meander_{}_{}WT_'.format(case,nWT)
        else:
            ds1 = xarray.open_dataset(os.path.join(datapath,'HubHeightWT1.nc_small'))
            ds2 = xarray.open_dataset(os.path.join(datapath,'HubHeightWT2.nc_small'))
            figbase = figsPath+'HubHeight_{}_{}WT_'.format(case,nWT)
        print('it',ds1.it.values[0],ds1.it.values[-1]  )

        y = ds1.y.values
        t = ds1.it.values * dt
        Itime = np.where(t>600)[0]
        t1=t[Itime[0]]
        t2=t[Itime[-1]]
        print('t1,t2', t1,t2)


        # --- Plot average wake deficit
        fig,axes = plt.subplots(1, 2, sharey=True, figsize=(12.4,5.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.06, right=0.98, top=0.88, bottom=0.11, hspace=0.20, wspace=0.20)
        for iP, x in enumerate(xPlanes):
            axes[0].plot(ds1.isel(x=iP, it=Itime).mean(dim='it', skipna=True).u, y/D, label='x= {}D'.format(int(x)))
            axes[1].plot(ds2.isel(x=iP, it=Itime).mean(dim='it', skipna=True).u, y/D, label='x= {}D'.format(int(x)))
            axes[0].set_xlabel('u [m/s]')
            axes[0].set_ylabel('y/D [-]')
        axes[0].legend()
        axes[0].set_title('WT 1')
        axes[1].set_title('WT 2')
        axes[1].legend()
        figname = figbase+ 'UTimeSeriesHH'
        fig.suptitle(figname)
        fig.savefig(figname+'.png')
        ds1.isel(x=iP, it=Itime).mean(dim='it', skipna=True).u
        import pdb; pdb.set_trace()


        # --- Plot mean at hub height as function of time
        fig,axes = plt.subplots(2, len(xPlanes), sharey=True, figsize=(12.8,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.06, right=0.98, top=0.88, bottom=0.11, hspace=0.20, wspace=0.20)
        for iP, x in enumerate(xPlanes):
            tArrival  = (x*D)    /(U0*(1-0.3))
            tArrival2 = ((x+7)*D)/(U0*(1-0.3))
            axes[0,iP].set_title('x= {}D'.format(int(x)))
            axes[0,iP].plot(t, ds1.isel(x=iP).mean(dim='y', skipna=True).u, label='WT1')
            axes[0,iP].axvline(tArrival, c='k', ls=':')

            axes[1,iP].plot(t, ds2.isel(x=iP).mean(dim='y', skipna=True).u, label='WT2')
            axes[1,iP].set_xlabel('Time [s]')
            axes[1,iP].axvline(tArrival,  c='k', ls=':')
            axes[1,iP].axvline(tArrival2, c='k', ls='--')
        axes[0,0].set_ylabel('Mean-y HH velocity [m/s]')
        axes[1,0].set_ylabel('Mean-y HH velocity [m/s]')
        figname = figbase+ 'MeanHoriProfilesHH'
        fig.suptitle(figname)
        fig.savefig(figname+'.png')



    # import pdb; pdb.set_trace()

plt.show()
