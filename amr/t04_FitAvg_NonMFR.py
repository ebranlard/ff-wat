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

figsPath='_figs/'
outPath='_out/'

# --- Parameters 
symmetric=False
# symmetric=True
removeBG='0WT'
removeBG='Outer'
# removeBG='Zero'

case = 'neutral'
# case = 'stable' 


label=case+'_rm'+removeBG
if symmetric:
    label+='_sym'


# --- Derived parameters
U0, dt, D, xyWT1, xyWT2, xPlanes = getSimParamsAMR(case)
figbase = figsPath+'KFit_'+label
outbase = outPath +'KFit_'+label
datapath0 = os.path.join('03-iea15-0WT/', case, 'processedData')
datapath2 = os.path.join('02-iea15-2WT/', case, 'processedData')

ds01= xarray.open_dataset(os.path.join(datapath0,'HubHeightWT1.nc_small')) # Wake behind turbine 1 - No turbine
ds21= xarray.open_dataset(os.path.join(datapath2,'HubHeightWT1.nc_small')) # Wake behind turbine 1
ds02= xarray.open_dataset(os.path.join(datapath0,'HubHeightWT2.nc_small')) # Wake behind turbine 2 - No turbine
ds22= xarray.open_dataset(os.path.join(datapath2,'HubHeightWT2.nc_small')) # Wake behind turbine 2
ds = ds21

# --- Intersection of time
ITime = common_itime(ds01, ds21, ds02, ds22)

KS = np.zeros((2, len(xPlanes), 2)) # iWT, xPlanes, (Kd, Kg)
# DU = np.zeros((len(xPlanes))) # iWT, xPlanes, (Kd, Kg)

fig,axes = plt.subplots(2, len(xPlanes), sharex=True, sharey=True, figsize=(12.8,4.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.06, right=0.98, top=0.88, bottom=0.11, hspace=0.20, wspace=0.20)
for iP, x in enumerate(xPlanes):
    axes[0,iP].set_title('x= {}D'.format(int(x)))
    var, u, k_fit, ks_fit, fitter = ds_fitK(ds21, iP, ITime, D, U0, symmetric=symmetric, ds0=ds01, removeBG=removeBG)
    KS[0,iP, :] = ks_fit
    axes[0,iP].plot(np.sqrt(var), var.y/D,        label='sqrt(Var)')
    axes[0,iP].plot(k_fit       , var.y/D,  'k--' ,label='fit')
    axes[0,iP].text(1,1, fitter.coeffsToString(sep='\n'), ha='center', va='center')

for iP, x in enumerate(xPlanes):
    var, u, k_fit, ks_fit, fitter = ds_fitK(ds22, iP, ITime, D, U0, symmetric=symmetric, ds0=ds02, removeBG=removeBG)
    KS[1,iP, :] = ks_fit
    axes[1,iP].plot(np.sqrt(var), var.y/D,        label=r'$\sqrt{\sigma^2-o}$')
    axes[1,iP].plot(k_fit       , var.y/D, 'k--' ,label='fit')
    axes[1,iP].text(1,1, fitter.coeffsToString(sep='\n'), ha='center', va='center')

    # Compute some kind of K factors based on first row
    print('---- iP', iP)
    kd =         KS[1,iP,0]    
    kg =         KS[1,iP,1]    
    print('kd={:.3f} kg={:.3f}'.format(kd, kg))
    kd = np.sqrt(KS[0,iP,0]**2 + KS[0,-1,0]**2)
    kg = np.sqrt(KS[0,iP,1]**2 + KS[0,-1,1]**2)
    print('kd={:.3f} kg={:.3f}'.format(kd, kg))
    K, KD, KG, du, gu = kWAT1D(var.y, U0, u, D, kd=kd, kg=kg)
    axes[1,iP].plot(K           , var.y/D,  ':',     label='superp')


#     axes[0,iP].plot(uG , y/D, 'k--', label='u_fit')
#     axes[1,iP].plot(KG**2, y/D, label='KG fit')
#     axes[1,iP].plot(varbg+y*0, y/D, label='Var')
    axes[0,iP].set_xlim([0.0,1.9])
    axes[1,iP].set_xlim([0.0,1.9])

    axes[1,iP].set_xlabel('k [-]')
axes[1,0].legend(fontsize=8, ncol=2, loc='upper center', borderpad=0)
axes[0,0].set_ylabel('y/D [-]')
axes[1,0].set_ylabel('y/D [-]')
figname = figbase
fig.suptitle(figname)
fig.savefig(figname+'.png')

# import pdb; pdb.set_trace()

pkl = PickleFile(data=KS)
# pkl['KS']=KS
pkl.write(outbase+'.pkl')




# "TKE"
# 2*tke_ds(ds01.isel(x=iP, it=ITime),   axis='it'),      
# 2*tke_ds(ds21.isel(x=iP, it=ITime),   axis='it'),      
# #          ds21.isel(x=iP, it=ITime).var(dim='it').u+tke0
# # 2*tke_ds(ds02.isel(x=iP, it=ITime),   axis='it'),      
# 2*tke_ds(ds22.isel(x=iP, it=ITime),   axis='it'),      
#          ds22.isel(x=iP, it=ITime).var(dim='it').u+tke0

# Var
# ds.isel(x=iP).var(dim='it').u




# axes[0,iP].set_title('x= {}D'.format(int(x)))
# #             axes[0,iP].plot(tke_ds(ds01.isel(x=iP, it=ITime), axis='it'), y/D, label='No turbine')
# #             axes[0,iP].plot(tke_ds(ds21.isel(x=iP, it=ITime), axis='it'), y/D, label='Turbine')
# #             axes[1,iP].plot(tke_ds(ds02.isel(x=iP, it=ITime), axis='it'), y/D, label='No turbine')
# #             axes[1,iP].plot(tke_ds(ds22.isel(x=iP, it=ITime), axis='it'), y/D, label='Turbine')
# varbg1 =  ds01.isel(x=iP, it=ITime).var(dim='it').u
# varbg2 =  ds02.isel(x=iP, it=ITime).var(dim='it').u
# 
# axes[0,iP].plot(ds01.isel(x=iP, it=ITime).var(dim='it').u - varbg1, y/D, label=None)
# axes[0,iP].plot(ds21.isel(x=iP, it=ITime).var(dim='it').u - varbg1, y/D, label='Turbine')
# axes[1,iP].plot(ds02.isel(x=iP, it=ITime).var(dim='it').u - varbg2, y/D, label=None)
# axes[1,iP].plot(ds22.isel(x=iP, it=ITime).var(dim='it').u - varbg2, y/D, label='Turbine')
# 
# 
# axes[0,iP].plot(du1*5, y/D, '--', c=python_colors(2) , label='Deficit')
# axes[1,iP].plot(du1*5, y/D, '--', c=python_colors(2) , label='Deficit')
# 
# axes[0,iP].plot(gu_r1*5, y/D, ':' , c=(0.5,0.5,0.5),  label='Gradient')
# axes[1,iP].plot(gu_r1*5, y/D, ':' , c=(0.5,0.5,0.5),  label='Gradient')
# 
# axes[0,iP].plot(k1**2, y/D, 'k--', label='k')
# axes[1,iP].plot(k2**2, y/D, 'k--', label='k')
# 
# 


plt.show()
