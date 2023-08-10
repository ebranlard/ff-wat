# 
# - Task 0a: Sanity check: FF output is as intended
#    - Uniform inflow (var(U0)=0) with WAT (W) and without (0), two turbines
#        - Check that mean (u_) of output planes are the same for (0) and (W)
#        - Check that Turbine 2 gets more TI with WAT
#        - Check that std of output planes of sim (W) is  (k1 u + k2 du)
#        - Compare wake position. Impact on meandering for increasing k.
# 
# - Task 0b: Get ready for meandering frame, and how much does var(U0) polute our signal
#    - Simulations with background turbulence (p_rerably from AMR wind but can use BTS), (W) and (0)
#         - Use wake tracker to extract outputs in MFOR   u(r, x) 
#         - Extract mean deficit u_ and gradient du_ in MFOR
#         - Subtract mean deficit u_ from wake time series in MFOR to get u_perturb
#         - Compute var of u_perturb(t,r,x). For case (0), this is var(U0)_0
#         - Check if var(u_perturb)(r,x) is close to (k1 u + k2 du)^2 - var(U0)_0
#         - Compute TKE profiles for the sake of it
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import weio
from weio.vtk_file import VTKFile
from helper_functions import *
import xarray as xr


# --- Parameters
folder ='steady/'
simbase_r = 'FF-NoWAT'
simbasewat = 'FF-WAT'
zHub=150
U0= 8
nPlanes=6
D=240
levels=np.linspace(2,10,20)/U0

kDef = 1
kGrad = 1



# --- Derived
outDir = os.path.join(folder, '_processedData')


dsr = xr.open_dataset(os.path.join(outDir,'planesREF.nc'))
dsw = xr.open_dataset(os.path.join(outDir,'planesWAT.nc'))

xPlanes = dsr.ix
y       = dsr.y
z       = dsr.z






# --- Plot average profile at hub height
def mygradient(y, x):
    # Forward differences;
    #d = np.diff(y)/np.diff(x)
    # Central differences
    z1  = np.hstack((y[0],  y[:-1]))
    z2  = np.hstack((y[1:], y[-1]))
    dx1 = np.hstack((0, np.diff(x)))
    dx2 = np.hstack((np.diff(x), 0))

    d = (z2-z1) / (dx2+dx1)
    return d

R=D/2


fig,axes = plt.subplots(8, nPlanes, sharey=True, figsize=(12.4,9.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.05, right=0.99, top=0.95, bottom=0.05, hspace=0.90, wspace=0.23)
for ix in xPlanes:
    Uh_r = dsr.u.sel(z=zHub, ix=ix).mean(dim='it').values
    Uh_w = dsw.u.sel(z=zHub, ix=ix).mean(dim='it').values

    DU_r = 1-Uh_r/U0
    DU_w = 1-Uh_w/U0
    gU_r  = np.gradient(Uh_r, y)*R/U0
    gU_w  = np.gradient(Uh_w, y)*R/U0
    gU_r2 = np.gradient(-DU_r, y/R)
    gU_w2 = np.gradient(-DU_w, y/R)

    kd_r = kDef * np.abs(DU_r)
    kd_w = kDef * np.abs(DU_w)

    kg_r = kGrad * np.abs(gU_r)
    kg_w = kGrad * np.abs(gU_w)

    kt_r = kd_r + kg_w
    kt_w = kd_w + kg_w

#  p%WAT_k_Def *  abs(1 - ((xd%Vx_wind_disk_filt(i)+y%Vx_wake2(iy,iz,i))/xd%Vx_wind_disk_filt(i)) ) & 
#  p%WAT_k_Grad/xd%Vx_wind_disk_filt(i) * R * ( abs(dvdr) + abs(dvdtheta_r) )

    vh_r = dsr.u.sel(z=zHub, ix=ix).std (dim='it').values**2
    vh_w = dsw.u.sel(z=zHub, ix=ix).std (dim='it').values**2

    j=0
    axes[j,ix].plot(Uh_r/U0, y/D, label='No WAT')
    axes[j,ix].plot(Uh_w/U0, y/D, label='WAT')
    axes[j,ix].set_title('x= {}D'.format(int(ix)))

    j+=1
    axes[j,ix].plot(DU_r, y/D, label='No WAT')
    axes[j,ix].plot(DU_w, y/D, label='WAT')
    if ix==0:
        axes[j,ix].set_title(r'Deficit $U_d=1-U/U_0$')

    j+=1
    axes[j,ix].plot(gU_r, y/D, label='No WAT')
    axes[j,ix].plot(gU_w, y/D, label='WAT')
    axes[j,ix].plot(gU_w2, y/D,'k--')
    if ix==0:
        axes[j,ix].set_title(r'Gradient: $dU/dy  R/U_0$')

    j+=1
    axes[j,ix].plot(kd_r, y/D, label='No WAT')
    axes[j,ix].plot(kd_w, y/D, label='WAT')
    if ix==0:
        axes[j,ix].set_title(r'$k_d |U_d|$')

    j+=1
    axes[j,ix].plot(kg_r, y/D, label='No WAT')
    axes[j,ix].plot(kg_w, y/D, label='WAT')
    if ix==0:
        axes[j,ix].set_title(r'$k_g |dU/dy|$')

    j+=1
    axes[j,ix].plot(kt_r, y/D, label='No WAT')
    axes[j,ix].plot(kt_w, y/D, label='WAT')
    if ix==0:
        axes[j,ix].set_title(r'$k$')


    j+=1
    axes[j,ix].plot(kt_r**2, y/D, label='No WAT')
    axes[j,ix].plot(kt_w**2, y/D, label='WAT')
    if ix==0:
        axes[j,ix].set_title(r'$k^2$')


    j+=1
    axes[j,ix].plot(vh_r, y/D, label='No WAT')
    axes[j,ix].plot(vh_w, y/D, label='WAT')
    if ix==0:
        axes[j,ix].set_title('Variance: var(U)')



#     #ax.pcolormesh(ds1.isel(x=iP).mean(dim='it')['u'].T)
#     ax.contour(ds1.isel(x=iP).mean(dim='it')['u'].T)
axes[0,ix].set_ylabel('y/D [-]')
axes[1,ix].set_ylabel('y/D [-]')
axes[2,ix].set_ylabel('y/D [-]')
# fig.tight_layout()



plt.show()
