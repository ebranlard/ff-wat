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
import welib.weio
from welib.weio.vtk_file import VTKFile
from welib.weio.fast_input_file import FASTInputFile
from helper_functions import *
import xarray as xr


# --- Parameters
kmin=-0.01
folder ='steady/'
SIMS =  []
# SIMS +=  ['FF-NoWAT']
SIMS +=  ['FF-WAT-kd00_kg00']; kmax=1.11
# SIMS +=  ['FF-WAT-kd00_kg00-Offset05']; kmax=0.5
# SIMS +=  ['FF-WAT-kd00_kg00-Offset20']; kmax=2
# SIMS +=  ['FF-WAT-kd00_kg00-Offset20-Zero']; kmax=2
# SIMS +=  ['FF-WAT-kd00_kg00-Offset20-Zero-Longer']; kmax=2
# SIMS +=  ['FF-WAT-kd00_kg10']; kmax=1.11
SIMS +=  ['FF-WAT-kd10_kg10']; kmax=1.11
# SIMS +=  ['FF-WAT-kd10_kg00']; kmax=0.35
# SIMS +=  ['FF-WAT-kd01_kg01']; kmax=0.25
# SIMS +=  ['FF-WAT-kd10_kg10-ModProj1']
# SIMS +=  ['FF-WAT-kd10_kg10-ModProj3']
# SIMS +=  ['FF-WAT-ModW-Cart']
# SIMS +=  ['FF-WAT-ModW-Polar']


zHub=150
U0= 8
nPlanes=6
D=240
levels=np.linspace(2,10,20)/U0



# --- Derived
outDir = os.path.join(folder, '_processedData')
LBLS=[l.replace('FF-','') for l in SIMS]

fstffile=os.path.join(folder, SIMS[1]+'.fstf')
fstf=FASTInputFile(fstffile)
kDef  =fstf['WAT_k_Def']
kGrad =fstf['WAT_k_Grad']
print('k',kDef, kGrad)
# kDef  = 0.5
# kGrad = 0.5

file1=os.path.join(outDir,'planes_{}.nc'.format(SIMS[0]))
file2=os.path.join(outDir,'planes_{}.nc'.format(SIMS[1]))

if not os.path.exists(file2):
    I, iMax, nDigits, sFmt = vtkIndices(folder, SIMS[1]+'.Low.DisYZ01')
    print(I[0],I[-1])

    I = np.asarray(I)
    iStart=100
    ISel = I[I>iStart]
    dsw = loadAllPlanes(folder, SIMS[1], nPlanes, ISel, outDir=outDir, verbose=False)
    dsr = xr.open_dataset(file1)
else:
    dsr = xr.open_dataset(file1)
    dsw = xr.open_dataset(file2)
print('Number of time Steps:', len(dsr.it), len(dsw.it))

xPlanes = dsr.ix
y       = dsr.y
z       = dsr.z



# --- Plot average profile at hub height

R=D/2

var_max = np.max(dsr.u.sel(z=zHub).var(dim='it').values)
print('>>> Var Max',var_max)

# fig,axes = plt.subplots(8, nPlanes, sharey=True, figsize=(12.4,9.8)) # (6.4,4.8)
fig,axes = plt.subplots(7, nPlanes, sharey=True, figsize=(12.4,9.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.05, right=0.99, top=0.93, bottom=0.05, hspace=0.90, wspace=0.23)
for i,ix in enumerate(xPlanes):
    j=-1
    Uh_r = dsr.u.sel(z=zHub, ix=ix).mean(dim='it').values
    Uh_w = dsw.u.sel(z=zHub, ix=ix).mean(dim='it').values
    vh_r = dsr.u.sel(z=zHub, ix=ix).var (dim='it').values
    vh_w = dsw.u.sel(z=zHub, ix=ix).var (dim='it').values

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

    kt_r = (kd_r + kg_r)
    kt_w = (kd_w + kg_w)

#  p%WAT_k_Def *  abs(1 - ((xd%Vx_wind_disk_filt(i)+y%Vx_wake2(iy,iz,i))/xd%Vx_wind_disk_filt(i)) ) & 
#  p%WAT_k_Grad/xd%Vx_wind_disk_filt(i) * R * ( abs(dvdr) + abs(dvdtheta_r) )


#     j+=1
#     axes[j,i].plot(Uh_r/U0, y/D, label=LBLS[0])
#     axes[j,i].plot(Uh_w/U0, y/D, label=LBLS[1])
#     if i==0:
#         axes[j,i].legend()
#     axes[j,i].set_title('x= {}D'.format(int(ix)))

    j+=1
    axes[j,i].plot(DU_r, y/D, label='No WAT')
    axes[j,i].plot(DU_w, y/D, label='WAT')
    if i==0:
        axes[j,i].set_title(r'Deficit $U_d=1-U/U_0$')
    else:
        axes[j,i].set_title('x= {}D'.format(int(ix)))

    # ---Gradient
    j+=1
    axes[j,i].plot(gU_r, y/D, label='No WAT')
    axes[j,i].plot(gU_w, y/D, label='WAT')
#     axes[j,i].plot(gU_w2, y/D,'k--')
    if i==0:
        axes[j,i].set_title(r'Gradient: $dU/dy  R/U_0$')

    # --- Ks
    j+=1
    axes[j,i].plot(kd_r, y/D, label='No WAT')
    axes[j,i].plot(kd_w, y/D, label='WAT')
    if i==0:
        axes[j,i].set_title(r'$k_d |U_d|$')
    axes[j,i].set_xlim([kmin,kmax+0.1])

    j+=1
    axes[j,i].plot(kg_r, y/D, label='No WAT')
    axes[j,i].plot(kg_w, y/D, label='WAT')
    if i==0:
        axes[j,i].set_title(r'$k_g |dU/dy|$')
    axes[j,i].set_xlim([kmin,kmax+0.1])

    j+=1
    axes[j,i].plot(kt_r, y/D, label='No WAT')
    axes[j,i].plot(kt_w, y/D, label='WAT')
    if i==0:
        axes[j,i].set_title(r'$k$')
    axes[j,i].set_xlim([kmin,kmax+0.1])


    j+=1
    axes[j,i].plot(kt_r**2, y/D, label='No WAT')
    axes[j,i].plot(kt_w**2, y/D, label='WAT')
    if i==0:
        axes[j,i].set_title(r'$k^2$')
    axes[j,i].set_xlim([kmin,kmax**2])


    j+=1
    axes[j,i].plot(vh_r, y/D, label='No WAT')
    axes[j,i].plot(vh_w, y/D, label='WAT')
    if i==0:
        axes[j,i].set_title('Variance: var(U)')
#     axes[j,i].axvline(kmax**2, ls='--', c='k')
    axes[j,i].set_xlim([kmin,kmax**2])

axes[0,i].set_ylabel('y/D [-]')
axes[1,i].set_ylabel('y/D [-]')
axes[2,i].set_ylabel('y/D [-]')
for ax in np.asarray(axes).ravel():
    ax.tick_params(direction='in', top=True, right=True, labelright=False, labeltop=False, which='both')
    ax.set_ylim([0,2])

# fig.tight_layout()


fig.suptitle('{} vs {}'.format(LBLS[0], LBLS[1]))
fig.savefig('_figs/Compare_{}{}.png'.format(LBLS[0], LBLS[1]))



plt.show()
