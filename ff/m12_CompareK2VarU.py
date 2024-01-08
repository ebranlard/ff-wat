import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import welib.weio as weio
from welib.weio.vtk_file import VTKFile
from welib.weio.fast_input_file import FASTInputFile
from helper_functions import *
import xarray as xr


# --- Parameters
kmin=-0.01
folder ='steady/'
# SIM +=  ['FF-NoWAT']
# SIM =  ['FF-WAT-kd00_kg00']; kmax=1.11
# SIM +=  ['FF-WAT-kd00_kg00-Offset05']; kmax=0.5
# SIM +=  ['FF-WAT-kd00_kg00-Offset20']; kmax=2
# SIM +=  ['FF-WAT-kd00_kg00-Offset20-Zero']; kmax=2
# SIM +=  ['FF-WAT-kd00_kg00-Offset20-Zero-Longer']; kmax=2
# SIM =  ['FF-WAT-kd00_kg10']; kmax=1.11
SIM =  'FF-WAT-kd10_kg10'; kmax=1.11
# SIM +=  ['FF-WAT-kd10_kg00']; kmax=0.35
# SIM +=  ['FF-WAT-kd01_kg01']; kmax=0.25
# SIM +=  ['FF-WAT-kd10_kg10-ModProj1']
# SIM +=  ['FF-WAT-kd10_kg10-ModProj3']
# SIM +=  ['FF-WAT-ModW-Cart']
# SIM +=  ['FF-WAT-ModW-Polar']


zHub=150
U0= 8
nPlanes=6
D=240
levels=np.linspace(2,10,20)/U0



# --- Derived
outDir = os.path.join(folder, '_processedData')
LBLS=SIM.replace('FF-','')

fstffile=os.path.join(folder, SIM+'.fstf')
fstf=FASTInputFile(fstffile)
kDef  =fstf['WAT_k_Def']
kGrad =fstf['WAT_k_Grad']
print('k',kDef, kGrad)
# kDef  = 0.5
# kGrad = 0.5

file1=os.path.join(outDir,'planes_{}.nc'.format(SIM))
dsr = xr.open_dataset(file1)
print('Number of time Steps:', len(dsr.it))

xPlanes = dsr.ix
y       = dsr.y
z       = dsr.z



# --- Plot average profile at hub height

R=D/2

var_max = np.max(dsr.u.sel(z=zHub).var(dim='it').values)
print('>>> Var Max',var_max)

# fig,axes = plt.subplots(8, nPlanes, sharey=True, figsize=(12.4,9.8)) # (6.4,4.8)
fig,axes = plt.subplots(3, nPlanes, sharey=True, figsize=(12.4,9.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.05, right=0.99, top=0.93, bottom=0.05, hspace=0.90, wspace=0.23)
for i,ix in enumerate(xPlanes):
    j=-1
    Uh_r = dsr.u.sel(z=zHub, ix=ix).mean(dim='it').values
    vh_r = dsr.u.sel(z=zHub, ix=ix).var (dim='it').values
    DU_r = 1-Uh_r/U0
    gU_r  = np.gradient(Uh_r, y)*R/U0
    gU_r2 = np.gradient(-DU_r, y/R)
    kd_r = kDef * np.abs(DU_r)
    kg_r = kGrad * np.abs(gU_r)
    kt_r = (kd_r + kg_r)
    j+=1
    axes[j,i].plot(DU_r, y/D, label='WAT')
    if i==0:
        axes[j,i].set_title(r'Deficit $U_d=1-U/U_0$')
    else:
        axes[j,i].set_title('x= {}D'.format(int(ix)))

    # ---Gradient
    j+=1
    axes[j,i].plot(gU_r, y/D, label='WAT')
    if i==0:
        axes[j,i].set_title(r'Gradient: $dU/dy  R/U_0$')

#     # --- Ks
#     j+=1
#     axes[j,i].plot(kd_r, y/D, label='WAT')
#     if i==0:
#         axes[j,i].set_title(r'$k_d |U_d|$')
#     axes[j,i].set_xlim([kmin,kmax+0.1])
# 
#     j+=1
#     axes[j,i].plot(kg_r, y/D, label='WAT')
#     if i==0:
#         axes[j,i].set_title(r'$k_g |dU/dy|$')
#     axes[j,i].set_xlim([kmin,kmax+0.1])

#     j+=1
#     axes[j,i].plot(kt_r, y/D, label='WAT')
#     if i==0:
#         axes[j,i].set_title(r'$k$')
#     axes[j,i].set_xlim([kmin,kmax+0.1])

    j+=1
    axes[j,i].plot(kt_r**2, y/D, label='k^2')
    axes[j,i].plot(vh_r, y/D, 'k--',label='Variance')
    if i==0:
        axes[j,i].set_title(r'$k^2$ and Variance')
#     axes[j,i].set_xlim([kmin,kmax**2])
    axes[j,i].legend()

axes[0,i].set_ylabel('y/D [-]')
axes[1,i].set_ylabel('y/D [-]')
axes[2,i].set_ylabel('y/D [-]')
for ax in np.asarray(axes).ravel():
    ax.tick_params(direction='in', top=True, right=True, labelright=False, labeltop=False, which='both')
#     ax.set_ylim([0,2])

# fig.tight_layout()


fig.suptitle('{} '.format(LBLS))
fig.savefig('_figs/CompareK2VarU_{}.png'.format(LBLS))



plt.show()
