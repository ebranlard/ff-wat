
# 
# - Task 0a: Sanity check: FF output is as intended
#    - Uniform inflow (var(U0)=0) with WAT (W) and without (0), two turbines
#        - Check that mean (u_) of output planes are the same for (0) and (W)
#        - Check that Turbine 2 gets more TI with WAT
#        - Check that std of output planes of sim (W) is  (k1 u + k2 du)
#        - Compare wake position. Impact on meandering for increasing k.
# 
# - Task 0b: Get ready for meandering frame, and how much does var(U0) polute our signal
#    - Simulations with background turbulence (preferably from AMR wind but can use BTS), (W) and (0)
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





folder ='steady/'
simbaseref = 'FF-NoWAT'
simbasewat = 'FF-WAT'

SIMS =  []
# SIMS +=  ['FF-NoWAT']
# SIMS +=  ['FF-WAT-kd10_kg10']
# SIMS +=  ['FF-WAT-kd00_kg00']
# SIMS +=  ['FF-WAT-kd10_kg00']
# SIMS +=  ['FF-WAT-kd00_kg10']
SIMS +=  ['FF-WAT-kd01_kg01']

iStart = 100
zHub=150
U0= 8
nPlanes=6
D=240
levels=np.linspace(2,10,20)/U0










# --- Derived
outDir = os.path.join(folder, '_processedData')
try:
    os.makedirs(outDir)
except:
    pass

vmin=np.min(levels)
vmax=np.max(levels)

I, iMax, nDigits, sFmt = vtkIndices(folder, simbaseref+'.Low.DisYZ01')
print(I)

I = np.asarray(I)
ISel = I[I>iStart]



# --- Load all and export
# U = np.zeros(len(y), len(z) len(Isel))
for simbase in SIMS:
    ds = loadAllPlanes(folder, simbase, nPlanes, ISel, outDir=outDir, verbose=False)






















# --- Load 2 compare

iPlane= 1
i = 110
y,z,Y,Z,Uref,V,W =  extractFFPlane('YZ', folder, simbaseref, iPlane=iPlane, iTime=i, verbose=True, sFmt=sFmt)
y,z,Y,Z,Uwat,V,W =  extractFFPlane('YZ', folder, simbasewat, iPlane=iPlane, iTime=i, verbose=True, sFmt=sFmt)

iH= np.argmin(np.abs(z-zHub))
print(Uref.shape)
print('len(z)', len(z), 'iH',iH, z[iH], zHub)


# --- Example of plot
fig,axes = plt.subplots(1, 2, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
img = axes[0].contourf(Y, Z, Uref.T/U0, levels=levels, vmin=vmin, vmax=vmax)
img = axes[1].contourf(Y, Z, Uwat.T/U0, levels=levels, vmin=vmin, vmax=vmax)
axes[0].set_aspect('equal')
axes[1].set_aspect('equal')
# ax.set_xlabel('')
# ax.set_ylabel('')
# ax.legend()
cb = colorbar(img, pad=0.1)




















plt.show()


if __name__ == '__main__':
    pass
