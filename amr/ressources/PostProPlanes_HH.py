import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import weio
from welib.essentials import *

import xarray

dt=1
xyWT1         = (706.25 , 641.25)   # Actuator.T1.base_position = 713.28 641.25 0    # hub at 701.25, 641.25, 0, considering overhang of 12.0313
xyWT2         = (2386.25, 641.25)   # Actuator.T2.base_position = 2393.28 641.25 0.  # hub at 2381.25, 641.25, 0, considering overhang of 12.0313
D             = 240



case = 'neutral'; U0=8.5 ; 
# case = 'stable';  U0=8.1;
LESpath = os.path.join('../02-iea15-2WT/', case)
datapath = os.path.join('../02-iea15-2WT/', case, 'processedData')

ds1 = xarray.open_dataset(os.path.join(datapath,'HubHeightWT1.nc'))
ds2 = xarray.open_dataset(os.path.join(datapath,'HubHeightWT2.nc'))
print('it',ds1.it.values)


y = ds1.y.values-xyWT1[1]
t = ds1.it.values * dt
Itime = np.where(t>600)[0]
t1=t[Itime[0]]
t2=t[Itime[-1]]

xPlanes = [0,1,2,3,4,5,6]


# --- Plot average wake deficit
fig,axes = plt.subplots(1, 2, sharey=True, figsize=(12.4,5.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.92, top=0.91, bottom=0.11, hspace=0.20, wspace=0.20)
for iP, x in enumerate(xPlanes):
    axes[0].plot(ds1.isel(x=iP, it=Itime).mean(dim='it').u, y/D, label='x= {}D'.format(int(x)))
    axes[1].plot(ds2.isel(x=iP, it=Itime).mean(dim='it').u, y/D, label='x= {}D'.format(int(x)))
    axes[0].set_xlabel('u [m/s]')
    axes[0].set_ylabel('y/D [-]')
axes[0].legend()
axes[0].set_title('WT 1')
axes[1].set_title('WT 2')
axes[1].legend()
fig.tight_layout()


plt.show()
