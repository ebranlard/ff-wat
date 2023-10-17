import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import weio
from welib.essentials import *

# import os, sys, shutil
# import xarray as xr
# import pickle
# from scipy.stats import gaussian_kde
from windtools.amrwind.post_processing  import Sampling
# from windtools.amrwind.post_processing  import ABLStatistics, Sampling, addDatetime
# from samwich.dataloaders import XarrayData
# from samwich.waketrackers import track, WakeTracker
# 
# from pyFAST.input_output import TurbSimFile, FASTOutputFile
# sys.path.append(os.path.abspath('/home/rthedin/utilities/'))
# from helper_fastfarm import readFFPlanes, readTurbineOutput, readTurbineOutputPar, readFFPlanesPar
# from helper import interpolate_to_heights, addLabels, calc_QOIs



# --- Turbine/simulation parameters
# case = 'neutral';nWT=2; U0=8.5; planeFileWT1 = 'planesT1129921.nc'; planeFileWT2 = 'planesT2129921.nc';
# case = 'stable'; nWT=2; U0=8.1; planeFileWT1 = 'planesT1121692.nc'; planeFileWT2 = 'planesT2121692.nc';
case = 'neutral';nWT=0; U0=8.5; planeFileWT1 = 'planesT1129921.nc'; planeFileWT2 = 'planesT2129921.nc';
case = 'stable'; nWT=0; U0=8.1; planeFileWT1 = 'planesT1121692.nc'; planeFileWT2 = 'planesT2121692.nc';


fixed_dt =   0.025
output_frequency = 40  # every 1 s
dt            = fixed_dt*output_frequency
R             = 120.97
D             = 2*R
OverHang      = -12.097571763912535
xyWT1         = (706.25 , 641.25)   # Actuator.T1.base_position = 713.28 641.25 0    # hub at 701.25, 641.25, 0, considering overhang of 12.0313
xyWT2         = (2386.25, 641.25)   # Actuator.T2.base_position = 2393.28 641.25 0.  # hub at 2381.25, 641.25, 0, considering overhang of 12.0313
D             = 240
ShftTilt      = -6.0 
Twr2Shft      = 4.349459414248071
TowerHt       = 144.386          


# --- derived parameters
if nWT==0:
    LESpath = os.path.join('03-iea15-{}WT/'.format(nWT), case); 
else:
    LESpath = os.path.join('02-iea15-{}WT/'.format(nWT), case); 
HubHeight  =  TowerHt + Twr2Shft + OverHang*np.sin(ShftTilt*np.pi/180)
outputPath   = os.path.join(LESpath, 'processedData')
planePathWT1 = os.path.join(LESpath, 'post_processing', planeFileWT1)
planePathWT2 = os.path.join(LESpath, 'post_processing', planeFileWT2)
print('>>> outputPath',outputPath)
try:
    os.makedirs(outputPath)
except:
    pass

# --- Read 
sp1 = Sampling(planePathWT1)
sp2 = Sampling(planePathWT2)
with Timer('Reading...'):
    # ds1 = sp1.read_single_group('pT1', itime=1000, ftime=1001, simCompleted=True, outputPath = os.path.join(outputPath, 'planesT1.zarr'))
    # ds2 = sp2.read_single_group('pT2', itime=1000, ftime=1001, simCompleted=True, outputPath = os.path.join(outputPath, 'planesT2.zarr'))
    ds1 = sp1.read_single_group('pT1', simCompleted=True).rename_dims({'samplingtimestep':'it'})
    ds2 = sp2.read_single_group('pT2', simCompleted=True).rename_dims({'samplingtimestep':'it'})
    # ds1.u # Shape ny x nz x nx x nt

# --- Find relevant indices
t = ds1.it.values * dt
z = ds1.z.values
y = ds1.y.values-xyWT1[1]
iH = np.argmin(np.abs(ds1.z.values-HubHeight))
iy0 = np.argmin(np.abs(y-0    ))
iym = np.argmin(np.abs(y-(-2*D)))
iyp = np.argmin(np.abs(y-(+2*D)))
print('z HubHeight {:9.3f}  Closest: {:9.3f} zIndex: {}'.format(HubHeight, ds1.z.values[iH] , iH))
print('y 0         {:9.3f}  Closest: {:9.3f} yIndex: {}'.format(xyWT1[1] , ds1.y.values[iy0], iy0))
# print('y Range    [{:8.3f}, {:8.3f}]'.format(y[iym], y[iyp]))
# Iy = list(np.arange(iym, iyp+1))
Itime = np.where(t>600)[0]
t1=t[Itime[0]]
t2=t[Itime[-1]]
print('Time range: [{:9.3f},{:9.3f}] T: {:9.3f}  Indices: {} {}'.format(t1, t2, t2-t1, Itime[0], Itime[-1]))

# planesT1.pT1.origin       = 706.25  141.25   11.25
# planesT1.pT1.offsets      = 0 240 480 720 960 1200 1440
# planesT2.pT2.origin       = 2386.25   141.25   11.25
# planesT2.pT2.offsets      = 0 240 480 720 960 1200 1440
xPlanes = [0,1,2,3,4,5,6]
print('xPlanes1:',np.around(np.asarray(ds1.x),2))
print('xPlanes2:',np.around(np.asarray(ds2.x),2))
print('xPlanes1:',np.around(np.asarray((ds1.x-xyWT1[0])/D),2), '[D]')
print('xPlanes2:',np.around(np.asarray((ds2.x-xyWT2[0])/D),2), '[D]')
nPlanes = len(xPlanes)


# --- Save only the data we wont
with Timer('Saving'):
    ds1.isel(z=iH).to_netcdf(os.path.join(outputPath, 'HubHeightWT1.nc'))
    ds2.isel(z=iH).to_netcdf(os.path.join(outputPath, 'HubHeightWT2.nc'))



