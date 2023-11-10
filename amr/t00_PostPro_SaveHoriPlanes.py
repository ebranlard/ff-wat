import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# External
from windtools.amrwind.post_processing  import Sampling
# Local 
import weio
from welib.essentials import *
# Super local
from helper_functions import *

# --- Turbine/simulation parameters
# case = 'neutral';nWT=2; planeFileWT1 = 'planesT1129921.nc'; planeFileWT2 = 'planesT2129921.nc';
# case = 'stable'; nWT=2; U0=8.1; planeFileWT1 = 'planesT1121692.nc'; planeFileWT2 = 'planesT2121692.nc';
# case = 'neutral';nWT=0; U0=8.5; planeFileWT1 = 'planesT1129921.nc'; planeFileWT2 = 'planesT2129921.nc';
# case = 'stable'; nWT=0; U0=8.1; planeFileWT1 = 'planesT1121692.nc'; planeFileWT2 = 'planesT2121692.nc';

Cases=[]
# Cases.append({'name':'neutral','nWT':2, 'planeFileWT':['planesT1129921.nc','planesT2129921.nc']})
# Cases.append({'name':'stable' ,'nWT':2, 'planeFileWT':['planesT1121692.nc','planesT2121692.nc']})
# Cases.append({'name':'neutral','nWT':0, 'planeFileWT':['planesT1129921.nc','planesT2129921.nc']})
# Cases.append({'name':'stable' ,'nWT':0, 'planeFileWT':['planesT1121692.nc','planesT2121692.nc']})
Cases.append({'name':'unstable' ,'nWT':2, 'planeFileWT':['planesT176826.nc','planesT276826.nc']})
                                                          



# --- Loop on all cases
for Case in Cases:

    case = Case['name']

    # ---
    U0, dt, D, xyWT1, xyWT2, xPlanes = getSimParamsAMR(case)
    HubHeight  =  150

    # --- derived parameters
    if Case['nWT']==0:
        LESpath = os.path.join('03-iea15-{}WT/'.format(Case['nWT']), case); 
    else:
        LESpath = os.path.join('02-iea15-{}WT/'.format(Case['nWT']), case); 
    # print('>>>> HUBHEIGHT', HubHeight)
    outputPath   = os.path.join(LESpath, 'processedData')
    planePathWT1 = os.path.join(LESpath, 'post_processing', Case['planeFileWT'][0])
    planePathWT2 = os.path.join(LESpath, 'post_processing', Case['planeFileWT'][1])
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
        ds1['y'] = ds1.y-xyWT1[1]
        ds2['y'] = ds2.y-xyWT1[1]
        ds1['x'] = ds1.x-xyWT1[0]
        ds2['x'] = ds2.x-xyWT1[0] # TODO WT2???

    # --- Find relevant indices
    t = ds1.it.values * dt
    z = ds1.z.values
    y = ds1.y.values
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
    print('xPlanes1:',np.around(np.asarray((ds1.x)/D),2), '[D]')
    print('xPlanes2:',np.around(np.asarray((ds2.x)/D),2), '[D]')
    nPlanes = len(xPlanes)


    # --- Save only the data we wont
    with Timer('Saving'):
        ds1.isel(z=iH).to_netcdf(os.path.join(outputPath, 'HubHeightWT1.nc_small'))
        ds2.isel(z=iH).to_netcdf(os.path.join(outputPath, 'HubHeightWT2.nc_small'))



