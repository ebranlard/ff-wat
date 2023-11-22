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


def saveHubHeightPlanes(caseName, Case, outDir='_out'):
    print(f'-----------------------------------------------------------------------')
    print(f'--- saveHubHeightPlanes - Case: {caseName}')
    print(f'-----------------------------------------------------------------------')

    outDir = os.path.join(outDir, 'planes')
    if not os.path.exists(outDir):
        os.makedirs(outDir)


    stability = Case['stability']
    LESpath   = Case['path']

    # ---
    U0, dt, D, xyWT, HubHeight, xPlanes = getSimParamsAMR(stability)
    planePathWT1 = os.path.join(LESpath, 'post_processing', Case['planeFileWT'][0])
    planePathWT2 = os.path.join(LESpath, 'post_processing', Case['planeFileWT'][1])
    outputPath   = os.path.join(outDir, caseName)
    try:
        os.makedirs(outputPath)
    except:
        pass

    # --- Read 
    try:
        print('Reading: ',planePathWT1)
        sp1 = Sampling(planePathWT1)
    except:
        FAIL('Error reading '+planePathWT1)
        raise Exception()
    try:
        print('Reading: ',planePathWT2)
        sp2 = Sampling(planePathWT2)
    except:
        FAIL('Error reading '+planePathWT2)
        raise Exception()
    with Timer('Reading...'):
        ds1 = sp1.read_single_group('pT1', simCompleted=True).rename_dims({'samplingtimestep':'it'})
        ds2 = sp2.read_single_group('pT2', simCompleted=True).rename_dims({'samplingtimestep':'it'})
        # ds1.u # Shape ny x nz x nx x nt
        ds1['y'] = ds1.y-xyWT[1][1]
        ds2['y'] = ds2.y-xyWT[1][1]
        ds1['x'] = ds1.x-xyWT[1][0]
        ds2['x'] = ds2.x-xyWT[1][0] # NOTE: we keep WT1 as reference

    # --- Find relevant indices
    t = ds1.samplingtimestep.values * dt
    z = ds1.z.values
    y = ds1.y.values
    iH = np.argmin(np.abs(ds1.z.values-HubHeight))
    iy0 = np.argmin(np.abs(y-0    ))
    iym = np.argmin(np.abs(y-(-2*D)))
    iyp = np.argmin(np.abs(y-(+2*D)))
    print('z HubHeight {:9.3f}  Closest: {:9.3f} zIndex: {}'.format(HubHeight, ds1.z.values[iH] , iH))
    print('y 0         {:9.3f}  Closest: {:9.3f} yIndex: {}'.format(xyWT[1][1] , ds1.y.values[iy0], iy0))
    # print('y Range    [{:8.3f}, {:8.3f}]'.format(y[iym], y[iyp]))
    # Iy = list(np.arange(iym, iyp+1))
    print('Time range: [{:9.3f},{:9.3f}] '.format(t[0], t[-1]))
    print('xPlanes1:',np.around(np.asarray(ds1.x),2))
    print('xPlanes2:',np.around(np.asarray(ds2.x),2))
    print('xPlanes1:',np.around(np.asarray((ds1.x)/D),2), '[D]')
    print('xPlanes2:',np.around(np.asarray((ds2.x)/D),2), '[D]')

    # --- Save only the data we wont
    print('Saving:',os.path.join(outputPath, 'HubHeightWT1.nc_small'))
    print('Saving:',os.path.join(outputPath, 'HubHeightWT2.nc_small'))
    with Timer('Saving'):
        ds1.isel(z=iH).to_netcdf(os.path.join(outputPath, 'HubHeightWT1.nc_small'))
        ds2.isel(z=iH).to_netcdf(os.path.join(outputPath, 'HubHeightWT2.nc_small'))

if __name__ == '__main__':
    caseNames = [] 
    caseNames += ['neutral2WT']
    caseNames += ['stable2WT'] 
    caseNames += ['neutral0WT']
    caseNames += ['stable0WT'] 
    caseNames += ['unstable2WT'] 
#     caseNames += ['unstable0WT'] 

    Cases = dict((k,v) for k,v in AllCases.items() if k in caseNames)
    print(Cases)

    for caseName, Case in Cases.items():
        try:
            saveHubHeightPlanes(caseName, Case, outDir='_out')
            OK(caseName)
        except:
            FAIL(caseName)

