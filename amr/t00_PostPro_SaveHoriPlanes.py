import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Super local
from helper_functions import *


def saveHubHeightLines(caseName, Case, outDir='_out'):
    print(f'-----------------------------------------------------------------------')
    print(f'--- saveHubHeightLines - Case: {caseName}')
    print(f'-----------------------------------------------------------------------')

    # --- derived parameters
    outDirHori = os.path.join(outDir, 'lines')
    if not os.path.exists(outDirHori):
        os.makedirs(outDirHori)
    U0, dt, D, xyWT, HubHeight, xPlanes = getSimParamsAMR(Case['stability'])

    for iWT in [1, 2]:
        # --- Read 
        planeBase = os.path.join(Case['path'], 'post_processing', 'planesT{}'.format(iWT))

        with Timer('Reading...'):
            ds = readPlanes(planeBase, group='pT{}'.format(iWT))
            # ds.u # Shape ny x nz x nx x nt 
            ds['x'] = ds.x-xyWT[1][0] # NOTE we use wxWT[1] as a reference for both j
            ds['y'] = ds.y-xyWT[1][1]

        # --- Find relevant indices
        t = ds.samplingtimestep.values * dt
        z = ds.z.values
        y = ds.y.values
        iH = np.argmin(np.abs(ds.z.values-HubHeight))
        iy0 = np.argmin(np.abs(y-0    ))
        print('y Range: [{:9.3f},{:9.3f}]'.format(np.min(y), np.max(y)))
        print('z Range: [{:9.3f},{:9.3f}]'.format(np.min(z), np.max(z)))
        print('i Range: [{:9d},{:9d}] '.format(ds.it.values[0], ds.it.values[-1]))
        print('t Range: [{:9.3f},{:9.3f}] '.format(t[0], t[-1]))
        print('xPlanes:',np.around(np.asarray((ds.x)/D),2), '[D]')
        print('z HubHeight {:9.3f}  Closest: {:9.3f} zIndex: {}'.format(HubHeight, ds.z.values[iH] , iH))
        print('y 0         {:9.3f}  Closest: {:9.3f} yIndex: {}'.format(xyWT[1][1] , ds.y.values[iy0], iy0))
        #print('xPlanes{}:'.format(iWT),np.around(np.asarray(ds.x),2))

        # --- Save data at hub height
        outFile = os.path.join(outDirHori, '{}_Inertial_WT{:d}.nc_small'.format(caseName, iWT))
        print('Writing: ', outFile)
        ds.isel(z=iH).to_netcdf(outFile)
        del ds

if __name__ == '__main__':
    caseNames = [] 
    caseNames += ['neutral2WT']
    caseNames += ['stable2WT'] 
    caseNames += ['neutral0WT']
    caseNames += ['stable0WT'] 
    caseNames += ['unstable2WT'] 
#     caseNames += ['unstable0WT'] 

    Cases = dict((k,v) for k,v in AllCases.items() if k in caseNames)

    for caseName, Case in Cases.items():
        with FailSafe(caseName, True):
            saveHubHeightLines(caseName, Case, outDir='_out')

