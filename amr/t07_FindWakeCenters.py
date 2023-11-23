""" 
Find wake center using Gaussian and Contour method, save as CSV
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from helper_functions import *
from windtools.amrwind.post_processing  import Sampling


def findWakeCenters(caseName, Case, plot=False, iTimeMin=500, outDir='_out'):
    print(f'-----------------------------------------------------------------------')
    print(f'--- findWakeCenters - Case: {caseName}')
    print(f'-----------------------------------------------------------------------')
    outDirTraj = os.path.join(outDir, 'trajectories')
    if not os.path.exists(outDirTraj):
        os.makedirs(outDirTraj)

    if Case['nWT']==0:
        WARN('Skipping case (no wake to be sought for if no turbine present)')
        return

    # ---
    U0, dt, D, xyWT, HubHeight, xPlanes = getSimParamsAMR(Case['stability'])

    for iWT in [1, 2]:
        # --- Read 
        planeBase = os.path.join(Case['path'], 'post_processing', 'planesT{}'.format(iWT))
        with Timer('Reading...'):
            ds = readPlanes(planeBase, Case['planeTimes'], group='pT{}'.format(iWT))
            ds['z'] = ds.z -(np.max(ds.z)+np.min(ds.z))/2
            ds['y'] = ds.y-xyWT[iWT][1]
            ds['x'] = np.around((ds.x-xyWT[iWT][0])/D)
        #print(ds)

        # --- Common variables useful for grid
        z = ds.z.values
        y = ds.y.values
        Z, Y = np.meshgrid(ds.z.values, ds.y.values)

        # Index coordinates of domain center (we take this as a reference
        IYBoundary = list(range(195,201))+list(range(0,5))
        ITime = np.arange(iTimeMin, np.max(ds.it))
        print('ITime : ', iTimeMin, ITime[-1]) 
        #ITime = np.arange(itimeMin, itimeMin+3)

        # --- Trajectors
        nP = len(ds.x)
        nt = len(ITime)
        WakeTrajectoriesC = np.zeros((nt, 2*nP))
        WakeTrajectoriesG = np.zeros((nt, 2*nP))

        # ---  Loop on planes
        print('Wake tracking...')
        for iP, x in enumerate(xPlanes):
            with Timer('Wake tracking, WT:{} Plane:{}'.format(iWT, iP)):
                # Mean shear over full simulation 
                shear = ds.isel(x=iP, y=IYBoundary).mean(dim=['it','y']).u.values
                #for it in np.arange(300,1500,500):
                # for it in [1300]:
                for it in range(0,iTimeMin):
                    ds['u'].loc[dict(x=iP, it=it)] = np.nan

                for it in ITime:
                    #if np.mod(it,100)==0:
                    #    print('iWT',iWT, 'iP',iP, 'it {}/{}'.format(it, ITime[-1]))
                    # --- Find wake center
                    cs= ds.isel(x=iP, it=it).u.values
                    U = ds.isel(x=iP, it=it).u.values
                    yc1, zc1, contour, ax = track_wake_center_plane(Y, Z, U, D, method='ConstantArea', shear=shear, plot=plot       )
                    yc2, zc2, contour, ax = track_wake_center_plane(Y, Z, U, D, method='Gaussian'    , shear=shear, plot=plot, ax=ax, col=fColrs(5), mk='d')

                    WakeTrajectoriesC[it-iTimeMin,2*iP:2*iP+2] = (yc1,zc1)
                    WakeTrajectoriesG[it-iTimeMin,2*iP:2*iP+2] = (yc2,zc2)

        # --- Save data
        cols = ['Time_[s]'] 
        for iP in xPlanes:
            cols += ['y{}'.format(iP) , 'z{}'.format(iP)]
        trajCfile = os.path.join(outDirTraj, '{}_MeanderWT{:d}_TrajectoriesC.csv'.format(caseName, iWT))
        trajGfile = os.path.join(outDirTraj, '{}_MeanderWT{:d}_TrajectoriesG.csv'.format(caseName, iWT))
        print('Writing: ', trajCfile)
        pd.DataFrame(data=np.column_stack((ITime, WakeTrajectoriesC)), columns=cols).to_csv(trajCfile, index=False)
        pd.DataFrame(data=np.column_stack((ITime, WakeTrajectoriesG)), columns=cols).to_csv(trajGfile, index=False)
                        

if __name__ == '__main__':
    caseNames =[]
#     caseNames += list(AllCases.keys())
    caseNames += ['neutral2WT']
    caseNames += ['stable2WT'] 
    caseNames += ['unstable2WT'] 
#     caseNames += ['neutral1WT']
#     caseNames += ['stable1WT'] 
#     caseNames += ['unstable1WT'] 
    Cases = dict((k,v) for k,v in AllCases.items() if k in caseNames)

    # --- Loop on all cases
    for caseName, Case in Cases.items():
        dfC, dfG = findWakeCenters(caseName, Case, plot=False)

    plt.show()
