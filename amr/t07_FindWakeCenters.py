""" 
Find wake center using Gaussian and Contour method, save as CSV
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import weio
from weio.pickle_file import PickleFile
from welib.essentials import *
from helper_functions import *
from windtools.amrwind.post_processing  import Sampling

itimeMin = 500 # time index before which we don't expect wakes to be present at all locations.
plot=False


Cases={}
Cases['neutral2'] ={'stability':'neutral','nWT':2, 'planeFileWT':['planesT1129921.nc','planesT2129921.nc']}
# Cases['stable2']  ={'stability':'stable' ,'nWT':2, 'planeFileWT':['planesT1121692.nc','planesT2121692.nc']}
# Cases['unstable2']={'stability':'unstable' ,'nWT':2, 'planeFileWT':['planesT176826.nc','planesT276826.nc']}
# Cases['neutral0'] ={'stability':'neutral','nWT':0, 'planeFileWT':['planesT1129921.nc','planesT2129921.nc']}
# Cases['stable0']  ={'stability':'stable' ,'nWT':0, 'planeFileWT':['planesT1121692.nc','planesT2121692.nc']}

# --- Loop on all cases
for casename, Case in Cases.items():

    case = Case['stability']

    # ---
    U0, dt, D, xyWT1, xyWT2, xPlanes = getSimParamsAMR(case)
    HubHeight  =  150
    refArea = np.pi*D**2/4
    LESpath    = os.path.join('02-iea15-{}WT/'.format(Case['nWT']), case);
    outputPath = os.path.join(LESpath, 'processedData')
    try:
        os.makedirs(outputPath)
    except:
        pass

    for iWT in [1, 2]:
        # --- Read 
        if iWT==1:
            xyWT=xyWT1
        else:
            xyWT=xyWT2
        planePathWT = os.path.join(LESpath, 'post_processing', Case['planeFileWT'][iWT-1])
        sp = Sampling(planePathWT)
        print('>>> Reading', planePathWT)
        with Timer('Reading...'):
            ds = sp.read_single_group('pT{}'.format(iWT), simCompleted=True).rename_dims({'samplingtimestep':'it'})
            ds['z'] = ds.z -(np.max(ds.z)+np.min(ds.z))/2
            ds['y'] = ds.y-xyWT[1]
            ds['x'] = np.around((ds.x-xyWT[0])/D)
        print(ds)


        # --- Common variables useful for grid
        z = ds.z.values
        y = ds.y.values
        Z, Y = np.meshgrid(ds.z.values, ds.y.values)

        # Index coordinates of domain center (we take this as a reference
        IYBoundary = list(range(195,201))+list(range(0,5))
        ITime = np.arange(itimeMin, np.max(ds.it))
        #ITime = np.arange(itimeMin, itimeMin+3)

        # --- Trajectors
        nP = len(ds.x)
        nt = len(ITime)
        WakeTrajectoriesC = np.zeros((nt, 2*nP))
        WakeTrajectoriesG = np.zeros((nt, 2*nP))

        def choose(y1,y2,yprev):
            if np.abs(y1-y2)>50:
                if np.abs(y1-yprev)<50:
                    y2=y1
                else:
                    y1=y2
            return y1, y2
        def prev_if_jump(y,yprev):
            if np.abs(y-yprev)>50:
                y=yprev
            return y

        # ---  Loop on planes
        for iP, x in enumerate(xPlanes):
            # Mean shear over full simulation 
            shear = ds.isel(x=iP, y=IYBoundary).mean(dim=['it','y']).u.values
            #for it in np.arange(300,1500,500):
            # for it in [1300]:
            for it in range(0,itimeMin):
                ds['u'].loc[dict(x=iP, it=it)] = np.nan

            for it in ITime:
                if np.mod(it,100)==0:
                    print('iWT',iWT, 'iP',iP, 'it',it)
                # --- Find wake center
                cs= ds.isel(x=iP, it=it).u.values
                U = ds.isel(x=iP, it=it).u.values
                yc1, zc1, contour, ax = track_wake_center_plane(Y, Z, U, D, method='ConstantArea', shear=shear, plot=plot       )
                yc2, zc2, contour, ax = track_wake_center_plane(Y, Z, U, D, method='Gaussian'    , shear=shear, plot=plot, ax=ax, col=fColrs(5), mk='d')

                WakeTrajectoriesC[it-itimeMin,2*iP:2*iP+2] = (yc1,zc1)
                WakeTrajectoriesG[it-itimeMin,2*iP:2*iP+2] = (yc2,zc2)

        # --- Save data
        cols = ['Time_[s]'] 
        for iP in xPlanes:
            cols += ['y{}'.format(iP) , 'z{}'.format(iP)]
        df = pd.DataFrame(data=np.column_stack((ITime, WakeTrajectoriesC)), columns=cols)
        df.to_csv(os.path.join(outputPath, 'MeanderWT{:d}_TrajectoriesC.csv'.format(iWT)), index=False)
        df = pd.DataFrame(data=np.column_stack((ITime, WakeTrajectoriesG)), columns=cols)
        df.to_csv(os.path.join(outputPath, 'MeanderWT{:d}_TrajectoriesG.csv'.format(iWT)), index=False)
                        

plt.show()
