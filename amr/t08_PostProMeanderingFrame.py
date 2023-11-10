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
check=False
plot=False


Cases={}
Cases['neutral2'] ={'stability':'neutral','nWT':2, 'planeFileWT':['planesT1129921.nc','planesT2129921.nc']}
# Cases['stable2']  ={'stability':'stable' ,'nWT':2, 'planeFileWT':['planesT1121692.nc','planesT2121692.nc']}
# Cases['unstable2']={'stability':'unstable' ,'nWT':2, 'planeFileWT':['planesT176826.nc','planesT276826.nc']}
# Cases['neutral0'] ={'stability':'neutral','nWT':0, 'planeFileWT':['planesT1129921.nc','planesT2129921.nc']}
# Cases['stable0']  ={'stability':'stable' ,'nWT':0, 'planeFileWT':['planesT1121692.nc','planesT2121692.nc']}

# --- Loop on all cases
for Case in Cases:

    case = Case['name']

    # ---
    U0, dt, D, xyWT1, xyWT2, xPlanes = getSimParamsAMR(case)
    HubHeight  =  150
    refArea = np.pi*D**2/4

    # --- derived parameters
    if Case['nWT']==0:
        LESpath = os.path.join('03-iea15-{}WT/'.format(Case['nWT']), case); 
    else:
        LESpath = os.path.join('02-iea15-{}WT/'.format(Case['nWT']), case); 
    outputPath   = os.path.join(LESpath, 'processedData')
    print('>>> outputPath',outputPath)
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
        # --- Reading trajectories
        dfTrajG = weio.read(os.path.join(outputPath, 'MeanderWT{:d}_TrajectoriesC.csv'.format(iWT))
        dfTrajC = weio.read(os.path.join(outputPath, 'MeanderWT{:d}_TrajectoriesC.csv'.format(iWT))

        planePathWT = os.path.join(LESpath, 'post_processing', Case['planeFileWT'][iWT-1])
        sp = Sampling(planePathWT)

        print('>>> Reading', planePathWT)
        with Timer('Reading...'):
            # ds1 = sp1.read_single_group('pT1', itime=1000, ftime=1001, simCompleted=True, outputPath = os.path.join(outputPath, 'planesT1.zarr'))
            # ds2 = sp2.read_single_group('pT2', itime=1000, ftime=1001, simCompleted=True, outputPath = os.path.join(outputPath, 'planesT2.zarr'))
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
        icz0= np.argmin(np.abs(z-0)) 
        icy0= np.argmin(np.abs(y-0))
        IY0 = np.arange(len(y))
        IZ0 = np.arange(len(z))
        iymax = len(y)
        izmax = len(z)
        dy = y[1]-y[0]
        dz = z[1]-z[0]
        IYBoundary = list(range(195,201))+list(range(0,5))
        ITime = np.arange(itimeMin, np.max(ds.it))
        #ITime = np.arange(itimeMin, itimeMin+3)

        # --- Trajectors
        nP = len(ds.x)
        nt = len(ITime)
        WakeTrajectoriesC = np.zeros((nt, 2*nP))
        WakeTrajectoriesG = np.zeros((nt, 2*nP))
        WakeTrajectoriesM = np.zeros((nt, 2*nP))

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
                yc1, zc1, contour, ax = track_wake_center_plane(Y, Z, U, D, method=method    , shear=shear, plot=plot       )
                yc2, zc2, contour, ax = track_wake_center_plane(Y, Z, U, D, method='Gaussian', shear=shear, plot=plot, ax=ax, col=fColrs(5), mk='d')
                yc1 = np.clip(yc1, -210, 210)
                yc2 = np.clip(yc2, -210, 210)
                zc1 = np.clip(zc1, -160, 50)
                zc2 = np.clip(zc2, -160, 50)
                if it>ITime[0]:
                    # filter simple jumps
                    #yc1 = prev_if_jump(yc1, yc1_prev)
                    #yc2 = prev_if_jump(yc2, yc2_prev)
                    #zc1 = prev_if_jump(zc1, zc1_prev)
                    #zc2 = prev_if_jump(zc2, zc2_prev)
                    # If difference between methods is too large use closest to previous value
                    yc1, yc2 = choose(yc1, yc2, yc_prev)
                    zc1, zc2 = choose(zc1, zc2, zc_prev)
                    yc = (1*yc1+3*yc2)/4
                    zc = (1*zc1+3*zc2)/4
                    # simple moving average
                    zc = (zc + zc_prev)/2
                    yc = (yc + yc_prev)/2
                else:
                    yc = (1*yc1+3*yc2)/4
                    zc = (1*zc1+3*zc2)/4

                yc1_prev = yc1
                yc2_prev = yc2
                zc1_prev = zc1
                zc2_prev = zc2
                yc_prev  = yc
                zc_prev  = zc


                WakeTrajectoriesC[it-itimeMin,2*iP:2*iP+2] = (yc1,zc1)
                WakeTrajectoriesG[it-itimeMin,2*iP:2*iP+2] = (yc2,zc2)
                WakeTrajectoriesM[it-itimeMin,2*iP:2*iP+2] = (yc ,zc )
                #     yc = yc1
                #     zc = zc1
                rz = np.array([0,0])
                rc = np.array([yc,zc])
                if check:
                    print('>>> CENTER BEFORE SHIFT ({:.3f}, {:.3f})'.format(*rc))
                    fig=ax.get_figure()
                    fig.savefig('_figs/Meandering_ip{:}_it{:}_before.png'.format(iP,it))

                # --- Move grid to meandering frame of reference (probably not done in the smartest way)
                icz = np.argmin(np.abs(z-zc))
                icy = np.argmin(np.abs(y-yc))
                # Difference between centers in index space
                diz = icz-icz0 
                diy = icy-icy0
                dg= np.array([diy*dz, diz*dz]) # in true dimension
                d = rc-rz
                y0 = y[icy]
                z0 = z[icz]
                dg = np.array([yc-y0, zc-z0]) # grid error
                # - Shift data
                U2 = np.zeros_like(U)*np.nan
                IY_in_old = IY0+diy
                IZ_in_old = IZ0+diz
                bYOK = np.logical_and(IY_in_old>=0, IY_in_old<iymax)
                bZOK = np.logical_and(IZ_in_old>=0, IZ_in_old<izmax)
                IY_in_new = IY0[bYOK]
                IZ_in_new = IZ0[bZOK]
                U2[np.ix_( IY_in_new, IZ_in_new ) ] = U[np.ix_( IY_in_old[bYOK], IZ_in_old[bZOK] ) ]

                if check:
                    print('>>> CENTER HAS MOVED    ({:.3f}, {:.3f})'.format(*dg))
                    print('>>> CENTER HAS MOVED    ({:.3f}, {:.3f})'.format(*d))
                    print('>>> GRID ERROR IS       ({:.3f}, {:.3f}) < ({:.3f}, {:.3f})'.format(dg[0], dg[1], dy, dz))
                    # --- Recompute wake center, should now be close to zero
                    shear2 = np.zeros_like(shear)*np.nan
                    shear2[IZ_in_new] = shear[IZ_in_old[bZOK]]
                    yc1, zc1, contour, ax = track_wake_center_plane(Y, Z, U2, D, method=method    , shear=shear2, plot=plot       )
                    yc2, zc2, contour, ax = track_wake_center_plane(Y, Z, U2, D, method='Gaussian', shear=shear2, plot=plot , ax=ax, col=fColrs(5), mk='d')
                    ax.set_xlim(np.min(Y.flatten()), np.max(Y.flatten()))
                    ax.set_ylim(np.min(Z.flatten()), np.max(Z.flatten()))
                    yc = (yc1+yc2)/2
                    zc = (zc1+zc2)/2
                    rcn = np.array([yc,zc])
                    rc2 = dg+rcn
                    fig=ax.get_figure()
                    fig.savefig('_figs/Meandering_ip{:}_it{:}_after.png'.format(iP,it))
                    print('>>> CENTER AFTER  SHIFT ({:.3f}, {:.3f})'.format(*rcn)) # Should match grid error
                    if np.abs(rcn[0])>dy or np.abs(rcn[1])>dz:
                        print('>>>>> Problem grid error')
                        import pdb; pdb.set_trace()

                # --- Replace data
                ds['u'].loc[dict(x=iP, it=it)] = U2

        # --- Save data
        cols = ['Time_[s]'] 
        for iP in xPlanes:
            cols += ['y{}'.format(iP) , 'z{}'.format(iP)]
        df = pd.DataFrame(data=np.column_stack((ITime, WakeTrajectoriesC)), columns=cols)
        df.to_csv(os.path.join(outputPath, 'MeanderWT{:d}_TrajectoriesC.csv'.format(iWT)), index=False)
        df = pd.DataFrame(data=np.column_stack((ITime, WakeTrajectoriesG)), columns=cols)
        df.to_csv(os.path.join(outputPath, 'MeanderWT{:d}_TrajectoriesG.csv'.format(iWT)), index=False)
        df = pd.DataFrame(data=np.column_stack((ITime, WakeTrajectoriesM)), columns=cols)
        df.to_csv(os.path.join(outputPath, 'MeanderWT{:d}_TrajectoriesM.csv'.format(iWT)), index=False)

        ds.isel(z=icz0).to_netcdf(os.path.join(outputPath, 'MeanderWT{:d}.nc_small'.format(iWT)))
#         ds2=ds.isel(z=icz0)
#         ds2.isel(x=iP, it=Itime).mean(dim='it', skipna=True)

        del ds
                        

plt.show()
