import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatch
# External
from windtools.amrwind.post_processing  import Sampling
# Local 
import welib.weio as weio
from helper_functions import *

# --- Loop on all cases
def saveMeanderingFrameLines(caseName, Case, plot=False, nFigsMax=3, outDir='_out', figDir='_fig', iTimeMin=500, returnFig=False):
    print(f'-----------------------------------------------------------------------')
    print(f'--- saveMeanderingFrameLines - Case: {caseName}')
    print(f'-----------------------------------------------------------------------')
    nWT  = Case['nWT']
    # --- derived parameters
    outDirHori = os.path.join(outDir  , 'lines')
    outDirTraj   = os.path.join(outDir, 'trajectories')
    figDirTrack  = os.path.join(figDir, '_figs_tracking')
    if not os.path.exists(outDirHori):
        os.makedirs(outDirHori)
    if not os.path.exists(figDirTrack):
        os.makedirs(figDirTrack)

    if nWT in [1,2]:
        caseNameForTraj = caseName
        for iWT in [1, 2]:
            fig = saveMeanderingFrameLines2(iWT, caseNameForTraj, caseName, Case, plot=plot, nFigsMax=nFigsMax, outDir=outDir, figDir=figDir, iTimeMin=iTimeMin, returnFig=returnFig)
            if returnFig:
                return fig
    else:
        # Using meandering frame from 2 WT case
        caseNameForTraj = caseName.replace('0WT','2WT')
        print('>>> Case 0WT using trajectories from:', caseNameForTraj)
        for iWT in [1, 2]:
            saveMeanderingFrameLines2(iWT, caseNameForTraj, caseName, Case, plot=plot, nFigsMax=nFigsMax, outDir=outDir, figDir=figDir, iTimeMin=iTimeMin, returnFig=returnFig)
        # Using meandering frame from 1 WT case
        caseNameForTraj = caseName.replace('0WT','1WT')
        print('>>> Case 0WT using trajectories from:', caseNameForTraj)
        for iWT in [1, 2]:
            saveMeanderingFrameLines2(iWT, caseNameForTraj, caseName, Case, plot=plot, nFigsMax=nFigsMax, outDir=outDir, figDir=figDir, iTimeMin=iTimeMin, returnFig=returnFig)

def saveMeanderingFrameLines2(iWT, caseNameForTraj, caseName, Case, plot=False, nFigsMax=3, outDir='_out', figDir='_fig', iTimeMin=500, returnFig=False):
    U0, dt, D, xyWT, HubHeight, xPlanes = getSimParamsAMR(Case['stability'])
    nWT  = Case['nWT']
    outDirHori = os.path.join(outDir  , 'lines')
    outDirTraj   = os.path.join(outDir, 'trajectories')
    figDirTrack  = os.path.join(figDir, '_figs_tracking')

    # --- Reading trajectories
    trajFile = os.path.join(outDirTraj, '{}_MeanderWT{:d}_TrajectoriesM.csv'.format(caseNameForTraj, iWT))
    print('Reading: ',trajFile)
    dfTrajM = weio.read(trajFile).toDataFrame()
    dfTrajM.index = dfTrajM['Time_[s]'].astype(int)

    # --- Read 
    planeBase = os.path.join(Case['path'], 'post_processing', 'planesT{}'.format(iWT))
    with Timer('Reading...'):
        ds = readPlanes(planeBase, group='pT{}'.format(iWT))
        ds['z'] = ds.z -(np.max(ds.z)+np.min(ds.z))/2
        ds['y'] = ds.y-xyWT[iWT][1]
        ds['x'] = np.around((ds.x-xyWT[iWT][0])/D)
    #print(ds)
    print('ITime :',ds.it.values[0], ds.it.values[-1])

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
    IYBoundary = list(range(195,201))+list(range(0,5)) # TODO TODO TODO Not Generic
    ITime = np.arange(iTimeMin, np.max(ds.it))
    if len(ITime)>len(dfTrajM):
        print('[WARN] THIS SIMULATION HAS MORE ITIME THAN TRAJECTORY', len(ITime), len(dfTrajM))
        ITime = ITime[:len(dfTrajM)]
    print('ITime :',ITime[0], ITime[-1])
    # ---  Loop on planes

    for iP, x in enumerate(xPlanes):
        with Timer('Shift to meandering frame, WT:{} Plane:{}'.format(iWT, iP)):
            nFigs=0
            # Mean shear over full simulation 
            shear = ds.isel(x=iP, y=IYBoundary).mean(dim=['it','y']).u.values
            #for it in np.arange(300,1500,500):
            # for it in [1300]:
            for it in range(0,iTimeMin):
                ds['u'].loc[dict(x=iP, it=it)] = np.nan

            # Wake centers
            Yc = dfTrajM['y{}'.format(iP)]
            Zc = dfTrajM['z{}'.format(iP)]
            for it in ITime:
                #if np.mod(it,100)==0:
                #    print('iWT',iWT, 'iP',iP, 'it',it)
                # --- Find wake center
                U = ds.isel(x=iP, it=it).u.values.copy()
                yc, zc = Yc[it], Zc[it]

                # --- Move grid to meandering frame of reference (probably not done in the smartest way)
                icz = np.argmin(np.abs(z-zc))
                icy = np.argmin(np.abs(y-yc))
                # Difference between centers in index space
                diz = icz-icz0 
                diy = icy-icy0
                # - Shift data
                U2 = np.zeros_like(U)*np.nan
                IY_in_old = IY0+diy
                IZ_in_old = IZ0+diz
                bYOK = np.logical_and(IY_in_old>=0, IY_in_old<iymax)
                bZOK = np.logical_and(IZ_in_old>=0, IZ_in_old<izmax)
                IY_in_new = IY0[bYOK]
                IZ_in_new = IZ0[bZOK]
                U2[np.ix_( IY_in_new, IZ_in_new ) ] = U[np.ix_( IY_in_old[bYOK], IZ_in_old[bZOK] ) ]

                # --- Replace data
                ds['u'].loc[dict(x=iP, it=it)] = U2

                # --- Plot
                if plot and nFigs<nFigsMax and np.mod(it,100)in [0,1,2,3,4,5]:
                    zOff= 250 # TODO TODO
                    u1    = 0
                    u2    = 10
                    mk    = 'o'
                    col   = 'w'
                    nFigs = nFigs+1
                    Z0 = Z.copy()+zOff

                    # --- Redo Wake center analysis for plot only
                    if nWT>0:
                        yc1, zc1, contour1, ax = track_wake_center_plane(Y, Z0, U, D, method='ConstantArea', shear=shear, plot=False       )
                        yc2, zc2, contour2, ax = track_wake_center_plane(Y, Z0, U, D, method='Gaussian'    , shear=shear, plot=False, ax=ax, col=fColrs(5), mk='d')
                        shear2 = np.zeros_like(shear)*np.nan
                        shear2[IZ_in_new] = shear[IZ_in_old[bZOK]]
                        yc1m, zc1m, contour1m, ax = track_wake_center_plane(Y, Z, U2, D, method='ConstantArea', shear=shear2, plot=False       )
                        yc2m, zc2m, contour2m, ax = track_wake_center_plane(Y, Z, U2, D, method='Gaussian'    , shear=shear2, plot=False , ax=ax, col=fColrs(5), mk='d')

                    if nWT==0:
                        figname ='{}__IN{}__{}D__t{}'.format(caseName, caseNameForTraj, ((iWT-1)*len(xPlanes))+iP, it)
                    else:
                        figname ='{}__{}D__t{}'.format(caseName, ((iWT-1)*len(xPlanes))+iP, it)
                    #print('>>> PLOTTING', nFigs, figname)
                    fig,axes = plt.subplots(2, 1, sharex=True, figsize=(6.4,6.8)) # (6.4,4.8)
                    fig.subplots_adjust(left=0.12, right=0.95, top=0.90, bottom=0.11, hspace=0.30, wspace=0.20)
                    clevels = np.linspace(u1, u2, 100)

                    # --- Plot in Inertial frame
                    ax=axes[0]
                    cf = ax.contourf(Y, Z0, U, clevels, cmap='viridis', extend='both')
                    #cb = fig.colorbar(cf, ticks=np.linspace(u1, u2, 11))
                    #cb.set_label(label=r'$U$ [m/s]',fontsize=14)
                    #cb.ax.tick_params(labelsize=12)
                    ax.set_xlim(np.min(Y.flatten()), np.max(Y.flatten()))
                    ax.set_ylim(np.min(Z0.flatten()), np.max(Z0.flatten()))
                    ax.axis('scaled')
                    ax.set_ylabel(r'$z$ [m]', fontsize=14)
                    Colrs = [fColrs(1), fColrs(4)]

                    if nWT>0:
                        ax.plot(yc1, zc1, marker='o' , c=Colrs[0], ms=8 , label = 'ConstantArea', alpha=0.6, markeredgewidth=1, markeredgecolor='k')
                        ax.plot(yc2, zc2, marker='d' , c=Colrs[1], ms=8 , label = 'Gaussian'    , alpha=0.6, markeredgewidth=1, markeredgecolor='k')
                    ax.plot(yc,  zc+zOff, marker='x'  , c='k'      , ms=10, label = 'Clean', alpha=1.0, markeredgewidth=1, markeredgecolor='k')
                    if nWT>0:
                        try:
                            ax.add_patch(mpatch.PathPatch(mpath.Path(contour1), lw=1,ls='-', facecolor='none', edgecolor=Colrs[0]))
                            ax.add_patch(mpatch.PathPatch(mpath.Path(contour2), lw=1,ls='-', facecolor='none', edgecolor=Colrs[1]))
                        except:
                            print('>>> FAIL')
                            pass
                    ax.legend()

                    # --- Plot in Inertial frame
                    ax=axes[1]
                    cf = ax.contourf(Y, Z, U2, clevels, cmap='viridis', extend='both')
                    #cb = fig.colorbar(cf, ticks=np.linspace(u1, u2, 11))
                    #cb.set_label(label=r'$U$ [m/s]',fontsize=14)
                    #cb.ax.tick_params(labelsize=12)
                    ax.set_xlim(np.min(Y.flatten()), np.max(Y.flatten()))
                    ax.set_ylim(np.min(Z.flatten()), np.max(Z.flatten()))
                    ax.axis('scaled')
                    ax.set_xlabel(r'$y$ [m]', fontsize=14)
                    ax.set_ylabel(r'$z$ [m]', fontsize=14)
                    if nWT>0:
                        ax.plot(yc1m, zc1m, marker='o' , c=Colrs[0]      , ms=8 , label = 'ConstantArea', alpha=0.6, markeredgewidth=1, markeredgecolor='k')
                        ax.plot(yc2m, zc2m, marker='d' , c=Colrs[1], ms=8 , label = 'Gaussian'    , alpha=0.6, markeredgewidth=1, markeredgecolor='k')
                    ax.plot(0   , 0   , marker='x' , c='k'      , ms=10, label = 'Clean'       , alpha=1.0, markeredgewidth=1, markeredgecolor='k')
                    if nWT>0:
                        try:
                            ax.add_patch(mpatch.PathPatch(mpath.Path(contour1m), lw=1,ls='-', facecolor='none', edgecolor=Colrs[0]))
                            ax.add_patch(mpatch.PathPatch(mpath.Path(contour2m), lw=1,ls='-', facecolor='none', edgecolor=Colrs[1]))
                        except:
                            print('>>> FAIL')
                            pass
                    for ax in np.array(axes).flatten():
                        ax.tick_params(direction='in', top=True, right=True, labelright=False, labeltop=False, which='both')
                    fig.suptitle(figname)
                    fig.savefig(os.path.join(figDirTrack,figname+'.png'))
                    if returnFig:
                        return fig
                    plt.close(fig)
                    # --- End of Plot
                    
    # --- Save data at mean wake center height
    if nWT==0:
        outFile = os.path.join(outDirHori, '{}_IN{}_Meander_WT{:d}.nc_small'.format(caseName, caseNameForTraj, iWT))
    else:
        outFile = os.path.join(outDirHori, '{}_Meander_WT{:d}.nc_small'.format(caseName, iWT))
    print('Writing: ', outFile)
    ds.isel(z=icz0).to_netcdf(outFile)
                        
if __name__ == '__main__':

    caseNames =[]
#     caseNames += list(AllCases.keys())

#     caseNames += ['neutral2WT']
#     caseNames += ['stable2WT'] 
    caseNames += ['unstable2WT'] 
#     caseNames += ['neutral1WT']
#     caseNames += ['stable1WT'] 
#     caseNames += ['unstable1WT'] 
#     caseNames += ['neutral0WT']
#     caseNames += ['stable0WT'] 
#     caseNames += ['unstable0WT'] 
    Cases = dict((k,v) for k,v in AllCases.items() if k in caseNames)

    for caseName, Case in Cases.items():
        with FailSafe(caseName, True):
            saveMeanderingFrameLines(caseName, Case, plot=True)

plt.show()
