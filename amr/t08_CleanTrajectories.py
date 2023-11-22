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
from welib.tools.signal_analysis import moving_average
from welib.tools.signal_analysis import moving_average_conv


def processTrajectories(outputPath, iWT, xPlanes):
    """ Read Gaussian and Contour trajectories, perform moving averaging on them """
    yCols = ['y{}'.format(iP) for iP in xPlanes]
    zCols = ['z{}'.format(iP) for iP in xPlanes]
    Cols = yCols+zCols
    # --- Reading trajectories
    dfG = weio.read(os.path.join(outputPath, 'MeanderWT{:d}_TrajectoriesC.csv'.format(iWT))).toDataFrame()
    dfC = weio.read(os.path.join(outputPath, 'MeanderWT{:d}_TrajectoriesG.csv'.format(iWT))).toDataFrame()
    dfG0 = dfG.copy()
    dfC0 = dfC.copy()
    dfM  = dfC.copy()
    dfM[Cols] = np.nan
    time = dfG['Time_[s]'].values
    def removeJumps(vv, col='y'):
        dd = np.diff(vv)
        I1 = np.where(dd>50)[0]
        I2 = np.where(dd>50)[0]+1
        vv[I2] = vv[I1]
        if col.find('y')>=0: 
            vv = np.clip(vv, -250, 250)
        else:
            vv = np.clip(vv, -160, 50)
        return vv
    for c in Cols:
        dfG[c] = moving_average_conv(removeJumps(dfG[c].values, c), n=20)
        dfC[c] = moving_average_conv(removeJumps(dfC[c].values, c), n=20)
        dfM[c] = (dfG[c]*3+1*dfC[c])/4.
        ma1 =  moving_average_conv(dfM[c], n=5)
#         ma2 =  moving_average(dfM[c], n=5)
        dfM[c] = ma1
#         dfM[c][0:20] = ma2[0:20]
#         dfM[c][-20:] = ma2[-20:]
    return dfG0, dfC0, dfM



# --- Loop on all cases
def cleanTrajectory(Case, plot=True):
    case = Case['stability']
    LESpath    = os.path.join('02-iea15-{}WT/'.format(Case['nWT']), case);
    outputPath = os.path.join(LESpath, 'processedData')
    U0, dt, D, xyWT, HubHeight, xPlanes = getSimParamsAMR(case)

    dfG01, dfC01, dfM1 = processTrajectories(outputPath, 1, xPlanes)
    dfG02, dfC02, dfM2 = processTrajectories(outputPath, 2, xPlanes)

    # --- Save data
    dfM1.to_csv(os.path.join(outputPath, 'MeanderWT{:d}_TrajectoriesM.csv'.format(1)), index=False)
    dfM2.to_csv(os.path.join(outputPath, 'MeanderWT{:d}_TrajectoriesM.csv'.format(2)), index=False)

    # --- Plot
    if plot:
        nC=10
        Colr1=[color_scales(nC, 'blue')[1], color_scales(nC, 'green')[1], 'k']
        Colrs=[color_scales(nC, 'blue')[3], color_scales(nC, 'green')[3], 'k']
        time = dfG01['Time_[s]'].values
        for iP in xPlanes:
            yCols = ['y{}'.format(iP)]
            zCols = ['z{}'.format(iP)]
            Cols = yCols+zCols

            fig,axes = plt.subplots(4, 1, sharex=True,  sharey=False, figsize=(12.8,10.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.10, right=0.99, top=0.96, bottom=0.06, hspace=0.10, wspace=0.10)
            def plotyz(iWT, dfG0, dfC0, dfM, legendSet=False):
                ax = axes[(iWT-1)+0]
                koff=0
                for c in yCols:
                    vG0= dfG0[c] - koff*dfM[c].mean()
                    vC0= dfC0[c] - koff*dfM[c].mean()
                    vM = dfM[c]  - koff*dfM[c].mean()
                    ax.plot(time, vG0, c=Colrs[0], label='Gaussian' if not legendSet else None, lw=0.8, ls='-')
                    ax.plot(time, vC0, c=Colrs[1], label='Contour'  if not legendSet else None, lw=0.8, ls='-')
                    ax.plot(time, vM, c=Colrs[2], label='Clean'    if not legendSet else None, lw=1.0, ls='-')
                ax.set_ylabel(r'$y_{Wake} - y_0$'+'\n Turbine {} [m]'.format(iWT))
                ax.set_ylim([-250,250])
                if not legendSet:
                    ax.legend()
                legendSet=True

                ax = axes[(iWT-1)+2]
                for c in zCols:
                    vG0= dfG0[c] - koff*dfM[c].mean()
                    vC0= dfC0[c] - koff*dfM[c].mean()
                    vM = dfM[c] - koff*dfM[c].mean()
                    ax.plot(time, vG0, c=Colrs[0], label='Gaussian' if not legendSet else None, lw=0.8)
                    ax.plot(time, vC0, c=Colrs[1], label='Contour'  if not legendSet else None, lw=0.8)
                    ax.plot(time, vM, c=Colrs[2], label='Clean'     if not legendSet else None, lw=1.0, ls='-')
                ax.set_ylabel(r'$z_{Wake} - z_0$'+'\n Turbine {} [m]'.format(iWT))
                ax.set_ylim([-200,100])
            plotyz(1, dfG01, dfC01, dfM1, False)
            plotyz(2, dfG02, dfC02, dfM2, True)
            for ax in np.array(axes).flatten():
                ax.tick_params(direction='in', top=True, right=True, labelright=False, labeltop=False, which='both')

            axes[3].set_xlabel('Time [s]')
            figname = '{}__{}D'.format(casename,iP)
            fig.suptitle(figname)
            fig.savefig('_figs_WakeTrajectories/' + figname+'.png')

    return dfM1, dfM2


if __name__ == '__main__':

    Cases={}
    Cases['neutral2WT'] ={'stability':'neutral','nWT':2, 'planeFileWT':['planesT1129921.nc','planesT2129921.nc']}
    Cases['stable2WT']  ={'stability':'stable' ,'nWT':2, 'planeFileWT':['planesT1121692.nc','planesT2121692.nc']}
    Cases['unstable2WT']={'stability':'unstable' ,'nWT':2, 'planeFileWT':['planesT176826.nc','planesT276826.nc']}
    # Cases['neutral0WT'] ={'stability':'neutral','nWT':0, 'planeFileWT':['planesT1129921.nc','planesT2129921.nc']}
    # Cases['stable0WT']  ={'stability':'stable' ,'nWT':0, 'planeFileWT':['planesT1121692.nc','planesT2121692.nc']}

    for casename, Case in Cases.items():
        cleanTrajectory(Case, plot=True)

plt.show()
