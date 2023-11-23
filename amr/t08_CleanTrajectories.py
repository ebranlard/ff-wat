import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import welib.weio as weio
from welib.essentials import *
from helper_functions import *
from windtools.amrwind.post_processing  import Sampling
from welib.tools.signal_analysis import moving_average
from welib.tools.signal_analysis import moving_average_conv


def processTrajectories(caseName, outDirTraj, iWT, xPlanes):
    """ Read Gaussian and Contour trajectories, perform moving averaging on them """
    yCols = ['y{}'.format(iP) for iP in xPlanes]
    zCols = ['z{}'.format(iP) for iP in xPlanes]
    Cols = yCols+zCols
    # --- Reading trajectories
    dfG = weio.read(os.path.join(outDirTraj, '{}_MeanderWT{:d}_TrajectoriesC.csv'.format(caseName, iWT))).toDataFrame()
    dfC = weio.read(os.path.join(outDirTraj, '{}_MeanderWT{:d}_TrajectoriesG.csv'.format(caseName, iWT))).toDataFrame()
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
def cleanTrajectory(caseName, Case, plot=True, outDir='_out', figDir='_figs'):
    print(f'-----------------------------------------------------------------------')
    print(f'--- cleanTrajectory - Case: {caseName}')
    print(f'-----------------------------------------------------------------------')

    outDirTraj = os.path.join(outDir, 'trajectories')
    figDirTraj = os.path.join(figDir, '_figs_trajectories')
    if not os.path.exists(figDirTraj):
        os.makedirs(figDirTraj)

    # --- Read trajectories
    U0, dt, D, xyWT, HubHeight, xPlanes = getSimParamsAMR(Case['stability'])
    dfG01, dfC01, dfM1 = processTrajectories(caseName, outDirTraj, 1, xPlanes)
    dfG02, dfC02, dfM2 = processTrajectories(caseName, outDirTraj, 2, xPlanes)

    # --- Save data
    dfM1.to_csv(os.path.join(outDirTraj, '{}_MeanderWT{:d}_TrajectoriesM.csv'.format(caseName, 1)), index=False)
    dfM2.to_csv(os.path.join(outDirTraj, '{}_MeanderWT{:d}_TrajectoriesM.csv'.format(caseName, 2)), index=False)

    if not plot:
        return dfM1, dfM2
    # --- Plot
    nC=10
    Colr1=[color_scales(nC, 'blue')[1], color_scales(nC, 'green')[1], 'k']
    Colrs=[color_scales(nC, 'blue')[3], color_scales(nC, 'green')[3], 'k']
    time = dfG01['Time_[s]'].values

    zRange = [11.250,  511.250]
    zMid = (zRange[1]+zRange[0])/2
    zHub = 150
    zOff = zMid - zHub
    #print('>>> zMid', zMid)
    #print('>>> zOff', zOff)

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
                vG0= (dfG0[c] - koff*dfM[c].mean()) + zOff
                vC0= (dfC0[c] - koff*dfM[c].mean()) + zOff
                vM = (dfM[c] - koff*dfM[c].mean() ) + zOff
                ax.plot(time, vG0, c=Colrs[0], label='Gaussian' if not legendSet else None, lw=0.8)
                ax.plot(time, vC0, c=Colrs[1], label='Contour'  if not legendSet else None, lw=0.8)
                ax.plot(time, vM, c=Colrs[2], label='Clean'     if not legendSet else None, lw=1.0, ls='-')
            ax.set_ylabel(r'$z_{Wake} - z_0$'+'\n Turbine {} [m]'.format(iWT))
            ax.set_ylim([-200+zOff,100+zOff])

        plotyz(1, dfG01, dfC01, dfM1, False)
        plotyz(2, dfG02, dfC02, dfM2, True)
        for ax in np.array(axes).flatten():
            ax.tick_params(direction='in', top=True, right=True, labelright=False, labeltop=False, which='both')

        axes[3].set_xlabel('Time [s]')
        figname = '{}__{}D'.format(caseName,iP)
        fig.suptitle(figname)
        fig.savefig(os.path.join(figDirTraj, figname+'.png'))
    return dfM1, dfM2



if __name__ == '__main__':

#     caseNames = AllCases.keys()
    caseNames =[]
    caseNames += ['neutral2WT']
#     caseNames += ['stable2WT'] 
#     caseNames += ['unstable2WT'] 
#     caseNames += ['neutral1WT']
#     caseNames += ['stable1WT'] 
#     caseNames += ['unstable1WT'] 
    Cases = dict((k,v) for k,v in AllCases.items() if k in caseNames)

    for caseName, Case in Cases.items():
        with FailSafe(caseName, True):
            cleanTrajectory(caseName, Case, plot=True)

plt.show()
