import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray
# Local 
import weio
from weio.pickle_file import PickleFile
from welib.essentials import *
from welib.tools.curve_fitting import model_fit
from helper_functions import *

outPath='_out_KFit/'
figsPath = '_figs_KFit/'


# ---- Read all
def loadallKFits(Meander, syms, cases, rms, smooth=0):
    D={}
    for case in cases:
        for symmetric in syms:
            for removeBG in rms:
                if Meander: 
                    label='KFit_Meander_'+case
                else:
                    label='KFit_Inertial_'+case
                if symmetric:
                    label+='_sym'
                label+='_rm'+removeBG+'_smooth{:d}'.format(smooth)

                filename = outPath+label+'.pkl'
                try:
                    pkl = PickleFile(filename)
                    D[label] = pkl['data']
                except:
                    print('>>> Missing', filename)
    return D

def plotKFits(Meander, D, smooth=0):
    for k,v in D.items():
        print('Key: ', k)
    xPlanes=np.arange(7)
    colrsNeut = color_scales(n=7, color='blue' )
    colrsStab = color_scales(n=7, color='green' )
    colrsUnst = color_scales(n=7, color='red' )

    fig,axes = plt.subplots(2, 2, sharex=True, figsize=(12.8,9.6)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    iNeut=0
    iStab=0
    iUnst=0
    for k,KS in D.items():
        label=k
        k = k.replace('KFit_','')
        k = k.replace('Meander','MF')
        k = k.replace('Inertial','IF')
        if 'neutral' in label:
            iNeut+=1
            col = colrsNeut[iNeut]
        elif 'unstable' in label:
            iUnst+=1
            col = colrsUnst[iUnst]
        elif 'stable' in label:
            iStab+=1
            col = colrsStab[iStab]
        for iRow in [0,1]:
            axes[iRow,0].plot(xPlanes, KS[iRow,:,0], c=col, label=k)
            axes[iRow,0].set_xlabel('')
            axes[iRow,0].set_ylabel('Kdef')

            axes[iRow,1].plot(xPlanes, KS[iRow,:,1], c=col, label=k)
            axes[iRow,1].set_xlabel('')
            axes[iRow,1].set_ylabel('Kgrad')

    if Meander:
        prefix='KFit_Meander_'
        sufix='Meander_'
    else:
        prefix='KFit_Inertial_'
        sufix='Inertial_'

#     KS = D[prefix+'neutral_sym_rm0WT']
#     axes[1,0].plot(xPlanes, np.sqrt(KS[0,-1,0]**2 + KS[0,:,0]**2) , c=colrsNeut[0], ls='--')
#     axes[1,1].plot(xPlanes, np.sqrt(KS[0,-1,1]**2 + KS[0,:,1]**2) , c=colrsNeut[0], ls='--')

#     KS = D[prefix+'stable_sym_rm0WT']
#     axes[1,0].plot(xPlanes, np.sqrt(KS[0,-1,0]**2 + KS[0,:,0]**2) , c=colrsStab[0], ls='--')
#     axes[1,1].plot(xPlanes, np.sqrt(KS[0,-1,1]**2 + KS[0,:,1]**2) , c=colrsStab[0], ls='--')

#     KS = D[prefix+'unstable_sym_rm0WT']
#     axes[1,0].plot(xPlanes, np.sqrt(KS[0,-1,0]**2 + KS[0,:,0]**2) , c=colrsUnst[0], ls='--')
#     axes[1,1].plot(xPlanes, np.sqrt(KS[0,-1,1]**2 + KS[0,:,1]**2) , c=colrsUnst[0], ls='--')


    for ax in np.array(axes).flatten():
        ax.set_ylim([0,6.5])
        ax.set_xlabel('x/D [-]')
        ax.tick_params(direction='in', top=True, right=True, labelright=False, labeltop=False, which='both')

    axes[0,0].legend(fontsize=9, ncol=2, loc='upper left')


    figname = 'Summary_Downstream_'+sufix+'_smooth'+str(smooth)
    fig.suptitle(figname)
    fig.savefig(figsPath+figname+'.png')




if __name__ == '__main__':
    syms = [True, False]
    cases=['neutral', 'stable','unstable']
    rms =['0WT', 'Outer']
    cases=['neutral']
    rms =['0WT']

    smooth=0

    D = loadallKFits(False, syms, cases, rms, smooth=smooth)
    plotKFits(False, D, smooth=smooth)

#     D = loadallKFits(True, syms, cases, rms, smooth=smooth)
#     plotKFits(True, D, smooth=smooth)

    plt.show()

    pass
