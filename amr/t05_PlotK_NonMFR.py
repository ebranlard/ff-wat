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

outPath='_out/'
figsPath = '_figs/'

syms = [True, False]
cases=['neutral', 'stable']
rms =['0WT', 'Outer']

# ---- Read all
D={}
for case in cases:
    for symmetric in syms:
        for removeBG in rms:
            label=case+'_rm'+removeBG
            if symmetric:
                label+='_sym'

            filename = outPath+'Kfit_'+label+'.pkl'
            pkl = PickleFile(filename)
            D[label] = pkl['data']

xPlanes=np.arange(7)
colrsNeut = color_scales(n=7, color='blue' )
colrsStab = color_scales(n=7, color='green' )

fig,axes = plt.subplots(2, 2, sharex=True, figsize=(12.8,9.6)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
iNeut=0
iStab=0
for k,KS in D.items():
    label=k
    if 'neutral' in label:
        iNeut+=1
        col = colrsNeut[iNeut]
    if 'stable' in label:
        iStab+=1
        col = colrsStab[iStab]
    for iRow in [0,1]:
        axes[iRow,0].plot(xPlanes, KS[iRow,:,0], c=col, label=k)
        axes[iRow,0].set_xlabel('')
        axes[iRow,0].set_ylabel('Kdef')

        axes[iRow,1].plot(xPlanes, KS[iRow,:,1], c=col, label=k)
        axes[iRow,1].set_xlabel('')
        axes[iRow,1].set_ylabel('Kgrad')

KS = D['neutral_rm0WT_sym']
axes[1,0].plot(xPlanes, np.sqrt(KS[0,-1,0]**2 + KS[0,:,0]**2) , c=colrsNeut[0], ls='--')
axes[1,1].plot(xPlanes, np.sqrt(KS[0,-1,1]**2 + KS[0,:,1]**2) , c=colrsNeut[0], ls='--')

KS = D['stable_rm0WT_sym']
axes[1,0].plot(xPlanes, np.sqrt(KS[0,-1,0]**2 + KS[0,:,0]**2) , c=colrsStab[0], ls='--')
axes[1,1].plot(xPlanes, np.sqrt(KS[0,-1,1]**2 + KS[0,:,1]**2) , c=colrsStab[0], ls='--')

for ax in np.array(axes).flatten():
    ax.set_ylim([0,5.5])
    ax.set_xlabel('x/D [-]')
    ax.tick_params(direction='in', top=True, right=True, labelright=False, labeltop=False, which='both')

axes[0,0].legend(fontsize=9, ncol=2, loc='upper left')


figname = figsPath+'KFit_KwrtDownstream'
fig.suptitle(figname)
fig.savefig(figname+'.png')

plt.show()

if __name__ == '__main__':
    pass
