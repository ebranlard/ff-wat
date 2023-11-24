import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import welib.weio as weio
from welib.essentials import *
from plothelper import *
setFigureFont(15)

export=False
export=True

outDirTraj='amr/_out_all/trajectories/'
df={}
for stab in stabilities:
    caseName = stab+'2WT' # TODO 
    df[stab] = weio.read(os.path.join(outDirTraj, '{}_MeanderWT{:d}_TrajectoriesM.csv'.format(caseName, 1))).toDataFrame()

zOff=111
D=240

stys=[':', '--', '-']

fig,axes = plt.subplots(2, 1, sharey=False, sharex=True, figsize=(6.4,3.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.145, right=0.97, top=0.982, bottom=0.14, hspace=0.08, wspace=0.20)
for i, stab in enumerate(stabilities):
    axes[0].plot(df[stab]['Time_[s]']-600, df[stab]['y6'], ls=stys[i], c=colrsStab[i], label=stab)
    axes[1].plot(df[stab]['Time_[s]']-600, df[stab]['z6']+zOff, ls=stys[i], c=colrsStab[i], label=stab)
axes[0].set_ylabel('y/D [-]')
axes[1].set_ylabel('z/D [-]')
axes[1].set_xlabel('Time [s]')
for ax in axes:
    ax.tick_params(direction='in', top=True, right=True, labelright=False, labeltop=False, which='both')
    ax.set_xlim([0,600])
axes[0].legend(fontsize=13)
fig.suptitle('WakeTracking')









if export:
    export2pdf()




plt.show()

if __name__ == '__main__':
    pass
