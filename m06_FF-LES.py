import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatch
# Local 
import welib.weio as weio
from welib.tools.fatigue import equivalent_load
from welib.weio.pickle_file import PickleFile
from welib.essentials import *
import sys
sys.path.append('amr')
from plothelper import *
from helper_functions import *
setFigureFont(13)



# --- Parameters
stabilities=['stable','neutral']
outFileBase = '_data/LES-FF'
colKeep= ['Time_[s]','RootMyc1_[kN-m]','TwrBsMyt_[kN-m]','TwrBsMxt_[kN-m]']
postPro=False
analyse=False
plot=True

tMin = 500
tMax = 1677 # TODO


# --- PostPro
if postPro:
    def readsafe(filename):
        print('Reading ', filename)
        if os.path.exists(filename):
            df = weio.read(filename).toDataFrame()
            try:
                df=df[colKeep]
                OK(filename)
                return df
            except:
                FAIL(filename+ ' >>> Bad columns:'+ df.columns)
                import pdb; pdb.set_trace()
        else:
            FAIL('Missing: '+ filename)
        return None

    dfs = {}
    for stab in stabilities:

        # --- FAST.Farm
        # amr/02-iea15-2WT/neutral/IEA-15-240-RWT-Monopile.T1.conv.outb
        dfs[stab+'_LES_T1'] = readsafe(os.path.join('amr','02-iea15-2WT',stab, 'IEA-15-240-RWT-Monopile.T1.conv.outb'))
        dfs[stab+'_LES_T2'] = readsafe(os.path.join('amr','02-iea15-2WT',stab, 'IEA-15-240-RWT-Monopile.T2.conv.outb'))
        # ff/stable/FF-NoWAT.T1.outb
        # ff/stable/FF-WAT-More.T1.outb
        # ff/stable/FF-WAT-More.T1.outb
        dfs[stab+'_NoWAT_T1'] = readsafe(os.path.join('ff',stab, 'FF-NoWAT.T1.outb'))
        dfs[stab+'_NoWAT_T2'] = readsafe(os.path.join('ff',stab, 'FF-NoWAT.T2.outb'))
        dfs[stab+'_WAT_T1'] = readsafe(os.path.join('ff',stab, 'FF-WAT-More.T1.outb'))
        dfs[stab+'_WAT_T2'] = readsafe(os.path.join('ff',stab, 'FF-WAT-More.T2.outb'))

    pkl = PickleFile(data=dfs)
    pkl.write(outFileBase+'.pkl')

if analyse:
    # --- Analyses
    pkl = PickleFile(outFileBase+'.pkl')
    print(pkl)

    # Restrict time
    for k,df in pkl.items():
        print('Time Range: ', df['Time_[s]'].min(),  df['Time_[s]'].max(), k)
        df = df[df['Time_[s]']>600]
        df = df[df['Time_[s]']<1800]
        pkl[k] = df

    # --- LEq / Sig
    Sig={}
    Leq={}
    for k,df in pkl.items():
        LeqBMy = equivalent_load(df['Time_[s]'].values, df['RootMyc1_[kN-m]'].values, m=10)
        LeqTMx = equivalent_load(df['Time_[s]'].values, df['TwrBsMxt_[kN-m]'].values, m=4)
        LeqTMy = equivalent_load(df['Time_[s]'].values, df['TwrBsMyt_[kN-m]'].values, m=4)
        sigBMy = df['RootMyc1_[kN-m]'].std()
        sigTMx = df['TwrBsMxt_[kN-m]'].std()
        sigTMy = df['TwrBsMyt_[kN-m]'].std()
        Leq[k] = {'BMy':LeqBMy, 'TMy':LeqTMy, 'TMx':LeqTMx}
        Sig[k] = {'BMy':sigBMy, 'TMy':sigTMy, 'TMx':sigTMx}
        print(Leq[k])
    # --- Ratios T1/T2
    SigRat={}
    LeqRat={}
    for k1 in Sig.keys():
        k = k1.replace('_T1','')
        if k1.find('T1')>1:
            k2 = k1.replace('T1','T2')
            s1 = Sig[k1]
            s2 = Sig[k2]
            l1 = Leq[k1]
            l2 = Leq[k2]
            SigRat[k]=dict( [(l, s2[l]/s1[l]) for l in s1.keys()  ] )
            LeqRat[k]=dict( [(l, l2[l]/l1[l]) for l in l1.keys()  ] )


    pkl = PickleFile()
    pkl['Leq']=Leq
    pkl['Sig']=Sig
    pkl['LeqRat']=LeqRat
    pkl['SigRat']=SigRat
    pkl.write(outFileBase+'SigLeq.pkl')



# --- Plot

stats = PickleFile(outFileBase+'SigLeq.pkl')
print(stats)
Cases=['NoWAT', 'WAT', 'LES']

var='BMy'
sstat='SigRat'
for stab in stabilities:

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    for c in Cases:
        k = stab+'_'+ c
        ax.plot(k, stats[sstat][var])
# 
#             off=1.2*width/2
#             off2=0
#             if iL==3:
#                 off2=-1.5*width
#             ax.bar(iL+off2-off, U25_AD, width, label=labels[0], ec=cVC  , fc=cAD)
#             ax.bar(iL+off2+off, U25_VC, width, label=labels[1], ec=cAD , fc=cVC)
#             eps=np.abs((U25_AD-U25_VC)/U25_AD)
#             ax.text(iL+off2,max(U25_AD,U25_VC)*1.0003,'{:.2f}%'.format(eps*100),fontsize=10,ha='center')


    ax.set_xlabel('')
    ax.set_ylabel(sstat+ ' '+ var)
    ax.legend()
plt.show()



















# import pdb; pdb.set_trace()


plt.show()

