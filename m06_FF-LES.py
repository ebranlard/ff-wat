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
figDir      = '_data/_figs_LES-FF'
# colKeep= ['Time_[s]','RootMyb1_[kN-m]','RootMxb1_[kN-m]','RootMxc1_[kN-m]','RootMyc1_[kN-m]','TwrBsMyt_[kN-m]','TwrBsMxt_[kN-m]']
colKeep= ['Time_[s]','RootMyb1_[kN-m]','RootMxb1_[kN-m]','RootMxc1_[kN-m]','RootMyc1_[kN-m]','TwrBsMyt_[kN-m]','TwrBsMxt_[kN-m]','YawBrMzp_[kN-m]']
# postPro=True
postPro=False
# plot=False
plot=True

tMin = 500
tMax = 1677 # TODO

# --- Derived params
if not os.path.exists(figDir):
    os.makedirs(figDir)
if not os.path.exists(os.path.dirname(outFileBase)):
    os.makedirs(os.path.dirname(outFileBase))

# --------------------------------------------------------------------------------}
# --- PostPro
# --------------------------------------------------------------------------------{
if postPro:
    def readsafe(filename):
        print('Reading ', filename)
        if os.path.exists(filename):
            df = weio.read(filename).toDataFrame()
            try:
                df=df[colKeep]
                dt = np.around(df['Time_[s]'].values[1]- df['Time_[s]'].values[0],4)
                if dt==0.005:
                    df = df.iloc[::4,:]
                elif dt==0.02:
                    df = df
                else:
                    FAIL('>>>>> NOT IMPLEMENTED'+ dt)
                    import pdb; pdb.set_trace()
                # 
                dt = np.around(df['Time_[s]'].values[1]- df['Time_[s]'].values[0],4)
                if dt!=0.02:
                    FAIL('Wrong dt')
                    import pdb; pdb.set_trace()
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
        dfs[stab+'_WATOpt_T1'] = readsafe(os.path.join('ff',stab, 'FF-WAT.T1.outb'))
        dfs[stab+'_WATOpt_T2'] = readsafe(os.path.join('ff',stab, 'FF-WAT.T2.outb'))
        dfs[stab+'_WATOptHR_T1'] = readsafe(os.path.join('ff',stab, 'FF-WAT-More.T1.outb'))
        dfs[stab+'_WATOptHR_T2'] = readsafe(os.path.join('ff',stab, 'FF-WAT-More.T2.outb'))
        dfs[stab+'_LES_T1'] = readsafe(os.path.join('amr','02-iea15-2WT',stab, 'IEA-15-240-RWT-Monopile.T1.conv.outb'))
        dfs[stab+'_LES_T2'] = readsafe(os.path.join('amr','02-iea15-2WT',stab, 'IEA-15-240-RWT-Monopile.T2.conv.outb'))
        # ff/stable/FF-NoWAT.T1.outb
        # ff/stable/FF-WAT-More.T1.outb
        # ff/stable/FF-WAT-More.T1.outb
        dfs[stab+'_NoWAT_T1'] = readsafe(os.path.join('ff',stab, 'FF-NoWAT.T1.outb'))
        dfs[stab+'_NoWAT_T2'] = readsafe(os.path.join('ff',stab, 'FF-NoWAT.T2.outb'))
        dfs[stab+'_WAT05_T1'] = readsafe(os.path.join('ff',stab, '_KEEP_ME', 'FF-WAT05.T1.outb'))
        dfs[stab+'_WAT05_T2'] = readsafe(os.path.join('ff',stab, '_KEEP_ME', 'FF-WAT05.T2.outb'))
        dfs[stab+'_WAT15_T1'] = readsafe(os.path.join('ff',stab, '_KEEP_ME', 'FF-WAT15.T1.outb'))
        dfs[stab+'_WAT15_T2'] = readsafe(os.path.join('ff',stab, '_KEEP_ME', 'FF-WAT15.T2.outb'))

    pkl = PickleFile(data=dfs)
    pkl.write(outFileBase+'.pkl')

    # --- Analyses
    pkl = PickleFile(outFileBase+'.pkl')
    print(pkl)

    # Restrict time
    for k,df in pkl.items():
        print('Time Range: ', df['Time_[s]'].min(),  df['Time_[s]'].max(), k)
        df = df[df['Time_[s]']>tMin]
        df = df[df['Time_[s]']<tMax]
        pkl[k] = df

    # --- LEq / Sig
    Sig={}
    Leq={}
    Muu={}
    for k,df in pkl.items():
        method ='fatpack'
        LeqBMxb = np.around(equivalent_load(df['Time_[s]'].values, df['RootMxb1_[kN-m]'].values, m=10, method=method),1)
        LeqBMyb = np.around(equivalent_load(df['Time_[s]'].values, df['RootMyb1_[kN-m]'].values, m=10, method=method),1)
        LeqBMx = np.around(equivalent_load(df['Time_[s]'].values, df['RootMxc1_[kN-m]'].values, m=10, method=method),1)
        LeqBMy = np.around(equivalent_load(df['Time_[s]'].values, df['RootMyc1_[kN-m]'].values, m=10, method=method),1)
        LeqTMx = np.around(equivalent_load(df['Time_[s]'].values, df['TwrBsMxt_[kN-m]'].values, m=4 , method=method),1)
        LeqTMy = np.around(equivalent_load(df['Time_[s]'].values, df['TwrBsMyt_[kN-m]'].values, m=4 , method=method),1)
        try:
            Mz = df['YawBrMzp_[kN-m]'].iloc[:,0]
        except:
            Mz = df['YawBrMzp_[kN-m]']

        LeqTMz = np.around(equivalent_load(df['Time_[s]'].values, Mz.values, m=4 , method=method),1)
        sigBMxb= np.around(df['RootMxb1_[kN-m]'].std(),1)
        sigBMyb= np.around(df['RootMyb1_[kN-m]'].std(),1)
        sigBMx = np.around(df['RootMxc1_[kN-m]'].std(),1)
        sigBMy = np.around(df['RootMyc1_[kN-m]'].std(),1)
        sigBMz = np.around(df['RootMyc1_[kN-m]'].std(),1)
        sigTMx = np.around(df['TwrBsMxt_[kN-m]'].std(),1)
        sigTMy = np.around(df['TwrBsMyt_[kN-m]'].std(),1)
        sigTMz = np.around(Mz.std(),1)

        muuBMxb= np.around(df['RootMxb1_[kN-m]'].mean(),2)
        muuBMyb= np.around(df['RootMyb1_[kN-m]'].mean(),2)
        muuBMx = np.around(df['RootMxc1_[kN-m]'].mean(),2)
        muuBMy = np.around(df['RootMyc1_[kN-m]'].mean(),2)
        muuTMx = np.around(df['TwrBsMxt_[kN-m]'].mean(),2)
        muuTMy = np.around(df['TwrBsMyt_[kN-m]'].mean(),2)
        muuTMz = np.around(Mz.mean(),4)

        Leq[k] = {'BMxb':LeqBMx, 'BMyb':LeqBMy,'BMxc':LeqBMx, 'BMyc':LeqBMy, 'TMy':LeqTMy, 'TMx':LeqTMx , 'TMz':LeqTMz}
        Sig[k] = {'BMxb':sigBMx, 'BMyb':sigBMy,'BMxc':sigBMx, 'BMyc':sigBMy, 'TMy':sigTMy, 'TMx':sigTMx , 'TMz':sigTMz}
        Muu[k] = {'BMxb':muuBMx, 'BMyb':muuBMy,'BMxc':muuBMx, 'BMyc':muuBMy, 'TMy':muuTMy, 'TMx':muuTMx , 'TMz':muuTMz}
        print(Leq[k], k)
    # --- Ratios T1/T2
    print('>>>> Ratios')
    MuuRat={}
    SigRat={}
    LeqRat={}
    for k1 in Sig.keys():
        k = k1.replace('_T1','')
        if k1.find('_T1')>1:
            k2 = k1.replace('_T1','_T2')
            s1 = Sig[k1]
            s2 = Sig[k2]
            l1 = Leq[k1]
            l2 = Leq[k2]
            m1 = Muu[k1]
            m2 = Muu[k2]
            SigRat[k]=dict( [(l, np.around(s2[l]/s1[l],3)) for l in s1.keys()  ] )
            LeqRat[k]=dict( [(l, np.around(l2[l]/l1[l],3)) for l in l1.keys()  ] )
            MuuRat[k]=dict( [(l, np.around(m2[l]/m1[l],3)) for l in m1.keys()  ] )
            print(k1, l1)
            print(k2, l2)
            print(k+ '_Rat', LeqRat[k])

    pkl = PickleFile()
    pkl['Muu']=Muu
    pkl['Sig']=Sig
    pkl['Leq']=Leq
    pkl['MuuRat']=MuuRat
    pkl['SigRat']=SigRat
    pkl['LeqRat']=LeqRat
    pkl.write(outFileBase+'SigLeq.pkl')



# --------------------------------------------------------------------------------}
# --- Plot 
# --------------------------------------------------------------------------------{
colrsStabLight=[ lighten_color(c, factor=1.1) for c in colrsStab]
prettyNames={}
prettyNames['BMxb']='Blade edge moment'
prettyNames['BMyb']='Blade flap moment'
prettyNames['BMxc']='Blade edge moment (c)'
prettyNames['BMyc']='Blade flap moment (c)'
prettyNames['TMx']='Tower side-side moment'
prettyNames['TMy']='Tower fore-aft moment'
prettyNames['TMz']='Tower yaw'
prettyStats={}
prettyStats['SigRat']=r'$\sigma$'
prettyStats['LeqRat']=r'$L_{eq}$'
prettyStats['MuuRat']=r'$\mu$'

# Cases=['NoWAT', 'WAT05', 'WAT15', 'WATOpt', 'WATOptHR', 'LES']
Cases=['NoWAT', 'WATOpt', 'LES']

if plot:
    stats = PickleFile(outFileBase+'SigLeq.pkl')
    print(stats)
    #for var in ['BMxb','BMyb','BMxc','BMyc', 'TMy', 'TMx']:
#     for var in ['BMyb', 'TMy', 'TMz']:
    for var in ['TMy','TMz']:
        for sstat in ['SigRat']: #,'LeqRat','MuuRat']:
    #for var in ['TMy']:
    #    for sstat in ['LeqRat']:
            fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
            width=0.8
            xticklabels=[]
            xticks=[]
            j=0
            for iStab, stab in enumerate(stabilities):
                for i,c in enumerate(Cases):
                    k = stab+'_'+ c
                    #ax.plot(k, stats[sstat][k][var])
                    ax.bar(j, stats[sstat][k][var], width, label=stab.capitalize() if (i==0) else None, ec=colrsStab[iStab], fc=colrsStabLight[iStab])
                    xticks.append(j)
                    xticklabels.append(c.replace('WAT',''))
                    j+=1
                j+=1
            if sstat.find('MuuRat')==0:
                ax.set_ylim([0.0,1.1])
                pass
            elif var.find('BMx')==0:
                ax.set_ylim([0.9,1.1])
            else:
                ax.set_ylim([0,4.0])
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels)
            ax.tick_params(direction='in', top=True, right=True, labelright=False, labeltop=False, which='both')
            # ax.set_title(stab.capitalize())
            figname = var+' ' + sstat
            ax.set_title(prettyStats[sstat] + ' ' + prettyNames[var])
            ax.set_xlabel('')
            ax.set_ylabel(r'Ratio $T_2/T_1$ - '+prettyStats[sstat] + ' ' + prettyNames[var]+ ' [-]')
            ax.legend()
            filename=os.path.join(figDir, figname.replace(' ','_')+'.png')
            fig.savefig(filename)
    plt.show()

















# 
#             off=1.2*width/2
#             off2=0
#             if iL==3:
#                 off2=-1.5*width
#             ax.bar(iL+off2-off, U25_AD, width, label=labels[0], ec=cVC  , fc=cAD)
#             ax.bar(iL+off2+off, U25_VC, width, label=labels[1], ec=cAD , fc=cVC)
#             eps=np.abs((U25_AD-U25_VC)/U25_AD)
#             ax.text(iL+off2,max(U25_AD,U25_VC)*1.0003,'{:.2f}%'.format(eps*100),fontsize=10,ha='center')
#                 for bar in bar2:
#                     bar.set_hatch('//')


# import pdb; pdb.set_trace()


plt.show()

