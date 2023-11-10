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


Cases={}
Cases['neutral2'] ={'stability':'neutral','nWT':2, 'planeFileWT':['planesT1129921.nc','planesT2129921.nc']}
# Cases['stable2']  ={'stability':'stable' ,'nWT':2, 'planeFileWT':['planesT1121692.nc','planesT2121692.nc']}
# Cases['unstable2']={'stability':'unstable' ,'nWT':2, 'planeFileWT':['planesT176826.nc','planesT276826.nc']}
# Cases['neutral0'] ={'stability':'neutral','nWT':0, 'planeFileWT':['planesT1129921.nc','planesT2129921.nc']}
# Cases['stable0']  ={'stability':'stable' ,'nWT':0, 'planeFileWT':['planesT1121692.nc','planesT2121692.nc']}

# --- Loop on all cases
for casename, Case in Cases.items():
    case = Case['stability']
    LESpath    = os.path.join('02-iea15-{}WT/'.format(Case['nWT']), case);
    outputPath = os.path.join(LESpath, 'processedData')
    U0, dt, D, xyWT1, xyWT2, xPlanes = getSimParamsAMR(case)
    for iWT in [1, 2]:
        # --- Reading trajectories
        dfG = weio.read(os.path.join(outputPath, 'MeanderWT{:d}_TrajectoriesC.csv'.format(iWT))).toDataFrame()
        dfC = weio.read(os.path.join(outputPath, 'MeanderWT{:d}_TrajectoriesC.csv'.format(iWT))).toDataFrame()
        dfG0 = dfG.copy()
        dfC0 = dfC.copy()

        yCols = ['y{}'.format(iP) for iP in xPlanes]
        zCols = ['z{}'.format(iP) for iP in xPlanes]
        Cols = yCols+zCols

#         DG=dfG[Cols].diff()
#         DC=dfC[Cols].diff()

        dfG[yCols] = dfG[yCols].clip(-210, 210)
        dfC[yCols] = dfC[yCols].clip(-210, 210)
        dfG[zCols] = dfG[zCols].clip(-160, 50 )
        dfC[zCols] = dfC[zCols].clip(-160, 50 )


#                 if it>ITime[0]:
#                     # If difference between methods is too large use closest to previous value
#                     yc1, yc2 = choose(yc1, yc2, yc_prev)
#                     zc1, zc2 = choose(zc1, zc2, zc_prev)
#                     yc = (1*yc1+3*yc2)/4
#                     zc = (1*zc1+3*zc2)/4
#                     # simple moving average
#                     zc = (zc + zc_prev)/2
#                     yc = (yc + yc_prev)/2
#                 else:
#                     yc = (1*yc1+3*yc2)/4
#                     zc = (1*zc1+3*zc2)/4
# 
#                 yc1_prev = yc1
#                 yc2_prev = yc2
#                 zc1_prev = zc1
#                 zc2_prev = zc2
#                 yc_prev  = yc
#                 zc_prev  = zc



        print(dfTrajG)



#         # --- Save data
#         cols = ['Time_[s]'] 
#         for iP in xPlanes:
#             cols += ['y{}'.format(iP) , 'z{}'.format(iP)]
#         df = pd.DataFrame(data=np.column_stack((ITime, WakeTrajectoriesM)), columns=cols)
#         df.to_csv(os.path.join(outputPath, 'MeanderWT{:d}_TrajectoriesM.csv'.format(iWT)), index=False)


plt.show()
