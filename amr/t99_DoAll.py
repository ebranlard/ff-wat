from helper_functions import *

from t00_PostPro_SaveHoriPlanes import saveHubHeightPlanes
from t01_Plot_HH import plotHH


if __name__ == '__main__':

    caseNames = AllCases.keys()
#     caseNames += ['neutral2WT']
#     caseNames += ['stable2WT'] 
#     caseNames += ['neutral0WT']
#     caseNames += ['stable0WT'] 
#     caseNames += ['unstable2WT'] 
# #     caseNames += ['unstable0WT'] 
# 
#     Cases = dict((k,v) for k,v in AllCases.items() if k in caseNames)
#     print(Cases)

    # --- Step 1 - Save planes at hubheight 
#     for caseName, Case in AllCases.items():
#         try:
#             saveHubHeightPlanes(caseName, Case)
#             OK(caseName)
#         except:
#             FAIL(caseName)


#     # --- Step  - Plot Mean deficit at hubheight and time series at hub height
#     for caseName, Case in AllCases.items():
#         try:
#             plotHH(caseName, Case, Meander=False)
#             OK(caseName)
#         except:
#             FAIL(caseName)
