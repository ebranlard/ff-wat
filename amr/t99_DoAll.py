from helper_functions import *
import glob

# Postpro
from t00_PostPro_SaveHoriPlanes import saveHubHeightLines
from t07_FindWakeCenters import findWakeCenters
from t08_CleanTrajectories import cleanTrajectory
from t09_PostProMeanderingFrame import saveMeanderingFrameLines
from t11_FitAvg_MFR import fitAvgProfiles
# Plot
from t01_Plot_HH import plotHH
from t02_Compare_Free_Waked import plotCompareFreeWake
from t10_CompareMean_MFR_Inertial import plotMFRInertialComp
from t12_PlotK import loadallKFits, plotKFits  


if __name__ == '__main__':
    # --- Parameters
    outDir = '_out_all_new'
    figDir = '_out_all_new'
    iTimeMin = 500
    raiseException=False

    # --- Select case names
    STEPS = [0, 1, 2, 3, 4, 5, 6]
    STEPS = [0]

    caseNames =[]
    caseNames += list(AllCases.keys())
#     caseNames += ['neutral2WT']
#     caseNames += ['stable2WT'] 
#     caseNames += ['unstable2WT'] 
#     caseNames += ['neutral0WT']
#     caseNames += ['stable0WT'] 
#     caseNames += ['unstable0WT'] 
#     caseNames += ['neutral1WT']
#     caseNames += ['stable1WT'] 
#     caseNames += ['unstable1WT'] 

    # --- Cases and stability from caseNames
    Cases = dict((k,v) for k,v in AllCases.items() if k in caseNames)
    stabilities = np.unique([Case['stability'] for _,Case in Cases.items()])
    nWTs = np.unique([Case['nWT'] for _,Case in Cases.items() if Case['nWT']>0])

    # --- Step 0 - Check files
    if 0 in STEPS:
        for caseName, Case in Cases.items():
            pattern = os.path.join(Case['path'], 'post_processing', 'planesT{}'.format(1)) +'*.nc'
            files = glob.glob(pattern)
            with FailSafe(caseName + ' Number of files: {}'.format(len(files)), raiseException=True):
                if len(files)==0:
                    raise Exception('No plane files found with pattern: {}'.format(pattern))

    # --- Step 1 - Save planes at hubheight in inertial frame
    if 1 in STEPS:
        for caseName, Case in Cases.items():
            with FailSafe(caseName, raiseException=raiseException):
                saveHubHeightLines(caseName, Case, outDir=outDir)

    # --- Step 2 - Find wake centers and clean trajectories (NOTE: only for simulations with 1 or 2 turbines)
    if 2 in STEPS:
        for caseName, Case in Cases.items():
            if Case['nWT']>0:
                with FailSafe(caseName, raiseException=raiseException):
                    findWakeCenters(caseName, Case, plot=False, outDir=outDir, iTimeMin=iTimeMin)

    # --- Step 3 - Clean trajectories
    if 3 in STEPS:
        for caseName, Case in Cases.items():
            if Case['nWT']>0:
                with FailSafe(caseName, raiseException=raiseException):
                    cleanTrajectory(caseName, Case, plot=True, outDir=outDir, figDir=figDir)

    # --- Step 4 - Save in meandering frame
    if 4 in STEPS:
        for caseName, Case in Cases.items():
            with FailSafe(caseName, raiseException=raiseException):
                saveMeanderingFrameLines(caseName, Case, plot=True, nFigsMax=10, outDir=outDir, figDir=figDir, iTimeMin=iTimeMin)

    # --- Step 5 - Fit K profile 
    if 5 in STEPS:
        BG=['Outer','0WT']
        for caseName, Case in Cases.items():
            if Case['nWT']>0:
                for smooth in [0,2,4]:
                    for Meander in [True, False]:
                        for symmetric in [True, False]:
                            for bg in BG:
                                label = getLabel(caseName, Meander, symmetric, bg, smooth)
                                with FailSafe(label, False):
                                    fitAvgProfiles(caseName, Case, Meander=Meander, symmetric=symmetric, removeBG=bg, smooth=smooth, outDir=outDir, figDir=figDir)



    # --------------------------------------------------------------------------------}
    # --- Plots 
    # --------------------------------------------------------------------------------{
    if 6 in STEPS:
        # --- Plot Mean deficit at hubheight and time series at hub height
        for caseName, Case in AllCases.items():
            with FailSafe(caseName, raiseException=raiseException):
                plotHH(caseName, Case, Meander=False, outDir=outDir, figDir=figDir)
        # --- Compare 0WT and 2WT - Inertial
        for stability in stabilities:
            with FailSafe(stability, raiseException=raiseException):
                plotCompareFreeWake(stability, Meander=False, outDir=outDir, figDir=figDir)

        # --- Compare Inertial MFR
        for caseName, Case in AllCases.items():
            with FailSafe(caseName, raiseException=raiseException):
                plotMFRInertialComp(caseName, Case, outDir=outDir, figDir=figDir, iTimeMin=iTimeMin)


        # --- Plot K Fits with downstream distance
        syms = [True, False]
        rms =['0WT', 'Outer']
        Meander=False
        Meander=True

        for nWT in nWTs:
            for Meander in [True, False]:
                for smooth in [0,2,4]:
                    with FailSafe('PlotKFit Meander{} nWT{} Smooth{}'.format(Meander, nWT, smooth), raiseException=raiseException):
                        D = loadallKFits(Meander, syms, Cases, rms, smooth=smooth, nWT=nWT, outDir=outDir)
                        plotKFits(Meander, D, smooth=smooth, nWT=nWT, figDir=figDir)


