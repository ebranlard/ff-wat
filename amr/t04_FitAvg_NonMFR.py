import matplotlib.pyplot as plt
from t11_FitAvg_MFR import fitAvgProfiles


if __name__ == '__main__':
    # --- Parameters 
    Meander=False
    cases=[]
    cases+=['neutral']
    cases+=['stable']
    cases+=['unstable']

    for case in cases:
        for symmetric in [True, False]:
            if case!='unstable':
                BG=['Outer','0WT']
            else:
                BG=['Outer']
            for bg in BG:
                fitAvgProfiles(case, Meander=Meander, symmetric=symmetric, removeBG=bg)


    plt.show()
