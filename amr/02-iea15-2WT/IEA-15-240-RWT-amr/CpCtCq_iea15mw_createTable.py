# -------------------------------------------------------------------------------------------------------------- #
# Read table of coefficients from the IEA 15 MW turbine and re-write in different format for the AeroDisk model
# Regis Thedin
# Sep 2022
# -------------------------------------------------------------------------------------------------------------- #

import pandas as pd
import numpy as np

# Full path to Cp, Ct, Cq file. Note this changes between branches.
cpctcqfile='/home/rthedin/repos/IEA-15-240-RWT/OpenFAST/IEA-15-240-RWT/Cp_Ct_Cq.IEA15MW.txt'

# Get pitch and TSR from first few lines
pitch = pd.read_csv(cpctcqfile, sep='  ', skiprows=4, nrows=1, header=None, engine='python').to_numpy()[0]
tsr = pd.read_csv(cpctcqfile, sep='   ', skiprows=6, nrows=1, header=None, engine='python').to_numpy()[0]


if len(pitch) < 2:
    raise valueError('There might have been an error reading the pitch line. Maybe the spacing between the numbers has changed? Try again.')
if len(tsr) < 2:
    raise valueError('There might have been an error reading the TSR line. Maybe the spacing between the numbers has changed? Try again.')
    
# Open the main tables    
cp = pd.read_csv(cpctcqfile, sep='  ', skiprows=12, nrows=len(tsr), header=None, engine='python').to_numpy()
ct = pd.read_csv(cpctcqfile, sep='  ', skiprows=42, nrows=len(tsr), header=None, engine='python').to_numpy()
cq = pd.read_csv(cpctcqfile, sep='  ', skiprows=72, nrows=len(tsr), header=None, engine='python').to_numpy()


with open('/projects/shellwind/rthedin/templateFiles/CpCtCq_iea15mw.csv', 'w') as file:
    file.write('# IEA 15 MW\n')
    file.write('TSR,   Pitch ,  C_Fx  ,   C_Fy   ,  C_Fz   ,  C_Mx   ,  C_My  ,   C_Mz\n')
    file.write('(-),   (deg) ,  (-)  ,   (-)   ,  (-)   ,  (-)   ,  (-)  ,   (-)\n')

    for t in range(len(tsr)):
        for p in range(len(pitch)):
            cfx = cp[t,p]
            cmx = cq[t,p]
            file.write(f"{tsr[t]:.4f},\t{pitch[p]:.6f},\t{cfx:.6f},\t{0:.4f},\t{0:.4f},\t{cmx:.4f},\t{0:.4f},\t{0:.4f}\n")


