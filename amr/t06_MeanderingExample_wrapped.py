
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import weio
from helper_functions import track_wake_center_plane


def u_gauss(y, z, U0, sig, udef=4,  z0=0, y0=0):
    """ Analytical Gaussian velocity profile u(y) = U - udef exp(-1/2 (r/sig)**2) """
    u = U0 - udef * np.exp(-1/2 * ((y-y0)**2 + (z-z0)**2)/sig**2)
    return u

U0=8
D=240
sig=0.5*D
HH = 150
y0 = 120
z0 = HH+20

y = np.arange(-500., 501., 5.) # 201
z = np.arange(0., 501., 5.) # 101

# --- Derived parameters
refArea = np.pi*D**2/4

# --- Create a dummy velocity profile with a wake at (y0,z0) and a power law shear
Z,Y = np.meshgrid(z,y)
U = u_gauss(Y, Z, U0, sig, z0=z0, y0=y0)
shear_in = U0*(z/HH)**0.1 # Power law
U = U + shear_in - U0 # Add shear

# --- Compute mean shear based on values at domain boundaries
shear = (np.mean(U[:5,:],axis=0) + np.mean(U[-5:,:],axis=0))/2

yc, zc, contour, ax = track_wake_center_plane(Y, Z, U, D, method='ConstantArea', shear=shear, plot=True)
yc, zc, contour, ax = track_wake_center_plane(Y, Z, U, D, method='Gaussian', shear=shear, plot=True)




plt.show()


