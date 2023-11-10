""" 
Illustrate how to track a wake center for a velocity profile at a given time.
A dummy velocity profile is created using a Gaussian wake and a power law shear.
Two wake tracking methods are illustrated.

"""
import numpy as np
import matplotlib.pyplot as plt
from samwich.waketrackers import track, WakeTracker
from samwich.dataloaders import PlanarData

# --- Inputs to generate a dummy velocity profile
U0  = 8     # Mean wind speed [m/s]
D   = 240   # Rotor diameter [m]
sig = 0.5*D # Wake extent for the dummy Gaussian wake [m]
HH  = 150   # Hub-height [m]
y0  = 120   # Horizontal wake center location [m]
z0  = HH+10 # Vertical wake center location [m]

# Grid where velocity is known in a given plane
y = np.arange(-500., 501., 5.) # shape: ny
z = np.arange(0., 501., 5.)    # shape: nz

# --- Create a dummy wake, centered (y0,z0) with a power law shear
Z, Y = np.meshgrid(z,y)                   # Z and Y have shape (ny x nz)
U = - (U0/2) * np.exp(-1/2 * ((Y-y0)**2 + (Z-z0)**2)/sig**2) # Gaussian wake deficit
u_shear = U0*(z/HH)**0.1 # Power law shear
U = U + u_shear  # Add shear

# --- Compute mean shear based on values at domain boundaries
shear_profile = (np.mean(U[:5,:],axis=0) + np.mean(U[-5:,:],axis=0))/2

# --- Format data for wake tracker
# NOTE: it's possible to provide multiple time steps, this is not done here
wakedata = PlanarData({'u':U, 'y':Y, 'z':Z}) # All inputs have shape (ny x nz)

# --- Use constant Area method to track the wake center
# Create a Constant Area wake tracker for that plane slice
wake = track(wakedata.sliceI(), method='ConstantArea', verbose=False)
# Remove shear to improve wake tracking
wake.remove_shear(wind_profile=shear_profile)
# Find the wake contour and center
yc, zc = wake.find_centers(ref_area=np.pi*D**2/4, weighted_center=lambda u: u**2) # y and z coordinates of the wake center
print('Calculated wake center:',yc, zc)
wake.plot_contour(vmin=-4, vmax=4, outline=True)

# --- Use Gaussian method to track the wake center
# Create a Gaussian wake tracker for that plane slice
wake = track(wakedata.sliceI(), method='Gaussian', verbose=False)
# Remove shear to improve wake tracking
wake.remove_shear(wind_profile=shear_profile)
# Find the wake contour and center
yc, zc = wake.find_centers(sigma=0.2*D, umin=None)
print('Calculated wake center:',yc, zc)
wake.plot_contour(vmin=-4, vmax=4, outline=True)

# Contour
contour = wake.paths[0]

plt.show()


