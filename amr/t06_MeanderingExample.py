import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatch
# Local 
import weio
from welib.essentials import *

from samwich.waketrackers import track, WakeTracker
from samwich.dataloaders import PlanarData


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
u1,u2 = U0-4, U0+4  # range for plot
s1,s2 = -4, 4       # range for plot when shear removed

# --- Create a dummy velocity profile with a wake at (y0,z0) and a power law shear
Z,Y = np.meshgrid(z,y)
U = u_gauss(Y, Z, U0, sig, z0=z0, y0=y0)
shear_in = U0*(z/HH)**0.1 # Power law
U = U + shear_in - U0 # Add shear

# --- Compute mean shear based on values at domain boundaries
shear_profile = (np.mean(U[:5,:],axis=0) + np.mean(U[-5:,:],axis=0))/2
dict4wake = {}
dict4wake['u'] = U
dict4wake['z'], dict4wake['y'] = np.meshgrid(z, y)
wakedata = PlanarData(dict4wake)

# --- Constant Area
ca_wake = track(wakedata.sliceI(), method='ConstantArea', verbose=False)
ca_wake.remove_shear(wind_profile=shear_profile)
ca_yc, ca_zc = ca_wake.find_centers(refArea, weighted_center=lambda u: u**2) # y and z coordinates of the wake center

# --- Gaussian
ga_wake = track(wakedata.sliceI(), method='Gaussian', verbose=False)
ga_wake.remove_shear(wind_profile=shear_profile)
ga_yc, ga_zc = ga_wake.find_centers(sigma=0.25*D, umin=None)


# --- Plot
colrs = [fColrs(4),fColrs(5)]
fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
clevels = np.linspace(u1, u2, 100)
cf = ax.contourf(Y, Z, U, clevels, cmap='viridis', extend='both')
ax.plot(y0,    z0,    'ko'            , ms=10, label='Input',             alpha=0.5, markeredgewidth=1, markeredgecolor='w')
ax.plot(ca_yc, ca_zc, 'd' , c=colrs[0], ms=8 , label = 'Contour method' , alpha=0.5, markeredgewidth=1, markeredgecolor='w')
ax.plot(ga_yc, ga_zc, 's' , c=colrs[1], ms=8 , label = 'Gaussian method', alpha=0.5, markeredgewidth=1, markeredgecolor='w')
ax.add_patch(mpatch.PathPatch(mpath.Path(ca_wake.paths[0]), lw=1,ls='-', facecolor='none', edgecolor=colrs[0]))
ax.add_patch(mpatch.PathPatch(mpath.Path(ga_wake.paths[0]), lw=1,ls='-', facecolor='none', edgecolor=colrs[1]))
ax.set_xlabel('y [m]')
ax.set_ylabel('z [m]')
ax.legend()
cb = fig.colorbar(cf, ticks=np.linspace(u1, u2, 11))
cb.set_label(label=r'$U$ [m/s]',fontsize=14)
cb.ax.tick_params(labelsize=12)
ax.set_xlim(np.min(y), np.max(y))
ax.set_ylim(np.min(z), np.max(z))
ax.axis('scaled')
ax.tick_params(axis='both', labelsize=12, size=10)
ax.set_xlabel(r'$y$ [m]', fontsize=14)
ax.set_ylabel(r'$z$ [m]', fontsize=14)

fig.savefig('./_figs/WakeTrackingExample.png')

# --- Plot shear
# fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
# fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
# ax.plot(shear_in, z, label='')
# ax.plot(shear_profile   , z, '--', label='')
# ax.set_xlabel('u')
# ax.set_ylabel('z')
# ax.legend()
# plt.show()



plt.show()


