import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatch
# Local 
import welib.weio as weio
from welib.essentials import *
import sys
sys.path.append('amr')
from helper_functions import *
from plothelper import *
# from t90_PostProMeanderingFrame import *

setFigureFont(12)

export=False
export=True


# --- Parameters
outDirTraj='amr/_out_all/trajectories/'
caseName = 'neutral2WT'

D=240

iWT = 1
iP  = 6
zOff= 250 # TODO TODO
u1    = 0
u2    = 11
cmap='viridis'
rasterized=True
mk    = 'o'
col   = 'w'

# it =600
it  = 812

iTimeMin=500
# --- derived parameters
Case=AllCases[caseName]
outDir = 'amr/_out_all/'
outDirHori = os.path.join(outDir  , 'lines')
outDirTraj   = os.path.join(outDir, 'trajectories')
U0, dt, D, xyWT, HubHeight, xPlanes = getSimParamsAMR(Case['stability'])
# --- Reading trajectories
trajFile = os.path.join(outDirTraj, '{}_MeanderWT{:d}_TrajectoriesM.csv'.format(caseName, iWT))
print('Reading: ',trajFile)
dfTrajM = weio.read(trajFile).toDataFrame()
dfTrajM.index = dfTrajM['Time_[s]'].astype(int)

# --- Read 
planeBase = os.path.join('amr', Case['path'], 'post_processing', 'planesT{}'.format(iWT))
with Timer('Reading...'):
    ds = readPlanes(planeBase, Case['planeTimes'], group='pT{}'.format(iWT))
    ds['z'] = ds.z -(np.max(ds.z)+np.min(ds.z))/2
    ds['y'] = ds.y-xyWT[iWT][1]
    ds['x'] = np.around((ds.x-xyWT[iWT][0])/D)
#print(ds)
print('ITime :',ds.it.values[0], ds.it.values[-1])

# --- Common variables useful for grid
z = ds.z.values
y = ds.y.values
Z, Y = np.meshgrid(ds.z.values, ds.y.values)

# Index coordinates of domain center (we take this as a reference
icz0= np.argmin(np.abs(z-0)) 
icy0= np.argmin(np.abs(y-0))
IY0 = np.arange(len(y))
IZ0 = np.arange(len(z))
iymax = len(y)
izmax = len(z)
dy = y[1]-y[0]
dz = z[1]-z[0]
IYBoundary = list(range(195,201))+list(range(0,5)) # TODO TODO TODO Not Generic
ITime = np.arange(iTimeMin, np.max(ds.it))
if len(ITime)>len(dfTrajM):
    print('[WARN] THIS SIMULATION HAS MORE ITIME THAN TRAJECTORY', len(ITime), len(dfTrajM))
    ITime = ITime[:len(dfTrajM)]
print('ITime :',ITime[0], ITime[-1])
# ---  Loop on planes

# Mean shear over full simulation 
shear = ds.isel(x=iP, y=IYBoundary).mean(dim=['it','y']).u.values
Yc = dfTrajM['y{}'.format(iP)]
Zc = dfTrajM['z{}'.format(iP)]

#if np.mod(it,100)==0:
#    print('iWT',iWT, 'iP',iP, 'it',it)
# --- Find wake center
U = ds.isel(x=iP, it=it).u.values.copy()
yc, zc = Yc[it], Zc[it]

# --- Move grid to meandering frame of reference (probably not done in the smartest way)
icz = np.argmin(np.abs(z-zc))
icy = np.argmin(np.abs(y-yc))
# Difference between centers in index space
diz = icz-icz0 
diy = icy-icy0
# - Shift data
U2 = np.zeros_like(U)*np.nan
IY_in_old = IY0+diy
IZ_in_old = IZ0+diz
bYOK = np.logical_and(IY_in_old>=0, IY_in_old<iymax)
bZOK = np.logical_and(IZ_in_old>=0, IZ_in_old<izmax)
IY_in_new = IY0[bYOK]
IZ_in_new = IZ0[bZOK]
U2[np.ix_( IY_in_new, IZ_in_new ) ] = U[np.ix_( IY_in_old[bYOK], IZ_in_old[bZOK] ) ]

# --- Replace data
ds['u'].loc[dict(x=iP, it=it)] = U2

# --- Plot
Z0 = Z.copy()+zOff

# --- Redo Wake center analysis for plot only
yc1, zc1, contour1, ax = track_wake_center_plane(Y, Z0, U, D, method='ConstantArea', shear=shear, plot=False       )
yc2, zc2, contour2, ax = track_wake_center_plane(Y, Z0, U, D, method='Gaussian'    , shear=shear, plot=False, ax=ax, col=fColrs(5), mk='d')
shear2 = np.zeros_like(shear)*np.nan
shear2[IZ_in_new] = shear[IZ_in_old[bZOK]]
yc1m, zc1m, contour1m, ax = track_wake_center_plane(Y, Z, U2, D, method='ConstantArea', shear=shear2, plot=False       )
yc2m, zc2m, contour2m, ax = track_wake_center_plane(Y, Z, U2, D, method='Gaussian'    , shear=shear2, plot=False , ax=ax, col=fColrs(5), mk='d')

# --------------------------------------------------------------------------------}
# --- FIGURE 
# --------------------------------------------------------------------------------{
fig,axes = plt.subplots(1, 2, sharex=True, figsize=(12.8,3.2)) # (6.4,4.8)
fig.subplots_adjust(left=0.070, right=0.92, top=0.95, bottom=0.12, hspace=0.30, wspace=0.20)
# fig.subplots_adjust(left=0.135, right=0.99, top=0.98, bottom=0.11, hspace=0.08, wspace=0.20)
clevels = np.linspace(u1, u2, 100)

# --- Plot in Inertial frame
ax=axes[0]
cf = ax.contourf(Y/D, Z0/D, U, levels=clevels, cmap=cmap)
for c in cf.collections:
    c.set_rasterized(rasterized)

ax.set_xlim(np.min(Y.flatten()) /D, np.max(Y.flatten()) /D)
ax.set_ylim(np.min(Z0.flatten())/D, np.max(Z0.flatten())/D)
ax.axis('scaled')
ax.set_xlabel(r'$y/D$ [-]')
ax.set_ylabel(r'$z/D$ [-]')
Colrs = [fColrs(1), fColrs(4)]

ax.plot(yc1/D,  zc1    /D, marker='o' , c=Colrs[0], ms=8 , label = 'Contour fit', alpha=0.6, markeredgewidth=1, markeredgecolor='k')
ax.plot(yc2/D,  zc2    /D, marker='d' , c=Colrs[1], ms=8 , label = 'Gaussian fit'    , alpha=0.6, markeredgewidth=1, markeredgecolor='k')
# ax.plot(yc /D, (zc+zOff)/D, marker='x'  , c='k'      , ms=10, label = 'Clean', alpha=1.0, markeredgewidth=1, markeredgecolor='k')
ax.add_patch(mpatch.PathPatch(mpath.Path(contour1/D), lw=1,ls='-', facecolor='none', edgecolor=Colrs[0]))
ax.add_patch(mpatch.PathPatch(mpath.Path(contour2/D), lw=1,ls='-', facecolor='none', edgecolor=Colrs[1]))
ax.legend()
ax.set_title('Global frame')

# --- Plot in MFR frame
ax=axes[1]
# ax.set_facecolor((0.8,0.8,0.8))

# Hatches and NA background
ax.set_xlim(np.min(Y.flatten())/D, np.max(Y.flatten())/D)
ax.set_ylim(np.min(Z.flatten())/D, np.max(Z.flatten())/D)
xlim = ax.get_xlim()
ylim = ax.get_ylim()
X2, Y2 = np.meshgrid(np.linspace(xlim[0], xlim[1], 10), np.linspace(ylim[0], ylim[1], 10))
ax.contourf(X2, Y2, np.ones((10, 10))*1.00, hatches="//",  cmap='bone', levels=np.linspace(0,1,10))


cf = ax.contourf(Y/D, Z/D, U2, levels=clevels, cmap=cmap)
for c in cf.collections:
    c.set_rasterized(rasterized)
ax.axis('scaled')
ax.set_xlabel(r'$y/D$ [-]')
ax.set_ylabel(r'$z/D$ [-]')
ax.plot(yc1m/D, zc1m/D, marker='o' , c=Colrs[0]      , ms=8 , label = 'Contour', alpha=0.6, markeredgewidth=1, markeredgecolor='k')
ax.plot(yc2m/D, zc2m/D, marker='d' , c=Colrs[1], ms=8 , label = 'Gaussian'    , alpha=0.6, markeredgewidth=1, markeredgecolor='k')
# ax.plot(0     , 0     , marker='x' , c='k'      , ms=10, label = 'Clean'       , alpha=1.0, markeredgewidth=1, markeredgecolor='k')
ax.add_patch(mpatch.PathPatch(mpath.Path(contour1m/D), lw=1,ls='-', facecolor='none', edgecolor=Colrs[0]))
ax.add_patch(mpatch.PathPatch(mpath.Path(contour2m/D), lw=1,ls='-', facecolor='none', edgecolor=Colrs[1]))
ax.set_title('Meandering frame')



for ax in np.array(axes).flatten():
    ax.tick_params(direction='in', top=True, right=True, labelright=False, labeltop=False, which='both')
    ax.set_xticks([-2,-1,0,1,2])


# --- COLORBAR
mappable = cf

pmin=[1,1,1,1]
pmax=[0,0,0,0]
for ax in axes:
    chartBox = ax.get_position()
    p    = np.around(ax.get_position().get_points().flatten(),3) # BBox: x1, y1,  x2 y2
    pmin = np.around([min(v1, v2) for v1,v2 in zip(pmin, p)] ,3)
    pmax = np.around([max(v1, v2) for v1,v2 in zip(pmax, p)] ,3)
pad = 0.015
w = 0.02
plots_height = pmax[3]-pmin[1]
plots_width  = pmax[2]-pmin[0]
# BoundingBox: Left (x), Bottom (y), Width, Height
BB = [pmax[2]+pad, pmin[1], w, plots_height]
cax = fig.add_axes(BB)


# from mpl_toolkits.axes_grid1 import make_axes_locatable
# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="6%", pad="2%")




cbar = fig.colorbar(mappable, cax=cax, ticks=np.arange(u1, u2+0.1, 2))
# cb = fig.colorbar(cf, )





#cb = fig.colorbar(cf, ticks=np.linspace(u1, u2, 11))
#cb.set_label(label=r'$U$ [m/s]',fontsize=14)
#cb.ax.tick_params(labelsize=12)


# fig.suptitle('WakeID')








if export:
    fig.savefig('../article_torque/figs/WakeID.pdf', dpi=100) #, dpi=300, rasterized=True)
#     export2pdf()


plt.show()
