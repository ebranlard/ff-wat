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
from plothelper import *
from helper_functions import *

setFigureFont(13)

# --- Parameters
export=False
export=True

cmap='viridis'
D=240
u1=2
u2=8.1
stime='00250'

# --- 
clevels = np.linspace(u1, u2, 100)

f1 =weio.read(f'ff/_data/vtk_ff/FF-WAT-kd00_kg00.Low.DisXY01.{stime}.vtk')
f2 =weio.read(f'ff/_data/vtk_ff/FF-WAT-kd10_kg10-BrkDwn.Low.DisXY01.{stime}.vtk')
print(f1)


fig,axes = plt.subplots(1, 2, sharex=True, figsize=(12.8,2.3)) # (6.4,4.8)
fig.subplots_adjust(left=0.005, right=0.95, top=0.97, bottom=0.210, hspace=0.20, wspace=0.005)

cf1=axes[0].contourf(f1.xp_grid/D, f1.yp_grid/D, f1.point_data_grid['Velocity'][:,:,0,0].T, levels=clevels, cmap=cmap)
cf2=axes[1].contourf(f2.xp_grid/D, f2.yp_grid/D, f2.point_data_grid['Velocity'][:,:,0,0].T, levels=clevels, cmap=cmap)

for c in cf1.collections:
    c.set_rasterized(True)
for c in cf2.collections:
    c.set_rasterized(True)


for ax in np.array(axes).flatten():
#     ax.set_ylim([0,4.0])
    ax.axis('scaled')
    ax.set_xlim([-0.5,10])
    ax.set_xlabel(r'$y/D$ [-]')
    ax.set_ylabel(r'$z/D$ [-]')
    ax.tick_params(direction='in', top=True, right=True, labelright=False, labeltop=False, which='both')





# --- Colorbar
mappable=cf1
pmin=[1,1,1,1]
pmax=[0,0,0,0]
for ax in axes:
    chartBox = ax.get_position()
    p    = np.around(ax.get_position().get_points().flatten(),3) # BBox: x1, y1,  x2 y2
    pmin = np.around([min(v1, v2) for v1,v2 in zip(pmin, p)] ,3)
    pmax = np.around([max(v1, v2) for v1,v2 in zip(pmax, p)] ,3)
pad = 0.038
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
# cax.set_title()
cbar.set_label('Axial velocity u [m/s]', fontsize=12, rotation=-90, labelpad=18)

if export:
    fig.savefig('../article_torque/figs/FF-WakeNoTI.pdf', dpi=100) #, dpi=300, rasterized=True)
    fig.suptitle('FF-WakeNoTI')
    export2png()

plt.show()
