import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
from weio.vtk_file import VTKFile

def colorbar(mappable, pad=0.1, size='10%'):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size=size, pad=pad)
    return fig.colorbar(mappable, cax=cax)
    
def vtkIndices(folder, base):
    """ """
    pattern = os.path.join(folder,'vtk_ff', base+'*.vtk')
    files = glob.glob(pattern)
    if len(files)==0:
        raise Exception('No file found with pattern',pattern)
    sIndices = np.array([ os.path.splitext(os.path.splitext(f)[0])[1][1:] for f in files ])
    nDigits = len(sIndices[0]) # Number of characters used for time stamp
    Indices = sIndices.astype(int)
    iMax= np.max(Indices)
    Indices = np.sort(Indices)
    sFmt = '{:0'+str(nDigits)+'d}'
    return Indices, iMax, nDigits, sFmt


def extractFFPlane(kind, folder, simbase, iPlane=1, iTime=None, removeBoundaries=None, U0=1, D=1, verbose=False, sFmt=None):
    """ 

    kind: 'XY', 'XZ', or 'YZ'

    OPTIONAL INPUTS:
     - U0: if provided, used to scale the velocity field
     - D:  if provided, used to scale the coordinates
     - removeBoundaries: number of index to remove on the boundaries

    """
    base = simbase + '.Low.Dis{}{:02d}.'.format(kind,iPlane)
    # Default argument
    if sFmt is None or iTime is None:
        I, iMax, nDigits, sFmt = vtkIndices(folder, base)
    if iTime is None:
        iTime= np.max(I)
    vtkFile = os.path.join(folder, 'vtk_ff', base +sFmt.format(iTime)+'.vtk')

    if verbose:
        print('Reading: ',vtkFile)
    vtkc = VTKFile(vtkFile)
    x  = vtkc.xp_grid 
    if kind=='YZ':
        y  = vtkc.yp_grid/D
        z  = vtkc.zp_grid/D
        Uc = vtkc.point_data_grid['Velocity']/U0
    else:
        raise NotImplementedError()
    U = Uc[0, :, :, 0] # ny x nz
    V = Uc[0, :, :, 1] # ny x nz
    W = Uc[0, :, :, 2] # ny x nz

#     if removeBoundaries is not None:
#         i=removeBoundaries
#         y  = y[i:-i]
#         z  = z[i:-i]
#         U = U[i:-i, i:-i]
#         V = V[i:-i, i:-i]
#         W = W[i:-i, i:-i]
    if kind=='YZ':
        Y,Z = np.meshgrid(y,z)
        return y,z,Y,Z,U,V,W



def loadAllPlanes(folder, simbase, nPlanes, ISel=None, kind='YZ', sFmt=None, outDir=None, verbose=False):
    print('>>> Processing', folder, simbase)
    if ISel is None or sFmt is None:
        I, iMax, nDigits, sFmt = vtkIndices(folder, simbase+'.Low.Dis{}01'.format(kind))
        if ISel is None:
            ISel=I

    # Read first plane for dimensions
    y,z,Y,Z,U,V,W =  extractFFPlane(kind, folder, simbase, iPlane=1, iTime=ISel[0], verbose=verbose, sFmt=sFmt)

    # 
    xPlanes = list(np.arange(nPlanes))
    dims=['y','z','ix','it']
    coords = {'y':y, 'z':z, 'ix':xPlanes, 'it': ISel}
    ds = xr.Dataset({'u': xr.DataArray( dims = dims, coords=coords)})
    for iPlane in range(nPlanes):
        for iTime in ISel:
            try:
                y,z,Y,Z,U,V,W =  extractFFPlane(kind, folder, simbase, iPlane=iPlane+1, iTime=iTime, verbose=verbose, sFmt=sFmt)
            except:
                print('FAIL:',iTime,iPlane)
            ds.u.loc[:,:, iPlane, iTime]=U

    if outDir is not None:
        fileOut = os.path.join(outDir, 'planes_'+simbase+'.nc')
        print('>>> fileOut', fileOut)
        ds.to_netcdf(fileOut)
    return ds
