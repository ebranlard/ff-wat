import os
import numpy as np
import pandas as pd
import glob
# Local 
import xarray
from windtools.amrwind.post_processing import Sampling


def readPlanes(outBase, group, simCompleted=False, ITimes=None, rate=40, verbose=False):
    print('> outBase:', outBase)
    if ITimes is None:
        # Figure out ITimes from list of *.nc files
        n = len(os.path.basename(outBase))
        filepaths = glob.glob(outBase+'*.nc')
        filenames=[os.path.splitext(os.path.basename(f))[0] for f in filepaths]
        ITimes = [int(f[n:]) for f in filenames]
        if verbose:
            print('> ITimes sorted                    :', ITimes) 
            print('> ITimes diff                      :',  np.array(np.diff(ITimes)))
            print('> ITimes diff divided by rate      :',  np.array(np.diff(ITimes))/rate, ', rate:',rate)
        if len(ITimes)==0:
            raise Exception('No files found.')
    ITimes.sort()
    ncPaths = [outBase+'{}.nc'.format(itime) for itime in ITimes ] 
    for ii, (itime, ncPath) in enumerate(zip(ITimes,ncPaths)):
        if verbose:
            print(ii, '---------------------------------------------------------------------------')
        if not os.path.exists(ncPath):
            raise Exception('File not found '+ncPath)

        # --- Reading
        sp = Sampling(ncPath)
        # Looking at next "ITime" to see how much time was spent between the two
        iLast0 = None
        DeltaIT = np.nan
        if ii+1<len(ITimes):
            DeltaIT = ITimes[ii+1]-ITimes[ii]
            #iLast0 = int(DeltaIT/rate)-1 # NOTE: uncomment to try using ftime
        if iLast0 is None:
            print(ii, 'Reading: {}        (all)'.format(os.path.basename(ncPath), 0, iLast0))
            ds = sp.read_single_group(group, simCompleted=True, var=['u'])
        else:
            print(ii, 'Reading: {}        (from {} to {})'.format(os.path.basename(ncPath), 0, iLast0))
            ds = sp.read_single_group(group, itime=0, ftime=iLast0+1, var=['u'])
        # --- DEBUG
        IT = ds.samplingtimestep.values
        if ii+1<len(ITimes):
            iLast0 = int(DeltaIT/rate)-1
        else:
            iLast0 = len(IT)-1
        umean = ds.u.isel(x=0).mean(dim=['y','z']).values
        INaN = np.where(np.isnan(umean))[0]
        I0   = np.where(umean==0)[0]
        if len(INaN)>0:
            iLast = I0[-1]
        else:
            iLast = len(IT)-1
        if verbose:
            print(ii, '> Wind speed at x=0, IT=[0, -1] : [{}, {}]'.format(umean[0], umean[-1]))
            print(ii, '> Shape of field u              :', ds.u.shape)
            print(ii, '> Number of NaN                 :', len(INaN))
            print(ii, '> Number of zeros               :', len(I0))
            print(ii, '> Last iTime index valid        :', iLast)
            print(ii, '> Difference with next ITime    :', iLast0+1, ', multiplied by rate:', DeltaIT)
            print(ii, '> Current IT range              : [{}, {}],  length:{}'.format(IT[0], IT[-1], len(IT)))
        # --- Truncating to keep only the relevant IT (note: we could also use ftime)
        # TODO DEBUG and uncomment me
        # import pdb; pdb.set_trace()
        if iLast0<len(IT)-1:
            if verbose:
                print(ii, '> Truncating current IT to      : [{}, {}] {}'.format(IT[0], IT[iLast0], sMsg))
            ds = ds.isel(samplingtimestep=IT[:iLast0+1])
            IT = ds.samplingtimestep.values
            if verbose:
                print(ii, '> Current IT range (after trunc): [{}, {}],  length:{}'.format(IT[0], IT[-1],len(IT)))
        else:
            if verbose:
                print(ii, '> Not truncating')
        if ii>0:
            # Add maximum sampling index
            OLD_IT  = ds0.samplingtimestep.values
            iOffset = OLD_IT[-1] +1
            THIS_IT = ds.samplingtimestep.values
            NEW_IT  = THIS_IT + iOffset
            print(ii, '> Changing Current IT           : from  [{} {}] to [{} {}]  (adding {})'.format(THIS_IT[0], THIS_IT[-1], NEW_IT[0], NEW_IT[-1], iOffset))
            ds['samplingtimestep'] = ds['samplingtimestep'] + iOffset
            IT = ds.samplingtimestep.values
            print(ii, '> Concatenate previous [{} {}], with current: [{} {}]...'.format(OLD_IT[0], OLD_IT[-1], IT[0], IT[-1]))
            ds = xarray.concat((ds0,ds), 'samplingtimestep')
            # import pdb; pdb.set_trace()
        # We store for future concatenation 
        ds0 = ds
    ds = ds.rename_dims({'samplingtimestep':'it'})
    IT = ds.samplingtimestep.values
    dI = np.unique(np.diff(IT))
    if verbose:
        print('Final IT range                   : [{}, {}],  length:{}'.format(IT[0], IT[-1],len(IT)))
    if len(dI)!=1:
        print('Unique diff is:', dI)
        raise Exception('IT is not continuous')
    return ds


if __name__ == '__main__':

    #rate = 40
    #print('Unstable 0WT:')
    #ITimes = np.array([41159, 65159, 105159, 149159])/rate
    #print(np.diff(ITimes), np.sum(np.diff(ITimes)))
    #print('Unstable 1WT:')
    #ITimes =np.array([41159, 97159, 153159])/rate
    #print(np.diff(ITimes), np.sum(np.diff(ITimes)))
    #print('---------------------------------------------------------------------------')

    #Case = {'path':'02-iea15-2WT/neutral/'}

    caseDir = '/projects/tcwnd/rthedin/amr_runs/02_2turbine_coherence_inflow/'
    Case = {'stability':'unstable','nWT':0, 'path':os.path.join(caseDir, '03_0turbine_unstable.W.8at150.20dTInv_0.05q_0.75z0_850zi_10.24x5.12x0.96km_res5m_2ref')}

    planeBase = os.path.join(Case['path'], 'post_processing', 'planesT1')
    ds = readPlanes(planeBase, group='pT1', verbose=True)

