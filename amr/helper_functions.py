import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import xarray
# External
from windtools.amrwind.post_processing  import Sampling
# Local
from welib.tools.curve_fitting import model_fit
from welib.essentials import *

AllCases={}

# caseDir='./'
# AllCases['stable2WT']   = {'stability':'stable'  ,'nWT':2, 'path':os.path.join(caseDir, '02-iea15-2WT/stable/')}
# AllCases['neutral2WT']  = {'stability':'neutral' ,'nWT':2, 'path':os.path.join(caseDir, '02-iea15-2WT/neutral/')}
# AllCases['unstable2WT'] = {'stability':'unstable','nWT':2, 'path':os.path.join(caseDir, '02-iea15-2WT/unstable')}
# 
# AllCases['stable1WT']   = {'stability':'stable'  ,'nWT':1, 'path':os.path.join(caseDir, '02-iea15-1WT/stable/')}
# AllCases['neutral1WT']  = {'stability':'neutral' ,'nWT':1, 'path':os.path.join(caseDir, '02-iea15-1WT/neutral/')}
# AllCases['unstable1WT'] = {'stability':'unstable','nWT':1, 'path':os.path.join(caseDir, '02-iea15-1WT/unstable/')}
# 
# AllCases['stable0WT']   = {'stability':'stable'  ,'nWT':0, 'path':os.path.join(caseDir, '02-iea15-0WT/stable/')}
# AllCases['neutral0WT']  = {'stability':'neutral' ,'nWT':0, 'path':os.path.join(caseDir, '02-iea15-0WT/neutral/')}
# AllCases['unstable0WT'] = {'stability':'unstable','nWT':0, 'path':os.path.join(caseDir, '02-iea15-0WT/unstable/')}


caseDir = '/projects/tcwnd/rthedin/amr_runs/02_2turbine_coherence_inflow/'
AllCases['stable2WT']   = {'stability':'stable'  ,'nWT':2, 'path':os.path.join(caseDir, '01_2turbine_stable.W.8at150.20dTInv_0.25cooling_0.05z0_450zi_3.84x1.28x0.9km_res2.5m_coriolis3days_1ref_postFix2/')}
AllCases['neutral2WT']  = {'stability':'neutral' ,'nWT':2, 'path':os.path.join(caseDir, '02_2turbine_neutral.W.8at150.10dTInv_0.75z0_750zi_3.84x1.28x0.9km_res2.5m_1ref_postFix2/')}
AllCases['unstable2WT'] = {'stability':'unstable','nWT':2, 'path':os.path.join(caseDir, '03_2turbine_unstable.W.8at150.20dTInv_0.05q_0.75z0_850zi_10.24x5.12x0.96km_res5m_2ref_postFix2/')}

AllCases['stable1WT']   = {'stability':'stable'  ,'nWT':1, 'path':os.path.join(caseDir, '01_1turbine_stable.W.8at150.20dTInv_0.25cooling_0.05z0_450zi_3.84x1.28x0.9km_res2.5m_coriol')}
AllCases['neutral1WT']  = {'stability':'neutral' ,'nWT':1, 'path':os.path.join(caseDir, '02_1turbine_neutral.W.8at150.10dTInv_0.75z0_750zi_3.84x1.28x0.9km_res2.5m_1ref_postFix2/')}
AllCases['unstable1WT'] = {'stability':'unstable','nWT':1, 'path':os.path.join(caseDir, '03_1turbine_unstable.W.8at150.20dTInv_0.05q_0.75z0_850zi_10.24x5.12x0.96km_res5m_2ref_postFix2/')}

AllCases['stable0WT']   = {'stability':'stable'  ,'nWT':0, 'path':os.path.join(caseDir, '01_0turbine_stable.W.8at150.20dTInv_0.25cooling_0.05z0_450zi_3.84x1.28x0.9km_res2.5m_coriol')}
AllCases['neutral0WT']  = {'stability':'neutral' ,'nWT':0, 'path':os.path.join(caseDir, '02_0turbine_neutral.W.8at150.10dTInv_0.75z0_750zi_3.84x1.28x0.9km_res2.5m_1ref_postFix2/')}
AllCases['unstable0WT'] = {'stability':'unstable','nWT':0, 'path':os.path.join(caseDir, '03_0turbine_unstable.W.8at150.20dTInv_0.05q_0.75z0_850zi_10.24x5.12x0.96km_res5m_2ref/')}


def getSimParamsAMR(stability):
    dt=1
    xyWT=dict()
    xyWT[1]         = (706.25 , 641.25)   # Actuator.T1.base_position = 713.28 641.25 0    # hub at 701.25, 641.25, 0, considering overhang of 12.0313
    xyWT[2]         = (2386.25, 641.25)   # Actuator.T2.base_position = 2393.28 641.25 0.  # hub at 2381.25, 641.25, 0, considering overhang of 12.0313
    D             = 240
    xPlanes = [0,1,2,3,4,5,6]
    if stability == 'neutral':
        U0=8.5;
    elif stability == 'stable': 
        U0=8.1;
    elif stability == 'unstable': 
        U0=7.23;
    # fixed_dt =   0.025
    # output_frequency = 40  # every 1 s
    # # R             = 120.97
    # OverHang      = -12.097571763912535
    # ShftTilt      = -6.0 
    # Twr2Shft      = 4.349459414248071
    # TowerHt       = 144.386          
    #HubHeight  =  TowerHt + Twr2Shft + OverHang*np.sin(ShftTilt*np.pi/180)
    HubHeight=150
    return U0, dt, D, xyWT, HubHeight, xPlanes


def eddyFilter(x_plane, Dmin, Dmax, Fmin, Exp):
    """ NOTE: x is x/D [-] """
    F = np.zeros_like(x_plane).astype(float)
    bLow = x_plane<= Dmin
    bHigh = x_plane >= Dmax
    F[bLow] = Fmin
    F[bHigh] = 1 
    b =  ~(np.logical_or(bLow, bHigh))
    F[b] =  Fmin + (1-Fmin)*( ( (x_plane[b]) - Dmin ) / (Dmax-Dmin) )**Exp
    return F

eddyFilter_Amb = {'Dmin':0, 'Dmax':1 , 'Fmin':1  , 'Exp':0.01}
eddyFilter_Shr = {'Dmin':3, 'Dmax':25, 'Fmin':0.2, 'Exp':0.1}



def readPlanes(outBase, group, simCompleted=False, ITimes=None):
    if ITimes is None:
        # Figure out ITimes from list of *.nc files
        n = len(os.path.basename(outBase))
        filepaths = glob.glob(outBase+'*.nc')
        filenames=[os.path.splitext(os.path.basename(f))[0] for f in filepaths]
        ITimes = [int(f[n:]) for f in filenames]
    ITimes.sort()
    ncPaths = [outBase+'{}.nc'.format(itime) for itime in ITimes ] 
    for ii, (itime, ncPath) in enumerate(zip(ITimes,ncPaths)):
        if not os.path.exists(ncPath):
            raise Exception('File not found '+ncPath)
        print('Reading: ',ncPath)
        sp = Sampling(ncPath)
        ds = sp.read_single_group(group, simCompleted=True, var=['u'])
        if ii>0:
            # Add maximum sampling index
            ds['samplingtimestep'] = ds['samplingtimestep'] + ds0.samplingtimestep.values[-1]+1
            print('Concat...(extremely slow)')
            ds = xarray.concat((ds0,ds), 'samplingtimestep')
        # We store for future concatenation 
        ds0 = ds
    ds = ds.rename_dims({'samplingtimestep':'it'})
    return ds


def readDataSet(filename):
    if not os.path.exists(filename):
        raise Exception('File not found:'+ filename)
    print('Reading:', filename)
    try:
        return xarray.open_dataset(filename)
    except Exception as e:
        FAIL('Failed to read:'+ filename)
        raise e

def getLabel(caseName, Meander, symmetric, removeBG, smooth):
    prefix = '_Meander_' if Meander else '_Inertial_'
    label='KFit_'+caseName+prefix
    if symmetric:
        label+='_sym'
    label+='_rm'+removeBG+'_smooth{:d}'.format(smooth)
    return label

class FailSafe(object):
    """ Run a command """
    def __init__(self, msg=None, raiseException=True):
        self.msg = msg
        self.raiseException = raiseException
    def __enter__(self):
        pass

    def __exit__(self, exc_type, exc_val, traceback):
        if exc_type is not None:
            # An exception occured
            FAIL(self.msg+': ' + exc_val.__repr__())
            return not self.raiseException
        else:
            OK(self.msg)
            return True



def FAIL(msg):
    HEADER = '\033[95m'
    RED = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    print(RED+'[FAIL] ' + msg + ENDC)

def WARN(msg):
    ORAN = '\033[93m'
    ENDC = '\033[0m'
    print(ORAN+'[WARN] ' + msg + ENDC)

def OK(msg):
    GREEN = '\033[92m'
    ENDC = '\033[0m'
    print(GREEN+'[ OK ] ' + msg + ENDC)

def track_wake_center_plane(Y, Z, U, D, method='ConstantArea', plot=False, removeShear=True, shear=None, 
        ax=None, u1=0, u2=10, mk='o', col='w' # Plot options
        ): 
    """ 
    INPUTS:
     - U shape: ny x nz
     - Y shape: ny x nz, e.g.: Z, Y  = np.meshgrid(z, y)
     - Z shape: ny x nz 
     - D: diameter
     - method: string in ['ConstantArea', 'Gaussian'] for now
     - shear: shape nz. 
    """
    from samwich.waketrackers import track, WakeTracker
    from samwich.dataloaders import PlanarData
    # --- Sanity
    ny, nz = U.shape
    assert(Y.shape == (ny, nz))
    assert(Z.shape == (ny, nz))

    # --- Compute shear profile
    if removeShear is None:
        if shear is None:
            # Compute mean shear based on values at domain boundaries
            shear = (np.mean(U[:5,:],axis=0) + np.mean(U[-5:,:],axis=0))/2
        assert(shear.shape == nz)


    # --- Format data for wake tracker
    # NOTE: it's possible to provide multiple time steps, this is not done here
    wakedata = PlanarData({'u':U, 'y':Y, 'z':Z}) # All inputs have shape (ny x nz)

    # Create a  wake tracker for that plane slice
    wake = track(wakedata.sliceI(), method=method, verbose=False)
    # Remove shear to improve wake tracking
    if removeShear:
        wake.remove_shear(wind_profile=shear)
    # Find the wake contour and center
    if method=='ConstantArea':
        yc, zc = wake.find_centers(ref_area=np.pi*D**2/4, weighted_center=lambda u: u**2) # y and z coordinates of the wake center
    elif method=='Gaussian':
        yc, zc = wake.find_centers(sigma=0.25*D, umin=None)
    else:
        raise NotImplementedError(method)

    yc=yc[0]
    zc=zc[0]
    contour = wake.paths[0]

    # --- Plot
    if plot:
        import matplotlib.path as mpath
        import matplotlib.patches as mpatch
        # wake.plot_contour(vmin=s1, vmax=s2, outline=True)
        # ga_wake.plot_contour(vmin=s1, vmax=s2, outline=True)
        newFig = ax is None
        if newFig:
            fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
            clevels = np.linspace(u1, u2, 100)
            cf = ax.contourf(Y, Z, U, clevels, cmap='viridis', extend='both')
            cb = fig.colorbar(cf, ticks=np.linspace(u1, u2, 11))
            cb.set_label(label=r'$U$ [m/s]',fontsize=14)
            cb.ax.tick_params(labelsize=12)

            ax.set_xlim(np.min(Y.flatten()), np.max(Y.flatten()))
            ax.set_ylim(np.min(Z.flatten()), np.max(Z.flatten()))
            ax.axis('scaled')
            ax.tick_params(axis='both', labelsize=12, size=10)
            ax.set_xlabel(r'$y$ [m]', fontsize=14)
            ax.set_ylabel(r'$z$ [m]', fontsize=14)

        ax.plot(yc, zc, marker=mk , c=col, ms=8 , label = method, alpha=0.5, markeredgewidth=1, markeredgecolor='w')
        ax.add_patch(mpatch.PathPatch(mpath.Path(contour), lw=1,ls='-', facecolor='none', edgecolor=col))
        #ax.legend()
        # --- Plot shear
        # fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        # fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        # ax.plot(shear_in, z, label='')
        # ax.plot(shear_profile   , z, '--', label='')
        # ax.set_xlabel('u')
        # ax.set_ylabel('z')
        # ax.legend()
        # plt.show()
    else:
        ax=None
    
    return yc, zc, contour, ax


def fitK(y, u, var, U0, D):
    # Precompute Gradient 
    K, KD, KG, du, gu = kWAT1D(y, U0, u, D, kd=0.1, kg=0.1)
#     uG, duG, guG, KG, KDG, KGG, fitter = fitGaussianAndComputeK(y, u, U0, D, kDef, kGrad)
    if np.isnan(y).sum()>0:
        raise Exception('>>>> NAN in y')
    if np.isnan(y).sum()>0:
        raise Exception('>>>> NAN in u')
    if np.isnan(var).sum()>0:
        raise Exception('>>>> NAN in var')


    def K_model(y, p): #, U0, u, D):
        kDef = p[0]
        kGrad = p[1]
        K, KD, KG =  kWAT(du, gu, U0, D, kDef, kGrad)
        return K

    k_fit, ks_fit, fitter = model_fit(K_model, y, np.sqrt(var), bounds={'kd':(0, 5),'kg':(0, 15)}, p0={'kd':1.0,'kg':1.0}) #, U0=U0, u=u, D=D)
    #k_fit, ks_fit, fitter = model_fit(K_model, y, np.sqrt(var), bounds={'kd':(0, 0.05),'kg':(0, 5)}, p0={'kd':0.0,'kg':1.0}) #, U0=U0, u=u, D=D)
    return k_fit, ks_fit, fitter


def ds_fitK(ds, iP, ITime, D, U0, symmetric=True, removeBG='None', ds0=None, smooth=0):
    if removeBG!='0WT':
        ds0=None
    var = compute_var(ds, iP, ITime, D, removeBG=removeBG, ysym=symmetric, BGlim=1.9, ds0=ds0)
    u   = compute_u  (ds, iP, ITime, ysym=symmetric)
    y   = var.y.values
    u   = u.values
    var = var.values
    bnan = np.isnan(var)
#     bnan2 = np.isnan(u)
#     nnan=np.sum(bnan)
#     nnan2=np.sum(bnan2)
    u   = u[~bnan]
    y   = y[~bnan]
    var = var[~bnan]

    from scipy.signal import savgol_filter

    # Apply the Savitzky-Golay filter
    if smooth==0:
        pass
    else:
        u_smoothed   = savgol_filter(u,   window_length=smooth*2, polyorder=1)
        var_smoothed = savgol_filter(var, window_length=smooth*2, polyorder=1)
        smoothLength = y[0:smooth*2]
        smoothLength = smoothLength[-1]-smoothLength[0]
        print('>>> Smooth Length', smoothLength, smoothLength/D)
#         fig=plt.figure()
#         plt.plot(u,            y, label='Original')
#         plt.plot(u_smoothed  , y, label='Smoothed')
#         plt.plot(var,          y, label='Original')
#         plt.plot(var_smoothed, y, label='Smoothed')
#         fig.suptitle('iP'+str(iP))
        u   = u_smoothed
        var = var_smoothed

    var[var<0]=0 # Smoothing might introduce negative values..
    assert(u.shape  ==y.shape)
    assert(var.shape==y.shape)

    if np.isnan(u).sum()>0:
        import pdb; pdb.set_trace()
    if np.isnan(var).sum()>0:
        import pdb; pdb.set_trace()
    if np.isnan(y).sum()>0:
        import pdb; pdb.set_trace()


    if len(u)==0:
        raise Exception('>>>>>>>>>>>>>>>>>>>> ALL VALUES ARE NAN for iP',iP, 'Check that varbg is not NAN')
#         print('>>>>>>>>>>>>>>>>>>>> ALL VALUES ARE NAN for iP',iP)
#         k_fit=[]
#         ks_fit=[np.nan, np.nan]
#         fitter=None
#     else:
    k_fit, ks_fit, fitter =  fitK(y, u, var, U0, D)
    return y, var, u, k_fit, ks_fit, fitter



def ds_ysym(ds, y0=0):
    """ Make ds symmetric about y"""
    y = ds.y.values-y0
    iy0 = np.argmin(np.abs(y-0    ))
    IP = np.arange(iy0,len(y))
    IN = np.arange(iy0,-1,-1)
    assert( all(y[IP]==-y[IN]) )
    ds2 = ds.isel(y=IP)
    ds2['y'] = y[IP]
    # NOTE: we concatenate and use a nan mean
    if isinstance(ds, xarray.DataArray):
        #mean1 = (ds.values[IP] + ds.values[IN])/2
        #mean2 = np.nanmean( np.stack( (ds.values[IP], ds.values[IN])), axis=0)
        ds2.loc[:] = np.nanmean( np.stack( (ds.values[IP], ds.values[IN])), axis=0)
    else:
        for k in list(ds.keys()):
            dims = ds[k].dims
            print('>>> TODO Deal with nan mean')
            import pdb; pdb.set_trace()
            ds2[k] = (   dims , (ds[k].values[IP]+ds[k].values[IN])/2  )
    return ds2

def u_gauss(y, U, ui, D, sig, kd=0.5, kg=0.5, use_numpy=True):
    """ 
    Analytical Gaussian velocity profile u(y) = U - ui exp(-1/2 (y/sig)**2)
    du = u-U # Deficit (negative)
    gu = gradient
    """
    if use_numpy:
        du = - ui * np.exp(-1/2 * (y/sig)**2)
    else:
        du = - ui *    exp(-1/2 * (y/sig)**2)
    u = U + du
    gu = - du * y / sig**2 
    KD = kd/U * abs(du)
    KG = kg*D/(2*U) * abs(gu)
    K = KD + KG
    K = (1/2*D*kg * abs(y) / sig**2 + kd)/U * abs(du)
    
    return u, du, gu, K, KD, KG


def kWAT(du, gu, U0, D, kd, kg):
    KD = kd/U0 * abs(du)
    KG = kg*D/(2*U0) * abs(gu)
    K = KD + KG
    return K, KD, KG


def kWAT1D(y, U, u, D, kd, kg, use_numpy=True):
    du = u-U
    if use_numpy:
        gu = np.gradient(du, y)
    else:
        gu = u.diff(y)        
    K, KD, KG =  kWAT(du, gu, U, D, kd, kg)
    return K, KD, KG, du, gu


def fitGaussianAndComputeK(y, u, U0, D, kDef, kGrad):
    bounds={'ui':(0.001,U0), 'sig':(0.001*D, 2*D)}
    guess ={'ui':U0/2, 'sig':D}
    u_fit, pfit, fitter = model_fit('eval: '+str(U0)+' - {ui}*np.exp(-1/2*(x/{sig})**2)', y, u, p0=guess, bounds=bounds) 
    c = fitter.model['coeffs']
    uG, duG, guG, KG, KDG, KGG = u_gauss(y, U0, ui=c['ui'], D=D, sig=c['sig'], kd=kDef, kg=kGrad)
    return  uG, duG, guG, KG, KDG, KGG, fitter

def fitGaussian(y, u, U0, D):
    y = np.asarray(y.copy())
    u = np.asarray(u.copy())
    bounds={'ui':(0.01,U0), 'sig':(0.001*D, 4*D)}
    guess ={'ui':U0/2, 'sig':D}
    u_fit, pfit, fitter = model_fit('eval: '+str(U0)+' - {ui}*np.exp(-1/2*(x/{sig})**2)', y, u, p0=guess, bounds=bounds) 
    c = fitter.model['coeffs']
    return u_fit, c['ui'], c['sig'], fitter

def fitGaussian01(y, u, U0, D):
    y = np.asarray(y.copy())
    u = np.asarray(u.copy())
    u[np.isnan(u)]=U0
    umin =np.min(u)
    umax =np.max(u)
    u= (u-umin)/(umax-umin)
    bounds={'sig':(0.001*D, 4*D)}
    guess ={'sig':D}
    u_fit, pfit, fitter = model_fit('eval: 1 - np.exp(-1/2*(x/{sig})**2)', y, u, p0=guess, bounds=bounds) 
    c = fitter.model['coeffs']
    return u_fit, c['sig'], fitter






def common_itime(*dss):
    ITime = set(dss[0].it.values)
    for ds in dss:
        ITime = ITime.intersection(ds.it.values)
    ITime = list(ITime)
    ITime.sort()
    return np.asarray(ITime)




def compute_var(ds, iP, ITime, D, removeBG='None', ysym=False, BGlim=1.2, ds0=None, zeroLim=None):

    # First compute variance
    var = ds.isel(x=iP, it=ITime).var(dim='it', skipna=True).u
    if ds0 is not None:
        var0 = ds0.isel(x=iP, it=ITime).var(dim='it', skipna=True).u

    # Then do symmetry (order important!)
    if ysym:
        var = ds_ysym(var)
        if ds0 is not None:
            var0 = ds_ysym(var0)
    if removeBG=='Zero':
        varbg =0
    elif removeBG=='0WT':
        varbg = var0
        varbg[np.isnan(varbg)] = np.nanmean(varbg)
    elif removeBG=='Outer':
        varbg = var.isel(y=np.abs(var.y)>BGlim*D).mean(skipna=True).values
    elif removeBG=='Min':
        varbg = np.min(var)
    else:
        raise NotImplementedError(removeBG)

    if np.asarray(np.isnan(varbg)).sum()>0:
        print('>>> VARBG is NaN')
        import pdb; pdb.set_trace()

    var -= varbg

    # Deal with negative values
    #var[var<0] = -var[var<0]
    var[var<0] = 0

    # ---
    if zeroLim is not None:
        var[abs(var.y)>zeroLim*D] = 0
    return var

def compute_u(ds, iP, ITime, ysym=True, zeroLim=None, D=None, U0=0):
    u = ds.isel(x = iP, it = ITime).mean(dim = 'it', skipna=True).u
    if ysym:
        u = ds_ysym(u)
    if zeroLim is not None:
        u[np.abs(u.y)>zeroLim*D] = U0
    return u


def compute_tke(u, v, w, axis=0):
    """
    Computes the turbulent kinetic energy (TKE) from the three velocity components.

    INPUTS:
    u (array-like): The wind component along the x-axis.
    v (array-like): The wind component along the y-axis.
    w (array-like): The wind component along the z-axis.

    Returns:
    float: The TKE.
    """
    tke = 0.5 * (np.var(u, axis=axis) + np.var(v, axis=axis) + np.var(w, axis=axis))

#     # Compute the mean velocity components.
#     u_bar = np.mean(u, axis=axis)
#     v_bar = np.mean(v, axis=axis)
#     w_bar = np.mean(w, axis=axis)
# 
#     # Compute the fluctuations from the mean.
#     u_prime = u - u_bar
#     v_prime = v - v_bar
#     w_prime = w - w_bar
# 
#     # Compute the TKE.
#     tke = 0.5 * (np.mean(u_prime ** 2, axis=axis) + np.mean(v_prime ** 2, axis=axis) + np.mean(w_prime ** 2, axis=axis))
    return tke

def tke_ds(ds, axis='it'):
    #Var = ((ds-ds.mean(dim='it'))**2).mean(dim='it')
    var = ds.var(dim=axis)
    tke = 0.5 * (var.u + var.v + var.w)
    return tke.values
