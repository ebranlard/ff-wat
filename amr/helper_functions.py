import numpy as np
import matplotlib.pyplot as plt
import xarray
from welib.tools.curve_fitting import model_fit

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
            ax.set_xlabel('y [m]')
            ax.set_ylabel('z [m]')
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

    def K_model(y, p): #, U0, u, D):
        kDef = p[0]
        kGrad = p[1]
        K, KD, KG =  kWAT(du, gu, U0, D, kDef, kGrad)
        return K

    k_fit, ks_fit, fitter = model_fit(K_model, var.y, np.sqrt(var), bounds={'kd':(0, 5),'kg':(0, 5)}, p0={'kd':1.0,'kg':1.0}) #, U0=U0, u=u, D=D)
    return k_fit, ks_fit, fitter


def ds_fitK(ds, iP, ITime, D, U0, symmetric=True, removeBG='None', ds0=None):
    var = compute_var(ds, iP, ITime, D, removeBG=removeBG, ysym=symmetric, BGlim=1.9, ds0=ds0)
    u   = compute_u  (ds, iP, ITime, ysym=symmetric)
    k_fit, ks_fit, fitter =  fitK(var.y, u, var, U0, D)
    return var, u, k_fit, ks_fit, fitter



def ds_ysym(ds, y0=0):
    """ Make ds symmetric about y"""
    y = ds.y.values-y0
    iy0 = np.argmin(np.abs(y-0    ))
    IP = np.arange(iy0,len(y))
    IN = np.arange(iy0,-1,-1)
    assert( all(y[IP]==-y[IN]) )
    ds2 = ds.isel(y=IP)
    ds2['y'] = y[IP]
    if isinstance(ds, xarray.DataArray):
        ds2.loc[:] = (ds.values[IP]+ds.values[IN])/2  
    else:
        for k in list(ds.keys()):
            dims = ds[k].dims
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
    #import pdb; pdb.set_trace()
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


def common_itime(*dss):
    ITime = set(dss[0].it.values)
    for ds in dss:
        ITime = ITime.intersection(ds.it.values)
    ITime = list(ITime)
    ITime.sort()
    return ITime



def getSimParamsAMR(case):
    dt=1
    xyWT1         = (706.25 , 641.25)   # Actuator.T1.base_position = 713.28 641.25 0    # hub at 701.25, 641.25, 0, considering overhang of 12.0313
    xyWT2         = (2386.25, 641.25)   # Actuator.T2.base_position = 2393.28 641.25 0.  # hub at 2381.25, 641.25, 0, considering overhang of 12.0313
    D             = 240
    xPlanes = [0,1,2,3,4,5,6]
    if case == 'neutral':
        U0=8.5;
    elif case == 'stable': 
        U0=8.1;
    elif case == 'unstable': 
        U0=7.23;
    # fixed_dt =   0.025
    # output_frequency = 40  # every 1 s
    # # R             = 120.97
    # OverHang      = -12.097571763912535
    # ShftTilt      = -6.0 
    # Twr2Shft      = 4.349459414248071
    # TowerHt       = 144.386          
    #HubHeight  =  TowerHt + Twr2Shft + OverHang*np.sin(ShftTilt*np.pi/180)
    hubHeight=150
    return U0, dt, D, xyWT1, xyWT2, xPlanes


def compute_var(ds, iP, ITime, D, removeBG='None', ysym=False, BGlim=1.9, ds0=None):

    # First compute variance
    var = ds.isel(x=iP, it=ITime).var(dim='it').u
    if ds0 is not None:
        var0 = ds0.isel(x=iP, it=ITime).var(dim='it').u

    # Then do symmetry (order important!)
    if ysym:
        var = ds_ysym(var)
        if ds0 is not None:
            var0 = ds_ysym(var0)
    if removeBG=='Zero':
        varbg =0
    elif removeBG=='0WT':
        varbg = var0
    elif removeBG=='Outer':
        varbg = var.isel(y=np.abs(var.y)/D>BGlim).mean().values
    elif removeBG=='Min':
        varbg = np.min(var)
    else:
        raise NotImplementedError(removeBG)

    var -= varbg

    # Deal with negative values
    #var[var<0] = -var[var<0]
    var[var<0] = 0
    return var

def compute_u(ds, iP, ITime, ysym=True):
    u = ds.isel(x = iP, it = ITime).mean(dim = 'it').u
    if ysym:
        u = ds_ysym(u)
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
