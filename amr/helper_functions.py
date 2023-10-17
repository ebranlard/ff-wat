
def getSimParamsAMR(case):
    dt=1
    xyWT1         = (706.25 , 641.25)   # Actuator.T1.base_position = 713.28 641.25 0    # hub at 701.25, 641.25, 0, considering overhang of 12.0313
    xyWT2         = (2386.25, 641.25)   # Actuator.T2.base_position = 2393.28 641.25 0.  # hub at 2381.25, 641.25, 0, considering overhang of 12.0313
    D             = 240
    xPlanes = [0,1,2,3,4,5,6]
    if case == 'neutral':
        U0=8.5;
    if case == 'stable': 
        U0=8.1;

    return U0, dt, D, xyWT1, xyWT2, xPlanes

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
