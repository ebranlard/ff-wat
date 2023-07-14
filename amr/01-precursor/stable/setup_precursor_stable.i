#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#          SIMULATION CONTROL           #
#.......................................#
time.stop_time               =   60000.0     # Max (simulated) time to evolve
time.max_step                =   -1         # Max number of time steps
time.fixed_dt         =   -1.0        # Use this constant dt if > 0 # changed this. feb 13
time.cfl              =   0.9         # CFL factor

incflo.verbose                           =   1          # incflo_level

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            INPUT AND OUTPUT           #
#.......................................#
time.plot_interval            =  20000       # Steps between plot files
time.checkpoint_interval      =  10000       # Steps between checkpoint files

ABL.bndry_file = "bndry_file.nc"
ABL.bndry_io_mode = 0
ABL.bndry_planes = xlo
ABL.bndry_output_start_time = 129999.0
ABL.bndry_output_format = native
ABL.bndry_var_names = velocity temperature tke

#io.restart_file                          = unused

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            GEOMETRY & BCs             #
#.......................................#
geometry.prob_lo        =   0. 0.    0.  # Lo corner coordinates
geometry.prob_hi        =   3840. 1280.  640.  # Hi corner coordinates
amr.n_cell              =  1536 512  256    # Grid cells at coarsest AMRlevel
amr.max_level           = 1           # Max AMR level in hierarchy 
geometry.is_periodic    =   1   1   0   # Periodicity x y z (0/1)
incflo.delp             =   0.  0.  0.  # Prescribed (cyclic) pressure gradient


# Boundary conditions
#xlo.type                                 = mass_inflow         
#xlo.density                              = 1.225               
#xlo.temperature                          = 290.0               
#xlo.tke                                  = 0.0
#xhi.type                                 = pressure_outflow    

#ylo.type                                 = mass_inflow         
#ylo.density                              = 1.225               
#ylo.temperature                          = 290.0               
#ylo.tke                                  = 0.0
#yhi.type                                 = pressure_outflow     

zlo.type                                 = wall_model
zhi.type                                 = slip_wall
zhi.temperature_type                     = fixed_gradient
zhi.temperature                          = 0.003

# MLMG settings
nodal_proj.mg_rtol = 1.0e-6
nodal_proj.mg_atol = 1.0e-12
mac_proj.mg_rtol = 1.0e-6
mac_proj.mg_atol = 1.0e-12
diffusion.mg_rtol = 1.0e-6
diffusion.mg_atol = 1.0e-12
temperature_diffusion.mg_rtol = 1.0e-10
temperature_diffusion.mg_atol = 1.0e-13

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#        Mesh refinement                #
#.......................................#


#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               PHYSICS                 #
#.......................................#
ICNS.source_terms = CoriolisForcing ABLForcing BoussinesqBuoyancy ABLMeanBoussinesq
incflo.velocity = 8.0 0.0 0.0

ABLForcing.abl_forcing_height= 150.

CoriolisForcing.latitude = 40.761  # Mayflower farm
CoriolisForcing.north_vector = 0.0 1.0 0.0
CoriolisForcing.east_vector = 1.0 0.0 0.0
CoriolisForcing.rotational_time_period = 432000

incflo.physics = ABL # Actuator

BoussinesqBuoyancy.reference_temperature = 300.0
ABL.reference_temperature = 300.0
ABL.temperature_heights                  = 0.0 400.0 500.0 2000.0
ABL.temperature_values                   = 300.0 300.0 310.0 324.5
ABL.perturb_temperature = false
ABL.cutoff_height = 50.0
ABL.perturb_velocity = true
ABL.perturb_ref_height = 200.0
ABL.Uperiods = 15.0
ABL.Vperiods = 10.0
ABL.deltaU = 0.25
ABL.deltaV = 0.25
ABL.kappa = .41
ABL.mo_gamma_m = 4.8
ABL.mo_gamma_h = 7.8
ABL.surface_roughness_z0 = 0.10
ABL.surface_temp_rate = -0.25
ABL.surface_temp_init = 300

incflo.use_godunov                       = 1
incflo.godunov_type                      = weno_z
incflo.gravity          =   0.  0. -9.81  # Gravitational force (3D)
incflo.density          = 1.225          # Reference density
transport.viscosity = 1.0e-5
transport.laminar_prandtl = 0.7
transport.turbulent_prandtl = 0.3333
turbulence.model = OneEqKsgsM84
#Smagorinsky_coeffs.Cs = 0.08
TKE.source_terms                         = KsgsM84Src
#TKE.interpolation                        = PiecewiseConstant          


#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#          POST-Processing              #
#.......................................#
#io.output_hdf5_plotfile                  = true
#io.hdf5_compression                      = "ZFP_ACCURACY@0.001"

incflo.post_processing                   = averaging



#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              AVERAGING                #
#.......................................#
averaging.type                           = TimeAveraging
averaging.labels                         = means stress

averaging.averaging_start_time           = 150000
averaging.averaging_window               = 600.0

averaging.means.fields                   = velocity
averaging.means.averaging_type           = ReAveraging

averaging.stress.fields                  = velocity
averaging.stress.averaging_type          = ReynoldsStress

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            MESH REFINEMENT            #
#.......................................#


#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               TURBINES                #
#.......................................#

