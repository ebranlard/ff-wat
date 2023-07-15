#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#          SIMULATION CONTROL           #
#.......................................#
time.stop_time               =   26800.0     # Max (simulated) time to evolve
time.max_step                =   -1         # Max number of time steps
time.fixed_dt         =   0.025        # Use this constant dt if > 0
time.cfl              =   0.95         # CFL factor

incflo.verbose                           =   1          # incflo_level

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            INPUT AND OUTPUT           #
#.......................................#
time.plot_interval            =  24000       # Steps between plot files, every 600 s
time.checkpoint_interval      =  8000       # Steps between checkpoint files, every 200 s

ABL.bndry_file = "/projects/car/rthedin/amr_runs/02_precursor_shell/stable.W.8at150.20dTinv_0.25cooling_0.1z0_450zi_3.84x1.28x0.9km_res2.5m_coriolis5days/bndry_file.nc"
ABL.bndry_io_mode = 1
ABL.bndry_planes = xlo
ABL.bndry_output_start_time = 24999.0
ABL.bndry_output_format = native
ABL.bndry_var_names = velocity temperature tke

io.KE_int = -1
io.line_plot_int = 1
io.restart_file                          = /projects/car/rthedin/amr_runs/02_precursor_shell/stable.W.8at150.20dTinv_0.25cooling_0.1z0_450zi_3.84x1.28x0.9km_res2.5m_coriolis5days/chk121692
io.outputs = actuator_src_term
io.derived_outputs = q_criterion

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            GEOMETRY & BCs             #
#.......................................#
geometry.prob_lo        =   0. 0.    0.  # Lo corner coordinates
geometry.prob_hi        =   3840. 1280.  640.  # Hi corner coordinates
amr.n_cell              =  1536 512  256    # Grid cells at coarsest AMRlevel
geometry.is_periodic    =   0   1   0   # Periodicity x y z (0/1)
incflo.delp             =   0.  0.  0.  # Prescribed (cyclic) pressure gradient


# Boundary conditions
xlo.type                                 = mass_inflow         
xlo.density                              = 1.225               
xlo.temperature                          = 0
xlo.tke                                  = 0.0
xhi.type                                 = pressure_outflow    

zlo.type                                 = wall_model
#zlo.temperature_type                     = fixed_gradient
#zlo.temperature                          = 0                       # !!!!!!!!!!!!!!!!!!!!!! does this matter for unstable/stable?!
#zlo.tke_type                             = zero_gradient

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
ICNS.source_terms = BoussinesqBuoyancy CoriolisForcing BodyForce ActuatorForcing ABLMeanBoussinesq
#---------------- Additions by calc_inflow_stats.py -----------------#
ABL.wall_shear_stress_type = "local"
ABL.inflow_outflow_mode    = true
ABL.wf_velocity            = 1.324928 0.091594
ABL.wf_vmag                = 1.340386632038509
ABL.wf_theta               = 298.4857241706984
BodyForce.magnitude        = 0.00011625217263954034 0.00015944554243391484 0.0
BoussinesqBuoyancy.read_temperature_profile = true
BoussinesqBuoyancy.tprofile_filename        = avg_theta.dat
#--------------------------------------------------------------------#
incflo.velocity = 8.0 0.0 0.0

ABLForcing.abl_forcing_height= 150.

CoriolisForcing.latitude = 40.761  # Mayflower farm
CoriolisForcing.north_vector = 0.0 1.0 0.0
CoriolisForcing.east_vector = 1.0 0.0 0.0

incflo.physics = ABL Actuator

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
TKE.interpolation                        = PiecewiseConstant          


#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#          POST-Processing              #
#.......................................#
#io.output_hdf5_plotfile                  = true
#io.hdf5_compression                      = "ZFP_ACCURACY@0.001"

incflo.post_processing                   = averaging samplingxy samplingxz planesT1 planesT2

# --- Sampling parameters ---
samplingxy.output_format                   = netcdf
samplingxy.output_frequency                = 2400  # every 60 s
samplingxy.fields                          = velocity #temperature
#---- sample defs ----
samplingxy.labels               = xyplane
samplingxy.xyplane.type         = PlaneSampler        
samplingxy.xyplane.num_points   = 1536 512
samplingxy.xyplane.origin       = 1.25  1.25   26.25
samplingxy.xyplane.axis1        = 3837.5    0.0   0.0      
samplingxy.xyplane.axis2        =   0.0  1277.5   0.0      
samplingxy.xyplane.normal       = 0.0 0.0 1.0         
samplingxy.xyplane.offsets      = 0 125 250


# --- Sampling parameters ---
samplingxz.output_format                   = netcdf
samplingxz.output_frequency                = 2400  # every 60 s 
samplingxz.fields                          = velocity #temperature
#---- sample defs ----
samplingxz.labels               = xzplane
samplingxz.xzplane.type         = PlaneSampler        
samplingxz.xzplane.num_points   = 1536 241
samplingxz.xzplane.origin       = 1.25  641.25   1.25
samplingxz.xzplane.axis1        = 3837.5 0.0   0.0     
samplingxz.xzplane.axis2        =   0.0  0.0   600.0
samplingxz.xzplane.normal       = 0.0 1.0 0.0         
samplingxz.xzplane.offsets      = 0


## --- Sampling parameters ---
#samplingyz.output_format                   = netcdf
#samplingyz.output_frequency                = 20  # every 0.5 s
#samplingyz.fields                          = velocity 
##---- sample defs ----
#samplingyz.labels               = yzplane
#samplingyz.yzplane.type         = PlaneSampler        
#samplingyz.yzplane.num_points   = 512  241
#samplingyz.yzplane.origin       = 1001.25  1.25   1.25
#samplingyz.yzplane.axis1        =   0.0 1277.5   0.0      
#samplingyz.yzplane.axis2        =   0.0  0.0     600.0
#samplingyz.yzplane.normal       = 1.0 0.0 0.0         
#samplingyz.yzplane.offsets      = 0 250 500 750 1000 1500 2000


# --- Sampling parameters ---
planesT1.output_format    = netcdf
planesT1.output_frequency = 40  # every 1 s
planesT1.fields           = velocity temperature
#---- sample defs ----
planesT1.labels           = pT1
planesT1.pT1.type         = PlaneSampler
planesT1.pT1.num_points   = 201 101
planesT1.pT1.origin       = 706.25  141.25   11.25
planesT1.pT1.axis1        = 0.0   1000   0.0
planesT1.pT1.axis2        = 0.0    0.0   500
planesT1.pT1.normal       = 1.0    0.0   0.0
planesT1.pT1.offsets      = 0 240 480 720 960 1200 1440

# --- Sampling parameters ---
planesT2.output_format    = netcdf
planesT2.output_frequency = 40  # every 1 s
planesT2.fields           = velocity temperature
#---- sample defs ----
planesT2.labels           = pT2
planesT2.pT2.type         = PlaneSampler
planesT2.pT2.num_points   = 201 101
planesT2.pT2.origin       = 2386.25   141.25   11.25
planesT2.pT2.axis1        = 0.0   1000   0.0
planesT2.pT2.axis2        = 0.0    0.0   500
planesT2.pT2.normal       = 1.0    0.0   0.0
planesT2.pT2.offsets      = 0 240 480 720 960 1200 1440


#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#              AVERAGING                #
#.......................................#
averaging.type                           = TimeAveraging
averaging.labels                         = means stress

averaging.averaging_start_time           = 25000
averaging.averaging_window               = 600.0

averaging.means.fields                   = velocity
averaging.means.averaging_type           = ReAveraging

averaging.stress.fields                  = velocity
averaging.stress.averaging_type          = ReynoldsStress

#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#            MESH REFINEMENT            #
#.......................................#
amr.max_level           =  1           # Max AMR level in hierarchy

#---- tagging defs ----
tagging.labels                           = T1_level_0_zone T2_level_0_zone

tagging.T1_level_0_zone.type             = GeometryRefinement  
tagging.T1_level_0_zone.shapes           = T1_level_0_zone     
tagging.T1_level_0_zone.level            = 0                   
tagging.T1_level_0_zone.T1_level_0_zone.type = box                 
tagging.T1_level_0_zone.T1_level_0_zone.origin = 580 490.0 10.0
tagging.T1_level_0_zone.T1_level_0_zone.xaxis = 240 0.0 0.0
tagging.T1_level_0_zone.T1_level_0_zone.yaxis = 0.0 300.0 0.0
tagging.T1_level_0_zone.T1_level_0_zone.zaxis = 0.0 0.0 300.0       
tagging.T2_level_0_zone.type             = GeometryRefinement  
tagging.T2_level_0_zone.shapes           = T2_level_0_zone     
tagging.T2_level_0_zone.level            = 0                   
tagging.T2_level_0_zone.T2_level_0_zone.type = box                 
tagging.T2_level_0_zone.T2_level_0_zone.origin = 2260.0 490.0 10.0  
tagging.T2_level_0_zone.T2_level_0_zone.xaxis = 240.0 0.0 0.0
tagging.T2_level_0_zone.T2_level_0_zone.yaxis = 0.0 300.0 0.0
tagging.T2_level_0_zone.T2_level_0_zone.zaxis = 0.0 0.0 300.0       


#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
#               TURBINES                #
#.......................................#


Actuator.labels = T1 T2
Actuator.type = TurbineFastLine

Actuator.TurbineFastLine.rotor_diameter = 240.0
Actuator.TurbineFastLine.hub_height = 150.
Actuator.TurbineFastLine.num_points_blade = 48
Actuator.TurbineFastLine.num_points_tower = 12
Actuator.TurbineFastLine.epsilon =  2.5 2.5 2.5
Actuator.TurbineFastLine.epsilon_tower = 2.5 2.5 2.5
Actuator.TurbineFastLine.openfast_start_time = 0.0
Actuator.TurbineFastLine.openfast_stop_time = 1801.0
Actuator.TurbineFastLine.nacelle_drag_coeff = 0.5
Actuator.TurbineFastLine.nacelle_area = 49.5
Actuator.TurbineFastLine.output_frequency = 10
Actuator.TurbineFastLine.density = 1.225

# hub at 701.25, 641.25, 0, considering overhang of 12.0313
Actuator.T1.base_position = 713.28 641.25 0
Actuator.T1.openfast_input_file = "turbine_iea15mw/IEA-15-240-RWT-Monopile.T1.fst"

# hub at 2381.25, 641.25, 0, considering overhang of 12.0313
Actuator.T2.base_position = 2393.28 641.25 0.
Actuator.T2.openfast_input_file = "turbine_iea15mw/IEA-15-240-RWT-Monopile.T2.fst"


