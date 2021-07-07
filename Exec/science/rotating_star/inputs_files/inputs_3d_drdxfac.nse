
# INITIAL MODEL
maestro.model_file = "initial_models/15m_aprox19.6400"
maestro.perturb_model = false

maestro.drdxfac = 5
#maestro.ppm_type = 1

# PROBLEM SIZE
geometry.prob_lo     =  0.0e0   0.0e0  0.0e0
geometry.prob_hi     =   5.e8   5.e8   5.e8

# BOUNDARY CONDITIONS
# 0 = Interior   3 = Symmetry
# 1 = Inflow     4 = Slipwall
# 2 = Outflow    5 = NoSlipWall
maestro.lo_bc = 3 3 3
maestro.hi_bc = 2 2 2
geometry.is_periodic =  0 0 0

# VERBOSITY
maestro.v              = 1       # verbosity

# DEBUG FOR NAN
amrex.fpe_trap_invalid = 1       # floating point exception

# GRIDDING AND REFINEMENT
#amr.n_cell             = 160 160 160
#amr.n_cell             = 64 64 64
amr.n_cell = 128 128 128
amr.max_grid_size      = 64
amr.max_level          = 0       # maximum level number allowed
maestro.regrid_int     = 10      # how often to regrid
amr.ref_ratio          = 2 2 2 2 2 2 # refinement ratio
amr.blocking_factor    = 8       # block factor in grid generation
amr.refine_grid_layout = 0       # chop grids up into smaller grids if nprocs > ngrids
maestro.octant = true

# TIME STEPPING
maestro.max_step  = 200
maestro.stop_time = 300000.
maestro.cfl       = 0.7    # cfl number for hyperbolic system
                           # In this test problem, the velocity is
		           # time-dependent.  We could use 0.9 in
		           # the 3D test, but need to use 0.7 in 2D
		           # to satisfy CFL condition.

# ALGORITHMIC OPTIONS
maestro.spherical = 1
maestro.evolve_base_state = false
maestro.fix_base_state = false
maestro.do_initial_projection = true
maestro.init_divu_iter        = 3
maestro.init_iter             = 1

maestro.anelastic_cutoff_density = 1.e3
maestro.base_cutoff_density = 1.e1

maestro.do_sponge = 0
maestro.sponge_center_density = 3.e6
maestro.sponge_start_factor = 3.333e0
maestro.sponge_kappa = 10.e0

maestro.init_shrink = 0.001
maestro.use_soundspeed_firstdt = true
maestro.use_divu_firstdt = true

maestro.use_tfromp = true
maestro.use_exact_base_state = false

# PLOTFILES
maestro.plot_base_name  = plt    # root name of plot file
maestro.plot_int   = 10      # number of timesteps between plot files
maestro.plot_deltat = -1.
amr.small_plot_vars = rho Pi p0 tpert circ_velocity radial_velocity

maestro.plot_aux = true
maestro.plot_auxdot = true

# CHECKPOINT
maestro.check_base_name = chk
maestro.chk_int         = 10

# tolerances for the initial projection
maestro.eps_init_proj_cart = 1.e-12
maestro.eps_init_proj_sph  = 1.e-10
# tolerances for the divu iterations
maestro.eps_divu_cart      = 1.e-12
maestro.eps_divu_sph       = 1.e-10
maestro.divu_iter_factor   = 100.
maestro.divu_level_factor  = 10.
# tolerances for the MAC projection
maestro.eps_mac            = 1.e-10
maestro.eps_mac_max        = 1.e-8
maestro.mac_level_factor   = 10.
maestro.eps_mac_bottom     = 1.e-3
# tolerances for the nodal projection
maestro.eps_hg             = 1.e-10
maestro.eps_hg_max         = 1.e-10
maestro.hg_level_factor    = 10.
maestro.eps_hg_bottom      = 1.e-4

# BURNING 
maestro.do_burning = true
maestro.small_dens = 1.e0
maestro.small_temp = 1.e5
maestro.burning_cutoff_density_lo = 5.e5
maestro.burning_cutoff_density_hi = 5.e6
maestro.do_heating = true
maestro.plot_Hext = true
maestro.heating_cutoff_density_lo = 5.e6

#SDC 
maestro.sdc_couple_velocity = true

# ROTATION
# For the 18M star, I think the Keplerian rot freq is about 1e-5 (if assume its radius is about that of the Sun?)
maestro.rotational_frequency = 1.e-3
maestro.rotation_radius = 2.e8

&probin

  ! override the default values of the probin namelist values here
  velpert_amplitude = 1.e7
  velpert_radius = 2.e7
  velpert_scale = 1.e7
  velpert_steep = 1.e5
  tag_density_1 = 1.e7
  use_analytic_heating = T

  !extern

  ! Note that some of the parameters in this
  ! namelist are specific to the default EOS,
  ! network, and/or integrator used in the
  ! makefile. If you try a different set of
  ! microphysics routines be sure to check that
  ! the parameters in here are consistent, as
  ! Fortran does not like seeing unknown variables
  ! in the namelist.

  use_eos_coulomb = T

  retry_burn = T

  call_eos_in_rhs = T
  jacobian = 1
  
  atol_spec = 1.e-6
  atol_temp = 1.e-6
  atol_enuc = 1.e-6

  rtol_spec = 1.e-6
  rtol_temp = 1.e-6
  rtol_enuc = 1.e-6

  renormalize_abundances = T

  rho_nse = 3.e8
  T_nse = 3.e9
/
