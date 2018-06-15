
!#########################################################
!#####          COMMON USER NAMELIST INPUTS          #####
!#########################################################
 GRID_FORMAT                        = -1
 GRIDFILE                           = ""
 PLOT3D_CUTFILE                     = ""
 PLOT3D_BCSFILE                     = ""
 GRID_SCALING_FACTOR                =  1.0
 PLOT3D_AGGLOMERATE_ORDER           = -1
 SOLUTION_ORDER                     =  3
 RUNGE_KUTTA_SCHEME                 =  3
 NUM_RK_STAGES                      =  3
 NUM_TIMESTEPS                      =  25000
 FINAL_TIME                         =  10.0
 CONSTANT_DT                        =  1.0E-6
 TIMESTEP_TYPE                      =  1
 MINIMUM_TIMESTEPS                  =  100
 CFL                                =  0.1
 CFL_BEG                            =  0.1
 CFL_END                            =  0.1
 CFL_CYCLES                         =  1
 ALL_ORDERS                         = -3
 ITESTCASE                          =  9
 LOAD_RESTART_FILE                  = .FALSE.
 RESTART_INTERVAL                   =  0
 RESTART_FILE                       = "out.rst"
 ITER_OUT_INTERVAL                  =  1
 RESULTS_INTERVAL                   =  0
 OUTPUT_DIR                         = "."
 OUTPUT_INTERVAL                    = -5
 GOVERNING_EQUATIONS                =  2
 MACHREF                            =  0.4
 REYREF                             = -1.0
 GAM                                =  1.4
 RHOREF                             = -1.0
 TREF                               =  300.0
 PREF                               =  100000.0
 RGASREF                            =  287.15
 ALPHA_AOAREF                       =  0.0
 BETA_AOAREF                        =  0.0
 PR                                 =  0.72
 SUTH_MUREF                         =  1.716E-5
 SUTH_TREF                          =  273.0
 SUTH_SREF                          =  110.4
 OUTPUT_TIME_AVERAGING              = .FALSE.
 TIME_AVE_FILE                      = "out.ave"
 TIME_SCALING_FACTOR                =  1.0
!#########################################################
!#####         ADVANCED USER NAMELIST INPUTS         #####
!#########################################################
 PLOT3D_I_STRIDE                    = -1
 PLOT3D_J_STRIDE                    = -1
 PLOT3D_K_STRIDE                    = -1
 PLOT3D_ALL_STRIDES                 = -1
 CORRECTION_FUNCTION                =  1
 LOC_SOLUTION_PTS                   =  1
 LOC_FLUX_PTS                       =  1
 INVS_FLUX_METHOD                   =  1
 VISC_FLUX_METHOD                   =  2
 PROJECTION_ORDER                   =  3
 QUADRATURE_ORDER                   =  3
 LOC_QUADRATURE_PTS                 =  1
 ERROR_ORDER                        = -3
 LOC_ERROR_PTS                      =  1
 ALL_POINTS                         = -1
 OUTPUT_ORDER                       = -3
 FILTER_OPTION                      =  0
 FILTER_ETAC                        =  0.0
 FILTER_ORDER                       =  16
 FILTER_ALPHA                       =  36.0
 FILTER_EPS                         =  0.0
 LIMITER_FILTER_OPTION              =  0
 LIMITER_FILTER_ETAC                =  0.0
 LIMITER_FILTER_ORDER               =  16
 LIMITER_FILTER_ALPHA               =  36.0
 LIMITER_FILTER_EPS                 =  0.0
 LIMITER_OPTION                     =  0
 FILTERING_INTERVAL                 =  0
 LIMITING_INTERVAL                  =  0
 WALLS_ARE_EXACT                    = .TRUE.
 SUB_INFLOW_METHOD                  =  1
 COMPLETELY_DISABLE_CGNS            = .FALSE.
 CONVERGENCE_ORDER_ABS_RES          = -16.0
 CONVERGENCE_ORDER_MAX_RES          = -16.0
 TIME_AVE_VEL_IS_AXISYMM            = .FALSE.
!#########################################################
!#####       BEST IF LEFT ALONE NAMELIST INPUTS      #####
!#########################################################
 LOC_INIT_PTS                       =  0
 VISC_VARIABLES_TO_AVE              =  1
 METIS_OPTION                       =  3
 RESTART_OUTPUT_USES_HOST_ROOTS     = .TRUE.
 LOC_OUTPUT_PTS                     =  0
 LUSTRE_STRIPE_COUNT                =  64
 LUSTRE_STRIPE_SIZE                 = "4m"
!#########################################################
!#####          BC CONTROL NAMELIST INPUTS           #####
!#########################################################
 BC_INPUT(1)%NAME         = ""
 BC_INPUT(1)%SET_BC_VALUE = -9999
 BC_INPUT(1)%T_STATIC     = -1.0
 BC_INPUT(1)%P_STATIC     = -1.0
 BC_INPUT(1)%RHO_STATIC   = -1.0
 BC_INPUT(1)%MACH         = -1.0
 BC_INPUT(1)%VX           =  0.0
 BC_INPUT(1)%VY           =  0.0
 BC_INPUT(1)%VZ           =  0.0
 BC_INPUT(1)%T_TOTAL      = -1.0
 BC_INPUT(1)%P_TOTAL      = -1.0
 BC_INPUT(1)%RHO_TOTAL    = -1.0
 BC_INPUT(1)%ALPHA_AOA    =  0.0
 BC_INPUT(1)%BETA_AOA     =  0.0
 BC_INPUT(1)%WALL_TEMP    = -1.0
 BC_INPUT(1)%RELAX%VALUE  =  1.0
 BC_INPUT(1)%RELAX%TIME   =  0.0
 BC_INPUT(1)%RELAX%ITER   =  0
!#########################################################
!#####      BUILT IN UTILITIES NAMELIST INPUTS       #####
!#########################################################
 CONVERT_RESTART_TO_CGNS            = .FALSE.
 TIME_AVERAGE_RESTART_FILES         = .FALSE.
 RESTART_LIST_FILE                  = ""
 OUTPUT_PLOT3D_TO_TECPLOT           = .FALSE.
!#########################################################
!#####         SPECIAL CASE NAMELIST INPUTS          #####
!#########################################################
 VORTEX_BIGR                        =  1.0
 VORTEX_BIGGAM                      =  5.0
 VORTEXGRID_LREF                    =  10.0
 VORTEXGRID_DOF                     =  48
 VORTEXGRID_RANDOM_PERTURB          = .FALSE.
 VORTEXGRID_PERTURB_PERCENT         =  30
 USE_BC_CONFLICTS                   = .TRUE.
!#########################################################
!#####           DEVELOPER NAMELIST INPUTS           #####
!#########################################################
 LOC_TRIANGLE_PTS                   =  2
 USE_OLD_RESTART_FORMAT             = .FALSE.
 CONTINUOUS_OUTPUT                  = .FALSE.
 CGNS_USE_QUEUE                     = .FALSE.
 INTERPOLATE_BEFORE_OUTPUT_VARIABLE = .TRUE.
 PROFILE_IO_ALL                     = .FALSE.
 PROFILE_IO_CGNS                    = .FALSE.
 PROFILE_IO_RESTART                 = .FALSE.
 USE_UNROLLED_DOT_PRODUCTS          = .TRUE.
 DUMP_MEMORY_USAGE                  = .FALSE.
 DUMP_FLUXES                        = .FALSE.
 CHECK_FLUX_POINT_COORDINATES       = .FALSE.
 POSTPROC_DEBUG_GRID                = .FALSE.
 OUTPUT_BFACE_ARRAY                 = .FALSE.
!#########################################################
!###   EXPERIMENTAL / USE-AT-OWN-RISK NAMELIST INPUTS  ###
!#########################################################
 MMS_OPT                            =  0
 MMS_INIT                           =  0.1
 MMS_OUTPUT_SOLUTION                = .FALSE.
 OUTPUT_BND_PROFILES                = .FALSE.
!#########################################################
!#####          NON-WORKING NAMELIST INPUTS          #####
!#########################################################
 FLUX_DIVERGENCE                    =  1
 TURBULENCE_MODEL                   =  0
 PTOTREF                            = -1.0
 TTOTREF                            = -1.0
 PTOT2P_RATIO                       = -1.0
 TINREF                             =  0.0
 TVRREF                             =  1.0E-6
 PRT                                =  0.72
 CHANNEL_BODY_FORCE                 =  0.0
 CHANNEL_INIT_FILE                  = ""
!#########################################################





























!###############################################################################
!###############################################################################
!###   DEFAULT NAMELIST WITH VALUES
!###############################################################################
!###############################################################################
 &INPUT

 GRID_FORMAT                        = -1
 GRIDFILE                           = ""
 PLOT3D_CUTFILE                     = ""
 PLOT3D_BCSFILE                     = ""
 GRID_SCALING_FACTOR                =  1.0

 PLOT3D_I_STRIDE                    = -1
 PLOT3D_J_STRIDE                    = -1
 PLOT3D_K_STRIDE                    = -1
 PLOT3D_ALL_STRIDES                 = -1
 PLOT3D_AGGLOMERATE_ORDER           = -1

 SOLUTION_ORDER                     =  3
 FLUX_DIVERGENCE                    =  1
 CORRECTION_FUNCTION                =  1

 LOC_SOLUTION_PTS                   =  1
 LOC_FLUX_PTS                       =  1
 LOC_TRIANGLE_PTS                   =  2

 LOC_INIT_PTS                       =  0

 RUNGE_KUTTA_SCHEME                 =  3
 NUM_RK_STAGES                      =  3

 INVS_FLUX_METHOD                   =  1
 VISC_FLUX_METHOD                   =  2
 VISC_VARIABLES_TO_AVE              =  1

 NUM_TIMESTEPS                      =  25000
 FINAL_TIME                         =  10.0
 CONSTANT_DT                        =  1.0E-6
 TIMESTEP_TYPE                      =  1
 MINIMUM_TIMESTEPS                  =  100

 CFL                                =  0.1
 CFL_BEG                            =  0.1
 CFL_END                            =  0.1
 CFL_CYCLES                         =  1

 PROJECTION_ORDER                   =  3
 QUADRATURE_ORDER                   =  3
 LOC_QUADRATURE_PTS                 =  1

 ERROR_ORDER                        = -3
 LOC_ERROR_PTS                      =  1

 ALL_ORDERS                         = -3
 ALL_POINTS                         = -1

 METIS_OPTION                       =  3

 MMS_OPT                            =  0
 MMS_INIT                           =  0.1
 MMS_OUTPUT_SOLUTION                = .FALSE.

 ITESTCASE                          =  9

 LOAD_RESTART_FILE                  = .FALSE.
 RESTART_INTERVAL                   =  0
 USE_OLD_RESTART_FORMAT             = .FALSE.
 RESTART_FILE                       = "out.rst"
 RESTART_OUTPUT_USES_HOST_ROOTS     = .TRUE.

 ITER_OUT_INTERVAL                  =  1
 RESULTS_INTERVAL                   =  0

 OUTPUT_DIR                         = "."
 OUTPUT_INTERVAL                    = -5
 CONTINUOUS_OUTPUT                  = .FALSE.
 OUTPUT_ORDER                       = -3
 LOC_OUTPUT_PTS                     =  0

 GOVERNING_EQUATIONS                =  2
 TURBULENCE_MODEL                   =  0

 MACHREF                            =  0.4
 REYREF                             = -1.0
 GAM                                =  1.4

 PTOTREF                            = -1.0
 TTOTREF                            = -1.0
 PTOT2P_RATIO                       = -1.0

 RHOREF                             = -1.0
 TREF                               =  300.0
 PREF                               =  100000.0
 RGASREF                            =  287.15

 TINREF                             =  0.0
 TVRREF                             =  1.0E-6

 ALPHA_AOAREF                       =  0.0
 BETA_AOAREF                        =  0.0

 PR                                 =  0.72
 PRT                                =  0.72

 SUTH_MUREF                         =  1.716E-5
 SUTH_TREF                          =  273.0
 SUTH_SREF                          =  110.4

 FILTER_OPTION                      =  0
 FILTER_ETAC                        =  0.0
 FILTER_ORDER                       =  16
 FILTER_ALPHA                       =  36.0
 FILTER_EPS                         =  0.0

 LIMITER_FILTER_OPTION              =  0
 LIMITER_FILTER_ETAC                =  0.0
 LIMITER_FILTER_ORDER               =  16
 LIMITER_FILTER_ALPHA               =  36.0
 LIMITER_FILTER_EPS                 =  0.0

 LIMITER_OPTION                     =  0

 FILTERING_INTERVAL                 =  0
 LIMITING_INTERVAL                  =  0

 VORTEX_BIGR                        =  1.0
 VORTEX_BIGGAM                      =  5.0
 VORTEXGRID_LREF                    =  10.0

 VORTEXGRID_DOF                     =  48
 VORTEXGRID_RANDOM_PERTURB          = .FALSE.
 VORTEXGRID_PERTURB_PERCENT         =  30

 CHANNEL_BODY_FORCE                 =  0.0
 CHANNEL_INIT_FILE                  = ""

 WALLS_ARE_EXACT                    = .TRUE.
 SUB_INFLOW_METHOD                  =  1
 USE_BC_CONFLICTS                   = .TRUE.

 CGNS_USE_QUEUE                     = .FALSE.
 INTERPOLATE_BEFORE_OUTPUT_VARIABLE = .TRUE.
 CONVERT_RESTART_TO_CGNS            = .FALSE.
 COMPLETELY_DISABLE_CGNS            = .FALSE.

 CONVERGENCE_ORDER_ABS_RES          = -16.0
 CONVERGENCE_ORDER_MAX_RES          = -16.0

 PROFILE_IO_ALL                     = .FALSE.
 PROFILE_IO_CGNS                    = .FALSE.
 PROFILE_IO_RESTART                 = .FALSE.

 OUTPUT_BND_PROFILES                = .FALSE.

 USE_UNROLLED_DOT_PRODUCTS          = .TRUE.

 LUSTRE_STRIPE_COUNT                =  64
 LUSTRE_STRIPE_SIZE                 = "4m"

 OUTPUT_TIME_AVERAGING              = .FALSE.
 TIME_AVE_FILE                      = "out.ave"
 TIME_SCALING_FACTOR                =  1.0
 TIME_AVE_VEL_IS_AXISYMM            = .FALSE.

 TIME_AVERAGE_RESTART_FILES         = .FALSE.
 RESTART_LIST_FILE                  = ""

 DUMP_MEMORY_USAGE                  = .FALSE.

 DUMP_FLUXES                        = .FALSE.
 CHECK_FLUX_POINT_COORDINATES       = .FALSE.

 POSTPROC_DEBUG_GRID                = .FALSE.

 OUTPUT_PLOT3D_TO_TECPLOT           = .FALSE.

 OUTPUT_BFACE_ARRAY                 = .FALSE.

 BC_INPUT(1)%NAME         = ""
 BC_INPUT(1)%SET_BC_VALUE = -9999
 BC_INPUT(1)%T_STATIC     = -1.0
 BC_INPUT(1)%P_STATIC     = -1.0
 BC_INPUT(1)%RHO_STATIC   = -1.0
 BC_INPUT(1)%MACH         = -1.0
 BC_INPUT(1)%VX           =  0.0
 BC_INPUT(1)%VY           =  0.0
 BC_INPUT(1)%VZ           =  0.0
 BC_INPUT(1)%T_TOTAL      = -1.0
 BC_INPUT(1)%P_TOTAL      = -1.0
 BC_INPUT(1)%RHO_TOTAL    = -1.0
 BC_INPUT(1)%ALPHA_AOA    =  0.0
 BC_INPUT(1)%BETA_AOA     =  0.0
 BC_INPUT(1)%WALL_TEMP    = -1.0
 BC_INPUT(1)%RELAX%VALUE  =  1.0
 BC_INPUT(1)%RELAX%TIME   =  0.0
 BC_INPUT(1)%RELAX%ITER   =  0

 /
!###############################################################################
!###############################################################################
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### GRID INFORMATION
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! GRID_FORMAT : the format of the input grid file
  !      -1 = Internally Generated Grid
  !       1 = Gambit Neutral
  !       2 = Gmsh Format
  !       3 = Plot3D
  !       4 = CGNS
  !
  integer :: grid_format = -1
  !
  ! GRIDFILE : the location of the input grid file relative to the
  !            current directory
  !
  character(len=150) :: gridfile = ""
  !
  ! PLOT3D_CUTFILE : the location of the file that specifies the
  !                  cut conditions for a multiblock Plot3D file
  !
  character(len=150) :: plot3d_cutfile = ""
  !
  ! PLOT3D_BCSFILE : the location of the file that specifies the
  !                  boundary conditions for block boundaries
  !
  character(len=150) :: plot3d_bcsfile = ""
  !
  ! GRID_SCALING_FACTOR : Multiplication factor (greater than 0) for the grid
  !                       coordinates to scale the grid to the desired units
  !                       (for example, to convert a grid from inches to meters
  !                       use 0.0254).
  !
  real :: grid_scaling_factor = 1.0
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### PLOT3D COARSENING OPTIONS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! PLOT3D_ALL_STRIDES : applies this coarsening stride to all directions
  !
  ! PLOT3D_I_STRIDE : applies this coarsening stride to only the I-direction
  ! PLOT3D_J_STRIDE : applies this coarsening stride to only the J-direction
  ! PLOT3D_K_STRIDE : applies this coarsening stride to only the K-direction
  !
  ! THESE STRIDES ARE APPLIED IN THE FOLLOWING ORDER
  !
  ! strides(:) = max(1,PLOT3D_ALL_STRIDES)
  !
  ! if (PLOT3D_I_STRIDE > 0) strides(1) = PLOT3D_I_STRIDE
  ! if (PLOT3D_J_STRIDE > 0) strides(2) = PLOT3D_J_STRIDE
  ! if (PLOT3D_K_STRIDE > 0) strides(3) = PLOT3D_K_STRIDE
  !
  integer :: plot3d_i_stride    = -1
  integer :: plot3d_j_stride    = -1
  integer :: plot3d_k_stride    = -1
  integer :: plot3d_all_strides = -1
  !
  ! PLOT3D_AGGLOMERATE_ORDER : agglomerates linear Plot3D cells into high-order
  !                            grid cells of this order
  !
  integer :: plot3d_agglomerate_order = -1
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### BASE FR OPTIONS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! SOLUTION_ORDER : specifies the order of accuracy of the solution
  !
  integer :: solution_order = 3
  !
  ! FLUX_DIVERGENCE : (NOT USED YET) the method used for computing the
  !                   divergence of the fluxes
  !      1 = using Lagrange polynomials
  !      2 = using chain rule
  !      3 = using chain rule with Jacobians
  !
  integer :: flux_divergence = 1
  !
  ! CORRECTION_FUNCTION : specifies the correction function to use
  !      1 = gDG    correction function (recovers DG method)
  !      2 = g2     correction function
  !      3 = gGauss correction function
  !
  integer :: correction_function = 1
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### LOCATION OF SOLUTION, FLUX, AND TRIANGULAR POINTS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! LOC_SOLUTION_PTS : specifies the location of the solution points
  !              used for tensor products
  !      1 = Legendre-Gauss         nodes
  !      2 = Legendre-Gauss-Lobatto nodes
  !
  integer :: loc_solution_pts = 1
  !
  ! LOC_FLUX_PTS : specifies the location of the flux points on cell faces
  !      1 = Legendre-Gauss         nodes
  !      2 = Legendre-Gauss-Lobatto nodes
  !
  integer :: loc_flux_pts = 1
  !
  ! LOC_TRIANGLE_PTS : specifies the location of the solution points
  !              used for triangles
  !      1 = Hesthaven-Warburton alpha-optimized nodes
  !      2 = Barycentric Lobatto nodes
  !
  integer :: loc_triangle_pts = 2
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### LOCATION OF INITIALIZATION POINTS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! LOC_INIT_PTS : nodal points used for initializing the flow, after which
  !                the initial solution is projected to the solution points
  !              1 = Legendre-Gauss         nodes
  !              2 = Legendre-Gauss-Lobatto nodes
  !  anything else = LOC_SOLUTION_PTS
  !
  integer :: loc_init_pts = 0
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### RUNGE-KUTTA OPTIONS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! RUNGE_KUTTA_SCHEME : specifies the Runge-Kutta method
  !      1 = classic n-stage RK method
  !      2 = 2-stage/2nd-order SSP-RK method
  !      3 = 3-stage/3rd-order SSP-RK method
  !      4 = 5-stage/4th-order SSP-RK method
  !      5 = 5-stage/4th-order Carpenter-Kennedy low storage RK method
  !
  integer :: Runge_Kutta_Scheme = 3
  !
  ! NUM_RK_STAGES : gives the number of Runge-Kutta stages to use if using the
  !              classic RK method (RUNGE_KUTTA_SCHEME = 1)
  !
  integer :: num_rk_stages = 3
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### COMMON FLUX OPTIONS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! INVS_FLUX_METHOD : specifies the approximate Riemann solver used to
  !                    compute the common inviscid fluxes at the interfaces
  !       1 = Roe flux with entropy fix
  !       2 = HLLC flux
  !       3 = LDFSS flux
  !       4 = Rotated-Roe-HLL flux
  !
  integer :: invs_flux_method = 1
  !
  ! VISC_FLUX_METHOD : specifies the method thats used for computing the
  !                    viscous fluxes at the interfaces
  !       1 = Bassi-Rebay 1  (I think this has stability problems)
  !       2 = Bassi-Rebay 2  (This seems to be the best in terms of stability
  !                           and largest CFL allowable)
  !       3 = Local DG (LDG) (Stable but seems to require a smaller CFL)
  !
  integer :: visc_flux_method = 2
  !
  ! VISC_VARIABLES_TO_AVE : determines which variables are averaged to the
  !                         interfaces when computing the common viscous fluxes
  !        1 = average the interface fluxes:
  !              Compute the common viscous flux by averaging the
  !              left and right viscous interface fluxes which are
  !              computed using the left and right interface
  !              solutions and gradients.
  !        2 = average the interface solutions and gradients:
  !              Compute common solutions and gradients by
  !              averaging the left and right interface
  !              solutions and gradients and use these
  !              values to compute the common viscous flux.
  !
  integer :: visc_variables_to_ave = 1
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### TIME STEPPING INFORMATION
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! NUM_TIMESTEPS : number of time steps for the simulation
  ! FINAL_TIME    : final simulation time to run the simulation to
  ! CONSTANT_DT   : global time step when using a constant global time step
  !
  !    if (FINAL_TIME > 0.0) NUM_TIMESTEPS = nint( FINAL_TIME / CONSTANT_DT )
  !
  integer :: num_timesteps = 25000
  real    :: Final_Time    = 10.0
  real    :: constant_dt   = 1.0e-6
  !
  ! TIMESTEP_TYPE : type of time stepping to use for the simulation
  !     0 = constant global time step
  !     1 = global time step from minimum computed cell time step
  !    -1 = local cell time step
  !     2 = global time step from minimum computed solution point time step
  !    -2 = local solution point time step
  !
  integer :: Timestep_Type = 1
  !
  ! MINIMUM_TIMESTEPS : the minimum number of time steps to complete
  !                     before checking for convergence
  !
  integer :: minimum_timesteps = 100
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### CFL CONTROL VARIABLES
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! CFL : base CFL number to use
  !
  real :: cfl = 0.1
  !
  ! CFL_BEG : beginning CFL number to use when CFL_CYCLES is greater than 1
  !
  real :: cfl_beg = 0.1
  !
  ! CFL_END : final CFL number after completing the number of time steps
  !           equal to CFL_CYCLES
  !
  real :: cfl_end = 0.1
  !
  ! CFL_CYCLES : the number of time steps used to ramp up the CFL number
  !              starting from CFL_BEG and ending at CFL_END
  !
  integer :: cfl_cycles = 1
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### OVER-INTEGRATION OPTIONS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! PROJECTION_ORDER : polynomial space (order of accuracy) to project the residual
  !           onto when over-integrating by projecting the residual
  !
  integer :: projection_order = 3
  !
  ! QUADRATURE_ORDER : quadrature order of accuracy when over-integrating by
  !           over-sampling at the quadrature points and projecting back
  !           down to the solution order of accuracy (SOLUTION_ORDER).
  !
  integer :: quadrature_order = 3
  !
  ! LOC_QUADRATURE_PTS : specifies the location of the quadrature points used for
  !              quadrature over-integration
  !      1 = Legendre-Gauss         nodes
  !      2 = Legendre-Gauss-Lobatto nodes
  !
  integer :: loc_quadrature_pts = 1
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### ERROR-INTEGRATION OPTIONS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! ERROR_ORDER : quadrature order to use when computing integrated
  !               error quantities
  !               ONLY USED IF > SOLUTION_ORDER
  !
  integer :: error_order = -3
  !
  ! LOC_ERROR_PTS : specifies the location of the quadrature points used for
  !                 computing integrated error quantities
  !                1 = Legendre-Gauss         nodes
  !                2 = Legendre-Gauss-Lobatto nodes
  !    anything else = NOT USED
  !
  integer :: loc_error_pts = 1
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### FORCE GLOBAL ORDER AND SOLUTION POINTS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! ALL_ORDERS : input parameter to easily set SOLUTION_ORDER, PROJECTION_ORDER, and QUADRATURE_ORDER
  !              all to the same value
  !        ONLY USED IF (all_orders >= 0)
  !
  integer :: all_orders = -3
  !
  ! ALL_POINTS : input parameter to easily set LOC_SOLUTION_PTS, LOC_FLUX_PTS, and
  !              LOC_QUADRATURE_PTS all to the same value.
  !                1 = Legendre-Gauss         nodes
  !                2 = Legendre-Gauss-Lobatto nodes
  !    anything else = NOT USED
  !
  integer :: all_points = -1
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### PARALLEL OPTIONS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! METIS_OPTION : specifies the METIS option used for partitioning the
  !                global grid
  !     1 = non-weighted partitions based on multilevel recursive bisection
  !         while minimizing the edge cut
  !     2 = non-weighted partitions based on multilevel     k-way algorithm
  !         while minimizing the edge cut
  !     3 = non-weighted partitions based on multilevel     k-way algorithm
  !         while minimizing the total communication volume
  !     4 =     weighted partitions based on multilevel     k-way algorithm
  !         while minimizing the total communication volume
  !
  integer :: metis_option = 3
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### METHOD OF MANUFACTURED SOLUTIONS (MMS) OPTIONS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! MMS_OPT : (DOES NOT WORK) specifies the MMS problem to use
  !
  integer :: mms_opt = 0
  !
  ! MMS_INIT : (DOES NOT WORK) initializes the solution to
  !            mms_init*(exact MMS solution)
  !
  real :: mms_init = 0.1
  !
  ! MMS_OUTPUT_SOLUTION : flag used to completely skip the solver and just
  !                       output the exact MMS solution
  !
  logical :: mms_output_solution = .FALSE.
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### INITIALIZATION / TEST CASE
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! ITESTCASE : specifies the test-case / problem to simulation affecting the
  !             initialization of the simulation
  !       1 = Generic_MeanFlow    (initialize to reference conditions)
  !       2 = Channel_Flow        (initialize to ...)
  !       4 = Laminar_BndLyr      (initialize to reference conditions)
  !       5 = Shu_Vortex          (diagonally propagating)
  !      -5 = Shu_Vortex          (stationary vortex)
  !       6 = Vortex_Transport    (propagate based on alpha_aoaref)
  !      -6 = Vortex_Transport    (propagate only in y)
  !       7 = Density_Transport   (diagonally propagating)
  !       8 = Entropy_Transport   (diagonally propagating)
  !       9 = Taylor_Green_Vortex
  !      10 = Nozzle_Jet initialization
  !      11 = Double Mach Reflection (not-working)
  !      12 = Infinite Cylinder
  !      13 = Freestream Preservation
  !
  integer :: itestcase = 9
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### RESTART OPTIONS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! LOAD_RESTART_FILE : flag used to enable/disable using the restart file
  !              .FALSE. = do not use restart file
  !              .TRUE.  = use restart file
  !
  logical :: load_restart_file = .FALSE.
  !
  ! RESTART_INTERVAL : number of time steps between writing a restart file
  !   NOTE: A RESTART FILE IS ALWAYS WRITTEN AFTER COMPLETING THE LAST
  !         TIME STEP ITERATION.
  !
  integer :: restart_interval = 0
  !
  ! USE_OLD_RESTART_FORMAT : flag for using the old restart file format
  !    .TRUE.  = use the old serial   restart file format
  !    .FALSE. = use the new parallel restart file format
  !
  logical :: use_old_restart_format = .FALSE.
  !
  ! RESTART_FILE : name of file to read restart solution from
  !
  character(len=150) :: restart_file = "out.rst"
  !
  ! RESTART_OUTPUT_USES_HOST_ROOTS : logical flag that whether all MPI processes
  !                                  are involved in writing to the restart file
  !                                  or if only the root processes of each host
  !                                  node write to the restart file
  !
  logical :: restart_output_uses_host_roots = .TRUE.
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### STATISTICAL OUTPUT OPTIONS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! ITER_OUT_INTERVAL : interval number of time steps between writing
  !                     iteration statistics to standard output
  !    if (ITER_OUT_INTERVAL == 0) then
  !      ITER_OUT_INTERVAL = 1
  !    if (ITER_OUT_INTERVAL < 0) then
  !      ITER_OUT_INTERVAL = max( 1 , NUM_TIMESTEPS / abs(ITER_OUT_INTERVAL) )
  !    end if
  !
  integer :: iter_out_interval = 1
  !
  ! RESULTS_INTERVAL : interval number of time steps between writing
  !                    convergence statistics to the results file
  !    if (RESULTS_INTERVAL == 0) then
  !      RESULTS_INTERVAL = NUM_TIMESTEPS / 500
  !    if (RESULTS_INTERVAL < 0) then
  !      RESULTS_INTERVAL = max( 1 , NUM_TIMESTEPS / abs(RESULTS_INTERVAL) )
  !    end if
  !
  integer :: results_interval = 0
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### SOLUTION OUTPUT OPTIONS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! OUTPUT_DIR : directory used for writing output files
  !
  character(len=150) :: output_dir = "."
  !
  ! OUTPUT_INTERVAL : interval number of time steps between writing a CGNS
  !                   solution file
  !    if (OUTPUT_INTERVAL == 0) then
  !      OUTPUT_INTERVAL = NUM_TIMESTEPS / 50
  !    if (OUTPUT_INTERVAL < 0) then
  !      OUTPUT_INTERVAL = max( 1 , NUM_TIMESTEPS / abs(OUTPUT_INTERVAL) )
  !    end if
  !
  integer :: output_interval = -5
  !
  ! CONTINUOUS_OUTPUT : enables/disables the post-processor that attempts to
  !                     create a smooth, continuous solution across all cells
  !                     before writing it to the CGNS solution file.
  !
  logical :: continuous_output = .FALSE. ! THIS IS HIGHLY EXPERIMENTAL
  !
  ! OUTPUT_ORDER : order of accuracy used when writing the CGNS solution files
  !                ONLY USED IF (OUTPUT_ORDER > SOLUTION_ORDER)
  !
  integer :: output_order = -3
  !
  ! LOC_OUTPUT_PTS : specifies the location of the solution points
  !                  used for the CGNS solution file
  !      0 = Equi-distant           nodes
  !      1 = Legendre-Gauss         nodes
  !      2 = Legendre-Gauss-Lobatto nodes
  !
  integer :: loc_output_pts = 0
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### GOVERNING EQUATIONS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! GOVERNING_EQUATIONS : specifies the governing equations to use
  !         1 = Euler equations
  !         2 = Navier-Stokes equations
  !
  integer :: governing_equations = 2
  !
  ! TURBULENCE_MODEL : (NOT USED YET) specifies the turbulence model to use
  !         0 = laminar / implicit LES
  !         1 = 1998 k-omega
  !        -1 = 2006 k-omega
  !
  integer :: turbulence_model = 0
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### NONDIMENSIONAL REFERENCE CONDITIONS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! MACHREF : reference Mach number
  ! REYREF  : reference Reynolds number
  ! GAM     : reference gamma (ratio of specific heats)
  !
  real :: machref =  0.4
  real :: reyref  = -1.0
  real :: gam     =  1.4
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### TOTAL REFERENCE CONDITIONS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! PTOTREF      : reference total pressure
  ! TTOTREF      : reference total temperature
  ! PTOT2P_RATIO : ratio of reference total pressure
  !                to reference static pressure
  !
  ! if (PTOT2P_RATIO > 0.0) then
  !   TTOTREF = TREF * (1.0 + 0.5*(GAM-1.0)*MACHREF**2)
  !   PTOTREF = PREF * (1.0 + 0.5*(GAM-1.0)*MACHREF**2)**(GAM/GAM-1.0)
  !   PTOT2P_RATIO = PTOTREF / PREF
  ! else
  !   TTOTREF = TREF * PTOT2P_RATIO**(GAM-1.0/GAM)
  !   PTOTREF = PREF * PTOT2P_RATIO
  ! end if
  !
  real :: ptotref      = -1.0
  real :: ttotref      = -1.0
  real :: ptot2p_ratio = -1.0
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### STATIC REFERENCE CONDITIONS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! RHOREF  : reference static density
  ! TREF    : reference static temperature
  ! PREF    : reference static pressure
  !
  ! SPECIFIY TWO OF THE THREE AND SET THE OTHER NEGATIVE. IF ALL THREE
  ! ARE SPECIFIED, A CHECK ON THE CONSISTENCY BETWEEN THEM WILL BE DONE.
  !
  real :: rhoref  =     -1.0
  real :: tref    =    300.0
  real :: pref    = 100000.0
  !
  ! RGASREF : reference gas constant
  !
  real :: rgasref =    287.15
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### REFERENCE TURBULENCE CONDITIONS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! TINREF : (NOT USED YET) reference turbulence intensity
  !
  real :: tinref = 0.0
  !
  ! TVRREF : (NOT USED YET) reference ratio of eddy viscosity
  !          to molecular viscosity
  !
  real :: tvrref = 1.0e-6
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### REFERENCE ANGLES OF ATTACK
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! ALPHA_AOAREF : reference angle of attack in the x-y plane
  !
  real :: alpha_aoaref = 0.0
  !
  ! BETA_AOAREF : reference angle of attack in the x-z plane
  !
  real :: beta_aoaref = 0.0
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### LAMINAR AND TURBULENT PRANDTL NUMBERS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! PR : laminar Prandtl number
  !
  real :: Pr = 0.72
  !
  ! PRT : (NOT USED YET) turbulent Prandtl number
  !
  real :: Prt = 0.72
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### SUTHERLANDS LAW COEFFICIENTS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! SUTH_MUREF : constant mu0 in Sutherlands law
  !
  real :: suth_muref = 1.716e-5
  !
  ! SUTH_TREF : constant T0 in Sutherlands law
  !
  real :: suth_Tref = 273.0
  !
  ! SUTH_SREF : constant S in Sutherlands law
  !
  real :: suth_Sref = 110.4
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### FILTER OPTIONS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! FILTER_OPTION : specifies the type of filter to use for stabilization
  !       0 = filter is not used
  !   FOR STRUCTURED ELEMENTS
  !       1 = tensor product of the 1D exponential filter where the
  !           1D filter is used separately in each direction and
  !               N = SOLUTION_ORDER
  !       2 = true tensor product of the 1D exponential filter where
  !               N = 2*SOLUTION_ORDER for 2D
  !               N = 3*SOLUTION_ORDER for 3D
  !   FOR UNSTRUCTURED ELEMENTS
  !       anything but 0 = exponential filter where N = SOLUTION_ORDER
  !
  ! Exponential filter is defined such that
  !
  ! etac = Filter_etac = real(Nc)/real(N)
  ! Nc = filter cutoff (0 <= Nc <= N)
  ! a = Filter_Alpha
  ! s = Filter_Order
  ! eps = Filter_eps
  !
  ! filter_diagonal(:) = 0.0
  ! do i = Nc,N
  !   eta = real(i)/real(N)
  !   filter_diagonal(i) = (1-eps) * exp( -a * ( (eta-etac)/(1.0-etac) )**s )
  ! end do
  !
  integer :: Filter_Option = 0
  !
  ! FILTER_ETAC : filter cutoff below which the low modes are left untouched
  !
  real :: Filter_etac = 0.0
  !
  ! FILTER_ORDER : order of the exponential filter
  !
  integer :: Filter_Order = 16
  !
  ! FILTER_ALPHA : maximum damping parameter
  !
  real :: Filter_Alpha = 36.0
  !
  ! FILTER_EPS : leading damping parameter
  !
  real :: Filter_eps = 0.0
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### LIMITER FILTER OPTIONS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! LIMITER_FILTER_OPTION : details are the same as for FILTER_OPTION, but this
  !                         is a second filter that is only used as a limiter
  !                         on any cells that have been marked as troubled
  !
  !  The idea is you can use FILTER_OPTION as a very weak filter that is
  !  constantly being applied regardless of whether or not the solution is
  !  oscillatory within the cell, and you can use LIMITER_FILTER_OPTION as
  !  a stronger filter that is only applied to cells that have been marked
  !  as troubled to prevent the oscillatory solutions in these cells from
  !  causing the simulation to diverge.
  !
  integer :: Limiter_Filter_Option = 0
  !
  ! LIMITER_FILTER_ETAC : same as FILTER_ETAC but for the limiter filter
  !
  real :: Limiter_Filter_etac = 0.0
  !
  ! LIMITER_FILTER_ORDER : same as FILTER_ORDER but for the limiter filter
  !
  integer :: Limiter_Filter_Order = 16
  !
  ! LIMITER_FILTER_ALPHA : same as FILTER_ALPHA but for the limiter filter
  !
  real :: Limiter_Filter_Alpha = 36.0
  !
  ! LIMITER_FILTER_EPS : same as FILTER_EPS but for the limiter filter
  !
  real :: Limiter_Filter_eps = 0.0
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### LIMITER OPTIONS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! LIMITER_OPTION : specifies what limiter to use on troubled cells
  !
  !       0 = do not use limiter
  !
  !      -1 = use the resolution indicator to mark troubled cells, and for each
  !           marked cell, project the solution to a lower order based on the
  !           value of the indicator for the cell
  !
  !     For values > 0:
  !
  !       1 = use the resolution indicator to mark troubled cells
  !       2 = use the jump indicator to mark troubled cells
  !
  !       if (LIMITER_FILTER_OPTION /= 0) then
  !         use the indicator to mark troubled cells
  !         use the limiter filter on all marked cells
  !       else
  !         use the indicator to mark troubled cells
  !         reverse loop from solution order to 1
  !           project marked cells to one order less than the current loop order
  !           use the indicator to recheck all troubled cells
  !           exit reverse loop if no more troubled cells exist
  !         end reverse loop
  !       end if
  !
  integer :: Limiter_Option = 0
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### CONTROL OVER FREQUENCY OF FILTERING/LIMITING
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! FILTERING_INTERVAL : frequency of applying the standard filter
  !
  !        0 = apply standard filter after each Runge-Kutta stage
  !      > 0 = apply standard filter after the last Runge-Kutta stage of every
  !            FILTERING_INTERVAL-th time step
  !      < 0 = standard filter will not be applied even if FILTER_OPTION /= 0
  !
  integer :: Filtering_Interval = 0
  !
  ! LIMITING_INTERVAL : frequency of applying the limiter filter
  !
  !        0 = apply limiter filter after each Runge-Kutta stage
  !      > 0 = apply limiter filter after the last Runge-Kutta stage of every
  !            LIMITING_INTERVAL-th time step
  !      < 0 = limiter filter will not be applied even if LIMITER_OPTION /= 0
  !
  integer :: Limiting_Interval = 0
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### ISENTROPIC VORTEX CONSTANTS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! VORTEX_BIGR : constant used for the Euler vortex problem
  !
  real :: vortex_bigR = 1.0
  !
  ! VORTEX_BIGGAM : constant used for the Euler vortex problem
  !
  real :: vortex_bigGam = 5.0
  !
  ! VORTEXGRID_LREF : reference length for the grid domain, grid is a square
  !                   with bounds of [-VortexGrid_Lref,+VortexGrid_Lref] in
  !                   both the X and Y directions
  !
  real :: VortexGrid_Lref = 10.0
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### TAYLOR-GREEN VORTEX CONSTANTS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! VORTEXGRID_DOF : number of solution points (degrees-of-freedom/DoF) in
  !                  each direction to use for the Taylor-Green vortex
  !                  problem
  !
  integer :: VortexGrid_DOF = 48
  !
  ! VORTEXGRID_RANDOM_PERTURB : flag for applying a random perturbation to
  !                             all the interior grid nodes for the
  !                             Taylor-Green vortex problem
  !
  logical :: VortexGrid_random_perturb = .FALSE.
  !
  ! VORTEXGRID_PERTURB_PERCENT : percent of the uniform grid cell length to
  !                              perturb the grid nodes
  !
  integer :: VortexGrid_perturb_percent = 30
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### CHANNEL FLOW PARAMETERS (NOT WORKING)
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! CHANNEL_BODY_FORCE : body force used to drive the channel flow problem
  !
  real :: channel_body_force = 0.0 ! NOT WORKING
  !
  ! CHANNEL_INIT_FILE : file used to initialize the solution
  !                     for the channel flow problem
  !
  character(len=150) :: channel_init_file = "" ! NOT WORKING
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### OPTIONS MODIFYING BC ALGORITHMS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! WALLS_ARE_EXACT : flag determining what solution is used on wall
  !                   boundary faces
  !    .TRUE.  = wall boundary solution is exact
  !        ([u_wall,v_wall,w_wall] = 0 for no-slip wall)
  !    .FALSE. = average of wall boundary solution and interior solution on
  !              the boundary face recovers exact wall solution
  !        ([u_wall,v_wall,w_wall] = -[u_int,v_int,w_int] for no-slip wall)
  !
  logical :: walls_are_exact = .TRUE.
  !
  ! SUB_INFLOW_METHOD : determines which subsonic inflow BC to use
  !
  !       1 = hold the total pressure and temperature constant at the inflow,
  !           use the outgoing characteristic from the interior to perform a
  !           Newton iteration to get the static temperature, speed of
  !           sound, and normal velocity magnitude for the exterior state
  !
  !       2 = hold inflow density and velocity constant,
  !           and use interior pressure
  !
  !       3 = hold total pressure and temperature constant at the inflow,
  !           and use interior velocity
  !
  !       4 = hold total pressure and temperature constant at the inflow,
  !           and use interior pressure
  !
  integer :: sub_inflow_method = 1
  !
  ! USE_BC_CONFLICTS : flag for determining whether or not to overwrite
  !                    boundary solutions in adjacent cells if co-located
  !                    solution points have different boundary conditions
  !                    NOTE: only relevant if using Lobatto points for
  !                          the solution points
  !
  !
  logical :: use_bc_conflicts = .TRUE.
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### CGNS OUTPUT OPTIONS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! CGNS_USE_QUEUE : flag enabling/disabling the CGNS queuing system when
  !                  writing parallel CGNS files
  !
  logical :: cgns_use_queue = .FALSE.
  !
  ! INTERPOLATE_BEFORE_OUTPUT_VARIABLE : flag determining whether the solution
  !                                      variables are computed before or after
  !                                      interpolating to the output nodes
  !
  logical :: interpolate_before_output_variable = .TRUE.
  !
  ! CONVERT_RESTART_TO_CGNS : flag used to read in a restart file and
  !                           immediately write a CGNS solution file
  !
  logical :: convert_restart_to_cgns = .FALSE.
  !
  ! COMPLETELY_DISABLE_CGNS : flag used to completely disable CGNS solutions
  !                           from being output at the output interval
  !
  logical :: completely_disable_cgns = .FALSE.
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### CONVERGENCE CRITERIA
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! CONVERGENCE_ORDER_ABS_RES : Orders of magnitude of the absolute residual
  !
  real :: convergence_order_abs_res = -16.0
  !
  ! CONVERGENCE_ORDER_MAX_RES : Orders of magnitude reduction of the residual
  !                             relative to the maximum residual
  !
  real :: convergence_order_max_res = -16.0
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### PROFILE PARALLEL IO
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! PROFILE_IO_ALL : if true, sets both profile_io_cgns
  !                  and profile_io_restart to true
  !
  logical :: profile_io_all = .FALSE.
  !
  ! PROFILE_IO_CGNS : flag to activate some simple timers that time how long
  !                   it takes to write each CGNS file, and output the results
  !                   to stdout.
  !
  logical :: profile_io_cgns = .FALSE.
  !
  ! PROFILE_IO_RESTART : flag to activate some simple times that time how long
  !                      it takes to write each restart or time-averaged restart
  !                      file, and output the results to stdout.
  !
  logical :: profile_io_restart = .FALSE.
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### OPTIONS TO OUTPUT BOUNDARY PROFILE
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! OUTPUT_BND_PROFILES : flag used to output the solution profiles on
  !                       certain boundary faces
  !
  logical :: output_bnd_profiles = .FALSE.
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### CONTROL DOT_PRODUCT PROCEDURES BOUND TO MATRIX DERIVED-TYPE
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! USE_UNROLLED_DOT_PRODUCTS : flag for enabling/disabling the manually
  !                             unrolled dot product functions
  !
  logical :: use_unrolled_dot_products = .TRUE.
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### LUSTRE FILE SYSTEM STRIPING PARAMTERS
  ! #####
  ! #####   ONLY APPLICABLE WHEN RUNNING GFR SO THAT
  ! #####   IT WRITES SOLUTION/RESTART FILES TO A
  ! #####   LUSTRE FILE SYSTEM
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! LUSTRE_STRIPE_COUNT : specifies the approximate number of stripes to use
  !                       when writing a solution or restart file
  !                       NOTE: it will take this value and find the closest
  !                             integer value satisfying
  !                                 mod(ncpu,lustre_stripe_count) == 0
  !
  integer :: lustre_stripe_count = 64
  !
  ! LUSTRE_STRIPE_SIZE : specifies the number of bytes to store on each stripe
  !                      before moving to the next stripe.
  !                      NOTE: must be a multiple of 65,536 bytes (64 KB) and
  !                            'k', 'm', or 'g' can be used to specify units of
  !                            KB, MB, or GB, respectively
  !                            for example: 4m => 4 MB => 4,194,304 bytes
  !
  !
  character(len=10) :: lustre_stripe_size = "4m"
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### TIME AVERAGING
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! OUTPUT_TIME_AVERAGING :
  !
  logical :: output_time_averaging = .FALSE.
  !
  ! TIME_AVE_FILE :
  !
  character(len=150) :: time_ave_file = "out.ave"
  !
  ! TIME_SCALING_FACTOR :
  !
  real :: time_scaling_factor = 1.0
  !
  ! TIME_AVE_VEL_IS_AXISYMM :
  !
  logical :: time_ave_vel_is_axisymm = .FALSE.
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### UTILITY TO TIME AVERAGE EXISTING RESTART FILES
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! TIME_AVERAGE_RESTART_FILES :
  !
  logical :: time_average_restart_files = .FALSE.
  !
  ! RESTART_LIST_FILE :
  !
  character(len=150) :: restart_list_file = ""
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### MEMORY USAGE
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! DUMP_MEMORY_USAGE : flag used to profile the memory usage
  !
  logical :: dump_memory_usage = .FALSE.
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### DUMP FLUX INFORMATION FOR DEBUGGING
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! DUMP_FLUXES : flag used to output the fluxes computed in the solver
  !
  logical :: dump_fluxes = .FALSE.
  !
  ! CHECK_FLUX_POINT_COORDINATES : flag used to:
  !
  ! 1) loop over each face
  !      loop over each flux point on this face
  !        interpolate left  cell coordinates to this flx-pt and write to file 1
  !        interpolate right cell coordinates to this flx-pt and write to file 2
  !      end loop over flux points
  !    end loop over faces
  !
  ! 2) loop over each cell
  !      loop over each cell face
  !        loop over each flux point on this face
  !          if cell is left  cell on this face
  !            interpolate cell coordinates to face_xyz(face_idx,fp_idx,1)
  !          if cell is right cell on this face
  !            interpolate cell coordinates to face_xyz(face_idx,fp_idx,2)
  !        end loop over flux points
  !      end loop over cell faces
  !    end loop over cells
  !    write face_xyz(:,:,1) to file 3
  !    write face_xyz(:,:,2) to file 4
  !
  ! If the coordinates match between files 1 and 2, and match between
  ! files 3 and 4, then the grid connectivity should be correct and
  ! co-located face quantities from the two cells on the face are
  ! consistent in terms of interpolation.
  !
  logical :: check_flux_point_coordinates = .FALSE.
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### ACTIVATE SPECIAL GRID FOR DEBUGGING THE POSTPROCESSOR
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! POSTPROC_DEBUG_GRID : flag to activate using an internally created grid
  !                       that is useful for debugging the post-processor
  !                       used for the CONTINUOUS_OUTPUT flag
  !
  logical :: postproc_debug_grid = .FALSE.
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### OPTION TO OUTPUT PLOT3D GRID TO TECPLOT IN UNSTRUCTURED FORMAT
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! OUTPUT_PLOT3D_TO_TECPLOT : flag to output a Tecplot file of the input
  !                            Plot3D grid file after reading it
  !
  logical :: output_plot3d_to_tecplot = .FALSE.
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### OPTION TO OUTPUT BOUNDARY CONDITIONS TO A TECPLOT FILE
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! OUTPUT_BFACE_ARRAY : flag to output the boundary conditions
  !
  logical :: output_bface_array = .FALSE.
  !
  ! ####################################################
  ! ####################################################
  ! #####
  ! ##### INPUT BOUNDARY CONDITIONS
  ! #####
  ! ####################################################
  ! ####################################################
  !
  ! BC_INPUT : derived type to specify boundary conditions
  !
  type(bc_input_t) :: bc_input(1:20)
  !
  type :: bc_input_t
    !
    ! BC_INPUT(i)%NAME : character string that matches the name of a boundary
    !                    group in the grid file
    !                    NOTE: case (upper/lower-case) does NOT matter
    !
    character(len=32) :: name = ""
    !
    ! BC_INPUT(*)%SET_BC_VALUE : overrides the boundary condition type specified
    !                            in the grid file with this type
    !
    !          Valid boundary condition values for bc_input(*)%set_bc_value
    !
    !              1 = generic inflow, subsonic or supersonic
    !              2 = subsonic inflow
    !              3 = supersonic inflow
    !             11 = generic outflow, subsonic or supersonic
    !             12 = subsonic outflow
    !             13 = supersonic outflow
    !             21 = characteristic inflow/outflow based on incoming/outgoing
    !                  Rieman invariants, can be subsonic or supersonic
    !             22 = generic inflow/outflow, subsonic or supersonic
    !             23 = fixed to freestream/reference flow conditions
    !             24 = fixed flow conditions, equivalent to supersonic inflow
    !             50 = slip wall
    !             52 = adiabatic no-slip wall
    !             53 = isothermal no-slip wall
    !             80 = symmetry boundary condition, equivalent to slip wall
    !
    integer :: set_bc_value = -9999
    !
    ! BC_INPUT(*)%T_STATIC     : static temperature
    ! BC_INPUT(*)%P_STATIC     : static pressure
    ! BC_INPUT(*)%RHO_STATIC   : static density
    ! BC_INPUT(*)%MACH         : Mach number
    ! BC_INPUT(*)%VX           : velocity in x-direction
    ! BC_INPUT(*)%VY           : velocity in y-direction
    ! BC_INPUT(*)%VZ           : velocity in z-direction
    ! BC_INPUT(*)%T_TOTAL      : total temperature
    ! BC_INPUT(*)%P_TOTAL      : total pressure
    ! BC_INPUT(*)%RHO_TOTAL    : total density
    ! BC_INPUT(*)%ALPHA_AOA    : angle of attack w.r.t. the x-y plane
    ! BC_INPUT(*)%BETA_AOA     : angle of attack w.r.t. the x-z plane
    ! BC_INPUT(*)%WALL_TEMP    : temperature for isothermal wall
    !
    !       An iterative loop will go through the specified flow conditions
    !       and try to compute the remaining unspecified flow conditions.
    !       A run-time error will occur if:
    !         A) the flow conditions are over-specified and inconsistent, or
    !         B) the flow conditions are under-specified and it is unable to
    !            compute the remaining unspecified flow conditions
    !
    !       An example of a well defined boundary condition using bc_input:
    !
    !        bc_input(1)%name = "INFLOWSUBSONIC"
    !        bc_input(1)%set_bc_value = 2 ! Subsonic Inflow BC
    !        bc_input(1)%p_static = 101325.0
    !        bc_input(1)%p_total = 121286.025
    !        bc_input(1)%t_total = 293.15
    !
    real :: t_static   = -1.0
    real :: p_static   = -1.0
    real :: rho_static = -1.0
    real :: mach       = -1.0
    real :: vx         =  0.0
    real :: vy         =  0.0
    real :: vz         =  0.0
    real :: t_total    = -1.0
    real :: p_total    = -1.0
    real :: rho_total  = -1.0
    real :: alpha_aoa  =  0.0
    real :: beta_aoa   =  0.0
    real :: wall_temp  = -1.0
    !
    ! BC_INPUT(*)%RELAX : experimental relaxation parameters that only
    !                     work for a few BCs
    !
    type(relaxation_t) :: relax
    !
  end type bc_input_t
  !
  type :: relaxation_t
    !
    ! BC_INPUT(*)%RELAX%VALUE  : relaxation value to limit the difference
    !                            between the interior and exterior states
    !                             0.0 < BC_INPUT(*)%RELAX%VALUE <= 1.0
    !
    real :: value = 1.0
    !
    ! BC_INPUT(*)%RELAX%TIME   : simulation time to end relaxation ( > 0.0)
    !
    real :: time = 0.0
    !
    ! BC_INPUT(*)%RELAX%ITER   : time-step to end relaxation ( > 0)
    !
    integer :: iter = 0
    !
    !       The value rlx is determined from the minimum relaxation value from the
    !       three relaxation parameters. With
    !          ITCUR = current iteration / time-step
    !          TIME  = current time within the simulation
    !
    !       rlx = 1.0
    !       rlx_iter = 1.0
    !       rlx_time = 1.0
    !       if (relax%iter > 0) then
    !         rlx_iter = max( 1.0 , real(itcur)/real(relax%iter) )
    !       end if
    !       if (relax%time > 0.0) then
    !         rlx_time = max( 1.0 , time/relax%time )
    !       end if
    !       rlx = max( 0.0 , min(relax%value,rlx_iter,rlx_time,1.0) )
    !
    !       bc-state = (1.0-rlx)*interior-state + rlx*exterior-state
    !
  end type relaxation_t
