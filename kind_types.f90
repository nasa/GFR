module rename_intrinsics
  !
  implicit none
  !
  public
  !
  intrinsic dot_product
  intrinsic matmul
  !
end module rename_intrinsics
!
!###############################################################################
!
module module_kind_types
  !
  !.. Use Statements ..
  use mpi
#ifdef USE_INTEL_MKL
  use blas95, only : dot
#else
  use rename_intrinsics, only : dot => dot_product
#endif
  !
  implicit none
  !
  public
  !
  ! ####################################
  ! ####################################
  ! #####   KIND TYPE PARAMETERS   #####
  ! ####################################
  ! ####################################
  !
  ! KIND definitions for single or double precision real variables (IEEE 754)
  !
  integer, parameter :: i4 = selected_int_kind(9)
  integer, parameter :: i8 = selected_int_kind(18)
  integer, parameter :: r4 = selected_real_kind(6,37)
  integer, parameter :: r8 = selected_real_kind(15,307)
  integer, parameter :: r16 = max( r8 , selected_real_kind(33, 4931) )
  !
  ! Aliases for kind type parameters
  !
  integer, parameter :: SP = r4
  integer, parameter :: DP = r8
  integer, parameter :: QP = r16
  integer, parameter :: WP = r8
  !
  integer, parameter :: WIP = i4
  integer, parameter :: INT_MPI = kind(MPI_INTEGER_KIND)
  integer, parameter :: INT_CGNS = i4
  integer, parameter :: INT_METIS = i4
  !
  ! ###########################################
  ! ###########################################
  ! #####   NUMERICAL/LOGICAL CONSTANTS   #####
  ! ###########################################
  ! ###########################################
  !
  ! Definitions for numerical contants
  !
  real(qp), parameter :: qzero   = 0.0_qp
  real(qp), parameter :: qhalf   = 0.5_qp
  real(qp), parameter :: qone    = 1.0_qp
  real(qp), parameter :: qtwo    = 2.0_qp
  real(qp), parameter :: qten    = 10.0_qp
  !
  real(wp), parameter :: zero    = 0.0_wp
  real(wp), parameter :: half    = 0.5_wp
  real(wp), parameter :: one     = 1.0_wp
  real(wp), parameter :: s1      = 1.2_wp
  real(wp), parameter :: s2      = 0.8_wp
  real(wp), parameter :: two     = 2.0_wp
  real(wp), parameter :: three   = 3.0_wp
  real(wp), parameter :: four    = 4.0_wp
  real(wp), parameter :: five    = 5.0_wp
  real(wp), parameter :: six     = 6.0_wp
  real(wp), parameter :: seven   = 7.0_wp
  real(wp), parameter :: eight   = 8.0_wp
  real(wp), parameter :: nine    = 9.0_wp
  real(wp), parameter :: ten     = 10.0_wp
  real(wp), parameter :: twelve  = 12.0_wp
  real(wp), parameter :: sixteen = 16.0_wp
  !
  real(wp), parameter :: one8   = one/eight
  real(wp), parameter :: one6   = one/six
  real(wp), parameter :: one5   = one/five
  real(wp), parameter :: one4   = one/four
  real(wp), parameter :: one3   = one/three
  real(wp), parameter :: two3   = two/three
  real(wp), parameter :: three4 = three/four
  !
  real(wp), parameter :: eps1  = 1.0e-01_wp
  real(wp), parameter :: eps2  = 1.0e-02_wp
  real(wp), parameter :: eps3  = 1.0e-03_wp
  real(wp), parameter :: eps6  = 1.0e-06_wp
  real(wp), parameter :: eps7  = 1.0e-07_wp
  real(wp), parameter :: eps9  = 1.0e-09_wp
  real(wp), parameter :: eps10 = 1.0e-10_wp
  real(wp), parameter :: eps11 = 1.0e-11_wp
  real(wp), parameter :: eps12 = 1.0e-12_wp
  real(wp), parameter :: eps14 = 1.0e-14_wp
  real(wp), parameter :: eps17 = 1.0e-17_wp
  real(wp), parameter :: small = eps17
  !
  real(wp), parameter :: pi = 3.14159265358979323846_wp
  !
  real(wp), parameter :: radians2degrees = two*nine*ten/pi
  real(wp), parameter :: degrees2radians = one/radians2degrees
  !
  real(wp), parameter :: rndoff  = epsilon(one)
  real(wp), parameter :: tny     = tiny(one)
  !
  logical, parameter :: true = .true.
  logical, parameter :: fals = .false.
  !
  ! #########################################
  ! #########################################
  ! #####   GENERAL GLOBAL PARAMETERS   #####
  ! #########################################
  ! #########################################
  !
  ! I/O constants
  !
  integer, parameter :: iout = 6
  integer, parameter :: inpt = 7
  integer, parameter :: iorst = 10
  !
  ! Boundary condition constants
  !
  ! Inflow BCs
  integer, parameter :: bc_generic_inflow = 1
  integer, parameter :: bc_sub_inflow = 2
  integer, parameter :: bc_sup_inflow = 3
  integer, parameter :: bc_mdot_inflow = 4
  ! Outflow BCs
  integer, parameter :: bc_generic_outflow = 11
  integer, parameter :: bc_sub_outflow = 12
  integer, parameter :: bc_sup_outflow = 13
  integer, parameter :: bc_mdot_outflow = 14
  ! Generic Inflow/Outflow BCs
  integer, parameter :: bc_characteristic = 21
  integer, parameter :: bc_generic_freeflow = 22
  integer, parameter :: bc_freestream = 23
  integer, parameter :: bc_fixed = 24
  ! Wall BCs
  integer, parameter :: bc_slip_wall = 50
  integer, parameter :: bc_euler_wall = 51
  integer, parameter :: bc_adiabatic_wall = 52
  integer, parameter :: bc_isothermal_wall = 53
  integer, parameter :: bc_default_wall = bc_adiabatic_wall
  ! MMS BCs
  integer, parameter :: bc_mms_dirichlet = 70
  ! Other BCs
  integer, parameter :: not_a_bc = 0
  integer, parameter :: bc_symmetry = 80
  integer, parameter :: bc_periodic = 81
  integer, parameter :: bc_cpu_bnd = 99
  integer, parameter :: bc_unknown = -9999
  ! BC Groups
  integer, parameter :: bc_comm(1:2) = [bc_periodic,bc_cpu_bnd]
  integer, parameter :: bc_walls(1:5) = [ bc_slip_wall, &
                                          bc_euler_wall, &
                                          bc_symmetry, &
                                          bc_adiabatic_wall, &
                                          bc_isothermal_wall ]
  !
  integer, parameter :: bc_heirarchy(1:17) = [ bc_isothermal_wall, &
                                               bc_adiabatic_wall, &
                                               bc_euler_wall, &
                                               bc_slip_wall, &
                                               bc_symmetry, &
                                               bc_fixed, &
                                               bc_freestream, &
                                               bc_sub_inflow, &
                                               bc_sup_inflow, &
                                               bc_mdot_inflow, &
                                               bc_generic_inflow, &
                                               bc_sub_outflow, &
                                               bc_sup_outflow, &
                                               bc_mdot_outflow, &
                                               bc_generic_outflow, &
                                               bc_characteristic, &
                                               bc_generic_freeflow ]
  integer, parameter :: bc_mms(1:1) = [bc_mms_dirichlet]
  !
  ! Parameters defining workshop test cases
  !
  integer, parameter :: Generic_MeanFlow = 1
  integer, parameter :: Channel_Flow     = 2
  integer, parameter :: NACA0012_Airfoil = 3
  integer, parameter :: Laminar_BndLyr   = 4
  integer, parameter :: Shu_Vortex       = 5
  integer, parameter :: Vortex_Transport = 6
  integer, parameter :: Density_Transport = 7
  integer, parameter :: Entropy_Transport = 8
  integer, parameter :: Taylor_Green_Vortex = 9
  integer, parameter :: Advec_Problems(1:2) = [ Density_Transport, &
                                                Entropy_Transport ]
  integer, parameter :: Vortex_Problems(1:2) = [ Shu_Vortex, &
                                                 Vortex_Transport ]
  integer, parameter :: Transport_Problems(1:4) = [ Shu_Vortex, &
                                                    Vortex_Transport, &
                                                    Advec_Problems ]
  !
  ! Location of 1D solution points
  integer, parameter :: Legendre_Gauss          = 1
  integer, parameter :: Legendre_Gauss_Lobatto  = 2
  integer, parameter :: Legendre_Nodes(1:2) = [Legendre_Gauss, &
                                               Legendre_Gauss_Lobatto]
  integer, parameter :: Chebyshev_Gauss         = 3
  integer, parameter :: Chebyshev_Gauss_Lobatto = 4
  integer, parameter :: Gauss_Lobatto(1:2) = [Legendre_Gauss_Lobatto, &
                                              Chebyshev_Gauss_Lobatto]
  !
  ! Method for creating triangle solution points
  integer, parameter :: DG_TriPoints = 1
  integer, parameter :: HT_TriPoints = 2
  !
  integer, parameter :: LagrangePolynomial = 1
  integer, parameter :: ChainRuleCV = 2
  integer, parameter :: ChainRuleUtd = 3
  integer, parameter :: DGWeakForm = 4
  integer, parameter :: DGStrongForm = 5
  integer, parameter :: TestQuadrature = 0
  integer, parameter :: TestViscous = 10
  integer, parameter :: OnPrimitiveVariables = 1
  integer, parameter :: OnConservativeVariables = 2
  integer, parameter :: OnFluxes = 3
  !
  ! Choices for correction function
  integer, parameter :: Correction_gDG = 1
  integer, parameter :: Correction_g2 = 2
  integer, parameter :: Correction_gGauss = 3
  !
  integer, parameter :: project_derivatives = 1
  integer, parameter :: derivative_of_proj_fluxes = 2
  !
  integer, parameter :: invs_flux_ROE   = 1
  integer, parameter :: invs_flux_HLLC  = 2
  integer, parameter :: invs_flux_LDFSS = 3
  !
  integer, parameter :: visc_flux_BR1 = 1
  integer, parameter :: visc_flux_BR2 = 2
  integer, parameter :: visc_flux_LDG = 3
  !
  ! Parameters giving what set of variables to
  ! average when computing the common viscous flux.
  !
  ! VISC_INTERFACE_FLUXES : Compute the common viscous flux by averaging the
  !                         left and right viscous interface fluxes which are
  !                         computed using the left and right interface
  !                         solutions and gradients.
  ! VISC_SOLUTION_AND_GRADIENTS : Compute common solutions and gradients by
  !                               averaging the left and right interface
  !                               solutions and gradients and use these
  !                               values to compute the common viscous flux.
  !
  integer, parameter :: visc_interface_fluxes = 1
  integer, parameter :: visc_solution_and_gradients = 2
  !
  ! Constants for timestepping method
  !
  integer, parameter :: Constant_Timestep = 0
  !
  integer, parameter :: Cell_Timestep = 1
  integer, parameter :: Global_Cell_Timestep = Cell_Timestep
  integer, parameter :: Local_Cell_Timestep = -Cell_Timestep
  !
  integer, parameter :: Point_Timestep = 2
  integer, parameter :: Global_Point_Timestep = Point_Timestep
  integer, parameter :: Local_Point_Timestep = -Point_Timestep
  !
  ! RK Parameters
  integer, parameter :: Classic_RK = 1
  integer, parameter :: TVD2_RK    = 2
  integer, parameter :: TVD3_RK    = 3
  integer, parameter :: TVD4_RK    = 4
  integer, parameter :: CK4_RK     = 5
  !
  ! Equation Types
  integer, parameter :: Euler_Eqns = 1
  integer, parameter :: NavierStokes_Eqns = 2
  !
  ! Turbulence Options
  integer, parameter :: laminar_flow = 0
  integer, parameter :: k_omega      = 1
  integer, parameter :: k_omega_1998 = +k_omega
  integer, parameter :: k_omega_2006 = -k_omega
  !
  ! Grid Types
  integer, parameter :: Gambit_Neutral = 1
  integer, parameter :: Gmsh_Format = 2
  integer, parameter :: Plot3D_Format = 3
  !
  ! Cell Geometry Types
  integer, parameter :: Geom_Node = 0
  integer, parameter :: Geom_Edge = 1
  integer, parameter :: Geom_Tria = 2
  integer, parameter :: Geom_Quad = 3
  integer, parameter :: Geom_Tetr = 4
  integer, parameter :: Geom_Pyra = 5
  integer, parameter :: Geom_Pris = 6
  integer, parameter :: Geom_Unknown = 7
  integer, parameter :: Geom_Hexa = 8
  !
  integer, parameter :: Geom_Min = min(Geom_Unknown, &
                                       Geom_Node,Geom_Edge, &
                                       Geom_Tria,Geom_Quad, &
                                       Geom_Tetr,Geom_Pyra,Geom_Pris,Geom_Hexa)
  !
  integer, parameter :: Geom_Max = max(Geom_Unknown, &
                                       Geom_Node,Geom_Edge, &
                                       Geom_Tria,Geom_Quad, &
                                       Geom_Tetr,Geom_Pyra,Geom_Pris,Geom_Hexa)
  !
  integer, parameter :: Geom_0D(1:1) = [Geom_Node]
  integer, parameter :: Geom_1D(1:1) = [Geom_Edge]
  integer, parameter :: Geom_2D(1:2) = [Geom_Tria,Geom_Quad]
  integer, parameter :: Geom_3D(1:4) = [Geom_Tetr,Geom_Pyra,Geom_Pris,Geom_Hexa]
  integer, parameter :: Geom_Valid(1:8) = [Geom_0D,Geom_1D,Geom_2D,Geom_3D]
  !
  integer, parameter :: geom_of_face(1:4) = [Geom_Node,Geom_Edge, &
                                             Geom_Tria,Geom_Quad]
  !
  integer, protected, save :: geom_dimen(Geom_Min:Geom_Max)
  integer, protected, save :: geom_nodes(Geom_Min:Geom_Max)
  integer, protected, save :: geom_edges(Geom_Min:Geom_Max)
  integer, protected, save :: geom_faces(Geom_Min:Geom_Max)
  integer, protected, save :: cell_face_geom(1:6,Geom_Min:Geom_Max)
  integer, protected, save :: num_face_nodes(1:6,Geom_Min:Geom_Max)
  integer, protected, save :: num_node_edges(1:8,Geom_Min:Geom_Max)
  integer, protected, save :: edge_nodes(1:2,1:12,Geom_Min:Geom_Max)
  integer, protected, save :: side_nodes(1:4,1:6,Geom_Min:Geom_Max)
  integer, protected, save :: node_edges(1:4,1:8,Geom_Min:Geom_Max)
  integer, protected, save :: node_nodes(1:4,1:8,Geom_Min:Geom_Max)
  !
  character(len=100), protected, save :: Geom_String(Geom_Min:Geom_Max)
  character(len=9),   protected, save :: Geom_Name(Geom_Min:Geom_Max)
  !
  integer, parameter :: Zmin3D_Face = 1
  integer, parameter :: Zmax3D_Face = 2
  integer, parameter :: Xmax3D_Face = 3
  integer, parameter :: Xmin3D_Face = 4
  integer, parameter :: Ymin3D_Face = 5
  integer, parameter :: Ymax3D_Face = 6
  !
  integer, parameter :: Ymin2D_Face = 1
  integer, parameter :: Xmax2D_Face = 2
  integer, parameter :: Ymax2D_Face = 3
  integer, parameter :: Xmin2D_Face = 4
  integer, parameter :: Zmin2D_Face = 0
  integer, parameter :: Zmax2D_Face = 0
  !
  integer, parameter :: Xmin1D = 0
  integer, parameter :: Xmax1D = 1
  !
  integer, parameter :: Xmin2D(1:2) = [3,0]
  integer, parameter :: Xmax2D(1:2) = [1,2]
  integer, parameter :: Ymin2D(1:2) = [0,1]
  integer, parameter :: Ymax2D(1:2) = [2,3]
  !
  integer, parameter :: Xmin3D(1:4) = [0,4,7,3]
  integer, parameter :: Xmax3D(1:4) = [1,2,6,5]
  integer, parameter :: Ymin3D(1:4) = [0,1,5,4]
  integer, parameter :: Ymax3D(1:4) = [3,7,6,2]
  integer, parameter :: Zmin3D(1:4) = [0,3,2,1]
  integer, parameter :: Zmax3D(1:4) = [4,5,6,7]
  !
  integer, parameter :: HostFace2D(1:6) = [Xmin2D_Face,Xmax2D_Face, &
                                           Ymin2D_Face,Ymax2D_Face, &
                                           Zmin2D_Face,Zmax2D_Face]
  integer, parameter :: HostFace3D(1:6) = [Xmin3D_Face,Xmax3D_Face, &
                                           Ymin3D_Face,Ymax3D_Face, &
                                           Zmin3D_Face,Zmax3D_Face]
  !
  ! Minimum/Maximum Order Allowed
  integer, parameter :: Order_Min = 0
  integer, parameter :: Order_Max = 25
  !
  ! tria_is_biunit : logical that changes the location of the reference
  !                  nodes in the definition of the reference triangle
  !         = true  -> defined by the three nodes (-1,-1) - ( 1,-1) - (-1, 1)
  !         = false -> defined by the three nodes ( 0, 0) - ( 1, 0) - ( 0, 1)
  !
  logical, parameter :: tria_is_biunit = fals
  !
  ! CL : coefficient that gives (2.0 / length of the 1D solution points).
  !      This is needed because of the difference between how the reference
  !      triangle and reference quad are defined.
  !      = 2.0 for triangles if tria_is_biunit is false
  !      = 1.0 for everything else
  !
  real(wp), protected, save :: cl(Geom_Min:Geom_Max)
  !
  ! General character variable for storing error messages
  !
  character(len=1000), save :: error_message
  !
  ! MPI Scalars 
  !
  integer(int_mpi), protected, save :: mypnum = 0_int_mpi
  integer(int_mpi), protected, save :: ncpu = 1_int_mpi
  integer(int_mpi),            save :: mpierr
  !
  integer(int_mpi), parameter :: abort = int( 8 , kind=kind(abort) )
  integer(int_mpi), parameter :: stop_mpi = int( 9 , kind=kind(stop_mpi) )
  !
  ! MPI datatypes to match Fortran intrinsic datatypes
  !
  ! working real precision
  integer(int_mpi), protected, save :: mpi_flttyp
  ! half working real precision
  integer(int_mpi), protected, save :: mpi_hlftyp
  ! double working real precision
  integer(int_mpi), protected, save :: mpi_dbltyp
  ! working integer precision
  integer(int_mpi), protected, save :: mpi_inttyp
  !
  logical, private, save :: random_is_initialized = fals
  !
  ! Definition of 1 Work Unit using TauBench Benchmark Test
  ! Using Intel 2013 compilers
  real(wp), parameter :: one_work_unit = 6.9205_wp
  ! Westmere nodes
 !real(wp), parameter :: one_work_unit = 25.551_wp/3.0_wp
  ! Sandy Bridge nodes
 !real(wp), parameter :: one_work_unit = 22.808_wp/3.0_wp
  ! Ivy Bridge nodes
 !real(wp), parameter :: one_work_unit = 19.096_wp/3.0_wp
  ! Using Intel version 10 compilers
 !real(wp), parameter :: one_work_unit = 7.026_wp
  !
  ! Global parameter to help clarify calls to alloc_error
  !
  logical, parameter :: skip_alloc_pause = fals
  !
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  ! ################                                             ###############
  ! ################   GENERIC INTERFACES FOR GLOBAL FUNCTIONS   ###############
  ! ################                                             ###############
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  !
  ! Generic interface for reallocating arrays of varying types and ranks
  !
  private :: reallocate1_i4,reallocate1_i8
  private :: reallocate1_r4,reallocate1_r8
  private :: reallocate1_c
  private :: reallocate2_i4,reallocate2_i8
  !
  interface reallocate
    module procedure reallocate1_i4, reallocate1_i8, &
                     reallocate1_r4, reallocate1_r8, &
                     reallocate2_i4, reallocate2_i8, &
                     reallocate2_r4, reallocate2_r8, &
                     reallocate1_c
  end interface reallocate
  !
  !
  !
  private :: cross_product_DP,cross_product_QP
  !
  interface cross_product
    module procedure cross_product_SP, cross_product_DP, cross_product_QP
  end interface cross_product
  !
  !
  !
  private :: deltafun_DP!,deltafun_QP
  !
  interface deltafun
    module procedure deltafun_DP!, deltafun_QP
  end interface deltafun
  !
  !
  !
  private :: selfdot_SP,selfdot_DP,selfdot_QP
  !
  interface selfdot
    module procedure selfdot_SP, selfdot_DP, selfdot_QP
  end interface selfdot
  !
  !
  !
  private :: rnorm2_SP,rnorm2_DP,rnorm2_QP
  !
  interface rnorm2
    module procedure rnorm2_SP, rnorm2_DP, rnorm2_QP
  end interface rnorm2
  !
  !
  !
  private :: unit_vector_SP,unit_vector_DP,unit_vector_QP
  !
  interface unit_vector
    module procedure unit_vector_SP, unit_vector_DP, unit_vector_QP
  end interface unit_vector
  !
  !
  !
  private :: outerprod_DP,outerprod_QP
  !
  interface outerprod
    module procedure outerprod_DP, outerprod_QP
  end interface outerprod
  !
  !
  !
  private :: reciprocal_SP, reciprocal_DP,reciprocal_QP
  !
  interface reciprocal
    module procedure reciprocal_SP, reciprocal_DP, reciprocal_QP
  end interface reciprocal
  !
  !
  !
  private :: chop_WP, chop_SP, chop_DP, chop_QP
  !
  interface chop
    module procedure chop_WP, chop_SP, chop_DP, chop_QP
  end interface chop
  !
  !
  !
  private :: reflect_DP,reflect_QP
  !
  interface reflect
    module procedure reflect_DP, reflect_QP
  end interface reflect
  !
  !
  !
  private :: rotate_2d_DP,rotate_2d_QP
  !
  interface rotate_2d
    module procedure rotate_2d_DP, rotate_2d_QP
  end interface rotate_2d
  !
  !
  !
  private :: rotate_DP,rotate_QP
  !
  interface rotate
    module procedure rotate_DP, rotate_QP
  end interface rotate
  !
  !
  !
  private :: node_center_SP,node_center_DP,node_center_QP
  !
  interface node_center
    module procedure node_center_SP, node_center_DP, node_center_QP
  end interface node_center
  !
  !
  !
  private :: face_normal_SP,face_normal_DP,face_normal_QP
  !
  interface face_normal
    module procedure face_normal_SP, face_normal_DP, face_normal_QP
  end interface face_normal
  !
  !
  !
  private :: compute_jacobian_SP,compute_jacobian_DP,compute_jacobian_QP
  !
  interface compute_jacobian
    module procedure compute_jacobian_SP, &
                     compute_jacobian_DP, &
                     compute_jacobian_QP
  end interface compute_jacobian
  !
  !
  !
  private :: reverse_metrics_SP,reverse_metrics_DP,reverse_metrics_QP
  !
  interface reverse_metrics
    module procedure reverse_metrics_SP, &
                     reverse_metrics_DP, &
                     reverse_metrics_QP
  end interface reverse_metrics
  !
  !
  !
  private :: invert_matrix_DP,invert_matrix_QP
  !
  interface invert_matrix
    module procedure invert_matrix_DP, invert_matrix_QP
  end interface invert_matrix
  !
  !
  !
  private :: invert_matrix90_DP,invert_matrix90_QP
  !
  interface invert_matrix90
    module procedure invert_matrix90_DP, invert_matrix90_QP
  end interface invert_matrix90
  !
  !
  !
  private :: ludcmp_DP,ludcmp_QP
  !
  interface ludcmp
    module procedure ludcmp_DP, ludcmp_QP
  end interface ludcmp
  !
  !
  !
  private :: lubksb_DP,lubksb_QP
  !
  interface lubksb
    module procedure lubksb_DP, lubksb_QP
  end interface lubksb
  !
  !
  !
  private :: determinant_DP,determinant_QP
  !
  interface determinant
    module procedure determinant_DP, determinant_QP
  end interface determinant
  !
  !
  !
  private :: random_i4,random_sp,random_dp,random_qp
  !
  interface random
    module procedure random_i4,random_sp,random_dp,random_qp
  end interface random
  !
  !
  !
  private :: findloc_i4,findloc_i8,findloc_r4,findloc_r8
  !
  interface findloc
    module procedure findloc_i4,findloc_i8,findloc_r4,findloc_r8
  end interface findloc
  !
  !
  !
  !
  private :: identity_matrix_dim_SP,identity_matrix_mat_SP
  private :: identity_matrix_dim_DP,identity_matrix_mat_DP
  private :: identity_matrix_dim_QP,identity_matrix_mat_QP
  !
  interface identity_matrix
    module procedure identity_matrix_dim_SP,identity_matrix_mat_SP, &
                     identity_matrix_dim_DP,identity_matrix_mat_DP, &
                     identity_matrix_dim_QP,identity_matrix_mat_QP
  end interface identity_matrix
  !
  !
  !
  private :: intseq_no_interval,intseq_with_interval
  !
  interface intseq
    module procedure intseq_no_interval,intseq_with_interval
  end interface intseq
  !
  !
  !
  private :: big_to_little_endian_i2
  private :: big_to_little_endian_i4
  private :: big_to_little_endian_i8
  private :: big_to_little_endian_r4
  private :: big_to_little_endian_r8
  !
  interface big_to_little_endian
    module procedure big_to_little_endian_i2, &
                     big_to_little_endian_i4, &
                     big_to_little_endian_i8, &
                     big_to_little_endian_r4, &
                     big_to_little_endian_r8
  end interface big_to_little_endian
  !
contains
!
!###############################################################################
!
subroutine initialize_runtime_parameters
  !
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  ! #####   I SHOULD PROBABLY MAKE SOME ASCII FIGURES HERE FOR EACH CELL   #####
  ! #####   TYPE TO HELP ILLUSTRATE THE DIFFERENT GEOMETRIES AND THE       #####
  ! #####   REQUIREMENTS FOR THE NODE ORDERINGS OF EACH GEOMETRY TYPE.     #####
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  !
  !.. Local Scalars ..
  integer :: this_geom,this_node,this_face,this_edge,i
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_runtime_parameters"
  !
continue
  !
  Geom_Name(Geom_Node) = "Geom_Node"
  Geom_Name(Geom_Edge) = "Geom_Edge"
  Geom_Name(Geom_Tria) = "Geom_Tria"
  Geom_Name(Geom_Quad) = "Geom_Quad"
  Geom_Name(Geom_Tetr) = "Geom_Tetr"
  Geom_Name(Geom_Pyra) = "Geom_Pyra"
  Geom_Name(Geom_Pris) = "Geom_Pris"
  Geom_Name(Geom_Unknown) = "Geom_Unkn"
  Geom_Name(Geom_Hexa) = "Geom_Hexa"
  !
  Geom_String(Geom_Node) = "Node (0-D)"
  Geom_String(Geom_Edge) = "Edge (1-D)"
  Geom_String(Geom_Tria) = "Triangle (2-D)"
  Geom_String(Geom_Quad) = "Quadrilateral (2-D)"
  Geom_String(Geom_Tetr) = "Tetrahedron (3-D)"
  Geom_String(Geom_Pyra) = "Pyramid (3-D)"
  Geom_String(Geom_Pris) = "Prism (3-D)"
  Geom_String(Geom_Unknown) = "Unknown Geometry"
  Geom_String(Geom_Hexa) = "Hexahedron (3-D)"
  !
  geom_dimen = 0
  geom_nodes = 0
  geom_edges = 0
  geom_faces = 0
  !
  num_face_nodes = 0
  num_node_edges = 0
  !
  edge_nodes = 0
  side_nodes = 0
  node_edges = 0
  node_nodes = 0
  !
  cl = one
  !
  geom_dimen(Geom_0D) = 0
  geom_dimen(Geom_1D) = 1
  geom_dimen(Geom_2D) = 2
  geom_dimen(Geom_3D) = 3
  !
  geom_nodes(Geom_Node) = 1
  geom_nodes(Geom_Edge) = 2
  geom_nodes(Geom_Tria) = 3
  geom_nodes(Geom_Quad) = 4
  geom_nodes(Geom_Tetr) = 4
  geom_nodes(Geom_Pyra) = 5
  geom_nodes(Geom_Pris) = 6
  geom_nodes(Geom_Hexa) = 8
  !
  geom_edges(Geom_Edge) = 1
  geom_edges(Geom_Tria) = 3
  geom_edges(Geom_Quad) = 4
  geom_edges(Geom_Tetr) = 6
  geom_edges(Geom_Pyra) = 8
  geom_edges(Geom_Pris) = 9
  geom_edges(Geom_Hexa) = 12
  !
  geom_faces(Geom_Edge) = 2
  geom_faces(Geom_Tria) = 3
  geom_faces(Geom_Quad) = 4
  geom_faces(Geom_Tetr) = 4
  geom_faces(Geom_Pyra) = 5
  geom_faces(Geom_Pris) = 5
  geom_faces(Geom_Hexa) = 6
  !
  ! Define the geometry type of each cell face
  !
  cell_face_geom(1:2,Geom_Edge) = Geom_Node
  cell_face_geom(1:3,Geom_Tria) = Geom_Edge
  cell_face_geom(1:4,Geom_Quad) = Geom_Edge
  cell_face_geom(1:4,Geom_Tetr) = Geom_Tria
  cell_face_geom(1:4,Geom_Pyra) = Geom_Tria
  cell_face_geom(5:5,Geom_Pyra) = Geom_Quad
  cell_face_geom(1:2,Geom_Pris) = Geom_Tria
  cell_face_geom(3:5,Geom_Pris) = Geom_Quad
  cell_face_geom(1:6,Geom_Hexa) = Geom_Quad
  !
  ! Store the number of nodes that define the faces of each geometry
  !
  do this_geom = Geom_Min,Geom_Max
    do this_face = 1,geom_faces(this_geom)
      num_face_nodes(this_face,this_geom) = &
                         geom_nodes(cell_face_geom(this_face,this_geom))
    end do
  end do
  !
  !
  ! Nodes defining the edges of each geometry
  !
  edge_nodes(1:2,1,Geom_Edge) = [0,1]
  !
  edge_nodes(1:2,1,Geom_Tria) = [0,1]
  edge_nodes(1:2,2,Geom_Tria) = [1,2]
  edge_nodes(1:2,3,Geom_Tria) = [2,0]
  !
  edge_nodes(1:2,1,Geom_Quad) = [0,1]
  edge_nodes(1:2,2,Geom_Quad) = [1,2]
  edge_nodes(1:2,3,Geom_Quad) = [2,3]
  edge_nodes(1:2,4,Geom_Quad) = [3,0]
  !
  edge_nodes(1:2,1,Geom_Tetr) = [0,1]
  edge_nodes(1:2,2,Geom_Tetr) = [1,2]
  edge_nodes(1:2,3,Geom_Tetr) = [2,0]
  edge_nodes(1:2,4,Geom_Tetr) = [0,3]
  edge_nodes(1:2,5,Geom_Tetr) = [1,3]
  edge_nodes(1:2,6,Geom_Tetr) = [2,3]
  !
  edge_nodes(1:2,1,Geom_Pyra) = [0,1]
  edge_nodes(1:2,2,Geom_Pyra) = [1,2]
  edge_nodes(1:2,3,Geom_Pyra) = [2,3]
  edge_nodes(1:2,4,Geom_Pyra) = [3,0]
  edge_nodes(1:2,5,Geom_Pyra) = [0,4]
  edge_nodes(1:2,6,Geom_Pyra) = [1,4]
  edge_nodes(1:2,7,Geom_Pyra) = [2,4]
  edge_nodes(1:2,8,Geom_Pyra) = [3,4]
  !
  edge_nodes(1:2,1,Geom_Pris) = [0,1]
  edge_nodes(1:2,2,Geom_Pris) = [1,2]
  edge_nodes(1:2,3,Geom_Pris) = [2,0]
  edge_nodes(1:2,4,Geom_Pris) = [0,3]
  edge_nodes(1:2,5,Geom_Pris) = [1,4]
  edge_nodes(1:2,6,Geom_Pris) = [2,5]
  edge_nodes(1:2,7,Geom_Pris) = [3,4]
  edge_nodes(1:2,8,Geom_Pris) = [4,5]
  edge_nodes(1:2,9,Geom_Pris) = [5,3]
  !
  edge_nodes(1:2, 1,Geom_Hexa) = [0,3]
  edge_nodes(1:2, 2,Geom_Hexa) = [3,2]
  edge_nodes(1:2, 3,Geom_Hexa) = [2,1]
  edge_nodes(1:2, 4,Geom_Hexa) = [1,0]
  edge_nodes(1:2, 5,Geom_Hexa) = [0,4]
  edge_nodes(1:2, 6,Geom_Hexa) = [1,5]
  edge_nodes(1:2, 7,Geom_Hexa) = [2,6]
  edge_nodes(1:2, 8,Geom_Hexa) = [3,7]
  edge_nodes(1:2, 9,Geom_Hexa) = [4,5]
  edge_nodes(1:2,10,Geom_Hexa) = [5,6]
  edge_nodes(1:2,11,Geom_Hexa) = [6,7]
  edge_nodes(1:2,12,Geom_Hexa) = [7,4]
  !
  ! Nodes defining the faces of each geometry
  !
  side_nodes(1:1,1,Geom_Edge) = Xmin1D
  side_nodes(1:1,2,Geom_Edge) = Xmax1D
  !
  side_nodes(1:2,1,Geom_Tria) = [0,1]
  side_nodes(1:2,2,Geom_Tria) = [1,2]
  side_nodes(1:2,3,Geom_Tria) = [2,0]
  !
 !integer, parameter :: Ymin2D_Face = 1
 !integer, parameter :: Xmax2D_Face = 2
 !integer, parameter :: Ymax2D_Face = 3
 !integer, parameter :: Xmin2D_Face = 4
 !integer, parameter :: Zmin2D_Face = 0
 !integer, parameter :: Zmax2D_Face = 0
  !
  side_nodes(1:2,Ymin2D_Face,Geom_Quad) = Ymin2D ! [0,1]
  side_nodes(1:2,Xmax2D_Face,Geom_Quad) = Xmax2D ! [1,2]
  side_nodes(1:2,Ymax2D_Face,Geom_Quad) = Ymax2D ! [2,3]
  side_nodes(1:2,Xmin2D_Face,Geom_Quad) = Xmin2D ! [3,0]
  !
  side_nodes(1:3,1,Geom_Tetr) = [1,2,3]
  side_nodes(1:3,2,Geom_Tetr) = [2,0,3]
  side_nodes(1:3,3,Geom_Tetr) = [3,0,1]
  side_nodes(1:3,4,Geom_Tetr) = [0,2,1]
  !
  side_nodes(1:3,1,Geom_Pyra) = [1,2,4]
  side_nodes(1:3,2,Geom_Pyra) = [2,3,4]
  side_nodes(1:3,3,Geom_Pyra) = [3,0,4]
  side_nodes(1:3,4,Geom_Pyra) = [0,1,4]
  side_nodes(1:4,5,Geom_Pyra) = [0,3,2,1]
  !
  side_nodes(1:3,1,Geom_Pris) = [0,2,1]
  side_nodes(1:3,2,Geom_Pris) = [3,4,5]
  side_nodes(1:4,3,Geom_Pris) = [0,1,4,3]
  side_nodes(1:4,4,Geom_Pris) = [1,2,5,4]
  side_nodes(1:4,5,Geom_Pris) = [3,5,2,0]
  !
 !integer, parameter :: Zmin3D_Face = 1
 !integer, parameter :: Zmax3D_Face = 2
 !integer, parameter :: Xmax3D_Face = 3
 !integer, parameter :: Xmin3D_Face = 4
 !integer, parameter :: Ymin3D_Face = 5
 !integer, parameter :: Ymax3D_Face = 6
  !
  side_nodes(1:4,Zmin3D_Face,Geom_Hexa) = Zmin3D ! [0,3,2,1]
  side_nodes(1:4,Zmax3D_Face,Geom_Hexa) = Zmax3D ! [4,5,6,7]
  side_nodes(1:4,Xmax3D_Face,Geom_Hexa) = Xmax3D ! [1,2,6,5]
  side_nodes(1:4,Xmin3D_Face,Geom_Hexa) = Xmin3D ! [0,4,7,3]
  side_nodes(1:4,Ymin3D_Face,Geom_Hexa) = Ymin3D ! [0,1,5,4]
  side_nodes(1:4,Ymax3D_Face,Geom_Hexa) = Ymax3D ! [3,7,6,2]
  !
  ! Create the connectivity information between the nodes and edges
  ! of each geometry type
  !
  do this_geom = Geom_Min,Geom_Max
    do this_node = 1,geom_nodes(this_geom)
      i = 0
      do this_edge = 1,geom_edges(this_geom)
        if (this_node-1 == edge_nodes(1,this_edge,this_geom)) then
          i = i + 1
          node_edges(i,this_node,this_geom) = +this_edge
          node_nodes(i,this_node,this_geom) = &
                                          edge_nodes(2,this_edge,this_geom)
        else if (this_node-1 == edge_nodes(2,this_edge,this_geom)) then
          i = i + 1
          node_edges(i,this_node,this_geom) = -this_edge
          node_nodes(i,this_node,this_geom) = &
                                          edge_nodes(1,this_edge,this_geom)
        end if
      end do
      num_node_edges(this_node,this_geom) = i
    end do
  end do
  !
  ! Set the value of cl for triangles depending on if the reference
  ! triangle is based on the unit or biunit triangle
  !
  if (.not. tria_is_biunit) cl(Geom_Tria) = two
  !
 !call stop_recon(stop_mpi,pname,error_message)
  !
end subroutine initialize_runtime_parameters
!
!###############################################################################
!
subroutine alloc_error(procedure_name,variable_name,ioptype,error_flag, &
                       iemessag,override_memory_pause)
  !
  ! Subroutine to generate an error message when
  ! an allocate or deallocate returns an error code
  !
  !.. Formal Arguments ..
  character(len=*),           intent(in) :: procedure_name
  character(len=*),           intent(in) :: variable_name
  integer,                    intent(in) :: ioptype
  integer,                    intent(in) :: error_flag
  character(len=*), optional, intent(in) :: iemessag
  logical,          optional, intent(in) :: override_memory_pause
  !
#ifdef TRACK_MEMORY_ALLOCATIONS
  !.. Local Scalars ..
  logical :: alloc_memory_pause
  character(len=50+len(procedure_name)+len(variable_name)) :: pause_message
#endif
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "alloc_error"
  !
continue
  !
  if (error_flag /= 0) then
    !
    write (iout,1)
    if (ioptype == 1) then
      write (iout,2) trim(adjustl(variable_name))
    else
      write (iout,3) trim(adjustl(variable_name))
    end if
    !
    write (iout,4) trim(adjustl(procedure_name)),error_flag
    if (present(iemessag)) write (iout,5) trim(adjustl(iemessag))
    write (iout,6)
    !
    call stop_recon(abort,pname)
    !
#ifdef TRACK_MEMORY_ALLOCATIONS
#define MEMORY_TEST_PAUSE
  else
    !
    alloc_memory_pause = true
    if (present(override_memory_pause)) then
      alloc_memory_pause = override_memory_pause
    end if
    !
    if (alloc_memory_pause) then
      write (pause_message,7) trim(merge("  ","de",ioptype == 1)), &
                              trim(adjustl(variable_name)), &
                              trim(adjustl(procedure_name))
      call memory_pause(pause_message)
    end if
    !
#endif
  end if
  !
  ! Format statements
  !
  1 format (/," ***************************************************** ")
  2 format (" UNABLE TO ALLOCATE MEMORY FOR VARIABLE ",a)
  3 format (" UNABLE TO DEALLOCATE MEMORY FOR VARIABLE ",a)
  4 format (" INSIDE PROCEDURE: ",a,/, &
            " ERROR INDICATOR = ",i0)
  5 format (" ERROR MESSAGE: ",a)
  6 format (" ***************************************************** ",//)
  7 format ("    After ",a,"allocating array '",a, &
            "' in subroutine '",a,"'")
  !
end subroutine alloc_error
!
!###############################################################################
!
subroutine memory_pause(message_1,message_2)
  !
  character(len=*),           intent(in) :: message_1
  character(len=*), optional, intent(in) :: message_2
  !
  character(len=1) :: keyval
  integer :: stop_value,ifmt
  character(len=300) :: msgfmt
  character(len=999) :: before_message
  character(len=999) :: after_message
  !
  character(len=*), parameter :: pname = "memory_pause"
  !
continue
  !
#ifdef MEMORY_TEST_PAUSE
  !
  if (present(message_2)) then
    after_message  = "    After: " // trim(adjustl(message_1))
    before_message = "   Before: " // trim(adjustl(message_2))
    write (msgfmt,1)
  else
    after_message  = trim(adjustl(message_1))
    before_message = ""
    write (msgfmt,2)
  end if
  !
  if (mypnum == 0) then
    !
    write (iout,msgfmt,advance='no') trim(adjustl(after_message)), &
                                     trim(adjustl(before_message))
    !
    flush (iout)
    read (*,*) keyval
    !
    stop_value = merge( 1 , 0 , uppercase(keyval)=="N" )
    !
  end if
  !
  if (ncpu > 1) then
    call mpi_bcast(stop_value,1,mpi_inttyp,0,mpi_comm_world,mpierr)
    call mpi_barrier(mpi_comm_world,mpierr)
  end if
  !
  if (stop_value == 1) then
    error_message = "Abort issused during memory test at: " // &
                    trim(adjustl(after_message))
    call stop_recon(stop_mpi,pname,error_message)
  end if
  !
  1 format ('(/," Memory test pause",/,a,/,a,/,', &
               '" Press ',"'N'",' to abort or any other key to continue:")')
  2 format ('(/," Memory test pause",/,4x,a,a,/,', &
               '" Press ',"'N'",' to abort or any other key to continue:")')
 !1 format (/," Memory test pause",/,a,/,a,/, &
 !            " Press 'N' to abort or any other key to continue:")
 !2 format (/," Memory test pause",/,4x,a,a,/, &
 !            " Press 'N' to abort or any other key to continue:")
  !
#endif
  !
end subroutine memory_pause
!
!###############################################################################
!
subroutine io_error(procedure_name,file_name,ioptype,iemessag,optmessag)
  !
  ! Subroutine to generate an error message when
  ! an open or close returns an error code
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: procedure_name
  character(len=*), intent(in) :: file_name
  integer,          intent(in) :: ioptype
  integer,          intent(in) :: iemessag
  !
  !.. Optional Arguments ..
  character(len=*), optional, intent(in) :: optmessag
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "io_error"
  !
continue
  !
  if (iemessag /= 0) then
    !
    write (iout,10)
    !
    ! Open/Close errors
    !
    if (ioptype == 1) then
      write (iout,1) trim(adjustl(file_name))
    else if (ioptype == 2) then
      write (iout,2) trim(adjustl(file_name))
    else if (ioptype == 3) then
      write (iout,3) trim(adjustl(file_name))
    else if (ioptype == 4) then
      write (iout,4) trim(adjustl(file_name))
      !
      ! Read/Write errors 
      ! Use the value of iemessag to determine the type of error
      !
    else if (ioptype == -1) then
      !
      if (iemessag > 0) then
        write (iout,5) trim(adjustl(file_name))
      else if (iemessag < 0) then
        write (iout,6) trim(adjustl(file_name))
      end if
      !
      ! Generic unknown error
      !
    else
      write (iout,99) trim(adjustl(file_name))
    end if
    !
    if (present(optmessag)) then
      write (iout,7) trim(adjustl(optmessag))
    end if
    !
    write (iout,11) trim(adjustl(procedure_name)),iemessag
    call stop_recon(abort,pname)
    !
  end if
  !
  ! Format Statements
  !
  1  format (" ERROR OPENING FILE '",a,"'")
  2  format (" ERROR CLOSING FILE '",a,"'")
  3  format (" ERROR REOPENING CGNS FILE '",a,"'")
  4  format (" UNABLE TO DETERMINE FORMAT OF PLOT3D FILE '",a,"'")
  5  format (" UNEXPECTED INPUT/OUTPUT ERROR WITH FILE '",a,"'")
  6  format (" UNEXPECTED END OF FILE/RECORD ERROR WITH FILE '",a,"'")
 !7  format (" ERROR OCCURED TRYING TO READ IN '",a,"'")
  7  format (" ERROR MESSAGE: '",a,"'")
  99 format (" UNKNOWN ERROR WITH FILE '",a,"'")
  10 format (/," ***************************************************** ")
  11 format (" INSIDE PROCEDURE: ",a,/, &
             " ERROR INDICATOR = ",i0,/, &
             " ***************************************************** ",//)
  !
end subroutine io_error
!
!###############################################################################
!
subroutine cgns_error(procedure_name,file_name,cgns_routine,iemessag, &
                      optmessag,optinteger)
  !
  !.. Include Statements ..
  include "cgnslib_f95.h"
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: procedure_name
  character(len=*), intent(in) :: file_name
  character(len=*), intent(in) :: cgns_routine
  integer,          intent(in) :: iemessag
  !
  !.. Optional Arguments ..
  character(len=*), optional, intent(in) :: optmessag
  integer,          optional, intent(in) :: optinteger
  !
  !.. Local Scalars ..
  integer :: ierror
  !
continue
  !
  if (iemessag == CG_OK) then
    ierror = 0
  else
    ierror = 1
  end if
  call mpi_allreduce(MPI_IN_PLACE,ierror,1,mpi_inttyp, &
                     MPI_MAX,mpi_comm_world,mpierr)
  if (ierror == 0) return
  !
  write (iout,10)
  !
  if (present(optmessag) .and. present(optinteger)) then
    write (iout,1) trim(adjustl(optmessag)),optinteger
  else if (present(optmessag)) then
    write (iout,2) trim(adjustl(optmessag))
  else if (present(optinteger)) then
    write (iout,3) optinteger
  end if
  !
  write (iout,11) trim(adjustl(procedure_name)),iemessag, &
                  trim(adjustl(cgns_routine))
  call cgp_error_exit_f
  !
  ! Format Statements
  !
  1  format (" ",a," = ",i0)
  2  format (" ",a)
  3  format (" ",i0)
  10 format (/," ***************************************************** ")
  11 format (" INSIDE PROCEDURE: ",a,/, &
             " ERROR INDICATOR = ",i0,/, &
             " ERROR RETURNED FROM CGNS PROCEDURE: ",a,/, &
             " ***************************************************** ",//)
  !
end subroutine cgns_error
!
!###############################################################################
!
subroutine initialize_mpi_environment
  !
  !.. Local Scalars ..
  integer :: r,p
  !
continue
  !
  ! Initialize the overall MPI session
  !
  call mpi_init(mpierr)
  call mpi_comm_size(mpi_comm_world,ncpu,mpierr)
  call mpi_comm_rank(mpi_comm_world,mypnum,mpierr)
  !
  ! MPI datatype for working integer precision
  !
 !r = range( 0 )
 !write (*,*) "integer precision"
 !write (*,*) r
 !write (*,*) selected_int_kind(r)
 !call mpi_type_create_f90_integer(r,mpi_inttyp,mpierr)
 !write (*,*) mpi_inttyp
  mpi_inttyp = merge(mpi_integer8,mpi_integer,wip==i8)
  !
  ! MPI datatype for half working real precision
  !
 !r = range( 0.0_sp )
 !p = precision( 0.0_sp )
 !write (*,*) "half working precision"
 !write (*,*) r,p
 !write (*,*) sp
 !call mpi_type_create_f90_real(p,r,mpi_hlftyp,mpierr)
 !write (*,*) mpi_hlftyp
  mpi_hlftyp = merge(mpi_real,mpi_double_precision,sp<r8)
  mpi_flttyp = merge(mpi_real,mpi_double_precision,wp<r8)
  !
  ! MPI datatype for working real precision
  !
 !r = range( 0.0_wp )
 !p = precision( 0.0_wp )
 !write (*,*) "working precision"
 !write (*,*) r,p
 !write (*,*) wp
 !call mpi_type_create_f90_real(p,r,mpi_flttyp,mpierr)
 !write (*,*) mpi_flttyp
 !mpi_flttyp = mpi_double_precision
  !
  ! MPI datatype for double working real precision
  ! NOTE: mpi_type_create_f90_real doesnt seem to support quad precision
  !
 !r = range( 0.0_qp )
 !p = precision ( 0.0_qp )
 !write (*,*) "double working precision"
 !write (*,*) r,p
 !write (*,*) qp
 !call mpi_type_create_f90_real(p,r,mpi_dbltyp,mpierr)
 !write (*,*) mpi_dbltyp
 !!
 !! Test to make sure mpi_in_place is Fortran compatible
 !!
 !if (is_in_place(mpi_in_place) /= 1) then
 !  write (*,*) "MPI_IN_PLACE does not seem to be compatible with Fortran!"
 !  call stop_recon(9,"initialize_mpi_environment")
 !end if
  !
end subroutine initialize_mpi_environment
!
!###############################################################################
!
subroutine finalize_mpi_environment
  !
continue
  !
  ! Finalize the overall MPI session
  !
  call mpi_finalize(mpierr)
  !
  ! Stop execution
  !
  flush(iout)
  !
  stop
  !
end subroutine finalize_mpi_environment
!
!###############################################################################
!
subroutine stop_recon(iend,procedure_name,optional_message)
  !
  !.. Formal Arguments ..
  integer, optional, intent(in) :: iend
  character(len=*), optional, intent(in) :: procedure_name
  character(len=*), optional, intent(in) :: optional_message
  !
  !.. Local Scalars ..
  integer(int_mpi) :: iflag
  character(len=10000) :: strng
  character(len=10000) :: output_message
  !
continue
  !
  iflag = abort
  if (present(iend)) iflag = int( iend , kind=kind(iflag) )
  !
  ! Reformat optional_message
  !
  if (present(optional_message)) then
    !
    if (scan(optional_message,new_line('a')) == 0) then
      !
      ! optional_message does NOT contain any new line characters so reformat
      ! the string so that it is multiple lines with a max of 50 characters
      ! per line.
      !
      output_message = reformat_string( optional_message )
      !
    else
      !
      ! optional_message already contains new line characters so just copy
      ! it straight into output_message
      !
      output_message = trim(adjustl(optional_message))
      !
    end if
    !
  end if
  !
  ! Write out the final information before stopping execution
  !
  if (mypnum == 0) write (iout,1)
  !
  if (ncpu > 0) then
    !
    ! First finalize the overall MPI session before stopping execution
    !
    if (ncpu == 1 .or. iflag == stop_mpi) then
      !
      if (mypnum == 0) then
        write (iout,2)
        if (present(procedure_name)) &
          write (iout,4) trim(adjustl(procedure_name))
        if (present(optional_message)) &
          write (iout,5) trim(adjustl(output_message))
        write (iout,6)
      end if
      !
      ! Stop all MPI processes
      !
      call mpi_finalize(mpierr)
      !
    else if (iflag == abort) then
      !
      ! Abort all MPI processes
      !
      write (iout,3) mypnum
      if (present(procedure_name)) &
        write (iout,4) trim(adjustl(procedure_name))
      if (present(optional_message)) &
        write (iout,5) trim(adjustl(output_message))
      write (iout,6)
      !
      call mpi_abort(mpi_comm_world,iflag,mpierr)
      !
    end if
    !
  end if
  !
  flush (iout)
  !
  stop
  !
  ! Format Statements
  !
  1 format (//,66("#"))
  2 format ("    CALLING PREMATURE MPI_FINALIZE ON ALL PROCESSORS! ")
  3 format ("           MPI ABORT ISSUED FROM PROCESSOR ",i0,"!")
  4 format (" STOP ISSUED FROM PROCEDURE '",a,"'")
  5 format (" ERROR MESSAGE: '",a,"'")
  6 format (66("#"),//)
  !
end subroutine stop_recon
!
!###############################################################################
!
! Include statements for functions that have a generic interface
!
include "Functions/reallocate.f90"
include "Functions/cross_product.f90"
include "Functions/deltafun.f90"
include "Functions/norms.f90"
include "Functions/outerprod.f90"
include "Functions/reciprocal.f90"
include "Functions/reflect.f90"
include "Functions/rotate.f90"
include "Functions/rotate_2d.f90"
include "Functions/matrix_ops.f90"
include "Functions/random.f90"
include "Functions/findloc.f90"
include "Functions/node_center.f90"
include "Functions/face_normal.f90"
include "Functions/jacobians.f90"
include "Functions/reverse_metrics.f90"
include "Functions/chop.f90"
include "Functions/intseq.f90"
include "Functions/big_to_little_endian.f90"
!
! Include statements for functions that DO NOT have a generic interface
!
include "Functions/uppercase.f90"
include "Functions/bc_integer_to_string.f90"
include "Functions/all_are_equal.f90"
include "Functions/reformat_string.f90"
include "Functions/ghost_nodes.f90"
!
!###############################################################################
!
end module module_kind_types
!
!###############################################################################
!
module module_memory
  !
  use module_kind_types, only : wp,zero,one
 !use module_kind_types, only : inttype => i8
  use module_kind_types, only : inttype => i4
  !
  implicit none
  !
  private
  !
  public :: memory
  public :: inttype
  !
  ! Derived type for working with memory information
  !
  type :: memory
    character(len=6) :: units = " Bytes"
    integer(inttype) :: mem = 0_inttype
    integer(inttype) :: siz = 0_inttype
    real(wp)         :: rmem = zero
  contains
    procedure :: get_units
    procedure :: memout
    procedure :: reset => reset_mem
    procedure, private :: reset_and_compute_memory_usage
    procedure, private :: reset_and_add_memory
    procedure, private :: compute_memory_usage
    procedure, private :: add_memory
    generic :: add => compute_memory_usage, add_memory
    generic :: set => reset_and_compute_memory_usage, &
                      reset_and_add_memory
    generic :: write(formatted) => memout
  end type memory
  !
  !.. Module Parameters ..
  real(wp), parameter :: prefix_multiple = 1024.0_wp
  !
contains
!
!###############################################################################
!
pure function get_units(this) result(return_value)
  !
  !.. Formal Arguments ..
  class(memory), intent(in) :: this
  !
  !.. Function Result ..
  type(memory) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: rmem
  character(len=6) :: c
  !
continue
  !
  rmem = this%rmem
  c = this%units
  !
  select case ( c(1:1) )
    !
    case ("M")
      !
      if (rmem/prefix_multiple > one) then
        rmem = rmem/prefix_multiple
        c(1:1) = "G"
      end if
      !
    case ("k")
      !
      if (rmem/prefix_multiple > one) then
        rmem = rmem/prefix_multiple
        c(1:1) = "M"
      end if
      !
      if (rmem/prefix_multiple > one) then
        rmem = rmem/prefix_multiple
        c(1:1) = "G"
      end if
      !
    case default
      !
      if (rmem/prefix_multiple > one) then
        rmem = rmem/prefix_multiple
        c(1:1) = "k"
      end if
      !
      if (rmem/prefix_multiple > one) then
        rmem = rmem/prefix_multiple
        c(1:1) = "M"
      end if
      !
      if (rmem/prefix_multiple > one) then
        rmem = rmem/prefix_multiple
        c(1:1) = "G"
      end if
      !
  end select
  !
  return_value%units = c
  return_value%mem = nint(rmem)
  return_value%siz = this%siz
  return_value%rmem = rmem
  !
end function get_units
!
!###############################################################################
!
subroutine memout(dtv,iunit,iotype,vlist,my_iostat,my_iomsg)
  !
  !.. Formal Arguments ..
  class(memory), intent(in) :: dtv
  integer, intent(in) :: iunit
  character(len=*), intent(in) :: iotype
  integer, intent(in) :: vlist(:)
  integer, intent(out) :: my_iostat
  character(len=*), intent(inout) :: my_iomsg
  !
continue
  !
  write (iunit,1,iostat=my_iostat) dtv%siz, dtv%rmem, trim(dtv%units)
  !
  1 format (" - size = ",i11,", memory = ",f8.2,1x,a)
  !
end subroutine memout
!
!###############################################################################
!
pure subroutine compute_memory_usage(this,obj_storage,obj_size)
  !
  class(memory), intent(inout) :: this
  integer(inttype), intent(in) :: obj_storage
  integer(inttype), optional, intent(in) :: obj_size
  !
  type(memory) :: obj
  !
continue
  !
  obj%siz = 1_inttype
  if (present(obj_size)) then
    obj%siz = obj_size
  end if
  !
  obj%mem = obj_storage / 8_inttype
  obj%rmem = real(obj%mem,kind=wp) * real(obj%siz,kind=wp)
  !
  call this%add(obj)
  !
end subroutine compute_memory_usage
!
!###############################################################################
!
pure subroutine reset_and_compute_memory_usage(this,obj_storage,obj_size)
  !
  class(memory), intent(inout) :: this
  integer(inttype), intent(in) :: obj_storage
  integer(inttype), optional, intent(in) :: obj_size
  !
  type(memory) :: obj
  !
continue
  !
  obj%siz = 1_inttype
  if (present(obj_size)) then
    obj%siz = obj_size
  end if
  !
  obj%mem = obj_storage / 8_inttype
  obj%rmem = real(obj%mem,kind=wp) * real(obj%siz,kind=wp)
  !
  call this%reset
  call this%add(obj)
  !
end subroutine reset_and_compute_memory_usage
!
!###############################################################################
!
pure subroutine add_memory(this,obj)
  !
  class(memory), intent(inout) :: this
  type(memory), intent(in) :: obj
  !
  type(memory) :: tmp_obj
  real(wp) :: rmem
  !
continue
  !
  select case ( this%units(1:1) )
    case ("k")
      this%rmem = this%rmem * prefix_multiple
    case ("M")
      this%rmem = this%rmem * (prefix_multiple**2)
    case ("G")
      this%rmem = this%rmem * (prefix_multiple**3)
  end select
  !
  rmem = obj%rmem
  select case ( obj%units(1:1) )
    case ("k")
      rmem = rmem * prefix_multiple
    case ("M")
      rmem = rmem * (prefix_multiple**2)
    case ("G")
      rmem = rmem * (prefix_multiple**3)
  end select
  !
  tmp_obj%rmem = this%rmem + rmem
  tmp_obj%siz  = this%siz  + obj%siz
  !
  tmp_obj = tmp_obj%get_units()
  !
  this%units = tmp_obj%units
  this%mem   = tmp_obj%mem
  this%siz   = tmp_obj%siz
  this%rmem  = tmp_obj%rmem
  !
end subroutine add_memory
!
!###############################################################################
!
pure subroutine reset_and_add_memory(this,obj)
  !
  class(memory), intent(inout) :: this
  type(memory), intent(in) :: obj
  !
  type(memory) :: tmp_obj
  real(wp) :: rmem
  !
continue
  !
  rmem = obj%rmem
  select case ( obj%units(1:1) )
    case ("k")
      rmem = rmem * prefix_multiple
    case ("M")
      rmem = rmem * prefix_multiple * prefix_multiple
    case ("G")
      rmem = rmem * prefix_multiple * prefix_multiple * prefix_multiple
  end select
  !
  tmp_obj%rmem = rmem
  tmp_obj%siz  = obj%siz
  !
  tmp_obj = tmp_obj%get_units()
  !
  call this%reset
  !
  this%units = tmp_obj%units
  this%mem   = tmp_obj%mem
  this%siz   = tmp_obj%siz
  this%rmem  = tmp_obj%rmem
  !
end subroutine reset_and_add_memory
!
!###############################################################################
!
elemental subroutine reset_mem(this)
  !
  !.. Formal Arguments ..
  class(memory), intent(inout) :: this
  !
continue
  !
  this%siz = 0_inttype
  this%mem = 0_inttype
  this%rmem = zero
  this%units = " Bytes"
  !
end subroutine reset_mem
!
!###############################################################################
!
end module module_memory
!
!###############################################################################
!
