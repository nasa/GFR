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
  real(qp), parameter :: qone    = 1.0_qp
  real(qp), parameter :: qtwo    = 2.0_qp
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
  real(wp), parameter :: eight   = 8.0_wp
  real(wp), parameter :: nine    = 9.0_wp
  real(wp), parameter :: ten     = 10.0_wp
  real(wp), parameter :: twelve  = 12.0_wp
  real(wp), parameter :: sixteen = 16.0_wp
  !
  real(wp), parameter :: one8   = one/eight
  real(wp), parameter :: one6   = one/six
  real(wp), parameter :: one4   = one/four
  real(wp), parameter :: one3   = one/three
  real(wp), parameter :: two3   = two/three
  real(wp), parameter :: three4 = three/four
  real(wp), parameter :: four3  = four/three
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
  !
  real(wp), parameter :: c0_03 = 0.03_wp
  real(wp), parameter :: c0_10 = 0.10_wp
  real(wp), parameter :: c0_15 = 0.15_wp
  real(wp), parameter :: c0_20 = 0.20_wp
  real(wp), parameter :: c0_30 = 0.30_wp
  real(wp), parameter :: c0_35 = 0.35_wp
  real(wp), parameter :: c0_45 = 0.45_wp
  real(wp), parameter :: c0_55 = 0.55_wp
  real(wp), parameter :: c0_65 = 0.65_wp
  real(wp), parameter :: c0_70 = 0.70_wp
  real(wp), parameter :: c0_80 = 0.80_wp
  real(wp), parameter :: c2_50 = 2.50_wp
  real(wp), parameter :: c3_50 = 3.50_wp
  real(wp), parameter :: c4_50 = 4.50_wp
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
  !
  ! Parameters defining workshop test cases
  !
  integer, parameter :: Generic_MeanFlow = 1
  integer, parameter :: Ringleb_Problem  = 2
  integer, parameter :: NACA0012_Airfoil = 3
  integer, parameter :: Laminar_BndLyr   = 4
  integer, parameter :: Shu_Vortex       = 5
  integer, parameter :: Vortex_Transport = 6
  integer, parameter :: Density_Transport = 7
  integer, parameter :: Entropy_Transport = 8
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
  integer, parameter :: project_derivatives = 1
  integer, parameter :: derivative_of_proj_fluxes = 2
  !
  integer, parameter :: visc_flux_BR1 = 1
  integer, parameter :: visc_flux_BR2 = 2
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
                                       Geom_Node,Geom_Edge,Geom_Tria,Geom_Quad, &
                                       Geom_Tetr,Geom_Pyra,Geom_Pris,Geom_Hexa)
  !
  integer, parameter :: Geom_Max = max(Geom_Unknown, &
                                       Geom_Node,Geom_Edge,Geom_Tria,Geom_Quad, &
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
  integer, protected, save, dimension(Geom_Min:Geom_Max)         :: geom_dimen
  integer, protected, save, dimension(Geom_Min:Geom_Max)         :: geom_nodes
  integer, protected, save, dimension(Geom_Min:Geom_Max)         :: geom_faces
  integer, protected, save, dimension(1:6,Geom_Min:Geom_Max)     :: num_face_nodes
  integer, protected, save, dimension(1:4,1:6,Geom_Min:Geom_Max) :: side_nodes
  !
  character(len=100), protected, save, dimension(Geom_Min:Geom_Max) :: Geom_String
  character(len=9), protected, save, dimension(Geom_Min:Geom_Max) :: Geom_Name
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
  ! CL : coefficient for the length of the face. This is needed because of
  !      the difference between how the reference triangle and reference
  !      quad are defined. 
  !      = 1.0 for triangles
  !      = 0.5 for quads
  !
  real(wp), protected, save :: cl(Geom_Min:Geom_Max)
  !
  ! MPI Scalars 
  !
  integer(int_mpi), protected, save :: mypnum
  integer(int_mpi), protected, save :: ncpu
  integer(int_mpi),            save :: mpierr
  !
  integer(int_mpi), parameter :: abort = int( 8 , kind=kind(abort) )
  integer(int_mpi), parameter :: stop_mpi = int( 9 , kind=kind(stop_mpi) )
  !
  ! MPI datatypes to match Fortran intrinsic datatypes
  !
  integer(int_mpi), protected, save :: mpi_flttyp ! working real precision
  integer(int_mpi), protected, save :: mpi_hlftyp ! half working real precision
  integer(int_mpi), protected, save :: mpi_dbltyp ! double working real precision
  integer(int_mpi), protected, save :: mpi_inttyp ! working integer precision
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
  geom_faces = 0
  num_face_nodes = 0
  side_nodes = 0
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
  geom_faces(Geom_Edge) = 2
  geom_faces(Geom_Tria) = 3
  geom_faces(Geom_Quad) = 4
  geom_faces(Geom_Tetr) = 4
  geom_faces(Geom_Pyra) = 5
  geom_faces(Geom_Pris) = 5
  geom_faces(Geom_Hexa) = 6
  !
  num_face_nodes(1:2,Geom_Edge) = geom_nodes(Geom_Node)
  num_face_nodes(1:3,Geom_Tria) = geom_nodes(Geom_Edge)
  num_face_nodes(1:4,Geom_Quad) = geom_nodes(Geom_Edge)
  num_face_nodes(1:4,Geom_Tetr) = geom_nodes(Geom_Tria)
  num_face_nodes(1:4,Geom_Pyra) = geom_nodes(Geom_Tria)
  num_face_nodes(5:5,Geom_Pyra) = geom_nodes(Geom_Quad)
  num_face_nodes(1:2,Geom_Pris) = geom_nodes(Geom_Tria)
  num_face_nodes(3:5,Geom_Pris) = geom_nodes(Geom_Quad)
  num_face_nodes(1:6,Geom_Hexa) = geom_nodes(Geom_Quad)
  !
  side_nodes(1:1,1,Geom_Edge) = 0
  side_nodes(1:1,2,Geom_Edge) = 1
  !
  side_nodes(1:2,1,Geom_Tria) = [0,1]
  side_nodes(1:2,2,Geom_Tria) = [1,2]
  side_nodes(1:2,3,Geom_Tria) = [2,0]
  !
  side_nodes(1:2,1,Geom_Quad) = [0,1]
  side_nodes(1:2,2,Geom_Quad) = [1,2]
  side_nodes(1:2,3,Geom_Quad) = [2,3]
  side_nodes(1:2,4,Geom_Quad) = [3,0]
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
  side_nodes(1:4,1,Geom_Hexa) = [0,3,2,1]
  side_nodes(1:4,2,Geom_Hexa) = [4,5,6,7]
  side_nodes(1:4,3,Geom_Hexa) = [1,2,6,5]
  side_nodes(1:4,4,Geom_Hexa) = [4,7,3,0]
  side_nodes(1:4,5,Geom_Hexa) = [0,1,5,4]
  side_nodes(1:4,6,Geom_Hexa) = [2,3,7,6]
  !
  if (.not. tria_is_biunit) cl(Geom_Tria) = half
  !
end subroutine initialize_runtime_parameters
!
!###############################################################################
!
!!subroutine reallocate1_i4(array,new_size)
!!  integer, intent(in) :: new_size
!!  integer(i4), allocatable, dimension(:), intent(inout) :: array
!!  integer(i4), allocatable, dimension(:) :: temp_array
!!  integer :: ierr,num_copy
!!continue
!!  ! Allocate the temporary array to store the old contents of array
!!  allocate ( temp_array(1:new_size) , source=0_i4 , stat=ierr )
!!  call alloc_error("reallocate1_i4","temp_array",1,ierr)
!!  ! If array is allocated, copy the old contents of array into the temporary
!!  ! array. If decreasing the size of array, make sure to limit the amount of
!!  ! information copied over from array to the size of temp_array.
!!  if (allocated(array)) then
!!    num_copy = min( size(temp_array) , size(array) )
!!    temp_array(1:num_copy) = array(1:num_copy)
!!  end if
!!  ! move_alloc does:
!!  !     1. deallocates array if it is allocated
!!  !     2. allocates array with the same type, type parameters,
!!  !        array bounds, and value as temp_array
!!  !     3. deallocates temp_array
!!  call move_alloc( from=temp_array , to=array )
!!end subroutine reallocate1_i4
!!!
!!!###############################################################################
!!!
!!subroutine reallocate1_i8(array,new_size)
!!  integer, intent(in) :: new_size
!!  integer(i8), allocatable, dimension(:), intent(inout) :: array
!!  integer(i8), allocatable, dimension(:) :: temp_array
!!  integer :: ierr,num_copy
!!continue
!!  allocate ( temp_array(1:new_size) , source=0_i8 , stat=ierr )
!!  call alloc_error("reallocate1_i8","temp_array",1,ierr)
!!  if (allocated(array)) then
!!    num_copy = min( size(temp_array) , size(array) )
!!    temp_array(1:num_copy) = array(1:num_copy)
!!  end if
!!  call move_alloc( from=temp_array , to=array )
!!end subroutine reallocate1_i8
!!!
!!!###############################################################################
!!!
!!subroutine reallocate1_r4(array,new_size)
!!  integer, intent(in) :: new_size
!!  real(r4), allocatable, dimension(:), intent(inout) :: array
!!  real(r4), allocatable, dimension(:) :: temp_array
!!  integer :: ierr,num_copy
!!continue
!!  allocate ( temp_array(1:new_size) , source=0.0_r4 , stat=ierr )
!!  call alloc_error("reallocate1_r4","temp_array",1,ierr)
!!  if (allocated(array)) then
!!    num_copy = min( size(temp_array) , size(array) )
!!    temp_array(1:num_copy) = array(1:num_copy)
!!  end if
!!  call move_alloc( from=temp_array , to=array )
!!end subroutine reallocate1_r4
!!!
!!!###############################################################################
!!!
!!subroutine reallocate1_r8(array,new_size)
!!  integer, intent(in) :: new_size
!!  real(r8), allocatable, dimension(:), intent(inout) :: array
!!  real(r8), allocatable, dimension(:) :: temp_array
!!  integer :: ierr,num_copy
!!continue
!!  allocate ( temp_array(1:new_size) , source=0.0_r8 , stat=ierr )
!!  call alloc_error("reallocate1_r8","temp_array",1,ierr)
!!  if (allocated(array)) then
!!    num_copy = min( size(temp_array) , size(array) )
!!    temp_array(1:num_copy) = array(1:num_copy)
!!  end if
!!  call move_alloc( from=temp_array , to=array )
!!end subroutine reallocate1_r8
!!!
!!!###############################################################################
!!!
!!subroutine reallocate1_c(array,new_size)
!!  integer, intent(in) :: new_size
!!  character(len=*), allocatable, dimension(:), intent(inout) :: array
!!  character(len=len(array)), allocatable, dimension(:) :: temp_array
!!  integer :: ierr,num_copy
!!continue
!! !write (iout,1) size(array),new_size
!!  allocate ( temp_array(1:new_size) , stat=ierr )
!!  call alloc_error("reallocate1_c","temp_array",1,ierr)
!!  if (allocated(array)) then
!!    num_copy = min( size(temp_array) , size(array) )
!!    temp_array(1:num_copy) = array(1:num_copy)
!!  end if
!!  call move_alloc( from=temp_array , to=array )
!!  1 format ("size(array) = ",i0,"; new_size = ",i0)
!!end subroutine reallocate1_c
!!!
!!!###############################################################################
!!!
!!subroutine reallocate2_i4(array,new_size1,new_size2)
!!  !
!!  !.. Formal Arguments ..
!!  integer, optional,                       intent(in) :: new_size1
!!  integer, optional,                       intent(in) :: new_size2
!!  integer(i4), allocatable, dimension(:,:), intent(inout) :: array
!!  !
!!  !.. Local Scalars ..
!!  integer :: temp_size1,temp_size2,num_copy1,num_copy2,ierr
!!  !
!!  !.. Local Allocatable Arrays ..
!!  integer(i4), allocatable, dimension(:,:) :: temp_array
!!  !
!!continue
!!  !
!!  ! Get the old dimensions of array
!!  !
!!  temp_size1 = 1
!!  temp_size2 = 1
!!  if (allocated(array)) then
!!    temp_size1 = size(array,1)
!!    temp_size2 = size(array,2)
!!  end if
!!  !
!!  if (present(new_size1)) temp_size1 = new_size1
!!  if (present(new_size2)) temp_size2 = new_size2
!!  !
!!  ! Allocate the temporary array to store the old contents of array
!!  !
!!  allocate ( temp_array(1:temp_size1,1:temp_size2) , source=0_i4 , stat=ierr )
!!  call alloc_error("reallocate2_i4","temp_array",1,ierr)
!!  !
!!  ! Copy the old contents of array into the temporary array
!!  !
!!  if (allocated(array)) then
!!    num_copy1 = min( size(temp_array,dim=1) , size(array,dim=1) )
!!    num_copy2 = min( size(temp_array,dim=2) , size(array,dim=2) )
!!    temp_array(1:num_copy1,1:num_copy2) = array(1:num_copy1,1:num_copy2)
!!  end if
!!  !
!!  ! move_alloc does:
!!  !     1. deallocates array
!!  !     2. allocates array with the same type, type parameters,
!!  !        array bounds, and value as temp_array
!!  !     3. deallocates temp_array
!!  !
!!  call move_alloc( from=temp_array , to=array )
!!  !
!!end subroutine reallocate2_i4
!!!
!!!###############################################################################
!!!
!!subroutine reallocate2_i8(array,new_size1,new_size2)
!!  !
!!  !.. Formal Arguments ..
!!  integer, optional,                       intent(in) :: new_size1
!!  integer, optional,                       intent(in) :: new_size2
!!  integer(i8), allocatable, dimension(:,:), intent(inout) :: array
!!  !
!!  !.. Local Scalars ..
!!  integer :: temp_size1,temp_size2,num_copy1,num_copy2,ierr
!!  !
!!  !.. Local Allocatable Arrays ..
!!  integer(i8), allocatable, dimension(:,:) :: temp_array
!!  !
!!continue
!!  !
!!  ! Get the old dimensions of array
!!  !
!!  temp_size1 = 1
!!  temp_size2 = 1
!!  if (allocated(array)) then
!!    temp_size1 = size(array,1)
!!    temp_size2 = size(array,2)
!!  end if
!!  !
!!  if (present(new_size1)) temp_size1 = new_size1
!!  if (present(new_size2)) temp_size2 = new_size2
!!  !
!!  ! Allocate the temporary array to store the old contents of array
!!  !
!!  allocate ( temp_array(1:temp_size1,1:temp_size2) , source=0_i8 , stat=ierr )
!!  call alloc_error("reallocate2_i8","temp_array",1,ierr)
!!  !
!!  ! Copy the old contents of array into the temporary array
!!  !
!!  if (allocated(array)) then
!!    num_copy1 = min( size(temp_array,dim=1) , size(array,dim=1) )
!!    num_copy2 = min( size(temp_array,dim=2) , size(array,dim=2) )
!!    temp_array(1:num_copy1,1:num_copy2) = array(1:num_copy1,1:num_copy2)
!!  end if
!!  !
!!  ! move_alloc does:
!!  !     1. deallocates array
!!  !     2. allocates array with the same type, type parameters,
!!  !        array bounds, and value as temp_array
!!  !     3. deallocates temp_array
!!  !
!!  call move_alloc( from=temp_array , to=array )
!!  !
!!end subroutine reallocate2_i8
!!!
!!!###############################################################################
!!!
!!subroutine reallocate2_r4(array,new_size1,new_size2)
!!  !
!!  !.. Formal Arguments ..
!!  integer,  optional,                       intent(in) :: new_size1
!!  integer,  optional,                       intent(in) :: new_size2
!!  real(r4), allocatable, dimension(:,:), intent(inout) :: array
!!  !
!!  !.. Local Scalars ..
!!  integer :: temp_size1,temp_size2,num_copy1,num_copy2,ierr
!!  !
!!  !.. Local Allocatable Arrays ..
!!  real(r4), allocatable, dimension(:,:) :: temp_array
!!  !
!!continue
!!  !
!!  ! Get the old dimensions of array
!!  !
!!  temp_size1 = 1
!!  temp_size2 = 1
!!  if (allocated(array)) then
!!    temp_size1 = size(array,1)
!!    temp_size2 = size(array,2)
!!  end if
!!  !
!!  if (present(new_size1)) temp_size1 = new_size1
!!  if (present(new_size2)) temp_size2 = new_size2
!!  !
!!  ! Allocate the temporary array to store the old contents of array
!!  !
!!  allocate ( temp_array(1:temp_size1,1:temp_size2) , source=0.0_r4 , stat=ierr )
!!  call alloc_error("reallocate2_r4","temp_array",1,ierr)
!!  !
!!  ! Copy the old contents of array into the temporary array
!!  !
!!  if (allocated(array)) then
!!    num_copy1 = min( size(temp_array,dim=1) , size(array,dim=1) )
!!    num_copy2 = min( size(temp_array,dim=2) , size(array,dim=2) )
!!    temp_array(1:num_copy1,1:num_copy2) = array(1:num_copy1,1:num_copy2)
!!  end if
!!  !
!!  ! move_alloc does:
!!  !     1. deallocates array
!!  !     2. allocates array with the same type, type parameters,
!!  !        array bounds, and value as temp_array
!!  !     3. deallocates temp_array
!!  !
!!  call move_alloc( from=temp_array , to=array )
!!  !
!!end subroutine reallocate2_r4
!!!
!!!###############################################################################
!!!
!!subroutine reallocate2_r8(array,new_size1,new_size2)
!!  !
!!  !.. Formal Arguments ..
!!  integer,  optional,                       intent(in) :: new_size1
!!  integer,  optional,                       intent(in) :: new_size2
!!  real(r8), allocatable, dimension(:,:), intent(inout) :: array
!!  !
!!  !.. Local Scalars ..
!!  integer :: temp_size1,temp_size2,num_copy1,num_copy2,ierr
!!  !
!!  !.. Local Allocatable Arrays ..
!!  real(r8), allocatable, dimension(:,:) :: temp_array
!!  !
!!continue
!!  !
!!  ! Get the old dimensions of array
!!  !
!!  temp_size1 = 1
!!  temp_size2 = 1
!!  if (allocated(array)) then
!!    temp_size1 = size(array,1)
!!    temp_size2 = size(array,2)
!!  end if
!!  !
!!  if (present(new_size1)) temp_size1 = new_size1
!!  if (present(new_size2)) temp_size2 = new_size2
!!  !
!!  ! Allocate the temporary array to store the old contents of array
!!  !
!!  allocate ( temp_array(1:temp_size1,1:temp_size2) , source=0.0_r8 , stat=ierr )
!!  call alloc_error("reallocate2_r8","temp_array",1,ierr)
!!  !
!!  ! Copy the old contents of array into the temporary array
!!  !
!!  if (allocated(array)) then
!!    num_copy1 = min( size(temp_array,dim=1) , size(array,dim=1) )
!!    num_copy2 = min( size(temp_array,dim=2) , size(array,dim=2) )
!!    temp_array(1:num_copy1,1:num_copy2) = array(1:num_copy1,1:num_copy2)
!!  end if
!!  !
!!  ! move_alloc does:
!!  !     1. deallocates array
!!  !     2. allocates array with the same type, type parameters,
!!  !        array bounds, and value as temp_array
!!  !     3. deallocates temp_array
!!  !
!!  call move_alloc( from=temp_array , to=array )
!!  !
!!end subroutine reallocate2_r8
!!!
!!!###############################################################################
!!!
subroutine alloc_error(procedure_name,variable_name,ioptype,error_flag,iemessag)
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
    call stop_recon(abort,"alloc_error")
    !
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
  !
end subroutine alloc_error
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
      if (present(optmessag)) then
        write (iout,7) trim(adjustl(optmessag))
      end if
      !
      ! Generic unknown error
      !
    else
      write (iout,99) trim(adjustl(file_name))
    end if
    !
    write (iout,11) trim(adjustl(procedure_name)),iemessag
    call stop_recon(abort,"io_error")
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
  7  format (" ERROR OCCURED TRYING TO READ IN '",a,"'")
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
continue
  !
  if (iemessag == CG_OK) return
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
  1 format (//,"##################################################################")
  2 format (   "    CALLING PREMATURE MPI_FINALIZE ON ALL PROCESSORS! ")
  3 format (   "           MPI ABORT ISSUED FROM PROCESSOR ",i0,"!")
  4 format (   " STOP ISSUED FROM PROCEDURE '",a,"'")
  5 format (   " ERROR MESSAGE: '",a,"'")
  6 format (   "##################################################################",//)
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
include "Functions/chop.f90"
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
  implicit none
  !
  private
  !
  ! Derived type for working with memory information
  !
  type, public :: memory
    character(len=6) :: units
    integer :: mem
    integer :: siz
  contains
    procedure :: get_units
    procedure :: memout
    generic :: write(formatted) => memout
  end type memory
  !
contains
!
!###############################################################################
!
pure function get_units(a) result(return_value)
  !
  !.. Formal Arguments ..
  class(memory), intent(in) :: a
  !
  !.. Function Result ..
  type(memory) :: return_value
  !
  !.. Local Scalars ..
  integer :: i
  character(len=6) :: c
  !
continue
  !
  i = a%mem / 8
  c = " Bytes"
  !
  if (i/1024 > 0) then
    i = i/1024
    c(1:1) = "k"
  end if
  !
  if (i/1024 > 0) then
    i = i/1024
    c(1:1) = "M"
  end if
  !
  if (i/1024 > 0) then
    i = i/1024
    c(1:1) = "G"
  end if
  !
  return_value%mem = i
  return_value%units = c
  return_value%siz = a%siz
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
  write (iunit,1,iostat=my_iostat) dtv%siz, dtv%mem, trim(dtv%units)
  !
  1 format (" - size = ",i8,", memory = ",i4,1x,a)
  !
end subroutine memout
!
!###############################################################################
!
end module module_memory
!
!###############################################################################
!
