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
module compatibility_pgi_mod
  !
  ! Interfaces to PGI built-in (non-standard) intrinsic procedures
  !
  ! NOTE: This is copied from the file 'lib3f.h' which can be found
  !       in the include directory of the PGI installation.
  !
  implicit none
  !
#ifdef SPECIAL_FOR_PGI
  public
  !
  interface
    !
    function system(command) result(ierr)
      character*(*), intent(in) :: command
      integer :: ierr
    end function system
    !
    function ierrno() result(ierr)
      integer :: ierr
    end function ierrno
    !
    subroutine gerror(string)
      character*(*), intent(out) :: string
    end subroutine gerror
    !
  end interface
  !
contains
!
!###############################################################################
!
subroutine execute_command_line(command,wait,exitstat,cmdstat,cmdmsg)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: command
  !
  !.. Optional Arguments ..
  logical,          optional,    intent(in) :: wait
  integer,          optional, intent(inout) :: exitstat
  integer,          optional, intent(inout) :: cmdstat
  character(len=*), optional, intent(inout) :: cmdmsg
  !
  !.. Local Scalars ..
  integer :: ierr
  character(len=len_trim(adjustl(command))) :: command_string
  !
continue
  !
  if (present(exitstat)) exitstat = 0
  if (present(cmdstat)) cmdstat = 0
  if (present(cmdmsg)) cmdmsg = ""
  !
  command_string = trim(adjustl(command))
  !
  ierr = system( command_string )
  !
  if (ierr /= 0) then
    !
    if (present(exitstat)) then
      exitstat = ierr
    end if
    !
    if (present(cmdstat)) then
      cmdstat = ierrno()
    end if
    !
    if (present(cmdmsg)) then
      call gerror( cmdmsg )
    end if
    !
  end if
  !
end subroutine execute_command_line
!
!###############################################################################
!
#endif
end module compatibility_pgi_mod
!
!###############################################################################
!
module compatibility_gcc_mod
  !
  ! Interfaces to GCC built-in (non-standard) intrinsic procedures
  !
  ! NOTE: These interfaces are defined based on the definitions of these
  !       procedures in the GNU Fortran Manual which can be found at
  !       https://gcc.gnu.org/onlinedocs/gcc-4.9.4/gfortran.pdf
  !
  ! NOTE: The C binding names are based on the result of the command
  !       `nm -D libgfortran.so | grep -i '\(system\|ierrno\|gerror\)'`
  !       which shows that these procedures are prefixed with '_gfortran'
  !       when compiled and linked into the libraries.
  !
  implicit none
  !
#ifdef SPECIAL_FOR_GCC
  private
  !
  public :: system
  public :: ierrno
  public :: gerror
  !
  interface
    !
    function gcc_system(command) result(ierr) bind(c,name="_gfortran_system")
      use, intrinsic :: iso_c_binding, only : c_char,c_int
      character(kind=c_char), intent(in) :: command(*)
      integer(c_int) :: ierr
    end function gcc_system
    !
    ! NOTE: There are two version of ierrno()
    !           _gfortran_ierrno_i4 which returns a 32-bit integer
    !           _gfortran_ierrno_i8 which returns a 64-bit integer
    !
    function ierrno() result(ierr) bind(c,name="_gfortran_ierrno_i4")
      use, intrinsic :: iso_c_binding, only : c_int
      integer(c_int) :: ierr
    end function ierrno
    !
    subroutine gcc_gerror(string) bind(c,name="_gfortran_gerror")
      use, intrinsic :: iso_c_binding, only : c_char
      character(kind=c_char), intent(out) :: string(*)
    end subroutine gcc_gerror
    !
  end interface
  !
contains
!
!###############################################################################
!
function system(command) result(ierr)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_char,c_int,c_null_char
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: command
  !
  !.. Function Result ..
  integer(c_int) :: ierr
  !
  !.. Local Scalars ..
  integer :: n,nlen
  character(len=len_trim(adjustl(command))) :: command_string
  character(len=100) :: frmt
  !
  !.. Local Arrays ..
  character(len=1,kind=c_char) :: c_command(1:len_trim(adjustl(command))+1)
  !
continue
  !
  command_string = trim(adjustl(command))
  write (*,100) "command_string",command_string
  100 format (a," = '",a,"'")
  !
  nlen = len(command_string)
  do n = 1,nlen
    c_command(n) = command_string(n:n)
  end do
  c_command(nlen+1) = c_null_char
  write (frmt,101) nlen+1
  101 format ('(a," = ',"'",'",',i0,'(a),"',"'",'")')
 !101 format (a," = '",<nlen+1>(a),"'")
  write (*,fmt=frmt) "c_command",(c_command(n),n=1,size(c_command))
  !
  ierr = gcc_system( c_command )
  !
end function system
!
!###############################################################################
!
subroutine gerror(string)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_char,c_int,c_null_char
  !
  !.. Formal Arguments ..
 !character(len=:), allocatable, intent(out) :: string
  character(len=*), intent(out) :: string
  !
  !.. Local Scalars ..
  integer :: n,nlen,ierr
  !
  !.. Local Arrays ..
  character(len=1,kind=c_char) :: c_string(1:301)
  !
continue
  !
  c_string(:) = c_null_char
  !
  call gcc_gerror( c_string )
  !
  do n = 1,size(c_string)
    if (c_string(n) == c_null_char) then
      nlen = n-1
      exit
    end if
  end do
  !
 !allocate ( character(len=nlen) :: string , stat=ierr )
  !
  do n = 1,min(nlen,len(string))
    string(n:n) = c_string(n)
  end do
  !
end subroutine gerror
!
!###############################################################################
!
#endif
end module compatibility_gcc_mod
!
!###############################################################################
!
module base_kind_types
  !
  !.. Use Statements ..
  use iso_fortran_env, only : int32
  use iso_fortran_env, only : int64
  use iso_fortran_env, only : real32
  use iso_fortran_env, only : real64
  use iso_fortran_env, only : real128
  use iso_fortran_env, only : output_unit
  use iso_fortran_env, only : logical_kinds
  !
  implicit none
  !
  private
  !
  ! KIND definitions for 32-bit or 64-bit integer variables (IEEE 754)
  !
  integer, public, parameter :: I4 = int32
  integer, public, parameter :: I8 = max( int32 , int64 )
  !
  ! KIND definitions for single or double precision real variables (IEEE 754)
  !
  integer, public, parameter :: R4 = real32
  integer, public, parameter :: R8 = real64
  integer, public, parameter :: R16 = max( real64 , real128 )
  !
  ! KIND definitions for logical variables
  !
#ifdef SPECIAL_FOR_PGI
  integer, public, parameter :: l1k = logical_kinds(1)
#else
#ifdef SPECIAL_FOR_GCC
  integer,         parameter :: lb_lk = lbound(logical_kinds,dim=1)
  integer,         parameter :: ub_lk = ubound(logical_kinds,dim=1)
  integer, public, parameter :: l1k = min(logical_kinds(min(lb_lk+0,ub_lk)), &
                                          logical_kinds(min(lb_lk+1,ub_lk)), &
                                          logical_kinds(min(lb_lk+2,ub_lk)), &
                                          logical_kinds(min(lb_lk+3,ub_lk)), &
                                          logical_kinds(min(lb_lk+4,ub_lk)))
#else
  integer, public, parameter :: l1k = minval(logical_kinds)
#endif
#endif
  integer, public, parameter :: ldk = kind(.true.)
  !
  ! Unit for I/O to standard output
  !
  integer, public, parameter :: iout = output_unit
  !
end module base_kind_types
!
!###############################################################################
!
module module_kind_types
  !
#include <mpi_defs.h>
  !
  !.. Use Statements ..
  use base_kind_types
  use _MPI_MODULE_
#ifdef USE_INTEL_MKL
  use blas95, only : dot
#else
  use rename_intrinsics, only : dot => dot_product
#endif
#ifndef SPECIAL_FOR_GCC
  use, intrinsic :: ieee_arithmetic, only : is_nan => ieee_is_nan
  use, intrinsic :: ieee_arithmetic, only : is_finite => ieee_is_finite
  use, intrinsic :: ieee_arithmetic, only : is_negative => ieee_is_negative
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
  ! Aliases for kind type parameters
  !
  integer, parameter :: LK = l1k
  integer, parameter :: SP = r4
  integer, parameter :: DP = r8
  integer, parameter :: QP = r16
  integer, parameter :: WP = r8
  !
  integer, parameter :: WIP = i4
  integer, parameter :: INT_MPI = kind(MPI_INTEGER_KIND)
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
  logical(lk), parameter :: true = .true.
  logical(lk), parameter :: fals = .false.
  !
  ! Constants for cshift and eoshift functions
  !
  integer, parameter :: to_end   = -1
  integer, parameter :: to_front =  1
  !
  ! #########################################
  ! #########################################
  ! #####   GENERAL GLOBAL PARAMETERS   #####
  ! #########################################
  ! #########################################
  !
  ! I/O constants
  !
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
  integer, parameter :: bc_unsupported = -9998
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
  integer, parameter :: Generic_MeanFlow        = 1
  integer, parameter :: Channel_Flow            = 2
  integer, parameter :: Inviscid_Gaussian_Bump  = 3
  integer, parameter :: Laminar_BndLyr          = 4
  integer, parameter :: Shu_Vortex              = 5
  integer, parameter :: Vortex_Transport        = 6
  integer, parameter :: Density_Transport       = 7
  integer, parameter :: Entropy_Transport       = 8
  integer, parameter :: Taylor_Green_Vortex     = 9
  integer, parameter :: Nozzle_Jet              = 10
  integer, parameter :: Double_Mach_Reflection  = 11
  integer, parameter :: Infinite_Cylinder       = 12
  integer, parameter :: Freestream_Preservation = 13
  integer, parameter :: Advec_Problems(1:2) = [ Density_Transport, &
                                                Entropy_Transport ]
  integer, parameter :: Vortex_Problems(1:2) = [ Shu_Vortex, &
                                                 Vortex_Transport ]
  integer, parameter :: Transport_Problems(1:4) = [ Shu_Vortex, &
                                                    Vortex_Transport, &
                                                    Advec_Problems ]
  !
  ! Location of 1D solution points
  integer, parameter :: Equi_Distant_Points     = 0
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
  integer, parameter :: AlphaOptimized_TriPoints = 1
  integer, parameter :: BarycentricLobatto_TriPoints = 2
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
  integer, parameter :: invs_flux_ROE          = 1
  integer, parameter :: invs_flux_HLLC         = 2
  integer, parameter :: invs_flux_LDFSS        = 3
  integer, parameter :: invs_flux_Rotated_RHLL = 4
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
  integer, parameter :: Internally_Generated_Grid = -1
  integer, parameter :: Gambit_Neutral = 1
  integer, parameter :: Gmsh_Format = 2
  integer, parameter :: Plot3D_Format = 3
  integer, parameter :: CGNS_Format = 4
  !
  ! CGNS constants
  integer, parameter :: reopen_CGNS = 2
  integer, parameter :: open_CGNS   =  1
  integer, parameter :: add_to_CGNS =  0
  integer, parameter :: close_CGNS  = -1
  integer, parameter :: dump_CGNS   = -2
  !
  ! Cell Geometry Family Types
  !
  type :: geom_family_t
    integer :: v = 1
  end type geom_family_t
  !
  type(geom_family_t), parameter :: Complete     = geom_family_t( v=1 )
  type(geom_family_t), parameter :: Edge_Defined = geom_family_t( v=2 )
  type(geom_family_t), parameter :: Face_Defined = geom_family_t( v=3 )
  !
  interface operator(==)
    module procedure geom_family_type_equals
  end interface
  !
  interface operator(/=)
    module procedure geom_family_type_not_equal
  end interface
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
  integer, protected, save :: nodes_to_geom(0:3,0:8)
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
  logical(lk), parameter :: tria_is_biunit = fals
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
  ! MPI Public Scalars
  !
  integer(int_mpi), protected, save :: mypnum = 0_int_mpi
  integer(int_mpi), protected, save :: ncpu = 1_int_mpi
  integer(int_mpi),            save :: mpierr
  !
  integer(int_mpi), protected, save :: mpi_comm_world_val
  integer(int_mpi), protected, save :: mpi_comm_self_val
  !
  integer(int_mpi), parameter :: abort = int( 8 , kind=kind(abort) )
  integer(int_mpi), parameter :: stop_mpi = int( 9 , kind=kind(stop_mpi) )
  !
  integer(int_MPI), parameter :: glb_root = 0_int_mpi
  integer(int_MPI), parameter :: host_root = 0_int_mpi
  integer(int_MPI), parameter :: cpu_root = 0_int_mpi
  !
  logical(lk), save :: i_am_cpu_root = true
  logical(lk), save :: i_am_host_root = true
  !
  _MPI_COMM_TYPE_, save :: my_host_comm
  _MPI_COMM_TYPE_, save :: my_cpu_comm
  _MPI_COMM_TYPE_, save :: host_roots_comm
  _MPI_COMM_TYPE_, save :: cpu_roots_comm
  !
  integer, save :: my_cpu_num = 0
  integer, save :: my_host_num = 0
  !
  character(len=MPI_MAX_PROCESSOR_NAME) :: my_host_name
  !
  integer(int_mpi), save :: my_cpu_rank = 0_int_mpi
  integer(int_mpi), save :: my_host_rank = 0_int_mpi
  !
  integer(int_mpi), save :: my_host_root = 0_int_mpi
  integer(int_mpi), save :: my_cpu_root = 0_int_mpi
  !
  integer(int_mpi), allocatable, save :: my_cpu_grp(:)
  integer(int_mpi), allocatable, save :: my_host_grp(:)
  !
  integer(int_mpi), allocatable, save :: cpu_roots_grp(:)
  integer(int_mpi), allocatable, save :: host_roots_grp(:)
  !
  ! MPI datatypes to match Fortran intrinsic datatypes
  !
  ! working real precision
  _MPI_DATA_TYPE_, protected, save :: mpi_flttyp
  ! half working real precision
  _MPI_DATA_TYPE_, protected, save :: mpi_hlftyp
  ! double working real precision
  _MPI_DATA_TYPE_, protected, save :: mpi_dbltyp
  ! working integer precision
  _MPI_DATA_TYPE_, protected, save :: mpi_inttyp
  !
  logical(lk), private, save :: random_is_initialized = fals
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
  logical(lk), parameter :: skip_alloc_pause = fals
  !
  ! Parameters for procedure debug_timer
  !
  integer, parameter :: start_timer = 0
  integer, parameter :: stop_timer = 1
  integer, parameter :: entering_procedure = 2
  integer, parameter :: leaving_procedure = 3
  integer, parameter :: start_accumulation = 4
  integer, parameter :: stop_accumulation = 5
  integer, parameter :: reset_accumulated_time = 6
  integer, parameter :: report_accumulated_time = 7
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
  private :: cross_product_SP,cross_product_DP,cross_product_QP
  !
  interface cross_product
    module procedure cross_product_SP
    module procedure cross_product_DP
#ifndef DISABLE_QP
    module procedure cross_product_QP
#endif
  end interface cross_product
  !
  !
  !
  private :: deltafun_DP_LDK,deltafun_DP_LK
  !
  interface deltafun
    module procedure deltafun_DP_LDK
    module procedure deltafun_DP_LK
  end interface deltafun
  !
  !
  !
  private :: deltafun_QP_LDK,deltafun_QP_LK
  !
  interface deltafun_QP
    module procedure deltafun_QP_LDK
    module procedure deltafun_QP_LK
  end interface deltafun_QP
  !
  !
  !
  private :: selfdot_SP,selfdot_DP,selfdot_QP
  !
  interface selfdot
    module procedure selfdot_SP
    module procedure selfdot_DP
#ifndef DISABLE_QP
    module procedure selfdot_QP
#endif
  end interface selfdot
  !
  !
  !
  private :: rnorm2_SP,rnorm2_DP,rnorm2_QP
  !
#ifdef SPECIAL_FOR_PGI
  interface norm2
    module procedure rnorm2_SP
    module procedure rnorm2_DP
#ifndef DISABLE_QP
    module procedure rnorm2_QP
#endif
  end interface norm2
#endif
  !
  !
  !
  private :: unit_vector_SP,unit_vector_DP,unit_vector_QP
  !
  interface unit_vector
    module procedure unit_vector_SP
    module procedure unit_vector_DP
#ifndef DISABLE_QP
    module procedure unit_vector_QP
#endif
  end interface unit_vector
  !
  !
  !
  private :: outerprod_DP,outerprod_QP
  !
  interface outerprod
    module procedure outerprod_DP
#ifndef DISABLE_QP
    module procedure outerprod_QP
#endif
  end interface outerprod
  !
  !
  !
  private :: reciprocal_SP,reciprocal_DP,reciprocal_QP
  !
  interface reciprocal
    module procedure reciprocal_SP
    module procedure reciprocal_DP
#ifndef DISABLE_QP
    module procedure reciprocal_QP
#endif
  end interface reciprocal
  !
  !
  !
  private :: chop_WP,chop_SP,chop_DP,chop_QP
  !
  interface chop
    module procedure chop_WP
    module procedure chop_SP
    module procedure chop_DP
#ifndef DISABLE_QP
    module procedure chop_QP
#endif
  end interface chop
  !
  !
  !
  private :: reflect_DP,reflect_QP
  !
  interface reflect
    module procedure reflect_DP
#ifndef DISABLE_QP
    module procedure reflect_QP
#endif
  end interface reflect
  !
  !
  !
  private :: rotate_2d_DP,rotate_2d_QP
  !
  interface rotate_2d
    module procedure rotate_2d_DP
#ifndef DISABLE_QP
    module procedure rotate_2d_QP
#endif
  end interface rotate_2d
  !
  !
  !
  private :: rotate_DP,rotate_QP
  !
  interface rotate
    module procedure rotate_DP
#ifndef DISABLE_QP
    module procedure rotate_QP
#endif
  end interface rotate
  !
  !
  !
  private :: node_center_SP,node_center_DP,node_center_QP
  !
  interface node_center
    module procedure node_center_SP
    module procedure node_center_DP
#ifndef DISABLE_QP
    module procedure node_center_QP
#endif
  end interface node_center
  !
  !
  !
  private :: face_normal_SP,face_normal_DP,face_normal_QP
  !
  interface face_normal
    module procedure face_normal_SP
    module procedure face_normal_DP
#ifndef DISABLE_QP
    module procedure face_normal_QP
#endif
  end interface face_normal
  !
  !
  !
  private :: compute_jacobian_SP,compute_jacobian_DP,compute_jacobian_QP
  !
  interface compute_jacobian
    module procedure compute_jacobian_SP
    module procedure compute_jacobian_DP
#ifndef DISABLE_QP
    module procedure compute_jacobian_QP
#endif
  end interface compute_jacobian
  !
  !
  !
  private :: reverse_metrics_SP,reverse_metrics_DP,reverse_metrics_QP
  !
  interface reverse_metrics
    module procedure reverse_metrics_SP
    module procedure reverse_metrics_DP
#ifndef DISABLE_QP
    module procedure reverse_metrics_QP
#endif
  end interface reverse_metrics
  !
  !
  !
  private :: invert_matrix_SP,invert_matrix_DP,invert_matrix_QP
  !
  interface invert_matrix
    module procedure invert_matrix_SP
    module procedure invert_matrix_DP
#ifndef DISABLE_QP
    module procedure invert_matrix_QP
#endif
  end interface invert_matrix
  !
  !
  !
  private :: invert_matrix90_SP,invert_matrix90_DP,invert_matrix90_QP
  !
  interface invert_matrix90
    module procedure invert_matrix90_SP
    module procedure invert_matrix90_DP
#ifndef DISABLE_QP
    module procedure invert_matrix90_QP
#endif
  end interface invert_matrix90
  !
  !
  !
  private :: ludcmp_SP,ludcmp_DP,ludcmp_QP
  !
  interface ludcmp
    module procedure ludcmp_SP
    module procedure ludcmp_DP
#ifndef DISABLE_QP
    module procedure ludcmp_QP
#endif
  end interface ludcmp
  !
  !
  !
  private :: lubksb_SP,lubksb_DP,lubksb_QP
  !
  interface lubksb
    module procedure lubksb_SP
    module procedure lubksb_DP
#ifndef DISABLE_QP
    module procedure lubksb_QP
#endif
  end interface lubksb
  !
  !
  !
  private :: determinant_SP,determinant_DP,determinant_QP
  !
  interface determinant
    module procedure determinant_SP
    module procedure determinant_DP
#ifndef DISABLE_QP
    module procedure determinant_QP
#endif
  end interface determinant
  !
  !
  !
  private :: random_i4,random_sp,random_dp,random_qp
  !
  interface random
    module procedure random_i4
    module procedure random_sp
    module procedure random_dp
#ifndef DISABLE_QP
    module procedure random_qp
#endif
  end interface random
  !
  !
  !
  private :: findloc_i4,recursive_findloc_i4
  private :: findloc_i8,recursive_findloc_i8
  private :: findloc_r4,recursive_findloc_r4
  private :: findloc_r8,recursive_findloc_r8
  !
  interface findloc
    module procedure findloc_i4,recursive_findloc_i4, &
                     findloc_i8,recursive_findloc_i8, &
                     findloc_r4,recursive_findloc_r4, &
                     findloc_r8,recursive_findloc_r8

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
    module procedure identity_matrix_dim_SP
    module procedure identity_matrix_mat_SP
    module procedure identity_matrix_dim_DP
    module procedure identity_matrix_mat_DP
#ifndef DISABLE_QP
    module procedure identity_matrix_dim_QP
    module procedure identity_matrix_mat_QP
#endif
  end interface identity_matrix
  !
  !
  !
  private :: intseq_no_interval_i1,intseq_with_interval_i1
  private :: intseq_no_interval_i2,intseq_with_interval_i2
  private :: intseq_no_interval_i4,intseq_with_interval_i4
  private :: intseq_no_interval_i8,intseq_with_interval_i8
  !
  interface intseq
    module procedure &
      intseq_no_interval_i1,intseq_with_interval_i1,intseq_rev_dir_i1, &
      intseq_no_interval_i2,intseq_with_interval_i2,intseq_rev_dir_i2, &
      intseq_no_interval_i4,intseq_with_interval_i4,intseq_rev_dir_i4, &
      intseq_no_interval_i8,intseq_with_interval_i8,intseq_rev_dir_i8
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
  !
  !
  private :: is_nan_SP,is_nan_DP,is_nan_QP
  !
  interface is_nan
#ifdef SPECIAL_FOR_GCC
    module procedure is_nan_SP
    module procedure is_nan_DP
#ifndef DISABLE_QP
    module procedure is_nan_QP
#endif
#endif
  end interface is_nan
  !
  !
  !
  private :: is_finite_SP,is_finite_DP,is_finite_QP
  !
  interface is_finite
#ifdef SPECIAL_FOR_GCC
    module procedure is_finite_SP
    module procedure is_finite_DP
#ifndef DISABLE_QP
    module procedure is_finite_QP
#endif
#endif
  end interface is_finite
  !
  !
  !
  private :: is_negative_SP,is_negative_DP,is_negative_QP
  !
  interface is_negative
#ifdef SPECIAL_FOR_GCC
    module procedure is_negative_SP
    module procedure is_negative_DP
#ifndef DISABLE_QP
    module procedure is_negative_QP
#endif
#endif
  end interface is_negative
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
  nodes_to_geom = 0
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
  ! Define the array to get the geometry type knowing
  ! the number of nodes and the dimension of the element
  !
  do this_geom = Geom_Min,Geom_Max
    nodes_to_geom( geom_dimen(this_geom) , &
                   geom_nodes(this_geom) ) = this_geom
  end do
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
  edge_nodes(1:2, 1,Geom_Hexa) = [0,1]
  edge_nodes(1:2, 2,Geom_Hexa) = [1,2]
  edge_nodes(1:2, 3,Geom_Hexa) = [2,3]
  edge_nodes(1:2, 4,Geom_Hexa) = [3,0]
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
 !call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  !
end subroutine initialize_runtime_parameters
!
!###############################################################################
!
subroutine alloc_error(procedure_name,variable_name,ioptype,line_number, &
                       file_name,error_flag,iemessag,override_memory_pause)
  !
  ! Subroutine to generate an error message when
  ! an allocate or deallocate returns an error code
  !
  !.. Formal Arguments ..
  character(len=*),           intent(in) :: procedure_name
  character(len=*),           intent(in) :: variable_name
  integer,                    intent(in) :: ioptype
  integer,                    intent(in) :: line_number
  character(len=*),           intent(in) :: file_name
  integer,                    intent(in) :: error_flag
  character(len=*), optional, intent(in) :: iemessag
  logical(lk),      optional, intent(in) :: override_memory_pause
  !
#ifdef TRACK_MEMORY_ALLOCATIONS
  !.. Local Scalars ..
  logical(lk) :: alloc_memory_pause
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
    write (iout,4) trim(adjustl(procedure_name)),line_number, &
                   trim(adjustl(file_name)),error_flag
    if (present(iemessag)) write (iout,5) trim(adjustl(iemessag))
    write (iout,6)
    !
    call stop_gfr(abort,pname,__LINE__,__FILE__)
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
            "         ON LINE # ",i0,/, &
            "      OF THE FILE: ",a,/, &
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
    call mpi_bcast(stop_value,1_int_mpi,mpi_inttyp, &
                   0_int_mpi,MPI_COMM_WORLD,mpierr)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
  if (stop_value == 1) then
    error_message = "Abort issused during memory test at: " // &
                    trim(adjustl(after_message))
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
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
subroutine io_error(procedure_name,io_file_name,ioptype,line_number, &
                    src_file_name,error_flag,optmessag)
  !
  ! Subroutine to generate an error message when
  ! an open or close returns an error code
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: procedure_name
  character(len=*), intent(in) :: io_file_name
  integer,          intent(in) :: ioptype
  integer,          intent(in) :: line_number
  character(len=*), intent(in) :: src_file_name
  integer,          intent(in) :: error_flag
  !
  !.. Optional Arguments ..
  character(len=*), optional, intent(in) :: optmessag
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "io_error"
  !
continue
  !
  if (error_flag /= 0) then
    !
    write (iout,10)
    !
    ! Open/Close errors
    !
    if (ioptype == 1) then
      write (iout,1) trim(adjustl(io_file_name))
    else if (ioptype == 2) then
      write (iout,2) trim(adjustl(io_file_name))
    else if (ioptype == 3) then
      write (iout,3) trim(adjustl(io_file_name))
    else if (ioptype == 4) then
      write (iout,4) trim(adjustl(io_file_name))
    else if (ioptype == 5) then
      write (iout,5) trim(adjustl(io_file_name))
      !
      ! Read/Write errors
      ! Use the value of error_flag to determine the type of error
      !
    else if (ioptype == -1) then
      !
      if (error_flag > 0) then
        write (iout,6) trim(adjustl(io_file_name))
      else if (error_flag < 0) then
        write (iout,7) trim(adjustl(io_file_name))
      end if
      !
      ! Generic unknown error
      !
    else
      write (iout,99) trim(adjustl(io_file_name))
    end if
    !
    if (present(optmessag)) then
      write (iout,8) trim(adjustl(optmessag))
    end if
    !
    write (iout,11) trim(adjustl(procedure_name)),error_flag
    call stop_gfr(abort,pname,line_number,src_file_name)
    !
  end if
  !
  ! Format Statements
  !
  1  format (" ERROR OPENING FILE '",a,"'")
  2  format (" ERROR CLOSING FILE '",a,"'")
  3  format (" ERROR REOPENING CGNS FILE '",a,"'")
  4  format (" UNABLE TO DETERMINE FORMAT OF PLOT3D FILE '",a,"'")
  5  format (" CGNS ERROR ENCOUNTERED WITH FILE '",a,"'")
  6  format (" UNEXPECTED INPUT/OUTPUT ERROR WITH FILE '",a,"'")
  7  format (" UNEXPECTED END OF FILE/RECORD ERROR WITH FILE '",a,"'")
 !8  format (" ERROR OCCURED TRYING TO READ IN '",a,"'")
  8  format (" ERROR MESSAGE: '",a,"'")
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
subroutine execute_command(command,wait_for_completion,exit_status, &
                           command_status,command_message,procedure_name)
  !
#define USE_EXECUTE_COMMAND_LINE
#ifdef USE_EXECUTE_COMMAND_LINE
#ifdef SPECIAL_FOR_PGI
  use compatibility_pgi_mod, only : execute_command_line
#endif
#else
#ifdef SPECIAL_FOR_PGI
  use compatibility_pgi_mod, only : system,ierrno,gerror
#else
#ifdef SPECIAL_FOR_GCC
  use compatibility_gcc_mod, only : system,ierrno,gerror
#else
  use ifport, only : system,ierrno
  use ifcore, only : gerror
#endif
#endif
#endif
  !.. Formal Arguments ..
  character(len=*),  intent(in) :: command
  logical(ldk),      intent(in) :: wait_for_completion
  integer,          intent(out) :: exit_status
  integer,          intent(out) :: command_status
  character(len=*), intent(out) :: command_message
  character(len=*),  intent(in) :: procedure_name
  !
  !.. Local Scalars ..
  integer :: n
  !
  !.. Allocatable Scalars ..
  character(len=:), allocatable :: command_string
 !character(len=len_trim(adjustl(command))) :: command_string
  !
  !.. Local Parameters ..
#ifdef USE_EXECUTE_COMMAND_LINE
  character(len=*), parameter :: command_line_procedure = "execute_command_line"
#else
  character(len=*), parameter :: command_line_procedure = "system"
#endif
  !
continue
  !
  command_string = trim(adjustl(command))
  !
  if (mypnum == 0) then
    !
    do n = 1,5
      !
      exit_status = 0
      command_status = 0
      command_message = ""
      !
#ifdef USE_EXECUTE_COMMAND_LINE
      call execute_command_line(command_string,wait_for_completion, &
                                exit_status,command_status,command_message )
#else
      exit_status = system( command_string )
      if (exit_status /= 0) then
        command_status = ierrno()
        call gerror(command_message)
      end if
#endif
      !
      call write_command_results(procedure_name,command,exit_status, &
                                 command_status,command_message, &
                                 command_line_procedure)
      !
      if (exit_status == 0) exit
      !
      flush (iout)
      !
    end do
    !
  end if
  !
end subroutine execute_command
!
!###############################################################################
!
subroutine write_command_results(procedure_name,command,exit_status, &
                                 command_status,command_message, &
                                 command_line_procedure)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: procedure_name
  character(len=*), intent(in) :: command
  integer,          intent(in) :: exit_status
  !
  !.. Optional Arguments ..
  integer,          optional, intent(in) :: command_status
  character(len=*), optional, intent(in) :: command_message
  character(len=*), optional, intent(in) :: command_line_procedure
  !
  !.. Local Scalars ..
  character(len=20) :: clp
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_command_results"
  !
continue
  !
  clp = "execute_command_line"
  if (present(command_line_procedure)) then
    clp = trim(adjustl(command_line_procedure))
  end if
  !
  if (mypnum == 0) then
    !
    write (iout,1) trim(adjustl(clp)), &
                   trim(adjustl(command)), &
                   trim(adjustl(procedure_name)), &
                   exit_status
    !
    if (present(command_message)) then
      write (iout,6) trim(adjustl(clp)),trim(adjustl(command_message))
    end if
    !
    if (present(command_status)) then
      !
      write (iout,5) trim(adjustl(clp)),command_status
      !
      if (command_status == -1) then
        write (iout,2)
      else if (command_status == -2) then
        write (iout,3)
      else if (command_status /= 0) then
        write (iout,4) trim(adjustl(clp))
      end if
      !
    end if
    !
    write (iout,*)
    flush (iout)
    !
  end if
  !
  ! Format Statements
  !
  1 format (/,2x,"Calling ",a," with the command",/, &
              5x,"'",a,"'",/, &
              2x,"from procedure '",a,"'",/, &
              5x,"status returned from command line = ",i0)
  2 format (/,2x,"ERROR RETURNED INDICATES THIS PROCESSOR ", &
                 "DOES NOT SUPPORT COMMAND LINE EXECUTION!!!",/)
  3 format (/,2x,"ERROR RETURNED INDICATES THIS PROCESSOR DOES ", &
                 "NOT SUPPORT ASYNCHRONOUS COMMAND EXECUTION!!!",/)
  4 format (/,2x,"AN ERROR OCCURRED CALLING ",a,"() !!!",/)
  5 format (5x,"status  returned from ",a,"() = ",i0)
  6 format (5x,"message returned from ",a,"() = '",a,"'")

  !
end subroutine write_command_results
!
!###############################################################################
!
subroutine initialize_mpi_environment
  !
  !.. Local Scalars ..
  integer(INT_MPI) :: r,p
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_mpi_environment"
continue
  !
  ! Initialize the overall MPI session
  !
  call mpi_init(mpierr)
  call mpi_comm_size(MPI_COMM_WORLD,ncpu,mpierr)
  call mpi_comm_rank(MPI_COMM_WORLD,mypnum,mpierr)
  !
  mpi_comm_world_val = _MPI_COMM_WORLD_VAL_
  mpi_comm_self_val = _MPI_COMM_SELF_VAL_
  !
  ! MPI datatype for working integer precision
  !
 !r = range( 0 )
 !write (*,*) "integer precision"
 !write (*,*) r
 !write (*,*) selected_int_kind(r)
 !call mpi_type_create_f90_integer(r,mpi_inttyp,mpierr)
 !write (*,*) mpi_inttyp
  if (WIP == I8) then
    mpi_inttyp = MPI_INTEGER8
  else
    mpi_inttyp = MPI_INTEGER
  end if
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
  if (SP < R8) then
    mpi_hlftyp = MPI_REAL
  else
    mpi_hlftyp = MPI_DOUBLE_PRECISION
  end if
  if (WP < R8) then
    mpi_flttyp = MPI_REAL
  else
    mpi_flttyp = MPI_DOUBLE_PRECISION
  end if
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
 !mpi_flttyp = MPI_DOUBLE_PRECISION
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
 !  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
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
subroutine stop_gfr(iend,procedure_name,line_number,file_name, &
                      optional_message)
  !
  !.. Formal Arguments ..
  integer,          intent(in) :: iend
  character(len=*), intent(in) :: procedure_name
  integer,          intent(in) :: line_number
  character(len=*), intent(in) :: file_name
  !
  !.. Optional Arguments ..
  character(len=*), optional, intent(in) :: optional_message
  !
  !.. Local Scalars ..
  integer(int_mpi) :: iflag
  character(len=10000) :: strng
  character(len=10000) :: output_message
  !
continue
  !
  iflag = int( iend , kind=kind(iflag) )
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
        write (iout,4) trim(adjustl(procedure_name)),line_number, &
                       trim(adjustl(file_name))
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
      write (iout,4) trim(adjustl(procedure_name)),line_number, &
                     trim(adjustl(file_name))
      if (present(optional_message)) &
        write (iout,5) trim(adjustl(output_message))
      write (iout,6)
      !
      call mpi_abort(MPI_COMM_WORLD,iflag,mpierr)
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
  2 format ("  CALLING PREMATURE MPI_FINALIZE ON ALL PROCESSORS! ")
  3 format ("  MPI ABORT ISSUED FROM PROCESSOR ",i0,"!")
  4 format ("         STOP ISSUED IN PROCEDURE '",a,"'",/, &
            "                        ON LINE #  ",i0,/, &
            "                      OF THE FILE '",a,"'")
  5 format (" ERROR MESSAGE: '",a,"'")
  6 format (66("#"),//)
  !
end subroutine stop_gfr
!
!###############################################################################
!
subroutine gfr_timer(current_wall_time)
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: current_wall_time
  !
  !.. Local Scalars ..
  integer(i8) :: clock_count
  integer(i8) :: clock_count_rate
  !
continue
  !
  call system_clock(count=clock_count,count_rate=clock_count_rate)
  current_wall_time = real(clock_count,kind=wp) / real(clock_count_rate,kind=wp)
  !
end subroutine gfr_timer
!
!###############################################################################
!
subroutine debug_timer(time_flag,message,processor,part_num,acc)
  !
  !.. Formal Arguments ..
  integer,          intent(in) :: time_flag
  !
  character(len=*), intent(in) :: message
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: processor
  integer, optional, intent(in) :: part_num
  integer, optional, intent(in) :: acc
  !
#ifdef PROFILE_ON
  !.. Local Scalars ..
  real(wp), save :: beg_wall_time
  real(wp)       :: end_wall_time
  real(wp)       :: wall_time
  real(wp)       :: temp_time
  !
  integer, save :: procedure_counter = 0
  integer       :: i,output_cpu
  !
  character(len=300) :: out_message
  !
  !.. Local Arrays ..
  real(wp), save :: beg_procedure_time(10)
  real(wp), save :: accumulated_time(10)
  real(wp), save :: accumulation_beg_time(10)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: leading_space = "                    "
  !
continue
  !
  output_cpu = 0
  if (present(processor)) output_cpu = processor
  !
  if (time_flag == start_timer) then
    !
    if (present(part_num)) then
      write (out_message,6) message,part_num
    else
      write (out_message,5) message
    end if
    !
    i = 3 + procedure_counter*2
    !
    call gfr_timer( beg_wall_time )
    !
    if (mypnum == output_cpu) then
      write (iout,1) leading_space(1:i),trim(adjustl(out_message))
    end if
    !
  else if (time_flag == stop_timer) then
    !
    call gfr_timer( end_wall_time )
    !
    wall_time = end_wall_time - beg_wall_time
    !
    if (present(part_num)) then
      write (out_message,6) message,part_num
    else
      write (out_message,5) message
    end if
    !
    i = 3 + procedure_counter*2
    !
    if (mypnum == output_cpu) then
      write (iout,2) leading_space(1:i),wall_time,trim(adjustl(out_message))
    end if
    !
    ! Reset beg_wall_time in case this procedure is called twice
    ! in a row with time_flag set to stop_timer
    !
    call gfr_timer( beg_wall_time )
    !
  else if (time_flag == reset_accumulated_time) then
    !
    if (present(acc)) then
      accumulated_time(acc) = zero
    end if
    !
  else if (time_flag == report_accumulated_time) then
    !
    if (present(part_num)) then
      write (out_message,6) message,part_num
    else
      write (out_message,5) message
    end if
    !
    i = 3 + procedure_counter*2
    !
    if (mypnum == output_cpu) then
      write (iout,2) leading_space(1:i),accumulated_time(acc), &
                     trim(adjustl(out_message))
    end if
    !
  else if (time_flag == start_accumulation) then
    !
    if (present(acc)) then
      call gfr_timer( accumulation_beg_time(acc) )
    end if
    !
  else if (time_flag == stop_accumulation) then
    !
    if (present(acc)) then
      !
      call gfr_timer( temp_time )
      !
      wall_time = temp_time - accumulation_beg_time(acc)
      !
      accumulated_time(acc) = accumulated_time(acc) + wall_time
      !
    end if
    !
  else if (time_flag == entering_procedure) then
    !
    i = 1 + procedure_counter*2
    !
    if (mypnum == output_cpu) then
      write (iout,3) leading_space(1:i),trim(adjustl(message))
    end if
    !
    procedure_counter = procedure_counter + 1
    !
    call gfr_timer( beg_procedure_time(procedure_counter) )
    !
  else if (time_flag == leaving_procedure) then
    !
    call gfr_timer( end_wall_time )
    !
    wall_time = end_wall_time - beg_procedure_time(procedure_counter)
    !
    procedure_counter = max( 0 , procedure_counter-1 )
    !
    i = 1 + procedure_counter*2
    !
    if (mypnum == output_cpu) then
      write (iout,4) leading_space(1:i),trim(adjustl(message)),wall_time
    end if
    !
  end if
  !
  flush (iout)
  !
  ! Format Statements
  !
  1 format (a,"## Starting timer for: ",a)
  2 format (a,"     -- TIME WAS",f12.6," seconds for ",a)
  3 format (a,"=> Entering procedure : ",a)
  4 format (a," <= Leaving procedure : ",a,"   after",f12.6," seconds")
  5 format (a)
  6 format (a," for partition # ",i0)
  !
#endif
  !
end subroutine debug_timer
!
!##############################################################################
!
elemental function geom_family_type_equals(a,b) result(return_value)
  type(geom_family_t), intent(in) :: a,b
  logical(lk) :: return_value
  return_value = (a%v == b%v)
end function geom_family_type_equals
!
!###############################################################################
!
elemental function geom_family_type_not_equal(a,b) result(return_value)
  type(geom_family_t), intent(in) :: a,b
  logical(lk) :: return_value
  return_value = (a%v /= b%v)
end function geom_family_type_not_equal
!
!###############################################################################
!
! Include statements for functions that have a generic interface
!
#include "Functions/reallocate.f90"
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
include "Functions/checks_for_bad_numbers.f90"
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
#ifndef DISABLE_DTIO
    generic :: write(formatted) => memout
#endif
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
