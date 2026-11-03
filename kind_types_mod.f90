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
#ifdef SPECIAL_FOR_NVIDIA
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
#ifdef SPECIAL_FOR_GNU
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
#ifdef SPECIAL_FOR_NVIDIA
  integer, public, parameter :: l1k = logical_kinds(1)
#else
#ifdef SPECIAL_FOR_GNU
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
module ieee_mod
  !
#ifndef SPECIAL_FOR_GNU
  !
  use, intrinsic :: ieee_arithmetic, only : is_nan => ieee_is_nan
  use, intrinsic :: ieee_arithmetic, only : is_finite => ieee_is_finite
  use, intrinsic :: ieee_arithmetic, only : is_negative => ieee_is_negative
  !
  use, intrinsic :: ieee_exceptions
  !
#endif
  !
  implicit none
  !
  public
  !
  character(len=100), save :: src_file_with_bad_number = ""
  integer,            save :: line_num_with_bad_number = 0
  !
#ifndef SPECIAL_FOR_GNU
  !
  type(ieee_status_type) :: default_ieee_status
  !
#else
  !
  integer, parameter :: ieee_usual(1:5) = 0
  !
  ! From Functions/checks_for_bad_numbers.f90
  !
  private :: is_nan_SP,is_nan_DP,is_nan_QP
  !
  interface is_nan
    module procedure is_nan_SP
    module procedure is_nan_DP
#ifdef ENABLE_QP
    module procedure is_nan_QP
#endif
  end interface is_nan
  !
  ! From Functions/checks_for_bad_numbers.f90
  !
  private :: is_finite_SP,is_finite_DP,is_finite_QP
  !
  interface is_finite
    module procedure is_finite_SP
    module procedure is_finite_DP
#ifdef ENABLE_QP
    module procedure is_finite_QP
#endif
  end interface is_finite
  !
  ! From Functions/checks_for_bad_numbers.f90
  !
  private :: is_negative_SP,is_negative_DP,is_negative_QP
  !
  interface is_negative
    module procedure is_negative_SP
    module procedure is_negative_DP
#ifdef ENABLE_QP
    module procedure is_negative_QP
#endif
  end interface is_negative
  !
#endif
  !
contains
!
!###############################################################################
!
subroutine get_default_ieee_status()
  ! This should only be called at code startup
continue
#ifndef SPECIAL_FOR_GNU
  call ieee_get_status(default_ieee_status)
#endif
end subroutine get_default_ieee_status
!
!###############################################################################
!
subroutine reset_ieee_status()
continue
#ifndef SPECIAL_FOR_GNU
  call ieee_set_status(default_ieee_status)
#endif
  src_file_with_bad_number = ""
  line_num_with_bad_number = 0
end subroutine reset_ieee_status
!
!###############################################################################
!
subroutine reset_fp_exceptions()
  !
  !.. Local Scalars ..
  integer :: n
  !
  !.. Local Arrays ..
  logical :: reset_exceptions(1:size(ieee_usual))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "reset_fp_exceptions"
  !
continue
  !
  do n = 1,size(reset_exceptions)
    reset_exceptions(n) = .false.
  end do
  !
#ifndef SPECIAL_FOR_GNU
  call ieee_set_flag(ieee_usual,reset_exceptions)
#endif
  !
  src_file_with_bad_number = ""
  line_num_with_bad_number = 0
  !
end subroutine reset_fp_exceptions
!
!###############################################################################
!
#include "Functions/checks_for_bad_numbers.f90"
!
!###############################################################################
!
end module ieee_mod
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
#ifdef USE_INTEL_MKL_BLAS
  use blas95, only : dot
#else
  use rename_intrinsics, only : dot => dot_product
#endif
  use ieee_mod
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
  !--------------------------------
  ! Aliases for logical kind types
  !--------------------------------
  !
  ! Logical kind of least precision (only needs true/false values)
  integer, parameter :: LK = l1k
  !
  !--------------------------------
  ! Aliases for real kind types
  !--------------------------------
  !
  ! 'Single Precision' (32-bit/4-byte) real
  integer, parameter :: SP = r4
  ! 'Double Precision' (64-bit/8-byte) real
  integer, parameter :: DP = r8
  ! 'Quadruple Precision' (128-bit/16-byte) real
  ! (or DP if ENABLE_QP is not active)
#ifdef ENABLE_QP
  integer, parameter :: QP = r16
#else
  integer, parameter :: QP = r8
#endif
  ! 'Working Precision' real: Used for most real scalars/arrays
  integer, parameter :: WP = r8
  ! 'Matrix Precision' real: Used for polynomial evals and creating matrices
  integer, parameter :: MP = QP
  !
  !--------------------------------
  ! Aliases for integer kind types
  !--------------------------------
  !
  ! 'Working Integer Precision': for now default integer precision (32-bit)
  integer, parameter :: WIP = i4
  ! 'MPI Integer Precision': determined by MPI constant MPI_INTEGER_KIND
  integer, parameter :: INT_MPI = kind(MPI_INTEGER_KIND)
  ! 'METIS Integer Precision': for now default integer precision (32-bit)
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
  real(wp), parameter :: eps4  = 1.0e-04_wp
  real(wp), parameter :: eps5  = 1.0e-05_wp
  real(wp), parameter :: eps6  = 1.0e-06_wp
  real(wp), parameter :: eps7  = 1.0e-07_wp
  real(wp), parameter :: eps8  = 1.0e-08_wp
  real(wp), parameter :: eps9  = 1.0e-09_wp
  real(wp), parameter :: eps10 = 1.0e-10_wp
  real(wp), parameter :: eps11 = 1.0e-11_wp
  real(wp), parameter :: eps12 = 1.0e-12_wp
  real(wp), parameter :: eps13 = 1.0e-13_wp
  real(wp), parameter :: eps14 = 1.0e-14_wp
  real(wp), parameter :: eps15 = 1.0e-15_wp
  real(wp), parameter :: eps16 = 1.0e-16_wp
  real(wp), parameter :: eps17 = 1.0e-17_wp
  real(wp), parameter :: small = eps17
  !
  real(wp), parameter :: ten2  = ten**2
  real(wp), parameter :: ten3  = ten**3
  !
  real(wp),  parameter :: pi = real(acos(-1.0_r16),kind=wp)
 !real(wp),  parameter :: pi = 3.14159265358979323846_wp
  !
  real(wp), parameter :: one80   = two*nine*ten
  real(wp), parameter :: three60 = two*one80
  real(wp), parameter :: radians2degrees = one80/pi
  real(wp), parameter :: degrees2radians = one/radians2degrees
  !
  real(wp), parameter :: inches_to_meters = 0.0254_wp
  real(wp), parameter :: meters_to_inches = one/inches_to_meters
  !
  real(wp), parameter :: rpm_to_mps = two*pi/(six*ten)
  !
  real(wp), parameter :: rndoff  = epsilon(one)
  real(wp), parameter :: tny     = tiny(one)
  real(wp), parameter :: exp_lim = real( int(log(huge(zero))) , kind=wp )
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
  integer, parameter :: bc_sem = 5
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
  integer, parameter :: bc_unsupported = -9998
  integer, parameter :: bc_unknown = -9999
  ! BC Groups
  integer, parameter :: bc_comm(1:2) = [bc_periodic,bc_cpu_bnd]
  integer, parameter :: bc_walls(1:5) = [ bc_slip_wall, &
                                          bc_euler_wall, &
                                          bc_symmetry, &
                                          bc_adiabatic_wall, &
                                          bc_isothermal_wall ]
  integer, parameter :: bc_physical_walls(1:4) = [ bc_slip_wall, &
                                                   bc_euler_wall, &
                                                   bc_adiabatic_wall, &
                                                   bc_isothermal_wall ]
  integer, parameter :: bc_noslip_walls(1:2) = [ bc_adiabatic_wall, &
                                                 bc_isothermal_wall ]
  integer, parameter :: bc_inflows(1:9) = [ bc_generic_inflow, &
                                            bc_sub_inflow, &
                                            bc_sup_inflow, &
                                            bc_mdot_inflow, &
                                            bc_sem, &
                                            bc_characteristic, &
                                            bc_generic_freeflow, &
                                            bc_freestream, &
                                            bc_fixed ]
  integer, parameter :: bc_outflows(1:6) = [ bc_generic_outflow, &
                                             bc_sub_outflow, &
                                             bc_sup_outflow, &
                                             bc_mdot_outflow, &
                                             bc_characteristic, &
                                             bc_generic_freeflow ]
  !$ACC DECLARE COPYIN(bc_walls)
  !
  integer, parameter :: bc_heirarchy(1:18) = [ bc_isothermal_wall, &
                                               bc_adiabatic_wall, &
                                               bc_euler_wall, &
                                               bc_slip_wall, &
                                               bc_symmetry, &
                                               bc_sem, &
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
  integer, parameter :: Generic_MeanFlow        = 1
  integer, parameter :: Channel_Flow            = 2
  integer, parameter :: Inviscid_Gaussian_Bump  = 3
  integer, parameter :: Laminar_BndLyr          = 4
  integer, parameter :: Shu_Vortex              = 5
  integer, parameter :: Shu_MultiVortex         = 55
  integer, parameter :: Vortex_Transport        = 6
  integer, parameter :: Density_Transport       = 7
  integer, parameter :: Entropy_Transport       = 8
  integer, parameter :: Taylor_Green_Vortex     = 9
  integer, parameter :: Subsonic_Jet            = 10
  integer, parameter :: Supersonic_Jet          = 11
  integer, parameter :: Infinite_Cylinder       = 12
  integer, parameter :: Freestream_Preservation = 13
  integer, parameter :: Reacting_Backstep       = 14
  integer, parameter :: H2_2D_Diffusion_Flame   = 15
  integer, parameter :: CH4_3D_Diffusion_Flame  = 16
  integer, parameter :: N2_O2_Diffusion         = 17
  integer, parameter :: Compression_Corner      = 18
  integer, parameter :: CH4_Flame_Speed_Test    = 19
  integer, parameter :: THX_Phase_3             = 20
  integer, parameter :: Taylor_Couette_Flow     = 21
  integer, parameter :: Freestream_TCF          = 22
  integer, parameter :: Sod_Shock_Tube          = 23
  integer, parameter :: Density_Pulse_Test      = 24
  integer, parameter :: Pipe_Channel_Flow       = 25
  integer, parameter :: Shu_Osher_Test          = 26
  integer, parameter :: Advection_JiangShu      = 27
  integer, parameter :: UIUC_Scramjet           = 70
  integer, parameter :: C2H4_Flame_Speed_Test    = 71
  integer, parameter :: Advec_Problems(1:2) = [ Density_Transport, &
                                                Entropy_Transport ]
  integer, parameter :: Vortex_Problems(1:3) = [ Shu_Vortex, &
                                                 Vortex_Transport, &
                                                 Shu_MultiVortex ]
  integer, parameter :: Transport_Problems(1:5) = [ Vortex_Problems, &
                                                    Advec_Problems ]
  integer, parameter :: Taylor_Problems(1:3) = [ Taylor_Green_Vortex, &
                                                 Taylor_Couette_Flow, &
                                                -Taylor_Couette_Flow ]
  integer, parameter :: Reacting_Flows(1:7) = [Reacting_Backstep, &
                                               H2_2D_Diffusion_Flame, &
                                               CH4_3D_Diffusion_Flame, &
                                               Compression_Corner, &
                                               CH4_Flame_Speed_Test, &
                                               UIUC_Scramjet, &
                                               C2H4_Flame_Speed_Test]
  integer, parameter :: Dual_Stream_Mixing_Flows(1:2) = &
                                          [H2_2D_Diffusion_Flame, &
                                           N2_O2_Diffusion]
  !
  ! Grouping the different sets of points with
  ! values that match across geometries
  !
  ! Equal distant points
  integer, parameter :: Equi_Distant_Points     = 0
  !
  ! Interior-only points (e.g. Gauss points)
  integer, parameter :: Legendre_Gauss          = 1
  integer, parameter :: WilliamsShunn_TriPoints = 1
  integer, parameter :: ShunnHam_TetPoints      = 1
  !
  ! Preferred sets that include points on the faces (e.g. Lobatto points)
  integer, parameter :: Legendre_Gauss_Lobatto = 2
  integer, parameter :: Barycentric_Points     = 2
  !
  ! Non-preferred sets that include points on the faces (e.g. alpha-opt points)
  integer, parameter :: AlphaOptimized_Points = -2
  !
  ! Other sets that arent used anywhere yet
  integer, parameter :: Chebyshev_Gauss         = 1003
  integer, parameter :: Chebyshev_Gauss_Lobatto = 1004
  !
  ! Arrays grouping different sets of points
  integer, parameter :: Legendre_Nodes(1:2) = [Legendre_Gauss, &
                                               Legendre_Gauss_Lobatto]
  integer, parameter :: Gauss_Lobatto(1:2) = [Legendre_Gauss_Lobatto, &
                                              Chebyshev_Gauss_Lobatto]
  integer, parameter :: TriTet_Lobatto(1:2) = [Barycentric_Points, &
                                               AlphaOptimized_Points]
!!!
!!! Location of 1D solution points
!!integer, parameter :: Equi_Distant_Points     = 0
!!integer, parameter :: Legendre_Gauss          = 1
!!integer, parameter :: Legendre_Gauss_Lobatto  = 2
!!integer, parameter :: Legendre_Nodes(1:2) = [Legendre_Gauss, &
!!                                             Legendre_Gauss_Lobatto]
!!integer, parameter :: Chebyshev_Gauss         = 3
!!integer, parameter :: Chebyshev_Gauss_Lobatto = 4
!!integer, parameter :: Gauss_Lobatto(1:2) = [Legendre_Gauss_Lobatto, &
!!                                            Chebyshev_Gauss_Lobatto]
!!!
!!! Method for creating triangle solution points
!!integer, parameter :: AlphaOptimized_TriPoints = 1
!!integer, parameter :: Barycentric_TriPoints = 2
!!integer, parameter :: WilliamsShunn_TriPoints = 3
!!!
!!! Method for creating tetrahedra solution points
!!integer, parameter :: AlphaOptimized_TetPoints = 1
!!integer, parameter :: Barycentric_TetPoints = 2
!!integer, parameter :: ShunnHam_TetPoints = 3
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
  integer, parameter :: invs_flux_AUSM_Plus    = 5
  integer, parameter :: invs_flux_Rusanov      = 6
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
  integer, parameter :: OERK       = 6
  integer, parameter :: EPERK      = 7
  !
  ! OERK Parameters
  integer, parameter :: OERK_5Stage_3rdOrder  = 53
  integer, parameter :: OERK_11Stage_3rdOrder = 113
  integer, parameter :: OERK_6Stage_4thOrder  = 64
  integer, parameter :: OERK_15Stage_4thOrder = 154
  integer, parameter :: OERK_8Stage_5thOrder  = 85
  integer, parameter :: OERK_16Stage_5thOrder = 165
  integer, parameter :: OERK_Schemes(1:6) = [OERK_5Stage_3rdOrder, &
                                             OERK_11Stage_3rdOrder, &
                                             OERK_6Stage_4thOrder, &
                                             OERK_15Stage_4thOrder, &
                                             OERK_8Stage_5thOrder, &
                                             OERK_16Stage_5thOrder]
  !
  ! Embedded-Pair ERK Parameters
  integer, parameter :: EPERK_DormandPrince_5thOrder = 54
  integer, parameter :: EPERK_DormandPrince_8thOrder = 853
  integer, parameter :: OEPERK_3rdOrder              = 532
  integer, parameter :: OEPERK_4thOrder              = 943
  integer, parameter :: OEPERK_5thOrder              = 1054
  integer, parameter :: EPERK_DP_Schemes(1:2) = [EPERK_DormandPrince_5thOrder, &
                                                 EPERK_DormandPrince_8thOrder]
  integer, parameter :: OEPERK_Schemes(1:3) = [OEPERK_3rdOrder, &
                                               OEPERK_4thOrder, &
                                               OEPERK_5thOrder]
  integer, parameter :: EPERK_Schemes(1:5) = [EPERK_DP_Schemes(:), &
                                              OEPERK_Schemes(:)]
  !
  ! All valid Runge-Kutta schemes (to verify input is valid)
  integer, parameter :: Valid_RK_Schemes(1:18) = [Classic_RK, &
                                                  TVD2_RK, TVD3_RK, &
                                                  TVD4_RK, CK4_RK, &
                                                  OERK, EPERK, &
                                                  OERK_Schemes(:), &
                                                  EPERK_Schemes(:)]
  !
  ! =========================
  ! Choices for Filter_Option
  ! =========================
  !
  !> Turns the global filter off.
  integer, parameter :: Filter_Option__Off = 0
  !
  !> Turns the global filter on, but filters quad/hex cells
  !> independently in each direction.
  !> If the 1d filter for P=3 provides values of [1,a,b,c]
  !> for modes [0:P] with 1>=a>=b>=c>=0, then for a quad
  !> the 2d filter is the tensor product of the 1d filters.
  !>         | 1   | a   | b   | c
  !>         | x^0 | x^1 | x^2 | x^3
  !>   -----------------------------
  !>   1 y^0 | 1*1 | a*1 | b*1 | c*1
  !>   a y^1 | 1*a | a*a | b*a | c*a
  !>   b y^2 | 1*b | a*b | b*b | c*b
  !>   c y^3 | 1*c | a*c | b*c | c*c
  integer, parameter :: Filter_Option__TensProd = 1
  !
  !> Turns the global filter on, but filters quad/hex cells
  !> by the sum of the modal exponents rather than cell order.
  !> If the 1d filter for P=3 provides values of [1,a,b,c,d,e,f]
  !> for modes [0:2P] with 1>=a>=b>=c>d>=e>f>==0, then for a quad
  !> the 2d filter varies along the diagonal based on the sum of the
  !> modal exponent (m+n) of x^m*y^n.
  !>         | 1   | a   | b   | c
  !>         | x^0 | x^1 | x^2 | x^3
  !>   -----------------------------
  !>   1 y^0 | 1   | a   | b   | c  
  !>   a y^1 | a   | b   | c   | d  
  !>   b y^2 | b   | c   | d   | e  
  !>   c y^3 | c   | d   | e   | f  
  !> This option is hidden from the user documentation.
  integer, parameter :: Filter_Option__ModalSum = 2
  !
  ! =================================
  ! Choices for Pyramid_Filter_Option
  ! =================================
  !
  !> Keep the cutoff order (Pc) the same and let the cutoff fraction (etac)
  !> vary for all pyramid levels (i=0:P). Local order, P_i=i.
  !>   Pc = etac*real(P).
  !>   etac_i = etac * real(P)/real(max(1,i)) = Pc / real(max(1,i)).
  !>   Pc_i = etac_i * real(P_i).
  !> Must use min(etac_i,1.0) to stay in valid range [0:1] for small i.
  !> For order P=4, etac=0.5:
  !>   i | P_i | Pc_i | etac_i
  !>   -------------------------------
  !>   4 | 4   | 2.0  | 2/4 = 1/2
  !>   3 | 3   | 2.0  | 2/3 = 2/3
  !>   2 | 2   | 2.0  | 2/2 = 1
  !>   1 | 1   | 2.0  | 2/1 = min(1,2)
  !>   0 | 0   | 2.0  | 2/1 = min(1,2)
  !> The base mode is not filtered to preserve the cell average value.
  integer, parameter :: Pyramid_Filter_Option__1 = 1
  !
  !> Keep the cutoff fraction (etac) the same and let the cutoff order (Pc)
  !> vary for all pyramid levels (i=0:P). Local order, P_i=i.
  !>   Pc = etac*real(i).
  !>   etac_i = etac.
  !>   Pc_i = etac_i * real(P_i).
  !> For order P=4, etac=0.5:
  !>   i | P_i | Pc_i | etac_i
  !>   -----------------------
  !>   4 | 4   | 2.0  | 0.5
  !>   3 | 3   | 1.5  | 0.5
  !>   2 | 2   | 1.0  | 0.5
  !>   1 | 1   | 0.5  | 0.5
  !>   0 | 0   | 0.0  | 0.5
  !> This results in a stronger filter compared to Pyramid_Filter_Option=1
  !> because it affects more modes.
  !> The base mode is not filtered to preserve the cell average value.
  integer, parameter :: Pyramid_Filter_Option__2 = 2
  !
  !> Apply the standard tensor product filter to pyramid cells.
  !> The base mode is not filtered to preserve the cell average value.
  integer, parameter :: Pyramid_Filter_Option__3 = 3
  !
  ! ============================
  ! Choices for Limiter_Option
  integer, parameter :: Limiter_Option__Off                 = 0
  integer, parameter :: Limiter_Option__Filter              = 1
  integer, parameter :: Limiter_Option__FilterScaled        = 2
  integer, parameter :: Limiter_Option__MLP                 = 3
  integer, parameter :: Limiter_Option__ProjectionRecursive = 9001
  integer, parameter :: Limiter_Option__ProjectionScaled    = 9002
  !
  ! Choices for Marker_Option
  integer, parameter :: Marker_Option__Off        = 0
  integer, parameter :: Marker_Option__Resolution = 1
  integer, parameter :: Marker_Option__Jump       = 2
  integer, parameter :: Marker_Option__MLP        = 3
  integer, parameter :: Marker_Option__Ducros     = 4
  !
  ! Choices for Marker_Variable
  integer, parameter :: Marker_Variable__Off    =  0
  integer, parameter :: Marker_Variable__Ducros = -2
  integer, parameter :: Marker_Variable__npres  = -1
  integer, parameter :: Marker_Variable__nec    =  1
  integer, parameter :: Marker_Variable__nmx    =  2
  integer, parameter :: Marker_Variable__nmy    =  3
  integer, parameter :: Marker_Variable__nmz    =  4
  integer, parameter :: Marker_Variable__nee    =  5
  integer, parameter :: Marker_Variable__all    =  6 ! nec,nmx,nmy,nmz,nee
  !
  ! Choices for thermofix%method
  integer, parameter :: TFM_OFF   = 0  !< No check/fix.
  integer, parameter :: TFM_LIMIT = 1  !< Limit HO modes via sponge filter.
  integer, parameter :: TFM_SCALE = 2  !< Scale HO modes by a user factor.
  integer, parameter :: TFM_LIMI2 = 3  !< Use Positivity Preservation
                                 !(Sheshadri Stanford Thesis, S0 -> -\infty)
  !> Valid choices for thermofix%method.
  integer, parameter :: Valid_TFM(1:4) = [TFM_OFF,TFM_LIMIT,TFM_SCALE,TFM_LIMI2]
  !
  !> Choices for type of method to modify the
  !> solution using with entropy limiter
  integer, parameter :: ENTLIM_OPT__OFF          = 0
  integer, parameter :: ENTLIM_OPT__Flat         = 1
  integer, parameter :: ENTLIM_OPT__ScaledFilter = 2
  !
  !> Choices for type of interpolant function to use to compute
  !> the damping ratio (badness) needed for the entropy limiter
  integer, parameter :: ENTLIM_FUN__LinearLvIhme    = 0
  integer, parameter :: ENTLIM_FUN__Linear          = 1
  integer, parameter :: ENTLIM_FUN__DensityWeighted = 2
  integer, parameter :: ENTLIM_FUN__LogScaled       = 3
  integer, parameter :: ENTLIM_FUN__LogScaledPade   = 4
  !
  !> Choices for positivity limiter options
  integer, parameter :: POSLIM_OPT__OFF = 0
  integer, parameter :: POSLIM_OPT__ON  = 1
  !
  ! Choices for ignition method
  integer, parameter :: IGN_OFF        = 0  ! No ignition source term
  ! IGN_DELTA_TEMP : delta temperature based source term added to energy eqn
  integer, parameter :: IGN_DELTA_TEMP = 1
  ! IGN_EREF_SPARK : percent of ref rho*E is added as source to energy eqn
  integer, parameter :: IGN_EREF_SPARK = 2
  ! IGN_REACT_TEMP : specified temperature used to compute reaction rates
  integer, parameter :: IGN_REACT_TEMP = 3
  ! IGN_REACT_PRES : react pressure used to scale species concentrations
  integer, parameter :: IGN_REACT_PRES = 4
  !
  ! Valid choices for ignition method
  integer, parameter :: Valid_IGN(1:5) = [IGN_OFF, &
                                          IGN_DELTA_TEMP, &
                                          IGN_EREF_SPARK, &
                                          IGN_REACT_TEMP, &
                                          IGN_REACT_PRES]
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
  ! Grid/Solution Types
  integer, parameter :: Unknown_Format = -99
  integer, parameter :: Internally_Generated_Grid = -1
  integer, parameter :: Restart_Format = 0
  integer, parameter :: Gambit_Neutral = 1
  integer, parameter :: Gmsh_Format = 2
  integer, parameter :: Plot3D_Format = 3
  integer, parameter :: CGNS_Format = 4
  !
  ! Parameterizations for Plot3D Agglomeration
  integer, parameter :: PmetAgglom__OFF          = 0
  integer, parameter :: PmetAgglom__Cartesian    = 1
  integer, parameter :: PmetAgglom__Axisymmetric = 2
  ! ANOTHER POSSIBLE SET OF AGGLOMERATION PARAMETERS TO ALLOW
  ! FOR AXISYMMETRIC GRIDS TO HAVE DIFFERENT AXES OF SYMMETRY
 !integer, parameter :: PmetAgglom__AxisymmX = 3
 !integer, parameter :: PmetAgglom__AxisymmY = 4
 !integer, parameter :: PmetAgglom__AxisymmZ = 5
  integer, parameter :: PmetAgglom__ON(1:2) = [PmetAgglom__Cartesian, &
                                               PmetAgglom__Axisymmetric]
  !
  ! CGNS constants
  integer, parameter :: reopen_CGNS = 2
  integer, parameter :: open_CGNS   =  1
  integer, parameter :: add_to_CGNS =  0
  integer, parameter :: close_CGNS  = -1
  integer, parameter :: dump_CGNS   = -2
  !
  ! Channel flow constants
  real(wp), parameter :: channel_grid_Lx = 0.036576_wp ! 0.12(ft)*0.3048(m/ft)
  real(wp), parameter :: channel_grid_Ly = 0.006096_wp ! 0.02(ft)*0.3048(m/ft)
  real(wp), parameter :: channel_grid_Lz = 0.018288_wp ! 0.06(ft)*0.3048(m/ft)
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
  character(len=15),  protected, save :: Geom_Full_Name(Geom_Min:Geom_Max)
  character(len=4),   protected, save :: Geom_Short_Name(Geom_Min:Geom_Max)
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
  integer, parameter :: Xmin3D(1:4) = [0,4,7,3] ! [1,5,8,4]
  integer, parameter :: Xmax3D(1:4) = [1,2,6,5] ! [2,3,7,6]
  integer, parameter :: Ymin3D(1:4) = [0,1,5,4] ! [1,2,6,5]
  integer, parameter :: Ymax3D(1:4) = [3,7,6,2] ! [4,8,7,3]
  integer, parameter :: Zmin3D(1:4) = [0,3,2,1] ! [1,4,3,2]
  integer, parameter :: Zmax3D(1:4) = [4,5,6,7] ! [5,6,7,8]
  !
  integer, parameter :: HostFace2D(1:6) = [Xmin2D_Face,Xmax2D_Face, &
                                           Ymin2D_Face,Ymax2D_Face, &
                                           Zmin2D_Face,Zmax2D_Face]
  ! HostFace3D is used to traditional values of 
  ! [imin,imax,jmin,jmax,kmin,kmax] = [1,2,3,4,5,6] to their
  ! corresponding GFR side/face values of
  ! [zmin,zmax,xmax,xmin,ymin,ymax] = [1,2,3,4,5,6] or
  ! HostFace3D[imin,imax,jmin,jmax,kmin,kmax] = [4,3,5,6,1,2] which is also
  ! HostFace3D[1,2,3,4,5,6] = [Xmin3D_Face,Xmax3D_Face, &
  !                            Ymin3D_Face,Ymax3D_Face, &
  !                            Zmin3D_Face,Zmax3D_Face]
  integer, parameter :: HostFace3D(1:6) = [Xmin3D_Face,Xmax3D_Face, &
                                           Ymin3D_Face,Ymax3D_Face, &
                                           Zmin3D_Face,Zmax3D_Face]
  !
  real(wp), protected, save :: ref_edge_nodes_xyz(1:1,1:2)
  real(wp), protected, save :: ref_tria_nodes_xyz(1:2,1:3)
  real(wp), protected, save :: ref_quad_nodes_xyz(1:2,1:4)
  real(wp), protected, save :: ref_tetr_nodes_xyz(1:3,1:4)
  real(wp), protected, save :: ref_pyra_nodes_xyz(1:3,1:5)
  real(wp), protected, save :: ref_pris_nodes_xyz(1:3,1:6)
  real(wp), protected, save :: ref_hexa_nodes_xyz(1:3,1:8)
  real(wp), protected, save :: ref_elem_nodes_xyz(1:3,1:8,Geom_Min:Geom_Max)
  !
  ! Minimum/Maximum Order Allowed
  integer, parameter :: Order_Min = 0
  integer, parameter :: Order_Max = 25
  !
  integer, save :: exclusive_element = -1
  !
  ! tria_is_biunit : logical that changes the location of the reference
  !                  nodes in the definition of the reference triangle
  !         = true  -> defined by the three nodes (-1,-1) - ( 1,-1) - (-1, 1)
  !         = false -> defined by the three nodes ( 0, 0) - ( 1, 0) - ( 0, 1)
  !
  logical(lk), save :: tria_is_biunit = true
 !logical(lk), parameter :: tria_is_biunit = true
 !logical(lk), parameter :: tria_is_biunit = fals
  !
  ! CL : coefficient that gives (2.0 / length of the 1D solution points).
  !      This is needed because of the difference between how the reference
  !      triangle and reference quad are defined.
  !      = 2.0 for triangles if tria_is_biunit is false
  !      = 1.0 for everything else
  !
  real(wp), protected, save :: cl(Geom_Min:Geom_Max)
  !$ACC DECLARE CREATE(cl)
  !
  ! General character variable for storing error messages
  !
  character(len=1000), save :: error_message
  character(len=1000), save :: emsg
  character(len=1000), save :: cerr
  !
  ! Variables for storing information about the HPC system
  ! that are defined by GFR environent variables.
  !
  character(len=:), allocatable, protected, save :: GFR_PrgEnv
  character(len=:), allocatable, protected, save :: GFR_Scheduler
  character(len=:), allocatable, protected, save :: GFR_Scheduler_Path
  character(len=:), allocatable, protected, save :: GFR_Scheduler_JobInfo_Cmnd
  character(len=:), allocatable, protected, save :: GFR_Scheduler_JobInfo_Flags
  integer,                       protected, save :: GFR_JobInfo_UsedTime
  integer,                       protected, save :: GFR_JobInfo_MaxTime
  character(len=:), allocatable, protected, save :: GFR_JobInfo_JobID
  !
  ! OpenMP Public Scalars
  !
  integer, save :: tid = 0
  integer, save :: nthreads = 1
  integer, save :: chunk_size = 1
  !
  character(len=*), parameter :: thread_report = &
                                    "(/,2x,'number of threads (',a,') = ',i0,/)"
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
  integer(int_mpi), parameter :: dump_stack = int( 7 , kind=kind(dump_stack) )
  integer(int_mpi), parameter :: abort      = int( 8 , kind=kind(abort     ) )
  integer(int_mpi), parameter :: stop_mpi   = int( 9 , kind=kind(stop_mpi  ) )
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
  _MPI_DATA_TYPE_, protected, save :: mpi_2flttyp
  ! half working real precision
  _MPI_DATA_TYPE_, protected, save :: mpi_hlftyp
  ! double working real precision
  _MPI_DATA_TYPE_, protected, save :: mpi_dbltyp
  ! working integer precision
  _MPI_DATA_TYPE_, protected, save :: mpi_inttyp
  !
  integer, public, save :: my_dbg_out = iout
  !
  character(len=*), parameter :: dbg_mpi_frmt = &
                              '("CPU-",i0," : ",a," at line ",i0," in file ",a)'
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
  logical(lk), public, save :: debug_timer_disabled = fals
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
  ! From Functions/reallocate.f90
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
  ! From Functions/init_array.f90
  !
  ! Versions of init_array with size inputs being 32-bit integers
  private :: &
      init_array_s4_i1_rank1, init_array_s4_i1_rank2, init_array_s4_i1_rank3, &
      init_array_s4_i2_rank1, init_array_s4_i2_rank2, init_array_s4_i2_rank3, &
      init_array_s4_i4_rank1, init_array_s4_i4_rank2, init_array_s4_i4_rank3, &
      init_array_s4_i8_rank1, init_array_s4_i8_rank2, init_array_s4_i8_rank3, &
      init_array_s4_r4_rank1, init_array_s4_r4_rank2, init_array_s4_r4_rank3, &
      init_array_s4_r8_rank1, init_array_s4_r8_rank2, init_array_s4_r8_rank3, &
      init_array_s4_r16_rank1,init_array_s4_r16_rank2,init_array_s4_r16_rank3, &
      init_array_s4_lk_rank1, init_array_s4_lk_rank2, init_array_s4_lk_rank3, &
      init_array_s4_l4_rank1, init_array_s4_l4_rank2, init_array_s4_l4_rank3
  !
  ! Versions of init_array with size inputs being 64-bit integers
  private :: &
      init_array_s8_i1_rank1, init_array_s8_i1_rank2, init_array_s8_i1_rank3, &
      init_array_s8_i2_rank1, init_array_s8_i2_rank2, init_array_s8_i2_rank3, &
      init_array_s8_i4_rank1, init_array_s8_i4_rank2, init_array_s8_i4_rank3, &
      init_array_s8_i8_rank1, init_array_s8_i8_rank2, init_array_s8_i8_rank3, &
      init_array_s8_r4_rank1, init_array_s8_r4_rank2, init_array_s8_r4_rank3, &
      init_array_s8_r8_rank1, init_array_s8_r8_rank2, init_array_s8_r8_rank3, &
      init_array_s8_r16_rank1,init_array_s8_r16_rank2,init_array_s8_r16_rank3, &
      init_array_s8_lk_rank1, init_array_s8_lk_rank2, init_array_s8_lk_rank3, &
      init_array_s8_l4_rank1, init_array_s8_l4_rank2, init_array_s8_l4_rank3
  !
  interface init_array
    module procedure &
       init_array_s4_i1_rank1, init_array_s4_i1_rank2, init_array_s4_i1_rank3, &
       init_array_s4_i2_rank1, init_array_s4_i2_rank2, init_array_s4_i2_rank3, &
       init_array_s4_i4_rank1, init_array_s4_i4_rank2, init_array_s4_i4_rank3, &
       init_array_s4_i8_rank1, init_array_s4_i8_rank2, init_array_s4_i8_rank3, &
       init_array_s4_r4_rank1, init_array_s4_r4_rank2, init_array_s4_r4_rank3, &
       init_array_s4_r8_rank1, init_array_s4_r8_rank2, init_array_s4_r8_rank3, &
       init_array_s4_lk_rank1, init_array_s4_lk_rank2, init_array_s4_lk_rank3, &
       init_array_s4_l4_rank1, init_array_s4_l4_rank2, init_array_s4_l4_rank3, &
       init_array_s8_i1_rank1, init_array_s8_i1_rank2, init_array_s8_i1_rank3, &
       init_array_s8_i2_rank1, init_array_s8_i2_rank2, init_array_s8_i2_rank3, &
       init_array_s8_i4_rank1, init_array_s8_i4_rank2, init_array_s8_i4_rank3, &
       init_array_s8_i8_rank1, init_array_s8_i8_rank2, init_array_s8_i8_rank3, &
       init_array_s8_r4_rank1, init_array_s8_r4_rank2, init_array_s8_r4_rank3, &
       init_array_s8_r8_rank1, init_array_s8_r8_rank2, init_array_s8_r8_rank3, &
       init_array_s8_lk_rank1, init_array_s8_lk_rank2, init_array_s8_lk_rank3, &
       init_array_s8_l4_rank1, init_array_s8_l4_rank2, init_array_s8_l4_rank3
#ifdef ENABLE_QP
    module procedure &
    init_array_s4_r16_rank1, init_array_s4_r16_rank2, init_array_s4_r16_rank3, &
    init_array_s8_r16_rank1, init_array_s8_r16_rank2, init_array_s8_r16_rank3
#endif
  end interface init_array
  !
  ! From Functions/spread.f90
  !
  private :: &
      spread_i4_rank0, spread_i4_rank1, spread_i4_rank2, spread_i4_rank3, &
      spread_i8_rank0, spread_i8_rank1, spread_i8_rank2, spread_i8_rank3, &
      spread_r4_rank0, spread_r4_rank1, spread_r4_rank2, spread_r4_rank3, &
      spread_r8_rank0, spread_r8_rank1, spread_r8_rank2, spread_r8_rank3, &
      spread_r16_rank0,spread_r16_rank1,spread_r16_rank2,spread_r16_rank3,&
      spread_lk_rank0, spread_lk_rank1, spread_lk_rank2, spread_lk_rank3, &
      spread_l4_rank0, spread_l4_rank1, spread_l4_rank2, spread_l4_rank3
  !
#ifdef SPECIAL_FOR_CRAY
  interface spread
    module procedure &
       spread_i4_rank0, spread_i4_rank1, spread_i4_rank2, spread_i4_rank3, &
       spread_i8_rank0, spread_i8_rank1, spread_i8_rank2, spread_i8_rank3, &
       spread_r4_rank0, spread_r4_rank1, spread_r4_rank2, spread_r4_rank3, &
       spread_r8_rank0, spread_r8_rank1, spread_r8_rank2, spread_r8_rank3, &
       spread_lk_rank0, spread_lk_rank1, spread_lk_rank2, spread_lk_rank3, &
       spread_l4_rank0, spread_l4_rank1, spread_l4_rank2, spread_l4_rank3
#ifdef ENABLE_QP
    module procedure &
       spread_r16_rank0, spread_r16_rank1, spread_r16_rank2, spread_r16_rank3
#endif
  end interface spread
#else
  intrinsic :: spread
#endif
  !
  ! From Functions/all_are_equal.f90
  !
  private :: all_are_equal_i4, all_are_equal_i8, &
             all_are_equal_r4, all_are_equal_r8, all_are_equal_r16
  !
  interface all_are_equal
    module procedure all_are_equal_i4, all_are_equal_i8, &
                     all_are_equal_r4, all_are_equal_r8
#ifdef ENABLE_QP
    module procedure all_are_equal_r16
#endif
  end interface all_are_equal
  !
  ! From Functions/all_are_equal.f90
  !
  private :: all_not_equal_i4, all_not_equal_i8, &
             all_not_equal_r4, all_not_equal_r8, all_not_equal_r16
  !
  interface all_not_equal
    module procedure all_not_equal_i4, all_not_equal_i8, &
                     all_not_equal_r4, all_not_equal_r8
#ifdef ENABLE_QP
    module procedure all_not_equal_r16
#endif
  end interface all_not_equal
  !
  ! From Functions/cross_product.f90
  !
  private :: cross_product_SP,cross_product_DP,cross_product_QP
  !
  interface cross_product
    module procedure cross_product_SP
    module procedure cross_product_DP
#ifdef ENABLE_QP
    module procedure cross_product_QP
#endif
  end interface cross_product
  !
  ! From Functions/deltafun.f90
  !
  private :: deltafun_DP_LDK,deltafun_DP_LK
  !
  interface deltafun
    module procedure deltafun_DP_LDK
    module procedure deltafun_DP_LK
  end interface deltafun
  !
  ! From Functions/deltafun.f90
  !
  private :: deltafun_QP_LDK,deltafun_QP_LK
  !
  interface deltafun_QP
    module procedure deltafun_QP_LDK
    module procedure deltafun_QP_LK
  end interface deltafun_QP
  !
  ! From Functions/norms.f90
  !
  private :: selfdot_SP,selfdot_DP,selfdot_QP
  !
  interface selfdot
    module procedure selfdot_SP
    module procedure selfdot_DP
#ifdef ENABLE_QP
    module procedure selfdot_QP
#endif
  end interface selfdot
  !
  ! From Functions/norms.f90
  !
  private :: rnorm2_SP,rnorm2_DP,rnorm2_QP
  !
#ifdef SPECIAL_FOR_NVIDIA
  interface norm2
    module procedure rnorm2_SP
    module procedure rnorm2_DP
#ifdef ENABLE_QP
    module procedure rnorm2_QP
#endif
  end interface norm2
#endif
  !
  ! From Functions/norms.f90
  !
  private :: unit_vector_SP,unit_vector_DP,unit_vector_QP
  !
  interface unit_vector
    module procedure unit_vector_SP
    module procedure unit_vector_DP
#ifdef ENABLE_QP
    module procedure unit_vector_QP
#endif
  end interface unit_vector
  !
  ! From Functions/reciprocal.f90
  !
  private :: reciprocal_SP,reciprocal_DP,reciprocal_QP
  !
  interface reciprocal
    module procedure reciprocal_SP
    module procedure reciprocal_DP
#ifdef ENABLE_QP
    module procedure reciprocal_QP
#endif
  end interface reciprocal
  !
  ! From Functions/chop.f90
  !
  private :: chop_WP,chop_SP,chop_DP,chop_QP
  !
  interface chop
    module procedure chop_WP
    module procedure chop_SP
    module procedure chop_DP
#ifdef ENABLE_QP
    module procedure chop_QP
#endif
  end interface chop
  !
  ! From Functions/face_normal.f90
  !
  private :: face_normal_SP,face_normal_DP,face_normal_QP
  !
  interface face_normal
    module procedure face_normal_SP
    module procedure face_normal_DP
#ifdef ENABLE_QP
    module procedure face_normal_QP
#endif
  end interface face_normal
  !
  ! From Functions/jacobians.f90
  !
  private :: compute_jacobian_SP,compute_jacobian_DP,compute_jacobian_QP
  !
  interface compute_jacobian
    module procedure compute_jacobian_SP
    module procedure compute_jacobian_DP
#ifdef ENABLE_QP
    module procedure compute_jacobian_QP
#endif
  end interface compute_jacobian
  !
  ! From Functions/reverse_metrics.f90
  !
  private :: reverse_metrics_SP,reverse_metrics_DP,reverse_metrics_QP
  !
  interface reverse_metrics
    module procedure reverse_metrics_SP
    module procedure reverse_metrics_DP
#ifdef ENABLE_QP
    module procedure reverse_metrics_QP
#endif
  end interface reverse_metrics
  !
  ! From Functions/matrix_ops.f90
  !
  private :: invert_matrix_SP,invert_matrix_DP,invert_matrix_QP
  !
  interface invert_matrix
    module procedure invert_matrix_SP
    module procedure invert_matrix_DP
#ifdef ENABLE_QP
    module procedure invert_matrix_QP
#endif
  end interface invert_matrix
  !
  ! From Functions/matrix_ops.f90
  !
  private :: ludcmp_SP,ludcmp_DP,ludcmp_QP
  !
  interface ludcmp
    module procedure ludcmp_SP
    module procedure ludcmp_DP
#ifdef ENABLE_QP
    module procedure ludcmp_QP
#endif
  end interface ludcmp
  !
  ! From Functions/matrix_ops.f90
  !
  private :: lubksb_SP,lubksb_DP,lubksb_QP
  !
  interface lubksb
    module procedure lubksb_SP
    module procedure lubksb_DP
#ifdef ENABLE_QP
    module procedure lubksb_QP
#endif
  end interface lubksb
  !
  ! From Functions/matrix_ops.f90
  !
  private :: determinant_SP,determinant_DP,determinant_QP
  !
  interface determinant
    module procedure determinant_SP
    module procedure determinant_DP
#ifdef ENABLE_QP
    module procedure determinant_QP
#endif
  end interface determinant
  !
  ! From Functions/random.f90
  !
  private :: random_i4,random_sp,random_dp,random_qp
  !
  interface random
    module procedure random_i4
    module procedure random_sp
    module procedure random_dp
#ifdef ENABLE_QP
    module procedure random_qp
#endif
  end interface random
  !
  ! From Functions/matrix_ops.f90
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
#ifdef ENABLE_QP
    module procedure identity_matrix_dim_QP
    module procedure identity_matrix_mat_QP
#endif
  end interface identity_matrix
  !
  ! From Functions/intseq.f90
  !
  private :: intseq_no_beg_i1, intseq_no_interval_i1
  private :: intseq_rev_dir_i1,intseq_with_interval_i1
  private :: intseq_no_beg_i2, intseq_no_interval_i2
  private :: intseq_rev_dir_i2,intseq_with_interval_i2
  private :: intseq_no_beg_i4, intseq_no_interval_i4
  private :: intseq_rev_dir_i4,intseq_with_interval_i4
  private :: intseq_no_beg_i8, intseq_no_interval_i8
  private :: intseq_rev_dir_i8,intseq_with_interval_i8
  !
  interface intseq
    module procedure &
      intseq_no_beg_i1, intseq_no_interval_i1,   &
      intseq_rev_dir_i1,intseq_with_interval_i1, &
      intseq_no_beg_i2, intseq_no_interval_i2,   &
      intseq_rev_dir_i2,intseq_with_interval_i2, &
      intseq_no_beg_i4, intseq_no_interval_i4,   &
      intseq_rev_dir_i4,intseq_with_interval_i4, &
      intseq_no_beg_i8, intseq_no_interval_i8,   &
      intseq_rev_dir_i8,intseq_with_interval_i8
  end interface intseq
  !
contains
!
!###############################################################################
!
subroutine initialize_gfr_environment
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_gfr_environment"
  !
continue
  !
  ! Save the default IEEE floating point status
  !
  call get_default_ieee_status
  !
  ! Initialize the system and job information variables using the
  ! GFR environment variables
  !
  call read_gfr_environment_variables
  !
  ! Initialize all the random geometry constants that
  ! must be initialized at run time
  !
  call initialize_runtime_parameters
  !
  ! Initialize the MPI environment
  !
  call initialize_mpi_environment
  !
end subroutine initialize_gfr_environment
!
!###############################################################################
!
subroutine read_gfr_environment_variables
  !
  !.. Local Scalars ..
  integer     :: n,ierr
  logical(lk) :: file_exists
  !
  character(len=16)    :: tmpfile
  character(len=10000) :: cvalue
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: cmd
  character(len=:), allocatable :: command
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_gfr_environment_variables"
  !
continue
  !
  ! Get the programming/compiler environment used to build GFR
  !
  call get_environment_variable("GFR_PrgEnv",cvalue,n,status=ierr)
  !
  GFR_PrgEnv = trim(adjustl(cvalue)) ! F2003 AUTO-REALLOC
  !
  ! Get the name of the batch scheduling software (BSS) used
  ! by the current system
  !
  call get_environment_variable("GFR_Scheduler",cvalue,n,status=ierr)
  !
  GFR_Scheduler = trim(adjustl(cvalue)) ! F2003 AUTO-REALLOC
  !
  ! Get the path to the BSS used by the current system
  !
  call get_environment_variable("GFR_Scheduler_Path",cvalue,n,status=ierr)
  !
  GFR_Scheduler_Path = trim(adjustl(cvalue)) ! F2003 AUTO-REALLOC
  !
  ! Set the name for the job-status command based on the BSS being used
  !
  select case ( uppercase(GFR_Scheduler) )
    !
    case ("PBS")
      !
      GFR_Scheduler_JobInfo_Cmnd  = "qstat" ! F2003 AUTO-REALLOC
      GFR_Scheduler_JobInfo_Flags = "-f"    ! F2003 AUTO-REALLOC
      !
      ! Get the PBS job ID
      call get_environment_variable("PBS_JOBID",cvalue,n,status=ierr)
      GFR_JobInfo_JobID = trim(adjustl(cvalue)) ! F2003 AUTO-REALLOC
      !
    case ("SLURM")
      !
      GFR_Scheduler_JobInfo_Cmnd  = "squeue" ! F2003 AUTO-REALLOC
      GFR_Scheduler_JobInfo_Flags = "-f"     ! F2003 AUTO-REALLOC
      !
      ! Get the Slurm job ID
      call get_environment_variable("SLURM_JOB_ID",cvalue,n,status=ierr)
      GFR_JobInfo_JobID = trim(adjustl(cvalue)) ! F2003 AUTO-REALLOC
      !
    case default
      !
      GFR_Scheduler = ""               ! F2003 AUTO-REALLOC
      GFR_Scheduler_Path = ""          ! F2003 AUTO-REALLOC
      GFR_Scheduler_JobInfo_Cmnd  = "" ! F2003 AUTO-REALLOC
      GFR_Scheduler_JobInfo_Flags = "" ! F2003 AUTO-REALLOC
      cmd = ""                         ! F2003 AUTO-REALLOC
      !
  end select
  !
  if (.not. allocated(cmd)) then
    ! F2003 AUTO-REALLOC
    cmd = GFR_Scheduler_Path // "/bin/" // GFR_Scheduler_JobInfo_Cmnd
  end if
  !
  inquire (file=cmd,exist=file_exists)
  !
  write (tmpfile,'("tmp.",i8.8,".tmp")') random(99999999)
  !
  ! F2003 AUTO-REALLOC
  command = cmd // &
            " " // GFR_Scheduler_JobInfo_Flags // &
            " " // GFR_JobInfo_JobID // &
            " |& /usr/bin/sed -n -e '/walltime/p' >| " // tmpfile
  !
  if (.not. file_exists) then
    GFR_Scheduler = ""
  end if
  !
end subroutine read_gfr_environment_variables
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
  Geom_Full_Name(Geom_Node) = "Node"
  Geom_Full_Name(Geom_Edge) = "Edge"
  Geom_Full_Name(Geom_Tria) = "Triangle"
  Geom_Full_Name(Geom_Quad) = "Quadrilateral"
  Geom_Full_Name(Geom_Tetr) = "Tetrahedra"
  Geom_Full_Name(Geom_Pyra) = "Pyramid"
  Geom_Full_Name(Geom_Pris) = "Prism"
  Geom_Full_Name(Geom_Unknown) = "Unknown"
  Geom_Full_Name(Geom_Hexa) = "Hexahedra"
  !
  Geom_Short_Name(Geom_Node) = "node"
  Geom_Short_Name(Geom_Edge) = "edge"
  Geom_Short_Name(Geom_Tria) = "tria"
  Geom_Short_Name(Geom_Quad) = "quad"
  Geom_Short_Name(Geom_Tetr) = "tetr"
  Geom_Short_Name(Geom_Pyra) = "pyra"
  Geom_Short_Name(Geom_Pris) = "pris"
  Geom_Short_Name(Geom_Unknown) = "unkn"
  Geom_Short_Name(Geom_Hexa) = "hexa"
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
  side_nodes(1:3,1,Geom_Tetr) = [0,1,3]
  side_nodes(1:3,2,Geom_Tetr) = [1,2,3]
  side_nodes(1:3,3,Geom_Tetr) = [2,0,3]
  side_nodes(1:3,4,Geom_Tetr) = [0,2,1]
 !! #### OLD PRE-UNSTRUCTURED GFR
 !side_nodes(1:3,1,Geom_Tetr) = [1,2,3]
 !side_nodes(1:3,2,Geom_Tetr) = [2,0,3]
 !side_nodes(1:3,3,Geom_Tetr) = [3,0,1]
 !side_nodes(1:3,4,Geom_Tetr) = [0,2,1]
  !
  side_nodes(1:3,1,Geom_Pyra) = [0,1,4]
  side_nodes(1:3,2,Geom_Pyra) = [1,2,4]
  side_nodes(1:3,3,Geom_Pyra) = [2,3,4]
  side_nodes(1:3,4,Geom_Pyra) = [3,0,4]
  side_nodes(1:4,5,Geom_Pyra) = [0,3,2,1]
 !! #### OLD PRE-UNSTRUCTURED GFR
 !side_nodes(1:3,1,Geom_Pyra) = [1,2,4]
 !side_nodes(1:3,2,Geom_Pyra) = [2,3,4]
 !side_nodes(1:3,3,Geom_Pyra) = [3,0,4]
 !side_nodes(1:3,4,Geom_Pyra) = [0,1,4]
 !side_nodes(1:4,5,Geom_Pyra) = [0,3,2,1]
  !
  side_nodes(1:3,1,Geom_Pris) = [0,2,1]
  side_nodes(1:3,2,Geom_Pris) = [3,4,5]
  side_nodes(1:4,3,Geom_Pris) = [0,1,4,3]
  side_nodes(1:4,4,Geom_Pris) = [1,2,5,4]
  side_nodes(1:4,5,Geom_Pris) = [0,3,5,2]
 !! #### OLD PRE-UNSTRUCTURED GFR
 !side_nodes(1:3,1,Geom_Pris) = [0,2,1]
 !side_nodes(1:3,2,Geom_Pris) = [3,4,5]
 !side_nodes(1:4,3,Geom_Pris) = [0,1,4,3]
 !side_nodes(1:4,4,Geom_Pris) = [1,2,5,4]
 !side_nodes(1:4,5,Geom_Pris) = [3,5,2,0]
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
  ! Define the node coordinates for each reference element
  !
  call define_ref_elem_node_coordinates
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
  !$ACC UPDATE DEVICE(cl)
  !
 !call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  !
end subroutine initialize_runtime_parameters
!
!###############################################################################
!
subroutine define_ref_elem_node_coordinates
  !
  ! Define the node coordinates for each reference element
  !
continue
  !
  ref_edge_nodes_xyz = reshape( [-one,one] , shape(ref_edge_nodes_xyz) )
  !
  ref_tria_nodes_xyz = reshape( [-one,-one, &
                                  one,-one, &
                                 -one, one] , shape(ref_tria_nodes_xyz) )
  if (.not. tria_is_biunit) ref_tria_nodes_xyz = half*(ref_tria_nodes_xyz+one)
  !
  ref_quad_nodes_xyz = reshape( [-one,-one, &
                                  one,-one, &
                                  one, one, &
                                 -one, one] , shape(ref_quad_nodes_xyz) )
  !
  ref_tetr_nodes_xyz = reshape( [-one,-one,-one, &
                                  one,-one,-one, &
                                 -one, one,-one, &
                                 -one,-one, one] , shape(ref_tetr_nodes_xyz) )
  !
  ref_pyra_nodes_xyz = reshape( [-one,-one,-one, &
                                  one,-one,-one, &
                                  one, one,-one, &
                                 -one, one,-one, &
                                 zero,zero, one] , shape(ref_pyra_nodes_xyz) )
  !
  ref_pris_nodes_xyz = reshape( [-one,-one,-one, &
                                  one,-one,-one, &
                                 -one, one,-one, &
                                 -one,-one, one, &
                                  one,-one, one, &
                                 -one, one, one] , shape(ref_pris_nodes_xyz) )
  !
  ref_hexa_nodes_xyz = reshape( [-one,-one,-one, &
                                  one,-one,-one, &
                                  one, one,-one, &
                                 -one, one,-one, &
                                 -one,-one, one, &
                                  one,-one, one, &
                                  one, one, one, &
                                 -one, one, one] , shape(ref_hexa_nodes_xyz) )
  !
  ! Initialize the parts of ref_elem_nodes_xyz that shouldnt be
  ! accessed to a number that will hopefully cause problems.
  !
  ref_elem_nodes_xyz = -huge(zero)
  !
  ! Gather the node coordinates for all reference element
  ! types to the combined array ref_elem_nodes_xyz
  !
  ref_elem_nodes_xyz(1:1,1:2,Geom_Edge) = ref_edge_nodes_xyz(1:1,1:2)
  ref_elem_nodes_xyz(1:2,1:3,Geom_Tria) = ref_tria_nodes_xyz(1:2,1:3)
  ref_elem_nodes_xyz(1:2,1:4,Geom_Quad) = ref_quad_nodes_xyz(1:2,1:4)
  ref_elem_nodes_xyz(1:3,1:4,Geom_Tetr) = ref_tetr_nodes_xyz(1:3,1:4)
  ref_elem_nodes_xyz(1:3,1:5,Geom_Pyra) = ref_pyra_nodes_xyz(1:3,1:5)
  ref_elem_nodes_xyz(1:3,1:6,Geom_Pris) = ref_pris_nodes_xyz(1:3,1:6)
  ref_elem_nodes_xyz(1:3,1:8,Geom_Hexa) = ref_hexa_nodes_xyz(1:3,1:8)
  !
end subroutine define_ref_elem_node_coordinates
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
#ifdef SPECIAL_FOR_NVIDIA
  use compatibility_pgi_mod, only : execute_command_line
#endif
#else
#ifdef SPECIAL_FOR_NVIDIA
  use compatibility_pgi_mod, only : system,ierrno,gerror
#else
#ifdef SPECIAL_FOR_GNU
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
  if (mypnum == glb_root) then
    write (iout,1)
    flush (iout)
  end if
  flush (iout)
  1 format (/,2x,"Starting GFR after initializing the MPI environment!",/)
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
    mpi_2flttyp = MPI_2REAL
  else
    mpi_flttyp = MPI_DOUBLE_PRECISION
    mpi_2flttyp = MPI_2DOUBLE_PRECISION
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
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
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
                    optional_message, optional_newline_char)
  !
#ifdef SPECIAL_FOR_INTEL
  !.. Use Statements ..
  use ifcore, only : tracebackqq
#endif
  !
  !.. Formal Arguments ..
  integer,          intent(in) :: iend
  character(len=*), intent(in) :: procedure_name
  integer,          intent(in) :: line_number
  character(len=*), intent(in) :: file_name
  !
  !.. Optional Arguments ..
  character(len=*), optional, intent(in) :: optional_message
  character(len=1), optional, intent(in) :: optional_newline_char
  !
  !.. Local Scalars ..
  real(wp) :: val
  integer(int_mpi) :: iflag,io
  character(len=10000) :: strng
  character(len=10000) :: output_message
  character(len=1)     :: newline_char
  !
continue
  !
  iflag = int( iend , kind=kind(iflag) )
  !
  ! Reformat optional_message
  !
  if (present(optional_message)) then
    output_message = format_message(optional_message,optional_newline_char)
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
    if (iflag == dump_stack) then
      !
      ! Abort by performing an divide by zero
      !
      write (iout,3) mypnum
      write (iout,4) trim(adjustl(procedure_name)),line_number, &
                     trim(adjustl(file_name))
      if (present(optional_message)) &
        write (iout,5) trim(adjustl(output_message))
      write (iout,6)
      flush (iout)
#ifdef SPECIAL_FOR_INTEL
      call tracebackqq()
#else
      open (newunit=io,file="ThisFileShouldNeverExist.EOF.read", &
                       status="old",form="unformatted")
      read (io) val
      val = one/val
#endif
      !
    else if (ncpu == 1 .or. iflag == stop_mpi) then
      !
      if (mypnum == 0) then
        write (iout,2)
        write (iout,4) trim(adjustl(procedure_name)),line_number, &
                       trim(adjustl(file_name))
        if (present(optional_message)) &
          write (iout,5) trim(adjustl(output_message))
        write (iout,6)
      end if
#ifdef ATHENA_PBS
      !
      ! If on NAS Athena (Turin) nodes, reset the FP exceptions before
      ! finalizing MPI
      !
      call reset_fp_exceptions
#endif
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
#include "Functions/debug_timer.f90"
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
function format_message(message,optional_newline_char) result(return_value)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: message
  !
  !.. Optional Arguments ..
  character(len=*), optional, intent(in) :: optional_newline_char
  !
  !.. Function Result ..
  character(len=:), allocatable :: return_value
  !
  !.. Local Scalars ..
  character(len=1) :: newline_char
  !
continue
  !
  ! Placeholder char in error_message to be converted to newline.
  !
  newline_char = "$"  ! default: $ becomes newline
  if ( present(optional_newline_char) ) newline_char = optional_newline_char
  !
  if (scan(message,new_line('a')) > 0) then
    !
    ! message already contains new line characters so just copy
    ! it straight into return_value.
    !
    return_value = trim(adjustl(message)) ! F2003 AUTO-REALLOC
    !
  elseif (scan(message,newline_char) > 0) then
    !
    ! Optional message contains a placeholder char for new_line.
    ! Replace placeholder with new_line and copy to return_value.
    !
    ! F2003 AUTO-REALLOC
    return_value = replace_char( message, newline_char, new_line('a') )
    !
    return_value = trim(adjustl(return_value)) ! F2003 AUTO-REALLOC
    !
  else
    !
    ! message does NOT contain any new line characters so reformat
    ! the string so that it is multiple lines with a max of 50 characters
    ! per line.
    !
    return_value = reformat_string( message ) ! F2003 AUTO-REALLOC
    !
  end if
  !
end function format_message
!
!###############################################################################
!> Replace all instances of one character in a string with another character.
function replace_char(string,old_char,new_char) result(new_string)
  character(len=*), intent(in) :: string
  character(len=1), intent(in) :: old_char
  character(len=1), intent(in) :: new_char
  character(len=len(string))   :: new_string
  integer n

  new_string = string
  do n = 1, len(new_string)
    if ( new_string(n:n) == old_char ) new_string(n:n) = new_char
  enddo
end function replace_char
!
!###############################################################################
!
! Include statements for functions that have a generic interface
!
#include "Functions/reallocate.f90"
include "Functions/init_array.f90"
include "Functions/spread.f90"
include "Functions/cross_product.f90"
include "Functions/deltafun.f90"
include "Functions/norms.f90"
include "Functions/reciprocal.f90"
#include "Functions/matrix_ops.f90"
include "Functions/random.f90"
include "Functions/face_normal.f90"
include "Functions/jacobians.f90"
include "Functions/reverse_metrics.f90"
include "Functions/chop.f90"
include "Functions/intseq.f90"
!
! Include statements for functions that DO NOT have a generic interface
!
include "Functions/lowercase.f90"
include "Functions/uppercase.f90"
include "Functions/bc_integer_to_string.f90"
include "Functions/all_are_equal.f90"
include "Functions/reformat_string.f90"
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
