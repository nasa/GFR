module module_kind_types
  !
  !.. Use Statements ..
  use mpi
  !
  implicit none
  !
  public
  !
 !include "mpif.h"
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
 !integer, parameter :: vhp = r4
 !integer, parameter :: vbp = r8
 !integer, parameter :: vdp = r16
  integer, parameter :: vhp = r4
  integer, parameter :: vbp = WP
  integer, parameter :: vdp = r16
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
  real(wp), parameter :: sixteen = 16.0_wp
  !
  real(wp), parameter :: one6   = one/six
  real(wp), parameter :: one4   = one/four
  real(wp), parameter :: one3   = one/three
  real(wp), parameter :: one2   = one/two
  real(wp), parameter :: two3   = two/three
  real(wp), parameter :: three4 = three/four
  real(wp), parameter :: four3  = four/three
  !
  real(wp), parameter :: eps1  = 1.0e-01_wp
  real(wp), parameter :: eps2  = 1.0e-02_wp
  real(wp), parameter :: eps3  = 1.0e-03_wp
  real(wp), parameter :: eps6  = 1.0e-06_wp
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
  real(wp), parameter :: c0_25 = 0.25_wp
  real(wp), parameter :: c0_30 = 0.30_wp
  real(wp), parameter :: c0_35 = 0.35_wp
  real(wp), parameter :: c0_45 = 0.45_wp
  real(wp), parameter :: c0_55 = 0.55_wp
  real(wp), parameter :: c0_65 = 0.65_wp
  real(wp), parameter :: c0_70 = 0.70_wp
  real(wp), parameter :: c0_75 = 0.75_wp
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
  ! Parameters controlled by pre-preprocessor directives
  !
#define p_order_fpp 2
#define n_order_fpp 4
#define q_order_fpp 4
#define iop_solpts_fpp 2
#define iop_flxpts_fpp 2
#define iop_qdtpts_fpp 2
#define iop_tripts_fpp 2
#define metis_opt_fpp 2
#define flux_div_fpp 0
#define iop_corr_fpp 3
#define dflux_opt_fpp 2
#define filt_opt_fpp 0
#define filt_ord_fpp 16
#define filt_alph_fpp 36.0
#define filt_etac_fpp 0.0
#define err_bnd_lim_fpp 2.0
#define err_val_fpp 0
  !
  integer, public, parameter :: p_order_ktmod = p_order_fpp
  integer, public, parameter :: n_order_ktmod = max(n_order_fpp,p_order_fpp)
  integer, public, parameter :: q_order_ktmod = max(n_order_fpp,q_order_fpp)
  !
  integer, public, parameter :: iop_solpts_ktmod = iop_solpts_fpp
  integer, public, parameter :: iop_flxpts_ktmod = iop_flxpts_fpp
  integer, public, parameter :: iop_qdtpts_ktmod = iop_qdtpts_fpp
  integer, public, parameter :: iop_tripts_ktmod = iop_tripts_fpp
  !
  integer, public, parameter :: metis_opt_ktmod = metis_opt_fpp
  !
  integer, public, parameter :: flux_div_ktmod = flux_div_fpp
  integer, public, parameter :: iop_corr_ktmod = iop_corr_fpp
  !
  integer, public, parameter :: dflux_opt_ktmod = dflux_opt_fpp
  !
  integer,  public, parameter :: Filter_Option_ktmod = filt_opt_fpp
  integer,  public, parameter :: Filter_Order_ktmod = filt_ord_fpp
  real(wp), public, parameter :: Filter_Alpha_ktmod = real(filt_alph_fpp,wp)
  real(wp), public, parameter :: Filter_etac_ktmod = real(filt_etac_fpp,wp)
 !real(wp), public, parameter :: ft1_ktmod = max(zero, min(one,Filter_etac_ktmod))
 !real(wp), public, parameter :: ft2_ktmod = real(n_order_ktmod,wp)
 !real(wp), public, parameter :: ft3_ktmod = anint(ft1_ktmod*ft2_ktmod)!+ eps6
 !integer,  public, parameter :: Filter_Cutoff_ktmod = nint(ft3_ktmod)
 !integer,  public, parameter :: Filter_Cutoff_ktmod2 = int(ft3_ktmod)
 !integer,  public, parameter :: Filter_Cutoff_ktmod3 = floor(ft3_ktmod)
  integer,  public, parameter :: &
    Filter_Cutoff_ktmod = int( max(zero,min(one,Filter_etac_ktmod))*real(n_order_ktmod,wp) + eps6 )
  !
  real(wp), public, parameter :: err_bnd_lim_ktmod = real(err_bnd_lim_fpp,wp)
  logical,  public, parameter :: bound_error_ktmod = (err_val_fpp /= 0)
  !
  logical, public, parameter :: &
    using_quadrature_ktmod = .not.(q_order_ktmod==n_order_ktmod .and. iop_solpts_ktmod==iop_qdtpts_ktmod)
  !
  ! I/O constants
  !
  integer, parameter :: iout = 6
  !
  ! Boundary condition constants
  !
  integer, parameter :: bc_inflow = 1
  integer, parameter :: bc_outflow = 2
  integer, parameter :: bc_wall = 3
  integer, parameter :: bc_periodic = 4
  integer, parameter :: bc_freestream = 5
  integer, parameter :: bc_cpu_bnd = 99
  integer, parameter :: bc_unknown = -9999
  integer, parameter :: bc_sub_inflow = 6
  integer, parameter :: bc_sub_outflow = 7
  integer, parameter :: bc_sup_inflow = 8
  integer, parameter :: bc_sup_outflow = 9
  integer, parameter :: bc_characteristic = 10
 !integer, parameter :: bc_freestream = 6
 !integer, parameter :: bc_wall = 7
 !integer, parameter :: bc_periodic = 8
  !
  integer, parameter :: bc_comm(1:2) = [bc_periodic,bc_cpu_bnd]
  !
  ! Parameters defining workshop test cases
  !
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
  integer, parameter :: Legendre_Gauss          = 1
  integer, parameter :: Legendre_Gauss_Lobatto  = 2
  integer, parameter :: Chebyshev_Gauss         = 3
  integer, parameter :: Chebyshev_Gauss_Lobatto = 4
  !
  integer, parameter :: LagrangePolynomial = 1
  integer, parameter :: ChainRuleCV = 2
  integer, parameter :: ChainRuleUtd = 3
  integer, parameter :: DGWeakForm = 4
  integer, parameter :: DGStrongForm = 5
  integer, parameter :: TestQuadrature = 0
  integer, parameter :: OnPrimitiveVariables = 1
  integer, parameter :: OnConservativeVariables = 2
  integer, parameter :: OnFluxes = 3
  !
  integer, parameter :: project_derivatives = 1
  integer, parameter :: derivative_of_proj_fluxes = 2
  !
  ! MPI Scalars 
  !
  integer, protected, save :: mypnum
  integer, protected, save :: ncpu
  integer,            save :: mpierr
  !
  integer, parameter :: abort = 8
  integer, parameter :: stop_mpi = 9
  !
  ! MPI datatypes to match Fortran intrinsic datatypes
  !
  integer, protected, save :: mpi_flttyp ! working real precision
  integer, protected, save :: mpi_hlftyp ! half working real precision
  integer, protected, save :: mpi_dbltyp ! double working real precision
  integer, protected, save :: mpi_inttyp ! working integer precision
  !
  logical, private, save :: random_is_initialized = fals
  !
  ! Generic interface for reallocating arrays of varying types and ranks
  !
  interface reallocate
    module procedure reallocate1_i4, reallocate1_i8, &
                     reallocate1_r4, reallocate1_r8, &
                     reallocate1_c,  reallocate2_i4
  end interface reallocate
  !
  interface cross_product
    module procedure cross_product_DP, cross_product_QP
  end interface cross_product
  !
  interface deltafun
    module procedure deltafun_DP!, deltafun_QP
  end interface deltafun
  !
  interface selfdot
    module procedure selfdot_DP, selfdot_QP
  end interface selfdot
  !
  interface rnorm2
    module procedure rnorm2_DP, rnorm2_QP
  end interface rnorm2
  !
  interface unit_vector
    module procedure unit_vector_DP, unit_vector_QP
  end interface unit_vector
  !
  interface outerprod
    module procedure outerprod_DP, outerprod_QP
  end interface outerprod
  !
  interface reciprocal
    module procedure reciprocal_DP, reciprocal_QP
  end interface reciprocal
  !
  interface reflect
    module procedure reflect_DP, reflect_QP
  end interface reflect
  !
  interface rotate_2d
    module procedure rotate_2d_DP, rotate_2d_QP
  end interface rotate_2d
  !
  interface rotate
    module procedure rotate_DP, rotate_QP
  end interface rotate
  !
  interface invert_matrix
    module procedure invert_matrix_DP, invert_matrix_QP
  end interface invert_matrix
  !
  interface invert_matrix90
    module procedure invert_matrix90_DP, invert_matrix90_QP
  end interface invert_matrix90
  !
  interface ludcmp
    module procedure ludcmp_DP, ludcmp_QP
  end interface ludcmp
  !
  interface lubksb
    module procedure lubksb_DP, lubksb_QP
  end interface lubksb
  !
  interface determinant
    module procedure determinant_DP, determinant_QP
  end interface determinant
  !
  interface random
    module procedure random_i4,random_sp,random_dp,random_qp
  end interface random
  !
  ! Make the non-generic functions private
  !
  private :: cross_product_DP
  private :: cross_product_QP
  private :: deltafun_DP
 !private :: deltafun_QP
  private :: selfdot_DP
  private :: selfdot_QP
  private :: rnorm2_DP
  private :: rnorm2_QP
  private :: unit_vector_DP
  private :: unit_vector_QP
  private :: outerprod_DP
  private :: outerprod_QP
  private :: reciprocal_DP
  private :: reciprocal_QP
  private :: reflect_DP
  private :: reflect_QP
  private :: rotate_DP
  private :: rotate_QP
  private :: rotate_2d_DP
  private :: rotate_2d_QP
  private :: invert_matrix_DP
  private :: invert_matrix_QP
  private :: invert_matrix90_DP
  private :: invert_matrix90_QP
  private :: ludcmp_DP
  private :: ludcmp_QP
  private :: lubksb_DP
  private :: lubksb_QP
  private :: determinant_DP
  private :: determinant_QP
  private :: random_i4
  private :: random_sp
  private :: random_dp
  private :: random_qp
  !
  ! Make the non-generic reallocation subroutines to private
  !
  private :: reallocate1_i4
  private :: reallocate1_i8
  private :: reallocate1_r4
  private :: reallocate1_r8
  private :: reallocate1_c
  private :: reallocate2_i4
  !
contains
!
!###############################################################################
!
subroutine reallocate1_i4(array,new_size)
  integer, intent(in) :: new_size
  integer(i4), allocatable, dimension(:), intent(inout) :: array
  integer(i4), allocatable, dimension(:) :: temp_array
  integer :: ierr
continue
  ! Allocate the temporary array to store the old contents of array
  allocate ( temp_array(1:new_size) , source=0_i4 , stat=ierr )
  call alloc_error("reallocate1_i4","temp_array",1,ierr)
  ! If array is allocated, copy the old contents
  ! of array into the temporary array
  if (allocated(array)) temp_array(1:size(array)) = array
  ! move_alloc does:
  !     1. deallocates array if it is allocated
  !     2. allocates array with the same type, type parameters,
  !        array bounds, and value as temp_array
  !     3. deallocates temp_array
  call move_alloc( from=temp_array , to=array )
end subroutine reallocate1_i4
!
!###############################################################################
!
subroutine reallocate1_i8(array,new_size)
  integer, intent(in) :: new_size
  integer(i8), allocatable, dimension(:), intent(inout) :: array
  integer(i8), allocatable, dimension(:) :: temp_array
  integer :: ierr
continue
  allocate ( temp_array(1:new_size) , source=0_i8 , stat=ierr )
  call alloc_error("reallocate1_i8","temp_array",1,ierr)
  if (allocated(array)) temp_array(1:size(array)) = array
  call move_alloc( from=temp_array , to=array )
end subroutine reallocate1_i8
!
!###############################################################################
!
subroutine reallocate1_r4(array,new_size)
  integer, intent(in) :: new_size
  real(r4), allocatable, dimension(:), intent(inout) :: array
  real(r4), allocatable, dimension(:) :: temp_array
  integer :: ierr
continue
  allocate ( temp_array(1:new_size) , source=0.0_r4 , stat=ierr )
  call alloc_error("reallocate1_r4","temp_array",1,ierr)
  if (allocated(array)) temp_array(1:size(array)) = array
  call move_alloc( from=temp_array , to=array )
end subroutine reallocate1_r4
!
!###############################################################################
!
subroutine reallocate1_r8(array,new_size)
  integer, intent(in) :: new_size
  real(r8), allocatable, dimension(:), intent(inout) :: array
  real(r8), allocatable, dimension(:) :: temp_array
  integer :: ierr
continue
  allocate ( temp_array(1:new_size) , source=0.0_r8 , stat=ierr )
  call alloc_error("reallocate1_r8","temp_array",1,ierr)
  if (allocated(array)) temp_array(1:size(array)) = array
  call move_alloc( from=temp_array , to=array )
end subroutine reallocate1_r8
!
!###############################################################################
!
subroutine reallocate1_c(array,new_size)
  integer, intent(in) :: new_size
  character(len=*), allocatable, dimension(:), intent(inout) :: array
  character(len=len(array)), allocatable, dimension(:) :: temp_array
  integer :: ierr
continue
  allocate ( temp_array(1:new_size) , stat=ierr )
  call alloc_error("reallocate1_c","temp_array",1,ierr)
  if (allocated(array)) temp_array(1:size(array)) = array
  call move_alloc( from=temp_array , to=array )
end subroutine reallocate1_c
!
!###############################################################################
!
subroutine reallocate2_i4(array,new_size1,new_size2)
  !
  !.. Formal Arguments ..
  integer, optional,                       intent(in) :: new_size1
  integer, optional,                       intent(in) :: new_size2
  integer(i4), allocatable, dimension(:,:), intent(inout) :: array
  !
  !.. Local Scalars ..
  integer :: old_size1,old_size2,ierr
  !
  !.. Local Allocatable Arrays ..
  integer(i4), allocatable, dimension(:,:) :: temp_array
  !
continue
  !
  ! Get the old dimensions of array
  !
  old_size1 = 1
  old_size2 = 1
  if (allocated(array)) then
    old_size1 = size(array,1)
    old_size2 = size(array,2)
  end if
  !
  ! Allocate the temporary array to store the old contents of array
  !
  if (present(new_size1) .and. present(new_size2)) then
    allocate ( temp_array(1:new_size1,1:new_size2) , source=0_i4 , stat=ierr )
  else if (present(new_size1)) then
    allocate ( temp_array(1:new_size1,1:old_size2) , source=0_i4 , stat=ierr )
  else if (present(new_size2)) then
    allocate ( temp_array(1:old_size1,1:new_size2) , source=0_i4 , stat=ierr )
  else
    allocate ( temp_array(1:old_size1,1:old_size2) , source=0_i4 , stat=ierr )
  end if
  call alloc_error("reallocate2_i4","temp_array",1,ierr)
  !
  ! Copy the old contents of array into the temporary array
  !
  if (allocated(array)) temp_array(1:old_size1,1:old_size2) = array
  !
  ! move_alloc does:
  !     1. deallocates array
  !     2. allocates array with the same type, type parameters,
  !        array bounds, and value as temp_array
  !     3. deallocates temp_array
  !
  call move_alloc( from=temp_array , to=array )
  !
end subroutine reallocate2_i4
!
!###############################################################################
!
subroutine alloc_error(procedure_name,variable_name,ioptype,iemessag)
  !
  ! Subroutine to generate an error message when
  ! an allocate or deallocate returns an error code
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: procedure_name
  character(len=*), intent(in) :: variable_name
  integer,          intent(in) :: ioptype
  integer,          intent(in) :: iemessag
  !
continue
  !
  if (iemessag /= 0) then
    !
    if (ioptype == 1) then
      write (iout,1) trim(adjustl(variable_name))
    else
      write (iout,2) trim(adjustl(variable_name))
    end if
    !
    write (iout,3) trim(adjustl(procedure_name)),iemessag
    call stop_recon(abort,"alloc_error")
    !
  end if
  !
  ! Format statements
  !
  1 format (/," ***************************************************** ",/, &
              " UNABLE TO ALLOCATE MEMORY FOR VARIABLE ",a)
  2 format (/," ***************************************************** ",/, &
              " UNABLE TO DEALLOCATE MEMORY FOR VARIABLE ",a)
  3 format (" INSIDE PROCEDURE: ",a,/, &
            " ERROR INDICATOR = ",i0,/, &
            " ***************************************************** ",//)
  !
end subroutine alloc_error
!
!###############################################################################
!
subroutine io_error(procedure_name,file_name,ioptype,iemessag)
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
continue
  !
  if (iemessag /= 0) then
    !
    write (iout,10)
    !
    if (ioptype == 1) then
      write (iout,1) trim(adjustl(file_name))
    else if (ioptype == 2) then
      write (iout,2) trim(adjustl(file_name))
    else if (ioptype == 3) then
      write (iout,3) trim(adjustl(file_name))
    else
      write (iout,4) trim(adjustl(file_name))
    end if
    !
    write (iout,11) trim(adjustl(procedure_name)),iemessag
    call stop_recon(abort,"io_error")
    !
  end if
  !
  ! Format Statements
  !
  1 format (" ERROR OPENING FILE '",a,"'")
  2 format (" ERROR CLOSING FILE '",a,"'")
  3 format (" ERROR REOPENING CGNS FILE '",a,"'")
  4 format (" UNKNOWN ERROR WITH FILE '",a,"'")
  10 format (/," ***************************************************** ")
  11 format (" INSIDE PROCEDURE: ",a,/, &
             " ERROR INDICATOR = ",i0,/, &
             " ***************************************************** ",//)
  !
end subroutine io_error
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
  mpi_inttyp = mpi_integer
  !
  ! MPI datatype for half working real precision
  !
 !r = range( 0.0_vhp )
 !p = precision( 0.0_vhp )
 !write (*,*) "half working precision"
 !write (*,*) r,p
 !write (*,*) vhp
 !call mpi_type_create_f90_real(p,r,mpi_hlftyp,mpierr)
 !write (*,*) mpi_hlftyp
  mpi_hlftyp = merge(mpi_real,mpi_double_precision,vhp<r8)
  mpi_flttyp = merge(mpi_real,mpi_double_precision,vbp<r8)
  !
  ! MPI datatype for working real precision
  !
 !r = range( 0.0_vbp )
 !p = precision( 0.0_vbp )
 !write (*,*) "working precision"
 !write (*,*) r,p
 !write (*,*) vbp
 !call mpi_type_create_f90_real(p,r,mpi_flttyp,mpierr)
 !write (*,*) mpi_flttyp
 !mpi_flttyp = mpi_double_precision
  !
  ! MPI datatype for double working real precision
  ! NOTE: mpi_type_create_f90_real doesnt seem to support quad precision
  !
 !r = range( 0.0_vdp )
 !p = precision ( 0.0_vdp )
 !write (*,*) "double working precision"
 !write (*,*) r,p
 !write (*,*) vdp
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
subroutine stop_recon(iend,procedure_name)
  !
  !.. Formal Arguments ..
  integer, optional, intent(in) :: iend
  character(len=*), optional, intent(in) :: procedure_name 
  !
  !.. Local Scalars ..
  integer :: iflag
  !
continue
  !
  iflag = abort
  if (present(iend)) iflag = iend
  !
  if (mypnum == 0) write (iout,1)
  !
  if (ncpu > 0) then
    !
    ! First finalize the overall MPI session before stoping execution
    !
    if (ncpu == 1 .or. iflag == stop_mpi) then
      !
      if (mypnum == 0) then
        write (iout,2)
        if (present(procedure_name)) &
          write (iout,4) trim(adjustl(procedure_name))
        write (iout,5)
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
      write (iout,5)
      !
      call mpi_abort(mpi_comm_world,iflag,mpierr)
      !
    end if
    !
  end if
  !
  flush(iout)
  !
  stop
  !
  ! Format Statements
  !
  1 format (//,"#########################################################")
  2 format (   "    CALLING PREMATURE MPI_FINALIZE ON ALL PROCESSORS! ")
  3 format (   "           MPI ABORT ISSUED FROM PROCESSOR ",i0,"!")
  4 format (   " STOP ISSUED FROM PROCEDURE '",a,"'")
  5 format (   "#########################################################",//)
  !
end subroutine stop_recon
!
!###############################################################################
!
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
include "Functions/uppercase.f90"
!
!###############################################################################
!
end module module_kind_types
