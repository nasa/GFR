module filter_mod
  !
!!!#define PROFILE_ON
  !
  !.. Use Statements ..
  use module_kind_types, global_debug_timer =>      debug_timer, &
                    pure_global_debug_timer => pure_debug_timer
  !
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq
  use eqn_idx, only : ncs
  !
  use matrix_types_mod, only : matrix
  use matrix_types_mod, only : create_matrix_sparsity
  use matrix_types_mod, only : output_sparsity
  !
  use ovar, only : Filter_Option,Limiter_Filter_Option
  use ovar, only : Filter_etac,  Limiter_Filter_etac,  Sponge_Filter_etac
  use ovar, only : Filter_Order, Limiter_Filter_Order, Sponge_Filter_Order
  use ovar, only : Filter_Alpha, Limiter_Filter_Alpha, Sponge_Filter_Alpha
  use ovar, only : Filter_eps,   Limiter_Filter_eps
  !
  implicit none
  !
  private
  !
  !
  ! Module Public Procedures
  !
  public :: init_filter_matrices
  public :: filter_usp
  public :: filter_memory_usage
  public :: move_filter_to_device
  !
  ! Module public arrays
  !
  !-----------------------------------------------------------------------------
  !
  ! filter_t : derived type that contains various types of filtering matrices
  !
  type, public :: filter_t
    ! exponential : matrix to filter a polynomial by multiplying the polynomial
    !               modes by an exponentially decreasing function
    type(matrix), allocatable :: exponential
    ! limiter : filter matrix to use as a limiter for troubled cells
    type(matrix), allocatable :: limiter
    ! cutoff : matrix to filter a polynomial by
    !          cutting-off (zeroing) higher modes
    type(matrix), allocatable :: cutoff(:)
    ! sponge : filter matrix to use as a stronger limiter
    !          for cells that are in a sponge region
    type(matrix), allocatable :: sponge
  end type filter_t
  !
  ! filter : filter matrices for all combinations of cell geometry and
  !          cell order that are possible in the current simulation
  !
  type(filter_t), public, save, allocatable :: filter(:,:)
  !
  !-----------------------------------------------------------------------------
  !
  ! filt_coef_t : derived type that stores 1D filter coefficients
  !
  type :: filt_coef_t
    ! coef : array that contains the 1D filter coefficients
    real(wp), allocatable :: coef(:)
  end type filt_coef_t
  !
  ! global_filt : array containing filter coefficients of varying P
  !               for the global filter, e.g. global_filt(i)%coef(1:i+1)
  !               are the 1D global filter coefficients for a polynomial
  !               of degree i.
  !
  type(filt_coef_t), allocatable :: global_filt(:)
  !
  ! limiter_filt : array containing filter coefficients of varying P
  !                for the limiter filter, e.g. limiter_filt(i)%coef(1:i+1)
  !                are the 1D limiter filter coefficients for a polynomial
  !                of degree i.
  !
  type(filt_coef_t), allocatable :: limiter_filt(:)
  !
  ! sponge_filt : array containing filter coefficients of varying P
  !               for the sponge filter, e.g. sponge_filt(i)%coef(1:i+1)
  !               are the 1D sponge filter coefficients for a polynomial
  !               of degree i.
  !
  type(filt_coef_t), allocatable :: sponge_filt(:)
  !
  !-----------------------------------------------------------------------------
  !
  ! Cutoff_Filter_Option : method used to compute the cutoff filter matrices
  !          = 1 : use a tensor product of the 1D filter with the filter cutoff
  !                based on the order of the 1D polynomial
  !          = 2 : the filter cutoff is based on total polynomial order,
  !                i.e., for the polynomial (x^a * y^b * z*c),
  !                      the total polynomial order = a+b+c
  !
  integer, parameter :: Cutoff_Filter_Option = 1
 !integer, parameter :: Cutoff_Filter_Option = 2
  !
  integer, save :: io
  !
contains
!
!###############################################################################
!
subroutine new_init_filter_matrices()
  !
  !.. Use Statements ..
  use ovar,            only : output_filter_diagonals
  use ovar,            only : use_sponge_filter
  use order_mod,       only : geom_solpts
  use order_mod,       only : n_min_geom,n_max_geom
  use order_mod,       only : n_min_order,n_max_order
  use order_mod,       only : geom_is_used
  use vandermonde_mod, only : vand
  !
  !.. Local Scalars ..
  integer :: ierr,gmin,gmax,omin,omax
  integer :: this_geom,this_order,npts
  integer :: cutoff_order,i,j
  character(len=200) :: array_name
  character(len=300) :: timer_message
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "new_init_filter_matrices"
  character(len=*), parameter :: fname = "filter_diagonals.dat"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (output_filter_diagonals .and. mypnum == glb_root) then
    open (newunit=io,file=fname,status="replace",action="write", &
                     iostat=ierr,iomsg=error_message)
    call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Allocate the filter array
  !
  gmin = n_min_geom
  gmax = n_max_geom
  omin = n_min_order
  omax = n_max_order
  !
  allocate ( filter(gmin:gmax,omin:omax) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"filter",1,__LINE__,__FILE__,ierr,error_message)
  !
  geom_loop: do this_geom = gmin,gmax
    !
    if (.not. geom_is_used(this_geom)) cycle geom_loop
    if (all(this_geom /= Geom_Valid)) cycle geom_loop
    if (this_geom == Geom_Node) cycle geom_loop
    !
    order_loop: do this_order = omin,omax
      !
      ! Get the number of solution points for this cell
      !
      npts = geom_solpts(this_geom,this_order)
      !
      ! ###################################
      ! #####  1. Exponential Filter  #####
      ! ###################################
      !
      if (Filter_Option > 0 .or. output_filter_diagonals) then
        !
        write (timer_message,200) trim(adjustl(Geom_Name(this_geom))),this_order
        200 format ("Exponential filter: Geometry = ",a,"  ;  Order = ",i0)
        call debug_timer(start_timer,timer_message)
        !
        ! Allocate the exponential filter matrix
        !
        allocate ( filter(this_geom,this_order)%exponential , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"exponential"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__, &
                         ierr,error_message)
        !
        ! Compute the exponential global filter matrix
        !
        call create_exp_filter_matrix(this_geom,this_order, &
                       Filter_Option,Filter_etac,Filter_Order,Filter_alpha, &
                       global_filt,filter(this_geom,this_order)%exponential%mat)
        !
        ! Write out the filter diagonals if requested
        !
        if (output_filter_diagonals .and. mypnum == glb_root) then
          write (io,800) trim(Geom_Full_Name(this_geom)), " Global Filter"
          if (Filter_Option == 0) write (io,802)
          do j = npts,1,-1
            write (io,801) (filter(this_geom,this_order)% &
                                             exponential%mat(i,j), i=1,npts)
          end do
        end if
        !
#ifdef SPECIAL_OUTPUT_SPARSITY
        ! Output the sparsity if requested
        !
        write (array_name,10) Geom_Short_Name(this_geom),this_order,"exp"
        call output_sparsity( filter(this_geom,this_order)%exponential , &
                              array_name )
#endif
        !
        call debug_timer(stop_timer,timer_message)
        !
      end if
      !
      ! ###############################
      ! #####  2. Limiter Filter  #####
      ! ###############################
      !
      if (Limiter_Filter_Option /= 0 .or. output_filter_diagonals) then
        !
        write (timer_message,201) trim(adjustl(Geom_Name(this_geom))),this_order
        201 format ("Limiter filter: Geometry = ",a,"  ;  Order = ",i0)
        call debug_timer(start_timer,timer_message)
        !
        ! Allocate the limiter filter matrix
        !
        allocate ( filter(this_geom,this_order)%limiter , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"limiter"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__, &
                         ierr,error_message)
        !
        ! Compute the exponential limiter filter matrix
        !
        call create_exp_filter_matrix(this_geom,this_order, &
                       Limiter_Filter_Option,Limiter_Filter_etac, &
                       Limiter_Filter_Order,Limiter_Filter_alpha, &
                       limiter_filt,filter(this_geom,this_order)%limiter%mat)
        !
        ! Write out the filter diagonals if requested
        !
        if (output_filter_diagonals .and. mypnum == glb_root) then
          !
          write (io,800) trim(Geom_Full_Name(this_geom)), "Limiter Filter"
          if (Limiter_Filter_Option == 0) write (io,802)
          do j = npts,1,-1
            write (io,801) (filter(this_geom,this_order)% &
                                                 limiter%mat(i,j), i=1,npts)
          end do
          !
        end if
        !
#ifdef SPECIAL_OUTPUT_SPARSITY
        ! Output the sparsity if requested
        !
        write (array_name,10) Geom_Short_Name(this_geom),this_order,"lim"
        call output_sparsity( filter(this_geom,this_order)%limiter , &
                              array_name )
#endif
        !
        call debug_timer(stop_timer,timer_message)
        !
      end if
      !
      ! ################################
      ! #####  3. Cut-Off Filters  #####
      ! ################################
      !
      ! Allocate the cutoff filter matrix
      !
      allocate ( filter(this_geom,this_order)%cutoff(1:this_order) , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"cutoff"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      cutoff_loop: do cutoff_order = 1,this_order
        write (timer_message,202) trim(adjustl(Geom_Name(this_geom))), &
                                  this_order,cutoff_order
        202 format ("Cutoff filter: Geometry = ",a,"  ;  Order = ",i0, &
                    "  ;  Cutoff Order = ",i0)
        call debug_timer(start_timer,timer_message)
        !
        ! Special cutoff filter for an edge element using cutoff_order+1
        ! NOTE: This will be replaced with the cutoff filter using cutoff_order
        !
        if (output_filter_diagonals .and. mypnum == glb_root) then
          if (cutoff_order == this_order .and. this_geom == Geom_Edge) then
            ! F2003 AUTO-REALLOC
            filter(this_geom,this_order)%cutoff(cutoff_order)%mat = &
                       create_cutoff_filter(this_geom,this_order,cutoff_order+1)
            write (io,800) trim(Geom_Full_Name(this_geom)), "Cutoff Filter",&
                           " for P-",this_order+1-(cutoff_order+1)
            do j = npts,1,-1
              write (io,801) (filter(this_geom,this_order)% &
                                      cutoff(cutoff_order)%mat(i,j), i=1,npts)
            end do
          end if
        end if
        !
        ! Compute the cutoff filter of order cutoff_order
        !
        ! F2003 AUTO-REALLOC
        filter(this_geom,this_order)%cutoff(cutoff_order)%mat = &
                         create_cutoff_filter(this_geom,this_order,cutoff_order)
        !
        ! Write out the filter diagonals if requested
        !
        if (output_filter_diagonals .and. mypnum == glb_root) then
          write (io,800) trim(Geom_Full_Name(this_geom)), "Cutoff Filter",&
                         " for P-",this_order+1-cutoff_order
          do j = npts,1,-1
            write (io,801) (filter(this_geom,this_order)% &
                                    cutoff(cutoff_order)%mat(i,j), i=1,npts)
          end do
        end if
        !
#ifdef SPECIAL_OUTPUT_SPARSITY
        ! Output the sparsity if requested
        !
        write (array_name,11) Geom_Short_Name(this_geom),this_order,cutoff_order
        call output_sparsity( filter(this_geom,this_order)% &
                                      cutoff(cutoff_order) , array_name )
#endif
        !
        call debug_timer(stop_timer,timer_message)
        !
      end do cutoff_loop
      !
      ! #################################
      ! #####  3. Sponge Filter  #####
      ! #################################
      !
      write (timer_message,203) trim(adjustl(Geom_Name(this_geom))),this_order
      203 format ("Sponge filter: Geometry = ",a,"  ;  Order = ",i0)
      call debug_timer(start_timer,timer_message)
      !
      ! Allocate the sponge filter matrix
      !
      allocate ( filter(this_geom,this_order)%sponge , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"sponge"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      ! Compute the exponential sponge filter matrix
      !
      call create_exp_filter_matrix(this_geom,this_order, &
                 1,Sponge_Filter_etac,Sponge_Filter_Order,Sponge_Filter_alpha, &
                 sponge_filt,filter(this_geom,this_order)%sponge%mat)
      !
      ! Write out the filter diagonals if requested
      !
      if (output_filter_diagonals .and. mypnum == glb_root) then
        !
        write (io,800) trim(Geom_Full_Name(this_geom)), "Sponge Filter"
        if (.not. use_sponge_filter) write (io,802)
        do j = npts,1,-1
          write (io,801) (filter(this_geom,this_order)% &
                                                sponge%mat(i,j), i=1,npts)
        end do
        !
      end if
      !
#ifdef SPECIAL_OUTPUT_SPARSITY
      ! Output the sparsity if requested
      !
      write (array_name,10) Geom_Short_Name(this_geom),this_order,"spg"
      call output_sparsity( filter(this_geom,this_order)%sponge , &
                            array_name )
#endif
      !
      call debug_timer(stop_timer,timer_message)
      !
    end do order_loop
    !
  end do geom_loop
  !
  if (output_filter_diagonals) then
    if (mypnum == glb_root) then
      flush (io)
      close (io,iostat=ierr,iomsg=error_message)
      call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
    end if
    write (error_message,803)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Check for sparsity
  !
  do this_geom = gmin,gmax
    do this_order = omin,omax
      if (allocated(filter(this_geom,this_order)%exponential)) then
        write (timer_message,300) trim(adjustl(Geom_Name(this_geom))),this_order
        300 format ("Exponential filter sparsity check:  Geometry = ", &
                    a,"  ;  Order = ",i0)
        call debug_timer(start_timer,timer_message)
        call create_matrix_sparsity( filter(this_geom,this_order)%exponential )
        call debug_timer(stop_timer,timer_message)
      end if
      if (allocated(filter(this_geom,this_order)%limiter)) then
        write (timer_message,301) trim(adjustl(Geom_Name(this_geom))),this_order
        301 format ("Limiter filter sparsity check:  Geometry = ", &
                    a,"  ;  Order = ",i0)
        call debug_timer(start_timer,timer_message)
        call create_matrix_sparsity( filter(this_geom,this_order)%limiter )
        call debug_timer(stop_timer,timer_message)
      end if
      if (allocated(filter(this_geom,this_order)%cutoff)) then
        write (timer_message,302) trim(adjustl(Geom_Name(this_geom))),this_order
        302 format ("Cutoff filter sparsity check:  Geometry = ", &
                    a,"  ;  Order = ",i0)
        call debug_timer(start_timer,timer_message)
        call create_matrix_sparsity( filter(this_geom,this_order)%cutoff )
        call debug_timer(stop_timer,timer_message)
      end if
      if (allocated(filter(this_geom,this_order)%sponge)) then
        write (timer_message,303) trim(adjustl(Geom_Name(this_geom))),this_order
        303 format ("Sponge filter sparsity check:  Geometry = ", &
                    a,"  ;  Order = ",i0)
        call debug_timer(start_timer,timer_message)
        call create_matrix_sparsity( filter(this_geom,this_order)%sponge )
        call debug_timer(stop_timer,timer_message)
      end if
    end do
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("filter(",a,",",i0,")%",a,:,"(",i0,")%mat")
  100 format (" A grid cell of an unknown geometry type was found!")
  !
  800 format (/,a,1x,a,:,a,i0)
  801 format (100(2x,es14.6))
  802 format ("  NOTE: CURRENTLY NOT USED BECAUSE ", &
              "IT IS TURNED OFF IN THE INPUT FILE!")
  803 format ("AFTER WRITING OUT THE FILTER DIAGONALS ", &
              "TO 'filter_diagonals.dat'.")
  !
  10 format (a4,".filter.p",i0,".",a3,".tec.dat")
  11 format (a4,".filter.p",i0,".cut-p",i0,".tec.dat")
  !
end subroutine new_init_filter_matrices
!
!###############################################################################
!
subroutine init_filter_matrices()
  !
  !.. Use Statements ..
  use ovar,            only : output_filter_diagonals
  use order_mod,       only : geom_solpts
  use order_mod,       only : n_min_geom,n_max_geom
  use order_mod,       only : n_min_order,n_max_order
  use order_mod,       only : geom_is_used
  use vandermonde_mod, only : vand
  !
  !.. Local Scalars ..
  integer :: ierr,gmin,gmax,omin,omax
  integer :: this_geom,this_order,npts
  integer :: cutoff_order,i,j
  character(len=200) :: array_name
  character(len=300) :: timer_message
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_filter_matrices"
  character(len=*), parameter :: fname = "filter_diagonals.dat"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (output_filter_diagonals .and. mypnum == glb_root) then
    open (newunit=io,file=fname,status="replace",action="write", &
                     iostat=ierr,iomsg=error_message)
    call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Allocate the filter array
  !
  gmin = n_min_geom
  gmax = n_max_geom
  omin = n_min_order
  omax = n_max_order
  !
  allocate ( filter(gmin:gmax,omin:omax) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"filter",1,__LINE__,__FILE__,ierr,error_message)
  !
  geom_loop: do this_geom = gmin,gmax
    !
    if (.not. geom_is_used(this_geom)) cycle geom_loop
    if (all(this_geom /= Geom_Valid)) cycle geom_loop
    if (this_geom == Geom_Node) cycle geom_loop
    !
    order_loop: do this_order = omin,omax
      !
      ! Get the number of solution points for this cell
      !
      npts = geom_solpts(this_geom,this_order)
      !
      ! ###################################
      ! #####  1. Exponential Filter  #####
      ! ###################################
      !
      if (Filter_Option > 0 .or. output_filter_diagonals) then
        !
        write (timer_message,200) trim(adjustl(Geom_Name(this_geom))),this_order
        200 format ("Exponential filter: Geometry = ",a,"  ;  Order = ",i0)
        call debug_timer(start_timer,timer_message)
        !
        ! Allocate the exponential filter matrix
        !
        allocate ( filter(this_geom,this_order)%exponential , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"exponential"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__, &
                         ierr,error_message)
        !
        ! Allocate the exponential filter matrix
        !
        allocate ( filter(this_geom,this_order)% &
                                    exponential%mat(1:npts,1:npts) , &
                   source=zero , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"exponential%mat"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__, &
                         ierr,error_message)
        !
        ! Compute the filter matrices depending on geometry type
        !
        select case (this_geom)
          case (Geom_Edge)
            !
            filter(this_geom,this_order)%exponential%mat = &
               FilterMatrix_Edge( vand(this_geom,this_order)%modal2nodal%mat , &
                                  vand(this_geom,this_order)%nodal2modal%mat )
            !
            if (output_filter_diagonals .and. mypnum == glb_root) then
              write (io,800) "Edge Standard Filter"
              if (Filter_Option == 0) write (io,802)
              do j = npts,1,-1
                write (io,801) (filter(this_geom,this_order)% &
                                                 exponential%mat(i,j), i=1,npts)
              end do
            end if
            !
          case (Geom_Tria)
            !
            filter(this_geom,this_order)%exponential%mat = &
               FilterMatrix_Tria( vand(this_geom,this_order)%modal2nodal%mat , &
                                  vand(this_geom,this_order)%nodal2modal%mat )
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,10) "tria",this_order,"exp"
            call output_sparsity( filter(this_geom,this_order)%exponential , &
                                  array_name )
#endif
            !
          case (Geom_Quad)
            !
            filter(this_geom,this_order)%exponential%mat = &
               FilterMatrix_Quad( vand(this_geom,this_order)%modal2nodal%mat , &
                                  vand(this_geom,this_order)%nodal2modal%mat )
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,10) "quad",this_order,"exp"
            call output_sparsity( filter(this_geom,this_order)%exponential , &
                                  array_name )
#endif
            !
          case (Geom_Tetr)
            !
            filter(this_geom,this_order)%exponential%mat = &
               FilterMatrix_Tetr( vand(this_geom,this_order)%modal2nodal%mat , &
                                  vand(this_geom,this_order)%nodal2modal%mat )
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,10) "tetr",this_order,"exp"
            call output_sparsity( filter(this_geom,this_order)%exponential , &
                                  array_name )
#endif
            !
          case (Geom_Pyra)
            !
            filter(this_geom,this_order)%exponential%mat = &
               FilterMatrix_Pyra( vand(this_geom,this_order)%modal2nodal%mat , &
                                  vand(this_geom,this_order)%nodal2modal%mat )
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,10) "pyra",this_order,"exp"
            call output_sparsity( filter(this_geom,this_order)%exponential , &
                                  array_name )
#endif
            !
          case (Geom_Pris)
            !
            filter(this_geom,this_order)%exponential%mat = &
               FilterMatrix_Pris( vand(this_geom,this_order)%modal2nodal%mat , &
                                  vand(this_geom,this_order)%nodal2modal%mat )
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,10) "pris",this_order,"exp"
            call output_sparsity( filter(this_geom,this_order)%exponential , &
                                  array_name )
#endif
            !
          case (Geom_Hexa)
            !
            filter(this_geom,this_order)%exponential%mat = &
               FilterMatrix_Hexa( vand(this_geom,this_order)%modal2nodal%mat , &
                                  vand(this_geom,this_order)%nodal2modal%mat )
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,10) "hexa",this_order,"exp"
            call output_sparsity( filter(this_geom,this_order)%exponential , &
                                  array_name )
#endif
            !
          case default
            !
            write (error_message,100)
            call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
            !
        end select
        !
        call debug_timer(stop_timer,timer_message)
        !
      end if
      !
      ! ###############################
      ! #####  2. Limiter Filter  #####
      ! ###############################
      !
      if (Limiter_Filter_Option /= 0 .or. output_filter_diagonals) then
        !
        write (timer_message,201) trim(adjustl(Geom_Name(this_geom))),this_order
        201 format ("Limiter filter: Geometry = ",a,"  ;  Order = ",i0)
        call debug_timer(start_timer,timer_message)
        !
        ! Allocate the limiter filter matrix
        !
        allocate ( filter(this_geom,this_order)%limiter , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"limiter"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__, &
                         ierr,error_message)
        !
        ! Allocate the limiter filter matrix
        !
        allocate ( filter(this_geom,this_order)%limiter%mat(1:npts,1:npts) , &
                   source=zero , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"limiter%mat"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__, &
                         ierr,error_message)
        !
        ! Compute the filter matrices depending on geometry type
        !
        select case (this_geom)
          case (Geom_Edge)
            !
            if (output_filter_diagonals .and. mypnum == glb_root) then
              !
              filter(this_geom,this_order)%limiter%mat = &
               FilterMatrix_Edge_with_output( vand(this_geom,this_order)% &
                                                             modal2nodal%mat , &
                                              vand(this_geom,this_order)% &
                                                             nodal2modal%mat , &
                                                   cutoff_order=-1 )
              write (io,800) "Edge Limiter Filter"
              if (Limiter_Filter_Option == 0) write (io,802)
              do j = npts,1,-1
                write (io,801) (filter(this_geom,this_order)% &
                                                     limiter%mat(i,j), i=1,npts)
              end do
              !
            else
              !
              filter(this_geom,this_order)%limiter%mat = &
               FilterMatrix_Edge( vand(this_geom,this_order)%modal2nodal%mat , &
                                  vand(this_geom,this_order)%nodal2modal%mat , &
                                  cutoff_order=-1 )
              !
            end if
            !
          case (Geom_Tria)
            !
            filter(this_geom,this_order)%limiter%mat = &
               FilterMatrix_Tria( vand(this_geom,this_order)%modal2nodal%mat , &
                                  vand(this_geom,this_order)%nodal2modal%mat , &
                                  cutoff_order=-1 )
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,10) "tria",this_order,"lim"
            call output_sparsity( filter(this_geom,this_order)%limiter , &
                                  array_name )
#endif
            !
          case (Geom_Quad)
            !
            filter(this_geom,this_order)%limiter%mat = &
               FilterMatrix_Quad( vand(this_geom,this_order)%modal2nodal%mat , &
                                  vand(this_geom,this_order)%nodal2modal%mat , &
                                  cutoff_order=-1 )
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,10) "quad",this_order,"lim"
            call output_sparsity( filter(this_geom,this_order)%limiter , &
                                  array_name )
#endif
            !
          case (Geom_Tetr)
            !
            filter(this_geom,this_order)%limiter%mat = &
               FilterMatrix_Tetr( vand(this_geom,this_order)%modal2nodal%mat , &
                                  vand(this_geom,this_order)%nodal2modal%mat , &
                                  cutoff_order=-1 )
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,10) "tetr",this_order,"lim"
            call output_sparsity( filter(this_geom,this_order)%limiter , &
                                  array_name )
#endif
            !
          case (Geom_Pyra)
            !
            filter(this_geom,this_order)%limiter%mat = &
               FilterMatrix_Pyra( vand(this_geom,this_order)%modal2nodal%mat , &
                                  vand(this_geom,this_order)%nodal2modal%mat , &
                                  cutoff_order=-1 )
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,10) "pyra",this_order,"lim"
            call output_sparsity( filter(this_geom,this_order)%limiter , &
                                  array_name )
#endif
            !
          case (Geom_Pris)
            !
            filter(this_geom,this_order)%limiter%mat = &
               FilterMatrix_Pris( vand(this_geom,this_order)%modal2nodal%mat , &
                                  vand(this_geom,this_order)%nodal2modal%mat , &
                                  cutoff_order=-1 )
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,10) "pris",this_order,"lim"
            call output_sparsity( filter(this_geom,this_order)%limiter , &
                                  array_name )
#endif
            !
          case (Geom_Hexa)
            !
            if (output_filter_diagonals .and. mypnum == glb_root) then
              !
              filter(this_geom,this_order)%limiter%mat = &
               FilterMatrix_Hexa_with_output( vand(this_geom,this_order)% &
                                                             modal2nodal%mat , &
                                              vand(this_geom,this_order)% &
                                                             nodal2modal%mat , &
                                                   cutoff_order=-1 )
              write (io,800) "Hexa Limiter Filter"
              if (Limiter_Filter_Option == 0) write (io,802)
              do j = npts,1,-1
                write (io,801) (filter(this_geom,this_order)% &
                                                     limiter%mat(i,j), i=1,npts)
              end do
              !
            else
              !
              filter(this_geom,this_order)%limiter%mat = &
                 FilterMatrix_Hexa( vand(this_geom,this_order)%modal2nodal%mat,&
                                    vand(this_geom,this_order)%nodal2modal%mat,&
                                    cutoff_order=-1 )
              !
            end if
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,10) "hexa",this_order,"lim"
            call output_sparsity( filter(this_geom,this_order)%limiter , &
                                  array_name )
#endif
            !
          case default
            !
            write (error_message,100)
            call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
            !
        end select
        !
        call debug_timer(stop_timer,timer_message)
        !
      end if
      !
      ! ################################
      ! #####  3. Cut-Off Filters  #####
      ! ################################
      !
      ! Allocate the cutoff filter matrix
      !
      allocate ( filter(this_geom,this_order)%cutoff(1:this_order) , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"cutoff"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      cutoff_loop: do cutoff_order = 1,this_order
        write (timer_message,202) trim(adjustl(Geom_Name(this_geom))), &
                                  this_order,cutoff_order
        202 format ("Cutoff filter: Geometry = ",a,"  ;  Order = ",i0, &
                    "  ;  Cutoff Order = ",i0)
        call debug_timer(start_timer,timer_message)
        !
        ! Allocate the cutoff filter matrix
        !
        allocate ( filter(this_geom,this_order)% &
                           cutoff(cutoff_order)%mat(1:npts,1:npts) , &
                   source=zero , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order, &
                             "cutoff",cutoff_order
        call alloc_error(pname,array_name,1,__LINE__,__FILE__, &
                         ierr,error_message)
        !
        if (output_filter_diagonals .and. mypnum == glb_root) then
          if (cutoff_order == this_order .and. this_geom == Geom_Edge) then
            filter(this_geom,this_order)%cutoff(cutoff_order)%mat = &
               FilterMatrix_Edge( vand(this_geom,this_order)%modal2nodal%mat , &
                                  vand(this_geom,this_order)%nodal2modal%mat , &
                                  cutoff_order+1 )
            write (io,800) "Edge Cutoff Filter",&
                           " for P-",this_order+1-(cutoff_order+1)
            do j = npts,1,-1
              write (io,801) (filter(this_geom,this_order)% &
                                      cutoff(cutoff_order)%mat(i,j), i=1,npts)
            end do
          end if
        end if
        !
        ! Compute the filter matrices depending on geometry type
        !
        select case (this_geom)
          case (Geom_Edge)
            !
            filter(this_geom,this_order)%cutoff(cutoff_order)%mat = &
               FilterMatrix_Edge( vand(this_geom,this_order)%modal2nodal%mat , &
                                  vand(this_geom,this_order)%nodal2modal%mat , &
                                  cutoff_order )
            !
            if (output_filter_diagonals .and. mypnum == glb_root) then
              write (io,800) "Edge Cutoff Filter",&
                             " for P-",this_order+1-cutoff_order
              do j = npts,1,-1
                write (io,801) (filter(this_geom,this_order)% &
                                        cutoff(cutoff_order)%mat(i,j), i=1,npts)
              end do
            end if
            !
          case (Geom_Tria)
            !
            filter(this_geom,this_order)%cutoff(cutoff_order)%mat = &
               FilterMatrix_Tria( vand(this_geom,this_order)%modal2nodal%mat , &
                                  vand(this_geom,this_order)%nodal2modal%mat , &
                                  cutoff_order )
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,11) "tria",this_order,cutoff_order
            call output_sparsity( filter(this_geom,this_order)% &
                                          cutoff(cutoff_order) , array_name )
#endif
            !
          case (Geom_Quad)
            !
            if (output_filter_diagonals .and. mypnum == glb_root) then
              !
              filter(this_geom,this_order)%cutoff(cutoff_order)%mat = &
                 FilterMatrix_Quad_with_output( vand(this_geom,this_order)%modal2nodal%mat , &
                                    vand(this_geom,this_order)%nodal2modal%mat , &
                                    cutoff_order )
              write (io,800) "Quad Cutoff Filter",&
                             " for P-",this_order+1-cutoff_order
              do j = npts,1,-1
                write (io,801) (filter(this_geom,this_order)% &
                                        cutoff(cutoff_order)%mat(i,j), i=1,npts)
              end do
              !
            else
              !
              filter(this_geom,this_order)%cutoff(cutoff_order)%mat = &
                 FilterMatrix_Quad( vand(this_geom,this_order)%modal2nodal%mat , &
                                    vand(this_geom,this_order)%nodal2modal%mat , &
                                    cutoff_order )
              !
            end if
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,11) "quad",this_order,cutoff_order
            call output_sparsity( filter(this_geom,this_order)% &
                                          cutoff(cutoff_order) , array_name )
#endif
            !
          case (Geom_Tetr)
            !
            filter(this_geom,this_order)%cutoff(cutoff_order)%mat = &
               FilterMatrix_Tetr( vand(this_geom,this_order)%modal2nodal%mat , &
                                  vand(this_geom,this_order)%nodal2modal%mat , &
                                  cutoff_order )
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,11) "tetr",this_order,cutoff_order
            call output_sparsity( filter(this_geom,this_order)% &
                                          cutoff(cutoff_order) , array_name )
#endif
            !
          case (Geom_Pyra)
            !
            filter(this_geom,this_order)%cutoff(cutoff_order)%mat = &
               FilterMatrix_Pyra( vand(this_geom,this_order)%modal2nodal%mat , &
                                  vand(this_geom,this_order)%nodal2modal%mat , &
                                  cutoff_order )
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,11) "pyra",this_order,cutoff_order
            call output_sparsity( filter(this_geom,this_order)% &
                                          cutoff(cutoff_order) , array_name )
#endif
            !
          case (Geom_Pris)
            !
            filter(this_geom,this_order)%cutoff(cutoff_order)%mat = &
               FilterMatrix_Pris( vand(this_geom,this_order)%modal2nodal%mat , &
                                  vand(this_geom,this_order)%nodal2modal%mat , &
                                  cutoff_order )
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,11) "pris",this_order,cutoff_order
            call output_sparsity( filter(this_geom,this_order)% &
                                          cutoff(cutoff_order) , array_name )
#endif
            !
          case (Geom_Hexa)
            !
            if (output_filter_diagonals .and. mypnum == glb_root) then
              !
              filter(this_geom,this_order)%cutoff(cutoff_order)%mat = &
                FilterMatrix_Hexa_with_output( vand(this_geom,this_order)%modal2nodal%mat , &
                                   vand(this_geom,this_order)%nodal2modal%mat , &
                                   cutoff_order )
              write (io,800) "Hexa Cutoff Filter",&
                             " for P-",this_order+1-cutoff_order
              do j = npts,1,-1
                write (io,801) (filter(this_geom,this_order)% &
                                        cutoff(cutoff_order)%mat(i,j), i=1,npts)
              end do
              !
            else
              !
              filter(this_geom,this_order)%cutoff(cutoff_order)%mat = &
                FilterMatrix_Hexa( vand(this_geom,this_order)%modal2nodal%mat , &
                                   vand(this_geom,this_order)%nodal2modal%mat , &
                                   cutoff_order )
              !
            end if
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,11) "hexa",this_order,cutoff_order
            call output_sparsity( filter(this_geom,this_order)% &
                                          cutoff(cutoff_order) , array_name )
#endif
            !
          case default
            !
            write (error_message,100)
            call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
            !
        end select
        !
        call debug_timer(stop_timer,timer_message)
        !
      end do cutoff_loop
      !
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
      ! #################################
      ! #####  3. Sponge Filter  #####
      ! #################################
      !
      write (timer_message,203) trim(adjustl(Geom_Name(this_geom))),this_order
      203 format ("Sponge filter: Geometry = ",a,"  ;  Order = ",i0)
      call debug_timer(start_timer,timer_message)
      !
      ! Allocate the sponge filter matrix
      !
      allocate ( filter(this_geom,this_order)%sponge , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"sponge"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      ! Allocate the sponge filter matrix
      !
      allocate ( filter(this_geom,this_order)%sponge%mat(1:npts,1:npts) , &
                 source=zero , stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"sponge%mat"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      ! Compute the filter matrices depending on geometry type
      !
      select case (this_geom)
        case (Geom_Edge)
          !
          filter(this_geom,this_order)%sponge%mat = &
           FilterMatrix_Edge( vand(this_geom,this_order)%modal2nodal%mat , &
                              vand(this_geom,this_order)%nodal2modal%mat , &
                              cutoff_order=-2 )
          !
        case (Geom_Tria)
          !
          filter(this_geom,this_order)%sponge%mat = &
             FilterMatrix_Tria( vand(this_geom,this_order)%modal2nodal%mat , &
                                vand(this_geom,this_order)%nodal2modal%mat , &
                                cutoff_order=-2 )
          !
        case (Geom_Quad)
          !
          filter(this_geom,this_order)%sponge%mat = &
             FilterMatrix_Quad( vand(this_geom,this_order)%modal2nodal%mat , &
                                vand(this_geom,this_order)%nodal2modal%mat , &
                                cutoff_order=-2 )
          !
        case (Geom_Tetr)
          !
          filter(this_geom,this_order)%sponge%mat = &
             FilterMatrix_Tetr( vand(this_geom,this_order)%modal2nodal%mat , &
                                vand(this_geom,this_order)%nodal2modal%mat , &
                                cutoff_order=-2 )
          !
        case (Geom_Pyra)
          !
          filter(this_geom,this_order)%sponge%mat = &
             FilterMatrix_Pyra( vand(this_geom,this_order)%modal2nodal%mat , &
                                vand(this_geom,this_order)%nodal2modal%mat , &
                                cutoff_order=-2 )
          !
        case (Geom_Pris)
          !
          filter(this_geom,this_order)%sponge%mat = &
             FilterMatrix_Pris( vand(this_geom,this_order)%modal2nodal%mat , &
                                vand(this_geom,this_order)%nodal2modal%mat , &
                                cutoff_order=-2 )
          !
        case (Geom_Hexa)
          !
          filter(this_geom,this_order)%sponge%mat = &
             FilterMatrix_Hexa( vand(this_geom,this_order)%modal2nodal%mat , &
                                vand(this_geom,this_order)%nodal2modal%mat , &
                                cutoff_order=-2 )
          !
        case default
          !
          write (error_message,100)
          call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
          !
      end select
      !
      call debug_timer(stop_timer,timer_message)
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
      !
    end do order_loop
    !
  end do geom_loop
  !
  if (output_filter_diagonals) then
    if (mypnum == glb_root) then
      flush (io)
      close (io,iostat=ierr,iomsg=error_message)
      call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
    end if
    write (error_message,803)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Check for sparsity
  !
  do this_geom = gmin,gmax
    do this_order = omin,omax
      if (allocated(filter(this_geom,this_order)%exponential)) then
        write (timer_message,300) trim(adjustl(Geom_Name(this_geom))),this_order
        300 format ("Exponential filter sparsity check:  Geometry = ", &
                    a,"  ;  Order = ",i0)
        call debug_timer(start_timer,timer_message)
        call create_matrix_sparsity( filter(this_geom,this_order)%exponential )
        call debug_timer(stop_timer,timer_message)
      end if
      if (allocated(filter(this_geom,this_order)%limiter)) then
        write (timer_message,301) trim(adjustl(Geom_Name(this_geom))),this_order
        301 format ("Limiter filter sparsity check:  Geometry = ", &
                    a,"  ;  Order = ",i0)
        call debug_timer(start_timer,timer_message)
        call create_matrix_sparsity( filter(this_geom,this_order)%limiter )
        call debug_timer(stop_timer,timer_message)
      end if
      if (allocated(filter(this_geom,this_order)%cutoff)) then
        write (timer_message,302) trim(adjustl(Geom_Name(this_geom))),this_order
        302 format ("Cutoff filter sparsity check:  Geometry = ", &
                    a,"  ;  Order = ",i0)
        call debug_timer(start_timer,timer_message)
        call create_matrix_sparsity( filter(this_geom,this_order)%cutoff )
        call debug_timer(stop_timer,timer_message)
      end if
      if (allocated(filter(this_geom,this_order)%sponge)) then
        write (timer_message,303) trim(adjustl(Geom_Name(this_geom))),this_order
        303 format ("Sponge filter sparsity check:  Geometry = ", &
                    a,"  ;  Order = ",i0)
        call debug_timer(start_timer,timer_message)
        call create_matrix_sparsity( filter(this_geom,this_order)%sponge )
        call debug_timer(stop_timer,timer_message)
      end if
    end do
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("filter(",a,",",i0,")%",a,:,"(",i0,")%mat")
  100 format (" A grid cell of an unknown geometry type was found!")
  !
  800 format (/,a,:,a,i0)
  801 format (100(2x,es14.6))
  802 format ("  NOTE: CURRENTLY NOT USED BECAUSE ", &
              "IT IS TURNED OFF IN THE INPUT FILE!")
  803 format ("AFTER WRITING OUT THE FILTER DIAGONALS ", &
              "TO 'filter_diagonals.dat'.")
  !
  10 format (a4,".filter.p",i0,".",a3,".tec.dat")
  11 format (a4,".filter.p",i0,".cut-p",i0,".tec.dat")
  !
end subroutine init_filter_matrices
!
!###############################################################################
!
subroutine filter_usp(update_gpsp)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use flowvar, only : usp,gpsp
  use ovar,    only : force_nd_subroutines
  !
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: update_gpsp
  !
  !.. Local Scalars ..
  integer     :: nr
  logical(lk) :: gpsp_update_requested
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "filter_usp"
  !
continue
  !
#ifdef DEBUG_ON
  call debug_timer(entering_procedure,pname)
#endif
  !
  gpsp_update_requested = fals
  if (present(update_gpsp)) then
    gpsp_update_requested = update_gpsp
  end if
  !
  nr = nme-nmb+1
  !
  if (force_nd_subroutines) then
    !
    call filter_usp_nd(1,nq,usp)
    !
  else
    !
    if (nr == 3) then
      call filter_usp_3d(usp)
    else if (nr == 2) then
      call filter_usp_2d(usp)
    else
      call filter_usp_1d(usp)
    end if
    !
    if (ncs > 0) then
      call filter_usp_nd(1,ncs,usp)
    end if
    !
  end if
  !
  if (gpsp_update_requested) then
    call update_thermo_properties(1,n_solpts,usp,gpsp, &
                                  loc="cell",caller_name=pname)
    call update_transport_properties(1,n_solpts,usp,gpsp)
  end if
  !
#ifdef DEBUG_ON
  call debug_timer(leaving_procedure,pname)
#endif
  !
end subroutine filter_usp
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine filter_usp_3d(usp)
  !
  !.. Use Statements ..
  use order_mod, only : maxSP
  use geovar,    only : ncell,cell
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,filt_ec,filt_mx,filt_my,filt_mz,filt_ee
  !
  !.. Local Arrays ..
  real(wp) :: unfiltered_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "filter_usp_3d"
  !
continue
  !
  ! Apply the filter to each cell
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,filter) &
  !$ACC&              PRIVATE(unfiltered_usp)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Save the current, unfiltered solution into the unfiltered_usp array.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      unfiltered_usp(n,nec) = usp(n0+n,nec)
      unfiltered_usp(n,nmx) = usp(n0+n,nmx)
      unfiltered_usp(n,nmy) = usp(n0+n,nmy)
      unfiltered_usp(n,nmz) = usp(n0+n,nmz)
      unfiltered_usp(n,nee) = usp(n0+n,nee)
    end do
    !
    ! Filter the solution and save it back into the usp array.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      filt_ec = zero
      filt_mx = zero
      filt_my = zero
      filt_mz = zero
      filt_ee = zero
      !
      do i = filter(this_geom,this_order)%exponential%col_ptr(n)+1, &
             filter(this_geom,this_order)%exponential%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%exponential%col_idx(i)
        a = filter(this_geom,this_order)%exponential%col_val(i)
        !
        filt_ec = filt_ec + a * unfiltered_usp(k,nec)
        filt_mx = filt_mx + a * unfiltered_usp(k,nmx)
        filt_my = filt_my + a * unfiltered_usp(k,nmy)
        filt_mz = filt_mz + a * unfiltered_usp(k,nmz)
        filt_ee = filt_ee + a * unfiltered_usp(k,nee)
        !
      end do
      !
      usp(n0+n,nec) = filt_ec
      usp(n0+n,nmx) = filt_mx
      usp(n0+n,nmy) = filt_my
      usp(n0+n,nmz) = filt_mz
      usp(n0+n,nee) = filt_ee
      !
    end do
    !
  end do
  !
end subroutine filter_usp_3d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine filter_usp_2d(usp)
  !
  !.. Use Statements ..
  use order_mod, only : maxSP
  use geovar,    only : ncell,cell
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,filt_ec,filt_mx,filt_my,filt_ee
  !
  !.. Local Arrays ..
  real(wp) :: unfiltered_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "filter_usp_2d"
  !
continue
  !
  ! Apply the filter to each cell
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,filter) &
  !$ACC&              PRIVATE(unfiltered_usp)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Save the current, unfiltered solution into the unfiltered_usp array.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      unfiltered_usp(n,nec) = usp(n0+n,nec)
      unfiltered_usp(n,nmx) = usp(n0+n,nmx)
      unfiltered_usp(n,nmy) = usp(n0+n,nmy)
      unfiltered_usp(n,nee) = usp(n0+n,nee)
    end do
    !
    ! Filter the solution and save it back into the usp array.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      filt_ec = zero
      filt_mx = zero
      filt_my = zero
      filt_ee = zero
      !
      do i = filter(this_geom,this_order)%exponential%col_ptr(n)+1, &
             filter(this_geom,this_order)%exponential%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%exponential%col_idx(i)
        a = filter(this_geom,this_order)%exponential%col_val(i)
        !
        filt_ec = filt_ec + a * unfiltered_usp(k,nec)
        filt_mx = filt_mx + a * unfiltered_usp(k,nmx)
        filt_my = filt_my + a * unfiltered_usp(k,nmy)
        filt_ee = filt_ee + a * unfiltered_usp(k,nee)
        !
      end do
      !
      usp(n0+n,nec) = filt_ec
      usp(n0+n,nmx) = filt_mx
      usp(n0+n,nmy) = filt_my
      usp(n0+n,nee) = filt_ee
      !
    end do
    !
  end do
  !
end subroutine filter_usp_2d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine filter_usp_1d(usp)
  !
  !.. Use Statements ..
  use order_mod, only : maxSP
  use geovar,    only : ncell,cell
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,filt_ec,filt_mx,filt_ee
  !
  !.. Local Arrays ..
  real(wp) :: unfiltered_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "filter_usp_1d"
  !
continue
  !
  ! Apply the filter to each cell
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,filter) &
  !$ACC&              PRIVATE(unfiltered_usp)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Save the current, unfiltered solution into the unfiltered_usp array.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      unfiltered_usp(n,nec) = usp(n0+n,nec)
      unfiltered_usp(n,nmx) = usp(n0+n,nmx)
      unfiltered_usp(n,nee) = usp(n0+n,nee)
    end do
    !
    ! Filter the solution and save it back into the usp array.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      filt_ec = zero
      filt_mx = zero
      filt_ee = zero
      !
      do i = filter(this_geom,this_order)%exponential%col_ptr(n)+1, &
             filter(this_geom,this_order)%exponential%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%exponential%col_idx(i)
        a = filter(this_geom,this_order)%exponential%col_val(i)
        !
        filt_ec = filt_ec + a * unfiltered_usp(k,nec)
        filt_mx = filt_mx + a * unfiltered_usp(k,nmx)
        filt_ee = filt_ee + a * unfiltered_usp(k,nee)
        !
      end do
      !
      usp(n0+n,nec) = filt_ec
      usp(n0+n,nmx) = filt_mx
      usp(n0+n,nee) = filt_ee
      !
    end do
    !
  end do
  !
end subroutine filter_usp_1d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine filter_usp_nd(mbeg,mend,usp)
  !
  !.. Use Statements ..
  use order_mod, only : maxSP
  use geovar,    only : ncell,cell
  !
  use thermo_mod, only : update_mixture_density
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: mbeg,mend
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,m,n,nc,np,n0,n1,n2
  integer  :: this_geom,this_order
  real(wp) :: a,val
  !
  !.. Local Arrays ..
  real(wp) :: unfiltered_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "filter_usp_nd"
  !
continue
  !
  ! Apply the filter to each cell
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,filter) &
  !$ACC&              PRIVATE(unfiltered_usp)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Save the current, unfiltered solution into the unfiltered_usp array.
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do m = mbeg,mend
      do n = 1,np
        unfiltered_usp(n,m) = usp(n0+n,m)
      end do
    end do
    !
    ! Filter the solution and save it back into the usp array.
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do m = mbeg,mend
      do n = 1,np
        !
        val = zero
        !
        do i = filter(this_geom,this_order)%exponential%col_ptr(n)+1, &
               filter(this_geom,this_order)%exponential%col_ptr(n+1)
          !
          k = filter(this_geom,this_order)%exponential%col_idx(i)
          a = filter(this_geom,this_order)%exponential%col_val(i)
          !
          val = val + a * unfiltered_usp(k,m)
          !
        end do
        !
        usp(n0+n,m) = val
        !
      end do
    end do
    !
  end do
  !
  ! Recompute the mixture density if this is a multi-species gas
  !
  if (ncs > 0) then
    !
    n1 = cell(1)%beg_sp
    n2 = cell(ncell)%end_sp
    !
    call update_mixture_density(n1,n2,usp)
    !
  end if
  !
end subroutine filter_usp_nd
!
!###############################################################################
!
function get_1D_exp_filter_coef(p,alpha,etac,s) result(return_value)
  !
  ! Routine to compute the 1D exponential filter coefficients.
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: p
  real(wp), intent(in) :: alpha
  real(wp), intent(in) :: etac
  integer,  intent(in) :: s
  !
  !.. Function Return Value ..
  real(wp) :: return_value(1:p+1)
  !
  !.. Local Scalars ..
  integer  :: i
  real(wp) :: rd,rn
  !
continue
  !
  ! Create the 1D exponential filter coefficients
  !
  if (p == 0) then
    !
    return_value(1) = one
    !
  else
    !
    do i = 0,p
      !
      rn = real(i,kind=wp)/real(p,kind=wp) - etac
      !
      if (rn > zero) then
        !
        rd = max( zero , one - etac )
        !
        if (rd == zero) then
          ! This protection from rd=0 shouldnt be needed. rd is only zero
          ! when etac>=1. However, when etac>=1, rn<=0, so we shouldnt be
          ! here when rd=0. Im keeping it here just in case.
          return_value(i+1) = exp(-alpha)
        else
          return_value(i+1) = exp(-alpha*(rn/rd)**s)
        end if
        !
      else
        !
        return_value(i+1) = one
        !
      end if
      !
    end do
    !
  end if
  !
end function get_1D_exp_filter_coef
!
!###############################################################################
!
function get_2P_quad_filter_coef(p,alpha,etac,s) result(return_value)
  !
  ! Routine to compute the exponential filter coefficients for a quad
  ! element based on the sum of the modal exponent, i.e., from [0:2*P]
  ! instead of [0:P].
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: p
  real(wp), intent(in) :: alpha
  real(wp), intent(in) :: etac
  integer,  intent(in) :: s
  !
  !.. Function Return Value ..
  real(wp) :: return_value(1:(p+1)**2)
  !
  !.. Local Scalars ..
  integer  :: i,j,l
  real(wp) :: rd,rn
  !
continue
  !
  ! Create the exponential filter coefficients
  !
  if (p == 0) then
    !
    return_value(1) = one
    !
  else
    !
    l = 0
    do j = 0,p
      do i = 0,p
        !
        l = l + 1
        !
        rn = real(i+j,kind=wp)/real(2*p,kind=wp) - etac
        !
        if (rn > zero) then
          !
          rd = max( zero , one - etac )
          !
          if (rd == zero) then
            ! This protection from rd=0 shouldnt be needed. rd is only zero
            ! when etac>=1. However, when etac>=1, rn<=0, so we shouldnt be
            ! here when rd=0. Im keeping it here just in case.
            return_value(l) = exp(-alpha)
          else
            return_value(l) = exp(-alpha*(rn/rd)**s)
          end if
          !
        else
          !
          return_value(l) = one
          !
        end if
        !
      end do
    end do
    !
  end if
  !
end function get_2P_quad_filter_coef
!
!###############################################################################
!
function get_3P_hex_filter_coef(p,alpha,etac,s) result(return_value)
  !
  ! Routine to compute the exponential filter coefficients for a hex
  ! element based on the sum of the modal exponent, i.e., from [0:3*P]
  ! instead of [0:P].
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: p
  real(wp), intent(in) :: alpha
  real(wp), intent(in) :: etac
  integer,  intent(in) :: s
  !
  !.. Function Return Value ..
  real(wp) :: return_value(1:(p+1)**3)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l
  real(wp) :: rd,rn
  !
continue
  !
  ! Create the exponential filter coefficients
  !
  if (p == 0) then
    !
    return_value(1) = one
    !
  else
    !
    l = 0
    do k = 0,p
      do j = 0,p
        do i = 0,p
          !
          l = l + 1
          !
          rn = real(i+j+k,kind=wp)/real(3*p,kind=wp) - etac
          !
          if (rn > zero) then
            !
            rd = max( zero , one - etac )
            !
            if (rd == zero) then
              ! This protection from rd=0 shouldnt be needed. rd is only zero
              ! when etac>=1. However, when etac>=1, rn<=0, so we shouldnt be
              ! here when rd=0. Im keeping it here just in case.
              return_value(l) = exp(-alpha)
            else
              return_value(l) = exp(-alpha*(rn/rd)**s)
            end if
            !
          else
            !
            return_value(l) = one
            !
          end if
          !
        end do
      end do
    end do
    !
  end if
  !
end function get_3P_hex_filter_coef
!
!###############################################################################
!
subroutine create_exp_filter_matrix(geom,p, &
                                    loc_filt_option,loc_filt_etac, &
                                    loc_filt_order,loc_filt_alpha, &
                                    filt,fmat)
  !
  ! Routine to compute an exponential filter matrix
  !
  !.. Use Statements ..
  use order_mod, only : geom_solpts
  use ovar,      only : Pyramid_Filter_Option
  !
  !.. Formal Arguments ..
  integer,                           intent(in) :: geom
  integer,                           intent(in) :: p
  integer,                           intent(in) :: loc_filt_option
  real(wp),                          intent(in) :: loc_filt_etac
  integer,                           intent(in) :: loc_filt_order
  real(wp),                          intent(in) :: loc_filt_alpha
  type(filt_coef_t), allocatable, intent(inout) :: filt(:)
  real(wp),          allocatable,   intent(out) :: fmat(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,ierr,nsp
  integer  :: s
  real(wp) :: alpha,etac
  !
  !.. Local Allocatable Arrays ..
  real(wp),          allocatable :: filt_diag(:)
  type(filt_coef_t), allocatable :: tmp_filt(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_exp_filter_matrix"
  !
continue
  !
  nsp = geom_solpts(geom,p)
  !
  etac = loc_filt_etac
  alpha = loc_filt_alpha
  s = loc_filt_order
  !
  ! Allocate the filt array if it has not yet been allocated,
  ! and compute the degree-P filter coefficients
  !
  if (.not. allocated(filt)) then
    !
    ! Allocate the derived type to store the filter coefficients
    !
    allocate ( filt(0:p) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"filt",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Compute the degree-P filter coefficients because these will
    ! always be needed. The rest of the filter coefficients will
    ! be special cases that will be computed as needed.
    !
    ! F2003 AUTO-REALLOC
    filt(p)%coef = get_1D_exp_filter_coef(p,alpha,etac,s)
    !
  end if
  !
  if (geom == Geom_Quad .and. loc_filt_option == 2) then
    !
    ! Special case of filter for quad element based on the sum of the modal
    ! exponents, i.e., for [0:2*P] instead of [0:P]
    !
    if (ubound(filt,dim=1) < 2*p) then
      !
      ! Expand filt
      !
      allocate ( tmp_filt(0:2*p) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"tmp_filt",1,__LINE__,__FILE__,ierr,error_message)
      !
      i = ubound(filt,dim=1)
      !
      tmp_filt(0:i) = filt(0:i)
      !
      call move_alloc( from=tmp_filt , to=filt )
      !
      ! F2003 AUTO-REALLOC
      filt(2*p)%coef = get_2P_quad_filter_coef(p,alpha,etac,s)
      !
    end if
    !
    filt_diag = filt(2*p)%coef(:) ! F2003 AUTO-REALLOC
    !
  else if (geom == Geom_Hexa .and. loc_filt_option == 2) then
    !
    ! Special case of filter for hex element based on the sum of the modal
    ! exponents, i.e., for [0:3*P] instead of [0:P]
    !
    if (ubound(filt,dim=1) < 3*p) then
      !
      ! Expand filt
      !
      allocate ( tmp_filt(0:3*p) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"tmp_filt",1,__LINE__,__FILE__,ierr,error_message)
      !
      i = ubound(filt,dim=1)
      !
      tmp_filt(0:i) = filt(0:i)
      !
      call move_alloc( from=tmp_filt , to=filt )
      !
      ! F2003 AUTO-REALLOC
      filt(3*p)%coef = get_3P_hex_filter_coef(p,alpha,etac,s)
      !
    end if
    !
    filt_diag = filt(3*p)%coef(:) ! F2003 AUTO-REALLOC
    !
  else if (geom == Geom_Pyra .and. Pyramid_Filter_Option /= 3) then
    !
    ! Special case of pyramid elements where the tensor products for k/=0
    ! layers use P-k 1D filter coefficients for the tensor product.
    !
    do i = 0,p-1
      if (.not. allocated(filt(i)%coef)) then
        if (Pyramid_Filter_Option == 1) then
          !
          ! This option keeps the cutoff order Pc relative to P for all i,
          ! where Pc = loc_filt_etac*real(P). To do this, we need to compute
          ! a new etac, etac_i, for i using
          !     etac_i = loc_filt_etac*real(p)/real(i)
          ! For smaller i, this may mean that Pc > i which would result in
          ! etac_i > 1, but etac_i should be in the interval [0:1]. For that
          ! reason, we add the min(1,max(0,...)) limiting on etac_i.
          !
          ! NOTE: Make sure we dont divide by zero when i=0. When i=0,
          !       we dont care about the value of etac because it doesnt
          !       affect the filter coefficient.
          !
          etac = loc_filt_etac*real(p,kind=wp)/real(max(1,i),kind=wp)
          etac = min( one , max(zero,etac) )
          !
          ! F2003 AUTO-REALLOC
          filt(i)%coef = get_1D_exp_filter_coef(i,alpha,etac,s)
          !
        else if (Pyramid_Filter_Option == 2) then
          !
          ! This option keeps etac at the same value loc_filt_etac for all i.
          ! This has the effect that the cutoff order Pc changes with each i.
          ! So if loc_filt_etac=0.5 and P=4, Pc takes the following values:
          !    i = P : Pc = etac*P = 0.5*4 = 2
          !    i = 3 : Pc = 1.5
          !    i = 2 : Pc = 1
          !    i = 1 : Pc = 0.5
          ! This results in a stronger filter compared to
          ! Pyramid_Filter_Option=1.
          !
          etac = loc_filt_etac
          !
          ! F2003 AUTO-REALLOC
          filt(i)%coef = get_1D_exp_filter_coef(i,alpha,etac,s)
          !
        end if
      end if
    end do
    !
    filt_diag = get_pyra_diag(filt(0:p)) ! F2003 AUTO-REALLOC
    !
  else
    !
    filt_diag = get_filt_diag(geom,filt(p)%coef) ! F2003 AUTO-REALLOC
    !
  end if
  ! 
  ! Finally, compute the filter matrix operator from the filter diagonal
  !
  fmat = compute_filter_matrix(geom,p,filt_diag) ! F2003 AUTO-REALLOC
  !
end subroutine create_exp_filter_matrix
!
!###############################################################################
!
function get_pyra_diag(filt) result(return_value)
  !
  !.. Use Statements ..
  use order_mod, only : geom_solpts
  !
  !.. Formal Arguments ..
  type(filt_coef_t), intent(in) :: filt(0:)
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,p,nsp
  !
continue
  !
  p = size(filt)-1
  nsp = geom_solpts(Geom_Pyra,p)
  !
  return_value = init_array(nsp,zero) ! F2003 AUTO-REALLOC
  !
  l = 0
  do k = 0,p
    do j = 0,p-k
      do i = 0,p-k
        l = l + 1
        return_value(l) = filt(p-k)%coef(i+1) * &
                          filt(p-k)%coef(j+1) * &
                          filt(p  )%coef(k+1)
      end do
    end do
  end do
  !
end function get_pyra_diag
!
!###############################################################################
!
function get_filt_diag(geom,filt_1D_coef) result(return_value)
  !
  !.. Use Statements ..
  use order_mod, only : geom_solpts
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: geom
  real(wp), intent(in) :: filt_1D_coef(:)
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,p,np,nsp
  !
continue
  !
  np = size(filt_1D_coef)
  p = np - 1
  nsp = geom_solpts(geom,p)
  !
  return_value = init_array(nsp,zero) ! F2003 AUTO-REALLOC
  !
  if (geom == Geom_Edge) then
    !
    do i = 1,np
      return_value(i) = filt_1D_coef(i)
    end do
    !
  else if (geom == Geom_Quad) then
    !
    l = 0
    do j = 1,np
      do i = 1,np
        l = l + 1
        return_value(l) = filt_1D_coef(i)*filt_1D_coef(j)
      end do
    end do
    !
  else if (geom == Geom_Tria) then
    !
    l = 0
    do j = 0,p
      do i = 0,p-j
        l = l + 1
        return_value(l) = filt_1D_coef(i+j+1)
      end do
    end do
    !
  else if (geom == Geom_Tetr) then
    !
    l = 0
    do k = 0,p
      do j = 0,p-k
        do i = 0,p-k-j
          l = l + 1
          return_value(l) = filt_1D_coef(i+j+k+1)
        end do
      end do
    end do
    !
  else if (geom == Geom_Pyra) then
    !
    ! NOTE: This is probably not the best way to create the pyramid
    !       filter coefficients. We should only get here if the input
    !       Pyramid_Filter_Option = 3, and it will be recommended to
    !       set Pyramid_Filter_Option to 1 (default) or 2.
    !
    l = 0
    do k = 0,p
      do j = 1,np-k
        do i = 1,np-k
          l = l + 1
          return_value(l) = filt_1D_coef(i)*filt_1D_coef(j)*filt_1D_coef(k+1)
        end do
      end do
    end do
    !
  else if (geom == Geom_Pris) then
    !
    l = 0
    do k = 1,np
      do j = 0,p
        do i = 0,p-j
          l = l + 1
          return_value(l) = filt_1D_coef(i+j+1)*filt_1D_coef(k)
        end do
      end do
    end do
    !
  else if (geom == Geom_Hexa) then
    !
    l = 0
    do k = 1,np
      do j = 1,np
        do i = 1,np
          l = l + 1
          return_value(l) = filt_1D_coef(i)*filt_1D_coef(j)*filt_1D_coef(k)
        end do
      end do
    end do
    !
  end if
  !
end function get_filt_diag
!
!###############################################################################
!
function compute_filter_matrix(geom,p,filt_diag) result(return_value)
  !
  ! Function to compute the final filter matrix operator
  ! provided the filter diagonal.
  !
  !.. Use Statements ..
  use vandermonde_mod, only : vand
#ifdef SPECIAL_FOR_INTEL
  use blas95, only : gemm
#endif
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: geom
  integer,  intent(in) :: p
  real(wp), intent(in) :: filt_diag(:)
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,nsp
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: fmat(:,:)
  !
continue
  !
  nsp = size(filt_diag,dim=1)
  !
#ifdef SPECIAL_FOR_INTEL
  !
  fmat = init_array(nsp,nsp,zero) ! F2003 AUTO-REALLOC
  return_value = init_array(nsp,nsp,zero) ! F2003 AUTO-REALLOC
  !
  ! Perform the matrix multiplication of (V x F) where F is the diagonal matrix
  ! with the filter coefficients on the diagonal.
  !
  do j = 1,nsp
    do i = 1,nsp
      fmat(i,j) = vand(geom,p)%modal2nodal%mat(i,j)*filt_diag(j)
    end do
  end do
  !
  ! Finally, multiply (V x F) x V_inv
  !
  call gemm(fmat,vand(geom,p)%nodal2modal%mat,return_value)
  !
#else
  !
  ! Perform the matrix multiplication of (V x F) where F is the diagonal matrix
  ! with the filter coefficients on the diagonal.
  !
  ! F2003 AUTO-REALLOC
  fmat = vand(geom,p)%modal2nodal%mat(:,:) * spread(filt_diag,dim=1,ncopies=nsp)
  !
  ! Finally, multiply (V x F) x V_inv
  !
  ! F2003 AUTO-REALLOC
  return_value = matmul( fmat(:,:) , vand(geom,p)%nodal2modal%mat(:,:) )
  !
#endif
  !
end function compute_filter_matrix
!
!###############################################################################
!
function create_cutoff_filter(geom,p,cutoff_order) result(return_value)
  !
  ! Routine to compute the cutoff filter matrix of degree P
  ! that cuts off all modes of degree cutoff_order and higher.
  !
  !.. Use Statements ..
  use order_mod, only : geom_solpts
  !
  !.. Formal Arguments ..
  integer, intent(in) :: geom
  integer, intent(in) :: p
  integer, intent(in) :: cutoff_order
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,nsp
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: filt_diag(:)
  !
continue
  !
  nsp = geom_solpts(geom,p)
  !
  ! Initialize the filter diagonal to all ones. We will go through and change
  ! any locations on the diagonal to zero that need to be cutoff.
  !
  filt_diag = init_array(nsp,one) ! F2003 AUTO-REALLOC
  !
  if (geom == Geom_Edge) then
    !
    do i = cutoff_order,p
      filt_diag(i+1) = zero
    end do
    !
  else if (geom == Geom_Quad) then
    !
    l = 0
    do j = 0,p
      do i = 0,p
        l = l + 1
        if (any([i,j] >= cutoff_order)) then
          filt_diag(l) = zero
        end if
      end do
    end do
    !
  else if (geom == Geom_Tria) then
    !
    l = 0
    do j = 0,p
      do i = 0,p-j
        l = l + 1
        if (i+j >= cutoff_order) then
          filt_diag(l) = zero
        end if
      end do
    end do
    !
  else if (geom == Geom_Tetr) then
    !
    l = 0
    do k = 0,p
      do j = 0,p-k
        do i = 0,p-k-j
          l = l + 1
          if (i+j+k >= cutoff_order) then
            filt_diag(l) = zero
          end if
        end do
      end do
    end do
    !
  else if (geom == Geom_Pyra) then
    !
    l = 0
    do k = 0,p
      do j = 0,p-k
        do i = 0,p-k
          l = l + 1
          if (any([i,j] >= cutoff_order-k .or. k >= cutoff_order)) then
            filt_diag(l) = zero
          end if
        end do
      end do
    end do
    !
  else if (geom == Geom_Pris) then
    !
    l = 0
    do k = 0,p
      do j = 0,p
        do i = 0,p-j
          l = l + 1
          if (any([i+j,k] >= cutoff_order)) then
            filt_diag(l) = zero
          end if
        end do
      end do
    end do
    !
  else if (geom == Geom_Hexa) then
    !
    l = 0
    do k = 0,p
      do j = 0,p
        do i = 0,p
          l = l + 1
          if (any([i,j,k] >= cutoff_order)) then
            filt_diag(l) = zero
          end if
        end do
      end do
    end do
    !
  end if
  ! 
  ! Finally, compute the filter matrix operator from the filter diagonal
  !
  return_value = compute_filter_matrix(geom,p,filt_diag) ! F2003 AUTO-REALLOC
  !
end function create_cutoff_filter
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function FilterMatrix_Edge(vand,vand_inv,cutoff_order) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Optional Arguments ..
  integer,  optional, intent(in) :: cutoff_order
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,nsp,order
  integer  :: Filter_Cutoff,ierr
  integer  :: local_Filter_Option
  integer  :: filt_order
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
  real(qp) :: feps_qp
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: filter_1D_qp(:)
  real(qp), allocatable :: fmat_qp(:,:)
  !
continue
  !
  ! Get the number of solution points
  !
  nsp = size(vand,dim=1)
  order = nsp-1
  !
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  !
  filt_order = Filter_Order
  !
  if (present(cutoff_order)) then
    if (cutoff_order >= 0) then
      local_Filter_Option = Cutoff_Filter_Option
      Filter_Cutoff = cutoff_order
      fetac_qp = real(cutoff_order,kind=qp) / real(order,kind=qp)
      feps_qp = qzero
    else if (cutoff_order == -2) then
      fetac_qp = real( Sponge_Filter_etac , kind=qp )
      Filter_Cutoff = int_cutoff_function( &
                        max(zero,min(one,Sponge_Filter_etac)) * &
                        real(order,kind=wp) )
      local_Filter_Option = 1
      feps_qp = qone
      falpha_qp = real( Sponge_Filter_Alpha , kind=qp )
      filt_order = Sponge_Filter_Order
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = int_cutoff_function( &
                        max(zero,min(one,Limiter_Filter_etac)) * &
                        real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp = real( Filter_etac , kind=qp )
    Filter_Cutoff = int_cutoff_function( &
                      max(zero,min(one,Filter_etac)) * &
                      real(order,kind=wp) )
    local_Filter_Option = Filter_Option
    feps_qp = real( one - Filter_eps , kind=qp )
  end if
  !
  ! Allocate the return array containing the working
  ! precision filter matrix
  !
  allocate ( return_value(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  ! Create the filter matrix using quad precision
  !
  allocate ( fmat_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  !
  ! Initialize the local filter matrices to identity matrices
  !
  do i = 1,nsp
    fmat_qp(i,i) = qone
  end do
  !
  ! Fill the diagonal of the filter matrix with the filter function
  !
  if (local_Filter_Option /= 0) then
    !
    do i = 0,order
      if (i >= Filter_Cutoff) then
        !
        rn_qp = real(i,kind=qp)/real(order,kind=qp) - fetac_qp
        rd_qp = max( qzero , qone - fetac_qp )
        !
        if (rd_qp == qzero) then
          fmat_qp(i+1,i+1) = feps_qp
        else
          fmat_qp(i+1,i+1) = feps_qp*exp(-falpha_qp*(rn_qp/rd_qp)**filt_order)
        end if
        !
      end if
    end do
    !
  end if
  !
  ! Compute the final filter matrix for quadrilaterals
  !
  fmat_qp = matmul( real(vand,kind=qp) , fmat_qp )
  fmat_qp = matmul( fmat_qp , real(vand_inv,kind=qp) )
  !
  ! Use the chop function to copy fmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( fmat_qp )
  !
end function FilterMatrix_Edge
!
!###############################################################################
!
function FilterMatrix_Edge_with_output(vand,vand_inv,cutoff_order) &
                                result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Optional Arguments ..
  integer,  optional, intent(in) :: cutoff_order
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,nsp,order
  integer  :: Filter_Cutoff,ierr
  integer  :: local_Filter_Option
  integer  :: filt_order
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
  real(qp) :: feps_qp
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: filter_1D_qp(:)
  real(qp), allocatable :: fmat_qp(:,:)
  !
continue
  !
  ! Get the number of solution points
  !
  nsp = size(vand,dim=1)
  order = nsp-1
  !
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  !
  filt_order = Filter_Order
  !
  if (present(cutoff_order)) then
    if (cutoff_order >= 0) then
      local_Filter_Option = Cutoff_Filter_Option
      Filter_Cutoff = cutoff_order
      fetac_qp = real(cutoff_order,kind=qp) / real(order,kind=qp)
      feps_qp = qzero
    else if (cutoff_order == -2) then
      fetac_qp = real( Sponge_Filter_etac , kind=qp )
      Filter_Cutoff = int_cutoff_function( &
                        max(zero,min(one,Sponge_Filter_etac)) * &
                        real(order,kind=wp) )
      local_Filter_Option = 1
      feps_qp = qone
      falpha_qp = real( Sponge_Filter_Alpha , kind=qp )
      filt_order = Sponge_Filter_Order
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = int_cutoff_function( &
                        max(zero,min(one,Limiter_Filter_etac)) * &
                        real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp = real( Filter_etac , kind=qp )
    Filter_Cutoff = int_cutoff_function( &
                      max(zero,min(one,Filter_etac)) * &
                      real(order,kind=wp) )
    local_Filter_Option = Filter_Option
    feps_qp = real( one - Filter_eps , kind=qp )
  end if
  !
  ! Allocate the return array containing the working
  ! precision filter matrix
  !
  allocate ( return_value(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  ! Create the filter matrix using quad precision
  !
  allocate ( fmat_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  !
  ! Initialize the local filter matrices to identity matrices
  !
  do i = 1,nsp
    fmat_qp(i,i) = qone
  end do
  !
  ! Fill the diagonal of the filter matrix with the filter function
  !
  if (local_Filter_Option /= 0) then
    !
    do i = 0,order
      if (i >= Filter_Cutoff) then
        !
        rn_qp = real(i,kind=qp)/real(order,kind=qp) - fetac_qp
        rd_qp = max( qzero , qone - fetac_qp )
        !
        if (rd_qp == qzero) then
          fmat_qp(i+1,i+1) = feps_qp
        else
          fmat_qp(i+1,i+1) = feps_qp*exp(-falpha_qp*(rn_qp/rd_qp)**filt_order)
        end if
        !
      end if
    end do
    !
    if (mypnum == glb_root) then
      if (present(cutoff_order)) then
        if (cutoff_order == -1) then
          write (io,102)
          102 format (/,"Diagonal for edge limiter filter")
        else
          write (io,101) order+1-(Filter_cutoff)
          101 format (/,"Diagonal for edge cutoff filter for P-",i0)
        end if
      else
        write (io,103)
        103 format (/,"Diagonal for edge standard filter")
      end if
      write (io,100) (fmat_qp(i,i),i=1,nsp)
      100 format (100(2x,es14.6))
    end if
    !
  end if
  !
  ! Compute the final filter matrix for quadrilaterals
  !
  fmat_qp = matmul( real(vand,kind=qp) , fmat_qp )
  fmat_qp = matmul( fmat_qp , real(vand_inv,kind=qp) )
  !
  ! Use the chop function to copy fmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( fmat_qp )
  !
end function FilterMatrix_Edge_with_output
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function FilterMatrix_Quad(vand,vand_inv,cutoff_order) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: cutoff_order
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,np,nsp,order
  integer  :: Filter_Cutoff,ierr
  integer  :: local_Filter_Option
  integer  :: filt_order
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
  real(qp) :: feps_qp
  !
  !.. Local Arrays ..
  real(qp), allocatable :: filter_1D_qp(:)
  real(qp), allocatable :: fmat_qp(:,:)
  !
continue
  !
  ! Get the number of solution points
  !
  nsp = size(vand,dim=1)
  order = solpts2order_quad(nsp)
  np = order+1
 !np = nint( sqrt(real(nsp,kind=wp)) )
 !order = np-1
  !
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  !
  filt_order = Filter_Order
  !
  if (present(cutoff_order)) then
    if (cutoff_order >= 0) then
      local_Filter_Option = Cutoff_Filter_Option
      i = merge(2,1,local_Filter_Option==2)
      Filter_Cutoff = i*cutoff_order
      fetac_qp = real(Filter_Cutoff,kind=qp) / real(i*order,kind=qp)
      feps_qp = qzero
    else if (cutoff_order == -2) then
      fetac_qp = real( Sponge_Filter_etac , kind=qp )
      Filter_Cutoff = int_cutoff_function( &
                        max(zero,min(one,Sponge_Filter_etac)) * &
                        real(order,kind=wp) )
      local_Filter_Option = 1
      feps_qp = qone
      falpha_qp = real( Sponge_Filter_Alpha , kind=qp )
      filt_order = Sponge_Filter_Order
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = int_cutoff_function( &
                        max(zero,min(one,Limiter_Filter_etac)) * &
                        real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp  = real( Filter_etac , kind=qp )
    Filter_Cutoff = int_cutoff_function( &
                      max(zero,min(one,Filter_etac)) * &
                      real(order,kind=wp) )
    local_Filter_Option = Filter_Option
    feps_qp = real( one - Filter_eps , kind=qp )
  end if
  !
  ! Allocate the return array containing the working
  ! precision filter matrix
  !
  allocate ( return_value(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  ! Create the filter matrix using quad precision
  !
  allocate ( fmat_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  !
  ! Initialize the local filter matrices to identity matrices
  !
  do i = 1,nsp
    fmat_qp(i,i) = qone
  end do
  !
  ! Compute the filter diagonal for quadrilaterals
  !
  if (local_Filter_Option == 1) then
    !
    ! Option 1 - Tensor product of the 1D filter
    !           (Equivalent to using a 1D filter in each direction)
    !
    ! Allocate the 1D filter vector and set all its values to one
    !
    allocate ( filter_1D_qp(1:np) , source=qone , stat=ierr )
    !
    ! First get the 1D filter diagonal
    !
    do i = 0,order
      if (i >= Filter_Cutoff) then
        !
        rn_qp = real(i,kind=qp)/real(order,kind=qp) - fetac_qp
        rd_qp = max( qzero , qone - fetac_qp )
        !
        if (rd_qp == qzero) then
          filter_1D_qp(i+1) = feps_qp
        else
          filter_1D_qp(i+1) = feps_qp * &
                              exp(-falpha_qp * (rn_qp/rd_qp)**filt_order)
        end if
        !
      end if
    end do
    !
    ! Now fill the filter diagonal for quads
    !
    do j = 0,order
      do i = 0,order
        !
        k = np*j + i + 1
        !
        fmat_qp(k,k) = filter_1D_qp(i+1)*filter_1D_qp(j+1)
        !
      end do
    end do
    !
  else if (local_Filter_Option == 2) then
    !
    ! Option 2 - The quad basis functions actually contains cross polynomial
    !            terms up to order 2*p (as opposed to just 1*p for triangles).
    !            This filter applies the cut-off frequency based on 2*p as
    !            opposed to basing it on just p in 1D like Option 1.
    !
    do i = 0,order
      do j = 0,order
        !
        if (i+j >= Filter_Cutoff) then
          !
          k = np*j + i + 1
          !
          rn_qp = real(i+j-Filter_Cutoff,kind=qp)
          rd_qp = max( qzero , real(2*order-Filter_Cutoff,kind=qp) )
          !
          if (rd_qp == qzero) then
            fmat_qp(k,k) = feps_qp
          else
            fmat_qp(k,k) = feps_qp*exp(-falpha_qp * (rn_qp/rd_qp)**filt_order)
          end if
          !
        end if
        !
      end do
    end do
    !
  end if
  !
  ! Compute the final filter matrix for quadrilaterals
  !
  fmat_qp = matmul( real(vand,kind=qp) , fmat_qp )
  fmat_qp = matmul( fmat_qp , real(vand_inv,kind=qp) )
  !
  ! Use the chop function to copy fmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( fmat_qp )
  !
end function FilterMatrix_Quad
!
!###############################################################################
!
function FilterMatrix_Quad_with_output(vand,vand_inv,cutoff_order) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: cutoff_order
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,np,nsp,order
  integer  :: Filter_Cutoff,ierr
  integer  :: local_Filter_Option
  integer  :: filt_order
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
  real(qp) :: feps_qp
  !
  !.. Local Arrays ..
  real(qp), allocatable :: filter_1D_qp(:)
  real(qp), allocatable :: fmat_qp(:,:)
  !
continue
  !
  ! Get the number of solution points
  !
  nsp = size(vand,dim=1)
  order = solpts2order_quad(nsp)
  np = order+1
 !np = nint( sqrt(real(nsp,kind=wp)) )
 !order = np-1
  !
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  !
  filt_order = Filter_Order
  !
  if (present(cutoff_order)) then
    if (cutoff_order >= 0) then
      local_Filter_Option = Cutoff_Filter_Option
      i = merge(2,1,local_Filter_Option==2)
      Filter_Cutoff = i*cutoff_order
      fetac_qp = real(Filter_Cutoff,kind=qp) / real(i*order,kind=qp)
      feps_qp = qzero
    else if (cutoff_order == -2) then
      fetac_qp = real( Sponge_Filter_etac , kind=qp )
      Filter_Cutoff = int_cutoff_function( &
                        max(zero,min(one,Sponge_Filter_etac)) * &
                        real(order,kind=wp) )
      local_Filter_Option = 1
      feps_qp = qone
      falpha_qp = real( Sponge_Filter_Alpha , kind=qp )
      filt_order = Sponge_Filter_Order
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = int_cutoff_function( &
                        max(zero,min(one,Limiter_Filter_etac)) * &
                        real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp  = real( Filter_etac , kind=qp )
    Filter_Cutoff = int_cutoff_function( &
                      max(zero,min(one,Filter_etac)) * &
                      real(order,kind=wp) )
    local_Filter_Option = Filter_Option
    feps_qp = real( one - Filter_eps , kind=qp )
  end if
  !
  ! Allocate the return array containing the working
  ! precision filter matrix
  !
  allocate ( return_value(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  ! Create the filter matrix using quad precision
  !
  allocate ( fmat_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  !
  ! Initialize the local filter matrices to identity matrices
  !
  do i = 1,nsp
    fmat_qp(i,i) = qone
  end do
  !
  ! Compute the filter diagonal for quadrilaterals
  !
  if (local_Filter_Option == 1) then
    !
    ! Option 1 - Tensor product of the 1D filter
    !           (Equivalent to using a 1D filter in each direction)
    !
    ! Allocate the 1D filter vector and set all its values to one
    !
    allocate ( filter_1D_qp(1:np) , source=qone , stat=ierr )
    !
    ! First get the 1D filter diagonal
    !
    do i = 0,order
      if (i >= Filter_Cutoff) then
        !
        rn_qp = real(i,kind=qp)/real(order,kind=qp) - fetac_qp
        rd_qp = max( qzero , qone - fetac_qp )
        !
        if (rd_qp == qzero) then
          filter_1D_qp(i+1) = feps_qp
        else
          filter_1D_qp(i+1) = feps_qp * &
                              exp(-falpha_qp * (rn_qp/rd_qp)**filt_order)
        end if
        !
      end if
    end do
    !
    ! Now fill the filter diagonal for quads
    !
    do j = 0,order
      do i = 0,order
        !
        k = np*j + i + 1
        !
        fmat_qp(k,k) = filter_1D_qp(i+1)*filter_1D_qp(j+1)
        !
      end do
    end do
    !
  else if (local_Filter_Option == 2) then
    !
    ! Option 2 - The quad basis functions actually contains cross polynomial
    !            terms up to order 2*p (as opposed to just 1*p for triangles).
    !            This filter applies the cut-off frequency based on 2*p as
    !            opposed to basing it on just p in 1D like Option 1.
    !
    do i = 0,order
      do j = 0,order
        !
        if (i+j >= Filter_Cutoff) then
          !
          k = np*j + i + 1
          !
          rn_qp = real(i+j-Filter_Cutoff,kind=qp)
          rd_qp = max( qzero , real(2*order-Filter_Cutoff,kind=qp) )
          !
          if (rd_qp == qzero) then
            fmat_qp(k,k) = feps_qp
          else
            fmat_qp(k,k) = feps_qp*exp(-falpha_qp * (rn_qp/rd_qp)**filt_order)
          end if
          !
        end if
        !
      end do
    end do
    !
  end if
  !
  if (local_Filter_Option /= 0) then
    !
    if (mypnum == glb_root) then
      if (present(cutoff_order)) then
        if (cutoff_order == -1) then
          write (io,102)
          102 format (/,"Diagonal for quad limiter filter")
        else
          write (io,104) order+1-(Filter_cutoff)
          104 format (/,"1D Diagonal for quad cutoff filter for P-",i0)
          write (io,105) (filter_1d_qp(i),i=1,order+1)
          105 format (100(2x,es14.6,:))
          write (io,101) order+1-(Filter_cutoff)
          101 format (/,"Diagonal for quad cutoff filter for P-",i0)
        end if
      else
        write (io,103)
        103 format (/,"Diagonal for quad standard filter")
      end if
      do j = 0,order
        do i = 0,order
          k = np*j + i + 1
          write (io,100) i,j,fmat_qp(k,k)
        end do
      end do
      100 format (2x,i2,1x,i2,4x,es14.6)
    end if
    !
  end if
  !
  ! Compute the final filter matrix for quadrilaterals
  !
  fmat_qp = matmul( real(vand,kind=qp) , fmat_qp )
  fmat_qp = matmul( fmat_qp , real(vand_inv,kind=qp) )
  !
  ! Use the chop function to copy fmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( fmat_qp )
  !
end function FilterMatrix_Quad_with_output
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function FilterMatrix_Tria(vand,vand_inv,cutoff_order) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: cutoff_order
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,np,nsp,order
  integer  :: Filter_Cutoff,ierr
  integer  :: local_Filter_Option
  integer  :: filt_order
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
  real(qp) :: feps_qp
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: fmat_qp(:,:)
  !
continue
  !
  ! Get the number of solution points
  !
  nsp = size(vand,dim=1)
  order = solpts2order_tria(nsp)
  np = order+1
  !
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  !
  filt_order = Filter_Order
  !
  if (present(cutoff_order)) then
    if (cutoff_order >= 0) then
      local_Filter_Option = 1
      Filter_Cutoff = cutoff_order
      fetac_qp = real(Filter_Cutoff,kind=qp) / real(order,kind=qp)
      feps_qp = qzero
    else if (cutoff_order == -2) then
      fetac_qp = real( Sponge_Filter_etac , kind=qp )
      Filter_Cutoff = int_cutoff_function( &
                        max(zero,min(one,Sponge_Filter_etac)) * &
                        real(order,kind=wp) )
      local_Filter_Option = 1
      feps_qp = qone
      falpha_qp = real( Sponge_Filter_Alpha , kind=qp )
      filt_order = Sponge_Filter_Order
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = int_cutoff_function( &
                        max(zero,min(one,Limiter_Filter_etac)) * &
                        real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp  = real( Filter_etac , kind=qp )
    Filter_Cutoff = int_cutoff_function( &
                      max(zero,min(one,Filter_etac)) * &
                      real(order,kind=wp) )
    local_Filter_Option = Filter_Option
    feps_qp = real( one - Filter_eps , kind=qp )
  end if
  !
  ! Allocate the return array containing the working
  ! precision filter matrix
  !
  allocate ( return_value(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  ! Create the filter matrix using quad precision
  !
  allocate ( fmat_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  !
  ! Initialize the local filter matrices to identity matrices
  !
  do i = 1,nsp
    fmat_qp(i,i) = qone
  end do
  !
  ! Compute the filter diagonal for triangles
  !
  if (local_Filter_Option /= 0 .and. Filter_Cutoff < order) then
    !
    k = 1
    do i = 0,order
      do j = 0,order-i
        !
        if (i+j >= Filter_Cutoff) then
          !
          rn_qp = real(i+j,kind=qp)/real(order,kind=qp) - fetac_qp
          rd_qp = max( qzero , qone - fetac_qp )
          !
          if (rd_qp == qzero) then
            fmat_qp(k,k) = feps_qp
          else
            fmat_qp(k,k) = feps_qp*exp(-falpha_qp * (rn_qp/rd_qp)**filt_order)
          end if
          !
        end if
        !
        k = k + 1
        !
      end do
    end do
    !
  end if
  !
  ! Compute the final filter matrix for triangles
  !
  fmat_qp = matmul( real(vand,kind=qp) , fmat_qp )
  fmat_qp = matmul( fmat_qp , real(vand_inv,kind=qp) )
  !
  ! Use the chop function to copy fmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( fmat_qp )
  !
end function FilterMatrix_Tria
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function FilterMatrix_Tetr(vand,vand_inv,cutoff_order) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: cutoff_order
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,np,nsp,order
  integer  :: Filter_Cutoff,ierr
  integer  :: local_Filter_Option
  integer  :: filt_order
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
  real(qp) :: feps_qp
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: fmat_qp(:,:)
  !
continue
  !
  ! Get the number of solution points
  !
  nsp = size(vand,dim=1)
  order = solpts2order_tetr(nsp)
  np = order+1
  !
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  !
  filt_order = Filter_Order
  !
  if (present(cutoff_order)) then
    if (cutoff_order >= 0) then
     !local_Filter_Option = Cutoff_Filter_Option
      local_Filter_Option = 1
      Filter_Cutoff = cutoff_order
      fetac_qp = real(Filter_Cutoff,kind=qp) / real(order,kind=qp)
      feps_qp = qzero
    else if (cutoff_order == -2) then
      fetac_qp = real( Sponge_Filter_etac , kind=qp )
      Filter_Cutoff = int_cutoff_function( &
                        max(zero,min(one,Sponge_Filter_etac)) * &
                        real(order,kind=wp) )
      local_Filter_Option = 1
      feps_qp = qone
      falpha_qp = real( Sponge_Filter_Alpha , kind=qp )
      filt_order = Sponge_Filter_Order
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = int_cutoff_function( &
                        max(zero,min(one,Limiter_Filter_etac)) * &
                        real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp  = real( Filter_etac , kind=qp )
    Filter_Cutoff = int_cutoff_function( &
                      max(zero,min(one,Filter_etac)) * &
                      real(order,kind=wp) )
    local_Filter_Option = Filter_Option
    feps_qp = real( one - Filter_eps , kind=qp )
  end if
  !
  ! Allocate the return array containing the working
  ! precision filter matrix
  !
  allocate ( return_value(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  ! Create the filter matrix using quad precision
  !
  allocate ( fmat_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  !
  ! Initialize the local filter matrices to identity matrices
  !
  do i = 1,nsp
    fmat_qp(i,i) = qone
  end do
  !
  ! Compute the filter diagonal for tetrahedra
  !
  if (local_Filter_Option /= 0 .and. Filter_Cutoff < order) then
    !
    l = 1
    do i = 0,order
      do j = 0,order-i
        do k = 0,order-i-j
          !
          if (i+j+k >= Filter_Cutoff) then
            !
            rn_qp = real(i+j+k,kind=qp)/real(order,kind=qp) - fetac_qp
            rd_qp = max( qzero , qone - fetac_qp )
            !
            if (rd_qp == qzero) then
              fmat_qp(l,l) = feps_qp
            else
              fmat_qp(l,l) = feps_qp*exp(-falpha_qp * (rn_qp/rd_qp)**filt_order)
            end if
            !
          end if
          !
          l = l + 1
          !
        end do
      end do
    end do
    !
  end if
  !
  ! Compute the final filter matrix for tetrahedra
  !
  fmat_qp = matmul( real(vand,kind=qp) , fmat_qp )
  fmat_qp = matmul( fmat_qp , real(vand_inv,kind=qp) )
  !
  ! Use the chop function to copy fmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( fmat_qp )
  !
end function FilterMatrix_Tetr
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function FilterMatrix_Pyra(vand,vand_inv,cutoff_order) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: cutoff_order
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,np,nsp,order
  integer  :: Filter_Cutoff,ierr
  integer  :: local_Filter_Option
  integer  :: filt_order
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
  real(qp) :: feps_qp
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: fmat_qp(:,:)
  !
continue
  !
  ! Get the number of solution points
  !
  nsp = size(vand,dim=1)
  order = solpts2order_pyra(nsp)
  np = order+1
  !
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  !
  filt_order = Filter_Order
  !
  if (present(cutoff_order)) then
    if (cutoff_order >= 0) then
     !local_Filter_Option = Cutoff_Filter_Option
      local_Filter_Option = 1
      Filter_Cutoff = cutoff_order
      fetac_qp = real(Filter_Cutoff,kind=qp) / real(order,kind=qp)
      feps_qp = qzero
    else if (cutoff_order == -2) then
      fetac_qp = real( Sponge_Filter_etac , kind=qp )
      Filter_Cutoff = int_cutoff_function( &
                        max(zero,min(one,Sponge_Filter_etac)) * &
                        real(order,kind=wp) )
      local_Filter_Option = 1
      feps_qp = qone
      falpha_qp = real( Sponge_Filter_Alpha , kind=qp )
      filt_order = Sponge_Filter_Order
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = int_cutoff_function( &
                        max(zero,min(one,Limiter_Filter_etac)) * &
                        real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp  = real( Filter_etac , kind=qp )
    Filter_Cutoff = int_cutoff_function( &
                      max(zero,min(one,Filter_etac)) * &
                      real(order,kind=wp) )
    local_Filter_Option = Filter_Option
    feps_qp = real( one - Filter_eps , kind=qp )
  end if
  !
  ! Allocate the return array containing the working
  ! precision filter matrix
  !
  allocate ( return_value(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  ! Create the filter matrix using quad precision
  !
  allocate ( fmat_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  !
  ! Initialize the local filter matrices to identity matrices
  !
  do i = 1,nsp
    fmat_qp(i,i) = qone
  end do
  !
  ! Compute the filter diagonal for pyramids
  !
  if (local_Filter_Option /= 0 .and. Filter_Cutoff < order) then
    !
    l = 1
    do k = 0,order
      do j = 0,order-k
        do i = 0,order-k
          !
          if (i+j+k >= Filter_Cutoff) then
            !
            rn_qp = real(i+j+k,kind=qp)/real(order,kind=qp) - fetac_qp
            rd_qp = max( qzero , qone - fetac_qp )
            !
            if (rd_qp == qzero) then
              fmat_qp(l,l) = feps_qp
            else
              fmat_qp(l,l) = feps_qp*exp(-falpha_qp * (rn_qp/rd_qp)**filt_order)
            end if
            !
          end if
          !
          l = l + 1
          !
        end do
      end do
    end do
    !
  end if
  !
  ! Compute the final filter matrix for pyramids
  !
  fmat_qp = matmul( real(vand,kind=qp) , fmat_qp )
  fmat_qp = matmul( fmat_qp , real(vand_inv,kind=qp) )
  !
  ! Use the chop function to copy fmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( fmat_qp )
  !
end function FilterMatrix_Pyra
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function FilterMatrix_Pris(vand,vand_inv,cutoff_order) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: cutoff_order
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,np,nsp,order
  integer  :: Filter_Cutoff,ierr
  integer  :: local_Filter_Option
  integer  :: filt_order
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
  real(qp) :: feps_qp
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: fmat_qp(:,:)
  !
continue
  !
  ! Get the number of solution points
  !
  nsp = size(vand,dim=1)
  order = solpts2order_pris(nsp)
  np = order+1
  !
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  !
  filt_order = Filter_Order
  !
  if (present(cutoff_order)) then
    if (cutoff_order >= 0) then
     !local_Filter_Option = Cutoff_Filter_Option
      local_Filter_Option = 1
      Filter_Cutoff = cutoff_order
      fetac_qp = real(Filter_Cutoff,kind=qp) / real(order,kind=qp)
      feps_qp = qzero
    else if (cutoff_order == -2) then
      fetac_qp = real( Sponge_Filter_etac , kind=qp )
      Filter_Cutoff = int_cutoff_function( &
                        max(zero,min(one,Sponge_Filter_etac)) * &
                        real(order,kind=wp) )
      local_Filter_Option = 1
      feps_qp = qone
      falpha_qp = real( Sponge_Filter_Alpha , kind=qp )
      filt_order = Sponge_Filter_Order
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = int_cutoff_function( &
                        max(zero,min(one,Limiter_Filter_etac)) * &
                        real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp  = real( Filter_etac , kind=qp )
    Filter_Cutoff = int_cutoff_function( &
                      max(zero,min(one,Filter_etac)) * &
                      real(order,kind=wp) )
    local_Filter_Option = Filter_Option
    feps_qp = real( one - Filter_eps , kind=qp )
  end if
  !
  ! Allocate the return array containing the working
  ! precision filter matrix
  !
  allocate ( return_value(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  ! Create the filter matrix using quad precision
  !
  allocate ( fmat_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  !
  ! Initialize the local filter matrices to identity matrices
  !
  do i = 1,nsp
    fmat_qp(i,i) = qone
  end do
  !
  ! Compute the filter diagonal for prisms
  !
  if (local_Filter_Option /= 0 .and. Filter_Cutoff < order) then
    !
    l = 1
    do k = 0,order
      do i = 0,order
        do j = 0,order-i
          !
          if (i+j+k >= Filter_Cutoff) then
            !
            rn_qp = real(i+j+k,kind=qp)/real(order,kind=qp) - fetac_qp
            rd_qp = max( qzero , qone - fetac_qp )
            !
            if (rd_qp == qzero) then
              fmat_qp(l,l) = feps_qp
            else
              fmat_qp(l,l) = feps_qp*exp(-falpha_qp * (rn_qp/rd_qp)**filt_order)
            end if
            !
          end if
          !
          l = l + 1
          !
        end do
      end do
    end do
    !
  end if
  !
  ! Compute the final filter matrix for prisms
  !
  fmat_qp = matmul( real(vand,kind=qp) , fmat_qp )
  fmat_qp = matmul( fmat_qp , real(vand_inv,kind=qp) )
  !
  ! Use the chop function to copy fmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( fmat_qp )
  !
end function FilterMatrix_Pris
!
!###############################################################################
!
!!!#define TIME_FILTER_MATMUL
#ifndef TIME_FILTER_MATMUL
pure &
#endif
function FilterMatrix_Hexa(vand,vand_inv,cutoff_order) result(return_value)
  !
#ifdef SPECIAL_FOR_INTEL
  !.. Use Statements ..
  use blas95, only : gemm
#endif
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: cutoff_order
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,np,nsp,order
  integer  :: Filter_Cutoff,ierr
  integer  :: local_Filter_Option
  integer  :: filt_order
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
  real(qp) :: feps_qp
  real(wp) :: a
  integer  :: num_zeros
  !
  !.. Local Arrays ..
  real(qp), allocatable :: filter_1D_qp(:)
  real(qp), allocatable :: filt_diag_qp(:)
  real(qp), allocatable :: fmat_qp(:,:)
  real(wp), allocatable :: filt_diag(:)
  real(wp), allocatable :: fmat(:,:)
  !
continue
  !
  ! Get the number of solution points
  !
  nsp = size(vand,dim=1)
  order = solpts2order_hexa(nsp)
  np = order+1
 !np = nint( real(nsp,kind=wp)**one3 )
 !order = np-1
  !
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  !
  filt_order = Filter_Order
  !
  if (present(cutoff_order)) then
    if (cutoff_order >= 0) then
      local_Filter_Option = Cutoff_Filter_Option
      i = merge(3,1,local_Filter_Option==2)
      Filter_Cutoff = i*cutoff_order
      fetac_qp = real(Filter_Cutoff,kind=qp) / real(i*order,kind=qp)
      feps_qp = qzero
    else if (cutoff_order == -2) then
      fetac_qp = real( Sponge_Filter_etac , kind=qp )
      Filter_Cutoff = int_cutoff_function( &
                        max(zero,min(one,Sponge_Filter_etac)) * &
                        real(order,kind=wp) )
      local_Filter_Option = 1
      feps_qp = qone
      falpha_qp = real( Sponge_Filter_Alpha , kind=qp )
      filt_order = Sponge_Filter_Order
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = int_cutoff_function( &
                        max(zero,min(one,Limiter_Filter_etac)) * &
                        real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp  = real( Filter_etac , kind=qp )
    Filter_Cutoff = int_cutoff_function( &
                      max(zero,min(one,Filter_etac)) * &
                      real(order,kind=wp) )
    local_Filter_Option = Filter_Option
    feps_qp = real( one - Filter_eps , kind=qp )
  end if
  !
  ! Allocate the return array containing the working
  ! precision filter matrix
  !
  allocate ( return_value(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  ! Create the filter matrix using quad precision
  !
  allocate ( filt_diag_qp(1:nsp) , source=qzero , stat=ierr )
  !
  ! Initialize the local filter matrices to identity matrices
  !
  do i = 1,nsp
    filt_diag_qp(i) = qone
  end do
  !
  ! Compute the filter diagonal for hexahedrals
  !
#ifdef TIME_FILTER_MATMUL
  call debug_timer(start_timer,"creating diagonal filter matrix")
#endif
  if (local_Filter_Option == 1) then
    !
    ! Option 1 - Tensor product of the 1D filter
    !           (Equivalent to using a 1D filter in each direction)
    !
    ! Allocate the 1D filter vector and set all its values to one
    !
    allocate ( filter_1D_qp(1:np) , source=qone , stat=ierr )
    !
    ! First get the 1D filter diagonal
    !
    do i = 0,order
      if (i >= Filter_Cutoff) then
        !
        rn_qp = real(i,kind=qp)/real(order,kind=qp) - fetac_qp
        rd_qp = max( qzero , qone - fetac_qp )
        !
        if (rd_qp == qzero) then
          filter_1D_qp(i+1) = feps_qp
        else
          filter_1D_qp(i+1) = feps_qp * &
                              exp(-falpha_qp * (rn_qp/rd_qp)**filt_order)
        end if
        !
      end if
    end do
    !
    ! Now fill the filter diagonal for hexes
    !
    do k = 0,order
      do j = 0,order
        do i = 0,order
          !
          l = np*np*k + np*j + i + 1
          !
          filt_diag_qp(l) = filter_1D_qp(i+1) * &
                            filter_1D_qp(j+1) * &
                            filter_1D_qp(k+1)
          !
        end do
      end do
    end do
    !
  else if (local_Filter_Option == 2) then
    !
    ! Option 2 - The hex basis functions actually contains cross polynomial
    !            terms up to order 3*p (as opposed to just 1*p for simplexes).
    !            This filter applies the cut-off frequency based on 3*p as
    !            opposed to basing it on just p in 1D like Option 1.
    !
    do k = 0,order
      do j = 0,order
        do i = 0,order
          !
          if (i+j+k >= Filter_Cutoff) then
            !
            l = np*np*k + np*j + i + 1
            !
            rn_qp = real(i+j+k-Filter_Cutoff,kind=qp)
            rd_qp = max( qzero , real(3*order-Filter_Cutoff,kind=qp) )
            !
            if (rd_qp == qzero) then
              filt_diag_qp(l) = feps_qp
            else
              filt_diag_qp(l) = feps_qp * &
                                exp(-falpha_qp*(rn_qp/rd_qp)**filt_order)
            end if
            !
          end if
          !
        end do
      end do
    end do
    !
  end if
#ifdef TIME_FILTER_MATMUL
  call debug_timer(stop_timer,"creating diagonal filter matrix")
#endif
  !
  ! Compute the final filter matrix for hexahedrals
  !
#ifdef SPECIAL_FOR_INTEL
  !
  allocate ( fmat(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  filt_diag = chop( filt_diag_qp(:) ) ! F2003 AUTO-REALLOCATION
  !
#ifdef TIME_FILTER_MATMUL
#ifdef DEBUG_ON
  num_zeros = count(abs(vand) < eps14)
  write (iout,100) "vand",num_zeros,size(vand), &
                   100.0_wp*real(num_zeros,kind=wp)/real(size(vand),kind=wp)
  num_zeros = count(abs(vand_inv) < eps14)
  write (iout,100) "vand_inv",num_zeros,size(vand_inv), &
                   100.0_wp*real(num_zeros,kind=wp)/real(size(vand_inv),kind=wp)
  100 format ("Sparsity of ",a8," is ",i9,"/",i9," (",f10.6," % )")
#endif
  call debug_timer(start_timer,"V x F")
#endif
  !
  do j = 1,nsp
    do i = 1,nsp
      fmat(i,j) = vand(i,j)*filt_diag(j)
    end do
  end do
  !
#ifdef TIME_FILTER_MATMUL
  call debug_timer(stop_timer,"V x F")
  !
  call debug_timer(start_timer,"V x F x V_inv")
#endif
  !
  call gemm(fmat,vand_inv,return_value)
  !
#ifdef TIME_FILTER_MATMUL
  call debug_timer(stop_timer,"V x F x V_inv")
#endif
  !
#else
  !
  fmat_qp = real(vand,kind=qp) * & ! F2003 AUTO-REALLOCATION
            spread(filt_diag_qp,dim=1,ncopies=size(vand,dim=1))
 !fmat_qp = matmul( real(vand,kind=qp) , fmat_qp )
  fmat_qp = matmul( fmat_qp , real(vand_inv,kind=qp) )
  !
  ! Use the chop function to copy fmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( fmat_qp )
  !
#endif
  !
end function FilterMatrix_Hexa
!
!###############################################################################
!
function FilterMatrix_Hexa_with_output(vand,vand_inv,cutoff_order) &
                                result(return_value)
  !
#ifdef SPECIAL_FOR_INTEL
  !.. Use Statements ..
  use blas95, only : gemm
#endif
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: cutoff_order
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,np,nsp,order
  integer  :: Filter_Cutoff,ierr
  integer  :: local_Filter_Option
  integer  :: filt_order
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
  real(qp) :: feps_qp
  real(wp) :: a
  integer  :: num_zeros
  !
  !.. Local Arrays ..
  real(qp), allocatable :: filter_1D_qp(:)
  real(qp), allocatable :: filt_diag_qp(:)
  real(qp), allocatable :: fmat_qp(:,:)
  real(wp), allocatable :: filt_diag(:)
  real(wp), allocatable :: fmat(:,:)
  !
continue
  !
  ! Get the number of solution points
  !
  nsp = size(vand,dim=1)
  order = solpts2order_hexa(nsp)
  np = order+1
 !np = nint( real(nsp,kind=wp)**one3 )
 !order = np-1
  !
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  !
  filt_order = Filter_Order
  !
  if (present(cutoff_order)) then
    if (cutoff_order >= 0) then
      local_Filter_Option = Cutoff_Filter_Option
      i = merge(3,1,local_Filter_Option==2)
      Filter_Cutoff = i*cutoff_order
      fetac_qp = real(Filter_Cutoff,kind=qp) / real(i*order,kind=qp)
      feps_qp = qzero
    else if (cutoff_order == -2) then
      fetac_qp = real( Sponge_Filter_etac , kind=qp )
      Filter_Cutoff = int_cutoff_function( &
                        max(zero,min(one,Sponge_Filter_etac)) * &
                        real(order,kind=wp) )
      local_Filter_Option = 1
      feps_qp = qone
      falpha_qp = real( Sponge_Filter_Alpha , kind=qp )
      filt_order = Sponge_Filter_Order
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = int_cutoff_function( &
                        max(zero,min(one,Limiter_Filter_etac)) * &
                        real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp  = real( Filter_etac , kind=qp )
    Filter_Cutoff = int_cutoff_function( &
                      max(zero,min(one,Filter_etac)) * &
                      real(order,kind=wp) )
    local_Filter_Option = Filter_Option
    feps_qp = real( one - Filter_eps , kind=qp )
  end if
  !
  ! Allocate the return array containing the working
  ! precision filter matrix
  !
  allocate ( return_value(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  ! Create the filter matrix using quad precision
  !
  allocate ( filt_diag_qp(1:nsp) , source=qzero , stat=ierr )
  !
  ! Initialize the local filter matrices to identity matrices
  !
  do i = 1,nsp
    filt_diag_qp(i) = qone
  end do
  !
  ! Compute the filter diagonal for hexahedrals
  !
  if (local_Filter_Option == 1) then
    !
    ! Option 1 - Tensor product of the 1D filter
    !           (Equivalent to using a 1D filter in each direction)
    !
    ! Allocate the 1D filter vector and set all its values to one
    !
    allocate ( filter_1D_qp(1:np) , source=qone , stat=ierr )
    !
    ! First get the 1D filter diagonal
    !
    do i = 0,order
      if (i >= Filter_Cutoff) then
        !
        rn_qp = real(i,kind=qp)/real(order,kind=qp) - fetac_qp
        rd_qp = max( qzero , qone - fetac_qp )
        !
        if (rd_qp == qzero) then
          filter_1D_qp(i+1) = feps_qp
        else
          filter_1D_qp(i+1) = feps_qp * &
                              exp(-falpha_qp * (rn_qp/rd_qp)**filt_order)
        end if
        !
      end if
    end do
    !
    ! Now fill the filter diagonal for hexes
    !
    do k = 0,order
      do j = 0,order
        do i = 0,order
          !
          l = np*np*k + np*j + i + 1
          !
          filt_diag_qp(l) = filter_1D_qp(i+1) * &
                            filter_1D_qp(j+1) * &
                            filter_1D_qp(k+1)
          !
        end do
      end do
    end do
    !
  else if (local_Filter_Option == 2) then
    !
    ! Option 2 - The hex basis functions actually contains cross polynomial
    !            terms up to order 3*p (as opposed to just 1*p for simplexes).
    !            This filter applies the cut-off frequency based on 3*p as
    !            opposed to basing it on just p in 1D like Option 1.
    !
    do k = 0,order
      do j = 0,order
        do i = 0,order
          !
          if (i+j+k >= Filter_Cutoff) then
            !
            l = np*np*k + np*j + i + 1
            !
            rn_qp = real(i+j+k-Filter_Cutoff,kind=qp)
            rd_qp = max( qzero , real(3*order-Filter_Cutoff,kind=qp) )
            !
            if (rd_qp == qzero) then
              filt_diag_qp(l) = feps_qp
            else
              filt_diag_qp(l) = feps_qp * &
                                exp(-falpha_qp*(rn_qp/rd_qp)**filt_order)
            end if
            !
          end if
          !
        end do
      end do
    end do
    !
  end if
  !
  if (local_Filter_Option /= 0) then
    !
    if (mypnum == glb_root) then
      if (present(cutoff_order)) then
        if (cutoff_order == -1) then
          write (io,102)
          102 format (/,"Diagonal for hexa limiter filter")
        else
          write (io,101) order+1-(Filter_cutoff)
          101 format (/,"Diagonal for hexa cutoff filter for P-",i0)
        end if
      else
        write (io,103)
        103 format (/,"Diagonal for hexa standard filter")
      end if
      do k = 0,order
        do j = 0,order
          do i = 0,order
            l = np*np+k + np*j + i + 1
            write (io,100) i,j,k,filt_diag_qp(l)
          end do
        end do
      end do
      100 format (2x,i2,1x,i2,1x,i2,4x,es14.6)
    end if
    !
  end if
  !
  ! Compute the final filter matrix for hexahedrals
  !
#ifdef SPECIAL_FOR_INTEL
  !
  allocate ( fmat(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  filt_diag = chop( filt_diag_qp(:) ) ! F2003 AUTO-REALLOCATION
  !
  do j = 1,nsp
    do i = 1,nsp
      fmat(i,j) = vand(i,j)*filt_diag(j)
    end do
  end do
  !
  call gemm(fmat,vand_inv,return_value)
  !
#else
  !
  fmat_qp = real(vand,kind=qp) * & ! F2003 AUTO-REALLOCATION
            spread(filt_diag_qp,dim=1,ncopies=size(vand,dim=1))
 !fmat_qp = matmul( real(vand,kind=qp) , fmat_qp )
  fmat_qp = matmul( fmat_qp , real(vand_inv,kind=qp) )
  !
  ! Use the chop function to copy fmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( fmat_qp )
  !
#endif
  !
end function FilterMatrix_Hexa_with_output
!
!###############################################################################
!
pure function int_cutoff_function(input_value) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: input_value
  !
  !.. Function Result ..
  integer :: return_value
  !
continue
  !
  ! Old integer conversion function
  ! This can cause an incorrect value for Filter_Cutoff.
  ! For example, take two very similar P3 filters:
  !   1) alpha = 36 ; etac = 0.75 ; order = 8
  !   2) alpha = 36 ; etac = 0.8  ; order = 8
  ! Analytically, the second filter is very slightly weaker than the first.
  ! Numerically, they should be the same when evaluated at finite integer
  ! orders (i=0,3). The filter cutoffs should be 2.25 and 2.4, so all integer
  ! orders less than these values should give a filter diagonal value of 1.
  ! However, using nint means the filter cutoffs will be 2, so the diagonal
  ! for order i=2 will not be 1. In fact the i=2 diagonal value for the first
  ! filter will be 0.9945281 whereas the value will be 0.2454495 for the
  ! second, which makes the second filter significantly weaker than the first
  ! when using the nint function! Both filters should have i=2 diagonal values
  ! of 1.
  !
 !return_value = nint( input_value )
  !
  ! Correct integer conversion function
  !
  return_value = ceiling( input_value )
  !
end function int_cutoff_function
!
!###############################################################################
!
subroutine output_all_filters()
  !
  !.. Use Statements ..
  use ovar, only : use_sponge_filter
  use ovar, only : Filter_Option
  use ovar, only : Limiter_Filter_Option
  !
  !.. Local Scalars ..
  integer :: geom,p,ierr,co,i,j
  integer :: io,iog,iol,ioc,ios
  character(len=200) :: array_name
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: filt_diag(:)
  real(wp), allocatable :: fmat(:,:)
  integer,  allocatable :: idx(:,:)
  !
  character(len=:), allocatable :: gfname,lfname,cfname,sfname
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "output_all_filters"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (mypnum == glb_root) then
    !
    ! Open a file for each filter type
    !
    gfname = filt_fname("global")  ! F2003 AUTO-REALLOC
    lfname = filt_fname("limiter") ! F2003 AUTO-REALLOC
    cfname = filt_fname("cutoff")  ! F2003 AUTO-REALLOC
    sfname = filt_fname("sponge")  ! F2003 AUTO-REALLOC
    !
    open (newunit=iog,file=gfname,status="replace",action="write", &
                      iostat=ierr,iomsg=error_message)
    call io_error(pname,gfname,1,__LINE__,__FILE__,ierr,error_message)
    !
    open (newunit=iol,file=lfname,status="replace",action="write", &
                      iostat=ierr,iomsg=error_message)
    call io_error(pname,lfname,1,__LINE__,__FILE__,ierr,error_message)
    !
    open (newunit=ioc,file=cfname,status="replace",action="write", &
                      iostat=ierr,iomsg=error_message)
    call io_error(pname,cfname,1,__LINE__,__FILE__,ierr,error_message)
    !
    open (newunit=ios,file=sfname,status="replace",action="write", &
                      iostat=ierr,iomsg=error_message)
    call io_error(pname,sfname,1,__LINE__,__FILE__,ierr,error_message)
    !
    !
    !
    geom_loop: &
    do geom = lbound(filter,dim=1),ubound(filter,dim=1)
      !
      write (iog,100) trim(Geom_Full_Name(geom))
      write (iol,100) trim(Geom_Full_Name(geom))
      write (ioc,100) trim(Geom_Full_Name(geom))
      write (ios,100) trim(Geom_Full_Name(geom))
      !
      order_loop: &
      do p = lbound(filter,dim=2),ubound(filter,dim=2)
        !
        ! Get the indexing for each diagonal term. This is the same
        ! regardless of the filter type
        !
        idx = get_geom_indices(geom,p) ! F2003 AUTO-REALLOC
        !
        !#################################
        !#################################
        !###  Global filter
        !#################################
        !#################################
        !
        io = iog
        !
        if (.not. allocated(filter(geom,p)%exponential)) then
          write (array_name,1) Geom_Name(geom),p,"exponential"
          write (io,900) trim(adjustl(array_name))
          goto 1001
        end if
        !
        if (.not. allocated(filter(geom,p)%exponential%mat)) then
          write (array_name,1) Geom_Name(geom),p,"exponential"
          write (io,901) trim(adjustl(array_name))
          goto 1001
        end if
        !
        write (io,101) "Global",p,trim(Geom_Full_Name(geom))
        !
        fmat = filter(geom,p)%exponential%mat(:,:) ! F2003 AUTO-REALLOC
        !
        filt_diag = extract_filter_diag(geom,p,fmat) ! F2003 AUTO-REALLOC
        !
        ! Write out the filter diagonal and indices
        !
        write (io,200)
        do i = 1,size(filt_diag)
          write (io,201) filt_diag(i), (idx(j,i), j=1,size(idx,dim=1))
        end do
        !
        ! Write out the actual filter matrix
        !
        write (io,300)
        if (Filter_Option == 0) write (io,302)
        do j = size(fmat,dim=2),1,-1
          write (io,301) (fmat(i,j), i=1,size(fmat,dim=1))
        end do
        !
        1001 continue
        !
        !#################################
        !#################################
        !###  Limiter filter
        !#################################
        !#################################
        !
        io = iol
        !
        if (.not. allocated(filter(geom,p)%limiter)) then
          write (array_name,1) Geom_Name(geom),p,"limiter"
          write (io,900) trim(adjustl(array_name))
          goto 1002
        end if
        !
        if (.not. allocated(filter(geom,p)%limiter%mat)) then
          write (array_name,1) Geom_Name(geom),p,"limiter"
          write (io,901) trim(adjustl(array_name))
          goto 1002
        end if
        !
        write (io,101) "Limiter",p,trim(Geom_Full_Name(geom))
        !
        fmat = filter(geom,p)%limiter%mat(:,:) ! F2003 AUTO-REALLOC
        !
        filt_diag = extract_filter_diag(geom,p,fmat) ! F2003 AUTO-REALLOC
        !
        ! Write out the filter diagonal and indices
        !
        write (io,200)
        do i = 1,size(filt_diag)
          write (io,201) filt_diag(i), (idx(j,i), j=1,size(idx,dim=1))
        end do
        !
        ! Write out the actual filter matrix
        !
        write (io,300)
        if (Limiter_Filter_Option == 0) write (io,302)
        do j = size(fmat,dim=2),1,-1
          write (io,301) (fmat(i,j), i=1,size(fmat,dim=1))
        end do
        !
        1002 continue
        !
        !#################################
        !#################################
        !###  Sponge filter
        !#################################
        !#################################
        !
        io = ios
        !
        if (.not. allocated(filter(geom,p)%sponge)) then
          write (array_name,1) Geom_Name(geom),p,"sponge"
          write (io,900) trim(adjustl(array_name))
          goto 1003
        end if
        !
        if (.not. allocated(filter(geom,p)%sponge%mat)) then
          write (array_name,1) Geom_Name(geom),p,"sponge"
          write (io,901) trim(adjustl(array_name))
          goto 1003
        end if
        !
        write (io,101) "Sponge",p,trim(Geom_Full_Name(geom))
        !
        fmat = filter(geom,p)%sponge%mat(:,:) ! F2003 AUTO-REALLOC
        !
        filt_diag = extract_filter_diag(geom,p,fmat) ! F2003 AUTO-REALLOC
        !
        ! Write out the filter diagonal and indices
        !
        write (io,200)
        do i = 1,size(filt_diag)
          write (io,201) filt_diag(i), (idx(j,i), j=1,size(idx,dim=1))
        end do
        !
        ! Write out the actual filter matrix
        !
        write (io,300)
        if (.not. use_sponge_filter) write (io,302)
        do j = size(fmat,dim=2),1,-1
          write (io,301) (fmat(i,j), i=1,size(fmat,dim=1))
        end do
        !
        1003 continue
        !
        !#################################
        !#################################
        !###  Cutoff filters
        !#################################
        !#################################
        !
        io = ioc
        !
        if (.not. allocated(filter(geom,p)%cutoff)) then
          write (array_name,1) Geom_Name(geom),p,"cutoff"
          write (io,900) trim(adjustl(array_name))
          goto 1004
        end if
        !
        co_loop: &
        do co = ubound(filter(geom,p)%cutoff,dim=1), &
                lbound(filter(geom,p)%cutoff,dim=1),-1
          !
          if (.not. allocated(filter(geom,p)%cutoff(co)%mat)) then
            write (array_name,1) Geom_Name(geom),p,"cutoff",co
            write (io,901) trim(adjustl(array_name))
            cycle co_loop
          end if
          !
          write (io,102) p+1-co,"Cutoff",p,trim(Geom_Full_Name(geom))
          !
          fmat = filter(geom,p)%cutoff(co)%mat(:,:) ! F2003 AUTO-REALLOC
          !
          filt_diag = extract_filter_diag(geom,p,fmat) ! F2003 AUTO-REALLOC
          !
          ! Write out the filter diagonal and indices
          !
          write (io,200)
          do i = 1,size(filt_diag)
            write (io,201) filt_diag(i), (idx(j,i), j=1,size(idx,dim=1))
          end do
          !
          ! Write out the actual filter matrix
          !
          write (io,300)
          do j = size(fmat,dim=2),1,-1
            write (io,301) (fmat(i,j), i=1,size(fmat,dim=1))
          end do
          !
        end do co_loop
        !
        1004 continue
        !
      end do order_loop
      !
    end do geom_loop
    !
    ! Close all of the files
    !
    flush (iog)
    close (iog,iostat=ierr,iomsg=error_message)
    call io_error(pname,gfname,2,__LINE__,__FILE__,ierr,error_message)
    !
    flush (iol)
    close (iol,iostat=ierr,iomsg=error_message)
    call io_error(pname,lfname,2,__LINE__,__FILE__,ierr,error_message)
    !
    flush (ioc)
    close (ioc,iostat=ierr,iomsg=error_message)
    call io_error(pname,cfname,2,__LINE__,__FILE__,ierr,error_message)
    !
    flush (ios)
    close (ios,iostat=ierr,iomsg=error_message)
    call io_error(pname,sfname,2,__LINE__,__FILE__,ierr,error_message)
    !
  end if
  !
  write (error_message,903)
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1   format ("filter(",a,",",i0,")%",a,:,"(",i0,")")
  100 format (//,3(80("#"),/), &
              4x,"Filter matrices for ",a," elements",/, &
              3(80("#"),/),/)
  101 format (/,80("="),/, &
              4x,a,1x,"filter for P",i0,1x,a," elements",/, &
              80("="),/)
  102 format (/,80("="),/, &
              4x,"P-",i0,1x,a,1x,"filter for P",i0,1x,a," elements",/, &
              80("="),/)
  200 format (2x,"Filter diagonal coefficients with modal indices",/)
  201 format (4x,es18.10,6x,i2,:,", ",i2,:,", ",i2)
  300 format (/,80("-"),//,2x,"Actual filter matrix operator",/)
  301 format (125(2x,es14.6))
  302 format ("  NOTE: CURRENTLY NOT USED BECAUSE ", &
              "IT IS TURNED OFF IN THE INPUT FILE!")
  !
  900 format ("The filter matrix ---  ",a,"  --- is not allocated")
  901 format ("The filter matrix ---  ",a,"  --- IS allocated,", &
              " the %mat array is NOT allocated!")
  903 format ("AFTER WRITING OUT THE FILTER DIAGONALS ", &
              "TO 'filter_diagonals.dat'.")
  !
end subroutine output_all_filters
!
!###############################################################################
!
function extract_filter_diag(geom,p,fmat) result(return_value)
  !
  ! Function to extract the filter diagonal from the filter matrix
  !
  !.. Use Statements ..
  use vandermonde_mod, only : vand
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: geom
  integer,  intent(in) :: p
  real(wp), intent(in) :: fmat(:,:)
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: i,nsp
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: filt(:,:)
  !
continue
  !
  nsp = size(fmat,dim=1)
  !
  return_value = init_array(nsp,zero) ! F2003 AUTO-REALLOC
  !
  filt = matmul(vand(geom,p)%nodal2modal%mat,fmat) ! F2003 AUTO-REALLOC
  filt(:,:) = matmul(filt,vand(geom,p)%modal2nodal%mat)
  !
  do i = 1,nsp
    return_value(i) = filt(i,i)
  end do
  !
end function extract_filter_diag
!
!###############################################################################
!
function filt_fname(filt_type) result(return_value)
  !
  !.. Use Statements ..
  use order_mod, only : geom_solpts
  use ovar,      only : Filter_Alpha,Limiter_Filter_Alpha,Sponge_Filter_Alpha
  use ovar,      only : Filter_etac, Limiter_Filter_etac, Sponge_Filter_etac
  use ovar,      only : Filter_order,Limiter_Filter_order,Sponge_Filter_order
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: filt_type
  !
  !.. Function Result ..
  character(len=:), allocatable :: return_value
  !
  !.. Local Scalars ..
  character(len=100) :: fname
  character(len=:), allocatable :: ftype
  !
continue
  !
  ftype = trim(adjustl(filt_type)) ! F2003 AUTO-REALLOC
  !
  if (uppercase(ftype) == "GLOBAL") then
    !
    write (fname,1) ftype,nint(Filter_Alpha),Filter_Order,Filter_etac
    !
  else if (uppercase(ftype) == "LIMITER") then
    !
    write (fname,1) ftype,nint(Limiter_Filter_Alpha), &
                    Limiter_Filter_Order,Limiter_Filter_etac
    !
  else if (uppercase(ftype) == "SPONGE") then
    !
    write (fname,1) ftype,nint(Sponge_Filter_Alpha), &
                    Sponge_Filter_Order,Sponge_Filter_etac
    !
  else if (uppercase(ftype) == "CUTOFF") then
    !
    fname = "filter_diagonal.cutoff.dat"
    !
  end if
  !
  return_value = trim(adjustl(fname)) ! F2003 AUTO-REALLOC
  !
  1 format ("filter_diagonal.",a,".a-",i0,"_s-",i0,"_eta-",f4.2,".dat")
  !
end function filt_fname
!
!###############################################################################
!
function get_geom_indices(geom,p) result(return_value)
  !
  !.. Use Statements ..
  use order_mod, only : geom_solpts
  !
  !.. Formal Arguments ..
  integer, intent(in) :: geom
  integer, intent(in) :: p
  !
  !.. Function Result ..
  integer, allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,nd,nsp
  !
continue
  !
  nsp = geom_solpts(geom,p)
  !
  nd = geom_dimen(geom)
  !
  return_value = init_array(nd,nsp,0) ! F2003 AUTO-REALLOC
  !
  if (geom == Geom_Edge) then
    !
    l = 0
    do i = 0,p
      l = l + 1
      return_value(1:nd,l) = [i]
    end do
    !
  else if (geom == Geom_Quad) then
    !
    l = 0
    do j = 0,p
      do i = 0,p
        l = l + 1
        return_value(1:nd,l) = [i,j]
      end do
    end do
    !
  else if (geom == Geom_Tria) then
    !
    l = 0
    do j = 0,p
      do i = 0,p-j
        l = l + 1
        return_value(1:nd,l) = [i,j]
      end do
    end do
    !
  else if (geom == Geom_Tetr) then
    !
    l = 0
    do k = 0,p
      do j = 0,p-k
        do i = 0,p-k-j
          l = l + 1
          return_value(1:nd,l) = [i,j,k]
        end do
      end do
    end do
    !
  else if (geom == Geom_Pyra) then
    !
    l = 0
    do k = 0,p
      do j = 0,p-k
        do i = 0,p-k
          l = l + 1
          return_value(1:nd,l) = [i,j,k]
        end do
      end do
    end do
    !
  else if (geom == Geom_Pris) then
    !
    l = 0
    do k = 0,p
      do j = 0,p
        do i = 0,p-j
          l = l + 1
          return_value(1:nd,l) = [i,j,k]
        end do
      end do
    end do
    !
  else if (geom == Geom_Hexa) then
    !
    l = 0
    do k = 0,p
      do j = 0,p
        do i = 0,p
          l = l + 1
          return_value(1:nd,l) = [i,j,k]
        end do
      end do
    end do
    !
  end if
  !
end function get_geom_indices
!
!###############################################################################
!
subroutine filter_memory_usage(iunit)
  !
  !.. Use Statements ..
  use module_memory
#ifndef DISABLE_DTIO
  use matrix_types_mod, only : matrix_memory_usage
#endif
  !
  !.. Formal Arguments ..
  integer, optional, intent(in) :: iunit
  !
  !.. Local Scalars ..
  integer :: iou,i,j,k
  character(len=36) :: array_name
  type(memory) :: array,dt
  type(memory) :: total
  !
#ifndef DISABLE_DTIO
continue
  !
  iou = iout
  if (present(iunit)) iou = iunit
  !
  call total%reset
  !
  write (iou,1)
  !
  ! Check the memory of the arrays within this module
  !
  if (allocated(filter)) then
    !
    call dt%set( storage_size(filter,kind=inttype) , &
                         size(filter,kind=inttype) )
    !
    do j = lbound(filter,dim=2),ubound(filter,dim=2)
      do i = lbound(filter,dim=1),ubound(filter,dim=1)
        !
        if (allocated(filter(i,j)%exponential)) then
          !
          call array%set( storage_size(filter(i,j)%exponential,kind=inttype) )
          !
          write (array_name,4) Geom_Name(i),j,"exponential"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          array = matrix_memory_usage( filter(i,j)%exponential )
          !
          write (array_name,4) Geom_Name(i),j,"exponential%mat"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
        end if
        !
        if (allocated(filter(i,j)%cutoff)) then
          !
          call array%set( storage_size(filter(i,j)%cutoff,kind=inttype) , &
                                  size(filter(i,j)%cutoff,kind=inttype) )
          !
          write (array_name,4) Geom_Name(i),j,"cutoff"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          do k = lbound(filter(i,j)%cutoff,dim=1), &
                 ubound(filter(i,j)%cutoff,dim=1)
            !
            array = matrix_memory_usage( filter(i,j)%cutoff(k) )
            !
            write (array_name,4) Geom_Name(i),j,"cutoff",k
            write (iou,2) array_name,array%get_units()
            !
            call dt%add(array)
            !
          end do
          !
        end if
        !
      end do
    end do
    !
    call total%add(dt)
    write (array_name,5) "filter"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "filter"
    !
  end if
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module filter_mod")
  2 format (" Array: ",a,dt)
  3 format (" Total memory = ",f8.2,1x,a)
  4 format ("filter(",a,",",i0,")",a,:,"(",i0,")%mat")
  5 format ("'",a,"'")
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  !
#endif
end subroutine filter_memory_usage
!
!###############################################################################
!
subroutine move_filter_to_device()
  !
  !.. Local Scalars ..
  integer :: i,i1,i2,j,j1,j2,k,k1,k2,m1,m2,n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "move_filter_to_device"
  !
continue
  !
  if (allocated(filter)) then
    !$ACC ENTER DATA COPYIN(filter)
    i1 = lbound(filter,dim=1)
    i2 = ubound(filter,dim=1)
    j1 = lbound(filter,dim=2)
    j2 = ubound(filter,dim=2)
    !$ACC ENTER DATA COPYIN(filter(i1:i2,j1:j2))
    do j = j1,j2
      do i = i1,i2
        !
        !---------------------------------------------------------------
        !
        if (allocated(filter(i,j)%exponential)) then
          !$ACC ENTER DATA COPYIN(filter(i,j)%exponential)
          if (allocated(filter(i,j)%exponential%mat)) then
            !$ACC ENTER DATA COPYIN(filter(i,j)%exponential%mat)
            m1 = lbound(filter(i,j)%exponential%mat,dim=1)
            m2 = ubound(filter(i,j)%exponential%mat,dim=1)
            n1 = lbound(filter(i,j)%exponential%mat,dim=2)
            n2 = ubound(filter(i,j)%exponential%mat,dim=2)
            !$ACC ENTER DATA COPYIN(filter(i,j)%exponential%mat(m1:m2,n1:n2))
          end if
          !
          if (allocated(filter(i,j)%exponential%row_ptr)) then
            !$ACC ENTER DATA COPYIN(filter(i,j)%exponential%row_ptr)
            m1 = lbound(filter(i,j)%exponential%row_ptr,dim=1)
            m2 = ubound(filter(i,j)%exponential%row_ptr,dim=1)
            !$ACC ENTER DATA COPYIN(filter(i,j)%exponential%row_ptr(m1:m2))
          end if
          !
          if (allocated(filter(i,j)%exponential%row_idx)) then
            !$ACC ENTER DATA COPYIN(filter(i,j)%exponential%row_idx)
            m1 = lbound(filter(i,j)%exponential%row_idx,dim=1)
            m2 = ubound(filter(i,j)%exponential%row_idx,dim=1)
            !$ACC ENTER DATA COPYIN(filter(i,j)%exponential%row_idx(m1:m2))
          end if
          !
          if (allocated(filter(i,j)%exponential%row_val)) then
            !$ACC ENTER DATA COPYIN(filter(i,j)%exponential%row_val)
            m1 = lbound(filter(i,j)%exponential%row_val,dim=1)
            m2 = ubound(filter(i,j)%exponential%row_val,dim=1)
            !$ACC ENTER DATA COPYIN(filter(i,j)%exponential%row_val(m1:m2))
          end if
          !
          if (allocated(filter(i,j)%exponential%col_ptr)) then
            !$ACC ENTER DATA COPYIN(filter(i,j)%exponential%col_ptr)
            m1 = lbound(filter(i,j)%exponential%col_ptr,dim=1)
            m2 = ubound(filter(i,j)%exponential%col_ptr,dim=1)
            !$ACC ENTER DATA COPYIN(filter(i,j)%exponential%col_ptr(m1:m2))
          end if
          !
          if (allocated(filter(i,j)%exponential%col_idx)) then
            !$ACC ENTER DATA COPYIN(filter(i,j)%exponential%col_idx)
            m1 = lbound(filter(i,j)%exponential%col_idx,dim=1)
            m2 = ubound(filter(i,j)%exponential%col_idx,dim=1)
            !$ACC ENTER DATA COPYIN(filter(i,j)%exponential%col_idx(m1:m2))
          end if
          !
          if (allocated(filter(i,j)%exponential%col_val)) then
            !$ACC ENTER DATA COPYIN(filter(i,j)%exponential%col_val)
            m1 = lbound(filter(i,j)%exponential%col_val,dim=1)
            m2 = ubound(filter(i,j)%exponential%col_val,dim=1)
            !$ACC ENTER DATA COPYIN(filter(i,j)%exponential%col_val(m1:m2))
          end if
        end if
        !
        !---------------------------------------------------------------
        !
        if (allocated(filter(i,j)%limiter)) then
          !$ACC ENTER DATA COPYIN(filter(i,j)%limiter)
          if (allocated(filter(i,j)%limiter%mat)) then
            !$ACC ENTER DATA COPYIN(filter(i,j)%limiter%mat)
            m1 = lbound(filter(i,j)%limiter%mat,dim=1)
            m2 = ubound(filter(i,j)%limiter%mat,dim=1)
            n1 = lbound(filter(i,j)%limiter%mat,dim=2)
            n2 = ubound(filter(i,j)%limiter%mat,dim=2)
            !$ACC ENTER DATA COPYIN(filter(i,j)%limiter%mat(m1:m2,n1:n2))
          end if
          !
          if (allocated(filter(i,j)%limiter%row_ptr)) then
            !$ACC ENTER DATA COPYIN(filter(i,j)%limiter%row_ptr)
            m1 = lbound(filter(i,j)%limiter%row_ptr,dim=1)
            m2 = ubound(filter(i,j)%limiter%row_ptr,dim=1)
            !$ACC ENTER DATA COPYIN(filter(i,j)%limiter%row_ptr(m1:m2))
          end if
          !
          if (allocated(filter(i,j)%limiter%row_idx)) then
            !$ACC ENTER DATA COPYIN(filter(i,j)%limiter%row_idx)
            m1 = lbound(filter(i,j)%limiter%row_idx,dim=1)
            m2 = ubound(filter(i,j)%limiter%row_idx,dim=1)
            !$ACC ENTER DATA COPYIN(filter(i,j)%limiter%row_idx(m1:m2))
          end if
          !
          if (allocated(filter(i,j)%limiter%row_val)) then
            !$ACC ENTER DATA COPYIN(filter(i,j)%limiter%row_val)
            m1 = lbound(filter(i,j)%limiter%row_val,dim=1)
            m2 = ubound(filter(i,j)%limiter%row_val,dim=1)
            !$ACC ENTER DATA COPYIN(filter(i,j)%limiter%row_val(m1:m2))
          end if
          !
          if (allocated(filter(i,j)%limiter%col_ptr)) then
            !$ACC ENTER DATA COPYIN(filter(i,j)%limiter%col_ptr)
            m1 = lbound(filter(i,j)%limiter%col_ptr,dim=1)
            m2 = ubound(filter(i,j)%limiter%col_ptr,dim=1)
            !$ACC ENTER DATA COPYIN(filter(i,j)%limiter%col_ptr(m1:m2))
          end if
          !
          if (allocated(filter(i,j)%limiter%col_idx)) then
            !$ACC ENTER DATA COPYIN(filter(i,j)%limiter%col_idx)
            m1 = lbound(filter(i,j)%limiter%col_idx,dim=1)
            m2 = ubound(filter(i,j)%limiter%col_idx,dim=1)
            !$ACC ENTER DATA COPYIN(filter(i,j)%limiter%col_idx(m1:m2))
          end if
          !
          if (allocated(filter(i,j)%limiter%col_val)) then
            !$ACC ENTER DATA COPYIN(filter(i,j)%limiter%col_val)
            m1 = lbound(filter(i,j)%limiter%col_val,dim=1)
            m2 = ubound(filter(i,j)%limiter%col_val,dim=1)
            !$ACC ENTER DATA COPYIN(filter(i,j)%limiter%col_val(m1:m2))
          end if
        end if
        !
        !---------------------------------------------------------------
        !
        if (allocated(filter(i,j)%cutoff)) then
          !$ACC ENTER DATA COPYIN(filter(i,j)%cutoff)
          k1 = lbound(filter(i,j)%cutoff,dim=1)
          k2 = ubound(filter(i,j)%cutoff,dim=1)
          !$ACC ENTER DATA COPYIN(filter(i,j)%cutoff(k1:k2))
          do k = k1,k2
            !
            if (allocated(filter(i,j)%cutoff(k)%mat)) then
              !$ACC ENTER DATA COPYIN(filter(i,j)%cutoff(k)%mat)
              m1 = lbound(filter(i,j)%cutoff(k)%mat,dim=1)
              m2 = ubound(filter(i,j)%cutoff(k)%mat,dim=1)
              n1 = lbound(filter(i,j)%cutoff(k)%mat,dim=2)
              n2 = ubound(filter(i,j)%cutoff(k)%mat,dim=2)
              !$ACC ENTER DATA COPYIN(filter(i,j)%cutoff(k)%mat(m1:m2,n1:n2))
            end if
            !
            if (allocated(filter(i,j)%cutoff(k)%row_ptr)) then
              !$ACC ENTER DATA COPYIN(filter(i,j)%cutoff(k)%row_ptr)
              m1 = lbound(filter(i,j)%cutoff(k)%row_ptr,dim=1)
              m2 = ubound(filter(i,j)%cutoff(k)%row_ptr,dim=1)
              !$ACC ENTER DATA COPYIN(filter(i,j)%cutoff(k)%row_ptr(m1:m2))
            end if
            !
            if (allocated(filter(i,j)%cutoff(k)%row_idx)) then
              !$ACC ENTER DATA COPYIN(filter(i,j)%cutoff(k)%row_idx)
              m1 = lbound(filter(i,j)%cutoff(k)%row_idx,dim=1)
              m2 = ubound(filter(i,j)%cutoff(k)%row_idx,dim=1)
              !$ACC ENTER DATA COPYIN(filter(i,j)%cutoff(k)%row_idx(m1:m2))
            end if
            !
            if (allocated(filter(i,j)%cutoff(k)%row_val)) then
              !$ACC ENTER DATA COPYIN(filter(i,j)%cutoff(k)%row_val)
              m1 = lbound(filter(i,j)%cutoff(k)%row_val,dim=1)
              m2 = ubound(filter(i,j)%cutoff(k)%row_val,dim=1)
              !$ACC ENTER DATA COPYIN(filter(i,j)%cutoff(k)%row_val(m1:m2))
            end if
            !
            if (allocated(filter(i,j)%cutoff(k)%col_ptr)) then
              !$ACC ENTER DATA COPYIN(filter(i,j)%cutoff(k)%col_ptr)
              m1 = lbound(filter(i,j)%cutoff(k)%col_ptr,dim=1)
              m2 = ubound(filter(i,j)%cutoff(k)%col_ptr,dim=1)
              !$ACC ENTER DATA COPYIN(filter(i,j)%cutoff(k)%col_ptr(m1:m2))
            end if
            !
            if (allocated(filter(i,j)%cutoff(k)%col_idx)) then
              !$ACC ENTER DATA COPYIN(filter(i,j)%cutoff(k)%col_idx)
              m1 = lbound(filter(i,j)%cutoff(k)%col_idx,dim=1)
              m2 = ubound(filter(i,j)%cutoff(k)%col_idx,dim=1)
              !$ACC ENTER DATA COPYIN(filter(i,j)%cutoff(k)%col_idx(m1:m2))
            end if
            !
            if (allocated(filter(i,j)%cutoff(k)%col_val)) then
              !$ACC ENTER DATA COPYIN(filter(i,j)%cutoff(k)%col_val)
              m1 = lbound(filter(i,j)%cutoff(k)%col_val,dim=1)
              m2 = ubound(filter(i,j)%cutoff(k)%col_val,dim=1)
              !$ACC ENTER DATA COPYIN(filter(i,j)%cutoff(k)%col_val(m1:m2))
            end if
            !
          end do
          !
        end if
        !
        !---------------------------------------------------------------
        !
        if (allocated(filter(i,j)%sponge)) then
          !$ACC ENTER DATA COPYIN(filter(i,j)%sponge)
          if (allocated(filter(i,j)%sponge%mat)) then
            !$ACC ENTER DATA COPYIN(filter(i,j)%sponge%mat)
            m1 = lbound(filter(i,j)%sponge%mat,dim=1)
            m2 = ubound(filter(i,j)%sponge%mat,dim=1)
            n1 = lbound(filter(i,j)%sponge%mat,dim=2)
            n2 = ubound(filter(i,j)%sponge%mat,dim=2)
            !$ACC ENTER DATA COPYIN(filter(i,j)%sponge%mat(m1:m2,n1:n2))
          end if
          !
          if (allocated(filter(i,j)%sponge%row_ptr)) then
            !$ACC ENTER DATA COPYIN(filter(i,j)%sponge%row_ptr)
            m1 = lbound(filter(i,j)%sponge%row_ptr,dim=1)
            m2 = ubound(filter(i,j)%sponge%row_ptr,dim=1)
            !$ACC ENTER DATA COPYIN(filter(i,j)%sponge%row_ptr(m1:m2))
          end if
          !
          if (allocated(filter(i,j)%sponge%row_idx)) then
            !$ACC ENTER DATA COPYIN(filter(i,j)%sponge%row_idx)
            m1 = lbound(filter(i,j)%sponge%row_idx,dim=1)
            m2 = ubound(filter(i,j)%sponge%row_idx,dim=1)
            !$ACC ENTER DATA COPYIN(filter(i,j)%sponge%row_idx(m1:m2))
          end if
          !
          if (allocated(filter(i,j)%sponge%row_val)) then
            !$ACC ENTER DATA COPYIN(filter(i,j)%sponge%row_val)
            m1 = lbound(filter(i,j)%sponge%row_val,dim=1)
            m2 = ubound(filter(i,j)%sponge%row_val,dim=1)
            !$ACC ENTER DATA COPYIN(filter(i,j)%sponge%row_val(m1:m2))
          end if
          !
          if (allocated(filter(i,j)%sponge%col_ptr)) then
            !$ACC ENTER DATA COPYIN(filter(i,j)%sponge%col_ptr)
            m1 = lbound(filter(i,j)%sponge%col_ptr,dim=1)
            m2 = ubound(filter(i,j)%sponge%col_ptr,dim=1)
            !$ACC ENTER DATA COPYIN(filter(i,j)%sponge%col_ptr(m1:m2))
          end if
          !
          if (allocated(filter(i,j)%sponge%col_idx)) then
            !$ACC ENTER DATA COPYIN(filter(i,j)%sponge%col_idx)
            m1 = lbound(filter(i,j)%sponge%col_idx,dim=1)
            m2 = ubound(filter(i,j)%sponge%col_idx,dim=1)
            !$ACC ENTER DATA COPYIN(filter(i,j)%sponge%col_idx(m1:m2))
          end if
          !
          if (allocated(filter(i,j)%sponge%col_val)) then
            !$ACC ENTER DATA COPYIN(filter(i,j)%sponge%col_val)
            m1 = lbound(filter(i,j)%sponge%col_val,dim=1)
            m2 = ubound(filter(i,j)%sponge%col_val,dim=1)
            !$ACC ENTER DATA COPYIN(filter(i,j)%sponge%col_val(m1:m2))
          end if
        end if
        !
        !---------------------------------------------------------------
        !
      end do
    end do
  end if
  !
end subroutine move_filter_to_device
!
!###############################################################################
!
#include "Functions/debug_timer.f90"
!
!###############################################################################
!
include "Functions/solpts2order.f90"
!
!###############################################################################
!
end module filter_mod
!
!###############################################################################
!
