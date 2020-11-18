module filter_mod
  !
!!!#define PROFILE_ON
  !
  !.. Use Statements ..
  use module_kind_types, global_debug_timer => debug_timer
  !
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq
  use eqn_idx, only : ncs
  !
  use generic_types_mod, only : matrix => matrix_gpu_t
  use generic_types_mod, only : create_matrix_sparsity
  use generic_types_mod, only : output_sparsity
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
  public :: output_filter_dt
  public :: input_filter_dt
  !
  ! Module public arrays
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
subroutine init_filter_matrices()
  !
  !.. Use Statements ..
  use ovar,      only : output_filter_diagonals
  use order_mod, only : geom_solpts
  use order_mod, only : n_min_geom,n_max_geom
  use order_mod, only : n_min_order,n_max_order
  use quadrature_mod, only : geom_is_used
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
            filter(this_geom,this_order)%limiter%mat = &
               FilterMatrix_Hexa( vand(this_geom,this_order)%modal2nodal%mat , &
                                  vand(this_geom,this_order)%nodal2modal%mat , &
                                  cutoff_order=-1 )
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
            filter(this_geom,this_order)%cutoff(cutoff_order)%mat = &
               FilterMatrix_Quad( vand(this_geom,this_order)%modal2nodal%mat , &
                                  vand(this_geom,this_order)%nodal2modal%mat , &
                                  cutoff_order )
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
            filter(this_geom,this_order)%cutoff(cutoff_order)%mat = &
              FilterMatrix_Hexa( vand(this_geom,this_order)%modal2nodal%mat , &
                                 vand(this_geom,this_order)%nodal2modal%mat , &
                                 cutoff_order )
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
subroutine filter_usp()
  !
  !.. Use Statements ..
  use flowvar, only : usp
  use ovar,    only : force_nd_subroutines
  !
  !.. Local Scalars ..
  integer :: nr
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
#ifdef DEBUG_ON
  call debug_timer(leaving_procedure,pname)
#endif
  !
end subroutine filter_usp
!
!###############################################################################
!
pure subroutine filter_usp_3d(usp)
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
pure subroutine filter_usp_2d(usp)
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
pure subroutine filter_usp_1d(usp)
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
pure subroutine filter_usp_nd(mbeg,mend,usp)
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
pure function FilterMatrix_Edge(vand,vand_inv,cutoff_order) result(return_value)
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
      Filter_Cutoff = nint( max(zero,min(one,Sponge_Filter_etac)) * &
                            real(order,kind=wp) )
      local_Filter_Option = 1
      feps_qp = qone
      falpha_qp = real( Sponge_Filter_Alpha , kind=qp )
      filt_order = Sponge_Filter_Order
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = nint( max(zero,min(one,Limiter_Filter_etac)) * &
                      real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp = real( Filter_etac , kind=qp )
    Filter_Cutoff = nint( max(zero,min(one,Filter_etac)) * real(order,kind=wp) )
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
      Filter_Cutoff = nint( max(zero,min(one,Sponge_Filter_etac)) * &
                            real(order,kind=wp) )
      local_Filter_Option = 1
      feps_qp = qone
      falpha_qp = real( Sponge_Filter_Alpha , kind=qp )
      filt_order = Sponge_Filter_Order
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = nint( max(zero,min(one,Limiter_Filter_etac)) * &
                      real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp = real( Filter_etac , kind=qp )
    Filter_Cutoff = nint( max(zero,min(one,Filter_etac)) * real(order,kind=wp) )
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
 pure &
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
      Filter_Cutoff = nint( max(zero,min(one,Sponge_Filter_etac)) * &
                            real(order,kind=wp) )
      local_Filter_Option = 1
      feps_qp = qone
      falpha_qp = real( Sponge_Filter_Alpha , kind=qp )
      filt_order = Sponge_Filter_Order
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = nint( max(zero,min(one,Limiter_Filter_etac)) * &
                      real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp  = real( Filter_etac , kind=qp )
    Filter_Cutoff = nint( max(zero,min(one,Filter_etac)) * real(order,kind=wp) )
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
pure function FilterMatrix_Tria(vand,vand_inv,cutoff_order) result(return_value)
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
      local_Filter_Option = Cutoff_Filter_Option
      Filter_Cutoff = cutoff_order
      fetac_qp = real(Filter_Cutoff,kind=qp) / real(order,kind=qp)
      feps_qp = qzero
    else if (cutoff_order == -2) then
      fetac_qp = real( Sponge_Filter_etac , kind=qp )
      Filter_Cutoff = nint( max(zero,min(one,Sponge_Filter_etac)) * &
                            real(order,kind=wp) )
      local_Filter_Option = 1
      feps_qp = qone
      falpha_qp = real( Sponge_Filter_Alpha , kind=qp )
      filt_order = Sponge_Filter_Order
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = nint( max(zero,min(one,Limiter_Filter_etac)) * &
                      real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp  = real( Filter_etac , kind=qp )
    Filter_Cutoff = nint( max(zero,min(one,Filter_etac)) * real(order,kind=wp) )
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
  if (local_Filter_Option /= 0) then
    !
    k = 1
    do i = 0,order
      do j = 0,order-i
        !
        if (i+j >= Filter_Cutoff) then
          !
          rn_qp = real(i,kind=qp)/real(order,kind=qp) - fetac_qp
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
pure function FilterMatrix_Tetr(vand,vand_inv,cutoff_order) result(return_value)
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
      local_Filter_Option = Cutoff_Filter_Option
      Filter_Cutoff = cutoff_order
      fetac_qp = real(Filter_Cutoff,kind=qp) / real(order,kind=qp)
      feps_qp = qzero
    else if (cutoff_order == -2) then
      fetac_qp = real( Sponge_Filter_etac , kind=qp )
      Filter_Cutoff = nint( max(zero,min(one,Sponge_Filter_etac)) * &
                            real(order,kind=wp) )
      local_Filter_Option = 1
      feps_qp = qone
      falpha_qp = real( Sponge_Filter_Alpha , kind=qp )
      filt_order = Sponge_Filter_Order
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = nint( max(zero,min(one,Limiter_Filter_etac)) * &
                      real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp  = real( Filter_etac , kind=qp )
    Filter_Cutoff = nint( max(zero,min(one,Filter_etac)) * real(order,kind=wp) )
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
  if (local_Filter_Option /= 0) then
    !
    k = 1
    do i = 0,order
      do j = 0,order-i
        !
        if (i+j >= Filter_Cutoff) then
          !
          rn_qp = real(i,kind=qp)/real(order,kind=qp) - fetac_qp
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
end function FilterMatrix_Tetr
!
!###############################################################################
!
pure function FilterMatrix_Pyra(vand,vand_inv,cutoff_order) result(return_value)
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
      local_Filter_Option = Cutoff_Filter_Option
      Filter_Cutoff = cutoff_order
      fetac_qp = real(Filter_Cutoff,kind=qp) / real(order,kind=qp)
      feps_qp = qzero
    else if (cutoff_order == -2) then
      fetac_qp = real( Sponge_Filter_etac , kind=qp )
      Filter_Cutoff = nint( max(zero,min(one,Sponge_Filter_etac)) * &
                            real(order,kind=wp) )
      local_Filter_Option = 1
      feps_qp = qone
      falpha_qp = real( Sponge_Filter_Alpha , kind=qp )
      filt_order = Sponge_Filter_Order
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = nint( max(zero,min(one,Limiter_Filter_etac)) * &
                      real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp  = real( Filter_etac , kind=qp )
    Filter_Cutoff = nint( max(zero,min(one,Filter_etac)) * real(order,kind=wp) )
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
  if (local_Filter_Option /= 0) then
    !
    k = 1
    do i = 0,order
      do j = 0,order-i
        !
        if (i+j >= Filter_Cutoff) then
          !
          rn_qp = real(i,kind=qp)/real(order,kind=qp) - fetac_qp
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
end function FilterMatrix_Pyra
!
!###############################################################################
!
pure function FilterMatrix_Pris(vand,vand_inv,cutoff_order) result(return_value)
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
      local_Filter_Option = Cutoff_Filter_Option
      Filter_Cutoff = cutoff_order
      fetac_qp = real(Filter_Cutoff,kind=qp) / real(order,kind=qp)
      feps_qp = qzero
    else if (cutoff_order == -2) then
      fetac_qp = real( Sponge_Filter_etac , kind=qp )
      Filter_Cutoff = nint( max(zero,min(one,Sponge_Filter_etac)) * &
                            real(order,kind=wp) )
      local_Filter_Option = 1
      feps_qp = qone
      falpha_qp = real( Sponge_Filter_Alpha , kind=qp )
      filt_order = Sponge_Filter_Order
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = nint( max(zero,min(one,Limiter_Filter_etac)) * &
                      real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp  = real( Filter_etac , kind=qp )
    Filter_Cutoff = nint( max(zero,min(one,Filter_etac)) * real(order,kind=wp) )
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
  if (local_Filter_Option /= 0) then
    !
    k = 1
    do i = 0,order
      do j = 0,order-i
        !
        if (i+j >= Filter_Cutoff) then
          !
          rn_qp = real(i,kind=qp)/real(order,kind=qp) - fetac_qp
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
      Filter_Cutoff = nint( max(zero,min(one,Sponge_Filter_etac)) * &
                            real(order,kind=wp) )
      local_Filter_Option = 1
      feps_qp = qone
      falpha_qp = real( Sponge_Filter_Alpha , kind=qp )
      filt_order = Sponge_Filter_Order
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = nint( max(zero,min(one,Limiter_Filter_etac)) * &
                      real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp  = real( Filter_etac , kind=qp )
    Filter_Cutoff = nint( max(zero,min(one,Filter_etac)) * real(order,kind=wp) )
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
subroutine filter_memory_usage(iunit)
  !
  !.. Use Statements ..
  use module_memory
#ifndef DISABLE_DTIO
  use generic_types_mod, only : matrix_memory_usage
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
subroutine output_filter_dt(io)
  !
  !.. Use Statements ..
  use generic_types_mod, only : output_matrix_dt
  !
  !.. Formal Arguments ..
  integer, intent(in) :: io
  !
  !.. Local Scalars ..
  integer :: i,i1,i2
  integer :: o,o1,o2
  integer :: g,g1,g2
  !
continue
  !
  !########################################################
 !type, public :: matrix
 !  real(wp), allocatable :: mat(:,:)
 !  integer,  allocatable :: row_ptr(:)
 !  integer,  allocatable :: row_idx(:)
 !  real(wp), allocatable :: row_val(:)
 !  integer,  allocatable :: col_ptr(:)
 !  integer,  allocatable :: col_idx(:)
 !  real(wp), allocatable :: col_val(:)
 !end type matrix
 !type, public :: filter_t
 !  type(matrix), allocatable :: exponential
 !  type(matrix), allocatable :: limiter
 !  type(matrix), allocatable :: cutoff(:)
 !end type filter_t
 !type(filter_t), public, save, allocatable :: filter(:,:)
  !########################################################
  !
  if (allocated(filter)) then
    !
    g1 = lbound(filter,dim=1)
    g2 = ubound(filter,dim=1)
    !
    o1 = lbound(filter,dim=2)
    o2 = ubound(filter,dim=2)
    !
    write (io) g1,g2,o1,o2
    !
    do o = o1,o2
      do g = g1,g2
        !
        if (allocated(filter(g,o)%exponential)) then
          write (io) 1
          call output_matrix_dt(io,filter(g,o)%exponential)
        else
          write (io) -1
        end if
        !
        if (allocated(filter(g,o)%limiter)) then
          write (io) 1
          call output_matrix_dt(io,filter(g,o)%limiter)
        else
          write (io) -1
        end if
        !
        if (allocated(filter(g,o)%cutoff)) then
          i1 = lbound(filter(g,o)%cutoff,dim=1)
          i2 = ubound(filter(g,o)%cutoff,dim=1)
          write (io) i1,i2
          do i = i1,i2
            call output_matrix_dt(io,filter(g,o)%cutoff(i))
          end do
        else
          write (io) -1,-1
        end if
        !
      end do
    end do
    !
  else
    !
    write (io) -1,-1,-1,-1
    !
  end if
  !
end subroutine output_filter_dt
!
!###############################################################################
!
subroutine input_filter_dt(io,success)
  !
  !.. Use Statements ..
  use generic_types_mod, only : input_matrix_dt
  !
  !.. Formal Arguments ..
  integer,      intent(in) :: io
  logical(lk), intent(out) :: success
  !
  !.. Local Scalars ..
  integer :: i,i1,i2
  integer :: o,o1,o2
  integer :: g,g1,g2
  logical(lk) :: local_success
  !
continue
  !
  !########################################################
 !type, public :: matrix
 !  real(wp), allocatable :: mat(:,:)
 !  integer,  allocatable :: row_ptr(:)
 !  integer,  allocatable :: row_idx(:)
 !  real(wp), allocatable :: row_val(:)
 !  integer,  allocatable :: col_ptr(:)
 !  integer,  allocatable :: col_idx(:)
 !  real(wp), allocatable :: col_val(:)
 !end type matrix
 !type, public :: filter_t
 !  type(matrix), allocatable :: exponential
 !  type(matrix), allocatable :: limiter
 !  type(matrix), allocatable :: cutoff(:)
 !end type filter_t
 !type(filter_t), public, save, allocatable :: filter(:,:)
  !########################################################
  !
  success = fals
  !
  read (io) g1,g2,o1,o2
  !
  if (all([g1,g2,o1,o2] == -1)) return
  if (g1 /= lbound(filter,dim=1)) return
  if (g2 /= ubound(filter,dim=1)) return
  if (o1 /= lbound(filter,dim=2)) return
  if (o2 /= ubound(filter,dim=2)) return
  !
  do o = o1,o2
    do g = g1,g2
      !
      read (io) i1
      if (i1 == -1 .and. allocated(filter(g,o)%exponential)) return
      if (i1 == 1) then
        if (.not. allocated(filter(g,o)%exponential)) return
        call input_matrix_dt(io,filter(g,o)%exponential,local_success)
        if (.not. local_success) return
      end if
      !
      read (io) i1
      if (i1 == -1 .and. allocated(filter(g,o)%limiter)) return
      if (i1 == 1) then
        if (.not. allocated(filter(g,o)%limiter)) return
        call input_matrix_dt(io,filter(g,o)%limiter,local_success)
        if (.not. local_success) return
      end if
      !
      read (io) i1,i2
      if (all([i1,i2] == -1)) then
        if (allocated(filter(g,o)%cutoff)) then
          return
        else
          cycle
        end if
      else
        if (.not. allocated(filter(g,o)%cutoff)) then
          return
        else
          if (i1 /= lbound(filter(g,o)%cutoff,dim=1)) return
          if (i2 /= ubound(filter(g,o)%cutoff,dim=1)) return
          do i = i1,i2
            call input_matrix_dt(io,filter(g,o)%cutoff(i),local_success)
            if (.not. local_success) return
          end do
        end if
      end if
      !
    end do
  end do
  !
  success = true
  !
end subroutine input_filter_dt
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
