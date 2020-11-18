module module_limiters
  !
  !.. Use Statements ..
  use module_kind_types, global_debug_timer => debug_timer
  !
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq
  use eqn_idx, only : ncs,npres
  !
  use ovar, only : force_nd_subroutines
  !
  implicit none
  !
  private
  !
  !
  !.. Public Procedures ..
  !
  public :: initialize_limiters
  public :: limit_by_projecting
  public :: limiters_memory_usage
  public :: resolution_indicator
  public :: jump_indicator
  !
  !
  !.. Module Local Scalars ..
  !
  logical(lk), save :: using_markers = fals
  logical(lk), save :: using_moment_limiter = fals
  logical(lk), save :: init_all_as_marked = true
  !
  !
  !.. Module Local Allocatable Arrays ..
  !
  real(wp), save, allocatable :: cell_ave_sol(:,:)
  real(wp), save, allocatable :: umodes(:,:)
  !
  logical(ldk), save, allocatable :: troubled_cell(:)
  !
  !
  !.. Module Public Protected Allocatable Arrays ..
  !
  real(wp), save, public, protected, allocatable :: marker(:)
  real(wp), save, public, protected, allocatable :: scaled_marker(:)
  !
  !
  !.. Module Parameters ..
  !
  integer, parameter :: init_troubled_cells = -11
  !
#ifdef TIME_MINIAPPS_FOR_GPU_COMPARISON
  public :: write_lim_timings
  integer,  save :: timings_counter = 0
  integer,  save :: timings_points(1:1000) = 0
  real(wp), save :: timings_beg(1:1000) = zero
  real(wp), save :: timings_end(1:1000) = zero
#endif
  !
contains
!
!###############################################################################
!
#ifdef TIME_MINIAPPS_FOR_GPU_COMPARISON
subroutine write_lim_timings()
  !
  integer :: i,l,n,io,ierr,ncalls
  character(len=100) :: fname
  !
  write (fname,2) mypnum
  open (newunit=io,file=fname,status="replace",action="write",iostat=ierr)
  l = maxval(timings_points(1:timings_counter))
  ncalls = count(timings_points(1:timings_counter) == l)
  write (io,10) mypnum+1,ncpu,ncalls,l
  i = 0
  do n = 1,timings_counter
    if (timings_points(n) == l) then
      i = i + 1
      write (io,1) i,timings_end(n)-timings_beg(n)
    end if
  end do
  close (io,iostat=ierr)
  !
  1 format (i4,2x,es23.15)
  2 format ("lim_timings.",i2.2,".dat")
  10 format ('VARIABLES = "Call #", "Wall Time"',/, &
             'ZONE T="lim: CPU ',i2,' of ',i2, &
             '", I=',i0,', DATAPACKING=POINT',/, &
             'AUXDATA loop_size="'i0,'"')
end subroutine write_lim_timings
#endif
!
!###############################################################################
!
subroutine initialize_limiters
  !
  !.. Use Statements ..
  use geovar, only : ncell,n_totpts
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_limiters"
  !
continue
  !
 !if (using_markers) then
    !
    ! Create the array for storing the cell averaged solutions
    !
 !  allocate ( cell_ave_sol(1:nq,1:ncell) , source=zero , &
 !             stat=ierr , errmsg=error_message )
 !  call alloc_error(pname,"cell_ave_sol",1,__LINE__,__FILE__,ierr, &
 !                   error_message)
    !
    !
    ! Create the array for identifying troubled cells
    !
    allocate ( troubled_cell(1:ncell) , source=.false. , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"troubled_cell",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !
 !end if
  !
  ! Create the array for storing the moments of each cell
  !
 !if (using_moment_limiter) then
 !  !
 !  allocate ( umodes(1:n_totpts,1:nq) , source=zero , &
 !             stat=ierr , errmsg=error_message )
 !  call alloc_error(pname,"umodes",1,__LINE__,__FILE__,ierr,error_message)
 !  !
 !end if
  !
end subroutine initialize_limiters
!
!###############################################################################
!
subroutine compute_cell_ave
  !
  !.. Use Statements ..
  use geovar,          only : ncell,cell
  use flowvar,         only : usp
  use vandermonde_mod, only : vand
  !
  !.. Local Scalars ..
  integer  :: m,n,nc,n0,np
  integer  :: this_geom,this_order
  real(wp) :: a,val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_cell_ave"
  !
continue
  !
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1
    np = cell(nc)%end_sp-n0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    do m = 1,nq
      !
      val = zero
      !
      do n = 1,np
        a = vand(this_geom,this_order)%nodal2modal%mat(1,n)
        val = val + a * usp(n0+n,m)
      end do
      !
      cell_ave_sol(m,nc) = half * val
      !
    end do
    !
  end do
  !
end subroutine compute_cell_ave
!
!###############################################################################
!
subroutine moment_limiter
  !
  !.. Use Statements ..
  use geovar,          only : ncell,cell
  use flowvar,         only : usp
  use vandermonde_mod, only : vand
  !
 !use parallel_mod, only : exchange_modes
  !
  !.. Local Scalars
  integer  :: i,k,m,n,nc,n0,np
  integer  :: this_geom,this_order
  real(wp) :: a,val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "moment_limiter"
  !
continue
  !
  if (init_all_as_marked) troubled_cell(:) = .true.
  !
  ! Get the solution modes in each cell of this partition
  !
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1
    np = cell(nc)%end_sp-n0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    do m = 1,nq
      do n = 1,np
        !
        val = zero
        !
        do i = vand(this_geom,this_order)%nodal2modal%col_ptr(n)+1, &
               vand(this_geom,this_order)%nodal2modal%col_ptr(n+1)
          !
          k = vand(this_geom,this_order)%nodal2modal%col_idx(i)
          a = vand(this_geom,this_order)%nodal2modal%col_val(i)
          !
          val = val + a * usp(n0+k,m)
          !
        end do
        !
        umodes(n,m) = val
        !
      end do
    end do
    !
  end do
  !
end subroutine moment_limiter
!
!###############################################################################
!
subroutine limit_by_projecting
  !
  !.. Use Statements ..
  use order_mod, only : n_max_order,n_order
  use flowvar,   only : usp
  use ovar,      only : Limiter_Filter_Option
  use ovar,      only : output_limiter_statistics,nlim
  !
  !.. Local Scalars ..
  integer     :: ntc,dropped_orders,ierr
  logical(lk) :: there_are_no_troubles
  logical(lk) :: first_output
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "limit_by_projecting"
  !
continue
  !
  first_output = true
  !
  ! Loop from the highest order to the lowest order and project
  ! any marked cells to one degree lower order until there are
  ! no more marked cells or we reach the lowest order.
  !
  ! Initialize the array troubled_cell to identify all troubled cells before
  ! attempting to use the projection limiter.
  !
  call identify_troubled_cells(troubled_cell,init_troubled_cells, &
                               there_are_no_troubles)
  !
  if (output_limiter_statistics) then
    !
    if (.not. allocated(nlim)) then
      allocate ( nlim(0:n_order) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"nlim",1,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    nlim(:) = 0
    !
  end if
  !
  ! If there are no troubled cells to begin with, just exit this subroutine
  !
  if (there_are_no_troubles) return
  !
  ! Loop from the maximum cell order down to first order, projecting each
  ! troubled cell down one order until no more troubled cells remain.
  !
  order_loop: do dropped_orders = 1,n_max_order
    !
    ! Loop through the cells applying the projection limiter to
    ! any troubled cells.
    !
    ! Count the number of troubled cells for this projection order
    ! if outputting the limiter statistics
    !
    if (output_limiter_statistics) then
      nlim( max(0,n_order-dropped_orders) ) = count(troubled_cell)
    end if
    !
#ifdef PROFILE_ON
    ntc = count(troubled_cell)
    if (ncpu > 1) then
      call mpi_allreduce(MPI_IN_PLACE,ntc,1_int_mpi,mpi_inttyp, &
                         MPI_SUM,MPI_COMM_WORLD,mpierr)
    end if
    if (mypnum == glb_root) then
      if (first_output) then
        write (iout,100) ntc,n_max_order-dropped_orders+1
        first_output = fals
      else
        write (iout,101) ntc,n_max_order-dropped_orders+1
      end if
    end if
    100 format (6x,"=> ",i8," troubled cells at P",i0)
    101 format (9x,i8," troubled cells at P",i0)
#endif
    !
    if (Limiter_Filter_Option /= 0) then
      call apply_filter_limiter_to_usp(troubled_cell,usp)
    else
      call apply_projct_limiter_to_usp(troubled_cell,usp,dropped_orders)
    end if
    !
    ! Recheck the troubled cells to see if there are no longer
    ! troubled due to the previous projection.
    !
    call identify_troubled_cells(troubled_cell,dropped_orders, &
                                 there_are_no_troubles)
    !
    ! Exit this loop if there are no more troubled cells remaining
    !
    if (there_are_no_troubles) then
#ifdef PROFILE_ON
      if (Limiter_Filter_Option == 0) then
        if (mypnum == glb_root) write (iout,200)
        200 format (8x,"No more troubled cells!")
      end if
#endif
      exit order_loop
    end if
    !
  end do order_loop
  !
end subroutine limit_by_projecting
!
!###############################################################################
!
subroutine apply_filter_limiter_to_usp(troubled_cell,usp)
  !
  !.. Formal Arguments ..
  logical(ldk), intent(inout) :: troubled_cell(:)
  real(wp),     intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "apply_filter_limiter_to_usp"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = nme-nmb+1
  !
  if (force_nd_subroutines) then
    !
    call apply_filter_limiter_to_usp_nd(1,nq,troubled_cell,usp)
    !
  else
    !
    if (nr == 3) then
      call apply_filter_limiter_to_usp_3d(troubled_cell,usp)
    else if (nr == 2) then
      call apply_filter_limiter_to_usp_2d(troubled_cell,usp)
    else
      call apply_filter_limiter_to_usp_1d(troubled_cell,usp)
    end if
    !
    if (ncs > 0) then
      call apply_filter_limiter_to_usp_nd(1,ncs,troubled_cell,usp)
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine apply_filter_limiter_to_usp
!
!###############################################################################
!
pure &
subroutine apply_filter_limiter_to_usp_3d(troubled_cell,usp)
  !
  !.. Use Statements ..
  use order_mod,  only : maxSP
  use geovar,     only : ncell,cell
  use filter_mod, only : filter
  !
  !.. Formal Arguments ..
  logical(ldk), intent(inout) :: troubled_cell(:)
  real(wp),     intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,filt_ec,filt_mx,filt_my,filt_mz,filt_ee
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "apply_filter_limiter_to_usp_3d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,troubled_cell,filter) &
  !$ACC&              PRIVATE(old_usp)
  filt_cell_loop: do nc = 1,ncell
    !
    ! If the current cell is not troubled, skip to the next one
    !
    if (.not.troubled_cell(nc)) then
      cycle filt_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      old_usp(n,nec) = usp(n0+n,nec)
      old_usp(n,nmx) = usp(n0+n,nmx)
      old_usp(n,nmy) = usp(n0+n,nmy)
      old_usp(n,nmz) = usp(n0+n,nmz)
      old_usp(n,nee) = usp(n0+n,nee)
    end do
    !
    ! Apply the filter limiter by multplying the filter
    ! matrix with the solution in the current cell
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
      do i = filter(this_geom,this_order)%limiter%col_ptr(n)+1, &
             filter(this_geom,this_order)%limiter%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%limiter%col_idx(i)
        a = filter(this_geom,this_order)%limiter%col_val(i)
        !
        filt_ec = filt_ec + a * old_usp(k,nec)
        filt_mx = filt_mx + a * old_usp(k,nmx)
        filt_my = filt_my + a * old_usp(k,nmy)
        filt_mz = filt_mz + a * old_usp(k,nmz)
        filt_ee = filt_ee + a * old_usp(k,nee)
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
    ! Unmark this cell as troubled because we only want to filter once
    !
    if (ncs == 0) then
      troubled_cell(nc) = .false.
    end if
    !
  end do filt_cell_loop
  !
end subroutine apply_filter_limiter_to_usp_3d
!
!###############################################################################
!
pure &
subroutine apply_filter_limiter_to_usp_2d(troubled_cell,usp)
  !
  !.. Use Statements ..
  use order_mod,  only : maxSP
  use geovar,     only : ncell,cell
  use filter_mod, only : filter
  !
  !.. Formal Arguments ..
  logical(ldk), intent(inout) :: troubled_cell(:)
  real(wp),     intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,filt_ec,filt_mx,filt_my,filt_ee
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "apply_filter_limiter_to_usp_2d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,troubled_cell,filter) &
  !$ACC&              PRIVATE(old_usp)
  filt_cell_loop: do nc = 1,ncell
    !
    ! If the current cell is not troubled, skip to the next one
    !
    if (.not.troubled_cell(nc)) then
      cycle filt_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      old_usp(n,nec) = usp(n0+n,nec)
      old_usp(n,nmx) = usp(n0+n,nmx)
      old_usp(n,nmy) = usp(n0+n,nmy)
      old_usp(n,nee) = usp(n0+n,nee)
    end do
    !
    ! Apply the filter limiter by multplying the filter
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      filt_ec = zero
      filt_mx = zero
      filt_my = zero
      filt_ee = zero
      !
      do i = filter(this_geom,this_order)%limiter%col_ptr(n)+1, &
             filter(this_geom,this_order)%limiter%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%limiter%col_idx(i)
        a = filter(this_geom,this_order)%limiter%col_val(i)
        !
        filt_ec = filt_ec + a * old_usp(k,nec)
        filt_mx = filt_mx + a * old_usp(k,nmx)
        filt_my = filt_my + a * old_usp(k,nmy)
        filt_ee = filt_ee + a * old_usp(k,nee)
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
    ! Unmark this cell as troubled because we only want to filter once
    !
    if (ncs == 0) then
      troubled_cell(nc) = .false.
    end if
    !
  end do filt_cell_loop
  !
end subroutine apply_filter_limiter_to_usp_2d
!
!###############################################################################
!
pure &
subroutine apply_filter_limiter_to_usp_1d(troubled_cell,usp)
  !
  !.. Use Statements ..
  use order_mod,  only : maxSP
  use geovar,     only : ncell,cell
  use filter_mod, only : filter
  !
  !.. Formal Arguments ..
  logical(ldk), intent(inout) :: troubled_cell(:)
  real(wp),     intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,filt_ec,filt_mx,filt_ee
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "apply_filter_limiter_to_usp_1d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,troubled_cell,filter) &
  !$ACC&              PRIVATE(old_usp)
  filt_cell_loop: do nc = 1,ncell
    !
    ! If the current cell is not troubled, skip to the next one
    !
    if (.not.troubled_cell(nc)) then
      cycle filt_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      old_usp(n,nec) = usp(n0+n,nec)
      old_usp(n,nmx) = usp(n0+n,nmx)
      old_usp(n,nee) = usp(n0+n,nee)
    end do
    !
    ! Apply the filter limiter by multplying the filter
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      filt_ec = zero
      filt_mx = zero
      filt_ee = zero
      !
      do i = filter(this_geom,this_order)%limiter%col_ptr(n)+1, &
             filter(this_geom,this_order)%limiter%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%limiter%col_idx(i)
        a = filter(this_geom,this_order)%limiter%col_val(i)
        !
        filt_ec = filt_ec + a * old_usp(k,nec)
        filt_mx = filt_mx + a * old_usp(k,nmx)
        filt_ee = filt_ee + a * old_usp(k,nee)
        !
      end do
      !
      usp(n0+n,nec) = filt_ec
      usp(n0+n,nmx) = filt_mx
      usp(n0+n,nee) = filt_ee
      !
    end do
    !
    ! Unmark this cell as troubled because we only want to filter once
    !
    if (ncs == 0) then
      troubled_cell(nc) = .false.
    end if
    !
  end do filt_cell_loop
  !
end subroutine apply_filter_limiter_to_usp_1d
!
!###############################################################################
!
pure &
subroutine apply_filter_limiter_to_usp_nd(mbeg,mend,troubled_cell,usp)
  !
  !.. Use Statements ..
  use order_mod,  only : maxSP
  use geovar,     only : ncell,cell
  use filter_mod, only : filter
  !
  !.. Formal Arguments ..
  integer,         intent(in) :: mbeg,mend
  logical(ldk), intent(inout) :: troubled_cell(:)
  real(wp),     intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,m,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,val
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "apply_filter_limiter_to_usp_nd"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,troubled_cell,filter) &
  !$ACC&              PRIVATE(old_usp)
  filt_cell_loop: do nc = 1,ncell
    !
    ! If the current cell is not troubled, skip to the next one
    !
    if (.not.troubled_cell(nc)) then
      cycle filt_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do m = mbeg,mend
      do n = 1,np
        old_usp(n,m) = usp(n0+n,m)
      end do
    end do
    !
    ! Apply the filter limiter by multplying the filter
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do m = mbeg,mend
      do n = 1,np
        !
        val = zero
        !
        do i = filter(this_geom,this_order)%limiter%col_ptr(n)+1, &
               filter(this_geom,this_order)%limiter%col_ptr(n+1)
          !
          k = filter(this_geom,this_order)%limiter%col_idx(i)
          a = filter(this_geom,this_order)%limiter%col_val(i)
          !
          val = val + a * old_usp(k,m)
          !
        end do
        !
        usp(n0+n,m) = val
        !
      end do
    end do
    !
    ! Unmark this cell as troubled because we only want to filter once
    !
    troubled_cell(nc) = .false.
    !
  end do filt_cell_loop
  !
end subroutine apply_filter_limiter_to_usp_nd
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine apply_projct_limiter_to_usp(troubled_cell,usp,dropped_orders)
  !
  !.. Formal Arguments ..
  integer,         intent(in) :: dropped_orders
  logical(ldk), intent(inout) :: troubled_cell(:)
  real(wp),     intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "apply_projct_limiter_to_usp"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = nme-nmb+1
  !
  if (force_nd_subroutines) then
    !
    call apply_projct_limiter_to_usp_nd(1,nq,troubled_cell,usp,dropped_orders)
    !
  else
    !
    if (nr == 3) then
      call apply_projct_limiter_to_usp_3d(troubled_cell,usp,dropped_orders)
    else if (nr == 2) then
      call apply_projct_limiter_to_usp_2d(troubled_cell,usp,dropped_orders)
    else
      call apply_projct_limiter_to_usp_1d(troubled_cell,usp,dropped_orders)
    end if
    !
    if (ncs > 0) then
      call apply_projct_limiter_to_usp_nd(1,ncs,troubled_cell,usp, &
                                          dropped_orders)
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine apply_projct_limiter_to_usp
!
!###############################################################################
!
pure &
subroutine apply_projct_limiter_to_usp_3d(troubled_cell,usp,dropped_orders)
  !
  !.. Use Statements ..
  use order_mod,      only : maxSP
  use geovar,         only : ncell,cell
  use projection_mod, only : projct
  !
  !.. Formal Arguments ..
  integer,         intent(in) :: dropped_orders
  logical(ldk), intent(inout) :: troubled_cell(:)
  real(wp),     intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  integer  :: proj_order
  real(wp) :: a,plim_ec,plim_mx,plim_my,plim_mz,plim_ee
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "apply_projct_limiter_to_usp_3d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,troubled_cell,projct) &
  !$ACC&              PRIVATE(old_usp)
  proj_cell_loop: do nc = 1,ncell
    !
    ! If the current cell is not troubled, skip to the next one
    !
    if (.not.troubled_cell(nc)) then
      cycle proj_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Logic in subroutine identify_troubled_cells
    ! prevents this from becoming less than 0
    !
    proj_order = this_order - dropped_orders
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      old_usp(n,nec) = usp(n0+n,nec)
      old_usp(n,nmx) = usp(n0+n,nmx)
      old_usp(n,nmy) = usp(n0+n,nmy)
      old_usp(n,nmz) = usp(n0+n,nmz)
      old_usp(n,nee) = usp(n0+n,nee)
    end do
    !
    ! Apply the projection limiter by multplying the projection
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      plim_ec = zero
      plim_mx = zero
      plim_my = zero
      plim_mz = zero
      plim_ee = zero
      !
      do i = projct(this_geom,this_order)%limiter(proj_order)%row_ptr(n)+1, &
             projct(this_geom,this_order)%limiter(proj_order)%row_ptr(n+1)
        !
        k = projct(this_geom,this_order)%limiter(proj_order)%row_idx(i)
        a = projct(this_geom,this_order)%limiter(proj_order)%row_val(i)
        !
        plim_ec = plim_ec + a * old_usp(k,nec)
        plim_mx = plim_mx + a * old_usp(k,nmx)
        plim_my = plim_my + a * old_usp(k,nmy)
        plim_mz = plim_mz + a * old_usp(k,nmz)
        plim_ee = plim_ee + a * old_usp(k,nee)
        !
      end do
      !
      usp(n0+n,nec) = plim_ec
      usp(n0+n,nmx) = plim_mx
      usp(n0+n,nmy) = plim_my
      usp(n0+n,nmz) = plim_mz
      usp(n0+n,nee) = plim_ee
      !
    end do
    !
    ! Unmark this cell as troubled if proj_order is 0
    !
    if (proj_order == 0 .and. ncs == 0) then
      troubled_cell(nc) = .false.
    end if
    !
  end do proj_cell_loop
  !
end subroutine apply_projct_limiter_to_usp_3d
!
!###############################################################################
!
pure &
subroutine apply_projct_limiter_to_usp_2d(troubled_cell,usp,dropped_orders)
  !
  !.. Use Statements ..
  use order_mod,      only : maxSP
  use geovar,         only : ncell,cell
  use projection_mod, only : projct
  !
  !.. Formal Arguments ..
  integer,         intent(in) :: dropped_orders
  logical(ldk), intent(inout) :: troubled_cell(:)
  real(wp),     intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  integer  :: proj_order
  real(wp) :: a,plim_ec,plim_mx,plim_my,plim_ee
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "apply_projct_limiter_to_usp_2d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,troubled_cell,projct) &
  !$ACC&              PRIVATE(old_usp)
  proj_cell_loop: do nc = 1,ncell
    !
    ! If the current cell is not troubled, skip to the next one
    !
    if (.not.troubled_cell(nc)) then
      cycle proj_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Logic in subroutine identify_troubled_cells
    ! prevents this from becoming less than 0
    !
    proj_order = this_order - dropped_orders
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      old_usp(n,nec) = usp(n0+n,nec)
      old_usp(n,nmx) = usp(n0+n,nmx)
      old_usp(n,nmy) = usp(n0+n,nmy)
      old_usp(n,nee) = usp(n0+n,nee)
    end do
    !
    ! Apply the projection limiter by multplying the projection
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      plim_ec = zero
      plim_mx = zero
      plim_my = zero
      plim_ee = zero
      !
      do i = projct(this_geom,this_order)%limiter(proj_order)%row_ptr(n)+1, &
             projct(this_geom,this_order)%limiter(proj_order)%row_ptr(n+1)
        !
        k = projct(this_geom,this_order)%limiter(proj_order)%row_idx(i)
        a = projct(this_geom,this_order)%limiter(proj_order)%row_val(i)
        !
        plim_ec = plim_ec + a * old_usp(k,nec)
        plim_mx = plim_mx + a * old_usp(k,nmx)
        plim_my = plim_my + a * old_usp(k,nmy)
        plim_ee = plim_ee + a * old_usp(k,nee)
        !
      end do
      !
      usp(n0+n,nec) = plim_ec
      usp(n0+n,nmx) = plim_mx
      usp(n0+n,nmy) = plim_my
      usp(n0+n,nee) = plim_ee
      !
    end do
    !
    ! Unmark this cell as troubled if proj_order is 0
    !
    if (proj_order == 0 .and. ncs == 0) then
      troubled_cell(nc) = .false.
    end if
    !
  end do proj_cell_loop
  !
end subroutine apply_projct_limiter_to_usp_2d
!
!###############################################################################
!
pure &
subroutine apply_projct_limiter_to_usp_1d(troubled_cell,usp,dropped_orders)
  !
  !.. Use Statements ..
  use order_mod,      only : maxSP
  use geovar,         only : ncell,cell
  use projection_mod, only : projct
  !
  !.. Formal Arguments ..
  integer,         intent(in) :: dropped_orders
  logical(ldk), intent(inout) :: troubled_cell(:)
  real(wp),     intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  integer  :: proj_order
  real(wp) :: a,plim_ec,plim_mx,plim_ee
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "apply_projct_limiter_to_usp_1d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,troubled_cell,projct) &
  !$ACC&              PRIVATE(old_usp)
  proj_cell_loop: do nc = 1,ncell
    !
    ! If the current cell is not troubled, skip to the next one
    !
    if (.not.troubled_cell(nc)) then
      cycle proj_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Logic in subroutine identify_troubled_cells
    ! prevents this from becoming less than 0
    !
    proj_order = this_order - dropped_orders
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      old_usp(n,nec) = usp(n0+n,nec)
      old_usp(n,nmx) = usp(n0+n,nmx)
      old_usp(n,nee) = usp(n0+n,nee)
    end do
    !
    ! Apply the projection limiter by multplying the projection
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      plim_ec = zero
      plim_mx = zero
      plim_ee = zero
      !
      do i = projct(this_geom,this_order)%limiter(proj_order)%row_ptr(n)+1, &
             projct(this_geom,this_order)%limiter(proj_order)%row_ptr(n+1)
        !
        k = projct(this_geom,this_order)%limiter(proj_order)%row_idx(i)
        a = projct(this_geom,this_order)%limiter(proj_order)%row_val(i)
        !
        plim_ec = plim_ec + a * old_usp(k,nec)
        plim_mx = plim_mx + a * old_usp(k,nmx)
        plim_ee = plim_ee + a * old_usp(k,nee)
        !
      end do
      !
      usp(n0+n,nec) = plim_ec
      usp(n0+n,nmx) = plim_mx
      usp(n0+n,nee) = plim_ee
      !
    end do
    !
    ! Unmark this cell as troubled if proj_order is 0
    !
    if (proj_order == 0 .and. ncs == 0) then
      troubled_cell(nc) = .false.
    end if
    !
  end do proj_cell_loop
  !
end subroutine apply_projct_limiter_to_usp_1d
!
!###############################################################################
!
pure &
subroutine apply_projct_limiter_to_usp_nd(mbeg,mend,troubled_cell,usp, &
                                          dropped_orders)
  !
  !.. Use Statements ..
  use order_mod,      only : maxSP
  use geovar,         only : ncell,cell
  use projection_mod, only : projct
  !
  !.. Formal Arguments ..
  integer,         intent(in) :: mbeg,mend
  integer,         intent(in) :: dropped_orders
  logical(ldk), intent(inout) :: troubled_cell(:)
  real(wp),     intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,m,n,nc,np,n0
  integer  :: this_geom,this_order
  integer  :: proj_order
  real(wp) :: a,val
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "apply_projct_limiter_to_usp_nd"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,troubled_cell,projct) &
  !$ACC&              PRIVATE(old_usp)
  proj_cell_loop: do nc = 1,ncell
    !
    ! If the current cell is not troubled, skip to the next one
    !
    if (.not.troubled_cell(nc)) then
      cycle proj_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Logic in subroutine identify_troubled_cells
    ! prevents this from becoming less than 0
    !
    proj_order = this_order - dropped_orders
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do m = mbeg,mend
      do n = 1,np
        old_usp(n,m) = usp(n0+n,m)
      end do
    end do
    !
    ! Apply the projection limiter by multplying the projection
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do m = mbeg,mend
      do n = 1,np
        !
        val = zero
        !
        do i = projct(this_geom,this_order)%limiter(proj_order)%row_ptr(n)+1, &
               projct(this_geom,this_order)%limiter(proj_order)%row_ptr(n+1)
          !
          k = projct(this_geom,this_order)%limiter(proj_order)%row_idx(i)
          a = projct(this_geom,this_order)%limiter(proj_order)%row_val(i)
          !
          val = val + a * old_usp(k,m)
          !
        end do
        !
        usp(n0+n,m) = val
        !
      end do
    end do
    !
    ! Unmark this cell as troubled if proj_order is 0
    !
    if (proj_order == 0) then
      troubled_cell(nc) = .false.
    end if
    !
  end do proj_cell_loop
  !
end subroutine apply_projct_limiter_to_usp_nd
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine identify_troubled_cells(troubled_cell,dropped_orders, &
                                   there_are_no_troubles)
  !
  !.. Use Statements ..
  use geovar,    only : ncell,cell
  use ovar,      only : Limiter_Option
  use order_mod, only : p_order,n_order
  !
  !.. Formal Arguments ..
  integer,         intent(in) :: dropped_orders
  logical(lk),    intent(out) :: there_are_no_troubles
  logical(ldk), intent(inout) :: troubled_cell(:)
  !
  !.. Local Scalars ..
  integer :: nc,proj_order
  integer :: orders_to_cutoff
  integer :: ord_mlt,ord_add
  integer :: actual_order,adjusted_order
  !
  logical(lk) :: there_are_troubles
  !
continue
  !
  proj_order = dropped_orders
  !
  ! Check if troubled_cell needs to be initialized. If it does, mark
  ! all cells as troubled before actually checking the cells. Also set
  ! proj_order to zero so there are no problems when checking if the
  ! cell solution has been projected down to a first order solution.
  !
  if (dropped_orders == init_troubled_cells) then
    !
    proj_order = 0
    !
    !$ACC PARALLEL LOOP PRESENT(troubled_cell)
    do nc = 1,ncell
      troubled_cell(nc) = .true.
    end do
    !
  end if
  !
  ! Compute the troubled marker
  !
  if (Limiter_Option == 1) then
   !orders_to_cutoff = proj_order+1
    orders_to_cutoff = proj_order
    call resolution_indicator(troubled_cell=troubled_cell, &
                              orders_to_cutoff=orders_to_cutoff)
  else
    call jump_indicator(troubled_cell=troubled_cell)
  end if
  !
  ord_mlt = 1
  ord_add = 0
  if (p_order < n_order) then
    ord_mlt = 0
    ord_add = p_order
  end if
  !
  ! Now loop through the cells and determine if there are any remaining
  ! troubled cells. Unmark any cells that the limiters will try to
  ! limit below P1, i.e. a 2nd order or linear polynomial function
  !
  there_are_troubles = .false.
  !
  !$ACC PARALLEL LOOP REDUCTION(.or.:there_are_troubles) &
  !$ACC&              PRESENT(cell,troubled_cell)
  do nc = 1,ncell
    !
    ! If the current cell is not troubled, skip to the next one
    !
    if (.not. troubled_cell(nc)) cycle
    !
    actual_order = ord_mlt*cell(nc)%order + ord_add
    !
    ! If the projection has reduced the solution to first order, mark
    ! this cell as no longer troubled and skip to the next cell
    !
    adjusted_order = actual_order - proj_order
    !
   !if (adjusted_order <= 1) then
    if (adjusted_order < 0 .or. &
        (adjusted_order <= 1 .and. cell(nc)%order > 1)) then
      troubled_cell(nc) = .false.
      cycle
    end if
    !
    ! If the current cell is determined to be troubled, change
    ! there_are_no_troubles to fals meaning that troubles cells exist.
    !
    if (troubled_cell(nc)) then
      there_are_troubles = true
    end if
    !
  end do
  !
  there_are_no_troubles = logical( .not. there_are_troubles , kind=lk )
  !
#ifdef PROFILE_ON
  nc = merge(1,0,there_are_no_troubles)
  call mpi_allreduce(MPI_IN_PLACE,nc,1_int_mpi,mpi_inttyp, &
                     MPI_MIN,MPI_COMM_WORLD,mpierr)
  there_are_no_troubles = (nc == 1)
#endif
  !
end subroutine identify_troubled_cells
!
!###############################################################################
!
subroutine resolution_indicator(apply_projection,troubled_cell, &
                                orders_to_cutoff)
  !
  !.. Use Statements ..
  use order_mod, only : p_order,n_order
  use geovar,    only : ncell,cell
  use flowvar,   only : usp
  use ovar,      only : output_limiter_statistics,nlim
  use ovar,      only : Limiter_Option
  !
  !.. Optional Arguments ..
  logical(lk),  optional,    intent(in) :: apply_projection
  logical(ldk), optional, intent(inout) :: troubled_cell(:)
  integer,      optional,    intent(in) :: orders_to_cutoff
  !
  !.. Local Scalars ..
  integer :: n,nc,nr,ord_mlt,ord_add
  integer :: this_order,limit_order
  integer :: num_limited_cells,ierr
  !
  !.. Local Allocatable Arrays ..
  logical(ldk), allocatable :: skip_cell(:)
  integer,      allocatable :: lmtr_order(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "resolution_indicator"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Allocate the module arrays marker and scaled_marker
  ! if they havent been allocated yet
  !
  if (.not. allocated(marker)) then
    allocate ( marker(1:ncell) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"marker",1,__LINE__,__FILE__,ierr,error_message)
    !$ACC UPDATE DEVICE(marker)
  end if
  !
  if (.not. allocated(scaled_marker)) then
    allocate ( scaled_marker(1:ncell) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"scaled_marker",1,__LINE__,__FILE__, &
                     ierr,error_message)
    !$ACC UPDATE DEVICE(scaled_marker)
  end if
  !
  ! Allocate the local worker arrays
  !
  allocate ( skip_cell(1:ncell) )
  allocate ( lmtr_order(1:ncell) )
  !$ACC DATA CREATE(skip_cell,lmtr_order)
  !
  ! Initialize skip_cell to the negation of troubled_cell if it is present.
  ! Otherwise, initialize to false to indicate that all cells will be
  ! checked by the resolution indicator.
  !
  if (present(troubled_cell)) then
    !$ACC PARALLEL LOOP PRESENT(troubled_cell,skip_cell)
    do nc = 1,ncell
      skip_cell(nc) = (.not. troubled_cell(nc))
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(skip_cell)
    do nc = 1,ncell
      skip_cell(nc) = .false.
    end do
  end if
  !
  ! Set the cutoff order for each cell
  !
  if (p_order < n_order) then
    ord_mlt = 0
    ord_add = p_order
  else
    ord_mlt = 1
   !ord_add = 0
    ord_add = -1
  end if
  !
  if (present(apply_projection)) then
   !if (apply_projection) then
      if (Limiter_Option == -2) then
        ord_add = 0
      end if
   !end if
    if (Limiter_Option == 1) then
      if (.not. present(orders_to_cutoff)) then
        ord_add = 0
      end if
    end if
  end if
  !
  !$ACC PARALLEL LOOP PRESENT(cell,lmtr_order)
  do nc = 1,ncell
    lmtr_order(nc) = ord_mlt*cell(nc)%order + ord_add
  end do
  !
  ! Modify the cutoff order if orders_to_cutoff is true
  !
  if (present(orders_to_cutoff)) then
    !$ACC PARALLEL LOOP PRESENT(lmtr_order)
    do nc = 1,ncell
      lmtr_order(nc) = lmtr_order(nc)-orders_to_cutoff+1
    end do
  end if
  !
  ! Set skip_cell to true for any cells that lmtr_order equals zero
  !
  !$ACC PARALLEL LOOP PRESENT(lmtr_order,skip_cell)
  do nc = 1,ncell
    if (lmtr_order(nc) < 0 .or. &
        (lmtr_order(nc) == 0 .and. cell(nc)%order > 1)) then
      skip_cell(nc) = .true.
    end if
  end do
  !
  ! Use the resolution indicator to compute marker and scaled_marker
  !
  call compute_resolution_marker(usp,skip_cell,lmtr_order,marker,scaled_marker)
  !
  ! Update the troubled_cell indicator if it is present
  !
  if (present(troubled_cell)) then
    !$ACC PARALLEL LOOP PRESENT(scaled_marker,troubled_cell,skip_cell)
    do nc = 1,ncell
      troubled_cell(nc) = (scaled_marker(nc) >= 1.0e-8_wp)
      skip_cell(nc) = .not. troubled_cell(nc)
    end do
  end if
  !
  ! Apply the projection limiter to usp if
  ! apply_projection is present and true
  !
  if (present(apply_projection)) then
    if (apply_projection) then
      !
      ! Initialize skip_cell to the negation of troubled_cell if it is present.
      ! Otherwise, initialize to false to indicate that all cells will be
      ! limited if needed.
      !
      if (present(troubled_cell)) then
        !$ACC PARALLEL LOOP PRESENT(troubled_cell,skip_cell)
        do nc = 1,ncell
          skip_cell(nc) = (.not. troubled_cell(nc))
        end do
      else
        !$ACC PARALLEL LOOP PRESENT(skip_cell)
        do nc = 1,ncell
          skip_cell(nc) = .false.
        end do
      end if
      !
      ! Find the limiter order to use for each cell based on scaled_marker.
      ! If the limiter order is the same as the cell order or troubled_cell
      ! has said not to limit cell, mark the cell as not to be limited.
      !
      !$ACC PARALLEL LOOP PRESENT(cell) &
      !$ACC&              PRESENT(scaled_marker) &
      !$ACC&              PRESENT(skip_cell,lmtr_order)
      do nc = 1,ncell
        !
        this_order = cell(nc)%order
        !
        if (skip_cell(nc)) then
          lmtr_order(nc) = this_order
        else
          limit_order = int( real(this_order,kind=wp)*(one-scaled_marker(nc)) )
          if (this_order <= 1 .or. Limiter_Option == -2) then
            lmtr_order(nc) = max(0,min(this_order,limit_order))
          else
            lmtr_order(nc) = max(1,min(this_order,limit_order))
          end if
        end if
        !
        if (lmtr_order(nc) == this_order) then
          skip_cell(nc) = .true.
        end if
        !
      end do
      !
      if (output_limiter_statistics) then
        !
        if (.not. allocated(nlim)) then
          allocate ( nlim(0:n_order) , stat=ierr , errmsg=error_message )
          call alloc_error(pname,"nlim",1,__LINE__,__FILE__,ierr,error_message)
        end if
        !
        nlim(:) = 0
        !
        do n = 0,n_order
          num_limited_cells = 0
          !$ACC PARALLEL LOOP REDUCTION(+:num_limited_cells) &
          !$ACC&              PRESENT(lmtr_order)
          do nc = 1,ncell
            if (.not. skip_cell(nc)) then
              if (lmtr_order(nc) == n) then
                num_limited_cells = num_limited_cells + 1
              end if
            end if
          end do
          nlim(n) = num_limited_cells
        end do
        !
      end if
      !
      ! Apply the projection limiter to the cells marked as needing limiting
      !
      nr = nme-nmb+1
      !
      if (force_nd_subroutines) then
        !
        call res_ind_apply_proj_nd(1,nq,usp,skip_cell,lmtr_order)
        !
      else
        !
        if (nr == 3) then
          call res_ind_apply_proj_3d(usp,skip_cell,lmtr_order)
        else if (nr == 2) then
          call res_ind_apply_proj_2d(usp,skip_cell,lmtr_order)
        else
          call res_ind_apply_proj_1d(usp,skip_cell,lmtr_order)
        end if
        !
        if (ncs > 0) then
          call res_ind_apply_proj_nd(1,ncs,usp,skip_cell,lmtr_order)
        end if
        !
      end if
      !
    end if
  end if
  !
  !$ACC END DATA
  deallocate ( lmtr_order )
  deallocate ( skip_cell )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine resolution_indicator
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
#ifndef TIME_MINIAPPS_FOR_GPU_COMPARISON
 pure &
#endif
subroutine compute_resolution_marker(usp,skip_cell,cutoff_order, &
                                     marker,scaled_marker)
  !
  !.. Use Statements ..
  use ovar, only : Resolution_Marker_Variable
  !
  !.. Formal Arguments ..
  real(wp),     intent(in) :: usp(:,:)
  logical(ldk), intent(in) :: skip_cell(:)
  integer,      intent(in) :: cutoff_order(:)
  real(wp),  intent(inout) :: marker(:)
  real(wp),  intent(inout) :: scaled_marker(:)
  !
  !.. Local Scalars ..
  integer :: nc,nr,nvar
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_resolution_marker"
  !
continue
  !
  nr = nme-nmb+1
  !
  !$ACC PARALLEL LOOP PRESENT(marker)
  do nc = 1,size(marker)
    marker(nc) = -huge(zero)
  end do
  !
  if (Resolution_Marker_Variable > 0) then
    nvar = Resolution_Marker_Variable
  else
    nvar = 0
  end if
  !
  if (nvar == 0) then
    !
    if (force_nd_subroutines) then
      !
      call compute_resolution_marker_nd(1,nq,usp,skip_cell,cutoff_order,marker)
      !
    else
      !
      if (nr == 3) then
        call compute_resolution_marker_3d(usp,skip_cell,cutoff_order,marker)
      else if (nr == 2) then
        call compute_resolution_marker_2d(usp,skip_cell,cutoff_order,marker)
      else
        call compute_resolution_marker_1d(usp,skip_cell,cutoff_order,marker)
      end if
      !
      if (ncs > 0) then
        call compute_resolution_marker_nd(1,ncs,usp,skip_cell, &
                                          cutoff_order,marker)
      end if
      !
    end if
    !
  else
    !
    call compute_resolution_marker_1var(usp(:,nvar),skip_cell, &
                                        cutoff_order,marker)
    !
  end if
  !
  ! Compute the scaled resolution marker
  !
  call compute_scaled_marker(skip_cell,marker,scaled_marker)
  !
end subroutine compute_resolution_marker
!
!###############################################################################
!
subroutine output_lim_data(usp,skip_cell,cutoff_order)
  !
  !.. Use Statements ..
  use geovar,         only : n_solpts
  use geovar,         only : output_cell_dt,input_cell_dt
  use filter_mod,     only : output_filter_dt,input_filter_dt
  use quadrature_mod, only : output_std_elem_dt,input_std_elem_dt
  !
  !.. Formal Arguments ..
  real(wp),     intent(in) :: usp(:,:)
  logical(ldk), intent(in) :: skip_cell(:)
  integer,      intent(in) :: cutoff_order(:)
  !
  !.. Local Scalars ..
  integer :: io,m,n,ierr
  logical(lk) :: success
  !
continue
  !
  open (newunit=io,file="lim_usp.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  write (io) n_solpts
  write (io) ((usp(n,m), n=1,n_solpts), m=1,nq)
  close (io)
  !
  open (newunit=io,file="lim_skip_cell.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  write (io) size(skip_cell)
  write (io) (skip_cell(n), n=1,size(skip_cell))
  close (io)
  !
  open (newunit=io,file="lim_cutoff_order.dat",status="replace",action="write",&
                   form="unformatted",iostat=ierr)
  write (io) size(cutoff_order)
  write (io) (cutoff_order(n), n=1,size(cutoff_order))
  close (io)
  !
  open (newunit=io,file="lim_cell.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  call output_cell_dt(io)
  close (io)
  !
  open (newunit=io,file="lim_filter.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  call output_filter_dt(io)
  close (io)
  !
  open (newunit=io,file="lim_std_elem.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  call output_std_elem_dt(io)
  close (io)
  !
  ! Check to make sure the derived types were output correctly
  !
  open (newunit=io,file="lim_cell.dat",status="old",action="read", &
                   form="unformatted",iostat=ierr)
  call input_cell_dt(io,success)
  close (io)
  if (success) then
    write (iout,*) "cell output correctly!"
  else
    write (iout,*) "ERROR: cell NOT output correctly"
  end if
  !
  open (newunit=io,file="lim_filter.dat",status="old",action="read", &
                   form="unformatted",iostat=ierr)
  call input_filter_dt(io,success)
  close (io)
  if (success) then
    write (iout,*) "filter output correctly!"
  else
    write (iout,*) "ERROR: filter NOT output correctly"
  end if
  !
  open (newunit=io,file="lim_std_elem.dat",status="old",action="read", &
                   form="unformatted",iostat=ierr)
  call input_std_elem_dt(io,success)
  close (io)
  if (success) then
    write (iout,*) "std_elem output correctly!"
  else
    write (iout,*) "ERROR: std_elem NOT output correctly"
  end if
  !
end subroutine output_lim_data
!
!###############################################################################
!
#ifndef TIME_MINIAPPS_FOR_GPU_COMPARISON
 pure &
#endif
subroutine compute_resolution_marker_3d(usp,skip_cell,cutoff_order,marker)
  !
  !.. Use Statements ..
  use order_mod,      only : maxSP
  use geovar,         only : n_solpts,ncell,cell
  use filter_mod,     only : filter
  use quadrature_mod, only : std_elem
  !
  !.. Formal Arguments ..
  real(wp),     intent(in) :: usp(:,:)
  logical(ldk), intent(in) :: skip_cell(:)
  integer,      intent(in) :: cutoff_order(:)
  real(wp),  intent(inout) :: marker(:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order,order
  real(wp) :: a,wt
  real(wp) :: val_ec,l2n_ec,l2d_ec,vn_ec,vd_ec,mark_ec
  real(wp) :: val_mx,l2n_mx,l2d_mx,vn_mx,vd_mx,mark_mx
  real(wp) :: val_my,l2n_my,l2d_my,vn_my,vd_my,mark_my
  real(wp) :: val_mz,l2n_mz,l2d_mz,vn_mz,vd_mz,mark_mz
  real(wp) :: val_ee,l2n_ee,l2d_ee,vn_ee,vd_ee,mark_ee
  !
  !.. Local Allocatable Arrays ..
  real(wp) :: var_p(1:maxSP,1:nq)
  real(wp) :: var_pm1(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_resolution_marker_3d"
  !
continue
  !
#ifdef TIME_MINIAPPS_FOR_GPU_COMPARISON
  timings_counter = timings_counter + 1
  timings_points(timings_counter) = ncell
  call gfr_timer( timings_beg(timings_counter) )
#endif
  !
  ! Compute u_hat, the degree-(P-1) solution of density at the sol-pts
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,filter,std_elem) &
  !$ACC&              PRESENT(usp,marker,skip_cell,cutoff_order) &
  !$ACC&              PRIVATE(var_p,var_pm1)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index of sol-pts for this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts for this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      var_p(n,nec) = usp(n0+n,nec)
      var_p(n,nmx) = usp(n0+n,nmx)
      var_p(n,nmy) = usp(n0+n,nmy)
      var_p(n,nmz) = usp(n0+n,nmz)
      var_p(n,nee) = usp(n0+n,nee)
    end do
    !
    if (skip_cell(nc)) then
      !
      !$ACC LOOP VECTOR
      do n = 1,np
        var_pm1(n,nec) = var_p(n,nec)
        var_pm1(n,nmx) = var_p(n,nmx)
        var_pm1(n,nmy) = var_p(n,nmy)
        var_pm1(n,nmz) = var_p(n,nmz)
        var_pm1(n,nee) = var_p(n,nee)
      end do
      !
    else
      !
      order = cutoff_order(nc) ! cutoff order for the filter
      !
      !$ACC LOOP VECTOR
      do n = 1,np
        !
        val_ec = zero
        val_mx = zero
        val_my = zero
        val_mz = zero
        val_ee = zero
        !
        do i = filter(this_geom,this_order)%cutoff(order)%col_ptr(n)+1, &
               filter(this_geom,this_order)%cutoff(order)%col_ptr(n+1)
          !
          k = filter(this_geom,this_order)%cutoff(order)%col_idx(i)
          a = filter(this_geom,this_order)%cutoff(order)%col_val(i)
          !
          val_ec = val_ec + a * var_p(k,nec)
          val_mx = val_mx + a * var_p(k,nmx)
          val_my = val_my + a * var_p(k,nmy)
          val_mz = val_mz + a * var_p(k,nmz)
          val_ee = val_ee + a * var_p(k,nee)
          !
        end do
        !
        var_pm1(n,nec) = val_ec
        var_pm1(n,nmx) = val_mx
        var_pm1(n,nmy) = val_my
        var_pm1(n,nmz) = val_mz
        var_pm1(n,nee) = val_ee
        !
      end do
      !
    end if
    !
    l2n_ec = zero ; l2d_ec = zero
    l2n_mx = zero ; l2d_mx = zero
    l2n_my = zero ; l2d_my = zero
    l2n_mz = zero ; l2d_mz = zero
    l2n_ee = zero ; l2d_ee = zero
    !
    !$ACC LOOP VECTOR REDUCTION(+:l2n_ec,l2d_ec) &
    !$ACC&            REDUCTION(+:l2n_mx,l2d_mx) &
    !$ACC&            REDUCTION(+:l2n_my,l2d_my) &
    !$ACC&            REDUCTION(+:l2n_mz,l2d_mz) &
    !$ACC&            REDUCTION(+:l2n_ee,l2d_ee)
    do n = 1,np
      !
      vd_ec = var_p(n,nec)             ! denominator is: var[P]
      vd_mx = var_p(n,nmx)             ! denominator is: var[P]
      vd_my = var_p(n,nmy)             ! denominator is: var[P]
      vd_mz = var_p(n,nmz)             ! denominator is: var[P]
      vd_ee = var_p(n,nee)             ! denominator is: var[P]
      !
      vn_ec = vd_ec - var_pm1(n,nec) ! numerator   is: var[P] - var[P-1]
      vn_mx = vd_mx - var_pm1(n,nmx) ! numerator   is: var[P] - var[P-1]
      vn_my = vd_my - var_pm1(n,nmy) ! numerator   is: var[P] - var[P-1]
      vn_mz = vd_mz - var_pm1(n,nmz) ! numerator   is: var[P] - var[P-1]
      vn_ee = vd_ee - var_pm1(n,nee) ! numerator   is: var[P] - var[P-1]
      !
      wt = std_elem(this_geom,this_order)%wts(n)
      !
      l2n_ec = l2n_ec + vn_ec*vn_ec*wt
      l2n_mx = l2n_mx + vn_mx*vn_mx*wt
      l2n_my = l2n_my + vn_my*vn_my*wt
      l2n_mz = l2n_mz + vn_mz*vn_mz*wt
      l2n_ee = l2n_ee + vn_ee*vn_ee*wt
      !
      l2d_ec = l2d_ec + vd_ec*vd_ec*wt
      l2d_mx = l2d_mx + vd_mx*vd_mx*wt
      l2d_my = l2d_my + vd_my*vd_my*wt
      l2d_mz = l2d_mz + vd_mz*vd_mz*wt
      l2d_ee = l2d_ee + vd_ee*vd_ee*wt
      !
    end do
    !
    mark_ec = l2n_ec / max(l2d_ec,rndoff)
    mark_mx = l2n_mx / max(l2d_mx,rndoff)
    mark_my = l2n_my / max(l2d_my,rndoff)
    mark_mz = l2n_mz / max(l2d_mz,rndoff)
    mark_ee = l2n_ee / max(l2d_ee,rndoff)
    !
    marker(nc) = max(mark_ec,mark_mx,mark_my,mark_mz,mark_ee)
    !
  end do
  !
#ifdef TIME_MINIAPPS_FOR_GPU_COMPARISON
  call gfr_timer( timings_end(timings_counter) )
#endif
  !
end subroutine compute_resolution_marker_3d
!
!###############################################################################
!
subroutine output_gold_lim_data(gold_marker)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: gold_marker(:)
  !
  !.. Local Scalars ..
  integer :: io,n,ierr
  !
continue
  !
  open (newunit=io,file="lim_gold_marker.dat",status="replace",action="write", &
                   form="unformatted",iostat=ierr)
  write (io) size(gold_marker)
  write (io) (gold_marker(n), n=1,size(gold_marker))
  close (io)
  !
end subroutine output_gold_lim_data
!
!###############################################################################
!
pure &
subroutine compute_resolution_marker_2d(usp,skip_cell,cutoff_order,marker)
  !
  !.. Use Statements ..
  use geovar,         only : n_solpts,ncell,cell
  use filter_mod,     only : filter
  use quadrature_mod, only : std_elem
  !
  !.. Formal Arguments ..
  real(wp),     intent(in) :: usp(:,:)
  logical(ldk), intent(in) :: skip_cell(:)
  integer,      intent(in) :: cutoff_order(:)
  real(wp),  intent(inout) :: marker(:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order,order
  real(wp) :: a,wt
  real(wp) :: val_ec,l2n_ec,l2d_ec,vn_ec,vd_ec,mark_ec
  real(wp) :: val_mx,l2n_mx,l2d_mx,vn_mx,vd_mx,mark_mx
  real(wp) :: val_my,l2n_my,l2d_my,vn_my,vd_my,mark_my
  real(wp) :: val_ee,l2n_ee,l2d_ee,vn_ee,vd_ee,mark_ee
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: var_pm1(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_resolution_marker_2d"
  !
continue
  !
  allocate ( var_pm1(1:n_solpts,1:nq) )
  !$ACC DATA CREATE(var_pm1)
  !
  ! Compute u_hat, the degree-(P-1) solution of density at the sol-pts
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,filter) &
  !$ACC&              PRESENT(usp,skip_cell,cutoff_order) &
  !$ACC&              PRESENT(var_pm1)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index of sol-pts for this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts for this cell
    !
    if (skip_cell(nc)) then
      !
      !$ACC LOOP VECTOR
      do n = 1,np
        var_pm1(n0+n,nec) = usp(n0+n,nec)
        var_pm1(n0+n,nmx) = usp(n0+n,nmx)
        var_pm1(n0+n,nmy) = usp(n0+n,nmy)
        var_pm1(n0+n,nee) = usp(n0+n,nee)
      end do
      !
    else
      !
      this_geom = cell(nc)%geom   ! geometry type of this cell
      this_order = cell(nc)%order ! solution order of this cell
      !
      order = cutoff_order(nc) ! cutoff order for the filter
      !
      !$ACC LOOP VECTOR
      do n = 1,np
        !
        val_ec = zero
        val_mx = zero
        val_my = zero
        val_ee = zero
        !
        do i = filter(this_geom,this_order)%cutoff(order)%col_ptr(n)+1, &
               filter(this_geom,this_order)%cutoff(order)%col_ptr(n+1)
          !
          k = filter(this_geom,this_order)%cutoff(order)%col_idx(i)
          a = filter(this_geom,this_order)%cutoff(order)%col_val(i)
          !
          val_ec = val_ec + a * usp(n0+k,nec)
          val_mx = val_mx + a * usp(n0+k,nmx)
          val_my = val_my + a * usp(n0+k,nmy)
          val_ee = val_ee + a * usp(n0+k,nee)
          !
        end do
        !
        var_pm1(n0+n,nec) = val_ec
        var_pm1(n0+n,nmx) = val_mx
        var_pm1(n0+n,nmy) = val_my
        var_pm1(n0+n,nee) = val_ee
        !
      end do
      !
    end if
    !
  end do
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,std_elem) &
  !$ACC&              PRESENT(usp,marker) &
  !$ACC&              PRESENT(var_pm1)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index of sol-pts for this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts for this cell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    l2n_ec = zero ; l2d_ec = zero
    l2n_mx = zero ; l2d_mx = zero
    l2n_my = zero ; l2d_my = zero
    l2n_ee = zero ; l2d_ee = zero
    !
    !$ACC LOOP VECTOR REDUCTION(+:l2n_ec,l2d_ec) &
    !$ACC&            REDUCTION(+:l2n_mx,l2d_mx) &
    !$ACC&            REDUCTION(+:l2n_my,l2d_my) &
    !$ACC&            REDUCTION(+:l2n_ee,l2d_ee)
    do n = 1,np
      !
      vd_ec = usp(n0+n,nec)             ! denominator is: var[P]
      vd_mx = usp(n0+n,nmx)             ! denominator is: var[P]
      vd_my = usp(n0+n,nmy)             ! denominator is: var[P]
      vd_ee = usp(n0+n,nee)             ! denominator is: var[P]
      !
      vn_ec = vd_ec - var_pm1(n0+n,nec) ! numerator   is: var[P] - var[P-1]
      vn_mx = vd_mx - var_pm1(n0+n,nmx) ! numerator   is: var[P] - var[P-1]
      vn_my = vd_my - var_pm1(n0+n,nmy) ! numerator   is: var[P] - var[P-1]
      vn_ee = vd_ee - var_pm1(n0+n,nee) ! numerator   is: var[P] - var[P-1]
      !
      wt = std_elem(this_geom,this_order)%wts(n)
      !
      l2n_ec = l2n_ec + vn_ec*vn_ec*wt
      l2n_mx = l2n_mx + vn_mx*vn_mx*wt
      l2n_my = l2n_my + vn_my*vn_my*wt
      l2n_ee = l2n_ee + vn_ee*vn_ee*wt
      !
      l2d_ec = l2d_ec + vd_ec*vd_ec*wt
      l2d_mx = l2d_mx + vd_mx*vd_mx*wt
      l2d_my = l2d_my + vd_my*vd_my*wt
      l2d_ee = l2d_ee + vd_ee*vd_ee*wt
      !
    end do
    !
    mark_ec = l2n_ec / max(l2d_ec,rndoff)
    mark_mx = l2n_mx / max(l2d_mx,rndoff)
    mark_my = l2n_my / max(l2d_my,rndoff)
    mark_ee = l2n_ee / max(l2d_ee,rndoff)
    !
    marker(nc) = max(mark_ec,mark_mx,mark_my,mark_ee)
    !
  end do
  !
  !$ACC END DATA
  deallocate ( var_pm1 )
  !
end subroutine compute_resolution_marker_2d
!
!###############################################################################
!
pure &
subroutine compute_resolution_marker_1d(usp,skip_cell,cutoff_order,marker)
  !
  !.. Use Statements ..
  use geovar,         only : n_solpts,ncell,cell
  use filter_mod,     only : filter
  use quadrature_mod, only : std_elem
  !
  !.. Formal Arguments ..
  real(wp),     intent(in) :: usp(:,:)
  logical(ldk), intent(in) :: skip_cell(:)
  integer,      intent(in) :: cutoff_order(:)
  real(wp),  intent(inout) :: marker(:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order,order
  real(wp) :: a,wt
  real(wp) :: val_ec,l2n_ec,l2d_ec,vn_ec,vd_ec,mark_ec
  real(wp) :: val_mx,l2n_mx,l2d_mx,vn_mx,vd_mx,mark_mx
  real(wp) :: val_ee,l2n_ee,l2d_ee,vn_ee,vd_ee,mark_ee
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: var_pm1(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_resolution_marker_2d"
  !
continue
  !
  allocate ( var_pm1(1:n_solpts,1:nq) )
  !$ACC DATA CREATE(var_pm1)
  !
  ! Compute u_hat, the degree-(P-1) solution of density at the sol-pts
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,filter) &
  !$ACC&              PRESENT(usp,skip_cell,cutoff_order) &
  !$ACC&              PRESENT(var_pm1)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index of sol-pts for this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts for this cell
    !
    if (skip_cell(nc)) then
      !
      !$ACC LOOP VECTOR
      do n = 1,np
        var_pm1(n0+n,nec) = usp(n0+n,nec)
        var_pm1(n0+n,nmx) = usp(n0+n,nmx)
        var_pm1(n0+n,nee) = usp(n0+n,nee)
      end do
      !
    else
      !
      this_geom = cell(nc)%geom   ! geometry type of this cell
      this_order = cell(nc)%order ! solution order of this cell
      !
      order = cutoff_order(nc) ! cutoff order for the filter
      !
      !$ACC LOOP VECTOR
      do n = 1,np
        !
        val_ec = zero
        val_mx = zero
        val_ee = zero
        !
        do i = filter(this_geom,this_order)%cutoff(order)%col_ptr(n)+1, &
               filter(this_geom,this_order)%cutoff(order)%col_ptr(n+1)
          !
          k = filter(this_geom,this_order)%cutoff(order)%col_idx(i)
          a = filter(this_geom,this_order)%cutoff(order)%col_val(i)
          !
          val_ec = val_ec + a * usp(n0+k,nec)
          val_mx = val_mx + a * usp(n0+k,nmx)
          val_ee = val_ee + a * usp(n0+k,nee)
          !
        end do
        !
        var_pm1(n0+n,nec) = val_ec
        var_pm1(n0+n,nmx) = val_mx
        var_pm1(n0+n,nee) = val_ee
        !
      end do
      !
    end if
    !
  end do
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,std_elem) &
  !$ACC&              PRESENT(usp,marker) &
  !$ACC&              PRESENT(var_pm1)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index of sol-pts for this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts for this cell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    l2n_ec = zero ; l2d_ec = zero
    l2n_mx = zero ; l2d_mx = zero
    l2n_ee = zero ; l2d_ee = zero
    !
    !$ACC LOOP VECTOR REDUCTION(+:l2n_ec,l2d_ec) &
    !$ACC&            REDUCTION(+:l2n_mx,l2d_mx) &
    !$ACC&            REDUCTION(+:l2n_ee,l2d_ee)
    do n = 1,np
      !
      vd_ec = usp(n0+n,nec)             ! denominator is: var[P]
      vd_mx = usp(n0+n,nmx)             ! denominator is: var[P]
      vd_ee = usp(n0+n,nee)             ! denominator is: var[P]
      !
      vn_ec = vd_ec - var_pm1(n0+n,nec) ! numerator   is: var[P] - var[P-1]
      vn_mx = vd_mx - var_pm1(n0+n,nmx) ! numerator   is: var[P] - var[P-1]
      vn_ee = vd_ee - var_pm1(n0+n,nee) ! numerator   is: var[P] - var[P-1]
      !
      wt = std_elem(this_geom,this_order)%wts(n)
      !
      l2n_ec = l2n_ec + vn_ec*vn_ec*wt
      l2n_mx = l2n_mx + vn_mx*vn_mx*wt
      l2n_ee = l2n_ee + vn_ee*vn_ee*wt
      !
      l2d_ec = l2d_ec + vd_ec*vd_ec*wt
      l2d_mx = l2d_mx + vd_mx*vd_mx*wt
      l2d_ee = l2d_ee + vd_ee*vd_ee*wt
      !
    end do
    !
    mark_ec = l2n_ec / max(l2d_ec,rndoff)
    mark_mx = l2n_mx / max(l2d_mx,rndoff)
    mark_ee = l2n_ee / max(l2d_ee,rndoff)
    !
    marker(nc) = max(mark_ec,mark_mx,mark_ee)
    !
  end do
  !
  !$ACC END DATA
  deallocate ( var_pm1 )
  !
end subroutine compute_resolution_marker_1d
!
!###############################################################################
!
pure &
subroutine compute_resolution_marker_nd(mbeg,mend,usp,skip_cell, &
                                        cutoff_order,marker)
  !
  !.. Use Statements ..
  use order_mod,      only : maxSP
  use geovar,         only : n_solpts,ncell,cell
  use filter_mod,     only : filter
  use quadrature_mod, only : std_elem
  !
  !.. Formal Arguments ..
  integer,      intent(in) :: mbeg,mend
  real(wp),     intent(in) :: usp(:,:)
  logical(ldk), intent(in) :: skip_cell(:)
  integer,      intent(in) :: cutoff_order(:)
  real(wp),  intent(inout) :: marker(:)
  !
  !.. Local Scalars ..
  integer  :: i,k,m,n,nc,np,n0
  integer  :: this_geom,this_order,order
  real(wp) :: a,wt,mark
  real(wp) :: val,l2_n,l2_d,vd,vn
  !
  !.. Local Arrays ..
  real(wp) :: val_n(1:maxSP,mbeg:mend)
  real(wp) :: val_d(1:maxSP,mbeg:mend)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: var_pm1(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_resolution_marker_nd"
  !
continue
  !
  allocate ( var_pm1(1:n_solpts,mbeg:mend) )
  !$ACC DATA CREATE(var_pm1)
  !
  ! Compute u_hat, the degree-(P-1) solution of density at the sol-pts
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,filter) &
  !$ACC&              PRESENT(usp,skip_cell,cutoff_order) &
  !$ACC&              PRESENT(var_pm1)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index of sol-pts for this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts for this cell
    !
    if (skip_cell(nc)) then
      !
      !$ACC LOOP VECTOR COLLAPSE(2)
      do m = mbeg,mend
        do n = 1,np
          var_pm1(n0+n,m) = usp(n0+n,m)
        end do
      end do
      !
    else
      !
      this_geom = cell(nc)%geom   ! geometry type of this cell
      this_order = cell(nc)%order ! solution order of this cell
      !
      order = cutoff_order(nc) ! cutoff order for the filter
      !
      !$ACC LOOP VECTOR COLLAPSE(2)
      do n = 1,np
        do m = mbeg,mend
          !
          val = zero
          !
          do i = filter(this_geom,this_order)%cutoff(order)%col_ptr(n)+1, &
                 filter(this_geom,this_order)%cutoff(order)%col_ptr(n+1)
            !
            k = filter(this_geom,this_order)%cutoff(order)%col_idx(i)
            a = filter(this_geom,this_order)%cutoff(order)%col_val(i)
            !
            val = val + a * usp(n0+k,m)
            !
          end do
          !
          var_pm1(n0+n,m) = val
          !
        end do
      end do
      !
    end if
    !
  end do
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,std_elem) &
  !$ACC&              PRESENT(usp,marker) &
  !$ACC&              PRESENT(var_pm1) &
  !$ACC&              PRIVATE(val_n,val_d)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index of sol-pts for this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts for this cell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do m = mbeg,mend
      do n = 1,np
        !
        val = usp(n0+n,m)
        !
        val_d(n,m) = val                   ! denominator is: var[P]
        val_n(n,m) = val - var_pm1(n0+n,m) ! numerator   is: var[P] - var[P-1]
        !
      end do
    end do
    !
    mark = -huge(zero)
    !
    do m = mbeg,mend
      !
      l2_n = zero
      l2_d = zero
      !
      !$ACC LOOP VECTOR REDUCTION(+:l2_n,l2_d)
      do n = 1,np
        !
        vn = val_n(n,m)
        vd = val_d(n,m)
        !
        wt = std_elem(this_geom,this_order)%wts(n)
        !
        l2_n = l2_n + vn*vn*wt
        l2_d = l2_d + vd*vd*wt
        !
      end do
      !
      mark = max( mark , l2_n/max(l2_d,rndoff) )
      !
    end do
    !
    marker(nc) = max(marker(nc),mark)
    !
  end do
  !
  !
  !$ACC END DATA
  deallocate ( var_pm1 )
  !
end subroutine compute_resolution_marker_nd
!
!###############################################################################
!
pure &
subroutine compute_resolution_marker_1var(var,skip_cell,cutoff_order,marker)
  !
  !.. Use Statements ..
  use order_mod,      only : maxSP
  use geovar,         only : n_solpts,ncell,cell
  use filter_mod,     only : filter
  use quadrature_mod, only : std_elem
  !
  !.. Formal Arguments ..
  real(wp),     intent(in) :: var(:)
  logical(ldk), intent(in) :: skip_cell(:)
  integer,      intent(in) :: cutoff_order(:)
  real(wp),  intent(inout) :: marker(:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order,order
  real(wp) :: l2_n,l2_d,var_d,var_n
  real(wp) :: a,val
  !
  !.. Local Allocatable Arrays ..
  real(wp) :: var_p(1:maxSP)
  real(wp) :: var_pm1(1:maxSP)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_resolution_marker_1var"
  !
continue
  !
  ! Compute u_hat, the degree-(P-1) solution of density at the sol-pts
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,filter,std_elem) &
  !$ACC&              PRESENT(var,marker,skip_cell,cutoff_order) &
  !$ACC&              PRIVATE(var_p,var_pm1)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index of sol-pts for this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts for this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      var_p(n) = var(n0+n)
    end do
    !
    if (skip_cell(nc)) then
      !
      !$ACC LOOP VECTOR
      do n = 1,np
        var_pm1(n) = var_p(n)
      end do
      !
    else
      !
      order = cutoff_order(nc) ! cutoff order for the filter
      !
      !$ACC LOOP VECTOR
      do n = 1,np
        !
        val = zero
        !
        do i = filter(this_geom,this_order)%cutoff(order)%col_ptr(n)+1, &
               filter(this_geom,this_order)%cutoff(order)%col_ptr(n+1)
          !
          k = filter(this_geom,this_order)%cutoff(order)%col_idx(i)
          a = filter(this_geom,this_order)%cutoff(order)%col_val(i)
          !
          val = val + a * var_p(k)
          !
        end do
        !
        var_pm1(n) = val
        !
      end do
      !
    end if
    !
    l2_n = zero ! numerator   for resolution indicator
    l2_d = zero ! denominator for resolution indicator
    !
    !$ACC LOOP VECTOR REDUCTION(+:l2_n,l2_d)
    do n = 1,np
      !
      var_d = var_p(n)           ! denominator var is: var[P]
      var_n = var_d - var_pm1(n) ! numerator   var is: var[P] - var[P-1]
      !
      l2_n = l2_n + var_n*var_n*std_elem(this_geom,this_order)%wts(n)
      l2_d = l2_d + var_d*var_d*std_elem(this_geom,this_order)%wts(n)
      !
    end do
    !
    marker(nc) = l2_n / max(l2_d,rndoff)
    !
  end do
  !
end subroutine compute_resolution_marker_1var
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
pure &
subroutine compute_scaled_marker(skip_cell,marker,scaled_marker)
  !
  !.. Use Statements ..
  use geovar, only : ncell,cell
  use ovar,   only : dpsi => marker_dpsi
  use ovar,   only : c1 => marker_c1
  use ovar,   only : c2 => marker_c2
  !
  !.. Formal Arguments ..
  logical(ldk), intent(in) :: skip_cell(:)
  real(wp),  intent(inout) :: marker(:)
  real(wp),  intent(inout) :: scaled_marker(:)
  !
  !.. Local Scalars ..
  integer  :: nc
  real(wp) :: l10_marker,psi0
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_scaled_marker"
  !
continue
  !
  ! Compute the log10 of marker and scaled_marker. For any skipped cells,
  ! set marker and scaled_marker to values that will prevent them from
  ! being marked as troubled cells.
  !
  !$ACC PARALLEL LOOP PRESENT(skip_cell,marker,scaled_marker)
  do nc = 1,ncell
    !
    if (skip_cell(nc)) then
      !
      marker(nc) = -two*ten
      scaled_marker(nc) = zero
      !
    else
      !
      l10_marker = log10( max( marker(nc) , rndoff ) )
      !
      psi0 = -one*(c1 + c2*log10(real(cell(nc)%order,kind=wp)))
      !
      if (l10_marker <= psi0-dpsi) then
        scaled_marker(nc) = zero
      else if (l10_marker >= psi0+dpsi) then
        scaled_marker(nc) = one
      else
        scaled_marker(nc) = half*(one + sin(half*pi*(l10_marker-psi0)/dpsi))
      end if
      !
      marker(nc) = l10_marker
      !
    end if
    !
  end do
  !
end subroutine compute_scaled_marker
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
pure &
subroutine res_ind_apply_proj_3d(usp,skip_cell,lmtr_order)
  !
  !.. Use Statements ..
  use order_mod,      only : maxSP
  use geovar,         only : ncell,cell
  use projection_mod, only : projct
  !
  !.. Formal Arguments ..
  logical(ldk), intent(in) :: skip_cell(:)
  integer,      intent(in) :: lmtr_order(:)
  real(wp),  intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order,order
  real(wp) :: a,plim_ec,plim_mx,plim_my,plim_mz,plim_ee
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "res_ind_apply_proj_3d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,skip_cell,lmtr_order,projct) &
  !$ACC&              PRIVATE(old_usp)
  limit_cell_loop: do nc = 1,ncell
    !
    if (skip_cell(nc)) then
      cycle limit_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index of sol-pts for this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts for this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    order = lmtr_order(nc) ! order to limit cell solution to
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      old_usp(n,nec) = usp(n0+n,nec)
      old_usp(n,nmx) = usp(n0+n,nmx)
      old_usp(n,nmy) = usp(n0+n,nmy)
      old_usp(n,nmz) = usp(n0+n,nmz)
      old_usp(n,nee) = usp(n0+n,nee)
    end do
    !
    ! Apply the projection limiter by multplying the projection
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      plim_ec = zero
      plim_mx = zero
      plim_my = zero
      plim_mz = zero
      plim_ee = zero
      !
      do i = projct(this_geom,this_order)%limiter(order)%row_ptr(n)+1, &
             projct(this_geom,this_order)%limiter(order)%row_ptr(n+1)
        !
        k = projct(this_geom,this_order)%limiter(order)%row_idx(i)
        a = projct(this_geom,this_order)%limiter(order)%row_val(i)
        !
        plim_ec = plim_ec + a * old_usp(k,nec)
        plim_mx = plim_mx + a * old_usp(k,nmx)
        plim_my = plim_my + a * old_usp(k,nmy)
        plim_mz = plim_mz + a * old_usp(k,nmz)
        plim_ee = plim_ee + a * old_usp(k,nee)
        !
      end do
      !
      usp(n0+n,nec) = plim_ec
      usp(n0+n,nmx) = plim_mx
      usp(n0+n,nmy) = plim_my
      usp(n0+n,nmz) = plim_mz
      usp(n0+n,nee) = plim_ee
      !
    end do
    !
  end do limit_cell_loop
  !
end subroutine res_ind_apply_proj_3d
!
!###############################################################################
!
pure &
subroutine res_ind_apply_proj_2d(usp,skip_cell,lmtr_order)
  !
  !.. Use Statements ..
  use order_mod,      only : maxSP
  use geovar,         only : ncell,cell
  use projection_mod, only : projct
  !
  !.. Formal Arguments ..
  logical(ldk), intent(in) :: skip_cell(:)
  integer,      intent(in) :: lmtr_order(:)
  real(wp),  intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order,order
  real(wp) :: a,plim_ec,plim_mx,plim_my,plim_ee
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "res_ind_apply_proj_2d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,skip_cell,lmtr_order,projct) &
  !$ACC&              PRIVATE(old_usp)
  limit_cell_loop: do nc = 1,ncell
    !
    if (skip_cell(nc)) then
      cycle limit_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index of sol-pts for this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts for this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    order = lmtr_order(nc) ! order to limit cell solution to
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      old_usp(n,nec) = usp(n0+n,nec)
      old_usp(n,nmx) = usp(n0+n,nmx)
      old_usp(n,nmy) = usp(n0+n,nmy)
      old_usp(n,nee) = usp(n0+n,nee)
    end do
    !
    ! Apply the projection limiter by multplying the projection
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      plim_ec = zero
      plim_mx = zero
      plim_my = zero
      plim_ee = zero
      !
      do i = projct(this_geom,this_order)%limiter(order)%row_ptr(n)+1, &
             projct(this_geom,this_order)%limiter(order)%row_ptr(n+1)
        !
        k = projct(this_geom,this_order)%limiter(order)%row_idx(i)
        a = projct(this_geom,this_order)%limiter(order)%row_val(i)
        !
        plim_ec = plim_ec + a * old_usp(k,nec)
        plim_mx = plim_mx + a * old_usp(k,nmx)
        plim_my = plim_my + a * old_usp(k,nmy)
        plim_ee = plim_ee + a * old_usp(k,nee)
        !
      end do
      !
      usp(n0+n,nec) = plim_ec
      usp(n0+n,nmx) = plim_mx
      usp(n0+n,nmy) = plim_my
      usp(n0+n,nee) = plim_ee
      !
    end do
    !
  end do limit_cell_loop
  !
end subroutine res_ind_apply_proj_2d
!
!###############################################################################
!
pure &
subroutine res_ind_apply_proj_1d(usp,skip_cell,lmtr_order)
  !
  !.. Use Statements ..
  use order_mod,      only : maxSP
  use geovar,         only : ncell,cell
  use projection_mod, only : projct
  !
  !.. Formal Arguments ..
  logical(ldk), intent(in) :: skip_cell(:)
  integer,      intent(in) :: lmtr_order(:)
  real(wp),  intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order,order
  real(wp) :: a,plim_ec,plim_mx,plim_ee
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "res_ind_apply_proj_1d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,skip_cell,lmtr_order,projct) &
  !$ACC&              PRIVATE(old_usp)
  limit_cell_loop: do nc = 1,ncell
    !
    if (skip_cell(nc)) then
      cycle limit_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index of sol-pts for this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts for this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    order = lmtr_order(nc) ! order to limit cell solution to
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      old_usp(n,nec) = usp(n0+n,nec)
      old_usp(n,nmx) = usp(n0+n,nmx)
      old_usp(n,nee) = usp(n0+n,nee)
    end do
    !
    ! Apply the projection limiter by multplying the projection
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      plim_ec = zero
      plim_mx = zero
      plim_ee = zero
      !
      do i = projct(this_geom,this_order)%limiter(order)%row_ptr(n)+1, &
             projct(this_geom,this_order)%limiter(order)%row_ptr(n+1)
        !
        k = projct(this_geom,this_order)%limiter(order)%row_idx(i)
        a = projct(this_geom,this_order)%limiter(order)%row_val(i)
        !
        plim_ec = plim_ec + a * old_usp(k,nec)
        plim_mx = plim_mx + a * old_usp(k,nmx)
        plim_ee = plim_ee + a * old_usp(k,nee)
        !
      end do
      !
      usp(n0+n,nec) = plim_ec
      usp(n0+n,nmx) = plim_mx
      usp(n0+n,nee) = plim_ee
      !
    end do
    !
  end do limit_cell_loop
  !
end subroutine res_ind_apply_proj_1d
!
!###############################################################################
!
pure &
subroutine res_ind_apply_proj_nd(mbeg,mend,usp,skip_cell,lmtr_order)
  !
  !.. Use Statements ..
  use order_mod,      only : maxSP
  use geovar,         only : ncell,cell
  use projection_mod, only : projct
  !
  !.. Formal Arguments ..
  integer,      intent(in) :: mbeg,mend
  logical(ldk), intent(in) :: skip_cell(:)
  integer,      intent(in) :: lmtr_order(:)
  real(wp),  intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,m,n,nc,np,n0
  integer  :: this_geom,this_order,order
  real(wp) :: a,val
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "res_ind_apply_proj_nd"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,skip_cell,lmtr_order,projct) &
  !$ACC&              PRIVATE(old_usp)
  limit_cell_loop: do nc = 1,ncell
    !
    if (skip_cell(nc)) then
      cycle limit_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index of sol-pts for this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts for this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    order = lmtr_order(nc) ! order to limit cell solution to
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do m = mbeg,mend
      do n = 1,np
        old_usp(n,m) = usp(n0+n,m)
      end do
    end do
    !
    ! Apply the projection limiter by multplying the projection
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do m = mbeg,mend
      do n = 1,np
        !
        val = zero
        !
        do i = projct(this_geom,this_order)%limiter(order)%row_ptr(n)+1, &
               projct(this_geom,this_order)%limiter(order)%row_ptr(n+1)
          !
          k = projct(this_geom,this_order)%limiter(order)%row_idx(i)
          a = projct(this_geom,this_order)%limiter(order)%row_val(i)
          !
          val = val + a * old_usp(k,m)
          !
        end do
        !
        usp(n0+n,m) = val
        !
      end do
    end do
    !
  end do limit_cell_loop
  !
end subroutine res_ind_apply_proj_nd
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine jump_indicator(troubled_cell)
  !
  !.. Use Statements ..
  use geovar,         only : nr,n_flxpts,ncell,cell,nface,face
  use quadrature_mod, only : std_elem
  use flowvar,        only : facegpsp
  use metrics_mod,    only : facenrm
  !
  !.. Optional Arguments ..
  logical(ldk), optional, intent(inout) :: troubled_cell(:)
  !
  !.. Local Scalars ..
  integer  :: nc,ncf,nf,k,nfp,f0,ierr
  integer  :: num_cell_faces
  integer  :: this_geom,face_geom,face_order
  real(wp) :: ds,rnx,rny,rnz,val,psi0
  real(wp) :: pres_l,pres_r,p_jump
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: pres_jump(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "jump_indicator"
  real(wp), parameter :: dpsi = half
  real(wp), parameter :: c1 = 2.25_wp
  real(wp), parameter :: c2 = three
  !
continue
  !
  if (.not. allocated(marker)) then
    allocate ( marker(1:ncell) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"marker",1,__LINE__,__FILE__,ierr,error_message)
    !$ACC UPDATE DEVICE(marker)
  end if
  !
  if (.not. allocated(scaled_marker)) then
    allocate ( scaled_marker(1:ncell) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"scaled_marker",1,__LINE__,__FILE__, &
                     ierr,error_message)
    !$ACC UPDATE DEVICE(scaled_marker)
  end if
  !
  allocate ( pres_jump(1:nface) )
  !$ACC DATA CREATE(pres_jump)
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP GANG &
    !$ACC&              VECTOR_LENGTH(32) &
    !$ACC&              PRESENT(face,facegpsp,pres_jump,std_elem)
    do nf = 1,nface
      !
      face_geom = face(nf)%geom  ! geometry of face
      face_order = face(nf)%order ! order of face
      !
      f0 = face(nf)%beg_fp - 1
      nfp = face(nf)%end_fp - f0
      !
      val = zero
      !
      !$ACC LOOP VECTOR REDUCTION(+:val)
      do k = 1,nfp
        !
        rnx = facenrm(1,f0+k)
        rny = facenrm(2,f0+k)
        rnz = facenrm(3,f0+k)
        ds = sqrt(rnx*rnx + rny*rny + rnz*rnz)
        !
        pres_l = facegpsp(f0+k,npres,1)
        pres_r = facegpsp(f0+k,npres,2)
        !
        p_jump = two * abs(pres_l - pres_r) / (pres_l + pres_r)
        !
        val = val + ds*p_jump*std_elem(face_geom,face_order)%wts(k)
        !
      end do
      !
      pres_jump(nf) = val
      !
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP GANG &
    !$ACC&              VECTOR_LENGTH(32) &
    !$ACC&              PRESENT(face,facegpsp,pres_jump,std_elem)
    do nf = 1,nface
      !
      face_geom = face(nf)%geom  ! geometry of face
      face_order = face(nf)%order ! order of face
      !
      f0 = face(nf)%beg_fp - 1
      nfp = face(nf)%end_fp - f0
      !
      val = zero
      !
      !$ACC LOOP VECTOR REDUCTION(+:val)
      do k = 1,nfp
        !
        rnx = facenrm(1,f0+k)
        rny = facenrm(2,f0+k)
        ds = sqrt(rnx*rnx + rny*rny)
        !
        pres_l = facegpsp(f0+k,npres,1)
        pres_r = facegpsp(f0+k,npres,2)
        !
        p_jump = two * abs(pres_l - pres_r) / (pres_l + pres_r)
        !
        val = val + ds*p_jump*std_elem(face_geom,face_order)%wts(k)
        !
      end do
      !
      pres_jump(nf) = val
      !
    end do
  else
    !$ACC PARALLEL LOOP GANG &
    !$ACC&              VECTOR_LENGTH(32) &
    !$ACC&              PRESENT(face,facegpsp,pres_jump,std_elem)
    do nf = 1,nface
      !
      face_geom = face(nf)%geom  ! geometry of face
      face_order = face(nf)%order ! order of face
      !
      f0 = face(nf)%beg_fp - 1
      nfp = face(nf)%end_fp - f0
      !
      val = zero
      !
      !$ACC LOOP VECTOR REDUCTION(+:val)
      do k = 1,nfp
        !
        ds = abs(facenrm(1,f0+k))
        !
        pres_l = facegpsp(f0+k,npres,1)
        pres_r = facegpsp(f0+k,npres,2)
        !
        p_jump = two * abs(pres_l - pres_r) / (pres_l + pres_r)
        !
        val = val + ds*p_jump*std_elem(face_geom,face_order)%wts(k)
        !
      end do
      !
      pres_jump(nf) = val
      !
    end do
  end if
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR &
  !$ACC&              PRESENT(marker,cell,pres_jump)
  do nc = 1,ncell
    !
    this_geom = cell(nc)%geom
    !
    num_cell_faces = size(cell(nc)%face)
    !
    !$ACC LOOP SEQ REDUCTION(+:val)
    do ncf = 1,num_cell_faces
      !
      nf = cell(nc)%face(ncf)%idx
      !
      val = val + cl(this_geom)*pres_jump(nf)
      !
    end do
    !
    marker(nc) = val
    !
  end do
  !
  !$ACC END DATA
  deallocate ( pres_jump )
  !
  !$ACC PARALLEL LOOP PRESENT(marker)
  do nc = 1,ncell
    val = max( marker(nc) , rndoff )
    marker(nc) = log10( val )
  end do
  !
  !$ACC PARALLEL LOOP PRESENT(cell,marker,scaled_marker)
  do nc = 1,ncell
    !
    psi0 = -one*(c1 + c2*log10(real(cell(nc)%order,kind=wp)))
    !
    if (marker(nc) <= psi0-dpsi) then
      scaled_marker(nc) = zero
    else if (marker(nc) >= psi0+dpsi) then
      scaled_marker(nc) = one
    else
      scaled_marker(nc) = half*(one + sin(half*pi*(marker(nc)-psi0)/dpsi))
    end if
    !
  end do
  !
  if (present(troubled_cell)) then
    !$ACC PARALLEL LOOP PRESENT(scaled_marker,troubled_cell)
    do nc = 1,ncell
      troubled_cell(nc) = (scaled_marker(nc) >= 1.0e-8_wp)
    end do
  end if
  !
end subroutine jump_indicator
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine limiters_memory_usage(iunit)
  !
  !.. Use Statements ..
  use module_memory
  !
  !.. Formal Arguments ..
  integer, optional, intent(in) :: iunit
  !
  !.. Local Scalars ..
  integer :: iou
  character(len=36) :: array_name
  type(memory) :: array
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
  if (allocated(cell_ave_sol)) then
    call array%set( storage_size(cell_ave_sol,kind=inttype) , &
                            size(cell_ave_sol,kind=inttype) )
    call total%add(array)
    write (array_name,5) "cell_ave_sol"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "cell_ave_sol"
  end if
  if (allocated(umodes)) then
    call array%set( storage_size(umodes,kind=inttype) , &
                            size(umodes,kind=inttype) )
    call total%add(array)
    write (array_name,5) "umodes"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "umodes"
  end if
  if (allocated(troubled_cell)) then
    call array%set( storage_size(troubled_cell,kind=inttype) , &
                            size(troubled_cell,kind=inttype) )
    call total%add(array)
    write (array_name,5) "troubled_cell"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "troubled_cell"
  end if
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module module_limiters")
  2 format (" Array: ",a7,dt)
  3 format (" Total memory = ",f8.2,1x,a)
  5 format ("'",a,"'")
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  !
#endif
end subroutine limiters_memory_usage
!
!###############################################################################
!
#include "Functions/debug_timer.f90"
!
!###############################################################################
!
end module module_limiters
