module module_limiters
  !
  !.. Use Statements ..
  use module_kind_types, global_debug_timer => debug_timer
  !
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,ntk,ntl,nq
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
contains
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
  use order_mod,  only : maxpts,n_max_order
  use geovar,     only : ncell,cell
  use flowvar,    only : usp
  use projection_mod, only : projct
  use filter_mod, only : filter
  use ovar, only : Limiter_Filter_Option
  !
  !.. Local Scalars ..
  integer :: k,m,n1,n2,nc,np,ntc
  integer :: this_geom,this_order
  integer :: dropped_orders,proj_order
  logical(lk) :: there_are_no_troubles
  logical(lk) :: first_output
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxpts,1:nq) :: tusp
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
      endif
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
  if (nr == 3) then
    call apply_filter_limiter_to_usp_3d(troubled_cell,usp)
  else if (nr == 2) then
    call apply_filter_limiter_to_usp_2d(troubled_cell,usp)
  else
    call apply_filter_limiter_to_usp_nd(troubled_cell,usp)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine apply_filter_limiter_to_usp
!
!###############################################################################
!
pure subroutine apply_filter_limiter_to_usp_3d(troubled_cell,usp)
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
    troubled_cell(nc) = .false.
    !
  end do filt_cell_loop
  !
end subroutine apply_filter_limiter_to_usp_3d
!
!###############################################################################
!
pure subroutine apply_filter_limiter_to_usp_2d(troubled_cell,usp)
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
    troubled_cell(nc) = .false.
    !
  end do filt_cell_loop
  !
end subroutine apply_filter_limiter_to_usp_2d
!
!###############################################################################
!
pure subroutine apply_filter_limiter_to_usp_nd(troubled_cell,usp)
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
    do m = 1,nq
      do n = 1,np
        old_usp(n,m) = usp(n0+n,m)
      end do
    end do
    !
    ! Apply the filter limiter by multplying the filter
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do m = 1,nq
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
  if (nr == 3) then
    call apply_projct_limiter_to_usp_3d(troubled_cell,usp,dropped_orders)
  else if (nr == 2) then
    call apply_projct_limiter_to_usp_2d(troubled_cell,usp,dropped_orders)
  else
    call apply_projct_limiter_to_usp_nd(troubled_cell,usp,dropped_orders)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine apply_projct_limiter_to_usp
!
!###############################################################################
!
pure subroutine apply_projct_limiter_to_usp_3d(troubled_cell,usp,dropped_orders)
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
    if (proj_order == 0) then
      troubled_cell(nc) = .false.
    end if
    !
  end do proj_cell_loop
  !
end subroutine apply_projct_limiter_to_usp_3d
!
!###############################################################################
!
pure subroutine apply_projct_limiter_to_usp_2d(troubled_cell,usp,dropped_orders)
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
    if (proj_order == 0) then
      troubled_cell(nc) = .false.
    end if
    !
  end do proj_cell_loop
  !
end subroutine apply_projct_limiter_to_usp_2d
!
!###############################################################################
!
pure subroutine apply_projct_limiter_to_usp_nd(troubled_cell,usp,dropped_orders)
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
    do m = 1,nq
      do n = 1,np
        old_usp(n,m) = usp(n0+n,m)
      end do
    end do
    !
    ! Apply the projection limiter by multplying the projection
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do m = 1,nq
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
  integer :: ord_mlt,ord_add,actual_order
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
      troubled_cell(nc) = .false.
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
    if (actual_order-proj_order <= 1) then
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
subroutine resolution_indicator(apply_projection,troubled_cell,orders_to_cutoff)
  !
  !.. Use Statements ..
  use order_mod, only : p_order,n_order
  use geovar,    only : ncell,cell
  use flowvar,   only : usp
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
  !.. Local Arrays ..
  integer :: nlim(0:n_order)
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
  ! Use the resolution indicator to compute marker and scaled_marker
  !
  call compute_resolution_marker(usp,skip_cell,lmtr_order,marker,scaled_marker)
  !
  ! Update the troubled_cell indicator if it is present
  !
  if (present(troubled_cell)) then
    !$ACC PARALLEL LOOP PRESENT(scaled_marker,troubled_cell)
    do nc = 1,ncell
      troubled_cell(nc) = (scaled_marker(nc) >= 1.0e-8_wp)
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
          lmtr_order(nc) = max(1,min(this_order,limit_order))
        end if
        !
        if (lmtr_order(nc) == this_order) then
          skip_cell(nc) = .true.
        end if
        !
      end do
      !
#ifdef DEBUG_ON
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
      if (ncpu > 1) then
        call mpi_allreduce(MPI_IN_PLACE,nlim,size(nlim,kind=int_mpi), &
                           mpi_inttyp,MPI_SUM,MPI_COMM_WORLD,mpierr)
      end if
      if (mypnum == glb_root) then
        do n = lbound(nlim,dim=1),ubound(nlim,dim=1)
          if (nlim(n) > 0) write (iout,1) nlim(n),n
        end do
      end if
#endif
      !
      ! Apply the projection limiter to the cells marked as needing limiting
      !
      nr = nme-nmb+1
      !
      if (nr == 3) then
        call res_ind_apply_proj_3d(usp,skip_cell,lmtr_order)
      else if (nr == 2) then
        call res_ind_apply_proj_2d(usp,skip_cell,lmtr_order)
      else
        call res_ind_apply_proj_nd(usp,skip_cell,lmtr_order)
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
  ! Format Statements
  !
  1 format (8x,i8," Cells limited to P",i0)
  !
end subroutine resolution_indicator
!
!###############################################################################
!
pure subroutine compute_resolution_marker(usp,skip_cell,cutoff_order, &
                                          marker,scaled_marker)
  !
  !.. Use Statements ..
  use order_mod,   only : maxSP
  use geovar,      only : n_solpts,ncell,cell
  use metrics_mod, only : wtdjac
  use filter_mod,  only : filter
  !
  !.. Formal Arguments ..
  real(wp),     intent(in) :: usp(:,:)
  logical(ldk), intent(in) :: skip_cell(:)
  integer,      intent(in) :: cutoff_order(:)
  real(wp),  intent(inout) :: marker(:)
  real(wp),  intent(inout) :: scaled_marker(:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0,n1,n2
  integer  :: this_geom,this_order,order
  real(wp) :: l2_n,l2_d,rho_d,rho_n
  real(wp) :: a,val,l10_marker,psi0
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: rho_pm1(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_resolution_marker"
  real(wp),         parameter :: dpsi = half
  real(wp),         parameter :: c1 = four
  real(wp),         parameter :: c2 = 4.25_wp
  !
continue
  !
  allocate ( rho_pm1(1:n_solpts) )
  !$ACC DATA CREATE(rho_pm1)
  !
  ! Compute u_hat, the degree-(P-1) solution of density at the sol-pts
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,filter) &
  !$ACC&              PRESENT(usp,skip_cell,cutoff_order) &
  !$ACC&              PRESENT(rho_pm1)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index of sol-pts for this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts for this cell
    !
    if (skip_cell(nc)) then
      !
      !$ACC LOOP VECTOR
      do n = 1,np
        rho_pm1(n0+n) = usp(n0+n,nec)
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
        val = zero
        !
        do i = filter(this_geom,this_order)%cutoff(order)%col_ptr(n)+1, &
               filter(this_geom,this_order)%cutoff(order)%col_ptr(n+1)
          !
          k = filter(this_geom,this_order)%cutoff(order)%col_idx(i)
          a = filter(this_geom,this_order)%cutoff(order)%col_val(i)
          !
          val = val + a * usp(n0+k,nec)
          !
        end do
        !
        rho_pm1(n0+n) = val
        !
      end do
      !
    end if
    !
  end do
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,wtdjac) &
  !$ACC&              PRESENT(usp,marker) &
  !$ACC&              PRESENT(rho_pm1)
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
    n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
    !
    l2_n = zero ! numerator   for resolution indicator
    l2_d = zero ! denominator for resolution indicator
    !
    !$ACC LOOP VECTOR REDUCTION(+:l2_n,l2_d)
    do n = n1,n2
      !
      rho_d = usp(n,nec)         ! denominator rho is: rho[P]
      rho_n = rho_d - rho_pm1(n) ! numerator   rho is: rho[P] - rho[P-1]
      !
      l2_n = l2_n + rho_n*rho_n*wtdjac(n)
      l2_d = l2_d + rho_d*rho_d*wtdjac(n)
      !
    end do
    !
    marker(nc) = max(l2_n,rndoff) / l2_d
    !
  end do
  !
  ! Compute the log10 of marker and scaled_marker. For any skipped cells,
  ! set marker and scaled_marker to values that will prevent them from
  ! being marked as troubled cells.
  !
  !$ACC PARALLEL LOOP PRESENT(skip_cell,cutoff_order,marker,scaled_marker)
  do nc = 1,ncell
    !
    if (skip_cell(nc)) then
      !
      marker(nc) = -two*ten
      scaled_marker(nc) = zero
      !
    else
      !
      l10_marker = log10( marker(nc) )
      !
      psi0 = -one*(c1 + c2*log10(real(cutoff_order(nc),kind=wp)))
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
  !$ACC END DATA
  deallocate ( rho_pm1 )
  !
end subroutine compute_resolution_marker
!
!###############################################################################
!
pure subroutine res_ind_apply_proj_3d(usp,skip_cell,lmtr_order)
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
pure subroutine res_ind_apply_proj_2d(usp,skip_cell,lmtr_order)
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
pure subroutine res_ind_apply_proj_nd(usp,skip_cell,lmtr_order)
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
    do m = 1,nq
      do n = 1,np
        old_usp(n,m) = usp(n0+n,m)
      end do
    end do
    !
    ! Apply the projection limiter by multplying the projection
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do m = 1,nq
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
!
subroutine resolution_indicator_old(apply_projection,troubled_cell, &
                                    orders_to_cutoff)
  !
  !.. Use Statements ..
  use geovar,         only : ncell,cell
  use order_mod,      only : maxSP,p_order,n_order
  use flowvar,        only : usp
  use metrics_mod,    only : wtdjac
  use filter_mod,     only : filter
  use projection_mod, only : projct
  !
  !.. Optional Arguments ..
  logical(lk),  optional,    intent(in) :: apply_projection
  logical(ldk), optional, intent(inout) :: troubled_cell(:)
  integer,      optional,    intent(in) :: orders_to_cutoff
  !
  !.. Local Scalars ..
  integer  :: i,k,m,n,nc,n1,n2,np,ierr
  integer  :: this_geom,this_order
  integer  :: ord_mlt,ord_add
  integer  :: cutoff_order,order
  real(wp) :: a,val,l2_p,l2_pm1,psi0
  !
  !.. Local Arrays ..
  real(wp) :: rho_p(1:maxSP)
  real(wp) :: rho_pm1(1:maxSP)
  real(wp) :: tusp(1:maxSP,1:nq)
  real(wp) :: usp_check(1:nq,1:maxSP)
  real(wp) :: usp_copy(1:nq,1:maxSP)
  integer  :: nlim(0:n_order)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "resolution_indicator_old"
  real(wp),         parameter :: dpsi = half
  real(wp),         parameter :: c1 = four
  real(wp),         parameter :: c2 = 4.25_wp
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
  if (p_order < n_order) then
    ord_mlt = 0
    ord_add = p_order
  else
    ord_mlt = 1
   !ord_add = 0
    ord_add = -1
  end if
  !
  do nc = 1,ncell
    !
    if (present(troubled_cell)) then
      if (.not. troubled_cell(nc)) then
        marker(nc) = -two*ten
        scaled_marker(nc) = zero
        cycle
      end if
    end if
    !
    n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
    n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
    np = n2-n1+1         ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    cutoff_order = ord_mlt*this_order + ord_add
    if (present(orders_to_cutoff)) then
      cutoff_order = cutoff_order-orders_to_cutoff+1
   !  write (iout,'(*(2x,i0))') this_order,cutoff_order,orders_to_cutoff
   !else
   !  write (iout,'(*(2x,i0))') this_order,cutoff_order
    end if
    !
    ! Copy u, the degree-P polynomial solution of density at the solution points
    !
    rho_p(1:np) = usp(n1:n2,nec)
    !
    ! Compute u_hat, the degree-(P-1) polynomial solution
    ! of density at the solution points
    !
    do n = 1,np
      !
      val = zero
      !
      do i = filter(this_geom,this_order)%cutoff(cutoff_order)%col_ptr(n)+1, &
             filter(this_geom,this_order)%cutoff(cutoff_order)%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%cutoff(cutoff_order)%col_idx(i)
        a = filter(this_geom,this_order)%cutoff(cutoff_order)%col_val(i)
        !
        val = val + a * usp(n1-1+k,nec)
        !
      end do
      !
      rho_pm1(n) = val
      !
    end do
    !
    ! Compute u - u_hat
    !
    rho_pm1(1:np) = rho_p(1:np) - rho_pm1(1:np)
    !
    ! Compute the L2 inner product of the two polynomials
    !
    l2_pm1 = sum( rho_pm1(1:np) * rho_pm1(1:np) * wtdjac(n1:n2) )
    l2_pm1 = max(l2_pm1,rndoff)
    !
    l2_p = sum( rho_p(1:np) * rho_p(1:np) * wtdjac(n1:n2) )
    !
    marker(nc) = log10( l2_pm1 / l2_p )
    !
    psi0 = -one*(c1 + c2*log10(real(cutoff_order,kind=wp)))
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
    troubled_cell(:) = (scaled_marker(:) >= 1.0e-8_wp)
  end if
  !
  if (present(apply_projection)) then
    if (apply_projection) then
      !
      nlim(:) = 0
      !
      do nc = 1,ncell
        !
        if (present(troubled_cell)) then
          if (.not. troubled_cell(nc)) cycle
        end if
        !
        n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
        n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
        np = n2-n1+1         ! number of sol-pts in this cell
        !
        this_geom = cell(nc)%geom   ! geometry type of this cell
        this_order = cell(nc)%order ! solution order of this cell
        !
        order = int( real(this_order,kind=wp)*(one-scaled_marker(nc)) )
        order = max(1,min(this_order,order))
        !
        if (order < this_order) then
          !
          nlim(order) = nlim(order) + 1
          !
          ! Store the transpose of usp for this cell to reduce cache misses
          !
          do k = 1,np
            do m = 1,nq
              tusp(k,m) = usp(n1-1+k,m)
            end do
          end do
          !
          ! Project the solution down to proj_order
          ! at each solution point of this cell
          !
          do m = 1,nq
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
                val = val + a * tusp(k,m)
                !
              end do
              !
              usp(n1-1+n,m) = val
              !
            end do
          end do
          !
        end if
        !
      end do
      !
#ifdef DEBUG_ON
      if (ncpu > 1) then
        call mpi_allreduce(MPI_IN_PLACE,nlim,size(nlim,kind=int_mpi), &
                           mpi_inttyp,MPI_SUM,MPI_COMM_WORLD,mpierr)
      end if
      if (mypnum == 0) then
        do n = lbound(nlim,dim=1),ubound(nlim,dim=1)
          if (nlim(n) > 0) write (iout,1) nlim(n),n
        end do
      end if
#endif
      !
    end if
  end if
  !
  1 format (8x,i8," Cells limited to P",i0)
  !
end subroutine resolution_indicator_old
!
!###############################################################################
!
subroutine jump_indicator(troubled_cell)
  !
  !.. Use Statements ..
  use geovar,         only : nr,n_flxpts,ncell,cell,nface,face
  use quadrature_mod, only : std_elem
  use flowvar,        only : faceusp_array
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
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: p_jump_fp(:)
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
  allocate ( p_jump_fp(1:n_flxpts) )
  allocate ( pres_jump(1:nface) )
  !$ACC DATA CREATE(p_jump_fp,pres_jump)
  !
  call compute_pressure_jump(faceusp_array,p_jump_fp)
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP GANG &
    !$ACC&              VECTOR_LENGTH(32) &
    !$ACC&              PRESENT(face,p_jump_fp,pres_jump,std_elem)
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
        val = val + ds*p_jump_fp(f0+k)*std_elem(face_geom,face_order)%wts(k)
        !
      end do
      !
      pres_jump(nf) = val
      !
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP GANG &
    !$ACC&              VECTOR_LENGTH(32) &
    !$ACC&              PRESENT(face,p_jump_fp,pres_jump,std_elem)
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
        val = val + ds*p_jump_fp(f0+k)*std_elem(face_geom,face_order)%wts(k)
        !
      end do
      !
      pres_jump(nf) = val
      !
    end do
  else
    !$ACC PARALLEL LOOP GANG &
    !$ACC&              VECTOR_LENGTH(32) &
    !$ACC&              PRESENT(face,p_jump_fp,pres_jump,std_elem)
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
        val = val + ds*p_jump_fp(f0+k)*std_elem(face_geom,face_order)%wts(k)
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
  deallocate ( p_jump_fp )
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
!
subroutine compute_pressure_jump(face_cv,pres_jump)
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp), intent(inout) :: pres_jump(:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_pressure_jump"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    call compute_pressure_jump_3d(face_cv,pres_jump)
  else if (nr == 2) then
    call compute_pressure_jump_2d(face_cv,pres_jump)
  else
    call compute_pressure_jump_1d(face_cv,pres_jump)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_pressure_jump
!
!###############################################################################
!
pure subroutine compute_pressure_jump_3d(face_cv,pres_jump)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp), intent(inout) :: pres_jump(:)
  !
  !.. Local Scalars ..
  integer  :: k,np
  real(wp) :: gm1
  real(wp) :: rho_l,rhou_l,rhov_l,rhow_l,rhoe_l,rhoq2_l,pres_l
  real(wp) :: rho_r,rhou_r,rhov_r,rhow_r,rhoe_r,rhoq2_r,pres_r
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_pressure_jump_3d"
  !
continue
  !
  np = size(pres_jump)
  !
  gm1 = gam-one
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,pres_jump)
  do k = 1,np
    !
    rho_l = face_cv(k,nec,1)
    rhou_l = face_cv(k,nmx,1)
    rhov_l = face_cv(k,nmy,1)
    rhow_l = face_cv(k,nmz,1)
    rhoe_l = face_cv(k,nee,1)
    !
    rhoq2_l = rhou_l*rhou_l + rhov_l*rhov_l + rhow_l*rhow_l
    !
    pres_l = gm1*(rhoe_l - half*rhoq2_l/rho_l)
    !
    rho_r = face_cv(k,nec,2)
    rhou_r = face_cv(k,nmx,2)
    rhov_r = face_cv(k,nmy,2)
    rhow_r = face_cv(k,nmz,2)
    rhoe_r = face_cv(k,nee,2)
    !
    rhoq2_r = rhou_r*rhou_r + rhov_r*rhov_r + rhow_r*rhow_r
    !
    pres_r = gm1*(rhoe_r - half*rhoq2_r/rho_r)
    !
    pres_jump(k) = two * abs(pres_l - pres_r) / (pres_l + pres_r)
    !
  end do
  !
end subroutine compute_pressure_jump_3d
!
!###############################################################################
!
pure subroutine compute_pressure_jump_2d(face_cv,pres_jump)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp), intent(inout) :: pres_jump(:)
  !
  !.. Local Scalars ..
  integer  :: k,np
  real(wp) :: gm1
  real(wp) :: rho_l,rhou_l,rhov_l,rhoe_l,rhoq2_l,pres_l
  real(wp) :: rho_r,rhou_r,rhov_r,rhoe_r,rhoq2_r,pres_r
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_pressure_jump_2d"
  !
continue
  !
  np = size(pres_jump)
  !
  gm1 = gam-one
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,pres_jump)
  do k = 1,np
    !
    rho_l = face_cv(k,nec,1)
    rhou_l = face_cv(k,nmx,1)
    rhov_l = face_cv(k,nmy,1)
    rhoe_l = face_cv(k,nee,1)
    !
    rhoq2_l = rhou_l*rhou_l + rhov_l*rhov_l
    !
    pres_l = gm1*(rhoe_l - half*rhoq2_l/rho_l)
    !
    rho_r = face_cv(k,nec,2)
    rhou_r = face_cv(k,nmx,2)
    rhov_r = face_cv(k,nmy,2)
    rhoe_r = face_cv(k,nee,2)
    !
    rhoq2_r = rhou_r*rhou_r + rhov_r*rhov_r
    !
    pres_r = gm1*(rhoe_r - half*rhoq2_r/rho_r)
    !
    pres_jump(k) = two * abs(pres_l - pres_r) / (pres_l + pres_r)
    !
  end do
  !
end subroutine compute_pressure_jump_2d
!
!###############################################################################
!
pure subroutine compute_pressure_jump_1d(face_cv,pres_jump)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp), intent(inout) :: pres_jump(:)
  !
  !.. Local Scalars ..
  integer  :: k,np
  real(wp) :: gm1
  real(wp) :: rho_l,rhou_l,rhoe_l,pres_l
  real(wp) :: rho_r,rhou_r,rhoe_r,pres_r
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_pressure_jump_1d"
  !
continue
  !
  np = size(pres_jump)
  !
  gm1 = gam-one
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,pres_jump)
  do k = 1,np
    !
    rho_l = face_cv(k,nec,1)
    rhou_l = face_cv(k,nmx,1)
    rhoe_l = face_cv(k,nee,1)
    !
    pres_l = gm1*(rhoe_l - half*rhou_l*rhou_l/rho_l)
    !
    rho_r = face_cv(k,nec,2)
    rhou_r = face_cv(k,nmx,2)
    rhoe_r = face_cv(k,nee,2)
    !
    pres_r = gm1*(rhoe_r - half*rhou_r*rhou_r/rho_r)
    !
    pres_jump(k) = two * abs(pres_l - pres_r) / (pres_l + pres_r)
    !
  end do
  !
end subroutine compute_pressure_jump_1d
!
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
include "Functions/pressure.f90"
!
!###############################################################################
!
end module module_limiters
