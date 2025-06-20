module limiters_mod
  !
  !.. Use Statements ..
  use module_kind_types, global_debug_timer =>      debug_timer, &
                    pure_global_debug_timer => pure_debug_timer
  !
  use eqn_idx_mod, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq
  use eqn_idx_mod, only : ncs,npres,ntemp
  !
  use ovar, only : force_nd_subroutines
  !
  ! Use nmarker, marker, and scaled_marker from their new temporary location
  ! in flowvar_mod to remove module circular dependency. This use statement is
  ! added to the module header since a lot of functions/subroutines in this
  ! module assume these entities exists in this module header and access them
  ! using "host association". This module will be cleaned up and reorganized in
  ! the near future, and this simplifies the temporary "fix" until that happens.
  use flowvar_mod, only : nmarker,marker,scaled_marker
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
  public :: ducros_indicator
  public :: cell_length_scale_max_box
  public :: apply_sponge_filter
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
  real(wp), save, allocatable :: umodes(:,:)
  !
  logical(ldk), save, allocatable :: troubled_cell(:)
  !$ACC DECLARE CREATE(troubled_cell)
  !
  logical(ldk), save, allocatable :: sponge_cell(:)
  !$ACC DECLARE CREATE(sponge_cell)
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
  use geomvar_mod, only : ncell,n_totpts
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
    !
    ! Create the array for identifying troubled cells
    !
    allocate ( troubled_cell(1:ncell) , source=.false. , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"troubled_cell",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !$ACC UPDATE DEVICE(troubled_cell)
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
subroutine moment_limiter
  !
  !.. Use Statements ..
  use geomvar_mod,     only : ncell,cell
  use flowvar_mod,     only : usp
  use vandermonde_mod, only : vand
  !
 !use parallel_mod, only : exchange_modes
  !
  !.. Local Scalars ..
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
subroutine limit_by_projecting(nst)
  !
  !.. Use Statements ..
  use order_mod,   only : n_max_order,n_order
  use flowvar_mod, only : usp
  use ovar,        only : Limiter_Filter_Option
  use ovar,        only : output_limiter_statistics,nlim
  !
  !.. Arguments ..
  integer, intent(in) :: nst  !< Runge-Kutta stage
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
  call identify_troubled_cells(nst,troubled_cell,init_troubled_cells, &
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
    101 format (9x,      i8," troubled cells at P",i0)
#endif
    !
    if (Limiter_Filter_Option /= 0) then
      call apply_filter_limiter_to_usp(troubled_cell,usp)
      exit order_loop
    else
      call apply_projct_limiter_to_usp(troubled_cell,usp,dropped_orders)
    end if
    !
    ! Recheck the troubled cells to see if there are no longer
    ! troubled due to the previous projection.
    !
    call identify_troubled_cells(nst,troubled_cell,dropped_orders, &
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
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine apply_sponge_filter(usp)
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "apply_sponge_filter"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (.not. allocated(sponge_cell)) then
    call mark_sponge_filter_cells
  end if
  !
  nr = nme-nmb+1
  !
  if (force_nd_subroutines) then
    !
    call apply_sponge_filter_nd(1,nq,usp)
    !
  else
    !
    if (nr == 3) then
      call apply_sponge_filter_3d(usp)
    else if (nr == 2) then
      call apply_sponge_filter_2d(usp)
    else
      call apply_sponge_filter_1d(usp)
    end if
    !
    if (ncs > 0) then
      call apply_sponge_filter_nd(1,ncs,usp)
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine apply_sponge_filter
!
!###############################################################################
!
pure &
subroutine apply_sponge_filter_3d(usp)
  !
  !.. Use Statements ..
  use order_mod,   only : maxSP
  use geomvar_mod, only : ncell,cell
  use filter_mod,  only : filter
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
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "apply_sponge_filter_3d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,sponge_cell,filter) &
  !$ACC&              PRIVATE(old_usp)
  filt_cell_loop: do nc = 1,ncell
    !
    ! If the current cell is not a sponge cell, skip to the next one
    !
    if (.not.sponge_cell(nc)) then
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
    ! Apply the sponge filter by multplying the filter
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
      do i = filter(this_geom,this_order)%sponge%col_ptr(n)+1, &
             filter(this_geom,this_order)%sponge%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%sponge%col_idx(i)
        a = filter(this_geom,this_order)%sponge%col_val(i)
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
  end do filt_cell_loop
  !
end subroutine apply_sponge_filter_3d
!
!###############################################################################
!
pure &
subroutine apply_sponge_filter_2d(usp)
  !
  !.. Use Statements ..
  use order_mod,   only : maxSP
  use geomvar_mod, only : ncell,cell
  use filter_mod,  only : filter
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
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "apply_sponge_filter_2d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,sponge_cell,filter) &
  !$ACC&              PRIVATE(old_usp)
  filt_cell_loop: do nc = 1,ncell
    !
    ! If the current cell is not a sponge cell, skip to the next one
    !
    if (.not.sponge_cell(nc)) then
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
    ! Apply the sponge filter by multplying the filter
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
      do i = filter(this_geom,this_order)%sponge%col_ptr(n)+1, &
             filter(this_geom,this_order)%sponge%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%sponge%col_idx(i)
        a = filter(this_geom,this_order)%sponge%col_val(i)
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
  end do filt_cell_loop
  !
end subroutine apply_sponge_filter_2d
!
!###############################################################################
!
pure &
subroutine apply_sponge_filter_1d(usp)
  !
  !.. Use Statements ..
  use order_mod,   only : maxSP
  use geomvar_mod, only : ncell,cell
  use filter_mod,  only : filter
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
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "apply_sponge_filter_1d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,sponge_cell,filter) &
  !$ACC&              PRIVATE(old_usp)
  filt_cell_loop: do nc = 1,ncell
    !
    ! If the current cell is not a sponge cell, skip to the next one
    !
    if (.not.sponge_cell(nc)) then
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
    ! Apply the sponge filter by multplying the filter
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      filt_ec = zero
      filt_mx = zero
      filt_ee = zero
      !
      do i = filter(this_geom,this_order)%sponge%col_ptr(n)+1, &
             filter(this_geom,this_order)%sponge%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%sponge%col_idx(i)
        a = filter(this_geom,this_order)%sponge%col_val(i)
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
  end do filt_cell_loop
  !
end subroutine apply_sponge_filter_1d
!
!###############################################################################
!
pure &
subroutine apply_sponge_filter_nd(mbeg,mend,usp)
  !
  !.. Use Statements ..
  use order_mod,   only : maxSP
  use geomvar_mod, only : ncell,cell
  use filter_mod,  only : filter
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: mbeg,mend
  real(wp), intent(inout) :: usp(:,:)
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
  character(len=*), parameter :: pname = "apply_sponge_filter_nd"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,sponge_cell,filter) &
  !$ACC&              PRIVATE(old_usp)
  filt_cell_loop: do nc = 1,ncell
    !
    ! If the current cell is not a sponge cell, skip to the next one
    !
    if (.not.sponge_cell(nc)) then
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
    ! Apply the sponge filter by multplying the filter
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do m = mbeg,mend
      do n = 1,np
        !
        val = zero
        !
        do i = filter(this_geom,this_order)%sponge%col_ptr(n)+1, &
               filter(this_geom,this_order)%sponge%col_ptr(n+1)
          !
          k = filter(this_geom,this_order)%sponge%col_idx(i)
          a = filter(this_geom,this_order)%sponge%col_val(i)
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
  end do filt_cell_loop
  !
end subroutine apply_sponge_filter_nd
!
!###############################################################################
!
subroutine mark_sponge_filter_cells()
  !
  !.. Use Statements ..
  use ovar,        only : Sponge_Beyond_X,grid_scaling_factor
  use ovar,        only : beta_aoaref,alpha_aoaref
  use geomvar_mod, only : ncell,cell,xyz
  !
  !.. Local Scalars ..
  integer  :: nc,n1,n2,i
  real(wp) :: xmax
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "mark_sponge_filter_cells"
  !
continue
  !
  ! Allocate the sponge_cell array
  !
  sponge_cell = init_array(ncell,.false.) ! F2003 AUTO-REALLOC
  !$ACC UPDATE DEVICE(sponge_cell)
  !
  ! Find all the cells in the sponge regions
  !
  xmax = Sponge_Beyond_X * grid_scaling_factor
  !
  if (abs(beta_aoaref-nine*ten) < eps6) then
    i = 3
  else if (abs(alpha_aoaref-nine*ten) < eps6) then
    i = 2
  else
    i = 1
  end if
  !$ACC PARALLEL LOOP PRESENT(cell,xyz,sponge_cell)
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    !
    if (all(xyz(i,n1:n2) > xmax)) sponge_cell(nc) = .true.
    !
  end do
  !
end subroutine mark_sponge_filter_cells
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
subroutine apply_filter_limiter_to_usp(troubled_cell,usp)
  !
  !.. Use Statements ..
  use flowvar_mod, only : cell_diss
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
    call apply_filter_limiter_to_usp_nd(1,nq,troubled_cell,usp,cell_diss)
    !
  else
    !
    if (nr == 3) then
      call apply_filter_limiter_to_usp_3d(troubled_cell,usp,cell_diss)
    else if (nr == 2) then
      call apply_filter_limiter_to_usp_2d(troubled_cell,usp,cell_diss)
    else
      call apply_filter_limiter_to_usp_1d(troubled_cell,usp,cell_diss)
    end if
    !
    if (ncs > 0) then
      call apply_filter_limiter_to_usp_nd(1,ncs,troubled_cell,usp,cell_diss)
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
subroutine apply_filter_limiter_to_usp_3d(troubled_cell,usp,cell_diss)
  !
  !.. Use Statements ..
  use order_mod,   only : maxSP
  use geomvar_mod, only : ncell,cell
  use flowvar_mod, only : get_marker_index
  use filter_mod,  only : filter
  use ovar,        only : nst => current_rk_stage
  use ovar,        only : Limiter_Use_Dynamic_Filter
  use ovar,        only : Track_Limiter_Filter_Dissipation
  use metrics_mod, only : wtdjac
  !
  !.. Formal Arguments ..
  logical(ldk), intent(inout) :: troubled_cell(:)
  real(wp),     intent(inout) :: usp(:,:)
  real(wp),     intent(inout) :: cell_diss(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0,midx
  integer  :: this_geom,this_order
  real(wp) :: a,filt_ec,filt_mx,filt_my,filt_mz,filt_ee,s,oms
  real(wp) :: diss_ec,diss_mx,diss_my,diss_mz,diss_ee
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "apply_filter_limiter_to_usp_3d"
  !
continue
  !
  midx = get_marker_index(nst)
  s = one
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
    ! Extract the scaled marker value for this cell if
    ! using the scaled marker for dynamic filtering
    ! NOTE: If not using dynamic filtering, s will remain unchanged from
    !       its initial value of 1, which will leave the diagonal values
    !       of the filter matrix unmodified.
    !
    if (Limiter_Use_Dynamic_Filter) then
      s = scaled_marker(nc,midx)
    end if
    oms = one-s
    !
    ! Apply the filter limiter by multplying the filter
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      filt_ec = oms*old_usp(n,nec)
      filt_mx = oms*old_usp(n,nmx)
      filt_my = oms*old_usp(n,nmy)
      filt_mz = oms*old_usp(n,nmz)
      filt_ee = oms*old_usp(n,nee)
      !
      do i = filter(this_geom,this_order)%limiter%col_ptr(n)+1, &
             filter(this_geom,this_order)%limiter%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%limiter%col_idx(i)
        a = filter(this_geom,this_order)%limiter%col_val(i)
        !
        ! Use the scaled marker to linearly scale the strength of the 
        ! filter. The scaling is applied to all elements of the base
        ! filter matrix, and a (1-scaling) term needs to be added to
        ! the diagonal elements of the matrix.
        ! NOTE: Fmat = Vmat.f.Vinv   
        !          where : Fmat = base filter matrix
        !                  f = diagonal matrix of the filter coefficients
        !       Scaled(Fmat) = SFmat = Vmat.[(1-s)*I + s*f).Vinv
        !          where : s = scaled marker (0=no filter, 1=full filter)
        !                  I = identity matrix
        !       SFmat = (1-s)*Vmat.I.Vinv + s*Vmat.f.Vinv
        !             = (1-s)*I + s*Fmat
        !       We can use the (1-s)*I part as the initialization of the
        !       filter sums because that value is independent of the actual
        !       filter matrix.
        !
       !a = a*s + merge(one-s,zero,n==k)
        !
        filt_ec = filt_ec + a * s * old_usp(k,nec)
        filt_mx = filt_mx + a * s * old_usp(k,nmx)
        filt_my = filt_my + a * s * old_usp(k,nmy)
        filt_mz = filt_mz + a * s * old_usp(k,nmz)
        filt_ee = filt_ee + a * s * old_usp(k,nee)
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
!#define DISS_FUNC( a , b ) (a-b)
!#define DISS_TOTAL( a ) (a)
!#define DISS_FUNC( a , b ) (abs(a-b))
#define DISS_FUNC( a , b ) ((a-b)**2)
#define DISS_TOTAL( a ) (sqrt(a))

    if (Track_Limiter_Filter_Dissipation) then
      diss_ec = zero
      diss_mx = zero
      diss_my = zero
      diss_mz = zero
      diss_ee = zero
      do n = 1,np
        diss_ec = diss_ec + DISS_FUNC(old_usp(n,nec),usp(n0+n,nec))*wtdjac(n0+n)
        diss_mx = diss_mx + DISS_FUNC(old_usp(n,nmx),usp(n0+n,nmx))*wtdjac(n0+n)
        diss_my = diss_my + DISS_FUNC(old_usp(n,nmy),usp(n0+n,nmy))*wtdjac(n0+n)
        diss_mz = diss_mz + DISS_FUNC(old_usp(n,nmz),usp(n0+n,nmz))*wtdjac(n0+n)
        diss_ee = diss_ee + DISS_FUNC(old_usp(n,nee),usp(n0+n,nee))*wtdjac(n0+n)
      end do
      cell_diss(nc,nec) = cell_diss(nc,nec) + DISS_TOTAL(diss_ec)
      cell_diss(nc,nmx) = cell_diss(nc,nmx) + DISS_TOTAL(diss_mx)
      cell_diss(nc,nmy) = cell_diss(nc,nmy) + DISS_TOTAL(diss_my)
      cell_diss(nc,nmz) = cell_diss(nc,nmz) + DISS_TOTAL(diss_mz)
      cell_diss(nc,nee) = cell_diss(nc,nee) + DISS_TOTAL(diss_ee)
    end if
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
subroutine apply_filter_limiter_to_usp_2d(troubled_cell,usp,cell_diss)
  !
  !.. Use Statements ..
  use order_mod,   only : maxSP
  use geomvar_mod, only : ncell,cell
  use flowvar_mod, only : get_marker_index
  use filter_mod,  only : filter
  use ovar,        only : nst => current_rk_stage
  use ovar,        only : Limiter_Use_Dynamic_Filter
  use ovar,        only : Track_Limiter_Filter_Dissipation
  use metrics_mod, only : wtdjac
  !
  !.. Formal Arguments ..
  logical(ldk), intent(inout) :: troubled_cell(:)
  real(wp),     intent(inout) :: usp(:,:)
  real(wp),     intent(inout) :: cell_diss(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0,midx
  integer  :: this_geom,this_order
  real(wp) :: a,filt_ec,filt_mx,filt_my,filt_ee,s,oms
  real(wp) :: diss_ec,diss_mx,diss_my,diss_ee
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "apply_filter_limiter_to_usp_2d"
  !
continue
  !
  midx = get_marker_index(nst)
  s = one
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
    ! Extract the scaled marker value for this cell if
    ! using the scaled marker for dynamic filtering
    ! NOTE: If not using dynamic filtering, s will remain unchanged from
    !       its initial value of 1, which will leave the diagonal values
    !       of the filter matrix unmodified.
    !
    if (Limiter_Use_Dynamic_Filter) then
      s = scaled_marker(nc,midx)
    end if
    oms = one-s
    !
    ! Apply the filter limiter by multplying the filter
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      filt_ec = oms*old_usp(n,nec)
      filt_mx = oms*old_usp(n,nmx)
      filt_my = oms*old_usp(n,nmy)
      filt_ee = oms*old_usp(n,nee)
      !
      do i = filter(this_geom,this_order)%limiter%col_ptr(n)+1, &
             filter(this_geom,this_order)%limiter%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%limiter%col_idx(i)
        a = filter(this_geom,this_order)%limiter%col_val(i)
        !
        ! Use the scaled marker to linearly scale the strength of the 
        ! filter. The scaling is applied to all elements of the base
        ! filter matrix, and a (1-scaling) term needs to be added to
        ! the diagonal elements of the matrix.
        ! NOTE: Fmat = Vmat.f.Vinv   
        !          where : Fmat = base filter matrix
        !                  f = diagonal matrix of the filter coefficients
        !       Scaled(Fmat) = SFmat = Vmat.[(1-s)*I + s*f).Vinv
        !          where : s = scaled marker (0=no filter, 1=full filter)
        !                  I = identity matrix
        !       SFmat = (1-s)*Vmat.I.Vinv + s*Vmat.f.Vinv
        !             = (1-s)*I + s*Fmat
        !       We can use the (1-s)*I part as the initialization of the
        !       filter sums because that value is independent of the actual
        !       filter matrix.
        !
       !a = a*s + merge(one-s,zero,n==k)
        !
        filt_ec = filt_ec + a * s * old_usp(k,nec)
        filt_mx = filt_mx + a * s * old_usp(k,nmx)
        filt_my = filt_my + a * s * old_usp(k,nmy)
        filt_ee = filt_ee + a * s * old_usp(k,nee)
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
    if (Track_Limiter_Filter_Dissipation) then
      diss_ec = zero
      diss_mx = zero
      diss_my = zero
      diss_ee = zero
      do n = 1,np
        diss_ec = diss_ec + DISS_FUNC(old_usp(n,nec),usp(n0+n,nec))*wtdjac(n0+n)
        diss_mx = diss_mx + DISS_FUNC(old_usp(n,nmx),usp(n0+n,nmx))*wtdjac(n0+n)
        diss_my = diss_my + DISS_FUNC(old_usp(n,nmy),usp(n0+n,nmy))*wtdjac(n0+n)
        diss_ee = diss_ee + DISS_FUNC(old_usp(n,nee),usp(n0+n,nee))*wtdjac(n0+n)
      end do
      cell_diss(nc,nec) = cell_diss(nc,nec) + DISS_TOTAL(diss_ec)
      cell_diss(nc,nmx) = cell_diss(nc,nmx) + DISS_TOTAL(diss_mx)
      cell_diss(nc,nmy) = cell_diss(nc,nmy) + DISS_TOTAL(diss_my)
      cell_diss(nc,nee) = cell_diss(nc,nee) + DISS_TOTAL(diss_ee)
    end if
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
subroutine apply_filter_limiter_to_usp_1d(troubled_cell,usp,cell_diss)
  !
  !.. Use Statements ..
  use order_mod,   only : maxSP
  use geomvar_mod, only : ncell,cell
  use flowvar_mod, only : get_marker_index
  use filter_mod,  only : filter
  use ovar,        only : nst => current_rk_stage
  use ovar,        only : Limiter_Use_Dynamic_Filter
  use ovar,        only : Track_Limiter_Filter_Dissipation
  use metrics_mod, only : wtdjac
  !
  !.. Formal Arguments ..
  logical(ldk), intent(inout) :: troubled_cell(:)
  real(wp),     intent(inout) :: usp(:,:)
  real(wp),     intent(inout) :: cell_diss(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0,midx
  integer  :: this_geom,this_order
  real(wp) :: a,filt_ec,filt_mx,filt_ee,s,oms
  real(wp) :: diss_ec,diss_mx,diss_ee
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "apply_filter_limiter_to_usp_1d"
  !
continue
  !
  midx = get_marker_index(nst)
  s = one
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
    ! Extract the scaled marker value for this cell if
    ! using the scaled marker for dynamic filtering
    ! NOTE: If not using dynamic filtering, s will remain unchanged from
    !       its initial value of 1, which will leave the diagonal values
    !       of the filter matrix unmodified.
    !
    if (Limiter_Use_Dynamic_Filter) then
      s = scaled_marker(nc,midx)
    end if
    oms = one-s
    !
    ! Apply the filter limiter by multplying the filter
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      filt_ec = oms*old_usp(n,nec)
      filt_mx = oms*old_usp(n,nmx)
      filt_ee = oms*old_usp(n,nee)
      !
      do i = filter(this_geom,this_order)%limiter%col_ptr(n)+1, &
             filter(this_geom,this_order)%limiter%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%limiter%col_idx(i)
        a = filter(this_geom,this_order)%limiter%col_val(i)
        !
        ! Use the scaled marker to linearly scale the strength of the 
        ! filter. The scaling is applied to all elements of the base
        ! filter matrix, and a (1-scaling) term needs to be added to
        ! the diagonal elements of the matrix.
        ! NOTE: Fmat = Vmat.f.Vinv   
        !          where : Fmat = base filter matrix
        !                  f = diagonal matrix of the filter coefficients
        !       Scaled(Fmat) = SFmat = Vmat.[(1-s)*I + s*f).Vinv
        !          where : s = scaled marker (0=no filter, 1=full filter)
        !                  I = identity matrix
        !       SFmat = (1-s)*Vmat.I.Vinv + s*Vmat.f.Vinv
        !             = (1-s)*I + s*Fmat
        !       We can use the (1-s)*I part as the initialization of the
        !       filter sums because that value is independent of the actual
        !       filter matrix.
        !
       !a = a*s + merge(one-s,zero,n==k)
        !
        filt_ec = filt_ec + a * s * old_usp(k,nec)
        filt_mx = filt_mx + a * s * old_usp(k,nmx)
        filt_ee = filt_ee + a * s * old_usp(k,nee)
        !
      end do
      !
      usp(n0+n,nec) = filt_ec
      usp(n0+n,nmx) = filt_mx
      usp(n0+n,nee) = filt_ee
      !
    end do
    !
    if (Track_Limiter_Filter_Dissipation) then
      diss_ec = zero
      diss_mx = zero
      diss_ee = zero
      do n = 1,np
        diss_ec = diss_ec + DISS_FUNC(old_usp(n,nec),usp(n0+n,nec))*wtdjac(n0+n)
        diss_mx = diss_mx + DISS_FUNC(old_usp(n,nmx),usp(n0+n,nmx))*wtdjac(n0+n)
        diss_ee = diss_ee + DISS_FUNC(old_usp(n,nee),usp(n0+n,nee))*wtdjac(n0+n)
      end do
      cell_diss(nc,nec) = cell_diss(nc,nec) + DISS_TOTAL(diss_ec)
      cell_diss(nc,nmx) = cell_diss(nc,nmx) + DISS_TOTAL(diss_mx)
      cell_diss(nc,nee) = cell_diss(nc,nee) + DISS_TOTAL(diss_ee)
    end if
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
subroutine apply_filter_limiter_to_usp_nd(mbeg,mend,troubled_cell,usp,cell_diss)
  !
  !.. Use Statements ..
  use order_mod,   only : maxSP
  use geomvar_mod, only : ncell,cell
  use flowvar_mod, only : get_marker_index
  use filter_mod,  only : filter
  use ovar,        only : nst => current_rk_stage
  use ovar,        only : Limiter_Use_Dynamic_Filter
  use ovar,        only : Track_Limiter_Filter_Dissipation
  use metrics_mod, only : wtdjac
  !
  !.. Formal Arguments ..
  integer,         intent(in) :: mbeg,mend
  logical(ldk), intent(inout) :: troubled_cell(:)
  real(wp),     intent(inout) :: usp(:,:)
  real(wp),     intent(inout) :: cell_diss(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,m,n,nc,np,n0,midx
  integer  :: this_geom,this_order
  real(wp) :: a,val,diss,s,oms
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "apply_filter_limiter_to_usp_nd"
  !
continue
  !
  midx = get_marker_index(nst)
  s = one
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
    ! Extract the scaled marker value for this cell if
    ! using the scaled marker for dynamic filtering
    ! NOTE: If not using dynamic filtering, s will remain unchanged from
    !       its initial value of 1, which will leave the diagonal values
    !       of the filter matrix unmodified.
    !
    if (Limiter_Use_Dynamic_Filter) then
      s = scaled_marker(nc,midx)
    end if
    oms = one-s
    !
    ! Apply the filter limiter by multplying the filter
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do m = mbeg,mend
      do n = 1,np
        !
        val = oms*old_usp(n,m)
        !
        do i = filter(this_geom,this_order)%limiter%col_ptr(n)+1, &
               filter(this_geom,this_order)%limiter%col_ptr(n+1)
          !
          k = filter(this_geom,this_order)%limiter%col_idx(i)
          a = filter(this_geom,this_order)%limiter%col_val(i)
          !
          ! Use the scaled marker to linearly scale the strength of the 
          ! filter. The scaling is applied to all elements of the base
          ! filter matrix, and a (1-scaling) term needs to be added to
          ! the diagonal elements of the matrix.
          ! NOTE: Fmat = Vmat.f.Vinv   
          !          where : Fmat = base filter matrix
          !                  f = diagonal matrix of the filter coefficients
          !       Scaled(Fmat) = SFmat = Vmat.[(1-s)*I + s*f).Vinv
          !          where : s = scaled marker (0=no filter, 1=full filter)
          !                  I = identity matrix
          !       SFmat = (1-s)*Vmat.I.Vinv + s*Vmat.f.Vinv
          !             = (1-s)*I + s*Fmat
          !       We can use the (1-s)*I part as the initialization of the
          !       filter sums because that value is independent of the actual
          !       filter matrix.
          !
         !a = a*s + merge(one-s,zero,n==k)
          !
          val = val + a * s * old_usp(k,m)
          !
        end do
        !
        usp(n0+n,m) = val
        !
      end do
    end do
    !
    if (Track_Limiter_Filter_Dissipation) then
      do m = mbeg,mend
        diss = zero
        do n = 1,np
          diss = diss + DISS_FUNC(old_usp(n,m),usp(n0+n,m))*wtdjac(n0+n)
        end do
        cell_diss(nc,m) = cell_diss(nc,m) + DISS_TOTAL(diss)
      end do
    end if
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
  use geomvar_mod,    only : ncell,cell
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
  use geomvar_mod,    only : ncell,cell
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
  use geomvar_mod,    only : ncell,cell
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
  use geomvar_mod,    only : ncell,cell
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
subroutine identify_troubled_cells(nst,troubled_cell,dropped_orders, &
                                   there_are_no_troubles)
  !
  !.. Use Statements ..
  use geomvar_mod, only : ncell,cell
  use ovar,        only : Limiter_Option
  use ovar,        only : Track_Limiter_Filter_Dissipation
  use order_mod,   only : p_order,n_order
  use flowvar_mod, only : cell_num_marks
  !
  !.. Formal Arguments ..
  integer,         intent(in) :: nst                   !< Runge-Kutta stage
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
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "identify_troubled_cells"
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
    if (allocated(sponge_cell)) then
      !
      ! If using the sponge cell limiter:
      !   1. initialize all non-sponge cells as troubled, and
      !   2. initialize all sponge cells as not troubled since these cells
      !      will have already been limited using the sponge filter
      !
      !$ACC PARALLEL LOOP PRESENT(sponge_cell,troubled_cell)
      do nc = 1,ncell
        troubled_cell(nc) = ( .not. sponge_cell(nc))
      end do
      !
    else
      !
      ! If not using the sponge cell limiter, initialize all cells as troubled
      !
      !$ACC PARALLEL LOOP PRESENT(troubled_cell)
      do nc = 1,ncell
        troubled_cell(nc) = .true.
      end do
      !
    end if
    !
  end if
  !
  ! Compute the troubled marker
  !
  if (Limiter_Option == 1) then
   !orders_to_cutoff = proj_order+1
    orders_to_cutoff = proj_order
    call resolution_indicator(nst=nst,troubled_cell=troubled_cell, &
                              orders_to_cutoff=orders_to_cutoff)
  else if (abs(Limiter_Option) == 2) then
    call jump_indicator(nst=nst,troubled_cell=troubled_cell)
  else if (abs(Limiter_Option) == 4) then
    call ducros_indicator(nst=nst,troubled_cell=troubled_cell)
  else
    ! Other values of Limiter_Option should never reach this routine.
    write (error_message,100) Limiter_Option
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
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
    ! Update the counter for how many times this cell has been marked
    !
    if (Track_Limiter_Filter_Dissipation) then
      cell_num_marks(nc) = cell_num_marks(nc) + one
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
100 format( "Unexpected value, Limiter_Option = ", i0 )
  !
end subroutine identify_troubled_cells
!
!###############################################################################
!
subroutine resolution_indicator(nst,apply_projection,troubled_cell, &
                                orders_to_cutoff)
  !
  !.. Use Statements ..
  use order_mod,   only : p_order,n_order
  use geomvar_mod, only : ncell,cell
  use flowvar_mod, only : usp
  use flowvar_mod, only : get_marker_index
  use ovar,        only : output_limiter_statistics,nlim
  use ovar,        only : Limiter_Option
  !
  !.. Optional Arguments ..
  integer,      optional,    intent(in) :: nst              !< Runge-Kutta stage
  logical(lk),  optional,    intent(in) :: apply_projection
  logical(ldk), optional, intent(inout) :: troubled_cell(:)
  integer,      optional,    intent(in) :: orders_to_cutoff
  !
  !.. Local Scalars ..
  integer :: n,nc,nr,ord_mlt,ord_add
  integer :: this_order,limit_order
  integer :: num_limited_cells,ierr
  integer :: imarker
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
    allocate ( marker(1:ncell,1:nmarker) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"marker",1,__LINE__,__FILE__,ierr,error_message)
    !$ACC UPDATE DEVICE(marker)
  end if
  !
  if (.not. allocated(scaled_marker)) then
    allocate ( scaled_marker(1:ncell,1:nmarker) , source=zero , &
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
      if (Limiter_Option < -2 .and. Limiter_Option /= -4) then
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
  ! Determine which marker position to fill based on Runge-Kutta stage.
  !
  imarker = get_marker_index(nst)
  !
  ! Use the resolution indicator to compute marker and scaled_marker
  !
  call compute_resolution_marker(usp,skip_cell,lmtr_order, &
                                 marker(:,imarker),scaled_marker(:,imarker))
  !
  ! Update the troubled_cell indicator if it is present
  !
  if (present(troubled_cell)) then
    !$ACC PARALLEL LOOP PRESENT(scaled_marker,troubled_cell,skip_cell)
    do nc = 1,ncell
      troubled_cell(nc) = (scaled_marker(nc,imarker) >= 1.0e-8_wp)
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
          limit_order = int( real(this_order,kind=wp) &
                           * (one-scaled_marker(nc,imarker)) )
          if (this_order <= 1 .or. &
              (Limiter_Option < -2 .and. Limiter_Option /= -4)) then
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
      call res_ind_apply_proj(usp,skip_cell,lmtr_order)
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
pure &
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
pure &
subroutine compute_resolution_marker_3d(usp,skip_cell,cutoff_order,marker)
  !
  !.. Use Statements ..
  use order_mod,        only : maxSP
  use geomvar_mod,      only : n_solpts,ncell,cell
  use filter_mod,       only : filter
  use element_type_mod, only : std_elem
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
end subroutine compute_resolution_marker_3d
!
!###############################################################################
!
pure &
subroutine compute_resolution_marker_2d(usp,skip_cell,cutoff_order,marker)
  !
  !.. Use Statements ..
  use geomvar_mod,      only : n_solpts,ncell,cell
  use filter_mod,       only : filter
  use element_type_mod, only : std_elem
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
  use geomvar_mod,      only : n_solpts,ncell,cell
  use filter_mod,       only : filter
  use element_type_mod, only : std_elem
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
  use order_mod,        only : maxSP
  use geomvar_mod,      only : n_solpts,ncell,cell
  use filter_mod,       only : filter
  use element_type_mod, only : std_elem
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
  use order_mod,        only : maxSP
  use geomvar_mod,      only : n_solpts,ncell,cell
  use filter_mod,       only : filter
  use element_type_mod, only : std_elem
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
  use geomvar_mod, only : ncell,cell
  use ovar,        only : dpsi => marker_dpsi
  use ovar,        only : c1 => marker_c1
  use ovar,        only : c2 => marker_c2
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
subroutine res_ind_apply_proj(usp,skip_cell,lmtr_order)
  !
  !.. Formal Arguments ..
  logical(ldk), intent(in) :: skip_cell(:)
  integer,      intent(in) :: lmtr_order(:)
  real(wp),  intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "res_ind_apply_proj"
  !
continue
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
end subroutine res_ind_apply_proj
!
!###############################################################################
!
pure &
subroutine res_ind_apply_proj_3d(usp,skip_cell,lmtr_order)
  !
  !.. Use Statements ..
  use order_mod,      only : maxSP
  use geomvar_mod,    only : ncell,cell
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
  use geomvar_mod,    only : ncell,cell
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
  use geomvar_mod,    only : ncell,cell
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
  use geomvar_mod,    only : ncell,cell
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
subroutine jump_indicator(nst,troubled_cell)
  !
  !> Compute the local marker using the jump indicator across adjacent cells.
  !> This is a wrapper subroutine to either call the original version
  !> or the modified version.
  !>
  !> Refs:
  !>  -(2008) Barter,      PhD MIT,        http://hdl.handle.net/1721.1/57606
  !>  -(2007) Barter,      AIAA-2007-3823, doi:10.2514/6.2007-3823
  !> See Also:
  !>  -(2004) Krivodonova, Appl.Num.Math., doi:10.1016/j.apnum.2003.11.002
  !>  -(2013) Bhatia,      AIAA-2013-0302, doi:10.2514/6.2013-302
  !>  -(2014) Bhatia,      PhD U.Florida,
  !>                       https://apps.dtic.mil/sti/citations/tr/ADA624182
  !
  !.. Use Statements ..
  use ovar, only : Limiter_Option
  !
  !.. Optional Arguments ..
  integer,      optional, intent(in)    :: nst              !< Runge-Kutta stage
  logical(ldk), optional, intent(inout) :: troubled_cell(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "jump_indicator"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (Limiter_Option > 0) then
    call jump_indicator_original(nst,troubled_cell)
  else if (Limiter_Option < 0) then
    call jump_indicator_modified(nst,troubled_cell)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine jump_indicator
!
!###############################################################################
!
subroutine jump_indicator_original(nst,troubled_cell)
  !
  !> Compute the local marker using the jump indicator across adjacent cells.
  !> This is the original version without any modifications
  !>
  !> Refs:
  !>  -(2008) Barter,      PhD MIT,        http://hdl.handle.net/1721.1/57606
  !>  -(2007) Barter,      AIAA-2007-3823, doi:10.2514/6.2007-3823
  !> See Also:
  !>  -(2004) Krivodonova, Appl.Num.Math., doi:10.1016/j.apnum.2003.11.002
  !>  -(2013) Bhatia,      AIAA-2013-0302, doi:10.2514/6.2013-302
  !>  -(2014) Bhatia,      PhD U.Florida,
  !>                       https://apps.dtic.mil/sti/citations/tr/ADA624182
  !
  !.. Use Statements ..
  use geomvar_mod,      only : nr,ncell,cell,face
  use element_type_mod, only : std_elem
  use flowvar_mod,      only : facegpsp,faceusp
  use flowvar_mod,      only : get_marker_index
  use metrics_mod,      only : facenrm
  use ovar,             only : dpsi => marker_dpsi
  use ovar,             only : marker_c1,marker_c2
  use ovar,             only : Default_Marker_c1_RI,Default_Marker_c2_RI
  use ovar,             only : Default_Marker_c1_JI,Default_Marker_c2_JI
  use ovar,             only : ji_gwt,ji_uwt,ji_gvar,ji_uvar
  !
  !.. Optional Arguments ..
  integer,      optional, intent(in)    :: nst              !< Runge-Kutta stage
  logical(ldk), optional, intent(inout) :: troubled_cell(:)
  !
  !.. Local Scalars ..
  integer  :: nc,ncf,nf,k,nfp,f0,ierr
  integer  :: sf,sl,sr
  integer  :: face_geom,face_order
  integer  :: imarker
  real(wp) :: ds,rnx,rny,rnz,val,psi0
  real(wp) :: vall,valr,area
  real(wp) :: numer,denom,wt
  real(wp) :: c1,c2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "jump_indicator_original"
  !
  real(wp), parameter :: face_sign(1:2) = [one,-one]
  integer,  parameter :: other_side(1:2) = [2,1]
  !
continue
  !
  ! If the values contained in marker_c1 and marker_c2 are respectively
  ! the same as Default_Marker_c1_RI and Default_Marker_c2_RI, that means
  ! the marker values are currently set to the default values for the
  ! resolution indicator. If this is the case, switch these over to the
  ! default values for the jump indicator.
  !
  if (abs(marker_c1-Default_Marker_c1_RI) <= eps6 .and. &
      abs(marker_c2-Default_Marker_c2_RI) <= eps6) then
    c1 = Default_Marker_c1_JI
    c2 = Default_Marker_c2_JI
  else
    c1 = marker_c1
    c2 = marker_c2
  end if
  !
  if (.not. allocated(marker)) then
    allocate ( marker(1:ncell,1:nmarker) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"marker",1,__LINE__,__FILE__,ierr,error_message)
    !$ACC UPDATE DEVICE(marker)
  end if
  !
  if (.not. allocated(scaled_marker)) then
    allocate ( scaled_marker(1:ncell,1:nmarker) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"scaled_marker",1,__LINE__,__FILE__, &
                     ierr,error_message)
    !$ACC UPDATE DEVICE(scaled_marker)
  end if
  !
  ! Determine which marker position to fill based on Runge-Kutta stage.
  !
  imarker = get_marker_index(nst)
  !
  if (nr == 3) then
    !
    do nc = 1,ncell
      !
      val = zero
      area = zero
      !
      do ncf = 1,size(cell(nc)%face)
        !
        nf = cell(nc)%face(ncf)%idx
        sf = cell(nc)%face(ncf)%side
        !
        sl = sf
        sr = other_side(sl)
        !
        face_geom = face(nf)%geom
        face_order = face(nf)%order
        !
        f0 = face(nf)%beg_fp - 1
        nfp = face(nf)%end_fp - f0
        !
        do k = 1,nfp
          !
          rnx = facenrm(1,f0+k)*face_sign(sl)
          rny = facenrm(2,f0+k)*face_sign(sl)
          rnz = facenrm(3,f0+k)*face_sign(sl)
          ds = sqrt(rnx*rnx + rny*rny + rnz*rnz)
          !
          vall = ji_gwt *      facegpsp(f0+k,ji_gvar,sl) &
               + ji_uwt * faceusp(f0+k,ji_uvar,sl)
          valr = ji_gwt *      facegpsp(f0+k,ji_gvar,sr) &
               + ji_uwt * faceusp(f0+k,ji_uvar,sr)
          !
          numer = abs( vall - valr ) * ds
          denom =    ( vall + valr ) * half
          !
          wt = std_elem(face_geom,face_order)%wts(k)
          !
          val = val + wt*numer/denom
          area = area + wt*ds
          !
        end do
        !
      end do
      !
      marker(nc,imarker) = log10( max( val/area , rndoff ) )
      !
    end do
    !
  else if (nr == 2) then
    !
    do nc = 1,ncell
      !
      val = zero
      area = zero
      !
      do ncf = 1,size(cell(nc)%face)
        !
        nf = cell(nc)%face(ncf)%idx
        sf = cell(nc)%face(ncf)%side
        !
        sl = sf
        sr = other_side(sl)
        !
        face_geom = face(nf)%geom
        face_order = face(nf)%order
        !
        f0 = face(nf)%beg_fp - 1
        nfp = face(nf)%end_fp - f0
        !
        do k = 1,nfp
          !
          rnx = facenrm(1,f0+k)*face_sign(sl)
          rny = facenrm(2,f0+k)*face_sign(sl)
          ds = sqrt(rnx*rnx + rny*rny)
          !
          vall = ji_gwt *      facegpsp(f0+k,ji_gvar,sl) &
               + ji_uwt * faceusp(f0+k,ji_uvar,sl)
          valr = ji_gwt *      facegpsp(f0+k,ji_gvar,sr) &
               + ji_uwt * faceusp(f0+k,ji_uvar,sr)
          !
          numer = abs( vall - valr ) * ds
          denom =    ( vall + valr ) * half
          !
          wt = std_elem(face_geom,face_order)%wts(k)
          !
          val = val + wt*numer/denom
          area = area + wt*ds
          !
        end do
        !
      end do
      !
      marker(nc,imarker) = log10( max( val/area , rndoff ) )
      !
    end do
    !
  else
    !
    do nc = 1,ncell
      !
      val = zero
      area = zero
      !
      do ncf = 1,size(cell(nc)%face)
        !
        nf = cell(nc)%face(ncf)%idx
        sf = cell(nc)%face(ncf)%side
        !
        sl = sf
        sr = other_side(sl)
        !
        face_geom = face(nf)%geom
        face_order = face(nf)%order
        !
        f0 = face(nf)%beg_fp - 1
        nfp = face(nf)%end_fp - f0
        !
        do k = 1,nfp
          !
          rnx = facenrm(1,f0+k)*face_sign(sl)
          ds = sqrt(rnx*rnx)
          !
          vall = ji_gwt *      facegpsp(f0+k,ji_gvar,sl) &
               + ji_uwt * faceusp(f0+k,ji_uvar,sl)
          valr = ji_gwt *      facegpsp(f0+k,ji_gvar,sr) &
               + ji_uwt * faceusp(f0+k,ji_uvar,sr)
          !
          numer = abs( vall - valr ) * ds
          denom =    ( vall + valr ) * half
          !
          wt = std_elem(face_geom,face_order)%wts(k)
          !
          val = val + wt*numer/denom
          area = area + wt*ds
          !
        end do
        !
      end do
      !
      marker(nc,imarker) = log10( max( val/area , rndoff ) )
      !
    end do
    !
  end if
  !
  !$ACC PARALLEL LOOP PRESENT(cell,marker,scaled_marker)
  do nc = 1,ncell
    !
    psi0 = -one*(c1 + c2*log10(real(cell(nc)%order,kind=wp)))
    !
    if      (marker(nc,imarker) <= psi0-dpsi) then
      scaled_marker(nc,imarker) = zero
    else if (marker(nc,imarker) >= psi0+dpsi) then
      scaled_marker(nc,imarker) = one
    else
      scaled_marker(nc,imarker) = half*( one &
                                + sin(half*pi*(marker(nc,imarker)-psi0)/dpsi) )
    end if
    !
  end do
  !
  if (present(troubled_cell)) then
    !$ACC PARALLEL LOOP PRESENT(scaled_marker,troubled_cell)
    do nc = 1,ncell
      troubled_cell(nc) = (scaled_marker(nc,imarker) >= 1.0e-8_wp)
    end do
  end if
  !
end subroutine jump_indicator_original
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
subroutine jump_indicator_modified(nst,troubled_cell)
  !
  !> Compute the local marker using the jump indicator across adjacent cells.
  !> This is a modified version of the jump indicator.
  !>
  !> Refs:
  !>  -(2008) Barter,      PhD MIT,        http://hdl.handle.net/1721.1/57606
  !>  -(2007) Barter,      AIAA-2007-3823, doi:10.2514/6.2007-3823
  !> See Also:
  !>  -(2004) Krivodonova, Appl.Num.Math., doi:10.1016/j.apnum.2003.11.002
  !>  -(2013) Bhatia,      AIAA-2013-0302, doi:10.2514/6.2013-302
  !>  -(2014) Bhatia,      PhD U.Florida,
  !>                       https://apps.dtic.mil/sti/citations/tr/ADA624182
  !
  !.. Use Statements ..
  use geomvar_mod,      only : nr,ncell,cell,face,nfbnd,bface
  use element_type_mod, only : std_elem
  use flowvar_mod,      only : facegpsp,faceusp
  use flowvar_mod,      only : get_marker_index
  use metrics_mod,      only : facenrm
  use ovar,             only : dpsi => marker_dpsi
  use ovar,             only : marker_c1,marker_c2
  use ovar,             only : Default_Marker_c1_RI,Default_Marker_c2_RI
  use ovar,             only : Default_Marker_c1_JI,Default_Marker_c2_JI
  use ovar,             only : ji_gwt,ji_uwt,ji_gvar,ji_uvar
  use cell_ave_mod,     only : compute_cell_ave_var
  !
  !.. Optional Arguments ..
  integer,      optional, intent(in)    :: nst              !< Runge-Kutta stage
  logical(ldk), optional, intent(inout) :: troubled_cell(:)
  !
  !.. Local Scalars ..
  integer  :: nc,ncf,nf,k,nfp,f0,ierr
  integer  :: sf,sl,sr
  integer  :: face_geom,face_order
  integer  :: imarker
  real(wp) :: ds,rnx,rny,rnz,val,psi0
  real(wp) :: vall,valr,area,valave,pointmx
  real(wp) :: numer,denom,wt
  real(wp) :: c1,c2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "jump_indicator_modified"
  !
  real(wp), parameter :: face_sign(1:2) = [one,-one]
  integer,  parameter :: other_side(1:2) = [2,1]
  !
  real(wp), allocatable :: cell_ave_p(:)
continue
  !
  !
  ! If the values contained in marker_c1 and marker_c2 are respectively
  ! the same as Default_Marker_c1_RI and Default_Marker_c2_RI, that means
  ! the marker values are currently set to the default values for the
  ! resolution indicator. If this is the case, switch these over to the
  ! default values for the jump indicator.
  !
  if (abs(marker_c1-Default_Marker_c1_RI) <= eps6 .and. &
      abs(marker_c2-Default_Marker_c2_RI) <= eps6) then
    c1 = Default_Marker_c1_JI
    c2 = Default_Marker_c2_JI
  else
    c1 = marker_c1
    c2 = marker_c2
  end if
  !
  if (.not. allocated(marker)) then
    allocate ( marker(1:ncell,1:nmarker) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"marker",1,__LINE__,__FILE__,ierr,error_message)
    !$ACC UPDATE DEVICE(marker)
  end if
  !
  if (.not. allocated(scaled_marker)) then
    allocate ( scaled_marker(1:ncell,1:nmarker) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"scaled_marker",1,__LINE__,__FILE__, &
                     ierr,error_message)
    !$ACC UPDATE DEVICE(scaled_marker)
  end if
  !
  ! Determine which marker position to fill based on Runge-Kutta stage.
  !
  imarker = get_marker_index(nst)
  !
  cell_ave_p = init_array(ncell,zero) ! f2003 auto-realloc
  call compute_cell_ave_var(nq+1,cell_ave_p)
  !
  if (nr == 3) then
    !
    do nc = 1,ncell
      !
      val      = zero
      area     = zero
      pointmx  = rndoff
      !
      do ncf = 1,size(cell(nc)%face)
        !
        nf = cell(nc)%face(ncf)%idx
        sf = cell(nc)%face(ncf)%side
        !
        if (nf<=nfbnd) then 
          if (any(bface(1,nf) == bc_walls)) cycle
        end if
        !
        sl = sf
        sr = other_side(sl)
        !
        face_geom = face(nf)%geom
        face_order = face(nf)%order
        !
        f0 = face(nf)%beg_fp - 1
        nfp = face(nf)%end_fp - f0
        !
        do k = 1,nfp
                    !
          rnx = facenrm(1,f0+k)*face_sign(sl)
          rny = facenrm(2,f0+k)*face_sign(sl)
          rnz = facenrm(3,f0+k)*face_sign(sl)
          ds = sqrt(rnx*rnx + rny*rny + rnz*rnz)
          !
          vall = ji_gwt *      facegpsp(f0+k,ji_gvar,sl) &
               + ji_uwt * faceusp(f0+k,ji_uvar,sl)
          valr = ji_gwt *      facegpsp(f0+k,ji_gvar,sr) &
               + ji_uwt * faceusp(f0+k,ji_uvar,sr)
          !
          numer = abs( vall - valr ) * ds
          denom = ( vall + valr ) * half
          !
          wt = std_elem(face_geom,face_order)%wts(k)
          !
          !!! if (ji_gwt >= 0.99) then
            !!!Determine if filtering would help this jump of not (i.e. if slope helps or hurts the jump)
            valave = cell_ave_p(nc)
            !if ( (vall-valave) * (valr-vall) >= zero) then
            !        numer = zero
            !endif
            !!!! endif
          !
          val = val + wt*numer/denom
          area = area + wt*ds
          pointmx = max(numer/(ds*denom) , pointmx)
          !
        end do
        !
      end do
      !
      marker(nc,imarker) = log10(pointmx)
      ! marker(nc,imarker) = log10( max( val/area , rndoff ) )
      !
    end do
    !
  else if (nr == 2) then
    !
    do nc = 1,ncell
      !
      val      = zero
      area     = zero
      pointmx  = rndoff
      !
      do ncf = 1,size(cell(nc)%face)
        !
        nf = cell(nc)%face(ncf)%idx
        sf = cell(nc)%face(ncf)%side
        !
        if (nf<=nfbnd) then 
          if (any(bface(1,nf) == bc_walls)) cycle
        end if
        !
        sl = sf
        sr = other_side(sl)
        !
        face_geom = face(nf)%geom
        face_order = face(nf)%order
        !
        f0 = face(nf)%beg_fp - 1
        nfp = face(nf)%end_fp - f0
        !
        do k = 1,nfp
          !
          rnx = facenrm(1,f0+k)*face_sign(sl)
          rny = facenrm(2,f0+k)*face_sign(sl)
          ds = sqrt(rnx*rnx + rny*rny)
          !
          vall = ji_gwt *      facegpsp(f0+k,ji_gvar,sl) &
               + ji_uwt * faceusp(f0+k,ji_uvar,sl)
          valr = ji_gwt *      facegpsp(f0+k,ji_gvar,sr) &
               + ji_uwt * faceusp(f0+k,ji_uvar,sr)
          !
          numer = abs( vall - valr ) * ds
          denom = ( vall + valr ) * half
          !
          wt = std_elem(face_geom,face_order)%wts(k)
          !
          !!! if (ji_gwt >= 0.99) then
          !!!Determine if filtering would help this jump of not (i.e. if slope helps or hurts the jump)
            !if ((vall-valave) * (valr-vall) >= zero) then
            !   numer = zero
            !endif
          !!!! endif
          !
          val = val + wt*numer/denom
          area = area + wt*ds
          pointmx = max( numer/(ds*denom) , pointmx)
          !
        end do
        !
      end do
      !
      marker(nc,imarker) = log10(pointmx)
      !marker(nc,imarker) = log10( max( val/area , rndoff ) )
      !
    end do
    !
  else
    !
    do nc = 1,ncell
      !
      val      = zero
      area     = zero
      pointmx  = rndoff
      !
      do ncf = 1,size(cell(nc)%face)
        !
        nf = cell(nc)%face(ncf)%idx
        sf = cell(nc)%face(ncf)%side
        !
        if (nf<=nfbnd) then 
          if (any(bface(1,nf) == bc_walls)) cycle
        end if
        !
        sl = sf
        sr = other_side(sl)
        !
        face_geom = face(nf)%geom
        face_order = face(nf)%order
        !
        f0 = face(nf)%beg_fp - 1
        nfp = face(nf)%end_fp - f0
        !
        do k = 1,nfp
          !
          rnx = facenrm(1,f0+k)*face_sign(sl)
          ds = sqrt(rnx*rnx)
          !
          vall = ji_gwt *      facegpsp(f0+k,ji_gvar,sl) &
               + ji_uwt * faceusp(f0+k,ji_uvar,sl)
          valr = ji_gwt *      facegpsp(f0+k,ji_gvar,sr) &
               + ji_uwt * faceusp(f0+k,ji_uvar,sr)
          !
          numer = abs( vall - valr ) * ds
          denom = ( vall + valr ) * half
          !
          wt = std_elem(face_geom,face_order)%wts(k)
          !
          !!! if (ji_gwt >= 0.99) then
            !!!Determine if filtering would help this jump of not (i.e. if slope helps or hurts the jump)
            !if ((vall-valave) * (valr-vall) >= zero) then
            !        numer = zero
            !endif
          !!! endif
          !
          val = val + wt*numer/denom
          area = area + wt*ds
          pointmx = max( numer/(ds*denom) , rndoff)
          !
        end do
        !
      end do
      !
      !marker(nc,imarker) = log10(pointmx)
      marker(nc,imarker) = log10( max( val/area , rndoff ) )
      !
    end do
    !
  end if
  !
  !$ACC PARALLEL LOOP PRESENT(cell,marker,scaled_marker)
  do nc = 1,ncell
    !
    psi0 = -one*(c1 + c2*log10(real(cell(nc)%order,kind=wp)))
    !
    if      (marker(nc,imarker) <= psi0-dpsi) then
      scaled_marker(nc,imarker) = zero
    else if (marker(nc,imarker) >= psi0+dpsi) then
      scaled_marker(nc,imarker) = one
    else
      scaled_marker(nc,imarker) = half*( one &
                                + sin(half*pi*(marker(nc,imarker)-psi0)/dpsi) )
    end if
    !
  end do
  !
  if (present(troubled_cell)) then
    !$ACC PARALLEL LOOP PRESENT(scaled_marker,troubled_cell)
    do nc = 1,ncell
      troubled_cell(nc) = (scaled_marker(nc,imarker) >= 1.0e-8_wp)
    end do
  end if
  !
end subroutine jump_indicator_modified
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
subroutine ducros_indicator(nst,troubled_cell)
  !
  !> Compute the local marker using a modified ducross shock sensor.
  !>
  !> Refs:                              ##### Reminder to find links for these #####
  !>  - Caleb's or other paper
  !> See Also:
  !>  - OG Ducros and common mods
  !>
  !
  !.. Use Statements ..


  use geomvar_mod,      only : nr,ncell,cell,face,nfbnd,bface,n_solpts,cell_order
  use element_type_mod, only : std_elem
  use eqn_idx_mod,      only : nmx,nmy,nmz,ntemp, nq
  use ovar,             only : gam,rgasref,machref,tref,pref,rhoref,aref,velref,Lref
  use gpu_mod,          only : grad_cv_to_grad_pv_gpu
  use flowvar_mod,      only : usp,gpsp,dusp
  use flowvar_mod,      only : get_marker_index
  use flux_mod,         only : discontinuous_gradient_at_solpts
  !
  !.. Optional Arguments ..
  integer,      optional, intent(in)    :: nst              !< Runge-Kutta stage
  logical(ldk), optional, intent(inout) :: troubled_cell(:)
  !
  !.. Local Scalars ..
  real(wp)      :: ducros_fs,eps_factor,char_len,val,cell_len,eps_d, &
                   dudx,dvdy,dwdz,divv,vortmags,f_duc,asp,velrefmag,S,valS, &
                   trm,dist,deld,musp,mumax,muscale,cbeta,cmax
  integer       :: ierr,imarker,nc,nsp,n1,n2
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: dpvdx(:,:,:), cell_lengths(:)
  real(wp) :: dvdxref
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "ducros_indicator"
  !
  !.. Constant Model Parameters ..
  ducros_fs  = 0.01_wp
  eps_factor = 1.0e-8_wp
  dist = 1e-4_wp
  deld = 0.2_wp * dist
  cbeta = 1
  cmax = 0.5
  !
continue
  !
  call jump_indicator(nst,troubled_cell) 
  !
  call discontinuous_gradient_at_solpts()
  !
  allocate(cell_lengths(ncell), source=zero, stat=ierr , errmsg=error_message)
  call alloc_error(pname,"cell_lengths",1,__LINE__,__FILE__,ierr,error_message)
  !
  if (.not. allocated(dpvdx)) then
     allocate(dpvdx(1:n_solpts,1:nr,1:nq), source=zero, stat=ierr , errmsg=error_message) 
     call alloc_error(pname,"dpvdx",1,__LINE__,__FILE__,ierr,error_message)
  endif
  !
  if (.not. allocated(marker)) then
    allocate ( marker(1:ncell,1:nmarker) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"marker",1,__LINE__,__FILE__,ierr,error_message)
    !$ACC UPDATE DEVICE(marker)
  end if
  !
  if (.not. allocated(scaled_marker)) then
    allocate ( scaled_marker(1:ncell,1:nmarker) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"scaled_marker",1,__LINE__,__FILE__, &
                     ierr,error_message)
    !$ACC UPDATE DEVICE(scaled_marker)
  end if
  !
  !  Determine which marker position to fill based on R-K stage
  imarker = get_marker_index(nst)
  !
  !Refference velocity magnitude
  !
  n1 = cell(1)%beg_sp
  n2 = cell(ncell)%end_sp
  call grad_cv_to_grad_pv_gpu(n1,n2,usp,gpsp,dusp,dpvdx)
  !
  call cell_length_scale_max_box(cell_lengths)
  !
  ! test outputs
  ! write(iout,*) 'velref(1,2,3): ', velref(1), velref(2), velref(3)
  ! write(iout,*) 'Lref: ', Lref
  ! write(iout,*) 'cell(1)%beg_sp, cell(nc)%end_sp, n_solpts', n1, n2, n_solpts
  ! write(iout,*) 'nmx, nmy: ', nmx, nmy
  ! write (iout,*) 'sum(dusp)', sum(abs(dusp))
  !
  do nc = 1,ncell
     !
     val = zero
     valS = zero
     !
     ! calculate epsilon normalization factor
     eps_d = (eps_factor * aref * machref )/ cell_lengths(nc)  
     !
     ! Now Loop Through Solution Points
     do nsp = cell(nc)%beg_sp, cell(nc)%end_sp
       !
       ! Calculate Divergence squared
       dudx = zero
       dvdy = zero
       dwdz = zero
       !
       dudx = aref * dpvdx(nsp,1,nmx)
       !
       if (nr >= 2) then
       dvdy = aref * dpvdx(nsp,2,nmy)
       end if
       !
       if (nr == 3) then
       dwdz = aref * dpvdx(nsp,3,nmz)
       end if
       !
       divv = dudx + dvdy
       !
       if (nr ==3) then
       divv = divv + dwdz
       endif
       !
       divv = sqrt(divv*divv)
       !
       ! Calculate Vorticity Magnitude Squared
       vortmags = zero
       !
       if (nr >= 2) then
       vortmags = vortmags + (aref*dpvdx(nsp,2,nmx) - aref*dpvdx(nsp,1,nmy))**2
       endif
       !
       if (nr == 3) then
       vortmags = vortmags + (aref*dpvdx(nsp,3,nmy) - aref*dpvdx(nsp,2,nmz))**2
       vortmags = vortmags + (aref*dpvdx(nsp,3,nmx) - aref*dpvdx(nsp,1,nmz))**2
       end if
       !
       ! Calculate the ducros sensor
       f_duc = (divv + eps_d) / (divv + vortmags + eps_d)
       f_duc = min( f_duc, one)
       !
       ! Calculate Sound Speed (modified to be tolerant to non-physical negative vals)
       asp = aref * sqrt(gam*max(gpsp(nsp,npres)/usp(nsp,nec), 1e-8_wp))
       !
       ! (don't) Calculate trm (requires wall distance: ydc)
       !trm = half * (one - tanh((ydc-dist)/deld))
       trm = zero
       !
       ! Calculate S quantity
       S = min((divv * cell_lengths(nc) * half) / (asp) , one) 
       S = S * S
       !
       val = max( val, min( one, max( S, trm, min(ducros_fs, f_duc))))
       valS = max( valS, (divv * cell_lengths(nc) * half) / (asp))
     end do
     !
     ! For now the same scaling as the version from REACTMB
     marker(nc,imarker) = valS
     scaled_marker(nc,imarker) = scaled_marker(nc,imarker)*val !!! can swap mult for max for very heavy marking
     !
  end do
  !
  if (present(troubled_cell)) then
    !$ACC PARALLEL LOOP PRESENT(scaled_marker,troubled_cell)
    do nc = 1,ncell
      troubled_cell(nc) = (scaled_marker(nc,imarker) >= 1.0e-8_wp)
    end do
  end if
  !
  deallocate(cell_lengths)
  !deallocate(dpvdx)
  !
end subroutine ducros_indicator
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
subroutine cell_length_scale_max_box(cell_lengths)
  !
  ! Finds the longest dimension of the x,y,z bounding box of the cell
  !
  use geomvar_mod, only : ncell,nodes_of_cell_ptr,nodes_of_cell,xyz_nodes, nr
  !
  !.. I/O Arrays ..
  real(wp),intent(out)     :: cell_lengths(:)
  !
  !.. Local Scalars ..
  real(wp)      :: xmin, xmax, ymin, ymax, zmin, zmax
  integer       :: ibeg,iend,jj,j,i,nc,nn,ii,maxnodes
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable      :: cell_nodes_xyz(:,:)
  !
  maxnodes = 8 !?
  !
  allocate(cell_nodes_xyz(nr,maxnodes))
  !
  do nc = 1,ncell
    jj = 0
    ibeg = nodes_of_cell_ptr(nc) + 1
    iend = nodes_of_cell_ptr(nc+1)
    do i = 1, (iend-ibeg)
      ii  = i + ibeg - 1
      nn = nodes_of_cell(ii)
      jj = jj+1
      cell_nodes_xyz(1:nr,jj) = xyz_nodes(1:nr,nn) 
    enddo
    !
    !  Initialize with coords of 1st node
    !
    xmax = cell_nodes_xyz(1,1)
    xmin = cell_nodes_xyz(1,1)
    !
    ymax = cell_nodes_xyz(2,1)
    ymin = cell_nodes_xyz(2,1)
    !
    if (nr ==3 ) then
      zmin = cell_nodes_xyz(3,1)
      zmax = cell_nodes_xyz(3,1)
    endif
    ! loop through nodes to find extrema
    do j = 2,jj
      !
      xmax = max(xmax, cell_nodes_xyz(1,j))
      xmin = min(xmin, cell_nodes_xyz(1,j))
      !
      ymax = max(ymax, cell_nodes_xyz(2,j))
      ymin = min(ymin, cell_nodes_xyz(2,j))
      !
      if (nr ==3 ) then
        zmax = max(zmax, cell_nodes_xyz(3,j))
        zmin = min(zmin, cell_nodes_xyz(3,j))
      endif
      !
    enddo
    !
    ! Find the maximum distance between dx,dy,dz
      if (nr >=3 ) then
        cell_lengths(nc) = max( xmax-xmin , ymax-ymin , zmax-zmin )
      else
        cell_lengths(nc) = max( xmax-xmin , ymax-ymin)
      endif
    !
  enddo
    deallocate(cell_nodes_xyz)
end subroutine cell_length_scale_max_box
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
  1 format (/,"Status of arrays in module limiters_mod")
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
end module limiters_mod
