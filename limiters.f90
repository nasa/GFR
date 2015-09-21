module module_limiters
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nq
  !
  implicit none
  !
  private
  !
  !
  !.. Public Procedures ..
  !
  public :: initialize_limiters
  public :: compute_cell_ave
  public :: limit_by_projecting
  public :: limiters_memory_usage
  !
  !
  !.. Module Local Scalars ..
  !
  logical, save :: using_markers = fals
  logical, save :: using_moment_limiter = fals
  logical, save :: init_all_as_marked = true
  !
  !
  !.. Module Local Allocatable Arrays ..
  !
  real(wp), save, allocatable, dimension(:,:) :: cell_ave_sol
  real(wp), save, allocatable, dimension(:,:) :: umodes
  !
  logical, save, allocatable, dimension(:) :: troubled_cell
  logical, save, allocatable, dimension(:) :: needs_limiting
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
  if (using_markers) then
    !
    ! Create the array for storing the cell averaged solutions
    !
    allocate ( cell_ave_sol(1:nq,1:ncell) , source=zero , stat=ierr )
    call alloc_error(pname,"cell_ave_sol",1,ierr)
    !
    !
    ! Create the array for identifying troubled cells
    !
    allocate ( troubled_cell(1:ncell) , source=fals , stat=ierr )
    call alloc_error(pname,"troubled_cell",1,ierr)
    !
  end if
  !
  ! Create the array for storing the moments of each cell
  !
  if (using_moment_limiter) then
    !
    allocate ( umodes(1:n_totpts,1:nq) , source=zero , stat=ierr )
    call alloc_error(pname,"umodes",1,ierr)
    !
    allocate ( needs_limiting(1:ncell) , stat=ierr )
    call alloc_error(pname,"needs_limiting",1,ierr)
    !
  end if
  !
end subroutine initialize_limiters
!
!###############################################################################
!
subroutine compute_cell_ave
  !
  !.. Use Statements ..
  use geovar, only : ncell,pts_in_cell_ptr
  use geovar, only : cell_geom,cell_order
  use order_mod, only : maxpts
  use flowvar, only : usp
  !
  use vandermonde_mod, only : vand
  !
  !.. Local Scalars ..
  integer :: nc,n1,n2,np,m
  integer :: this_geom,this_order
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxpts) :: const_mode
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_cell_ave"
  !
continue
  !
  do nc = 1,ncell
    !
    n1 = pts_in_cell_ptr(nc)+1
    n2 = pts_in_cell_ptr(nc+1)
    np = n2-n1+1
    !
    this_geom = cell_geom(nc)
    this_order = cell_order(nc)
    !
    const_mode(1:np) = vand(this_geom,this_order)%nodal2modal%mat(1,1:np)
    !
    do m = 1,nq
      cell_ave_sol(m,nc) = half * dot( const_mode(1:np) , usp(m,n1:n2) )
    end do
    !
  end do
  !
  call stop_recon(stop_mpi,pname)
  !
end subroutine compute_cell_ave
!
!###############################################################################
!
subroutine moment_limiter
  !
  !.. Use Statements ..
  use geovar, only : ncell,cell_geom,cell_order
  use geovar, only : pts_in_cell_ptr
  use flowvar, only : usp
  !
  use vandermonde_mod, only : vand
  !
 !use parallel_mod, only : exchange_modes
  !
  !.. Local Scalars
  integer :: nc,n1,n2,np,m
  integer :: this_geom,this_order
  !
#ifndef USE_MATRIX_PROCEDURES
  !.. Local Pointers ..
  real(wp), pointer, contiguous :: nodal2modal(:,:)
#endif
  !
continue
  !
#ifndef USE_MATRIX_PROCEDURES
  nodal2modal => null()
#endif
  !
  if (init_all_as_marked) needs_limiting = true
  !
  ! Get the solution modes in each cell of this partition
  !
  do nc = 1,ncell
    !
    n1 = pts_in_cell_ptr(nc)+1
    n2 = pts_in_cell_ptr(nc+1)
    np = n2-n1+1
    !
    this_geom = cell_geom(nc)
    this_order = cell_order(nc)
    !
#ifdef USE_MATRIX_PROCEDURES
    do m = 1,nq
      umodes(n1:n2,m) = &
        vand(this_geom,this_order)%nodal2modal%mv_mult( usp(m,n1:n2) )
    end do
#else
    nodal2modal => vand(this_geom,this_order)%nodal2modal%mat
    !
    do m = 1,nq
      umodes(n1:n2,m) = matmul( nodal2modal , usp(m,n1:n2) )
    end do
    !
    if (associated(nodal2modal)) nodal2modal => null()
#endif
    !
  end do
  !
  ! Exchange modes of cells on partition boundaries
  !
 !if (ncpu > 1) call exchange_modes(umodes)
  !
end subroutine moment_limiter
!
!###############################################################################
!
subroutine limit_by_projecting
  !
  !.. Use Statements ..
  use order_mod,  only : maxpts,n_max_order
  use geovar,     only : ncell,pts_in_cell_ptr
  use geovar,     only : cell_geom,cell_order
  use flowvar,    only : usp
  use projection_mod, only : projct
  !
  !.. Local Scalars ..
  integer :: k,m,n1,n2,nc,nsp
  integer :: this_geom,this_order
  integer :: dropped_orders,proj_order
  logical :: there_are_no_troubles
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxpts,1:nq) :: tusp
  !
#ifndef USE_MATRIX_PROCEDURES
  !.. Local Pointers ..
  real(wp), pointer, contiguous :: plim(:,:)
#endif
  !
continue
  !
  ! Loop from the highest order to the lowest order and project
  ! any marked cells to one degree lower order until there are
  ! no more marked cells or we reach the lowest order.
  !
  ! Initialize the local pointers to disassociated
  !
#ifndef USE_MATRIX_PROCEDURES
  plim => null()
#endif
  !
  ! Initialize the array troubled_cell to identify all troubled cells before
  ! attempting to use the projection limiter.
  !
  call identify_troubled_cells(init_troubled_cells,there_are_no_troubles)
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
    cell_loop: do nc = 1,ncell
      !
      ! If the current cell is not troubled, skip to the next one
      !
      if (.not.troubled_cell(nc)) cycle cell_loop
      !
      n1 = pts_in_cell_ptr(nc)+1
      n2 = pts_in_cell_ptr(nc+1)
      nsp = n2-n1+1
      !
      this_geom = cell_geom(nc)
      this_order = cell_order(nc)
      !
      proj_order = this_order - dropped_orders
      !
      ! Store a copy of usp for this cell because the individual
      ! solutions at the solution points within the cell will be
      ! altered during the projection computation. While we are
      ! doing this, store the copy as its transpose to reduce
      ! cache misses.
      !
      do k = 1,nsp
        do m = 1,nq
          tusp(k,m) = usp(m,n1-1+k)
        end do
      end do
      !
      ! Apply the limiter by multplying the projection matrix with
      ! the solution in the current cell
      !
#ifdef USE_MATRIX_PROCEDURES
      do m = 1,nq
        usp(m,n1:n2) = &
          projct(this_geom,this_order)%toCell(proj_order)%mv_mult( tusp(1:nsp,m) )
      end do
#else
      !
      ! Assign local pointer to projection matrix as an alias to simplify the code
      ! 
      plim => projct(this_geom,this_order)%toCell(proj_order)%mat
      !
      do m = 1,nq
        usp(m,n1:n2) = matmul( plim , tusp(1:nsp,m) )
      end do
      !
      ! Disassociate plim before continuing to the next cell
      !
      if (associated(plim)) plim => null()
#endif
      !
    end do cell_loop
    !
    ! Recheck the troubled cells to see if there are no longer
    ! troubled due to the previous projection.
    !
    call identify_troubled_cells(dropped_orders,there_are_no_troubles)
    !
    ! Exit this loop if there are no more troubled cells remaining
    !
    if (there_are_no_troubles) exit order_loop
    !
  end do order_loop
  !
end subroutine limit_by_projecting
!
!###############################################################################
!
subroutine identify_troubled_cells(dropped_orders,there_are_no_troubles)
  !
  !.. Use Statements ..
  use geovar,  only : ncell,pts_in_cell_ptr
  use geovar,  only : cell_geom,cell_order
  use flowvar, only : usp
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: dropped_orders
  logical, intent(out) :: there_are_no_troubles
  !
  !.. Local Scalars ..
  integer :: nc,n1,n2,proj_order
  integer :: this_geom,this_order
  !
continue
  !
  proj_order = dropped_orders
  !
  ! Initialize there_are_no_troubles to true meaning that no troubled cells
  ! exist. If one is found to exist in the following loop, this will be
  ! changed to false.
  !
  there_are_no_troubles = true
  !
  ! Check if troubled_cell needs to be initialized. If it does, mark
  ! all cells as troubled before actually checking the cells. Also set
  ! proj_order to zero so there are no problems when checking if the
  ! cell solution has been projected down to a first order solution.
  !
  if (dropped_orders == init_troubled_cells) then
    troubled_cell = true
    proj_order = 0
  end if
  !
  ! Now loop through the 
  !
  do nc = 1,ncell
    !
    ! If the current cell is not troubled, skip to the next one
    !
    if (.not. troubled_cell(nc)) cycle
    !
    n1 = pts_in_cell_ptr(nc)+1
    n2 = pts_in_cell_ptr(nc+1)
    !
    this_geom = cell_geom(nc)
    this_order = cell_order(nc)
    !
    ! If the projection has reduced the solution to first order, mark
    ! this cell as no longer troubled and skip to the next cell
    !
    if (this_order-proj_order <= 0) then
      troubled_cell(nc) = fals
      cycle
    end if
    !
    ! Recheck if this cell is still troubled after the previous projection
    !
    ! ###########################################################
    ! ###########################################################
    ! ###########################################################
    ! #####   NEED TO FIND WAY TO IDENTIFY TROUBLED CELLS   #####
    ! ###########################################################
    ! ###########################################################
    ! ###########################################################
    !
    ! If the current cell is determined to be troubled, change
    ! there_are_no_troubles to fals meaning that troubles cells exist.
    !
    if (troubled_cell(nc)) there_are_no_troubles = fals
    !
  end do
  !
end subroutine identify_troubled_cells
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
  if (allocated(needs_limiting)) then
    call array%set( storage_size(needs_limiting,kind=inttype) , &
                            size(needs_limiting,kind=inttype) )
    call total%add(array)
    write (array_name,5) "needs_limiting"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "needs_limiting"
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
end subroutine limiters_memory_usage
!
!###############################################################################
!
end module module_limiters
