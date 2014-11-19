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
  public :: initialize_limiters
  public :: compute_cell_ave
  !
  logical, save :: using_markers = fals
  logical, save :: using_moment_limiter = fals
  logical, save :: init_all_as_marked = true
  !
  real(wp), save, allocatable, dimension(:,:) :: cell_ave_sol
  real(wp), save, allocatable, dimension(:,:) :: umodes
  !
  logical, save, allocatable, dimension(:) :: troubled_cell
  logical, save, allocatable, dimension(:) :: needs_limiting
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
continue
  !
  if (using_markers) then
    !
    ! Create the array for storing the cell averaged solutions
    !
    allocate ( cell_ave_sol(1:nq,1:ncell) , source=zero , stat=ierr )
    call alloc_error("initialize_limiters","cell_ave_sol",1,ierr)
    !
    !
    ! Create the array for identifying troubled cells
    !
    allocate ( troubled_cell(1:ncell) , source=fals , stat=ierr )
    call alloc_error("initialize_limiters","troubled_cell",1,ierr)
    !
  end if
  !
  ! Create the array for storing the moments of each cell
  !
  if (using_moment_limiter) then
    !
    allocate ( umodes(1:n_totpts,1:nq) , source=zero , stat=ierr )
    call alloc_error("initialize_limiters","umodes",1,ierr)
    !
    allocate ( needs_limiting(1:ncell) , stat=ierr )
    call alloc_error("initialize_limiters","needs_limiting",1,ierr)
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
  use geovar, only : vand_inv
  use geovar, only : n_sp2dt,n_sp2dq
  use flowvar, only : usp
  !
  !.. Local Scalars ..
  integer :: i,j,nc,n1,n2,np,nw,m
  !
  !.. Local Arrays ..
  real(wp), dimension(1:n_sp2dq) :: inv_vand
  !
continue
  !
  do nc = 1,ncell
    !
    n1 = pts_in_cell_ptr(nc)+1
    n2 = pts_in_cell_ptr(nc+1)
    np = n2-n1+1
    nw = merge(1,2,np==n_sp2dt)
    !
    inv_vand(1:np) = vand_inv(1,1:np,nw)
    !
    do m = 1,nq
      cell_ave_sol(m,nc) = half * dot( inv_vand(1:np) , usp(m,n1:n2) )
    end do
    !
  end do
  !
  call stop_recon(stop_mpi,"compute_cell_ave")
  !
end subroutine compute_cell_ave
!
!###############################################################################
!
subroutine moment_limiter
  !
  !.. Use Statements ..
  use geovar, only : ncell,n_sp2dt
  use geovar, only : pts_in_cell_ptr,vand_inv
  use flowvar, only : usp
  !
  use parallel_mod, only : exchange_modes
  !
  !.. Local Scalars
  integer :: nc,n1,n2,np,nw,m
  !
continue
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
    nw = merge(1,2,np==n_sp2dt)
    !
    do m = 1,nq
      umodes(n1:n2,m) = matmul( vand_inv(1:np,1:np,nw) , usp(m,n1:n2) )
    end do
    !
  end do
  !
  ! Exchange modes of cells on partition boundaries
  !
  if (ncpu > 1) call exchange_modes(umodes)
  !
  
  !
end subroutine moment_limiter
!
!###############################################################################
!
subroutine limit_by_projecting
  !
  !.. Use Statements ..
  use geovar,  only : ncell,n_order,n_sp2dt,n_sp2dq
  use geovar,  only : pts_in_cell_ptr,projct_limiter
  use flowvar, only : usp
  !
  !.. Local Scalars ..
  integer :: k,m,n,n1,n2,nc,nsp,nw
  !
  !.. Local Arrays ..
  real(wp), dimension(1:n_sp2dq,1:nq) :: tusp
  !
continue
  !
  ! Loop from the highest order to the lowest order and project
  ! any marked cells to one degree lower order until there are
  ! no more marked cells or we reach the lowest order.
  ! NOTE: The projection operators are stored in the array projct_limiter such
  !       that when looping from n = 1 -> n_order, the array gives:
  !        projct_limiter(:,:,:,n) => projects cell solution to (n_order-n) order
  !        i.e., projct_limiter(:,:,:,1)         => projects cell solution to
  !                                                 (n_order-1) order
  !              projct_limiter(:,:,:,n_order-1) => projects cell solution to
  !                                                 a linear solution (P1 order)
  !              projct_limiter(:,:,:,n_order)   => projects cell solution to
  !                                                 a constant solution (P0 order)
  !
  moment_loop: do n = 1,n_order
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
      nw = merge(1,2,nsp==n_sp2dt)
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
      do k = 1,nsp
        do m = 1,nq
          usp(m,n1-1+k) = dot( projct_limiter(1:nsp,k,nw,n) , &
                                       tusp(1:nsp,m) )
        end do
      end do
      !
    end do cell_loop
    !
    ! Recheck the troubled cells to see if they are no longer
    ! troubled due to the previous projection.
    !
    call recheck_troubled_cells
    !
    ! If there are no more troubled cells, we no longer need to do
    ! any projections so exit this outer moment loop.
    !
    if (.not. any(troubled_cell)) exit moment_loop
    ! 
  end do moment_loop
  !
  !
  !
end subroutine limit_by_projecting
!
!###############################################################################
!
subroutine recheck_troubled_cells
  !
  !.. Use Statements ..
  use geovar,  only : ncell,n_sp2dt,pts_in_cell_ptr
  use flowvar, only : usp
  !
  !.. Local Scalars ..
  integer :: nc,n1,n2,nsp,nw
  !
continue
  !
  do nc = 1,ncell
    !
    n1 = pts_in_cell_ptr(nc)+1
    n2 = pts_in_cell_ptr(nc+1)
    nsp = n2-n1+1
    nw = merge(1,2,nsp==n_sp2dt)
    !
  end do
  !
end subroutine recheck_troubled_cells
!
!###############################################################################
!
end module module_limiters
