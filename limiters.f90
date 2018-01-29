module module_limiters
  !
  !.. Use Statements ..
  use module_kind_types
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
  public :: compute_cell_ave
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
  logical(lk), save, allocatable :: troubled_cell(:)
  logical(lk), save, allocatable :: needs_limiting(:)
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
    allocate ( troubled_cell(1:ncell) , source=fals , &
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
 !  allocate ( needs_limiting(1:ncell) , source=fals , &
 !             stat=ierr , errmsg=error_message )
 !  call alloc_error(pname,"needs_limiting",1,__LINE__,__FILE__,ierr, &
 !                   error_message)
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
  use geovar, only : ncell,cell
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
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    np = n2-n1+1
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    const_mode(1:np) = vand(this_geom,this_order)%nodal2modal%mat(1,1:np)
    !
    do m = 1,nq
      cell_ave_sol(m,nc) = half * dot( const_mode(1:np) , usp(m,n1:n2) )
    end do
    !
  end do
  !
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  !
end subroutine compute_cell_ave
!
!###############################################################################
!
subroutine moment_limiter
  !
  !.. Use Statements ..
  use geovar, only : ncell,cell
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
continue
  !
  if (init_all_as_marked) needs_limiting = true
  !
  ! Get the solution modes in each cell of this partition
  !
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    np = n2-n1+1
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    do m = 1,nq
      umodes(n1:n2,m) = &
        vand(this_geom,this_order)%nodal2modal%mv_mult( usp(m,n1:n2) )
    end do
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
    cell_loop: do nc = 1,ncell
      !
      ! If the current cell is not troubled, skip to the next one
      !
      if (.not.troubled_cell(nc)) cycle cell_loop
      !
      n1 = cell(nc)%beg_sp
      n2 = cell(nc)%end_sp
      np = n2-n1+1
      !
      this_geom = cell(nc)%geom
      this_order = cell(nc)%order
      !
      ! Logic in subroutine identify_troubled_cells
      ! prevents this from becoming less than 0
      !
      proj_order = this_order - dropped_orders
      !
      ! Store a copy of usp for this cell because the individual
      ! solutions at the solution points within the cell will be
      ! altered during the projection computation. While we are
      ! doing this, store the copy as its transpose to reduce
      ! cache misses.
      !
      do k = 1,np
        do m = 1,nq
          tusp(k,m) = usp(m,n1-1+k)
        end do
      end do
      !
      ! Apply the limiter by multplying the projection matrix with
      ! the solution in the current cell
      !
      if (Limiter_Filter_Option /= 0) then
        do m = 1,nq
          usp(m,n1:n2) = filter(this_geom,this_order)% &
                                              limiter% &
                         mv_mult( tusp(1:np,m) )
        end do
        troubled_cell(nc) = fals ! only filter once
      else
        do m = 1,nq
          usp(m,n1:n2) = projct(this_geom,this_order)% &
                                  limiter(proj_order)% &
                         vm_mult( tusp(1:np,m) )
        end do
      end if
      !
      ! Unmark this cell as troubled if proj_order is 0
      !
      if (proj_order == 0) troubled_cell(nc) = fals
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
subroutine resolution_indicator(apply_projection,troubled_cell,orders_to_cutoff)
  !
  !.. Use Statements ..
  use geovar, only : ncell,cell
  use order_mod, only : maxSP,p_order,n_order
  use quadrature_mod, only : std_elem
  use flowvar, only : usp
  use metrics_mod, only : metrics_dt
  use filter_mod, only : filter
  use projection_mod, only : projct
  !
  !.. Optional Arguments ..
  logical(lk), optional,    intent(in) :: apply_projection
  logical(lk), optional, intent(inout) :: troubled_cell(:)
  integer,     optional,    intent(in) :: orders_to_cutoff
  !
  !.. Local Scalars ..
  integer  :: k,m,n,nc,n1,n2,np,ierr
  integer  :: this_geom,this_order
  integer  :: ord_mlt,ord_add
  integer  :: cutoff_order,lmtr_order
  real(wp) :: l2_p,l2_pm1,psi0,fk
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
  character(len=*), parameter :: pname = "resolution_indicator"
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
  end if
  !
  if (.not. allocated(scaled_marker)) then
    allocate ( scaled_marker(1:ncell) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"scaled_marker",1,__LINE__,__FILE__, &
                     ierr,error_message)
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
    rho_p(1:np) = usp(nec,n1:n2)
    !
    ! Compute u_hat, the degree-(P-1) polynomial solution
    ! of density at the solution points
    !
    rho_pm1(1:np) = filter(this_geom,this_order)% &
                            cutoff(cutoff_order)% &
                    mv_mult( rho_p(1:np) )
    !
    ! Compute u - u_hat
    !
    rho_pm1(1:np) = rho_p(1:np) - rho_pm1(1:np)
    !
    ! Compute the L2 inner product of the two polynomials
    !
    l2_pm1 = sum( rho_pm1(1:np) * rho_pm1(1:np) * metrics_dt(nc)%jac(1:np) * &
                  std_elem(this_geom,this_order)%wts(1:np) )
    l2_pm1 = max(l2_pm1,rndoff)
    !
    l2_p = sum( rho_p(1:np) * rho_p(1:np) * metrics_dt(nc)%jac(1:np) * &
                std_elem(this_geom,this_order)%wts(1:np) )
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
        lmtr_order = int( real(this_order,kind=wp)*(one-scaled_marker(nc)) )
        lmtr_order = max(1,min(this_order,lmtr_order))
        !
        if (lmtr_order < this_order) then
          !
          nlim(lmtr_order) = nlim(lmtr_order) + 1
          !
          ! Store the transpose of usp for this cell to reduce cache misses
          !
          do k = 1,np
            do m = 1,nq
              tusp(k,m) = usp(m,n1-1+k)
            end do
          end do
          !
          ! Project the solution down to proj_order
          ! at each solution point of this cell
          !
          do m = 1,nq
            usp(m,n1:n2) = projct(this_geom,this_order)% &
                                    limiter(lmtr_order)% &
                           vm_mult( tusp(1:np,m) )
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
end subroutine resolution_indicator
!
!###############################################################################
!
subroutine jump_indicator(troubled_cell)
  !
  !.. Use Statements ..
  use geovar, only : nr,ncell,cell,nface,face
  use order_mod, only : maxFP,geom_solpts
  use quadrature_mod, only : std_elem
  use flowvar, only : faceusp
  use metrics_mod, only : geofa
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(inout) :: troubled_cell(:)
  !
  !.. Local Scalars ..
  integer :: nc,nf,k,nfp,ierr
  integer :: left_geom,right_geom,face_geom
  integer :: left_cell,right_cell,face_order
  real(wp) :: ds,dsl,dsr,p_ave,p_diff,p_jump,psi0
  !
  !.. Local Arrays ..
  real(wp) :: unrm(1:nr)
  real(wp) :: pl_fp(1:maxFP)
  real(wp) :: pr_fp(1:maxFP)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: area_sum(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "jump_indicator"
  real(wp), parameter :: dpsi = half
  real(wp), parameter :: c1 = 2.25_wp
 !real(wp), parameter :: c1 = three
  real(wp), parameter :: c2 = three
 !real(wp), parameter :: c2 = three
  !
continue
  !
  if (.not. allocated(marker)) then
    allocate ( marker(1:ncell) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"marker",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (.not. allocated(scaled_marker)) then
    allocate ( scaled_marker(1:ncell) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"scaled_marker",1,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
 !allocate ( area_sum(1:ncell) , source=zero , &
 !           stat=ierr , errmsg=error_message )
 !call alloc_error(pname,"area_sum",1,__LINE__,__FILE__,ierr,error_message)
  !
  marker(:) = zero
  !
  do nf = 1,nface
    !
    face_geom = face(nf)%geom  ! geometry of face
    face_order = face(nf)%order ! order of face
    !
    nfp = geom_solpts(face_geom,face_order) ! number of flux points on this face
    !
    left_cell = face(nf)%left%cell ! left  cell on face
    right_cell = face(nf)%right%cell ! right cell on face
    !
    left_geom = cell(left_cell)%geom
    if (right_cell <= ncell) then
      right_geom = cell(right_cell)%geom
    else
      right_geom = face_geom
    end if
    !
    ! For all the flux points on the current face, extract the conservative
    ! conservative variables for the left and right cells from faceusp
    !
    pl_fp(1:nfp) = pressure_cv( faceusp(nf)%v(1:nq,1:nfp,1) )
    pr_fp(1:nfp) = pressure_cv( faceusp(nf)%v(1:nq,1:nfp,2) )
    !
    ! Correct the gradients of the conservative variables using the
    ! solution jumps at the interface and the correction function
    !
    ! Use the LDG flux formulation for computing the common interface solutions
    !   u_com = half*(u_left + u_right) - &
    !                (u_left - u_right)*dot(LDG_beta,unit_normal)
    !
    do k = 1,nfp
     !!
     !! Get the unit normal at this flux point
     !!
     !unrm(1:nr) = unit_vector( geofa(nf)%v(1:nr,k) )
      !
      ! Get the face length for the current face
      !
      ds = norm2( geofa(nf)%v(1:nr,k) )
      dsl = cl(left_geom)*ds
      dsr = cl(right_geom)*ds
      !
      ! Compute the jump value at this flux point
      !
     !p_ave = half * (pl_fp(k) + pr_fp(k))
     !p_diff = abs(pl_fp(k) - pr_fp(k))
     !!
     !p_jump = p_diff / p_ave
      p_jump = two * abs(pl_fp(k) - pr_fp(k)) / (pl_fp(k) + pr_fp(k))
      !
      ! Multiply the jump value by the quadrature weight for this flux point
      !
      p_jump = p_jump * std_elem(face_geom,face_order)%wts(k)
      !
      ! Add the jump*weight to the integrated jump indicator for both cells
      !
      marker(left_cell) = marker(left_cell) + dsl*p_jump
     !area_sum(left_cell) = area_sum(left_cell) + dsl
      !
      if (right_cell <= ncell) then
        marker(right_cell) = marker(right_cell) + dsr*p_jump
       !area_sum(right_cell) = area_sum(right_cell) + dsr
      end if
      !
    end do
    !
  end do
  !
 !marker(1:ncell) = marker(1:ncell) / area_sum(1:ncell)
  !
  where (marker(:) > zero)
    marker(:) = log10( marker(:) )
  else where
    marker(:) = log10( rndoff )
  end where
  !
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
    troubled_cell(:) = (scaled_marker(:) >= 1.0e-8_wp)
  end if
  !
  if (allocated(area_sum)) then
    deallocate ( area_sum , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"area_sum",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
end subroutine jump_indicator
!
!###############################################################################
!
subroutine identify_troubled_cells(dropped_orders,there_are_no_troubles)
  !
  !.. Use Statements ..
  use geovar,  only : ncell,cell
  use flowvar, only : usp
  use ovar,    only : Limiter_Option
  use order_mod, only : p_order,n_order
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: dropped_orders
  logical(lk), intent(out) :: there_are_no_troubles
  !
  !.. Local Scalars ..
  integer :: nc,n1,n2,proj_order
  integer :: this_geom,this_order
  integer :: orders_to_cutoff
  integer :: ord_mlt,ord_add,actual_order
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
  ! Compute the troubled marker
  !
  if (Limiter_Option == 1) then
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
  ! Now loop through the
  !
  do nc = 1,ncell
    !
    ! If the current cell is not troubled, skip to the next one
    !
    if (.not. troubled_cell(nc)) cycle
    !
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    actual_order = ord_mlt*this_order + ord_add
    !
    ! If the projection has reduced the solution to first order, mark
    ! this cell as no longer troubled and skip to the next cell
    !
    if (actual_order-proj_order <= 1) then
      troubled_cell(nc) = fals
      cycle
    end if
    !
    ! If the current cell is determined to be troubled, change
    ! there_are_no_troubles to fals meaning that troubles cells exist.
    !
    if (troubled_cell(nc)) there_are_no_troubles = fals
    !
  end do
  !
#ifdef PROFILE_ON
  n1 = merge(1,0,there_are_no_troubles)
  call mpi_allreduce(MPI_IN_PLACE,n1,1_int_mpi,mpi_inttyp, &
                     MPI_MIN,MPI_COMM_WORLD,mpierr)
  there_are_no_troubles = (n1 == 1)
#endif
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
#endif
end subroutine limiters_memory_usage
!
!###############################################################################
!
include "Functions/pressure.f90"
!
!###############################################################################
!
end module module_limiters
