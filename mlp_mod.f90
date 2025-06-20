!###############################################################################
!> Heiarchical multidimensional limiting process (hMLP).
!> Ref:
!>   Park & Kim (2016), J.Comp.Physics, Vol.308, pp.57-80.
!>   http://dx.doi.org/10.1016/j.jcp.2015.12.020
!###############################################################################
module mlp_mod
  !
  !.. Use Statements ..
  use module_kind_types,      global_debug_timer =>      debug_timer, &
                         pure_global_debug_timer => pure_debug_timer
  use eqn_idx_mod, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq
  use eqn_idx_mod, only : ncs
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
  public :: interp_cv_to_nodepts
  public :: limit_by_mlp
  !
  !
  !.. Private Procedures ..
  !
  private :: allocate_mlp
  private :: deallocate_mlp
  private :: compute_mlp_node_data
  private :: get_proj_to_usp
  private :: get_proj_to_usp_3d
  private :: get_proj_to_usp_2d
  private :: get_proj_to_usp_1d
  private :: get_proj_to_usp_nd
  private :: interp_cv_to_nodepts_3d
  private :: interp_cv_to_nodepts_2d
  private :: interp_cv_to_nodepts_1d
  private :: interp_cv_to_nodepts_nd
  !
  !.. Module Local Scalars ..
  !
  !> MLP marker variable. Usually set to density (nec).
  integer :: nmk
  !
  !> MLP relaxation factor in the range [0.0:1.0], where:
  !> 0 = fully unlimited (no limiting of solution), and
  !> 1 = fully limited solution.
  !> Intermediate values permit blending of the limited/unlimited solutions.
  integer :: mlp_relaxation = one
  !
  !> Flag indicating whether mlp module arrays have been allocated.
  logical(lk) :: mlp_is_allocated = fals
  !
  !.. Module Local Allocatable Arrays ..
  !
  !> MLP node data. Contains the min/max of the cell averaged value 
  !> for all cells that share the given node.  This value is only for
  !> the MLP marker target variable (ie, density).
  !> The array is sized for the nodes on the local partiton.
  real(wp), allocatable :: mlp_node_data(:,:)  ! (1:2,1:nnode)
  !
contains
  !
!###############################################################################
!> Allocate arrays and initialize values before using the mlp routines.
!###############################################################################
subroutine allocate_mlp
  !
  !.. Use Statements ..
  use geomvar_mod,  only : ncell,nnode
  use ovar,         only : output_marker_variables,MLP_Marker_Variable
  use flowvar_mod,  only : marker,scaled_marker,nmarker
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "allocate_mlp"
  !
  continue
  !
  ! Allocate Arrays
  !
  if (.not. allocated(mlp_node_data)) then
    allocate( mlp_node_data(1:2,1:nnode), source=zero, &
              stat=ierr , errmsg=error_message )
    call alloc_error(pname,"mlp_node_data",1,__LINE__,__FILE__,ierr, &
                       error_message)
  end if
  !
  ! Because MLP computes the markers as part of the limiting process,
  ! only need to store it if user requests it for output.
  !
  if ( output_marker_variables ) then
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
  end if
  !
  ! Initialize values
  !
  ! Set the marker variable.
  nmk = MLP_Marker_Variable
  !
  ! Indicate initialization has been done.
  mlp_is_allocated = true
  !
end subroutine allocate_mlp
!###############################################################################
!> Deallocate arrays after using the mlp routines.
!###############################################################################
subroutine deallocate_mlp
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "deallocate_mlp"
  !
  continue
  !
  if ( allocated(mlp_node_data) ) then
    deallocate( mlp_node_data , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"mlp_node_data",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  mlp_is_allocated = fals
  !
end subroutine deallocate_mlp
!###############################################################################
!> Compute min|max cell averaged values from all cells sharing a node.
!>
!> Store the min|max in the node-based module arrays.
!> Only stores the one conserved variable chosen as the marker.
!> Values based only on the local partition, no parallel communication.
!> Operates on all cell nodes.
!###############################################################################
subroutine compute_mlp_node_data(cell_ave_sol)
  !
  use geomvar_mod,  only : ncell,cell, &
                           nodes_of_cell_ptr, nodes_of_cell
! use geomvar_mod,  only : cells_with_node_ptr, cells_with_node
  use parallel_mod, only : exch_mlp_data
  !
  !.. Arguments ..
  !> Cell averaged solution variables for all cells.
  real(wp), intent(in) :: cell_ave_sol(:,:)  ! (1:nq,1:ncell)
  !
  !.. Local Scalars ..
  integer  :: nc,n,nn
  real(wp) :: var
  !
continue
  !
  ! Compute min/max value from all cells on local partition with this
  ! local node.
  !
  mlp_node_data(1,:) =  huge(one)  ! largest  possible value
  mlp_node_data(2,:) = -huge(one)  ! smallest possible value
  !
  ! loop cell
  do nc = 1, ncell
    !
    ! loop node of cell
    do n = 1,size(cell(nc)%node)
      nn = cell(nc)%node(n)%idx  ! index in nodeusp array
      !
      var = cell_ave_sol(nmk,nc)
      mlp_node_data(1,nn) = min( mlp_node_data(1,nn), var )
      mlp_node_data(2,nn) = max( mlp_node_data(2,nn), var )
    end do ! loop n
  end do ! loop nc
  !
  ! Exchange values both ways with other partitions to obtain
  ! min/max value from all cells on all partitions with this local node.
  !
  if ( ncpu > 1 ) then
    call exch_mlp_data( mlp_node_data )
  end if
  !
end subroutine compute_mlp_node_data
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
!###############################################################################
!> Project the nodal solution for the given cell to a lower order.
!>
!> Returns the projected solution on the original points without modifying
!> the original solution array.  Operates on a single cell.
!###############################################################################
subroutine get_proj_to_usp(cell,proj_order,cell_usp_proj)
  !
  !.. Use Statements ..
  use geomvar_mod, only : cell_t
  use flowvar_mod, only : usp
  use filter_mod,  only : filter
  !
  !.. Formal Arguments ..
  !> Cell information.
  type(cell_t), intent(in)  :: cell
  !
  !> Projection order.
  integer,      intent(in)  :: proj_order
  !
  !> Projected Pm nodal solution point values for this cell, where Pm < Pn.
  !> Working range is (1:np,1:nq), but array may be sized larger.
  real(wp),     intent(out) :: cell_usp_proj(:,:)
  !
  !.. Local Scalars ..
  integer  :: nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_proj_to_usp"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = nme-nmb+1
  !
  if (force_nd_subroutines) then
    !
    call get_proj_to_usp_nd(cell,proj_order,cell_usp_proj)
    !
  else
    !
    if (nr == 3) then
      call get_proj_to_usp_3d(cell,proj_order,cell_usp_proj)
    else if (nr == 2) then
      call get_proj_to_usp_2d(cell,proj_order,cell_usp_proj)
    else
      call get_proj_to_usp_1d(cell,proj_order,cell_usp_proj)
    end if
    !
    if (ncs > 0) then
      call get_proj_to_usp_nd(cell,proj_order,cell_usp_proj)
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_proj_to_usp
!###############################################################################
!> Project the nodal solution for the given cell to a lower order.
!>
!> Optimized routine for 3d.
!###############################################################################
subroutine get_proj_to_usp_3d(cell,proj_order,cell_usp_proj)
  !
  !.. Use Statements ..
  use geomvar_mod, only : cell_t
  use flowvar_mod, only : usp
  use filter_mod,  only : filter
  !
  !.. Formal Arguments ..
  !> Cell information.
  type(cell_t), intent(in)  :: cell
  !
  !> Projection order.
  integer,      intent(in)  :: proj_order
  !
  !> Projected Pm nodal solution point values for this cell, where Pm < Pn.
  !> Working range is (1:np,1:nq), but array may be sized larger.
  real(wp),     intent(out) :: cell_usp_proj(:,:)
  !
  !.. Local Scalars ..
  integer  :: this_geom, this_order
  integer  :: i,k,n,n0,np
  real(wp) :: a, var_ec, var_mx, var_my, var_mz, var_ee
  integer  :: cutoff_order
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_proj_to_usp_3d"
  !
continue
  !
  n0 = cell%beg_sp - 1     ! offset index for sol-pts in this cell
  np = cell%end_sp - n0    ! number of sol-pts in this cell
  !
  this_geom  = cell%geom   ! this cell geometry type
  this_order = cell%order  ! this cell solution order
  !
  if ( proj_order == this_order ) then
    !
    ! Project to same order, simply copy values.
    !
    do n = 1,np  ! solution points
      cell_usp_proj(n,nec) = usp(n0+n,nec)
      cell_usp_proj(n,nmx) = usp(n0+n,nmx)
      cell_usp_proj(n,nmy) = usp(n0+n,nmy)
      cell_usp_proj(n,nmz) = usp(n0+n,nmz)
      cell_usp_proj(n,nee) = usp(n0+n,nee)
    end do ! loop n
  else
    !
    ! Project by multplying the cutoff filter matrix with the
    ! solution in the current cell. Done this way because the
    ! projection filter matrix may not be defined for all lower orders.
    !
    ! To project to a given order, must cutoff all higher modes.
    cutoff_order = proj_order + 1   ! cut modes at/above this order
    !
    do n = 1,np  ! solution points
      !
      var_ec = zero
      var_mx = zero
      var_my = zero
      var_mz = zero
      var_ee = zero
      !
      do i = filter(this_geom,this_order)%cutoff(cutoff_order)%col_ptr(n)+1, &
             filter(this_geom,this_order)%cutoff(cutoff_order)%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%cutoff(cutoff_order)%col_idx(i)
        a = filter(this_geom,this_order)%cutoff(cutoff_order)%col_val(i)
        !
        var_ec = var_ec + a * usp(n0+k,nec)
        var_mx = var_mx + a * usp(n0+k,nmx)
        var_my = var_my + a * usp(n0+k,nmy)
        var_mz = var_mz + a * usp(n0+k,nmz)
        var_ee = var_ee + a * usp(n0+k,nee)
        !
      end do ! loop i
      !
      ! cell_usp_proj(1:np,1:nq) for 1 cell.  No n0 offset to first index.
      ! Array size may be larger than working range.
      cell_usp_proj(n,nec) = var_ec
      cell_usp_proj(n,nmx) = var_mx
      cell_usp_proj(n,nmy) = var_my
      cell_usp_proj(n,nmz) = var_mz
      cell_usp_proj(n,nee) = var_ee
      !
    end do ! loop n
    !
  end if
  !
end subroutine get_proj_to_usp_3d
!###############################################################################
!> Project the nodal solution for the given cell to a lower order.
!>
!> Optimized routine for 2d.
!###############################################################################
subroutine get_proj_to_usp_2d(cell,proj_order,cell_usp_proj)
  !
  !.. Use Statements ..
  use geomvar_mod, only : cell_t
  use flowvar_mod, only : usp
  use filter_mod,  only : filter
  !
  !.. Formal Arguments ..
  !> Cell information.
  type(cell_t), intent(in)  :: cell
  !
  !> Projection order.
  integer,      intent(in)  :: proj_order
  !
  !> Projected Pm nodal solution point values for this cell, where Pm < Pn.
  !> Working range is (1:np,1:nq), but array may be sized larger.
  real(wp),     intent(out) :: cell_usp_proj(:,:)
  !
  !.. Local Scalars ..
  integer  :: this_geom, this_order
  integer  :: i,k,n,n0,np
  real(wp) :: a, var_ec, var_mx, var_my, var_ee
  integer  :: cutoff_order
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_proj_to_usp_2d"
  !
continue
  !
  n0 = cell%beg_sp - 1     ! offset index for sol-pts in this cell
  np = cell%end_sp - n0    ! number of sol-pts in this cell
  !
  this_geom  = cell%geom   ! this cell geometry type
  this_order = cell%order  ! this cell solution order
  !
  if ( proj_order == this_order ) then
    !
    ! Project to same order, simply copy values.
    !
    do n = 1,np  ! solution points
      cell_usp_proj(n,nec) = usp(n0+n,nec)
      cell_usp_proj(n,nmx) = usp(n0+n,nmx)
      cell_usp_proj(n,nmy) = usp(n0+n,nmy)
      cell_usp_proj(n,nee) = usp(n0+n,nee)
    end do ! loop n
  else
    !
    ! Project by multplying the cutoff filter matrix with the
    ! solution in the current cell. Done this way because the
    ! projection filter matrix may not be defined for all lower orders.
    !
    ! To project to a given order, must cutoff all higher modes.
    cutoff_order = proj_order + 1   ! cut modes at/above this order
    !
    do n = 1,np  ! solution points
      !
      var_ec = zero
      var_mx = zero
      var_my = zero
      var_ee = zero
      !
      do i = filter(this_geom,this_order)%cutoff(cutoff_order)%col_ptr(n)+1, &
             filter(this_geom,this_order)%cutoff(cutoff_order)%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%cutoff(cutoff_order)%col_idx(i)
        a = filter(this_geom,this_order)%cutoff(cutoff_order)%col_val(i)
        !
        var_ec = var_ec + a * usp(n0+k,nec)
        var_mx = var_mx + a * usp(n0+k,nmx)
        var_my = var_my + a * usp(n0+k,nmy)
        var_ee = var_ee + a * usp(n0+k,nee)
        !
      end do ! loop i
      !
      ! cell_usp_proj(1:np,1:nq) for 1 cell.  No n0 offset to first index.
      ! Array size may be larger than working range.
      cell_usp_proj(n,nec) = var_ec
      cell_usp_proj(n,nmx) = var_mx
      cell_usp_proj(n,nmy) = var_my
      cell_usp_proj(n,nee) = var_ee
      !
    end do ! loop n
    !
  end if
  !
end subroutine get_proj_to_usp_2d
!###############################################################################
!> Project the nodal solution for the given cell to a lower order.
!>
!> Optimized routine for 1d.
!###############################################################################
subroutine get_proj_to_usp_1d(cell,proj_order,cell_usp_proj)
  !
  !.. Use Statements ..
  use geomvar_mod, only : cell_t
  use flowvar_mod, only : usp
  use filter_mod,  only : filter
  !
  !.. Formal Arguments ..
  !> Cell information.
  type(cell_t), intent(in)  :: cell
  !
  !> Projection order.
  integer,      intent(in)  :: proj_order
  !
  !> Projected Pm nodal solution point values for this cell, where Pm < Pn.
  !> Working range is (1:np,1:nq), but array may be sized larger.
  real(wp),     intent(out) :: cell_usp_proj(:,:)
  !
  !.. Local Scalars ..
  integer  :: this_geom, this_order
  integer  :: i,k,n,n0,np
  real(wp) :: a, var_ec, var_mx, var_ee
  integer  :: cutoff_order
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_proj_to_usp_1d"
  !
continue
  !
  n0 = cell%beg_sp - 1     ! offset index for sol-pts in this cell
  np = cell%end_sp - n0    ! number of sol-pts in this cell
  !
  this_geom  = cell%geom   ! this cell geometry type
  this_order = cell%order  ! this cell solution order
  !
  if ( proj_order == this_order ) then
    !
    ! Project to same order, simply copy values.
    !
    do n = 1,np  ! solution points
      cell_usp_proj(n,nec) = usp(n0+n,nec)
      cell_usp_proj(n,nmx) = usp(n0+n,nmx)
      cell_usp_proj(n,nee) = usp(n0+n,nee)
    end do ! loop n
  else
    !
    ! Project by multplying the cutoff filter matrix with the
    ! solution in the current cell. Done this way because the
    ! projection filter matrix may not be defined for all lower orders.
    !
    ! To project to a given order, must cutoff all higher modes.
    cutoff_order = proj_order + 1   ! cut modes at/above this order
    !
    do n = 1,np  ! solution points
      !
      var_ec = zero
      var_mx = zero
      var_ee = zero
      !
      do i = filter(this_geom,this_order)%cutoff(cutoff_order)%col_ptr(n)+1, &
             filter(this_geom,this_order)%cutoff(cutoff_order)%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%cutoff(cutoff_order)%col_idx(i)
        a = filter(this_geom,this_order)%cutoff(cutoff_order)%col_val(i)
        !
        var_ec = var_ec + a * usp(n0+k,nec)
        var_mx = var_mx + a * usp(n0+k,nmx)
        var_ee = var_ee + a * usp(n0+k,nee)
        !
      end do ! loop i
      !
      ! cell_usp_proj(1:np,1:nq) for 1 cell.  No n0 offset to first index.
      ! Array size may be larger than working range.
      cell_usp_proj(n,nec) = var_ec
      cell_usp_proj(n,nmx) = var_mx
      cell_usp_proj(n,nee) = var_ee
      !
    end do ! loop n
    !
  end if
  !
end subroutine get_proj_to_usp_1d
!###############################################################################
!> Project the nodal solution for the given cell to a lower order.
!>
!> General routine for nd.
!###############################################################################
subroutine get_proj_to_usp_nd(cell,proj_order,cell_usp_proj)
  !
  !.. Use Statements ..
  use geomvar_mod, only : cell_t
  use flowvar_mod, only : usp
  use filter_mod,  only : filter
  !
  !.. Formal Arguments ..
  !> Cell information.
  type(cell_t), intent(in)  :: cell
  !
  !> Projection order.
  integer,      intent(in)  :: proj_order
  !
  !> Projected Pm nodal solution point values for this cell, where Pm < Pn.
  !> Working range is (1:np,1:nq), but array may be sized larger.
  real(wp),     intent(out) :: cell_usp_proj(:,:)
  !
  !.. Local Scalars ..
  integer  :: this_geom, this_order
  integer  :: i,k,m,n,n0,np
  real(wp) :: a, var
  integer  :: cutoff_order
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_proj_to_usp_nd"
  !
continue
  !
  n0 = cell%beg_sp - 1     ! offset index for sol-pts in this cell
  np = cell%end_sp - n0    ! number of sol-pts in this cell
  !
  this_geom  = cell%geom   ! this cell geometry
  this_order = cell%order  ! this cell order
  !
  if ( proj_order == this_order ) then
    !
    ! Project to same order, simply copy values.
    !
    do m = 1,nq  ! solution variables
      do n = 1,np  ! solution points
        cell_usp_proj(n,m) = usp(n0+n,m)
      end do ! loop n
    end do ! loop m
  else
    !
    ! Project by multplying the cutoff filter matrix with the
    ! solution in the current cell. Done this way because the
    ! projection filter matrix may not be defined for all lower orders.
    !
    ! To project to a given order, must cutoff all higher modes.
    cutoff_order = proj_order + 1   ! cut modes at/above this order
    !
    do m = 1,nq  ! solution variables
      do n = 1,np  ! solution points
        var = zero
        do i = filter(this_geom,this_order)%cutoff(cutoff_order)%col_ptr(n)+1, &
               filter(this_geom,this_order)%cutoff(cutoff_order)%col_ptr(n+1)
          !
          k = filter(this_geom,this_order)%cutoff(cutoff_order)%col_idx(i)
          a = filter(this_geom,this_order)%cutoff(cutoff_order)%col_val(i)
          !
          var = var + a * usp(n0+k,m)
          !
        end do ! loop i
        !
        ! cell_usp_proj(1:np,1:nq) for 1 cell.  No n0 offset to first index.
        ! Array size may be larger than working range.
        cell_usp_proj(n,m) = var
        !
      end do ! loop n
    end do ! loop m
    !
  end if
  !
end subroutine get_proj_to_usp_nd
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
!###############################################################################
!> Interpolate the solution conserved variables to the nodes of the given cell.
!> Operates on a single cell.
!###############################################################################
subroutine interp_cv_to_nodepts(cell,n0,usp,cell_nodeusp)
  !
  !.. Use Statements ..
  use geomvar_mod,       only : cell_t
  use interpolation_mod, only : interp
  !
  !.. Formal Arguments ..
  !> Cell information.
  type(cell_t), intent(in)  :: cell
  !
  !> Cell index solution point offset for usp.
  !> Use cell%beg_sp-1 if usp contains data for multiple cells.
  !> Use 0             if usp contains data for a single cell.
  integer,      intent(in)  :: n0
  !
  !> Conserved variables at cell solution points.
  !> This could be flowvar_mod:usp which contains data at all cells,
  !> or a local usp array sized for only a single cell.
  !> Working range is (n0+1:n0+np,1:nq), but array may be sized larger.
  real(wp),     intent(in)  :: usp(:,:)
  !
  !> Conserved variables at node points of this single cell.
  !> Working range is (1:np,1:nq), but array may be sized larger.
  real(wp),     intent(out) :: cell_nodeusp(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_cv_to_nodepts"
  !
  !.. Local Scalars ..
  integer  :: nr
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = nme-nmb+1
  !
  if (force_nd_subroutines) then
    !
    call interp_cv_to_nodepts_nd(cell,n0,usp,cell_nodeusp)
    !
  else
    !
    if (nr == 3) then
      call interp_cv_to_nodepts_3d(cell,n0,usp,cell_nodeusp)
    else if (nr == 2) then
      call interp_cv_to_nodepts_2d(cell,n0,usp,cell_nodeusp)
    else
      call interp_cv_to_nodepts_1d(cell,n0,usp,cell_nodeusp)
    end if
    !
    if (ncs > 0) then
      call interp_cv_to_nodepts_nd(cell,n0,usp,cell_nodeusp)
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine interp_cv_to_nodepts
!
!###############################################################################
!> Interpolate the solution conserved variables to the nodes of the given cell.
!> Operates on a single cell.
!>
!> Optimized routine for 3d.
!###############################################################################
pure &
subroutine interp_cv_to_nodepts_3d(cell,n0,usp,cell_nodeusp)
  !
  !.. Use Statements ..
  use geomvar_mod,       only : cell_t
  use interpolation_mod, only : interp
  !
  !.. Formal Arguments ..
  !> Cell information.
  type(cell_t), intent(in)  :: cell
  !
  !> Cell index solution point offset for usp.
  !> Use cell%beg_sp-1 if usp contains data for multiple cells.
  !> Use 0             if usp contains data for a single cell.
  integer,      intent(in)  :: n0
  !
  !> Conserved variables at cell solution points.
  !> This could be flowvar_mod:usp which contains data at all cells,
  !> or a local usp array sized for only a single cell.
  !> Working range is (n0+1:n0+np,1:nq), but array may be sized larger.
  real(wp),     intent(in)  :: usp(:,:)
  !
  !> Conserved variables at node points of this single cell.
  !> Working range is (1:np,1:nq), but array may be sized larger.
  real(wp),     intent(out) :: cell_nodeusp(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_cv_to_nodepts_3d"
  !
  !.. Local Scalars ..
  integer  :: this_geom,this_order
  integer  :: i,j,lf,n
  real(wp) :: a, var_ec, var_mx, var_my, var_mz, var_ee
  !
continue
  !
  this_geom = cell%geom    ! this cell geometry
  this_order = cell%order  ! this cell order
  !
  cell_nodes: do n = 1,size(cell%node)
    !
    lf = n
    !
    var_ec = zero
    var_mx = zero
    var_my = zero
    var_mz = zero
    var_ee = zero
    !
    do i = interp(this_geom,this_order)%toNode%row_ptr(lf)+1, &
           interp(this_geom,this_order)%toNode%row_ptr(lf+1)
      !
      j = interp(this_geom,this_order)%toNode%row_idx(i)
      a = interp(this_geom,this_order)%toNode%row_val(i)
      !
      var_ec = var_ec + a * usp(n0+j,nec)
      var_mx = var_mx + a * usp(n0+j,nmx)
      var_my = var_my + a * usp(n0+j,nmy)
      var_mz = var_mz + a * usp(n0+j,nmz)
      var_ee = var_ee + a * usp(n0+j,nee)
    end do ! loop i
    !
    cell_nodeusp(n,nec) = var_ec
    cell_nodeusp(n,nmx) = var_mx
    cell_nodeusp(n,nmy) = var_my
    cell_nodeusp(n,nmz) = var_mz
    cell_nodeusp(n,nee) = var_ee
    !
  end do cell_nodes
  !
end subroutine interp_cv_to_nodepts_3d
!
!###############################################################################
!> Interpolate the solution conserved variables to the nodes of the given cell.
!> Operates on a single cell.
!>
!> Optimized routine for 2d.
!###############################################################################
pure &
subroutine interp_cv_to_nodepts_2d(cell,n0,usp,cell_nodeusp)
  !
  !.. Use Statements ..
  use geomvar_mod,       only : cell_t
  use interpolation_mod, only : interp
  !
  !.. Formal Arguments ..
  !> Cell information.
  type(cell_t), intent(in)  :: cell
  !
  !> Cell index solution point offset for usp.
  !> Use cell%beg_sp-1 if usp contains data for multiple cells.
  !> Use 0             if usp contains data for a single cell.
  integer,      intent(in)  :: n0
  !
  !> Conserved variables at cell solution points.
  !> This could be flowvar_mod:usp which contains data at all cells,
  !> or a local usp array sized for only a single cell.
  !> Working range is (n0+1:n0+np,1:nq), but array may be sized larger.
  real(wp),     intent(in)  :: usp(:,:)
  !
  !> Conserved variables at node points of this single cell.
  !> Working range is (1:np,1:nq), but array may be sized larger.
  real(wp),     intent(out) :: cell_nodeusp(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_cv_to_nodepts_2d"
  !
  !.. Local Scalars ..
  integer  :: this_geom,this_order
  integer  :: i,j,lf,n
  real(wp) :: a, var_ec, var_mx, var_my, var_ee
  !
continue
  !
  this_geom = cell%geom    ! this cell geometry
  this_order = cell%order  ! this cell order
  !
  cell_nodes: do n = 1,size(cell%node)
    !
    lf = n
    !
    var_ec = zero
    var_mx = zero
    var_my = zero
    var_ee = zero
    !
    do i = interp(this_geom,this_order)%toNode%row_ptr(lf)+1, &
           interp(this_geom,this_order)%toNode%row_ptr(lf+1)
      !
      j = interp(this_geom,this_order)%toNode%row_idx(i)
      a = interp(this_geom,this_order)%toNode%row_val(i)
      !
      var_ec = var_ec + a * usp(n0+j,nec)
      var_mx = var_mx + a * usp(n0+j,nmx)
      var_my = var_my + a * usp(n0+j,nmy)
      var_ee = var_ee + a * usp(n0+j,nee)
    end do ! loop i
    !
    cell_nodeusp(n,nec) = var_ec
    cell_nodeusp(n,nmx) = var_mx
    cell_nodeusp(n,nmy) = var_my
    cell_nodeusp(n,nee) = var_ee
    !
  end do cell_nodes
  !
end subroutine interp_cv_to_nodepts_2d
!
!###############################################################################
!> Interpolate the solution conserved variables to the nodes of the given cell.
!> Operates on a single cell.
!>
!> Optimized routine for 1d.
!###############################################################################
pure &
subroutine interp_cv_to_nodepts_1d(cell,n0,usp,cell_nodeusp)
  !
  !.. Use Statements ..
  use geomvar_mod,       only : cell_t
  use interpolation_mod, only : interp
  !
  !.. Formal Arguments ..
  !> Cell information.
  type(cell_t), intent(in)  :: cell
  !
  !> Cell index solution point offset for usp.
  !> Use cell%beg_sp-1 if usp contains data for multiple cells.
  !> Use 0             if usp contains data for a single cell.
  integer,      intent(in)  :: n0
  !
  !> Conserved variables at cell solution points.
  !> This could be flowvar_mod:usp which contains data at all cells,
  !> or a local usp array sized for only a single cell.
  !> Working range is (n0+1:n0+np,1:nq), but array may be sized larger.
  real(wp),     intent(in)  :: usp(:,:)
  !
  !> Conserved variables at node points of this single cell.
  !> Working range is (1:np,1:nq), but array may be sized larger.
  real(wp),     intent(out) :: cell_nodeusp(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_cv_to_nodepts_1d"
  !
  !.. Local Scalars ..
  integer  :: this_geom,this_order
  integer  :: i,j,lf,n
  real(wp) :: a, var_ec, var_mx, var_ee
  !
continue
  !
  this_geom = cell%geom    ! this cell geometry
  this_order = cell%order  ! this cell order
  !
  cell_nodes: do n = 1,size(cell%node)
    !
    lf = n
    !
    var_ec = zero
    var_mx = zero
    var_ee = zero
    !
    do i = interp(this_geom,this_order)%toNode%row_ptr(lf)+1, &
           interp(this_geom,this_order)%toNode%row_ptr(lf+1)
      !
      j = interp(this_geom,this_order)%toNode%row_idx(i)
      a = interp(this_geom,this_order)%toNode%row_val(i)
      !
      var_ec = var_ec + a * usp(n0+j,nec)
      var_mx = var_mx + a * usp(n0+j,nmx)
      var_ee = var_ee + a * usp(n0+j,nee)
    end do ! loop i
    !
    cell_nodeusp(n,nec) = var_ec
    cell_nodeusp(n,nmx) = var_mx
    cell_nodeusp(n,nee) = var_ee
    !
  end do cell_nodes
  !
end subroutine interp_cv_to_nodepts_1d
!
!###############################################################################
!> Interpolate the solution conserved variables to the nodes of the given cell.
!> Operates on a single cell.
!>
!> General routine for nd.
!###############################################################################
pure &
subroutine interp_cv_to_nodepts_nd(cell,n0,usp,cell_nodeusp)
  !
  !.. Use Statements ..
  use geomvar_mod,       only : cell_t
  use interpolation_mod, only : interp
  !
  !.. Formal Arguments ..
  !> Cell information.
  type(cell_t), intent(in)  :: cell
  !
  !> Cell index solution point offset for usp.
  !> Use cell%beg_sp-1 if usp contains data for multiple cells.
  !> Use 0             if usp contains data for a single cell.
  integer,      intent(in)  :: n0
  !
  !> Conserved variables at cell solution points.
  !> This could be flowvar_mod:usp which contains data at all cells,
  !> or a local usp array sized for only a single cell.
  !> Working range is (n0+1:n0+np,1:nq), but array may be sized larger.
  real(wp),     intent(in)  :: usp(:,:)
  !
  !> Conserved variables at node points of this single cell.
  !> Working range is (1:np,1:nq), but array may be sized larger.
  real(wp),     intent(out) :: cell_nodeusp(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_cv_to_nodepts_nd"
  !
  !.. Local Scalars ..
  integer  :: this_geom,this_order
  integer  :: i,j,lf,m,n
  real(wp) :: a, var
  !
continue
  !
  this_geom = cell%geom    ! this cell geometry
  this_order = cell%order  ! this cell order
  !
  cell_nodes: do n = 1,size(cell%node)
    !
    lf = n
    !
    do m = 1,nq  ! solution variables
      var = zero
      do i = interp(this_geom,this_order)%toNode%row_ptr(lf)+1, &
             interp(this_geom,this_order)%toNode%row_ptr(lf+1)
        !
        j = interp(this_geom,this_order)%toNode%row_idx(i)
        a = interp(this_geom,this_order)%toNode%row_val(i)
        !
        var = var + a * usp(n0+j,m)
      end do ! loop i
      !
      cell_nodeusp(n,m) = var
    end do ! loop m
    !
  end do cell_nodes
  !
end subroutine interp_cv_to_nodepts_nd
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
!###############################################################################
!> Compute the limiter for hMLP and apply to the solution polynomial.
!> Also update the troubled cell markers if requested for output.
!> 
!> Note: Technically, this is the hMLP-P1 algorithm.
!> Ref:
!>   Park and Kim (2016),
!>   "Heirarchical Multi-dimensional Limiting Strategy for Correction
!>    Procedure via Reconstruction," J.Comp.Phys, V.308,pp.57-80.
!###############################################################################
subroutine limit_by_mlp(nst,relax)
  !
  !.. Use Statements ..
  use geomvar_mod,     only : ncell,cell,face,bface,nfbnd
  use flowvar_mod,     only : usp
  use cell_ave_mod,    only : allocate_cell_ave,deallocate_cell_ave, &
                              compute_cell_ave_sol,cell_ave_sol
  use flowvar_mod,     only : marker,scaled_marker,get_marker_index
  use order_mod,       only : maxSP
  use ovar,            only : output_marker_variables
  !
  !.. Arguments ..
  integer,            intent(in) :: nst    !< Runge-Kutta stage.
  real(wp), optional, intent(in) :: relax  !< Relaxation factor, between [0,1].
  !
  !.. Local Parameters ..
  integer, parameter :: maxNP = 8          !  max nodes per cell
  character(len=*), parameter :: pname = "limit_by_mlp"
  !
  !.. Local Arrays ..
  ! MLP test at each node of cell.
  logical, dimension(maxNP) :: mlp_test
  logical, dimension(maxNP) :: node_on_bface
  !
  ! P1, Pm projected solution at each node of single cell.
  real(wp), dimension(1:maxNP,1:nq) :: cell_nodeusp_p1, cell_nodeusp_pm
  !
  ! P1, Pn, Pm projected solution at Pn solution points of single cell.
  real(wp), dimension(1:maxSP,1:nq) :: cell_usp_p1, cell_usp_pm
  !
  !.. Local Scalars ..
  integer  :: m_order,m,n,nc,n0,np,nn
  integer  :: ncf,nf,fq,nfnodes
  integer  :: this_geom,this_order,proj_order
  integer  :: num_cell_faces,num_cell_nodes
  integer  :: imarker
  !
  !> Tolerance for non-zero P1 mode. If the P1 mode is below this tolerance
  !> limiting of that mode will not be done (ie, it is already zero).
  !> Prevents divide-by-zero.
  real(wp) :: eps_p1mode = tny
  !
  !> Tolerance for deactivating MLP in a cell. Variations of the solution
  !> within a cell below this fraction of the cell average will result in
  !> the cell being marked as smooth.
  real(wp) :: eps_deact  = 0.001_wp
  !
  !> Deactivation test Eq.33 of the original paper includes a volume term
  !> intended to account for coarse grids, but the units are inconsistent
  !> with the other terms.  Set =0 for now to ignore it.
  real(wp) :: vol = zero     ! cell volume ToDo: FixMe
  !
  real(wp) :: q,qave,qave_min,qave_max,qmin,qmax,q_lin,q_hot
  real(wp) :: p1mode,term1,term2,r
  real(wp) :: cell_limiter   ! Continuous [0.,1.] = [Bad,Good] cell.
  real(wp) :: marker_val,smarker_val
  real(wp) :: relaxation     ! local variable for optional arg relax
  !
  logical(lk) :: local_min_test,local_max_test,deactivate_test,smooth_test
  logical(lk) :: good_node,good_cell,apply_slope_limiter
  !
  continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Set default value for optional argument.
  relaxation = mlp_relaxation
  if ( present(relax) ) relaxation = relax
  !
  ! Determine which marker position to fill based on Runge-Kutta stage.
  imarker = get_marker_index(nst)
  !
  ! Allocate arrays.
  call allocate_mlp()
  call allocate_cell_ave()
  !
  ! Compute cell ave solution values for each cell.
  call compute_cell_ave_sol()
  !
  ! Compute min|max of cell ave solution values for cells sharing each node.
  ! For target marker variable only.
  call compute_mlp_node_data(cell_ave_sol)
  !
  !
  loop_cells: &
  do nc = 1, ncell
    !
    ! cell info
    this_geom  = cell(nc)%geom  ! this cell geometry
    this_order = cell(nc)%order ! this cell order
    n0 = cell(nc)%beg_sp - 1    ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0   ! number of sol-pts in this cell
    num_cell_faces = size(cell(nc)%face)
    num_cell_nodes = size(cell(nc)%node)
    !
    qave = cell_ave_sol(nmk,nc)  ! cell ave value for marker variable
    !
    ! cell_usp_pm(1:np,1:nq) is Pm projection of usp for this cell.
    ! Pm=Pn to start, so copy.
    do m = 1,nq  ! solution variables
      do n = 1,np  ! solution points
        cell_usp_pm(n,m) = usp(n0+n,m)
      end do
    end do
    !
    ! cell_usp_p1(1:np,1:nq) is P1 projection of usp for this cell.
    proj_order = 1
    call get_proj_to_usp(cell(nc),proj_order,cell_usp_p1)
    !
    ! cell_nodeusp_p1(1:nn,1:nq) is interpolation of cell_usp_p1 to nodes.
    call interp_cv_to_nodepts(cell(nc),0,cell_usp_p1,cell_nodeusp_p1)
    !
    ! loop n = node of cell
    do n = 1,num_cell_nodes
      !
      ! Compute P1 MLP test.
      !
      nn = cell(nc)%node(n)%idx   ! index in nodeusp
      ! min|max cell ave values from all cells sharing this node, eq.15
      qave_min = mlp_node_data(1,nn)
      qave_max = mlp_node_data(2,nn)
      ! P1 projection at the node location
      qmin = cell_nodeusp_p1(n,nmk)
      qmax = cell_nodeusp_p1(n,nmk)
      ! P1 mlp test for a good node, eq.29: [T,F]=[Good,Bad]
      mlp_test(n) = ( ( qave_min .le. qmin ) .and. &
                      ( qmax .le. qave_max ) )
      !
      ! Determine if node is on a boundary face.
      ! Treat all boundary nodes as good, per Deneb code.
      !
      node_on_bface(n) = fals
      ! loop ncf = face of cell
      do ncf = 1,num_cell_faces
        nf = cell(nc)%face(ncf)%idx  ! global index of face
        nfnodes = num_face_nodes(ncf,this_geom)  ! num nodes on face
        ! skip faces that do not contain this node
        if ( all( nn /= face(nf)%nodes(1:nfnodes) ) ) cycle ! ncf
        ! skip interior faces
        if ( nf > nfbnd ) cycle ! ncf
        ! skip boundary communication faces
        if ( any(bface(1,nf) == bc_comm(:)) ) cycle ! ncf
        ! node is on this boundary face
        node_on_bface(n) = true
        exit ! ncf
      end do
      !
    end do ! loop n
    !
    apply_slope_limiter = true  ! in case:  this_order < 2
    !
    ! Note: If this_order<2, then m_order=this_order but loop not executed.
    loop_modes: &
    do m_order = this_order, 2, -1
      !
      ! cell_nodeusp_pm(1:nn,1:nq) is interpolation of cell_usp_pm to nodes.
      call interp_cv_to_nodepts(cell(nc),0,cell_usp_pm,cell_nodeusp_pm)
      !
      ! loop n = node of cell
      do n = 1,num_cell_nodes
        nn  = cell(nc)%node(n)%idx   ! index in nodeusp
        !
        ! good_node begins with mlp_test, then adjusts for smooth regions
        good_node = ( mlp_test(n) .or. node_on_bface(n) )  ! [T,F]=[Good,Bad]
        if ( .not. good_node ) then  ! test for smooth regions
          ! cell_nodeusp_p1(n,nmk) = P1 solution at the node location.
          ! cell_nodeusp_pm(n,nmk) = Pm solution at the node location.
          ! Separate solution q = qave + q_lin + q_hot via eq.30
          q     = cell_nodeusp_pm(n,nmk)
          q_lin = cell_nodeusp_p1(n,nmk) - qave          ! Pn_projected_slope
          q_hot = cell_nodeusp_pm(n,nmk) &
                - cell_nodeusp_p1(n,nmk)                 ! P1_filtered_Pn
          qave_min = mlp_node_data(1,nn)                 ! min cell ave
          qave_max = mlp_node_data(2,nn)                 ! max cell ave
          local_max_test = ( q_lin .gt. zero .and. &
                             q_hot .lt. zero .and. &
                             q     .gt. qave_min   )  ! C1 test, eq.32
          local_min_test = ( q_lin .lt. zero .and. &
                             q_hot .gt. zero .and. &
                             q     .lt. qave_max   )  ! C2 test, eq.32
          deactivate_test = ( abs( q - qave ) .le. &
                              max( eps_deact * abs(qave), vol ) )  ! eq.33
          smooth_test = ( local_max_test .or. local_min_test .or. &
                          deactivate_test )
          good_node = smooth_test
          ! any bad node marks a bad cell, so stop testing
          if ( .not. good_node ) exit  ! break loop n
        end if
      end do ! n
      !
      good_cell = good_node
      apply_slope_limiter = fals
      !
      if ( good_cell ) then
        exit loop_modes ! Break loop m_order. m_order will keep current value.
      else
        if ( m_order .gt. 2 ) then
          ! Project the solution down to order P(m-1).
          ! cell_usp_pm(1:np,1:nq) is P(m-1) projection for this cell.
          proj_order = m_order - 1
          call get_proj_to_usp(cell(nc),proj_order,cell_usp_pm)
          ! Naturally cycle loop_modes.
        else if ( m_order .eq. 2 ) then
          ! Project the solution down to order P1.
          ! Only constant and linear modes remain. Apply slope limiter below.
          apply_slope_limiter = true
          ! Naturally exit loop_modes. m_order will have value 1.
        end if ! check m_order
      end if ! check good_cell
      !
    end do loop_modes ! loop m_order
    !
    ! Compute the P1 slope limiter.
    !
    cell_limiter = one  ! Default is no limiting
    if ( this_order >= 1 .and. apply_slope_limiter ) then
      ! Apply the mlp-u1 slope limiter via eq.20-21, algorithm 1.
      !
      ! loop n = node of cell
      do n = 1,size(cell(nc)%node)
        nn  = cell(nc)%node(n)%idx   ! index in nodeusp
        ! cell_nodeusp_p1(n,nmk) = P1 solution at the node location.
        ! Compute the P1 mode at the node location.
        p1mode = cell_nodeusp_p1(n,nmk) - qave
        if ( abs(p1mode) .ge. eps_p1mode ) then
          qave_min = mlp_node_data(1,nn)    ! min cell ave
          qave_max = mlp_node_data(2,nn)    ! max cell ave
          term1 = ( qave_min - qave ) / p1mode
          term2 = ( qave_max - qave ) / p1mode
          r = max( term1, term2 )
          cell_limiter = min( cell_limiter, min(r,one) )
        end if
      end do ! loop n
      ! Reconstruct the Pm solution from the average and limited p1 mode.
      ! q(x) = qave + cell_limiter * p1mode(x)
      do m = 1,nq  ! solution variables
        qave = cell_ave_sol(m,nc)
        do n = 1,np  ! solution points
          p1mode =  cell_usp_p1(n,m) - qave  ! P1 mode at solution point.
          cell_usp_pm(n,m) = qave + cell_limiter * p1mode
        end do ! loop n
      end do ! loop m
    else if ( this_order == 0 ) then
      ! Pn=P0 solution, no P1 mode to begin with.
      ! Set cell_limiter so marker=0 (ie, order P0).
      cell_limiter = zero
    else
      ! Pm solution with m<=n, but no P1 limiting.
      cell_limiter = one
    end if
    !
    ! cell_usp_pm(1:np,1:nq) now contains the limited solution.
    !
    ! Raw marker indicates the polynomial order of the limited solution,
    ! with possible real values of [0:1,2,3,...,Pn=this_order]
    marker_val  = real(m_order-1 ,wp) + cell_limiter
    !
    ! Scaled marker is a binary flag to mark troubled cells,
    ! with possible real values of [0,1]=[Good,Bad] cell.
    smarker_val = merge(zero,one,marker_val==real(this_order,wp))
    !
    ! Set troubled cell markers, if requested for output.
    !
    if ( output_marker_variables ) then
      !
      marker(nc,imarker) = marker_val
      scaled_marker(nc,imarker) = smarker_val
      !
    end if
    !
    if ( this_order > 0 ) then
      !
      ! Replace usp with the limited soluton for this cell in cell_usp_pm.
      ! Relaxation factor permits blending the unlimited/limited solutions.
      ! Typically has value 1. so blending not used.
      !
      do m = 1,nq  ! solution variables
        do n = 1,np  ! solution points
          usp(n0+n,m) = (    relaxation) * cell_usp_pm(n,m) &
                      + (one-relaxation) * usp(n0+n,m)
        end do
      end do
      !
    end if
    !
  end do loop_cells
  !
  ! Deallocate arrays.
  call deallocate_mlp()
  call deallocate_cell_ave()
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine limit_by_mlp
!###############################################################################
!
#include "Functions/debug_timer.f90"
!
!###############################################################################
!
end module mlp_mod
