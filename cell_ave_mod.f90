!###############################################################################
!> Module for computing the cell averaged solution values.
!###############################################################################
module cell_ave_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx_mod, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq,ncs
  !
  implicit none
  !
  private
  !
  !
  !.. Public Procedures ..
  !
  public :: allocate_cell_ave
  public :: deallocate_cell_ave
  public :: compute_cell_ave_sol
  public :: compute_cell_ave_var
  public :: cell_ave_memory_usage
  !
  !.. Private Procedures ..
  !
  private :: set_cell_ave_coef
  private :: compute_cell_ave_sol_3d
  private :: compute_cell_ave_sol_2d
  private :: compute_cell_ave_sol_1d
  private :: compute_cell_ave_sol_nd
  !
  !.. Module Public Local Scalars ..
  !
  !> Status flag indicating whether module arrays have been allocated.
  !> (0)  arrays not allocated yet.
  !> (1+) arrays are allocated.
  !> Each call to allocate_cell_ave/deallocate_cell_ave will 
  !> increment/decrement this counter and allocate/deallocate
  !> arrays only on first/last use.
  integer, public, protected :: cell_ave_is_allocated = 0
  !
  !.. Module Public Allocatable Arrays ..
  !
  !> Conserved variables cell average solution values for each cell.
  !> Dim (1:nq,1:ncell).
  real(wp), public, protected, allocatable :: cell_ave_sol(:,:)
  !
  !.. Module Local Arrays ..
  !
  !> Cell averaging coefficient for various geometries.
  !> Dim (Geom_Min:Geom_Max).
  !> Default is zero so that cell average will be zero unless the coefficient
  !> is set for that cell geometry. Call the set routine before use.
  real(wp) :: cell_ave_coef(Geom_Min:Geom_Max) = zero
  !
contains
  !
!###############################################################################
!> Allocate module arrays for cell averaging.
!>
!> This acts as a push request for memory, incrementing the stack counter and
!> allocating if needed.
!###############################################################################
subroutine allocate_cell_ave
  !
  !.. Use Statements ..
  use geomvar_mod, only : ncell
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "allocate_cell_ave"
  !
continue
  !
  ! Allocate array for storing the cell averaged solution
  !
  if (.not. allocated(cell_ave_sol)) then
    allocate( cell_ave_sol(1:nq,1:ncell), source=zero, &
              stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cell_ave_sol",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !
    ! Set values of the geometry coefficients used for cell averaging.
    call set_cell_ave_coef
    !
  end if
  !
  cell_ave_is_allocated = max(cell_ave_is_allocated + 1, 1)
  !
end subroutine allocate_cell_ave
!###############################################################################
!> Deallocate module arrays for cell averaging.
!>
!> This acts as a pop request for memory, decrementing the stack counter and
!> deallocating when no longer needed.
!###############################################################################
subroutine deallocate_cell_ave(force_dealloc)
  !
  !.. Formal Arguments ..
  !> Force the deallocation to occur and reset the stack counter.
  logical(lk), optional, intent(in) :: force_dealloc
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "deallocate_cell_ave"
  !
  continue
  !
  cell_ave_is_allocated = max(cell_ave_is_allocated - 1, 0)
  !
  if ( present(force_dealloc) ) then
    if (force_dealloc) cell_ave_is_allocated = 0
  end if
  !
  if ( cell_ave_is_allocated == 0 ) then
    !
    if ( allocated(cell_ave_sol) ) then
      deallocate( cell_ave_sol , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"cell_ave_sol",2,__LINE__,__FILE__, &
                       ierr,error_message)
    end if
    !
  end if
  !
end subroutine deallocate_cell_ave
!###############################################################################
!> Set values of the geometry coefficients used for cell averaging.
!###############################################################################
subroutine set_cell_ave_coef()
  !
  cell_ave_coef(Geom_Edge) = half * sqrt(two)
  cell_ave_coef(Geom_Tria) = half * sqrt(two)
  cell_ave_coef(Geom_Quad) = half
  cell_ave_coef(Geom_Tetr) = half * sqrt(three)
  cell_ave_coef(Geom_Pyra) = half * sqrt(three) / sqrt(two)
  cell_ave_coef(Geom_Pris) = half
  cell_ave_coef(Geom_Hexa) = sqrt(two)/four
  !
end subroutine set_cell_ave_coef
!###############################################################################
!> Compute cell averaged values of conserved solution variables for each cell.
!> Wrapper routine for calling optimized versions.
!###############################################################################
subroutine compute_cell_ave_sol
  !
  !.. Use Statements ..
  use ovar,      only : force_nd_subroutines
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_cell_ave_sol"
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if ( cell_ave_is_allocated == 0 ) then
    write (error_message,100)
    100 format("cell_ave arrays used before allocated.")
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  if (force_nd_subroutines) then
    !
    call compute_cell_ave_sol_nd(1,nq)
    !
  else
    !
    else if (nr == 3) then
      call compute_cell_ave_sol_3d
    else if (nr == 2) then
      call compute_cell_ave_sol_2d
    else if (nr == 1) then
      call compute_cell_ave_sol_1d
    end if
    !
    if (ncs > 0) then
      call compute_cell_ave_sol_nd(1,ncs)
      ! Mathematically, the cell average of the mixture density polynomial
      ! computed in the *_[123]d subroutines above should be equivalent to
      ! summing the cell averages of the species densities.  If you are
      ! concerned that there is enough numerical roundoff to make these
      ! inconsistent, enabling the following call to update_cell_ave_density
      ! will recompute the cell average of the mixture density by summing
      ! the averages of the species densities to make sure the cell averages
      ! are consistent.
     !call update_cell_ave_density
    end if
    !
  end if
  !
end subroutine compute_cell_ave_sol
!###############################################################################
!> Compute the conserved variables cell average solution values for each cell.
!> Optimized for 3d cells.
!###############################################################################
subroutine compute_cell_ave_sol_3d
  !
  !.. Use Statements ..
  use geomvar_mod,     only : ncell,cell
  use flowvar_mod,     only : usp
  use vandermonde_mod, only : vand
  !
  !.. Local Scalars ..
  integer  :: n,nc,n0,np
  integer  :: this_geom,this_order
  real(wp) :: a,val_ec,val_mx,val_my,val_mz,val_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_cell_ave_sol_3d"
  !
continue
  !
  ! loop cells
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1
    np = cell(nc)%end_sp-n0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    val_ec = zero
    val_mx = zero
    val_my = zero
    val_mz = zero
    val_ee = zero
    !
    ! loop solution points
    do n = 1,np
      !
      a = vand(this_geom,this_order)%nodal2modal%mat(1,n)
      val_ec = val_ec + a * usp(n0+n,nec)
      val_mx = val_mx + a * usp(n0+n,nmx)
      val_my = val_my + a * usp(n0+n,nmy)
      val_mz = val_mz + a * usp(n0+n,nmz)
      val_ee = val_ee + a * usp(n0+n,nee)
      !
    end do ! loop n
    !
    ! Fill the module valiable.
    cell_ave_sol(nec,nc) = val_ec * cell_ave_coef(this_geom)
    cell_ave_sol(nmx,nc) = val_mx * cell_ave_coef(this_geom)
    cell_ave_sol(nmy,nc) = val_my * cell_ave_coef(this_geom)
    cell_ave_sol(nmz,nc) = val_mz * cell_ave_coef(this_geom)
    cell_ave_sol(nee,nc) = val_ee * cell_ave_coef(this_geom)
    !
  end do ! loop nc
  !
end subroutine compute_cell_ave_sol_3d
!###############################################################################
!> Compute the conserved variables cell average solution values for each cell.
!> Optimized for 2d cells.
!###############################################################################
subroutine compute_cell_ave_sol_2d
  !
  !.. Use Statements ..
  use geomvar_mod,     only : ncell,cell
  use flowvar_mod,     only : usp
  use vandermonde_mod, only : vand
  !
  !.. Local Scalars ..
  integer  :: n,nc,n0,np
  integer  :: this_geom,this_order
  real(wp) :: a,val_ec,val_mx,val_my,val_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_cell_ave_sol_2d"
  !
continue
  !
  ! loop cells
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1
    np = cell(nc)%end_sp-n0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    val_ec = zero
    val_mx = zero
    val_my = zero
    val_ee = zero
    !
    ! loop solution points
    do n = 1,np
      !
      a = vand(this_geom,this_order)%nodal2modal%mat(1,n)
      val_ec = val_ec + a * usp(n0+n,nec)
      val_mx = val_mx + a * usp(n0+n,nmx)
      val_my = val_my + a * usp(n0+n,nmy)
      val_ee = val_ee + a * usp(n0+n,nee)
      !
    end do ! loop n
    !
    ! Fill the module valiable.
    cell_ave_sol(nec,nc) = val_ec * cell_ave_coef(this_geom)
    cell_ave_sol(nmx,nc) = val_mx * cell_ave_coef(this_geom)
    cell_ave_sol(nmy,nc) = val_my * cell_ave_coef(this_geom)
    cell_ave_sol(nee,nc) = val_ee * cell_ave_coef(this_geom)
    !
  end do ! loop nc
  !
end subroutine compute_cell_ave_sol_2d
!###############################################################################
!> Compute the conserved variables cell average solution values for each cell.
!> Optimized for 1d cells.
!###############################################################################
subroutine compute_cell_ave_sol_1d
  !
  !.. Use Statements ..
  use geomvar_mod,     only : ncell,cell
  use flowvar_mod,     only : usp
  use vandermonde_mod, only : vand
  !
  !.. Local Scalars ..
  integer  :: n,nc,n0,np
  integer  :: this_geom,this_order
  real(wp) :: a,val_ec,val_mx,val_my,val_mz,val_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_cell_ave_sol_1d"
  !
continue
  !
  ! loop cells
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1
    np = cell(nc)%end_sp-n0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    val_ec = zero
    val_mx = zero
    val_ee = zero
    !
    ! loop solution points
    do n = 1,np
      !
      a = vand(this_geom,this_order)%nodal2modal%mat(1,n)
      val_ec = val_ec + a * usp(n0+n,nec)
      val_mx = val_mx + a * usp(n0+n,nmx)
      val_ee = val_ee + a * usp(n0+n,nee)
      !
    end do ! loop n
    !
    ! Fill the module valiable.
    cell_ave_sol(nec,nc) = val_ec * cell_ave_coef(this_geom)
    cell_ave_sol(nmx,nc) = val_mx * cell_ave_coef(this_geom)
    cell_ave_sol(nee,nc) = val_ee * cell_ave_coef(this_geom)
    !
  end do ! loop nc
  !
end subroutine compute_cell_ave_sol_1d
!###############################################################################
!> Compute the conserved variables cell average solution values for each cell.
!> General routine for nd cells.
!###############################################################################
subroutine compute_cell_ave_sol_nd(mbeg,mend)
  !
  !.. Use Statements ..
  use geomvar_mod,     only : ncell,cell
  use flowvar_mod,     only : usp
  use vandermonde_mod, only : vand
  !
  !.. Formal Arguments ..
  integer, intent(in) :: mbeg,mend
  !
  !.. Local Scalars ..
  integer  :: m,n,nc,n0,np
  integer  :: this_geom,this_order
  real(wp) :: a,val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_cell_ave_sol_nd"
  !
continue
  !
  ! loop cells
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1
    np = cell(nc)%end_sp-n0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! loop solution variables
    do m = mbeg,mend
      !
      val = zero
      !
      ! loop solution points
      do n = 1,np
        a = vand(this_geom,this_order)%nodal2modal%mat(1,n)
        val = val + a * usp(n0+n,m)
      end do ! loop n
      !
      ! Fill the module variable.
      cell_ave_sol(m,nc) = val * cell_ave_coef(this_geom)
      !
    end do ! loop m
    !
  end do ! loop nc
  !
end subroutine compute_cell_ave_sol_nd
!###############################################################################
!> Compute the cell average mixture density from the cell averages of the
!> species densities. This is only used for a thermally perfect gas mixture.
!###############################################################################
subroutine update_cell_ave_density
  !
  !.. Use Statements ..
  use geomvar_mod, only : ncell
  !
  !.. Local Scalars ..
  integer :: nc
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_cell_ave_density"
  !
continue
  !
  do nc = 1,ncell
    cell_ave_sol(nec,nc) = sum(cell_ave_sol(1:ncs,nc))
  end do
  !
end subroutine update_cell_ave_density
!###############################################################################
!> Compute cell average value of one flow variable for each cell.
!###############################################################################
subroutine compute_cell_ave_var(meq,cell_ave)
  !
  !.. Use Statements ..
  use ovar,            only : gam,rgasref_nd
  use geomvar_mod,     only : ncell,cell
  use flowvar_mod,     only : usp
  use vandermonde_mod, only : vand
  !
  !.. Formal Arguments ..
  !
  !> Index of the variable to average.
  !> * (1:nq)  Conserved variable.
  !> * (nq+1)  Pressure.
  !> * (nq+2)  Temperature.
  integer,     intent(in) :: meq
  !
  !> The cell averaged variable for each cell.
  real(wp), intent(inout) :: cell_ave(:)
  !
  !.. Local Scalars ..
  integer  :: n,nc,nr,n0,np
  integer  :: this_geom,this_order
  real(wp) :: a,val,val_ec,val_mx,val_my,val_mz,val_ee
  real(wp) :: rho,rhou,rhov,rhow,rhoe,rhoq2
  real(wp) :: gm1,rgasref_nd_inv,p
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_cell_ave_var"
  !
continue
  !
  nr = nme-nmb+1
  !
  gm1 = gam - one
  rgasref_nd_inv = one/rgasref_nd
  !
  call set_cell_ave_coef
  !
  check_meq: &
  if (meq <= nq) then
    !
    ! Average the conserved variable meq
    !
    ! loop cells
    do nc = 1,ncell
      !
      n0 = cell(nc)%beg_sp-1
      np = cell(nc)%end_sp-n0
      !
      this_geom = cell(nc)%geom
      this_order = cell(nc)%order
      !
      val = zero
      !
      ! loop solution points
      do n = 1,np
        a = vand(this_geom,this_order)%nodal2modal%mat(1,n)
        val = val + a * usp(n0+n,meq)
      end do
      !
      cell_ave(nc) = val * cell_ave_coef(this_geom)
      !
    end do  ! loop nc
    !
  else check_meq
    !
    ! Average pressure (meq=neq+1) or temperature (meq=neq+2).
    !
    check_nr: &
    if (nr == 3) then  ! 3d
      !
      ! loop cells
      do nc = 1,ncell
        !
        n0 = cell(nc)%beg_sp-1
        np = cell(nc)%end_sp-n0
        !
        this_geom = cell(nc)%geom
        this_order = cell(nc)%order
        !
        val_ec = zero
        val_mx = zero
        val_my = zero
        val_mz = zero
        val_ee = zero
        !
        ! loop solution points
        do n = 1,np
          a = vand(this_geom,this_order)%nodal2modal%mat(1,n)
          val_ec = val_ec + a * usp(n0+n,nec)
          val_mx = val_mx + a * usp(n0+n,nmx)
          val_my = val_my + a * usp(n0+n,nmy)
          val_mz = val_mz + a * usp(n0+n,nmz)
          val_ee = val_ee + a * usp(n0+n,nee)
        end do
        !
        rho  = val_ec * cell_ave_coef(this_geom)
        rhou = val_mx * cell_ave_coef(this_geom)
        rhov = val_my * cell_ave_coef(this_geom)
        rhow = val_mz * cell_ave_coef(this_geom)
        rhoe = val_ee * cell_ave_coef(this_geom)
        !
        rhoq2 = rhou*rhou + rhov*rhov + rhow*rhow
        !
        p = gm1*(rhoe - half*rhoq2/rho)
        !
        cell_ave(nc) = merge(p,p*rgasref_nd_inv/rho,meq==nq+1)
        !
      end do  ! loop nc
      !
    else if (nr == 2) then  ! 2d
      !
      ! loop cells
      do nc = 1,ncell
        !
        n0 = cell(nc)%beg_sp-1
        np = cell(nc)%end_sp-n0
        !
        this_geom = cell(nc)%geom
        this_order = cell(nc)%order
        !
        val_ec = zero
        val_mx = zero
        val_my = zero
        val_ee = zero
        !
        ! loop solution points
        do n = 1,np
          a = vand(this_geom,this_order)%nodal2modal%mat(1,n)
          val_ec = val_ec + a * usp(n0+n,nec)
          val_mx = val_mx + a * usp(n0+n,nmx)
          val_my = val_my + a * usp(n0+n,nmy)
          val_ee = val_ee + a * usp(n0+n,nee)
        end do
        !
        rho  = val_ec * cell_ave_coef(this_geom)
        rhou = val_mx * cell_ave_coef(this_geom)
        rhov = val_my * cell_ave_coef(this_geom)
        rhoe = val_ee * cell_ave_coef(this_geom)
        !
        rhoq2 = rhou*rhou + rhov*rhov
        !
        p = gm1*(rhoe - half*rhoq2/rho)
        !
        cell_ave(nc) = merge(p,p*rgasref_nd_inv/rho,meq==nq+1)
        !
      end do  ! loop nc
      !
    else check_nr ! 1d
      !
      ! loop cells
      do nc = 1,ncell
        !
        n0 = cell(nc)%beg_sp-1
        np = cell(nc)%end_sp-n0
        !
        this_geom = cell(nc)%geom
        this_order = cell(nc)%order
        !
        val_ec = zero
        val_mx = zero
        val_ee = zero
        !
        ! loop solution points
        do n = 1,np
          a = vand(this_geom,this_order)%nodal2modal%mat(1,n)
          val_ec = val_ec + a * usp(n0+n,nec)
          val_mx = val_mx + a * usp(n0+n,nmx)
          val_ee = val_ee + a * usp(n0+n,nee)
        end do
        !
        rho  = val_ec * cell_ave_coef(this_geom)
        rhou = val_mx * cell_ave_coef(this_geom)
        rhoe = val_ee * cell_ave_coef(this_geom)
        !
        rhoq2 = rhou*rhou
        !
        p = gm1*(rhoe - half*rhoq2/rho)
        !
        cell_ave(nc) = merge(p,p*rgasref_nd_inv/rho,meq==nq+1)
        !
      end do  ! loop nc
      !
    end if check_nr
    !
  end if check_meq
  !
end subroutine compute_cell_ave_var
!###############################################################################
!> Output memory usage for this module.
!###############################################################################
subroutine cell_ave_memory_usage(iunit)
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
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module cell_ave_mod")
  2 format (" Array: ",a7,dt)
  3 format (" Total memory = ",f8.2,1x,a)
  5 format ("'",a,"'")
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  !
#endif
end subroutine cell_ave_memory_usage
!###############################################################################
end module cell_ave_mod
