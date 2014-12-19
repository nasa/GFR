module filter_mod
  !
  !.. Use Statements ..
  use module_kind_types
  !
  use generic_types_mod, only : matrix
  !
  use ovar, only : Filter_Option,Filter_etac
  use ovar, only : Filter_Order,Filter_Alpha
  !
  implicit none
  !
  private
  !
  !
  ! Module Public Procedures
  !
  public :: initialize_filter_matrices
  !
  ! Module public arrays
  !
  ! filter : 
  !
  type(matrix), public, save, target, allocatable :: filter(:,:)
  !
  ! Generic interfaces for local module functions
  !
  interface get_quad_filter
    module procedure get_quad_filter_DP, get_quad_filter_QP
  end interface get_quad_filter
  !
  interface get_tria_filter
    module procedure get_tria_filter_DP, get_tria_filter_QP
  end interface get_tria_filter
  !
contains
!
!###############################################################################
!
subroutine initialize_filter_matrices()
  !
  !.. Use Statements ..
  use order_mod, only : n_order,geom_solpts
  use vandermonde_mod, only : vand
  !
  !.. Local Scalars ..
  integer :: ierr,gmin,gmax,omin,omax
  integer :: this_geom,this_order,npts
  character(len=200) :: error_message
  character(len=200) :: array_name
  !
  !.. Local Pointers ..
  real(wp), pointer, contiguous :: fmat(:,:)
  real(wp), pointer, contiguous :: vmat(:,:)
  real(wp), pointer, contiguous :: vmat_inv(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_filter_matrices"
  !
continue
  !
  return
  !
  fmat => null()
  vmat => null()
  vmat_inv => null()
  !
  ! Allocate the filter array
  !
  gmin = Geom_Edge
  gmax = Geom_Quad
  omin = n_order
  omax = n_order
  !
  allocate ( filter(gmin:gmax,omin:omax) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"filter",1,ierr,error_message)
  !
  ! ###########################################################################
  ! ###########################################################################
  ! ###########################################################################
  ! ##    NEED TO FIGURE OUT HOW TO CREATE THE FILTERS IN QUAD PRECISION     ##
  ! ## AND THEN CHOP TO WORKING PRECISION BEFORE STORING IN THE FILTER ARRAY ##
  ! ###########################################################################
  ! ###########################################################################
  ! ###########################################################################
  !
  ! Initialize the Edge filter to zero
  !
  this_geom = Geom_Edge
  this_order = n_order
  !
  npts = geom_solpts(this_geom,this_order)
  !
  allocate ( filter(this_geom,this_order)%mat(1:npts,1:npts) , source=zero , &
             stat=ierr , errmsg=error_message )
  write (array_name,1) Geom_Name(this_geom),this_order
  call alloc_error(pname,array_name,1,ierr,error_message)
  !
  ! Create the filter for triangle elements
  !
  this_geom = Geom_Tria
  this_order = n_order
  !
  npts = geom_solpts(this_geom,this_order)
  !
  allocate ( filter(this_geom,this_order)%mat(1:npts,1:npts) , source=zero , &
             stat=ierr , errmsg=error_message )
  write (array_name,1) Geom_Name(this_geom),this_order
  call alloc_error(pname,array_name,1,ierr,error_message)
  !
 !associate ( fmat => filter(this_geom,this_order)%mat , &
 !            vmat => vand(this_geom,this_order)%modal2nodal%mat , &
 !            vmat_inv => vand(this_geom,this_order)%nodal2modal%mat )
 !  fmat = get_tria_filter(vmat,vmat_inv)
 !end associate
  !
  fmat => filter(this_geom,this_order)%mat
  vmat => vand(this_geom,this_order)%modal2nodal%mat
  vmat_inv => vand(this_geom,this_order)%nodal2modal%mat
  !
  fmat = get_tria_filter(vmat,vmat_inv)
  !
  if (associated(fmat)) fmat => null()
  if (associated(vmat)) vmat => null()
  if (associated(vmat_inv)) vmat_inv => null()
  !
  ! Create the filter for quadrilateral elements
  !
  this_geom = Geom_Quad
  this_order = n_order
  !
  npts = geom_solpts(this_geom,this_order)
  !
  allocate ( filter(this_geom,this_order)%mat(1:npts,1:npts) , source=zero , &
             stat=ierr , errmsg=error_message )
  write (array_name,1) Geom_Name(this_geom),this_order
  call alloc_error(pname,array_name,1,ierr,error_message)
  !
  fmat => filter(this_geom,this_order)%mat
  vmat => vand(this_geom,this_order)%modal2nodal%mat
  vmat_inv => vand(this_geom,this_order)%nodal2modal%mat
  !
  fmat = get_quad_filter(vmat,vmat_inv)
  !
  if (associated(fmat)) fmat => null()
  if (associated(vmat)) vmat => null()
  if (associated(vmat_inv)) vmat_inv => null()
  !
  ! Format Statements
  !
  1 format ("filter(",a,",",i0,")%mat")
  !
end subroutine initialize_filter_matrices
!
!###############################################################################
!
pure function get_tria_filter_DP(vand,vand_inv) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Use Statements ..
  use order_mod, only : n_order,n_sp2dt
  !
  !.. Local Parameters ..
  integer, parameter :: lp = dp
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one = 1.0_lp
  !
  !.. Formal Arguments ..
  real(lp), dimension(1:n_sp2dt,1:n_sp2dt), intent(in) :: vand
  real(lp), dimension(1:n_sp2dt,1:n_sp2dt), intent(in) :: vand_inv
  !
  !.. Function Return Value ..
  real(lp), dimension(1:n_sp2dt,1:n_sp2dt) :: return_value
  !
  !.. Local Scalars ..
  integer  :: i,j,k
  integer  :: Filter_Cutoff
  real(lp) :: rd,rn,falpha,fetac
  !
continue
  !
  ! Compute the filter cutoff level
  !
  Filter_Cutoff = nint( max(zero, min(one,Filter_etac) ) * real(n_order,kind=lp) )
  !
  ! Covert Filter_Alpha to local precision
  !
  falpha = real( Filter_Alpha , kind=lp )
  fetac  = real( Filter_etac , kind=lp )
  !
  ! Initialize the local filter matrices to identity matrices
  !
  return_value(:,:) = zero
  do i = 1,n_sp2dt
    return_value(i,i) = one
  end do
  !
  ! Compute the filter diagonal for triangles
  !
  if (Filter_Option /= 0) then
    !
    k = 1
    do i = 0,n_order
      do j = 0,n_order-i
        !
        if (i+j >= Filter_Cutoff) then
          !
          rn = real(i,kind=lp)/real(n_order,kind=lp) - fetac
          rd = one - fetac
          !
          return_value(k,k) = exp(-falpha * (rn/rd)**Filter_Order)
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
  return_value = matmul( vand , return_value )
  return_value = matmul( return_value , vand_inv )
  !
end function get_tria_filter_DP
!
!###############################################################################
!
pure function get_tria_filter_QP(vand,vand_inv) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Use Statements ..
  use order_mod, only : n_order,n_sp2dt
  !
  !.. Local Parameters ..
  integer, parameter :: lp = qp
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one = 1.0_lp
  !
  !.. Formal Arguments ..
  real(lp), dimension(1:n_sp2dt,1:n_sp2dt), intent(in) :: vand
  real(lp), dimension(1:n_sp2dt,1:n_sp2dt), intent(in) :: vand_inv
  !
  !.. Function Return Value ..
  real(lp), dimension(1:n_sp2dt,1:n_sp2dt) :: return_value
  !
  !.. Local Scalars ..
  integer  :: i,j,k
  integer  :: Filter_Cutoff
  real(lp) :: rd,rn,falpha,fetac
  !
continue
  !
  ! Compute the filter cutoff level
  !
  Filter_Cutoff = nint( max(zero, min(one,Filter_etac) ) * real(n_order,kind=lp) )
  !
  ! Covert Filter_Alpha to local precision
  !
  falpha = real( Filter_Alpha , kind=lp )
  fetac  = real( Filter_etac , kind=lp )
  !
  ! Initialize the local filter matrices to identity matrices
  !
  return_value(:,:) = zero
  do i = 1,n_sp2dt
    return_value(i,i) = one
  end do
  !
  ! Compute the filter diagonal for triangles
  !
  if (Filter_Option /= 0) then
    !
    k = 1
    do i = 0,n_order
      do j = 0,n_order-i
        !
        if (i+j >= Filter_Cutoff) then
          !
          rn = real(i,kind=lp)/real(n_order,kind=lp) - fetac
          rd = one - fetac
          !
          return_value(k,k) = exp(-falpha * (rn/rd)**Filter_Order)
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
  return_value = matmul( vand , return_value )
  return_value = matmul( return_value , vand_inv )
  !
end function get_tria_filter_QP
!
!###############################################################################
!
pure function get_quad_filter_DP(vand,vand_inv) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Use Statements ..
  use order_mod, only : n_order,n_sp1de,n_sp2dq
  !
  !.. Local Parameters ..
  integer, parameter :: lp = dp
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one = 1.0_lp
  !
  !.. Formal Arguments ..
  real(lp), dimension(1:n_sp2dq,1:n_sp2dq), intent(in) :: vand
  real(lp), dimension(1:n_sp2dq,1:n_sp2dq), intent(in) :: vand_inv
  !
  !.. Function Return Value ..
  real(lp), dimension(1:n_sp2dq,1:n_sp2dq) :: return_value
  !
  !.. Local Scalars ..
  integer  :: i,j,k
  integer  :: Filter_Cutoff
  real(lp) :: rd,rn,falpha,fetac
  !
  !.. Local Arrays ..
  real(lp), dimension(1:n_sp1de) :: filter_1D
  !
continue
  !
  ! Compute the filter cutoff level
  !
  Filter_Cutoff = nint( max(zero, min(one,Filter_etac) ) * real(n_order,kind=lp) )
  !
  ! Covert Filter_Alpha to local precision
  !
  falpha = real( Filter_Alpha , kind=lp )
  fetac  = real( Filter_etac , kind=lp )
  !
  ! Initialize the local filter matrices to identity matrices
  !
  filter_1D(:) = one
  !
  return_value(:,:) = zero
  do i = 1,n_sp2dq
    return_value(i,i) = one
  end do
  !
  ! Compute the filter diagonal for quadrilaterals
  !
  if (Filter_Option == 1) then
    !
    ! Option 1 - Tensor product of the 1D filter
    !           (Equivalent to using a 1D filter in each direction)
    !
    ! First get the 1D filter diagonal
    !
    do i = 0,n_order
      if (i >= Filter_Cutoff) then
        !
        rn = real(i,kind=lp)/real(n_order,kind=lp) - fetac
        rd = one - fetac
        !
        filter_1D(i+1) = exp(-falpha * (rn/rd)**Filter_Order)
        !
      end if
    end do
    !
    ! Now fill the filter diagonal for quads
    !
    do j = 0,n_order
      do i = 0,n_order
        !
        k = n_sp1de*j + i + 1
        !
        return_value(k,k) = filter_1D(i+1)*filter_1D(j+1)
        !
      end do
    end do
    !
  else if (Filter_Option == 2) then
    !
    ! Option 2 - The quad basis functions actually contains cross polynomial
    !            terms up to order 2*p (as opposed to just 1*p for triangles).
    !            This filter applies the cut-off frequency based on 2*p as
    !            opposed to basing it on just p in 1D like Option 1.
    !
    do i = 0,n_order
      do j = 0,n_order
        !
        if (i+j >= Filter_Cutoff) then
          !
          k = n_sp1de*j + i + 1
          !
          rn = real(i+j-Filter_Cutoff,kind=lp)
          rd = real(2*n_order-Filter_Cutoff,kind=lp)
          !
          return_value(k,k) = exp(-falpha * (rn/rd)**Filter_Order)
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
  return_value = matmul( vand , return_value )
  return_value = matmul( return_value , vand_inv )
  !
end function get_quad_filter_DP
!
!###############################################################################
!
pure function get_quad_filter_QP(vand,vand_inv) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Use Statements ..
  use order_mod, only : n_order,n_sp1de,n_sp2dq
  !
  !.. Local Parameters ..
  integer, parameter :: lp = qp
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one = 1.0_lp
  !
  !.. Formal Arguments ..
  real(lp), dimension(1:n_sp2dq,1:n_sp2dq), intent(in) :: vand
  real(lp), dimension(1:n_sp2dq,1:n_sp2dq), intent(in) :: vand_inv
  !
  !.. Function Return Value ..
  real(lp), dimension(1:n_sp2dq,1:n_sp2dq) :: return_value
  !
  !.. Local Scalars ..
  integer  :: i,j,k
  integer  :: Filter_Cutoff
  real(lp) :: rd,rn,falpha,fetac
  !
  !.. Local Arrays ..
  real(lp), dimension(1:n_sp1de) :: filter_1D
  !
continue
  !
  ! Compute the filter cutoff level
  !
  Filter_Cutoff = nint( max(zero, min(one,Filter_etac) ) * real(n_order,kind=lp) )
  !
  ! Covert Filter_Alpha to local precision
  !
  falpha = real( Filter_Alpha , kind=lp )
  fetac  = real( Filter_etac , kind=lp )
  !
  ! Initialize the local filter matrices to identity matrices
  !
  filter_1D(:) = one
  !
  return_value(:,:) = zero
  do i = 1,n_sp2dq
    return_value(i,i) = one
  end do
  !
  ! Compute the filter diagonal for quadrilaterals
  !
  if (Filter_Option == 1) then
    !
    ! Option 1 - Tensor product of the 1D filter
    !           (Equivalent to using a 1D filter in each direction)
    !
    ! First get the 1D filter diagonal
    !
    do i = 0,n_order
      if (i >= Filter_Cutoff) then
        !
        rn = real(i,kind=lp)/real(n_order,kind=lp) - fetac
        rd = one - fetac
        !
        filter_1D(i+1) = exp(-falpha * (rn/rd)**Filter_Order)
        !
      end if
    end do
    !
    ! Now fill the filter diagonal for quads
    !
    do j = 0,n_order
      do i = 0,n_order
        !
        k = n_sp1de*j + i + 1
        !
        return_value(k,k) = filter_1D(i+1)*filter_1D(j+1)
        !
      end do
    end do
    !
  else if (Filter_Option == 2) then
    !
    ! Option 2 - The quad basis functions actually contains cross polynomial
    !            terms up to order 2*p (as opposed to just 1*p for triangles).
    !            This filter applies the cut-off frequency based on 2*p as
    !            opposed to basing it on just p in 1D like Option 1.
    !
    do i = 0,n_order
      do j = 0,n_order
        !
        if (i+j >= Filter_Cutoff) then
          !
          k = n_sp1de*j + i + 1
          !
          rn = real(i+j-Filter_Cutoff,kind=lp)
          rd = real(2*n_order-Filter_Cutoff,kind=lp)
          !
          return_value(k,k) = exp(-falpha * (rn/rd)**Filter_Order)
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
  return_value = matmul( vand , return_value )
  return_value = matmul( return_value , vand_inv )
  !
end function get_quad_filter_QP
!
!###############################################################################
!
!
!###############################################################################
!
end module filter_mod
