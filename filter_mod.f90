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
  public :: init_filter_matrices
  public :: filter_usp
  public :: filter_memory_usage
  !
  ! Module public arrays
  !
  ! filter : 
  !
  type(matrix), public, save, target, allocatable :: filter(:,:)
  !
contains
!
!###############################################################################
!
subroutine init_filter_matrices()
  !
  !.. Use Statements ..
  use order_mod, only : geom_solpts
  use order_mod, only : n_min_geom,n_max_geom
  use order_mod, only : n_min_order,n_max_order
  use quadrature_mod, only : geom_is_used
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
  character(len=*), parameter :: pname = "init_filter_matrices"
  !
continue
  !
  ! Initialize the local pointers to disassociated
  !
  fmat => null()
  vmat => null()
  vmat_inv => null()
  !
  ! Allocate the filter array
  !
  gmin = n_min_geom
  gmax = n_max_geom
  omin = n_min_order
  omax = n_max_order
  !
  allocate ( filter(gmin:gmax,omin:omax) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"filter",1,ierr,error_message)
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
      ! Allocate the correction matrix
      !
      allocate ( filter(this_geom,this_order)%mat(1:npts,1:npts) , &
                 source=zero , stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order
      call alloc_error(pname,array_name,1,ierr,error_message)
      !
      ! Assign local pointers as aliases to simplify the code
      !
      fmat => filter(this_geom,this_order)%mat
      vmat => vand(this_geom,this_order)%modal2nodal%mat
      vmat_inv => vand(this_geom,this_order)%nodal2modal%mat
      !
      ! Compute the correction matrices depending on geometry type
      !
      select case (this_geom)
        case (Geom_Edge)
          !
          fmat = FilterMatrix_Edge( vmat , vmat_inv )
          !
        case (Geom_Tria)
          !
          fmat = FilterMatrix_Tria( vmat , vmat_inv )
          !
        case (Geom_Quad)
          !
          fmat = FilterMatrix_Quad( vmat , vmat_inv )
          !
        case (Geom_Tetr)
          !
          fmat = FilterMatrix_Tetr( vmat , vmat_inv )
          !
        case (Geom_Pyra)
          !
          fmat = FilterMatrix_Pyra( vmat , vmat_inv )
          !
        case (Geom_Pris)
          !
          fmat = FilterMatrix_Pris( vmat , vmat_inv )
          !
        case (Geom_Hexa)
          !
          fmat = FilterMatrix_Hexa( vmat , vmat_inv )
          !
        case default
          !
          write (error_message,2)
          call stop_recon(stop_mpi,pname,error_message)
          !
      end select
      !
      ! Remove all associations of the local pointers before continuing
      !
      if (associated(fmat)) fmat => null()
      if (associated(vmat)) vmat => null()
      if (associated(vmat_inv)) vmat_inv => null()
      !
    end do order_loop
    !
  end do geom_loop
  !
  ! Check for sparsity
  !
  call filter%check_for_sparsity
  !
  ! Format Statements
  !
  1 format ("filter(",a,",",i0,")%mat")
  2 format (" A grid cell of an unknown geometry type was found!")
  !
end subroutine init_filter_matrices
!
!###############################################################################
!
subroutine filter_usp
  !
  !.. Use Statements ..
  use geovar,     only : ncell,pts_in_cell_ptr
  use geovar,     only : cell_geom,cell_order
  use flowvar,    only : usp
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,m
  integer :: this_geom,this_order
  !
continue
  !
  ! Apply the filter to each cell
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1
    n2 = pts_in_cell_ptr(n+1)
    !
    this_geom = cell_geom(n)
    this_order = cell_order(n)
    !
    do m = 1,size(usp,dim=1)
#ifdef USE_MATRIX_PROCEDURES
      usp(m,n1:n2) = filter(this_geom,this_order)%mv_mult( usp(m,n1:n2) )
#else
      usp(m,n1:n2) = matmul( filter(this_geom,this_order)%mat , usp(m,n1:n2) )
#endif
    end do
    !
  end do
  !
end subroutine filter_usp
!
!###############################################################################
!
pure function FilterMatrix_Edge(vand,vand_inv) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,nsp,order
  integer  :: Filter_Cutoff,ierr
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
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
  ! Compute the filter cutoff level
  !
  Filter_Cutoff = nint( max(zero, min(one,Filter_etac) ) * real(order,kind=wp) )
  !
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  fetac_qp  = real( Filter_etac , kind=qp )
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
  do i = 0,order
    if (i >= Filter_Cutoff) then
      !
      rn_qp = real(i,kind=qp)/real(order,kind=qp) - fetac_qp
      rd_qp = qone - fetac_qp
      !
      fmat_qp(i+1,i+1) = exp(-falpha_qp * (rn_qp/rd_qp)**Filter_Order)
      !
    end if
  end do
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
pure function FilterMatrix_Quad(vand,vand_inv) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,np,nsp,order
  integer  :: Filter_Cutoff,ierr
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
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
  ! Compute the filter cutoff level
  !
  Filter_Cutoff = nint( max(zero, min(one,Filter_etac) ) * real(order,kind=wp) )
  !
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  fetac_qp  = real( Filter_etac , kind=qp )
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
  if (Filter_Option == 1) then
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
        rd_qp = qone - fetac_qp
        !
        filter_1D_qp(i+1) = exp(-falpha_qp * (rn_qp/rd_qp)**Filter_Order)
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
  else if (Filter_Option == 2) then
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
          rd_qp = real(2*order-Filter_Cutoff,kind=qp)
          !
          fmat_qp(k,k) = exp(-falpha_qp * (rn_qp/rd_qp)**Filter_Order)
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
pure function FilterMatrix_Tria(vand,vand_inv) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,np,nsp,order
  integer  :: Filter_Cutoff,ierr
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
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
  ! Compute the filter cutoff level
  !
  Filter_Cutoff = nint( max(zero, min(one,Filter_etac) ) * real(order,kind=wp) )
  !
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  fetac_qp  = real( Filter_etac , kind=qp )
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
  fmat_qp(:,:) = qzero
  do i = 1,nsp
    fmat_qp(i,i) = qone
  end do
  !
  ! Compute the filter diagonal for triangles
  !
  if (Filter_Option /= 0) then
    !
    k = 1
    do i = 0,order
      do j = 0,order-i
        !
        if (i+j >= Filter_Cutoff) then
          !
          rn_qp = real(i,kind=qp)/real(order,kind=qp) - fetac_qp
          rd_qp = qone - fetac_qp
          !
          fmat_qp(k,k) = exp(-falpha_qp * (rn_qp/rd_qp)**Filter_Order)
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
pure function FilterMatrix_Tetr(vand,vand_inv) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,np,nsp,order
  integer  :: Filter_Cutoff,ierr
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
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
  ! Compute the filter cutoff level
  !
  Filter_Cutoff = nint( max(zero, min(one,Filter_etac) ) * real(order,kind=wp) )
  !
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  fetac_qp  = real( Filter_etac , kind=qp )
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
  fmat_qp(:,:) = qzero
  do i = 1,nsp
    fmat_qp(i,i) = qone
  end do
  !
  ! Compute the filter diagonal for triangles
  !
  if (Filter_Option /= 0) then
    !
    k = 1
    do i = 0,order
      do j = 0,order-i
        !
        if (i+j >= Filter_Cutoff) then
          !
          rn_qp = real(i,kind=qp)/real(order,kind=qp) - fetac_qp
          rd_qp = qone - fetac_qp
          !
          fmat_qp(k,k) = exp(-falpha_qp * (rn_qp/rd_qp)**Filter_Order)
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
pure function FilterMatrix_Pyra(vand,vand_inv) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,np,nsp,order
  integer  :: Filter_Cutoff,ierr
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
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
  ! Compute the filter cutoff level
  !
  Filter_Cutoff = nint( max(zero, min(one,Filter_etac) ) * real(order,kind=wp) )
  !
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  fetac_qp  = real( Filter_etac , kind=qp )
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
  fmat_qp(:,:) = qzero
  do i = 1,nsp
    fmat_qp(i,i) = qone
  end do
  !
  ! Compute the filter diagonal for triangles
  !
  if (Filter_Option /= 0) then
    !
    k = 1
    do i = 0,order
      do j = 0,order-i
        !
        if (i+j >= Filter_Cutoff) then
          !
          rn_qp = real(i,kind=qp)/real(order,kind=qp) - fetac_qp
          rd_qp = qone - fetac_qp
          !
          fmat_qp(k,k) = exp(-falpha_qp * (rn_qp/rd_qp)**Filter_Order)
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
pure function FilterMatrix_Pris(vand,vand_inv) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,np,nsp,order
  integer  :: Filter_Cutoff,ierr
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
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
  ! Compute the filter cutoff level
  !
  Filter_Cutoff = nint( max(zero, min(one,Filter_etac) ) * real(order,kind=wp) )
  !
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  fetac_qp  = real( Filter_etac , kind=qp )
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
  fmat_qp(:,:) = qzero
  do i = 1,nsp
    fmat_qp(i,i) = qone
  end do
  !
  ! Compute the filter diagonal for triangles
  !
  if (Filter_Option /= 0) then
    !
    k = 1
    do i = 0,order
      do j = 0,order-i
        !
        if (i+j >= Filter_Cutoff) then
          !
          rn_qp = real(i,kind=qp)/real(order,kind=qp) - fetac_qp
          rd_qp = qone - fetac_qp
          !
          fmat_qp(k,k) = exp(-falpha_qp * (rn_qp/rd_qp)**Filter_Order)
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
pure function FilterMatrix_Hexa(vand,vand_inv) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,np,nsp,order
  integer  :: Filter_Cutoff,ierr
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
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
  order = solpts2order_hexa(nsp)
  np = order+1
 !np = nint( real(nsp,kind=wp)**one3 )
 !order = np-1
  !
  ! Compute the filter cutoff level
  !
  Filter_Cutoff = nint( max(zero, min(one,Filter_etac) ) * real(order,kind=wp) )
  !
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  fetac_qp  = real( Filter_etac , kind=qp )
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
  ! Compute the filter diagonal for hexahedrals
  !
  if (Filter_Option == 1) then
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
        rd_qp = qone - fetac_qp
        !
        filter_1D_qp(i+1) = exp(-falpha_qp * (rn_qp/rd_qp)**Filter_Order)
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
          fmat_qp(l,l) = filter_1D_qp(i+1)*filter_1D_qp(j+1)*filter_1D_qp(k+1)
          !
        end do
      end do
    end do
    !
  else if (Filter_Option == 2) then
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
            rd_qp = real(3*order-Filter_Cutoff,kind=qp)
            !
            fmat_qp(l,l) = exp(-falpha_qp * (rn_qp/rd_qp)**Filter_Order)
            !
          end if
          !
        end do
      end do
    end do
    !
  end if
  !
  ! Compute the final filter matrix for hexahedrals
  !
  fmat_qp = matmul( real(vand,kind=qp) , fmat_qp )
  fmat_qp = matmul( fmat_qp , real(vand_inv,kind=qp) )
  !
  ! Use the chop function to copy fmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( fmat_qp )
  !
end function FilterMatrix_Hexa
!
!###############################################################################
!
subroutine filter_memory_usage(iunit)
  !
  !.. Use Statements ..
  use module_memory
  !
  !.. Formal Arguments ..
  integer, optional, intent(in) :: iunit
  !
  !.. Local Scalars ..
  integer :: iou,i,j
  character(len=36) :: array_name
  type(memory) :: array,dt
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
        array = filter(i,j)%memory_usage()
        !
        write (array_name,4) Geom_Name(i),j
        write (iou,2) array_name,array%get_units()
        !
        call dt%add(array)
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
  4 format ("filter(",a,",",i0,")%mat")
  5 format ("'",a,"'")
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  !
end subroutine filter_memory_usage
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
