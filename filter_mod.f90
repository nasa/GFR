module filter_mod
  !
  !.. Use Statements ..
  use module_kind_types
  !
  use generic_types_mod, only : matrix
  !
  use ovar, only : Filter_Option,Limiter_Filter_Option
  use ovar, only : Filter_etac,  Limiter_Filter_etac
  use ovar, only : Filter_Order, Limiter_Filter_Order
  use ovar, only : Filter_Alpha, Limiter_Filter_Alpha
  use ovar, only : Filter_eps,   Limiter_Filter_eps
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
  ! filter_t : derived type that contains various types of filtering matrices
  !
  type, public :: filter_t
    ! exponential : matrix to filter a polynomial by multiplying the polynomial
    !               modes by an exponentially decreasing function
    type(matrix), allocatable :: exponential
    ! limiter : filter matrix to use as a limiter for troubled cells
    type(matrix), allocatable :: limiter
    ! cutoff : matrix to filter a polynomial by
    !          cutting-off (zeroing) higher modes
    type(matrix), allocatable :: cutoff(:)
  end type filter_t
  !
  ! filter : filter matrices for all combinations of cell geometry and
  !          cell order that are possible in the current simulation
  !
  type(filter_t), public, save, target, allocatable :: filter(:,:)
  !
  ! Cutoff_Filter_Option : method used to compute the cutoff filter matrices
  !          = 1 : use a tensor product of the 1D filter with the filter cutoff
  !                based on the order of the 1D polynomial
  !          = 2 : the filter cutoff is based on total polynomial order,
  !                i.e., for the polynomial (x^a * y^b * z*c),
  !                      the total polynomial order = a+b+c
  !
  integer, parameter :: Cutoff_Filter_Option = 1
 !integer, parameter :: Cutoff_Filter_Option = 2
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
  integer :: cutoff_order
  character(len=200) :: array_name
  !
  !.. Local Pointers ..
  type(matrix), pointer :: filt
 !real(wp), pointer, contiguous :: fmat(:,:)
  real(wp), pointer :: fmat(:,:)
 !real(wp), pointer, contiguous :: vmat(:,:)
  real(wp), pointer :: vmat(:,:)
 !real(wp), pointer, contiguous :: vmat_inv(:,:)
  real(wp), pointer :: vmat_inv(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_filter_matrices"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Initialize the local pointers to disassociated
  !
  filt => null()
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
  call alloc_error(pname,"filter",1,__LINE__,__FILE__,ierr,error_message)
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
      ! Assign local pointers as aliases to simplify the code
      !
      vmat => vand(this_geom,this_order)%modal2nodal%mat
      vmat_inv => vand(this_geom,this_order)%nodal2modal%mat
      !
      ! ###################################
      ! #####  1. Exponential Filter  #####
      ! ###################################
      !
      ! Allocate the exponential filter matrix
      !
      allocate ( filter(this_geom,this_order)%exponential , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"exponential"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Assign local pointers as aliases to simplify the code
      !
      filt => filter(this_geom,this_order)%exponential
      !
      ! Allocate the exponential filter matrix
      !
      allocate ( filt%mat(1:npts,1:npts) , &
                 source=zero , stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"exponential%mat"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Assign local pointers as aliases to simplify the code
      !
      fmat => filt%mat
      !
      ! Compute the filter matrices depending on geometry type
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
          write (error_message,100)
          call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
          !
      end select
      !
      ! Remove all associations of the local pointers before continuing
      !
      if (associated(fmat)) fmat => null()
      if (associated(filt)) filt => null()
      !
      ! ###############################
      ! #####  2. Limiter Filter  #####
      ! ###############################
      !
      ! Allocate the limiter filter matrix
      !
      allocate ( filter(this_geom,this_order)%limiter , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"limiter"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Assign local pointers as aliases to simplify the code
      !
      filt => filter(this_geom,this_order)%limiter
      !
      ! Allocate the limiter filter matrix
      !
      allocate ( filt%mat(1:npts,1:npts) , &
                 source=zero , stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"limiter%mat"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Assign local pointers as aliases to simplify the code
      !
      fmat => filt%mat
      !
      ! Compute the filter matrices depending on geometry type
      !
      select case (this_geom)
        case (Geom_Edge)
          !
          fmat = FilterMatrix_Edge( vmat , vmat_inv , cutoff_order=-1 )
          !
        case (Geom_Tria)
          !
          fmat = FilterMatrix_Tria( vmat , vmat_inv , cutoff_order=-1 )
          !
        case (Geom_Quad)
          !
          fmat = FilterMatrix_Quad( vmat , vmat_inv , cutoff_order=-1 )
          !
        case (Geom_Tetr)
          !
          fmat = FilterMatrix_Tetr( vmat , vmat_inv , cutoff_order=-1 )
          !
        case (Geom_Pyra)
          !
          fmat = FilterMatrix_Pyra( vmat , vmat_inv , cutoff_order=-1 )
          !
        case (Geom_Pris)
          !
          fmat = FilterMatrix_Pris( vmat , vmat_inv , cutoff_order=-1 )
          !
        case (Geom_Hexa)
          !
          fmat = FilterMatrix_Hexa( vmat , vmat_inv , cutoff_order=-1 )
          !
        case default
          !
          write (error_message,100)
          call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
          !
      end select
      !
      ! Remove all associations of the local pointers before continuing
      !
      if (associated(fmat)) fmat => null()
      if (associated(filt)) filt => null()
      !
      ! ################################
      ! #####  3. Cut-Off Filters  #####
      ! ################################
      !
      ! Allocate the cutoff filter matrix
      !
      allocate ( filter(this_geom,this_order)%cutoff(1:this_order) , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"cutoff"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      cutoff_loop: do cutoff_order = 1,this_order
        !
        ! Assign local pointers as aliases to simplify the code
        !
        filt => filter(this_geom,this_order)%cutoff(cutoff_order)
        !
        ! Allocate the cutoff filter matrix
        !
        allocate ( filt%mat(1:npts,1:npts) , &
                   source=zero , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order, &
                             "cutoff",cutoff_order
        call alloc_error(pname,array_name,1,__LINE__,__FILE__, &
                         ierr,error_message)
        !
        ! Assign local pointers as aliases to simplify the code
        !
        fmat => filt%mat
        !
        ! Compute the filter matrices depending on geometry type
        !
        select case (this_geom)
          case (Geom_Edge)
            !
            fmat = FilterMatrix_Edge( vmat , vmat_inv , cutoff_order )
            !
          case (Geom_Tria)
            !
            fmat = FilterMatrix_Tria( vmat , vmat_inv , cutoff_order )
            !
          case (Geom_Quad)
            !
            fmat = FilterMatrix_Quad( vmat , vmat_inv , cutoff_order )
            !
          case (Geom_Tetr)
            !
            fmat = FilterMatrix_Tetr( vmat , vmat_inv , cutoff_order )
            !
          case (Geom_Pyra)
            !
            fmat = FilterMatrix_Pyra( vmat , vmat_inv , cutoff_order )
            !
          case (Geom_Pris)
            !
            fmat = FilterMatrix_Pris( vmat , vmat_inv , cutoff_order )
            !
          case (Geom_Hexa)
            !
            fmat = FilterMatrix_Hexa( vmat , vmat_inv , cutoff_order )
            !
          case default
            !
            write (error_message,100)
            call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
            !
        end select
        !
        ! Remove all associations of the local pointers before continuing
        !
        if (associated(fmat)) fmat => null()
        if (associated(filt)) filt => null()
        !
      end do cutoff_loop
      !
      ! Remove all associations of the local pointers before continuing
      !
      if (associated(fmat)) fmat => null()
      if (associated(filt)) filt => null()
      if (associated(vmat)) vmat => null()
      if (associated(vmat_inv)) vmat_inv => null()
      !
    end do order_loop
    !
  end do geom_loop
  !
  ! Check for sparsity
  !
  do this_geom = gmin,gmax
    do this_order = omin,omax
      if (allocated(filter(this_geom,this_order)%exponential)) then
        call filter(this_geom,this_order)%exponential%check_for_sparsity
      end if
      if (allocated(filter(this_geom,this_order)%limiter)) then
        call filter(this_geom,this_order)%limiter%check_for_sparsity
      end if
      if (allocated(filter(this_geom,this_order)%cutoff)) then
        call filter(this_geom,this_order)%cutoff%check_for_sparsity
      end if
    end do
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("filter(",a,",",i0,")%",a,:,"(",i0,")%mat")
  100 format (" A grid cell of an unknown geometry type was found!")
  !
end subroutine init_filter_matrices
!
!###############################################################################
!
subroutine filter_usp
  !
  !.. Use Statements ..
  use geovar,     only : ncell,cell
  use flowvar,    only : usp
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,m
  integer :: this_geom,this_order
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "filter_usp"
  !
continue
  !
 !call debug_timer(entering_procedure,pname)
  !
  ! Apply the filter to each cell
  !
  do n = 1,ncell
    !
    n1 = cell(n)%beg_sp
    n2 = cell(n)%end_sp
    !
    this_geom = cell(n)%geom
    this_order = cell(n)%order
    !
    do m = 1,size(usp,dim=1)
      usp(m,n1:n2) = filter(this_geom,this_order)% &
                                      exponential% &
                     mv_mult( usp(m,n1:n2) )
    end do
    !
  end do
  !
 !call debug_timer(leaving_procedure,pname)
  !
end subroutine filter_usp
!
!###############################################################################
!
pure function FilterMatrix_Edge(vand,vand_inv,cutoff_order) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Optional Arguments ..
  integer,  optional, intent(in) :: cutoff_order
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,nsp,order
  integer  :: Filter_Cutoff,ierr
  integer  :: local_Filter_Option
  integer  :: filt_order
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
  real(qp) :: feps_qp
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
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  !
  filt_order = Filter_Order
  !
  if (present(cutoff_order)) then
    if (cutoff_order >= 0) then
      local_Filter_Option = Cutoff_Filter_Option
      Filter_Cutoff = cutoff_order
      fetac_qp = real(cutoff_order,kind=qp) / real(order,kind=qp)
      feps_qp = qzero
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = nint( max(zero,min(one,Limiter_Filter_etac)) * &
                      real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp = real( Filter_etac , kind=qp )
    Filter_Cutoff = nint( max(zero,min(one,Filter_etac)) * real(order,kind=wp) )
    local_Filter_Option = Filter_Option
    feps_qp = real( one - Filter_eps , kind=qp )
  end if
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
  if (local_Filter_Option /= 0) then
    !
    do i = 0,order
      if (i >= Filter_Cutoff) then
        !
        rn_qp = real(i,kind=qp)/real(order,kind=qp) - fetac_qp
        rd_qp = max( qzero , qone - fetac_qp )
        !
        if (rd_qp == qzero) then
          fmat_qp(i+1,i+1) = feps_qp
        else
          fmat_qp(i+1,i+1) = feps_qp*exp(-falpha_qp*(rn_qp/rd_qp)**filt_order)
        end if
        !
      end if
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
end function FilterMatrix_Edge
!
!###############################################################################
!
 pure &
function FilterMatrix_Quad(vand,vand_inv,cutoff_order) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: cutoff_order
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,np,nsp,order
  integer  :: Filter_Cutoff,ierr
  integer  :: local_Filter_Option
  integer  :: filt_order
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
  real(qp) :: feps_qp
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
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  !
  filt_order = Filter_Order
  !
  if (present(cutoff_order)) then
    if (cutoff_order >= 0) then
      local_Filter_Option = Cutoff_Filter_Option
      i = merge(2,1,local_Filter_Option==2)
      Filter_Cutoff = i*cutoff_order
      fetac_qp = real(Filter_Cutoff,kind=qp) / real(i*order,kind=qp)
      feps_qp = qzero
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = nint( max(zero,min(one,Limiter_Filter_etac)) * &
                      real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp  = real( Filter_etac , kind=qp )
    Filter_Cutoff = nint( max(zero,min(one,Filter_etac)) * real(order,kind=wp) )
    local_Filter_Option = Filter_Option
    feps_qp = real( one - Filter_eps , kind=qp )
  end if
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
  if (local_Filter_Option == 1) then
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
        rd_qp = max( qzero , qone - fetac_qp )
        !
        if (rd_qp == qzero) then
          filter_1D_qp(i+1) = feps_qp
        else
          filter_1D_qp(i+1) = feps_qp * &
                              exp(-falpha_qp * (rn_qp/rd_qp)**filt_order)
        end if
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
  else if (local_Filter_Option == 2) then
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
          rd_qp = max( qzero , real(2*order-Filter_Cutoff,kind=qp) )
          !
          if (rd_qp == qzero) then
            fmat_qp(k,k) = feps_qp
          else
            fmat_qp(k,k) = feps_qp*exp(-falpha_qp * (rn_qp/rd_qp)**filt_order)
          end if
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
pure function FilterMatrix_Tria(vand,vand_inv,cutoff_order) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: cutoff_order
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,np,nsp,order
  integer  :: Filter_Cutoff,ierr
  integer  :: local_Filter_Option
  integer  :: filt_order
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
  real(qp) :: feps_qp
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
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  !
  filt_order = Filter_Order
  !
  if (present(cutoff_order)) then
    if (cutoff_order >= 0) then
      local_Filter_Option = Cutoff_Filter_Option
      Filter_Cutoff = cutoff_order
      fetac_qp = real(Filter_Cutoff,kind=qp) / real(order,kind=qp)
      feps_qp = qzero
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = nint( max(zero,min(one,Limiter_Filter_etac)) * &
                      real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp  = real( Filter_etac , kind=qp )
    Filter_Cutoff = nint( max(zero,min(one,Filter_etac)) * real(order,kind=wp) )
    local_Filter_Option = Filter_Option
    feps_qp = real( one - Filter_eps , kind=qp )
  end if
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
  ! Compute the filter diagonal for triangles
  !
  if (local_Filter_Option /= 0) then
    !
    k = 1
    do i = 0,order
      do j = 0,order-i
        !
        if (i+j >= Filter_Cutoff) then
          !
          rn_qp = real(i,kind=qp)/real(order,kind=qp) - fetac_qp
          rd_qp = max( qzero , qone - fetac_qp )
          !
          if (rd_qp == qzero) then
            fmat_qp(k,k) = feps_qp
          else
            fmat_qp(k,k) = feps_qp*exp(-falpha_qp * (rn_qp/rd_qp)**filt_order)
          end if
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
pure function FilterMatrix_Tetr(vand,vand_inv,cutoff_order) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: cutoff_order
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,np,nsp,order
  integer  :: Filter_Cutoff,ierr
  integer  :: local_Filter_Option
  integer  :: filt_order
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
  real(qp) :: feps_qp
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
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  !
  filt_order = Filter_Order
  !
  if (present(cutoff_order)) then
    if (cutoff_order >= 0) then
      local_Filter_Option = Cutoff_Filter_Option
      Filter_Cutoff = cutoff_order
      fetac_qp = real(Filter_Cutoff,kind=qp) / real(order,kind=qp)
      feps_qp = qzero
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = nint( max(zero,min(one,Limiter_Filter_etac)) * &
                      real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp  = real( Filter_etac , kind=qp )
    Filter_Cutoff = nint( max(zero,min(one,Filter_etac)) * real(order,kind=wp) )
    local_Filter_Option = Filter_Option
    feps_qp = real( one - Filter_eps , kind=qp )
  end if
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
  ! Compute the filter diagonal for triangles
  !
  if (local_Filter_Option /= 0) then
    !
    k = 1
    do i = 0,order
      do j = 0,order-i
        !
        if (i+j >= Filter_Cutoff) then
          !
          rn_qp = real(i,kind=qp)/real(order,kind=qp) - fetac_qp
          rd_qp = max( qzero , qone - fetac_qp )
          !
          if (rd_qp == qzero) then
            fmat_qp(k,k) = feps_qp
          else
            fmat_qp(k,k) = feps_qp*exp(-falpha_qp * (rn_qp/rd_qp)**filt_order)
          end if
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
pure function FilterMatrix_Pyra(vand,vand_inv,cutoff_order) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: cutoff_order
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,np,nsp,order
  integer  :: Filter_Cutoff,ierr
  integer  :: local_Filter_Option
  integer  :: filt_order
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
  real(qp) :: feps_qp
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
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  !
  filt_order = Filter_Order
  !
  if (present(cutoff_order)) then
    if (cutoff_order >= 0) then
      local_Filter_Option = Cutoff_Filter_Option
      Filter_Cutoff = cutoff_order
      fetac_qp = real(Filter_Cutoff,kind=qp) / real(order,kind=qp)
      feps_qp = qzero
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = nint( max(zero,min(one,Limiter_Filter_etac)) * &
                      real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp  = real( Filter_etac , kind=qp )
    Filter_Cutoff = nint( max(zero,min(one,Filter_etac)) * real(order,kind=wp) )
    local_Filter_Option = Filter_Option
    feps_qp = real( one - Filter_eps , kind=qp )
  end if
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
  ! Compute the filter diagonal for triangles
  !
  if (local_Filter_Option /= 0) then
    !
    k = 1
    do i = 0,order
      do j = 0,order-i
        !
        if (i+j >= Filter_Cutoff) then
          !
          rn_qp = real(i,kind=qp)/real(order,kind=qp) - fetac_qp
          rd_qp = max( qzero , qone - fetac_qp )
          !
          if (rd_qp == qzero) then
            fmat_qp(k,k) = feps_qp
          else
            fmat_qp(k,k) = feps_qp*exp(-falpha_qp * (rn_qp/rd_qp)**filt_order)
          end if
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
pure function FilterMatrix_Pris(vand,vand_inv,cutoff_order) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: cutoff_order
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,np,nsp,order
  integer  :: Filter_Cutoff,ierr
  integer  :: local_Filter_Option
  integer  :: filt_order
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
  real(qp) :: feps_qp
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
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  !
  filt_order = Filter_Order
  !
  if (present(cutoff_order)) then
    if (cutoff_order >= 0) then
      local_Filter_Option = Cutoff_Filter_Option
      Filter_Cutoff = cutoff_order
      fetac_qp = real(Filter_Cutoff,kind=qp) / real(order,kind=qp)
      feps_qp = qzero
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = nint( max(zero,min(one,Limiter_Filter_etac)) * &
                      real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp  = real( Filter_etac , kind=qp )
    Filter_Cutoff = nint( max(zero,min(one,Filter_etac)) * real(order,kind=wp) )
    local_Filter_Option = Filter_Option
    feps_qp = real( one - Filter_eps , kind=qp )
  end if
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
  ! Compute the filter diagonal for triangles
  !
  if (local_Filter_Option /= 0) then
    !
    k = 1
    do i = 0,order
      do j = 0,order-i
        !
        if (i+j >= Filter_Cutoff) then
          !
          rn_qp = real(i,kind=qp)/real(order,kind=qp) - fetac_qp
          rd_qp = max( qzero , qone - fetac_qp )
          !
          if (rd_qp == qzero) then
            fmat_qp(k,k) = feps_qp
          else
            fmat_qp(k,k) = feps_qp*exp(-falpha_qp * (rn_qp/rd_qp)**filt_order)
          end if
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
pure function FilterMatrix_Hexa(vand,vand_inv,cutoff_order) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vand
  real(wp), dimension(:,:), intent(in) :: vand_inv
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: cutoff_order
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,np,nsp,order
  integer  :: Filter_Cutoff,ierr
  integer  :: local_Filter_Option
  integer  :: filt_order
  real(qp) :: rd_qp,fetac_qp
  real(qp) :: rn_qp,falpha_qp
  real(qp) :: feps_qp
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
  ! Covert Filter_Alpha and Filter_etac to quad precision
  !
  falpha_qp = real( Filter_Alpha , kind=qp )
  !
  filt_order = Filter_Order
  !
  if (present(cutoff_order)) then
    if (cutoff_order >= 0) then
      local_Filter_Option = Cutoff_Filter_Option
      i = merge(3,1,local_Filter_Option==2)
      Filter_Cutoff = i*cutoff_order
      fetac_qp = real(Filter_Cutoff,kind=qp) / real(i*order,kind=qp)
      feps_qp = qzero
    else
      fetac_qp = real( Limiter_Filter_etac , kind=qp )
      Filter_Cutoff = nint( max(zero,min(one,Limiter_Filter_etac)) * &
                      real(order,kind=wp) )
      local_Filter_Option = Limiter_Filter_Option
      feps_qp = real( one - Limiter_Filter_eps , kind=qp )
      falpha_qp = real( Limiter_Filter_Alpha , kind=qp )
      filt_order = Limiter_Filter_Order
    end if
  else
    fetac_qp  = real( Filter_etac , kind=qp )
    Filter_Cutoff = nint( max(zero,min(one,Filter_etac)) * real(order,kind=wp) )
    local_Filter_Option = Filter_Option
    feps_qp = real( one - Filter_eps , kind=qp )
  end if
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
  if (local_Filter_Option == 1) then
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
        rd_qp = max( qzero , qone - fetac_qp )
        !
        if (rd_qp == qzero) then
          filter_1D_qp(i+1) = feps_qp
        else
          filter_1D_qp(i+1) = feps_qp * &
                              exp(-falpha_qp * (rn_qp/rd_qp)**filt_order)
        end if
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
  else if (local_Filter_Option == 2) then
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
            rd_qp = max( qzero , real(3*order-Filter_Cutoff,kind=qp) )
            !
            if (rd_qp == qzero) then
              fmat_qp(l,l) = feps_qp
            else
              fmat_qp(l,l) = feps_qp*exp(-falpha_qp*(rn_qp/rd_qp)**filt_order)
            end if
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
  integer :: iou,i,j,k
  character(len=36) :: array_name
  type(memory) :: array,dt
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
        if (allocated(filter(i,j)%exponential)) then
          !
          call array%set( storage_size(filter(i,j)%exponential,kind=inttype) )
          !
          write (array_name,4) Geom_Name(i),j,"exponential"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          array = filter(i,j)%exponential%memory_usage()
          !
          write (array_name,4) Geom_Name(i),j,"exponential%mat"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
        end if
        !
        if (allocated(filter(i,j)%cutoff)) then
          !
          call array%set( storage_size(filter(i,j)%cutoff,kind=inttype) , &
                                  size(filter(i,j)%cutoff,kind=inttype) )
          !
          write (array_name,4) Geom_Name(i),j,"cutoff"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          do k = lbound(filter(i,j)%cutoff,dim=1), &
                 ubound(filter(i,j)%cutoff,dim=1)
            !
            array = filter(i,j)%cutoff(k)%memory_usage()
            !
            write (array_name,4) Geom_Name(i),j,"cutoff",k
            write (iou,2) array_name,array%get_units()
            !
            call dt%add(array)
            !
          end do
          !
        end if
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
  4 format ("filter(",a,",",i0,")",a,:,"(",i0,")%mat")
  5 format ("'",a,"'")
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  !
#endif
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
