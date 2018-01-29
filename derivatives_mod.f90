module derivatives_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use generic_types_mod, only : matrix
  !
  implicit none
  !
  private
  !
  !
  !.. Public Procedures ..
  !
  public :: init_derivative_matrices
  public :: derivatives_memory_usage
  !
  !
  !
  ! derivative_matrices : derived type that stores the derivative matrices of
  !                       each local coordinate direction for a given
  !                       combination of cell geometry and cell order
  !
  type, public :: derivative_matrices
    type(matrix), allocatable :: d(:)
  contains
    private
    procedure :: grad_p1 => gradient_for_point_rank1
    procedure :: grad_p2 => gradient_for_point_rank2
    procedure :: grad_c1 => gradient_for_cell_rank1
    procedure :: grad_c2 => gradient_for_cell_rank2
    procedure, public  :: div => divergence_derivative_matrices
    generic,   public  :: grad => grad_p1,grad_p2,grad_c1,grad_c2
  end type derivative_matrices
  !
  ! deriv : derivative matrix for all combinations of cell geometry and
  !        and cell order that are possible in the current simulation
  !
  !     For cell 'n' :
  !
  !         The derivative matrix (giving the contribution of each solution
  !         point within cell 'n' to the solution point 'k' that is currently
  !         of interest) is given by the expression
  !
  !         deriv(this_geom,this_order)%d(l)%mat(:,k)
  !
  !     Where :
  !
  !        this_geom = cell(n)%geom
  !        this_order = cell(n)%order
  !        l = the local coordinate direction within the reference cell,
  !            i.e., (1,2,3) for (xi,eta,zeta) or (r,s,t)
  !        k = the solution point within cell 'n' that is currently of interest
  !
  type(derivative_matrices), public, save, target, allocatable :: deriv(:,:)
  !
  ! dqmat : same as deriv but gives the derivative matrices for the quadrature
  !         points
  !
  type(derivative_matrices), public, save, target, allocatable :: dqmat(:,:)
  !
contains
!
!###############################################################################
!
pure function divergence_derivative_matrices(this,f) result(return_value)
  !
  !.. Formal Arguments ..
  class(derivative_matrices), intent(in) :: this
  real(wp),                   intent(in) :: f(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(f,dim=1))
  !
  !.. Local Scalars ..
  integer :: n,l
  !
continue
  !
  return_value = zero
  !
  do l = 1,size(this%d)
    do n = 1,size(f,dim=1)
      ! Inner loop over n means:
      !    we reuse f(:,l) and d(l) as much as possible
      !    once they are loaded into cache
      ! instead of:
      !    for n, loading f(:,l) and d(l) into cache and using it, loading
      !    f(:,l+1) and d(l+1) into cache and using it, and repeating again
      !    if size(this%d) == 3; start over for n+1 and continue until
      !    n == size(f,dim=1).
      return_value(n) = return_value(n) + this%d(l)%dot( n , f(:,l) )
    end do
  end do
  !
end function divergence_derivative_matrices
!
!###############################################################################
!
pure function gradient_for_point_rank1(this,n,f) result(return_value)
  !
  !.. Formal Arguments ..
  class(derivative_matrices), intent(in) :: this
  integer,                    intent(in) :: n
  real(wp),                   intent(in) :: f(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(this%d))
  !
  !.. Local Scalars ..
  integer :: l
  !
continue
  !
  return_value = zero
  !
  do l = 1,size(this%d)
    return_value(l) = this%d(l)%dot( n , f(:) )
  end do
  !
end function gradient_for_point_rank1
!
!###############################################################################
!
pure function gradient_for_point_rank2(this,n,f) result(return_value)
  !
  !.. Formal Arguments ..
  class(derivative_matrices), intent(in) :: this
  integer,                    intent(in) :: n
  real(wp),                   intent(in) :: f(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(this%d),1:size(f,dim=2))
  !
  !.. Local Scalars ..
  integer :: l,m
  !
continue
  !
  return_value = zero
  !
  do m = 1,size(f,dim=2)
    do l = 1,size(this%d)
      return_value(l,m) = this%d(l)%dot( n , f(:,m) )
    end do
  end do
  !
end function gradient_for_point_rank2
!
!###############################################################################
!
pure function gradient_for_cell_rank1(this,f) result(return_value)
  !
  !.. Formal Arguments ..
  class(derivative_matrices), intent(in) :: this
  real(wp),                   intent(in) :: f(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(f,dim=1),1:size(this%d))
  !
  !.. Local Scalars ..
  integer :: l
  !
continue
  !
  return_value = zero
  !
  do l = 1,size(this%d) ! Derivative reference direction
    return_value(:,l) = this%d(l)%vm_mult( f(:) )
  end do
  !
end function gradient_for_cell_rank1
!
!###############################################################################
!
pure function gradient_for_cell_rank2(this,f) result(return_value)
  !
  !.. Formal Arguments ..
  class(derivative_matrices), intent(in) :: this
  real(wp),                   intent(in) :: f(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(f,dim=1),1:size(f,dim=2),1:size(this%d))
  !
  !.. Local Scalars ..
  integer :: l,m
  !
continue
  !
  return_value = zero
  !
  do l = 1,size(this%d) ! Derivative reference direction
    do m = 1,size(f,dim=2) ! variable to differentiate in direction 'l'
      return_value(:,m,l) = this%d(l)%vm_mult( f(:,m) )
    end do
  end do
  !
end function gradient_for_cell_rank2
!
!###############################################################################
!
subroutine init_derivative_matrices()
  !
  !.. Use Statements ..
  use order_mod,      only : n_min_geom,n_max_geom
  use order_mod,      only : n_min_order,n_max_order
  use order_mod,      only : n_order
  use ovar,           only : loc_solution_pts
  use quadrature_mod, only : geom_is_used
  use quadrature_mod, only : std_elem
  !
  !.. Local Scalars ..
  integer :: gmin,gmax,omin,omax,ierr
  integer :: this_geom,this_order,ndim
  integer :: n,n1,n2,i
  character(len=200) :: array_name
  !
  !.. Local Pointers ..
  type(derivative_matrices), pointer :: this_dmat
 !real(wp), pointer, contiguous :: solpts(:,:)
  real(wp), pointer :: solpts(:,:)
 !real(wp), pointer, contiguous :: d1_mat(:,:)
  real(wp), pointer :: d1_mat(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_derivative_matrices"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Initialize the local pointers to disassociated
  !
  this_dmat => null()
  solpts    => null()
  d1_mat    => null()
  !
  gmin = n_min_geom
  gmax = n_max_geom
  omin = n_min_order
  omax = n_max_order
  !
  ! Allocate the vand array
  !
  allocate ( deriv(gmin:gmax,omin:omax) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"deriv",1,__LINE__,__FILE__,ierr,error_message)
  !
  geom_loop: do this_geom = gmin,gmax
    !
    if (.not. geom_is_used(this_geom)) cycle geom_loop
    if (all(this_geom /= Geom_Valid)) cycle geom_loop
    if (this_geom == Geom_Node) cycle geom_loop
    !
    ! Get the number of dimensions for this cell type
    !
    ndim = geom_dimen(this_geom)
    !
    order_loop: do this_order = omin,omax
      !
      ! Assign local pointer as alias to simplify the code
      !
      this_dmat => deriv(this_geom,this_order)
      !
      ! Allocate the number of dimensional derivative matrices
      !
      allocate ( this_dmat%d(1:ndim) , stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                       error_message,skip_alloc_pause)
      !
      ! Compute the derivative matrices depending on geometry type
      !
      select case (this_geom)
        case (Geom_Edge)
          !
          this_dmat%d = Dmatrix_Edge( this_order , loc_solution_pts )
          !
        case (Geom_Tria)
          !
          solpts => std_elem(this_geom,this_order)%pts
          !
          this_dmat%d = Dmatrix_Tria( solpts )
          !
        case (Geom_Quad)
          !
          d1_mat => deriv(Geom_Edge,this_order)%d(1)%mat
          !
          this_dmat%d = Dmatrix_Quad( d1_mat )
          !
        case (Geom_Tetr)
          !
          solpts => std_elem(this_geom,this_order)%pts
          !
          this_dmat%d = Dmatrix_Tetr( solpts )
          !
        case (Geom_Pyra)
          !
          solpts => std_elem(this_geom,this_order)%pts
          !
          this_dmat%d = Dmatrix_Pyra( solpts )
          !
        case (Geom_Pris)
          !
          solpts => std_elem(this_geom,this_order)%pts
          !
          this_dmat%d = Dmatrix_Pris( solpts )
          !
        case (Geom_Hexa)
          !
          d1_mat => deriv(Geom_Edge,this_order)%d(1)%mat
          !
          this_dmat%d = Dmatrix_Hexa( d1_mat )
          !
        case default
          !
          write (error_message,2)
          call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
          !
      end select
      !
      ! Remove all associations of the local pointers before continuing
      !
      if (associated(this_dmat)) this_dmat => null()
      if (associated(solpts   )) solpts    => null()
      if (associated(d1_mat   )) d1_mat    => null()
      !
    end do order_loop
    !
  end do geom_loop
  !
  ! Check for sparsity
  !
  do this_geom = gmin,gmax
    do this_order = omin,omax
      if (allocated(deriv(this_geom,this_order)%d)) then
        call deriv(this_geom,this_order)%d%check_for_sparsity
      end if
    end do
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("deriv(",a,",",i0,")%d")
  2 format (" A grid cell of an unknown geometry type was found!")
  !
end subroutine init_derivative_matrices
!
!###############################################################################
!
pure function Dmatrix_Edge(n,location) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : derivative_matrix_legendre_gauss
  use polynomial_mod, only : derivative_matrix_legendre_gauss_lobatto
  !
  !.. Formal Arguments ..
  integer, intent(in) :: n
  integer, intent(in) :: location
  !
  !.. Function Result ..
  type(matrix), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: nsp,ierr
  !
  !.. Local Arrays ..
  real(qp), dimension(1:n+1,1:n+1) :: dmat_qp
  !
continue
  !
  nsp = n+1
  !
  ! Allocate the dimensions of the derivative matrix, return_value
  !
  allocate ( return_value(1:1) , stat=ierr )
  !
  ! Allocate each individual derivative matrix within return_value
  !
  allocate ( return_value(1)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  ! Compute the derivative matrix using quad precision
  !
  select case (location)
    !
    case (Legendre_Gauss)
      !
      ! Derivative matrix for Legendre Gauss solution points
      !
      dmat_qp = derivative_matrix_legendre_gauss(n,qzero)
      !
    case (Legendre_Gauss_Lobatto)
      !
      ! Derivative matrix for Legendre Gauss Lobatto solution points
      !
      dmat_qp = derivative_matrix_legendre_gauss_lobatto(n,qzero)
      !
    case default
      !
      ! Default to Legendre Gauss solution points
      !
      dmat_qp = derivative_matrix_legendre_gauss(n,qzero)
      !
  end select
  !
  !!!!!
  !!!!! Alternative Method if this is ever needed
  !!!!!
  !!!!return_value = invert_matrix(V)
  !!!!return_value = matmul(Vr,return_value)
  !
  ! We need to transpose this derivative matrix to make its use more efficient
  !
  dmat_qp = transpose( dmat_qp )
  !
  ! Use the chop function to copy dmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value(1)%mat = chop( dmat_qp )
  !
end function Dmatrix_Edge
!
!###############################################################################
!
pure function Dmatrix_Quad(d1_mat) result(return_value)
  !
  ! Create the derivative matrices for quad elements.
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: d1_mat
  !
  !.. Function Result ..
  type(matrix), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: i,j,np,nsp,ierr
  !
  !.. Local Arrays ..
  integer, dimension(1:size(d1_mat,dim=1)) :: idx
  !
continue
  !
  ! ###########################################################################
  ! ###########################################################################
  ! ####  NOTE: No need to use chopping function in this procedure. The    ####
  ! ####        quad derivative matrix is really just the 1D derivative    ####
  ! ####        matrix applied in each direction of the standard quad      ####
  ! ####        element. The 1D derivative matrix is obtained by the       ####
  ! ####        function Dmatrix_Edge which applies the chopping function  ####
  ! ####        to its result so there is no need to repeat it again.      ####
  ! ###########################################################################
  ! ###########################################################################
  !
  np = size(d1_mat,dim=1)
  nsp = np*np
  !
  ! Allocate the dimensions of the derivative matrix, return_value
  !
  allocate ( return_value(1:2) , stat=ierr )
  !
  ! Allocate each individual derivative matrix within return_value
  !
  allocate ( return_value(1)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  allocate ( return_value(2)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  ! Create the quadrilateral derivative matrix in the x/xi direction
  !
  do j = 1,np
    !
    ! Create the indices within the derivative matrix in the x/xi direction
    ! that correspond to the locations within the 1D derivative matrix
    !
    idx = [ ( (j-1)*np+i , i=1,np ) ]
    !
    ! Copy the 1D derivative matrix to the corresponding locations given by idx
    ! NOTE: We need to copy the transpose due to the
    !       way the derivative matrices are stored
    !
   !return_value(1)%mat(idx,idx) = transpose( d1_mat )
    return_value(1)%mat(idx,idx) = d1_mat
    !
  end do
  !
  ! Create the quadrilateral derivative matrix in the y/eta direction
  !
  do i = 1,np
    !
    ! Create the indices within the derivative matrix in the y/eta direction
    ! that correspond to the locations within the 1D derivative matrix
    !
    idx = [ ( (j-1)*np+i , j=1,np ) ]
    !
    ! Copy the 1D derivative matrix to the corresponding locations given by idx
    ! NOTE: We need to copy the transpose due to the
    !       way the derivative matrices are stored
    !
   !return_value(2)%mat(idx,idx) = transpose( d1_mat )
    return_value(2)%mat(idx,idx) = d1_mat
    !
  end do
  !
end function Dmatrix_Quad
!
!###############################################################################
!
pure function Dmatrix_Tria(rs) result(return_value)
  !
  !.. Use Statements ..
  use vandermonde_mod, only : vand
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Function Result ..
  type(matrix), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: nsp,this_order,ierr
  !
  !.. Local Allocatable Arrays ..
  type(matrix), allocatable :: grad_vand(:)
  real(qp), allocatable :: dmat_qp(:,:)
  real(qp), allocatable :: grad_vand_qp(:,:)
  real(qp), allocatable :: vand_inv_qp(:,:)
  !
continue
  !
  nsp = size(rs,dim=2)
  this_order = np2n(nsp)
  !
  ! Allocate the dimensions of the derivative matrix, return_value
  !
  allocate ( return_value(1:2) , stat=ierr )
  !
  ! Allocate each individual derivative matrix within return_value
  !
  allocate ( return_value(1)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  allocate ( return_value(2)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  ! Compute the derivative of the Vandermonde matrix
  !
  allocate ( grad_vand(1:2) , stat=ierr )
  allocate ( grad_vand(1)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  allocate ( grad_vand(2)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  grad_vand(1:2) = GradVandermonde_Tria(rs)
  !
  ! Allocate the quad precision arrays needed to finish
  ! computing the derivative matrices
  !
  allocate ( dmat_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  allocate ( grad_vand_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  allocate ( vand_inv_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  !
  ! Create a quad precision copy of the inverse of the Vandermonde matrix
  !
  vand_inv_qp = real( vand(Geom_Tria,this_order)%nodal2modal%mat , kind=qp )
  !
  ! Compute the derivative matrix in the r/x/xi direction
  !
  grad_vand_qp = real( grad_vand(1)%mat , kind=qp )
  dmat_qp = qtwo * matmul( grad_vand_qp , vand_inv_qp )
  !
  ! Use the chop function to copy dmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  ! NOTE: Transpose this matrix to get it in the right ordering
  !
  return_value(1)%mat = transpose( chop( dmat_qp ) )
  !
  ! Now repeat to compute the derivative matrix in the s/y/eta direction
  !
  grad_vand_qp = real( grad_vand(2)%mat , kind=qp )
  dmat_qp = qtwo * matmul( grad_vand_qp , vand_inv_qp )
  !
  ! Use the chop function to copy dmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  ! NOTE: Transpose this matrix to get it in the right ordering
  !
  return_value(2)%mat = transpose( chop( dmat_qp ) )
  !
end function Dmatrix_Tria
!
!###############################################################################
!
pure function Dmatrix_Tetr(rs) result(return_value)
  !
  !.. Use Statements ..
  use vandermonde_mod, only : vand
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Function Result ..
  type(matrix), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: nsp,this_order,ierr
  !
  !.. Local Allocatable Arrays ..
  type(matrix), allocatable :: grad_vand(:)
  real(qp), allocatable :: dmat_qp(:,:)
  real(qp), allocatable :: grad_vand_qp(:,:)
  real(qp), allocatable :: vand_inv_qp(:,:)
  !
continue
  !
  nsp = size(rs,dim=2)
  this_order = np2n(nsp)
  !
  ! Allocate the dimensions of the derivative matrix, return_value
  !
  allocate ( return_value(1:3) , stat=ierr )
  !
  ! Allocate each individual derivative matrix within return_value
  !
  allocate ( return_value(1)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  allocate ( return_value(2)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  allocate ( return_value(3)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  ! Compute the derivative of the Vandermonde matrix
  !
  allocate ( grad_vand(1:2) , stat=ierr )
  allocate ( grad_vand(1)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  allocate ( grad_vand(2)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  grad_vand(1:2) = GradVandermonde_Tria(rs)
  !
  ! Allocate the quad precision arrays needed to finish
  ! computing the derivative matrices
  !
  allocate ( dmat_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  allocate ( grad_vand_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  allocate ( vand_inv_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  !
  ! Create a quad precision copy of the inverse of the Vandermonde matrix
  !
  vand_inv_qp = real( vand(Geom_Tria,this_order)%nodal2modal%mat , kind=qp )
  !
  ! Compute the derivative matrix in the r/x/xi direction
  !
  grad_vand_qp = real( grad_vand(1)%mat , kind=qp )
  dmat_qp = qtwo * matmul( grad_vand_qp , vand_inv_qp )
  !
  ! Use the chop function to copy dmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  ! NOTE: Transpose this matrix to get it in the right ordering
  !
  return_value(1)%mat = transpose( chop( dmat_qp ) )
  !
  ! Now repeat to compute the derivative matrix in the s/y/eta direction
  !
  grad_vand_qp = real( grad_vand(2)%mat , kind=qp )
  dmat_qp = qtwo * matmul( grad_vand_qp , vand_inv_qp )
  !
  ! Use the chop function to copy dmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  ! NOTE: Transpose this matrix to get it in the right ordering
  !
  return_value(2)%mat = transpose( chop( dmat_qp ) )
  return_value(3)%mat = transpose( chop( dmat_qp ) )
  !
end function Dmatrix_Tetr
!
!###############################################################################
!
pure function Dmatrix_Pyra(rs) result(return_value)
  !
  !.. Use Statements ..
  use vandermonde_mod, only : vand
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Function Result ..
  type(matrix), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: nsp,this_order,ierr
  !
  !.. Local Allocatable Arrays ..
  type(matrix), allocatable :: grad_vand(:)
  real(qp), allocatable :: dmat_qp(:,:)
  real(qp), allocatable :: grad_vand_qp(:,:)
  real(qp), allocatable :: vand_inv_qp(:,:)
  !
continue
  !
  nsp = size(rs,dim=2)
  this_order = np2n(nsp)
  !
  ! Allocate the dimensions of the derivative matrix, return_value
  !
  allocate ( return_value(1:3) , stat=ierr )
  !
  ! Allocate each individual derivative matrix within return_value
  !
  allocate ( return_value(1)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  allocate ( return_value(2)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  allocate ( return_value(3)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  ! Compute the derivative of the Vandermonde matrix
  !
  allocate ( grad_vand(1:2) , stat=ierr )
  allocate ( grad_vand(1)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  allocate ( grad_vand(2)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  grad_vand(1:2) = GradVandermonde_Tria(rs)
  !
  ! Allocate the quad precision arrays needed to finish
  ! computing the derivative matrices
  !
  allocate ( dmat_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  allocate ( grad_vand_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  allocate ( vand_inv_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  !
  ! Create a quad precision copy of the inverse of the Vandermonde matrix
  !
  vand_inv_qp = real( vand(Geom_Tria,this_order)%nodal2modal%mat , kind=qp )
  !
  ! Compute the derivative matrix in the r/x/xi direction
  !
  grad_vand_qp = real( grad_vand(1)%mat , kind=qp )
  dmat_qp = qtwo * matmul( grad_vand_qp , vand_inv_qp )
  !
  ! Use the chop function to copy dmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  ! NOTE: Transpose this matrix to get it in the right ordering
  !
  return_value(1)%mat = transpose( chop( dmat_qp ) )
  !
  ! Now repeat to compute the derivative matrix in the s/y/eta direction
  !
  grad_vand_qp = real( grad_vand(2)%mat , kind=qp )
  dmat_qp = qtwo * matmul( grad_vand_qp , vand_inv_qp )
  !
  ! Use the chop function to copy dmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  ! NOTE: Transpose this matrix to get it in the right ordering
  !
  return_value(2)%mat = transpose( chop( dmat_qp ) )
  return_value(3)%mat = transpose( chop( dmat_qp ) )
  !
end function Dmatrix_Pyra
!
!###############################################################################
!
pure function Dmatrix_Pris(rs) result(return_value)
  !
  !.. Use Statements ..
  use vandermonde_mod, only : vand
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Function Result ..
  type(matrix), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: nsp,this_order,ierr
  !
  !.. Local Allocatable Arrays ..
  type(matrix), allocatable :: grad_vand(:)
  real(qp), allocatable :: dmat_qp(:,:)
  real(qp), allocatable :: grad_vand_qp(:,:)
  real(qp), allocatable :: vand_inv_qp(:,:)
  !
continue
  !
  nsp = size(rs,dim=2)
  this_order = np2n(nsp)
  !
  ! Allocate the dimensions of the derivative matrix, return_value
  !
  allocate ( return_value(1:3) , stat=ierr )
  !
  ! Allocate each individual derivative matrix within return_value
  !
  allocate ( return_value(1)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  allocate ( return_value(2)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  allocate ( return_value(3)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  ! Compute the derivative of the Vandermonde matrix
  !
  allocate ( grad_vand(1:2) , stat=ierr )
  allocate ( grad_vand(1)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  allocate ( grad_vand(2)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  grad_vand(1:2) = GradVandermonde_Tria(rs)
  !
  ! Allocate the quad precision arrays needed to finish
  ! computing the derivative matrices
  !
  allocate ( dmat_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  allocate ( grad_vand_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  allocate ( vand_inv_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  !
  ! Create a quad precision copy of the inverse of the Vandermonde matrix
  !
  vand_inv_qp = real( vand(Geom_Tria,this_order)%nodal2modal%mat , kind=qp )
  !
  ! Compute the derivative matrix in the r/x/xi direction
  !
  grad_vand_qp = real( grad_vand(1)%mat , kind=qp )
  dmat_qp = qtwo * matmul( grad_vand_qp , vand_inv_qp )
  !
  ! Use the chop function to copy dmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  ! NOTE: Transpose this matrix to get it in the right ordering
  !
  return_value(1)%mat = transpose( chop( dmat_qp ) )
  !
  ! Now repeat to compute the derivative matrix in the s/y/eta direction
  !
  grad_vand_qp = real( grad_vand(2)%mat , kind=qp )
  dmat_qp = qtwo * matmul( grad_vand_qp , vand_inv_qp )
  !
  ! Use the chop function to copy dmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  ! NOTE: Transpose this matrix to get it in the right ordering
  !
  return_value(2)%mat = transpose( chop( dmat_qp ) )
  return_value(3)%mat = transpose( chop( dmat_qp ) )
  !
end function Dmatrix_Pris
!
!###############################################################################
!
pure function Dmatrix_Hexa(d1_mat) result(return_value)
  !
  ! Create the derivative matrices for hexahedral elements.
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: d1_mat
  !
  !.. Function Result ..
  type(matrix), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: i,j,k,np,nsp,ierr
  !
  !.. Local Arrays ..
  integer, dimension(1:size(d1_mat,dim=1)) :: idx
  !
continue
  !
  ! ###########################################################################
  ! ###########################################################################
  ! ####  NOTE: No need to use chopping function in this procedure. The    ####
  ! ####        hexa derivative matrix is really just the 1D derivative    ####
  ! ####        matrix applied in each direction of the standard hexa      ####
  ! ####        element. The 1D derivative matrix is obtained by the       ####
  ! ####        function Dmatrix_Edge which applies the chopping function  ####
  ! ####        to its result so there is no need to repeat it again.      ####
  ! ###########################################################################
  ! ###########################################################################
  !
  np = size(d1_mat,dim=1)
  nsp = np*np*np
  !
  ! Allocate the dimensions of the derivative matrix, return_value
  !
  allocate ( return_value(1:3) , stat=ierr )
  !
  ! Allocate each individual derivative matrix within return_value
  !
  allocate ( return_value(1)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  allocate ( return_value(2)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  allocate ( return_value(3)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  ! Create the hexahedral derivative matrix in the x/xi direction
  !
  do k = 1,np
    do j = 1,np
      !
      ! Create the indices within the derivative matrix in the x/xi direction
      ! that correspond to the locations within the 1D derivative matrix
      !
      idx = [ ( (k-1)*np*np + (j-1)*np + i , i=1,np ) ]
      !
      ! Copy the 1D derivative matrix to the
      ! corresponding locations given by idx
      ! NOTE: We need to copy the transpose due to the
      !       way the derivative matrices are stored
      !
      return_value(1)%mat(idx,idx) = d1_mat
      !
    end do
  end do
  !
  ! Create the hexahedral derivative matrix in the y/eta direction
  !
  do k = 1,np
    do i = 1,np
      !
      ! Create the indices within the derivative matrix in the y/eta direction
      ! that correspond to the locations within the 1D derivative matrix
      !
      idx = [ ( (k-1)*np*np + (j-1)*np + i , j=1,np ) ]
      !
      ! Copy the 1D derivative matrix to the
      ! corresponding locations given by idx
      ! NOTE: We need to copy the transpose due to the
      !       way the derivative matrices are stored
      !
      return_value(2)%mat(idx,idx) = d1_mat
      !
    end do
  end do
  !
  ! Create the hexahedral derivative matrix in the z/zeta direction
  !
  do j = 1,np
    do i = 1,np
      !
      ! Create the indices within the derivative matrix in the z/zeta direction
      ! that correspond to the locations within the 1D derivative matrix
      !
      idx = [ ( (k-1)*np*np + (j-1)*np + i , k=1,np ) ]
      !
      ! Copy the 1D derivative matrix to the
      ! corresponding locations given by idx
      ! NOTE: We need to copy the transpose due to the
      !       way the derivative matrices are stored
      !
      return_value(3)%mat(idx,idx) = d1_mat
      !
    end do
  end do
  !
end function Dmatrix_Hexa
!
!###############################################################################
!
pure function GradVandermonde_Edge(xi) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : normalized_jacobi_poly
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xi
  !
  !.. Function Result ..
  real(wp), dimension(1:size(xi),1:size(xi)) :: return_value
  !
  !.. Local Scalars ..
  integer  :: i,j
  real(wp) :: rj
  !
continue
  !
  return_value(:,:) = zero
  !
  do j = 2,size(xi)
    do i = 1,size(xi)
      rj = real(j-1,kind=wp)
      return_value(i,j) = sqrt(rj*(rj+one)) * &
                          normalized_jacobi_poly(j-2,one,one,xi(i))
    end do
  end do
  !
end function GradVandermonde_Edge
!
!###############################################################################
!
pure function GradVandermonde_Tria(rs) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : normalized_jacobi_poly
  use polynomial_mod, only : grad_normalized_jacobi_poly
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Function Result ..
  type(matrix), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer  :: i,im1,j,k,m,this_order,nsp,ierr
  real(qp) :: r,s,ri,twoiph,alpha,alphb,beta
  real(qp) :: a,b,fa,gb,dfa,dgb,h1mb,dmdr,dmds
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: vr_qp(:,:)
  real(qp), allocatable :: vs_qp(:,:)
  !
  !.. Local Constants ..
  real(qp), parameter :: qsqrt2 = sqrt(qtwo)
  !
continue
  !
  nsp = size(rs,dim=2)
  this_order = np2n(nsp)
  !
  ! Allocate the dimensions of the derivative of the Vandermonde matrices
  !
  allocate ( return_value(1:2) , stat=ierr )
  !
  ! Allocate each individual Vandermonde matrix derivative within return_value
  !
  allocate ( return_value(1)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  allocate ( return_value(2)%mat(1:nsp,1:nsp) , source=zero , stat=ierr )
  !
  ! Create the derivatives of the Vandermonde matrices using quad precision
  !
  allocate ( vr_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  allocate ( vs_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  !
  m = 1
  do i = 0,this_order
    !
    ri   = real( i , kind=qp )
    im1  = max(i-1,0)
    !
    twoiph = qsqrt2 * qtwo**i
    !
    alpha = qzero
    alphb = real( 2*i+1 , kind=qp )
    beta  = qzero
    !
    do j = 0,this_order-i
      !
      do k = 1,nsp
        !
        ! Get a,b from r,s
        !
        r = real( rs(1,k) , kind=qp )
        s = real( rs(2,k) , kind=qp )
        !
        a = qtwo * (qone+r) * reciprocal(qone-s) - qone
        b = s
        !
        fa = normalized_jacobi_poly(i,alpha,beta,a)
        gb = normalized_jacobi_poly(j,alphb,beta,b)
        !
        dfa = grad_normalized_jacobi_poly(i,alpha,beta,a)
        dgb = grad_normalized_jacobi_poly(j,alphb,beta,b)
        !
        h1mb = qhalf*(qone-b)
        !
        dmdr = dfa * gb * h1mb**im1
        dmds = qhalf*dmdr*(qone+a) + fa*(dgb*h1mb**i - qhalf*ri*gb*h1mb**im1)
        !
        vr_qp(k,m) = twoiph * dmdr
        vs_qp(k,m) = twoiph * dmds
        !
      end do
      !
      m = m + 1
      !
    end do
  end do
  !
  ! Use the chop function to copy vr_qp and vs_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value(1)%mat = chop( vr_qp )
  return_value(2)%mat = chop( vs_qp )
  !
end function GradVandermonde_Tria
!
!###############################################################################
!
subroutine derivatives_memory_usage(iunit)
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
  if (allocated(deriv)) then
    !
    call dt%set( storage_size(deriv,kind=inttype) , &
                         size(deriv,kind=inttype) )
    !
    do j = lbound(deriv,dim=2),ubound(deriv,dim=2)
      do i = lbound(deriv,dim=1),ubound(deriv,dim=1)
        !
        if (allocated(deriv(i,j)%d)) then
          !
          call array%set( storage_size(deriv(i,j)%d,kind=inttype) , &
                                  size(deriv(i,j)%d,kind=inttype) )
          !
          write (array_name,4) "deriv",Geom_Name(i),j
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          do k = lbound(deriv(i,j)%d,dim=1),ubound(deriv(i,j)%d,dim=1)
            !
            array = deriv(i,j)%d(k)%memory_usage()
            !
            write (array_name,4) "deriv",Geom_Name(i),j,k
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
    write (array_name,5) "deriv"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "deriv"
    !
  end if
  !
  if (allocated(dqmat)) then
    !
    call dt%set( storage_size(dqmat,kind=inttype) , &
                         size(dqmat,kind=inttype) )
    !
    do j = lbound(dqmat,dim=2),ubound(dqmat,dim=2)
      do i = lbound(dqmat,dim=1),ubound(dqmat,dim=1)
        !
        if (allocated(dqmat(i,j)%d)) then
          !
          call array%set( storage_size(dqmat(i,j)%d,kind=inttype) , &
                                  size(dqmat(i,j)%d,kind=inttype) )
          !
          write (array_name,4) "dqmat",Geom_Name(i),j
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          do k = lbound(dqmat(i,j)%d,dim=1),ubound(dqmat(i,j)%d,dim=1)
            !
            array = dqmat(i,j)%d(k)%memory_usage()
            !
            write (array_name,4) "dqmat",Geom_Name(i),j,k
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
    write (array_name,5) "dqmat"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "dqmat"
    !
  end if
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module derivatives_mod")
  2 format (" Array: ",a,dt)
  3 format (" Total memory = ",f8.2,1x,a)
  4 format (a,"(",a,",",i0,")%d",:,"(",i0,")%mat")
  5 format ("'",a,"'")
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  !
#endif
end subroutine derivatives_memory_usage
!
!###############################################################################
!
include "Functions/np2n.f90"
!
!###############################################################################
!
end module derivatives_mod
!
!###############################################################################
!
