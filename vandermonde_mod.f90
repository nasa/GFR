module vandermonde_mod
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
  public :: init_vandermonde_matrices
  public :: Vandermonde_Edge
  public :: Vandermonde_Tria
  public :: vandermonde_memory_usage
  public :: Compute_Basis_Matrix
  !
  !
  !
  ! vandermonde_matrix : derived type that stores the matrices for converting
  !                      the polynomial solution within a cell between the nodal
  !                      and modal forms (two matrices, one for converting in
  !                      each of the two directions) for a given combination of
  !                      cell geometry and cell order.
  !
  ! vandermonde_matrix%nodal2modal : converts a nodal solution to modal form
  ! vandermonde_matrix%modal2nodal : converts a modal solution to nodal form
  !
  ! NOTE: the modal2nodal matrix is the 'actual' generalized Vandermonde matrix
  !       whereas the nodal2modal matrix is the inverse of the generalized
  !       Vandermonde matrix
  !
  type, public :: vandermonde_matrix
    type(matrix) :: nodal2modal
    type(matrix) :: modal2nodal
  end type vandermonde_matrix
  !
  ! vand : projection matrix for all combinations of cell geometry and
  !        and cell order that are possible in the current simulation
  !
  !     For cell 'n' :
  !
  !         The expression for converting the nodal solution in cell 'n' to
  !         its modal form is
  !
  !     matmul( vand(this_geom,this_order)%nodal2modal%mat , nodal_solution )
  !
  !         The expression for converting the modal solution in cell 'n' to
  !         its nodal form is
  !
  !     matmul( vand(this_geom,this_order)%modal2nodal%mat , modal_solution )
  !
  !     Where :
  !        nsp = # of solution points within this cell geometry determined
  !              by the order of the cell
  !        this_geom = cell(n)%geom
  !        this_order = cell(n)%order
  !        mat has dimensions and extents mat(1:nsp,1:nsp)
  !        nodal_solution(1:nsp) : nodal solution vector for cell 'n'
  !        modal_solution(1:nsp) : modal solution vector for cell 'n'
  !
  type(vandermonde_matrix), public, save, target, allocatable :: vand(:,:)
  !
  ! monomials : generalized Vandermonde matrix to convert a polynomial between
  !             the nodal coefficients and the monomial modal coefficients.
  !
  type(vandermonde_matrix), public, save, target, allocatable :: monomials(:,:)
  !
  !
contains
!
!###############################################################################
!
subroutine init_vandermonde_matrices()
  !
  !.. Use Statements ..
  use order_mod,      only : geom_solpts
  use order_mod,      only : n_min_geom,n_max_geom
  use order_mod,      only : n_min_order,n_max_order
  use quadrature_mod, only : geom_is_used
  use quadrature_mod, only : std_elem
  !
  !.. Local Scalars ..
  integer :: gmin,gmax,omin,omax,ierr
  integer :: this_geom,this_order,npts
  character(len=200) :: array_name
  !
  !.. Local Pointers ..
  type(vandermonde_matrix), pointer :: this_vand
 !real(wp), pointer, contiguous :: solpts(:,:)
 !real(wp), pointer :: solpts(:,:)
 !real(wp), pointer, contiguous :: vmat(:,:)
  real(wp), pointer :: vmat(:,:)
 !real(wp), pointer, contiguous :: vmat_inv(:,:)
  real(wp), pointer :: vmat_inv(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_vandermonde_matrices"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Initialize the local pointers to disassociated
  !
  this_vand => null()
 !solpts => null()
  vmat => null()
  vmat_inv => null()
  !
  gmin = n_min_geom
  gmax = n_max_geom
  omin = n_min_order
  omax = n_max_order
  !
  ! Allocate the vand array
  !
  allocate ( vand(gmin:gmax,omin:omax) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"vand",1,__LINE__,__FILE__,ierr,error_message)
  !
  geom_loop: do this_geom = gmin,gmax
    !
    if (.not. geom_is_used(this_geom)) cycle geom_loop
    if (all(this_geom /= Geom_Valid)) cycle geom_loop
    if (this_geom == Geom_Node) cycle geom_loop
    !
    order_loop: do this_order = omin,omax
      !
      ! Assign local pointer as alias to simplify the code
      !
      this_vand => vand(this_geom,this_order)
      !
      ! Get the number of solution points for this cell
      !
      npts = geom_solpts(this_geom,this_order)
      !
      ! Allocate the modal2nodal / vandermonde matrix
      !
      allocate ( this_vand%modal2nodal%mat(1:npts,1:npts) , &
                 source=zero , stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"modal2nodal%mat"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Allocate the nodal2modal / vandermonde inverse matrix
      !
      allocate ( this_vand%nodal2modal%mat(1:npts,1:npts) , &
                 source=zero , stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"nodal2modal%mat"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Assign local pointers to these matrices as aliases to simplify the code
      !
      vmat => this_vand%modal2nodal%mat
      vmat_inv => this_vand%nodal2modal%mat
      !
      ! Compute the Vandermonde matrix
      !
      vmat(:,:) = Compute_Basis_Matrix(this_geom,this_order)
      !
      ! Compute the inverse of the Vandermonde matrix that was just computed
      !
      vmat_inv(:,:) = Vandermonde_Inverse( vmat )
      !
      ! Remove all associations of the local pointers before continuing
      !
      if (associated(this_vand)) this_vand => null()
      if (associated(vmat     )) vmat      => null()
      if (associated(vmat_inv )) vmat_inv  => null()
      !
    end do order_loop
    !
  end do geom_loop
  !
  call vand%modal2nodal%check_for_sparsity
  call vand%nodal2modal%check_for_sparsity
 !do this_geom = gmin,gmax
 !  do this_order = omin,omax
 !    call vand(this_geom,this_order)%modal2nodal%check_for_sparsity
 !    call vand(this_geom,this_order)%nodal2modal%check_for_sparsity
 !  end do
 !end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("vand(",a,",",i0,")%",a)
  2 format (" A grid cell of an unknown geometry type was found!")
  !
end subroutine init_vandermonde_matrices
!
!###############################################################################
!
function Compute_Basis_Matrix(elem_geom,elem_order,basis_order,nodal_points) &
                       result(return_value)
  !
  !.. Use Statements ..
  use order_mod,      only : geom_solpts
  use quadrature_mod, only : std_elem
  !
  !.. Formal Arguments ..
  integer, intent(in) :: elem_geom
  integer, intent(in) :: elem_order
  !
  !.. Optional Arguments ..
  integer,  optional,         intent(in) :: basis_order
  real(wp), optional, target, intent(in) :: nodal_points(:,:)
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: npts,nmodes,ierr
  !
  !.. Local Pointers ..
 !real(wp), pointer, contiguous :: solpts(:,:)
  real(wp), pointer :: solpts(:,:)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: edgpts(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "Compute_Basis_Matrix"
  !
continue
  !
  solpts => null()
  !
  npts = geom_solpts(elem_geom,elem_order)
  if (present(basis_order)) then
    nmodes = geom_solpts(elem_geom,basis_order)
  else
    nmodes = npts
  end if
  !
  allocate ( return_value(1:npts,1:nmodes) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"return_value",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Compute the basis matrix depending on the geometry type
  !
  select case (elem_geom)
    case (Geom_Edge)
      !
      if (present(nodal_points)) then
        edgpts = nodal_points(1,:) ! F2003 AUTO-REALLOCATION
       !solpts => nodal_points
      else
        edgpts = std_elem(Geom_Edge,elem_order)%pts(1,:) ! F2003 AUTO-REALLOC
       !solpts => std_elem(elem_geom,elem_order)%pts
      end if
      !
      return_value(:,:) = Vandermonde_Edge( edgpts , basis_order )
     !return_value(:,:) = Vandermonde_Edge( solpts(1,:) , basis_order )
      !
    case (Geom_Tria)
      !
      if (present(nodal_points)) then
        solpts => nodal_points
      else
        solpts => std_elem(elem_geom,elem_order)%pts
      end if
      !
      return_value(:,:) = Vandermonde_Tria( solpts , basis_order )
      !
    case (Geom_Quad)
      !
      if (present(nodal_points)) then
        edgpts = nodal_points(1,:) ! F2003 AUTO-REALLOCATION
       !solpts => nodal_points
      else
        edgpts = std_elem(Geom_Edge,elem_order)%pts(1,:) ! F2003 AUTO-REALLOC
       !solpts => std_elem(Geom_Edge,elem_order)%pts
      end if
      !
      return_value(:,:) = Vandermonde_Quad( edgpts , basis_order )
     !return_value(:,:) = Vandermonde_Quad( solpts(1,:) , basis_order )
      !
    case (Geom_Tetr)
      !
      if (present(nodal_points)) then
        solpts => nodal_points
      else
        solpts => std_elem(elem_geom,elem_order)%pts
      end if
      !
      return_value(:,:) = Vandermonde_Tetr( solpts , basis_order )
      !
    case (Geom_Pyra)
      !
      if (present(nodal_points)) then
        solpts => nodal_points
      else
        solpts => std_elem(elem_geom,elem_order)%pts
      end if
      !
      return_value(:,:) = Vandermonde_Pyra( solpts , basis_order )
      !
    case (Geom_Pris)
      !
      if (present(nodal_points)) then
        solpts => nodal_points
      else
        solpts => std_elem(elem_geom,elem_order)%pts
      end if
      !
      return_value(:,:) = Vandermonde_Pris( solpts , basis_order )
      !
    case (Geom_Hexa)
      !
      if (present(nodal_points)) then
        edgpts = nodal_points(1,:) ! F2003 AUTO-REALLOCATION
       !solpts => nodal_points
      else
        edgpts = std_elem(Geom_Edge,elem_order)%pts(1,:) ! F2003 AUTO-REALLOC
       !solpts => std_elem(Geom_Edge,elem_order)%pts
      end if
      !
      return_value(:,:) = Vandermonde_Hexa( edgpts , basis_order )
     !return_value(:,:) = Vandermonde_Hexa( solpts(1,:) , basis_order )
      !
    case default
      !
      write (error_message,1)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      !
  end select
  !
  if (associated(solpts)) solpts => null()
  !
  ! Format Statements
  !
  1 format (" An error occured by attempting to create a Basis", &
            " Matrix for a grid cell of unknown geometry type!")
  !
end function Compute_Basis_Matrix
!
!###############################################################################
!
pure function Vandermonde_Inverse(vand_mat) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: vand_mat(:,:)
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: n1,n2,ierr
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: vmat_inv_qp(:,:)
  !
continue
  !
  n1 = size(vand_mat,dim=1)
  n2 = size(vand_mat,dim=2)
  !
  ! Allocate the return array containing the working
  ! precision inverse of the Vandermonde matrix
  !
  allocate ( return_value(1:n1,1:n2) , source=zero , stat=ierr )
  !
  ! Create the Vandermonde matrix and its inverse using quad precision
  !
  allocate ( vmat_inv_qp(1:n1,1:n2) , source=qzero , stat=ierr )
  !
  ! Convert the Vandermonde matrix to quad precision and temporarily
  ! store in the array for the inverse Vandermonde matrix
  !
  vmat_inv_qp = real( vand_mat , kind=qp )
  !
  ! Compute the inverse of the Vandermonde matrix
  !
  vmat_inv_qp = invert_matrix( vmat_inv_qp )
  !
  ! Use the chop function to copy vmat_inv_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( vmat_inv_qp )
  !
end function Vandermonde_Inverse
!
!###############################################################################
!
pure function Vandermonde_Edge(xi,number_of_modes) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : normalized_jacobi_poly
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xi(:)
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: number_of_modes
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars
  integer  :: i,j,nqp,nqo,nsp,nso,ierr
  real(qp) :: x
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: vmat_qp(:,:)
  !
continue
  !
  nqp = size(xi)
  nqo = nqp-1
  !
  if (present(number_of_modes)) then
    nso = number_of_modes
  else
    nso = nqo
  end if
  nsp = nso+1
  !
  ! Allocate the return array containing the
  ! working precision Vandermonde matrix
  !
  allocate ( return_value(1:nqp,1:nsp) , source=zero , stat=ierr )
  !
  ! Create the Vandermonde matrix using quad precision
  !
  allocate ( vmat_qp(1:nqp,1:nsp) , source=qzero , stat=ierr )
  !
  do j = 1,nsp
    do i = 1,nqp
      !
      x = real( xi(i) , kind=qp )
      !
      vmat_qp(i,j) = normalized_jacobi_poly(j-1,qzero,qzero,x)
      !
    end do
  end do
  !
  ! Use the chop function to copy vmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( vmat_qp )
  !
end function Vandermonde_Edge
!
!###############################################################################
!
pure function Vandermonde_Quad(xi,number_of_modes) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : normalized_jacobi_poly
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xi(:)
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: number_of_modes
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,ni,nj,m,mi,mj,ierr
  integer  :: nsp,nso,nsq,nqp,nqo,nqq
  real(qp) :: x,y,hx,hy
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: vmat_qp(:,:)
  !
continue
  !
  nqp = size(xi)
  nqo = nqp-1
  nqq = nqp*nqp
  !
  if (present(number_of_modes)) then
    nso = number_of_modes
  else
    nso = nqo
  end if
  nsp = nso+1
  nsq = nsp*nsp
  !
  ! Allocate the return array containing the
  ! working precision Vandermonde matrix
  !
  allocate ( return_value(1:nqq,1:nsq) , source=zero , stat=ierr )
  !
  ! Create the Vandermonde matrix using quad precision
  !
  allocate ( vmat_qp(1:nqq,1:nsq) , source=qzero , stat=ierr )
  !
  do nj = 1,nqp
    do ni = 1,nqp
      !
      ! n is the 1D index of the solution point at
      ! which each basis function is being evaluated
      !
      n = nqp*(nj-1) + ni
      !
      x = real( xi(ni) , kind=qp )
      y = real( xi(nj) , kind=qp )
      !
      do mj = 0,nso
        do mi = 0,nso
          !
          ! m is the index of the current basis function being evaluated
          !
          m = nsp*mj + mi + 1
          !
          hx = normalized_jacobi_poly(mi,qzero,qzero,x)
          hy = normalized_jacobi_poly(mj,qzero,qzero,y)
          !
          vmat_qp(n,m) = hx*hy
          !
        end do
      end do
      !
    end do
  end do
  !
  ! Use the chop function to copy vmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( vmat_qp )
  !
end function Vandermonde_Quad
!
!###############################################################################
!
pure function Vandermonde_Tria(rs,number_of_modes) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : normalized_jacobi_poly
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: rs(:,:)
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: number_of_modes
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,n,np,m,ierr
  real(qp) :: r,s,a,b,h1,h2,alpha,beta
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: vmat_qp(:,:)
  !
continue
  !
  np = size(rs,2)
  n  = np2n(np)
  !
  ! Allocate the return array containing the
  ! working precision Vandermonde matrix
  !
  allocate ( return_value(1:np,1:np) , source=zero , stat=ierr )
  !
  ! Create the Vandermonde matrix using quad precision
  !
  allocate ( vmat_qp(1:np,1:np) , source=qzero , stat=ierr )
  !
  m = 1
  do i = 0,n
    do j = 0,n-i
      !
      ! Evaluate the basis function, psi(m), at the np solution points
      ! where m = j + i*(n+1) + 1 - i*(i-1)/2
      !
      do k = 1,np
        !
        ! Get a,b from r,s
        !
        r = real( rs(1,k) , kind=qp )
        s = real( rs(2,k) , kind=qp )
        !
        a = qtwo * (qone+r) * reciprocal(qone-s) - qone
        b = s
        !
        ! Get the value of the Jacobi polynomial for a and b
        !
        alpha = qzero
        beta  = qzero
        h1 = normalized_jacobi_poly(i,alpha,beta,a)
        !
        alpha = real( 2*i+1 , kind=qp )
        beta  = qzero
        h2 = normalized_jacobi_poly(j,alpha,beta,b)
        !
        ! Evaluate the basis function psi(m) at node k
        !
        vmat_qp(k,m) = sqrt(qtwo) * h1 * h2 * (qone-b)**i
        !
      end do
      !
      m = m + 1
      !
    end do
  end do
  !
  ! Use the chop function to copy vmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( vmat_qp )
  !
end function Vandermonde_Tria
!
!###############################################################################
!
pure function Vandermonde_Tetr(rs,number_of_modes) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : normalized_jacobi_poly
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: rs(:,:)
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: number_of_modes
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,n,np,m,ierr
  real(qp) :: r,s,t,a,b,c,h1,h2,h3,alpha,beta
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: vmat_qp(:,:)
  !
continue
  !
  np = size(rs,dim=2)
  !
  ! Find the order for this tetrahedra
  !
  n = -1
  do i = 0,256
    if (Tetr_Solpts(i) == np) then
      n = i
      exit
    end if
  end do
  !
  ! Allocate the return array containing the
  ! working precision Vandermonde matrix
  !
  allocate ( return_value(1:np,1:np) , source=zero , stat=ierr )
  !
  ! If n is still -1, return with the Vandermonde matrix still being zero
  !
  if (n == -1) return
  !
  ! Create the Vandermonde matrix using quad precision
  !
  allocate ( vmat_qp(1:np,1:np) , source=qzero , stat=ierr )
  !
  m = 1
  do i = 0,n
    do j = 0,n-i
      do k = 0,n-i-j
        !
        ! Evaluate the basis function, psi(m), at the np solution points
        ! where m = j + i*(n+1) + 1 - i*(i-1)/2
        !
        do l = 1,np
          !
          ! Get a,b from r,s
          !
          r = real( rs(1,l) , kind=qp )
          s = real( rs(2,l) , kind=qp )
          t = real( rs(3,l) , kind=qp )
          !
          a = qtwo * (qone+r) * reciprocal(-s-t) - qone
          b = qtwo * (qone+s) * reciprocal(qone-t) - qone
          c = t
          !
          ! Get the value of the Jacobi polynomial for a and b
          !
          alpha = qzero
          beta  = qzero
          h1 = normalized_jacobi_poly(i,alpha,beta,a)
          !
          alpha = real( 2*i+1 , kind=qp )
          beta  = qzero
          h2 = normalized_jacobi_poly(j,alpha,beta,b)
          !
          alpha = real( 2*(i+j)+2 , kind=qp )
          beta  = qzero
          h3 = normalized_jacobi_poly(k,alpha,beta,c)
          !
          ! Evaluate the basis function psi(m) at node l
          !
          vmat_qp(l,m) = qtwo*sqrt(qtwo) * h1 * h2 * h3 * (qone-b)**(i  ) &
                                                        * (qone-c)**(i+j)
          !
        end do
        !
        m = m + 1
        !
      end do
    end do
  end do
  !
  ! Use the chop function to copy vmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( vmat_qp )
  !
end function Vandermonde_Tetr
!
!###############################################################################
!
pure function Vandermonde_Pyra(rs,number_of_modes) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : normalized_jacobi_poly
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: rs(:,:)
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: number_of_modes
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,n,np,m,ierr
  real(qp) :: r,s,a,b,h1,h2,alpha,beta
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: vmat_qp(:,:)
  !
continue
  !
  np = size(rs,2)
  n  = np2n(np)
  !
  ! Allocate the return array containing the
  ! working precision Vandermonde matrix
  !
  allocate ( return_value(1:np,1:np) , source=zero , stat=ierr )
  !
end function Vandermonde_Pyra
!
!###############################################################################
!
pure function Vandermonde_Pris(rs,number_of_modes) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : normalized_jacobi_poly
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: rs(:,:)
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: number_of_modes
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,n,np,m,ierr
  real(qp) :: r,s,a,b,h1,h2,alpha,beta
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: vmat_qp(:,:)
  !
continue
  !
  np = size(rs,2)
  n  = np2n(np)
  !
  ! Allocate the return array containing the
  ! working precision Vandermonde matrix
  !
  allocate ( return_value(1:np,1:np) , source=zero , stat=ierr )
  !
end function Vandermonde_Pris
!
!###############################################################################
!
pure function Vandermonde_Hexa(xi,number_of_modes) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : normalized_jacobi_poly
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xi(:)
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: number_of_modes
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,ni,nj,nk,m,mi,mj,mk,ierr
  integer  :: nsp,nso,nsq,nqp,nqo,nqq
  real(qp) :: x,y,z,hx,hy,hz
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: vmat_qp(:,:)
  !
continue
  !
  nqp = size(xi)
  nqo = nqp-1
  nqq = nqp*nqp*nqp
  !
  if (present(number_of_modes)) then
    nso = number_of_modes
  else
    nso = nqo
  end if
  nsp = nso+1
  nsq = nsp*nsp*nsp
  !
  ! Allocate the return array containing the
  ! working precision Vandermonde matrix
  !
  allocate ( return_value(1:nqq,1:nsq) , source=zero , stat=ierr )
  !
  ! Create the Vandermonde matrix using quad precision
  !
  allocate ( vmat_qp(1:nqq,1:nsq) , source=qzero , stat=ierr )
  !
  do nk = 1,nqp
    do nj = 1,nqp
      do ni = 1,nqp
        !
        ! n is the 1D index of the solution point at
        ! which each basis function is being evaluated
        !
        n = nqp*nqp*(nk-1) + nqp*(nj-1) + ni
        !
        x = real( xi(ni) , kind=qp )
        y = real( xi(nj) , kind=qp )
        z = real( xi(nk) , kind=qp )
        !
        do mk = 0,nso
          do mj = 0,nso
            do mi = 0,nso
              !
              ! m is the index of the current basis function being evaluated
              !
              m = nsp*nsp*mk + nsp*mj + mi + 1
              !
              hx = normalized_jacobi_poly(mi,qzero,qzero,x)
              hy = normalized_jacobi_poly(mj,qzero,qzero,y)
              hz = normalized_jacobi_poly(mk,qzero,qzero,z)
              !
              vmat_qp(n,m) = hx*hy*hz
              !
            end do
          end do
        end do
        !
      end do
    end do
  end do
  !
  ! Use the chop function to copy vmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( vmat_qp )
  !
end function Vandermonde_Hexa
!
!###############################################################################
!
subroutine vandermonde_memory_usage(iunit)
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
  ! Memory usage of the vand array of type vandermonde_matrix
  !
  if (allocated(vand)) then
    !
    call dt%set( storage_size(vand,kind=inttype) , &
                         size(vand,kind=inttype) )
    !
    do j = lbound(vand,dim=2),ubound(vand,dim=2)
      do i = lbound(vand,dim=1),ubound(vand,dim=1)
        !
        if (allocated(vand(i,j)%nodal2modal%mat)) then
          !
          call array%set( storage_size(vand(i,j)%nodal2modal,kind=inttype) )
          !
          write (array_name,4) "vand",Geom_Name(i),j,"nodal2modal"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          array = vand(i,j)%nodal2modal%memory_usage()
          !
          write (array_name,4) "vand",Geom_Name(i),j,"nodal2modal%mat"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
        end if
        !
        if (allocated(vand(i,j)%modal2nodal%mat)) then
          !
          call array%set( storage_size(vand(i,j)%modal2nodal,kind=inttype) )
          !
          write (array_name,4) "vand",Geom_Name(i),j,"modal2nodal"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          array = vand(i,j)%modal2nodal%memory_usage()
          !
          write (array_name,4) "vand",Geom_Name(i),j,"modal2nodal%mat"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
        end if
        !
      end do
    end do
    !
    call total%add(dt)
    write (array_name,5) "vand"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "vand"
    !
  end if
  !
  ! Memory usage of the monomials array of type vandermonde_matrix
  !
  if (allocated(monomials)) then
    !
    call dt%set( storage_size(monomials,kind=inttype) , &
                         size(monomials,kind=inttype) )
    !
    do j = lbound(monomials,dim=2),ubound(monomials,dim=2)
      do i = lbound(monomials,dim=1),ubound(monomials,dim=1)
        !
        if (allocated(monomials(i,j)%nodal2modal%mat)) then
          !
          call array%set(storage_size(monomials(i,j)%nodal2modal,kind=inttype))
          !
          write (array_name,4) "monomials",Geom_Name(i),j,"nodal2modal"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          array = monomials(i,j)%nodal2modal%memory_usage()
          !
          write (array_name,4) "monomials",Geom_Name(i),j,"nodal2modal%mat"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
        end if
        !
        if (allocated(monomials(i,j)%modal2nodal%mat)) then
          !
          call array%set(storage_size(monomials(i,j)%modal2nodal,kind=inttype))
          !
          write (array_name,4) "monomials",Geom_Name(i),j,"modal2nodal"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          array = monomials(i,j)%modal2nodal%memory_usage()
          !
          write (array_name,4) "monomials",Geom_Name(i),j,"modal2nodal%mat"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
        end if
        !
      end do
    end do
    !
    call total%add(dt)
    write (array_name,5) "monomials"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "monomials"
    !
  end if
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module vandermonde_mod")
  2 format (" Array: ",a,dt)
  3 format (" Total memory = ",f8.2,1x,a)
  4 format (a,"(",a,",",i0,")%",a)
  5 format ("'",a,"'")
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  !
#endif
end subroutine vandermonde_memory_usage
!
!###############################################################################
!
include "Functions/np2n.f90"
include "Functions/cell_solpts.f90"
!
!###############################################################################
!
end module vandermonde_mod
!
!###############################################################################
!
