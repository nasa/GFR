module vandermonde_mod
  !
!!!#define PROFILE_ON
  !
  !.. Use Statements ..
  use module_kind_types, global_debug_timer =>      debug_timer, &
                    pure_global_debug_timer => pure_debug_timer
  !
  use matrix_types_mod, only : matrix
  use matrix_types_mod, only : create_matrix_sparsity
  use matrix_types_mod, only : output_matrix
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
  public :: move_vand_to_device
  public :: move_monomials_to_device
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
  use order_mod,        only : geom_solpts
  use order_mod,        only : n_min_geom,n_max_geom
  use order_mod,        only : n_min_order,n_max_order
  use quadrature_mod,   only : geom_is_used
  use element_type_mod, only : std_elem
  !
  !.. Local Scalars ..
  integer :: gmin,gmax,omin,omax,ierr
  integer :: this_geom,this_order,npts
  character(len=200) :: array_name
  character(len=300) :: timer_message
  !
  !.. Local Pointers ..
  type(vandermonde_matrix), pointer :: this_vand
! real(wp), pointer :: solpts(:,:)
  real(wp), pointer :: vmat(:,:)
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
      write (timer_message,100) trim(adjustl(Geom_Name(this_geom))),this_order
      100 format ("modal2nodal Geometry = ",a,"  ;  Order = ",i0)
      call debug_timer(start_timer,timer_message)
      vmat(:,:) = Compute_Basis_Matrix(this_geom,this_order)
      call debug_timer(stop_timer,timer_message)
      !
      ! Compute the inverse of the Vandermonde matrix that was just computed
      !
      write (timer_message,101) trim(adjustl(Geom_Name(this_geom))),this_order
      101 format ("nodal2modal Geometry = ",a,"  ;  Order = ",i0)
      call debug_timer(start_timer,timer_message)
      vmat_inv(:,:) = Vandermonde_Inverse( vmat )
      call debug_timer(stop_timer,timer_message)
      !
      error_message = "Modal-to-Nodal Vandermonde"
      call output_matrix(this_vand%modal2nodal, &
                         error_message,this_geom,this_order)
      !
      error_message = "Nodal-to-Modal Vandermonde"
      call output_matrix(this_vand%nodal2modal, &
                         error_message,this_geom,this_order)
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
  do this_geom = gmin,gmax
    do this_order = omin,omax
      if (allocated(vand(this_geom,this_order)%modal2nodal%mat)) then
        write (timer_message,200) trim(adjustl(Geom_Name(this_geom))),this_order
        200 format ("Sparsity check modal2nodal:  Geometry = ",a, &
                    "  ;  Order = ",i0)
        call debug_timer(start_timer,timer_message)
        call create_matrix_sparsity( vand(this_geom,this_order)%modal2nodal )
        call debug_timer(stop_timer,timer_message)
      end if
      if (allocated(vand(this_geom,this_order)%nodal2modal%mat)) then
        write (timer_message,201) trim(adjustl(Geom_Name(this_geom))),this_order
        201 format ("Sparsity check nodal2modal:  Geometry = ",a, &
                    "  ;  Order = ",i0)
        call debug_timer(start_timer,timer_message)
        call create_matrix_sparsity( vand(this_geom,this_order)%nodal2modal )
        call debug_timer(stop_timer,timer_message)
      end if
    end do
  end do
 !call create_matrix_sparsity( vand%modal2nodal )
 !call create_matrix_sparsity( vand%nodal2modal )
 !do this_geom = gmin,gmax
 !  do this_order = omin,omax
 !    call create_matrix_sparsity( vand(this_geom,this_order)%modal2nodal )
 !    call create_matrix_sparsity( vand(this_geom,this_order)%nodal2modal )
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
  use order_mod,        only : geom_solpts
  use element_type_mod, only : std_elem
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
#ifndef DISABLE_PURE
pure &
#endif
function Vandermonde_Inverse(vand_mat) result(return_value)
  !
  use lapack95, only : getrf,getri
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
  real(qp), allocatable :: vmat_inv(:,:)
  integer,  allocatable :: ipiv(:)
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
  ! Convert the Vandermonde matrix to quad precision and temporarily
  ! store in the array for the inverse Vandermonde matrix
  !
  vmat_inv = vand_mat ! USING F2003 AUTO-REALLOCATION
  !
  ! Compute the inverse of the Vandermonde matrix
  !
#ifdef SPECIAL_FOR_INTEL
  vmat_inv(:,:) = invert_matrix_lapack95( vmat_inv )
#else
 !vmat_inv(:,:) = invert_matrix( vmat_inv_qp )
  vmat_inv(:,:) = invert_matrix90( vmat_inv )
#endif
  !
  ! Use the chop function to copy vmat_inv_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value(:,:) = chop( vmat_inv )
  !
! allocate ( vmat_inv_qp(1:n1,1:n2) , source=qzero , stat=ierr )
! !
! ! Convert the Vandermonde matrix to quad precision and temporarily
! ! store in the array for the inverse Vandermonde matrix
! !
! vmat_inv_qp(:,:) = real( vand_mat , kind=qp )
! !
! ! Compute the inverse of the Vandermonde matrix
! !
!!vmat_inv_qp(:,:) = invert_matrix( vmat_inv_qp )
! vmat_inv_qp(:,:) = invert_matrix90( vmat_inv_qp )
! !
! ! Use the chop function to copy vmat_inv_qp into return_value while also
! ! setting any values that are less than working machine epsilon to zero.
! !
! return_value(:,:) = chop( vmat_inv_qp )
  !
end function Vandermonde_Inverse
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function Vandermonde_Edge(xi,number_of_modes) result(return_value)
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
#ifndef DISABLE_PURE
pure &
#endif
function Vandermonde_Quad(xi,number_of_modes) result(return_value)
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
#ifndef DISABLE_PURE
pure &
#endif
function Vandermonde_Tria(rs,number_of_modes) result(return_value)
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
#ifndef DISABLE_PURE
pure &
#endif
function Vandermonde_Tetr(rs,number_of_modes) result(return_value)
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
  integer  :: i,j,k,l,n,np,nb,m,ierr
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
  if (present(number_of_modes)) then
    n = number_of_modes
  else
    n = solpts2order(Geom_Tetr,np)
  end if
  nb = cell_solpts(Geom_Tetr,n)
  !
  ! If n is still -1, return with the Vandermonde matrix still being zero
  !
  if (n == -1) then
    return_value = init_array(np,np,zero)
    return
  end if
  !
  ! Create the Vandermonde matrix using quad precision
  !
  vmat_qp = init_array(np,nb,qzero) ! F2003 AUTO-REALLOC
  !
  m = 0
  do i = 0,n
    do j = 0,n-i
      do k = 0,n-i-j
        !
        m = m + 1
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
          alpha = real( 2*i + 1 , kind=qp )
          beta  = qzero
          h2 = normalized_jacobi_poly(j,alpha,beta,b)
          !
          alpha = real( 2*i + 2*j + 2 , kind=qp )
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
      end do
    end do
  end do
  !
  ! Use the chop function to copy vmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( vmat_qp ) ! F2003 AUTO-REALLOC
  !
end function Vandermonde_Tetr
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function Vandermonde_Pyra(rs,number_of_modes) result(return_value)
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
  integer  :: i,j,k,l,n,np,nb,m,ierr
  real(qp) :: r,s,t,a,b,c,h1,h2,h3,alpha,beta
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: vmat_qp(:,:)
  !
continue
  !
  np = size(rs,2)
  !
  ! Find the order for this pyramid
  !
  if (present(number_of_modes)) then
    n = number_of_modes
  else
    n = solpts2order(Geom_Pyra,np)
  end if
  nb = cell_solpts(Geom_Pyra,n)
  !
  ! If n is still -1, return with the Vandermonde matrix still being zero
  !
  if (n == -1) then
    return_value = init_array(np,np,zero)
    return
  end if
  !
  ! Create the Vandermonde matrix using quad precision
  !
  vmat_qp = init_array(np,nb,qzero) ! F2003 AUTO-REALLOC
  !
  m = 0
  do k = 0,n
    do j = 0,n-k
      do i = 0,n-k
        !
        m = m + 1
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
          a = qtwo * r * reciprocal(qone-t)
          b = qtwo * s * reciprocal(qone-t)
          c = t
          !
          ! Get the value of the Jacobi polynomial for a and b
          !
          alpha = qzero
          beta  = qzero
          h1 = normalized_jacobi_poly(i,alpha,beta,a)
          !
          alpha = qzero
          beta  = qzero
          h2 = normalized_jacobi_poly(j,alpha,beta,b)
          !
          alpha = real( 2*i + 2*j + 2 , kind=qp )
          beta  = qzero
          h3 = normalized_jacobi_poly(k,alpha,beta,c)
          !
          ! Evaluate the basis function psi(m) at node l
          !
          vmat_qp(l,m) = qtwo * h1 * h2 * h3 * (qone-c)**(i+j)
          !
        end do
        !
      end do
    end do
  end do
  !
  ! Use the chop function to copy vmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( vmat_qp ) ! F2003 AUTO-REALLOC
  !
end function Vandermonde_Pyra
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function Vandermonde_Pris(rs,number_of_modes) result(return_value)
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
  integer  :: i,j,k,l,n,np,nb,m,ierr
  real(qp) :: r,s,t,a,b,c,h1,h2,h3,alpha,beta
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: vmat_qp(:,:)
  !
continue
  !
  np = size(rs,2)
  !
  ! Find the order for this prism
  !
  if (present(number_of_modes)) then
    n = number_of_modes
  else
    n = solpts2order(Geom_Pris,np)
  end if
  nb = cell_solpts(Geom_Pris,n)
  !
  ! If n is still -1, return with the Vandermonde matrix still being zero
  !
  if (n == -1) then
    return_value = init_array(np,np,zero)
    return
  end if
  !
  ! Create the Vandermonde matrix using quad precision
  !
  vmat_qp = init_array(np,nb,qzero) ! F2003 AUTO-REALLOC
  !
  m = 0
  do k = 0,n
    do i = 0,n
      do j = 0,n-i
        !
        m = m + 1
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
          a = qtwo * (qone+r) * reciprocal(qone-s) - qone
          b = s
          c = t
          !
          ! Get the value of the Jacobi polynomial for a and b
          !
          alpha = qzero
          beta  = qzero
          h1 = normalized_jacobi_poly(i,alpha,beta,a)
          !
          alpha = real( 2*i + 1 , kind=qp )
          beta  = qzero
          h2 = normalized_jacobi_poly(j,alpha,beta,b)
          !
          alpha = qzero
          beta  = qzero
          h3 = normalized_jacobi_poly(k,alpha,beta,c)
          !
          ! Evaluate the basis function psi(m) at node l
          !
          vmat_qp(l,m) = sqrt(qtwo) * h1 * h2 * h3 * (qone-b)**i
          !
        end do
        !
      end do
    end do
  end do
  !
  ! Use the chop function to copy vmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( vmat_qp ) ! F2003 AUTO-REALLOC
  !
end function Vandermonde_Pris
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function Vandermonde_Hexa(xi,number_of_modes) result(return_value)
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
  real(qp), allocatable :: vnjp_qp(:,:)
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
  ! Create the Vandermonde matrix using quad precision
  !
  allocate ( vmat_qp(1:nqq,1:nsq) , source=qzero , stat=ierr )
  !
  allocate ( vnjp_qp(0:nso,1:nqp) , source=qzero , stat=ierr )
  !
  do n = 1,nqp
    do m = 0,nso
      x = real( xi(n) , kind=qp )
      vnjp_qp(m,n) = normalized_jacobi_poly(m,qzero,qzero,x)
    end do
  end do
  !
  n = 0
  do nk = 1,nqp
    do nj = 1,nqp
      do ni = 1,nqp
        !
        ! n is the 1D index of the solution point at
        ! which each basis function is being evaluated
        !
       !n = nqp*nqp*(nk-1) + nqp*(nj-1) + ni
        n = n + 1
        !
        m = 0
        do mk = 0,nso
          do mj = 0,nso
            do mi = 0,nso
              !
              ! m is the index of the current basis function being evaluated
              !
             !m = nsp*nsp*mk + nsp*mj + mi + 1
             !!
             !hx = vnjp_qp(mi,ni)
             !hy = vnjp_qp(mj,nj)
             !hz = vnjp_qp(mk,nk)
             !vmat_qp(n,m) = hx*hy*hz
              !
              m = m + 1
              vmat_qp(n,m) = vnjp_qp(mi,ni)*vnjp_qp(mj,nj)*vnjp_qp(mk,nk)
              !
            end do
          end do
        end do
        !
      end do
    end do
  end do
  !
! do nk = 1,nqp
!   do nj = 1,nqp
!     do ni = 1,nqp
!       !
!       ! n is the 1D index of the solution point at
!       ! which each basis function is being evaluated
!       !
!       n = nqp*nqp*(nk-1) + nqp*(nj-1) + ni
!       !
!       x = real( xi(ni) , kind=qp )
!       y = real( xi(nj) , kind=qp )
!       z = real( xi(nk) , kind=qp )
!       !
!       do mk = 0,nso
!         do mj = 0,nso
!           do mi = 0,nso
!             !
!             ! m is the index of the current basis function being evaluated
!             !
!             m = nsp*nsp*mk + nsp*mj + mi + 1
!             !
!             hx = normalized_jacobi_poly(mi,qzero,qzero,x)
!             hy = normalized_jacobi_poly(mj,qzero,qzero,y)
!             hz = normalized_jacobi_poly(mk,qzero,qzero,z)
!             !
!             vmat_qp(n,m) = hx*hy*hz
!             !
!           end do
!         end do
!       end do
!       !
!     end do
!   end do
! end do
  !
  ! Use the chop function to copy vmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( vmat_qp ) ! F2003 AUTO-REALLOC
  !
end function Vandermonde_Hexa
!
!###############################################################################
!
subroutine vandermonde_memory_usage(iunit)
  !
  !.. Use Statements ..
  use module_memory
#ifndef DISABLE_DTIO
  use matrix_types_mod, only : matrix_memory_usage
#endif
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
          array = matrix_memory_usage( vand(i,j)%nodal2modal )
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
          array = matrix_memory_usage( vand(i,j)%modal2nodal )
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
          array = matrix_memory_usage( monomials(i,j)%nodal2modal )
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
          array = matrix_memory_usage( monomials(i,j)%modal2nodal )
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
subroutine move_vand_to_device()
  !
  !.. Local Scalars ..
  integer :: i,i1,i2,j,j1,j2,m1,m2,n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "move_vand_to_device"
  !
continue
  !
  if (allocated(vand)) then
    !$ACC ENTER DATA COPYIN(vand)
    i1 = lbound(vand,dim=1)
    i2 = ubound(vand,dim=1)
    j1 = lbound(vand,dim=2)
    j2 = ubound(vand,dim=2)
    !$ACC ENTER DATA COPYIN(vand(i1:i2,j1:j2))
    do j = j1,j2
      do i = i1,i2
        !
        !---------------------------------------------------------------
        !
        if (allocated(vand(i,j)%nodal2modal%mat)) then
          !$ACC ENTER DATA COPYIN(vand(i,j)%nodal2modal%mat)
          m1 = lbound(vand(i,j)%nodal2modal%mat,dim=1)
          m2 = ubound(vand(i,j)%nodal2modal%mat,dim=1)
          n1 = lbound(vand(i,j)%nodal2modal%mat,dim=2)
          n2 = ubound(vand(i,j)%nodal2modal%mat,dim=2)
          !$ACC ENTER DATA COPYIN(vand(i,j)%nodal2modal%mat(m1:m2,n1:n2))
        end if
        !
        if (allocated(vand(i,j)%nodal2modal%row_ptr)) then
          !$ACC ENTER DATA COPYIN(vand(i,j)%nodal2modal%row_ptr)
          m1 = lbound(vand(i,j)%nodal2modal%row_ptr,dim=1)
          m2 = ubound(vand(i,j)%nodal2modal%row_ptr,dim=1)
          !$ACC ENTER DATA COPYIN(vand(i,j)%nodal2modal%row_ptr(m1:m2))
        end if
        !
        if (allocated(vand(i,j)%nodal2modal%row_idx)) then
          !$ACC ENTER DATA COPYIN(vand(i,j)%nodal2modal%row_idx)
          m1 = lbound(vand(i,j)%nodal2modal%row_idx,dim=1)
          m2 = ubound(vand(i,j)%nodal2modal%row_idx,dim=1)
          !$ACC ENTER DATA COPYIN(vand(i,j)%nodal2modal%row_idx(m1:m2))
        end if
        !
        if (allocated(vand(i,j)%nodal2modal%row_val)) then
          !$ACC ENTER DATA COPYIN(vand(i,j)%nodal2modal%row_val)
          m1 = lbound(vand(i,j)%nodal2modal%row_val,dim=1)
          m2 = ubound(vand(i,j)%nodal2modal%row_val,dim=1)
          !$ACC ENTER DATA COPYIN(vand(i,j)%nodal2modal%row_val(m1:m2))
        end if
        !
        if (allocated(vand(i,j)%nodal2modal%col_ptr)) then
          !$ACC ENTER DATA COPYIN(vand(i,j)%nodal2modal%col_ptr)
          m1 = lbound(vand(i,j)%nodal2modal%col_ptr,dim=1)
          m2 = ubound(vand(i,j)%nodal2modal%col_ptr,dim=1)
          !$ACC ENTER DATA COPYIN(vand(i,j)%nodal2modal%col_ptr(m1:m2))
        end if
        !
        if (allocated(vand(i,j)%nodal2modal%col_idx)) then
          !$ACC ENTER DATA COPYIN(vand(i,j)%nodal2modal%col_idx)
          m1 = lbound(vand(i,j)%nodal2modal%col_idx,dim=1)
          m2 = ubound(vand(i,j)%nodal2modal%col_idx,dim=1)
          !$ACC ENTER DATA COPYIN(vand(i,j)%nodal2modal%col_idx(m1:m2))
        end if
        !
        if (allocated(vand(i,j)%nodal2modal%col_val)) then
          !$ACC ENTER DATA COPYIN(vand(i,j)%nodal2modal%col_val)
          m1 = lbound(vand(i,j)%nodal2modal%col_val,dim=1)
          m2 = ubound(vand(i,j)%nodal2modal%col_val,dim=1)
          !$ACC ENTER DATA COPYIN(vand(i,j)%nodal2modal%col_val(m1:m2))
        end if
        !
        !---------------------------------------------------------------
        !
        if (allocated(vand(i,j)%modal2nodal%mat)) then
          !$ACC ENTER DATA COPYIN(vand(i,j)%modal2nodal%mat)
          m1 = lbound(vand(i,j)%modal2nodal%mat,dim=1)
          m2 = ubound(vand(i,j)%modal2nodal%mat,dim=1)
          n1 = lbound(vand(i,j)%modal2nodal%mat,dim=2)
          n2 = ubound(vand(i,j)%modal2nodal%mat,dim=2)
          !$ACC ENTER DATA COPYIN(vand(i,j)%modal2nodal%mat(m1:m2,n1:n2))
        end if
        !
        if (allocated(vand(i,j)%modal2nodal%row_ptr)) then
          !$ACC ENTER DATA COPYIN(vand(i,j)%modal2nodal%row_ptr)
          m1 = lbound(vand(i,j)%modal2nodal%row_ptr,dim=1)
          m2 = ubound(vand(i,j)%modal2nodal%row_ptr,dim=1)
          !$ACC ENTER DATA COPYIN(vand(i,j)%modal2nodal%row_ptr(m1:m2))
        end if
        !
        if (allocated(vand(i,j)%modal2nodal%row_idx)) then
          !$ACC ENTER DATA COPYIN(vand(i,j)%modal2nodal%row_idx)
          m1 = lbound(vand(i,j)%modal2nodal%row_idx,dim=1)
          m2 = ubound(vand(i,j)%modal2nodal%row_idx,dim=1)
          !$ACC ENTER DATA COPYIN(vand(i,j)%modal2nodal%row_idx(m1:m2))
        end if
        !
        if (allocated(vand(i,j)%modal2nodal%row_val)) then
          !$ACC ENTER DATA COPYIN(vand(i,j)%modal2nodal%row_val)
          m1 = lbound(vand(i,j)%modal2nodal%row_val,dim=1)
          m2 = ubound(vand(i,j)%modal2nodal%row_val,dim=1)
          !$ACC ENTER DATA COPYIN(vand(i,j)%modal2nodal%row_val(m1:m2))
        end if
        !
        if (allocated(vand(i,j)%modal2nodal%col_ptr)) then
          !$ACC ENTER DATA COPYIN(vand(i,j)%modal2nodal%col_ptr)
          m1 = lbound(vand(i,j)%modal2nodal%col_ptr,dim=1)
          m2 = ubound(vand(i,j)%modal2nodal%col_ptr,dim=1)
          !$ACC ENTER DATA COPYIN(vand(i,j)%modal2nodal%col_ptr(m1:m2))
        end if
        !
        if (allocated(vand(i,j)%modal2nodal%col_idx)) then
          !$ACC ENTER DATA COPYIN(vand(i,j)%modal2nodal%col_idx)
          m1 = lbound(vand(i,j)%modal2nodal%col_idx,dim=1)
          m2 = ubound(vand(i,j)%modal2nodal%col_idx,dim=1)
          !$ACC ENTER DATA COPYIN(vand(i,j)%modal2nodal%col_idx(m1:m2))
        end if
        !
        if (allocated(vand(i,j)%modal2nodal%col_val)) then
          !$ACC ENTER DATA COPYIN(vand(i,j)%modal2nodal%col_val)
          m1 = lbound(vand(i,j)%modal2nodal%col_val,dim=1)
          m2 = ubound(vand(i,j)%modal2nodal%col_val,dim=1)
          !$ACC ENTER DATA COPYIN(vand(i,j)%modal2nodal%col_val(m1:m2))
        end if
        !
        !---------------------------------------------------------------
        !
      end do
    end do
  end if
  !
end subroutine move_vand_to_device
!
!###############################################################################
!
subroutine move_monomials_to_device()
  !
  !.. Local Scalars ..
  integer :: i,i1,i2,j,j1,j2,m1,m2,n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "move_monomials_to_device"
  !
continue
  !
  if (allocated(monomials)) then
    !$ACC ENTER DATA COPYIN(monomials)
    i1 = lbound(monomials,dim=1)
    i2 = ubound(monomials,dim=1)
    j1 = lbound(monomials,dim=2)
    j2 = ubound(monomials,dim=2)
    !$ACC ENTER DATA COPYIN(monomials(i1:i2,j1:j2))
    do j = j1,j2
      do i = i1,i2
        !
        !---------------------------------------------------------------
        !
        !$ACC ENTER DATA COPYIN(monomials(i,j)%nodal2modal)
        if (allocated(monomials(i,j)%nodal2modal%mat)) then
          !$ACC ENTER DATA COPYIN(monomials(i,j)%nodal2modal%mat)
          m1 = lbound(monomials(i,j)%nodal2modal%mat,dim=1)
          m2 = ubound(monomials(i,j)%nodal2modal%mat,dim=1)
          n1 = lbound(monomials(i,j)%nodal2modal%mat,dim=2)
          n2 = ubound(monomials(i,j)%nodal2modal%mat,dim=2)
          !$ACC ENTER DATA COPYIN(monomials(i,j)%nodal2modal%mat(m1:m2,n1:n2))
        end if
        !
        if (allocated(monomials(i,j)%nodal2modal%row_ptr)) then
          !$ACC ENTER DATA COPYIN(monomials(i,j)%nodal2modal%row_ptr)
          m1 = lbound(monomials(i,j)%nodal2modal%row_ptr,dim=1)
          m2 = ubound(monomials(i,j)%nodal2modal%row_ptr,dim=1)
          !$ACC ENTER DATA COPYIN(monomials(i,j)%nodal2modal%row_ptr(m1:m2))
        end if
        !
        if (allocated(monomials(i,j)%nodal2modal%row_idx)) then
          !$ACC ENTER DATA COPYIN(monomials(i,j)%nodal2modal%row_idx)
          m1 = lbound(monomials(i,j)%nodal2modal%row_idx,dim=1)
          m2 = ubound(monomials(i,j)%nodal2modal%row_idx,dim=1)
          !$ACC ENTER DATA COPYIN(monomials(i,j)%nodal2modal%row_idx(m1:m2))
        end if
        !
        if (allocated(monomials(i,j)%nodal2modal%row_val)) then
          !$ACC ENTER DATA COPYIN(monomials(i,j)%nodal2modal%row_val)
          m1 = lbound(monomials(i,j)%nodal2modal%row_val,dim=1)
          m2 = ubound(monomials(i,j)%nodal2modal%row_val,dim=1)
          !$ACC ENTER DATA COPYIN(monomials(i,j)%nodal2modal%row_val(m1:m2))
        end if
        !
        if (allocated(monomials(i,j)%nodal2modal%col_ptr)) then
          !$ACC ENTER DATA COPYIN(monomials(i,j)%nodal2modal%col_ptr)
          m1 = lbound(monomials(i,j)%nodal2modal%col_ptr,dim=1)
          m2 = ubound(monomials(i,j)%nodal2modal%col_ptr,dim=1)
          !$ACC ENTER DATA COPYIN(monomials(i,j)%nodal2modal%col_ptr(m1:m2))
        end if
        !
        if (allocated(monomials(i,j)%nodal2modal%col_idx)) then
          !$ACC ENTER DATA COPYIN(monomials(i,j)%nodal2modal%col_idx)
          m1 = lbound(monomials(i,j)%nodal2modal%col_idx,dim=1)
          m2 = ubound(monomials(i,j)%nodal2modal%col_idx,dim=1)
          !$ACC ENTER DATA COPYIN(monomials(i,j)%nodal2modal%col_idx(m1:m2))
        end if
        !
        if (allocated(monomials(i,j)%nodal2modal%col_val)) then
          !$ACC ENTER DATA COPYIN(monomials(i,j)%nodal2modal%col_val)
          m1 = lbound(monomials(i,j)%nodal2modal%col_val,dim=1)
          m2 = ubound(monomials(i,j)%nodal2modal%col_val,dim=1)
          !$ACC ENTER DATA COPYIN(monomials(i,j)%nodal2modal%col_val(m1:m2))
        end if
        !
        !---------------------------------------------------------------
        !
        !$ACC ENTER DATA COPYIN(monomials(i,j)%modal2nodal)
        if (allocated(monomials(i,j)%modal2nodal%mat)) then
          !$ACC ENTER DATA COPYIN(monomials(i,j)%modal2nodal%mat)
          m1 = lbound(monomials(i,j)%modal2nodal%mat,dim=1)
          m2 = ubound(monomials(i,j)%modal2nodal%mat,dim=1)
          n1 = lbound(monomials(i,j)%modal2nodal%mat,dim=2)
          n2 = ubound(monomials(i,j)%modal2nodal%mat,dim=2)
          !$ACC ENTER DATA COPYIN(monomials(i,j)%modal2nodal%mat(m1:m2,n1:n2))
        end if
        !
        if (allocated(monomials(i,j)%modal2nodal%row_ptr)) then
          !$ACC ENTER DATA COPYIN(monomials(i,j)%modal2nodal%row_ptr)
          m1 = lbound(monomials(i,j)%modal2nodal%row_ptr,dim=1)
          m2 = ubound(monomials(i,j)%modal2nodal%row_ptr,dim=1)
          !$ACC ENTER DATA COPYIN(monomials(i,j)%modal2nodal%row_ptr(m1:m2))
        end if
        !
        if (allocated(monomials(i,j)%modal2nodal%row_idx)) then
          !$ACC ENTER DATA COPYIN(monomials(i,j)%modal2nodal%row_idx)
          m1 = lbound(monomials(i,j)%modal2nodal%row_idx,dim=1)
          m2 = ubound(monomials(i,j)%modal2nodal%row_idx,dim=1)
          !$ACC ENTER DATA COPYIN(monomials(i,j)%modal2nodal%row_idx(m1:m2))
        end if
        !
        if (allocated(monomials(i,j)%modal2nodal%row_val)) then
          !$ACC ENTER DATA COPYIN(monomials(i,j)%modal2nodal%row_val)
          m1 = lbound(monomials(i,j)%modal2nodal%row_val,dim=1)
          m2 = ubound(monomials(i,j)%modal2nodal%row_val,dim=1)
          !$ACC ENTER DATA COPYIN(monomials(i,j)%modal2nodal%row_val(m1:m2))
        end if
        !
        if (allocated(monomials(i,j)%modal2nodal%col_ptr)) then
          !$ACC ENTER DATA COPYIN(monomials(i,j)%modal2nodal%col_ptr)
          m1 = lbound(monomials(i,j)%modal2nodal%col_ptr,dim=1)
          m2 = ubound(monomials(i,j)%modal2nodal%col_ptr,dim=1)
          !$ACC ENTER DATA COPYIN(monomials(i,j)%modal2nodal%col_ptr(m1:m2))
        end if
        !
        if (allocated(monomials(i,j)%modal2nodal%col_idx)) then
          !$ACC ENTER DATA COPYIN(monomials(i,j)%modal2nodal%col_idx)
          m1 = lbound(monomials(i,j)%modal2nodal%col_idx,dim=1)
          m2 = ubound(monomials(i,j)%modal2nodal%col_idx,dim=1)
          !$ACC ENTER DATA COPYIN(monomials(i,j)%modal2nodal%col_idx(m1:m2))
        end if
        !
        if (allocated(monomials(i,j)%modal2nodal%col_val)) then
          !$ACC ENTER DATA COPYIN(monomials(i,j)%modal2nodal%col_val)
          m1 = lbound(monomials(i,j)%modal2nodal%col_val,dim=1)
          m2 = ubound(monomials(i,j)%modal2nodal%col_val,dim=1)
          !$ACC ENTER DATA COPYIN(monomials(i,j)%modal2nodal%col_val(m1:m2))
        end if
        !
        !---------------------------------------------------------------
        !
      end do
    end do
  end if
  !
end subroutine move_monomials_to_device
!
!###############################################################################
!
include "Functions/solpts2order.f90"
include "Functions/np2n.f90"
include "Functions/cell_solpts.f90"
#include "Functions/debug_timer.f90"
!
!###############################################################################
!
end module vandermonde_mod
!
!###############################################################################
!
