module correction_mod
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
  public :: init_correction_matrices
  public :: correction_memory_usage
  !
  !
  !
  ! correct : correction matrix for all combinations of cell geometry and
  !           and cell order that are possible in the current simulation
  !
  !     For cell 'n' :
  !
  !         The correction matrix (giving the contribution to each solution
  !         point within cell 'n' from the interface flux point 'k' that is
  !         currently of interest) is given by the expression
  !
  !         correct(this_geom,this_order)%mat(:,k)
  !
  !     Where :
  !
  !        this_geom = cell(n)%geom
  !        this_order = cell(n)%order
  !        k = the flx-pt on a face of cell 'n' that is currently of interest
  !
  type(matrix), public, save, target, allocatable :: correct(:,:)
  !
contains
!
!###############################################################################
!
subroutine init_correction_matrices()
  !
  !.. Use Statements ..
  use order_mod,      only : geom_solpts,geom_flxpts
  use order_mod,      only : n_min_geom,n_max_geom
  use order_mod,      only : n_min_order,n_max_order
  use order_mod,      only : n_order
  use quadrature_mod, only : geom_is_used
  use quadrature_mod, only : std_elem
  !
  !.. Local Scalars ..
  integer :: gmin,gmax,omin,omax
  integer :: this_geom,this_order
  integer :: nsolpts,nflxpts,ierr
  integer :: n,n1,n2,i
  character(len=200) :: array_name
  !
  !.. Local Pointers ..
  type(matrix), pointer :: this_correct
 !real(wp), pointer, contiguous :: solpts(:,:)
  real(wp), pointer :: solpts(:,:)
 !real(wp), pointer, contiguous :: cmat(:,:)
  real(wp), pointer :: cmat(:,:)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: edgpts(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_correction_matrices"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Initialize the local pointers to disassociated
  !
  this_correct => null()
  solpts => null()
  cmat => null()
  !
  gmin = n_min_geom
  gmax = n_max_geom
  omin = n_min_order
  omax = n_max_order
  !
  ! Allocate the correct array
  !
  allocate ( correct(gmin:gmax,omin:omax) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"correct",1,__LINE__,__FILE__,ierr,error_message)
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
      this_correct => correct(this_geom,this_order)
      !
      ! Get the number of solution points for this cell
      !
      nsolpts = geom_solpts(this_geom,this_order)
      !
      ! Get the number of face/flux points for this cell
      !
      nflxpts = geom_flxpts(this_geom,this_order)
      !
      ! Allocate the correction matrix
      !
      allocate ( this_correct%mat(1:nsolpts,1:nflxpts) , &
                 source=zero , stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Assign local pointers to these matrices as aliases to simplify the code
      !
      cmat => this_correct%mat
      !
      ! Compute the correction matrices depending on geometry type
      !
      select case (this_geom)
        case (Geom_Edge)
          !
          edgpts = std_elem(this_geom,this_order)%pts(1,:) ! F2003 AUTO-REALLOC
         !solpts => std_elem(this_geom,this_order)%pts
          !
          cmat(:,1) = CorrectionMatrix_Edge( edgpts )
         !cmat(:,1) = CorrectionMatrix_Edge( solpts(1,:) )
          cmat(:,2) = cmat(nsolpts:1:-1,1)
          !
        case (Geom_Tria)
          !
          solpts => std_elem(this_geom,this_order)%pts
          !
          cmat = CorrectionMatrix_Tria( solpts )
          !
        case (Geom_Quad)
          !
          edgpts = std_elem(Geom_Edge,this_order)%pts(1,:) ! F2003 AUTO-REALLOC
         !solpts => std_elem(Geom_Edge,this_order)%pts
          !
          cmat = CorrectionMatrix_Quad( edgpts )
         !cmat = CorrectionMatrix_Quad( solpts(1,:) )
          !
        case (Geom_Tetr)
          !
          solpts => std_elem(this_geom,this_order)%pts
          !
          cmat = CorrectionMatrix_Tetr( solpts )
          !
        case (Geom_Pyra)
          !
          solpts => std_elem(this_geom,this_order)%pts
          !
          cmat = CorrectionMatrix_Pyra( solpts )
          !
        case (Geom_Pris)
          !
          solpts => std_elem(this_geom,this_order)%pts
          !
          cmat = CorrectionMatrix_Pris( solpts )
          !
        case (Geom_Hexa)
          !
          edgpts = std_elem(Geom_Edge,this_order)%pts(1,:) ! F2003 AUTO-REALLOC
         !solpts => std_elem(Geom_Edge,this_order)%pts
          !
          cmat = CorrectionMatrix_Hexa( edgpts )
         !cmat = CorrectionMatrix_Hexa( solpts(1,:) )
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
      if (associated(this_correct)) this_correct => null()
      if (associated(solpts   ))    solpts       => null()
      if (associated(cmat     ))    cmat         => null()
      !
    end do order_loop
    !
  end do geom_loop
  !
  ! Check for sparsity within each of the correction matrices
  !
  call correct%check_for_sparsity
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("correct(",a,",",i0,")%mat")
  2 format (" A grid cell of an unknown geometry type was found!")
  !
end subroutine init_correction_matrices
!
!###############################################################################
!
pure function CorrectionMatrix_Edge(xi,eval_deriv) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : correction_function
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xi
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: eval_deriv
  !
  !.. Function Result ..
  real(wp), dimension(1:size(xi)) :: return_value
  !
  !.. Local Scalars ..
  integer  :: n,k,j,npts,ndeg
  logical(lk) :: compute_derivative
  real(qp) :: cn,cnm1,cnm2
  !
  !.. Local Arrays ..
  real(qp), dimension(1:size(xi)) :: x
  real(qp), dimension(1:size(xi)) :: cmat
  real(qp), dimension(0:size(xi)) :: coef_c
  real(qp), dimension(0:size(xi)) :: coef_dc
  !
  real(qp), dimension( 0:size(xi), 1:size(xi)) :: xexp
  real(qp), dimension(-1:size(xi),-1:size(xi)) :: coef_Le
  real(qp), dimension( 0:size(xi), 0:size(xi)) :: coef_RR
  !
continue
  !
  ! #########################################################################
  !
  !       THE LOOP INDICES IN THIS FUNCTION ARE USED FOR THE FOLLOWING
  !       SPECIFIC TASKS
  !
  !    n : index used for looping through polynomial degrees / exponents
  !    k : index used for looping over the solution points
  !    j : index used for looping through the coefficients of a polynomial
  !
  ! #########################################################################
  !
  npts = size(xi) ! number of solution points
  ndeg = size(xi) ! polynomial degree of the correction function
  !
  compute_derivative = true
  if (present(eval_deriv)) then
    compute_derivative = eval_deriv
  end if
  !
  ! Initialize the function result
  !
  return_value(:) = zero
  !
  ! Initialize the local arrays
  !
  x(:) = real( xi(:) , kind=qp )
  !
  cmat(:)      = qzero
  coef_c(:)    = qzero
  coef_dc(:)   = qzero
  coef_Le(:,:) = qzero
  coef_RR(:,:) = qzero
  !
  ! Initialize specific coefficients of Legendre polynomial
  !
  coef_Le(0,0) = qone
  coef_Le(1,1) = qone
  !
  ! Get the value of each nodal coordinate raised to all possible powers
  !   using xk = x(k), xexp(0:neg,k) contains the values of
  !        [ xk^0, xk^1, xk^2, ..., xk^(ndeg-1), xk^ndeg]
  !
  xexp(0,:) = qone
  do k = 1,npts
    do n = 1,ndeg
      xexp(n,k) = xexp(n-1,k) * x(k)
    end do
  end do
  !
  ! Get the coefficients for every Legendre polynomial up to degree ndeg
  !
  do n = 2,ndeg
    !
    cnm1 = real( 2*n-1 , kind=qp ) / real( n , kind=qp )
    cnm2 = real(   n-1 , kind=qp ) / real( n , kind=qp )
    !
    ! Loop to get the coefficients for Legendre polynomial of degree n
    !
    do j = 0,n
      coef_Le(j,n) = cnm1*coef_Le(j-1,n-1) - cnm2*coef_Le(j,n-2)
    end do
    !
  end do
  !
  ! Get the coefficients for every right Radau polynomial up to degree ndeg
  !
  do n = 1,ndeg
    !
    cn = qhalf * real( (-1)**n , kind=qp )
    !
    ! Loop to get the coefficients for right Radau polynomial of degree n
    !
    do j = 0,n
      coef_RR(j,n) = cn * (coef_Le(j,n) - coef_Le(j,n-1))
    end do
    !
  end do
  !
  ! Get the coefficients for the correction function
  !
  if (correction_function == Correction_g2) then
    !
    ! Use the g2 correction function
    !
    cn   = real( ndeg-1 , kind=qp ) / real( 2*ndeg-1 , kind=qp )
    cnm1 = real( ndeg   , kind=qp ) / real( 2*ndeg-1 , kind=qp )
    !
    do j = 0,ndeg
      coef_c(j) = cn*coef_RR(j,ndeg) + cnm1*coef_RR(j,ndeg-1)
    end do
    !
  else if (correction_function == Correction_gGauss) then
    !
    ! Use the gGauss correction function
    !
    cn   = real( ndeg-1 , kind=qp ) / real( 2*ndeg-1 , kind=qp )
    cnm1 = real( ndeg   , kind=qp ) / real( 2*ndeg-1 , kind=qp )
    !
    do j = 0,ndeg
      coef_c(j) = cn*coef_RR(j,ndeg) + cnm1*coef_RR(j,ndeg-1)
    end do
    !
  else
    !
    ! Use the right Radau polynomial for the correction function
    !
    do j = 0,ndeg
      coef_c(j) = coef_RR(j,ndeg)
    end do
    !
  end if
  !
  if (compute_derivative) then
    !
    ! Compute the coefficients for the derivative of the correction function
    !
    do j = 1,ndeg
      coef_dc(j-1) = real( j , kind=qp ) * coef_c(j)
    end do
    !
    ! Finally, get the value for the derivative of the
    ! correction function at each 1D solution point
    !
    do k = 1,npts
      cmat(k) = -sum( xexp(0:ndeg-1,k) * coef_dc(0:ndeg-1) )
    end do
    !
  else
    !
    ! Finally get the value of the correction function
    ! at each 1D solution point
    !
    do k = 1,npts
      cmat(k) = sum( xexp(0:ndeg,k) * coef_c(0:ndeg) )
    end do
    !
  end if
  !
  ! Use the chop function to copy cmat into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( cmat )
  !
end function CorrectionMatrix_Edge
!
!###############################################################################
!
pure function CorrectionMatrix_Quad_DG(pts1d,wts1d,eval_deriv) &
                                result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly1D
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: pts1d
  real(wp), dimension(:), intent(in) :: wts1d
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: eval_deriv
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,k,nf,np,nsp,nfp,ierr
  logical(lk) :: compute_derivative
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: cmat_qp(:,:)
  real(qp), allocatable :: pts1d_qp(:)
  real(qp), allocatable :: wts1d_qp(:)
  !
continue
  !
  np = size(pts1d)
  nsp = np*np
  nfp = 4*np
  !
  compute_derivative = true
  if (present(eval_deriv)) then
    compute_derivative = eval_deriv
  end if
  !
  ! Allocate the return array containing the working
  ! precision correction matrix
  !
  allocate ( return_value(1:nsp,1:nfp) , source=zero , stat=ierr )
  !
  ! Create the correction matrix using quad precision and make quad precision
  ! copies of the 1D solution points and their corresponding quadrature weights
  !
  allocate ( cmat_qp(1:nsp,1:nfp) , source=qzero , stat=ierr )
  allocate ( pts1d_qp(1:np) , source=qzero , stat=ierr )
  allocate ( wts1d_qp(1:np) , source=qzero , stat=ierr )
  !
  ! Face 1 - bottom/south edge of quad
  !
  nf = 1
  do i = 1,np
    do j = 1,np
      k = (j-1)*np + i
      cmat_qp(k,nf) = eval_LagrangePoly1D(j,-qone,pts1d_qp) / wts1d_qp(j)
    end do
    nf = nf + 1
  end do
  !
  ! Face 2 - right/east edge of quad
  !
  nf = np+1
  do j = 1,np
    do i = np,1,-1
      k = (j-1)*np + i
      cmat_qp(k,nf) = eval_LagrangePoly1D(i,qone,pts1d_qp) / wts1d_qp(i)
    end do
    nf = nf + 1
  end do
  !
  ! Face 3 - top/north edge of quad
  !
  nf = 2*np+1
  do i = np,1,-1
    do j = np,1,-1
      k = (j-1)*np + i
      cmat_qp(k,nf) = eval_LagrangePoly1D(j,qone,pts1d_qp) / wts1d_qp(j)
    end do
    nf = nf + 1
  end do
  !
  ! Face 4 - left/west edge of quad
  !
  nf = 3*np+1
  do j = np,1,-1
    do i = 1,np
      k = (j-1)*np + i
      cmat_qp(k,nf) = eval_LagrangePoly1D(i,-qone,pts1d_qp) / wts1d_qp(i)
    end do
    nf = nf + 1
  end do
  !
  ! Use the chop function to copy cmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( cmat_qp )
  !
end function CorrectionMatrix_Quad_DG
!
!###############################################################################
!
pure function CorrectionMatrix_Quad(xi,eval_deriv) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xi
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: eval_deriv
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,k,np,nf,nsp,nfp,ierr
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: cmat1D(:)
  !
continue
  !
  ! ###########################################################################
  ! ###########################################################################
  ! ####  NOTE: No need to use chopping function in this procedure. The    ####
  ! ####        quad correction matrix is really just the 1D correction    ####
  ! ####        matrix applied to each face/flux point of the standard     ####
  ! ####        quad element. The 1D correction matrix is obtained by the  ####
  ! ####        function CorrectionMatrix_Edge which applies the chopping  ####
  ! ####        function to its result so there is no need to repeat this  ####
  ! ####        operation again.                                           ####
  ! ###########################################################################
  ! ###########################################################################
  !
  np = size(xi)
  nsp = np*np
  nfp = 4*np
  !
  ! Allocate the return array containing the working
  ! precision correction matrix
  !
  allocate ( return_value(1:nsp,1:nfp) , source=zero , stat=ierr )
  !
  ! Create the 1D correction matrix
  !
  allocate ( cmat1D(1:np) , source=zero , stat=ierr )
  !
  ! 1D Correction matrix
  !
  cmat1D = CorrectionMatrix_Edge(xi,eval_deriv)
  !
  ! Face 1 - y/eta min face of quad (bottom/south)
  !
  nf = 1
  do j = 1,np
    do i = 1,np
      k = j + (i-1)*np
      return_value(k,nf) = cmat1D(i)
    end do
    nf = nf + 1
  end do
  !
  ! Face 2 - x/xi max face of quad (right/east)
  !
  nf = np+1
  do j = 1,np
    do i = 1,np
      k = (j-1)*np + np - i + 1
      return_value(k,nf) = cmat1D(i)
    end do
    nf = nf + 1
  end do
  !
  ! Face 3 - y/eta max face of quad (top/north)
  !
  nf = 2*np+1
  do j = 1,np
    do i = 1,np
      k = (np-i)*np + np - j + 1
      return_value(k,nf) = cmat1D(i)
    end do
    nf = nf + 1
  end do
  !
  ! Face 4 - x/xi min face of quad (left/west)
  !
  nf = 3*np+1
  do j = 1,np
    do i = 1,np
      k = (np-j)*np + i
      return_value(k,nf) = cmat1D(i)
    end do
    nf = nf + 1
  end do
  !
end function CorrectionMatrix_Quad
!
!###############################################################################
!
pure function CorrectionMatrix_Tria(rs,eval_deriv) result(return_value)
  !
  ! RS : r and s coordinates for the standard
  !      triangle element on the interval [-1,1]
  !
  !.. Use Statements ..
  use vandermonde_mod, only : vand
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: eval_deriv
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,n1,n2,np,nsp,nfp
  integer :: this_order,ierr
  logical(lk) :: compute_derivative
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: Fmask(:,:)
  real(qp), allocatable :: cmat_qp(:,:)
  real(qp), allocatable :: mmat1D_qp(:,:)
  real(qp), allocatable :: vedge_qp(:,:)
  real(qp), allocatable :: vtria_qp(:,:)
  !
continue
  !
  nsp = size(rs,dim=2)
  this_order = np2n(nsp)
  np = this_order+1
  nfp = 3*np
  !
  compute_derivative = true
  if (present(eval_deriv)) then
    compute_derivative = eval_deriv
  end if
  !
  ! Allocate the return array containing the working
  ! precision correction matrix
  !
  allocate ( return_value(1:nsp,1:nfp) , source=zero , stat=ierr )
  !
  ! Create the correction matrix using quad precision and allocate the
  ! inverse mass matrix and Fmask arrays.
  !
  allocate ( cmat_qp(1:nsp,1:nfp) , source=qzero , stat=ierr )
  allocate ( mmat1D_qp(1:np,1:np) , source=qzero , stat=ierr )
  allocate ( Fmask(1:np,1:3) , source=0 , stat=ierr )
  !
  ! Make quad precision copies of the Edge and Tria Vandermonde matrices.
  !
  allocate ( vedge_qp(1:np,1:np) , source=qzero , stat=ierr )
  allocate ( vtria_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  !
  vedge_qp = real( vand(Geom_Edge,this_order)%modal2nodal%mat , kind=qp )
  vtria_qp = real( vand(Geom_Tria,this_order)%modal2nodal%mat , kind=qp )
  !
  ! Face 1 - bottom edge of triangle
  !
  do i = 1,np
    Fmask(i,1) = i
  end do
  !
  ! Face 2 - diagonal edge of triangle
  !
  j = 0
  do i = 1,np
    j = j + np + 1 - i
    Fmask(i,2) = j
  end do
  !
  ! Face 3 - vertical edge of triangle
  !
  j = nsp + 1
  do i = 1,np
    j = j - i
    Fmask(i,3) = j
  end do
  !
  ! DG Book's method for getting Fmask
  !
 !Fmask(:,1) = pack( (/(i,i=1,nsp)/) , abs(rs(2,:) + one    ) < eps12 )
 !Fmask(:,2) = pack( (/(i,i=1,nsp)/) , abs(rs(1,:) + rs(2,:)) < eps12 )
 !Fmask(:,3) = pack( (/(i,i=1,nsp)/) , abs(rs(1,:) + one    ) < eps12 )
  !
  ! Compute the inverse mass matrix along an edge
  !
  mmat1D_qp = invert_matrix( matmul( vedge_qp , transpose(vedge_qp) ) )
  !
  ! Face 1
  !
  n1 = 1  ! index of first node on this face
  n2 = np ! index of last  node on this face
  !
  cmat_qp(Fmask(:,1),n1:n2) = mmat1D_qp
  !
  ! Face 2
  !
  n1 = np+1 ! index of first node on this face
  n2 = 2*np ! index of last  node on this face
  !
  cmat_qp(Fmask(:,2),n1:n2) = mmat1D_qp
  !
  ! Face 3
  !
  n1 = 2*np+1 ! index of first node on this face
  n2 = 3*np   ! index of last  node on this face
  !
  cmat_qp(Fmask(:,3),n1:n2) = mmat1D_qp
  !
  ! Finally, compute the correction matrix
  !
  cmat_qp = matmul( transpose(vtria_qp) , cmat_qp )
  cmat_qp = qtwo*matmul( vtria_qp , cmat_qp )
  !
  ! Use the chop function to copy cmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( cmat_qp )
  !
end function CorrectionMatrix_Tria
!
!###############################################################################
!
pure function CorrectionMatrix_Tetr(rs,eval_deriv) result(return_value)
  !
  ! RS : r and s coordinates for the standard
  !      triangle element on the interval [-1,1]
  !
  !.. Use Statements ..
  use vandermonde_mod, only : vand
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: eval_deriv
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,n1,n2,np,nsp,nfp
  integer :: this_order,ierr
  logical(lk) :: compute_derivative
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: Fmask(:,:)
  real(qp), allocatable :: cmat_qp(:,:)
  real(qp), allocatable :: mmat1D_qp(:,:)
  real(qp), allocatable :: vedge_qp(:,:)
  real(qp), allocatable :: vtria_qp(:,:)
  !
continue
  !
  nsp = size(rs,dim=2)
  this_order = np2n(nsp)
  np = this_order+1
  nfp = 3*np
  !
  compute_derivative = true
  if (present(eval_deriv)) then
    compute_derivative = eval_deriv
  end if
  !
  ! Allocate the return array containing the working
  ! precision correction matrix
  !
  allocate ( return_value(1:nsp,1:nfp) , source=zero , stat=ierr )
  !
  ! Create the correction matrix using quad precision and allocate the
  ! inverse mass matrix and Fmask arrays.
  !
  allocate ( cmat_qp(1:nsp,1:nfp) , source=qzero , stat=ierr )
  allocate ( mmat1D_qp(1:np,1:np) , source=qzero , stat=ierr )
  allocate ( Fmask(1:np,1:3) , source=0 , stat=ierr )
  !
  ! Make quad precision copies of the Edge and Tria Vandermonde matrices.
  !
  allocate ( vedge_qp(1:np,1:np) , source=qzero , stat=ierr )
  allocate ( vtria_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  !
  vedge_qp = real( vand(Geom_Edge,this_order)%modal2nodal%mat , kind=qp )
  vtria_qp = real( vand(Geom_Tria,this_order)%modal2nodal%mat , kind=qp )
  !
  ! Face 1 - bottom edge of triangle
  !
  do i = 1,np
    Fmask(i,1) = i
  end do
  !
  ! Face 2 - diagonal edge of triangle
  !
  j = 0
  do i = 1,np
    j = j + np + 1 - i
    Fmask(i,2) = j
  end do
  !
  ! Face 3 - vertical edge of triangle
  !
  j = nsp + 1
  do i = 1,np
    j = j - i
    Fmask(i,3) = j
  end do
  !
  ! DG Book's method for getting Fmask
  !
 !Fmask(:,1) = pack( (/(i,i=1,nsp)/) , abs(rs(2,:) + one    ) < eps12 )
 !Fmask(:,2) = pack( (/(i,i=1,nsp)/) , abs(rs(1,:) + rs(2,:)) < eps12 )
 !Fmask(:,3) = pack( (/(i,i=1,nsp)/) , abs(rs(1,:) + one    ) < eps12 )
  !
  ! Compute the inverse mass matrix along an edge
  !
  mmat1D_qp = invert_matrix( matmul( vedge_qp , transpose(vedge_qp) ) )
  !
  ! Face 1
  !
  n1 = 1  ! index of first node on this face
  n2 = np ! index of last  node on this face
  !
  cmat_qp(Fmask(:,1),n1:n2) = mmat1D_qp
  !
  ! Face 2
  !
  n1 = np+1 ! index of first node on this face
  n2 = 2*np ! index of last  node on this face
  !
  cmat_qp(Fmask(:,2),n1:n2) = mmat1D_qp
  !
  ! Face 3
  !
  n1 = 2*np+1 ! index of first node on this face
  n2 = 3*np   ! index of last  node on this face
  !
  cmat_qp(Fmask(:,3),n1:n2) = mmat1D_qp
  !
  ! Finally, compute the correction matrix
  !
  cmat_qp = matmul( transpose(vtria_qp) , cmat_qp )
  cmat_qp = qtwo*matmul( vtria_qp , cmat_qp )
  !
  ! Use the chop function to copy cmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( cmat_qp )
  !
end function CorrectionMatrix_Tetr
!
!###############################################################################
!
pure function CorrectionMatrix_Pyra(rs,eval_deriv) result(return_value)
  !
  ! RS : r and s coordinates for the standard
  !      triangle element on the interval [-1,1]
  !
  !.. Use Statements ..
  use vandermonde_mod, only : vand
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: eval_deriv
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,n1,n2,np,nsp,nfp
  integer :: this_order,ierr
  logical(lk) :: compute_derivative
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: Fmask(:,:)
  real(qp), allocatable :: cmat_qp(:,:)
  real(qp), allocatable :: mmat1D_qp(:,:)
  real(qp), allocatable :: vedge_qp(:,:)
  real(qp), allocatable :: vtria_qp(:,:)
  !
continue
  !
  nsp = size(rs,dim=2)
  this_order = np2n(nsp)
  np = this_order+1
  nfp = 3*np
  !
  compute_derivative = true
  if (present(eval_deriv)) then
    compute_derivative = eval_deriv
  end if
  !
  ! Allocate the return array containing the working
  ! precision correction matrix
  !
  allocate ( return_value(1:nsp,1:nfp) , source=zero , stat=ierr )
  !
  ! Create the correction matrix using quad precision and allocate the
  ! inverse mass matrix and Fmask arrays.
  !
  allocate ( cmat_qp(1:nsp,1:nfp) , source=qzero , stat=ierr )
  allocate ( mmat1D_qp(1:np,1:np) , source=qzero , stat=ierr )
  allocate ( Fmask(1:np,1:3) , source=0 , stat=ierr )
  !
  ! Make quad precision copies of the Edge and Tria Vandermonde matrices.
  !
  allocate ( vedge_qp(1:np,1:np) , source=qzero , stat=ierr )
  allocate ( vtria_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  !
  vedge_qp = real( vand(Geom_Edge,this_order)%modal2nodal%mat , kind=qp )
  vtria_qp = real( vand(Geom_Tria,this_order)%modal2nodal%mat , kind=qp )
  !
  ! Face 1 - bottom edge of triangle
  !
  do i = 1,np
    Fmask(i,1) = i
  end do
  !
  ! Face 2 - diagonal edge of triangle
  !
  j = 0
  do i = 1,np
    j = j + np + 1 - i
    Fmask(i,2) = j
  end do
  !
  ! Face 3 - vertical edge of triangle
  !
  j = nsp + 1
  do i = 1,np
    j = j - i
    Fmask(i,3) = j
  end do
  !
  ! DG Book's method for getting Fmask
  !
 !Fmask(:,1) = pack( (/(i,i=1,nsp)/) , abs(rs(2,:) + one    ) < eps12 )
 !Fmask(:,2) = pack( (/(i,i=1,nsp)/) , abs(rs(1,:) + rs(2,:)) < eps12 )
 !Fmask(:,3) = pack( (/(i,i=1,nsp)/) , abs(rs(1,:) + one    ) < eps12 )
  !
  ! Compute the inverse mass matrix along an edge
  !
  mmat1D_qp = invert_matrix( matmul( vedge_qp , transpose(vedge_qp) ) )
  !
  ! Face 1
  !
  n1 = 1  ! index of first node on this face
  n2 = np ! index of last  node on this face
  !
  cmat_qp(Fmask(:,1),n1:n2) = mmat1D_qp
  !
  ! Face 2
  !
  n1 = np+1 ! index of first node on this face
  n2 = 2*np ! index of last  node on this face
  !
  cmat_qp(Fmask(:,2),n1:n2) = mmat1D_qp
  !
  ! Face 3
  !
  n1 = 2*np+1 ! index of first node on this face
  n2 = 3*np   ! index of last  node on this face
  !
  cmat_qp(Fmask(:,3),n1:n2) = mmat1D_qp
  !
  ! Finally, compute the correction matrix
  !
  cmat_qp = matmul( transpose(vtria_qp) , cmat_qp )
  cmat_qp = qtwo*matmul( vtria_qp , cmat_qp )
  !
  ! Use the chop function to copy cmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( cmat_qp )
  !
end function CorrectionMatrix_Pyra
!
!###############################################################################
!
pure function CorrectionMatrix_Pris(rs,eval_deriv) result(return_value)
  !
  ! RS : r and s coordinates for the standard
  !      triangle element on the interval [-1,1]
  !
  !.. Use Statements ..
  use vandermonde_mod, only : vand
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: eval_deriv
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,n1,n2,np,nsp,nfp
  integer :: this_order,ierr
  logical(lk) :: compute_derivative
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: Fmask(:,:)
  real(qp), allocatable :: cmat_qp(:,:)
  real(qp), allocatable :: mmat1D_qp(:,:)
  real(qp), allocatable :: vedge_qp(:,:)
  real(qp), allocatable :: vtria_qp(:,:)
  !
continue
  !
  nsp = size(rs,dim=2)
  this_order = np2n(nsp)
  np = this_order+1
  nfp = 3*np
  !
  compute_derivative = true
  if (present(eval_deriv)) then
    compute_derivative = eval_deriv
  end if
  !
  ! Allocate the return array containing the working
  ! precision correction matrix
  !
  allocate ( return_value(1:nsp,1:nfp) , source=zero , stat=ierr )
  !
  ! Create the correction matrix using quad precision and allocate the
  ! inverse mass matrix and Fmask arrays.
  !
  allocate ( cmat_qp(1:nsp,1:nfp) , source=qzero , stat=ierr )
  allocate ( mmat1D_qp(1:np,1:np) , source=qzero , stat=ierr )
  allocate ( Fmask(1:np,1:3) , source=0 , stat=ierr )
  !
  ! Make quad precision copies of the Edge and Tria Vandermonde matrices.
  !
  allocate ( vedge_qp(1:np,1:np) , source=qzero , stat=ierr )
  allocate ( vtria_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  !
  vedge_qp = real( vand(Geom_Edge,this_order)%modal2nodal%mat , kind=qp )
  vtria_qp = real( vand(Geom_Tria,this_order)%modal2nodal%mat , kind=qp )
  !
  ! Face 1 - bottom edge of triangle
  !
  do i = 1,np
    Fmask(i,1) = i
  end do
  !
  ! Face 2 - diagonal edge of triangle
  !
  j = 0
  do i = 1,np
    j = j + np + 1 - i
    Fmask(i,2) = j
  end do
  !
  ! Face 3 - vertical edge of triangle
  !
  j = nsp + 1
  do i = 1,np
    j = j - i
    Fmask(i,3) = j
  end do
  !
  ! DG Book's method for getting Fmask
  !
 !Fmask(:,1) = pack( (/(i,i=1,nsp)/) , abs(rs(2,:) + one    ) < eps12 )
 !Fmask(:,2) = pack( (/(i,i=1,nsp)/) , abs(rs(1,:) + rs(2,:)) < eps12 )
 !Fmask(:,3) = pack( (/(i,i=1,nsp)/) , abs(rs(1,:) + one    ) < eps12 )
  !
  ! Compute the inverse mass matrix along an edge
  !
  mmat1D_qp = invert_matrix( matmul( vedge_qp , transpose(vedge_qp) ) )
  !
  ! Face 1
  !
  n1 = 1  ! index of first node on this face
  n2 = np ! index of last  node on this face
  !
  cmat_qp(Fmask(:,1),n1:n2) = mmat1D_qp
  !
  ! Face 2
  !
  n1 = np+1 ! index of first node on this face
  n2 = 2*np ! index of last  node on this face
  !
  cmat_qp(Fmask(:,2),n1:n2) = mmat1D_qp
  !
  ! Face 3
  !
  n1 = 2*np+1 ! index of first node on this face
  n2 = 3*np   ! index of last  node on this face
  !
  cmat_qp(Fmask(:,3),n1:n2) = mmat1D_qp
  !
  ! Finally, compute the correction matrix
  !
  cmat_qp = matmul( transpose(vtria_qp) , cmat_qp )
  cmat_qp = qtwo*matmul( vtria_qp , cmat_qp )
  !
  ! Use the chop function to copy cmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( cmat_qp )
  !
end function CorrectionMatrix_Pris
!
!###############################################################################
!
pure function CorrectionMatrix_Hexa(xi,eval_deriv) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xi
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: eval_deriv
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: nf,nep,nfp,nbp,nsp
  integer :: i,j,k,l,ierr
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: cmat1D(:)
  !
continue
  !
  ! #########################################################################
  ! #########################################################################
  ! #########################################################################
  ! #########################################################################
  ! #####  NOTE: THE ORDER OF THE FLUX POINTS IS SUCH THAT THE LOWEST   #####
  ! #####        COORDINATE DIRECTION IS THE INNER MOST LOOP, i.e.,     #####
  ! #####          FOR AN X-FACE, INNER LOOP OVER Y, OUTER LOOP OVER Z  #####
  ! #####          FOR A  Y-FACE, INNER LOOP OVER X, OUTER LOOP OVER Z  #####
  ! #####          FOR A  Z-FACE, INNER LOOP OVER X, OUTER LOOP OVER Y  #####
  ! #########################################################################
  ! #########################################################################
  ! #########################################################################
  ! #########################################################################
  !
  ! ###########################################################################
  ! ###########################################################################
  ! ####  NOTE: No need to use chopping function in this procedure. The    ####
  ! ####        hexa correction matrix is really just the 1D correction    ####
  ! ####        matrix applied to each face/flux point of the standard     ####
  ! ####        hexa element. The 1D correction matrix is obtained by the  ####
  ! ####        function CorrectionMatrix_Edge which applies the chopping  ####
  ! ####        function to its result so there is no need to repeat this  ####
  ! ####        operation again.                                           ####
  ! ###########################################################################
  ! ###########################################################################
  !
  nep = size(xi)
  nfp = nep*nep
  nbp = 6*nfp
  nsp = nep*nep*nep
  !
  ! Allocate the return array containing the working
  ! precision correction matrix
  !
  allocate ( return_value(1:nsp,1:nbp) , source=zero , stat=ierr )
  !
  ! Create the 1D correction matrix
  !
  allocate ( cmat1D(1:nep) , source=zero , stat=ierr )
  !
  ! 1D Correction matrix
  !
  cmat1D = CorrectionMatrix_Edge(xi,eval_deriv)
  !
  ! ############################################################
  ! ############################################################
  ! #####  Face 1 - z/zeta min face of quad (bottom/down)  #####
  ! ############################################################
  ! ############################################################
  !
  nf = 0*nfp+1 ! First point on this face
  !
  ! Loops over points on this face
  !
  do j = 1,nep
    do i = 1,nep
      !
      ! Loop over interior solution points in z/zeta direction
      !
      do k = 1,nep
        l = i + (j-1)*nep + (k-1)*nfp
        return_value(l,nf) = cmat1D(k)
      end do
      !
      nf = nf + 1 ! Proceed to next point on this face
      !
    end do
  end do
  !
  ! #######################################################
  ! #######################################################
  ! #####  Face 2 - z/zeta max face of quad (top/up)  #####
  ! #######################################################
  ! #######################################################
  !
  nf = 1*nfp+1 ! First point on this face
  !
  ! Loops over points on this face
  !
  do j = 1,nep
    do i = 1,nep
      !
      ! Loop over interior solution points in z/zeta direction
      !
      do k = 1,nep
        l = i + (j-1)*nep + nfp*(nep-k)
        return_value(l,nf) = cmat1D(k)
      end do
      !
      nf = nf + 1 ! Proceed to next point on this face
      !
    end do
  end do
  !
  ! #########################################################
  ! #########################################################
  ! #####  Face 3 - x/xi max face of quad (right/east)  #####
  ! #########################################################
  ! #########################################################
  !
  nf = 2*nfp+1 ! First point on this face
  !
  ! Loops over points on this face
  !
  do k = 1,nep
    do j = 1,nep
      !
      ! Loop over interior solution points in x/xi direction
      !
      do i = 1,nep
        l = j*nep + (k-1)*nfp - i + 1
        return_value(l,nf) = cmat1D(i)
      end do
      !
      nf = nf + 1 ! Proceed to next point on this face
      !
    end do
  end do
  !
  ! ########################################################
  ! ########################################################
  ! #####  Face 4 - x/xi min face of quad (left/west)  #####
  ! ########################################################
  ! ########################################################
  !
  nf = 3*nfp+1 ! First point on this face
  !
  ! Loops over points on this face
  !
  do k = 1,nep
    do j = 1,nep
      !
      ! Loop over interior solution points in x/xi direction
      !
      do i = 1,nep
        l = (j-1)*nep + (k-1)*nfp + i
        return_value(l,nf) = cmat1D(i)
      end do
      !
      nf = nf + 1 ! Proceed to next point on this face
      !
    end do
  end do
  !
  ! ###########################################################
  ! ###########################################################
  ! #####  Face 5 - y/eta min face of quad (front/south)  #####
  ! ###########################################################
  ! ###########################################################
  !
  nf = 4*nfp+1 ! First point on this face
  !
  ! Loops over points on this face
  !
  do k = 1,nep
    do i = 1,nep
      !
      ! Loop over interior solution points in y/eta direction
      !
      do j = 1,nep
        l = (k-1)*nfp + i + (j-1)*nep
        return_value(l,nf) = cmat1D(j)
      end do
      !
      nf = nf + 1 ! Proceed to next point on this face
      !
    end do
  end do
  !
  ! ##########################################################
  ! ##########################################################
  ! #####  Face 6 - y/eta max face of quad (back/north)  #####
  ! ##########################################################
  ! ##########################################################
  !
  nf = 5*nfp+1 ! First point on this face
  !
  ! Loops over points on this face
  !
  do k = 1,nep
    do i = 1,nep
      !
      ! Loop over interior solution points in y/eta direction
      !
      do j = 1,nep
        l = (k-1)*nfp + i + (nep-j)*nep
        return_value(l,nf) = cmat1D(j)
      end do
      !
      nf = nf + 1 ! Proceed to next point on this face
      !
    end do
  end do
  !
end function CorrectionMatrix_Hexa
!
!###############################################################################
!
subroutine correction_memory_usage(iunit)
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
  ! Memory usage for the correct array of type matrix
  !
  if (allocated(correct)) then
    !
    call dt%set( storage_size(correct,kind=inttype) , &
                         size(correct,kind=inttype) )
    !
    do j = lbound(correct,dim=2),ubound(correct,dim=2)
      do i = lbound(correct,dim=1),ubound(correct,dim=1)
        !
        array = correct(i,j)%memory_usage()
        !
        write (array_name,4) "correct",Geom_Name(i),j
        write (iou,2) array_name,array%get_units()
        !
        call dt%add(array)
        !
      end do
    end do
    !
    call total%add(dt)
    write (array_name,5) "correct"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "correct"
    !
  end if
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module correction_mod")
  2 format (" Array: ",a,dt)
  3 format (" Total memory = ",f8.2,1x,a)
  4 format (a,"(",a,",",i0,")%mat")
  5 format ("'",a,"'")
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  !
#endif
end subroutine correction_memory_usage
!
!###############################################################################
!
include "Functions/np2n.f90"
!
!###############################################################################
!
end module correction_mod
!
!###############################################################################
!
