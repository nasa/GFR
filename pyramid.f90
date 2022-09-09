module pyramid_mod
  !
  ! Module that contains various functions and routines related to
  ! creating the prism elements and the solution nodes within
  !
  use module_kind_types, only : MP,LK,true,fals
  use module_kind_types, only : Legendre_Gauss_Lobatto
  use module_kind_types, only : Equi_Distant_Points
#ifdef SPECIAL_FOR_INTEL
  use module_kind_types, only : invert_matrix_lapack95
#else
  use module_kind_types, only : invert_matrix
#endif
  !
  implicit none
  !
  !.. Module Precision ..
  integer, parameter :: lp = MP
  !
  private
  !
  public :: PyrNodes3D
  public :: get_pyramid_weights_at_solpts
  !
  !.. Module Parameter Constants ..
  real(lp), parameter :: zero   = 0.0_lp
  real(lp), parameter :: half   = 0.5_lp
  real(lp), parameter :: one    = 1.0_lp
  real(lp), parameter :: eight3 = 8.0_lp/3.0_lp
  real(lp), parameter :: eps10  = 1.0e-10_lp
  !
contains
!
!###############################################################################
!
pure integer function pyr_pts(p)
  integer, intent(in) :: p
  pyr_pts = (p+1)*(p+2)*(2*p+3)/6
end function pyr_pts
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function PyrNodes3D(p,location) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : nodes_legendre_gauss
  use polynomial_mod, only : nodes_legendre_gauss_lobatto
  use polynomial_mod, only : equal_distant_quadrature
  !
  !.. Formal Arguments ..
  integer, intent(in) :: p
  integer, intent(in) :: location
  !
  !.. Function Result..
  real(lp) :: return_value(1:3,1:pyr_pts(p))
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,n
  !
  !.. Local Arrays ..
  real(lp) :: zeta(1:p+1)
  !
  !.. Local Allocatable Arrays ..
  real(lp), allocatable :: xi(:)
  !
continue
  !
  if (p == 0) then
    return_value(:,:) = reshape( [zero,zero,-half], [3,1] )
    return
  end if
  !
  ! Get the 1D Legendre-Gauss-Lobatto nodes
  !
  if (location == Legendre_Gauss_Lobatto) then
    call nodes_legendre_gauss_lobatto(p,zeta)
  else if (location == Equi_Distant_Points) then
    call equal_distant_quadrature(zeta)
  else
    call nodes_legendre_gauss(p+1,zeta)
  end if
  !
  n = 0
  do k = 1,size(zeta)-1
    !
    xi = spread(zero,dim=1,ncopies=p-k+2) ! F2003 AUTO-REALLOC
    !
    if (location == Legendre_Gauss_Lobatto) then
      call nodes_legendre_gauss_lobatto(p-k+1,xi)
    else if (location == Equi_Distant_Points) then
      call equal_distant_quadrature(xi)
    else
      call nodes_legendre_gauss(p-k+2,xi)
    end if
    !
    ! Scale xi to the current zeta value
    !
    xi(:) = half*(one-zeta(k)) * xi(:)
    !
    do j = 1,size(xi)
      do i = 1,size(xi)
        n = n + 1
        return_value(1,n) = xi(i)
        return_value(2,n) = xi(j)
        return_value(3,n) = zeta(k)
      end do
    end do
    !
  end do
  !
  n = n + 1
  return_value(1,n) = zero
  return_value(2,n) = zero
  return_value(3,n) = zeta(size(zeta))
  !
end function PyrNodes3D
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_pyramid_weights_at_solpts(rst,nquadrature) result(return_value)
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! NOTE: THIS FUNCTION IS KNOWN TO COMPUTE INCORRECT QUADRATURE WEIGHTS!!!
  !       USE FUNCTION get_generic_quadrature_weights FROM THE MODULE
  !       quadrature_mod INSTEAD!!! 
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  !.. Use Statements ..
  use polynomial_mod, only : legendre_gauss_quadrature
  !
  !.. Formal Arguments ..
  real(lp), intent(in) :: rst(:,:)
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: nquadrature
  !
  !.. Function Result ..
  real(lp) :: return_value(1:size(rst,2))
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,m,n,p,nquad
  integer  :: ix,iy,iz
  integer  :: nspts,nfpts
  real(lp) :: x,y,z,wx,wy,wz,a,b,c
  !
  !.. Local Arrays ..
  real(lp) :: V(1:size(rst,2),1:size(rst,2))
  real(lp) :: xyz(1:3,1:size(rst,2))
  !
  !.. Local Allocatable Arrays ..
  real(lp), allocatable :: xi(:)
  real(lp), allocatable :: wts(:)
  !
continue
  !
  ! On entry, r and s need to be defined on the interval [-1,1]
  !
  nspts = size(rst,2)
  do n = 0,nspts
    if (pyr_pts(n) >= nspts) then
      p = n
      exit
    end if
  end do
  nfpts = p+1
  !
  nquad = 2*nfpts - 1 ! quadtrature order is 2*P
  !
  if (present(nquadrature)) then
    nquad = nquadrature
  end if
  !
  xi = spread(zero,dim=1,ncopies=nquad) ! F2003 AUTO-REALLOC
  wts = spread(zero,dim=1,ncopies=nquad) ! F2003 AUTO-REALLOC
  !
  ! Compute the coefficients of the Vandermonde matrix for
  ! the Lagrange polynomials using monomial basis functions
  !
  xyz(:,:) = half * (rst(:,:) + one)
  !
  do n = 1,nspts
    l = 0
    x = xyz(1,n)
    y = xyz(2,n)
    z = xyz(3,n)
    do i = 0,p
      do j = 0,i
        do k = 0,j
          l = l + 1
          V(l,n) = x**(i-j) * y**(j-k) * z**(k)
        end do
      end do
    end do
    do m = p+1,2*p
      do i = p,m-p,-1
        do j = m-i,m-p,-1
          l = l + 1
          k = m-i-j
          V(l,n) = x**(i) * y**(j) * z**(k)
        end do
      end do
    end do
  end do
  !
  ! Now get the inverse of this matrix for use later
  !
#ifdef SPECIAL_FOR_INTEL
  V = invert_matrix_lapack95(V)
#else
  V = invert_matrix(V)
#endif
  !
  ! Get the Gauss quadrature points on the interval [-1,1] and
  ! the corresponding weights necessary to perform the integration
  !
  call legendre_gauss_quadrature(xi,wts)
  !
  ! Reset the quadrature points to the interval [0,1]
  !
  xi(:) = half*(one+xi(:))
  !
  ! Initialize the pyramid weights to zero
  !
  return_value(:) = zero
  !
  ! Loop over each solution point performing a double integral on the
  ! corresponding Lagrange polynomial over the domain of the pyramid
  ! i.e., n    : current solution point
  !       L(n) : Lagrange polynomial for n
  !
  !   /\ 1      /\ 1-zeta/2 /\ 1-zeta/2
  !   |         |           |
  !   | d(zeta) | d(eta)    | d(xi) [L(n)]
  !   |         |           |
  !  \/ 0      \/ zeta/2   \/ zeta/2
  !
  do n = 1,nspts
    !
    ! Outer integral loop
    !
    z_integral: &
    do iz = 1,nquad
      !
      z = xi(iz)
      wz = wts(iz)
      !
      ! Middle integral loop
      !
      a = zero
      y_integral: &
      do iy = 1,nquad
        !
        ! Multiply eta by (1-z), then add z/2 to it,
        ! to scale eta from [0,1] to [zeta/2,1-zeta/2]
        !
        ! Multiply wts by (1-z) to scale it as well
        !
        y = (one-z) * xi(iy) + half*z
        wy = (one-z) * wts(iy)
        !
        ! Inner integral loop
        !
        b = zero
        x_integral: &
        do ix = 1,nquad
          !
          ! Multiply xi by (1-z), then add z/2 to it,
          ! to scale xi from [0,1] to [zeta/2,1-zeta/2]
          !
          ! Multiply wts by (1-z) to scale it as well
          !
          x = (one-z) * xi(ix) + half*z
          wx = (one-z) * wts(ix)
          !
          ! Evaluate the Lagrange polynomial at the point x,y,z
          !
          l = 0
          c = zero
          do i = 0,p
            do j = 0,i
              do k = 0,j
                l = l + 1
                c = c + V(n,l) * x**(i-j) * y**(j-k) * z**(k)
              end do
            end do
          end do
          do m = p+1,2*p
            do i = p,m-p,-1
              do j = m-i,m-p,-1
                l = l + 1
                k = m-i-j
                c = c + V(n,l) * x**(i) * y**(j) * z**(k)
              end do
            end do
          end do
          !
          ! Add the weighted polynomial value to the inner integral result
          !
          b = b + wx*c
          !
        end do x_integral
        !
        a = a + wy*b
        !
      end do y_integral
      !
      ! Add the weighted inner integral result to the outer integral result
      !
      return_value(n) = return_value(n) + wz*a
      !
    end do z_integral
    !
  end do
  !
  ! Zero out any weights that are smaller than eps10
  !
  return_value(:) = merge(return_value(:),zero,abs(return_value(:))>eps10)
  !
  ! Make sure the weights sum to 8/3
  !
 !return_value(:) = eight3 * return_value(:) / sum(return_value(:))
  !
end function get_pyramid_weights_at_solpts
!
!###############################################################################
!
end module pyramid_mod
!
!###############################################################################
!
