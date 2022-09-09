module prism_mod
  !
  ! Module that contains various functions and routines related to
  ! creating the prism elements and the solution nodes within
  !
  use module_kind_types, only : MP
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
  public :: PriNodes3D
  public :: get_prism_weights_at_solpts
  !
  !.. Module Parameter Constants ..
  real(lp), parameter :: zero  = 0.0_lp
  real(lp), parameter :: one3  = 1.0_lp/3.0_lp
  real(lp), parameter :: half  = 0.5_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: four  = 4.0_lp
  real(lp), parameter :: eps10 = 1.0e-10_lp
  !
contains
!
!###############################################################################
!
pure integer function pri_pts(p)
  integer, intent(in) :: p
  pri_pts = (p+1)*(p+1)*(p+2)/2
end function pri_pts
!
!###############################################################################
!
pure integer function tri_pts(p)
  integer, intent(in) :: p
  tri_pts = (p+1)*(p+2)/2
end function tri_pts
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function PriNodes3D(p,tria_loc,edge_loc,tria_edge_loc) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : nodes_legendre_gauss
  use polynomial_mod, only : nodes_legendre_gauss_lobatto
  use polynomial_mod, only : equal_distant_quadrature
  use triangle_mod,   only : TriNodes2D
  !
  !.. Formal Arguments ..
  integer, intent(in) :: p
  integer, intent(in) :: tria_loc
  integer, intent(in) :: edge_loc
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: tria_edge_loc
  !
  !.. Function Result..
  real(lp) :: return_value(1:3,1:pri_pts(p))
  !
  !.. Local Scalars ..
  integer :: i,j,n
  !
  !.. Local Arrays ..
  real(lp) :: xi(1:p+1)
  real(lp) :: rs(1:2,1:tri_pts(p))
  !
continue
  !
  if (p == 0) then
    return_value(:,:) = reshape( [-one3,-one3,zero], [3,1] )
    return
  end if
  !
  ! Get the 1D Legendre-Gauss-Lobatto nodes
  !
  if (edge_loc == Legendre_Gauss_Lobatto) then
    call nodes_legendre_gauss_lobatto(p,xi)
  else if (edge_loc == Equi_Distant_Points) then
    call equal_distant_quadrature(xi)
  else
    call nodes_legendre_gauss(p+1,xi)
  end if
  !
  ! Get the 2D Barycentric Gauss triangle nodes
  !
  rs(:,:) = TriNodes2D(p,tria_loc,tria_edge_loc)
  !
  ! Create a tensor product of these two sets of points
  !
  n = 0
  do j = 1,size(xi)
    do i = 1,size(rs,dim=2)
      n = n + 1
      return_value(1,n) = rs(1,i)
      return_value(2,n) = rs(2,i)
      return_value(3,n) = xi(j)
    end do
  end do
  !
end function PriNodes3D
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_prism_weights_at_solpts(rst) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : legendre_gauss_quadrature
  !
  !.. Formal Arguments ..
  real(lp), intent(in) :: rst(:,:)
  !
  !.. Function Result ..
  real(lp) :: return_value(1:size(rst,2))
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,n,nquad
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
    if (pri_pts(n) >= nspts) then
      nfpts = n+1
      exit
    end if
  end do
  !
  nquad = 2*nfpts - 1 ! quadtrature order is 2*P
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
    do i = 1,nfpts
      do j = 1,i
        do k = 1,nfpts
          l = l + 1
          V(l,n) = x**(i-j) * y**(j-1) * z**(k-1)
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
  ! Initialize the prism weights to zero
  !
  return_value(:) = zero
  !
  ! Loop over each solution point performing a double integral on the
  ! corresponding Lagrange polynomial over the domain of the prism
  ! i.e., n    : current solution point
  !       L(n) : Lagrange polynomial for n
  !
  !   /\ 1    /\ 1-xi  /\ 1
  !   |       |        |
  !   | d(xi) | d(eta) | d(zeta) [L(n)]
  !   |       |        |
  !  \/ 0    \/ 0     \/ 0
  !
  do n = 1,nspts
    !
    ! Outer integral loop
    !
    x_integral: &
    do ix = 1,nquad
      !
      x = xi(ix)
      wx = wts(ix)
      !
      ! Middle integral loop
      !
      a = zero
      y_integral: &
      do iy = 1,nquad
        !
        ! Multiply eta and wts by (1-x) to scale eta from [0,1] to [0,1-xi]
        !
        y = (one-x) * xi(iy)
        wy = (one-x) * wts(iy)
        !
        ! Inner integral loop
        !
        b = zero
        z_integral: &
        do iz = 1,nquad
          !
          z = xi(iz)
          wz = wts(iz)
          !
          ! Evaluate the Lagrange polynomial at the point x,y,z
          !
          l = 0
          c = zero
          do i = 1,nfpts
            do j = 1,i
              do k = 1,nfpts
                l = l + 1
                c = c + V(n,l) * x**(i-j) * y**(j-1) * z**(k-1)
              end do
            end do
          end do
          !
          ! Add the weighted polynomial value to the inner integral result
          !
          b = b + wz*c
          !
        end do z_integral
        !
        a = a + wy*b
        !
      end do y_integral
      !
      ! Add the weighted inner integral result to the outer integral result
      !
      return_value(n) = return_value(n) + wx*a
      !
    end do x_integral
    !
  end do
  !
  ! Zero out any weights that are smaller than eps10
  !
  return_value(:) = merge(return_value(:),zero,abs(return_value(:))>eps10)
  !
  ! Make sure the weights sum to 4
  !
  return_value(:) = four * return_value(:) / sum(return_value(:))
  !
end function get_prism_weights_at_solpts
!
!###############################################################################
!
end module prism_mod
!
!###############################################################################
!
