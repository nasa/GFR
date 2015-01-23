module polynomial_mod
  !
  ! Module that contains various functions and routines
  ! related to Jacobi, Legendre, and other polynomials
  !
  use module_kind_types
  !
  implicit none
  !
  private
  !
  public :: legendre_gauss_quadrature
  public :: legendre_gauss_lobatto_quadrature
  public :: JacobiGL
  public :: eval_jacobi_poly
  public :: eval_legendre_poly
  public :: eval_chebyshev_poly
  public :: eval_normalized_jacobi_poly
  public :: normalized_jacobi_poly
  public :: grad_normalized_jacobi_poly
  public :: nodes_jacobi_gauss
  public :: nodes_legendre_gauss
  public :: nodes_right_radau
  public :: nodes_chebyshev_gauss
  public :: nodes_jacobi_gauss_lobatto
  public :: nodes_legendre_gauss_lobatto
  public :: nodes_chebyshev_gauss_lobatto
  public :: weights_jacobi_gauss
  public :: weights_legendre_gauss
  public :: weights_chebyshev_gauss
  public :: weights_jacobi_gauss_lobatto
  public :: weights_legendre_gauss_lobatto
  public :: weights_right_radau
  public :: weights_chebyshev_gauss_lobatto
  public :: weights_dgbook_jacobi_gauss_lobatto
  public :: derivative_matrix_legendre_gauss
  public :: derivative_matrix_legendre_gauss_lobatto
  public :: derivative_matrix_jacobi_gauss_lobatto
  public :: derivative_matrix_chebyshev_gauss_lobatto
  public :: eval_LagrangePoly
  public :: eval_LagrangePoly1D
  public :: eval_LagrangePoly2D
  public :: eval_D2LagrangeDx2
  public :: eval_DLagrangeDx
  !
  interface JacobiGL
    module procedure JacobiGL_DP, & 
                     JacobiGL_QP
  end interface JacobiGL
  !
  interface eval_jacobi_poly
    module procedure eval_jacobi_poly_DP, & 
                     eval_jacobi_poly_QP
  end interface eval_jacobi_poly
  !
  interface eval_legendre_poly
    module procedure eval_legendre_poly_DP, & 
                     eval_legendre_poly_QP
  end interface eval_legendre_poly
  !
  interface eval_chebyshev_poly
    module procedure eval_chebyshev_poly_DP, & 
                     eval_chebyshev_poly_QP
  end interface eval_chebyshev_poly
  !
  interface eval_normalized_jacobi_poly
    module procedure eval_normalized_jacobi_poly_DP, & 
                     eval_normalized_jacobi_poly_QP
  end interface eval_normalized_jacobi_poly
  !
  interface normalized_jacobi_poly
    module procedure normalized_jacobi_poly_DP, & 
                     normalized_jacobi_poly_QP
  end interface normalized_jacobi_poly
  !
  interface grad_normalized_jacobi_poly
    module procedure grad_normalized_jacobi_poly_DP, & 
                     grad_normalized_jacobi_poly_QP
  end interface grad_normalized_jacobi_poly
  !
  interface nodes_jacobi_gauss
    module procedure nodes_jacobi_gauss_DP, & 
                     nodes_jacobi_gauss_QP
  end interface nodes_jacobi_gauss
  !
  interface nodes_legendre_gauss
    module procedure nodes_legendre_gauss_DP, & 
                     nodes_legendre_gauss_QP
  end interface nodes_legendre_gauss
  !
  interface nodes_right_radau
    module procedure nodes_right_radau_DP, & 
                     nodes_right_radau_QP
  end interface nodes_right_radau
  !
  interface nodes_chebyshev_gauss
    module procedure nodes_chebyshev_gauss_DP, & 
                     nodes_chebyshev_gauss_QP
  end interface nodes_chebyshev_gauss
  !
  interface nodes_jacobi_gauss_lobatto
    module procedure nodes_jacobi_gauss_lobatto_DP, & 
                     nodes_jacobi_gauss_lobatto_QP
  end interface nodes_jacobi_gauss_lobatto
  !
  interface nodes_legendre_gauss_lobatto
    module procedure nodes_legendre_gauss_lobatto_DP, & 
                     nodes_legendre_gauss_lobatto_QP
  end interface nodes_legendre_gauss_lobatto
  !
  interface nodes_chebyshev_gauss_lobatto
    module procedure nodes_chebyshev_gauss_lobatto_DP, & 
                     nodes_chebyshev_gauss_lobatto_QP
  end interface nodes_chebyshev_gauss_lobatto
  !
  interface weights_jacobi_gauss
    module procedure weights_jacobi_gauss_DP, & 
                     weights_jacobi_gauss_QP
  end interface weights_jacobi_gauss
  !
  interface weights_legendre_gauss
    module procedure weights_legendre_gauss_DP, & 
                     weights_legendre_gauss_QP
  end interface weights_legendre_gauss
  !
  interface weights_chebyshev_gauss
    module procedure weights_chebyshev_gauss_DP, & 
                     weights_chebyshev_gauss_QP
  end interface weights_chebyshev_gauss
  !
  interface weights_jacobi_gauss_lobatto
    module procedure weights_jacobi_gauss_lobatto_DP, & 
                     weights_jacobi_gauss_lobatto_QP
  end interface weights_jacobi_gauss_lobatto
  !
  interface weights_legendre_gauss_lobatto
    module procedure weights_legendre_gauss_lobatto_DP, & 
                     weights_legendre_gauss_lobatto_QP
  end interface weights_legendre_gauss_lobatto
  !
  interface weights_right_radau
    module procedure weights_right_radau_DP, & 
                     weights_right_radau_QP
  end interface weights_right_radau
  !
  interface weights_chebyshev_gauss_lobatto
    module procedure weights_chebyshev_gauss_lobatto_DP, & 
                     weights_chebyshev_gauss_lobatto_QP
  end interface weights_chebyshev_gauss_lobatto
  !
  interface weights_dgbook_jacobi_gauss_lobatto
    module procedure weights_dgbook_jacobi_gauss_lobatto_DP, & 
                     weights_dgbook_jacobi_gauss_lobatto_QP
  end interface weights_dgbook_jacobi_gauss_lobatto
  !
  interface derivative_matrix_legendre_gauss
    module procedure derivative_matrix_legendre_gauss_DP, & 
                     derivative_matrix_legendre_gauss_QP
  end interface derivative_matrix_legendre_gauss
  !
  interface derivative_matrix_legendre_gauss_lobatto
    module procedure derivative_matrix_legendre_gauss_lobatto_DP, & 
                     derivative_matrix_legendre_gauss_lobatto_QP
  end interface derivative_matrix_legendre_gauss_lobatto
  !
  interface derivative_matrix_jacobi_gauss_lobatto
    module procedure derivative_matrix_jacobi_gauss_lobatto_DP, & 
                     derivative_matrix_jacobi_gauss_lobatto_QP
  end interface derivative_matrix_jacobi_gauss_lobatto
  !
  interface derivative_matrix_chebyshev_gauss_lobatto
    module procedure derivative_matrix_chebyshev_gauss_lobatto_DP, & 
                     derivative_matrix_chebyshev_gauss_lobatto_QP
  end interface derivative_matrix_chebyshev_gauss_lobatto
  !
  interface gammafun
    module procedure gammafun_DP, &
                     gammafun_QP
  end interface gammafun
  !
  interface log_gammafun
    module procedure log_gammafun_DP, &
                     log_gammafun_QP
  end interface log_gammafun
  !
  interface eval_LagrangePoly
    module procedure eval_LagrangePoly_DP, &
                     eval_LagrangePoly_QP
  end interface eval_LagrangePoly
  !
  interface eval_LagrangePoly1D
    module procedure eval_LagrangePoly1D_DP, &
                     eval_LagrangePoly1D_QP
  end interface eval_LagrangePoly1D
  !
  interface eval_LagrangePoly2D
    module procedure eval_LagrangePoly2D_DP, &
                     eval_LagrangePoly2D_QP
  end interface eval_LagrangePoly2D
  !
  interface eval_D2LagrangeDx2
    module procedure eval_D2LagrangeDx2_DP, &
                     eval_D2LagrangeDx2_QP
  end interface eval_D2LagrangeDx2
  !
  interface eval_DLagrangeDx
    module procedure eval_DLagrangeDx_DP, &
                     eval_DLagrangeDx_QP
  end interface eval_DLagrangeDx
  !
  real(DP), parameter :: zero_DP  = 0.0_DP
  real(DP), parameter :: half_DP  = 0.5_DP
  real(DP), parameter :: one_DP   = 1.0_DP
  real(DP), parameter :: two_DP   = 2.0_DP
  real(DP), parameter :: three_DP = 3.0_DP
  real(DP), parameter :: four_DP  = 4.0_DP
  real(DP), parameter :: six_DP   = 6.0_DP
  real(DP), parameter :: c0_25_DP = 0.25_DP
  real(DP), parameter :: c4_50_DP = 4.5_DP
  real(DP), parameter :: eps14_DP = 1.0e-14_DP
  real(DP), parameter :: eps17_DP = 1.0e-17_DP
  real(DP), parameter :: pi_DP    = 3.14159265358979323846_DP
  !
  real(QP), parameter :: zero_QP  = 0.0_QP
  real(QP), parameter :: half_QP  = 0.5_QP
  real(QP), parameter :: one_QP   = 1.0_QP
  real(QP), parameter :: two_QP   = 2.0_QP
  real(QP), parameter :: three_QP = 3.0_QP
  real(QP), parameter :: four_QP  = 4.0_QP
  real(QP), parameter :: six_QP   = 6.0_QP
  real(QP), parameter :: c0_25_QP = 0.25_QP
  real(QP), parameter :: c4_50_QP = 4.5_QP
  real(QP), parameter :: eps14_QP = 1.0e-14_QP
  real(QP), parameter :: eps17_QP = 1.0e-17_QP
  real(QP), parameter :: pi_QP    = 3.14159265358979323846_QP
  !
contains
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
!  WORKING PRECISION VERSIONS OF THE POLYNOMIAL PROCEDURES
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
pure subroutine legendre_gauss_quadrature(xi,weights)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(inout) :: xi
  real(wp), dimension(:), intent(inout) :: weights
  !
  !.. Local Scalars ..
  integer :: npts
  !
  !.. Local Arrays ..
  real(qp) :: pts_qp(1:size(xi))
  real(qp) :: wts_qp(1:size(xi))
  real(qp) :: val_qp(1:size(xi))
  !
continue
  !
  npts = size(xi)
  !
  call nodes_legendre_gauss(npts,pts_qp,val_qp)
  call weights_legendre_gauss(npts,pts_qp,val_qp,wts_qp)
  !
  xi = chop( pts_qp )
  weights = chop( wts_qp )
  !
end subroutine legendre_gauss_quadrature
!
!###############################################################################
!
pure subroutine legendre_gauss_lobatto_quadrature(xi,weights)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(inout) :: xi
  real(wp), dimension(:), intent(inout) :: weights
  !
  !.. Local Scalars ..
  integer :: npts
  !
  !.. Local Arrays ..
  real(qp) :: pts_qp(1:size(xi))
  real(qp) :: wts_qp(1:size(xi)) 
  real(qp) :: val_qp(1:size(xi))
  !
continue
  !
  npts = size(xi)
  !
  call nodes_legendre_gauss_lobatto(npts-1,pts_qp,val_qp)
  call weights_legendre_gauss_lobatto(npts-1,val_qp,wts_qp)
  !
  xi = chop( pts_qp )
  weights = chop( wts_qp )
  !
end subroutine legendre_gauss_lobatto_quadrature
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
!  DOUBLE PRECISION VERSIONS OF THE POLYNOMIAL PROCEDURES
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
pure subroutine JacobiGL_DP(alpha,beta,n,xi,weights)
  !
  integer,                      intent(in) :: n
  !
  real(DP),                    intent(in) :: alpha
  real(DP),                    intent(in) :: beta
  !
  real(DP), dimension(1:n+1), intent(out) :: xi
  real(DP), dimension(1:n+1), intent(out) :: weights
  !
  real(DP), dimension(1:n+1) :: vn
  !
continue
  !
  xi = zero_DP
  weights = zero_DP
  !
  call nodes_jacobi_gauss_lobatto(n,alpha,beta,xi,vn)
  !
  call weights_dgbook_jacobi_gauss_lobatto(n,alpha,beta,xi,weights)
  !
end subroutine JacobiGL_DP
!
!###############################################################################
!
pure subroutine eval_legendre_poly_DP(n,x,y,dy,d2y)
  !
  integer,    intent(in) :: n
  !
  real(DP),  intent(in) :: x
  real(DP), intent(out) :: y
  real(DP), intent(out) :: dy
  real(DP), intent(out) :: d2y
  !
  real(DP) :: alpha,beta
  !
continue
  !
  alpha = zero_DP
  beta  = zero_DP
  !
  call eval_jacobi_poly(n,alpha,beta,x,y,dy,d2y)
  !
end subroutine eval_legendre_poly_DP
!
!###############################################################################
!
pure function normalized_jacobi_poly_DP(N,alpha,beta,x) result(return_value)
  !
  integer,   intent(in) :: N
  !
  real(DP), intent(in) :: alpha
  real(DP), intent(in) :: beta
  real(DP), intent(in) :: x
  !
  real(DP) :: return_value
  !
  integer   :: i
  !
  real(DP) :: ab,ab1,ab2,ab3,amb,a1,b1
  real(DP) :: gamma0,gamma1
  real(DP) :: aold,anew,bnew,h1,ri
  !
  real(DP) :: PL(N+1)
  !
continue
  !
  ab = alpha + beta
  amb = alpha - beta
  a1 = alpha + one_DP
  b1 = beta + one_DP
  ab1 = ab + one_DP
  ab2 = ab + two_DP
  ab3 = ab + three_DP
  !
  gamma0 = two_DP**ab1 / ab1 * gammafun(a1) * gammafun(b1) / gammafun(ab1)
  !
  PL(1) = one_DP/sqrt(gamma0)
  !
  if (N == 0) then
    !
    return_value = PL(1)
    !
  else
    !
    gamma1 = a1 * b1 / ab3 * gamma0
    !
    PL(2) = half_DP*( ab2*x + amb ) / sqrt(gamma1)
    !
    aold = two_DP / ab2 * sqrt( a1*b1/ab3 )
    !
    do i = 1,N-1
      ri = real(i,kind=DP)
      h1 = two_DP*ri + ab
      anew = two_DP / (h1+two_DP) * sqrt( (ri+one_DP)*(ri+ab1)*(ri+a1)*(ri+b1)/(h1+one_DP)/(h1+three_DP) )
      bnew = (beta*beta - alpha*alpha) / h1 / (h1+two_DP)
      PL(i+2) = one_DP/anew*( -aold*PL(i) + (x-bnew)*PL(i+1) )
      aold = anew
    end do
    !
    return_value = PL(N+1)
    !
  end if
  !
end function normalized_jacobi_poly_DP
!
!###############################################################################
!
pure function grad_normalized_jacobi_poly_DP(n,alpha,beta,x) result(return_value)
  !
  integer,  parameter :: lp = dp
  !
  integer,   intent(in) :: n
  real(lp), intent(in) :: alpha
  real(lp), intent(in) :: beta
  real(lp), intent(in) :: x
  !
  real(lp) :: return_value
  !
  !.. Local Scalars ..
  real(lp) :: rn,a1,b1,ab1
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !
  if (n == 0) then
    return_value = zero
  else
    !
    rn  = real(n,kind=lp)
    a1  = alpha + one
    b1  = beta  + one
    ab1 = alpha + b1
    !
    return_value = sqrt(rn*(rn+ab1)) * normalized_jacobi_poly(n-1,a1,b1,x)
    !
  end if
  !
end function grad_normalized_jacobi_poly_DP
!
!###############################################################################
!
pure subroutine eval_normalized_jacobi_poly_DP(n,alpha,beta,x,y,dy)
  !
  integer,    intent(in) :: n
  !
  real(DP),  intent(in) :: alpha
  real(DP),  intent(in) :: beta
  real(DP),  intent(in) :: x
  real(DP), intent(out) :: y
  real(DP), intent(out) :: dy
  !
  integer :: n1
  !
  real(DP) :: a1,b1,ab1,rn
  !
continue
  !
  a1  = alpha + one_DP
  b1  = beta  + one_DP
  ab1 = alpha + beta + one_DP
  rn  = real(n,kind=DP)
  n1  = n - 1
  !
  y = normalized_jacobi_poly(n,alpha,beta,x)
  !
  dy = sqrt(rn*(rn+ab1))*normalized_jacobi_poly(n1,a1,b1,x)
  !
end subroutine eval_normalized_jacobi_poly_DP
!
!###############################################################################
!
pure subroutine eval_jacobi_poly_DP(n,alpha,beta,x,y,dy,d2y)
  !
  integer,    intent(in) :: n
  !
  real(DP),  intent(in) :: alpha
  real(DP),  intent(in) :: beta
  real(DP),  intent(in) :: x
  real(DP), intent(out) :: y
  real(DP), intent(out) :: dy
  real(DP), intent(out) :: d2y
  !
  integer :: i
  !
  real(DP) :: ab,di,c0,c1,c2,c3,c4
  real(DP) :: ym,yp,dym,dyp,d2ym,d2yp
  !
continue
  !
  y   = one_DP
  dy  = zero_DP
  d2y = zero_DP
  !
  if (n > 0) then
    !
    ab  = alpha + beta
    !
    y   = half_DP*(ab + two_DP)*x + half_DP*(alpha - beta)
    dy  = half_DP*(ab + two_DP)
    d2y = zero_DP
    !
  end if
  !
  if (n > 1) then
    !
    yp   = one_DP
    dyp  = zero_DP
    d2yp = zero_DP
    !
    do i = 2,n
      !
      di = real(i,kind=DP)
      c0 = two_DP*di + ab
      c1 = two_DP * di * (di+ab) * (c0-two_DP)
      c2 = (c0-one_DP) * (c0-two_DP) * c0
      c3 = (c0-one_DP) * (alpha-beta) * ab
      c4 = two_DP * c0 * (di+alpha-one_DP) * (di+beta-one_DP)
      !
      ym = y
      y  = ((c2*x + c3)*y - c4*yp)/c1
      yp = ym
      !
      dym = dy
      dy  = ((c2*x + c3)*dy - c4*dyp + c2*yp)/c1
      dyp = dym
      !
      d2ym = d2y
      d2y  = ((c2*x + c3)*d2y - c4*d2yp + two_DP*c2*dyp)/c1
      d2yp = d2ym
      !
    end do
    !
  end if
  !
end subroutine eval_jacobi_poly_DP
!
!###############################################################################
!
pure subroutine eval_chebyshev_poly_DP(n,x,y,dy,d2y)
  !
  integer,    intent(in) :: n
  !
  real(DP),  intent(in) :: x
  real(DP), intent(out) :: y
  real(DP), intent(out) :: dy
  real(DP), intent(out) :: d2y
  !
  integer :: k
  !
  real(DP) :: yp,ym,dyp,dym,d2yp,d2ym
  !
continue
  !
  y   = one_DP
  dy  = zero_DP
  d2y = zero_DP
  !
  if (n > 0) then
    y   = x
    dy  = one_DP
    d2y = zero_DP
  end if
  !
  if (n > 1) then
    !
    yp   = one_DP
    dyp  = zero_DP
    d2yp = zero_DP
    !
    do k = 2,n
      !
      ym = y
      y  = two_DP*x*y - yp
      yp = ym
      !
      dym = dy
      dy  = two_DP*x*dy + two_DP*yp - dyp
      dyp = dym
      !
      d2ym = d2y
      d2y  = two_DP*x*d2y + four_DP*dyp - d2yp
      d2yp = d2ym
      !
    end do
    !
  end if
  !
end subroutine eval_chebyshev_poly_DP
!
!###############################################################################
!
pure subroutine nodes_jacobi_gauss_DP(n,alpha,beta,xi,dxi)
  !
  integer,                    intent(in) :: n
  !
  real(DP),                  intent(in) :: alpha
  real(DP),                  intent(in) :: beta
  real(DP), dimension(1:n), intent(out) :: xi
  real(DP), dimension(1:n), intent(out) :: dxi
  !
  integer :: i,it
  !
  real(DP) :: const,ri,x,y,dy,d2y
  !
continue
  !
  if (n > 0) then
    xi(1)  = (beta - alpha) / (alpha + beta + two_DP)
    dxi(1) = half_DP*(alpha + beta) + one_DP
  end if
  !
  if (n > 1) then
    !
    const = half_DP*pi_DP / (two_DP*real(n,kind=DP) + alpha + beta + one_DP)
    !
    do i = 1,n
      !
      ri = real(i,kind=DP)
      x  = -cos( const * (four_DP*ri + alpha + beta - one_DP) )
      !
      do it = 1,8
        call eval_jacobi_poly(n,alpha,beta,x,y,dy,d2y)
        if (abs(y) <= eps17_DP) exit
        x = x - y/dy
      end do
      !
      if (abs(x) <= eps17_DP) x = zero_DP
      !
      xi(i)  = x
      dxi(i) = dy
      !
    end do
    !
    do i = 1,n/2
      xi (n-i+1) = -xi (i)
     !dxi(n-i+1) = -dxi(i)
    end do
    !
  end if
  !
end subroutine nodes_jacobi_gauss_DP
!
!###############################################################################
!
pure subroutine nodes_legendre_gauss_DP(n,xi,dxi)
  !
  integer,                    intent(in) :: n
  !
  real(DP), dimension(1:n), intent(out) :: xi
  real(DP), dimension(1:n), intent(out) :: dxi
  !
  real(DP) :: alpha,beta
  !
continue
  !
  call nodes_jacobi_gauss(n,zero_DP,zero_DP,xi,dxi)
  !
end subroutine nodes_legendre_gauss_DP
!
!###############################################################################
!
pure subroutine nodes_right_radau_DP(n,xi)
  !
  integer,                        intent(in) :: n
  !
  real(DP), dimension(1:n+1), intent(inout) :: xi
  !
  integer :: i,it
  real(DP) :: x,y,dy
  !
  real(DP), dimension(0:n+1) :: cfrr
  real(DP), dimension(0:n) :: cdrr
  !
continue
  !
  ! Set the final zero_DP which will always be unity
  !
  xi(n+1) = one_DP
  !
  if (n > 0) then
    !
    ! Get the coefficients of the right Radau polynomial
    !
    cfrr(:) = coefrr_DP(n+1)
    !
    do i = 0,n
      cdrr(i) = real(i+1,kind=DP) * cfrr(i+1)
    end do
    !
    xi(n+1) = one_DP
    !
    do i = 1,n
      !
      x  = -cos( pi_DP * real(i,kind=DP) / real(n+1,kind=DP) )
      !
      do it = 1,8
        !
        y = sum( (x**intseq(0,n+1)) * cfrr(:) )
        dy = sum( (x**intseq(0,n)) * cdrr(:) )
        !
        x = x - y/dy
        !
      end do
      !
      xi(i) = x
      !
    end do
    !
  end if
  !
end subroutine nodes_right_radau_DP
!
!###############################################################################
!
pure function coefrr_DP(n) result(return_value)
  !
  integer, intent(in) :: n
  !
  real(DP), dimension(0:n) :: return_value
  !
  integer :: k,ip
  !
  real(DP), dimension(0:n)   :: cfk
  real(DP), dimension(0:n)   :: cfkmt
  real(DP), dimension(0:n-1) :: cfkm1
  !
continue
  !
  cfk = coefle_DP(n)
  cfkm1 = coefle_DP(n-1)
  !
  do k = 0,n-1
    cfkmt(k) = cfkm1(k)
  end do
  !
  cfkmt(n) = zero_DP
  ip = (-1)**n
  !
  do k = 0,n
    return_value(k) = half_DP * real(ip,kind=DP) * (cfk(k) - cfkmt(k))
  end do
  !
end function coefrr_DP
!
!###############################################################################
!
pure function coefle_DP(n) result(return_value)
  !
  integer, intent(in) :: n
  !
  real(DP), dimension(0:n) :: return_value
  !
  integer :: i,j,k,m
  !
  real(DP), dimension(-1:n,-1:n) :: cle
  !
continue
  !
  cle(:,:) = zero_DP
  !
  cle(0,0) = one_DP
  cle(1,1) = one_DP
  !
  do k = 2,n
    do m = 0,k
      cle(k,m) = (real(2*k-1,kind=DP)/real(k,kind=DP))*cle(k-1,m-1) - &
                 (real(  k-1,kind=DP)/real(k,kind=DP))*cle(k-2,m  )
    end do
  end do
  !
  do j = 0,n
    return_value(j) = cle(n,j)
  end do
  !
end function coefle_DP
!
!###############################################################################
!
pure subroutine nodes_chebyshev_gauss_DP(n,xi)
  !
  integer,                    intent(in) :: n
  !
  real(DP), dimension(1:n), intent(out) :: xi
  !
  real(DP) :: alpha,beta
  real(DP), dimension(1:n) :: dxi
  !
continue
  !
  alpha = -half_DP
  beta  = -half_DP
  !
  call nodes_jacobi_gauss(n,alpha,beta,xi,dxi)
  !
end subroutine nodes_chebyshev_gauss_DP
!
!###############################################################################
!
pure subroutine nodes_jacobi_gauss_lobatto_DP(n,alpha,beta,xi,vn)
  !
  integer,                     intent(in) :: n
  !
  real(DP),                   intent(in) :: alpha
  real(DP),                   intent(in) :: beta
  real(DP), dimension(1:n+1), intent(out) :: xi
  real(DP), dimension(1:n+1), intent(out) :: vn
  !
  integer :: i,it,n1
  !
  real(DP) :: x,y,dy,d2y,ri,rn,ab,a1,b1,const
  !
continue
  !
  if (n == 0) then
    !
    x = zero_DP
    call eval_jacobi_poly(n,alpha,beta,x,y,dy,d2y)
    xi(1) = x
    vn(1) = y
    !
  else
    !
    x = -one_DP
    call eval_jacobi_poly(n,alpha,beta,x,y,dy,d2y)
    xi(1) = x
    vn(1) = y
    !
    x = one_DP
    call eval_jacobi_poly(n,alpha,beta,x,y,dy,d2y)
    xi(n+1) = x
    vn(n+1) = y
    !
  end if
  !
  if (n > 1) then
    !
    rn = real(n,kind=DP)
    ab = alpha + beta
    a1 = alpha + one_DP
    b1 = beta  + one_DP
    n1 = n-1
    const = half_DP * pi_DP / (two_DP*rn + ab + one_DP)
    !
    do i = 2,n
      !
      ri = real(i-1,kind=DP)
      x  = -cos( const * (four_DP*ri + ab + one_DP) )
      !
      do it = 1,8
        call eval_jacobi_poly(n1,a1,b1,x,y,dy,d2y)
        if (abs(y) <= eps17_DP) exit
        x = x - y/dy
      end do
      !
      if (abs(x) <= eps17_DP) x = zero_DP
      !
      xi(i) = x
      vn(i) = -half_DP * dy * (one_DP - x*x) / rn
      !
    end do
    !
  end if
  !
end subroutine nodes_jacobi_gauss_lobatto_DP
!
!###############################################################################
!
pure subroutine nodes_legendre_gauss_lobatto_DP(n,xi,vn)
  !
  integer,                      intent(in) :: n
  real(DP), dimension(1:n+1), intent(out) :: xi
  real(DP), dimension(1:n+1), intent(out) :: vn
  !
  integer :: i,it,n2
  !
  real(DP) :: sn,x,y,dy,d2y
  !
continue
  !
  if (n == 0) then
    xi(1) = zero_DP
    vn(1) = one_DP
  else
    !
    n2 = (n-1)/2
    sn = real(2*n - 4*n2 - 3,kind=DP)
    !
    xi(1) = -one_DP
    vn(1) =  sn
    !
    xi(n+1) = one_DP
    vn(n+1) = one_DP
    !
  end if
  !
  if (n > 1) then
    x = zero_DP
    call eval_legendre_poly(n,x,y,dy,d2y)
    xi(n2+2) = x
    vn(n2+2) = y
  end if
  !
  if (n > 2) then
    !
    do i = 1,n2
      !
      x = cos( pi_DP * real(i,kind=DP) / real(n,kind=DP) )
      !
      do it = 1,8
        call eval_legendre_poly(n,x,y,dy,d2y)
        x = x - dy/d2y
      end do
      !
      xi(i+1)   = -x
      xi(n-i+1) =  x
      !
      vn(i+1)   = y*sn
      vn(n-i+1) = y
      !
    end do
    !
  end if
  !
end subroutine nodes_legendre_gauss_lobatto_DP
!
!###############################################################################
!
pure subroutine nodes_chebyshev_gauss_lobatto_DP(n,xi)
  !
  integer,                      intent(in) :: n
  !
  real(DP), dimension(1:n+1), intent(out) :: xi
  !
  integer :: i,n2
  !
continue
  !
  if (n == 0) then
    xi(1) = zero_DP
  else
    xi(1)   = -one_DP
    xi(n+1) =  one_DP
  end if
  !
  if (n > 1) then
    n2 = (n-1)/2
    xi(n2+2) = zero_DP
  end if
  !
  if (n > 2) then
    do i = 1,n2
      xi(i+1)   = -cos( pi_DP * real(i,kind=DP) / real(n,kind=DP) )
      xi(n-i+1) =  cos( pi_DP * real(i,kind=DP) / real(n,kind=DP) )
    end do
  end if
  !
end subroutine nodes_chebyshev_gauss_lobatto_DP
!
!###############################################################################
!
elemental function gammafun_DP(x) result(return_value)
  !
  real(DP), intent(in) :: x
  !
  real(DP) :: return_value
  !
  integer   :: ind,k
  real(DP) :: xx,pr,s,g
  !
  real(DP), parameter :: s_start = 0.426401432711220868_DP
  real(DP), parameter :: gamma_constant(1:11) = &
             [ -0.524741987629368444e+00_DP,  0.116154405493589130e+00_DP, &
               -0.765978624506602380e-02_DP,  0.899719449391378898e-04_DP, &
               -0.194536980009534621e-07_DP,  0.199382839513630987e-10_DP, &
               -0.204209590209541319e-11_DP,  0.863896817907000175e-13_DP, &
                0.152237501608472336e-13_DP, -0.825725175277719950e-14_DP, &
                0.299734782205224610e-14_DP ]
continue
  !
  xx = x
  return_value = one_DP
  !
  do
    if (abs(xx-one_DP) < eps14_DP) return
    if (xx < one_DP) exit
    xx = xx - one_DP
    return_value = return_value*xx
  end do
  !
  ind = 0
  !
  if (xx < half_DP) then
    ind = 1
    return_value = return_value*pi_DP/sin(pi_DP*xx)
    xx = one_DP - xx
  end if
  !
  pr = one_DP
  s = s_start
  !
  do k = 1,11
    pr = pr * (xx-real(k,kind=DP)) / (xx+real(k-1,kind=DP))
    s  = s + gamma_constant(k)*pr
  end do
  !
  g = s * exp(one_DP-xx) * (xx+c4_50_DP)**(xx-half_DP)
  !
  if (ind == 1) then
    return_value = return_value/g
  else
    return_value = return_value*g
  end if
  !
end function gammafun_DP
!
!###############################################################################
!
elemental function log_gammafun_DP(x) result(return_value)
  real(DP), intent(in) :: x
  real(DP) :: return_value
  return_value = log( gammafun(x) )
end function log_gammafun_DP
!
!###############################################################################
!
pure subroutine weights_jacobi_gauss_DP(n,alpha,beta,xi,dy,weights)
  !***************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE JACOBI GAUSS FORMULA
  !   N  = ORDER OF THE FORMULA
  !   alpha  = PARAMETER > -1
  !   beta  = PARAMETER > -1
  !   CS = ZEROES OF THE JACOBI POLYNOMIAL, CS(I), I=1,N
  !   DZ = VECTOR OF THE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
  !   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
  !***************************************************************
  integer,                    intent(in) :: n
  !
  real(DP),                  intent(in) :: alpha
  real(DP),                  intent(in) :: beta
  !
  real(DP), dimension(1:n),  intent(in) :: xi
  real(DP), dimension(1:n),  intent(in) :: dy
  real(DP), dimension(1:n), intent(out) :: weights
  !
  integer :: i,m
  !
  real(DP) :: ab,a2,b2,c,ri
  !
continue
  !
  weights(:) = one_DP
  !
  ab = alpha + beta + two_DP
  a2 = alpha + two_DP
  b2 = beta  + two_DP
  !
  c = half_DP * (two_DP**ab) * gammafun(a2) * gammafun(b2) / gammafun(ab)
  !
  do i = 2,n
    ri = real(i,kind=DP)
    c  = c * (ri+alpha) * (ri+beta) / (ri*(ri+alpha+beta))
  end do
  !
  do i = 1,n
    weights(i) = c / ((one_DP - xi(i)*xi(i)) * dy(i) * dy(i))
  end do
  !
end subroutine weights_jacobi_gauss_DP
!
!################################################################################
!
pure subroutine weights_legendre_gauss_DP(n,xi,dy,weights)
  !****************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS FORMULA
  !   N  = ORDER OF THE FORMULA
  !   CS = ZEROES OF THE LEGENDRE POLYNOMIAL, CS(I), I=1,N
  !   DZ = VECTOR OF THE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
  !   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
  !****************************************************************
  integer,                    intent(in) :: n
  !
  real(DP), dimension(1:n),  intent(in) :: xi
  real(DP), dimension(1:n),  intent(in) :: dy
  real(DP), dimension(1:n), intent(out) :: weights
  !
  integer :: i
  !
continue
  !
  weights(:) = two_DP
  !
  do i = 1,n
    weights(i) = two_DP / ((one_DP - xi(i)*xi(i)) * dy(i) * dy(i))
  end do
  !
end subroutine weights_legendre_gauss_DP
!
!################################################################################
!
pure subroutine weights_chebyshev_gauss_DP(n,weights)
  !****************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE CHEBYSHEV GAUSS FORMULA
  !   N  = ORDER OF THE FORMULA
  !   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
  !****************************************************************
  integer,                    intent(in) :: n
  !
  real(DP), dimension(1:n), intent(out) :: weights
  !
  integer :: i
  !
continue
  !
  weights(:) = pi_DP
  !
  do i = 1,n
    weights(i) = pi_DP/real(n,kind=DP)
  end do
  !
end subroutine weights_chebyshev_gauss_DP
!
!################################################################################
!
pure subroutine weights_jacobi_gauss_lobatto_DP(n,alpha,beta,xi,weights)
  !********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE JACOBI GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   alpha  = PARAMETER > -1
  !   beta  = PARAMETER > -1
  !   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !********************************************************************
  integer,                      intent(in) :: n
  !
  real(DP),                    intent(in) :: alpha
  real(DP),                    intent(in) :: beta
  !
  real(DP), dimension(1:n+1),  intent(in) :: xi
  real(DP), dimension(1:n+1), intent(out) :: weights
  !
  integer :: i
  !
  real(DP) :: a1,a2,b1,b2,ab,ab1,ab2
  real(DP) :: c,c1,c2,c3,ri,rn,su
  real(DP) :: y,z,dy,d2y
  !
continue
  !
  a1  = alpha + one_DP
  a2  = alpha + two_DP
  b1  = beta  + one_DP
  b2  = beta  + two_DP
  ab  = alpha + beta
  ab1 = alpha + beta + one_DP
  ab2 = alpha + beta + two_DP
  !
  c = (two_DP**ab) * gammafun(a1) * gammafun(b1) / gammafun(ab2)
  !
  weights(:) = two_DP*c
  !
  if (n == 1) then
    weights(1)   = two_DP * c * a1 / ab2
    weights(n+1) = two_DP * c * b1 / ab2
  end if
  !
  if (n > 1) then
    !
    rn = real(n,kind=DP)
    c  = c * (two_DP*rn + ab) / (rn + ab1)
    c1 = c * a1 / (b2 * ab2)
    c2 = c * b1 / (a2 * ab2)
    c3 = half_DP * c * a1 * b1
    !
    do i = 2,n-1
      ri = real(i-1,kind=DP)
      c1 = c1*(ri+a1)*ri / ((ri+ab2)*(ri+b2))
      c2 = c2*(ri+b1)*ri / ((ri+ab2)*(ri+a2))
      c3 = c3*(ri+a1)*(ri+b1) / ((ri+two_DP)*(ri+ab1))
    end do
    !
    su = zero_DP
    do i = 2,n
      su = su + xi(i)
    end do
    !
    weights(1)   = c1*(rn-one_DP-su)
    weights(n+1) = c2*(rn-one_DP+su)
    !
    do i = 2,n
      call eval_jacobi_poly(n  ,alpha,beta,xi(i),y,dy,d2y)
      call eval_jacobi_poly(n-1,alpha,beta,xi(i),z,dy,d2y)
      weights(i) = -c3/y/dy
    end do
    !
  end if
  !
end subroutine weights_jacobi_gauss_lobatto_DP
!
!################################################################################
!
pure subroutine weights_legendre_gauss_lobatto_DP(n,vn,weights)
  !**********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N
  !   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !**********************************************************************
  integer,                      intent(in) :: n
  !
  real(DP), dimension(1:n+1),  intent(in) :: vn
  real(DP), dimension(1:n+1), intent(out) :: weights
  !
  integer :: i
  !
  real(DP) :: c
  !
continue
  !
  c = two_DP
  if (n > 0) then
    c = c / (real(n,kind=DP)*(real(n,kind=DP)+one_DP))
  end if
  !
  do i = 1,n+1
    weights(i) = c / (vn(i)*vn(i))
  end do
  !
end subroutine weights_legendre_gauss_lobatto_DP
!
!################################################################################
!
pure subroutine weights_right_radau_DP(n,xi,weights)
  !**********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   XI = RIGHT RADAU NODES, XI(I), I=0,N
  !   VN = VALUES OF THE RIGHT RADAU POLYNOMIAL AT THE NODES, VN(I), I=0,N
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !**********************************************************************
  integer,                      intent(in) :: n
  !
  real(DP), dimension(1:n+1),  intent(in) :: xi
  real(DP), dimension(1:n+1), intent(out) :: weights
  !
  integer :: i
  !
  real(DP) :: c,y,dy,d2y
  !
continue
  !
  c = two_DP
  !
  if (n > 0) then
    !
    c = one_DP / (real(n+1,kind=DP)*real(n+1,kind=DP))
    !
    do i = 1,n+1
      call eval_legendre_poly(n,xi(i),y,dy,d2y)
      weights(i) = c * (one_DP+xi(i)) / (y*y)
    end do
    !
  end if
  !
end subroutine weights_right_radau_DP
!
!################################################################################
!
pure subroutine weights_chebyshev_gauss_lobatto_DP(n,weights)
  !***********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE CHEBYSHEV GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !***********************************************************************
  integer,                      intent(in) :: n
  !
  real(DP), dimension(1:n+1), intent(out) :: weights
  !
  integer :: i
  !
continue
  !
  weights(:) = pi_DP
  !
  if (n > 0) then
    weights(1)   = half_DP*pi_DP/real(n,kind=DP)
    weights(n+1) = half_DP*pi_DP/real(n,kind=DP)
  end if
  !
  do i = 2,n
    weights(i) = pi_DP/real(n,kind=DP)
  end do
  !
end subroutine weights_chebyshev_gauss_lobatto_DP
!
!###############################################################################
!
pure subroutine weights_dgbook_jacobi_gauss_lobatto_DP(n,alpha,beta,xi,weights)
  !********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE JACOBI GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   alpha  = PARAMETER > -1
  !   beta  = PARAMETER > -1
  !   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !********************************************************************
  integer,                      intent(in) :: n
  !
  real(DP),                    intent(in) :: alpha
  real(DP),                    intent(in) :: beta
  !
  real(DP), dimension(1:n+1),  intent(in) :: xi
  real(DP), dimension(1:n+1), intent(out) :: weights
  !
  integer :: i
  !
  real(DP) :: a1,a2,b1,b2,ab,ab1,ab2
  real(DP) :: c,c1,c2,c3,ri,rn,su
  real(DP) :: y,z,dy,d2y,scl
  !
continue
  !
  weights(:) = zero_DP
  !
  if (n > 1) then
    !
    a1  = alpha + one_DP
    a2  = alpha + two_DP
    b1  = beta  + one_DP
    b2  = beta  + two_DP
    ab  = alpha + beta
    ab1 = alpha + beta + one_DP
    ab2 = alpha + beta + two_DP
    !
    c = (two_DP**ab) * gammafun(a1) * gammafun(b1) / gammafun(ab2)
    !
    rn = real(n,kind=DP)
    c  = c * (two_DP*rn + ab) / (rn + ab1)
    c1 = c * a1 / (b2 * ab2)
    c2 = c * b1 / (a2 * ab2)
    c3 = half_DP * c * a1 * b1
    !
    do i = 2,n-1
      ri = real(i-1,kind=DP)
      c1 = c1*(ri+a1)*ri / ((ri+ab2)*(ri+b2))
      c2 = c2*(ri+b1)*ri / ((ri+ab2)*(ri+a2))
      c3 = c3*(ri+a1)*(ri+b1) / ((ri+two_DP)*(ri+ab1))
    end do
    !
    su = zero_DP
    do i = 2,n
      su = su + xi(i)
    end do
    !
    scl = ( (one_DP-xi(1))**a1 ) * ( (one_DP+xi(1))**b1 )
    weights(1)   = scl*c1*(rn-one_DP-su)
    !
    scl = ( (one_DP-xi(n+1))**a1 ) * ( (one_DP+xi(n+1))**b1 )
    weights(n+1) = scl*c2*(rn-one_DP+su)
    !
    do i = 2,n
      call eval_jacobi_poly(n  ,alpha,beta,xi(i),y,dy,d2y)
      call eval_jacobi_poly(n-1,alpha,beta,xi(i),z,dy,d2y)
      scl = ( (one_DP-xi(i))**a1 ) * ( (one_DP+xi(i))**b1 )
      weights(i) = -scl*c3/y/dy
    end do
    !
  end if
  !
end subroutine weights_dgbook_jacobi_gauss_lobatto_DP
!
!################################################################################
!
pure subroutine old_derivative_matrix_legendre_gauss_DP(n,cs,dz,dma)
  !
  !*************************************************************
  !*************************************************************
  !********************                       ******************
  !********************     NOT IN MANUAL     ******************
  !********************                       ******************
  !*************************************************************
  !*************************************************************
  !
  !***********************************************************************
  !  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE
  !  LEGENDRE GAUSS NODES
  !  N   = NUMBER OF LEGENDRE GAUSS NODES
  !  CS  = VECTOR OF THE ZEROES, CS(I), I=1,N
  !  DZ  = LEGENDRE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
  !  DMA = DERIVATIVE MATRIX, DMA(I,J), I=1,N  J=1,N
  !***********************************************************************
  !
  !.. Formal Arguments ..
  integer,                        intent(in) :: n
  real(DP), dimension(1:n),      intent(in) :: cs
  real(DP), dimension(1:n),      intent(in) :: dz
  real(DP), dimension(1:n,1:n), intent(out) :: dma
  !
  !.. Local Scalars ..
  integer   :: i,j
  real(DP) :: vi,zi,vj,zj
  !
continue
  !
  dma(1,1) = zero_DP
  !
  if (n <= 1) return
  !
  do i = 1,n
    vi = dz(i)
    zi = cs(i)
    do j = 1,n
      if (i /= j) then
        vj = dz(j)
        zj = cs(j)
        dma(i,j) = vi/(vj*(zi-zj))
      else
        dma(i,i) = zi/(one_DP-zi*zi)
      end if
    end do
  end do
  !
end subroutine old_derivative_matrix_legendre_gauss_DP
!
!################################################################################
!
pure subroutine derivative_matrix_jacobi_gauss_lobatto_DP(n,alpha,beta,et,vn,dma)
  !
  !***********************************************************************
  !  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE
  !  JACOBI GAUSS-LOBATTO NODES
  !  N   = PARAMETER RELATIVE TO THE DIMENSION OF THE MATRIX
  !  ALPHA   = PARAMETER >-1
  !  BETA   = PARAMETER >-1
  !  ET  = VECTOR OF THE NODES, ET(I), I=0,N
  !  VN  = VALUES OF THE JACOBI POLYNOMIAL AT THE NODES, VN(I), I=0,N
  !  DMA = DERIVATIVE MATRIX, DMA(I,J), I=0,N  J=0,N
  !***********************************************************************
  !
  !.. Formal Arguments ..
  integer,                        intent(in) :: n
  real(DP),                      intent(in) :: alpha
  real(DP),                      intent(in) :: beta
  real(DP), dimension(1:n),      intent(in) :: et
  real(DP), dimension(1:n),      intent(in) :: vn
  real(DP), dimension(1:n,1:n), intent(out) :: dma
  !
  !.. Local Scalars ..
  integer   :: i,j
  real(DP) :: a1,b1,ab,dn,ei,ej
  real(DP) :: c1,c2,c3,c4,vi,vj
  !
continue
  !
  dma(1,1) = zero_DP
  if (n == 1) return
  !
  a1 = alpha+one_DP
  b1 = beta+one_DP
  ab = alpha+beta
  dn = real(n-1,kind=DP)
  c1 = dn*(dn+ab+one_DP)
  c2 = a1*vn(1)/(b1*vn(n))
  dma(1,1) = half_DP*(alpha-c1)/(beta+two_DP)
  dma(n,n) = half_DP*(c1-beta)/(alpha+two_DP)
  dma(1,n) = -half_DP*c2
  dma(n,1) = half_DP/c2
  if (n == 2) return
  !
  c3 = vn(1)/b1
  c4 = vn(n)/a1
  do j = 2,n-1
    vj = vn(j)
    ej = et(j)
    dma(1,j) = -c3/(vj*(one_DP+ej))
    dma(n,j) = c4/(vj*(one_DP-ej))
    dma(j,1) = vj/(c3*(one_DP+ej))
    dma(j,n) = -vj/(c4*(one_DP-ej))
  end do
  !
  do i = 2,n-1
    vi = vn(i)
    ei = et(i)
    do j = 2,n-1
      if (i /= j) then
        vj = vn(j)
        ej = et(j)
        dma(i,j) = vi/(vj*(ei-ej))
      else
        dma(i,i) = half_DP*(ab*ei+alpha-beta)/(one_DP-ei*ei)
      end if
    end do
  end do
  !
end subroutine derivative_matrix_jacobi_gauss_lobatto_DP
!
!################################################################################
!
pure subroutine old_derivative_matrix_legendre_gauss_lobatto_DP(n,et,vn,dma)
  !
  !***********************************************************************
  !  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE
  !  LEGENDRE GAUSS-LOBATTO NODES
  !  N   = NUMBER OF LEGENDRE GAUSS-LOBATTO NODES
  !  ET  = VECTOR OF THE NODES, ET(I), I=1,N
  !  VN  = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=1,N
  !  DMA = DERIVATIVE MATRIX, DMA(I,J), I=1,N  J=1,N
  !***********************************************************************
  !
  !.. Formal Arguments ..
  integer,                        intent(in) :: n
  real(DP), dimension(1:n),      intent(in) :: et
  real(DP), dimension(1:n),      intent(in) :: vn
  real(DP), dimension(1:n,1:n), intent(out) :: dma
  !
  !.. Local Scalars ..
  integer   :: i,j
  real(DP) :: ei,ej,vi,vj,dn,c
  !
continue
  !
  do i = 1,n
    vi = vn(i)
    ei = et(i)
    do j = 1,n
      if (i /= j) then
        vj = vn(j)
        ej = et(j)
        dma(i,j) = vi/(vj*(ei-ej))
      else
        dma(i,i) = zero_DP
      end if
    end do
  end do
  !
  if (n > 1) then
    dn = real(n-1,kind=DP)
    c  = c0_25_DP*dn*(dn+one_DP)
    dma(1,1) = -c
    dma(n,n) = c
  end if
  !
end subroutine old_derivative_matrix_legendre_gauss_lobatto_DP
!
!################################################################################
!
pure subroutine derivative_matrix_chebyshev_gauss_lobatto_DP(n,et,dma)
  !
  !***********************************************************************
  !  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE
  !  CHEBYSHEV GAUSS-LOBATTO NODES
  !  N   = PARAMETER RELATIVE TO THE DIMENSION OF THE MATRIX
  !  ET  = VECTOR OF THE NODES, ET(I), I=0,N
  !  DMA = DERIVATIVE MATRIX, DMA(I,J), I=0,N  J=0,N
  !***********************************************************************
  !
  !.. Formal Arguments ..
  integer,                        intent(in) :: n
  real(DP), dimension(1:n),      intent(in) :: et
  real(DP), dimension(1:n,1:n), intent(out) :: dma
  !
  !.. Local Scalars ..
  integer   :: i,j
  real(DP) :: dn,cn,ej,ei
  real(DP) :: sn,sgni,sgnj
  !
continue
  !
  dma(1,1) = zero_DP
  if (n == 1) return
  !
  dn = real((n-1),kind=DP)
  cn = (two_DP*dn*dn+one_DP)/six_DP
  sn = real(1+4*((n-1)/2)-2*(n-1),kind=DP)
  dma(1,1) = -cn
  dma(n,n) = cn
  dma(1,n) = -half_DP*sn
  dma(n,1) = half_DP*sn
  if (n == 2) return
  !
  sgnj = -one_DP
  do j = 2,n-1
    ej = et(j)
    dma(1,j) = -two_DP*sgnj/(one_DP+ej)
    dma(n,j) = two_DP*sgnj*sn/(one_DP-ej)
    dma(j,1) = half_DP*sgnj/(one_DP+ej)
    dma(j,n) = -half_DP*sgnj*sn/(one_DP-ej)
    sgnj = -sgnj
  end do
  !
  sgni = -one_DP
  do i = 2,n-1
    ei = et(i)
    sgnj = -one_DP
    do j = 2,n-1
      if (i /= j) then
        ej = et(j)
        dma(i,j) = sgni*sgnj/(ei-ej)
      else
        dma(i,i) = -half_DP*ei/(one_DP-ei*ei)
      end if
      sgnj = -sgnj
    end do
    sgni = -sgni
  end do
  !
end subroutine derivative_matrix_chebyshev_gauss_lobatto_DP
!
!###############################################################################
!
pure function derivative_matrix_legendre_gauss_DP(n,a) result(return_value)
  !
  !*************************************************************
  !*************************************************************
  !********************                       ******************
  !********************     NOT IN MANUAL     ******************
  !********************                       ******************
  !*************************************************************
  !*************************************************************
  !
  !***********************************************************************
  !  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE
  !  LEGENDRE GAUSS NODES
  !***********************************************************************
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: n ! ORDER OF LEGENDRE GAUSS POLYNOMIAL
  real(lp), intent(in) :: a ! REAL VALUE TO DISTINGUISH BETWEEN DP AND QP VERSIONS
  !
  !.. Function Result ..
  real(lp), dimension(1:n+1,1:n+1) :: return_value ! DERIVATIVE MATRIX
  !
  !.. Local Scalars ..
  integer  :: i,j,np
  real(lp) :: vi,zi,vj,zj
  !
  !.. Local Arrays ..
  real(lp), dimension(1:n+1) :: cs ! VECTOR OF THE ZEROES
  real(lp), dimension(1:n+1) :: dz ! LEGENDRE DERIVATIVES AT THE ZEROES
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  np = n + 1
  !
  return_value(1,1) = zero
  !
  if (n <= 1) return
  !
  ! Compute the zeros of the Legendre Gauss polynomial of
  ! order 'n' and the value of the derivative at each zero
  !
  call nodes_legendre_gauss(np,cs,dz)
  !
  do i = 1,np
    vi = dz(i)
    zi = cs(i)
    do j = 1,np
      if (i /= j) then
        vj = dz(j)
        zj = cs(j)
        return_value(i,j) = vi/(vj*(zi-zj))
      else
        return_value(i,i) = zi/(one-zi*zi)
      end if
    end do
  end do
  !
end function derivative_matrix_legendre_gauss_DP
!
!################################################################################
!
pure function derivative_matrix_legendre_gauss_lobatto_DP(n,a) result(return_value)
  !
  !***********************************************************************
  !  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE
  !  LEGENDRE GAUSS-LOBATTO NODES
  !***********************************************************************
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: n ! ORDER OF LEGENDRE GAUSS-LOBATTO POLYNOMIAL
  real(lp), intent(in) :: a ! REAL VALUE TO DISTINGUISH BETWEEN DP AND QP VERSIONS
  !
  !.. Function Result ..
  real(lp), dimension(1:n+1,1:n+1) :: return_value ! DERIVATIVE MATRIX
  !
  !.. Local Scalars ..
  integer  :: i,j,np
  real(lp) :: ei,ej,vi,vj,dn,c
  !
  !.. Local Arrays ..
  real(lp), dimension(1:n+1) :: et ! VECTOR OF THE NODES
  real(lp), dimension(1:n+1) :: vn ! VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one4 = 0.25_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  np = n + 1
  !
  ! Compute the zeros of the Legendre Gauss Lobatto polynomial of
  ! order 'n' and the value of the Legendre Gauss polynomial of the
  ! same order at each zero
  !
  call nodes_legendre_gauss_lobatto(n,et,vn)
  !
  do i = 1,np
    vi = vn(i)
    ei = et(i)
    do j = 1,np
      if (i /= j) then
        vj = vn(j)
        ej = et(j)
        return_value(i,j) = vi/(vj*(ei-ej))
      else
        return_value(i,i) = zero
      end if
    end do
  end do
  !
  if (np > 1) then
    dn = real(np-1,kind=lp)
    c  = one4 * dn * (dn+one)
    return_value(1,1) = -c
    return_value(np,np) = c
  end if
  !
end function derivative_matrix_legendre_gauss_lobatto_DP
!
!################################################################################
!
pure function eval_LagrangePoly1D_DP(k,x,xi) result(return_value)
  !
  ! Get the value of the 1D Lagrange polynomial at point x
  !
  ! K  : The value of the Lagrange polynomial is 1 at xi(k)
  ! X  : Location at which to evaluate the k-th Lagrange polynomial
  ! XI : Array of all the interpolation points
  !
  ! This evaluates the Lagrange polynomial at the location x.
  ! The Lagrange polynomial corresponds to a function that
  ! is zero at all iterpolation points xi(i) for i=1,nfp
  ! except at xi(k) where the function is equal to 1.
  !
  !.. Formal Arguments ..
  integer,                 intent(in) :: k
  real(dp),               intent(in) :: x
  real(dp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(dp) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
  !
  integer, parameter :: lp = dp
  !
continue
  !
  nfp = size(xi)
  imask(:) = (/ (i,i=1,nfp) /)
  !
  return_value = product( x    -xi(1:nfp) , imask /= k ) / &
                 product( xi(k)-xi(1:nfp) , imask /= k )
  !
end function eval_LagrangePoly1D_DP
!
!###############################################################################
!
pure function eval_LagrangePoly2D_DP(k,l,x,y,xi,eta) result(return_value)
  !
  ! Get value of the 2D Lagrange polynomial at point x,y
  !
  ! K   : The value of the Lagrange polynomial in coordinate XI  is 1 at xi(k)
  ! L   : The value of the Lagrange polynomial in coordinate ETA is 1 at eta(l)
  ! X   : Location at which to evaluate the k-th Lagrange polynomial in coordinate XI
  ! Y   : Location at which to evaluate the l-th Lagrange polynomial in coordinate ETA
  ! XI  : Array of all the interpolation points in XI  coordinate
  ! ETA : Array of all the interpolation points in ETA coordinate
  !
  ! This evaluates the 2D Lagrange polynomial at the location x,y.
  ! The 2D Lagrange polynomial corresponds to a function that
  ! is zero at all iterpolation points [xi(i),eta(j)] for [i=1,ifp;j=1,jfp]
  ! except at [xi(k),eta(l)] where the function is equal to 1.
  !
  !.. Formal Arguments ..
  integer,                 intent(in) :: k
  integer,                 intent(in) :: l
  real(dp),               intent(in) :: x
  real(dp),               intent(in) :: y
  real(dp), dimension(:), intent(in) :: xi
  real(dp), dimension(:), intent(in) :: eta
  !
  !.. Function Return Value ..
  real(dp) :: return_value
  !
  integer, parameter :: lp = dp
  !
continue
  !
  return_value = eval_LagrangePoly1D(k,x,xi) * eval_LagrangePoly1D(l,y,eta)
  !
end function eval_LagrangePoly2D_DP
!
!###############################################################################
!
pure function eval_LagrangePoly_DP(ijk,xyz,xi) result(return_value)
  !
  ! Get the value of a tensor product of Lagrange polynomials of arbitrary dimension
  !
  ! IJK : The value of the Lagrange polynomial tensor product is 1 at ijk
  ! XYZ : Location at which to evaluate the Lagrange polynomial
  ! XI  : Array of the 1D solution points
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,  dimension(:), intent(in) :: ijk
  real(lp), dimension(:), intent(in) :: xyz
  real(lp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(lp) :: return_value
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Constants ..
  real(lp), parameter :: one = 1.0_lp
  !
continue
  !
  return_value = one
  do i = 1,size(ijk)
    return_value = return_value * eval_LagrangePoly1D(ijk(i),xyz(i),xi)
  end do
  !
end function eval_LagrangePoly_DP
!
!###############################################################################
!
pure function eval_DLagrangeDx_DP(m,x,xi) result(return_value)
  !
  ! Get the value of the derivative of the
  ! m-th Lagrange polynomial at the point x.
  !
  !                  nfp                nfp
  !   d              ___               _____
  !  -- [L_m (x)] =  \         1        | |   x   - x(j)
  !  dx              /__  -----------   | |  -----------
  !                       x(m) - x(l)   j=1  x(m) - x(j)
  !                  l=1               j/=l
  !                 l/=m               j/=m
  !
  ! m  : index where the value at xi(m) is 1
  ! X  : x location at which to evaluate the derivative
  !      of the m-th Lagrange polynomial 
  ! XI : array of all the interpolation points
  !
  !.. Formal Arguments ..
  integer,                 intent(in) :: m
  real(dp),               intent(in) :: x
  real(dp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(dp) :: return_value
  !
  !.. Local Scalars ..
  integer :: l,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
  logical, dimension(1:size(xi)) :: lmask
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  nfp = size(xi)
  !
  imask(:) = (/ (l,l=1,nfp) /)
  !
  return_value = zero
  !
  if (nfp > 1) then
    !
    do l = 1,nfp
      if (l == m) cycle
      lmask = (imask/=m .and. imask/=l)
      return_value = return_value + product( x    -xi(1:nfp) , lmask ) / &
                                    product( xi(m)-xi(1:nfp) , lmask ) / &
                                    (xi(m) - xi(l))
    end do
    !
  end if
  !
end function eval_DLagrangeDx_DP
!
!###############################################################################
!
pure function eval_D2LagrangeDx2_DP(m,x,xi) result(return_value)
  !
  ! Get value of second derivative of
  ! m-th Lagrange polynomial at point x.
  !
  !                   nfp                nfp                nfp
  !  d2               ___                ___               _____
  ! --- [L_m (x)]  =  \         1        \         1        | |   x   - x(j)
  ! dx2               /__  -----------   /__  -----------   | |  -----------
  !                        x(m) - x(l)        x(m) - x(n)   j=1  x(m) - x(j)
  !                   l=1                n=1               j/=l
  !                  l/=m               n/=l               j/=n
  !                                     n/=m               j/=m
  !
  ! m  : index where the value at xi(m) is 1
  ! X  : x location at which to evaluate the second
  !      derivative of the m-th Lagrange polynomial 
  ! XI : array of all the interpolation points
  !
  !.. Formal Arguments ..
  integer,                 intent(in) :: m
  real(dp),               intent(in) :: x
  real(dp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(dp) :: return_value
  !
  !.. Local Scalars ..
  integer :: l,n,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
  logical, dimension(1:size(xi)) :: lmask
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  nfp = size(xi)
  !
  imask(:) = (/ (l,l=1,nfp) /)
  !
  return_value = zero
  !
  if (nfp > 2) then
    !
    do l = 1,nfp
      if (l == m) cycle
      do n = 1,nfp
        if ( (n == m) .or. (n == l) ) cycle
        lmask = (imask/=m .and. imask/=l .and. imask/=n)
        return_value = return_value + product( x    -xi(1:nfp) , lmask ) / &
                                      product( xi(m)-xi(1:nfp) , lmask ) / &
                                      (xi(m) - xi(n))
      end do
      return_value = return_value / (xi(m) - xi(l))
    end do
    !
  end if
  !
end function eval_D2LagrangeDx2_DP
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
!  QUAD PRECISION VERSIONS OF THE POLYNOMIAL PROCEDURES
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
pure subroutine JacobiGL_QP(alpha,beta,n,xi,weights)
  !
  integer,                      intent(in) :: n
  !
  real(QP),                    intent(in) :: alpha
  real(QP),                    intent(in) :: beta
  !
  real(QP), dimension(1:n+1), intent(out) :: xi
  real(QP), dimension(1:n+1), intent(out) :: weights
  !
  real(QP), dimension(1:n+1) :: vn
  !
continue
  !
  xi = zero_QP
  weights = zero_QP
  !
  call nodes_jacobi_gauss_lobatto(n,alpha,beta,xi,vn)
  !
  call weights_dgbook_jacobi_gauss_lobatto(n,alpha,beta,xi,weights)
  !
end subroutine JacobiGL_QP
!
!###############################################################################
!
pure subroutine eval_legendre_poly_QP(n,x,y,dy,d2y)
  !
  integer,    intent(in) :: n
  !
  real(QP),  intent(in) :: x
  real(QP), intent(out) :: y
  real(QP), intent(out) :: dy
  real(QP), intent(out) :: d2y
  !
  real(QP) :: alpha,beta
  !
continue
  !
  alpha = zero_QP
  beta  = zero_QP
  !
  call eval_jacobi_poly(n,alpha,beta,x,y,dy,d2y)
  !
end subroutine eval_legendre_poly_QP
!
!###############################################################################
!
pure function normalized_jacobi_poly_QP(N,alpha,beta,x) result(return_value)
  !
  integer,   intent(in) :: N
  !
  real(QP), intent(in) :: alpha
  real(QP), intent(in) :: beta
  real(QP), intent(in) :: x
  !
  real(QP) :: return_value
  !
  integer   :: i
  !
  real(QP) :: ab,ab1,ab2,ab3,amb,a1,b1
  real(QP) :: gamma0,gamma1
  real(QP) :: aold,anew,bnew,h1,ri
  !
  real(QP) :: PL(N+1)
  !
continue
  !
  ab = alpha + beta
  amb = alpha - beta
  a1 = alpha + one_QP
  b1 = beta + one_QP
  ab1 = ab + one_QP
  ab2 = ab + two_QP
  ab3 = ab + three_QP
  !
  gamma0 = two_QP**ab1 / ab1 * gammafun(a1) * gammafun(b1) / gammafun(ab1)
  !
  PL(1) = one_QP/sqrt(gamma0)
  !
  if (N == 0) then
    !
    return_value = PL(1)
    !
  else
    !
    gamma1 = a1 * b1 / ab3 * gamma0
    !
    PL(2) = half_QP*( ab2*x + amb ) / sqrt(gamma1)
    !
    aold = two_QP / ab2 * sqrt( a1*b1/ab3 )
    !
    do i = 1,N-1
      ri = real(i,kind=QP)
      h1 = two_QP*ri + ab
      anew = two_QP / (h1+two_QP) * sqrt( (ri+one_QP)*(ri+ab1)*(ri+a1)*(ri+b1)/(h1+one_QP)/(h1+three_QP) )
      bnew = (beta*beta - alpha*alpha) / h1 / (h1+two_QP)
      PL(i+2) = one_QP/anew*( -aold*PL(i) + (x-bnew)*PL(i+1) )
      aold = anew
    end do
    !
    return_value = PL(N+1)
    !
  end if
  !
end function normalized_jacobi_poly_QP
!
!###############################################################################
!
pure function grad_normalized_jacobi_poly_QP(n,alpha,beta,x) result(return_value)
  !
  integer, parameter :: lp = qp
  !
  integer,  intent(in) :: n
  real(lp), intent(in) :: alpha
  real(lp), intent(in) :: beta
  real(lp), intent(in) :: x
  !
  real(lp) :: return_value
  !
  !.. Local Scalars ..
  real(lp) :: rn,a1,b1,ab1
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !
  if (n == 0) then
    return_value = zero
  else
    !
    rn  = real(n,kind=lp)
    a1  = alpha + one
    b1  = beta  + one
    ab1 = alpha + b1
    !
    return_value = sqrt(rn*(rn+ab1)) * normalized_jacobi_poly(n-1,a1,b1,x)
    !
  end if
  !
end function grad_normalized_jacobi_poly_QP
!
!###############################################################################
!
pure subroutine eval_normalized_jacobi_poly_QP(n,alpha,beta,x,y,dy)
  !
  integer,    intent(in) :: n
  !
  real(QP),  intent(in) :: alpha
  real(QP),  intent(in) :: beta
  real(QP),  intent(in) :: x
  real(QP), intent(out) :: y
  real(QP), intent(out) :: dy
  !
  integer :: n1
  !
  real(QP) :: a1,b1,ab1,rn
  !
continue
  !
  a1  = alpha + one_QP
  b1  = beta  + one_QP
  ab1 = alpha + beta + one_QP
  rn  = real(n,kind=QP)
  n1  = n - 1
  !
  y = normalized_jacobi_poly(n,alpha,beta,x)
  !
  dy = sqrt(rn*(rn+ab1))*normalized_jacobi_poly(n1,a1,b1,x)
  !
end subroutine eval_normalized_jacobi_poly_QP
!
!###############################################################################
!
pure subroutine eval_jacobi_poly_QP(n,alpha,beta,x,y,dy,d2y)
  !
  integer,    intent(in) :: n
  !
  real(QP),  intent(in) :: alpha
  real(QP),  intent(in) :: beta
  real(QP),  intent(in) :: x
  real(QP), intent(out) :: y
  real(QP), intent(out) :: dy
  real(QP), intent(out) :: d2y
  !
  integer :: i
  !
  real(QP) :: ab,di,c0,c1,c2,c3,c4
  real(QP) :: ym,yp,dym,dyp,d2ym,d2yp
  !
continue
  !
  y   = one_QP
  dy  = zero_QP
  d2y = zero_QP
  !
  if (n > 0) then
    !
    ab  = alpha + beta
    !
    y   = half_QP*(ab + two_QP)*x + half_QP*(alpha - beta)
    dy  = half_QP*(ab + two_QP)
    d2y = zero_QP
    !
  end if
  !
  if (n > 1) then
    !
    yp   = one_QP
    dyp  = zero_QP
    d2yp = zero_QP
    !
    do i = 2,n
      !
      di = real(i,kind=QP)
      c0 = two_QP*di + ab
      c1 = two_QP * di * (di+ab) * (c0-two_QP)
      c2 = (c0-one_QP) * (c0-two_QP) * c0
      c3 = (c0-one_QP) * (alpha-beta) * ab
      c4 = two_QP * c0 * (di+alpha-one_QP) * (di+beta-one_QP)
      !
      ym = y
      y  = ((c2*x + c3)*y - c4*yp)/c1
      yp = ym
      !
      dym = dy
      dy  = ((c2*x + c3)*dy - c4*dyp + c2*yp)/c1
      dyp = dym
      !
      d2ym = d2y
      d2y  = ((c2*x + c3)*d2y - c4*d2yp + two_QP*c2*dyp)/c1
      d2yp = d2ym
      !
    end do
    !
  end if
  !
end subroutine eval_jacobi_poly_QP
!
!###############################################################################
!
pure subroutine eval_chebyshev_poly_QP(n,x,y,dy,d2y)
  !
  integer,    intent(in) :: n
  !
  real(QP),  intent(in) :: x
  real(QP), intent(out) :: y
  real(QP), intent(out) :: dy
  real(QP), intent(out) :: d2y
  !
  integer :: k
  !
  real(QP) :: yp,ym,dyp,dym,d2yp,d2ym
  !
continue
  !
  y   = one_QP
  dy  = zero_QP
  d2y = zero_QP
  !
  if (n > 0) then
    y   = x
    dy  = one_QP
    d2y = zero_QP
  end if
  !
  if (n > 1) then
    !
    yp   = one_QP
    dyp  = zero_QP
    d2yp = zero_QP
    !
    do k = 2,n
      !
      ym = y
      y  = two_QP*x*y - yp
      yp = ym
      !
      dym = dy
      dy  = two_QP*x*dy + two_QP*yp - dyp
      dyp = dym
      !
      d2ym = d2y
      d2y  = two_QP*x*d2y + four_QP*dyp - d2yp
      d2yp = d2ym
      !
    end do
    !
  end if
  !
end subroutine eval_chebyshev_poly_QP
!
!###############################################################################
!
pure subroutine nodes_jacobi_gauss_QP(n,alpha,beta,xi,dxi)
  !
  integer,                    intent(in) :: n
  !
  real(QP),                  intent(in) :: alpha
  real(QP),                  intent(in) :: beta
  real(QP), dimension(1:n), intent(out) :: xi
  real(QP), dimension(1:n), intent(out) :: dxi
  !
  integer :: i,it
  !
  real(QP) :: const,ri,x,y,dy,d2y
  !
continue
  !
  if (n > 0) then
    xi(1)  = (beta - alpha) / (alpha + beta + two_QP)
    dxi(1) = half_QP*(alpha + beta) + one_QP
  end if
  !
  if (n > 1) then
    !
    const = half_QP*pi_QP / (two_QP*real(n,kind=QP) + alpha + beta + one_QP)
    !
    do i = 1,n
      !
      ri = real(i,kind=QP)
      x  = -cos( const * (four_QP*ri + alpha + beta - one_QP) )
      !
      do it = 1,8
        call eval_jacobi_poly(n,alpha,beta,x,y,dy,d2y)
        if (abs(y) <= eps17_QP) exit
        x = x - y/dy
      end do
      !
      if (abs(x) <= eps17_QP) x = zero_QP
      !
      xi(i)  = x
      dxi(i) = dy
      !
    end do
    !
    do i = 1,n/2
      xi (n-i+1) = -xi (i)
     !dxi(n-i+1) = -dxi(i)
    end do
    !
  end if
  !
end subroutine nodes_jacobi_gauss_QP
!
!###############################################################################
!
pure subroutine nodes_legendre_gauss_QP(n,xi,dxi)
  !
  integer,                    intent(in) :: n
  !
  real(QP), dimension(1:n), intent(out) :: xi
  real(QP), dimension(1:n), intent(out) :: dxi
  !
  real(QP) :: alpha,beta
  !
continue
  !
  call nodes_jacobi_gauss(n,zero_QP,zero_QP,xi,dxi)
  !
end subroutine nodes_legendre_gauss_QP
!
!###############################################################################
!
pure subroutine nodes_right_radau_QP(n,xi)
  !
  integer,                        intent(in) :: n
  !
  real(QP), dimension(1:n+1), intent(inout) :: xi
  !
  integer :: i,it
  real(QP) :: x,y,dy
  !
  real(QP), dimension(0:n+1) :: cfrr
  real(QP), dimension(0:n) :: cdrr
  !
continue
  !
  ! Set the final zero_QP which will always be unity
  !
  xi(n+1) = one_QP
  !
  if (n > 0) then
    !
    ! Get the coefficients of the right Radau polynomial
    !
    cfrr(:) = coefrr_QP(n+1)
    !
    do i = 0,n
      cdrr(i) = real(i+1,kind=QP) * cfrr(i+1)
    end do
    !
    xi(n+1) = one_QP
    !
    do i = 1,n
      !
      x  = -cos( pi_QP * real(i,kind=QP) / real(n+1,kind=QP) )
      !
      do it = 1,8
        !
        y = sum( (x**intseq(0,n+1)) * cfrr(:) )
        dy = sum( (x**intseq(0,n)) * cdrr(:) )
        !
        x = x - y/dy
        !
      end do
      !
      xi(i) = x
      !
    end do
    !
  end if
  !
end subroutine nodes_right_radau_QP
!
!###############################################################################
!
pure function coefrr_QP(n) result(return_value)
  !
  integer, intent(in) :: n
  !
  real(QP), dimension(0:n) :: return_value
  !
  integer :: k,ip
  !
  real(QP), dimension(0:n)   :: cfk
  real(QP), dimension(0:n)   :: cfkmt
  real(QP), dimension(0:n-1) :: cfkm1
  !
continue
  !
  cfk = coefle_QP(n)
  cfkm1 = coefle_QP(n-1)
  !
  do k = 0,n-1
    cfkmt(k) = cfkm1(k)
  end do
  !
  cfkmt(n) = zero_QP
  ip = (-1)**n
  !
  do k = 0,n
    return_value(k) = half_QP * real(ip,kind=QP) * (cfk(k) - cfkmt(k))
  end do
  !
end function coefrr_QP
!
!###############################################################################
!
pure function coefle_QP(n) result(return_value)
  !
  integer, intent(in) :: n
  !
  real(QP), dimension(0:n) :: return_value
  !
  integer :: i,j,k,m
  !
  real(QP), dimension(-1:n,-1:n) :: cle
  !
continue
  !
  cle(:,:) = zero_QP
  !
  cle(0,0) = one_QP
  cle(1,1) = one_QP
  !
  do k = 2,n
    do m = 0,k
      cle(k,m) = (real(2*k-1,kind=QP)/real(k,kind=QP))*cle(k-1,m-1) - &
                 (real(  k-1,kind=QP)/real(k,kind=QP))*cle(k-2,m  )
    end do
  end do
  !
  do j = 0,n
    return_value(j) = cle(n,j)
  end do
  !
end function coefle_QP
!
!###############################################################################
!
pure subroutine nodes_chebyshev_gauss_QP(n,xi)
  !
  integer,                    intent(in) :: n
  !
  real(QP), dimension(1:n), intent(out) :: xi
  !
  real(QP) :: alpha,beta
  real(QP), dimension(1:n) :: dxi
  !
continue
  !
  alpha = -half_QP
  beta  = -half_QP
  !
  call nodes_jacobi_gauss(n,alpha,beta,xi,dxi)
  !
end subroutine nodes_chebyshev_gauss_QP
!
!###############################################################################
!
pure subroutine nodes_jacobi_gauss_lobatto_QP(n,alpha,beta,xi,vn)
  !
  integer,                     intent(in) :: n
  !
  real(QP),                   intent(in) :: alpha
  real(QP),                   intent(in) :: beta
  real(QP), dimension(1:n+1), intent(out) :: xi
  real(QP), dimension(1:n+1), intent(out) :: vn
  !
  integer :: i,it,n1
  !
  real(QP) :: x,y,dy,d2y,ri,rn,ab,a1,b1,const
  !
continue
  !
  if (n == 0) then
    !
    x = zero_QP
    call eval_jacobi_poly(n,alpha,beta,x,y,dy,d2y)
    xi(1) = x
    vn(1) = y
    !
  else
    !
    x = -one_QP
    call eval_jacobi_poly(n,alpha,beta,x,y,dy,d2y)
    xi(1) = x
    vn(1) = y
    !
    x = one_QP
    call eval_jacobi_poly(n,alpha,beta,x,y,dy,d2y)
    xi(n+1) = x
    vn(n+1) = y
    !
  end if
  !
  if (n > 1) then
    !
    rn = real(n,kind=QP)
    ab = alpha + beta
    a1 = alpha + one_QP
    b1 = beta  + one_QP
    n1 = n-1
    const = half_QP * pi_QP / (two_QP*rn + ab + one_QP)
    !
    do i = 2,n
      !
      ri = real(i-1,kind=QP)
      x  = -cos( const * (four_QP*ri + ab + one_QP) )
      !
      do it = 1,8
        call eval_jacobi_poly(n1,a1,b1,x,y,dy,d2y)
        if (abs(y) <= eps17_QP) exit
        x = x - y/dy
      end do
      !
      if (abs(x) <= eps17_QP) x = zero_QP
      !
      xi(i) = x
      vn(i) = -half_QP * dy * (one_QP - x*x) / rn
      !
    end do
    !
  end if
  !
end subroutine nodes_jacobi_gauss_lobatto_QP
!
!###############################################################################
!
pure subroutine nodes_legendre_gauss_lobatto_QP(n,xi,vn)
  !
  integer,                      intent(in) :: n
  real(QP), dimension(1:n+1), intent(out) :: xi
  real(QP), dimension(1:n+1), intent(out) :: vn
  !
  integer :: i,it,n2
  !
  real(QP) :: sn,x,y,dy,d2y
  !
continue
  !
  if (n == 0) then
    xi(1) = zero_QP
    vn(1) = one_QP
  else
    !
    n2 = (n-1)/2
    sn = real(2*n - 4*n2 - 3,kind=QP)
    !
    xi(1) = -one_QP
    vn(1) =  sn
    !
    xi(n+1) = one_QP
    vn(n+1) = one_QP
    !
  end if
  !
  if (n > 1) then
    x = zero_QP
    call eval_legendre_poly(n,x,y,dy,d2y)
    xi(n2+2) = x
    vn(n2+2) = y
  end if
  !
  if (n > 2) then
    !
    do i = 1,n2
      !
      x = cos( pi_QP * real(i,kind=QP) / real(n,kind=QP) )
      !
      do it = 1,8
        call eval_legendre_poly(n,x,y,dy,d2y)
        x = x - dy/d2y
      end do
      !
      xi(i+1)   = -x
      xi(n-i+1) =  x
      !
      vn(i+1)   = y*sn
      vn(n-i+1) = y
      !
    end do
    !
  end if
  !
end subroutine nodes_legendre_gauss_lobatto_QP
!
!###############################################################################
!
pure subroutine nodes_chebyshev_gauss_lobatto_QP(n,xi)
  !
  integer,                      intent(in) :: n
  !
  real(QP), dimension(1:n+1), intent(out) :: xi
  !
  integer :: i,n2
  !
continue
  !
  if (n == 0) then
    xi(1) = zero_QP
  else
    xi(1)   = -one_QP
    xi(n+1) =  one_QP
  end if
  !
  if (n > 1) then
    n2 = (n-1)/2
    xi(n2+2) = zero_QP
  end if
  !
  if (n > 2) then
    do i = 1,n2
      xi(i+1)   = -cos( pi_QP * real(i,kind=QP) / real(n,kind=QP) )
      xi(n-i+1) =  cos( pi_QP * real(i,kind=QP) / real(n,kind=QP) )
    end do
  end if
  !
end subroutine nodes_chebyshev_gauss_lobatto_QP
!
!###############################################################################
!
elemental function gammafun_QP(x) result(return_value)
  !
  real(QP), intent(in) :: x
  !
  real(QP) :: return_value
  !
  integer   :: ind,k
  real(QP) :: xx,pr,s,g
  !
  real(QP), parameter :: s_start = 0.426401432711220868_QP
  real(QP), parameter :: gamma_constant(1:11) = &
             [ -0.524741987629368444e+00_QP,  0.116154405493589130e+00_QP, &
               -0.765978624506602380e-02_QP,  0.899719449391378898e-04_QP, &
               -0.194536980009534621e-07_QP,  0.199382839513630987e-10_QP, &
               -0.204209590209541319e-11_QP,  0.863896817907000175e-13_QP, &
                0.152237501608472336e-13_QP, -0.825725175277719950e-14_QP, &
                0.299734782205224610e-14_QP ]
  !
continue
  !
  xx = x
  return_value = one_QP
  !
  do
    if (abs(xx-one_QP) < eps14_QP) return
    if (xx < one_QP) exit
    xx = xx - one_QP
    return_value = return_value*xx
  end do
  !
  ind = 0
  !
  if (xx < half_QP) then
    ind = 1
    return_value = return_value*pi_QP/sin(pi_QP*xx)
    xx = one_QP - xx
  end if
  !
  pr = one_QP
  s = s_start
  !
  do k = 1,11
    pr = pr * (xx-real(k,kind=QP)) / (xx+real(k-1,kind=QP))
    s  = s + gamma_constant(k)*pr
  end do
  !
  g = s * exp(one_QP-xx) * (xx+c4_50_QP)**(xx-half_QP)
  !
  if (ind == 1) then
    return_value = return_value/g
  else
    return_value = return_value*g
  end if
  !
end function gammafun_QP
!
!###############################################################################
!
elemental function log_gammafun_QP(x) result(return_value)
  real(QP), intent(in) :: x
  real(QP) :: return_value
  return_value = log( gammafun(x) )
end function log_gammafun_QP
!
!###############################################################################
!
pure subroutine weights_jacobi_gauss_QP(n,alpha,beta,xi,dy,weights)
  !***************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE JACOBI GAUSS FORMULA
  !   N  = ORDER OF THE FORMULA
  !   alpha  = PARAMETER > -1
  !   beta  = PARAMETER > -1
  !   CS = ZEROES OF THE JACOBI POLYNOMIAL, CS(I), I=1,N
  !   DZ = VECTOR OF THE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
  !   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
  !***************************************************************
  integer,                    intent(in) :: n
  !
  real(QP),                  intent(in) :: alpha
  real(QP),                  intent(in) :: beta
  !
  real(QP), dimension(1:n),  intent(in) :: xi
  real(QP), dimension(1:n),  intent(in) :: dy
  real(QP), dimension(1:n), intent(out) :: weights
  !
  integer :: i,m
  !
  real(QP) :: ab,a2,b2,c,ri
  !
continue
  !
  weights(:) = one_QP
  !
  ab = alpha + beta + two_QP
  a2 = alpha + two_QP
  b2 = beta  + two_QP
  !
  c = half_QP * (two_QP**ab) * gammafun(a2) * gammafun(b2) / gammafun(ab)
  !
  do i = 2,n
    ri = real(i,kind=QP)
    c  = c * (ri+alpha) * (ri+beta) / (ri*(ri+alpha+beta))
  end do
  !
  do i = 1,n
    weights(i) = c / ((one_QP - xi(i)*xi(i)) * dy(i) * dy(i))
  end do
  !
end subroutine weights_jacobi_gauss_QP
!
!################################################################################
!
pure subroutine weights_legendre_gauss_QP(n,xi,dy,weights)
  !****************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS FORMULA
  !   N  = ORDER OF THE FORMULA
  !   CS = ZEROES OF THE LEGENDRE POLYNOMIAL, CS(I), I=1,N
  !   DZ = VECTOR OF THE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
  !   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
  !****************************************************************
  integer,                    intent(in) :: n
  !
  real(QP), dimension(1:n),  intent(in) :: xi
  real(QP), dimension(1:n),  intent(in) :: dy
  real(QP), dimension(1:n), intent(out) :: weights
  !
  integer :: i
  !
continue
  !
  weights(:) = two_QP
  !
  do i = 1,n
    weights(i) = two_QP / ((one_QP - xi(i)*xi(i)) * dy(i) * dy(i))
  end do
  !
end subroutine weights_legendre_gauss_QP
!
!################################################################################
!
pure subroutine weights_chebyshev_gauss_QP(n,weights)
  !****************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE CHEBYSHEV GAUSS FORMULA
  !   N  = ORDER OF THE FORMULA
  !   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
  !****************************************************************
  integer,                    intent(in) :: n
  !
  real(QP), dimension(1:n), intent(out) :: weights
  !
  integer :: i
  !
continue
  !
  weights(:) = pi_QP
  !
  do i = 1,n
    weights(i) = pi_QP/real(n,kind=QP)
  end do
  !
end subroutine weights_chebyshev_gauss_QP
!
!################################################################################
!
pure subroutine weights_jacobi_gauss_lobatto_QP(n,alpha,beta,xi,weights)
  !********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE JACOBI GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   alpha  = PARAMETER > -1
  !   beta  = PARAMETER > -1
  !   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !********************************************************************
  integer,                      intent(in) :: n
  !
  real(QP),                    intent(in) :: alpha
  real(QP),                    intent(in) :: beta
  !
  real(QP), dimension(1:n+1),  intent(in) :: xi
  real(QP), dimension(1:n+1), intent(out) :: weights
  !
  integer :: i
  !
  real(QP) :: a1,a2,b1,b2,ab,ab1,ab2
  real(QP) :: c,c1,c2,c3,ri,rn,su
  real(QP) :: y,z,dy,d2y
  !
continue
  !
  a1  = alpha + one_QP
  a2  = alpha + two_QP
  b1  = beta  + one_QP
  b2  = beta  + two_QP
  ab  = alpha + beta
  ab1 = alpha + beta + one_QP
  ab2 = alpha + beta + two_QP
  !
  c = (two_QP**ab) * gammafun(a1) * gammafun(b1) / gammafun(ab2)
  !
  weights(:) = two_QP*c
  !
  if (n == 1) then
    weights(1)   = two_QP * c * a1 / ab2
    weights(n+1) = two_QP * c * b1 / ab2
  end if
  !
  if (n > 1) then
    !
    rn = real(n,kind=QP)
    c  = c * (two_QP*rn + ab) / (rn + ab1)
    c1 = c * a1 / (b2 * ab2)
    c2 = c * b1 / (a2 * ab2)
    c3 = half_QP * c * a1 * b1
    !
    do i = 2,n-1
      ri = real(i-1,kind=QP)
      c1 = c1*(ri+a1)*ri / ((ri+ab2)*(ri+b2))
      c2 = c2*(ri+b1)*ri / ((ri+ab2)*(ri+a2))
      c3 = c3*(ri+a1)*(ri+b1) / ((ri+two_QP)*(ri+ab1))
    end do
    !
    su = zero_QP
    do i = 2,n
      su = su + xi(i)
    end do
    !
    weights(1)   = c1*(rn-one_QP-su)
    weights(n+1) = c2*(rn-one_QP+su)
    !
    do i = 2,n
      call eval_jacobi_poly(n  ,alpha,beta,xi(i),y,dy,d2y)
      call eval_jacobi_poly(n-1,alpha,beta,xi(i),z,dy,d2y)
      weights(i) = -c3/y/dy
    end do
    !
  end if
  !
end subroutine weights_jacobi_gauss_lobatto_QP
!
!################################################################################
!
pure subroutine weights_legendre_gauss_lobatto_QP(n,vn,weights)
  !**********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N
  !   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !**********************************************************************
  integer,                      intent(in) :: n
  !
  real(QP), dimension(1:n+1),  intent(in) :: vn
  real(QP), dimension(1:n+1), intent(out) :: weights
  !
  integer :: i
  !
  real(QP) :: c
  !
continue
  !
  c = two_QP
  if (n > 0) then
    c = c / (real(n,kind=QP)*(real(n,kind=QP)+one_QP))
  end if
  !
  do i = 1,n+1
    weights(i) = c / (vn(i)*vn(i))
  end do
  !
end subroutine weights_legendre_gauss_lobatto_QP
!
!################################################################################
!
pure subroutine weights_right_radau_QP(n,xi,weights)
  !**********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   XI = RIGHT RADAU NODES, XI(I), I=0,N
  !   VN = VALUES OF THE RIGHT RADAU POLYNOMIAL AT THE NODES, VN(I), I=0,N
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !**********************************************************************
  integer,                      intent(in) :: n
  !
  real(QP), dimension(1:n+1),  intent(in) :: xi
  real(QP), dimension(1:n+1), intent(out) :: weights
  !
  integer :: i
  !
  real(QP) :: c,y,dy,d2y
  !
continue
  !
  c = two_QP
  !
  if (n > 0) then
    !
    c = one_QP / (real(n+1,kind=QP)*real(n+1,kind=QP))
    !
    do i = 1,n+1
      call eval_legendre_poly(n,xi(i),y,dy,d2y)
      weights(i) = c * (one_QP+xi(i)) / (y*y)
    end do
    !
  end if
  !
end subroutine weights_right_radau_QP
!
!################################################################################
!
pure subroutine weights_chebyshev_gauss_lobatto_QP(n,weights)
  !***********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE CHEBYSHEV GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !***********************************************************************
  integer,                      intent(in) :: n
  !
  real(QP), dimension(1:n+1), intent(out) :: weights
  !
  integer :: i
  !
continue
  !
  weights(:) = pi_QP
  !
  if (n > 0) then
    weights(1)   = half_QP*pi_QP/real(n,kind=QP)
    weights(n+1) = half_QP*pi_QP/real(n,kind=QP)
  end if
  !
  do i = 2,n
    weights(i) = pi_QP/real(n,kind=QP)
  end do
  !
end subroutine weights_chebyshev_gauss_lobatto_QP
!
!###############################################################################
!
pure subroutine weights_dgbook_jacobi_gauss_lobatto_QP(n,alpha,beta,xi,weights)
  !********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE JACOBI GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   alpha  = PARAMETER > -1
  !   beta  = PARAMETER > -1
  !   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !********************************************************************
  integer,                      intent(in) :: n
  !
  real(QP),                    intent(in) :: alpha
  real(QP),                    intent(in) :: beta
  !
  real(QP), dimension(1:n+1),  intent(in) :: xi
  real(QP), dimension(1:n+1), intent(out) :: weights
  !
  integer :: i
  !
  real(QP) :: a1,a2,b1,b2,ab,ab1,ab2
  real(QP) :: c,c1,c2,c3,ri,rn,su
  real(QP) :: y,z,dy,d2y,scl
  !
continue
  !
  weights(:) = zero_QP
  !
  if (n > 1) then
    !
    a1  = alpha + one_QP
    a2  = alpha + two_QP
    b1  = beta  + one_QP
    b2  = beta  + two_QP
    ab  = alpha + beta
    ab1 = alpha + beta + one_QP
    ab2 = alpha + beta + two_QP
    !
    c = (two_QP**ab) * gammafun(a1) * gammafun(b1) / gammafun(ab2)
    !
    rn = real(n,kind=QP)
    c  = c * (two_QP*rn + ab) / (rn + ab1)
    c1 = c * a1 / (b2 * ab2)
    c2 = c * b1 / (a2 * ab2)
    c3 = half_QP * c * a1 * b1
    !
    do i = 2,n-1
      ri = real(i-1,kind=QP)
      c1 = c1*(ri+a1)*ri / ((ri+ab2)*(ri+b2))
      c2 = c2*(ri+b1)*ri / ((ri+ab2)*(ri+a2))
      c3 = c3*(ri+a1)*(ri+b1) / ((ri+two_QP)*(ri+ab1))
    end do
    !
    su = zero_QP
    do i = 2,n
      su = su + xi(i)
    end do
    !
    scl = ( (one_QP-xi(1))**a1 ) * ( (one_QP+xi(1))**b1 )
    weights(1)   = scl*c1*(rn-one_QP-su)
    !
    scl = ( (one_QP-xi(n+1))**a1 ) * ( (one_QP+xi(n+1))**b1 )
    weights(n+1) = scl*c2*(rn-one_QP+su)
    !
    do i = 2,n
      call eval_jacobi_poly(n  ,alpha,beta,xi(i),y,dy,d2y)
      call eval_jacobi_poly(n-1,alpha,beta,xi(i),z,dy,d2y)
      scl = ( (one_QP-xi(i))**a1 ) * ( (one_QP+xi(i))**b1 )
      weights(i) = -scl*c3/y/dy
    end do
    !
  end if
  !
end subroutine weights_dgbook_jacobi_gauss_lobatto_QP
!
!################################################################################
!
pure subroutine old_derivative_matrix_legendre_gauss_QP(n,cs,dz,dma)
  !
  !*************************************************************
  !*************************************************************
  !********************                       ******************
  !********************     NOT IN MANUAL     ******************
  !********************                       ******************
  !*************************************************************
  !*************************************************************
  !
  !***********************************************************************
  !  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE
  !  LEGENDRE GAUSS NODES
  !  N   = NUMBER OF LEGENDRE GAUSS NODES
  !  CS  = VECTOR OF THE ZEROES, CS(I), I=1,N
  !  DZ  = LEGENDRE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
  !  DMA = DERIVATIVE MATRIX, DMA(I,J), I=1,N  J=1,N
  !***********************************************************************
  !
  !.. Formal Arguments ..
  integer,                        intent(in) :: n
  real(QP), dimension(1:n),      intent(in) :: cs
  real(QP), dimension(1:n),      intent(in) :: dz
  real(QP), dimension(1:n,1:n), intent(out) :: dma
  !
  !.. Local Scalars ..
  integer   :: i,j
  real(QP) :: vi,zi,vj,zj
  !
continue
  !
  dma(1,1) = zero_QP
  !
  if (n <= 1) return
  !
  do i = 1,n
    vi = dz(i)
    zi = cs(i)
    do j = 1,n
      if (i /= j) then
        vj = dz(j)
        zj = cs(j)
        dma(i,j) = vi/(vj*(zi-zj))
      else
        dma(i,i) = zi/(one_QP-zi*zi)
      end if
    end do
  end do
  !
end subroutine old_derivative_matrix_legendre_gauss_QP
!
!################################################################################
!
pure subroutine derivative_matrix_jacobi_gauss_lobatto_QP(n,alpha,beta,et,vn,dma)
  !
  !***********************************************************************
  !  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE
  !  JACOBI GAUSS-LOBATTO NODES
  !  N   = PARAMETER RELATIVE TO THE DIMENSION OF THE MATRIX
  !  ALPHA   = PARAMETER >-1
  !  BETA   = PARAMETER >-1
  !  ET  = VECTOR OF THE NODES, ET(I), I=0,N
  !  VN  = VALUES OF THE JACOBI POLYNOMIAL AT THE NODES, VN(I), I=0,N
  !  DMA = DERIVATIVE MATRIX, DMA(I,J), I=0,N  J=0,N
  !***********************************************************************
  !
  !.. Formal Arguments ..
  integer,                        intent(in) :: n
  real(QP),                      intent(in) :: alpha
  real(QP),                      intent(in) :: beta
  real(QP), dimension(1:n),      intent(in) :: et
  real(QP), dimension(1:n),      intent(in) :: vn
  real(QP), dimension(1:n,1:n), intent(out) :: dma
  !
  !.. Local Scalars ..
  integer   :: i,j
  real(QP) :: a1,b1,ab,dn,ei,ej
  real(QP) :: c1,c2,c3,c4,vi,vj
  !
continue
  !
  dma(1,1) = zero_QP
  if (n == 1) return
  !
  a1 = alpha+one_QP
  b1 = beta+one_QP
  ab = alpha+beta
  dn = real(n-1,kind=QP)
  c1 = dn*(dn+ab+one_QP)
  c2 = a1*vn(1)/(b1*vn(n))
  dma(1,1) = half_QP*(alpha-c1)/(beta+two_QP)
  dma(n,n) = half_QP*(c1-beta)/(alpha+two_QP)
  dma(1,n) = -half_QP*c2
  dma(n,1) = half_QP/c2
  if (n == 2) return
  !
  c3 = vn(1)/b1
  c4 = vn(n)/a1
  do j = 2,n-1
    vj = vn(j)
    ej = et(j)
    dma(1,j) = -c3/(vj*(one_QP+ej))
    dma(n,j) = c4/(vj*(one_QP-ej))
    dma(j,1) = vj/(c3*(one_QP+ej))
    dma(j,n) = -vj/(c4*(one_QP-ej))
  end do
  !
  do i = 2,n-1
    vi = vn(i)
    ei = et(i)
    do j = 2,n-1
      if (i /= j) then
        vj = vn(j)
        ej = et(j)
        dma(i,j) = vi/(vj*(ei-ej))
      else
        dma(i,i) = half_QP*(ab*ei+alpha-beta)/(one_QP-ei*ei)
      end if
    end do
  end do
  !
end subroutine derivative_matrix_jacobi_gauss_lobatto_QP
!
!################################################################################
!
pure subroutine old_derivative_matrix_legendre_gauss_lobatto_QP(n,et,vn,dma)
  !
  !***********************************************************************
  !  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE
  !  LEGENDRE GAUSS-LOBATTO NODES
  !  N   = NUMBER OF LEGENDRE GAUSS-LOBATTO NODES
  !  ET  = VECTOR OF THE NODES, ET(I), I=1,N
  !  VN  = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=1,N
  !  DMA = DERIVATIVE MATRIX, DMA(I,J), I=1,N  J=1,N
  !***********************************************************************
  !
  !.. Formal Arguments ..
  integer,                        intent(in) :: n
  real(QP), dimension(1:n),      intent(in) :: et
  real(QP), dimension(1:n),      intent(in) :: vn
  real(QP), dimension(1:n,1:n), intent(out) :: dma
  !
  !.. Local Scalars ..
  integer   :: i,j
  real(QP) :: ei,ej,vi,vj,dn,c
  !
continue
  !
  do i = 1,n
    vi = vn(i)
    ei = et(i)
    do j = 1,n
      if (i /= j) then
        vj = vn(j)
        ej = et(j)
        dma(i,j) = vi/(vj*(ei-ej))
      else
        dma(i,i) = zero_QP
      end if
    end do
  end do
  !
  if (n > 1) then
    dn = real(n-1,kind=QP)
    c  = c0_25_QP*dn*(dn+one_QP)
    dma(1,1) = -c
    dma(n,n) = c
  end if
  !
end subroutine old_derivative_matrix_legendre_gauss_lobatto_QP
!
!################################################################################
!
pure subroutine derivative_matrix_chebyshev_gauss_lobatto_QP(n,et,dma)
  !
  !***********************************************************************
  !  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE
  !  CHEBYSHEV GAUSS-LOBATTO NODES
  !  N   = PARAMETER RELATIVE TO THE DIMENSION OF THE MATRIX
  !  ET  = VECTOR OF THE NODES, ET(I), I=0,N
  !  DMA = DERIVATIVE MATRIX, DMA(I,J), I=0,N  J=0,N
  !***********************************************************************
  !
  !.. Formal Arguments ..
  integer,                        intent(in) :: n
  real(QP), dimension(1:n),      intent(in) :: et
  real(QP), dimension(1:n,1:n), intent(out) :: dma
  !
  !.. Local Scalars ..
  integer   :: i,j
  real(QP) :: dn,cn,ej,ei
  real(QP) :: sn,sgni,sgnj
  !
continue
  !
  dma(1,1) = zero_QP
  if (n == 1) return
  !
  dn = real((n-1),kind=QP)
  cn = (two_QP*dn*dn+one_QP)/six_QP
  sn = real(1+4*((n-1)/2)-2*(n-1),kind=QP)
  dma(1,1) = -cn
  dma(n,n) = cn
  dma(1,n) = -half_QP*sn
  dma(n,1) = half_QP*sn
  if (n == 2) return
  !
  sgnj = -one_QP
  do j = 2,n-1
    ej = et(j)
    dma(1,j) = -two_QP*sgnj/(one_QP+ej)
    dma(n,j) = two_QP*sgnj*sn/(one_QP-ej)
    dma(j,1) = half_QP*sgnj/(one_QP+ej)
    dma(j,n) = -half_QP*sgnj*sn/(one_QP-ej)
    sgnj = -sgnj
  end do
  !
  sgni = -one_QP
  do i = 2,n-1
    ei = et(i)
    sgnj = -one_QP
    do j = 2,n-1
      if (i /= j) then
        ej = et(j)
        dma(i,j) = sgni*sgnj/(ei-ej)
      else
        dma(i,i) = -half_QP*ei/(one_QP-ei*ei)
      end if
      sgnj = -sgnj
    end do
    sgni = -sgni
  end do
  !
end subroutine derivative_matrix_chebyshev_gauss_lobatto_QP
!
!################################################################################
!
pure function derivative_matrix_legendre_gauss_QP(n,a) result(return_value)
  !
  !*************************************************************
  !*************************************************************
  !********************                       ******************
  !********************     NOT IN MANUAL     ******************
  !********************                       ******************
  !*************************************************************
  !*************************************************************
  !
  !***********************************************************************
  !  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE
  !  LEGENDRE GAUSS NODES
  !***********************************************************************
  !
  !.. Function Precision ..
  integer, parameter :: lp = qp
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: n ! ORDER OF LEGENDRE GAUSS POLYNOMIAL
  real(lp), intent(in) :: a ! REAL VALUE TO DISTINGUISH BETWEEN DP AND QP VERSIONS
  !
  !.. Function Result ..
  real(lp), dimension(1:n+1,1:n+1) :: return_value ! DERIVATIVE MATRIX
  !
  !.. Local Scalars ..
  integer  :: i,j,np
  real(lp) :: vi,zi,vj,zj
  !
  !.. Local Arrays ..
  real(lp), dimension(1:n+1) :: cs ! VECTOR OF THE ZEROES
  real(lp), dimension(1:n+1) :: dz ! LEGENDRE DERIVATIVES AT THE ZEROES
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  np = n + 1
  !
  return_value(1,1) = zero
  !
  if (n <= 1) return
  !
  ! Compute the zeros of the Legendre Gauss polynomial of
  ! order 'n' and the value of the derivative at each zero
  !
  call nodes_legendre_gauss(np,cs,dz)
  !
  do i = 1,np
    vi = dz(i)
    zi = cs(i)
    do j = 1,np
      if (i /= j) then
        vj = dz(j)
        zj = cs(j)
        return_value(i,j) = vi/(vj*(zi-zj))
      else
        return_value(i,i) = zi/(one-zi*zi)
      end if
    end do
  end do
  !
end function derivative_matrix_legendre_gauss_QP
!
!################################################################################
!
pure function derivative_matrix_legendre_gauss_lobatto_QP(n,a) result(return_value)
  !
  !***********************************************************************
  !  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE
  !  LEGENDRE GAUSS-LOBATTO NODES
  !***********************************************************************
  !
  !.. Function Precision ..
  integer, parameter :: lp = qp
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: n ! ORDER OF LEGENDRE GAUSS-LOBATTO POLYNOMIAL
  real(lp), intent(in) :: a ! REAL VALUE TO DISTINGUISH BETWEEN DP AND QP VERSIONS
  !
  !.. Function Result ..
  real(lp), dimension(1:n+1,1:n+1) :: return_value ! DERIVATIVE MATRIX
  !
  !.. Local Scalars ..
  integer  :: i,j,np
  real(lp) :: ei,ej,vi,vj,dn,c
  !
  !.. Local Arrays ..
  real(lp), dimension(1:n+1) :: et ! VECTOR OF THE NODES
  real(lp), dimension(1:n+1) :: vn ! VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one4 = 0.25_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  np = n + 1
  !
  ! Compute the zeros of the Legendre Gauss Lobatto polynomial of
  ! order 'n' and the value of the Legendre Gauss polynomial of the
  ! same order at each zero
  !
  call nodes_legendre_gauss_lobatto(n,et,vn)
  !
  do i = 1,np
    vi = vn(i)
    ei = et(i)
    do j = 1,np
      if (i /= j) then
        vj = vn(j)
        ej = et(j)
        return_value(i,j) = vi/(vj*(ei-ej))
      else
        return_value(i,i) = zero
      end if
    end do
  end do
  !
  if (np > 1) then
    dn = real(np-1,kind=lp)
    c  = one4 * dn * (dn+one)
    return_value(1,1) = -c
    return_value(np,np) = c
  end if
  !
end function derivative_matrix_legendre_gauss_lobatto_QP
!
!################################################################################
!
pure function eval_LagrangePoly1D_QP(k,x,xi) result(return_value)
  !
  ! Get the value of the 1D Lagrange polynomial at point x
  !
  ! K  : The value of the Lagrange polynomial is 1 at xi(k)
  ! X  : Location at which to evaluate the k-th Lagrange polynomial
  ! XI : Array of all the interpolation points
  !
  ! This evaluates the Lagrange polynomial at the location x.
  ! The Lagrange polynomial corresponds to a function that
  ! is zero at all iterpolation points xi(i) for i=1,nfp
  ! except at xi(k) where the function is equal to 1.
  !
  !.. Formal Arguments ..
  integer,                 intent(in) :: k
  real(qp),               intent(in) :: x
  real(qp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(qp) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
  !
  integer, parameter :: lp = qp
  !
continue
  !
  nfp = size(xi)
  imask(:) = (/ (i,i=1,nfp) /)
  !
  return_value = product( x    -xi(1:nfp) , imask /= k ) / &
                 product( xi(k)-xi(1:nfp) , imask /= k )
  !
end function eval_LagrangePoly1D_QP
!
!###############################################################################
!
pure function eval_LagrangePoly2D_QP(k,l,x,y,xi,eta) result(return_value)
  !
  ! Get value of the 2D Lagrange polynomial at point x,y
  !
  ! K   : The value of the Lagrange polynomial in coordinate XI  is 1 at xi(k)
  ! L   : The value of the Lagrange polynomial in coordinate ETA is 1 at eta(l)
  ! X   : Location at which to evaluate the k-th Lagrange polynomial in coordinate XI
  ! Y   : Location at which to evaluate the l-th Lagrange polynomial in coordinate ETA
  ! XI  : Array of all the interpolation points in XI  coordinate
  ! ETA : Array of all the interpolation points in ETA coordinate
  !
  ! This evaluates the 2D Lagrange polynomial at the location x,y.
  ! The 2D Lagrange polynomial corresponds to a function that
  ! is zero at all iterpolation points [xi(i),eta(j)] for [i=1,ifp;j=1,jfp]
  ! except at [xi(k),eta(l)] where the function is equal to 1.
  !
  !.. Formal Arguments ..
  integer,                 intent(in) :: k
  integer,                 intent(in) :: l
  real(qp),               intent(in) :: x
  real(qp),               intent(in) :: y
  real(qp), dimension(:), intent(in) :: xi
  real(qp), dimension(:), intent(in) :: eta
  !
  !.. Function Return Value ..
  real(qp) :: return_value
  !
  integer, parameter :: lp = qp
  !
continue
  !
  return_value = eval_LagrangePoly1D(k,x,xi) * eval_LagrangePoly1D(l,y,eta)
  !
end function eval_LagrangePoly2D_QP
!
!###############################################################################
!
pure function eval_LagrangePoly_QP(ijk,xyz,xi) result(return_value)
  !
  ! Get the value of a tensor product of Lagrange polynomials of arbitrary dimension
  !
  ! IJK : The value of the Lagrange polynomial tensor product is 1 at ijk
  ! XYZ : Location at which to evaluate the Lagrange polynomial
  ! XI  : Array of the 1D solution points
  !
  !.. Function Precision ..
  integer, parameter :: lp = qp
  !
  !.. Formal Arguments ..
  integer,  dimension(:), intent(in) :: ijk
  real(lp), dimension(:), intent(in) :: xyz
  real(lp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(lp) :: return_value
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Constants ..
  real(lp), parameter :: one = 1.0_lp
  !
continue
  !
  return_value = one
  do i = 1,size(ijk)
    return_value = return_value * eval_LagrangePoly1D(ijk(i),xyz(i),xi)
  end do
  !
end function eval_LagrangePoly_QP
!
!###############################################################################
!
pure function eval_DLagrangeDx_QP(m,x,xi) result(return_value)
  !
  ! Get the value of the derivative of the
  ! m-th Lagrange polynomial at the point x.
  !
  !                  nfp                nfp
  !   d              ___               _____
  !  -- [L_m (x)] =  \         1        | |   x   - x(j)
  !  dx              /__  -----------   | |  -----------
  !                       x(m) - x(l)   j=1  x(m) - x(j)
  !                  l=1               j/=l
  !                 l/=m               j/=m
  !
  ! m  : index where the value at xi(m) is 1
  ! X  : x location at which to evaluate the derivative
  !      of the m-th Lagrange polynomial 
  ! XI : array of all the interpolation points
  !
  !.. Formal Arguments ..
  integer,                 intent(in) :: m
  real(qp),               intent(in) :: x
  real(qp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(qp) :: return_value
  !
  !.. Local Scalars ..
  integer :: l,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
  logical, dimension(1:size(xi)) :: lmask
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  nfp = size(xi)
  !
  imask(:) = (/ (l,l=1,nfp) /)
  !
  return_value = zero
  !
  if (nfp > 1) then
    !
    do l = 1,nfp
      if (l == m) cycle
      lmask = (imask/=m .and. imask/=l)
      return_value = return_value + product( x    -xi(1:nfp) , lmask ) / &
                                    product( xi(m)-xi(1:nfp) , lmask ) / &
                                    (xi(m) - xi(l))
    end do
    !
  end if
  !
end function eval_DLagrangeDx_QP
!
!###############################################################################
!
pure function eval_D2LagrangeDx2_QP(m,x,xi) result(return_value)
  !
  ! Get value of second derivative of
  ! m-th Lagrange polynomial at point x.
  !
  !                   nfp                nfp                nfp
  !  d2               ___                ___               _____
  ! --- [L_m (x)]  =  \         1        \         1        | |   x   - x(j)
  ! dx2               /__  -----------   /__  -----------   | |  -----------
  !                        x(m) - x(l)        x(m) - x(n)   j=1  x(m) - x(j)
  !                   l=1                n=1               j/=l
  !                  l/=m               n/=l               j/=n
  !                                     n/=m               j/=m
  !
  ! m  : index where the value at xi(m) is 1
  ! X  : x location at which to evaluate the second
  !      derivative of the m-th Lagrange polynomial 
  ! XI : array of all the interpolation points
  !
  !.. Formal Arguments ..
  integer,                 intent(in) :: m
  real(qp),               intent(in) :: x
  real(qp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(qp) :: return_value
  !
  !.. Local Scalars ..
  integer :: l,n,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
  logical, dimension(1:size(xi)) :: lmask
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  nfp = size(xi)
  !
  imask(:) = (/ (l,l=1,nfp) /)
  !
  return_value = zero
  !
  if (nfp > 2) then
    !
    do l = 1,nfp
      if (l == m) cycle
      do n = 1,nfp
        if ( (n == m) .or. (n == l) ) cycle
        lmask = (imask/=m .and. imask/=l .and. imask/=n)
        return_value = return_value + product( x    -xi(1:nfp) , lmask ) / &
                                      product( xi(m)-xi(1:nfp) , lmask ) / &
                                      (xi(m) - xi(n))
      end do
      return_value = return_value / (xi(m) - xi(l))
    end do
    !
  end if
  !
end function eval_D2LagrangeDx2_QP
!
!###############################################################################
!
end module polynomial_mod
