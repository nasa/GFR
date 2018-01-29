module polynomial_mod
  !
  ! Module that contains various functions and routines
  ! related to Jacobi, Legendre, and other polynomials
  !
  use module_kind_types
  !
  implicit none
#ifdef DISABLE_QP
  integer, parameter :: local_qp = SP
#else
  integer, parameter :: local_qp = QP
#endif
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
  public :: eval_LagrangePoly3D
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
  interface weights_right_radau
    module procedure weights_right_radau_DP, &
                     weights_right_radau_QP
  end interface weights_right_radau
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
  interface coefrr
    module procedure coefrr_DP, &
                     coefrr_QP
  end interface coefrr
  !
  interface coefle
    module procedure coefle_DP, &
                     coefle_QP
  end interface coefle
  !
  interface eval_LagrangePoly
    module procedure eval_LagrangePoly_DP, &
                     eval_LagrangePoly_QP
  end interface eval_LagrangePoly
  !
  interface eval_LagrangePoly1D
    module procedure eval_LagrangePoly1D_DP, &
                     eval_LagrangePoly1D_QP, &
                     eval_LagrangePoly1D_arrays_DP, &
                     eval_LagrangePoly1D_arrays_QP
  end interface eval_LagrangePoly1D
  !
  interface eval_LagrangePoly2D
    module procedure eval_LagrangePoly2D_DP, &
                     eval_LagrangePoly2D_QP, &
                     eval_LagrangePoly2D_arrays_DP, &
                     eval_LagrangePoly2D_arrays_QP
  end interface eval_LagrangePoly2D
  !
  interface eval_LagrangePoly3D
    module procedure eval_LagrangePoly3D_DP, &
                     eval_LagrangePoly3D_QP, &
                     eval_LagrangePoly3D_arrays_DP, &
                     eval_LagrangePoly3D_arrays_QP
  end interface eval_LagrangePoly3D
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
  real(wp), intent(inout) :: xi(:)
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(inout) :: weights(:)
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
  !
  if (present(weights)) then
    weights = chop( wts_qp )
  end if
  !
end subroutine legendre_gauss_quadrature
!
!###############################################################################
!
pure subroutine legendre_gauss_lobatto_quadrature(xi,weights)
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: xi(:)
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(inout) :: weights(:)
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
  !
  if (present(weights)) then
    weights = chop( wts_qp )
  end if
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
! PROCEDURES OF VARYING PRECISION FOR THE MODULE GENERIC INTERFACES
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
pure subroutine JacobiGL_DP(alpha,beta,n,xi,weights)
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: n
  real(lp),                    intent(in) :: alpha
  real(lp),                    intent(in) :: beta
  real(lp), dimension(1:n+1), intent(out) :: xi
  real(lp), dimension(1:n+1), intent(out) :: weights
  !
  !.. Local Arrays ..
  real(lp), dimension(1:n+1) :: vn
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  xi = zero
  weights = zero
  !
  call nodes_jacobi_gauss_lobatto(n,alpha,beta,xi,vn)
  !
  call weights_dgbook_jacobi_gauss_lobatto(n,alpha,beta,xi,weights)
  !
end subroutine JacobiGL_DP
!
!###############################################################################
!
pure subroutine JacobiGL_QP(alpha,beta,n,xi,weights)
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: n
  real(lp),                    intent(in) :: alpha
  real(lp),                    intent(in) :: beta
  real(lp), dimension(1:n+1), intent(out) :: xi
  real(lp), dimension(1:n+1), intent(out) :: weights
  !
  !.. Local Arrays ..
  real(lp), dimension(1:n+1) :: vn
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  xi = zero
  weights = zero
  !
  call nodes_jacobi_gauss_lobatto(n,alpha,beta,xi,vn)
  !
  call weights_dgbook_jacobi_gauss_lobatto(n,alpha,beta,xi,weights)
  !
end subroutine JacobiGL_QP
!
!###############################################################################
!
pure subroutine eval_jacobi_poly_DP(n,alpha,beta,x,y,dy,d2y)
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp),  intent(in) :: alpha
  real(lp),  intent(in) :: beta
  real(lp),  intent(in) :: x
  real(lp), intent(out) :: y
  real(lp), intent(out) :: dy
  real(lp), intent(out) :: d2y
  !
  !.. Local Scalars ..
  integer  :: i
  real(lp) :: ab,di,c0,c1,c2,c3,c4
  real(lp) :: ym,yp,dym,dyp,d2ym,d2yp
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: half = 0.5_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: two  = 2.0_lp
  !
continue
  !
  y   = one
  dy  = zero
  d2y = zero
  !
  if (n > 0) then
    !
    ab  = alpha + beta
    !
    y   = half*(ab + two)*x + half*(alpha - beta)
    dy  = half*(ab + two)
    d2y = zero
    !
  end if
  !
  if (n > 1) then
    !
    yp   = one
    dyp  = zero
    d2yp = zero
    !
    do i = 2,n
      !
      di = real(i,kind=lp)
      c0 = two*di + ab
      c1 = two * di * (di+ab) * (c0-two)
      c2 = (c0-one) * (c0-two) * c0
      c3 = (c0-one) * (alpha-beta) * ab
      c4 = two * c0 * (di+alpha-one) * (di+beta-one)
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
      d2y  = ((c2*x + c3)*d2y - c4*d2yp + two*c2*dyp)/c1
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
pure subroutine eval_jacobi_poly_QP(n,alpha,beta,x,y,dy,d2y)
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp),  intent(in) :: alpha
  real(lp),  intent(in) :: beta
  real(lp),  intent(in) :: x
  real(lp), intent(out) :: y
  real(lp), intent(out) :: dy
  real(lp), intent(out) :: d2y
  !
  !.. Local Scalars ..
  integer  :: i
  real(lp) :: ab,di,c0,c1,c2,c3,c4
  real(lp) :: ym,yp,dym,dyp,d2ym,d2yp
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: half = 0.5_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: two  = 2.0_lp
  !
continue
  !
  y   = one
  dy  = zero
  d2y = zero
  !
  if (n > 0) then
    !
    ab  = alpha + beta
    !
    y   = half*(ab + two)*x + half*(alpha - beta)
    dy  = half*(ab + two)
    d2y = zero
    !
  end if
  !
  if (n > 1) then
    !
    yp   = one
    dyp  = zero
    d2yp = zero
    !
    do i = 2,n
      !
      di = real(i,kind=lp)
      c0 = two*di + ab
      c1 = two * di * (di+ab) * (c0-two)
      c2 = (c0-one) * (c0-two) * c0
      c3 = (c0-one) * (alpha-beta) * ab
      c4 = two * c0 * (di+alpha-one) * (di+beta-one)
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
      d2y  = ((c2*x + c3)*d2y - c4*d2yp + two*c2*dyp)/c1
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
pure subroutine eval_legendre_poly_DP(n,x,y,dy,d2y)
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp),  intent(in) :: x
  real(lp), intent(out) :: y
  real(lp), intent(out) :: dy
  real(lp), intent(out) :: d2y
  !
  !.. Local Scalars ..
  real(lp) :: alpha,beta
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  alpha = zero
  beta  = zero
  !
  call eval_jacobi_poly(n,alpha,beta,x,y,dy,d2y)
  !
end subroutine eval_legendre_poly_DP
!
!###############################################################################
!
pure subroutine eval_legendre_poly_QP(n,x,y,dy,d2y)
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp),  intent(in) :: x
  real(lp), intent(out) :: y
  real(lp), intent(out) :: dy
  real(lp), intent(out) :: d2y
  !
  !.. Local Scalars ..
  real(lp) :: alpha,beta
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  alpha = zero
  beta  = zero
  !
  call eval_jacobi_poly(n,alpha,beta,x,y,dy,d2y)
  !
end subroutine eval_legendre_poly_QP
!
!###############################################################################
!
pure subroutine eval_chebyshev_poly_DP(n,x,y,dy,d2y)
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp),  intent(in) :: x
  real(lp), intent(out) :: y
  real(lp), intent(out) :: dy
  real(lp), intent(out) :: d2y
  !
  !.. Local Scalars ..
  integer  :: k
  real(lp) :: yp,ym,dyp,dym,d2yp,d2ym
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: two  = 2.0_lp
  real(lp), parameter :: four = 4.0_lp
  !
continue
  !
  y   = one
  dy  = zero
  d2y = zero
  !
  if (n > 0) then
    y   = x
    dy  = one
    d2y = zero
  end if
  !
  if (n > 1) then
    !
    yp   = one
    dyp  = zero
    d2yp = zero
    !
    do k = 2,n
      !
      ym = y
      y  = two*x*y - yp
      yp = ym
      !
      dym = dy
      dy  = two*x*dy + two*yp - dyp
      dyp = dym
      !
      d2ym = d2y
      d2y  = two*x*d2y + four*dyp - d2yp
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
pure subroutine eval_chebyshev_poly_QP(n,x,y,dy,d2y)
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp),  intent(in) :: x
  real(lp), intent(out) :: y
  real(lp), intent(out) :: dy
  real(lp), intent(out) :: d2y
  !
  !.. Local Scalars ..
  integer  :: k
  real(lp) :: yp,ym,dyp,dym,d2yp,d2ym
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: two  = 2.0_lp
  real(lp), parameter :: four = 4.0_lp
  !
continue
  !
  y   = one
  dy  = zero
  d2y = zero
  !
  if (n > 0) then
    y   = x
    dy  = one
    d2y = zero
  end if
  !
  if (n > 1) then
    !
    yp   = one
    dyp  = zero
    d2yp = zero
    !
    do k = 2,n
      !
      ym = y
      y  = two*x*y - yp
      yp = ym
      !
      dym = dy
      dy  = two*x*dy + two*yp - dyp
      dyp = dym
      !
      d2ym = d2y
      d2y  = two*x*d2y + four*dyp - d2yp
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
pure subroutine eval_normalized_jacobi_poly_DP(n,alpha,beta,x,y,dy)
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp),  intent(in) :: alpha
  real(lp),  intent(in) :: beta
  real(lp),  intent(in) :: x
  real(lp), intent(out) :: y
  real(lp), intent(out) :: dy
  !
  !.. Local Scalars ..
  integer  :: n1
  real(lp) :: a1,b1,ab1,rn
  !
  !.. Local Constants ..
  real(lp), parameter :: one = 1.0_lp
  !
continue
  !
  a1  = alpha + one
  b1  = beta  + one
  ab1 = alpha + beta + one
  rn  = real(n,kind=lp)
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
pure subroutine eval_normalized_jacobi_poly_QP(n,alpha,beta,x,y,dy)
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp),  intent(in) :: alpha
  real(lp),  intent(in) :: beta
  real(lp),  intent(in) :: x
  real(lp), intent(out) :: y
  real(lp), intent(out) :: dy
  !
  !.. Local Scalars ..
  integer  :: n1
  real(lp) :: a1,b1,ab1,rn
  !
  !.. Local Constants ..
  real(lp), parameter :: one = 1.0_lp
  !
continue
  !
  a1  = alpha + one
  b1  = beta  + one
  ab1 = alpha + beta + one
  rn  = real(n,kind=lp)
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
pure function normalized_jacobi_poly_DP(N,alpha,beta,x) result(return_value)
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: N
  real(lp), intent(in) :: alpha
  real(lp), intent(in) :: beta
  real(lp), intent(in) :: x
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
  !.. Local Scalars ..
  integer  :: i
  real(lp) :: ab,ab1,ab2,ab3,amb,a1,b1
  real(lp) :: gamma0,gamma1
  real(lp) :: aold,anew,bnew,h1,ri
  real(lp) :: sqrt_numer,sqrt_denom
  !
  !.. Local Arrays ..
  real(lp) :: PL(N+1)
  !
  !.. Local Constants ..
  real(lp), parameter :: half  = 0.5_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  real(lp), parameter :: three = 3.0_lp
  !
continue
  !
  ab = alpha + beta
  amb = alpha - beta
  a1 = alpha + one
  b1 = beta + one
  ab1 = ab + one
  ab2 = ab + two
  ab3 = ab + three
  !
  gamma0 = two**ab1 / ab1 * gammafun(a1) * gammafun(b1) / gammafun(ab1)
  !
  PL(1) = one/sqrt(gamma0)
  !
  if (N == 0) then
    !
    return_value = PL(1)
    !
  else
    !
    gamma1 = a1 * b1 / ab3 * gamma0
    !
    PL(2) = half*( ab2*x + amb ) / sqrt(gamma1)
    !
    aold = two / ab2 * sqrt( a1*b1/ab3 )
    !
    do i = 1,N-1
      ri = real(i,kind=lp)
      h1 = two*ri + ab
      sqrt_numer = (ri+one)*(ri+ab1)*(ri+a1)*(ri+b1)
      sqrt_denom = (h1+one)*(h1+three)
      anew = sqrt( sqrt_numer / sqrt_denom ) * two / (h1+two)
      bnew = (beta*beta - alpha*alpha) / h1 / (h1+two)
      PL(i+2) = one/anew*( -aold*PL(i) + (x-bnew)*PL(i+1) )
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
pure function normalized_jacobi_poly_QP(N,alpha,beta,x) result(return_value)
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: N
  real(lp), intent(in) :: alpha
  real(lp), intent(in) :: beta
  real(lp), intent(in) :: x
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
  !.. Local Scalars ..
  integer  :: i
  real(lp) :: ab,ab1,ab2,ab3,amb,a1,b1
  real(lp) :: gamma0,gamma1
  real(lp) :: aold,anew,bnew,h1,ri
  real(lp) :: sqrt_numer,sqrt_denom
  !
  !.. Local Arrays ..
  real(lp) :: PL(N+1)
  !
  !.. Local Constants ..
  real(lp), parameter :: half  = 0.5_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  real(lp), parameter :: three = 3.0_lp
  !
continue
  !
  ab = alpha + beta
  amb = alpha - beta
  a1 = alpha + one
  b1 = beta + one
  ab1 = ab + one
  ab2 = ab + two
  ab3 = ab + three
  !
  gamma0 = two**ab1 / ab1 * gammafun(a1) * gammafun(b1) / gammafun(ab1)
  !
  PL(1) = one/sqrt(gamma0)
  !
  if (N == 0) then
    !
    return_value = PL(1)
    !
  else
    !
    gamma1 = a1 * b1 / ab3 * gamma0
    !
    PL(2) = half*( ab2*x + amb ) / sqrt(gamma1)
    !
    aold = two / ab2 * sqrt( a1*b1/ab3 )
    !
    do i = 1,N-1
      ri = real(i,kind=lp)
      h1 = two*ri + ab
      sqrt_numer = (ri+one)*(ri+ab1)*(ri+a1)*(ri+b1)
      sqrt_denom = (h1+one)*(h1+three)
      anew = sqrt( sqrt_numer / sqrt_denom ) * two / (h1+two)
      bnew = (beta*beta - alpha*alpha) / h1 / (h1+two)
      PL(i+2) = one/anew*( -aold*PL(i) + (x-bnew)*PL(i+1) )
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
pure function grad_normalized_jacobi_poly_DP(n,alpha,beta,x) &
                                      result(return_value)
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: n
  real(lp), intent(in) :: alpha
  real(lp), intent(in) :: beta
  real(lp), intent(in) :: x
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
  !.. Local Scalars ..
  real(lp) :: rn,a1,b1,ab1
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  if (n == 0) then
    !
    return_value = zero
    !
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
pure function grad_normalized_jacobi_poly_QP(n,alpha,beta,x) &
                                      result(return_value)
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: n
  real(lp), intent(in) :: alpha
  real(lp), intent(in) :: beta
  real(lp), intent(in) :: x
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
  !.. Local Scalars ..
  real(lp) :: rn,a1,b1,ab1
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  if (n == 0) then
    !
    return_value = zero
    !
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
pure subroutine nodes_jacobi_gauss_DP(n,alpha,beta,xi,dxi)
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                   intent(in) :: n
  real(lp),                  intent(in) :: alpha
  real(lp),                  intent(in) :: beta
  real(lp), dimension(1:n), intent(out) :: xi
  real(lp), dimension(1:n), intent(out) :: dxi
  !
  !.. Local Scalars ..
  integer  :: i,it
  real(lp) :: const,ri,x,y,dy,d2y
  !
  !.. Local Constants ..
  real(lp), parameter :: zero  = 0.0_lp
  real(lp), parameter :: half  = 0.5_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  real(lp), parameter :: four  = 4.0_lp
  real(lp), parameter :: eps17 = 1.0e-17_lp
  real(lp), parameter :: pi    = 3.14159265358979323846_lp
  !
continue
  !
  if (n > 0) then
    xi(1)  = (beta - alpha) / (alpha + beta + two)
    dxi(1) = half*(alpha + beta) + one
  end if
  !
  if (n > 1) then
    !
    const = half*pi / (two*real(n,kind=lp) + alpha + beta + one)
    !
    do i = 1,n
      !
      ri = real(i,kind=lp)
      x  = -cos( const * (four*ri + alpha + beta - one) )
      !
      do it = 1,8
        call eval_jacobi_poly(n,alpha,beta,x,y,dy,d2y)
        if (abs(y) <= eps17) exit
        x = x - y/dy
      end do
      !
      if (abs(x) <= eps17) x = zero
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
pure subroutine nodes_jacobi_gauss_QP(n,alpha,beta,xi,dxi)
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                   intent(in) :: n
  real(lp),                  intent(in) :: alpha
  real(lp),                  intent(in) :: beta
  real(lp), dimension(1:n), intent(out) :: xi
  real(lp), dimension(1:n), intent(out) :: dxi
  !
  !.. Local Scalars ..
  integer  :: i,it
  real(lp) :: const,ri,x,y,dy,d2y
  !
  !.. Local Constants ..
  real(lp), parameter :: zero  = 0.0_lp
  real(lp), parameter :: half  = 0.5_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  real(lp), parameter :: four  = 4.0_lp
  real(lp), parameter :: eps17 = 1.0e-17_lp
  real(lp), parameter :: pi    = 3.14159265358979323846_lp
  !
continue
  !
  if (n > 0) then
    xi(1)  = (beta - alpha) / (alpha + beta + two)
    dxi(1) = half*(alpha + beta) + one
  end if
  !
  if (n > 1) then
    !
    const = half*pi / (two*real(n,kind=lp) + alpha + beta + one)
    !
    do i = 1,n
      !
      ri = real(i,kind=lp)
      x  = -cos( const * (four*ri + alpha + beta - one) )
      !
      do it = 1,8
        call eval_jacobi_poly(n,alpha,beta,x,y,dy,d2y)
        if (abs(y) <= eps17) exit
        x = x - y/dy
      end do
      !
      if (abs(x) <= eps17) x = zero
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
pure subroutine nodes_legendre_gauss_DP(n,xi,dxi)
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                   intent(in) :: n
  real(lp), dimension(1:n), intent(out) :: xi
  real(lp), dimension(1:n), intent(out) :: dxi
  !
  !.. Local Scalars ..
  real(lp) :: alpha,beta
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  call nodes_jacobi_gauss(n,zero,zero,xi,dxi)
  !
end subroutine nodes_legendre_gauss_DP
!
!###############################################################################
!
pure subroutine nodes_legendre_gauss_QP(n,xi,dxi)
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                   intent(in) :: n
  real(lp), dimension(1:n), intent(out) :: xi
  real(lp), dimension(1:n), intent(out) :: dxi
  !
  !.. Local Scalars ..
  real(lp) :: alpha,beta
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  call nodes_jacobi_gauss(n,zero,zero,xi,dxi)
  !
end subroutine nodes_legendre_gauss_QP
!
!###############################################################################
!
pure subroutine nodes_right_radau_DP(n,xi)
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                       intent(in) :: n
  real(lp), dimension(1:n+1), intent(inout) :: xi
  !
  !.. Local Scalars ..
  integer  :: i,it
  real(lp) :: x,y,dy
  !
  !.. Local Arrays ..
  real(lp), dimension(0:n+1) :: cfrr
  real(lp), dimension(0:n) :: cdrr
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: pi   = 3.14159265358979323846_lp
  !
continue
  !
  ! Set the final zero which will always be unity
  !
  xi(n+1) = one
  !
  if (n > 0) then
    !
    ! Get the coefficients of the right Radau polynomial
    !
    cfrr(:) = coefrr(n+1,zero)
    !
    do i = 0,n
      cdrr(i) = real(i+1,kind=lp) * cfrr(i+1)
    end do
    !
    xi(n+1) = one
    !
    do i = 1,n
      !
      x  = -cos( pi * real(i,kind=lp) / real(n+1,kind=lp) )
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
pure subroutine nodes_right_radau_QP(n,xi)
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                       intent(in) :: n
  real(lp), dimension(1:n+1), intent(inout) :: xi
  !
  !.. Local Scalars ..
  integer  :: i,it
  real(lp) :: x,y,dy
  !
  !.. Local Arrays ..
  real(lp), dimension(0:n+1) :: cfrr
  real(lp), dimension(0:n) :: cdrr
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: pi   = 3.14159265358979323846_lp
  !
continue
  !
  ! Set the final zero which will always be unity
  !
  xi(n+1) = one
  !
  if (n > 0) then
    !
    ! Get the coefficients of the right Radau polynomial
    !
    cfrr(:) = coefrr(n+1,zero)
    !
    do i = 0,n
      cdrr(i) = real(i+1,kind=lp) * cfrr(i+1)
    end do
    !
    xi(n+1) = one
    !
    do i = 1,n
      !
      x  = -cos( pi * real(i,kind=lp) / real(n+1,kind=lp) )
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
pure subroutine nodes_chebyshev_gauss_DP(n,xi)
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                   intent(in) :: n
  real(lp), dimension(1:n), intent(out) :: xi
  !
  !.. Local Scalars ..
  real(lp) :: alpha,beta
  !
  !.. Local Arrays ..
  real(lp), dimension(1:n) :: dxi
  !
  !.. Local Constants ..
  real(lp), parameter :: half = 0.5_lp
  !
continue
  !
  alpha = -half
  beta  = -half
  !
  call nodes_jacobi_gauss(n,alpha,beta,xi,dxi)
  !
end subroutine nodes_chebyshev_gauss_DP
!
!###############################################################################
!
pure subroutine nodes_chebyshev_gauss_QP(n,xi)
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                   intent(in) :: n
  real(lp), dimension(1:n), intent(out) :: xi
  !
  !.. Local Scalars ..
  real(lp) :: alpha,beta
  !
  !.. Local Arrays ..
  real(lp), dimension(1:n) :: dxi
  !
  !.. Local Constants ..
  real(lp), parameter :: half = 0.5_lp
  !
continue
  !
  alpha = -half
  beta  = -half
  !
  call nodes_jacobi_gauss(n,alpha,beta,xi,dxi)
  !
end subroutine nodes_chebyshev_gauss_QP
!
!###############################################################################
!
pure subroutine nodes_jacobi_gauss_lobatto_DP(n,alpha,beta,xi,vn)
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: n
  real(lp),                    intent(in) :: alpha
  real(lp),                    intent(in) :: beta
  real(lp), dimension(1:n+1), intent(out) :: xi
  real(lp), dimension(1:n+1), intent(out) :: vn
  !
  !.. Local Scalars ..
  integer  :: i,it,n1
  real(lp) :: x,y,dy,d2y,ri,rn,ab,a1,b1,const
  !
  !.. Local Constants ..
  real(lp), parameter :: zero  = 0.0_lp
  real(lp), parameter :: half  = 0.5_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  real(lp), parameter :: four  = 4.0_lp
  real(lp), parameter :: eps17 = 1.0e-17_lp
  real(lp), parameter :: pi    = 3.14159265358979323846_lp
  !
continue
  !
  if (n == 0) then
    !
    x = zero
    call eval_jacobi_poly(n,alpha,beta,x,y,dy,d2y)
    xi(1) = x
    vn(1) = y
    !
  else
    !
    x = -one
    call eval_jacobi_poly(n,alpha,beta,x,y,dy,d2y)
    xi(1) = x
    vn(1) = y
    !
    x = one
    call eval_jacobi_poly(n,alpha,beta,x,y,dy,d2y)
    xi(n+1) = x
    vn(n+1) = y
    !
  end if
  !
  if (n > 1) then
    !
    rn = real(n,kind=lp)
    ab = alpha + beta
    a1 = alpha + one
    b1 = beta  + one
    n1 = n-1
    const = half * pi / (two*rn + ab + one)
    !
    do i = 2,n
      !
      ri = real(i-1,kind=lp)
      x  = -cos( const * (four*ri + ab + one) )
      !
      do it = 1,8
        call eval_jacobi_poly(n1,a1,b1,x,y,dy,d2y)
        if (abs(y) <= eps17) exit
        x = x - y/dy
      end do
      !
      if (abs(x) <= eps17) x = zero
      !
      xi(i) = x
      vn(i) = -half * dy * (one - x*x) / rn
      !
    end do
    !
  end if
  !
end subroutine nodes_jacobi_gauss_lobatto_DP
!
!###############################################################################
!
pure subroutine nodes_jacobi_gauss_lobatto_QP(n,alpha,beta,xi,vn)
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: n
  real(lp),                    intent(in) :: alpha
  real(lp),                    intent(in) :: beta
  real(lp), dimension(1:n+1), intent(out) :: xi
  real(lp), dimension(1:n+1), intent(out) :: vn
  !
  !.. Local Scalars ..
  integer  :: i,it,n1
  real(lp) :: x,y,dy,d2y,ri,rn,ab,a1,b1,const
  !
  !.. Local Constants ..
  real(lp), parameter :: zero  = 0.0_lp
  real(lp), parameter :: half  = 0.5_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  real(lp), parameter :: four  = 4.0_lp
  real(lp), parameter :: eps17 = 1.0e-17_lp
  real(lp), parameter :: pi    = 3.14159265358979323846_lp
  !
continue
  !
  if (n == 0) then
    !
    x = zero
    call eval_jacobi_poly(n,alpha,beta,x,y,dy,d2y)
    xi(1) = x
    vn(1) = y
    !
  else
    !
    x = -one
    call eval_jacobi_poly(n,alpha,beta,x,y,dy,d2y)
    xi(1) = x
    vn(1) = y
    !
    x = one
    call eval_jacobi_poly(n,alpha,beta,x,y,dy,d2y)
    xi(n+1) = x
    vn(n+1) = y
    !
  end if
  !
  if (n > 1) then
    !
    rn = real(n,kind=lp)
    ab = alpha + beta
    a1 = alpha + one
    b1 = beta  + one
    n1 = n-1
    const = half * pi / (two*rn + ab + one)
    !
    do i = 2,n
      !
      ri = real(i-1,kind=lp)
      x  = -cos( const * (four*ri + ab + one) )
      !
      do it = 1,8
        call eval_jacobi_poly(n1,a1,b1,x,y,dy,d2y)
        if (abs(y) <= eps17) exit
        x = x - y/dy
      end do
      !
      if (abs(x) <= eps17) x = zero
      !
      xi(i) = x
      vn(i) = -half * dy * (one - x*x) / rn
      !
    end do
    !
  end if
  !
end subroutine nodes_jacobi_gauss_lobatto_QP
!
!###############################################################################
!
pure subroutine nodes_legendre_gauss_lobatto_DP(n,xi,vn)
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: n
  real(lp), dimension(1:n+1), intent(out) :: xi
  real(lp), dimension(1:n+1), intent(out) :: vn
  !
  !.. Local Scalars ..
  integer  :: i,it,n2
  real(lp) :: sn,x,y,dy,d2y
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: pi   = 3.14159265358979323846_lp
  !
continue
  !
  if (n == 0) then
    xi(1) = zero
    vn(1) = one
  else
    !
    n2 = (n-1)/2
    sn = real(2*n - 4*n2 - 3,kind=lp)
    !
    xi(1) = -one
    vn(1) =  sn
    !
    xi(n+1) = one
    vn(n+1) = one
    !
  end if
  !
  if (n > 1) then
    x = zero
    call eval_legendre_poly(n,x,y,dy,d2y)
    xi(n2+2) = x
    vn(n2+2) = y
  end if
  !
  if (n > 2) then
    !
    do i = 1,n2
      !
      x = cos( pi * real(i,kind=lp) / real(n,kind=lp) )
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
pure subroutine nodes_legendre_gauss_lobatto_QP(n,xi,vn)
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: n
  real(lp), dimension(1:n+1), intent(out) :: xi
  real(lp), dimension(1:n+1), intent(out) :: vn
  !
  !.. Local Scalars ..
  integer  :: i,it,n2
  real(lp) :: sn,x,y,dy,d2y
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: pi   = 3.14159265358979323846_lp
  !
continue
  !
  if (n == 0) then
    xi(1) = zero
    vn(1) = one
  else
    !
    n2 = (n-1)/2
    sn = real(2*n - 4*n2 - 3,kind=lp)
    !
    xi(1) = -one
    vn(1) =  sn
    !
    xi(n+1) = one
    vn(n+1) = one
    !
  end if
  !
  if (n > 1) then
    x = zero
    call eval_legendre_poly(n,x,y,dy,d2y)
    xi(n2+2) = x
    vn(n2+2) = y
  end if
  !
  if (n > 2) then
    !
    do i = 1,n2
      !
      x = cos( pi * real(i,kind=lp) / real(n,kind=lp) )
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
pure subroutine nodes_chebyshev_gauss_lobatto_DP(n,xi)
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: n
  real(lp), dimension(1:n+1), intent(out) :: xi
  !
  !.. Local Scalars ..
  integer :: i,n2
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: pi   = 3.14159265358979323846_lp
  !
continue
  !
  if (n == 0) then
    xi(1) = zero
  else
    xi(1)   = -one
    xi(n+1) =  one
  end if
  !
  if (n > 1) then
    n2 = (n-1)/2
    xi(n2+2) = zero
  end if
  !
  if (n > 2) then
    do i = 1,n2
      xi(i+1)   = -cos( pi * real(i,kind=lp) / real(n,kind=lp) )
      xi(n-i+1) =  cos( pi * real(i,kind=lp) / real(n,kind=lp) )
    end do
  end if
  !
end subroutine nodes_chebyshev_gauss_lobatto_DP
!
!###############################################################################
!
pure subroutine nodes_chebyshev_gauss_lobatto_QP(n,xi)
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: n
  real(lp), dimension(1:n+1), intent(out) :: xi
  !
  !.. Local Scalars ..
  integer :: i,n2
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: pi   = 3.14159265358979323846_lp
  !
continue
  !
  if (n == 0) then
    xi(1) = zero
  else
    xi(1)   = -one
    xi(n+1) =  one
  end if
  !
  if (n > 1) then
    n2 = (n-1)/2
    xi(n2+2) = zero
  end if
  !
  if (n > 2) then
    do i = 1,n2
      xi(i+1)   = -cos( pi * real(i,kind=lp) / real(n,kind=lp) )
      xi(n-i+1) =  cos( pi * real(i,kind=lp) / real(n,kind=lp) )
    end do
  end if
  !
end subroutine nodes_chebyshev_gauss_lobatto_QP
!
!###############################################################################
!
pure subroutine weights_jacobi_gauss_DP(n,alpha,beta,xi,dy,weights)
  !
  !***************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE JACOBI GAUSS FORMULA
  !   N  = ORDER OF THE FORMULA
  !   alpha  = PARAMETER > -1
  !   beta  = PARAMETER > -1
  !   CS = ZEROES OF THE JACOBI POLYNOMIAL, CS(I), I=1,N
  !   DZ = VECTOR OF THE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
  !   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
  !***************************************************************
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                   intent(in) :: n
  real(lp),                  intent(in) :: alpha
  real(lp),                  intent(in) :: beta
  real(lp), dimension(1:n),  intent(in) :: xi
  real(lp), dimension(1:n),  intent(in) :: dy
  real(lp), dimension(1:n), intent(out) :: weights
  !
  !.. Local Scalars ..
  integer  :: i,m
  real(lp) :: ab,a2,b2,c,ri
  !
  !.. Local Constants ..
  real(lp), parameter :: half = 0.5_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: two  = 2.0_lp
  !
continue
  !
  weights(:) = one
  !
  ab = alpha + beta + two
  a2 = alpha + two
  b2 = beta  + two
  !
  c = half * (two**ab) * gammafun(a2) * gammafun(b2) / gammafun(ab)
  !
  do i = 2,n
    ri = real(i,kind=lp)
    c  = c * (ri+alpha) * (ri+beta) / (ri*(ri+alpha+beta))
  end do
  !
  do i = 1,n
    weights(i) = c / ((one - xi(i)*xi(i)) * dy(i) * dy(i))
  end do
  !
end subroutine weights_jacobi_gauss_DP
!
!###############################################################################
!
pure subroutine weights_jacobi_gauss_QP(n,alpha,beta,xi,dy,weights)
  !
  !***************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE JACOBI GAUSS FORMULA
  !   N  = ORDER OF THE FORMULA
  !   alpha  = PARAMETER > -1
  !   beta  = PARAMETER > -1
  !   CS = ZEROES OF THE JACOBI POLYNOMIAL, CS(I), I=1,N
  !   DZ = VECTOR OF THE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
  !   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
  !***************************************************************
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                   intent(in) :: n
  real(lp),                  intent(in) :: alpha
  real(lp),                  intent(in) :: beta
  real(lp), dimension(1:n),  intent(in) :: xi
  real(lp), dimension(1:n),  intent(in) :: dy
  real(lp), dimension(1:n), intent(out) :: weights
  !
  !.. Local Scalars ..
  integer  :: i,m
  real(lp) :: ab,a2,b2,c,ri
  !
  !.. Local Constants ..
  real(lp), parameter :: half = 0.5_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: two  = 2.0_lp
  !
continue
  !
  weights(:) = one
  !
  ab = alpha + beta + two
  a2 = alpha + two
  b2 = beta  + two
  !
  c = half * (two**ab) * gammafun(a2) * gammafun(b2) / gammafun(ab)
  !
  do i = 2,n
    ri = real(i,kind=lp)
    c  = c * (ri+alpha) * (ri+beta) / (ri*(ri+alpha+beta))
  end do
  !
  do i = 1,n
    weights(i) = c / ((one - xi(i)*xi(i)) * dy(i) * dy(i))
  end do
  !
end subroutine weights_jacobi_gauss_QP
!
!###############################################################################
!
pure subroutine weights_legendre_gauss_DP(n,xi,dy,weights)
  !
  !****************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS FORMULA
  !   N  = ORDER OF THE FORMULA
  !   CS = ZEROES OF THE LEGENDRE POLYNOMIAL, CS(I), I=1,N
  !   DZ = VECTOR OF THE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
  !   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
  !****************************************************************
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                   intent(in) :: n
  real(lp), dimension(1:n),  intent(in) :: xi
  real(lp), dimension(1:n),  intent(in) :: dy
  real(lp), dimension(1:n), intent(out) :: weights
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Constants ..
  real(lp), parameter :: one = 1.0_lp
  real(lp), parameter :: two = 2.0_lp
  !
continue
  !
  weights(:) = two
  !
  do i = 1,n
    weights(i) = two / ((one - xi(i)*xi(i)) * dy(i) * dy(i))
  end do
  !
end subroutine weights_legendre_gauss_DP
!
!###############################################################################
!
pure subroutine weights_legendre_gauss_QP(n,xi,dy,weights)
  !
  !****************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS FORMULA
  !   N  = ORDER OF THE FORMULA
  !   CS = ZEROES OF THE LEGENDRE POLYNOMIAL, CS(I), I=1,N
  !   DZ = VECTOR OF THE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
  !   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
  !****************************************************************
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                   intent(in) :: n
  real(lp), dimension(1:n),  intent(in) :: xi
  real(lp), dimension(1:n),  intent(in) :: dy
  real(lp), dimension(1:n), intent(out) :: weights
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Constants ..
  real(lp), parameter :: one = 1.0_lp
  real(lp), parameter :: two = 2.0_lp
  !
continue
  !
  weights(:) = two
  !
  do i = 1,n
    weights(i) = two / ((one - xi(i)*xi(i)) * dy(i) * dy(i))
  end do
  !
end subroutine weights_legendre_gauss_QP
!
!###############################################################################
!
pure subroutine weights_right_radau_DP(n,xi,weights)
  !
  !**********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   XI = RIGHT RADAU NODES, XI(I), I=0,N
  !   VN = VALUES OF THE RIGHT RADAU POLYNOMIAL AT THE NODES, VN(I), I=0,N
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !**********************************************************************
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: n
  real(lp), dimension(1:n+1),  intent(in) :: xi
  real(lp), dimension(1:n+1), intent(out) :: weights
  !
  !.. Local Scalars ..
  integer  :: i
  real(lp) :: c,y,dy,d2y
  !
  !.. Local Constants ..
  real(lp), parameter :: one = 1.0_lp
  real(lp), parameter :: two = 2.0_lp
  !
continue
  !
  c = two
  !
  if (n > 0) then
    !
    c = one / (real(n+1,kind=lp)*real(n+1,kind=lp))
    !
    do i = 1,n+1
      call eval_legendre_poly(n,xi(i),y,dy,d2y)
      weights(i) = c * (one+xi(i)) / (y*y)
    end do
    !
  end if
  !
end subroutine weights_right_radau_DP
!
!###############################################################################
!
pure subroutine weights_right_radau_QP(n,xi,weights)
  !
  !**********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   XI = RIGHT RADAU NODES, XI(I), I=0,N
  !   VN = VALUES OF THE RIGHT RADAU POLYNOMIAL AT THE NODES, VN(I), I=0,N
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !**********************************************************************
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: n
  real(lp), dimension(1:n+1),  intent(in) :: xi
  real(lp), dimension(1:n+1), intent(out) :: weights
  !
  !.. Local Scalars ..
  integer  :: i
  real(lp) :: c,y,dy,d2y
  !
  !.. Local Constants ..
  real(lp), parameter :: one = 1.0_lp
  real(lp), parameter :: two = 2.0_lp
  !
continue
  !
  c = two
  !
  if (n > 0) then
    !
    c = one / (real(n+1,kind=lp)*real(n+1,kind=lp))
    !
    do i = 1,n+1
      call eval_legendre_poly(n,xi(i),y,dy,d2y)
      weights(i) = c * (one+xi(i)) / (y*y)
    end do
    !
  end if
  !
end subroutine weights_right_radau_QP
!
!###############################################################################
!
pure subroutine weights_chebyshev_gauss_DP(n,weights)
  !
  !****************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE CHEBYSHEV GAUSS FORMULA
  !   N  = ORDER OF THE FORMULA
  !   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
  !****************************************************************
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                   intent(in) :: n
  real(lp), dimension(1:n), intent(out) :: weights
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Constants ..
  real(lp), parameter :: pi = 3.14159265358979323846_lp
  !
continue
  !
  weights(:) = pi
  !
  do i = 1,n
    weights(i) = pi/real(n,kind=lp)
  end do
  !
end subroutine weights_chebyshev_gauss_DP
!
!###############################################################################
!
pure subroutine weights_chebyshev_gauss_QP(n,weights)
  !
  !****************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE CHEBYSHEV GAUSS FORMULA
  !   N  = ORDER OF THE FORMULA
  !   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
  !****************************************************************
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                   intent(in) :: n
  real(lp), dimension(1:n), intent(out) :: weights
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Constants ..
  real(lp), parameter :: pi = 3.14159265358979323846_lp
  !
continue
  !
  weights(:) = pi
  !
  do i = 1,n
    weights(i) = pi/real(n,kind=lp)
  end do
  !
end subroutine weights_chebyshev_gauss_QP
!
!###############################################################################
!
pure subroutine weights_jacobi_gauss_lobatto_DP(n,alpha,beta,xi,weights)
  !
  !********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE JACOBI GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   alpha  = PARAMETER > -1
  !   beta  = PARAMETER > -1
  !   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !********************************************************************
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: n
  real(lp),                    intent(in) :: alpha
  real(lp),                    intent(in) :: beta
  real(lp), dimension(1:n+1),  intent(in) :: xi
  real(lp), dimension(1:n+1), intent(out) :: weights
  !
  !.. Local Scalars ..
  integer  :: i
  real(lp) :: a1,a2,b1,b2,ab,ab1,ab2
  real(lp) :: c,c1,c2,c3,ri,rn,su
  real(lp) :: y,z,dy,d2y
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: half = 0.5_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: two  = 2.0_lp
  !
continue
  !
  a1  = alpha + one
  a2  = alpha + two
  b1  = beta  + one
  b2  = beta  + two
  ab  = alpha + beta
  ab1 = alpha + beta + one
  ab2 = alpha + beta + two
  !
  c = (two**ab) * gammafun(a1) * gammafun(b1) / gammafun(ab2)
  !
  weights(:) = two*c
  !
  if (n == 1) then
    weights(1)   = two * c * a1 / ab2
    weights(n+1) = two * c * b1 / ab2
  end if
  !
  if (n > 1) then
    !
    rn = real(n,kind=lp)
    c  = c * (two*rn + ab) / (rn + ab1)
    c1 = c * a1 / (b2 * ab2)
    c2 = c * b1 / (a2 * ab2)
    c3 = half * c * a1 * b1
    !
    do i = 2,n-1
      ri = real(i-1,kind=lp)
      c1 = c1*(ri+a1)*ri / ((ri+ab2)*(ri+b2))
      c2 = c2*(ri+b1)*ri / ((ri+ab2)*(ri+a2))
      c3 = c3*(ri+a1)*(ri+b1) / ((ri+two)*(ri+ab1))
    end do
    !
    su = zero
    do i = 2,n
      su = su + xi(i)
    end do
    !
    weights(1)   = c1*(rn-one-su)
    weights(n+1) = c2*(rn-one+su)
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
!###############################################################################
!
pure subroutine weights_jacobi_gauss_lobatto_QP(n,alpha,beta,xi,weights)
  !
  !********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE JACOBI GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   alpha  = PARAMETER > -1
  !   beta  = PARAMETER > -1
  !   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !********************************************************************
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: n
  real(lp),                    intent(in) :: alpha
  real(lp),                    intent(in) :: beta
  real(lp), dimension(1:n+1),  intent(in) :: xi
  real(lp), dimension(1:n+1), intent(out) :: weights
  !
  !.. Local Scalars ..
  integer  :: i
  real(lp) :: a1,a2,b1,b2,ab,ab1,ab2
  real(lp) :: c,c1,c2,c3,ri,rn,su
  real(lp) :: y,z,dy,d2y
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: half = 0.5_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: two  = 2.0_lp
  !
continue
  !
  a1  = alpha + one
  a2  = alpha + two
  b1  = beta  + one
  b2  = beta  + two
  ab  = alpha + beta
  ab1 = alpha + beta + one
  ab2 = alpha + beta + two
  !
  c = (two**ab) * gammafun(a1) * gammafun(b1) / gammafun(ab2)
  !
  weights(:) = two*c
  !
  if (n == 1) then
    weights(1)   = two * c * a1 / ab2
    weights(n+1) = two * c * b1 / ab2
  end if
  !
  if (n > 1) then
    !
    rn = real(n,kind=lp)
    c  = c * (two*rn + ab) / (rn + ab1)
    c1 = c * a1 / (b2 * ab2)
    c2 = c * b1 / (a2 * ab2)
    c3 = half * c * a1 * b1
    !
    do i = 2,n-1
      ri = real(i-1,kind=lp)
      c1 = c1*(ri+a1)*ri / ((ri+ab2)*(ri+b2))
      c2 = c2*(ri+b1)*ri / ((ri+ab2)*(ri+a2))
      c3 = c3*(ri+a1)*(ri+b1) / ((ri+two)*(ri+ab1))
    end do
    !
    su = zero
    do i = 2,n
      su = su + xi(i)
    end do
    !
    weights(1)   = c1*(rn-one-su)
    weights(n+1) = c2*(rn-one+su)
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
!###############################################################################
!
pure subroutine weights_legendre_gauss_lobatto_DP(n,vn,weights)
  !
  !**********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N
  !   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !**********************************************************************
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: n
  real(lp), dimension(1:n+1),  intent(in) :: vn
  real(lp), dimension(1:n+1), intent(out) :: weights
  !
  !.. Local Scalars ..
  integer  :: i
  real(lp) :: c
  !
  !.. Local Constants ..
  real(lp), parameter :: one = 1.0_lp
  real(lp), parameter :: two = 2.0_lp
  !
continue
  !
  c = two
  if (n > 0) then
    c = c / (real(n,kind=lp)*(real(n,kind=lp)+one))
  end if
  !
  do i = 1,n+1
    weights(i) = c / (vn(i)*vn(i))
  end do
  !
end subroutine weights_legendre_gauss_lobatto_DP
!
!###############################################################################
!
pure subroutine weights_legendre_gauss_lobatto_QP(n,vn,weights)
  !
  !**********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N
  !   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !**********************************************************************
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: n
  real(lp), dimension(1:n+1),  intent(in) :: vn
  real(lp), dimension(1:n+1), intent(out) :: weights
  !
  !.. Local Scalars ..
  integer  :: i
  real(lp) :: c
  !
  !.. Local Constants ..
  real(lp), parameter :: one = 1.0_lp
  real(lp), parameter :: two = 2.0_lp
  !
continue
  !
  c = two
  if (n > 0) then
    c = c / (real(n,kind=lp)*(real(n,kind=lp)+one))
  end if
  !
  do i = 1,n+1
    weights(i) = c / (vn(i)*vn(i))
  end do
  !
end subroutine weights_legendre_gauss_lobatto_QP
!
!###############################################################################
!
pure subroutine weights_chebyshev_gauss_lobatto_DP(n,weights)
  !
  !***********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE CHEBYSHEV GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !***********************************************************************
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: n
  real(lp), dimension(1:n+1), intent(out) :: weights
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Constants ..
  real(lp), parameter :: half = 0.5_lp
  real(lp), parameter :: pi   = 3.14159265358979323846_lp
  !
continue
  !
  weights(:) = pi
  !
  if (n > 0) then
    weights(1)   = half*pi/real(n,kind=lp)
    weights(n+1) = half*pi/real(n,kind=lp)
  end if
  !
  do i = 2,n
    weights(i) = pi/real(n,kind=lp)
  end do
  !
end subroutine weights_chebyshev_gauss_lobatto_DP
!
!###############################################################################
!
pure subroutine weights_chebyshev_gauss_lobatto_QP(n,weights)
  !
  !***********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE CHEBYSHEV GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !***********************************************************************
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: n
  real(lp), dimension(1:n+1), intent(out) :: weights
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Constants ..
  real(lp), parameter :: half = 0.5_lp
  real(lp), parameter :: pi   = 3.14159265358979323846_lp
  !
continue
  !
  weights(:) = pi
  !
  if (n > 0) then
    weights(1)   = half*pi/real(n,kind=lp)
    weights(n+1) = half*pi/real(n,kind=lp)
  end if
  !
  do i = 2,n
    weights(i) = pi/real(n,kind=lp)
  end do
  !
end subroutine weights_chebyshev_gauss_lobatto_QP
!
!###############################################################################
!
pure subroutine weights_dgbook_jacobi_gauss_lobatto_DP(n,alpha,beta,xi,weights)
  !
  !********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE JACOBI GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   alpha  = PARAMETER > -1
  !   beta  = PARAMETER > -1
  !   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !********************************************************************
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: n
  real(lp),                    intent(in) :: alpha
  real(lp),                    intent(in) :: beta
  real(lp), dimension(1:n+1),  intent(in) :: xi
  real(lp), dimension(1:n+1), intent(out) :: weights
  !
  !.. Local Scalars ..
  integer  :: i
  real(lp) :: a1,a2,b1,b2,ab,ab1,ab2
  real(lp) :: c,c1,c2,c3,ri,rn,su
  real(lp) :: y,z,dy,d2y,scl
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: half = 0.5_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: two  = 2.0_lp
  !
continue
  !
  weights(:) = zero
  !
  if (n > 1) then
    !
    a1  = alpha + one
    a2  = alpha + two
    b1  = beta  + one
    b2  = beta  + two
    ab  = alpha + beta
    ab1 = alpha + beta + one
    ab2 = alpha + beta + two
    !
    c = (two**ab) * gammafun(a1) * gammafun(b1) / gammafun(ab2)
    !
    rn = real(n,kind=lp)
    c  = c * (two*rn + ab) / (rn + ab1)
    c1 = c * a1 / (b2 * ab2)
    c2 = c * b1 / (a2 * ab2)
    c3 = half * c * a1 * b1
    !
    do i = 2,n-1
      ri = real(i-1,kind=lp)
      c1 = c1*(ri+a1)*ri / ((ri+ab2)*(ri+b2))
      c2 = c2*(ri+b1)*ri / ((ri+ab2)*(ri+a2))
      c3 = c3*(ri+a1)*(ri+b1) / ((ri+two)*(ri+ab1))
    end do
    !
    su = zero
    do i = 2,n
      su = su + xi(i)
    end do
    !
    scl = ( (one-xi(1))**a1 ) * ( (one+xi(1))**b1 )
    weights(1)   = scl*c1*(rn-one-su)
    !
    scl = ( (one-xi(n+1))**a1 ) * ( (one+xi(n+1))**b1 )
    weights(n+1) = scl*c2*(rn-one+su)
    !
    do i = 2,n
      call eval_jacobi_poly(n  ,alpha,beta,xi(i),y,dy,d2y)
      call eval_jacobi_poly(n-1,alpha,beta,xi(i),z,dy,d2y)
      scl = ( (one-xi(i))**a1 ) * ( (one+xi(i))**b1 )
      weights(i) = -scl*c3/y/dy
    end do
    !
  end if
  !
end subroutine weights_dgbook_jacobi_gauss_lobatto_DP
!
!###############################################################################
!
pure subroutine weights_dgbook_jacobi_gauss_lobatto_QP(n,alpha,beta,xi,weights)
  !
  !********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE JACOBI GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   alpha  = PARAMETER > -1
  !   beta  = PARAMETER > -1
  !   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !********************************************************************
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: n
  real(lp),                    intent(in) :: alpha
  real(lp),                    intent(in) :: beta
  real(lp), dimension(1:n+1),  intent(in) :: xi
  real(lp), dimension(1:n+1), intent(out) :: weights
  !
  !.. Local Scalars ..
  integer  :: i
  real(lp) :: a1,a2,b1,b2,ab,ab1,ab2
  real(lp) :: c,c1,c2,c3,ri,rn,su
  real(lp) :: y,z,dy,d2y,scl
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: half = 0.5_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: two  = 2.0_lp
  !
continue
  !
  weights(:) = zero
  !
  if (n > 1) then
    !
    a1  = alpha + one
    a2  = alpha + two
    b1  = beta  + one
    b2  = beta  + two
    ab  = alpha + beta
    ab1 = alpha + beta + one
    ab2 = alpha + beta + two
    !
    c = (two**ab) * gammafun(a1) * gammafun(b1) / gammafun(ab2)
    !
    rn = real(n,kind=lp)
    c  = c * (two*rn + ab) / (rn + ab1)
    c1 = c * a1 / (b2 * ab2)
    c2 = c * b1 / (a2 * ab2)
    c3 = half * c * a1 * b1
    !
    do i = 2,n-1
      ri = real(i-1,kind=lp)
      c1 = c1*(ri+a1)*ri / ((ri+ab2)*(ri+b2))
      c2 = c2*(ri+b1)*ri / ((ri+ab2)*(ri+a2))
      c3 = c3*(ri+a1)*(ri+b1) / ((ri+two)*(ri+ab1))
    end do
    !
    su = zero
    do i = 2,n
      su = su + xi(i)
    end do
    !
    scl = ( (one-xi(1))**a1 ) * ( (one+xi(1))**b1 )
    weights(1)   = scl*c1*(rn-one-su)
    !
    scl = ( (one-xi(n+1))**a1 ) * ( (one+xi(n+1))**b1 )
    weights(n+1) = scl*c2*(rn-one+su)
    !
    do i = 2,n
      call eval_jacobi_poly(n  ,alpha,beta,xi(i),y,dy,d2y)
      call eval_jacobi_poly(n-1,alpha,beta,xi(i),z,dy,d2y)
      scl = ( (one-xi(i))**a1 ) * ( (one+xi(i))**b1 )
      weights(i) = -scl*c3/y/dy
    end do
    !
  end if
  !
end subroutine weights_dgbook_jacobi_gauss_lobatto_QP
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
  real(lp), intent(in) :: a ! REAL VALUE TO DISTINGUISH
  !                           BETWEEN DP AND QP VERSIONS
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
  if (np <= 1) return
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
!###############################################################################
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
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: n ! ORDER OF LEGENDRE GAUSS POLYNOMIAL
  real(lp), intent(in) :: a ! REAL VALUE TO DISTINGUISH
  !                           BETWEEN DP AND QP VERSIONS
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
  if (np <= 1) return
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
!###############################################################################
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
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                       intent(in) :: n
  real(lp), dimension(1:n),      intent(in) :: cs
  real(lp), dimension(1:n),      intent(in) :: dz
  real(lp), dimension(1:n,1:n), intent(out) :: dma
  !
  !.. Local Scalars ..
  integer  :: i,j
  real(lp) :: vi,zi,vj,zj
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  dma(1,1) = zero
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
        dma(i,i) = zi/(one-zi*zi)
      end if
    end do
  end do
  !
end subroutine old_derivative_matrix_legendre_gauss_DP
!
!###############################################################################
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
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                       intent(in) :: n
  real(lp), dimension(1:n),      intent(in) :: cs
  real(lp), dimension(1:n),      intent(in) :: dz
  real(lp), dimension(1:n,1:n), intent(out) :: dma
  !
  !.. Local Scalars ..
  integer  :: i,j
  real(lp) :: vi,zi,vj,zj
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  dma(1,1) = zero
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
        dma(i,i) = zi/(one-zi*zi)
      end if
    end do
  end do
  !
end subroutine old_derivative_matrix_legendre_gauss_QP
!
!###############################################################################
!
pure function derivative_matrix_legendre_gauss_lobatto_DP(n,a) &
                                                   result(return_value)
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
  real(lp), intent(in) :: a ! REAL VALUE TO DISTINGUISH
  !                           BETWEEN DP AND QP VERSIONS
  !.. Function Result ..
  real(lp), dimension(1:n+1,1:n+1) :: return_value ! DERIVATIVE MATRIX
  !
  !.. Local Scalars ..
  integer  :: i,j,np
  real(lp) :: ei,ej,vi,vj,dn,c
  !
  !.. Local Arrays ..
  real(lp), dimension(1:n+1) :: et ! VECTOR OF THE NODES
  real(lp), dimension(1:n+1) :: vn ! VALUES OF THE LEGENDRE POLYNOMIAL
  !                                  AT THE NODES
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
!###############################################################################
!
pure function derivative_matrix_legendre_gauss_lobatto_QP(n,a) &
                                                   result(return_value)
  !
  !***********************************************************************
  !  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE
  !  LEGENDRE GAUSS-LOBATTO NODES
  !***********************************************************************
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: n ! ORDER OF LEGENDRE GAUSS-LOBATTO POLYNOMIAL
  real(lp), intent(in) :: a ! REAL VALUE TO DISTINGUISH
  !                           BETWEEN DP AND QP VERSIONS
  !.. Function Result ..
  real(lp), dimension(1:n+1,1:n+1) :: return_value ! DERIVATIVE MATRIX
  !
  !.. Local Scalars ..
  integer  :: i,j,np
  real(lp) :: ei,ej,vi,vj,dn,c
  !
  !.. Local Arrays ..
  real(lp), dimension(1:n+1) :: et ! VECTOR OF THE NODES
  real(lp), dimension(1:n+1) :: vn ! VALUES OF THE LEGENDRE POLYNOMIAL
  !                                  AT THE NODES
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
!###############################################################################
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
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                       intent(in) :: n
  real(lp), dimension(1:n),      intent(in) :: et
  real(lp), dimension(1:n),      intent(in) :: vn
  real(lp), dimension(1:n,1:n), intent(out) :: dma
  !
  !.. Local Scalars ..
  integer  :: i,j
  real(lp) :: ei,ej,vi,vj,dn,c
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one4 = 0.25_lp
  real(lp), parameter :: one  = 1.0_lp
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
        dma(i,i) = zero
      end if
    end do
  end do
  !
  if (n > 1) then
    dn = real(n-1,kind=lp)
    c  = one4*dn*(dn+one)
    dma(1,1) = -c
    dma(n,n) = c
  end if
  !
end subroutine old_derivative_matrix_legendre_gauss_lobatto_DP
!
!###############################################################################
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
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                       intent(in) :: n
  real(lp), dimension(1:n),      intent(in) :: et
  real(lp), dimension(1:n),      intent(in) :: vn
  real(lp), dimension(1:n,1:n), intent(out) :: dma
  !
  !.. Local Scalars ..
  integer  :: i,j
  real(lp) :: ei,ej,vi,vj,dn,c
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one4 = 0.25_lp
  real(lp), parameter :: one  = 1.0_lp
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
        dma(i,i) = zero
      end if
    end do
  end do
  !
  if (n > 1) then
    dn = real(n-1,kind=lp)
    c  = one4*dn*(dn+one)
    dma(1,1) = -c
    dma(n,n) = c
  end if
  !
end subroutine old_derivative_matrix_legendre_gauss_lobatto_QP
!
!###############################################################################
!
pure subroutine derivative_matrix_jacobi_gauss_lobatto_DP(n,alpha,beta, &
                                                          et,vn,dma)
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
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                       intent(in) :: n
  real(lp),                      intent(in) :: alpha
  real(lp),                      intent(in) :: beta
  real(lp), dimension(1:n),      intent(in) :: et
  real(lp), dimension(1:n),      intent(in) :: vn
  real(lp), dimension(1:n,1:n), intent(out) :: dma
  !
  !.. Local Scalars ..
  integer  :: i,j
  real(lp) :: a1,b1,ab,dn,ei,ej
  real(lp) :: c1,c2,c3,c4,vi,vj
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: half = 0.5_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: two  = 2.0_lp
  !
continue
  !
  dma(1,1) = zero
  if (n == 1) return
  !
  a1 = alpha+one
  b1 = beta+one
  ab = alpha+beta
  dn = real(n-1,kind=lp)
  c1 = dn*(dn+ab+one)
  c2 = a1*vn(1)/(b1*vn(n))
  dma(1,1) = half*(alpha-c1)/(beta+two)
  dma(n,n) = half*(c1-beta)/(alpha+two)
  dma(1,n) = -half*c2
  dma(n,1) = half/c2
  if (n == 2) return
  !
  c3 = vn(1)/b1
  c4 = vn(n)/a1
  do j = 2,n-1
    vj = vn(j)
    ej = et(j)
    dma(1,j) = -c3/(vj*(one+ej))
    dma(n,j) = c4/(vj*(one-ej))
    dma(j,1) = vj/(c3*(one+ej))
    dma(j,n) = -vj/(c4*(one-ej))
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
        dma(i,i) = half*(ab*ei+alpha-beta)/(one-ei*ei)
      end if
    end do
  end do
  !
end subroutine derivative_matrix_jacobi_gauss_lobatto_DP
!
!###############################################################################
!
pure subroutine derivative_matrix_jacobi_gauss_lobatto_QP(n,alpha,beta, &
                                                          et,vn,dma)
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
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                       intent(in) :: n
  real(lp),                      intent(in) :: alpha
  real(lp),                      intent(in) :: beta
  real(lp), dimension(1:n),      intent(in) :: et
  real(lp), dimension(1:n),      intent(in) :: vn
  real(lp), dimension(1:n,1:n), intent(out) :: dma
  !
  !.. Local Scalars ..
  integer  :: i,j
  real(lp) :: a1,b1,ab,dn,ei,ej
  real(lp) :: c1,c2,c3,c4,vi,vj
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: half = 0.5_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: two  = 2.0_lp
  !
continue
  !
  dma(1,1) = zero
  if (n == 1) return
  !
  a1 = alpha+one
  b1 = beta+one
  ab = alpha+beta
  dn = real(n-1,kind=lp)
  c1 = dn*(dn+ab+one)
  c2 = a1*vn(1)/(b1*vn(n))
  dma(1,1) = half*(alpha-c1)/(beta+two)
  dma(n,n) = half*(c1-beta)/(alpha+two)
  dma(1,n) = -half*c2
  dma(n,1) = half/c2
  if (n == 2) return
  !
  c3 = vn(1)/b1
  c4 = vn(n)/a1
  do j = 2,n-1
    vj = vn(j)
    ej = et(j)
    dma(1,j) = -c3/(vj*(one+ej))
    dma(n,j) = c4/(vj*(one-ej))
    dma(j,1) = vj/(c3*(one+ej))
    dma(j,n) = -vj/(c4*(one-ej))
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
        dma(i,i) = half*(ab*ei+alpha-beta)/(one-ei*ei)
      end if
    end do
  end do
  !
end subroutine derivative_matrix_jacobi_gauss_lobatto_QP
!
!###############################################################################
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
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                       intent(in) :: n
  real(lp), dimension(1:n),      intent(in) :: et
  real(lp), dimension(1:n,1:n), intent(out) :: dma
  !
  !.. Local Scalars ..
  integer  :: i,j
  real(lp) :: dn,cn,ej,ei
  real(lp) :: sn,sgni,sgnj
  !
  !.. Local Constants ..
  real(lp), parameter :: zero  = 0.0_lp
  real(lp), parameter :: half  = 0.5_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  real(lp), parameter :: six   = 6.0_lp
  !
continue
  !
  dma(1,1) = zero
  if (n == 1) return
  !
  dn = real((n-1),kind=lp)
  cn = (two*dn*dn+one)/six
  sn = real(1+4*((n-1)/2)-2*(n-1),kind=lp)
  dma(1,1) = -cn
  dma(n,n) = cn
  dma(1,n) = -half*sn
  dma(n,1) = half*sn
  if (n == 2) return
  !
  sgnj = -one
  do j = 2,n-1
    ej = et(j)
    dma(1,j) = -two*sgnj/(one+ej)
    dma(n,j) = two*sgnj*sn/(one-ej)
    dma(j,1) = half*sgnj/(one+ej)
    dma(j,n) = -half*sgnj*sn/(one-ej)
    sgnj = -sgnj
  end do
  !
  sgni = -one
  do i = 2,n-1
    ei = et(i)
    sgnj = -one
    do j = 2,n-1
      if (i /= j) then
        ej = et(j)
        dma(i,j) = sgni*sgnj/(ei-ej)
      else
        dma(i,i) = -half*ei/(one-ei*ei)
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
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                       intent(in) :: n
  real(lp), dimension(1:n),      intent(in) :: et
  real(lp), dimension(1:n,1:n), intent(out) :: dma
  !
  !.. Local Scalars ..
  integer  :: i,j
  real(lp) :: dn,cn,ej,ei
  real(lp) :: sn,sgni,sgnj
  !
  !.. Local Constants ..
  real(lp), parameter :: zero  = 0.0_lp
  real(lp), parameter :: half  = 0.5_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  real(lp), parameter :: six   = 6.0_lp
  !
continue
  !
  dma(1,1) = zero
  if (n == 1) return
  !
  dn = real((n-1),kind=lp)
  cn = (two*dn*dn+one)/six
  sn = real(1+4*((n-1)/2)-2*(n-1),kind=lp)
  dma(1,1) = -cn
  dma(n,n) = cn
  dma(1,n) = -half*sn
  dma(n,1) = half*sn
  if (n == 2) return
  !
  sgnj = -one
  do j = 2,n-1
    ej = et(j)
    dma(1,j) = -two*sgnj/(one+ej)
    dma(n,j) = two*sgnj*sn/(one-ej)
    dma(j,1) = half*sgnj/(one+ej)
    dma(j,n) = -half*sgnj*sn/(one-ej)
    sgnj = -sgnj
  end do
  !
  sgni = -one
  do i = 2,n-1
    ei = et(i)
    sgnj = -one
    do j = 2,n-1
      if (i /= j) then
        ej = et(j)
        dma(i,j) = sgni*sgnj/(ei-ej)
      else
        dma(i,i) = -half*ei/(one-ei*ei)
      end if
      sgnj = -sgnj
    end do
    sgni = -sgni
  end do
  !
end subroutine derivative_matrix_chebyshev_gauss_lobatto_QP
!
!###############################################################################
!
elemental function gammafun_DP(x) result(return_value)
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  real(lp), intent(in) :: x
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
  !.. Local Scalars ..
  integer  :: ind,k
  real(lp) :: xx,pr,s,g
  !
  !.. Local Constants ..
  real(lp), parameter :: half  = 0.5_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: nine2 = 4.5_lp
  real(lp), parameter :: eps14 = 1.0e-14_lp
  real(lp), parameter :: pi    = 3.14159265358979323846_lp
  !
  real(lp), parameter :: s_start = 0.426401432711220868_lp
  real(lp), parameter :: gamma_constant(1:11) = &
             [ -0.524741987629368444e+00_lp,  0.116154405493589130e+00_lp, &
               -0.765978624506602380e-02_lp,  0.899719449391378898e-04_lp, &
               -0.194536980009534621e-07_lp,  0.199382839513630987e-10_lp, &
               -0.204209590209541319e-11_lp,  0.863896817907000175e-13_lp, &
                0.152237501608472336e-13_lp, -0.825725175277719950e-14_lp, &
                0.299734782205224610e-14_lp ]
  !
continue
  !
  xx = x
  return_value = one
  !
  do
    if (abs(xx-one) < eps14) return
    if (xx < one) exit
    xx = xx - one
    return_value = return_value*xx
  end do
  !
  ind = 0
  !
  if (xx < half) then
    ind = 1
    return_value = return_value*pi/sin(pi*xx)
    xx = one - xx
  end if
  !
  pr = one
  s = s_start
  !
  do k = 1,11
    pr = pr * (xx-real(k,kind=lp)) / (xx+real(k-1,kind=lp))
    s  = s + gamma_constant(k)*pr
  end do
  !
  g = s * exp(one-xx) * (xx+nine2)**(xx-half)
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
elemental function gammafun_QP(x) result(return_value)
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  real(lp), intent(in) :: x
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
  !.. Local Scalars ..
  integer  :: ind,k
  real(lp) :: xx,pr,s,g
  !
  !.. Local Constants ..
  real(lp), parameter :: half  = 0.5_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: nine2 = 4.5_lp
  real(lp), parameter :: eps14 = 1.0e-14_lp
  real(lp), parameter :: pi    = 3.14159265358979323846_lp
  !
  real(lp), parameter :: s_start = 0.426401432711220868_lp
  real(lp), parameter :: gamma_constant(1:11) = &
             [ -0.524741987629368444e+00_lp,  0.116154405493589130e+00_lp, &
               -0.765978624506602380e-02_lp,  0.899719449391378898e-04_lp, &
               -0.194536980009534621e-07_lp,  0.199382839513630987e-10_lp, &
               -0.204209590209541319e-11_lp,  0.863896817907000175e-13_lp, &
                0.152237501608472336e-13_lp, -0.825725175277719950e-14_lp, &
                0.299734782205224610e-14_lp ]
  !
continue
  !
  xx = x
  return_value = one
  !
  do
    if (abs(xx-one) < eps14) return
    if (xx < one) exit
    xx = xx - one
    return_value = return_value*xx
  end do
  !
  ind = 0
  !
  if (xx < half) then
    ind = 1
    return_value = return_value*pi/sin(pi*xx)
    xx = one - xx
  end if
  !
  pr = one
  s = s_start
  !
  do k = 1,11
    pr = pr * (xx-real(k,kind=lp)) / (xx+real(k-1,kind=lp))
    s  = s + gamma_constant(k)*pr
  end do
  !
  g = s * exp(one-xx) * (xx+nine2)**(xx-half)
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
elemental function log_gammafun_DP(x) result(return_value)
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  real(lp), intent(in) :: x
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
continue
  !
  return_value = log( gammafun(x) )
  !
end function log_gammafun_DP
!
!###############################################################################
!
elemental function log_gammafun_QP(x) result(return_value)
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  real(lp), intent(in) :: x
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
continue
  !
  return_value = log( gammafun(x) )
  !
end function log_gammafun_QP
!
!###############################################################################
!
pure function coefrr_DP(n,a) result(return_value)
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: n
  real(lp), intent(in) :: a ! REAL VALUE TO DISTINGUISH
  !                           BETWEEN DP AND QP VERSIONS
  !.. Function Result ..
  real(lp), dimension(0:n) :: return_value
  !
  !.. Local Scalars ..
  integer :: k,ip
  !
  !.. Local Arrays ..
  real(lp), dimension(0:n)   :: cfk
  real(lp), dimension(0:n)   :: cfkmt
  real(lp), dimension(0:n-1) :: cfkm1
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: half = 0.5_lp
  !
continue
  !
  cfk = coefle(n,zero)
  cfkm1 = coefle(n-1,zero)
  !
  do k = 0,n-1
    cfkmt(k) = cfkm1(k)
  end do
  !
  cfkmt(n) = zero
  ip = (-1)**n
  !
  do k = 0,n
    return_value(k) = half * real(ip,kind=lp) * (cfk(k) - cfkmt(k))
  end do
  !
end function coefrr_DP
!
!###############################################################################
!
pure function coefrr_QP(n,a) result(return_value)
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: n
  real(lp), intent(in) :: a ! REAL VALUE TO DISTINGUISH
  !                           BETWEEN DP AND QP VERSIONS
  !.. Function Result ..
  real(lp), dimension(0:n) :: return_value
  !
  !.. Local Scalars ..
  integer :: k,ip
  !
  !.. Local Arrays ..
  real(lp), dimension(0:n)   :: cfk
  real(lp), dimension(0:n)   :: cfkmt
  real(lp), dimension(0:n-1) :: cfkm1
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: half = 0.5_lp
  !
continue
  !
  cfk = coefle(n,zero)
  cfkm1 = coefle(n-1,zero)
  !
  do k = 0,n-1
    cfkmt(k) = cfkm1(k)
  end do
  !
  cfkmt(n) = zero
  ip = (-1)**n
  !
  do k = 0,n
    return_value(k) = half * real(ip,kind=lp) * (cfk(k) - cfkmt(k))
  end do
  !
end function coefrr_QP
!
!###############################################################################
!
pure function coefle_DP(n,a) result(return_value)
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: n
  real(lp), intent(in) :: a ! REAL VALUE TO DISTINGUISH
  !                           BETWEEN DP AND QP VERSIONS
  !.. Function Result ..
  real(lp), dimension(0:n) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,k,m
  !
  !.. Local Arrays ..
  real(lp), dimension(-1:n,-1:n) :: cle
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  cle(:,:) = zero
  !
  cle(0,0) = one
  cle(1,1) = one
  !
  do k = 2,n
    do m = 0,k
      cle(k,m) = (real(2*k-1,kind=lp)/real(k,kind=lp))*cle(k-1,m-1) - &
                 (real(  k-1,kind=lp)/real(k,kind=lp))*cle(k-2,m  )
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
pure function coefle_QP(n,a) result(return_value)
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: n
  real(lp), intent(in) :: a ! REAL VALUE TO DISTINGUISH
  !                           BETWEEN DP AND QP VERSIONS
  !.. Function Result ..
  real(lp), dimension(0:n) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,k,m
  !
  !.. Local Arrays ..
  real(lp), dimension(-1:n,-1:n) :: cle
  !
  !.. Local Constants ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  cle(:,:) = zero
  !
  cle(0,0) = one
  cle(1,1) = one
  !
  do k = 2,n
    do m = 0,k
      cle(k,m) = (real(2*k-1,kind=lp)/real(k,kind=lp))*cle(k-1,m-1) - &
                 (real(  k-1,kind=lp)/real(k,kind=lp))*cle(k-2,m  )
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
pure function eval_LagrangePoly_DP(ijk,xyz,xi) result(return_value)
  !
  ! Get the value of a tensor product of Lagrange polynomials
  ! of arbitrary dimension
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
  !.. Function Result ..
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
pure function eval_LagrangePoly_QP(ijk,xyz,xi) result(return_value)
  !
  ! Get the value of a tensor product of Lagrange polynomials
  ! of arbitrary dimension
  !
  ! IJK : The value of the Lagrange polynomial tensor product is 1 at ijk
  ! XYZ : Location at which to evaluate the Lagrange polynomial
  ! XI  : Array of the 1D solution points
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,  dimension(:), intent(in) :: ijk
  real(lp), dimension(:), intent(in) :: xyz
  real(lp), dimension(:), intent(in) :: xi
  !
  !.. Function Result ..
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
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                intent(in) :: k
  real(lp),               intent(in) :: x
  real(lp), dimension(:), intent(in) :: xi
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
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
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                intent(in) :: k
  real(lp),               intent(in) :: x
  real(lp), dimension(:), intent(in) :: xi
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
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
pure function eval_LagrangePoly1D_arrays_DP(x,xi) result(return_value)
  !
  ! Get the values of the 1D Lagrange polynomials at point x
  !
  ! X  : Location at which to evaluate the k-th Lagrange polynomial
  ! XI : Array of all the interpolation points
  !
  ! This evaluates the Lagrange polynomial at the location x.
  ! The Lagrange polynomial corresponds to a function that
  ! is zero at all iterpolation points xi(i) for i=1,nfp
  ! except at xi(k) where the function is equal to 1.
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  real(lp),               intent(in) :: x
  real(lp), dimension(:), intent(in) :: xi
  !
  !.. Function Result ..
  real(lp) :: return_value(1:size(xi))
  !
  !.. Local Scalars ..
  integer :: i,k,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
  !
continue
  !
  nfp = size(xi)
  imask(:) = intseq(1,nfp)
  !
  do k = 1,nfp
    return_value(k) = product( x    -xi(1:nfp) , imask /= k ) / &
                      product( xi(k)-xi(1:nfp) , imask /= k )
  end do
  !
end function eval_LagrangePoly1D_arrays_DP
!
!###############################################################################
!
pure function eval_LagrangePoly1D_arrays_QP(x,xi) result(return_value)
  !
  ! Get the values of the 1D Lagrange polynomials at point x
  !
  ! X  : Location at which to evaluate the k-th Lagrange polynomial
  ! XI : Array of all the interpolation points
  !
  ! This evaluates the Lagrange polynomial at the location x.
  ! The Lagrange polynomial corresponds to a function that
  ! is zero at all iterpolation points xi(i) for i=1,nfp
  ! except at xi(k) where the function is equal to 1.
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  real(lp),               intent(in) :: x
  real(lp), dimension(:), intent(in) :: xi
  !
  !.. Function Result ..
  real(lp) :: return_value(1:size(xi))
  !
  !.. Local Scalars ..
  integer :: i,k,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
  !
continue
  !
  nfp = size(xi)
  imask(:) = intseq(1,nfp)
  !
  do k = 1,nfp
    return_value(k) = product( x    -xi(1:nfp) , imask /= k ) / &
                      product( xi(k)-xi(1:nfp) , imask /= k )
  end do
  !
end function eval_LagrangePoly1D_arrays_QP
!
!###############################################################################
!
pure function eval_LagrangePoly2D_arrays_DP(ij,xy,xi,eta) result(return_value)
  !
  ! Get value of the 2D Lagrange polynomial at point x,y
  !
  ! K   : The value of the Lagrange polynomial in coordinate XI  is 1 at xi(k)
  ! L   : The value of the Lagrange polynomial in coordinate ETA is 1 at eta(l)
  ! X   : Location at which to evaluate the k-th Lagrange polynomial
  !       in coordinate XI
  ! Y   : Location at which to evaluate the l-th Lagrange polynomial in
  !       coordinate ETA
  ! XI  : Array of all the interpolation points in XI  coordinate
  ! ETA : Array of all the interpolation points in ETA coordinate
  !
  ! This evaluates the 2D Lagrange polynomial at the location x,y.
  ! The 2D Lagrange polynomial corresponds to a function that
  ! is zero at all iterpolation points [xi(i),eta(j)] for [i=1,ifp;j=1,jfp]
  ! except at [xi(k),eta(l)] where the function is equal to 1.
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,  dimension(:), intent(in) :: ij
  real(lp), dimension(:), intent(in) :: xy
  real(lp), dimension(:), intent(in) :: xi
  real(lp), dimension(:), intent(in) :: eta
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
continue
  !
  return_value = eval_LagrangePoly1D(ij(1),xy(1),xi) * &
                 eval_LagrangePoly1D(ij(2),xy(2),eta)
  !
end function eval_LagrangePoly2D_arrays_DP
!
!###############################################################################
!
pure function eval_LagrangePoly2D_arrays_QP(ij,xy,xi,eta) result(return_value)
  !
  ! Get value of the 2D Lagrange polynomial at point x,y
  !
  ! K   : The value of the Lagrange polynomial in coordinate XI  is 1 at xi(k)
  ! L   : The value of the Lagrange polynomial in coordinate ETA is 1 at eta(l)
  ! X   : Location at which to evaluate the k-th Lagrange polynomial
  !       in coordinate XI
  ! Y   : Location at which to evaluate the l-th Lagrange polynomial in
  !       coordinate ETA
  ! XI  : Array of all the interpolation points in XI  coordinate
  ! ETA : Array of all the interpolation points in ETA coordinate
  !
  ! This evaluates the 2D Lagrange polynomial at the location x,y.
  ! The 2D Lagrange polynomial corresponds to a function that
  ! is zero at all iterpolation points [xi(i),eta(j)] for [i=1,ifp;j=1,jfp]
  ! except at [xi(k),eta(l)] where the function is equal to 1.
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,  dimension(:), intent(in) :: ij
  real(lp), dimension(:), intent(in) :: xy
  real(lp), dimension(:), intent(in) :: xi
  real(lp), dimension(:), intent(in) :: eta
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
continue
  !
  return_value = eval_LagrangePoly1D(ij(1),xy(1),xi) * &
                 eval_LagrangePoly1D(ij(2),xy(2),eta)
  !
end function eval_LagrangePoly2D_arrays_QP
!
!###############################################################################
!
pure function eval_LagrangePoly2D_DP(k,l,x,y,xi,eta) result(return_value)
  !
  ! Get value of the 2D Lagrange polynomial at point x,y
  !
  ! K   : The value of the Lagrange polynomial in coordinate XI  is 1 at xi(k)
  ! L   : The value of the Lagrange polynomial in coordinate ETA is 1 at eta(l)
  ! X   : Location at which to evaluate the k-th Lagrange polynomial
  !       in coordinate XI
  ! Y   : Location at which to evaluate the l-th Lagrange polynomial in
  !       coordinate ETA
  ! XI  : Array of all the interpolation points in XI  coordinate
  ! ETA : Array of all the interpolation points in ETA coordinate
  !
  ! This evaluates the 2D Lagrange polynomial at the location x,y.
  ! The 2D Lagrange polynomial corresponds to a function that
  ! is zero at all iterpolation points [xi(i),eta(j)] for [i=1,ifp;j=1,jfp]
  ! except at [xi(k),eta(l)] where the function is equal to 1.
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                intent(in) :: k
  integer,                intent(in) :: l
  real(lp),               intent(in) :: x
  real(lp),               intent(in) :: y
  real(lp), dimension(:), intent(in) :: xi
  real(lp), dimension(:), intent(in) :: eta
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
continue
  !
  return_value = eval_LagrangePoly1D(k,x,xi) * eval_LagrangePoly1D(l,y,eta)
  !
end function eval_LagrangePoly2D_DP
!
!###############################################################################
!
pure function eval_LagrangePoly2D_QP(k,l,x,y,xi,eta) result(return_value)
  !
  ! Get value of the 2D Lagrange polynomial at point x,y
  !
  ! K   : The value of the Lagrange polynomial in coordinate XI  is 1 at xi(k)
  ! L   : The value of the Lagrange polynomial in coordinate ETA is 1 at eta(l)
  ! X   : Location at which to evaluate the k-th Lagrange polynomial
  !       in coordinate XI
  ! Y   : Location at which to evaluate the l-th Lagrange polynomial in
  !       coordinate ETA
  ! XI  : Array of all the interpolation points in XI  coordinate
  ! ETA : Array of all the interpolation points in ETA coordinate
  !
  ! This evaluates the 2D Lagrange polynomial at the location x,y.
  ! The 2D Lagrange polynomial corresponds to a function that
  ! is zero at all iterpolation points [xi(i),eta(j)] for [i=1,ifp;j=1,jfp]
  ! except at [xi(k),eta(l)] where the function is equal to 1.
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                intent(in) :: k
  integer,                intent(in) :: l
  real(lp),               intent(in) :: x
  real(lp),               intent(in) :: y
  real(lp), dimension(:), intent(in) :: xi
  real(lp), dimension(:), intent(in) :: eta
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
continue
  !
  return_value = eval_LagrangePoly1D(k,x,xi) * eval_LagrangePoly1D(l,y,eta)
  !
end function eval_LagrangePoly2D_QP
!
!###############################################################################
!
pure function eval_LagrangePoly3D_arrays_DP(ijk,xyz,xi,eta,zeta) &
                                     result(return_value)
  !
  ! Get value of the 3D Lagrange polynomial at point x,y,z
  !
  ! I    : The value of the Lagrange polynomial in coord XI   is 1 at xi(i)
  ! J    : The value of the Lagrange polynomial in coord ETA  is 1 at eta(j)
  ! K    : The value of the Lagrange polynomial in coord ZETA is 1 at zeta(k)
  ! X    : Location at which to evaluate the i-th Lagrange polynomial
  !        in coordinate XI
  ! Y    : Location at which to evaluate the j-th Lagrange polynomial in
  !        coordinate ETA
  ! Z    : Location at which to evaluate the k-th Lagrange polynomial in
  !        coordinate ZETA
  ! XI   : Array of all the interpolation points in XI   coordinate
  ! ETA  : Array of all the interpolation points in ETA  coordinate
  ! ZETA : Array of all the interpolation points in ZETA coordinate
  !
  ! This evaluates the 3D Lagrange polynomial at the location x,y,z.
  ! The 3D Lagrange polynomial corresponds to a function that
  ! is zero at all iterpolation points [xi(a),eta(b),zeta(c)] for
  ! [a=1,size(xi);b=1,size(eta),c=1,size(zeta)] except at
  ! [xi(i),eta(j),zeta(k)] where the function is equal to 1.
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,  dimension(:), intent(in) :: ijk
  real(lp), dimension(:), intent(in) :: xyz
  real(lp), dimension(:), intent(in) :: xi
  real(lp), dimension(:), intent(in) :: eta
  real(lp), dimension(:), intent(in) :: zeta
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
continue
  !
  return_value = eval_LagrangePoly1D(ijk(1),xyz(1),xi) * &
                 eval_LagrangePoly1D(ijk(2),xyz(2),eta) * &
                 eval_LagrangePoly1D(ijk(3),xyz(3),zeta)
  !
end function eval_LagrangePoly3D_arrays_DP
!
!###############################################################################
!
pure function eval_LagrangePoly3D_arrays_QP(ijk,xyz,xi,eta,zeta) &
                                     result(return_value)
  !
  ! Get value of the 3D Lagrange polynomial at point x,y,z
  !
  ! I    : The value of the Lagrange polynomial in coord XI   is 1 at xi(i)
  ! J    : The value of the Lagrange polynomial in coord ETA  is 1 at eta(j)
  ! K    : The value of the Lagrange polynomial in coord ZETA is 1 at zeta(k)
  ! X    : Location at which to evaluate the i-th Lagrange polynomial
  !        in coordinate XI
  ! Y    : Location at which to evaluate the j-th Lagrange polynomial in
  !        coordinate ETA
  ! Z    : Location at which to evaluate the k-th Lagrange polynomial in
  !        coordinate ZETA
  ! XI   : Array of all the interpolation points in XI   coordinate
  ! ETA  : Array of all the interpolation points in ETA  coordinate
  ! ZETA : Array of all the interpolation points in ZETA coordinate
  !
  ! This evaluates the 3D Lagrange polynomial at the location x,y,z.
  ! The 3D Lagrange polynomial corresponds to a function that
  ! is zero at all iterpolation points [xi(a),eta(b),zeta(c)] for
  ! [a=1,size(xi);b=1,size(eta),c=1,size(zeta)] except at
  ! [xi(i),eta(j),zeta(k)] where the function is equal to 1.
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,  dimension(:), intent(in) :: ijk
  real(lp), dimension(:), intent(in) :: xyz
  real(lp), dimension(:), intent(in) :: xi
  real(lp), dimension(:), intent(in) :: eta
  real(lp), dimension(:), intent(in) :: zeta
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
continue
  !
  return_value = eval_LagrangePoly1D(ijk(1),xyz(1),xi) * &
                 eval_LagrangePoly1D(ijk(2),xyz(2),eta) * &
                 eval_LagrangePoly1D(ijk(3),xyz(3),zeta)
  !
end function eval_LagrangePoly3D_arrays_QP
!
!###############################################################################
!
pure function eval_LagrangePoly3D_DP(i,j,k,x,y,z,xi,eta,zeta) &
                              result(return_value)
  !
  ! Get value of the 3D Lagrange polynomial at point x,y,z
  !
  ! I    : The value of the Lagrange polynomial in coord XI   is 1 at xi(i)
  ! J    : The value of the Lagrange polynomial in coord ETA  is 1 at eta(j)
  ! K    : The value of the Lagrange polynomial in coord ZETA is 1 at zeta(k)
  ! X    : Location at which to evaluate the i-th Lagrange polynomial
  !        in coordinate XI
  ! Y    : Location at which to evaluate the j-th Lagrange polynomial in
  !        coordinate ETA
  ! Z    : Location at which to evaluate the k-th Lagrange polynomial in
  !        coordinate ZETA
  ! XI   : Array of all the interpolation points in XI   coordinate
  ! ETA  : Array of all the interpolation points in ETA  coordinate
  ! ZETA : Array of all the interpolation points in ZETA coordinate
  !
  ! This evaluates the 3D Lagrange polynomial at the location x,y,z.
  ! The 3D Lagrange polynomial corresponds to a function that
  ! is zero at all iterpolation points [xi(a),eta(b),zeta(c)] for
  ! [a=1,size(xi);b=1,size(eta),c=1,size(zeta)] except at
  ! [xi(i),eta(j),zeta(k)] where the function is equal to 1.
  !
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                intent(in) :: i,j,k
  real(lp),               intent(in) :: x,y,z
  real(lp), dimension(:), intent(in) :: xi
  real(lp), dimension(:), intent(in) :: eta
  real(lp), dimension(:), intent(in) :: zeta
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
continue
  !
  return_value = eval_LagrangePoly1D(i,x,xi) * &
                 eval_LagrangePoly1D(j,y,eta) * &
                 eval_LagrangePoly1D(k,z,zeta)
  !
end function eval_LagrangePoly3D_DP
!
!###############################################################################
!
pure function eval_LagrangePoly3D_QP(i,j,k,x,y,z,xi,eta,zeta) &
                              result(return_value)
  !
  ! Get value of the 3D Lagrange polynomial at point x,y,z
  !
  ! I    : The value of the Lagrange polynomial in coord XI   is 1 at xi(i)
  ! J    : The value of the Lagrange polynomial in coord ETA  is 1 at eta(j)
  ! K    : The value of the Lagrange polynomial in coord ZETA is 1 at zeta(k)
  ! X    : Location at which to evaluate the i-th Lagrange polynomial
  !        in coordinate XI
  ! Y    : Location at which to evaluate the j-th Lagrange polynomial in
  !        coordinate ETA
  ! Z    : Location at which to evaluate the k-th Lagrange polynomial in
  !        coordinate ZETA
  ! XI   : Array of all the interpolation points in XI   coordinate
  ! ETA  : Array of all the interpolation points in ETA  coordinate
  ! ZETA : Array of all the interpolation points in ZETA coordinate
  !
  ! This evaluates the 3D Lagrange polynomial at the location x,y,z.
  ! The 3D Lagrange polynomial corresponds to a function that
  ! is zero at all iterpolation points [xi(a),eta(b),zeta(c)] for
  ! [a=1,size(xi);b=1,size(eta),c=1,size(zeta)] except at
  ! [xi(i),eta(j),zeta(k)] where the function is equal to 1.
  !
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                intent(in) :: i,j,k
  real(lp),               intent(in) :: x,y,z
  real(lp), dimension(:), intent(in) :: xi
  real(lp), dimension(:), intent(in) :: eta
  real(lp), dimension(:), intent(in) :: zeta
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
continue
  !
  return_value = eval_LagrangePoly1D(i,x,xi) * &
                 eval_LagrangePoly1D(j,y,eta) * &
                 eval_LagrangePoly1D(k,z,zeta)
  !
end function eval_LagrangePoly3D_QP
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
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                intent(in) :: m
  real(lp),               intent(in) :: x
  real(lp), dimension(:), intent(in) :: xi
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
  !.. Local Scalars ..
  integer :: l,n,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
  logical(lk), dimension(1:size(xi)) :: lmask
  !
  !.. Local Constants ..
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
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                intent(in) :: m
  real(lp),               intent(in) :: x
  real(lp), dimension(:), intent(in) :: xi
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
  !.. Local Scalars ..
  integer :: l,n,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
  logical(lk), dimension(1:size(xi)) :: lmask
  !
  !.. Local Constants ..
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
  !.. Function Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,                intent(in) :: m
  real(lp),               intent(in) :: x
  real(lp), dimension(:), intent(in) :: xi
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
  !.. Local Scalars ..
  integer :: l,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
  logical(lk), dimension(1:size(xi)) :: lmask
  !
  !.. Local Constants ..
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
  !.. Function Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,                intent(in) :: m
  real(lp),               intent(in) :: x
  real(lp), dimension(:), intent(in) :: xi
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
  !.. Local Scalars ..
  integer :: l,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
  logical(lk), dimension(1:size(xi)) :: lmask
  !
  !.. Local Constants ..
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
end module polynomial_mod
