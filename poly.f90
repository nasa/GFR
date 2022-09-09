module polynomial_mod
  !
!!!#define DISABLE_PURE
  !
  ! Module that contains various functions and routines
  ! related to Jacobi, Legendre, and other polynomials
  !
  use module_kind_types, only : MP,LK
  use module_kind_types, only : chop
  use module_kind_types, only : intseq
  !
  implicit none
  !
  !.. Module Precision ..
  integer, parameter :: lp = MP
  !
  private
  !
  public :: legendre_gauss_quadrature
  public :: legendre_gauss_lobatto_quadrature
  public :: equal_distant_quadrature
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
  interface eval_LagrangePoly1D
    module procedure eval_single_LagrangePoly1D, &
                     eval_all_LagrangePoly1D
  end interface eval_LagrangePoly1D
  !
  interface eval_LagrangePoly2D
    module procedure eval_single_LagrangePoly2D, &
                     eval_all_LagrangePoly2D
  end interface eval_LagrangePoly2D
  !
  interface eval_LagrangePoly3D
    module procedure eval_single_LagrangePoly3D, &
                     eval_all_LagrangePoly3D
  end interface eval_LagrangePoly3D
  !
  !.. Module Parameter Constants ..
  real(lp), parameter :: zero  = 0.0_lp
  real(lp), parameter :: one4  = 0.25_lp
  real(lp), parameter :: half  = 0.5_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  real(lp), parameter :: three = 3.0_lp
  real(lp), parameter :: four  = 4.0_lp
  real(lp), parameter :: six   = 6.0_lp
  !
  real(lp), parameter :: eps17 = 1.0e-17_lp
  !
  real(lp), parameter :: pi    = 3.14159265358979323846_lp
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
#ifndef DISABLE_PURE
pure &
#endif
subroutine legendre_gauss_quadrature(xi,weights)
  !
  !.. Formal Arguments ..
  real(lp), intent(inout) :: xi(:)
  !
  !.. Optional Arguments ..
  real(lp), optional, intent(inout) :: weights(:)
  !
  !.. Local Scalars ..
  integer :: npts
  !
  !.. Local Arrays ..
  real(lp) :: pts(1:size(xi))
  real(lp) :: wts(1:size(xi))
  real(lp) :: val(1:size(xi))
  !
continue
  !
  npts = size(xi)
  !
  call nodes_legendre_gauss(npts,pts,val)
  call weights_legendre_gauss(npts,pts,val,wts)
  !
  xi = real( chop(pts) , kind=lp )
  !
  if (present(weights)) then
    weights = real( chop(wts) , kind=lp )
  end if
  !
end subroutine legendre_gauss_quadrature
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine legendre_gauss_lobatto_quadrature(xi,weights)
  !
  !.. Formal Arguments ..
  real(lp), intent(inout) :: xi(:)
  !
  !.. Optional Arguments ..
  real(lp), optional, intent(inout) :: weights(:)
  !
  !.. Local Scalars ..
  integer :: npts
  !
  !.. Local Arrays ..
  real(lp) :: pts(1:size(xi))
  real(lp) :: wts(1:size(xi))
  real(lp) :: val(1:size(xi))
  !
continue
  !
  npts = size(xi)
  !
  call nodes_legendre_gauss_lobatto(npts-1,pts,val)
  call weights_legendre_gauss_lobatto(npts-1,val,wts)
  !
  xi = real( chop(pts) , kind=lp )
  !
  if (present(weights)) then
    weights = real( chop(wts) , kind=lp )
  end if
  !
end subroutine legendre_gauss_lobatto_quadrature
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine equal_distant_quadrature(xi,weights)
  !
  !.. Formal Arguments ..
  real(lp), intent(inout) :: xi(:)
  !
  !.. Optional Arguments ..
  real(lp), optional, intent(inout) :: weights(:)
  !
  !.. Local Scalars ..
  integer :: p,k
  !
  !.. Local Arrays ..
  real(lp) :: pts(1:size(xi))
  real(lp) :: wts(1:size(xi))
  !
  real(lp) :: ga_pts(1:size(xi)+1)
  real(lp) :: ga_wts(1:size(xi)+1)
  !
continue
  !
  p = size(xi)-1
  !
  pts(:) = real(intseq(-p,p,2),kind=lp) / real(max(p,1),kind=lp)
  !
  call legendre_gauss_quadrature(ga_pts,ga_wts)
  !
  ! Each quadrature weight is equal to the integral of the Lagrange
  ! polynomial for the respective quadrature point.
  !
  ! Integrate the Lagrange polynomials for the equal distant points
  ! using a Gauss quadrature of order p+1. Each Lagrange polynomial is
  ! of order p, and integration using Gauss quadrature is exact to 2p+1,
  ! so p+1 Gauss quadrature is exact to 2p+3. Using p+1 is just being
  ! extra. For 1D, this extra-ness isnt computationally expensive.
  !
  wts(:) = zero
  do k = 1,size(ga_pts)
    wts(:) = wts(:) + ga_wts(k)*eval_all_LagrangePoly1D(ga_pts(k),pts)
  end do
  !
  xi = real( chop(pts) , kind=lp )
  !
  if (present(weights)) then
    weights = real( chop(wts) , kind=lp )
  end if
  !
end subroutine equal_distant_quadrature
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
#ifndef DISABLE_PURE
pure &
#endif
subroutine JacobiGL(alpha,beta,n,xi,weights)
  !
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp),  intent(in) :: alpha
  real(lp),  intent(in) :: beta
  real(lp), intent(out) :: xi(1:n+1)
  real(lp), intent(out) :: weights(1:n+1)
  !
  !.. Local Arrays ..
  real(lp) :: vn(1:n+1)
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
end subroutine JacobiGL
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine eval_jacobi_poly(n,alpha,beta,x,y,dy,d2y)
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
end subroutine eval_jacobi_poly
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine eval_legendre_poly(n,x,y,dy,d2y)
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
continue
  !
  alpha = zero
  beta  = zero
  !
  call eval_jacobi_poly(n,alpha,beta,x,y,dy,d2y)
  !
end subroutine eval_legendre_poly
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine eval_chebyshev_poly(n,x,y,dy,d2y)
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
end subroutine eval_chebyshev_poly
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine eval_normalized_jacobi_poly(n,alpha,beta,x,y,dy)
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
end subroutine eval_normalized_jacobi_poly
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function normalized_jacobi_poly(N,alpha,beta,x) result(return_value)
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
  gamma0 = two**ab1 / ab1 * gamma(a1) * gamma(b1) / gamma(ab1)
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
end function normalized_jacobi_poly
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function grad_normalized_jacobi_poly(n,alpha,beta,x) &
                                      result(return_value)
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
end function grad_normalized_jacobi_poly
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine nodes_jacobi_gauss(n,alpha,beta,xi,dxi)
  !
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp),  intent(in) :: alpha
  real(lp),  intent(in) :: beta
  real(lp), intent(out) :: xi(1:n)
  !
  !.. Optional Arguments ..
  real(lp), optional, intent(out) :: dxi(1:n)
  !
  !.. Local Scalars ..
  integer  :: i,it
  real(lp) :: const,ri,x,y,dy,d2y
  !
  !.. Local Arrays ..
  real(lp) :: val(1:n)
  !
continue
  !
  if (n > 0) then
    xi(1)  = (beta - alpha) / (alpha + beta + two)
    val(1) = half*(alpha + beta) + one
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
      val(i) = dy
      !
    end do
    !
    do i = 1,n/2
      xi (n-i+1) = -xi (i)
     !val(n-i+1) = -val(i)
    end do
    !
  end if
  !
  if (present(dxi)) then
    dxi(:) = val(:)
  end if
  !
end subroutine nodes_jacobi_gauss
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine nodes_legendre_gauss(n,xi,dxi)
  !
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp), intent(out) :: xi(1:n)
  !
  !.. Optional Arguments ..
  real(lp), optional, intent(out) :: dxi(1:n)
  !
  !.. Local Scalars ..
  real(lp) :: alpha,beta
  !
continue
  !
  call nodes_jacobi_gauss(n,zero,zero,xi,dxi)
#ifdef DISABLE_PURE
  call check_for_bad_number(__LINE__,__FILE__)
#endif
  !
end subroutine nodes_legendre_gauss
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine nodes_right_radau(n,xi)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n
  real(lp), intent(inout) :: xi(1:n+1)
  !
  !.. Local Scalars ..
  integer  :: i,it
  real(lp) :: x,y,dy
  !
  !.. Local Arrays ..
  real(lp) :: cfrr(0:n+1)
  real(lp) :: cdrr(0:n)
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
    cfrr(:) = coefrr(n+1)
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
end subroutine nodes_right_radau
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine nodes_chebyshev_gauss(n,xi)
  !
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp), intent(out) :: xi(1:n)
  !
  !.. Local Scalars ..
  real(lp) :: alpha,beta
  !
  !.. Local Arrays ..
  real(lp) :: dxi(1:n)
  !
continue
  !
  alpha = -half
  beta  = -half
  !
  call nodes_jacobi_gauss(n,alpha,beta,xi,dxi)
  !
end subroutine nodes_chebyshev_gauss
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine nodes_jacobi_gauss_lobatto(n,alpha,beta,xi,vn)
  !
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp),  intent(in) :: alpha
  real(lp),  intent(in) :: beta
  real(lp), intent(out) :: xi(1:n+1)
  real(lp), intent(out) :: vn(1:n+1)
  !
  !.. Local Scalars ..
  integer  :: i,it,n1
  real(lp) :: x,y,dy,d2y,ri,rn,ab,a1,b1,const
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
end subroutine nodes_jacobi_gauss_lobatto
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine nodes_legendre_gauss_lobatto(n,xi,vn)
  !
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp), intent(out) :: xi(1:n+1)
  !
  !.. Optional Arguments ..
  real(lp), optional, intent(out) :: vn(1:n+1)
  !
  !.. Local Scalars ..
  integer  :: i,it,n2
  real(lp) :: sn,x,y,dy,d2y
  !
  !.. Local Arrays ..
  real(lp) :: val(1:n+1)
  !
continue
  !
  if (n == 0) then
    xi(1) = zero
    val(1) = one
  else
    !
    n2 = (n-1)/2
    sn = real(2*n - 4*n2 - 3,kind=lp)
    !
    xi(1) = -one
    val(1) =  sn
    !
    xi(n+1) = one
    val(n+1) = one
    !
  end if
  !
  if (n > 1) then
    x = zero
    call eval_legendre_poly(n,x,y,dy,d2y)
    xi(n2+2) = x
    val(n2+2) = y
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
      val(i+1)   = y*sn
      val(n-i+1) = y
      !
    end do
    !
  end if
  !
  if (present(vn)) then
    vn(:) = val(:)
  end if
  !
end subroutine nodes_legendre_gauss_lobatto
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine nodes_chebyshev_gauss_lobatto(n,xi)
  !
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp), intent(out) :: xi(1:n+1)
  !
  !.. Local Scalars ..
  integer :: i,n2
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
end subroutine nodes_chebyshev_gauss_lobatto
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine weights_jacobi_gauss(n,alpha,beta,xi,dy,weights)
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
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp),  intent(in) :: alpha
  real(lp),  intent(in) :: beta
  real(lp),  intent(in) :: xi(1:n)
  real(lp),  intent(in) :: dy(1:n)
  real(lp), intent(out) :: weights(1:n)
  !
  !.. Local Scalars ..
  integer  :: i,m
  real(lp) :: ab,a2,b2,c,ri
  !
continue
  !
  weights(:) = one
  !
  ab = alpha + beta + two
  a2 = alpha + two
  b2 = beta  + two
  !
  c = half * (two**ab) * gamma(a2) * gamma(b2) / gamma(ab)
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
end subroutine weights_jacobi_gauss
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine weights_legendre_gauss(n,xi,dy,weights)
  !
  !****************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS FORMULA
  !   N  = ORDER OF THE FORMULA
  !   CS = ZEROES OF THE LEGENDRE POLYNOMIAL, CS(I), I=1,N
  !   DZ = VECTOR OF THE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
  !   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
  !****************************************************************
  !
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp),  intent(in) :: xi(1:n)
  real(lp),  intent(in) :: dy(1:n)
  real(lp), intent(out) :: weights(1:n)
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Arrays ..
  real(lp) :: dy2_inv(1:n)
  real(lp) :: dy2(1:n)
  !
continue
  !
  weights(:) = two
  !
  do i = 1,n
    dy2(i) = dy(i)*dy(i)
  end do
#ifdef DISABLE_PURE
  call check_for_bad_number(__LINE__,__FILE__)
#endif
  !
  do i = 1,n
    dy2_inv(i) = one / dy2(i)
  end do
#ifdef DISABLE_PURE
  call check_for_bad_number(__LINE__,__FILE__)
#endif
  !
  do i = 1,n
    weights(i) = two * dy2_inv(i) / (one - xi(i)*xi(i))
  end do
#ifdef DISABLE_PURE
  call check_for_bad_number(__LINE__,__FILE__)
#endif
  !
! do i = 1,n
!   weights(i) = two / ((one - (xi(i)**2)) * (dy(i)**2) )
!#ifdef DISABLE_PURE
!   if (mypnum == glb_root) then
!     write (3000+mypnum,100) i, xi(i), dy(i), &
!                      one-xi(i)*xi(i), &
!                     (one-xi(i)*xi(i))*dy(i)*dy(i), &
!                     weights(i)
!     100 format ("i = ",i0," ; xi = ",es14.6," ; dy = ",es14.6, &
!                 " ; 1-xi*xi = ",es14.6," ; (1-xi*xi)*dy*dy = ",es14.6, &
!                 " ; weight = ",es14.6)
!     flush (3000+mypnum)
!   end if
!#endif
! end do
!#ifdef DISABLE_PURE
! if (mypnum == glb_root) then
!   write (3000+mypnum,*)
!   flush (3000+mypnum)
! end if
! call mpi_barrier(MPI_COMM_WORLD,mpierr)
! call check_for_bad_number(__LINE__,__FILE__)
!#endif
  !
end subroutine weights_legendre_gauss
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine weights_right_radau(n,xi,weights)
  !
  !**********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   XI = RIGHT RADAU NODES, XI(I), I=0,N
  !   VN = VALUES OF THE RIGHT RADAU POLYNOMIAL AT THE NODES, VN(I), I=0,N
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !**********************************************************************
  !
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp),  intent(in) :: xi(1:n+1)
  real(lp), intent(out) :: weights(1:n+1)
  !
  !.. Local Scalars ..
  integer  :: i
  real(lp) :: c,y,dy,d2y
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
end subroutine weights_right_radau
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine weights_chebyshev_gauss(n,weights)
  !
  !****************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE CHEBYSHEV GAUSS FORMULA
  !   N  = ORDER OF THE FORMULA
  !   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
  !****************************************************************
  !
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp), intent(out) :: weights(1:n)
  !
  !.. Local Scalars ..
  integer :: i
  !
continue
  !
  weights(:) = pi
  !
  do i = 1,n
    weights(i) = pi/real(n,kind=lp)
  end do
  !
end subroutine weights_chebyshev_gauss
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine weights_jacobi_gauss_lobatto(n,alpha,beta,xi,weights)
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
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp),  intent(in) :: alpha
  real(lp),  intent(in) :: beta
  real(lp),  intent(in) :: xi(1:n+1)
  real(lp), intent(out) :: weights(1:n+1)
  !
  !.. Local Scalars ..
  integer  :: i
  real(lp) :: a1,a2,b1,b2,ab,ab1,ab2
  real(lp) :: c,c1,c2,c3,ri,rn,su
  real(lp) :: y,z,dy,d2y
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
  c = (two**ab) * gamma(a1) * gamma(b1) / gamma(ab2)
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
end subroutine weights_jacobi_gauss_lobatto
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine weights_legendre_gauss_lobatto(n,vn,weights)
  !
  !**********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N
  !   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !**********************************************************************
  !
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp),  intent(in) :: vn(1:n+1)
  real(lp), intent(out) :: weights(1:n+1)
  !
  !.. Local Scalars ..
  integer  :: i
  real(lp) :: c
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
end subroutine weights_legendre_gauss_lobatto
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine weights_chebyshev_gauss_lobatto(n,weights)
  !
  !***********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE CHEBYSHEV GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !***********************************************************************
  !
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp), intent(out) :: weights(1:n+1)
  !
  !.. Local Scalars ..
  integer :: i
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
end subroutine weights_chebyshev_gauss_lobatto
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine weights_dgbook_jacobi_gauss_lobatto(n,alpha,beta,xi,weights)
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
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp),  intent(in) :: alpha
  real(lp),  intent(in) :: beta
  real(lp),  intent(in) :: xi(1:n+1)
  real(lp), intent(out) :: weights(1:n+1)
  !
  !.. Local Scalars ..
  integer  :: i
  real(lp) :: a1,a2,b1,b2,ab,ab1,ab2
  real(lp) :: c,c1,c2,c3,ri,rn,su
  real(lp) :: y,z,dy,d2y,scl
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
    c = (two**ab) * gamma(a1) * gamma(b1) / gamma(ab2)
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
end subroutine weights_dgbook_jacobi_gauss_lobatto
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function derivative_matrix_legendre_gauss(n) result(return_value)
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
  !.. Formal Arguments ..
  integer, intent(in) :: n ! ORDER OF LEGENDRE GAUSS POLYNOMIAL
  !
  !.. Function Result ..
  real(lp) :: return_value(1:n+1,1:n+1) ! DERIVATIVE MATRIX
  !
  !.. Local Scalars ..
  integer  :: i,j,np
  !
  !.. Local Arrays ..
  real(lp) :: cs(1:n+1) ! VECTOR OF THE ZEROES
  real(lp) :: dz(1:n+1) ! LEGENDRE DERIVATIVES AT THE ZEROES
  !
  real(lp) :: val(1:n+1,1:n+1)
  real(lp) :: val_inv(1:n+1,1:n+1)
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
    do j = 1,np
      if (i /= j) then
        val(i,j) = dz(j)*(cs(i)-cs(j))
      else
        val(i,j) = one-cs(i)*cs(i)
      end if
    end do
  end do
#ifdef DISABLE_PURE
  call check_for_bad_number(__LINE__,__FILE__)
#endif
  !
  do i = 1,np
    do j = 1,np
      val_inv(i,j) = one / val(i,j)
    end do
  end do
#ifdef DISABLE_PURE
  call check_for_bad_number(__LINE__,__FILE__)
#endif
  !
  do i = 1,np
    do j = 1,np
      if (i /= j) then
        return_value(i,j) = dz(i) * val_inv(i,j)
      else
        return_value(i,j) = cs(i) * val_inv(i,j)
      end if
    end do
  end do
#ifdef DISABLE_PURE
  call check_for_bad_number(__LINE__,__FILE__)
#endif
  !
 !do i = 1,np
 !  do j = 1,np
 !    if (i /= j) then
 !      return_value(i,j) = dz(i)/(dz(j)*(cs(i)-cs(j)))
 !    else
 !      return_value(i,j) = cs(i)/(one-cs(i)*cs(i))
 !    end if
 !  end do
 !end do
!#ifdef DISABLE_PURE
 !call check_for_bad_number(__LINE__,__FILE__)
!#endif
  !
end function derivative_matrix_legendre_gauss
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function derivative_matrix_legendre_gauss_lobatto(n) result(return_value)
  !
  !***********************************************************************
  !  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE
  !  LEGENDRE GAUSS-LOBATTO NODES
  !***********************************************************************
  !
  !.. Formal Arguments ..
  integer, intent(in) :: n ! ORDER OF LEGENDRE GAUSS-LOBATTO POLYNOMIAL
  !
  !.. Function Result ..
  real(lp) :: return_value(1:n+1,1:n+1) ! DERIVATIVE MATRIX
  !
  !.. Local Scalars ..
  integer  :: i,j,np
  real(lp) :: dn,c
  !
  !.. Local Arrays ..
  real(lp) :: et(1:n+1) ! VECTOR OF THE NODES
  real(lp) :: vn(1:n+1) ! VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES
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
    do j = 1,np
      if (i /= j) then
        return_value(i,j) = vn(i)/(vn(j)*(et(i)-et(j)))
      else
        return_value(i,j) = zero
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
end function derivative_matrix_legendre_gauss_lobatto
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine derivative_matrix_jacobi_gauss_lobatto(n,alpha,beta, &
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
  !.. Formal Arguments ..
  integer,   intent(in) :: n
  real(lp),  intent(in) :: alpha
  real(lp),  intent(in) :: beta
  real(lp),  intent(in) :: et(1:n)
  real(lp),  intent(in) :: vn(1:n)
  real(lp), intent(out) :: dma(1:n,1:n)
  !
  !.. Local Scalars ..
  integer  :: i,j
  real(lp) :: a1,b1,ab,dn
  real(lp) :: c1,c2,c3,c4
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
    dma(1,j) = -c3/(vn(j)*(one+et(j)))
    dma(n,j) = c4/(vn(j)*(one-et(j)))
    dma(j,1) = vn(j)/(c3*(one+et(j)))
    dma(j,n) = -vn(j)/(c4*(one-et(j)))
  end do
  !
  do i = 2,n-1
    do j = 2,n-1
      if (i /= j) then
        dma(i,j) = vn(i)/(vn(j)*(et(i)-et(j)))
      else
        dma(i,j) = half*(ab*et(i)+alpha-beta)/(one-et(i)*et(i))
      end if
    end do
  end do
  !
end subroutine derivative_matrix_jacobi_gauss_lobatto
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine derivative_matrix_chebyshev_gauss_lobatto(n,et,dma)
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
  integer,   intent(in) :: n
  real(lp),  intent(in) :: et(1:n)
  real(lp), intent(out) :: dma(1:n,1:n)
  !
  !.. Local Scalars ..
  integer  :: i,j
  real(lp) :: dn,cn
  real(lp) :: sn,sgni,sgnj
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
    dma(1,j) = -two*sgnj/(one+et(j))
    dma(n,j) = two*sgnj*sn/(one-et(j))
    dma(j,1) = half*sgnj/(one+et(j))
    dma(j,n) = -half*sgnj*sn/(one-et(j))
    sgnj = -sgnj
  end do
  !
  sgni = -one
  do i = 2,n-1
    sgnj = -one
    do j = 2,n-1
      if (i /= j) then
        dma(i,j) = sgni*sgnj/(et(i)-et(j))
      else
        dma(i,j) = -half*et(i)/(one-et(i)*et(i))
      end if
      sgnj = -sgnj
    end do
    sgni = -sgni
  end do
  !
end subroutine derivative_matrix_chebyshev_gauss_lobatto
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function coefrr(n) result(return_value)
  !
  ! Function to compute the coefficients of
  ! the right Radau polynomial of degree n
  !
  !.. Formal Arguments ..
  integer, intent(in) :: n
  !
  !.. Function Result ..
  real(lp) :: return_value(0:n)
  !
  !.. Local Scalars ..
  integer :: k,ip
  !
  !.. Local Arrays ..
  real(lp)   :: cfk(0:n)
  real(lp)   :: cfkmt(0:n)
  real(lp) :: cfkm1(0:n-1)
  !
continue
  !
  cfk = coefle(n)
  cfkm1 = coefle(n-1)
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
end function coefrr
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function coefle(n) result(return_value)
  !
  ! Function to compute the coefficients of
  ! the Legendre polynomial of degree n
  !
  !.. Formal Arguments ..
  integer, intent(in) :: n
  !
  !.. Function Result ..
  real(lp) :: return_value(0:n)
  !
  !.. Local Scalars ..
  integer :: i,j,k,m
  !
  !.. Local Arrays ..
  real(lp) :: cle(-1:n,-1:n)
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
end function coefle
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function eval_LagrangePoly(ijk,xyz,xi) result(return_value)
  !
  ! Get the value of a tensor product of Lagrange polynomials
  ! of arbitrary dimension
  !
  ! IJK : The value of the Lagrange polynomial tensor product is 1 at ijk
  ! XYZ : Location at which to evaluate the Lagrange polynomial
  ! XI  : Array of the 1D solution points
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: ijk(:)
  real(lp), intent(in) :: xyz(:)
  real(lp), intent(in) :: xi(:)
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
  !.. Local Scalars ..
  integer :: i
  !
continue
  !
  return_value = one
  do i = 1,size(ijk)
    return_value = return_value * eval_LagrangePoly1D(ijk(i),xyz(i),xi)
  end do
  !
end function eval_LagrangePoly
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function eval_single_LagrangePoly1D(k,x,xi) result(return_value)
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
  integer,  intent(in) :: k
  real(lp), intent(in) :: x
  real(lp), intent(in) :: xi(:)
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,nfp
  !
  !.. Local Arrays ..
  integer :: imask(1:size(xi))
  !
continue
  !
  nfp = size(xi)
  imask(:) = (/ (i,i=1,nfp) /)
  !
  return_value = product( x    -xi(1:nfp) , imask /= k ) / &
                 product( xi(k)-xi(1:nfp) , imask /= k )
  !
end function eval_single_LagrangePoly1D
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function eval_all_LagrangePoly1D(x,xi) result(return_value)
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
  !.. Formal Arguments ..
  real(lp), intent(in) :: x
  real(lp), intent(in) :: xi(:)
  !
  !.. Function Result ..
  real(lp) :: return_value(1:size(xi))
  !
  !.. Local Scalars ..
  integer :: i,k,nfp
  !
  !.. Local Arrays ..
  integer :: imask(1:size(xi))
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
end function eval_all_LagrangePoly1D
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function eval_all_LagrangePoly2D(ij,xy,xi,eta) result(return_value)
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
  !.. Formal Arguments ..
  integer,  intent(in) :: ij(:)
  real(lp), intent(in) :: xy(:)
  real(lp), intent(in) :: xi(:)
  real(lp), intent(in) :: eta(:)
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
continue
  !
  return_value = eval_LagrangePoly1D(ij(1),xy(1),xi) * &
                 eval_LagrangePoly1D(ij(2),xy(2),eta)
  !
end function eval_all_LagrangePoly2D
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function eval_single_LagrangePoly2D(k,l,x,y,xi,eta) result(return_value)
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
  !.. Formal Arguments ..
  integer,  intent(in) :: k
  integer,  intent(in) :: l
  real(lp), intent(in) :: x
  real(lp), intent(in) :: y
  real(lp), intent(in) :: xi(:)
  real(lp), intent(in) :: eta(:)
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
continue
  !
  return_value = eval_LagrangePoly1D(k,x,xi) * eval_LagrangePoly1D(l,y,eta)
  !
end function eval_single_LagrangePoly2D
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function eval_all_LagrangePoly3D(ijk,xyz,xi,eta,zeta) &
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
  !.. Formal Arguments ..
  integer,  intent(in) :: ijk(:)
  real(lp), intent(in) :: xyz(:)
  real(lp), intent(in) :: xi(:)
  real(lp), intent(in) :: eta(:)
  real(lp), intent(in) :: zeta(:)
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
end function eval_all_LagrangePoly3D
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function eval_single_LagrangePoly3D(i,j,k,x,y,z,xi,eta,zeta) &
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
  !.. Formal Arguments ..
  integer,  intent(in) :: i,j,k
  real(lp), intent(in) :: x,y,z
  real(lp), intent(in) :: xi(:)
  real(lp), intent(in) :: eta(:)
  real(lp), intent(in) :: zeta(:)
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
end function eval_single_LagrangePoly3D
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function eval_D2LagrangeDx2(m,x,xi) result(return_value)
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
  integer,  intent(in) :: m
  real(lp), intent(in) :: x
  real(lp), intent(in) :: xi(:)
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
  !.. Local Scalars ..
  integer :: l,n,nfp
  !
  !.. Local Arrays ..
  integer :: imask(1:size(xi))
  logical(lk) :: lmask(1:size(xi))
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
end function eval_D2LagrangeDx2
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function eval_DLagrangeDx(m,x,xi) result(return_value)
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
  integer,  intent(in) :: m
  real(lp), intent(in) :: x
  real(lp), intent(in) :: xi(:)
  !
  !.. Function Result ..
  real(lp) :: return_value
  !
  !.. Local Scalars ..
  integer :: l,nfp
  !
  !.. Local Arrays ..
  integer :: imask(1:size(xi))
  logical(lk) :: lmask(1:size(xi))
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
end function eval_DLagrangeDx
!
!###############################################################################
!
#include "Functions/check_for_bad_number.f90"
!
!###############################################################################
!
end module polynomial_mod
