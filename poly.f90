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
 !include "kind_defs.h"
 !include "constants.h"
  !
  public :: JacobiGL
  !
  public :: eval_jacobi_poly
  public :: eval_legendre_poly
  public :: eval_chebyshev_poly
  public :: eval_normalized_jacobi_poly
  public :: normalized_jacobi_poly
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
  !
  real(vbp), parameter :: s_start = 0.426401432711220868_vbp
  real(vbp), parameter :: gamma_constant(1:11) = &
             [ -0.524741987629368444e+00_vbp,  0.116154405493589130e+00_vbp, &
               -0.765978624506602380e-02_vbp,  0.899719449391378898e-04_vbp, &
               -0.194536980009534621e-07_vbp,  0.199382839513630987e-10_vbp, &
               -0.204209590209541319e-11_vbp,  0.863896817907000175e-13_vbp, &
                0.152237501608472336e-13_vbp, -0.825725175277719950e-14_vbp, &
                0.299734782205224610e-14_vbp ]
  !
contains
!
!###############################################################################
!
subroutine JacobiGL(alpha,beta,n,xi,weights)
  !
  integer,                      intent(in) :: n
  !
  real(vbp),                    intent(in) :: alpha
  real(vbp),                    intent(in) :: beta
  !
  real(vbp), dimension(1:n+1), intent(out) :: xi
  real(vbp), dimension(1:n+1), intent(out) :: weights
  !
  real(vbp), dimension(1:n+1) :: vn
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
pure subroutine eval_legendre_poly(n,x,y,dy,d2y)
  !
  integer,    intent(in) :: n
  !
  real(vbp),  intent(in) :: x
  real(vbp), intent(out) :: y
  real(vbp), intent(out) :: dy
  real(vbp), intent(out) :: d2y
  !
  real(vbp) :: alpha,beta
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
pure function normalized_jacobi_poly(N,alpha,beta,x) result(return_value)
  !
  integer,   intent(in) :: N
  !
  real(vbp), intent(in) :: alpha
  real(vbp), intent(in) :: beta
  real(vbp), intent(in) :: x
  !
  real(vbp) :: return_value
  !
  integer   :: i
  !
  real(vbp) :: ab,ab1,ab2,ab3,amb,a1,b1
  real(vbp) :: gamma0,gamma1
  real(vbp) :: aold,anew,bnew,h1,ri
  !
  real(vbp) :: PL(N+1)
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
      ri = real(i,vbp)
      h1 = two*ri + ab
      anew = two / (h1+two) * sqrt( (ri+one)*(ri+ab1)*(ri+a1)*(ri+b1)/(h1+one)/(h1+three) )
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
subroutine eval_normalized_jacobi_poly(n,alpha,beta,x,y,dy)
  !
  integer,    intent(in) :: n
  !
  real(vbp),  intent(in) :: alpha
  real(vbp),  intent(in) :: beta
  real(vbp),  intent(in) :: x
  real(vbp), intent(out) :: y
  real(vbp), intent(out) :: dy
  !
  integer :: n1
  !
  real(vbp) :: a1,b1,ab1,rn
  !
continue
  !
  a1  = alpha + one
  b1  = beta  + one
  ab1 = alpha + beta + one
  rn  = real(n,vbp)
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
pure subroutine eval_jacobi_poly(n,alpha,beta,x,y,dy,d2y)
  !
  integer,    intent(in) :: n
  !
  real(vbp),  intent(in) :: alpha
  real(vbp),  intent(in) :: beta
  real(vbp),  intent(in) :: x
  real(vbp), intent(out) :: y
  real(vbp), intent(out) :: dy
  real(vbp), intent(out) :: d2y
  !
  integer :: i
  !
  real(vbp) :: ab,di,c0,c1,c2,c3,c4
  real(vbp) :: ym,yp,dym,dyp,d2ym,d2yp
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
      di = real(i,vbp)
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
subroutine eval_chebyshev_poly(n,x,y,dy,d2y)
  !
  integer,    intent(in) :: n
  !
  real(vbp),  intent(in) :: x
  real(vbp), intent(out) :: y
  real(vbp), intent(out) :: dy
  real(vbp), intent(out) :: d2y
  !
  integer :: k
  !
  real(vbp) :: yp,ym,dyp,dym,d2yp,d2ym
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
pure subroutine nodes_jacobi_gauss(n,alpha,beta,xi,dxi)
  !
  integer,                    intent(in) :: n
  !
  real(vbp),                  intent(in) :: alpha
  real(vbp),                  intent(in) :: beta
  real(vbp), dimension(1:n), intent(out) :: xi
  real(vbp), dimension(1:n), intent(out) :: dxi
  !
  integer :: i,it
  !
  real(vbp) :: const,ri,x,y,dy,d2y
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
    const = half*pi / (two*real(n,vbp) + alpha + beta + one)
    !
    do i = 1,n
      !
      ri = real(i,vbp)
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
      dxi(n-i+1) = -dxi(i)
    end do
    !
  end if
  !
end subroutine nodes_jacobi_gauss
!
!###############################################################################
!
pure subroutine nodes_legendre_gauss(n,xi,dxi)
  !
  integer,                    intent(in) :: n
  !
  real(vbp), dimension(1:n), intent(out) :: xi
  real(vbp), dimension(1:n), intent(out) :: dxi
  !
  real(vbp) :: alpha,beta
  !
continue
  !
  call nodes_jacobi_gauss(n,zero,zero,xi,dxi)
  !
end subroutine nodes_legendre_gauss
!
!###############################################################################
!
subroutine nodes_right_radau(n,xi)
  !
  integer,                        intent(in) :: n
  !
  real(vbp), dimension(1:n+1), intent(inout) :: xi
  !
  integer :: i,it
  real(vbp) :: x,y,dy
  !
  real(vbp), dimension(0:n+1) :: cfrr
  real(vbp), dimension(0:n) :: cdrr
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
      cdrr(i) = real(i+1,vbp) * cfrr(i+1)
    end do
    !
    xi(n+1) = one
    !
    do i = 1,n
      !
      x  = -cos( pi * real(i,vbp) / real(n+1,vbp) )
      !
      do it = 1,8
        !
        y = sum( (x**[0:n+1]) * cfrr(:) )
        dy = sum( (x**[0:n]) * cdrr(:) )
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
pure function coefrr(n) result(return_value)
  !
  integer, intent(in) :: n
  !
  real(vbp), dimension(0:n) :: return_value
  !
  integer :: k,ip
  !
  real(vbp), dimension(0:n)   :: cfk
  real(vbp), dimension(0:n)   :: cfkmt
  real(vbp), dimension(0:n-1) :: cfkm1
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
    return_value(k) = half * real(ip,vbp) * (cfk(k) - cfkmt(k))
  end do
  !
end function coefrr
!
!###############################################################################
!
pure function coefle(n) result(return_value)
  !
  integer, intent(in) :: n
  !
  real(vbp), dimension(0:n) :: return_value
  !
  integer :: i,j,k,m
  !
  real(vbp), dimension(-1:n,-1:n) :: cle
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
      cle(k,m) = (real(2*k-1,vbp)/real(k,vbp))*cle(k-1,m-1) - &
                 (real(  k-1,vbp)/real(k,vbp))*cle(k-2,m  )
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
subroutine nodes_chebyshev_gauss(n,xi)
  !
  integer,                    intent(in) :: n
  !
  real(vbp), dimension(1:n), intent(out) :: xi
  !
  real(vbp) :: alpha,beta
  real(vbp), dimension(1:n) :: dxi
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
subroutine nodes_jacobi_gauss_lobatto(n,alpha,beta,xi,vn)
  !
  integer,                     intent(in) :: n
  !
  real(vbp),                   intent(in) :: alpha
  real(vbp),                   intent(in) :: beta
  real(vbp), dimension(1:n+1), intent(out) :: xi
  real(vbp), dimension(1:n+1), intent(out) :: vn
  !
  integer :: i,it,n1
  !
  real(vbp) :: x,y,dy,d2y,ri,rn,ab,a1,b1,const
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
    rn = real(n,vbp)
    ab = alpha + beta
    a1 = alpha + one
    b1 = beta  + one
    n1 = n-1
    const = half * pi / (two*rn + ab + one)
    !
    do i = 2,n
      !
      ri = real(i-1,vbp)
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
pure subroutine nodes_legendre_gauss_lobatto(n,xi,vn)
  !
  integer,                      intent(in) :: n
  real(vbp), dimension(1:n+1), intent(out) :: xi
  real(vbp), dimension(1:n+1), intent(out) :: vn
  !
  integer :: i,it,n2
  !
  real(vbp) :: sn,x,y,dy,d2y
  !
continue
  !
  if (n == 0) then
    xi(1) = zero
    vn(1) = one
  else
    !
    n2 = (n-1)/2
    sn = real(2*n - 4*n2 - 3,vbp)
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
      x = cos( pi * real(i,vbp) / real(n,vbp) )
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
end subroutine nodes_legendre_gauss_lobatto
!
!###############################################################################
!
subroutine nodes_chebyshev_gauss_lobatto(n,xi)
  !
  integer,                      intent(in) :: n
  !
  real(vbp), dimension(1:n+1), intent(out) :: xi
  !
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
      xi(i+1)   = -cos( pi * real(i,vbp) / real(n,vbp) )
      xi(n-i+1) =  cos( pi * real(i,vbp) / real(n,vbp) )
    end do
  end if
  !
end subroutine nodes_chebyshev_gauss_lobatto
!
!###############################################################################
!
elemental function gammafun(x) result(return_value)
  !
  real(vbp), intent(in) :: x
  !
  real(vbp) :: return_value
  !
  integer   :: ind,k
  real(vbp) :: xx,pr,s,g
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
    pr = pr * (xx-real(k,vbp)) / (xx+real(k-1,vbp))
    s  = s + gamma_constant(k)*pr
  end do
  !
  g = s * exp(one-xx) * (xx+c4_50)**(xx-half)
  !
  if (ind == 1) then
    return_value = return_value/g
  else
    return_value = return_value*g
  end if
  !
end function gammafun
!
!###############################################################################
!
elemental function log_gammafun(x) result(return_value)
  real(vbp), intent(in) :: x
  real(vbp) :: return_value
  return_value = log( gammafun(x) )
end function log_gammafun
!
!###############################################################################
!
subroutine weights_jacobi_gauss(n,alpha,beta,xi,dy,weights)
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
  real(vbp),                  intent(in) :: alpha
  real(vbp),                  intent(in) :: beta
  !
  real(vbp), dimension(1:n),  intent(in) :: xi
  real(vbp), dimension(1:n),  intent(in) :: dy
  real(vbp), dimension(1:n), intent(out) :: weights
  !
  integer :: i,m
  !
  real(vbp) :: ab,a2,b2,c,ri
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
    ri = real(i,vbp)
    c  = c * (ri+alpha) * (ri+beta) / (ri*(ri+alpha+beta))
  end do
  !
  do i = 1,n
    weights(i) = c / ((one - xi(i)*xi(i)) * dy(i) * dy(i))
  end do
  !
end subroutine weights_jacobi_gauss
!
!################################################################################
!
pure subroutine weights_legendre_gauss(n,xi,dy,weights)
  !****************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS FORMULA
  !   N  = ORDER OF THE FORMULA
  !   CS = ZEROES OF THE LEGENDRE POLYNOMIAL, CS(I), I=1,N
  !   DZ = VECTOR OF THE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
  !   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
  !****************************************************************
  integer,                    intent(in) :: n
  !
  real(vbp), dimension(1:n),  intent(in) :: xi
  real(vbp), dimension(1:n),  intent(in) :: dy
  real(vbp), dimension(1:n), intent(out) :: weights
  !
  integer :: i
  !
continue
  !
  weights(:) = two
  !
  do i = 1,n
    weights(i) = two / ((one - xi(i)*xi(i)) * dy(i) * dy(i))
  end do
  !
end subroutine weights_legendre_gauss
!
!################################################################################
!
subroutine weights_chebyshev_gauss(n,weights)
  !****************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE CHEBYSHEV GAUSS FORMULA
  !   N  = ORDER OF THE FORMULA
  !   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
  !****************************************************************
  integer,                    intent(in) :: n
  !
  real(vbp), dimension(1:n), intent(out) :: weights
  !
  integer :: i
  !
continue
  !
  weights(:) = pi
  !
  do i = 1,n
    weights(i) = pi/real(n,vbp)
  end do
  !
end subroutine weights_chebyshev_gauss
!
!################################################################################
!
subroutine weights_jacobi_gauss_lobatto(n,alpha,beta,xi,weights)
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
  real(vbp),                    intent(in) :: alpha
  real(vbp),                    intent(in) :: beta
  !
  real(vbp), dimension(1:n+1),  intent(in) :: xi
  real(vbp), dimension(1:n+1), intent(out) :: weights
  !
  integer :: i
  !
  real(vbp) :: a1,a2,b1,b2,ab,ab1,ab2
  real(vbp) :: c,c1,c2,c3,ri,rn,su
  real(vbp) :: y,z,dy,d2y
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
    rn = real(n,vbp)
    c  = c * (two*rn + ab) / (rn + ab1)
    c1 = c * a1 / (b2 * ab2)
    c2 = c * b1 / (a2 * ab2)
    c3 = half * c * a1 * b1
    !
    do i = 2,n-1
      ri = real(i-1,vbp)
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
!################################################################################
!
subroutine weights_legendre_gauss_lobatto(n,vn,weights)
  !**********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N
  !   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !**********************************************************************
  integer,                      intent(in) :: n
  !
  real(vbp), dimension(1:n+1),  intent(in) :: vn
  real(vbp), dimension(1:n+1), intent(out) :: weights
  !
  integer :: i
  !
  real(vbp) :: c
  !
continue
  !
  c = two
  if (n > 0) then
    c = c / (real(n,vbp)*(real(n,vbp)+one))
  end if
  !
  do i = 1,n+1
    weights(i) = c / (vn(i)*vn(i))
  end do
  !
end subroutine weights_legendre_gauss_lobatto
!
!################################################################################
!
subroutine weights_right_radau(n,xi,weights)
  !**********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   XI = RIGHT RADAU NODES, XI(I), I=0,N
  !   VN = VALUES OF THE RIGHT RADAU POLYNOMIAL AT THE NODES, VN(I), I=0,N
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !**********************************************************************
  integer,                      intent(in) :: n
  !
  real(vbp), dimension(1:n+1),  intent(in) :: xi
  real(vbp), dimension(1:n+1), intent(out) :: weights
  !
  integer :: i
  !
  real(vbp) :: c,y,dy,d2y
  !
continue
  !
  c = two
  !
  if (n > 0) then
    !
    c = one / (real(n+1,vbp)*real(n+1,vbp))
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
!################################################################################
!
subroutine weights_chebyshev_gauss_lobatto(n,weights)
  !***********************************************************************
  !   COMPUTES THE WEIGHTS RELATIVE TO THE CHEBYSHEV GAUSS-LOBATTO FORMULA
  !   N  = ORDER OF THE FORMULA
  !   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
  !***********************************************************************
  integer,                      intent(in) :: n
  !
  real(vbp), dimension(1:n+1), intent(out) :: weights
  !
  integer :: i
  !
continue
  !
  weights(:) = pi
  !
  if (n > 0) then
    weights(1)   = half*pi/real(n,vbp)
    weights(n+1) = half*pi/real(n,vbp)
  end if
  !
  do i = 2,n
    weights(i) = pi/real(n,vbp)
  end do
  !
end subroutine weights_chebyshev_gauss_lobatto
!
!###############################################################################
!
subroutine weights_dgbook_jacobi_gauss_lobatto(n,alpha,beta,xi,weights)
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
  real(vbp),                    intent(in) :: alpha
  real(vbp),                    intent(in) :: beta
  !
  real(vbp), dimension(1:n+1),  intent(in) :: xi
  real(vbp), dimension(1:n+1), intent(out) :: weights
  !
  integer :: i
  !
  real(vbp) :: a1,a2,b1,b2,ab,ab1,ab2
  real(vbp) :: c,c1,c2,c3,ri,rn,su
  real(vbp) :: y,z,dy,d2y,scl
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
    rn = real(n,vbp)
    c  = c * (two*rn + ab) / (rn + ab1)
    c1 = c * a1 / (b2 * ab2)
    c2 = c * b1 / (a2 * ab2)
    c3 = half * c * a1 * b1
    !
    do i = 2,n-1
      ri = real(i-1,vbp)
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
!################################################################################
!
subroutine derivative_matrix_legendre_gauss(n,cs,dz,dma)
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
  real(vbp), dimension(1:n),      intent(in) :: cs
  real(vbp), dimension(1:n),      intent(in) :: dz
  real(vbp), dimension(1:n,1:n), intent(out) :: dma
  !
  !.. Local Scalars ..
  integer   :: i,j
  real(vbp) :: vi,zi,vj,zj
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
end subroutine derivative_matrix_legendre_gauss
!
!################################################################################
!
subroutine derivative_matrix_jacobi_gauss_lobatto(n,alpha,beta,et,vn,dma)
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
  real(vbp),                      intent(in) :: alpha
  real(vbp),                      intent(in) :: beta
  real(vbp), dimension(1:n),      intent(in) :: et
  real(vbp), dimension(1:n),      intent(in) :: vn
  real(vbp), dimension(1:n,1:n), intent(out) :: dma
  !
  !.. Local Scalars ..
  integer   :: i,j
  real(vbp) :: a1,b1,ab,dn,ei,ej
  real(vbp) :: c1,c2,c3,c4,vi,vj
  !
continue
  !
  dma(1,1) = zero
  if (n == 1) return
  !
  a1 = alpha+one
  b1 = beta+one
  ab = alpha+beta
  dn = real(n-1,vbp)
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
end subroutine derivative_matrix_jacobi_gauss_lobatto
!
!################################################################################
!
subroutine derivative_matrix_legendre_gauss_lobatto(n,et,vn,dma)
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
  real(vbp), dimension(1:n),      intent(in) :: et
  real(vbp), dimension(1:n),      intent(in) :: vn
  real(vbp), dimension(1:n,1:n), intent(out) :: dma
  !
  !.. Local Scalars ..
  integer   :: i,j
  real(vbp) :: ei,ej,vi,vj,dn,c
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
    dn = real(n-1,vbp)
    c  = c0_25*dn*(dn+one)
    dma(1,1) = -c
    dma(n,n) = c
  end if
  !
end subroutine derivative_matrix_legendre_gauss_lobatto
!
!################################################################################
!
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
  integer,                        intent(in) :: n
  real(vbp), dimension(1:n),      intent(in) :: et
  real(vbp), dimension(1:n,1:n), intent(out) :: dma
  !
  !.. Local Scalars ..
  integer   :: i,j
  real(vbp) :: dn,cn,ej,ei
  real(vbp) :: sn,sgni,sgnj
  !
continue
  !
  dma(1,1) = zero
  if (n == 1) return
  !
  dn = real((n-1),vbp)
  cn = (two*dn*dn+one)/six
  sn = real(1+4*((n-1)/2)-2*(n-1),vbp)
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
end subroutine derivative_matrix_chebyshev_gauss_lobatto
!
!################################################################################
!
!
!################################################################################
!
!
!################################################################################
!
!
!################################################################################
!
!
!################################################################################
!
!
!################################################################################
!
!
!################################################################################
!
end module polynomial_mod
