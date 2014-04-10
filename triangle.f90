module triangle_mod
  !
  ! Module that contains various functions and routines related to
  ! creating the triangular elements and the solution nodes within
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
  public :: TriNodes2D_DG
  public :: TriNodes2D_HT
  !
  real(vbp), parameter :: csq3  = sqrt(three)
  real(vbp), parameter :: csq3i = one/csq3
  real(vbp), parameter :: cxt   = cos( two*pi/three)
  real(vbp), parameter :: cxr   = cos(four*pi/three)
  real(vbp), parameter :: cyt   = sin( two*pi/three)
  real(vbp), parameter :: cyr   = sin(four*pi/three)
  !
  real(vbp), dimension(1:15), parameter :: &
      alpopt = (/ 0.0000_vbp, 0.0000_vbp, 1.4152_vbp, 0.1001_vbp, 0.2751_vbp, &
                  0.9800_vbp, 1.0999_vbp, 1.2832_vbp, 1.3648_vbp, 1.4773_vbp, &
                  1.4959_vbp, 1.5743_vbp, 1.5770_vbp, 1.6223_vbp, 1.6258_vbp  /)
  !
contains
!
!###############################################################################
!
pure function EquiTriNodes2D(N) result(return_value)
  !
  integer, intent(in) :: N
  !
  real(vbp), dimension(1:2,1:(N+1)*(N+2)/2) :: return_value
  !
  integer :: i,j,k,np
  !
  real(vbp) :: rN,ri,rj,r,s,t
  !
continue
  !
  rN = max(real(N,vbp),eps10)
  np = (N+1)*(N+2)/2
  !
  k = 1
  do i = 1,N+1
    !
    ri = real(i,vbp)
    !
    do j = 1,N+2-i
      !
      rj = real(j,vbp)
      !
      r = (rj-one)/rN
      s = (ri-one)/rN
      t = one - r - s
      !
      return_value(1,k) = r - t
      return_value(2,k) = (two*s - r - t)*csq3i
      !
      k = k + 1
      !
    end do
    !
  end do
  !
end function EquiTriNodes2D
!
!###############################################################################
!
pure function Warpfactor(N,r) result(return_value)
  !
  use startup_mod,    only : invert_matrix, &
                             Vandermonde1D
  use polynomial_mod, only : nodes_legendre_gauss_lobatto, &
                             JacobiP => normalized_jacobi_poly
  !
  integer, intent(in) :: N
  real(vbp), dimension(1:(N+1)*(N+2)/2), intent(in) :: r
  !
  real(vbp), dimension(1:(N+1)*(N+2)/2) :: return_value
  !
  !.. Local Scalars ..
  integer   :: i,j,np,nn
  real(vbp) :: sf
  !
  !.. Local Arrays ..
  real(vbp), dimension(1:N+1) :: LGL_nodes
  real(vbp), dimension(1:N+1) :: EQD_nodes
  real(vbp), dimension(1:(N+1)*(N+2)/2) :: warp
  real(vbp), dimension(1:N+1,1:N+1) :: V_EQD
  real(vbp), dimension(1:N+1,1:(N+1)*(N+2)/2) :: Pmat
  real(vbp), dimension(1:N+1,1:(N+1)*(N+2)/2) :: Lmat
  !
continue
  !
  np = N+1
  nn = (N+1)*(N+2)/2
  !
  ! Get the 1D Legendre-Gauss-Lobatto nodes
  !
  call nodes_legendre_gauss_lobatto(N,LGL_nodes,EQD_nodes)
  !
  ! Get the 1D Equidistant nodes
  !
  do i = 1,np
    EQD_nodes(i) = two*real(i-1,vbp)/real(N,vbp) - one
  end do
  !
  ! Compute the transpose of the Vandermonde matrix based on EQD_nodes
  !
  V_EQD = transpose(Vandermonde1D(EQD_nodes))
  !
  ! Evaluate Lagrange polynomial at r
  !
  do i = 1,np
    do j = 1,nn
      Pmat(i,j) = JacobiP(i-1,zero,zero,r(j))
    end do
  end do
  !
  ! Compute the Lagrange polynomials based on equidistant nodes
  ! solving the equation in the middle of page 47
  !
  Lmat = matmul(invert_matrix(V_EQD),Pmat)
  !
  ! Compute warp factor, Eq 6.7 on page 175
  !
  warp = matmul(transpose(Lmat),(LGL_nodes-EQD_nodes))
  !
  ! Scale factor, equation at bottom of page 176
  !
  do i = 1,nn
    sf = one - ( deltafun(abs(r(i))<one-eps10) * r(i) )**2
    return_value(i) = warp(i)/sf + warp(i)*( deltafun(abs(r(i))<one-eps10) - one )
  end do
  !
end function Warpfactor
!
!###############################################################################
!
pure function Nodes2D(N) result(return_value)
  !
  integer, intent(in) :: N
  real(vbp), dimension(1:2,1:(N+1)*(N+2)/2) :: return_value
  !
  !.. Local Scalars ..
  real(vbp) :: alpha
  !
  !.. Local Arrays ..
  real(vbp), dimension(1:(N+1)*(N+2)/2) :: r,warp_r
  real(vbp), dimension(1:(N+1)*(N+2)/2) :: s,warp_s
  real(vbp), dimension(1:(N+1)*(N+2)/2) :: t,warp_t
  !
continue
  !
  ! Set optimized parameter, alpha, depending on order N
  !
  if (N<16 .and. N>0) then
    alpha = alpopt(N)
  else
    alpha = five/three
  end if
  !
  return_value(:,:) = EquiTriNodes2D(N)
  !
  r(:) = ( three*return_value(1,:) - csq3*return_value(2,:) + two)/six
  s(:) = (                           csq3*return_value(2,:) + one)/three
  t(:) = (-three*return_value(1,:) - csq3*return_value(2,:) + two)/six
  !
  warp_r = Warpfactor(N,t-s)
  warp_s = Warpfactor(N,r-t)
  warp_t = Warpfactor(N,s-r)
  !
  warp_r = (four*s*t) * warp_r * (one + (alpha*r)**2)
  warp_s = (four*r*t) * warp_s * (one + (alpha*s)**2)
  warp_t = (four*r*s) * warp_t * (one + (alpha*t)**2)
  !
  return_value(1,:) = return_value(1,:) + cxr*warp_r + cxt*warp_t + warp_s
  return_value(2,:) = return_value(2,:) + cyr*warp_r + cyt*warp_t
  !
end function Nodes2D
!
!###############################################################################
!
pure function xytors(xy) result(return_value)
  !
  real(vbp), dimension(:,:), intent(in) :: xy
  real(vbp), dimension(size(xy,1),size(xy,2)) :: return_value
  !
  !.. Local Scalars ..
  integer   :: i
  real(vbp) :: r,s,t
  !
continue
  !
  do i = 1,size(xy,2)
    r = ( three*xy(1,i) - csq3*xy(2,i) + two)/six
    s = (sqrt(three)*xy(2,i) + one)/three
    t = (-three*xy(1,i) - csq3*xy(2,i) + two)/six
    return_value(1,i) = r - s - t
    return_value(2,i) = s - r - t
  end do
  !
end function xytors
!
!###############################################################################
!
pure function TriNodes2D_DG(N) result(return_value)
  !
  integer, intent(in) :: N
  !
  real(vbp), dimension(1:2,1:(N+1)*(N+2)/2) :: return_value
  !
continue
  !
  return_value(:,:) = Nodes2D(N)
  return_value(:,:) = xytors(return_value)
  !
  ! Adjust the triangle nodes to lie in the unit right triangle
  ! defined by the vertices: (0,0) - (1,0) - (0,1)
  !
  return_value(:,:) = half * (return_value(:,:) + one)
  !
end function TriNodes2D_DG
!
!###############################################################################
!
function TriNodes2D_HT(N) result(return_value)
  !
  use polynomial_mod, only : nodes_legendre_gauss_lobatto, &
                             nodes_legendre_gauss
  !
  integer, intent(in) :: N
  !
  real(vbp), dimension(1:2,1:(N+1)*(N+2)/2) :: return_value
  !
  !.. Local Scalars ..
  integer   :: i,j,k,isp,np1d,np2d,m
  real(vbp) :: r,s,t,sm
  !
  !.. Local Arrays ..
  real(vbp), dimension(1:N+1) :: xi
  real(vbp), dimension(1:N+1) :: vp
  !
continue
  !
  np1d = N+1
  np2d = (N+1)*(N+2)/2
  !
  ! Get the 1D Legendre-Gauss-Lobatto nodes
  !
  call nodes_legendre_gauss_lobatto(N,xi,vp)
 !call nodes_legendre_gauss(np1d,xi,vp)
  !
  ! Adjust the 1D LGL nodes to lie on the interval [0,1]
  !
  xi = (xi + one)/two
  !
  ! Get the triangle nodes in r,s space and defined
  ! by the vertices: (0,0) - (1,0) - (0,1)
  !
  isp = 0
  do j = 1,np1d
    do i = 1,np1d+1-j
      !
      isp = isp + 1
      k = np1d + 2 - i - j
      !
      r = xi(i)
      s = xi(j)
      t = xi(k)
      sm = r + s + t
      !
      return_value(1,isp) = r/sm
      return_value(2,isp) = s/sm
      !
    end do
  end do
  !
  ! Adjust the triangle nodes to lie in the bi-unit right
  ! triangle defined by the vertices: (-1,-1) - (1,-1) - (-1,1)
  !
 !return_value = return_value*two - one
  !
end function TriNodes2D_HT
!
!###############################################################################
!
include "Functions/deltafun.f90"
!
!###############################################################################
!
end module triangle_mod
