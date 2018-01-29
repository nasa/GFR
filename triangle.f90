module triangle_mod
  !
  ! Module that contains various functions and routines related to
  ! creating the triangular elements and the solution nodes within
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
  public :: TriNodes2D_AlphaOptimized
  public :: TriNodes2D_BarycentricLobatto
  !
  interface EquiTriNodes2D
    module procedure EquiTriNodes2D_DP, &
                     EquiTriNodes2D_QP
  end interface EquiTriNodes2D
  !
  interface Warpfactor
    module procedure Warpfactor_DP, &
                     Warpfactor_QP
  end interface Warpfactor
  !
  interface Nodes2D
    module procedure Nodes2D_DP, &
                     Nodes2D_QP
  end interface Nodes2D
  !
  interface xytors
    module procedure xytors_DP, &
                     xytors_QP
  end interface xytors
  !
  interface TriNodes2D_AlphaOptimized
    module procedure TriNodes2D_AlphaOptimized_DP, &
                     TriNodes2D_AlphaOptimized_QP
  end interface TriNodes2D_AlphaOptimized
  !
  interface TriNodes2D_BarycentricLobatto
    module procedure TriNodes2D_BarycentricLobatto_DP, &
                     TriNodes2D_BarycentricLobatto_QP
  end interface TriNodes2D_BarycentricLobatto
  !
contains
!
!###############################################################################
!
pure function EquiTriNodes2D_DP(N,val) result(return_value)
  !
  integer, parameter :: lp = dp
  !
  integer,  intent(in) :: N
  real(lp), intent(in) :: val
  !
  real(lp), dimension(1:2,1:(N+1)*(N+2)/2) :: return_value
  !
  integer :: i,j,k,np
  !
  real(lp) :: rN,ri,rj,r,s,t
  !
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  real(lp), parameter :: three = 3.0_lp
  real(lp), parameter :: csq3  = sqrt(three)
  real(lp), parameter :: csq3i = one/csq3
  real(lp), parameter :: eps10 = 1.0e-10_lp
  !
continue
  !
  rN = max(real(N,kind=lp),eps10)
  np = (N+1)*(N+2)/2
  !
  k = 1
  do i = 1,N+1
    !
    ri = real(i,kind=lp)
    !
    do j = 1,N+2-i
      !
      rj = real(j,kind=lp)
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
end function EquiTriNodes2D_DP
!
!###############################################################################
!
pure function EquiTriNodes2D_QP(N,val) result(return_value)
  !
  integer, parameter :: lp = local_qp
  !
  integer,  intent(in) :: N
  real(lp), intent(in) :: val
  !
  real(lp), dimension(1:2,1:(N+1)*(N+2)/2) :: return_value
  !
  integer :: i,j,k,np
  !
  real(lp) :: rN,ri,rj,r,s,t
  !
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  real(lp), parameter :: three = 3.0_lp
  real(lp), parameter :: csq3  = sqrt(three)
  real(lp), parameter :: csq3i = one/csq3
  real(lp), parameter :: eps10 = 1.0e-10_lp
  !
continue
  !
  rN = max(real(N,kind=lp),eps10)
  np = (N+1)*(N+2)/2
  !
  k = 1
  do i = 1,N+1
    !
    ri = real(i,kind=lp)
    !
    do j = 1,N+2-i
      !
      rj = real(j,kind=lp)
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
end function EquiTriNodes2D_QP
!
!###############################################################################
!
pure function Warpfactor_DP(N,r) result(return_value)
  !
  use polynomial_mod, only : nodes_legendre_gauss_lobatto
  use polynomial_mod, only : JacobiP => normalized_jacobi_poly
  !
  integer, parameter :: lp = dp
  !
  integer,  intent(in) :: N
  real(lp), dimension(1:(N+1)*(N+2)/2), intent(in) :: r
  !
  real(lp), dimension(1:(N+1)*(N+2)/2) :: return_value
  !
  !.. Local Scalars ..
  integer  :: i,j,np,nn
  real(lp) :: sf
  !
  !.. Local Arrays ..
  real(lp), dimension(1:N+1) :: LGL_nodes
  real(lp), dimension(1:N+1) :: EQD_nodes
  real(lp), dimension(1:(N+1)*(N+2)/2) :: warp
  real(lp), dimension(1:N+1,1:N+1) :: V_EQD
  real(lp), dimension(1:N+1,1:(N+1)*(N+2)/2) :: Pmat
  real(lp), dimension(1:N+1,1:(N+1)*(N+2)/2) :: Lmat
  !
  real(lp), parameter :: zero  = 0.0_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  real(lp), parameter :: eps10 = 1.0e-10_lp
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
    EQD_nodes(i) = two*real(i-1,kind=lp)/real(N,kind=lp) - one
  end do
  !
  ! Compute the transpose of the Vandermonde matrix based on EQD_nodes
  !
  do j = 1,np
    do i = 1,np
      V_EQD(i,j) = JacobiP(j-1,zero,zero,EQD_nodes(i))
    end do
  end do
  V_EQD = transpose(V_EQD)
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
    return_value(i) = warp(i)/sf + &
                      warp(i)*( deltafun(abs(r(i))<one-eps10) - one )
  end do
  !
end function Warpfactor_DP
!
!###############################################################################
!
pure function Warpfactor_QP(N,r) result(return_value)
  !
  use polynomial_mod, only : nodes_legendre_gauss_lobatto
  use polynomial_mod, only : JacobiP => normalized_jacobi_poly
  !
  integer, parameter :: lp = local_qp
  !
  integer,  intent(in) :: N
  real(lp), dimension(1:(N+1)*(N+2)/2), intent(in) :: r
  !
  real(lp), dimension(1:(N+1)*(N+2)/2) :: return_value
  !
  !.. Local Scalars ..
  integer  :: i,j,np,nn
  real(lp) :: sf
  !
  !.. Local Arrays ..
  real(lp), dimension(1:N+1) :: LGL_nodes
  real(lp), dimension(1:N+1) :: EQD_nodes
  real(lp), dimension(1:(N+1)*(N+2)/2) :: warp
  real(lp), dimension(1:N+1,1:N+1) :: V_EQD
  real(lp), dimension(1:N+1,1:(N+1)*(N+2)/2) :: Pmat
  real(lp), dimension(1:N+1,1:(N+1)*(N+2)/2) :: Lmat
  !
  real(lp), parameter :: zero  = 0.0_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  real(lp), parameter :: eps10 = 1.0e-10_lp
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
    EQD_nodes(i) = two*real(i-1,kind=lp)/real(N,kind=lp) - one
  end do
  !
  ! Compute the transpose of the Vandermonde matrix based on EQD_nodes
  !
  do j = 1,np
    do i = 1,np
      V_EQD(i,j) = JacobiP(j-1,zero,zero,EQD_nodes(i))
    end do
  end do
  V_EQD = transpose(V_EQD)
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
    sf = one - ( deltafun_QP(abs(r(i))<one-eps10) * r(i) )**2
    return_value(i) = warp(i)/sf + &
                      warp(i)*( deltafun_QP(abs(r(i))<one-eps10) - one )
  end do
  !
end function Warpfactor_QP
!
!###############################################################################
!
pure function Nodes2D_DP(N,val) result(return_value)
  !
  integer, parameter :: lp = dp
  !
  integer,  intent(in) :: N
  real(lp), intent(in) :: val
  real(lp), dimension(1:2,1:(N+1)*(N+2)/2) :: return_value
  !
  !.. Local Scalars ..
  real(lp) :: alpha
  !
  !.. Local Arrays ..
  real(lp), dimension(1:(N+1)*(N+2)/2) :: r,warp_r
  real(lp), dimension(1:(N+1)*(N+2)/2) :: s,warp_s
  real(lp), dimension(1:(N+1)*(N+2)/2) :: t,warp_t
  !
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  real(lp), parameter :: three = 3.0_lp
  real(lp), parameter :: four  = 4.0_lp
  real(lp), parameter :: five  = 5.0_lp
  real(lp), parameter :: six   = 6.0_lp
  real(lp), parameter :: pi    = 3.14159265358979323846_lp
  !
  real(lp), parameter :: csq3  = sqrt(three)
  real(lp), parameter :: csq3i = one/csq3
  real(lp), parameter :: cxt   = cos( two*pi/three)
  real(lp), parameter :: cxr   = cos(four*pi/three)
  real(lp), parameter :: cyt   = sin( two*pi/three)
  real(lp), parameter :: cyr   = sin(four*pi/three)
  !
  real(lp), dimension(1:15), parameter :: &
         alpopt = (/ 0.0000_lp, 0.0000_lp, 1.4152_lp, 0.1001_lp, 0.2751_lp, &
                     0.9800_lp, 1.0999_lp, 1.2832_lp, 1.3648_lp, 1.4773_lp, &
                     1.4959_lp, 1.5743_lp, 1.5770_lp, 1.6223_lp, 1.6258_lp  /)
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
  return_value(:,:) = EquiTriNodes2D(N,val)
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
end function Nodes2D_DP
!
!###############################################################################
!
pure function Nodes2D_QP(N,val) result(return_value)
  !
  integer, parameter :: lp = local_qp
  !
  integer,  intent(in) :: N
  real(lp), intent(in) :: val
  real(lp), dimension(1:2,1:(N+1)*(N+2)/2) :: return_value
  !
  !.. Local Scalars ..
  real(lp) :: alpha
  !
  !.. Local Arrays ..
  real(lp), dimension(1:(N+1)*(N+2)/2) :: r,warp_r
  real(lp), dimension(1:(N+1)*(N+2)/2) :: s,warp_s
  real(lp), dimension(1:(N+1)*(N+2)/2) :: t,warp_t
  !
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  real(lp), parameter :: three = 3.0_lp
  real(lp), parameter :: four  = 4.0_lp
  real(lp), parameter :: five  = 5.0_lp
  real(lp), parameter :: six   = 6.0_lp
  real(lp), parameter :: pi    = 3.14159265358979323846_lp
  !
  real(lp), parameter :: csq3  = sqrt(three)
  real(lp), parameter :: csq3i = one/csq3
  real(lp), parameter :: cxt   = cos( two*pi/three)
  real(lp), parameter :: cxr   = cos(four*pi/three)
  real(lp), parameter :: cyt   = sin( two*pi/three)
  real(lp), parameter :: cyr   = sin(four*pi/three)
  !
  real(lp), dimension(1:15), parameter :: &
         alpopt = (/ 0.0000_lp, 0.0000_lp, 1.4152_lp, 0.1001_lp, 0.2751_lp, &
                     0.9800_lp, 1.0999_lp, 1.2832_lp, 1.3648_lp, 1.4773_lp, &
                     1.4959_lp, 1.5743_lp, 1.5770_lp, 1.6223_lp, 1.6258_lp  /)
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
  return_value(:,:) = EquiTriNodes2D(N,val)
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
end function Nodes2D_QP
!
!###############################################################################
!
pure function xytors_DP(xy) result(return_value)
  !
  integer, parameter :: lp = dp
  !
  real(lp), dimension(:,:), intent(in) :: xy
  real(lp), dimension(size(xy,1),size(xy,2)) :: return_value
  !
  !.. Local Scalars ..
  integer   :: i
  real(lp) :: r,s,t
  !
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  real(lp), parameter :: three = 3.0_lp
  real(lp), parameter :: six   = 6.0_lp
  real(lp), parameter :: csq3  = sqrt(three)
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
end function xytors_DP
!
!###############################################################################
!
pure function xytors_QP(xy) result(return_value)
  !
  integer, parameter :: lp = local_qp
  !
  real(lp), dimension(:,:), intent(in) :: xy
  real(lp), dimension(size(xy,1),size(xy,2)) :: return_value
  !
  !.. Local Scalars ..
  integer   :: i
  real(lp) :: r,s,t
  !
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  real(lp), parameter :: three = 3.0_lp
  real(lp), parameter :: six   = 6.0_lp
  real(lp), parameter :: csq3  = sqrt(three)
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
end function xytors_QP
!
!###############################################################################
!
pure function TriNodes2D_AlphaOptimized_DP(N,val) result(return_value)
  !
  integer, parameter :: lp = dp
  !
  integer,  intent(in) :: N
  real(lp), intent(in) :: val
  !
  real(lp), dimension(1:2,1:(N+1)*(N+2)/2) :: return_value
  !
  real(lp), parameter :: half = 0.5_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  return_value(:,:) = Nodes2D(N,val)
  return_value(:,:) = xytors(return_value)
  !
  ! Adjust the triangle nodes to lie in the unit right triangle
  ! defined by the vertices: (0,0) - (1,0) - (0,1)
  !
  return_value(:,:) = half * (return_value(:,:) + one)
  !
end function TriNodes2D_AlphaOptimized_DP
!
!###############################################################################
!
pure function TriNodes2D_AlphaOptimized_QP(N,val) result(return_value)
  !
  integer, parameter :: lp = local_qp
  !
  integer,  intent(in) :: N
  real(lp), intent(in) :: val
  !
  real(lp), dimension(1:2,1:(N+1)*(N+2)/2) :: return_value
  !
  real(lp), parameter :: half = 0.5_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  return_value(:,:) = Nodes2D(N,val)
  return_value(:,:) = xytors(return_value)
  !
  ! Adjust the triangle nodes to lie in the unit right triangle
  ! defined by the vertices: (0,0) - (1,0) - (0,1)
  !
  return_value(:,:) = half * (return_value(:,:) + one)
  !
end function TriNodes2D_AlphaOptimized_QP
!
!###############################################################################
!
pure function TriNodes2D_BarycentricLobatto_DP(N,val) result(return_value)
  !
  use polynomial_mod, only : nodes_legendre_gauss_lobatto
  !
  integer, parameter :: lp = dp
  !
  integer,  intent(in) :: N
  real(lp), intent(in) :: val
  !
  real(lp), dimension(1:2,1:(N+1)*(N+2)/2) :: return_value
  !
  !.. Local Scalars ..
  integer  :: i,j,k,isp,np1d,np2d,m
  real(lp) :: r,s,t,sm
  !
  !.. Local Arrays ..
  real(lp), dimension(1:N+1) :: xi
  real(lp), dimension(1:N+1) :: vp
  !
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  !
continue
  !
  np1d = N+1
  np2d = (N+1)*(N+2)/2
  !
  ! Get the 1D Legendre-Gauss-Lobatto nodes
  !
  call nodes_legendre_gauss_lobatto(N,xi,vp)
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
end function TriNodes2D_BarycentricLobatto_DP
!
!###############################################################################
!
pure function TriNodes2D_BarycentricLobatto_QP(N,val) result(return_value)
  !
  use polynomial_mod, only : nodes_legendre_gauss_lobatto
  !
  integer, parameter :: lp = local_qp
  !
  integer,  intent(in) :: N
  real(lp), intent(in) :: val
  !
  real(lp), dimension(1:2,1:(N+1)*(N+2)/2) :: return_value
  !
  !.. Local Scalars ..
  integer  :: i,j,k,isp,np1d,np2d,m
  real(lp) :: r,s,t,sm
  !
  !.. Local Arrays ..
  real(lp), dimension(1:N+1) :: xi
  real(lp), dimension(1:N+1) :: vp
  !
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  !
continue
  !
  np1d = N+1
  np2d = (N+1)*(N+2)/2
  !
  ! Get the 1D Legendre-Gauss-Lobatto nodes
  !
  call nodes_legendre_gauss_lobatto(N,xi,vp)
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
end function TriNodes2D_BarycentricLobatto_QP
!
!###############################################################################
!
end module triangle_mod
