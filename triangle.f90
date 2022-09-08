module triangle_mod
  !
  ! Module that contains various functions and routines related to
  ! creating the triangular elements and the solution nodes within
  !
  use module_kind_types, only : MP,LK,true,fals
  use module_kind_types, only : Legendre_Gauss
  use module_kind_types, only : Legendre_Gauss_Lobatto
  use module_kind_types, only : AlphaOptimized_Points
  use module_kind_types, only : Barycentric_Points
  use module_kind_types, only : Equi_Distant_Points
  use module_kind_types, only : intseq
#ifdef SPECIAL_FOR_INTEL
  use module_kind_types, only : invert_matrix_lapack95
#else
  use module_kind_types, only : invert_matrix
#endif
  use module_kind_types, only : tria_is_biunit
  !
  implicit none
  !
  !.. Module Precision ..
  integer, parameter :: lp = MP
  !
  private
  !
  public :: TriNodes2D
  public :: get_triangle_weights_at_solpts
  public :: old_get_triangle_weights_at_solpts
  !
  !.. Module Parameter Constants ..
  real(lp), parameter :: zero  = 0.0_lp
  real(lp), parameter :: half  = 0.5_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  real(lp), parameter :: three = 3.0_lp
  real(lp), parameter :: four  = 4.0_lp
  real(lp), parameter :: five  = 5.0_lp
  real(lp), parameter :: six   = 6.0_lp
  !
  real(lp), parameter :: one3  = one/three
  real(lp), parameter :: two3  = two/three
  !
  real(lp), parameter :: eps10 = 1.0e-10_lp
  real(lp), parameter :: eps12 = 1.0e-12_lp
  !
  !
contains
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
function TriNodes2D(p,location,edge_loc) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: p
  integer, intent(in) :: location
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: edge_loc
  !
  !.. Function Result ..
  real(lp) :: return_value(1:2,1:tri_pts(p))
  !
  !.. Local Scalars ..
  logical(lk) :: use_Gauss
  !
continue
  !
  use_Gauss = fals
  if (present(edge_loc)) then
    if (edge_loc == Legendre_Gauss) then
      use_Gauss = true
    end if
  end if
  !
  if (location == AlphaOptimized_Points) then
    !
    return_value(:,:) = TriNodes2D_AlphaOptimized(p,use_Gauss)
    !
  else if (location == Barycentric_Points) then
    !
    return_value(:,:) = TriNodes2D_Barycentric(p,use_Gauss)
    !
  else if (location == Equi_Distant_Points) then
    !
    return_value(:,:) = TriNodes2D_EqualDistant(p)
    !
  else
    !
    return_value(:,:) = TriNodes2D_WilliamsShunn(p)
    !
  end if
  !
end function TriNodes2D
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function TriNodes2D_EqualDistant(p) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: p
  !
  !.. Function Result ..
  real(lp) :: return_value(1:2,1:tri_pts(p))
  !
  !.. Local Scalars ..
  integer :: i,j,n,np
  !
  !.. Local Allocatable Arrays ..
  real(lp), allocatable :: eq_pts(:)
  !
continue
  !
  np = tri_pts(p)
  !
  if (p == 0) then
    return_value(:,:) = -one3
    return
  end if
  !
  ! Create the 1D equal distant points
  !
  ! F2003 AUTO-REALLOC
  eq_pts = real(intseq(-p,p,2),kind=lp) / real(max(p,1),kind=lp)
  !
  ! Create the equidistant nodes in the equilateral triangle
  !
  n = 0
  do j = 1,p+1
    do i = 1,p+2-j
      !
      n = n + 1
      !
      return_value(1,n) = eq_pts(i)
      return_value(2,n) = eq_pts(j)
      !
    end do
  end do
  !
end function TriNodes2D_EqualDistant
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function TriNodes2D_AlphaOptimized(p,use_Gauss) result(return_value)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: p
  logical(lk), intent(in) :: use_Gauss
  !
  !.. Function Result ..
  real(lp) :: return_value(1:2,1:tri_pts(p))
  !
  !.. Local Scalars ..
  integer  :: i,j,n,np
  real(lp) :: a,b,c,xp,yp,alpha
  !
  !.. Local Arrays ..
  real(lp), dimension(1:tri_pts(p)) :: r,warp_r,rmt,x,xeq
  real(lp), dimension(1:tri_pts(p)) :: s,warp_s,smr,y,yeq
  real(lp), dimension(1:tri_pts(p)) :: t,warp_t,tms
  !
  !.. Local Parameters ..
  real(lp), parameter :: csq3  = sqrt(three)
  real(lp), parameter :: csq3i = one/csq3
  !
  real(lp), parameter :: crx = one
  real(lp), parameter :: cry = -csq3i     ! -1/sqrt(3)
  real(lp), parameter :: csy = two*csq3i  !  2/sqrt(3)
  !
  real(lp), parameter :: cxt = -half      ! cos( two*pi/three)
  real(lp), parameter :: cxr = -half      ! cos(four*pi/three)
  real(lp), parameter :: cyt =  half*csq3 ! sin( two*pi/three)
  real(lp), parameter :: cyr = -half*csq3 ! sin(four*pi/three)
  !
  real(lp), parameter :: &
       alpopt(1:15) = [ zero,      zero,      1.4152_lp, 0.1001_lp, 0.2751_lp, &
                        0.9800_lp, 1.0999_lp, 1.2832_lp, 1.3648_lp, 1.4773_lp, &
                        1.4959_lp, 1.5743_lp, 1.5770_lp, 1.6223_lp, 1.6258_lp]
  !
continue
  !
  np = tri_pts(p)
  !
  if (p == 0) then
    return_value(:,:) = -half
    return
  end if
  !
  ! Set optimized parameter, alpha, depending on order P
  !
  alpha = five/three
  if (p <= 15) then
    alpha = alpopt(p)
  end if
  !
  ! Create the equidistant nodes in the equilateral triangle
  !
  n = 0
  do j = 1,p+1
    do i = 1,p+2-j
      !
      n = n + 1
      !
      a = (real(i,kind=lp) - one)/real(p,kind=lp)
      b = (real(j,kind=lp) - one)/real(p,kind=lp)
      c = one - a - b
      !
      xeq(n) = a - c
      yeq(n) = (two*b - a - c)*csq3i
      !
    end do
  end do
  !
  ! Compute the barycentric coordinates for the equilateral triangle
  !
 !r(:) = -half*csq3i*(yeq(:)+csq3i) + half*(xeq(:)+one)
 !s(:) =       csq3i*(yeq(:)+csq3i)
 !t(:) = one - r(:) - s(:)
  !
  r(:) = ( three*xeq(:) - csq3*yeq(:) + two)/six
  s(:) = (                csq3*yeq(:) + one)/three
  t(:) = (-three*xeq(:) - csq3*yeq(:) + two)/six
  !
  ! Get the vectors for each face of the triangle
  !
  tms(:) = t(:) - s(:)
  rmt(:) = r(:) - t(:)
  smr(:) = s(:) - r(:)
  !
  ! Compute the warp shifts due to each face of the equilateral triangle
  !
  warp_r(:) = Warpfactor(p,tms,use_Gauss)
  warp_s(:) = Warpfactor(p,rmt,use_Gauss)
  warp_t(:) = Warpfactor(p,smr,use_Gauss)
  !
  warp_r(:) = (four*s(:)*t(:)) * warp_r(:) * (one + (alpha*r(:))**2)
  warp_s(:) = (four*r(:)*t(:)) * warp_s(:) * (one + (alpha*s(:))**2)
  warp_t(:) = (four*r(:)*s(:)) * warp_t(:) * (one + (alpha*t(:))**2)
  !
  ! Apply the face shifts to the equidistant nodes on the equilateral triangle
  !
  x(:) = xeq(:) + cxr*warp_r(:) + cxt*warp_t(:) + warp_s(:)
  y(:) = yeq(:) + cyr*warp_r(:) + cyt*warp_t(:)
  !
  ! Convert node points from (x,y) in equilateral triangle
  ! to (r,s) coordinates in the standard triangle
  !
  do n = 1,np
    !
    xp = x(n) - half
    yp = y(n) - half*csq3i
    !
    return_value(1,n) = cry*yp + crx*xp
    return_value(2,n) = csy*yp
    !
  end do
  !
  return_value(:,:) = merge(return_value(:,:),zero,abs(return_value(:,:))>eps12)
  !
end function TriNodes2D_AlphaOptimized
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function Warpfactor(p,r,use_Gauss) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : nodes_legendre_gauss
  use polynomial_mod, only : nodes_legendre_gauss_lobatto
  use polynomial_mod, only : JacobiP => normalized_jacobi_poly
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: p
  real(lp),    intent(in) :: r(1:tri_pts(p))
  logical(lk), intent(in) :: use_Gauss
  !
  !.. Function Result ..
  real(lp) :: return_value(1:tri_pts(p))
  !
  !.. Local Scalars ..
  integer  :: i,j,np,nn
  real(lp) :: sf,val
  !
  !.. Local Arrays ..
  real(lp) :: LGL_nodes(1:p+1)
  real(lp) :: EQD_nodes(1:p+1)
  real(lp) :: warp(1:tri_pts(p))
  real(lp) :: V_EQD(1:p+1,1:p+1)
  real(lp) :: Pmat(1:p+1,1:tri_pts(p))
  real(lp) :: Lmat(1:p+1,1:tri_pts(p))
  !
continue
  !
  np = p+1
  nn = tri_pts(p)
  !
  ! Get the 1D Legendre-Gauss-Lobatto nodes
  !
  if (use_Gauss) then
    call nodes_legendre_gauss(p+1,LGL_nodes)
  else
    call nodes_legendre_gauss_lobatto(p,LGL_nodes)
  end if
  !
  ! Get the 1D Equidistant nodes
  !
  do i = 1,np
    EQD_nodes(i) = two*real(i-1,kind=lp)/real(p,kind=lp) - one
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
#ifdef SPECIAL_FOR_INTEL
  Lmat = matmul(invert_matrix_lapack95(V_EQD),Pmat)
#else
  Lmat = matmul(invert_matrix(V_EQD),Pmat)
#endif
  !
  ! Compute warp factor, Eq 6.7 on page 175
  !
  warp = matmul(transpose(Lmat),(LGL_nodes-EQD_nodes))
  !
  ! Scale factor, equation at bottom of page 176
  !
  do i = 1,nn
    val = merge(one, zero, abs(r(i)) < one-eps10)
    sf = one - (val*r(i))**2
    return_value(i) = warp(i)/sf + warp(i)*(val - one)
  end do
  !
end function Warpfactor
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function TriNodes2D_Barycentric(p,use_Gauss) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : nodes_legendre_gauss
  use polynomial_mod, only : nodes_legendre_gauss_lobatto
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: p
  logical(lk), intent(in) :: use_Gauss
  !
  !.. Function Result ..
  real(lp) :: return_value(1:2,1:tri_pts(p))
  !
  !.. Local Scalars ..
  integer  :: i,j,k,n,np1d
  real(lp) :: r,s,t,sm
  !
  !.. Local Arrays ..
  real(lp) :: xi(1:p+1)
  !
continue
  !
  np1d = p+1
  !
  ! Get the 1D Legendre-Gauss-Lobatto nodes
  !
  if (use_Gauss) then
    call nodes_legendre_gauss(p+1,xi)
  else
    call nodes_legendre_gauss_lobatto(p,xi)
  end if
  !
  ! Adjust the 1D LGL nodes to lie on the interval [0,1]
  !
  xi(:) = half*(xi(:) + one)
  !
  ! Get the triangle nodes in r,s space and defined
  ! by the vertices: (0,0) - (1,0) - (0,1)
  !
  n = 0
  do j = 1,np1d
    do i = 1,np1d+1-j
      !
      n = n + 1
      k = np1d + 2 - i - j
      !
      r = xi(i)
      s = xi(j)
      t = xi(k)
      sm = r + s + t
      !
      return_value(1,n) = r/sm
      return_value(2,n) = s/sm
      !
    end do
  end do
  !
  ! Adjust the triangle nodes to lie in the bi-unit right
  ! triangle defined by the vertices: (-1,-1) - (1,-1) - (-1,1)
  !
  if (tria_is_biunit) then
    return_value(:,:) = return_value(:,:)*two - one
  end if
  !
  return_value(:,:) = merge(return_value(:,:),zero,abs(return_value(:,:))>eps12)
  !
end function TriNodes2D_Barycentric
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_triangle_weights_at_solpts(rs) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : legendre_gauss_quadrature
  !
  !.. Formal Arguments ..
  real(lp), intent(in) :: rs(:,:)
  !
  !.. Function Result ..
  real(lp) :: return_value(1:size(rs,2))
  !
  !.. Local Scalars ..
  integer  :: i,j,l,n,nquad
  integer  :: ix,iy
  integer  :: nspts,nfpts
  real(lp) :: x,y,wx,wy,a,b
  !
  !.. Local Arrays ..
  real(lp) :: V(1:size(rs,2),1:size(rs,2))
  real(lp) :: xy(1:2,1:size(rs,2))
  !
  !.. Local Allocatable Arrays ..
  real(lp), allocatable :: xi(:)
  real(lp), allocatable :: wts(:)
  !
continue
  !
  ! On entry, r and s need to be defined on the interval [-1,1]
  !
  nspts = size(rs,2)
  nfpts = ( (nint(sqrt(real(9-8*(1-nspts),kind=lp)))-3)/2 ) + 1
  !
  nquad = 2*nfpts - 1 ! quadtrature order is 2*P
  !
  xi = spread(zero,dim=1,ncopies=nquad) ! F2003 AUTO-REALLOC
  wts = spread(zero,dim=1,ncopies=nquad) ! F2003 AUTO-REALLOC
  !
  ! Compute the coefficients of the Vandermonde matrix for
  ! the Lagrange polynomials using monomial basis functions
  !
  xy(:,:) = rs(1:2,:)
  if (tria_is_biunit) then
    xy(:,:) = half * (rs(1:2,:) + one)
  end if
  !
  do n = 1,nspts
    l = 0
    x = xy(1,n)
    y = xy(2,n)
    do i = 1,nfpts
      do j = 1,i
        l = l + 1
        V(l,n) = x**(i-j) * y**(j-1)
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
  ! Reset the quadrature points and weights to the interval [0,1]
  !
  xi(:) = half*(one+xi(:))
  wts(:) = wts(:)/sum(wts(:))
  !
  ! Initialize the triangle weights to zero
  !
  return_value(:) = zero
  !
  ! Loop over each solution point performing a double integral on the
  ! corresponding Lagrange polynomial over the domain of the triangle
  ! i.e., n    : current solution point
  !       L(n) : Lagrange polynomial for n
  !
  !   /\ 1    /\ 1-xi
  !   |       |
  !   | d(xi) | d(eta) [L(n)]
  !   |       |
  !  \/ 0    \/ 0
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
      ! Inner integral loop
      !
      a = zero
      y_integral: &
      do iy = 1,nquad
        !
        y = (one-x) * xi(iy)
        wy = (one-x) * wts(iy)
        !
        ! Evaluate the Lagrange polynomial at the point x,y
        !
        l = 0
        b = zero
        do i = 1,nfpts
          do j = 1,i
            l = l + 1
            b = b + V(n,l) * x**(i-j) * y**(j-1)
          end do
        end do
        !
        ! Add the weighted polynomial value to the inner integral result
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
  ! Make sure the weights sum to the right amount
  !
  if (tria_is_biunit) then
    ! Make sure the weights sum to 2
    return_value(:) = two * return_value(:) / sum(return_value(:))
  else
    ! Make sure the weights sum to 0.5
    return_value(:) = half * return_value(:) / sum(return_value(:))
  end if
  !
end function get_triangle_weights_at_solpts
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function old_get_triangle_weights_at_solpts(rs) result(return_value)
  !
  ! NOTE : This seems to fail for large n_order (>15)
  !
  !.. Use Statements ..
  use polynomial_mod, only : nodes_legendre_gauss
  use polynomial_mod, only : weights_legendre_gauss
  !
  !.. Formal Arguments ..
  real(lp), intent(in) :: rs(:,:)
  !
  !.. Function Result ..
  real(lp) :: return_value(1:size(rs,2))
  !
  !.. Local Scalars ..
  integer  :: i,j,k,n,m,np
  integer  :: nspts,nfpts
  real(lp) :: g,h,x,y
  !
  !.. Local Parameters ..
  integer, parameter :: nquad = 20
  !
  !.. Local Arrays ..
  real(lp) :: xi(1:nquad)
  real(lp) :: dxi(1:nquad)
  real(lp) :: weights(1:nquad)
  real(lp) :: V(1:size(rs,2),1:size(rs,2))
  !
continue
  !
  ! On entry, r and s need to be defined on the interval [0,1]
  !
  nspts = size(rs,2)
  nfpts = ( (nint(sqrt(real(9-8*(1-nspts),kind=lp)))-3)/2 ) + 1
  !
  ! Compute the coefficients of the Vandermonde matrix for
  ! the Lagrange polynomials using monomial basis functions
  !
  do n = 1,nspts
    k = 0
    do i = 1,nfpts
      do j = 1,i
        k = k + 1
        V(k,n) = rs(1,n)**(i-j) * rs(2,n)**(j-1)
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
  call nodes_legendre_gauss(nquad,xi,dxi)
  call weights_legendre_gauss(nquad,xi,dxi,weights)
  !
  ! Reset the quadrature points and weights to the interval [0,1]
  !
  xi(:) = half*(one+xi(:))
  weights(:) = half*weights(:)
  !
  ! Initialize the triangle weights to zero
  !
  return_value(:) = zero
  !
  ! Loop over each solution point performing a double integral on the
  ! corresponding Lagrange polynomial over the domain of the triangle
  ! i.e., ipt    : current solution point
  !       L(ipt) : Lagrange polynomial for ipt
  !
  !   /\ 1    /\ 1-xi
  !   |       |
  !   | d(xi) | d(eta) [L(ipt)]
  !   |       |
  !  \/ 0    \/ 0
  !
  do np = 1,nspts
    !
    ! Outer integral loop
    !
    outer_integral: do n = 1,nquad
      !
      ! Inner integral loop
      !
      h = zero
      inner_integral: do m = 1,nquad
        !
        x = xi(n)
        y = (one-xi(n))*xi(m)
        !
        ! Evaluate the polynomial at the point x,y
        !
        k = 0
        g = zero
        do i = 1,nfpts
          do j = 1,i
            k = k + 1
            g = g + V(np,k) * x**(i-j) * y**(j-1)
          end do
        end do
        !
        ! Add the weighted polynomial value to the inner integral result
        !
        h = h + (one-xi(n)) * weights(m) * g
        !
      end do inner_integral
      !
      ! Add the weighted inner integral result to the outer integral result
      !
      return_value(np) = return_value(np) + weights(n) * h
      !
    end do outer_integral
    !
  end do
  !
end function old_get_triangle_weights_at_solpts
!
!###############################################################################
!
pure &
function TriNodes2D_WilliamsShunn(p) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: p
  !
  !.. Function Result
  real(lp), allocatable :: return_value(:,:)
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: idx(:)
  !
  !.. Local Parameters ..
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  ! P0 triangle with a single solution point at the center of the triangle.
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  real(lp), parameter :: P0_pts(2,1) = -one3
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  ! P1 triangle solution points using Williams-Shunn triangle points.
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  real(lp), parameter :: P1_pts(2,3) = reshape( [-two3,  one3, &
                                                  one3, -two3, &
                                                 -two3, -two3], shape(P1_pts) )
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  ! P2 triangle solution points using Williams-Shunn triangle points.
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  real(lp), parameter :: P2_pts(2,6) = reshape( [ &
        -0.81684757298044_lp,   0.63369514596088_lp, &
         0.63369514596088_lp,  -0.81684757298044_lp, &
        -0.81684757298044_lp,  -0.81684757298044_lp, &
        -0.108103018168072_lp, -0.783793963663856_lp, &
        -0.783793963663856_lp, -0.108103018168072_lp, &
        -0.108103018168072_lp, -0.108103018168072_lp], &
         shape(P2_pts) )
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  ! P3 triangle solution points using Williams-Shunn triangle points.
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  real(lp), parameter :: P3_pts(2,10) = reshape( [ &
        -one3,                 -one3, &
        -0.888871894660414_lp,  0.777743789320828_lp, &
         0.777743789320828_lp, -0.888871894660414_lp, &
        -0.888871894660414_lp, -0.888871894660414_lp, &
         0.268421495491446_lp, -0.408932576528214_lp, &
        -0.859488918963232_lp, -0.408932576528214_lp, &
        -0.408932576528214_lp, -0.859488918963232_lp, &
         0.268421495491446_lp, -0.859488918963232_lp, &
        -0.859488918963232_lp,  0.268421495491446_lp, &
        -0.408932576528214_lp,  0.268421495491446_lp], &
         shape(P3_pts) )
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  ! P4 triangle solution points using Williams-Shunn triangle points.
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  real(lp), parameter :: P4_pts(2,15) = reshape( [ &
        -0.928258244608532_lp,  0.856516489217064_lp, &
         0.856516489217064_lp, -0.928258244608532_lp, &
        -0.928258244608532_lp, -0.928258244608532_lp, &
        -0.516541208464066_lp,  0.033082416928132_lp, &
         0.033082416928132_lp, -0.516541208464066_lp, &
        -0.516541208464066_lp, -0.516541208464066_lp, &
        -0.051382424445842_lp, -0.897235151108316_lp, &
        -0.897235151108316_lp, -0.051382424445842_lp, &
        -0.051382424445842_lp, -0.051382424445842_lp, &
         0.502367262212968_lp, -0.5969922362364_lp, &
        -0.905375025976568_lp, -0.5969922362364_lp, &
        -0.5969922362364_lp,   -0.905375025976568_lp, &
         0.502367262212968_lp, -0.905375025976568_lp, &
        -0.905375025976568_lp,  0.502367262212968_lp, &
        -0.5969922362364_lp,    0.502367262212968_lp], &
         shape(P4_pts) )
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  ! P5 triangle solution points using Williams-Shunn triangle points.
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  real(lp), parameter :: P5_pts(2,21) = reshape( [ &
        -0.943774095634672_lp,  0.887548191269344_lp, &
         0.887548191269344_lp, -0.943774095634672_lp, &
        -0.943774095634672_lp, -0.943774095634672_lp, &
        -0.645721803061366_lp,  0.291443606122732_lp, &
         0.291443606122732_lp, -0.645721803061366_lp, &
        -0.645721803061366_lp, -0.645721803061366_lp, &
        -0.188982808265134_lp, -0.622034383469732_lp, &
        -0.622034383469732_lp, -0.188982808265134_lp, &
        -0.188982808265134_lp, -0.188982808265134_lp, &
         0.635801960056998_lp, -0.702868375458226_lp, &
        -0.932933584598772_lp, -0.702868375458226_lp, &
        -0.702868375458226_lp, -0.932933584598772_lp, &
         0.635801960056998_lp, -0.932933584598772_lp, &
        -0.932933584598772_lp,  0.635801960056998_lp, &
        -0.702868375458226_lp,  0.635801960056998_lp, &
         0.209957823550266_lp, -0.285607402768638_lp, &
        -0.924350420781628_lp, -0.285607402768638_lp, &
        -0.285607402768638_lp, -0.924350420781628_lp, &
         0.209957823550266_lp, -0.924350420781628_lp, &
        -0.924350420781628_lp,  0.209957823550266_lp, &
        -0.285607402768638_lp,  0.209957823550266_lp], &
         shape(P5_pts) )
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  ! P6 triangle solution points using Williams-Shunn triangle points.
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  real(lp), parameter :: P6_pts(2,28) = reshape( [ &
        -one3,                 -one3, &
        -0.960045625755614_lp,  0.920091251511228_lp, &
         0.920091251511228_lp, -0.960045625755614_lp, &
        -0.960045625755614_lp, -0.960045625755614_lp, &
        -0.736556464940004_lp,  0.473112929880008_lp, &
         0.473112929880008_lp, -0.736556464940004_lp, &
        -0.736556464940004_lp, -0.736556464940004_lp, &
        -0.029729306413078_lp, -0.940541387173844_lp, &
        -0.940541387173844_lp, -0.029729306413078_lp, &
        -0.029729306413078_lp, -0.029729306413078_lp, &
         0.7358224202359_lp,   -0.784096036307978_lp, &
        -0.951726383927922_lp, -0.784096036307978_lp, &
        -0.784096036307978_lp, -0.951726383927922_lp, &
         0.7358224202359_lp,   -0.951726383927922_lp, &
        -0.951726383927922_lp,  0.7358224202359_lp, &
        -0.784096036307978_lp,  0.7358224202359_lp, &
         0.401745140761446_lp, -0.458318454156866_lp, &
        -0.94342668660458_lp,  -0.458318454156866_lp, &
        -0.458318454156866_lp, -0.94342668660458_lp, &
         0.401745140761446_lp, -0.94342668660458_lp, &
        -0.94342668660458_lp,   0.401745140761446_lp, &
        -0.458318454156866_lp,  0.401745140761446_lp, &
         0.073309368412276_lp, -0.366900802310766_lp, &
        -0.70640856610151_lp,  -0.366900802310766_lp, &
        -0.366900802310766_lp, -0.70640856610151_lp, &
         0.073309368412276_lp, -0.70640856610151_lp, &
        -0.70640856610151_lp,   0.073309368412276_lp, &
        -0.366900802310766_lp,  0.073309368412276_lp], &
         shape(P6_pts) )
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  ! P7 triangle solution points using Williams-Shunn triangle points.
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  real(lp), parameter :: P7_pts(2,36) = reshape( [ &
        -0.95765715444107_lp,   0.91531430888214_lp, &
         0.91531430888214_lp,  -0.95765715444107_lp, &
        -0.95765715444107_lp,  -0.95765715444107_lp, &
        -0.798831205208224_lp,  0.597662410416448_lp, &
         0.597662410416448_lp, -0.798831205208224_lp, &
        -0.798831205208224_lp, -0.798831205208224_lp, &
        -0.457923384576136_lp, -0.084153230847728_lp, &
        -0.084153230847728_lp, -0.457923384576136_lp, &
        -0.457923384576136_lp, -0.457923384576136_lp, &
        -0.119617483192336_lp, -0.760765033615328_lp, &
        -0.760765033615328_lp, -0.119617483192336_lp, &
        -0.119617483192336_lp, -0.119617483192336_lp, &
         0.759959282854462_lp, -0.963486641850504_lp, &
        -0.796472641003958_lp, -0.963486641850504_lp, &
        -0.963486641850504_lp, -0.796472641003958_lp, &
         0.759959282854462_lp, -0.796472641003958_lp, &
        -0.796472641003958_lp,  0.759959282854462_lp, &
        -0.963486641850504_lp,  0.759959282854462_lp, &
         0.165124045727344_lp, -0.211933456660026_lp, &
        -0.953190589067318_lp, -0.211933456660026_lp, &
        -0.211933456660026_lp, -0.953190589067318_lp, &
         0.165124045727344_lp, -0.953190589067318_lp, &
        -0.953190589067318_lp,  0.165124045727344_lp, &
        -0.211933456660026_lp,  0.165124045727344_lp, &
         0.503061229085564_lp, -0.547508938181542_lp, &
        -0.955552290904022_lp, -0.547508938181542_lp, &
        -0.547508938181542_lp, -0.955552290904022_lp, &
         0.503061229085564_lp, -0.955552290904022_lp, &
        -0.955552290904022_lp,  0.503061229085564_lp, &
        -0.547508938181542_lp,  0.503061229085564_lp, &
         0.27147436652621_lp,  -0.501841544757336_lp, &
        -0.769632821768874_lp, -0.501841544757336_lp, &
        -0.501841544757336_lp, -0.769632821768874_lp, &
         0.27147436652621_lp,  -0.769632821768874_lp, &
        -0.769632821768874_lp,  0.27147436652621_lp, &
        -0.501841544757336_lp,  0.27147436652621_lp], &
         shape(P7_pts) )
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  ! P8 triangle solution points using Williams-Shunn triangle points.
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  real(lp), parameter :: P8_pts(2,45) = reshape( [ &
        -0.895845867331362_lp,  0.791691734662724_lp, &
         0.791691734662724_lp, -0.895845867331362_lp, &
        -0.895845867331362_lp, -0.895845867331362_lp, &
        -0.772586942778008_lp,  0.545173885556016_lp, &
         0.545173885556016_lp, -0.772586942778008_lp, &
        -0.772586942778008_lp, -0.772586942778008_lp, &
        -0.565662308452164_lp,  0.131324616904328_lp, &
         0.131324616904328_lp, -0.565662308452164_lp, &
        -0.565662308452164_lp, -0.565662308452164_lp, &
        -0.227854451277354_lp, -0.544291097445292_lp, &
        -0.544291097445292_lp, -0.227854451277354_lp, &
        -0.227854451277354_lp, -0.227854451277354_lp, &
        -0.02014813417346_lp,  -0.95970373165308_lp, &
        -0.95970373165308_lp,  -0.02014813417346_lp, &
        -0.02014813417346_lp,  -0.02014813417346_lp, &
         0.669376181231556_lp, -0.703503756777164_lp, &
        -0.965872424454392_lp, -0.703503756777164_lp, &
        -0.703503756777164_lp, -0.965872424454392_lp, &
         0.669376181231556_lp, -0.965872424454392_lp, &
        -0.965872424454392_lp,  0.669376181231556_lp, &
        -0.703503756777164_lp,  0.669376181231556_lp, &
         0.414495985308496_lp, -0.57957674232286_lp, &
        -0.834919242985636_lp, -0.57957674232286_lp, &
        -0.57957674232286_lp,  -0.834919242985636_lp, &
         0.414495985308496_lp, -0.834919242985636_lp, &
        -0.834919242985636_lp,  0.414495985308496_lp, &
        -0.57957674232286_lp,   0.414495985308496_lp, &
         0.925803833285858_lp, -0.92588151521622_lp, &
        -0.999922318069638_lp, -0.92588151521622_lp, &
        -0.92588151521622_lp,  -0.999922318069638_lp, &
         0.925803833285858_lp, -0.999922318069638_lp, &
        -0.999922318069638_lp,  0.925803833285858_lp, &
        -0.92588151521622_lp,   0.925803833285858_lp, &
         0.348012222825796_lp, -0.38207975901364_lp, &
        -0.965932463812156_lp, -0.38207975901364_lp, &
        -0.38207975901364_lp,  -0.965932463812156_lp, &
         0.348012222825796_lp, -0.965932463812156_lp, &
        -0.965932463812156_lp,  0.348012222825796_lp, &
        -0.38207975901364_lp,   0.348012222825796_lp, &
         0.08389791764566_lp,  -0.279351646321358_lp, &
        -0.804546271324302_lp, -0.279351646321358_lp, &
        -0.279351646321358_lp, -0.804546271324302_lp, &
         0.08389791764566_lp,  -0.804546271324302_lp, &
        -0.804546271324302_lp,  0.08389791764566_lp, &
        -0.279351646321358_lp,  0.08389791764566_lp], &
         shape(P8_pts) )
  !
continue
  !
  if (p <= 8) then
    idx = TriNodes2D_Reorder(p) ! F2003 AUTO-REALLOC
  end if
  !
  if (p == 0) then
    return_value = P0_pts(:,idx(:)) ! F2003 AUTO-REALLOC
  else if (p == 1) then
    return_value = P1_pts(:,idx(:)) ! F2003 AUTO-REALLOC
  else if (p == 2) then
    return_value = P2_pts(:,idx(:)) ! F2003 AUTO-REALLOC
  else if (p == 3) then
    return_value = P3_pts(:,idx(:)) ! F2003 AUTO-REALLOC
  else if (p == 4) then
    return_value = P4_pts(:,idx(:)) ! F2003 AUTO-REALLOC
  else if (p == 5) then
    return_value = P5_pts(:,idx(:)) ! F2003 AUTO-REALLOC
  else if (p == 6) then
    return_value = P6_pts(:,idx(:)) ! F2003 AUTO-REALLOC
  else if (p == 7) then
    return_value = P7_pts(:,idx(:)) ! F2003 AUTO-REALLOC
  else if (p == 8) then
    return_value = P8_pts(:,idx(:)) ! F2003 AUTO-REALLOC
  else
    return_value = TriNodes2D_Barycentric(p,fals) ! F2003 AUTO-REALLOC
  end if
  !
end function TriNodes2D_WilliamsShunn
!
!###############################################################################
!
pure &
function TriNodes2D_Reorder(p) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: p
  !
  !.. Function Result
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Parameters ..
  integer, parameter :: P0_idx(1) = [1]
  !
  integer, parameter :: P1_idx(3) = [3,2, &
                                       1]
  !
  integer, parameter :: P2_idx(6) = [3,4,2, &
                                       5,6, &
                                         1]
  !
  integer, parameter :: P3_idx(10) = [ 4, 7, 8, 3, &
                                       6, 1, 5,    &
                                       9,10,       &
                                       2]
  !
  integer, parameter :: P4_idx(15) = [ 3,12, 7,13, 2, &
                                      11, 6, 5,10,    &
                                       8, 4, 9,       &
                                      14,15,          &
                                       1]
  !
  integer, parameter :: P5_idx(21) = [ 3,12,18,19,13, 2, &
                                      11, 6, 7, 5,10,    &
                                      17, 8, 9,16,       &
                                      20, 4,21,          &
                                      14,15,             &
                                       1]
  !
  integer, parameter :: P6_idx(28) = [ 4,13,19, 8,20,14, 3, &
                                      12, 7,25,26, 6,11,    &
                                      18,24, 1,23,17,       &
                                       9,27,28,10,          &
                                      21, 5,22,             &
                                      15,16,                &
                                       2]
  !
  integer, parameter :: P7_idx(36) = [ 3,14,27,21,22,28,13, 2, &
                                      15, 6,33,10,34, 5,16,    &
                                      26,32, 9, 8,31,25,       &
                                      20,11, 7,12,19,          &
                                      23,35,36,24,             &
                                      29, 4,30,                &
                                      18,17,                   &
                                       1]
  !
  integer, parameter :: P8_idx(45) = [ 3,30,18,36,13,37,19,31, 2, &
                                      29, 6,24,42,43,25, 5,28,    &
                                      17,23, 9,10, 8,22,16,       &
                                      35,41,11,12,40,34,          &
                                      14,44, 7,45,15,             &
                                      38,26,27,39,                &
                                      20, 4,21,                   &
                                      32,33,                      &
                                       1]
  !
continue
  !
  if (p == 0) then
    return_value = P0_idx(:) ! F2003 AUTO-REALLOC
  else if (p == 1) then
    return_value = P1_idx(:) ! F2003 AUTO-REALLOC
  else if (p == 2) then
    return_value = P2_idx(:) ! F2003 AUTO-REALLOC
  else if (p == 3) then
    return_value = P3_idx(:) ! F2003 AUTO-REALLOC
  else if (p == 4) then
    return_value = P4_idx(:) ! F2003 AUTO-REALLOC
  else if (p == 5) then
    return_value = P5_idx(:) ! F2003 AUTO-REALLOC
  else if (p == 6) then
    return_value = P6_idx(:) ! F2003 AUTO-REALLOC
  else if (p == 7) then
    return_value = P7_idx(:) ! F2003 AUTO-REALLOC
  else if (p == 8) then
    return_value = P8_idx(:) ! F2003 AUTO-REALLOC
  else
    return_value = [(i, i=1,tri_pts(p))] ! F2003 AUTO-REALLOC
  end if
  !
end function TriNodes2D_Reorder
!
!###############################################################################
!
end module triangle_mod
!
!###############################################################################
!
