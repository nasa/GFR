module tetrahedra_mod
  !
  ! Module that contains various functions and routines related to
  ! creating the tetrahedra elements and the solution nodes within
  !
  use module_kind_types
  !
  implicit none
  !
#ifdef DISABLE_QP
  integer, parameter :: local_qp = SP
#else
  integer, parameter :: local_qp = QP
#endif
  !
  private
  !
  public :: TetNodes3D
  public :: get_tetrahedron_weights_at_solpts
  !
  interface TetNodes3D
    module procedure TetNodes3D_DP, &
                     TetNodes3D_QP
  end interface TetNodes3D
  !
  interface TetNodes3D_AlphaOptimized
    module procedure TetNodes3D_AlphaOptimized_DP, &
                     TetNodes3D_AlphaOptimized_QP
  end interface TetNodes3D_AlphaOptimized
  !
  interface TetNodes3D_Barycentric
    module procedure TetNodes3D_Barycentric_DP, &
                     TetNodes3D_Barycentric_QP
  end interface TetNodes3D_Barycentric
  !
  interface get_tetrahedron_weights_at_solpts
    module procedure get_tetrahedron_weights_at_solpts_DP, &
                     get_tetrahedron_weights_at_solpts_QP
  end interface get_tetrahedron_weights_at_solpts
  !
  interface GetFaceShift3D
    module procedure GetFaceShift3D_DP, &
                     GetFaceShift3D_QP
  end interface GetFaceShift3D
  !
  interface WarpShiftFace3D
    module procedure WarpShiftFace3D_DP, &
                     WarpShiftFace3D_QP
  end interface WarpShiftFace3D
  !
  interface evalWarp
    module procedure evalWarp_DP, &
                     evalWarp_QP
  end interface evalWarp
  !
contains
!
!###############################################################################
!
pure integer function tet_pts(p)
  integer, intent(in) :: p
  tet_pts = (p+1)*(p+2)*(p+3)/6
end function tet_pts
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function TetNodes3D_DP(p,location,val,use_Gauss) result(return_value)
  !
  !.. Local Declared Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: p
  real(lp), intent(in) :: val
  integer,  intent(in) :: location
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: use_Gauss
  !
  !.. Function Result..
  real(lp) :: return_value(1:3,1:tet_pts(p))
  !
  !.. Local Scalars ..
  integer :: bary_loc
  !
  !.. Local Parameters ..
  real(lp), parameter :: half = 0.5_lp
  !
continue
  !
  if (p == 0) then
    return_value(:,:) = -half
    return
  end if
  !
  if (location == AlphaOptimized_TetPoints) then
    !
    return_value(:,:) = TetNodes3D_AlphaOptimized(p,val,use_Gauss)
    !
  else if (location == Barycentric_TetPoints) then
    !
    bary_loc = Legendre_Gauss_Lobatto
    !
    if (present(use_Gauss)) then
      if (use_Gauss) then
        bary_loc = Legendre_Gauss
      end if
    end if
    !
    return_value(:,:) = TetNodes3D_Barycentric(p,val,bary_loc)
    !
  else
    !
    return_value(:,:) = real(TetNodes3D_ShunnHam(p),kind=lp)
    !
  end if
  !
end function TetNodes3D_DP
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function TetNodes3D_QP(p,location,val,use_Gauss) result(return_value)
  !
  !.. Local Declared Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: p
  real(lp), intent(in) :: val
  integer,  intent(in) :: location
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: use_Gauss
  !
  !.. Function Result..
  real(lp) :: return_value(1:3,1:tet_pts(p))
  !
  !.. Local Scalars ..
  integer :: bary_loc
  !
  !.. Local Parameters ..
  real(lp), parameter :: half = 0.5_lp
  !
continue
  !
  if (p == 0) then
    return_value(:,:) = -half
    return
  end if
  !
  if (location == AlphaOptimized_TetPoints) then
    !
    return_value(:,:) = TetNodes3D_AlphaOptimized(p,val,use_Gauss)
    !
  else if (location == Barycentric_TetPoints) then
    !
    bary_loc = Legendre_Gauss_Lobatto
    !
    if (present(use_Gauss)) then
      if (use_Gauss) then
        bary_loc = Legendre_Gauss
      end if
    end if
    !
    return_value(:,:) = TetNodes3D_Barycentric(p,val,bary_loc)
    !
  else
    !
    return_value(:,:) = real(TetNodes3D_ShunnHam(p),kind=lp)
    !
  end if
  !
end function TetNodes3D_QP
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function TetNodes3D_AlphaOptimized_DP(p,val,use_Gauss) result(return_value)
  !
  !.. Local Declared Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: p
  real(lp), intent(in) :: val
  !
  !.. Formal Arguments ..
  logical(lk), optional, intent(in) :: use_Gauss
  !
  !.. Function Result ..
  real(lp) :: return_value(1:3,1:tet_pts(p))
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,np
  real(lp) :: alpha,xp,yp,zp
  !
  !.. Local Arrays ..
  real(lp) :: t1(1:3)
  real(lp) :: t2(1:3)
  !
  real(lp), dimension(1:tet_pts(p)) :: x,xeq,shiftx
  real(lp), dimension(1:tet_pts(p)) :: y,yeq,shifty
  real(lp), dimension(1:tet_pts(p)) :: z,zeq,shiftz
  real(lp), dimension(1:tet_pts(p)) :: l1,l2,l3,l4
  !
  !.. Local Parameters ..
  real(lp), parameter :: zero  = 0.0_lp
  real(lp), parameter :: half  = 0.5_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  real(lp), parameter :: three = 3.0_lp
  real(lp), parameter :: six   = 6.0_lp
  !
  real(lp), parameter :: one6  = one/six
  real(lp), parameter :: one3  = one/three
  real(lp), parameter :: two3  = two/three
  !
  real(lp), parameter :: eps12 = 1.0e-12_lp
  !
  real(lp), parameter :: csq3  = sqrt(three)
  real(lp), parameter :: csq3i = one/csq3
  real(lp), parameter :: csq6i = one/sqrt(six)
  !
  real(lp), parameter :: crx = one
  real(lp), parameter :: cry = -csq3i
  real(lp), parameter :: crz = -csq6i
  real(lp), parameter :: csy = two*csq3i
  real(lp), parameter :: csz = crz
  real(lp), parameter :: ctz = three*csq6i
  !
  real(lp), parameter :: &
       alpopt(1:15) = [ zero,      zero,      zero,      0.1002_lp, 1.1332_lp, &
                        1.5608_lp, 1.3413_lp, 1.2577_lp, 1.1603_lp, 1.0153_lp, &
                        0.6080_lp, 0.4523_lp, 0.8856_lp, 0.8717_lp, 0.9655_lp]
  !
continue
  !
  np = tet_pts(p)
  !
  ! Set optimized parameter, alpha, depending on order P
  !
  if (p == 0) then
    return_value(:,:) = -half
    return
  end if
  !
  alpha = one
  if (p <= 16) alpha = alpopt(p)
  !
  ! Create the equidistant nodes in the reference tetrahedron
  !
  l = 1
  do k = 1,p+1
    do j = 1,p+2-k
      do i = 1,p+3-j-k
        !
        xeq(l) = -one + real(2*(i-1),kind=lp)/real(p,kind=lp)
        yeq(l) = -one + real(2*(j-1),kind=lp)/real(p,kind=lp)
        zeq(l) = -one + real(2*(k-1),kind=lp)/real(p,kind=lp)
        !
        l = l + 1
        !
      end do
    end do
  end do
  !
  ! Compute the barycentric coordinate vectors (lamda1,lamda2,lamda3,lamda4)
  !
  l1(:) =  half*(one + zeq(:))
  l2(:) =  half*(one + yeq(:))
  l3(:) = -half*(one + xeq(:) + yeq(:) + zeq(:))
  l4(:) =  half*(one + xeq(:))
  !
  ! Compute the equidistant nodes on the equilateral tetrahedron
  ! 
  x(:) =   -one*l3(:) +   one*l4(:)!+      zero*l2(:) +        zero*l1(:)
  y(:) = -csq3i*l3(:) - csq3i*l4(:) + two*csq3i*l2(:)!+        zero*l1(:)
  z(:) = -csq6i*l3(:) - csq6i*l4(:) -     csq6i*l2(:) + three*csq6i*l1(:)
  !
  !###########################################################################
  !###########################################################################
  !###########################################################################
  !
  !  Compute the warp shifts due to each face of the equilateral tetrahedron
  !
  !###########################################################################
  !###########################################################################
  !###########################################################################
  !
  shiftx(:) = zero
  shifty(:) = zero
  shiftz(:) = zero
  !
  !##############
  !##  Face 1  ##
  !##############
  !
  ! Get the two orthogonal axis tangents for face 1
  !
  t1(:) = [one,  zero, zero]
  t2(:) = [zero, one,  zero]
  !
  call GetFaceShift3D(p,alpha,t1,t2, &
                      la=l1, lb=l2, lc=l3, ld=l4, &
                      shiftx=shiftx,shifty=shifty,shiftz=shiftz, &
                      use_Gauss=use_Gauss)
  !
  !##############
  !##  Face 2  ##
  !##############
  !
  ! Get the two orthogonal axis tangents for face 2
  !
  t1(:) = [one,  zero, zero          ]
  t2(:) = [zero, one3, two3*sqrt(two)]
  !
  call GetFaceShift3D(p,alpha,t1,t2, &
                      la=l2, lb=l1, lc=l3, ld=l4, &
                      shiftx=shiftx,shifty=shifty,shiftz=shiftz, &
                      use_Gauss=use_Gauss)
  !
  !##############
  !##  Face 3  ##
  !##############
  !
  ! Get the two orthogonal axis tangents for face 3
  !
  t1(:) = [-half,        half*csq3, zero          ]
  t2(:) = [-half*csq3i, -one6,      two3*sqrt(two)]
  !
  call GetFaceShift3D(p,alpha,t1,t2, &
                      la=l3, lb=l1, lc=l4, ld=l2, &
                      shiftx=shiftx,shifty=shifty,shiftz=shiftz, &
                      use_Gauss=use_Gauss)
  !
  !##############
  !##  Face 4  ##
  !##############
  !
  ! Get the two orthogonal axis tangents for face 4
  !
  t1(:) = [half,        half*csq3, zero          ]
  t2(:) = [half*csq3i, -one6,      two3*sqrt(two)]
  !
  call GetFaceShift3D(p,alpha,t1,t2, &
                      la=l4, lb=l1, lc=l3, ld=l2, &
                      shiftx=shiftx,shifty=shifty,shiftz=shiftz, &
                      use_Gauss=use_Gauss)
  !
  !#############################################################
  !#############################################################
  !#############################################################
  !
  ! Apply the face shifts to the equidistant nodes (x,y,z)
  ! on the equilateral tetrahedron
  !
  !#############################################################
  !#############################################################
  !#############################################################
  !
  x(:) = x(:) + shiftx(:)
  y(:) = y(:) + shifty(:)
  z(:) = z(:) + shiftz(:)
  !
  ! Convert node points from (x,y,z) in equilateral tetrahedron
  ! to (r,s,t) coordinates in the standard tetrahedron
  !
  do l = 1,np
    !
    xp = x(l) - one
    yp = y(l) - csq3i
    zp = z(l) - csq6i
    !
    return_value(1,l) = crz*zp + cry*yp + crx*xp
    return_value(2,l) = csz*zp + csy*yp
    return_value(3,l) = ctz*zp
    !
  end do
  !
  return_value(:,:) = merge(return_value(:,:),zero,abs(return_value(:,:))>eps12)
  !
end function TetNodes3D_AlphaOptimized_DP
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function TetNodes3D_AlphaOptimized_QP(p,val,use_Gauss) result(return_value)
  !
  !.. Local Declared Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: p
  real(lp), intent(in) :: val
  !
  !.. Formal Arguments ..
  logical(lk), optional, intent(in) :: use_Gauss
  !
  !.. Function Result ..
  real(lp) :: return_value(1:3,1:tet_pts(p))
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,np
  real(lp) :: alpha,xp,yp,zp
  !
  !.. Local Arrays ..
  real(lp) :: t1(1:3)
  real(lp) :: t2(1:3)
  !
  real(lp), dimension(1:tet_pts(p)) :: x,xeq,shiftx
  real(lp), dimension(1:tet_pts(p)) :: y,yeq,shifty
  real(lp), dimension(1:tet_pts(p)) :: z,zeq,shiftz
  real(lp), dimension(1:tet_pts(p)) :: l1,l2,l3,l4
  !
  !.. Local Parameters ..
  real(lp), parameter :: zero  = 0.0_lp
  real(lp), parameter :: half  = 0.5_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  real(lp), parameter :: three = 3.0_lp
  real(lp), parameter :: six   = 6.0_lp
  !
  real(lp), parameter :: one6  = one/six
  real(lp), parameter :: one3  = one/three
  real(lp), parameter :: two3  = two/three
  !
  real(lp), parameter :: eps12 = 1.0e-12_lp
  !
  real(lp), parameter :: csq3  = sqrt(three)
  real(lp), parameter :: csq3i = one/csq3
  real(lp), parameter :: csq6i = one/sqrt(six)
  !
  real(lp), parameter :: crx = one
  real(lp), parameter :: cry = -csq3i
  real(lp), parameter :: crz = -csq6i
  real(lp), parameter :: csy = two*csq3i
  real(lp), parameter :: csz = crz
  real(lp), parameter :: ctz = three*csq6i
  !
  real(lp), parameter :: &
       alpopt(1:15) = [ zero,      zero,      zero,      0.1002_lp, 1.1332_lp, &
                        1.5608_lp, 1.3413_lp, 1.2577_lp, 1.1603_lp, 1.0153_lp, &
                        0.6080_lp, 0.4523_lp, 0.8856_lp, 0.8717_lp, 0.9655_lp]
  !
continue
  !
  np = tet_pts(p)
  !
  ! Set optimized parameter, alpha, depending on order P
  !
  if (p == 0) then
    return_value(:,:) = -half
    return
  end if
  !
  alpha = one
  if (p <= 16) alpha = alpopt(p)
  !
  ! Create the equidistant nodes in the reference tetrahedron
  !
  l = 1
  do k = 1,p+1
    do j = 1,p+2-k
      do i = 1,p+3-j-k
        !
        xeq(l) = -one + real(2*(i-1),kind=lp)/real(p,kind=lp)
        yeq(l) = -one + real(2*(j-1),kind=lp)/real(p,kind=lp)
        zeq(l) = -one + real(2*(k-1),kind=lp)/real(p,kind=lp)
        !
        l = l + 1
        !
      end do
    end do
  end do
  !
  ! Compute the barycentric coordinate vectors (lamda1,lamda2,lamda3,lamda4)
  !
  l1(:) =  half*(one + zeq(:))
  l2(:) =  half*(one + yeq(:))
  l3(:) = -half*(one + xeq(:) + yeq(:) + zeq(:))
  l4(:) =  half*(one + xeq(:))
  !
  ! Compute the equidistant nodes on the equilateral tetrahedron
  ! 
  x(:) =   -one*l3(:) +   one*l4(:)!+      zero*l2(:) +        zero*l1(:)
  y(:) = -csq3i*l3(:) - csq3i*l4(:) + two*csq3i*l2(:)!+        zero*l1(:)
  z(:) = -csq6i*l3(:) - csq6i*l4(:) -     csq6i*l2(:) + three*csq6i*l1(:)
  !
  !###########################################################################
  !###########################################################################
  !###########################################################################
  !
  !  Compute the warp shifts due to each face of the equilateral tetrahedron
  !
  !###########################################################################
  !###########################################################################
  !###########################################################################
  !
  shiftx(:) = zero
  shifty(:) = zero
  shiftz(:) = zero
  !
  !##############
  !##  Face 1  ##
  !##############
  !
  ! Get the two orthogonal axis tangents for face 1
  !
  t1(:) = [one,  zero, zero]
  t2(:) = [zero, one,  zero]
  !
  call GetFaceShift3D(p,alpha,t1,t2, &
                      la=l1, lb=l2, lc=l3, ld=l4, &
                      shiftx=shiftx,shifty=shifty,shiftz=shiftz, &
                      use_Gauss=use_Gauss)
  !
  !##############
  !##  Face 2  ##
  !##############
  !
  ! Get the two orthogonal axis tangents for face 2
  !
  t1(:) = [one,  zero, zero          ]
  t2(:) = [zero, one3, two3*sqrt(two)]
  !
  call GetFaceShift3D(p,alpha,t1,t2, &
                      la=l2, lb=l1, lc=l3, ld=l4, &
                      shiftx=shiftx,shifty=shifty,shiftz=shiftz, &
                      use_Gauss=use_Gauss)
  !
  !##############
  !##  Face 3  ##
  !##############
  !
  ! Get the two orthogonal axis tangents for face 3
  !
  t1(:) = [-half,        half*csq3, zero          ]
  t2(:) = [-half*csq3i, -one6,      two3*sqrt(two)]
  !
  call GetFaceShift3D(p,alpha,t1,t2, &
                      la=l3, lb=l1, lc=l4, ld=l2, &
                      shiftx=shiftx,shifty=shifty,shiftz=shiftz, &
                      use_Gauss=use_Gauss)
  !
  !##############
  !##  Face 4  ##
  !##############
  !
  ! Get the two orthogonal axis tangents for face 4
  !
  t1(:) = [half,        half*csq3, zero          ]
  t2(:) = [half*csq3i, -one6,      two3*sqrt(two)]
  !
  call GetFaceShift3D(p,alpha,t1,t2, &
                      la=l4, lb=l1, lc=l3, ld=l2, &
                      shiftx=shiftx,shifty=shifty,shiftz=shiftz, &
                      use_Gauss=use_Gauss)
  !
  !#############################################################
  !#############################################################
  !#############################################################
  !
  ! Apply the face shifts to the equidistant nodes (x,y,z)
  ! on the equilateral tetrahedron
  !
  !#############################################################
  !#############################################################
  !#############################################################
  !
  x(:) = x(:) + shiftx(:)
  y(:) = y(:) + shifty(:)
  z(:) = z(:) + shiftz(:)
  !
  ! Convert node points from (x,y,z) in equilateral tetrahedron
  ! to (r,s,t) coordinates in the standard tetrahedron
  !
  do l = 1,np
    !
    xp = x(l) - one
    yp = y(l) - csq3i
    zp = z(l) - csq6i
    !
    return_value(1,l) = crz*zp + cry*yp + crx*xp
    return_value(2,l) = csz*zp + csy*yp
    return_value(3,l) = ctz*zp
    !
  end do
  !
  return_value(:,:) = merge(return_value(:,:),zero,abs(return_value(:,:))>eps12)
  !
end function TetNodes3D_AlphaOptimized_QP
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine GetFaceShift3D_DP(p,alpha,t1,t2,la,lb,lc,ld,shiftx,shifty,shiftz, &
                          use_Gauss)
  !
  !.. Local Declared Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: p
  real(lp),    intent(in) :: alpha
  real(lp),    intent(in) :: t1(:)
  real(lp),    intent(in) :: t2(:)
  real(lp),    intent(in) :: la(:)
  real(lp),    intent(in) :: lb(:)
  real(lp),    intent(in) :: lc(:)
  real(lp),    intent(in) :: ld(:)
  real(lp), intent(inout) :: shiftx(:)
  real(lp), intent(inout) :: shifty(:)
  real(lp), intent(inout) :: shiftz(:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: use_Gauss
  !
  !.. Local Scalars ..
  integer :: i,np
  !
  !.. Local Arrays ..
  real(lp) :: blend(1:size(la))
  real(lp) :: denom(1:size(la))
  real(lp) :: warp1(1:size(la))
  real(lp) :: warp2(1:size(la))
  !
  !.. Local Parameters ..
  real(lp), parameter :: half  = 0.5_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: eps10 = 1.0e-10_lp
  !
continue
  !
  np = size(la)
  !
  ! Compute warp tangential to face
  !
  call WarpShiftFace3D(p,alpha,lb,lc,ld,warp1,warp2,use_Gauss)
  !
  blend(:) = lb(:)*lc(:)*ld(:)
  denom(:) = (lb(:) + half*la(:)) * (lc(:) + half*la(:)) * (ld(:) + half*la(:))
  !
  do i = 1,np
    if (denom(i) > eps10) then
      blend(i) = (one + (alpha*la(i))**2) * blend(i) / denom(i)
    end if
  end do
  !
  ! Compute warp and blend
  !
  shiftx(:) = shiftx(:) + t1(1)*blend(:)*warp1(:) + t2(1)*blend(:)*warp2(:)
  shifty(:) = shifty(:) + t1(2)*blend(:)*warp1(:) + t2(2)*blend(:)*warp2(:)
  shiftz(:) = shiftz(:) + t1(3)*blend(:)*warp1(:) + t2(3)*blend(:)*warp2(:)
  !
  ! Fix face warp
  !
  do i = 1,np
    if ( la(i) < eps10 .and. any([lb(i),lc(i),ld(i)] < eps10) ) then
      shiftx(i) = t1(1)*warp1(i) + t2(1)*warp2(i)
      shifty(i) = t1(2)*warp1(i) + t2(2)*warp2(i)
      shiftz(i) = t1(3)*warp1(i) + t2(3)*warp2(i)
    end if
  end do
  !
end subroutine GetFaceShift3D_DP
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine GetFaceShift3D_QP(p,alpha,t1,t2,la,lb,lc,ld,shiftx,shifty,shiftz, &
                          use_Gauss)
  !
  !.. Local Declared Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: p
  real(lp),    intent(in) :: alpha
  real(lp),    intent(in) :: t1(:)
  real(lp),    intent(in) :: t2(:)
  real(lp),    intent(in) :: la(:)
  real(lp),    intent(in) :: lb(:)
  real(lp),    intent(in) :: lc(:)
  real(lp),    intent(in) :: ld(:)
  real(lp), intent(inout) :: shiftx(:)
  real(lp), intent(inout) :: shifty(:)
  real(lp), intent(inout) :: shiftz(:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: use_Gauss
  !
  !.. Local Scalars ..
  integer :: i,np
  !
  !.. Local Arrays ..
  real(lp) :: blend(1:size(la))
  real(lp) :: denom(1:size(la))
  real(lp) :: warp1(1:size(la))
  real(lp) :: warp2(1:size(la))
  !
  !.. Local Parameters ..
  real(lp), parameter :: half  = 0.5_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: eps10 = 1.0e-10_lp
  !
continue
  !
  np = size(la)
  !
  ! Compute warp tangential to face
  !
  call WarpShiftFace3D(p,alpha,lb,lc,ld,warp1,warp2,use_Gauss)
  !
  blend(:) = lb(:)*lc(:)*ld(:)
  denom(:) = (lb(:) + half*la(:)) * (lc(:) + half*la(:)) * (ld(:) + half*la(:))
  !
  do i = 1,np
    if (denom(i) > eps10) then
      blend(i) = (one + (alpha*la(i))**2) * blend(i) / denom(i)
    end if
  end do
  !
  ! Compute warp and blend
  !
  shiftx(:) = shiftx(:) + t1(1)*blend(:)*warp1(:) + t2(1)*blend(:)*warp2(:)
  shifty(:) = shifty(:) + t1(2)*blend(:)*warp1(:) + t2(2)*blend(:)*warp2(:)
  shiftz(:) = shiftz(:) + t1(3)*blend(:)*warp1(:) + t2(3)*blend(:)*warp2(:)
  !
  ! Fix face warp
  !
  do i = 1,np
    if ( la(i) < eps10 .and. any([lb(i),lc(i),ld(i)] < eps10) ) then
      shiftx(i) = t1(1)*warp1(i) + t2(1)*warp2(i)
      shifty(i) = t1(2)*warp1(i) + t2(2)*warp2(i)
      shiftz(i) = t1(3)*warp1(i) + t2(3)*warp2(i)
    end if
  end do
  !
end subroutine GetFaceShift3D_QP
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine WarpShiftFace3D_DP(p,alpha,l1,l2,l3,warpx,warpy,use_Gauss)
  !
  !.. Use Statements ..
  use polynomial_mod, only : nodes_legendre_gauss_lobatto
  use polynomial_mod, only : nodes_legendre_gauss
  !
  !.. Local Declared Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: p
  real(lp),    intent(in) :: alpha
  real(lp),    intent(in) :: l1(:)
  real(lp),    intent(in) :: l2(:)
  real(lp),    intent(in) :: l3(:)
  real(lp), intent(inout) :: warpx(:)
  real(lp), intent(inout) :: warpy(:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: use_Gauss
  !
  !.. Local Arrays ..
  real(lp) :: LGL_nodes(1:p+1)
  !
  real(lp), dimension(1:size(l1)) :: b1,l3m2,wf1,w1
  real(lp), dimension(1:size(l1)) :: b2,l1m3,wf2,w2
  real(lp), dimension(1:size(l1)) :: b3,l2m1,wf3,w3
  !
  !.. Local Parameters ..
  real(lp), parameter :: half  = 0.5_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: three = 3.0_lp
  real(lp), parameter :: four  = 4.0_lp
  real(lp), parameter :: cxt   = -half             ! cos( two*pi/three)
  real(lp), parameter :: cxr   = -half             ! cos(four*pi/three)
  real(lp), parameter :: cyt   =  half*sqrt(three) ! sin( two*pi/three)
  real(lp), parameter :: cyr   = -half*sqrt(three) ! sin(four*pi/three)
  !
continue
  !
  ! Get the Lobatto nodes
  !
  call nodes_legendre_gauss_lobatto(p,LGL_nodes)
  if (present(use_Gauss)) then
    if (use_Gauss) then
      call nodes_legendre_gauss(p+1,LGL_nodes)
    end if
  end if
  !
  LGL_nodes(:) = -LGL_nodes(:)
  !
  ! Compute blending function at each node for each edge
  !
  b1(:) = l2(:)*l3(:)
  b2(:) = l1(:)*l3(:)
  b3(:) = l1(:)*l2(:)
  !
  ! Compute the amount of warp for each node, for each edge
  !
  l3m2(:) = l3(:)-l2(:)
  l1m3(:) = l1(:)-l3(:)
  l2m1(:) = l2(:)-l1(:)
  !
  wf1(:) = four*evalWarp(p,LGL_nodes,l3m2)
  wf2(:) = four*evalWarp(p,LGL_nodes,l1m3)
  wf3(:) = four*evalWarp(p,LGL_nodes,l2m1)
  !
  ! Combine the blending and warp
  !
  w1(:) = b1(:)*wf1(:)*(one + (alpha*l1(:))**2)
  w2(:) = b2(:)*wf2(:)*(one + (alpha*l2(:))**2)
  w3(:) = b3(:)*wf3(:)*(one + (alpha*l3(:))**2)
  !
  ! Evaluate the shift in the equilateral triangle face
  !
  warpx(:) = cxt*w2(:) + cxr*w3(:) + w1(:)
  warpy(:) = cyt*w2(:) + cyr*w3(:)
  !
end subroutine WarpShiftFace3D_DP
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine WarpShiftFace3D_QP(p,alpha,l1,l2,l3,warpx,warpy,use_Gauss)
  !
  !.. Use Statements ..
  use polynomial_mod, only : nodes_legendre_gauss_lobatto
  use polynomial_mod, only : nodes_legendre_gauss
  !
  !.. Local Declared Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: p
  real(lp),    intent(in) :: alpha
  real(lp),    intent(in) :: l1(:)
  real(lp),    intent(in) :: l2(:)
  real(lp),    intent(in) :: l3(:)
  real(lp), intent(inout) :: warpx(:)
  real(lp), intent(inout) :: warpy(:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: use_Gauss
  !
  !.. Local Arrays ..
  real(lp) :: LGL_nodes(1:p+1)
  !
  real(lp), dimension(1:size(l1)) :: b1,l3m2,wf1,w1
  real(lp), dimension(1:size(l1)) :: b2,l1m3,wf2,w2
  real(lp), dimension(1:size(l1)) :: b3,l2m1,wf3,w3
  !
  !.. Local Parameters ..
  real(lp), parameter :: half  = 0.5_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: three = 3.0_lp
  real(lp), parameter :: four  = 4.0_lp
  real(lp), parameter :: cxt   = -half             ! cos( two*pi/three)
  real(lp), parameter :: cxr   = -half             ! cos(four*pi/three)
  real(lp), parameter :: cyt   =  half*sqrt(three) ! sin( two*pi/three)
  real(lp), parameter :: cyr   = -half*sqrt(three) ! sin(four*pi/three)
  !
continue
  !
  ! Get the Lobatto nodes
  !
  call nodes_legendre_gauss_lobatto(p,LGL_nodes)
  if (present(use_Gauss)) then
    if (use_Gauss) then
      call nodes_legendre_gauss(p+1,LGL_nodes)
    end if
  end if
  !
  LGL_nodes(:) = -LGL_nodes(:)
  !
  ! Compute blending function at each node for each edge
  !
  b1(:) = l2(:)*l3(:)
  b2(:) = l1(:)*l3(:)
  b3(:) = l1(:)*l2(:)
  !
  ! Compute the amount of warp for each node, for each edge
  !
  l3m2(:) = l3(:)-l2(:)
  l1m3(:) = l1(:)-l3(:)
  l2m1(:) = l2(:)-l1(:)
  !
  wf1(:) = four*evalWarp(p,LGL_nodes,l3m2)
  wf2(:) = four*evalWarp(p,LGL_nodes,l1m3)
  wf3(:) = four*evalWarp(p,LGL_nodes,l2m1)
  !
  ! Combine the blending and warp
  !
  w1(:) = b1(:)*wf1(:)*(one + (alpha*l1(:))**2)
  w2(:) = b2(:)*wf2(:)*(one + (alpha*l2(:))**2)
  w3(:) = b3(:)*wf3(:)*(one + (alpha*l3(:))**2)
  !
  ! Evaluate the shift in the equilateral triangle face
  !
  warpx(:) = cxt*w2(:) + cxr*w3(:) + w1(:)
  warpy(:) = cyt*w2(:) + cyr*w3(:)
  !
end subroutine WarpShiftFace3D_QP
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function evalWarp_DP(p,xi,xout) result(return_value)
  !
  !.. Local Declared Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: p
  real(lp), intent(in) :: xi(:)
  real(lp), intent(in) :: xout(:)
  !
  !.. Function Result..
  real(lp) :: return_value(1:size(xout))
  !
  !.. Local Scalars ..
  integer :: i,j,k
  !
  !.. Local Arrays ..
  real(lp) :: xeq(1:p+1)
  real(lp) :: d(1:size(xout))
  !
  !.. Local Parameters ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  return_value(:) = zero
  !
  do i = 1,p+1
    xeq(i) = -one + real(2*(p+1-i),kind=lp)/real(p,kind=lp)
  end do
  !
  do i = 1,p+1
    !
    d(:) = xi(i) - xeq(i)
    !
    do j = 2,p
      if (i /= j) then
        d(:) = d(:) * (xout(:)-xeq(j)) / (xeq(i)-xeq(j))
      end if
    end do
    !
    if (i /= 1) then
      d(:) = -d(:) / (xeq(i)-xeq(1))
    end if
    !
    if (i /= p+1) then
      d(:) = d(:) / (xeq(i)-xeq(p+1))
    end if
    !
    return_value(:) = return_value(:) + d(:)
    !
  end do
  !
end function evalWarp_DP
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function evalWarp_QP(p,xi,xout) result(return_value)
  !
  !.. Local Declared Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: p
  real(lp), intent(in) :: xi(:)
  real(lp), intent(in) :: xout(:)
  !
  !.. Function Result..
  real(lp) :: return_value(1:size(xout))
  !
  !.. Local Scalars ..
  integer :: i,j,k
  !
  !.. Local Arrays ..
  real(lp) :: xeq(1:p+1)
  real(lp) :: d(1:size(xout))
  !
  !.. Local Parameters ..
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  return_value(:) = zero
  !
  do i = 1,p+1
    xeq(i) = -one + real(2*(p+1-i),kind=lp)/real(p,kind=lp)
  end do
  !
  do i = 1,p+1
    !
    d(:) = xi(i) - xeq(i)
    !
    do j = 2,p
      if (i /= j) then
        d(:) = d(:) * (xout(:)-xeq(j)) / (xeq(i)-xeq(j))
      end if
    end do
    !
    if (i /= 1) then
      d(:) = -d(:) / (xeq(i)-xeq(1))
    end if
    !
    if (i /= p+1) then
      d(:) = d(:) / (xeq(i)-xeq(p+1))
    end if
    !
    return_value(:) = return_value(:) + d(:)
    !
  end do
  !
end function evalWarp_QP
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function TetNodes3D_Barycentric_DP(p,val,location) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : nodes_legendre_gauss
  use polynomial_mod, only : nodes_legendre_gauss_lobatto
  !
  !.. Local Declared Precision ..
  integer, parameter :: lp = dp
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: p
  real(lp), intent(in) :: val
  integer,  intent(in) :: location
  !
  !.. Function Result..
  real(lp) :: return_value(1:3,1:tet_pts(p))
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,n,np
  real(lp) :: a,b,c,d,sm
  !
  !.. Local Arrays ..
  real(lp) :: xi(1:p+1)
  !
  !.. Local Parameters ..
  real(lp), parameter :: half = 0.5_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: two  = 2.0_lp
  !
continue
  !
  if (p == 0) then
    return_value(:,:) = -half
    return
  end if
  !
  np = p+1
  !
  ! Get the 1D Legendre-Gauss-Lobatto nodes
  !
  if (location == Legendre_Gauss) then
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
  do k = 1,np
    do j = 1,np+1-k
      do i = 1,np+2-k-j
        !
        n = n + 1
        !
        l = np + 3 - i - j - k
        !
        a = xi(i)
        b = xi(j)
        c = xi(k)
        d = xi(l)
        !
        sm = a + b + c + d
        !
        return_value(1,n) = a/sm
        return_value(2,n) = b/sm
        return_value(3,n) = c/sm
        !
      end do
    end do
  end do
  !
  ! Adjust the triangle nodes to lie in the bi-unit right
  ! tetrahedron defined by the vertices:     (-1,-1,1)
  !                                          /   |   \
  !                                  (-1,-1,-1)--|----(-1,1,-1)
  !                                           \  |   /
  !                                           (1,-1,-1)
  ! (-1,-1,-1) - (1,-1,-1) - (-1,1,-1) - (-1,-1,1)
  !
  return_value(:,:) = return_value(:,:)*two - one
  !
end function TetNodes3D_Barycentric_DP
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function TetNodes3D_Barycentric_QP(p,val,location) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : nodes_legendre_gauss
  use polynomial_mod, only : nodes_legendre_gauss_lobatto
  !
  !.. Local Declared Precision ..
  integer, parameter :: lp = local_qp
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: p
  real(lp), intent(in) :: val
  integer,  intent(in) :: location
  !
  !.. Function Result..
  real(lp) :: return_value(1:3,1:tet_pts(p))
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,n,np
  real(lp) :: a,b,c,d,sm
  !
  !.. Local Arrays ..
  real(lp) :: xi(1:p+1)
  !
  !.. Local Parameters ..
  real(lp), parameter :: half = 0.5_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: two  = 2.0_lp
  !
continue
  !
  if (p == 0) then
    return_value(:,:) = -half
    return
  end if
  !
  np = p+1
  !
  ! Get the 1D Legendre-Gauss-Lobatto nodes
  !
  if (location == Legendre_Gauss) then
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
  do k = 1,np
    do j = 1,np+1-k
      do i = 1,np+2-k-j
        !
        n = n + 1
        !
        l = np + 3 - i - j - k
        !
        a = xi(i)
        b = xi(j)
        c = xi(k)
        d = xi(l)
        !
        sm = a + b + c + d
        !
        return_value(1,n) = a/sm
        return_value(2,n) = b/sm
        return_value(3,n) = c/sm
        !
      end do
    end do
  end do
  !
  ! Adjust the triangle nodes to lie in the bi-unit right
  ! tetrahedron defined by the vertices:     (-1,-1,1)
  !                                          /   |   \
  !                                  (-1,-1,-1)--|----(-1,1,-1)
  !                                           \  |   /
  !                                           (1,-1,-1)
  ! (-1,-1,-1) - (1,-1,-1) - (-1,1,-1) - (-1,-1,1)
  !
  return_value(:,:) = return_value(:,:)*two - one
  !
end function TetNodes3D_Barycentric_QP
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_tetrahedron_weights_at_solpts_DP(rst) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : legendre_gauss_quadrature
  !
  !.. Local Declared Precision ..
  integer, parameter :: lp = dp
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
  !.. Local Parameters ..
  real(lp), parameter :: zero  = 0.0_lp
  real(lp), parameter :: half  = 0.5_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: four3 = 4.0_lp/3.0_lp
  real(lp), parameter :: eps10 = 1.0e-10_lp
  !
continue
  !
  ! On entry, r and s need to be defined on the interval [-1,1]
  !
  nspts = size(rst,2)
  do n = 0,nspts
    if (tet_pts(n) >= nspts) then
      nfpts = n+1
      exit
    end if
  end do
  !
  nquad = 2*nfpts - 1 ! quadtrature order is 2*P
  !
  xi = init_array(nquad,zero) ! F2003 AUTO-REALLOC
  wts = init_array(nquad,zero) ! F2003 AUTO-REALLOC
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
        do k = 1,j
          l = l + 1
          V(l,n) = x**(i-j) * y**(j-k) * z**(k-1)
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
  ! Initialize the tetrahedron weights to zero
  !
  return_value(:) = zero
  !
  ! Loop over each solution point performing a double integral on the
  ! corresponding Lagrange polynomial over the domain of the tetrahedron
  ! i.e., n    : current solution point
  !       L(n) : Lagrange polynomial for n
  !
  !   /\ 1    /\ 1-xi  /\ 1-xi-eta
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
          ! Multipliy zeta and wts by (1-x-y) to scale
          ! zeta from [0-1] to [0,1-xi-eta]
          !
          z = (one-x-y) * xi(iz)
          wz = (one-x-y) * wts(iz)
          !
          ! Evaluate the Lagrange polynomial at the point x,y,z
          !
          l = 0
          c = zero
          do i = 1,nfpts
            do j = 1,i
              do k = 1,j
                l = l + 1
                c = c + V(n,l) * x**(i-j) * y**(j-k) * z**(k-1)
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
  ! Make sure the weights sum to 4/3
  !
  return_value(:) = four3 * return_value(:) / sum(return_value(:))
  !
end function get_tetrahedron_weights_at_solpts_DP
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_tetrahedron_weights_at_solpts_QP(rst) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : legendre_gauss_quadrature
  !
  !.. Local Declared Precision ..
  integer, parameter :: lp = local_qp
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
  !.. Local Parameters ..
  real(lp), parameter :: zero  = 0.0_lp
  real(lp), parameter :: half  = 0.5_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: four3 = 4.0_lp/3.0_lp
  real(lp), parameter :: eps10 = 1.0e-10_lp
  !
continue
  !
  ! On entry, r and s need to be defined on the interval [-1,1]
  !
  nspts = size(rst,2)
  do n = 0,nspts
    if (tet_pts(n) >= nspts) then
      nfpts = n+1
      exit
    end if
  end do
  !
  nquad = 2*nfpts - 1 ! quadtrature order is 2*P
  !
  xi = init_array(nquad,zero) ! F2003 AUTO-REALLOC
  wts = init_array(nquad,zero) ! F2003 AUTO-REALLOC
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
        do k = 1,j
          l = l + 1
          V(l,n) = x**(i-j) * y**(j-k) * z**(k-1)
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
  ! Initialize the tetrahedron weights to zero
  !
  return_value(:) = zero
  !
  ! Loop over each solution point performing a double integral on the
  ! corresponding Lagrange polynomial over the domain of the tetrahedron
  ! i.e., n    : current solution point
  !       L(n) : Lagrange polynomial for n
  !
  !   /\ 1    /\ 1-xi  /\ 1-xi-eta
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
          ! Multipliy zeta and wts by (1-x-y) to scale
          ! zeta from [0-1] to [0,1-xi-eta]
          !
          z = (one-x-y) * xi(iz)
          wz = (one-x-y) * wts(iz)
          !
          ! Evaluate the Lagrange polynomial at the point x,y,z
          !
          l = 0
          c = zero
          do i = 1,nfpts
            do j = 1,i
              do k = 1,j
                l = l + 1
                c = c + V(n,l) * x**(i-j) * y**(j-k) * z**(k-1)
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
  ! Make sure the weights sum to 4/3
  !
  return_value(:) = four3 * return_value(:) / sum(return_value(:))
  !
end function get_tetrahedron_weights_at_solpts_QP
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function TetNodes3D_ShunnHam(p) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: p
  !
  !.. Function Results ..
  real(wp), allocatable :: return_value(:,:)
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
  ! P0 tetrahedron with a single solution point at the center of the tetrahedron
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  real(wp), parameter :: P0_pts(3,1) = -half
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  ! P1 tetrahedron solution points using Shunn-Ham tetrahedron points.
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  real(wp), parameter :: P1_pts(3,4) = reshape( [ &
        -0.723606797749978_wp, -0.723606797749978_wp,  0.170820393249934_wp, &
         0.170820393249934_wp, -0.723606797749978_wp, -0.723606797749978_wp, &
        -0.723606797749978_wp,  0.170820393249934_wp, -0.723606797749978_wp, &
        -0.723606797749978_wp, -0.723606797749978_wp, -0.723606797749978_wp], &
         [3,4] )
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  ! P2 tetrahedron solution points using Shunn-Ham tetrahedron points.
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  real(wp), parameter :: P2_pts(3,10) = reshape( [ &
       -0.8523301965475532_wp, -0.8523301965475532_wp,  0.5569905896426596_wp, &
        0.5569905896426596_wp, -0.8523301965475532_wp, -0.8523301965475532_wp, &
       -0.8523301965475532_wp,  0.5569905896426596_wp, -0.8523301965475532_wp, &
       -0.8523301965475532_wp, -0.8523301965475532_wp, -0.8523301965475532_wp, &
       -0.187511312231898_wp,  -0.812488687768102_wp,  -0.187511312231898_wp, &
       -0.812488687768102_wp,  -0.187511312231898_wp,  -0.187511312231898_wp, &
       -0.187511312231898_wp,  -0.812488687768102_wp,  -0.812488687768102_wp, &
       -0.812488687768102_wp,  -0.187511312231898_wp,  -0.812488687768102_wp, &
       -0.812488687768102_wp,  -0.812488687768102_wp,  -0.187511312231898_wp, &
       -0.187511312231898_wp,  -0.187511312231898_wp,  -0.812488687768102_wp], &
        [3,10] )
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  ! P3 tetrahedron solution points using Shunn-Ham tetrahedron points.
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  real(wp), parameter :: P3_pts(3,20) = reshape( [ &
      -0.9352948105455122_wp, -0.9352948105455122_wp,  0.8058844316365366_wp, &
       0.8058844316365366_wp, -0.9352948105455122_wp, -0.9352948105455122_wp, &
      -0.9352948105455122_wp,  0.8058844316365366_wp, -0.9352948105455122_wp, &
      -0.9352948105455122_wp, -0.9352948105455122_wp, -0.9352948105455122_wp, &
      -0.380461391454276_wp,  -0.380461391454276_wp,  -0.858615825637172_wp, &
      -0.858615825637172_wp,  -0.380461391454276_wp,  -0.380461391454276_wp, &
      -0.380461391454276_wp,  -0.858615825637172_wp,  -0.380461391454276_wp, &
      -0.380461391454276_wp,  -0.380461391454276_wp,  -0.380461391454276_wp, &
       0.2331930661238736_wp, -0.8792791169497158_wp, -0.474634832224442_wp, &
       0.2331930661238736_wp, -0.8792791169497158_wp, -0.8792791169497158_wp, &
      -0.8792791169497158_wp, -0.8792791169497158_wp,  0.2331930661238736_wp, &
      -0.474634832224442_wp,   0.2331930661238736_wp, -0.8792791169497158_wp, &
      -0.8792791169497158_wp, -0.474634832224442_wp,   0.2331930661238736_wp, &
      -0.8792791169497158_wp,  0.2331930661238736_wp, -0.8792791169497158_wp, &
      -0.474634832224442_wp,  -0.8792791169497158_wp,  0.2331930661238736_wp, &
      -0.8792791169497158_wp, -0.474634832224442_wp,  -0.8792791169497158_wp, &
      -0.8792791169497158_wp, -0.8792791169497158_wp, -0.474634832224442_wp, &
      -0.8792791169497158_wp,  0.2331930661238736_wp, -0.474634832224442_wp, &
      -0.474634832224442_wp,  -0.8792791169497158_wp, -0.8792791169497158_wp, &
       0.2331930661238736_wp, -0.474634832224442_wp,  -0.8792791169497158_wp], &
       [3,20] )
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  ! P4 tetrahedron solution points using Shunn-Ham tetrahedron points.
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  real(wp), parameter :: P4_pts(3,35) = reshape( [ &
       -half,                  -half,                  -half, &
       -0.946526448891253_wp,  -0.946526448891253_wp,   0.839579346673759_wp, &
        0.839579346673759_wp,  -0.946526448891253_wp,  -0.946526448891253_wp, &
       -0.946526448891253_wp,   0.839579346673759_wp,  -0.946526448891253_wp, &
       -0.946526448891253_wp,  -0.946526448891253_wp,  -0.946526448891253_wp, &
       -0.9095091999689656_wp, -0.0904908000310344_wp, -0.9095091999689656_wp, &
       -0.0904908000310344_wp, -0.9095091999689656_wp, -0.9095091999689656_wp, &
       -0.9095091999689656_wp, -0.0904908000310344_wp, -0.0904908000310344_wp, &
       -0.0904908000310344_wp, -0.9095091999689656_wp, -0.0904908000310344_wp, &
       -0.0904908000310344_wp, -0.0904908000310344_wp, -0.9095091999689656_wp, &
       -0.9095091999689656_wp, -0.9095091999689656_wp, -0.0904908000310344_wp, &
        0.4955197769636168_wp, -0.9217955187287024_wp, -0.651928739506212_wp, &
        0.4955197769636168_wp, -0.9217955187287024_wp, -0.9217955187287024_wp, &
       -0.9217955187287024_wp, -0.9217955187287024_wp,  0.4955197769636168_wp, &
       -0.651928739506212_wp,   0.4955197769636168_wp, -0.9217955187287024_wp, &
       -0.9217955187287024_wp, -0.651928739506212_wp,   0.4955197769636168_wp, &
       -0.9217955187287024_wp,  0.4955197769636168_wp, -0.9217955187287024_wp, &
       -0.651928739506212_wp,  -0.9217955187287024_wp,  0.4955197769636168_wp, &
       -0.9217955187287024_wp, -0.651928739506212_wp,  -0.9217955187287024_wp, &
       -0.9217955187287024_wp, -0.9217955187287024_wp, -0.651928739506212_wp, &
       -0.9217955187287024_wp,  0.4955197769636168_wp, -0.651928739506212_wp, &
       -0.651928739506212_wp,  -0.9217955187287024_wp, -0.9217955187287024_wp, &
        0.4955197769636168_wp, -0.651928739506212_wp,  -0.9217955187287024_wp, &
        0.006237290029196_wp,  -0.55359792407537_wp,   -0.899041441878456_wp, &
        0.006237290029196_wp,  -0.55359792407537_wp,   -0.55359792407537_wp, &
       -0.55359792407537_wp,   -0.55359792407537_wp,    0.006237290029196_wp, &
       -0.899041441878456_wp,   0.006237290029196_wp,  -0.55359792407537_wp, &
       -0.55359792407537_wp,   -0.899041441878456_wp,   0.006237290029196_wp, &
       -0.55359792407537_wp,    0.006237290029196_wp,  -0.55359792407537_wp, &
       -0.899041441878456_wp,  -0.55359792407537_wp,    0.006237290029196_wp, &
       -0.55359792407537_wp,   -0.899041441878456_wp,  -0.55359792407537_wp, &
       -0.55359792407537_wp,   -0.55359792407537_wp,   -0.899041441878456_wp, &
       -0.55359792407537_wp,    0.006237290029196_wp,  -0.899041441878456_wp, &
       -0.899041441878456_wp,  -0.55359792407537_wp,   -0.55359792407537_wp, &
        0.006237290029196_wp,  -0.899041441878456_wp,  -0.55359792407537_wp], &
        [3,35] )
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  ! P5 tetrahedron solution points using Shunn-Ham tetrahedron points.
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  real(wp), parameter :: P5_pts(3,56) = reshape( [ &
       -0.9700958696938816_wp, -0.9700958696938816_wp,  0.9102876090816448_wp, &
        0.9102876090816448_wp, -0.9700958696938816_wp, -0.9700958696938816_wp, &
       -0.9700958696938816_wp,  0.9102876090816448_wp, -0.9700958696938816_wp, &
       -0.9700958696938816_wp, -0.9700958696938816_wp, -0.9700958696938816_wp, &
       -0.731043330414012_wp,  -0.731043330414012_wp,   0.193129991242036_wp, &
        0.193129991242036_wp,  -0.731043330414012_wp,  -0.731043330414012_wp, &
       -0.731043330414012_wp,   0.193129991242036_wp,  -0.731043330414012_wp, &
       -0.731043330414012_wp,  -0.731043330414012_wp,  -0.731043330414012_wp, &
        0.55995201688308_wp,   -0.931807957607477_wp,  -0.696336101668126_wp, &
        0.55995201688308_wp,   -0.931807957607477_wp,  -0.931807957607477_wp, &
       -0.931807957607477_wp,  -0.931807957607477_wp,   0.55995201688308_wp, &
       -0.696336101668126_wp,   0.55995201688308_wp,   -0.931807957607477_wp, &
       -0.931807957607477_wp,  -0.696336101668126_wp,   0.55995201688308_wp, &
       -0.931807957607477_wp,   0.55995201688308_wp,   -0.931807957607477_wp, &
       -0.696336101668126_wp,  -0.931807957607477_wp,   0.55995201688308_wp, &
       -0.931807957607477_wp,  -0.696336101668126_wp,  -0.931807957607477_wp, &
       -0.931807957607477_wp,  -0.931807957607477_wp,  -0.696336101668126_wp, &
       -0.931807957607477_wp,   0.55995201688308_wp,   -0.696336101668126_wp, &
       -0.696336101668126_wp,  -0.931807957607477_wp,  -0.931807957607477_wp, &
        0.55995201688308_wp,   -0.696336101668126_wp,  -0.931807957607477_wp, &
        0.1053112862120352_wp, -0.9075896991699966_wp, -0.290131887872042_wp, &
        0.1053112862120352_wp, -0.9075896991699966_wp, -0.9075896991699966_wp, &
       -0.9075896991699966_wp, -0.9075896991699966_wp,  0.1053112862120352_wp, &
       -0.290131887872042_wp,   0.1053112862120352_wp, -0.9075896991699966_wp, &
       -0.9075896991699966_wp, -0.290131887872042_wp,   0.1053112862120352_wp, &
       -0.9075896991699966_wp,  0.1053112862120352_wp, -0.9075896991699966_wp, &
       -0.290131887872042_wp,  -0.9075896991699966_wp,  0.1053112862120352_wp, &
       -0.9075896991699966_wp, -0.290131887872042_wp,  -0.9075896991699966_wp, &
       -0.9075896991699966_wp, -0.9075896991699966_wp, -0.290131887872042_wp, &
       -0.9075896991699966_wp,  0.1053112862120352_wp, -0.290131887872042_wp, &
       -0.290131887872042_wp,  -0.9075896991699966_wp, -0.9075896991699966_wp, &
        0.1053112862120352_wp, -0.290131887872042_wp,  -0.9075896991699966_wp, &
        0.076208645776001_wp,  -0.543619077862478_wp,  -0.988970490051045_wp, &
        0.076208645776001_wp,  -0.543619077862478_wp,  -0.543619077862478_wp, &
       -0.543619077862478_wp,  -0.543619077862478_wp,   0.076208645776001_wp, &
       -0.988970490051045_wp,   0.076208645776001_wp,  -0.543619077862478_wp, &
       -0.543619077862478_wp,  -0.988970490051045_wp,   0.076208645776001_wp, &
       -0.543619077862478_wp,   0.076208645776001_wp,  -0.543619077862478_wp, &
       -0.988970490051045_wp,  -0.543619077862478_wp,   0.076208645776001_wp, &
       -0.543619077862478_wp,  -0.988970490051045_wp,  -0.543619077862478_wp, &
       -0.543619077862478_wp,  -0.543619077862478_wp,  -0.988970490051045_wp, &
       -0.543619077862478_wp,   0.076208645776001_wp,  -0.988970490051045_wp, &
       -0.988970490051045_wp,  -0.543619077862478_wp,  -0.543619077862478_wp, &
        0.076208645776001_wp,  -0.988970490051045_wp,  -0.543619077862478_wp, &
       -0.607632480850882_wp,  -0.295389479824012_wp,  -0.801588559501094_wp, &
       -0.607632480850882_wp,  -0.295389479824012_wp,  -0.295389479824012_wp, &
       -0.295389479824012_wp,  -0.295389479824012_wp,  -0.607632480850882_wp, &
       -0.801588559501094_wp,  -0.607632480850882_wp,  -0.295389479824012_wp, &
       -0.295389479824012_wp,  -0.801588559501094_wp,  -0.607632480850882_wp, &
       -0.295389479824012_wp,  -0.607632480850882_wp,  -0.295389479824012_wp, &
       -0.801588559501094_wp,  -0.295389479824012_wp,  -0.607632480850882_wp, &
       -0.295389479824012_wp,  -0.801588559501094_wp,  -0.295389479824012_wp, &
       -0.295389479824012_wp,  -0.295389479824012_wp,  -0.801588559501094_wp, &
       -0.295389479824012_wp,  -0.607632480850882_wp,  -0.801588559501094_wp, &
       -0.801588559501094_wp,  -0.295389479824012_wp,  -0.295389479824012_wp, &
       -0.607632480850882_wp,  -0.801588559501094_wp,  -0.295389479824012_wp], &
        [3,56] )
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  ! P6 tetrahedron solution points using Shunn-Ham tetrahedron points.
  !
  !=============================================================================
  !=============================================================================
  !=============================================================================
  !
  real(wp), parameter :: P6_pts(3,84) = reshape( [ &
        -0.946243051170366_wp, -0.946243051170366_wp,  0.838729153511098_wp, &
         0.838729153511098_wp, -0.946243051170366_wp, -0.946243051170366_wp, &
        -0.946243051170366_wp,  0.838729153511098_wp, -0.946243051170366_wp, &
        -0.946243051170366_wp, -0.946243051170366_wp, -0.946243051170366_wp, &
        -0.62571864839306_wp,  -0.62571864839306_wp,  -0.12284405482082_wp, &
        -0.12284405482082_wp,  -0.62571864839306_wp,  -0.62571864839306_wp, &
        -0.62571864839306_wp,  -0.12284405482082_wp,  -0.62571864839306_wp, &
        -0.62571864839306_wp,  -0.62571864839306_wp,  -0.62571864839306_wp, &
        -0.355777136338286_wp, -0.355777136338286_wp, -0.932668590985142_wp, &
        -0.932668590985142_wp, -0.355777136338286_wp, -0.355777136338286_wp, &
        -0.355777136338286_wp, -0.932668590985142_wp, -0.355777136338286_wp, &
        -0.355777136338286_wp, -0.355777136338286_wp, -0.355777136338286_wp, &
        -0.947151670255874_wp, -0.052848329744126_wp, -0.947151670255874_wp, &
        -0.052848329744126_wp, -0.947151670255874_wp, -0.947151670255874_wp, &
        -0.947151670255874_wp, -0.052848329744126_wp, -0.052848329744126_wp, &
        -0.052848329744126_wp, -0.947151670255874_wp, -0.052848329744126_wp, &
        -0.052848329744126_wp, -0.052848329744126_wp, -0.947151670255874_wp, &
        -0.947151670255874_wp, -0.947151670255874_wp, -0.052848329744126_wp, &
        -0.704090524054712_wp, -0.295909475945288_wp, -0.704090524054712_wp, &
        -0.295909475945288_wp, -0.704090524054712_wp, -0.704090524054712_wp, &
        -0.704090524054712_wp, -0.295909475945288_wp, -0.295909475945288_wp, &
        -0.295909475945288_wp, -0.704090524054712_wp, -0.295909475945288_wp, &
        -0.295909475945288_wp, -0.295909475945288_wp, -0.704090524054712_wp, &
        -0.704090524054712_wp, -0.704090524054712_wp, -0.295909475945288_wp, &
         0.464619819385896_wp, -0.958093115559888_wp, -0.54843358826612_wp, &
         0.464619819385896_wp, -0.958093115559888_wp, -0.958093115559888_wp, &
        -0.958093115559888_wp, -0.958093115559888_wp,  0.464619819385896_wp, &
        -0.54843358826612_wp,   0.464619819385896_wp, -0.958093115559888_wp, &
        -0.958093115559888_wp, -0.54843358826612_wp,   0.464619819385896_wp, &
        -0.958093115559888_wp,  0.464619819385896_wp, -0.958093115559888_wp, &
        -0.54843358826612_wp,  -0.958093115559888_wp,  0.464619819385896_wp, &
        -0.958093115559888_wp, -0.54843358826612_wp,  -0.958093115559888_wp, &
        -0.958093115559888_wp, -0.958093115559888_wp, -0.54843358826612_wp, &
        -0.958093115559888_wp,  0.464619819385896_wp, -0.54843358826612_wp, &
        -0.54843358826612_wp,  -0.958093115559888_wp, -0.958093115559888_wp, &
         0.464619819385896_wp, -0.54843358826612_wp,  -0.958093115559888_wp, &
         0.295115188173952_wp, -0.806020533753068_wp, -0.683074120667816_wp, &
         0.295115188173952_wp, -0.806020533753068_wp, -0.806020533753068_wp, &
        -0.806020533753068_wp, -0.806020533753068_wp,  0.295115188173952_wp, &
        -0.683074120667816_wp,  0.295115188173952_wp, -0.806020533753068_wp, &
        -0.806020533753068_wp, -0.683074120667816_wp,  0.295115188173952_wp, &
        -0.806020533753068_wp,  0.295115188173952_wp, -0.806020533753068_wp, &
        -0.683074120667816_wp, -0.806020533753068_wp,  0.295115188173952_wp, &
        -0.806020533753068_wp, -0.683074120667816_wp, -0.806020533753068_wp, &
        -0.806020533753068_wp, -0.806020533753068_wp, -0.683074120667816_wp, &
        -0.806020533753068_wp,  0.295115188173952_wp, -0.683074120667816_wp, &
        -0.683074120667816_wp, -0.806020533753068_wp, -0.806020533753068_wp, &
         0.295115188173952_wp, -0.683074120667816_wp, -0.806020533753068_wp, &
         0.585878512939236_wp, -0.804783674219116_wp, -0.976311164501004_wp, &
         0.585878512939236_wp, -0.804783674219116_wp, -0.804783674219116_wp, &
        -0.804783674219116_wp, -0.804783674219116_wp,  0.585878512939236_wp, &
        -0.976311164501004_wp,  0.585878512939236_wp, -0.804783674219116_wp, &
        -0.804783674219116_wp, -0.976311164501004_wp,  0.585878512939236_wp, &
        -0.804783674219116_wp,  0.585878512939236_wp, -0.804783674219116_wp, &
        -0.976311164501004_wp, -0.804783674219116_wp,  0.585878512939236_wp, &
        -0.804783674219116_wp, -0.976311164501004_wp, -0.804783674219116_wp, &
        -0.804783674219116_wp, -0.804783674219116_wp, -0.976311164501004_wp, &
        -0.804783674219116_wp,  0.585878512939236_wp, -0.976311164501004_wp, &
        -0.976311164501004_wp, -0.804783674219116_wp, -0.804783674219116_wp, &
         0.585878512939236_wp, -0.976311164501004_wp, -0.804783674219116_wp, &
         0.082368825600474_wp, -0.732883678592864_wp, -0.406997958913752_wp, &
        -0.732883678592864_wp, -0.406997958913752_wp,  0.082368825600474_wp, &
        -0.732883678592864_wp,  0.082368825600474_wp, -0.406997958913752_wp, &
        -0.406997958913752_wp,  0.082368825600474_wp, -0.942487188093858_wp, &
        -0.406997958913752_wp,  0.082368825600474_wp, -0.732883678592864_wp, &
        -0.942487188093858_wp,  0.082368825600474_wp, -0.406997958913752_wp, &
        -0.406997958913752_wp, -0.732883678592864_wp, -0.942487188093858_wp, &
        -0.406997958913752_wp, -0.942487188093858_wp,  0.082368825600474_wp, &
        -0.942487188093858_wp, -0.732883678592864_wp, -0.406997958913752_wp, &
         0.082368825600474_wp, -0.732883678592864_wp, -0.942487188093858_wp, &
        -0.942487188093858_wp,  0.082368825600474_wp, -0.732883678592864_wp, &
        -0.406997958913752_wp, -0.732883678592864_wp,  0.082368825600474_wp, &
        -0.942487188093858_wp, -0.732883678592864_wp,  0.082368825600474_wp, &
         0.082368825600474_wp, -0.942487188093858_wp, -0.406997958913752_wp, &
        -0.406997958913752_wp, -0.942487188093858_wp, -0.732883678592864_wp, &
        -0.732883678592864_wp, -0.942487188093858_wp, -0.406997958913752_wp, &
         0.082368825600474_wp, -0.942487188093858_wp, -0.732883678592864_wp, &
        -0.732883678592864_wp, -0.942487188093858_wp,  0.082368825600474_wp, &
        -0.732883678592864_wp,  0.082368825600474_wp, -0.942487188093858_wp, &
        -0.942487188093858_wp, -0.406997958913752_wp,  0.082368825600474_wp, &
         0.082368825600474_wp, -0.406997958913752_wp, -0.942487188093858_wp, &
         0.082368825600474_wp, -0.406997958913752_wp, -0.732883678592864_wp, &
        -0.732883678592864_wp, -0.406997958913752_wp, -0.942487188093858_wp, &
        -0.942487188093858_wp, -0.406997958913752_wp, -0.732883678592864_wp], &
         [3,84] )
  !
continue
  !
  if (p <= 6) then
    idx = TetNodes3D_Reorder(p) ! F2003 AUTO-REALLOC
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
  else
    ! F2003 AUTO-REALLOC
    return_value = TetNodes3D_Barycentric(p,zero,Legendre_Gauss_Lobatto)
  end if
  !
end function TetNodes3D_ShunnHam
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function TetNodes3D_Reorder(p) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: p
  !
  !.. Function Result ..
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Parameters ..
  integer, parameter :: P0_idx(1) = [1]
  !
  integer, parameter :: P1_idx(4) = [4,2, &
                                       3, &
                                       1]
  !
  integer, parameter :: P2_idx(10) = [ 4, 7, 2, &
                                          8,10, &
                                             3, &
                                          9, 5, &
                                             6, &
                                             1]
  !
  integer, parameter :: P3_idx(20) = [ 4,19,10, 2, &
                                         16, 5,20, &
                                            14,12, &
                                                3, &
                                         17, 7, 9, &
                                             6, 8, &
                                               18, &
                                            11,15, &
                                               13, &
                                                1]
  !
  integer, parameter :: P4_idx(35) = [ 5,22, 7,13, 3, &
                                         19,32,24,23, &
                                             6,33,10, &
                                               17,15, &
                                                   4, &
                                         20,31,35,12, &
                                            34, 1,25, &
                                               27,29, &
                                                  21, &
                                            11,28, 9, &
                                               30,26, &
                                                   8, &
                                               14,18, &
                                                  16, &
                                                   2]
  !
  integer, parameter :: P5_idx(56) = [ 4,19,31,22,10, 2, &
                                         16,41,54,33,20, &
                                            28,45,53,32, &
                                               26,42,24, &
                                                  14,12, &
                                                      3, &
                                         17,40,49,44, 9, &
                                            43, 8, 6,34, &
                                               51, 7,47, &
                                                  36,38, &
                                                     18, &
                                            29,56,52,21, &
                                               48, 5,50, &
                                                  55,46, &
                                                     30, &
                                               23,37,27, &
                                                  39,35, &
                                                     25, &
                                                  11,15, &
                                                     13, &
                                                      1]
  !
  integer, parameter :: P6_idx(84) = [ 4,47,35,14,26,38, 2, &
                                         44,57,67,70,49,48, &
                                            32,83, 9,81,36, &
                                               13,79,64,17, &
                                                  30,58,28, &
                                                     42,40, &
                                                         3, &
                                         45,56,75,77,60,37, &
                                            59, 8,20, 6,50, &
                                               84,19,23,82, &
                                                  71, 7,65, &
                                                     52,54, &
                                                        46, &
                                            33,76,11,74,25, &
                                               69,24,22,61, &
                                                  10,21,12, &
                                                     66,63, &
                                                        34, &
                                               18,78,68,16, &
                                                  73, 5,72, &
                                                     80,62, &
                                                        15, &
                                                  27,53,31, &
                                                     55,51, &
                                                        29, &
                                                     39,43, &
                                                        41, &
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
  else
    return_value = [(i, i=1,tet_pts(p))] ! F2003 AUTO-REALLOC
  end if
  !
end function TetNodes3D_Reorder
!
!###############################################################################
!
end module tetrahedra_mod
!
!###############################################################################
!
