module bc_slip_wall_mod
  !
  !.. Use Statements
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nee,nq
  use eqn_idx, only : ncs
  !
  implicit none
  !
  private
  !
  public :: get_bc_slip_wall_solution
  public :: get_bc_slip_wall_gradient
  !
contains
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine get_bc_slip_wall_solution(n1,n2,facenrm,face_cv)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Slip wall / symmetry boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conserved variables at the flux points
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  if (nr == 3) then
    call get_bc_slip_wall_solution_3d(n1,n2,facenrm,face_cv)
  else if (nr == 2) then
    call get_bc_slip_wall_solution_2d(n1,n2,facenrm,face_cv)
  else
    call get_bc_slip_wall_solution_1d(n1,n2,facenrm,face_cv)
  end if
  !
end subroutine get_bc_slip_wall_solution
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine get_bc_slip_wall_solution_3d(n1,n2,facenrm,face_cv)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Slip wall / symmetry boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conserved variables at the flux points
  !
  !.. Use Statements ..
  use ovar, only : walls_are_exact
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: rho,rho_inv,vx,vy,vz,vn,rhoe
  real(wp) :: velocity_const,energy_const
  !
continue
  !
  velocity_const = merge(one,two,walls_are_exact)
  energy_const = merge(-one,zero,walls_are_exact)
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    rnz = facenrm(3,n)
    ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
    unx = rnx/ds
    uny = rny/ds
    unz = rnz/ds
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    vz = face_cv(n,nmz,1)*rho_inv
    rhoe = face_cv(n,nee,1)
    !
    vn = vx*unx + vy*uny + vz*unz
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rho*(vx - velocity_const*vn*unx)
    face_cv(n,nmy,2) = rho*(vy - velocity_const*vn*uny)
    face_cv(n,nmz,2) = rho*(vz - velocity_const*vn*unz)
    face_cv(n,nee,2) = rhoe + energy_const*half*rho*vn*vn
    !
  end do
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv)
    do m = 1,ncs
      do n = n1,n2
        face_cv(n,m,2) = face_cv(n,m,1)
      end do
    end do
  end if
  !
end subroutine get_bc_slip_wall_solution_3d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine get_bc_slip_wall_solution_2d(n1,n2,facenrm,face_cv)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Slip wall / symmetry boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conserved variables at the flux points
  !
  !.. Use Statements ..
  use ovar, only : walls_are_exact
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: ds,rnx,rny,unx,uny
  real(wp) :: rho,rho_inv,vx,vy,vn,rhoe
  real(wp) :: velocity_const,energy_const
  !
continue
  !
  velocity_const = merge(one,two,walls_are_exact)
  energy_const = merge(-one,zero,walls_are_exact)
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds = sqrt( rnx*rnx + rny*rny )
    unx = rnx/ds
    uny = rny/ds
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    rhoe = face_cv(n,nee,1)
    !
    vn = vx*unx + vy*uny
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rho*(vx - velocity_const*vn*unx)
    face_cv(n,nmy,2) = rho*(vy - velocity_const*vn*uny)
    face_cv(n,nee,2) = rhoe + energy_const*half*rho*vn*vn
    !
  end do
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv)
    do m = 1,ncs
      do n = n1,n2
        face_cv(n,m,2) = face_cv(n,m,1)
      end do
    end do
  end if
  !
end subroutine get_bc_slip_wall_solution_2d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine get_bc_slip_wall_solution_1d(n1,n2,facenrm,face_cv)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Slip wall / symmetry boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conserved variables at the flux points
  !
  !.. Use Statements ..
  use ovar, only : walls_are_exact
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: rho,vx,rhoe
  real(wp) :: velocity_const,energy_const
  !
continue
  !
  velocity_const = merge(one,two,walls_are_exact)
  energy_const = merge(-one,zero,walls_are_exact)
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv)
  do n = n1,n2
    !
    rho = face_cv(n,nec,1)
    vx = face_cv(n,nmx,1)/rho
    rhoe = face_cv(n,nee,1)
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rho*(one - velocity_const)*vx
    face_cv(n,nee,2) = rhoe + energy_const*half*rho*vx*vx
    !
  end do
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv)
    do m = 1,ncs
      do n = n1,n2
        face_cv(n,m,2) = face_cv(n,m,1)
      end do
    end do
  end if
  !
end subroutine get_bc_slip_wall_solution_1d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine get_bc_slip_wall_gradient(n1,n2,facenrm,face_cv,face_gp,face_dcvdx)
  !
  ! ###########################################################
  ! ###########################################################
  ! ##  Slip wall / symmetry boundary condition for gradients
  ! ###########################################################
  ! ###########################################################
  !
  ! N1         : Beginning face point for this boundary condition
  ! N2         : Ending    face point for this boundary condition
  ! FACENRM    : Face outward normal vector at the flux points
  ! FACE_CV    : Conserved    variables at the flux points
  ! FACE_GP    : Gas property variables at the flux points
  ! FACE_DCVDX : Gradients of the conserved variables at the flux points
  !
  !.. Use Statements ..
  use gpu_mod, only : grad_cv_to_grad_pv_gpu
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:)
  real(wp),    intent(in) :: face_gp(:,:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr,np
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: dpvdx(:,:,:)
  !
continue
  !
  np = size(face_dcvdx,dim=1)
  nr = size(face_dcvdx,dim=2)
  !
  allocate ( dpvdx(1:np,1:nr,1:nq) )
  !$ACC DATA CREATE(dpvdx)
  !
  ! Convert the interior host cell conservative gradients to primitive
  !
  call grad_cv_to_grad_pv_gpu(n1,n2,face_cv,face_gp,face_dcvdx(:,:,:,1),dpvdx)
  !
  ! Get the primitive variable gradients at the ghost cell flux points
  !
  if (nr == 3) then
    call get_bc_slip_wall_gradient_3d(n1,n2,facenrm,dpvdx,face_dcvdx(:,:,:,2))
  else if (nr == 2) then
    call get_bc_slip_wall_gradient_2d(n1,n2,facenrm,dpvdx,face_dcvdx(:,:,:,2))
  else
    call get_bc_slip_wall_gradient_1d(n1,n2,facenrm,dpvdx,face_dcvdx(:,:,:,2))
  end if
  !
  !$ACC END DATA
  deallocate ( dpvdx )
  !
end subroutine get_bc_slip_wall_gradient
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine get_bc_slip_wall_gradient_3d(n1,n2,facenrm,dpvdx,face_dcvdx)
  !
  ! ###########################################################
  ! ###########################################################
  ! ##  Slip wall / symmetry boundary condition for gradients
  ! ###########################################################
  ! ###########################################################
  !
  ! N1         : Beginning face point for this boundary condition
  ! N2         : Ending    face point for this boundary condition
  ! FACENRM    : Face outward normal vector at the flux points
  ! FACE_CV    : Conserved variables at the flux points
  ! FACE_DCVDX : Gradients of the conserved variables at the flux points
  !
  !.. Use Statements ..
  use ovar, only : walls_are_exact
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: dpvdx(:,:,:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: drdn,dudn,dvdn,dwdn,dtdn,dsdn
  real(wp) :: reflect_const
  !
continue
  !
  reflect_const = merge(one,two,walls_are_exact)
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_dcvdx,dpvdx)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    rnz = facenrm(3,n)
    ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
    unx = rnx/ds
    uny = rny/ds
    unz = rnz/ds
    !
    drdn = dpvdx(n,1,nec)*unx + dpvdx(n,2,nec)*uny + dpvdx(n,3,nec)*unz
    dudn = dpvdx(n,1,nmx)*unx + dpvdx(n,2,nmx)*uny + dpvdx(n,3,nmx)*unz
    dvdn = dpvdx(n,1,nmy)*unx + dpvdx(n,2,nmy)*uny + dpvdx(n,3,nmy)*unz
    dwdn = dpvdx(n,1,nmz)*unx + dpvdx(n,2,nmz)*uny + dpvdx(n,3,nmz)*unz
    dtdn = dpvdx(n,1,nee)*unx + dpvdx(n,2,nee)*uny + dpvdx(n,3,nee)*unz
    !
    face_dcvdx(n,1,nec) = dpvdx(n,1,nec) - reflect_const*unx*drdn
    face_dcvdx(n,2,nec) = dpvdx(n,2,nec) - reflect_const*uny*drdn
    face_dcvdx(n,3,nec) = dpvdx(n,3,nec) - reflect_const*unz*drdn
    !
    face_dcvdx(n,1,nmx) = dpvdx(n,1,nmx) - reflect_const*unx*dudn
    face_dcvdx(n,2,nmx) = dpvdx(n,2,nmx) - reflect_const*uny*dudn
    face_dcvdx(n,3,nmx) = dpvdx(n,3,nmx) - reflect_const*unz*dudn
    !
    face_dcvdx(n,1,nmy) = dpvdx(n,1,nmy) - reflect_const*unx*dvdn
    face_dcvdx(n,2,nmy) = dpvdx(n,2,nmy) - reflect_const*uny*dvdn
    face_dcvdx(n,3,nmy) = dpvdx(n,3,nmy) - reflect_const*unz*dvdn
    !
    face_dcvdx(n,1,nmz) = dpvdx(n,1,nmz) - reflect_const*unx*dwdn
    face_dcvdx(n,2,nmz) = dpvdx(n,2,nmz) - reflect_const*uny*dwdn
    face_dcvdx(n,3,nmz) = dpvdx(n,3,nmz) - reflect_const*unz*dwdn
    !
    face_dcvdx(n,1,nee) = dpvdx(n,1,nee) - reflect_const*unx*dtdn
    face_dcvdx(n,2,nee) = dpvdx(n,2,nee) - reflect_const*uny*dtdn
    face_dcvdx(n,3,nee) = dpvdx(n,3,nee) - reflect_const*unz*dtdn
    !
  end do
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(facenrm,face_dcvdx,dpvdx)
    do m = 1,ncs
      do n = n1,n2
        !
        rnx = facenrm(1,n)
        rny = facenrm(2,n)
        rnz = facenrm(3,n)
        ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
        unx = rnx/ds
        uny = rny/ds
        unz = rnz/ds
        !
        dsdn = dpvdx(n,1,m)*unx + dpvdx(n,2,m)*uny + dpvdx(n,3,m)*unz
        !
        face_dcvdx(n,1,m) = dpvdx(n,1,m) - reflect_const*unx*dsdn
        face_dcvdx(n,2,m) = dpvdx(n,2,m) - reflect_const*uny*dsdn
        face_dcvdx(n,3,m) = dpvdx(n,3,m) - reflect_const*unz*dsdn
        !
      end do
    end do
  end if
  !
end subroutine get_bc_slip_wall_gradient_3d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine get_bc_slip_wall_gradient_2d(n1,n2,facenrm,dpvdx,face_dcvdx)
  !
  ! ###########################################################
  ! ###########################################################
  ! ##  Slip wall / symmetry boundary condition for gradients
  ! ###########################################################
  ! ###########################################################
  !
  ! N1         : Beginning face point for this boundary condition
  ! N2         : Ending    face point for this boundary condition
  ! FACENRM    : Face outward normal vector at the flux points
  ! FACE_CV    : Conserved variables at the flux points
  ! FACE_DCVDX : Gradients of the conserved variables at the flux points
  !
  !.. Use Statements ..
  use ovar, only : walls_are_exact
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: dpvdx(:,:,:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: ds,rnx,rny,unx,uny
  real(wp) :: drdn,dudn,dvdn,dtdn,dsdn
  real(wp) :: reflect_const
  !
continue
  !
  reflect_const = merge(one,two,walls_are_exact)
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_dcvdx,dpvdx)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds = sqrt( rnx*rnx + rny*rny )
    unx = rnx/ds
    uny = rny/ds
    !
    drdn = dpvdx(n,1,nec)*unx + dpvdx(n,2,nec)*uny
    dudn = dpvdx(n,1,nmx)*unx + dpvdx(n,2,nmx)*uny
    dvdn = dpvdx(n,1,nmy)*unx + dpvdx(n,2,nmy)*uny
    dtdn = dpvdx(n,1,nee)*unx + dpvdx(n,2,nee)*uny
    !
    face_dcvdx(n,1,nec) = dpvdx(n,1,nec) - reflect_const*unx*drdn
    face_dcvdx(n,2,nec) = dpvdx(n,2,nec) - reflect_const*uny*drdn
    !
    face_dcvdx(n,1,nmx) = dpvdx(n,1,nmx) - reflect_const*unx*dudn
    face_dcvdx(n,2,nmx) = dpvdx(n,2,nmx) - reflect_const*uny*dudn
    !
    face_dcvdx(n,1,nmy) = dpvdx(n,1,nmy) - reflect_const*unx*dvdn
    face_dcvdx(n,2,nmy) = dpvdx(n,2,nmy) - reflect_const*uny*dvdn
    !
    face_dcvdx(n,1,nee) = dpvdx(n,1,nee) - reflect_const*unx*dtdn
    face_dcvdx(n,2,nee) = dpvdx(n,2,nee) - reflect_const*uny*dtdn
    !
  end do
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(facenrm,face_dcvdx,dpvdx)
    do m = 1,ncs
      do n = n1,n2
        !
        rnx = facenrm(1,n)
        rny = facenrm(2,n)
        ds = sqrt( rnx*rnx + rny*rny )
        unx = rnx/ds
        uny = rny/ds
        !
        dsdn = dpvdx(n,1,m)*unx + dpvdx(n,2,m)*uny
        !
        face_dcvdx(n,1,m) = dpvdx(n,1,m) - reflect_const*unx*dsdn
        face_dcvdx(n,2,m) = dpvdx(n,2,m) - reflect_const*uny*dsdn
        !
      end do
    end do
  end if
  !
end subroutine get_bc_slip_wall_gradient_2d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine get_bc_slip_wall_gradient_1d(n1,n2,facenrm,dpvdx,face_dcvdx)
  !
  ! ###########################################################
  ! ###########################################################
  ! ##  Slip wall / symmetry boundary condition for gradients
  ! ###########################################################
  ! ###########################################################
  !
  ! N1         : Beginning face point for this boundary condition
  ! N2         : Ending    face point for this boundary condition
  ! FACENRM    : Face outward normal vector at the flux points
  ! FACE_CV    : Conserved variables at the flux points
  ! FACE_DCVDX : Gradients of the conserved variables at the flux points
  !
  !.. Use Statements ..
  use ovar, only : walls_are_exact
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: dpvdx(:,:,:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: reflect_const
  !
continue
  !
 !reflect_const = merge(one,two,walls_are_exact)
  reflect_const = merge(zero,-one,walls_are_exact)
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_dcvdx,dpvdx)
  do m = 1,nq
    do n = n1,n2
     !face_dcvdx(n,1,m) = (one - reflect_const)*dpvdx(n,1,m)
      face_dcvdx(n,1,m) = reflect_const*dpvdx(n,1,m)
    end do
  end do
  !
end subroutine get_bc_slip_wall_gradient_1d
!
!###############################################################################
!
end module bc_slip_wall_mod
