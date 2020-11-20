module bc_mms_mod
  !
  !.. Use Statements
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nee
  !
  implicit none
  !
  private
  !
  public :: get_bc_mms_dirichlet_solution
  !
contains
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine get_bc_mms_dirichlet_solution(n1,n2,facenrm,facexyz,face_cv)
  !
  ! #############################################################
  ! #############################################################
  ! ##  MMS Dirichlet (generic freeflow) Boundary Condition
  ! #############################################################
  ! #############################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conserved variables at the flux points
  ! FACEXYZ : XYZ coordinates for all the flux points
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: facexyz(:,:)
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
    call get_bc_mms_dirichlet_solution_3d(n1,n2,facenrm,facexyz,face_cv)
  else if (nr == 2) then
    call get_bc_mms_dirichlet_solution_2d(n1,n2,facenrm,facexyz,face_cv)
  else
    call get_bc_mms_dirichlet_solution_1d(n1,n2,facenrm,facexyz,face_cv)
  end if
  !
end subroutine get_bc_mms_dirichlet_solution
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine get_bc_mms_dirichlet_solution_3d(n1,n2,facenrm,facexyz,face_cv)
  !
  ! #########################################################
  ! #########################################################
  ! ##  MMS Dirichlet (generic freeflow) Boundary Condition
  ! ##    - At each flux point, use the normal velocity
  ! ##      Mach number to determine if it is an inflow
  ! ##      or outflow boundary condition and whether
  ! ##      it is subsonic or supersonic.
  ! #########################################################
  ! #########################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conserved variables at the flux points
  ! FACEXYZ : XYZ coordinates for all the flux points
  !
  !.. Use Statements ..
  use ovar,    only : gam
  use mms_mod, only : mms_solution_at_xyz_gpu
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: facexyz(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,gm1_inv
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: rho,rhou,rhov,rhow,rhoe,rhoq2,ske,pres,mach
  real(wp) :: rho_mms,rhou_mms,rhov_mms,rhow_mms,rhoe_mms
  real(wp) :: rhoq2_mms,ske_mms,pres_mms
  !
continue
  !
  gm1 = gam - one
  gm1_inv = one/gm1
  !
  call mms_solution_at_xyz_gpu(n1,n2,facexyz,face_cv(:,:,2), &
                               return_nondim=true, &
                               return_cv=true)
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
  do n = n1,n2
    !
    ! Unit normal at this flux point
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    rnz = facenrm(3,n)
    ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
    unx = rnx/ds
    uny = rny/ds
    unz = rnz/ds
    !
    rho_mms = face_cv(n,nec,2)
    rhou_mms = face_cv(n,nmx,2)
    rhov_mms = face_cv(n,nmy,2)
    rhow_mms = face_cv(n,nmz,2)
    rhoe_mms = face_cv(n,nee,2)
    rhoq2_mms = rhou_mms*rhou_mms + rhov_mms*rhov_mms + rhow_mms*rhow_mms
    ske_mms = half*rhoq2_mms/rho_mms
    pres_mms = gm1*(rhoe_mms - ske_mms)
    !
    ! Compute the Mach number of the normal velocity at this flux point
    !
    rho = face_cv(n,nec,1)
    rhou = face_cv(n,nmx,1)
    rhov = face_cv(n,nmy,1)
    rhow = face_cv(n,nmz,1)
    rhoe = face_cv(n,nee,1)
    rhoq2 = rhou*rhou + rhov*rhov + rhow*rhow
    ske = half*rhoq2/rho
    pres = gm1*(rhoe - ske)
    !
    ! NOTE: The speed of sound (denominator) is gam*pres*rho instead of
    !       gam*pres/rho because the normal velocity (numerator) is rho*vn
    !       instead of just vn.
    !
    mach = (rhou*unx + rhov*uny + rhow*unz) / sqrt( max(gam*pres*rho,rndoff) )
    !
    ! If mach >= one, do nothing because in that case we want to extrapolate
    ! everything from the interior and we already have those values stored
    ! in rho, rhou, rhov, rhow, and rhoe.
    !
    if (mach <= -one) then
      !
      ! Supersonic inflow
      !
      ! Set all variables to free stream
      !
      rho = rho_mms
      rhou = rhou_mms
      rhov = rhov_mms
      rhow = rhow_mms
      rhoe = rhoe_mms
      !
    else if (mach > -one .and. mach <= zero) then
      !
      ! Subsonic inflow
      !
      ! Set all variables but pressure to free stream
      !
      rho = rho_mms
      rhou = rhou_mms
      rhov = rhov_mms
      rhow = rhow_mms
      rhoe = pres*gm1_inv + ske_mms
      !
    else if (mach < one .and. mach > zero) then
      !
      ! Subsonic outflow
      !
      ! Extrapolate interior density and velocities to boundary
      ! Compute the total energy using free stream pressure
      ! and interior specific kinetic energy
      !
      rhoe = pres_mms*gm1_inv + ske
      !
    end if
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rhou
    face_cv(n,nmy,2) = rhov
    face_cv(n,nmz,2) = rhow
    face_cv(n,nee,2) = rhoe
    !
  end do
  !
end subroutine get_bc_mms_dirichlet_solution_3d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine get_bc_mms_dirichlet_solution_2d(n1,n2,facenrm,facexyz,face_cv)
  !
  ! #########################################################
  ! #########################################################
  ! ##  MMS Dirichlet (generic freeflow) Boundary Condition
  ! ##    - At each flux point, use the normal velocity
  ! ##      Mach number to determine if it is an inflow
  ! ##      or outflow boundary condition and whether
  ! ##      it is subsonic or supersonic.
  ! #########################################################
  ! #########################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conserved variables at the flux points
  ! FACEXYZ : XYZ coordinates for all the flux points
  !
  !.. Use Statements ..
  use ovar,    only : gam
  use mms_mod, only : mms_solution_at_xyz_gpu
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: facexyz(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,gm1_inv
  real(wp) :: ds,rnx,rny,unx,uny
  real(wp) :: rho,rhou,rhov,rhoe,rhoq2,ske,pres,mach
  real(wp) :: rho_mms,rhou_mms,rhov_mms,rhoe_mms
  real(wp) :: rhoq2_mms,ske_mms,pres_mms
  !
continue
  !
  gm1 = gam - one
  gm1_inv = one/gm1
  !
  call mms_solution_at_xyz_gpu(n1,n2,facexyz,face_cv(:,:,2), &
                               return_nondim=true, &
                               return_cv=true)
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
  do n = n1,n2
    !
    ! Unit normal at this flux point
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds = sqrt( rnx*rnx + rny*rny )
    unx = rnx/ds
    uny = rny/ds
    !
    rho_mms = face_cv(n,nec,2)
    rhou_mms = face_cv(n,nmx,2)
    rhov_mms = face_cv(n,nmy,2)
    rhoe_mms = face_cv(n,nee,2)
    rhoq2_mms = rhou_mms*rhou_mms + rhov_mms*rhov_mms
    ske_mms = half*rhoq2_mms/rho_mms
    pres_mms = gm1*(rhoe_mms - ske_mms)
    !
    ! Compute the Mach number of the normal velocity at this flux point
    !
    rho = face_cv(n,nec,1)
    rhou = face_cv(n,nmx,1)
    rhov = face_cv(n,nmy,1)
    rhoe = face_cv(n,nee,1)
    rhoq2 = rhou*rhou + rhov*rhov
    ske = half*rhoq2/rho
    pres = gm1*(rhoe - ske)
    !
    ! NOTE: The speed of sound (denominator) is gam*pres*rho instead of
    !       gam*pres/rho because the normal velocity (numerator) is rho*vn
    !       instead of just vn.
    !
    mach = (rhou*unx + rhov*uny) / sqrt( max(gam*pres*rho,rndoff) )
    !
    ! If mach >= one, do nothing because in that case we want to extrapolate
    ! everything from the interior and we already have those values stored
    ! in rho, rhou, rhov, rhow, and rhoe.
    !
    if (mach <= -one) then
      !
      ! Supersonic inflow
      !
      ! Set all variables to free stream
      !
      rho = rho_mms
      rhou = rhou_mms
      rhov = rhov_mms
      rhoe = rhoe_mms
      !
    else if (mach > -one .and. mach <= zero) then
      !
      ! Subsonic inflow
      !
      ! Set all variables but pressure to free stream
      !
      rho = rho_mms
      rhou = rhou_mms
      rhov = rhov_mms
      rhoe = pres*gm1_inv + ske_mms
      !
    else if (mach < one .and. mach > zero) then
      !
      ! Subsonic outflow
      !
      ! Extrapolate interior density and velocities to boundary
      ! Compute the total energy using free stream pressure
      ! and interior specific kinetic energy
      !
      rhoe = pres_mms*gm1_inv + ske
      !
    end if
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rhou
    face_cv(n,nmy,2) = rhov
    face_cv(n,nee,2) = rhoe
    !
  end do
  !
end subroutine get_bc_mms_dirichlet_solution_2d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine get_bc_mms_dirichlet_solution_1d(n1,n2,facenrm,facexyz,face_cv)
  !
  ! #########################################################
  ! #########################################################
  ! ##  MMS Dirichlet (generic freeflow) Boundary Condition
  ! ##    - At each flux point, use the normal velocity
  ! ##      Mach number to determine if it is an inflow
  ! ##      or outflow boundary condition and whether
  ! ##      it is subsonic or supersonic.
  ! #########################################################
  ! #########################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conserved variables at the flux points
  ! FACEXYZ : XYZ coordinates for all the flux points
  !
  !.. Use Statements ..
  use ovar,    only : gam
  use mms_mod, only : mms_solution_at_xyz_gpu
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: facexyz(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,gm1_inv
  real(wp) :: ds,rnx,unx
  real(wp) :: rho,rhou,rhoe,rhoq2,ske,pres,mach
  real(wp) :: rho_mms,rhou_mms,rhoe_mms
  real(wp) :: rhoq2_mms,ske_mms,pres_mms
  !
continue
  !
  gm1 = gam - one
  gm1_inv = one/gm1
  !
  call mms_solution_at_xyz_gpu(n1,n2,facexyz,face_cv(:,:,2), &
                               return_nondim=true, &
                               return_cv=true)
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
  do n = n1,n2
    !
    ! Unit normal at this flux point
    !
    unx = facenrm(1,n)/abs(facenrm(1,n))
    !
    rho_mms = face_cv(n,nec,2)
    rhou_mms = face_cv(n,nmx,2)
    rhoe_mms = face_cv(n,nee,2)
    ske_mms = half*rhou_mms*rhou_mms/rho_mms
    pres_mms = gm1*(rhoe_mms - ske_mms)
    !
    ! Compute the Mach number of the normal velocity at this flux point
    !
    rho = face_cv(n,nec,1)
    rhou = face_cv(n,nmx,1)
    rhoe = face_cv(n,nee,1)
    ske = half*rhou*rhou/rho
    pres = gm1*(rhoe - ske)
    !
    ! NOTE: The speed of sound (denominator) is gam*pres*rho instead of
    !       gam*pres/rho because the normal velocity (numerator) is rho*vn
    !       instead of just vn.
    !
    mach = rhou*unx / sqrt( max(gam*pres*rho,rndoff) )
    !
    ! If mach >= one, do nothing because in that case we want to extrapolate
    ! everything from the interior and we already have those values stored
    ! in rho, rhou, and rhoe.
    !
    if (mach <= -one) then
      !
      ! Supersonic inflow
      !
      ! Set all variables to free stream
      !
      rho = rho_mms
      rhou = rhou_mms
      rhoe = rhoe_mms
      !
    else if (mach > -one .and. mach <= zero) then
      !
      ! Subsonic inflow
      !
      ! Set all variables but pressure to free stream
      !
      rho = rho_mms
      rhou = rhou_mms
      rhoe = pres*gm1_inv + ske_mms
      !
    else if (mach < one .and. mach > zero) then
      !
      ! Subsonic outflow
      !
      ! Extrapolate interior density and velocities to boundary
      ! Compute the total energy using free stream pressure
      ! and interior specific kinetic energy
      !
      rhoe = pres_mms*gm1_inv + ske
      !
    end if
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rhou
    face_cv(n,nee,2) = rhoe
    !
  end do
  !
end subroutine get_bc_mms_dirichlet_solution_1d
!
!###############################################################################
!
end module bc_mms_mod
