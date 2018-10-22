module bc_riemann_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nee
  use eqn_idx, only : ncs,npres
  !
  implicit none
  !
  private
  !
  public :: get_bc_characteristic_solution
  !
contains
!
!###############################################################################
!
#ifdef DISABLE_FOR_GPU
pure &
#endif
subroutine get_bc_characteristic_solution(n1,n2,facenrm,face_cv,face_gp, &
                                                pv_in,pv_fs,relax)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Characteristic boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! CV_IN   : Boundary condition flow variables in conservative form
  ! RELAX   : Boundary condition relaxation parameters
  !
  !.. Use Statements ..
  use ovar, only : relaxation_t
  !
  use bc_relax_mod, only : get_relaxed_pv_bc
  !
  !.. Formal Arguments ..
  integer,            intent(in) :: n1,n2
  real(wp),           intent(in) :: facenrm(:,:)
  real(wp),           intent(in) :: pv_in(:)
  real(wp),           intent(in) :: pv_fs(:)
  type(relaxation_t), intent(in) :: relax
  real(wp),        intent(inout) :: face_cv(:,:,:)
  real(wp),           intent(in) :: face_gp(:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
  !.. Local Arrays ..
  real(wp) :: pv_bc(1:size(pv_in))
  !
  !.. Local Parameters ..
  logical(lk), parameter :: veln_is_int = true
 !logical(lk), parameter :: veln_is_int = fals
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  ! Determine the relaxation coefficient
  !
  call get_relaxed_pv_bc(relax,pv_fs,pv_in,pv_rlx=pv_bc)
  !
  if (nr == 3) then
    call get_bc_characteristic_solution_3d(n1,n2,facenrm,face_cv,face_gp, &
                                                 pv_bc,veln_is_int)
  else if (nr == 2) then
    call get_bc_characteristic_solution_2d(n1,n2,facenrm,face_cv,face_gp, &
                                                 pv_bc,veln_is_int)
  else
    call get_bc_characteristic_solution_1d(n1,n2,facenrm,face_cv,face_gp, &
                                                 pv_bc,veln_is_int)
  end if
  !
end subroutine get_bc_characteristic_solution
!
!###############################################################################
!
pure &
subroutine get_bc_characteristic_solution_3d(n1,n2,facenrm,face_cv,face_gp, &
                                                   pv_bc,veln_is_int)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Characteristic boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_BC   : Boundary condition flow variables in conservative form
  ! RELAX   : Boundary condition relaxation parameters
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  logical(lk), intent(in) :: veln_is_int
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: gm1,one_over_gm1,two_over_gm1
  real(wp) :: rho,vx,vy,vz,vn,pres,sc,rp,rm,unc,ac,vel
  real(wp) :: rho_inv,velx,vely,velz,q2
  real(wp) :: rho_int,vx_int,vy_int,vz_int,pres_int,aspd_int,vn_int
  real(wp) :: rho_bc,vx_bc,vy_bc,vz_bc,pres_bc,aspd_bc,vn_bc
  !
continue
  !
  ! Precompute gamma-1 and its inverse
  !
  gm1 = gam - one
  one_over_gm1 = one / gm1
  two_over_gm1 = two / gm1
  !
  rho_bc = pv_bc(nec)
  vx_bc = pv_bc(nmx)
  vy_bc = pv_bc(nmy)
  vz_bc = pv_bc(nmz)
  pres_bc = pv_bc(nee)
  !
  aspd_bc = sqrt( max(gam*pres_bc/rho_bc,rndoff) ) ! Ext. speed of sound
  !
  ! Compute the boundary solution at each flux point
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp)
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
    ! Get the host cell variables at this flux point
    !
    rho_int = face_cv(n,nec,1)
    rho_inv = one/rho_int
    vx_int = face_cv(n,nmx,1)*rho_inv
    vy_int = face_cv(n,nmy,1)*rho_inv
    vz_int = face_cv(n,nmz,1)*rho_inv
    pres_int = face_gp(n,npres,1)
    !
    ! Compute the speed of sound in the host cell
    !
    aspd_int = sqrt( max(gam*pres_int/rho_int,rndoff) ) ! Int. speed of sound
    !
    ! Compute the normal velocity at this flux point
    !
    ! Interior normal velocity
    vn_int = vx_int*unx + vy_int*uny + vz_int*unz
    ! Exterior normal velocity
    vn_bc = vx_bc*unx + vy_bc*uny + vz_bc*unz
    !
    ! Apply boundary conditions based on the
    ! normal Mach number to the boundary face
    !
    rp = vn_int + two_over_gm1*aspd_int ! Outgoing Riemann invariant
    rm = vn_bc  - two_over_gm1*aspd_bc  ! Incoming Reimann invariant
    !
    ! Supersonic inflow boundary
    !    Set R+ so that 0.5*(R+ + R-) = vn_bc
    !
    if (vn_bc <= -aspd_bc) then
      rp = vn_bc + two_over_gm1*aspd_bc
    end if
    !
    ! Supersonic outflow boundary
    !    Set R- so that 0.5*(R+ + R-) = vn_int
    !
    if (vn_int >= aspd_int) then
      rm = vn_int - two_over_gm1*aspd_int
    end if
    !
    unc = half*(rp+rm)
    ac  = one4*gm1*(rp-rm)
    !
    vel = merge( vn_int , unc , veln_is_int )
    if (vel > zero) then ! outflow boundary
      vx = vx_int
      vy = vy_int
      vz = vz_int
      vn = vn_int
      sc = pres_int / (rho_int**gam)
    else                 ! inflow boundary
      vx = vx_bc
      vy = vy_bc
      vz = vz_bc
      vn = vn_bc
      sc = pres_bc / (rho_bc**gam)
    end if
    !
    rho = ( ac * ac / (gam*sc) )**one_over_gm1
    pres = sc * (rho**gam)
    !
    velx = vx + (unc - vn)*unx
    vely = vy + (unc - vn)*uny
    velz = vz + (unc - vn)*unz
    !
   !q2 = vx*vx + vy*vy + vz*vz
    q2 = velx*velx + vely*vely + velz*velz
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rho*velx
    face_cv(n,nmy,2) = rho*vely
    face_cv(n,nmz,2) = rho*velz
    face_cv(n,nee,2) = pres*one_over_gm1 + half*rho*q2
    !
  end do
  !
  ! Fix the species mass fractions at the boundary.
  ! NOTE: The characteristic boundary condition is not really valid for gas
  !       mixtures because the Riemann invariants assume a frozen gamma.
  ! NOTE: We need to output a WARNING of some kind before the simulation
  !       starts to warn the user that the characteristic boundary condition
  !       is not really valid for gas mixtures.
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv,pv_bc)
    do m = 1,ncs
      do n = n1,n2
        face_cv(n,m,2) = pv_bc(m) * face_cv(n,nec,2)
      end do
    end do
  end if
  !
end subroutine get_bc_characteristic_solution_3d
!
!###############################################################################
!
pure &
subroutine get_bc_characteristic_solution_2d(n1,n2,facenrm,face_cv,face_gp, &
                                                   pv_bc,veln_is_int)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Characteristic boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_BC   : Boundary condition flow variables in conservative form
  ! RELAX   : Boundary condition relaxation parameters
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  logical(lk), intent(in) :: veln_is_int
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: ds,rnx,rny,unx,uny
  real(wp) :: gm1,one_over_gm1,two_over_gm1
  real(wp) :: rho,vx,vy,vn,pres,sc,rp,rm,unc,ac,vel
  real(wp) :: rho_inv,velx,vely,q2
  real(wp) :: rho_int,vx_int,vy_int,pres_int,aspd_int,vn_int
  real(wp) :: rho_bc,vx_bc,vy_bc,pres_bc,aspd_bc,vn_bc
  !
continue
  !
  ! Precompute gamma-1 and its inverse
  !
  gm1 = gam - one
  one_over_gm1 = one / gm1
  two_over_gm1 = two / gm1
  !
  rho_bc = pv_bc(nec)
  vx_bc = pv_bc(nmx)
  vy_bc = pv_bc(nmy)
  pres_bc = pv_bc(nee)
  !
  aspd_bc = sqrt( max(gam*pres_bc/rho_bc,rndoff) ) ! Ext. speed of sound
  !
  ! Compute the boundary solution at each flux point
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp)
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
    ! Get the host cell variables at this flux point
    !
    rho_int = face_cv(n,nec,1)
    rho_inv = one/rho_int
    vx_int = face_cv(n,nmx,1)*rho_inv
    vy_int = face_cv(n,nmy,1)*rho_inv
    pres_int = face_gp(n,npres,1)
    !
    ! Compute the speed of sound in the host cell
    !
    aspd_int = sqrt( max(gam*pres_int/rho_int,rndoff) ) ! Int. speed of sound
    !
    ! Compute the normal velocity at this flux point
    !
    ! Interior normal velocity
    vn_int = vx_int*unx + vy_int*uny
    ! Exterior normal velocity
    vn_bc = vx_bc*unx + vy_bc*uny
    !
    ! Apply boundary conditions based on the
    ! normal Mach number to the boundary face
    !
    rp = vn_int + two_over_gm1*aspd_int ! Outgoing Riemann invariant
    rm = vn_bc  - two_over_gm1*aspd_bc  ! Incoming Reimann invariant
    !
    ! Supersonic inflow boundary
    !    Set R+ so that 0.5*(R+ + R-) = vn_bc
    !
    if (vn_bc <= -aspd_bc) then
      rp = vn_bc + two_over_gm1*aspd_bc
    end if
    !
    ! Supersonic outflow boundary
    !    Set R- so that 0.5*(R+ + R-) = vn_int
    !
    if (vn_int >= aspd_int) then
      rm = vn_int - two_over_gm1*aspd_int
    end if
    !
    unc = half*(rp+rm)
    ac  = one4*gm1*(rp-rm)
    !
    vel = merge( vn_int , unc , veln_is_int )
    if (vel > zero) then ! outflow boundary
      vx = vx_int
      vy = vy_int
      vn = vn_int
      sc = pres_int / (rho_int**gam)
    else                 ! inflow boundary
      vx = vx_bc
      vy = vy_bc
      vn = vn_bc
      sc = pres_bc / (rho_bc**gam)
    end if
    !
    rho = ( ac * ac / (gam*sc) )**one_over_gm1
    pres = sc * (rho**gam)
    !
    velx = vx + (unc - vn)*unx
    vely = vy + (unc - vn)*uny
    !
   !q2 = vx*vx + vy*vy
    q2 = velx*velx + vely*vely
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rho*velx
    face_cv(n,nmy,2) = rho*vely
    face_cv(n,nee,2) = pres*one_over_gm1 + half*rho*q2
    !
  end do
  !
  ! Fix the species mass fractions at the boundary.
  ! NOTE: The characteristic boundary condition is not really valid for gas
  !       mixtures because the Riemann invariants assume a frozen gamma.
  ! NOTE: We need to output a WARNING of some kind before the simulation
  !       starts to warn the user that the characteristic boundary condition
  !       is not really valid for gas mixtures.
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv,pv_bc)
    do m = 1,ncs
      do n = n1,n2
        face_cv(n,m,2) = pv_bc(m) * face_cv(n,nec,2)
      end do
    end do
  end if
  !
end subroutine get_bc_characteristic_solution_2d
!
!###############################################################################
!
pure &
subroutine get_bc_characteristic_solution_1d(n1,n2,facenrm,face_cv,face_gp, &
                                                   pv_bc,veln_is_int)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Characteristic boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_BC   : Boundary condition flow variables in conservative form
  ! RELAX   : Boundary condition relaxation parameters
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  logical(lk), intent(in) :: veln_is_int
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: ds,rnx,unx
  real(wp) :: gm1,one_over_gm1,two_over_gm1
  real(wp) :: rho,vx,vn,pres,sc,rp,rm,unc,ac,vel
  real(wp) :: rho_inv,velx,q2
  real(wp) :: rho_int,vx_int,pres_int,aspd_int,vn_int
  real(wp) :: rho_bc,vx_bc,pres_bc,aspd_bc,vn_bc
  !
continue
  !
  ! Precompute gamma-1 and its inverse
  !
  gm1 = gam - one
  one_over_gm1 = one / gm1
  two_over_gm1 = two / gm1
  !
  rho_bc = pv_bc(nec)
  vx_bc = pv_bc(nmx)
  pres_bc = pv_bc(nee)
  !
  aspd_bc = sqrt( max(gam*pres_bc/rho_bc,rndoff) ) ! Ext. speed of sound
  !
  ! Compute the boundary solution at each flux point
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp)
  do n = n1,n2
    !
    ! Unit normal at this flux point
    !
    rnx = facenrm(1,n)
    ds = sqrt( rnx*rnx )
    unx = rnx/ds
    !
    ! Get the host cell variables at this flux point
    !
    rho_int = face_cv(n,nec,1)
    rho_inv = one/rho_int
    vx_int = face_cv(n,nmx,1)*rho_inv
    pres_int = face_gp(n,npres,1)
    !
    ! Compute the speed of sound in the host cell
    !
    aspd_int = sqrt( max(gam*pres_int/rho_int,rndoff) ) ! Int. speed of sound
    !
    ! Compute the normal velocity at this flux point
    !
    ! Interior normal velocity
    vn_int = vx_int*unx
    ! Exterior normal velocity
    vn_bc = vx_bc*unx
    !
    ! Apply boundary conditions based on the
    ! normal Mach number to the boundary face
    !
    rp = vn_int + two_over_gm1*aspd_int ! Outgoing Riemann invariant
    rm = vn_bc  - two_over_gm1*aspd_bc  ! Incoming Reimann invariant
    !
    ! Supersonic inflow boundary
    !    Set R+ so that 0.5*(R+ + R-) = vn_bc
    !
    if (vn_bc <= -aspd_bc) then
      rp = vn_bc + two_over_gm1*aspd_bc
    end if
    !
    ! Supersonic outflow boundary
    !    Set R- so that 0.5*(R+ + R-) = vn_int
    !
    if (vn_int >= aspd_int) then
      rm = vn_int - two_over_gm1*aspd_int
    end if
    !
    unc = half*(rp+rm)
    ac  = one4*gm1*(rp-rm)
    !
    vel = merge( vn_int , unc , veln_is_int )
    if (vel > zero) then ! outflow boundary
      vx = vx_int
      vn = vn_int
      sc = pres_int / (rho_int**gam)
    else                 ! inflow boundary
      vx = vx_bc
      vn = vn_bc
      sc = pres_bc / (rho_bc**gam)
    end if
    !
    rho = ( ac * ac / (gam*sc) )**one_over_gm1
    pres = sc * (rho**gam)
    !
    velx = vx + (unc - vn)*unx
    !
   !q2 = vx*vx
    q2 = velx*velx
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rho*velx
    face_cv(n,nee,2) = pres*one_over_gm1 + half*rho*q2
    !
  end do
  !
  ! Fix the species mass fractions at the boundary.
  ! NOTE: The characteristic boundary condition is not really valid for gas
  !       mixtures because the Riemann invariants assume a frozen gamma.
  ! NOTE: We need to output a WARNING of some kind before the simulation
  !       starts to warn the user that the characteristic boundary condition
  !       is not really valid for gas mixtures.
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv,pv_bc)
    do m = 1,ncs
      do n = n1,n2
        face_cv(n,m,2) = pv_bc(m) * face_cv(n,nec,2)
      end do
    end do
  end if
  !
end subroutine get_bc_characteristic_solution_1d
!
!###############################################################################
!
end module bc_riemann_mod
