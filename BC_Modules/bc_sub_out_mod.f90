module bc_sub_out_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx_mod, only : nec,nmx,nmy,nmz,nee
  use eqn_idx_mod, only : ncs,npres,ngamm,nenth
  !
  implicit none
  !
  private
  !
  public :: get_bc_sub_outflow_solution
  !
contains
!
!###############################################################################
!
#ifdef DISABLE_FOR_GPU
pure &
#endif
subroutine get_bc_sub_outflow_solution(n1,n2,facenrm,face_cv,face_gp, &
                                             pv_in,gp_in,relax)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Subsonic outflow boundary condition
  ! ##    - Extrapolate density
  ! ##    - Extrapolate velocities
  ! ##    - Fix pressure to free stream
  ! #####################################################
  ! #####################################################
  !
  !.. Use Statements ..
  use ovar, only : gam
  use ovar, only : relaxation_t
  !
  use bc_relax_mod, only : get_relax_coefficient
  !
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  !
  !.. Formal Arguments ..
  integer,            intent(in) :: n1,n2
  real(wp),           intent(in) :: facenrm(:,:)
  real(wp),           intent(in) :: pv_in(:)
  real(wp),           intent(in) :: gp_in(:)
  real(wp),        intent(inout) :: face_cv(:,:,:)
  real(wp),        intent(inout) :: face_gp(:,:,:)
  type(relaxation_t), intent(in) :: relax
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: rlx,gam_over_gm1
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  call get_relax_coefficient(relax,rlx)
  !
  ! Load enthalpy into the gas properties auxiliary variable location for
  ! the interior state, and load gamma into the exterior state
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP PRESENT(face_gp)
    do n = n1,n2
      face_gp(n,ngamm,2) = face_gp(n,ngamm,1)
    end do
  else
    gam_over_gm1 = gam / (gam - one)
    !$ACC PARALLEL LOOP PRESENT(face_cv,face_gp)
    do n = n1,n2
      face_gp(n,nenth,1) = gam_over_gm1 * face_gp(n,npres,1) / face_cv(n,nec,1)
      ! Load gamma into the right state because (nenth == ngamm) for a
      ! calorically perfect gas.
      face_gp(n,ngamm,2) = gam
    end do
  end if
  !
  if (nr == 3) then
    call get_bc_sub_outflow_solution_3d(n1,n2,facenrm,face_cv,face_gp, &
                                              pv_in,gp_in,rlx)
  else if (nr == 2) then
    call get_bc_sub_outflow_solution_2d(n1,n2,facenrm,face_cv,face_gp, &
                                              pv_in,gp_in,rlx)
  else
    call get_bc_sub_outflow_solution_1d(n1,n2,facenrm,face_cv,face_gp, &
                                              pv_in,gp_in,rlx)
  end if
  !
end subroutine get_bc_sub_outflow_solution
!
!###############################################################################
!
#ifdef DISABLE_FOR_GPU
pure &
#endif
subroutine get_bc_sub_outflow_solution_3d(n1,n2,facenrm,face_cv,face_gp, &
                                                pv_bc,gp_bc,rlx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Subsonic outflow boundary condition
  ! ##    - Extrapolate density
  ! ##    - Extrapolate velocities
  ! ##    - Fix pressure to free stream
  ! #####################################################
  ! #####################################################
  !
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_BC   : Boundary condition flow variables in primitive form
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp),    intent(in) :: gp_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp),    intent(in) :: rlx
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: rho,rho_inv,vx,vy,vz,vn,q2,aspd
  real(wp) :: pres,gamm,enth,mach,ptot,total_cond
  real(wp) :: aspd_bc,mach_bc,pres_bc,ptot_bc
  real(wp) :: rho_ext,pres_ext
  !
continue
  !
  if (ncs > 0) then
    gamm = gp_bc(ngamm)
  else
    gamm = gam
  end if
  !
  pres_bc = pv_bc(nee)
  aspd_bc = sqrt( max( gamm*pres_bc/pv_bc(nec) , rndoff ) )
  mach_bc = sqrt( pv_bc(nmx)*pv_bc(nmx) + &
                  pv_bc(nmy)*pv_bc(nmy) + &
                  pv_bc(nmz)*pv_bc(nmz) ) / aspd_bc
  total_cond = one + half*(gamm-one)*mach_bc*mach_bc
  ptot_bc = pres_bc * total_cond**(gamm/(gamm-one))
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    rnz = facenrm(3,n)
    ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
    unx = rnx/ds ! These point outside the domain
    uny = rny/ds
    unz = rnz/ds
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    vz = face_cv(n,nmz,1)*rho_inv
    !
    vn = vx*unx + vy*uny + vz*unz
    q2 = vx*vx + vy*vy + vz*vz
    !
    pres = face_gp(n,npres,1)
    gamm = face_gp(n,ngamm,2) ! this was loaded into the right state
    enth = face_gp(n,nenth,1)
    !
    aspd = sqrt(max(gamm*pres*rho_inv,rndoff))
    mach = abs(vn)/aspd
    !
    ! If the normal velocity from the interior says this is an inflow
    ! boundary, set the ghost velocity to the magnitude of the interior
    ! velocity multiplied by the unit normal
    ! NOTE: q2 will remain the same
    !
    if (vn < zero) then
      vx = sqrt(q2)*unx
      vy = sqrt(q2)*uny
      vz = sqrt(q2)*unz
    end if
    !
    ! Set the boundary pressure to the reference pressure if the normal
    ! velocity from the interior is subsonic, otherwise set the boundary
    ! pressure to the reference total pressure
    ! NOTE: This computes a boundary density using the assumption that
    !       the boundary temperature is equal to the interior temperature.
    !
    pres_ext = (one-rlx)*pres + rlx*pres_bc
    rho_ext = rho*pres_ext/pres
   !enth_ext = gamm*pres_ext / ((gamm-one)*rho_ext)
    !
    if (mach >= one) then
      total_cond = one + half*(gamm-one)*mach*mach
      ptot = pres * total_cond**(gamm/(gamm-one))
      pres_ext = (one-rlx)*ptot + rlx*ptot_bc
      rho_ext = rho*pres_ext/ptot
    end if
    !
    ! THESE COMPUTE THE BOUNDARY STATES USING
    !     INTERIOR : enthalpy (temperature) and velocities
    !     EXTERIOR : pressure
    !
    face_cv(n,nec,2) = rho_ext
    face_cv(n,nmx,2) = rho_ext*vx
    face_cv(n,nmy,2) = rho_ext*vy
    face_cv(n,nmz,2) = rho_ext*vz
    face_cv(n,nee,2) = rho_ext*enth - pres_ext + half*rho_ext*q2
    !
    ! THESE COMMENTED OUT LINES COMPUTE THE BOUNDARY STATES USING
    !     INTERIOR : rho and velocities
    !     EXTERIOR : pressure
    ! NOTE: I THINK THE ENTHALPY/TEMPERATURE BASED STATES WORK BETTER???
    !
   !face_cv(n,nec,2) = rho
   !face_cv(n,nmx,2) = rho*vx
   !face_cv(n,nmy,2) = rho*vy
   !face_cv(n,nmz,2) = rho*vz
   !face_cv(n,nee,2) = pres_ext/(gamm-one) + half*rho*q2
    !
  end do
  !
  ! Extrapolate the interior mass fractions to the exterior.
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv)
    do m = 1,ncs
      do n = n1,n2
        face_cv(n,m,2) = face_cv(n,nec,2)*face_cv(n,m,1)/face_cv(n,nec,1)
      end do
    end do
  end if
  !
end subroutine get_bc_sub_outflow_solution_3d
!
!###############################################################################
!
#ifdef DISABLE_FOR_GPU
pure &
#endif
subroutine get_bc_sub_outflow_solution_2d(n1,n2,facenrm,face_cv,face_gp, &
                                                pv_bc,gp_bc,rlx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Subsonic outflow boundary condition
  ! ##    - Extrapolate density
  ! ##    - Extrapolate velocities
  ! ##    - Fix pressure to free stream
  ! #####################################################
  ! #####################################################
  !
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_BC   : Boundary condition flow variables in primitive form
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp),    intent(in) :: gp_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp),    intent(in) :: rlx
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: ds,rnx,rny,unx,uny
  real(wp) :: rho,rho_inv,vx,vy,vn,q2,aspd
  real(wp) :: pres,gamm,enth,mach,ptot,total_cond
  real(wp) :: aspd_bc,mach_bc,pres_bc,ptot_bc
  real(wp) :: rho_ext,pres_ext
  !
continue
  !
  if (ncs > 0) then
    gamm = gp_bc(ngamm)
  else
    gamm = gam
  end if
  !
  pres_bc = pv_bc(nee)
  aspd_bc = sqrt( max( gamm*pres_bc/pv_bc(nec) , rndoff ) )
  mach_bc = sqrt( pv_bc(nmx)*pv_bc(nmx) + &
                  pv_bc(nmy)*pv_bc(nmy) ) / aspd_bc
  total_cond = one + half*(gamm-one)*mach_bc*mach_bc
  ptot_bc = pres_bc * total_cond**(gamm/(gamm-one))
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp)
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
    !
    vn = vx*unx + vy*uny
    q2 = vx*vx + vy*vy
    !
    pres = face_gp(n,npres,1)
    gamm = face_gp(n,ngamm,2) ! this was loaded into the right state
    enth = face_gp(n,nenth,1)
    !
    aspd = sqrt(max(gamm*pres*rho_inv,rndoff))
    mach = abs(vn)/aspd
    !
    ! If the normal velocity from the interior says this is an inflow
    ! boundary, set the ghost velocity to the magnitude of the interior
    ! velocity multiplied by the unit normal
    ! NOTE: q2 will remain the same
    !
    if (vn < zero) then
      vx = sqrt(q2)*unx
      vy = sqrt(q2)*uny
    end if
    !
    ! Set the boundary pressure to the reference pressure if the normal
    ! velocity from the interior is subsonic, otherwise set the boundary
    ! pressure to the reference total pressure
    !
    pres_ext = (one-rlx)*pres + rlx*pres_bc
    rho_ext = rho*pres_ext/pres
    !
    if (mach >= one) then
      total_cond = one + half*(gamm-one)*mach*mach
      ptot = pres * total_cond**(gamm/(gamm-one))
      pres_ext = (one-rlx)*ptot + rlx*ptot_bc
      rho_ext = rho*pres_ext/ptot
    end if
    !
    face_cv(n,nec,2) = rho_ext
    face_cv(n,nmx,2) = rho_ext*vx
    face_cv(n,nmy,2) = rho_ext*vy
    face_cv(n,nee,2) = rho_ext*enth - pres_ext + half*rho_ext*q2
    !
  end do
  !
  ! Extrapolate the interior mass fractions to the exterior.
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv)
    do m = 1,ncs
      do n = n1,n2
        face_cv(n,m,2) = face_cv(n,nec,2)*face_cv(n,m,1)/face_cv(n,nec,1)
      end do
    end do
  end if
  !
end subroutine get_bc_sub_outflow_solution_2d
!
!###############################################################################
!
#ifdef DISABLE_FOR_GPU
pure &
#endif
subroutine get_bc_sub_outflow_solution_1d(n1,n2,facenrm,face_cv,face_gp, &
                                                pv_bc,gp_bc,rlx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Subsonic outflow boundary condition
  ! ##    - Extrapolate density
  ! ##    - Extrapolate velocities
  ! ##    - Fix pressure to free stream
  ! #####################################################
  ! #####################################################
  !
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_BC   : Boundary condition flow variables in primitive form
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp),    intent(in) :: gp_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp),    intent(in) :: rlx
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: ds,rnx,unx
  real(wp) :: rho,rho_inv,vx,vn,q2,aspd
  real(wp) :: pres,gamm,enth,mach,ptot,total_cond
  real(wp) :: aspd_bc,mach_bc,pres_bc,ptot_bc
  real(wp) :: rho_ext,pres_ext
  !
continue
  !
  if (ncs > 0) then
    gamm = gp_bc(ngamm)
  else
    gamm = gam
  end if
  !
  pres_bc = pv_bc(nee)
  aspd_bc = sqrt( max( gamm*pres_bc/pv_bc(nec) , rndoff ) )
  mach_bc = sqrt( pv_bc(nmx)*pv_bc(nmx) ) / aspd_bc
  total_cond = one + half*(gamm-one)*mach_bc*mach_bc
  ptot_bc = pres_bc * total_cond**(gamm/(gamm-one))
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    ds = sqrt( rnx*rnx )
    unx = rnx/ds
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    !
    vn = vx*unx
    q2 = vx*vx
    !
    pres = face_gp(n,npres,1)
    gamm = face_gp(n,ngamm,2) ! this was loaded into the right state
    enth = face_gp(n,nenth,1)
    !
    aspd = sqrt(max(gamm*pres*rho_inv,rndoff))
    mach = abs(vn)/aspd
    !
    ! If the normal velocity from the interior says this is an inflow
    ! boundary, set the ghost velocity to the magnitude of the interior
    ! velocity multiplied by the unit normal
    ! NOTE: q2 will remain the same
    !
    if (vn < zero) then
      vx = sqrt(q2)*unx
    end if
    !
    ! Set the boundary pressure to the reference pressure if the normal
    ! velocity from the interior is subsonic, otherwise set the boundary
    ! pressure to the reference total pressure
    !
    pres_ext = (one-rlx)*pres + rlx*pres_bc
    rho_ext = rho*pres_ext/pres
    !
    if (mach >= one) then
      total_cond = one + half*(gamm-one)*mach*mach
      ptot = pres * total_cond**(gamm/(gamm-one))
      pres_ext = (one-rlx)*ptot + rlx*ptot_bc
      rho_ext = rho*pres_ext/ptot
    end if
    !
    face_cv(n,nec,2) = rho_ext
    face_cv(n,nmx,2) = rho_ext*vx
    face_cv(n,nee,2) = rho_ext*enth - pres_ext + half*rho_ext*q2
    !
  end do
  !
  ! Extrapolate the interior mass fractions to the exterior.
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv)
    do m = 1,ncs
      do n = n1,n2
        face_cv(n,m,2) = face_cv(n,nec,2)*face_cv(n,m,1)/face_cv(n,nec,1)
      end do
    end do
  end if
  !
end subroutine get_bc_sub_outflow_solution_1d
!
!###############################################################################
!
end module bc_sub_out_mod
