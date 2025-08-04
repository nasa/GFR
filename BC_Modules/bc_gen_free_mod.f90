module bc_gen_free_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx_mod, only : nec,nmx,nmy,nmz,nee
  use eqn_idx_mod, only : ncs,npres,ntemp,ncp,ngamm,nenth
  !
  implicit none
  !
  private
  !
  public :: get_bc_generic_freeflow_solution
  !
contains
!
!###############################################################################
!
subroutine get_bc_generic_freeflow_solution(n1,n2,facenrm,face_cv,face_gp, &
                                            pv_in,gp_in,pv_fs,gp_fs, &
                                            relax,newton_failures)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Generic freeflow boundary condition
  ! ##    - At each flux point, use the normal velocity
  ! ##      Mach number to determine if it is an inflow
  ! ##      or outflow boundary condition and whether
  ! ##      it is subsonic or supersonic.
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conserved variables at the flux points
  ! FACE_GP : Gas property variables at the flux points
  ! PV_IN   : Primitive    variables specified for this BC
  ! GP_IN   : Gas property variables specified for this BC
  ! PV_FS   : Primitive    variables for the freestream conditions
  ! GP_FS   : Gas property variables for the freestream conditions
  ! NEWTON_FAILURES : Number of Newton iterations that failed to converge
  !
  !.. Use Statements ..
  use ovar, only : relaxation_t
  !
  use bc_relax_mod, only : get_relaxed_total_bc
  !
  !.. Formal Arguments ..
  integer,            intent(in) :: n1,n2
  real(wp),           intent(in) :: facenrm(:,:)
  real(wp),           intent(in) :: pv_in(:)
  real(wp),           intent(in) :: gp_in(:)
  real(wp),           intent(in) :: pv_fs(:)
  real(wp),           intent(in) :: gp_fs(:)
  type(relaxation_t), intent(in) :: relax
  real(wp),        intent(inout) :: face_cv(:,:,:)
  real(wp),        intent(inout) :: face_gp(:,:,:)
  integer,         intent(inout) :: newton_failures
  !
  !.. Local Scalars ..
  integer  :: nr
  real(wp) :: rlx
  !
  !.. Local Arrays ..
  real(wp) :: cv_bc(1:size(pv_in))
  real(wp) :: pv_bc(1:size(pv_in))
  real(wp) :: gp_bc(1:size(gp_in))
  !
continue
  !
  newton_failures = 0
  !
  nr = size(facenrm,dim=1)
  !
  ! NEED TO GET PV_BC, WITH:
  !           PV_BC(1:NCS)   = MASS_FRACTIONS
  !           PV_BC(NEC)     = TOTAL DENSITY
  !           PV_BC(NMB:NME) = VELOCITIES
  !           PV_BC(NEE)     = TOTAL PRESSURE
  !
  call get_relaxed_total_bc(relax,pv_fs,gp_fs,pv_in,gp_in,cv_rlx=cv_bc, &
                                                          pv_rlx=pv_bc, &
                                                          gp_rlx=gp_bc, &
                                                          return_rlx=rlx)
  !
  if (ncs > 0) then
    !
    ! This inflow condition is very similar to get_bc_sub_inflow_sol_tpg_[123]d
    ! except that it checks the interior velocity normal to the inflow plane,
    ! and switches to a fixed inflow boundary if the interior normal velocity is
    ! determined to be supersonic.
    !
    if (nr == 3) then
      call get_bc_gen_freeflow_sol_tpg_3d(n1,n2,facenrm,face_cv,face_gp, &
                                          cv_bc,pv_bc,gp_bc,rlx,newton_failures)
    else if (nr == 2) then
      call get_bc_gen_freeflow_sol_tpg_2d(n1,n2,facenrm,face_cv,face_gp, &
                                          cv_bc,pv_bc,gp_bc,rlx,newton_failures)
    else
      call get_bc_gen_freeflow_sol_tpg_1d(n1,n2,facenrm,face_cv,face_gp, &
                                          cv_bc,pv_bc,gp_bc,rlx,newton_failures)
    end if
    !
  else
    !
    ! This inflow condition is very similar to get_bc_sub_inflow_sol_cpg_[123]d
    ! except that it checks the interior velocity normal to the inflow plane,
    ! and switches to a fixed inflow boundary if the interior normal velocity is
    ! determined to be supersonic.
    !
    if (nr == 3) then
      call get_bc_gen_freeflow_sol_cpg_3d(n1,n2,facenrm,face_cv,face_gp, &
                                          cv_bc,pv_bc,gp_bc,rlx)
    else if (nr == 2) then
      call get_bc_gen_freeflow_sol_cpg_2d(n1,n2,facenrm,face_cv,face_gp, &
                                          cv_bc,pv_bc,gp_bc,rlx)
    else
      call get_bc_gen_freeflow_sol_cpg_1d(n1,n2,facenrm,face_cv,face_gp, &
                                          cv_bc,pv_bc,gp_bc,rlx)
    end if
    !
  end if
  !
end subroutine get_bc_generic_freeflow_solution
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine get_bc_gen_freeflow_sol_cpg_3d(n1,n2,facenrm,face_cv,face_gp, &
                                          cv_bc,pv_bc,gp_bc,rlx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Generic freeflow boundary condition
  ! ##    - At each flux point, use the normal velocity
  ! ##      Mach number to determine if it is an inflow
  ! ##      or outflow boundary condition and whether
  ! ##      it is subsonic or supersonic.
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conserved variables at the flux points
  ! FACE_GP : Gas property variables at the flux points
  ! CV_BC   : Conserved    variables specified for this BC
  ! PV_BC   : Primitive    variables specified for this BC
  ! GP_BC   : Gas property variables specified for this BC
  ! RLX     : Current relaxation value for this BC
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: cv_bc(:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp),    intent(in) :: gp_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp),    intent(in) :: rlx
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: rho,vel,pres,total_cond
  real(wp) :: pres_bc,rtot_bc,ptot_bc
  real(wp) :: vmag,dir_cosx,dir_cosy,dir_cosz
  real(wp) :: rho_int,rho_inv,vx,vy,vz,q2
  real(wp) :: vn_int,pres_int,aspd_int,mach_int
  real(wp) :: gm1,one_over_gm1,one_over_1mg,gam_over_1mg
  !
continue
  !
  gm1 = gam - one
  one_over_gm1 = one / gm1
  one_over_1mg = -one_over_gm1
  gam_over_1mg = gam * one_over_1mg
  !
  ! Get the boundary static pressure
  !
  pres_bc = gp_bc(npres)
  !
  ! Get the boundary total density and total pressure
  ! NOTE: These total conditions have been saved into static
  !       density and static pressure locations of pv_bc.
  !
  rtot_bc = pv_bc(nec)
  ptot_bc = pv_bc(nee)
  !
  ! Get the directional cosines of the specified boundary condition velocity
  !
  vmag = sqrt( pv_bc(nmx)*pv_bc(nmx) + &
               pv_bc(nmy)*pv_bc(nmy) + &
               pv_bc(nmz)*pv_bc(nmz) )
  !
  dir_cosx = pv_bc(nmx) / max(vmag,rndoff)
  dir_cosy = pv_bc(nmy) / max(vmag,rndoff)
  dir_cosz = pv_bc(nmz) / max(vmag,rndoff)
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp,cv_bc)
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
    rho_int = face_cv(n,nec,1)
    rho_inv = one/rho_int
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    vz = face_cv(n,nmz,1)*rho_inv
    !
    vn_int = vx*unx + vy*uny + vz*unz
    q2 = vx*vx + vy*vy + vz*vz
    !
    pres_int = face_gp(n,npres,1)
    !
    aspd_int = sqrt(max(gam*pres_int*rho_inv,rndoff))
    mach_int = vn_int/aspd_int
    !
    if (mach_int >= one) then
      !
      ! This is a supersonic outflow boundary.
      ! Extrapolate all variables to the boundary.
      !
      face_cv(n,nec,2) = face_cv(n,nec,1)
      face_cv(n,nmx,2) = face_cv(n,nmx,1)
      face_cv(n,nmy,2) = face_cv(n,nmy,1)
      face_cv(n,nmz,2) = face_cv(n,nmz,1)
      face_cv(n,nee,2) = face_cv(n,nee,1)
      !
    else if (mach_int <= -one) then
      !
      ! This is a supersonic inflow boundary. Fix all
      ! variables to the specified boundary conditions.
      !
      face_cv(n,nec,2) = cv_bc(nec)
      face_cv(n,nmx,2) = cv_bc(nmx)
      face_cv(n,nmy,2) = cv_bc(nmy)
      face_cv(n,nmz,2) = cv_bc(nmz)
      face_cv(n,nee,2) = cv_bc(nee)
      !
    else if (mach_int >= zero) then
      !
      ! This is a subsonic outflow boundary. Set the pressure
      ! to the boundary condition pressure, and extrapolate all
      ! other variables from the interior to the boundary.
      !
      pres = (one-rlx)*pres_int + rlx*pres_bc
      rho = rho_int*pres/pres_int
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*vx
      face_cv(n,nmy,2) = rho*vy
      face_cv(n,nmz,2) = rho*vz
      face_cv(n,nee,2) = pres*one_over_gm1 + half*rho*q2
      !
    else
      !
      ! This is a subsonic inflow boundary. Fix total pressure, total
      ! density (and thus total temperature), and angle of attach, and
      ! use the interior normal Mach number.
      !
      total_cond = one + half*gm1*mach_int*mach_int
      !
      ! Compute the density and pressure using the inflow total density,
      ! total pressure, and interior normal Mach number
      !
      rho  = rtot_bc * total_cond**one_over_1mg
      pres = ptot_bc * total_cond**gam_over_1mg
      !
      vel = abs(mach_int) * sqrt( max( gam*pres/rho , rndoff ) )
      !
      ! Set the boundary velocity in the same direction
      ! as the specified boundary condition velocity
      !
      vx = vel*dir_cosx
      vy = vel*dir_cosy
      vz = vel*dir_cosz
      !
      ! Set the boundary velocity normal to the boundary face
      !
     !vx = vel*unx
     !vy = vel*uny
     !vz = vel*unz
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*vx
      face_cv(n,nmy,2) = rho*vy
      face_cv(n,nmz,2) = rho*vz
      face_cv(n,nee,2) = pres*one_over_gm1 + half*rho*vel*vel
      !
    end if
    !
  end do
  !
end subroutine get_bc_gen_freeflow_sol_cpg_3d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine get_bc_gen_freeflow_sol_cpg_2d(n1,n2,facenrm,face_cv,face_gp, &
                                          cv_bc,pv_bc,gp_bc,rlx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Generic freeflow boundary condition
  ! ##    - At each flux point, use the normal velocity
  ! ##      Mach number to determine if it is an inflow
  ! ##      or outflow boundary condition and whether
  ! ##      it is subsonic or supersonic.
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conserved variables at the flux points
  ! FACE_GP : Gas property variables at the flux points
  ! CV_BC   : Conserved    variables specified for this BC
  ! PV_BC   : Primitive    variables specified for this BC
  ! GP_BC   : Gas property variables specified for this BC
  ! RLX     : Current relaxation value for this BC
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: cv_bc(:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp),    intent(in) :: gp_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp),    intent(in) :: rlx
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,rny,unx,uny
  real(wp) :: rho,vel,pres,total_cond
  real(wp) :: pres_bc,rtot_bc,ptot_bc
  real(wp) :: vmag,dir_cosx,dir_cosy
  real(wp) :: rho_int,rho_inv,vx,vy,q2
  real(wp) :: vn_int,pres_int,aspd_int,mach_int
  real(wp) :: gm1,one_over_gm1,one_over_1mg,gam_over_1mg
  !
continue
  !
  gm1 = gam - one
  one_over_gm1 = one / gm1
  one_over_1mg = -one_over_gm1
  gam_over_1mg = gam * one_over_1mg
  !
  ! Get the boundary static pressure
  !
  pres_bc = gp_bc(npres)
  !
  ! Get the boundary total density and total pressure
  ! NOTE: These total conditions have been saved into static
  !       density and static pressure locations of pv_bc.
  !
  rtot_bc = pv_bc(nec)
  ptot_bc = pv_bc(nee)
  !
  ! Get the directional cosines of the specified boundary condition velocity
  !
  vmag = sqrt( pv_bc(nmx)*pv_bc(nmx) + &
               pv_bc(nmy)*pv_bc(nmy) )
  !
  dir_cosx = pv_bc(nmx) / max(vmag,rndoff)
  dir_cosy = pv_bc(nmy) / max(vmag,rndoff)
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp,cv_bc)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds = sqrt( rnx*rnx + rny*rny )
    unx = rnx/ds
    uny = rny/ds
    !
    rho_int = face_cv(n,nec,1)
    rho_inv = one/rho_int
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    !
    vn_int = vx*unx + vy*uny
    q2 = vx*vx + vy*vy
    !
    pres_int = face_gp(n,npres,1)
    !
    aspd_int = sqrt(max(gam*pres_int*rho_inv,rndoff))
    mach_int = vn_int/aspd_int
    !
    if (mach_int >= one) then
      !
      ! This is a supersonic outflow boundary.
      ! Extrapolate all variables to the boundary.
      !
      face_cv(n,nec,2) = face_cv(n,nec,1)
      face_cv(n,nmx,2) = face_cv(n,nmx,1)
      face_cv(n,nmy,2) = face_cv(n,nmy,1)
      face_cv(n,nee,2) = face_cv(n,nee,1)
      !
    else if (mach_int <= -one) then
      !
      ! This is a supersonic inflow boundary. Fix all
      ! variables to the specified boundary conditions.
      !
      face_cv(n,nec,2) = cv_bc(nec)
      face_cv(n,nmx,2) = cv_bc(nmx)
      face_cv(n,nmy,2) = cv_bc(nmy)
      face_cv(n,nee,2) = cv_bc(nee)
      !
    else if (mach_int >= zero) then
      !
      ! This is a subsonic outflow boundary. Set the pressure
      ! to the boundary condition pressure, and extrapolate all
      ! other variables from the interior to the boundary.
      !
      pres = (one-rlx)*pres_int + rlx*pres_bc
      rho = rho_int*pres/pres_int
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*vx
      face_cv(n,nmy,2) = rho*vy
      face_cv(n,nee,2) = pres*one_over_gm1 + half*rho*q2
      !
    else
      !
      ! This is a subsonic inflow boundary. Fix total pressure, total
      ! density (and thus total temperature), and angle of attach, and
      ! use the interior normal Mach number.
      !
      total_cond = one + half*gm1*mach_int*mach_int
      !
      ! Compute the density and pressure using the inflow total density,
      ! total pressure, and interior normal Mach number
      !
      rho  = rtot_bc * total_cond**one_over_1mg
      pres = ptot_bc * total_cond**gam_over_1mg
      !
      vel = abs(mach_int) * sqrt( max( gam*pres/rho , rndoff ) )
      !
      ! Set the boundary velocity in the same direction
      ! as the specified boundary condition velocity
      !
      vx = vel*dir_cosx
      vy = vel*dir_cosy
      !
      ! Set the boundary velocity normal to the boundary face
      !
     !vx = vel*unx
     !vy = vel*uny
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*vx
      face_cv(n,nmy,2) = rho*vy
      face_cv(n,nee,2) = pres*one_over_gm1 + half*rho*vel*vel
      !
    end if
    !
  end do
  !
end subroutine get_bc_gen_freeflow_sol_cpg_2d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine get_bc_gen_freeflow_sol_cpg_1d(n1,n2,facenrm,face_cv,face_gp, &
                                          cv_bc,pv_bc,gp_bc,rlx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Generic freeflow boundary condition
  ! ##    - At each flux point, use the normal velocity
  ! ##      Mach number to determine if it is an inflow
  ! ##      or outflow boundary condition and whether
  ! ##      it is subsonic or supersonic.
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conserved variables at the flux points
  ! FACE_GP : Gas property variables at the flux points
  ! CV_BC   : Conserved    variables specified for this BC
  ! PV_BC   : Primitive    variables specified for this BC
  ! GP_BC   : Gas property variables specified for this BC
  ! RLX     : Current relaxation value for this BC
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: cv_bc(:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp),    intent(in) :: gp_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp),    intent(in) :: rlx
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,unx
  real(wp) :: rho,vel,pres,total_cond
  real(wp) :: pres_bc,rtot_bc,ptot_bc
  real(wp) :: vmag,dir_cosx
  real(wp) :: rho_int,rho_inv,vx,q2
  real(wp) :: vn_int,pres_int,aspd_int,mach_int
  real(wp) :: gm1,one_over_gm1,one_over_1mg,gam_over_1mg
  !
continue
  !
  gm1 = gam - one
  one_over_gm1 = one / gm1
  one_over_1mg = -one_over_gm1
  gam_over_1mg = gam * one_over_1mg
  !
  ! Get the boundary static pressure
  !
  pres_bc = gp_bc(npres)
  !
  ! Get the boundary total density and total pressure
  ! NOTE: These total conditions have been saved into static
  !       density and static pressure locations of pv_bc.
  !
  rtot_bc = pv_bc(nec)
  ptot_bc = pv_bc(nee)
  !
  ! Get the directional cosines of the specified boundary condition velocity
  !
  vmag = sqrt( pv_bc(nmx)*pv_bc(nmx) )
  !
  dir_cosx = pv_bc(nmx) / max(vmag,rndoff)
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp,cv_bc)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    ds = sqrt( rnx*rnx )
    unx = rnx/ds
    !
    rho_int = face_cv(n,nec,1)
    rho_inv = one/rho_int
    vx = face_cv(n,nmx,1)*rho_inv
    !
    vn_int = vx*unx
    q2 = vx*vx
    !
    pres_int = face_gp(n,npres,1)
    !
    aspd_int = sqrt(max(gam*pres_int*rho_inv,rndoff))
    mach_int = vn_int/aspd_int
    !
    if (mach_int >= one) then
      !
      ! This is a supersonic outflow boundary.
      ! Extrapolate all variables to the boundary.
      !
      face_cv(n,nec,2) = face_cv(n,nec,1)
      face_cv(n,nmx,2) = face_cv(n,nmx,1)
      face_cv(n,nee,2) = face_cv(n,nee,1)
      !
    else if (mach_int <= -one) then
      !
      ! This is a supersonic inflow boundary. Fix all
      ! variables to the specified boundary conditions.
      !
      face_cv(n,nec,2) = cv_bc(nec)
      face_cv(n,nmx,2) = cv_bc(nmx)
      face_cv(n,nee,2) = cv_bc(nee)
      !
    else if (mach_int >= zero) then
      !
      ! This is a subsonic outflow boundary. Set the pressure
      ! to the boundary condition pressure, and extrapolate all
      ! other variables from the interior to the boundary.
      !
      pres = (one-rlx)*pres_int + rlx*pres_bc
      rho = rho_int*pres/pres_int
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*vx
      face_cv(n,nee,2) = pres*one_over_gm1 + half*rho*q2
      !
    else
      !
      ! This is a subsonic inflow boundary. Fix total pressure, total
      ! density (and thus total temperature), and angle of attach, and
      ! use the interior normal Mach number.
      !
      total_cond = one + half*gm1*mach_int*mach_int
      !
      ! Compute the density and pressure using the inflow total density,
      ! total pressure, and interior normal Mach number
      !
      rho  = rtot_bc * total_cond**one_over_1mg
      pres = ptot_bc * total_cond**gam_over_1mg
      !
      vel = abs(mach_int) * sqrt( max( gam*pres/rho , rndoff ) )
      !
      ! Set the boundary velocity in the same direction
      ! as the specified boundary condition velocity
      !
      vx = vel*dir_cosx
      !
      ! Set the boundary velocity normal to the boundary face
      !
     !vx = vel*unx
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*vx
      face_cv(n,nee,2) = pres*one_over_gm1 + half*rho*vel*vel
      !
    end if
    !
  end do
  !
end subroutine get_bc_gen_freeflow_sol_cpg_1d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine get_bc_gen_freeflow_sol_tpg_3d(n1,n2,facenrm,face_cv,face_gp, &
                                          cv_bc,pv_bc,gp_bc,rlx,newton_failures)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Generic freeflow boundary condition
  ! ##    - At each flux point, use the normal velocity
  ! ##      Mach number to determine if it is an inflow
  ! ##      or outflow boundary condition and whether
  ! ##      it is subsonic or supersonic.
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conserved variables at the flux points
  ! FACE_GP : Gas property variables at the flux points
  ! CV_BC   : Conserved    variables specified for this BC
  ! PV_BC   : Primitive    variables specified for this BC
  ! GP_BC   : Gas property variables specified for this BC
  ! RLX     : Current relaxation value for this BC
  ! NEWTON_FAILURES : Number of Newton iterations that failed to converge
  !
  !.. Use Statements ..
  use ovar, only : tmin,tmax,tref,conv_temp
  !
  use thermo_mod, only : compute_mixture_gamma
  use thermo_mod, only : compute_mixture_rgas
  use thermo_mod, only : compute_mixture_enthalpy
  use thermo_mod, only : compute_mixture_cp
  use thermo_mod, only : compute_cp_over_temp_integral
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: cv_bc(:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp),    intent(in) :: gp_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp), intent(inout) :: face_gp(:,:,:)
  real(wp),    intent(in) :: rlx
  integer,  intent(inout) :: newton_failures
  !
  !.. Local Scalars ..
  integer  :: k,m,n
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: vmag,dir_cosx,dir_cosy,dir_cosz
  real(wp) :: pres_bc,rtot_bc,ptot_bc,ttot_bc,htot_bc,rmix_bc
  real(wp) :: rho_int,rho_inv,vx,vy,vz,vn_int,q2
  real(wp) :: pres_int,temp_int,gamm_int,enth_int,aspd_int,mach_int
  real(wp) :: rho,vel,pres,temp,enth,cp,gamm,htot,rmix,val
  real(wp) :: tmin_nd,tmax_nd,dtmx,func,dfunc,dtemp,dterr
  !
  logical(lk) :: convergence_failed
  !
  !.. Local Arrays ..
  real(wp) :: mass_frac(1:ncs)
  !
  !.. Local Parameters ..
  integer,  parameter :: newton_iter = 100
 !real(wp), parameter :: conv_temp = eps6
  !
continue
  !
  tmin_nd = tmin / tref
  tmax_nd = tmax / tref
  !
  newton_failures = 0
  !
  ! Get the boundary static pressure
  !
  pres_bc = gp_bc(npres)
  !
  ! Get the boundary total density and total pressure
  ! NOTE: These total conditions have been saved into static
  !       density and static pressure locations of pv_bc.
  !
  rtot_bc = pv_bc(nec)
  ptot_bc = pv_bc(nee)
  !
  ! Get the directional cosines of the specified boundary condition velocity
  !
  vmag = sqrt( pv_bc(nmx)*pv_bc(nmx) + &
               pv_bc(nmy)*pv_bc(nmy) + &
               pv_bc(nmz)*pv_bc(nmz) )
  !
  dir_cosx = pv_bc(nmx) / max(vmag,rndoff)
  dir_cosy = pv_bc(nmy) / max(vmag,rndoff)
  dir_cosz = pv_bc(nmz) / max(vmag,rndoff)
  !
  ! Compute the mixture gas constant at the inflow plane
  !
  call compute_mixture_rgas(pv_bc(1:ncs),rmix_bc)
  !
  ! Compute the total temperature (dimensional) at the inflow plane
  !
  ttot_bc = ptot_bc / (rtot_bc*rmix_bc)
  !
  ! Compute the total enthalpy at the inflow plane
  !
  call compute_mixture_enthalpy(pv_bc(1:ncs),ttot_bc,htot_bc)
  !
  dtmx = zero
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp,cv_bc,pv_bc,gp_bc) &
  !$ACC&              PRIVATE(mass_frac) &
  !$ACC&              REDUCTION(+:newton_failures) &
  !$ACC&              REDUCTION(max:dtmx)
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
    rho_int = face_cv(n,nec,1)
    rho_inv = one/rho_int
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    vz = face_cv(n,nmz,1)*rho_inv
    !
    vn_int = vx*unx + vy*uny + vz*unz
    q2 = vx*vx + vy*vy + vz*vz
    !
    pres_int = face_gp(n,npres,1)
    temp_int = face_gp(n,ntemp,1)
    gamm_int = face_gp(n,ngamm,1)
    enth_int = face_gp(n,nenth,1)
    !
    aspd_int = sqrt(max(gamm_int*pres_int*rho_inv,rndoff))
    mach_int = vn_int/aspd_int
    !
    if (mach_int >= one) then
      !
      ! This is a supersonic outflow boundary.
      ! Extrapolate all variables to the boundary.
      !
      !$ACC LOOP SEQ
      do m = 1,size(cv_bc)
        face_cv(n,m,2) = face_cv(n,m,1)
      end do
      !
      !$ACC LOOP SEQ
      do m = 1,size(face_gp,dim=2)
        face_gp(n,m,2) = face_gp(n,m,1)
      end do
      !
    else if (mach_int <= -one) then
      !
      ! This is a supersonic inflow boundary. Fix all
      ! variables to the specified boundary conditions.
      !
      !$ACC LOOP SEQ
      do m = 1,size(cv_bc)
        face_cv(n,m,2) = cv_bc(m)
      end do
      !
      !$ACC LOOP SEQ
      do m = 1,size(face_gp,dim=2)
        face_gp(n,m,2) = gp_bc(m)
      end do
      !
    else if (mach_int >= zero) then
      !
      ! This is a subsonic outflow boundary. Set the pressure
      ! to the boundary condition pressure, and extrapolate all
      ! other variables from the interior to the boundary.
      !
      !$ACC LOOP SEQ
      do m = 1,ncs
        mass_frac(m) = face_cv(n,m,1)*rho_inv
      end do
      !
      pres = (one-rlx)*pres_int + rlx*pres_bc
      rho = rho_int*pres/pres_int
      !
      !$ACC LOOP SEQ
      do m = 1,ncs
        face_cv(n,m,2) = rho*mass_frac(m)
      end do
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*vx
      face_cv(n,nmy,2) = rho*vy
      face_cv(n,nmz,2) = rho*vz
      face_cv(n,nee,2) = rho*enth_int - pres + half*rho*q2
      !
    else
      !
      ! This is a subsonic inflow boundary. Fix total pressure, total
      ! density (and thus total temperature), and angle of attach, and
      ! use the interior normal Mach number.
      !
      ! Compute static temperature using a Newton-Raphson method
      !
      temp = temp_int
      !
      !$ACC LOOP SEQ
      do m = 1,ncs
        mass_frac(m) = pv_bc(m)
      end do
      !
      convergence_failed = true
      !$ACC LOOP SEQ
      newton_iteration: &
      do k = 1,newton_iter
        !
        call compute_mixture_enthalpy(mass_frac,temp,enth)
        !
        call compute_mixture_cp(mass_frac,temp,cp)
        !
        gamm = rmix_bc
        !
        call compute_mixture_gamma(cp,gamm)
        !
        htot = enth + half * gamm*rmix_bc*temp * mach_int*mach_int
        !
        func = htot - htot_bc
        dfunc = cp + half*gamm*rmix_bc * mach_int*mach_int
        !
        dtemp = func/dfunc
       !func = (htot - htot_bc) / (cp + half*gamm*rmix_bc*mach_int**2)
        !
        temp = min( tmax_nd , max(tmin_nd,temp-dtemp) )
        !
        dterr = abs(dtemp)/temp
        !
        if (dterr <= conv_temp) then
          convergence_failed = fals
          exit newton_iteration
        end if
        !
      end do newton_iteration
      !
      if (dterr > dtmx) dtmx = dterr
      !
      if (convergence_failed) then
        newton_failures = newton_failures + 1
      end if
      !
      ! Compute the integral of Cp/T for each species and sum
      !
      call compute_cp_over_temp_integral(mass_frac,val,temp,ttot_bc)
      !
      ! Compute static pressure
      !
      pres = ptot_bc / exp(val/rmix_bc)
      !
      ! Compute thermodynamic properties given pressure, temperature, and
      ! mass fractions
      !
      ! Compute gas constant of the mixture
      ! NOTE: Rmix is stored in face_gp(:,ngamm,2) temporarily
      !
      call compute_mixture_rgas(mass_frac,rmix)
      !
      ! Use the ideal gas law to compute the updated density now that we
      ! have the updated mixture gas constant, temperature, and pressure.
      !
      rho = pres / (rmix*temp)
      !
      ! Compute the mixture static enthalpy
      !
      call compute_mixture_enthalpy(mass_frac,temp,enth)
      !
      ! Compute the mixture specific heat at constant pressure, Cp
      !
      call compute_mixture_cp(mass_frac,temp,cp)
      !
      ! Compute the ratio of the specific heats for the gas mixture
      ! NOTE: Upon entry to compute_mixture_gamma:
      !           face_gp(:,ngamm,2) = Rgas_mix (the mixture gas constant)
      !       After returning from compute_mixture_gamma:
      !           face_gp(:,ngamm,2) = gamma_mix (the mixture gamma)
      !
      gamm = rmix
      call compute_mixture_gamma(cp,gamm)
      !
      ! Compute the velocity at the boundary using a fixed Mach number.
      !
      vel = abs(mach_int) * sqrt(max(gamm*pres/rho,rndoff))
      !
      ! Compute the conserved variables for this boundary point
      !
      !$ACC LOOP SEQ
      do m = 1,ncs
        face_cv(n,m,2) = mass_frac(m)*rho
      end do
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*vel*dir_cosx
      face_cv(n,nmy,2) = rho*vel*dir_cosy
      face_cv(n,nmz,2) = rho*vel*dir_cosz
     !face_cv(n,nee,2) = rho*enth - pres + half*rho*vel*vel
      face_cv(n,nee,2) = rho*htot_bc - rho*rmix_bc*temp
      !
      face_gp(n,npres,2) = pres
      face_gp(n,ntemp,2) = temp
      face_gp(n,ngamm,2) = gamm
      face_gp(n,nenth,2) = enth
      face_gp(n,ncp  ,2) = cp
      !
      !
    end if
    !
  end do
  !
end subroutine get_bc_gen_freeflow_sol_tpg_3d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine get_bc_gen_freeflow_sol_tpg_2d(n1,n2,facenrm,face_cv,face_gp, &
                                          cv_bc,pv_bc,gp_bc,rlx,newton_failures)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Generic freeflow boundary condition
  ! ##    - At each flux point, use the normal velocity
  ! ##      Mach number to determine if it is an inflow
  ! ##      or outflow boundary condition and whether
  ! ##      it is subsonic or supersonic.
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conserved variables at the flux points
  ! FACE_GP : Gas property variables at the flux points
  ! CV_BC   : Conserved    variables specified for this BC
  ! PV_BC   : Primitive    variables specified for this BC
  ! GP_BC   : Gas property variables specified for this BC
  ! RLX     : Current relaxation value for this BC
  ! NEWTON_FAILURES : Number of Newton iterations that failed to converge
  !
  !.. Use Statements ..
  use ovar, only : tmin,tmax,tref,conv_temp
  !
  use thermo_mod, only : compute_mixture_gamma
  use thermo_mod, only : compute_mixture_rgas
  use thermo_mod, only : compute_mixture_enthalpy
  use thermo_mod, only : compute_mixture_cp
  use thermo_mod, only : compute_cp_over_temp_integral
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: cv_bc(:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp),    intent(in) :: gp_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp), intent(inout) :: face_gp(:,:,:)
  real(wp),    intent(in) :: rlx
  integer,  intent(inout) :: newton_failures
  !
  !.. Local Scalars ..
  integer  :: k,m,n
  real(wp) :: ds,rnx,rny,unx,uny
  real(wp) :: vmag,dir_cosx,dir_cosy
  real(wp) :: pres_bc,rtot_bc,ptot_bc,ttot_bc,htot_bc,rmix_bc
  real(wp) :: rho_int,rho_inv,vx,vy,vn_int,q2
  real(wp) :: pres_int,temp_int,gamm_int,enth_int,aspd_int,mach_int
  real(wp) :: rho,vel,pres,temp,enth,cp,gamm,htot,rmix,val
  real(wp) :: tmin_nd,tmax_nd,dtmx,func,dfunc,dtemp,dterr
  !
  logical(lk) :: convergence_failed
  !
  !.. Local Arrays ..
  real(wp) :: mass_frac(1:ncs)
  !
  !.. Local Parameters ..
  integer,  parameter :: newton_iter = 100
 !real(wp), parameter :: conv_temp = eps6
  !
continue
  !
  tmin_nd = tmin / tref
  tmax_nd = tmax / tref
  !
  newton_failures = 0
  !
  ! Get the boundary static pressure
  !
  pres_bc = gp_bc(npres)
  !
  ! Get the boundary total density and total pressure
  ! NOTE: These total conditions have been saved into static
  !       density and static pressure locations of pv_bc.
  !
  rtot_bc = pv_bc(nec)
  ptot_bc = pv_bc(nee)
  !
  ! Get the directional cosines of the specified boundary condition velocity
  !
  vmag = sqrt( pv_bc(nmx)*pv_bc(nmx) + &
               pv_bc(nmy)*pv_bc(nmy) )
  !
  dir_cosx = pv_bc(nmx) / max(vmag,rndoff)
  dir_cosy = pv_bc(nmy) / max(vmag,rndoff)
  !
  ! Compute the mixture gas constant at the inflow plane
  !
  call compute_mixture_rgas(pv_bc(1:ncs),rmix_bc)
  !
  ! Compute the total temperature (dimensional) at the inflow plane
  !
  ttot_bc = ptot_bc / (rtot_bc*rmix_bc)
  !
  ! Compute the total enthalpy at the inflow plane
  !
  call compute_mixture_enthalpy(pv_bc(1:ncs),ttot_bc,htot_bc)
  !
  dtmx = zero
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp,cv_bc,pv_bc,gp_bc) &
  !$ACC&              PRIVATE(mass_frac) &
  !$ACC&              REDUCTION(+:newton_failures) &
  !$ACC&              REDUCTION(max:dtmx)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds = sqrt( rnx*rnx + rny*rny )
    unx = rnx/ds
    uny = rny/ds
    !
    rho_int = face_cv(n,nec,1)
    rho_inv = one/rho_int
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    !
    vn_int = vx*unx + vy*uny
    q2 = vx*vx + vy*vy
    !
    pres_int = face_gp(n,npres,1)
    temp_int = face_gp(n,ntemp,1)
    gamm_int = face_gp(n,ngamm,1)
    enth_int = face_gp(n,nenth,1)
    !
    aspd_int = sqrt(max(gamm_int*pres_int*rho_inv,rndoff))
    mach_int = vn_int/aspd_int
    !
    if (mach_int >= one) then
      !
      ! This is a supersonic outflow boundary.
      ! Extrapolate all variables to the boundary.
      !
      !$ACC LOOP SEQ
      do m = 1,size(cv_bc)
        face_cv(n,m,2) = face_cv(n,m,1)
      end do
      !
      !$ACC LOOP SEQ
      do m = 1,size(face_gp,dim=2)
        face_gp(n,m,2) = face_gp(n,m,1)
      end do
      !
    else if (mach_int <= -one) then
      !
      ! This is a supersonic inflow boundary. Fix all
      ! variables to the specified boundary conditions.
      !
      !$ACC LOOP SEQ
      do m = 1,size(cv_bc)
        face_cv(n,m,2) = cv_bc(m)
      end do
      !
      !$ACC LOOP SEQ
      do m = 1,size(face_gp,dim=2)
        face_gp(n,m,2) = gp_bc(m)
      end do
      !
    else if (mach_int >= zero) then
      !
      ! This is a subsonic outflow boundary. Set the pressure
      ! to the boundary condition pressure, and extrapolate all
      ! other variables from the interior to the boundary.
      !
      !$ACC LOOP SEQ
      do m = 1,ncs
        mass_frac(m) = face_cv(n,m,1)*rho_inv
      end do
      !
      pres = (one-rlx)*pres_int + rlx*pres_bc
      rho = rho_int*pres/pres_int
      !
      !$ACC LOOP SEQ
      do m = 1,ncs
        face_cv(n,m,2) = rho*mass_frac(m)
      end do
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*vx
      face_cv(n,nmy,2) = rho*vy
      face_cv(n,nee,2) = rho*enth_int - pres + half*rho*q2
      !
    else
      !
      ! This is a subsonic inflow boundary. Fix total pressure, total
      ! density (and thus total temperature), and angle of attach, and
      ! use the interior normal Mach number.
      !
      ! Compute static temperature using a Newton-Raphson method
      !
      temp = temp_int
      !
      !$ACC LOOP SEQ
      do m = 1,ncs
        mass_frac(m) = pv_bc(m)
      end do
      !
      convergence_failed = true
      !$ACC LOOP SEQ
      newton_iteration: &
      do k = 1,newton_iter
        !
        call compute_mixture_enthalpy(mass_frac,temp,enth)
        !
        call compute_mixture_cp(mass_frac,temp,cp)
        !
        gamm = rmix_bc
        !
        call compute_mixture_gamma(cp,gamm)
        !
        htot = enth + half * gamm*rmix_bc*temp * mach_int*mach_int
        !
        func = htot - htot_bc
        dfunc = cp + half*gamm*rmix_bc * mach_int*mach_int
        !
        dtemp = func/dfunc
       !func = (htot - htot_bc) / (cp + half*gamm*rmix_bc*mach_int**2)
        !
        temp = min( tmax_nd , max(tmin_nd,temp-dtemp) )
        !
        dterr = abs(dtemp)/temp
        !
        if (dterr <= conv_temp) then
          convergence_failed = fals
          exit newton_iteration
        end if
        !
      end do newton_iteration
      !
      if (dterr > dtmx) dtmx = dterr
      !
      if (convergence_failed) then
        newton_failures = newton_failures + 1
      end if
      !
      ! Compute the integral of Cp/T for each species and sum
      !
      call compute_cp_over_temp_integral(mass_frac,val,temp,ttot_bc)
      !
      ! Compute static pressure
      !
      pres = ptot_bc / exp(val/rmix_bc)
      !
      ! Compute thermodynamic properties given pressure, temperature, and
      ! mass fractions
      !
      ! Compute gas constant of the mixture
      ! NOTE: Rmix is stored in face_gp(:,ngamm,2) temporarily
      !
      call compute_mixture_rgas(mass_frac,rmix)
      !
      ! Use the ideal gas law to compute the updated density now that we
      ! have the updated mixture gas constant, temperature, and pressure.
      !
      rho = pres / (rmix*temp)
      !
      ! Compute the mixture static enthalpy
      !
      call compute_mixture_enthalpy(mass_frac,temp,enth)
      !
      ! Compute the mixture specific heat at constant pressure, Cp
      !
      call compute_mixture_cp(mass_frac,temp,cp)
      !
      ! Compute the ratio of the specific heats for the gas mixture
      ! NOTE: Upon entry to compute_mixture_gamma:
      !           face_gp(:,ngamm,2) = Rgas_mix (the mixture gas constant)
      !       After returning from compute_mixture_gamma:
      !           face_gp(:,ngamm,2) = gamma_mix (the mixture gamma)
      !
      gamm = rmix
      call compute_mixture_gamma(cp,gamm)
      !
      ! Compute the velocity at the boundary using a fixed Mach number.
      !
      vel = abs(mach_int) * sqrt(max(gamm*pres/rho,rndoff))
      !
      ! Compute the conserved variables for this boundary point
      !
      !$ACC LOOP SEQ
      do m = 1,ncs
        face_cv(n,m,2) = mass_frac(m)*rho
      end do
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*vel*dir_cosx
      face_cv(n,nmy,2) = rho*vel*dir_cosy
     !face_cv(n,nee,2) = rho*enth - pres + half*rho*vel*vel
      face_cv(n,nee,2) = rho*htot_bc - rho*rmix_bc*temp
      !
      face_gp(n,npres,2) = pres
      face_gp(n,ntemp,2) = temp
      face_gp(n,ngamm,2) = gamm
      face_gp(n,nenth,2) = enth
      face_gp(n,ncp  ,2) = cp
      !
      !
    end if
    !
  end do
  !
end subroutine get_bc_gen_freeflow_sol_tpg_2d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine get_bc_gen_freeflow_sol_tpg_1d(n1,n2,facenrm,face_cv,face_gp, &
                                          cv_bc,pv_bc,gp_bc,rlx,newton_failures)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Generic freeflow boundary condition
  ! ##    - At each flux point, use the normal velocity
  ! ##      Mach number to determine if it is an inflow
  ! ##      or outflow boundary condition and whether
  ! ##      it is subsonic or supersonic.
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conserved variables at the flux points
  ! FACE_GP : Gas property variables at the flux points
  ! CV_BC   : Conserved    variables specified for this BC
  ! PV_BC   : Primitive    variables specified for this BC
  ! GP_BC   : Gas property variables specified for this BC
  ! RLX     : Current relaxation value for this BC
  ! NEWTON_FAILURES : Number of Newton iterations that failed to converge
  !
  !.. Use Statements ..
  use ovar, only : tmin,tmax,tref,conv_temp
  !
  use thermo_mod, only : compute_mixture_gamma
  use thermo_mod, only : compute_mixture_rgas
  use thermo_mod, only : compute_mixture_enthalpy
  use thermo_mod, only : compute_mixture_cp
  use thermo_mod, only : compute_cp_over_temp_integral
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: cv_bc(:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp),    intent(in) :: gp_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp), intent(inout) :: face_gp(:,:,:)
  real(wp),    intent(in) :: rlx
  integer,  intent(inout) :: newton_failures
  !
  !.. Local Scalars ..
  integer  :: k,m,n
  real(wp) :: ds,rnx,unx
  real(wp) :: vmag,dir_cosx
  real(wp) :: pres_bc,rtot_bc,ptot_bc,ttot_bc,htot_bc,rmix_bc
  real(wp) :: rho_int,rho_inv,vx,vn_int,q2
  real(wp) :: pres_int,temp_int,gamm_int,enth_int,aspd_int,mach_int
  real(wp) :: rho,vel,pres,temp,enth,cp,gamm,htot,rmix,val
  real(wp) :: tmin_nd,tmax_nd,dtmx,func,dfunc,dtemp,dterr
  !
  logical(lk) :: convergence_failed
  !
  !.. Local Arrays ..
  real(wp) :: mass_frac(1:ncs)
  !
  !.. Local Parameters ..
  integer,  parameter :: newton_iter = 100
 !real(wp), parameter :: conv_temp = eps6
  !
continue
  !
  tmin_nd = tmin / tref
  tmax_nd = tmax / tref
  !
  newton_failures = 0
  !
  ! Get the boundary static pressure
  !
  pres_bc = gp_bc(npres)
  !
  ! Get the boundary total density and total pressure
  ! NOTE: These total conditions have been saved into static
  !       density and static pressure locations of pv_bc.
  !
  rtot_bc = pv_bc(nec)
  ptot_bc = pv_bc(nee)
  !
  ! Get the directional cosines of the specified boundary condition velocity
  !
  vmag = sqrt( pv_bc(nmx)*pv_bc(nmx) )
  !
  dir_cosx = pv_bc(nmx) / max(vmag,rndoff)
  !
  ! Compute the mixture gas constant at the inflow plane
  !
  call compute_mixture_rgas(pv_bc(1:ncs),rmix_bc)
  !
  ! Compute the total temperature (dimensional) at the inflow plane
  !
  ttot_bc = ptot_bc / (rtot_bc*rmix_bc)
  !
  ! Compute the total enthalpy at the inflow plane
  !
  call compute_mixture_enthalpy(pv_bc(1:ncs),ttot_bc,htot_bc)
  !
  dtmx = zero
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp,cv_bc,pv_bc,gp_bc) &
  !$ACC&              PRIVATE(mass_frac) &
  !$ACC&              REDUCTION(+:newton_failures) &
  !$ACC&              REDUCTION(max:dtmx)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    ds = sqrt( rnx*rnx )
    unx = rnx/ds
    !
    rho_int = face_cv(n,nec,1)
    rho_inv = one/rho_int
    vx = face_cv(n,nmx,1)*rho_inv
    !
    vn_int = vx*unx
    q2 = vx*vx
    !
    pres_int = face_gp(n,npres,1)
    temp_int = face_gp(n,ntemp,1)
    gamm_int = face_gp(n,ngamm,1)
    enth_int = face_gp(n,nenth,1)
    !
    aspd_int = sqrt(max(gamm_int*pres_int*rho_inv,rndoff))
    mach_int = vn_int/aspd_int
    !
    if (mach_int >= one) then
      !
      ! This is a supersonic outflow boundary.
      ! Extrapolate all variables to the boundary.
      !
      !$ACC LOOP SEQ
      do m = 1,size(cv_bc)
        face_cv(n,m,2) = face_cv(n,m,1)
      end do
      !
      !$ACC LOOP SEQ
      do m = 1,size(face_gp,dim=2)
        face_gp(n,m,2) = face_gp(n,m,1)
      end do
      !
    else if (mach_int <= -one) then
      !
      ! This is a supersonic inflow boundary. Fix all
      ! variables to the specified boundary conditions.
      !
      !$ACC LOOP SEQ
      do m = 1,size(cv_bc)
        face_cv(n,m,2) = cv_bc(m)
      end do
      !
      !$ACC LOOP SEQ
      do m = 1,size(face_gp,dim=2)
        face_gp(n,m,2) = gp_bc(m)
      end do
      !
    else if (mach_int >= zero) then
      !
      ! This is a subsonic outflow boundary. Set the pressure
      ! to the boundary condition pressure, and extrapolate all
      ! other variables from the interior to the boundary.
      !
      !$ACC LOOP SEQ
      do m = 1,ncs
        mass_frac(m) = face_cv(n,m,1)*rho_inv
      end do
      !
      pres = (one-rlx)*pres_int + rlx*pres_bc
      rho = rho_int*pres/pres_int
      !
      !$ACC LOOP SEQ
      do m = 1,ncs
        face_cv(n,m,2) = rho*mass_frac(m)
      end do
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*vx
      face_cv(n,nee,2) = rho*enth_int - pres + half*rho*q2
      !
    else
      !
      ! This is a subsonic inflow boundary. Fix total pressure, total
      ! density (and thus total temperature), and angle of attach, and
      ! use the interior normal Mach number.
      !
      ! Compute static temperature using a Newton-Raphson method
      !
      temp = temp_int
      !
      !$ACC LOOP SEQ
      do m = 1,ncs
        mass_frac(m) = pv_bc(m)
      end do
      !
      convergence_failed = true
      !$ACC LOOP SEQ
      newton_iteration: &
      do k = 1,newton_iter
        !
        call compute_mixture_enthalpy(mass_frac,temp,enth)
        !
        call compute_mixture_cp(mass_frac,temp,cp)
        !
        gamm = rmix_bc
        !
        call compute_mixture_gamma(cp,gamm)
        !
        htot = enth + half * gamm*rmix_bc*temp * mach_int*mach_int
        !
        func = htot - htot_bc
        dfunc = cp + half*gamm*rmix_bc * mach_int*mach_int
        !
        dtemp = func/dfunc
       !func = (htot - htot_bc) / (cp + half*gamm*rmix_bc*mach_int**2)
        !
        temp = min( tmax_nd , max(tmin_nd,temp-dtemp) )
        !
        dterr = abs(dtemp)/temp
        !
        if (dterr <= conv_temp) then
          convergence_failed = fals
          exit newton_iteration
        end if
        !
      end do newton_iteration
      !
      if (dterr > dtmx) dtmx = dterr
      !
      if (convergence_failed) then
        newton_failures = newton_failures + 1
      end if
      !
      ! Compute the integral of Cp/T for each species and sum
      !
      call compute_cp_over_temp_integral(mass_frac,val,temp,ttot_bc)
      !
      ! Compute static pressure
      !
      pres = ptot_bc / exp(val/rmix_bc)
      !
      ! Compute thermodynamic properties given pressure, temperature, and
      ! mass fractions
      !
      ! Compute gas constant of the mixture
      ! NOTE: Rmix is stored in face_gp(:,ngamm,2) temporarily
      !
      call compute_mixture_rgas(mass_frac,rmix)
      !
      ! Use the ideal gas law to compute the updated density now that we
      ! have the updated mixture gas constant, temperature, and pressure.
      !
      rho = pres / (rmix*temp)
      !
      ! Compute the mixture static enthalpy
      !
      call compute_mixture_enthalpy(mass_frac,temp,enth)
      !
      ! Compute the mixture specific heat at constant pressure, Cp
      !
      call compute_mixture_cp(mass_frac,temp,cp)
      !
      ! Compute the ratio of the specific heats for the gas mixture
      ! NOTE: Upon entry to compute_mixture_gamma:
      !           face_gp(:,ngamm,2) = Rgas_mix (the mixture gas constant)
      !       After returning from compute_mixture_gamma:
      !           face_gp(:,ngamm,2) = gamma_mix (the mixture gamma)
      !
      gamm = rmix
      call compute_mixture_gamma(cp,gamm)
      !
      ! Compute the velocity at the boundary using a fixed Mach number.
      !
      vel = abs(mach_int) * sqrt(max(gamm*pres/rho,rndoff))
      !
      ! Compute the conserved variables for this boundary point
      !
      !$ACC LOOP SEQ
      do m = 1,ncs
        face_cv(n,m,2) = mass_frac(m)*rho
      end do
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*vel*dir_cosx
     !face_cv(n,nee,2) = rho*enth - pres + half*rho*vel*vel
      face_cv(n,nee,2) = rho*htot_bc - rho*rmix_bc*temp
      !
      face_gp(n,npres,2) = pres
      face_gp(n,ntemp,2) = temp
      face_gp(n,ngamm,2) = gamm
      face_gp(n,nenth,2) = enth
      face_gp(n,ncp  ,2) = cp
      !
      !
    end if
    !
  end do
  !
end subroutine get_bc_gen_freeflow_sol_tpg_1d
!
!###############################################################################
!
end module bc_gen_free_mod
