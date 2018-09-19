module bc_gen_in_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nee
  use eqn_idx, only : ncs,npres,ntemp,ncp,ngamm,nenth
  !
  implicit none
  !
  private
  !
  public :: get_bc_generic_inflow_solution
  !
contains
!
!###############################################################################
!
pure &
subroutine get_bc_generic_inflow_solution(n1,n2,facenrm,face_cv,face_gp, &
                                          pv_in,gp_in,pv_fs,gp_fs, &
                                          relax,newton_failures)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Supersonic inflow boundary condition
  ! ##    - Fix all variables to the inflow condition
  ! #####################################################
  ! #####################################################
  !
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
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
  integer :: nr
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
                                                          gp_rlx=gp_bc)
  !
  if (ncs > 0) then
    !
    ! This inflow condition is very similar to get_bc_sub_inflow_sol_tpg_[123]d
    ! except that it checks the interior velocity normal to the inflow plane,
    ! and switches to a fixed inflow boundary if the interior normal velocity is
    ! determined to be supersonic.
    !
    if (nr == 3) then
      call get_bc_gen_inflow_sol_tpg_3d(n1,n2,facenrm,face_cv,face_gp, &
                                        cv_bc,pv_bc,gp_bc,newton_failures)
    else if (nr == 2) then
      call get_bc_gen_inflow_sol_tpg_2d(n1,n2,facenrm,face_cv,face_gp, &
                                        cv_bc,pv_bc,gp_bc,newton_failures)
    else
      call get_bc_gen_inflow_sol_tpg_1d(n1,n2,facenrm,face_cv,face_gp, &
                                        cv_bc,pv_bc,gp_bc,newton_failures)
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
      call get_bc_gen_inflow_sol_cpg_3d(n1,n2,facenrm,face_cv,face_gp, &
                                        cv_bc,pv_bc)
    else if (nr == 2) then
      call get_bc_gen_inflow_sol_cpg_2d(n1,n2,facenrm,face_cv,face_gp, &
                                        cv_bc,pv_bc)
    else
      call get_bc_gen_inflow_sol_cpg_1d(n1,n2,facenrm,face_cv,face_gp, &
                                        cv_bc,pv_bc)
    end if
    !
  end if
  !
end subroutine get_bc_generic_inflow_solution
!
!###############################################################################
!
pure &
subroutine get_bc_gen_inflow_sol_cpg_3d(n1,n2,facenrm,face_cv,face_gp, &
                                        cv_bc,pv_bc)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: cv_bc(:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
 !real(wp) :: aoa,tan_aoa
  real(wp) :: vmag,dir_cosx,dir_cosy,dir_cosz
  real(wp) :: total_cond,rtot_bc,ptot_bc
  real(wp) :: rho,vx,vy,vz,vel,pres
  real(wp) :: rho_inv,aspd_int,mach_int,veln_int
  real(wp) :: gm1,one_over_gm1,one_over_1mg,gam_over_1mg
  !
continue
  !
  gm1 = gam - one
  one_over_gm1 = one / gm1
  one_over_1mg = -one_over_gm1
  gam_over_1mg = gam * one_over_1mg
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
  dir_cosx = pv_bc(nmx) / vmag
  dir_cosy = pv_bc(nmy) / vmag
  dir_cosz = pv_bc(nmz) / vmag
 !!
 !aoa = atan2(pv_bc(nmy),pv_bc(nmx))
 !tan_aoa = tan( aoa )
 !!
 !dir_cosx = one / sqrt(one + tan_aoa**2)
 !dir_cosy = tan_aoa * dir_cosx
 !dir_cosz = zero
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
    ! Get the host cell primitive variables at this flux point
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    vz = face_cv(n,nmz,1)*rho_inv
    !
    pres = face_gp(n,npres,1)
    !
    ! Normal velocity of the interior/host cell at this flux point
    !
    veln_int = vx*unx + vy*uny + vz*unz
    !
    ! Compute the speed of sound at this flux point
    !
    aspd_int = sqrt( max( gam*pres*rho_inv , rndoff ) )
    mach_int = veln_int / aspd_int
    !
    if (mach_int <= -one) then
      !
      ! This is a supersonic inflow so fix the flow conditions
      ! at the boundary face to the specified boundary conditions.
      !
      face_cv(n,nec,2) = cv_bc(nec)
      face_cv(n,nmx,2) = cv_bc(nmx)
      face_cv(n,nmy,2) = cv_bc(nmy)
      face_cv(n,nmz,2) = cv_bc(nmz)
      face_cv(n,nee,2) = cv_bc(nee)
      !
    else
      !
      ! This is a subsonic boundary. Fix total pressure, total density
      ! (and thus total temperature), and angle of attach, and use the
      ! interior normal Mach number.
      !
      total_cond = one + half*gm1*mach_int*mach_int
      !
      ! Compute the density and pressure using the inflow total density,
      ! total pressure, and interior normal Mach number
      !
      rho  = rtot_bc * total_cond**one_over_1mg
      pres = ptot_bc * total_cond**gam_over_1mg
      !
      vel = mach_int * sqrt( max( gam*pres/rho , rndoff ) )
      !
      ! Set the boundary velocity in the same direction
      ! as the specified boundary condition velocity
      !
      vx = abs(vel)*dir_cosx
      vy = abs(vel)*dir_cosy
      vz = abs(vel)*dir_cosz
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
end subroutine get_bc_gen_inflow_sol_cpg_3d
!
!###############################################################################
!
pure &
subroutine get_bc_gen_inflow_sol_cpg_2d(n1,n2,facenrm,face_cv,face_gp, &
                                        cv_bc,pv_bc)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: cv_bc(:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,rny,unx,uny
 !real(wp) :: aoa,tan_aoa
  real(wp) :: vmag,dir_cosx,dir_cosy
  real(wp) :: total_cond,rtot_bc,ptot_bc
  real(wp) :: rho,vx,vy,vel,pres
  real(wp) :: rho_inv,aspd_int,mach_int,veln_int
  real(wp) :: gm1,one_over_gm1,one_over_1mg,gam_over_1mg
  !
continue
  !
  gm1 = gam - one
  one_over_gm1 = one / gm1
  one_over_1mg = -one_over_gm1
  gam_over_1mg = gam * one_over_1mg
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
  dir_cosx = pv_bc(nmx) / vmag
  dir_cosy = pv_bc(nmy) / vmag
 !!
 !aoa = atan2(pv_bc(nmy),pv_bc(nmx))
 !tan_aoa = tan( aoa )
 !!
 !dir_cosx = one / sqrt(one + tan_aoa**2)
 !dir_cosy = tan_aoa * dir_cosx
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
    ! Get the host cell primitive variables at this flux point
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    !
    pres = face_gp(n,npres,1)
    !
    ! Normal velocity of the interior/host cell at this flux point
    !
    veln_int = vx*unx + vy*uny
    !
    ! Compute the speed of sound at this flux point
    !
    aspd_int = sqrt( max( gam*pres*rho_inv , rndoff ) )
    mach_int = veln_int / aspd_int
    !
    if (mach_int <= -one) then
      !
      ! This is a supersonic inflow so fix the flow conditions
      ! at the boundary face to the specified boundary conditions.
      !
      face_cv(n,nec,2) = cv_bc(nec)
      face_cv(n,nmx,2) = cv_bc(nmx)
      face_cv(n,nmy,2) = cv_bc(nmy)
      face_cv(n,nee,2) = cv_bc(nee)
      !
    else
      !
      ! This is a subsonic boundary. Fix total pressure, total density
      ! (and thus total temperature), and angle of attach, and use the
      ! interior normal Mach number.
      !
      total_cond = one + half*gm1*mach_int*mach_int
      !
      ! Compute the density and pressure using the inflow total density,
      ! total pressure, and interior normal Mach number
      !
      rho  = rtot_bc * total_cond**one_over_1mg
      pres = ptot_bc * total_cond**gam_over_1mg
      !
      vel = mach_int * sqrt( max( gam*pres/rho , rndoff ) )
      !
      ! Set the boundary velocity in the same direction
      ! as the specified boundary condition velocity
      !
      vx = abs(vel)*dir_cosx
      vy = abs(vel)*dir_cosy
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
end subroutine get_bc_gen_inflow_sol_cpg_2d
!
!###############################################################################
!
pure &
subroutine get_bc_gen_inflow_sol_cpg_1d(n1,n2,facenrm,face_cv,face_gp, &
                                        cv_bc,pv_bc)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: cv_bc(:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,unx
 !real(wp) :: aoa,tan_aoa
  real(wp) :: vmag,dir_cosx
  real(wp) :: total_cond,rtot_bc,ptot_bc
  real(wp) :: rho,vx,vel,pres
  real(wp) :: rho_inv,aspd_int,mach_int,veln_int
  real(wp) :: gm1,one_over_gm1,one_over_1mg,gam_over_1mg
  !
continue
  !
  gm1 = gam - one
  one_over_gm1 = one / gm1
  one_over_1mg = -one_over_gm1
  gam_over_1mg = gam * one_over_1mg
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
  dir_cosx = pv_bc(nmx) / vmag
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
    ! Get the host cell primitive variables at this flux point
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    !
    pres = face_gp(n,npres,1)
    !
    ! Normal velocity of the interior/host cell at this flux point
    !
    veln_int = vx*unx
    !
    ! Compute the speed of sound at this flux point
    !
    aspd_int = sqrt( max( gam*pres*rho_inv , rndoff ) )
    mach_int = veln_int / aspd_int
    !
    if (mach_int <= -one) then
      !
      ! This is a supersonic inflow so fix the flow conditions
      ! at the boundary face to the specified boundary conditions.
      !
      face_cv(n,nec,2) = cv_bc(nec)
      face_cv(n,nmx,2) = cv_bc(nmx)
      face_cv(n,nee,2) = cv_bc(nee)
      !
    else
      !
      ! This is a subsonic boundary. Fix total pressure, total density
      ! (and thus total temperature), and angle of attach, and use the
      ! interior normal Mach number.
      !
      total_cond = one + half*gm1*mach_int*mach_int
      !
      ! Compute the density and pressure using the inflow total density,
      ! total pressure, and interior normal Mach number
      !
      rho  = rtot_bc * total_cond**one_over_1mg
      pres = ptot_bc * total_cond**gam_over_1mg
      !
      vel = mach_int * sqrt( max( gam*pres/rho , rndoff ) )
      !
      ! Set the boundary velocity in the same direction
      ! as the specified boundary condition velocity
      !
      vx = abs(vel)*dir_cosx
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
end subroutine get_bc_gen_inflow_sol_cpg_1d
!
!###############################################################################
!
pure &
subroutine get_bc_gen_inflow_sol_tpg_3d(n1,n2,facenrm,face_cv,face_gp, &
                                        cv_bc,pv_bc,gp_bc,newton_failures)
  !
  ! Routine to set subsonic inflow boundary conditions on any subset of a
  ! block face
  !
  ! face_cv     : conserved variable array to be set
  ! face_gp    : gas variable array to be set
  ! facenrm  : metrics in the b.c. coordinate system
  ! pv_bc   : array containing primitive variables for this boundary condition
  !           NOTE: For this subroutine, the pv_bc array should contain:
  !                 pv_bc(1:ncs)   = relaxed species mass fractions
  !                 pv_bc(nec)     = relaxed total density
  !                 pv_bc(nmb:nme) = specified inflow velocities
  !                 pv_bc(nee)     = relaxed total pressure
  !
  !.. Use Statements ..
  use ovar, only : tref,tmin,tmax
  !
  use thermo_mod, only : compute_mixture_gamma
  use thermo_mod, only : compute_mixture_rgas
  use thermo_mod, only : compute_mixture_enthalpy
  use thermo_mod, only : compute_mixture_cp
  use thermo_mod, only : compute_cp_over_temp_integral
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv_bc(:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp),    intent(in) :: gp_bc(:)
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp), intent(inout) :: face_gp(:,:,:)
  integer,  intent(inout) :: newton_failures
  !
  !.. Local Scalars ..
  integer  :: k,m,n
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: vmag,dir_cosx,dir_cosy,dir_cosz
  real(wp) :: rho,rho_inv,vx,vy,vz,vn,vel,pres
  real(wp) :: gamm,htot,mach,temp,enth,cp,rmix,val
  real(wp) :: rmix_bc,rtot_bc,ptot_bc,ttot_bc,htot_bc
  real(wp) :: dtemp,dtmx,dterr,func,dfunc,tmin_nd,tmax_nd
  !
  logical(lk) :: convergence_reached
  !
  !.. Local Arrays ..
  real(wp) :: mass_frac(1:ncs)
  !
  !.. Local Parameters ..
  integer,  parameter :: newton_iter = 100
  real(wp), parameter :: conv_temp = eps6
  !
continue
  !
  tmin_nd = tmin / tref
  tmax_nd = tmax / tref
  !
  newton_failures = 0
  !
  ! Get the inflow total pressure and total density
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
  dir_cosx = pv_bc(nmx) / vmag
  dir_cosy = pv_bc(nmy) / vmag
  dir_cosz = pv_bc(nmz) / vmag
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
  ! Stagnation pressure and temperature specified, Mach number extrapolated
  !
  dtmx = zero
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp,cv_bc,pv_bc,gp_bc) &
  !$ACC&              PRIVATE(mass_frac) &
  !$ACC&              REDUCTION(+:newton_failures) &
  !$ACC&              REDUCTION(max:dtmx)
  do n = n1,n2
    !
    ! Face normal at this boundary point
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    rnz = facenrm(3,n)
    ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
    unx = rnx/ds
    uny = rny/ds
    unz = rnz/ds
    !
    !$ACC LOOP SEQ
    do m = 1,ncs
      mass_frac(m) = pv_bc(m)
    end do
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    !
    ! Extract velocity from the interior face point
    !
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    vz = face_cv(n,nmz,1)*rho_inv
    !
    vn = vx*unx + vy*uny + vz*unz
    !
    ! Extract gas properties from the interior face point
    !
    pres = face_gp(n,npres,1)
    temp = face_gp(n,ntemp,1)
    gamm = face_gp(n,ngamm,1)
    !
    mach = vn / sqrt(max(gamm*pres*rho_inv,rndoff))
    !
    if (mach <= -one) then
      !
      ! This is a supersonic inflow so fix the flow conditions
      ! at the boundary face to the specified boundary conditions.
      !
      !$ACC LOOP SEQ
      do m = 1,size(cv_bc)
        face_cv(n,m,2) = cv_bc(m)
      end do
      !
      !$ACC LOOP SEQ
      do m = 1,size(gp_bc)
        face_gp(n,m,2) = gp_bc(m)
      end do
      !
    else
      !
      ! This is a subsonic boundary. Fix total pressure, total density
      ! (and thus total temperature), and angle of attach, and use the
      ! interior normal Mach number.
      !
      ! Compute static temperature using a Newton-Raphson method
      !
      convergence_reached = fals
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
        htot = enth + half * gamm*rmix_bc*temp * mach*mach
        !
        func = htot - htot_bc
        dfunc = cp + half*gamm*rmix_bc * mach*mach
        !
        dtemp = func/dfunc
       !func = (htot - htot_bc) / (cp + half*gamm*rmix_bc*mach**2)
        !
        temp = min( tmax_nd , max(tmin_nd,temp-dtemp) )
        !
        dterr = abs(dtemp)/temp
        !
        if (dterr <= conv_temp) then
          convergence_reached = true
          exit newton_iteration
        end if
        !
      end do newton_iteration
      !
      if (dterr > dtmx) dtmx = dterr
      !
      if (.not. convergence_reached) then
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
      vel = abs(mach) * sqrt(max(gamm*pres/rho,rndoff))
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
    end if
    !
  end do
  !
  ! Finished setting subsonic inflow boundary conditions
  !
end subroutine get_bc_gen_inflow_sol_tpg_3d
!
!###############################################################################
!
pure &
subroutine get_bc_gen_inflow_sol_tpg_2d(n1,n2,facenrm,face_cv,face_gp, &
                                        cv_bc,pv_bc,gp_bc,newton_failures)
  !
  ! Routine to set subsonic inflow boundary conditions on any subset of a
  ! block face
  !
  ! face_cv     : conserved variable array to be set
  ! face_gp    : gas variable array to be set
  ! facenrm  : metrics in the b.c. coordinate system
  ! pv_bc   : array containing primitive variables for this boundary condition
  !           NOTE: For this subroutine, the pv_bc array should contain:
  !                 pv_bc(1:ncs)   = relaxed species mass fractions
  !                 pv_bc(nec)     = relaxed total density
  !                 pv_bc(nmb:nme) = specified inflow velocities
  !                 pv_bc(nee)     = relaxed total pressure
  !
  !.. Use Statements ..
  use ovar, only : tref,tmin,tmax
  !
  use thermo_mod, only : compute_mixture_gamma
  use thermo_mod, only : compute_mixture_rgas
  use thermo_mod, only : compute_mixture_enthalpy
  use thermo_mod, only : compute_mixture_cp
  use thermo_mod, only : compute_cp_over_temp_integral
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv_bc(:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp),    intent(in) :: gp_bc(:)
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp), intent(inout) :: face_gp(:,:,:)
  integer,  intent(inout) :: newton_failures
  !
  !.. Local Scalars ..
  integer  :: k,m,n
  real(wp) :: ds,rnx,rny,unx,uny
  real(wp) :: vmag,dir_cosx,dir_cosy
  real(wp) :: rho,rho_inv,vx,vy,vn,vel,pres
  real(wp) :: gamm,htot,mach,temp,enth,cp,rmix,val
  real(wp) :: rmix_bc,rtot_bc,ptot_bc,ttot_bc,htot_bc
  real(wp) :: dtemp,dtmx,dterr,func,dfunc,tmin_nd,tmax_nd
  !
  logical(lk) :: convergence_reached
  !
  !.. Local Arrays ..
  real(wp) :: mass_frac(1:ncs)
  !
  !.. Local Parameters ..
  integer,  parameter :: newton_iter = 100
  real(wp), parameter :: conv_temp = eps6
  !
continue
  !
  tmin_nd = tmin / tref
  tmax_nd = tmax / tref
  !
  newton_failures = 0
  !
  ! Get the inflow total pressure and total density
  !
  rtot_bc = pv_bc(nec)
  ptot_bc = pv_bc(nee)
  !
  ! Get the directional cosines of the specified boundary condition velocity
  !
  vmag = sqrt( pv_bc(nmx)*pv_bc(nmx) + &
               pv_bc(nmy)*pv_bc(nmy) )
  !
  dir_cosx = pv_bc(nmx) / vmag
  dir_cosy = pv_bc(nmy) / vmag
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
  ! Stagnation pressure and temperature specified, Mach number extrapolated
  !
  dtmx = zero
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp,cv_bc,pv_bc,gp_bc) &
  !$ACC&              PRIVATE(mass_frac) &
  !$ACC&              REDUCTION(+:newton_failures) &
  !$ACC&              REDUCTION(max:dtmx)
  do n = n1,n2
    !
    ! Face normal at this boundary point
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds = sqrt( rnx*rnx + rny*rny )
    unx = rnx/ds
    uny = rny/ds
    !
    !$ACC LOOP SEQ
    do m = 1,ncs
      mass_frac(m) = pv_bc(m)
    end do
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    !
    ! Extract velocity from the interior face point
    !
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    !
    vn = vx*unx + vy*uny
    !
    ! Extract gas properties from the interior face point
    !
    pres = face_gp(n,npres,1)
    temp = face_gp(n,ntemp,1)
    gamm = face_gp(n,ngamm,1)
    !
    mach = vn / sqrt(max(gamm*pres*rho_inv,rndoff))
    !
    if (mach <= -one) then
      !
      ! This is a supersonic inflow so fix the flow conditions
      ! at the boundary face to the specified boundary conditions.
      !
      !$ACC LOOP SEQ
      do m = 1,size(cv_bc)
        face_cv(n,m,2) = cv_bc(m)
      end do
      !
      !$ACC LOOP SEQ
      do m = 1,size(gp_bc)
        face_gp(n,m,2) = gp_bc(m)
      end do
      !
    else
      !
      ! This is a subsonic boundary. Fix total pressure, total density
      ! (and thus total temperature), and angle of attach, and use the
      ! interior normal Mach number.
      !
      ! Compute static temperature using a Newton-Raphson method
      !
      convergence_reached = fals
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
        htot = enth + half * gamm*rmix_bc*temp * mach*mach
        !
        func = htot - htot_bc
        dfunc = cp + half*gamm*rmix_bc * mach*mach
        !
        dtemp = func/dfunc
       !func = (htot - htot_bc) / (cp + half*gamm*rmix_bc*mach**2)
        !
        temp = min( tmax_nd , max(tmin_nd,temp-dtemp) )
        !
        dterr = abs(dtemp)/temp
        !
        if (dterr <= conv_temp) then
          convergence_reached = true
          exit newton_iteration
        end if
        !
      end do newton_iteration
      !
      if (dterr > dtmx) dtmx = dterr
      !
      if (.not. convergence_reached) then
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
      vel = abs(mach) * sqrt(max(gamm*pres/rho,rndoff))
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
    end if
    !
  end do
  !
  ! Finished setting subsonic inflow boundary conditions
  !
end subroutine get_bc_gen_inflow_sol_tpg_2d
!
!###############################################################################
!
pure &
subroutine get_bc_gen_inflow_sol_tpg_1d(n1,n2,facenrm,face_cv,face_gp, &
                                        cv_bc,pv_bc,gp_bc,newton_failures)
  !
  ! Routine to set subsonic inflow boundary conditions on any subset of a
  ! block face
  !
  ! face_cv     : conserved variable array to be set
  ! face_gp    : gas variable array to be set
  ! facenrm  : metrics in the b.c. coordinate system
  ! pv_bc   : array containing primitive variables for this boundary condition
  !           NOTE: For this subroutine, the pv_bc array should contain:
  !                 pv_bc(1:ncs)   = relaxed species mass fractions
  !                 pv_bc(nec)     = relaxed total density
  !                 pv_bc(nmb:nme) = specified inflow velocities
  !                 pv_bc(nee)     = relaxed total pressure
  !
  !.. Use Statements ..
  use ovar, only : tref,tmin,tmax
  !
  use thermo_mod, only : compute_mixture_gamma
  use thermo_mod, only : compute_mixture_rgas
  use thermo_mod, only : compute_mixture_enthalpy
  use thermo_mod, only : compute_mixture_cp
  use thermo_mod, only : compute_cp_over_temp_integral
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv_bc(:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp),    intent(in) :: gp_bc(:)
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp), intent(inout) :: face_gp(:,:,:)
  integer,  intent(inout) :: newton_failures
  !
  !.. Local Scalars ..
  integer  :: k,m,n
  real(wp) :: ds,rnx,unx
  real(wp) :: vmag,dir_cosx
  real(wp) :: rho,rho_inv,vx,vn,vel,pres
  real(wp) :: gamm,htot,mach,temp,enth,cp,rmix,val
  real(wp) :: rmix_bc,rtot_bc,ptot_bc,ttot_bc,htot_bc
  real(wp) :: dtemp,dtmx,dterr,func,dfunc,tmin_nd,tmax_nd
  !
  logical(lk) :: convergence_reached
  !
  !.. Local Arrays ..
  real(wp) :: mass_frac(1:ncs)
  !
  !.. Local Parameters ..
  integer,  parameter :: newton_iter = 100
  real(wp), parameter :: conv_temp = eps6
  !
continue
  !
  tmin_nd = tmin / tref
  tmax_nd = tmax / tref
  !
  newton_failures = 0
  !
  ! Get the inflow total pressure and total density
  !
  rtot_bc = pv_bc(nec)
  ptot_bc = pv_bc(nee)
  !
  ! Get the directional cosines of the specified boundary condition velocity
  !
  vmag = sqrt( pv_bc(nmx)*pv_bc(nmx) )
  !
  dir_cosx = pv_bc(nmx) / vmag
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
  ! Stagnation pressure and temperature specified, Mach number extrapolated
  !
  dtmx = zero
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp,cv_bc,pv_bc,gp_bc) &
  !$ACC&              PRIVATE(mass_frac) &
  !$ACC&              REDUCTION(+:newton_failures) &
  !$ACC&              REDUCTION(max:dtmx)
  do n = n1,n2
    !
    ! Face normal at this boundary point
    !
    rnx = facenrm(1,n)
    ds = sqrt( rnx*rnx )
    unx = rnx/ds
    !
    !$ACC LOOP SEQ
    do m = 1,ncs
      mass_frac(m) = pv_bc(m)
    end do
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    !
    ! Extract velocity from the interior face point
    !
    vx = face_cv(n,nmx,1)*rho_inv
    !
    vn = vx*unx
    !
    ! Extract gas properties from the interior face point
    !
    pres = face_gp(n,npres,1)
    temp = face_gp(n,ntemp,1)
    gamm = face_gp(n,ngamm,1)
    !
    mach = vn / sqrt(max(gamm*pres*rho_inv,rndoff))
    !
    if (mach <= -one) then
      !
      ! This is a supersonic inflow so fix the flow conditions
      ! at the boundary face to the specified boundary conditions.
      !
      !$ACC LOOP SEQ
      do m = 1,size(cv_bc)
        face_cv(n,m,2) = cv_bc(m)
      end do
      !
      !$ACC LOOP SEQ
      do m = 1,size(gp_bc)
        face_gp(n,m,2) = gp_bc(m)
      end do
      !
    else
      !
      ! This is a subsonic boundary. Fix total pressure, total density
      ! (and thus total temperature), and angle of attach, and use the
      ! interior normal Mach number.
      !
      ! Compute static temperature using a Newton-Raphson method
      !
      convergence_reached = fals
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
        htot = enth + half * gamm*rmix_bc*temp * mach*mach
        !
        func = htot - htot_bc
        dfunc = cp + half*gamm*rmix_bc * mach*mach
        !
        dtemp = func/dfunc
       !func = (htot - htot_bc) / (cp + half*gamm*rmix_bc*mach**2)
        !
        temp = min( tmax_nd , max(tmin_nd,temp-dtemp) )
        !
        dterr = abs(dtemp)/temp
        !
        if (dterr <= conv_temp) then
          convergence_reached = true
          exit newton_iteration
        end if
        !
      end do newton_iteration
      !
      if (dterr > dtmx) dtmx = dterr
      !
      if (.not. convergence_reached) then
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
      vel = abs(mach) * sqrt(max(gamm*pres/rho,rndoff))
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
    end if
    !
  end do
  !
  ! Finished setting subsonic inflow boundary conditions
  !
end subroutine get_bc_gen_inflow_sol_tpg_1d
!
!###############################################################################
!
end module bc_gen_in_mod
