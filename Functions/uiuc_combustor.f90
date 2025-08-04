#ifndef DISABLE_PURE
pure &
#endif
function uiuc_scramjet_combustor_init(xyz,bc,ref) result(return_value)
  !
  !.. Use Statements ..
  use eqn_idx_mod, only : nec,nmx,nmy,nmz,nee,nq,ncs
  use eqn_idx_mod, only : npres,ntemp,ngamm,nenth
  use eqn_idx_mod, only : NCD_Axial,NCD_Vert,NCD_Span
  use eqn_idx_mod, only : NM_Axial,NM_Vert,NM_Span
  use ovar,        only : bc_in_t
  use ovar,        only : grid_scaling_factor
  use ovar,        only : alpha_aoaref,beta_aoaref
  !
  use thermo_mod,  only : compute_mixture_enthalpy
  use thermo_mod,  only : compute_total_temperature
  use thermo_mod,  only : compute_total_pressure
  !
  !.. Formal Arguments ..
  real(wp),      intent(in) :: xyz(:,:)
  type(bc_in_t), intent(in) :: bc
  type(bc_in_t), intent(in) :: ref
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(xyz,dim=2),1:nq)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,m,np
  real(wp) :: rho,p,t,h,rgas,q2,p0,t0,h0
  real(wp) :: vmag,theta
  !
  logical(lk) :: inside_injector,inside_cavity
  !
  !.. Local Automatic Arrays ..
  real(wp) :: cv(1:nq)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: x(:),y(:),z(:),r(:)
  !
  !.. Local Parameters ..
  real(wp), parameter :: inj_angle = 45.0_wp*degrees2radians
  !
  real(wp), parameter :: inj_inflow_axial_min = 296.8_wp ! [mm]
  !
  real(wp), parameter :: inj_outflow_axial_min = 300.91_wp ! [mm]
  real(wp), parameter :: inj_outflow_axial_max = 301.97_wp ! [mm]
  !
  real(wp), parameter :: cavity_axial_min = 343.4_wp ! [mm]
  real(wp), parameter :: cavity_axial_max = 406.0_wp ! [mm]
  !
  ! Radius of the outer wall 
  real(wp), parameter :: wall_radius = 17.45_wp ! [mm]
  ! Radius that would definitely be an inside cell
  real(wp), parameter :: max_inside_radius = 17.25_wp ! [mm]
  real(wp), parameter :: radius_sensitivity = 0.030_wp ! [mm]
 !real(wp), parameter :: radius_sensitivity = 0.005_wp ! [mm]
  !
  ! Injector Mach to impose non-zero velocity
  real(wp), parameter :: inj_mach = 0.08_wp
  !
continue
  !
  np = size(xyz,dim=2)
  !
  inside_injector = fals
  inside_cavity = fals
  !
  !-----------------------------------------------------------------------------
  !
  ! Save the axial, vertical, and spanwise coordinates into arrays with
  ! names that make sense for calculations below while also removing the
  ! grid scaling. Also precompute the radial coordinate
  !
  x = xyz(NCD_Axial,:) / grid_scaling_factor ! F2003 AUTO-REALLOC
  y = xyz(NCD_Vert ,:) / grid_scaling_factor ! F2003 AUTO-REALLOC
  z = xyz(NCD_Span ,:) / grid_scaling_factor ! F2003 AUTO-REALLOC
  r = init_array(np,zero)                    ! F2003 AUTO-REALLOC
  do k = 1,np
    r(k) = sqrt( y(k)*y(k) + z(k)*z(k) )
  end do
  !
  !-----------------------------------------------------------------------------
  !
  ! Check the coordinates of this cell to determine if
  ! it is in the injector, cavity, or somewhere else.
  !
  one_pass_loop: &
  do i = 1,1
    !
    ! Eliminate all the easy cells that are definitely not
    ! in the injector or cavity.
    !
    if ( all(r(:) < max_inside_radius) ) exit one_pass_loop
    if ( all(x(:) < inj_inflow_axial_min) ) exit one_pass_loop
    if ( all(x(:) > cavity_axial_max) ) exit one_pass_loop
    !
    ! Check to see if this cell is inside the injector
    !
    injector_test_loop: &
    do j = 1,1
      !
      do k = 1,np
        !
        if (x(k) > inj_outflow_axial_max) exit injector_test_loop
        !
        if (x(k) < inj_outflow_axial_min) then
          if (r(k) < wall_radius) exit injector_test_loop
        end if
        !
      end do
      !
      ! If we have made it here, this cell is most likely in the injector if
      ! any of the solution points are larger or very close to the wall radius.
      !
      if ( any(r(:) >= wall_radius-radius_sensitivity) ) then
        inside_injector = true
        exit one_pass_loop
      end if
      !
    end do injector_test_loop
    !
    ! Check to see if this cell is inside the cavity
    !
    cavity_test_loop: &
    do j = 1,1
      !
      if ( all(x(:) < cavity_axial_min) ) exit cavity_test_loop
      !
      ! If we have made it here, this cell is most likely in the cavity if
      ! any of the solution points are larger or very close to the wall radius.
      !
      if ( any(r(:) > wall_radius-radius_sensitivity) ) then
        inside_cavity = true
        exit one_pass_loop
      end if
      !
    end do cavity_test_loop
    !
  end do one_pass_loop
  !
  !-----------------------------------------------------------------------------
  !
  ! Set the initial conditions to use for this cell depending
  ! on whether it is in the injector, cavity, or somewhere else.
  !
  if (inside_injector) then
    !
    ! Special initial conditions for a cell inside the injector
    !
    ! Use temperature and species mass fractions
    ! from injector boundary condition
    !
    rho = bc%cv(nec)
    p = bc%gp(npres)
    t = bc%gp(ntemp)
    h = bc%gp(nenth)
    !
    ! Mixture gas constant for the boundary conditions
    !
    rgas = p/(rho*t)
    !
    ! Compute total pressure from the boundary condition which
    ! will be used as the static pressure for the initial conditions
    ! inside the injector
    !
    call compute_total_pressure(bc%cv,bc%gp,p0,dimensional=fals)
    !
    ! Update the mixture static density initial condition using
    ! the total pressure as the static pressure
    !
    rho = p0/(rgas*t)
    !
    do k = 1,np
      !
      theta = atan2(z(k),y(k))
      !
      return_value(k,1:ncs) = rho*bc%pv(1:ncs)
      return_value(k,nec) = rho
      !
      ! vaxial = sqrt(2)*vmag
      ! vr = -vaxial
      ! vvert = vr*cos(theta)
      ! vspan = vr*sin(theta)
      !
      return_value(k,NM_Axial) =  sqrt(half)*rho*vmag
      return_value(k,NM_Vert ) = -sqrt(half)*rho*vmag*cos(theta)
      return_value(k,NM_Span ) = -sqrt(half)*rho*vmag*sin(theta)
      !
      return_value(k,nee) = rho*h - p0 + half*rho*q2
      !
    end do
    !
  else if (inside_cavity) then
    !
    ! Special initial conditions for a cell inside the cavity
    !
    ! Static thermo variables for the reference conditions
    !
    rho = ref%cv(nec)
    p = ref%gp(npres)
    t = ref%gp(ntemp)
    h = ref%gp(nenth)
    !
    ! Mixture gas constant for the reference conditions
    !
    rgas = p/(rho*t)
    !
    ! Square of the velocity magnitude for the reference conditions
    !
    q2 = ref%pv(nmx)*2 + ref%pv(nmy)*2 + ref%pv(nmz)*2
    !
    ! Compute total enthalpy for the reference conditions
    !
    h0 = h + half*q2
    !
    ! Compute total temperature which will be used for
    ! the static temperature initial condition
    !
    call compute_total_temperature(ref%cv,ref%gp,t0,dimensional=fals)
    !
    ! Update the mixture static density initial condition using
    ! the total temperature as the static temperature
    !
    rho = p/(rgas*t0)
    !
    do k = 1,np
      !
      return_value(k,1:ncs) = rho*ref%pv(1:ncs)
      !
      return_value(k,nec) = rho
      return_value(k,nmx) = zero
      return_value(k,nmy) = zero
      return_value(k,nmz) = zero
      return_value(k,nee) = rho*h0 - p
      !
    end do
    !
  else
    !
    ! All other cells use reference conditions as initial conditions
    !
    do m = 1,nq
      do k = 1,np
        return_value(k,m) = ref%cv(m)
      end do
    end do
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
end function uiuc_scramjet_combustor_init
