#ifndef DISABLE_PURE
pure &
#endif
function uiuc_scramjet_combustor_init(xyz,bc,ref) result(return_value)
  !
  !.. Use Statements ..
  use eqn_idx_mod, only : nec,nmx,nmy,nmz,nee,nq,ncs
  use eqn_idx_mod, only : npres,ntemp,ngamm,nenth
  use eqn_idx_mod, only : NCD_Vert,NCD_Span
  use eqn_idx_mod, only : NM_Axial,NM_Vert,NM_Span
  use ovar,        only : bc_in_t
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
  integer  :: k,m,np,scramjet_location
  real(wp) :: rho,p,t,h,gam,rgas,q2,a2,m2,p0,t0,h0
  real(wp) :: vmag,theta
  !
  logical(lk) :: inside_injector,inside_cavity
  !
  !.. Local Automatic Arrays ..
  real(wp) :: cv(1:nq)
  !
  !.. Local Parameters ..
  ! Injector Mach to impose non-zero velocity
  real(wp), parameter :: inj_mach = 0.08_wp
  !
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
  ! Check the coordinates of this cell to determine if
  ! it is in the injector, cavity, or somewhere else.
  !
  scramjet_location = get_uiuc_scramjet_location(xyz)
  !
  if (scramjet_location == 1) then
    inside_injector = true
  else if (abs(scramjet_location) == 2) then
    inside_cavity = true
  end if
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
    gam = bc%gp(ngamm)
    h = bc%gp(nenth)
    !
    ! Mixture gas constant for the boundary conditions
    !
    rgas = p/(rho*t)
    !
    a2 = gam*rgas*t
    m2 = inj_mach*inj_mach
    q2 = m2*a2
    vmag = sqrt(q2)
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
     !theta = atan2(z(k),y(k))
      theta = atan2(xyz(NCD_Span,k),xyz(NCD_Vert,k))
      !
      return_value(k,1:ncs) = rho*bc%pv(1:ncs)
      return_value(k,nec) = rho
      !
      ! vaxial = sqrt(1/2)*vmag
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
    call compute_total_temperature(ref%cv,ref%gp,t0)
    !
    t = t0
    !
    ! Update the static enthalpy with the static temperature initial condition
    !
    call compute_mixture_enthalpy(ref%pv(1:ncs),t,h)
    !
    ! Update the mixture static density initial condition using
    ! the total temperature as the static temperature
    !
    rho = p/(rgas*t)
    !
    do k = 1,np
      !
      return_value(k,1:ncs) = rho*ref%pv(1:ncs)
      !
      return_value(k,nec) = rho
      return_value(k,nmx) = zero
      return_value(k,nmy) = zero
      return_value(k,nmz) = zero
      return_value(k,nee) = rho*h - p
     !return_value(k,nee) = rho*h0 - p
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
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_uiuc_scramjet_location(xyz,wall_radius_sensitivity) &
                             result(return_value)
  !
  ! This function will test to see if the coordinates specified by xyz are
  ! located within the injector or cavity regions of the UIUC scramjet
  ! combustor geometry. The returned value from this function will be:
  !
  ! return_value = 
  !     0 : the xyz coordinates are neither inside the injector or cavity
  !     1 : the xyz coordinates are inside the injector
  !     2 : the xyz coordinates are inside the cavity
  !
  ! Optional Arguments:
  !   wall_radius_sensitivity
  !       this is a way to override the default 'fudge' factor to identify cells
  !       that might not be fully inside the injector or cavity, but are enough
  !       of the cell is to warrant it being identified as such. So far, this
  !       is only used to identify cells within the injector to apply the
  !       sponge filter to, since we want to be less generous in that situation
  !       to mark what cells are inside the injector.
  !
  !.. Use Statements ..
  use eqn_idx_mod, only : NCD_Axial,NCD_Vert,NCD_Span
  use geomvar_mod, only : n_global_cell
  use ovar,        only : grid_scaling_factor
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xyz(:,:)
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: wall_radius_sensitivity
  !
  !.. Function Result ..
  integer :: return_value
  !
  !.. Local Scalars ..
  integer  :: i,j,k,np
  real(wp) :: radius_sensitivity
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: x(:),y(:),z(:),r(:)
  !
  !.. Local Parameters ..
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
  real(wp), parameter :: cavity_radius = 28.48926_wp ! [mm]
  ! Radius that would definitely be an inside cell
  real(wp), parameter :: max_inside_radius = 17.25_wp ! [mm]
  real(wp), parameter :: default_radius_sensitivity = 0.030_wp ! [mm]
 !real(wp), parameter :: default_radius_sensitivity = 0.005_wp ! [mm]
  !
  real(wp), parameter :: cavity_sponge_axial  = cavity_axial_min + 0.3_wp ! [mm]
  real(wp), parameter :: cavity_sponge_radius = cavity_radius - 0.045_wp  ! [mm]
  !
continue
  !
  np = size(xyz,dim=2)
  !
  ! Initialze the function result as being neither inside the
  ! injector or cavity. This will only be changed if the following
  ! loops determine these xyz coordinates are inside either the
  ! injector or cavity.
  !
  return_value = 0
  !
  radius_sensitivity = default_radius_sensitivity
  if (present(wall_radius_sensitivity)) then
    radius_sensitivity = wall_radius_sensitivity
  end if
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
        return_value = 1
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
        !
        return_value = 2
        !
        ! This next set of conditionals is a hack for the initial full-length
        ! 90-degree sector grid created with GridPro. This is to add a few cells
        ! cavity cells to the sponge filter. This should only catch one layer of
        ! cells off the cavity outer wall that are adjacent and just downstream
        ! of the cavity 'backstep' wall. This hack is attempting to prevent
        ! these particular cells from causing simulations to crash, which they
        ! have done quite regularly.
        !
        if (n_global_cell == 367344) then
          if (any(r(:) > cavity_sponge_radius)) then
            if (any(x(:) < cavity_sponge_axial)) then
              return_value = -2
            end if
          end if
        end if
        !
        exit one_pass_loop
        !
      end if
      !
    end do cavity_test_loop
    !
  end do one_pass_loop
  !
end function get_uiuc_scramjet_location
