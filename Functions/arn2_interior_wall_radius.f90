pure function arn2_interior_wall_radius(x,eps_1) result(return_value)
  !
  ! NOTE: This function assumes the incoming xyz is in meters
  !       since all the parameter coordinates are in meters.
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: x
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: eps_1
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: e1
  !
  !.. Local Parameters ..
  real(wp), parameter :: pa(2) = [-0.179705_wp,0.0762_wp]
  real(wp), parameter :: pb(2) = [-0.006350_wp,0.0254_wp]
  !
  real(wp), parameter :: c0 =  0.0256084_wp
  real(wp), parameter :: c1 =  0.0667531_wp
  real(wp), parameter :: c2 =  5.4427200_wp
  real(wp), parameter :: c3 = 19.5023000_wp
  !
continue
  !
  e1 = zero
  if (present(eps_1)) e1 = eps_1
  !
  if (x <= pa(1)) then
    return_value = pa(2)
  else if (x >= pb(1)) then
    return_value = pb(2)
  else
    return_value = c0 + c1*x + c2*x*x + c3*x*x*x + e1
  end if
  !
end function arn2_interior_wall_radius
!
!###############################################################################
!
pure function arn2_interior_wall_slope(x) result(return_value)
  !
  ! NOTE: This function assumes the incoming xyz is in meters
  !       since all the parameter coordinates are in meters.
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: x
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Parameters ..
  real(wp), parameter :: pa(2) = [-0.179705_wp,0.0762_wp]
  real(wp), parameter :: pb(2) = [-0.006350_wp,0.0254_wp]
  !
  real(wp), parameter :: c0 =  0.0256084_wp
  real(wp), parameter :: c1 =  0.0667531_wp
  real(wp), parameter :: c2 =  5.4427200_wp
  real(wp), parameter :: c3 = 19.5023000_wp
  !
continue
  !
  if (x <= pa(1)) then
    return_value = zero
  else if (x >= pb(1)) then
    return_value = zero
  else
    return_value = c1 + two*c2*x + three*c3*x*x
  end if
  !
end function arn2_interior_wall_slope
!
!###############################################################################
!
pure function point_is_in_arn2_nozzle(xyz) result(return_value)
  !
  ! NOTE: This function assumes the incoming xyz is in meters
  !       since all the parameter coordinates are in meters.
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xyz(:)
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: x,r,radius
  !
  !.. Local Parameters ..
  real(wp), parameter :: x_in   = -0.196596_wp
  real(wp), parameter :: y_in   =  0.0762_wp
  real(wp), parameter :: x_exit =  0.0_wp
  !
  real(wp), parameter :: eps_1 = eps3
  !
continue
  !
  return_value = fals
  !
  x = xyz(1)
  !
  if (x >= x_in-eps_1 .and. x < x_exit+eps_1) then
    !
    r = norm2(xyz(2:))
    !
    if (r <= y_in+eps_1) then
      radius = arn2_interior_wall_radius(x,eps_1)
      return_value = ( r <= radius )
    end if
    !
  end if
  !
end function point_is_in_arn2_nozzle
!
!###############################################################################
!
pure function point_is_in_arn2_ext_inflow(xyz) result(return_value)
  !
  ! NOTE: This function assumes the incoming xyz is in meters
  !       since all the parameter coordinates are in meters.
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xyz(:)
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: x
  !
  !.. Local Parameters ..
  real(wp), parameter :: y_ext            =  4.25_wp*inches_to_meters
 !real(wp), parameter :: x_end_ext_inflow = -0.17_wp
  real(wp), parameter :: x_end_ext_inflow = -2.5_wp*inches_to_meters
  !
  real(wp), parameter :: eps_1 = eps3
  !
continue
  !
  return_value = fals
  !
  x = xyz(1)
  !
  if (x < x_end_ext_inflow+eps_1) then
    !
    return_value = (.not. point_is_in_arn2_nozzle(xyz))
   !return_value = (norm2(xyz(2:)) >= y_ext-eps_1)
    !
  end if
  !
end function point_is_in_arn2_ext_inflow
!
!###############################################################################
!
pure function point_is_in_arn2_plume(xyz) result(return_value)
  !
  ! NOTE: This function assumes the incoming xyz is in meters
  !       since all the parameter coordinates are in meters.
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xyz(:)
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: r
  !
  !.. Local Parameters ..
  real(wp), parameter :: r_outer_lip  = 1.03_wp*inches_to_meters
  real(wp), parameter :: r_plume_exit =  7.0_wp*inches_to_meters
  real(wp), parameter :: x_plume_exit = 60.0_wp*inches_to_meters
  !
  real(wp), parameter :: plume_slope = (r_plume_exit-r_outer_lip)/x_plume_exit
  !
continue
  !
  return_value = fals
  !
  if (xyz(1) > zero) then
    !
    r = sqrt(xyz(2)*xyz(2) + xyz(3)*xyz(3))
    !
    return_value = ( r <= (plume_slope*xyz(1) + r_outer_lip) )
    !
  end if
  !
end function point_is_in_arn2_plume
!
!###############################################################################
!
subroutine modify_arn2_farfield_mach()
  !
  ! NOTE: This function assumes the incoming xyz is in meters
  !       since all the parameter coordinates are in meters.
  !
  !.. Use Statements ..
  use eqn_idx_mod, only : nec,nmx,nmy,nmz,nee
  use geomvar_mod, only : n_solpts,xyz,ncell,cell
  use ovar,        only : gam,machref
  use flowvar_mod, only : usp
  !
  !.. Local Scalars ..
  integer  :: n,nc,n1,n2,ierr
  real(wp) :: r,rho,vx,vy,vz,q,pres,aspd,mach,qref,gm1
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: modify_cell(:)
  !
  !.. Local Parameters ..
  real(wp), parameter :: r_outer_lip = 1.03_wp*inches_to_meters
  real(wp), parameter :: r_outer_plume = 5.0_wp*inches_to_meters
  real(wp), parameter :: eps_1 = eps4
  !
continue
  !
  gm1 = gam-one
  !
  allocate ( modify_cell(1:ncell) , source=0 , stat=ierr )
  !
  cell_loop: &
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    !
    do n = n1,n2
      !
      r = sqrt(xyz(2,n)*xyz(2,n) + xyz(3,n)*xyz(3,n))
      !
      if (r >= r_outer_plume) then
        modify_cell(nc) = 2
        cycle cell_loop
      else if (xyz(1,n) <= zero) then
        if (.not. point_is_in_arn2_nozzle(xyz(:,n))) then
          modify_cell(nc) = 1
          cycle cell_loop
        end if
      else 
        if (.not. point_is_in_arn2_plume(xyz(:,n))) then
          if (abs(machref-0.005_wp) < eps12) then
            modify_cell(nc) = 2
          else
            modify_cell(nc) = 1
          end if
          cycle cell_loop
        end if
      end if
      !
    end do
    !
  end do cell_loop
  !
  modify_loop: &
  do nc = 1,ncell
    !
    if (modify_cell(nc) == 0) cycle modify_loop
    !
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    !
    do n = n1,n2
      !
      rho = usp(n,nec)
      vx = usp(n,nmx)/rho
      vy = usp(n,nmy)/rho
      vz = usp(n,nmz)/rho
      q = sqrt(vx*vx + vy*vy + vz*vz)
      pres = gm1*(usp(n,nee)-half*rho*q*q)
      aspd = sqrt(max(eps6,gam*pres/rho))
      mach = q/aspd
      !
      qref = machref*aspd
      !
      if (modify_cell(nc) == 2) then
        usp(n,nmx) = qref*rho
        usp(n,nmy) = zero
        usp(n,nmz) = zero
        usp(n,nee) = pres/gm1 + half*rho*qref*qref
      else if (modify_cell(nc) == 1 .and. mach < machref) then
        usp(n,nmx) = vx/q * qref*rho
        usp(n,nmy) = vy/q * qref*rho
        usp(n,nmz) = vz/q * qref*rho
        usp(n,nee) = pres/gm1 + half*rho*qref*qref
      end if
      !
    end do
    !
  end do modify_loop
  !
 !do n = 1,n_solpts
 !  !
 !  rho = usp(n,nec)
 !  vx = usp(n,nmx)/rho
 !  vy = usp(n,nmy)/rho
 !  vz = usp(n,nmz)/rho
 !  q = sqrt(vx*vx + vy*vy + vz*vz)
 !  pres = gm1*(usp(n,nee)-half*rho*q*q)
 !  aspd = sqrt(max(eps6,gam*pres/rho))
 !  mach = q/aspd
 !  !
 !  modify_point = fals
 !  !
 !  if (mach < machref) then
 !    if (xyz(1,n) > zero) then
 !      r = sqrt(xyz(2,n)*xyz(2,n) + xyz(3,n)*xyz(3,n))
 !      if (r > r_outer_lip+eps_1) then
 !        modify_point = true
 !      end if
 !    else if (.not. point_is_in_arn2_nozzle(xyz(:,n))) then
 !      modify_point = true
 !    end if
 !  end if
 !  !
 !  if (modify_point) then
 !    qref = machref*aspd
 !    usp(n,nmx) = vx/q * qref*rho
 !    usp(n,nmy) = vy/q * qref*rho
 !    usp(n,nmz) = vz/q * qref*rho
 !    usp(n,nee) = pres/gm1 + half*rho*qref*qref
 !  end if
 !  !
 !end do
  !
  if (allocated(modify_cell)) then
    deallocate ( modify_cell , stat=ierr )
  end if
  !
end subroutine modify_arn2_farfield_mach
!
!###############################################################################
!
