!==========================================================================
!> Initialization for density pulse convection problem.
!>
!> Baseline state is freestream, but pulse has increased density.
!> Pulse is defined by its center and width.
!> Pulse can be repeated with a gap between pulses.
!==========================================================================
!
#ifndef DISABLE_PURE
pure &
#endif
function pulse_test_init(xyz) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : velref,aref,uref,Lref
  use ovar, only : gam,rhoref,tref,pref
  use ovar, only : rgasref
  use ovar, only : pulse_test_form,pulse_test_form_dim
  use ovar, only : pulse_test_center,pulse_test_width,pulse_test_gap
  use ovar, only : pulse_test_repeat
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xyz
  !
  !.. Function Return Value ..
  real(wp), dimension(1:2+size(xyz)) :: return_value
  !
  !.. Local Scalars ..
  integer  :: nc,mb,me,ne,nx,nr,i,n
  real(wp) :: gm1,rhoa,rhoa2
  real(wp) :: rho,pres
  real(wp) :: velref_mag,pulse
  real(wp) :: alpha,beta,delta,gap,w
  real(wp) :: is_combo,s_off,s_inc
  !
  !.. Local Arrays ..
  real(wp), dimension(3)   :: vel, flow_dir
  real(wp), dimension(3,4) :: c
  real(wp), dimension(4)   :: dist, wt
  !
continue
  !
  nx = size(xyz,dim=1)                ! =(1:3) grid dim:    1d,2d,3d
  nr = min(nx,pulse_test_form_dim,2)  ! =(1:2) pulse form:  1d,2d
  nc = 1
  mb = nc + 1
  me = nc + nx
  ne = me + 1
  !
  ! Some needed constants
  !
  gm1 = gam - one
  rhoa = rhoref * aref
  rhoa2 = rhoa * aref
  !
  w   = pulse_test_width
  gap = pulse_test_gap
  delta = w / (four*ten)
  alpha = two / w
  !
  !
  ! Determine the unit flow direction based on freestream angle.
  ! If grid < 3d, zero out those components not represented.
  !
  velref_mag = sqrt( sum( velref(1:nx) * velref(1:nx) ) )
  flow_dir(1:3)  = zero  ! in case nx<3
  flow_dir(1:nx) = velref(1:nx) / velref_mag
  !
  !
  ! -------------------------
  ! Baseline flow properties.
  ! -------------------------
  !
  pres = pref
  rho  = pres / (rgasref*tref)
  vel(1:3) = flow_dir(1:3) * uref
  !
  !
  ! ----------------------
  ! Pulse flow properties.
  ! ----------------------
  !
  ! Set whether this is a combo pulse or not.
  is_combo = 0 ; if ( pulse_test_form == 0 ) is_combo = 1  ! Is combo.
  !
  ! Combo pulse must offset each sub-pulse.  Else offset is 0.
  s_off = ( w + gap ) * is_combo
  !
  ! Set the increment between pulse centers of the same type.
  s_inc = ( w + gap ) * (1 + 3*is_combo)
  !
  ! Set the center of each initial pulse.
  do i = 1, 4
    c(:,i) = pulse_test_center(:) + flow_dir(:) * s_off*max(i-1,0)
  end do
  !
  !
  pulse = zero
  loop_pulse_test_repeat : &
  do n = 0, pulse_test_repeat  ! number of repeated pulses
    !
    ! Compute distance from point xyz to pulse center.
    !
    if ( pulse_test_form_dim == 1 ) then  ! 1d pulse form
      do i = 1, 4  ! pulse types
        ! 1d Streamwise distance to pulse center.
        dist(i) =  abs( sum( ( xyz(1:nx) - c(1:nx,i) )*flow_dir(1:nx) ) )
      end do
      wt(1:4) = [one,one,four,six]  ! weights for Gauss & ellipse
      beta  = log10(two) / ((six*delta)**2)  ! Gauss "slope"
    else  ! 2d pulse form
      do i = 1, 4  ! pulse types
        ! 2d Radial distance to pulse center.
        dist(i) = sqrt( sum( ( xyz(1:nx) - c(1:nx,i) )**2 ) )
      end do
      wt(1:4) = [one,one,four,six]  ! weights for Gauss & ellipse
      beta  = log10(two) / ((six*delta)**2)  ! Gauss "slope"
      !
      ! Park & Kim (2014) use this 3d form for Gauss & ellipse.
      ! http://dx.doi.org/10.1016/j.compfluid.2013.11.030
      ! Different weights and slope.
!     wt(1:4) = [1.,1.,1.,1.]  ! weights for Gauss & ellipse
!     beta  = log10(2.) / (360.*delta**2)  ! Gauss "slope"
    end if  ! pulse_test_form
    !
    ! Compute pulse amplitude at point xyz. Peak value is generally 1.
    !
    if      ( any(pulse_test_form==[0,1]) .and. dist(1) <= w/two ) then
      ! Gauss
      pulse = ( fun_gauss1d( dist(1), -delta, beta )*wt(1) &
              + fun_gauss1d( dist(1),  delta, beta )*wt(2) &
              + fun_gauss1d( dist(1),   zero, beta )*wt(3) ) / wt(4)
      exit loop_pulse_test_repeat
    else if ( any(pulse_test_form==[0,2]) .and. dist(2) <= w/two ) then
      ! Square
      pulse = one
      exit loop_pulse_test_repeat
    else if ( any(pulse_test_form==[0,3]) .and. dist(3) <= w/two ) then
      ! Triangle
      pulse = max( one - dist(3) / (w/two) , zero )
      exit loop_pulse_test_repeat
    else if ( any(pulse_test_form==[0,4]) .and. dist(4) <= w/two ) then
      ! Ellipse
      pulse = ( fun_ellip1d( dist(4), -delta, alpha )*wt(1) &
              + fun_ellip1d( dist(4),  delta, alpha )*wt(2) &
              + fun_ellip1d( dist(4),   zero, alpha )*wt(3) ) / wt(4)
      exit loop_pulse_test_repeat
    end if
    !
    ! Translate the center of repeated pulse.
    !
    do i = 1, 4  ! pulse types
      c(1:nr,i) = c(1:nr,i) - flow_dir(1:nr)*s_inc
    end do
    !
  end do loop_pulse_test_repeat
  !
  ! Output to debug the pulse profile.
  ! To output from here, must compile with DISABLE_PURE.
  !
!   write(*,911) '[DEBUG] x,dist(1:4),pulse = ',xyz(1),dist(1:4),pulse
! 911 format(a,99e12.4)
  !
  ! Add the pulse to density.  Velocity and pressure the stay fixed.
  !
  rho = rho + pulse
  !
  ! Evaluate the dimensional initial conditions using conserved variables
  !
  return_value(nc)    = rho
  return_value(mb:me) = rho*vel(1:nx)
  return_value(ne)    = pres/gm1 + half*rho*selfdot(vel)
  !
  ! Non-dimensionalize the initial conditions
  !
  return_value(nc)    = return_value(nc)    / rhoref
  return_value(mb:me) = return_value(mb:me) / rhoa
  return_value(ne)    = return_value(ne)    / rhoa2
  !
end function pulse_test_init
!=====================================================================
!> Baseline 1D Gauss function with unit amplitude.
!=====================================================================
!
#ifndef DISABLE_PURE
pure &
#endif
function fun_gauss1d(x,x0,beta) result(return_value)

  real(wp) :: return_value
  real(wp) ,intent(in) :: x, x0, beta

  return_value = exp(-beta*(x-x0)**2)

end function fun_gauss1d
!=====================================================================
!> Baseline 1D half ellipse function with unit amplitude.
!=====================================================================
!
#ifndef DISABLE_PURE
pure &
#endif
function fun_ellip1d(x,x0,alpha) result(return_value)

  real(wp) :: return_value
  real(wp) ,intent(in) :: x, x0, alpha

  return_value = sqrt( max( one - alpha**2 * (x-x0)**2, zero ) )

end function fun_ellip1d
!=====================================================================
