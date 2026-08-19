!###############################################################################
!> Routines for advection problem of Jiang-Shu.
!> Ref:
!> * Jiang, G-S and Shu, C-W,
!>   NASA-CR-198228 / ICASE-95-73, Oct 1995, see p.17.
!> * Jiang, G-S and Shu, C-W,
!>   J.Comp.Phys., Vol.126, pp.202-228, 1996, see p.213.
!> Note:
!> * They solve a scalar advection equation with rho=[0.0:1.0].
!> * We solve the Euler equations with rho=[1.0:2.0] to avoid negative density.
!###############################################################################
module advection_js_mod
  !
  use module_kind_types
  !
  implicit none
  !
  private
  !
  public  :: advection_js_init_problem
  public  :: advection_js_init_flow
  public  :: advection_js_exact_solution
  private :: fun_gauss1d
  private :: fun_ellip1d
  !
  ! Grid Domain
  real(wp), parameter, public :: xmin = -one  ! [m]
  real(wp), parameter, public :: xmax =  one  ! [m]
  !
  ! Flow Conditions
  real(wp), parameter :: default_density   = 1.0_wp  ! [kg/m3]
  real(wp), parameter :: default_velocityx = 1.0_wp  ! [m/s]
  real(wp), parameter :: default_velocityy = zero    ! [m/s]
  real(wp), parameter :: default_velocityz = zero    ! [m/s]
  real(wp), parameter :: default_pressure  = 1.0_wp  ! [Pa]
  !
contains
!
!###############################################################################
!> Initialize parameters for the Jiang-Shu advection problem.
!>
!> Sets the correct freestream reference conditions.
!> Sets the correct inflow boundary condition inputs for bc_input.
!> Does not initialize the solution field.
!###############################################################################
subroutine advection_js_init_problem()
  !
  !.. Use Statements ..
  use ovar, only : rgasref,gam
  use ovar, only : aoa_alpharef,aoa_betaref
  use ovar, only : rhoref,pref,tref,uref,aref,machref,reyref
  use thermo_mod, only : rgas_cpg
  !
continue
  !
  ! Override reference conditions
  !
  if (gam <= zero) gam = 1.4_wp
  if (rgasref <= zero) rgasref = rgas_cpg
  !
  aoa_alpharef = zero
  aoa_betaref  = zero
  rhoref       = default_density
  pref         = default_pressure
  tref         = pref / (rhoref*rgasref)
  uref         = default_velocityx
  aref         = sqrt( max( gam*rgasref*tref , rndoff ) )
  machref      = uref / aref
  reyref       = -one
  !
end subroutine advection_js_init_problem
!###############################################################################
!> Solution initialization at a point for the Jiang-Shu advection problem.
!###############################################################################
pure function advection_js_init_flow(xyz) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xyz  !< Solution point [m].
  !
  !.. Function Return Value ..
  !> Conserved solution variables [nondim] at the solution point.
  real(wp), dimension(1:2+size(xyz)) :: return_value
  !
continue
  !
  return_value = advection_js_exact_solution(xyz,time=zero,nondim=true)
  !
end function advection_js_init_flow
!###############################################################################
!> Exact solution at a point for the Jiang-Shu advection problem.
!###############################################################################
pure &
function advection_js_exact_solution(xyz,time,nondim) result(return_value)
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr
  use eqn_idx_mod, only : nec,nmx,nmy,nmz,nee,nq,nmb,nme
  use ovar,        only : velref,aref,uref
  use ovar,        only : gam,rhoref,tref,pref
  use ovar,        only : rgasref
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: xyz(:)  !< Solution point [m].
  real(wp),    intent(in) :: time    !< Solution time [s].
  logical(lk), intent(in) :: nondim  !< Nondim return values.
  !
  !.. Function Return Value ..
  !> Conserved solution variables at the solution point.
  !> Units depend on the nondim input flag.
  real(wp), dimension(1:nq) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: xloc,xc,halfwidth
  real(wp) :: gm1,rhoa,rhoa2
  real(wp) :: rho,pres
  real(wp) :: alpha  ! controls "slope" of semi-ellipse
  real(wp) :: beta   ! controls "halfwidth" of gauss
  real(wp) :: delta  ! x-offset [m] of overset pulses
  real(wp) :: pulse
  !
  !.. Local Arrays ..
  real(wp), dimension(3) :: vel
  !
continue
  !
  ! ---------
  ! Constants
  ! ---------
  !
  delta = 0.005_wp                   ! x-offset [m] of overset pulses
  beta  = log10(two)/(six*delta)**2  ! controls "halfwidth" of gauss
  !
  gm1   = gam - one
  rhoa  = rhoref * aref
  rhoa2 = rhoa   * aref
  !
  !
  ! -------------------------------------
  ! Dimensional baseline flow properties.
  ! -------------------------------------
  !
  rho    = default_density
  vel(1) = default_velocityx
  vel(2) = default_velocityy
  vel(3) = default_velocityz
  pres   = default_pressure
  !
  !
  ! ----------------------------------
  ! Dimensional density pulse profile.
  ! ----------------------------------
  !
  ! Shift the physical coordinate based on convected time.
  ! Keep coordinate within the physical domain.
  !
  xloc = xyz(1)                    ! local copy
  if ( xloc == xmin ) xloc = xmax  ! force periodic bc
  xloc = xloc - uref * time        ! time shift
  do while ( xloc <= xmin )        ! keep within domain
    xloc = xloc + (xmax-xmin)
  end do
  !
  ! Baseline pulse
  pulse = zero
  !
  ! Overlay Gauss pulse
  ! Gauss is asymptotic, so halfwidth is infinite.
  xc    = -0.7_wp  ! center of pulse
  pulse = pulse + (        fun_gauss1d(xloc,xc-delta,beta) &
                  +        fun_gauss1d(xloc,xc+delta,beta) &
                  + four * fun_gauss1d(xloc,xc      ,beta) ) / six
  !
  ! Overlay Square pulse
  xc        = -0.3_wp  ! center of pulse
  halfwidth =  0.1_wp
  if ( abs(xloc-xc) <= halfwidth ) then
    pulse = pulse + one
  end if
  !
  ! Overlay Triangle pulse
  xc        =  0.1_wp  ! center of pulse
  halfwidth =  0.1_wp
  if ( abs(xloc-xc) <= halfwidth ) then
    pulse = pulse + one - abs(xloc-xc)/halfwidth
  end if
  !
  ! Overlay Semi-Ellipse pulse
  xc        = 0.5_wp  ! center of pulse
  halfwidth = 0.1_wp
  alpha     = one / halfwidth  ! controls "halfwidth" of semi-ellipse
  if ( abs(xloc-xc) <= halfwidth ) then
    pulse = pulse + (        fun_ellip1d(xloc,xc-delta,alpha) &
                    +        fun_ellip1d(xloc,xc+delta,alpha) &
                    + four * fun_ellip1d(xloc,xc      ,alpha) ) / six
  end if
  !
  !
  ! ----------------------
  ! Pulse flow properties.
  ! ----------------------
  !
  rho = rho + pulse
  !
  ! Dimensional conserved variables
  !
  return_value(nec)     = rho
  return_value(nmb:nme) = rho*vel(1:nr)
  return_value(nee)     = pres/gm1 + half*rho*selfdot(vel(1:nr))
  !
  ! Non-dimensional conserved variables
  !
  if ( nondim ) then
    return_value(nec)     = return_value(nec)     / rhoref
    return_value(nmb:nme) = return_value(nmb:nme) / rhoa
    return_value(nee)     = return_value(nee)     / rhoa2
  end if
  !
end function advection_js_exact_solution
!###############################################################################
!> Baseline 1D Gauss function with unit amplitude.
!###############################################################################
!
pure &
function fun_gauss1d(x,x0,beta) result(return_value)
  !
  real(wp) :: return_value
  real(wp) ,intent(in) :: x, x0, beta
  !
  return_value = exp(-beta*(x-x0)**2)
  !
end function fun_gauss1d
!###############################################################################
!> Baseline 1D half ellipse function with unit amplitude.
!###############################################################################
!
pure &
function fun_ellip1d(x,x0,alpha) result(return_value)
  !
  real(wp) :: return_value
  real(wp) ,intent(in) :: x, x0, alpha
  !
  return_value = sqrt( max( one - alpha**2 * (x-x0)**2, zero ) )
  !
end function fun_ellip1d
!###############################################################################
end module advection_js_mod
