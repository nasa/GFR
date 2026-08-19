!###############################################################################
!> Routines for Sod shock tube problem.
!> Ref:
!> * Sod, G.A., J.Comp.Phys., Vol.27, pp.1-31, 1978, see p.14.
!> * John, J., Gas Dynamics, 2006, 3rd Ed, Section 5.4, pp.173-184.
!>
!> Note:
!> *  We use the notation of John (2006).
!> *  We use extrapolation boundaries at xmin,xmax to avoid reflections.
!>    Otherwise, the exact solution is more complicated.
!> *  p4>p1
!>
!>   --------------------------------------------
!>   |                     |                    |
!>   |    p4, r4, u4       |     p1, r1, u1     |   time = 0.0
!>   |                     |                    |
!>   --------------------------------------------
!>   xmin                  xd                   xmax
!>
!>   --------------------------------------------
!>   |     |      |           | Uc       | Us   |
!>   |  4  |<-----|    3    --|->  2   --|--> 1 |   time > 0.0
!>   |     |      |           |          |      |
!>   --------------------------------------------
!>   xmin  expansion          slip     shock    xmax
!>         xh     xt                     xs
!###############################################################################
module sod_shock_tube_mod
  !
  use module_kind_types
  !
  implicit none
  !
  private
  !
  public  :: sod_shock_tube_init_problem
  public  :: sod_shock_tube_init_flow
  public  :: sod_shock_tube_exact_solution
  private :: get_p2p1
  !
  ! Grid Domain
  real(wp), parameter, public :: xmin = -half              ! left  bound [m]
  real(wp), parameter, public :: xmax =  half              ! right bound [m]
  real(wp), parameter, public :: xd   =  half*(xmin+xmax)  ! diaphram    [m]
  !
  ! Flow Conditions
  real(wp) :: p1,r1,u1,a1,m1,t1,s1,g1,rgas1       ! Dim-al primitive variables
  real(wp) :: p2,r2,u2,a2,m2,t2,s2,g2,rgas2       ! Dim-al primitive variables
  real(wp) :: p3,r3,u3,a3,m3,t3,s3,g3,rgas3       ! Dim-al primitive variables
  real(wp) :: p4,r4,u4,a4,m4,t4,s4,g4,rgas4,g4m1  ! Dim-al primitive variables
  real(wp) :: p2p1,r2r1,t2t1,a2a1  ! state 2 / state 1
  real(wp) :: p3p4,r3r4,t3t4,a3a4  ! state 3 / state 4
  real(wp) :: uc                   ! u[m/s] of contact discontinuity
  real(wp) :: us                   ! u[m/s] of shock
  real(wp) :: xh,xt                ! x[m] of expansion head, tail
  real(wp) :: xc                   ! x[m] of contact_discontinuity
  real(wp) :: xs                   ! x[m] of shock
  !
contains
!
!###############################################################################
!> Initialize parameters for the Sod shock tube problem.
!>
!> Sets the correct problem constants.
!> Sets the correct freestream reference conditions.
!> Sets the correct inflow boundary condition inputs for bc_input.
!> Does not initialize the solution field.
!###############################################################################
subroutine sod_shock_tube_init_problem()
  !
  !.. Use Statements ..
  use ovar, only : rgasref,gam
  use ovar, only : aoa_alpharef,aoa_betaref
  use ovar, only : rhoref,pref,tref,uref,aref,machref,reyref
  use thermo_mod, only : rgas_cpg
  !
  !.. Local Scalars ..
  integer  :: stat
  real(wp) :: tol
  real(wp) :: g1m1,g1p1,cp1,cp2
  !
continue
  !
  ! Set initial states.
  !
  ! Left State
  r4    = one             ! [kg/m3]
  p4    = one             ! [Pa=N/m2]
  u4    = zero            ! [m/s]
  g4    = 1.4_wp
  rgas4 = rgas_cpg        ! [m^2/(s^2*K)] = [J/(kg*K)]
  !
  ! Right State
  r1    = r4 / eight      ! [kg/m3]
  p1    = p4 / ten        ! [Pa=N/m2]
  u1    = zero            ! [m/s]
  g1    = g4    
  rgas1 = rgas4           ! [m^2/(s^2*K)] = [J/(kg*K)]
  !
  ! Constants
  rgas3 = rgas4
  rgas2 = rgas1
  g3    = g4
  g2    = g1
  g4m1  = g4 - one
  g1m1  = g1 - one
  g1p1  = g1 + one
  !
  ! Compute the final states.
  ! Done once here rather than within the exact solution routine
  ! which is called for every solution point.
  !
  ! Step 1:  Speed of sound from perfect gas relations.
  a4 = sqrt(g4*p4/r4)     ! [m/s]
  a1 = sqrt(g1*p1/r1)     ! [m/s]
  !
  ! Step 2:  Newton-secant iteration for p2/p1, Eq. 5.39.
  tol = 1.0E-05  ! tolerance for p2/p1 iteration convergence
  call get_p2p1(g1,p1,a1,g4,p4,a4,tol, p2p1,stat)
  p2 = p2p1 * p1
  !
  ! Step 3:  Shock velocity, Eq. 5.32.
  us = a1 * sqrt( (g1p1*p2p1 + g1m1) / (two*g1) )
  !
  ! Step 4:  Flow velocity, Eq. 5.35.
  uc = (a1/g1) * (p2p1-one) * sqrt( (two*g1/g1p1) / (p2p1 + g1m1/g1p1) )
  u2 = uc
  u3 = u2
  !
  ! Step 5:  T2/T1, Eq. 5.37.
  t2t1 = ( (g1m1/g1p1)*p2p1 + one ) / ( (g1m1/g1p1)/p2p1 + one )
  a2a1 = sqrt(t2t1)
  a2   = a2a1 * a1
  !
  ! Step 6:  r2/r1, Eq. 5.34.
  r2r1 = ( (g1p1/g1m1)*p2p1 + one ) / ( g1p1/g1m1 + p2p1 )
  r2r1 = p2p1 / t2t1
  r2   = r2r1 * r1
  !
  ! Step 7:  T3/T4, Eq. 5.38.
  p3p4 = p2p1 / (p4/p1)
  t3t4 = p3p4**(g4m1/g4)
  a3a4 = sqrt(t3t4)
  !
  p3   = p3p4 * p4
  a3   = a3a4 * a4
  !
  ! Step 8:  r3/r4, Eq. 5.30.
  r3r4 = p3p4**(one/g4)
! r3r4 = p3p4 / t3t4    !! = p3p4**(one-g4m1/g4) = p3p4**(one/g4)
  r3   = r3r4 * r4
  !
  ! Mach Numbers
  m1 = u1 / a1
  m2 = u2 / a2
  m3 = u3 / a3
  m4 = u4 / a4
  !
  ! Temperature [K]
  t4 = p4 / (r4*rgas4)
  t3 = p3 / (r3*rgas3)
  t2 = p2 / (r2*rgas2)
  t1 = p1 / (r1*rgas1)

  ! Entropy, assuming Cp=constant (ie, g1=g4).
  cp2 = (g2/(g2-one)) * rgas2
  cp1 = (g1/(g1-one)) * rgas1
  s4 = zero  ! left state is ref entropy state
  s3 = s4    ! isentropic expansion
  s2 = s3 + cp2 * log(t2/t3) - rgas2 * log(p2/p3)
  s1 = s2 + cp1 * log(t1/t2) - rgas1 * log(p1/p2)
! s4 = s4 / rgas4  ! For plotting vs Hirsch.
! s3 = s3 / rgas3
! s2 = s2 / rgas2
! s1 = s1 / rgas1
  !
  !
  ! Override reference conditions
  !
  gam          = g4
  rgasref      = rgas4
  aoa_alpharef = zero
  aoa_betaref  = zero
  rhoref       = r4
  pref         = p4
  tref         = t4
  uref         = u4
  aref         = a4
  machref      = m4
  reyref       = -one
  !
  !
  ! Output flow states.
  !
  if ( true ) then
    write (*,*) '    '
    write (*,*) '#########################################'
    write (*,*) '#####   Testcase:  Sod Shock Tube   #####'
    write (*,*) '#########################################'
    write (*,*) '    '
    write (*,*) '--------------------------------------------'
    write (*,*) '|     |      |           | Uc       | Us   |'
    write (*,*) '|  4  |<-----|    3    --|->  2   --|--> 1 |   tm > 0.0'
    write (*,*) '|     |      |           |          |      |'
    write (*,*) '--------------------------------------------'
    write (*,*) 'xmin  expansion          slip     shock    xmax'
    write (*,*) '      xh     xt                     xs'
    write (*,*) '    '
    write (*,*) 'gamma_1           = ', g1
    write (*,*) 'gamma_4           = ', g4
    write (*,*) 'xmin     location = ', xmin, ' [m]'
    write (*,*) 'diaphram location = ', xd  , ' [m]'
    write (*,*) 'xmax     location = ', xmax, ' [m]'
    write (*,*) '    '
    write (*,*) 'State 1:  Right'
    write (*,*) '    p1 = ', p1, ' [Pa]'
    write (*,*) '    r1 = ', r1, ' [kg/m3]'
    write (*,*) '    u1 = ', u1, ' [m/s]'
    write (*,*) '    a1 = ', a1, ' [m/s]'
    write (*,*) '    T1 = ', t1, ' [K]'
    write (*,*) '    m1 = ', m1
    write (*,*) '    '
    write (*,*) 'Shock speed  = ', us, ' [m/s]'
    write (*,*) '    '
    write (*,*) 'State 2:'
    write (*,*) '    p2 = ', p2, ' [Pa]'
    write (*,*) '    r2 = ', r2, ' [kg/m3]'
    write (*,*) '    u2 = ', u2, ' [m/s]'
    write (*,*) '    a2 = ', a2, ' [m/s]'
    write (*,*) '    T2 = ', t2, ' [K]'
    write (*,*) '    m2 = ', m2
    write (*,*) '    '
    write (*,*) 'Contact discontinuity speed  = ', uc, ' [m/s]'
    write (*,*) '    '
    write (*,*) 'State 3:'
    write (*,*) '    p3 = ', p3, ' [Pa]'
    write (*,*) '    r3 = ', r3, ' [kg/m3]'
    write (*,*) '    u3 = ', u3, ' [m/s]'
    write (*,*) '    a3 = ', a3, ' [m/s]'
    write (*,*) '    T3 = ', t3, ' [K]'
    write (*,*) '    m3 = ', m3
    write (*,*) '    '
    write (*,*) 'State 4:  Left'
    write (*,*) '    p4 = ', p4, ' [Pa]'
    write (*,*) '    r4 = ', r4, ' [kg/m3]'
    write (*,*) '    u4 = ', u4, ' [m/s]'
    write (*,*) '    a4 = ', a4, ' [m/s]'
    write (*,*) '    T4 = ', t4, ' [K]'
    write (*,*) '    m4 = ', m4
    write (*,*) '    '
    write (*,*) 'p2p1 = ', p2p1
    write (*,*) '    '
  end if

end subroutine sod_shock_tube_init_problem
!###############################################################################
!> Solution initialization at a point for the Sod shock tube problem.
!###############################################################################
pure function sod_shock_tube_init_flow(xyz) result(return_value)
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
  return_value = sod_shock_tube_exact_solution(xyz,time=zero,nondim=true)
  !
end function sod_shock_tube_init_flow
!###############################################################################
!> Exact solution at a point for the Sod shock tube problem.
!###############################################################################
pure function sod_shock_tube_exact_solution(xyz,time,nondim) result(return_value)
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr
  use eqn_idx_mod, only : nec,nmx,nmy,nmz,nee,nq,nmb,nme
  use ovar,        only : rhoref,aref
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
  real(wp) :: rhoa,rhoa2
  real(wp) :: rho,pres
  real(wp) :: rx,px,ux
  real(wp) :: x,xh,xt,xc,xs
  real(wp) :: gm1
  !
  !.. Local Arrays ..
  real(wp), dimension(3) :: vel
  !
continue
  !
  ! Constants
  !
  rhoa  = rhoref * aref
  rhoa2 = rhoa   * aref
  !
  ! Location of important flow regions.
  !
  xh = xd + ( u4 - a4 ) * time  ! Loc. of expansion region head.
  xt = xd + ( u3 - a3 ) * time  ! Loc. of expansion region tail.
  xc = xd + uc * time           ! Loc. of contact discontinuity.
  xs = xd + us * time           ! Loc. of shock.
  x = xyz(1)                    ! Loc. of solution point.
  !
  ! Dimensional primitive variables.
  !
  if ( x < xh ) then       ! left of expansion head
    rho      = r4
    pres     = p4
    vel(1:3) = [u4,zero,zero]
    gm1      = g4 - one
  else if ( x < xt ) then  ! left of expansion tail
    ux = u4 + u3 * ( x - xh ) / ( xt - xh )
    rx = r4 * ( one - half * g4m1 * ( ux / a4 ) )**( two      / g4m1 )
    px = p4 * ( one - half * g4m1 * ( ux / a4 ) )**( two * g4 / g4m1 )
    rho      = rx
    pres     = px
    vel(1:3) = [ux,zero,zero]
    gm1      = g4 - one
  else if ( x < xc ) then  ! left of contact discontinuity
    rho      = r3
    pres     = p3
    vel(1:3) = [u3,zero,zero]
    gm1      = g3 - one
  else if ( x < xs ) then  ! left of shock
    rho      = r2
    pres     = p2
    vel(1:3) = [u2,zero,zero]
    gm1      = g2 - one
  else                     ! otherwise, right of shock
    rho      = r1
    pres     = p1
    vel(1:3) = [u1,zero,zero]
    gm1      = g1 - one
  end if
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
end function sod_shock_tube_exact_solution
!###############################################################################
!> Newton-Raphson method for iterating on p2/p1.  Eq.5.39.
!###############################################################################
subroutine get_p2p1(g1,p1,a1,g4,p4,a4,tol, p2p1,stat)
  !
  !.. Formal Arguments ..
  real(wp)  ,intent(in)  :: g1,p1,a1
  real(wp)  ,intent(in)  :: g4,p4,a4
  real(wp)  ,intent(in)  :: tol       !> tolerance for p2/p1
  real(wp)  ,intent(out) :: p2p1
  integer   ,intent(out) :: stat      !> convergence status, (0=good,1=fail)
  !
  !.. Local Scalars ..
  real :: c1,c2,c3
  real :: p,dp
  real :: tn,td,t, dtndp,dtddp,dtdp
  real :: f, dfdp
  !
  integer :: iter, itermax = 50
  !
continue
  !
  ! Constants
  c1 = ( (g4-one) / (two*g1) ) * (a1/a4)
  c2 = ( (g1+one) / (two*g1) )
  c3 = two*g4 / (g4-one)
  !
  ! Initial guess for p=p2/p1.
  p = (nine/ten) * p4 / p1
  !
  ! Newton iteration loop
  stat = 1  ! iter has not yet converged
  do iter = 1, itermax
    !
    ! tn,td,t are temporary terms dependent on p.
    tn    = c1*(p-one)
    td    = sqrt(one+c2*(p-one))
    t     = one - tn / td
    !
    ! Derivatives wrt p.
    dtndp = c1
    dtddp = -( half / td ) * c2
    dtdp  = -( td * dtndp - tn *dtddp ) / td**2
    !
    ! f = f(p) Eq.5.39, where p = p2/p1.  dfdp is the derivative.
    f     =      (p4/p1) * t**(c3)            - p
    dfdp  = c3 * (p4/p1) * t**(c3-one) * dtdp - one
    !
    ! Newton-Raphson update.
    dp = -f / dfdp
    p  = p + dp
    !
    ! Check for convergece.
    if ( abs(dp) < tol ) then
      stat = 0  ! iter has converged
      exit
    endif
    !
  enddo
  !
  p2p1 = p
  !
end subroutine get_p2p1
!###############################################################################
end module sod_shock_tube_mod
