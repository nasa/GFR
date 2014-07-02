module initialization_mod
  !
  ! Module containing the functions and routines to
  ! initialize the solution
  !
  use module_kind_types
 !use mpivar
  !
  implicit none
  !
  private
  !
 !include "kind_defs.h"
 !include "constants.h"
  !
 !public :: initialize_solution
  public :: initialize_flow
  !
  integer, parameter :: AdvectionEquation = 1
  integer, parameter :: BurgersEquation   = 2
  integer, parameter :: EulerEquations    = 3
  !
  integer, parameter :: SineWave          =  1
  integer, parameter :: CosineWave        =  2
  integer, parameter :: SinePlusOneWave   =  3
  integer, parameter :: CosinePlusOneWave =  4
  integer, parameter :: CosineSquaredWave =  5
  integer, parameter :: SquareWave        =  6
  integer, parameter :: GaussianWave      =  7
  integer, parameter :: EllipseWave       =  8
  integer, parameter :: ComboWave         =  9
  integer, parameter :: ZJsComboWave      = 10
  integer, parameter :: AliasingWave      = 11
  !
  integer, parameter :: SodsProblem     = 1
  integer, parameter :: LaxsProblem     = 2
  integer, parameter :: ShuOsherProblem = 3
  !
contains
!
!###############################################################################
!
subroutine initialize_flow()
  !
  ! Routine to initialize the flow
  !
  !.. Use Statements ..
  use geovar
  use metvar
  use ovar
  use flowvar
  use io_mod,      only : read_restart_file
  use generic_mod, only : usp2v,project_to_porder
  !
  !.. Local Scalars ..
  integer  :: n1,n2,nc,np,ns,m,nf,ierr
  real(wp) :: c1,c3,c4,gm1,rt2
  real(wp) :: rho_inf,vx_inf,vy_inf,p_inf,T_inf,a_inf,mach_inf
  real(wp) :: Rgas,bigR,xc,yc,Lx,Ly,cp,xt,yt,rt,du,dv,dT,beta
  real(wp) :: vx,vy,Temp,rho,pres,rho_nd,vx_nd,vy_nd,pres_nd
  real(wp) :: bigGam,func
  real(wp) :: rhoa,rhoa2
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable, dimension(:,:) :: uspstart
  !
continue
  !
  ! Some constants from the free-stream conditions
  !
  gm1 = gam - one
  rhoa = rhoref * aref
  rhoa2 = rhoa * aref
  !
  ! Allocate the solution arrays
  !
  allocate ( usp(1:neqns,1:n_totpts) , source=zero , stat=ierr )
  call alloc_error("initialize_flow","usp",1,ierr)
  allocate ( vsp(1:neqns,1:n_totpts) , source=zero , stat=ierr )
  call alloc_error("initialize_flow","vsp",1,ierr)
  allocate ( residual(1:neqns,1:n_solpts) , source=zero , stat=ierr )
  call alloc_error("initialize_flow","residual",1,ierr)
  allocate ( uoldsp(1:neqns,1:n_solpts) , source=zero , stat=ierr )
  call alloc_error("initialize_flow","uoldsp",1,ierr)
  allocate ( faceusp(1:neqns,1:n_fp1de,1:2,1:nface) , source=zero , stat=ierr )
  call alloc_error("initialize_flow","faceusp",1,ierr)
  allocate ( faceflx(1:neqns,1:n_fp1de,1:2,1:nface) , source=zero , stat=ierr )
  call alloc_error("initialize_flow","faceusp",1,ierr)
  !
  allocate ( uspstart(1:neqns,1:n_solpts) , source=zero , stat=ierr )
  call alloc_error("initialize_flow","uspstart",1,ierr)
  !
  ! Allocate the storage for the Runge-Kutta methods
  !
  ns = 0
  if (iRKtype == Classic_RK) then
    ns = num_stages
  else if (iRKtype == CK4_RK .or. iRKtype == TVD4_RK) then
    ns = 1
    if (iRKtype == TVD4_RK) then
      allocate ( usp_sum(1:neqns,1:n_solpts) , source=zero , stat=ierr )
      call alloc_error("initialize_flow","usp_sum",1,ierr)
    end if
  end if
  !
  allocate ( rssprk(1:neqns,1:n_solpts,1:ns) , source=zero , stat=ierr )
  call alloc_error("initialize_flow","rssprk",1,ierr)
  !
  ! Initialize uspstart depending on the test case requested
  !
  if (any(itestcase == Transport_Problems)) then
    !
    ! Initialize the inflow boundary conditions to the reference conditions
    !
    bc_in(1) = rhoref
    bc_in(2) = vxref
    bc_in(3) = vyref
    bc_in(4) = pref
    !
    ! Nondimensionalize the inflow boundary conditions
    !
    bc_in(1) = bc_in(1) / rhoref
    bc_in(2) = bc_in(2) / aref
    bc_in(3) = bc_in(3) / aref
    bc_in(4) = bc_in(4) / rhoa2
    !
   !!
   !! Vortex parameters
   !!
   !bigGam = 13.5_wp
   !bigR = 1.5_wp
   !!
   !! Switch between vortex transport (2D Euler equations problem)
   !!           and density transport (2D advection equation problem)
   !!
   !c4 = deltafun(itestcase == Vortex_Transport)
   !!
   !! For Vincent 2011 vortex, pressure has a leading term of 1/(gam*mach**2)
   !! This term acts to nondimensionalize pressure by
   !! p_nondim = p_dim/(rhoref*uref**2)
   !! Changing the leading term to 1/gam changes
   !! the pressure nondimensionalization to
   !! p_nondim = p_dim/(rhoref*aref**2)
   !!
   !c3 = half*machref*machref*gm1
   !!
   !do nc = 1,ncell
   !  !
   !  n1 = pts_in_cell_ptr(nc)+1
   !  n2 = pts_in_cell_ptr(nc+1)
   !  !
   !  do np = n1,n2
   !    !
   !    ! Vincent 2011 (very similar to workshop) vortex
   !    !
   !    xt = xyz(1,np)
   !    yt = xyz(2,np)
   !    !
   !    func = (one - xt*xt - yt*yt)/(two*bigR*bigR)
   !    !
   !    c1 = bigGam*exp(func)/(two*pi)
   !    !
   !    ! Compute the dimensional fluxuations du, dv, and dt
   !    !
   !    du = -c1*yt/bigR * uref
   !    dv =  c1*xt/bigR * uref
   !    dT = -c3*c1*c1 * tref
   !    !
   !    ! Compute the dimensional velocities and temperature
   !    !
   !    vx   = du + vxref
   !    vy   = dv + vyref
   !    temp = dt + tref
   !    !
   !    ! Compute the remaining dimensional primitive variables
   !    !
   !    rho = rhoref * (temp/tref)**(one/gm1)
   !    pres = pref  * (temp/tref)**(gam/gm1)
   !    !
   !    ! Compute the dimensional conserved variables
   !    !
   !    uspstart(1,np) = rho
   !    uspstart(2,np) = rho*vx
   !    uspstart(3,np) = rho*vy
   !    uspstart(4,np) = pres/gm1*c4 + half*rho*(vx*vx + vy*vy)
   !    !
   !  end do
   !  !
   !end do
   !!
   !! Nondimensionalize the conserved variables
   !!
   !uspstart(1,1:n_solpts) = uspstart(1,1:n_solpts) / rhoref
   !uspstart(2,1:n_solpts) = uspstart(2,1:n_solpts) / rhoa
   !uspstart(3,1:n_solpts) = uspstart(3,1:n_solpts) / rhoa
   !uspstart(4,1:n_solpts) = uspstart(4,1:n_solpts) / rhoa2
   !!
    !
    if (init_loc == 1) then
      !
      write (*,4) "Gauss"
      do nc = 1,ncell
        n1 = pts_in_cell_ptr(nc)+1
        n2 = pts_in_cell_ptr(nc+1)
        do np = n1,n2
          uspstart(:,np) = vortex2D( xyz_gauss(1,np) , xyz_gauss(2,np) )
        end do
        do m = 1,neqns
          uspstart(m,n1:n2) = matmul( gauss2numer , uspstart(m,n1:n2) )
        end do
      end do
      !
    else if (init_loc == 2) then
      !
      write (*,4) "Lobatto"
      do nc = 1,ncell
        n1 = pts_in_cell_ptr(nc)+1
        n2 = pts_in_cell_ptr(nc+1)
        do np = n1,n2
          uspstart(:,np) = vortex2D( xyz_lobatto(1,np) , xyz_lobatto(2,np) )
        end do
        do m = 1,neqns
          uspstart(m,n1:n2) = matmul( lobatto2numer , uspstart(m,n1:n2) )
        end do
      end do
      !
    else if (init_loc == -1) then
      !
      write (*,4) "right Radau"
      do nc = 1,ncell
        n1 = pts_in_cell_ptr(nc)+1
        n2 = pts_in_cell_ptr(nc+1)
        do np = n1,n2
          uspstart(:,np) = vortex2D( xyz_radau(1,np) , xyz_radau(2,np) )
        end do
        do m = 1,neqns
          uspstart(m,n1:n2) = matmul( radau2numer , uspstart(m,n1:n2) )
        end do
      end do
      !
    else
      !
      do np = 1,n_solpts
        uspstart(:,np) = vortex2D( xyz(1,np) , xyz(2,np) )
      end do
      !
    end if
    !
    ! Lets go through the boundary faces and set the velocity to zero.
    ! This might prevent a very small shear layer at the periodic boundaries
    !
   !do nf = 1,nfbnd
   !  !
   !  if (bface(1,nf) == bc_periodic) then
   !    nc = cells_on_face(1,nf)
   !    n1 = pts_in_cell_ptr(nc)+1
   !    n2 = pts_in_cell_ptr(nc+1)
   !    uspstart(2:neqns-1,n1:n2) = zero
   !  end if
   !  !
   !end do
    !
  else
    !
    write (iout,5)
    call stop_recon(8,"initialize_flow")
    !
  end if
  !
  ! Check for a restart job, otherwise initialize usp to uspstart
  !
  if (irestart == 1) then
    !
    ! Read in usp and the current iteration, time step,
    ! and elapsed time from the restart file
    !
    call read_restart_file
    !
    ! Initialize the iteration counter, time step,
    ! and elapsed time to that in the restart file
    !
    d_t   = dtrst
    itcur = itrst
    time  = tmrst
    !
  else
    !
    ! Project the initial solution to p_order
    !
    call project_to_porder(ncell,pts_in_cell_ptr,uspstart)
    !
    ! Initialize usp to uspstart
    !
    usp(:,1:n_solpts) = uspstart(:,1:n_solpts)
    !
    ! Initialize the iteration counter, time step, and elapsed time
    !
    d_t   = dt_in
    itcur = 0
    time  = zero
    !
  end if
  !
  ! Initialize the backup conserved variables solution arrays
  !
  uoldsp(:,1:n_solpts) = usp(:,1:n_solpts)
  !
  ! Initialize the primitive variables solution array
  !
  vsp(:,1:n_solpts) = usp2v( usp(:,1:n_solpts) )
  !
  ! Output the free-stream conditions
  !
  if (mypnum == 0) then
    if (itestcase /= Density_Transport .and. itestcase /= Entropy_Transport) then
      !
      write (iout,1) rhoref,pref,tref,machref,aref, &
                     vxref,vyref,aoaref,gam,rgasref
      write (iout,2) rhoref/rhoref,pref/rhoa2,tref/tref,aref/aref, &
                     vxref/aref,vyref/aref,rgasref*tref/aref/aref
      !
    else
      !
      write (iout,3) rhoref,pref,tref,machref,uref/machref, &
                     vxref,vyref,aoaref,gam,rgasref
      !
    end if
    !
  end if
  !
  ! Deallocate uspstart
  !
  deallocate ( uspstart , stat=ierr )
  call alloc_error("initialize_flow","uspstart",2,ierr)
  !
  ! Format Statements
  !
  1 format (/," Dimensional Free-stream conditions",/, &
              " Density         = ",es14.6," kg/m^3",/, &
              " Pressure        = ",es14.6," N/m^2",/, &
              " Temperature     = ",es14.6," K",/, &
              " Mach Number     = ",es14.6,/, &
              " Speed of Sound  = ",es14.6," m/s",/, &
              " X-Velocity      = ",es14.6," m/s",/, &
              " Y-Velocity      = ",es14.6," m/s",/, &
              " Angle of Attack = ",es14.6," degrees"/, &
              " Gamma           = ",es14.6,/, &
              " Gas Constant    = ",es14.6," J/(kg*K)",/)
  2 format (/," Non-dimensional Free-stream conditions",/, &
              " Density         = ",es14.6,/, &
              " Pressure        = ",es14.6,/, &
              " Temperature     = ",es14.6,/, &
              " Speed of Sound  = ",es14.6,/, &
              " X-Velocity      = ",es14.6,/, &
              " Y-Velocity      = ",es14.6,/, &
              " Gas Constant    = ",es14.6,/)
  3 format (/," The density vortex transport problem is being",/, &
              " used to model the linear 2D advection equation.",//, &
              " The Non-dimensional Free-stream conditions are",/, &
              " Density         = ",es14.6,/, &
              " Pressure        = ",es14.6,/, &
              " Temperature     = ",es14.6,/, &
              " Mach Number     = ",es14.6,/, &
              " Speed of Sound  = ",es14.6,/, &
              " X-Velocity      = ",es14.6,/, &
              " Y-Velocity      = ",es14.6,/, &
              " Angle of Attack = ",es14.6,/, &
              " Gamma           = ",es14.6,/, &
              " Gas Constant    = ",es14.6,/)
  4 format (/," Computing the initial solution at the ",a," points.",/)
  5 format (/," ***************** ERROR! ****************",/, &
              "  THIS TEST CASE IS CURRENTLY NOT SET UP!",/, &
              " ***************** ERROR! ****************",/)
  !
end subroutine initialize_flow
!
!###############################################################################
!
!subroutine initialize_solution(EquationType,option,nsolpt,nq,x,q)
!
!  integer,  intent(in) :: EquationType
!  integer,  intent(in) :: option
!  integer,  intent(in) :: nsolpt
!  integer, intent(out) :: nq
!
!  real(wp),              dimension(1:nsolpt), intent(in) :: x
!  real(wp), allocatable, dimension(:,:),   intent(inout) :: q
!
!  real(wp) :: xmin,xmax,xlen,xcen,xhw
!
!continue
!
!  xmin = x(1)
!  xmax = x(nsolpt)
!  xlen = xmax-xmin
!  xcen = xmin + half*xlen
!  xhw  = xlen/four
!
!  if (EquationType == AdvectionEquation) then
!    nq = 1
!  else if (EquationType == BurgersEquation) then
!    nq = 1
!  else if (EquationType == EulerEquations) then
!    nq = 3
!  end if
!
!  if (allocated(q)) deallocate(q)
!  allocate ( q(1:nsolpt,1:nq) )
!  q(:,:) = zero
!
!  if (EquationType == AdvectionEquation .or. &
!      EquationType == BurgersEquation ) then
!
!    if (option == SineWave) then
!      q(:,1) =  sine_wave(x)
!    else if (option == CosineWave) then
!      q(:,1) =  cosine_wave(x)
!    else if (option == SinePlusOneWave) then
!      q(:,1) =  sine_wave(x) + one
!    else if (option == CosinePlusOneWave) then
!      q(:,1) =  cosine_wave(x) + one
!    else if (option == CosineSquaredWave) then
!      q(:,1) = -cosine_wave(x)*half + half
!    else if (option == SquareWave) then
!      q(:,1) = square_wave(x,xcen,xhw)
!    else if (option == GaussianWave) then
!      q(:,1) = gaussian_wave(x,xcen,c0_10*c0_75)*half
!    else if (option == EllipseWave) then
!      q(:,1) = ellipse_wave(x,xcen,xhw)
!    else if (option == ComboWave) then
!      q(:,1) = gaussian_wave(x,xmin+c0_15*xlen,c0_03*xlen) + &
!               square_wave  (x,xmin+c0_45*xlen,c0_10*xlen) + &
!               ellipse_wave (x,xmin+c0_80*xlen,c0_10*xlen)
!    else if (option == ZJsComboWave) then
!      q(:,1) = gaussian_wave(x,xmin+c0_15*xlen,c0_03/two*xlen) + &
!               square_wave  (x,xmin+c0_35*xlen,c0_10/two*xlen) + &
!               triangle_wave(x,xmin+c0_55*xlen,c0_10/two*xlen) + &
!               ellipse_wave (x,xmin+c0_75*xlen,c0_10/two*xlen)
!    else if (option == AliasingWave) then
!      q(:,1) = sine_wave(x,four)
!    end if
!
!! else if (EquationType == EulerEquations) then
!
!!   if (option == SodsProblem) then
!!     
!!   else if (option == LaxsProblem) then
!!     
!!   else if (option == ShuOsherProblem) then
!!     
!!   end if
!
!  end if
!
!end subroutine initialize_solution
!!
!!###############################################################################
!!
!pure function sine_wave(x,freq) result(return_value)
!
!  real(wp), optional,     intent(in) :: freq
!  real(wp), dimension(:), intent(in) :: x
!
!  real(wp), dimension(size(x)) :: return_value
!
!  real(wp) :: xmin,xmax,xhw
!
!continue
!
!  xmin = x(lbound(x,1))
!  xmax = x(ubound(x,1))
!  xhw  = half*(xmax-xmin)
!  if (present(freq)) xhw = xhw/freq
!
!  return_value = sin(pi*(x-xmin)/xhw)
!
!end function sine_wave
!!
!!###############################################################################
!!
!pure function cosine_wave(x,freq) result(return_value)
!
!  real(wp), optional,     intent(in) :: freq
!  real(wp), dimension(:), intent(in) :: x
!
!  real(wp), dimension(size(x)) :: return_value
!
!  real(wp) :: xmin,xmax,xhw
!
!continue
!
!  xmin = x(lbound(x,1))
!  xmax = x(ubound(x,1))
!  xhw  = half*(xmax-xmin)
!  if (present(freq)) xhw = xhw/freq
!
!  return_value = cos(pi*(x-xmin)/xhw)
!
!end function cosine_wave
!!
!!###############################################################################
!!
!pure function square_wave(x,xcen,hw) result(return_value)
!
!  real(wp),               intent(in) :: xcen
!  real(wp),               intent(in) :: hw
!  real(wp), dimension(:), intent(in) :: x
!
!  real(wp), dimension(size(x)) :: return_value
!
!  integer :: i
!
!continue
!
!  return_value = zero
!  forall (i = 1:size(x), abs(x(i)-xcen) < hw)
!    return_value(i) = one
!  end forall
!
!end function square_wave
!!
!!###############################################################################
!!
!pure function gaussian_wave(x,xcen,hw) result(return_value)
!
!  real(wp),               intent(in) :: xcen
!  real(wp),               intent(in) :: hw
!  real(wp), dimension(:), intent(in) :: x
!
!  real(wp), dimension(size(x)) :: return_value
!
!continue
!
!  return_value = exp( -log(two) * ((x-xcen)/hw)**2 )
!
!end function gaussian_wave
!!
!!###############################################################################
!!
!pure function ellipse_wave(x,xcen,hw) result(return_value)
!
!  real(wp),               intent(in) :: xcen
!  real(wp),               intent(in) :: hw
!  real(wp), dimension(:), intent(in) :: x
!
!  real(wp), dimension(size(x)) :: return_value
!
!  integer :: i
!
!continue
!
!  return_value = zero
!  forall (i = 1:size(x), abs(x(i)-xcen) < hw)
!    return_value(i) = sqrt( max( zero, one-((x(i)-xcen)/hw)**2 ) )
!  end forall
!
!end function ellipse_wave
!!
!!###############################################################################
!!
!pure function triangle_wave(x,xcen,hw) result(return_value)
!
!  real(wp),               intent(in) :: xcen
!  real(wp),               intent(in) :: hw
!  real(wp), dimension(:), intent(in) :: x
!
!  real(wp), dimension(size(x)) :: return_value
!
!  integer :: i
!
!continue
!
!  return_value = zero
!  forall (i = 1:size(x), abs(x(i)-xcen) < hw)
!    return_value(i) = one - abs((x(i)-xcen)/hw)
!  end forall
!
!end function triangle_wave
!
!###############################################################################
!
include "Functions/vortex.f90"
!
!###############################################################################
!
end module initialization_mod
