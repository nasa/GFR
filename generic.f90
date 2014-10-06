module generic_mod
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
  public :: getbc
  public :: rgetdtmax
  public :: get_dtsp
  public :: rresi
  public :: rupdate_Classic_RK
  public :: rupdate_CK4_RK
  public :: rupdate_TVD2_RK
  public :: rupdate_TVD3_RK
  public :: rupdate_TVD4_RK
  public :: rgeterr
  public :: getmaxent
  public :: compute_final_error
  public :: compute_final_advec_error
  public :: compute_final_vortex_error
  public :: test_error
  public :: project_to_porder
  public :: init_generic_mod
  !
  real(wp), parameter :: o = zero
  !
  ! Classic N stage Runge-Kutta coefficients
  !
  real(wp), parameter, dimension(4,4,4) :: &
  alphaRK = reshape( &                                                  ! Coefficients
       (/ONE ,o,o,o,o, o  ,o,o, o  ,o,  o   ,o, o  , o  , o  , o  ,   & ! 1-Stage
         HALF,o,o,o,o,ONE ,o,o, o  ,o,  o   ,o, o  , o  , o  , o  ,   & ! 2-Stage
         ONE3,o,o,o,o,TWO3,o,o,ONE4,o,THREE4,o, o  , o  , o  , o  ,   & ! 3-Stage
         HALF,o,o,o,o,HALF,o,o, o  ,o, ONE  ,o,ONE6,ONE3,ONE3,ONE6/), & ! 4-Stage
       (/4,4,4/) )
  !
  ! TVD 4th-order/5-stage Runge-Kutta coefficients
  !
  real(wp), parameter, dimension(1:2,1:5) :: &
    betaRK_TVD = reshape( (/ 0.00000000000000_wp,0.39175222700395_wp, &
                             0.00000000000000_wp,0.36841059262959_wp, &
                             0.00000000000000_wp,0.25189177424738_wp, &
                             0.00000000000000_wp,0.54497475021237_wp, &
                             0.08460416338212_wp,0.22600748319395_wp  /), &
                          (/2,5/) )
  !
  real(wp), parameter, dimension(1:3,1:5) :: &
    alphaRK_TVD = reshape( (/ 1.00000000000000_wp,0.00000000000000_wp,0.0_wp, &
                              0.44437049406734_wp,0.55562950593266_wp,0.0_wp, &
                              0.62010185138540_wp,0.37989814861460_wp,0.0_wp, &
                              0.17807995410773_wp,0.82192004589227_wp,0.0_wp, &
                              0.00683325884039_wp,0.34833675773694_wp,1.0_wp  /), &
                           (/3,5/) )
  !
  real(wp), parameter, dimension(1:5) :: gammaRK_TVD = (/ 0.00000000000000_wp, &
                                                           0.51723167208978_wp, &
                                                           0.12759831133288_wp, &
                                                           0.00000000000000_wp, &
                                                           0.00000000000000_wp /)
  !
  character(len=100), save :: iter_out_fmt
  !
contains
!
!###############################################################################
!
subroutine rgetdtmax
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use ovar,    only : gam,Final_Time,CFL,time,dt_in,d_t
  use metvar,  only : jacobian,dx,dy
  use flowvar, only : usp,dtsp
  !
  !.. Local Scalars ..
  integer   :: n
  real(wp) :: dxdr,dydr,dxds,dyds
  real(wp) :: rho,vx,vy,pres,aspd
  real(wp) :: alfa,beta,gama,dttemp
  real(wp) :: dt_global
  !
continue
  !
  !... time step for cfl = 1
  !
  if (abs(dt_in) <= eps12) then
    !
    dt_global = huge(zero)
    !
    do n = 1,n_solpts
      !
      dxdr = dx(n)%dr
      dydr = dy(n)%dr
      dxds = dx(n)%ds
      dyds = dy(n)%ds
      !
      rho  = usp(1,n)
      vx   = usp(2,n)/usp(1,n)
      vy   = usp(3,n)/usp(1,n)
      pres = pressure_cv_sp( usp(:,n) )
      !
      aspd = sqrt(gam*pres/rho)
      !
      alfa = abs( vx*dydr - vy*dxdr )
      beta = abs( vx*dyds - vy*dxds )
      !
      gama = aspd*sqrt(dxdr**2 + dydr**2 + dxds**2 + dyds**2)
      !
      dttemp = two * jacobian(n) / (alfa + beta + gama)
      dt_global = min(dt_global,dttemp)
      !
    end do
    !
    dttemp = dt_global
    if (ncpu > 1) then
      call mpi_allreduce(dttemp,dt_global,1,mpi_flttyp, &
                         mpi_min,mpi_comm_world,mpierr)
    end if
    !
    d_t = dt_global*CFL
    !
  else
    !
    d_t = dt_in
    !
  end if
  !
  if (Final_Time > zero) then
    if (time + d_t > Final_Time + eps9) then
      d_t = Final_Time - time
    end if
  end if
  !
  dtsp(:) = d_t
  !
  time = time + d_t
  !
end subroutine rgetdtmax
!
!###############################################################################
!
subroutine get_dtsp(dtmin)
  !
  !.. Use Statements ..
  use ovar, only : Timestep_Type,time,d_t
  !
  !.. rmal Arguments ..
  real(wp), intent(out) :: dtmin
  !
continue
  !
  if (Timestep_Type == Constant_Timestep) then
    !
    call get_constant_timestep(dtmin)
    !
  else if (abs(Timestep_Type) == Cell_Timestep) then
    !
    call get_cell_timestep(dtmin)
    !
  else if (abs(Timestep_Type) == Point_Timestep) then
    !
    call get_point_timestep(dtmin)
    !
  end if
  !
  time = time + dtmin
  d_t = dtmin
  !
end subroutine get_dtsp
!
!###############################################################################
!
subroutine get_constant_timestep(dtmin)
  !
  !.. Use Statements ..
  use ovar,    only : dt_in,Final_Time,time,CFL
  use flowvar, only : dtsp
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: dtmin
  !
continue
  !
  dtmin = CFL * dt_in
  !
  if (Final_Time > zero) then
    if (time + dtmin > Final_Time + eps9) then
      dtmin = Final_Time - time
    end if
  end if
  !
  dtsp(:) = dtmin
  !
end subroutine get_constant_timestep
!
!###############################################################################
!
subroutine get_cell_timestep(dtmin)
  !
  !.. Use Statements ..
  use geovar,  only : nr,n_sp2dt,nfbnd,nface,interp_s2f
  use geovar,  only : cells_on_face,pts_in_cell_ptr
  use ovar,    only : gam,Timestep_Type,Final_Time,time,CFL
  use metvar,  only : geofa,jacobian
  use flowvar, only : dtsp,faceusp
  use geovar,  only : nq => neqns
  use geovar,  only : nsp => n_sp1de
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: dtmin
  !
  !.. Local Scalars ..
  integer  :: nf,k
  integer  :: jl,kbl,l1,l2,npl,wl,kl
  integer  :: jr,kbr,r1,r2,npr,wr,kr
  real(wp) :: ds,eigni,eignv,dttemp,pp1sq
  real(wp) :: dsl,vell,al,jacl,dtl,mul,hl
  real(wp) :: dsr,velr,ar,jacr,dtr,mur,hr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nq) :: vl,vr
  !
continue
  !
  pp1sq = real( nsp**2 , kind=wp )
  !
  dtmin = huge(zero)
  dtsp(:) = dtmin
  !
  do nf = 1,nfbnd
    !
    jl = cells_on_face(1,nf) ! left cell on face
    !
    kbl = cells_on_face(5,nf) ! beginning flux point for left cell
    !
    l1 = pts_in_cell_ptr(jl)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(jl+1) ! ending    index for sol-pts in left cell
    !
    npl = l2-l1+1 ! number of solution points in left cell
    !
    wl = merge( 1 , 2 , npl == n_sp2dt ) ! check if jl is a tria or quad
    !
    do k = 1,nsp
      !
      kl = kbl-1 + k ! index for left cell
      !
      vl(1:nq) = usp2v_sp( faceusp(1:nq,k,1,nf) )
      !
      mul = viscosity_cv_sp( faceusp(1:nq,k,1,nf) )
      !
      ds = rnorm2( geofa(1:nr,k,nf) )
      dsl = ds * cl(wl) ! face length/area for left cell
      !
      ! Normal velocity at this flux point
      !
     !vell = dot( vl(2:nq-1) , unit_vector(geofa(1:nr,k,nf)) )
      vell = selfdot( vl(2:nq-1) )
      !
      ! Speed of sound at this face
      !
      al = max( gam*vl(nq)/vl(1) , rndoff )
      !
      ! Interpolate the jacobian to the flux point
      !
      jacl = dot( interp_s2f(1:npl,kl,wl) , jacobian(l1:l2) )
      !
      hl = half * dsl / jacl
      !
      eigni = abs(vell) + al
      eignv = pp1sq * hl * mul
      !
      ! Compute the timestep for this flux point
      !
      dtl = min( dtsp(l1) , CFL / (pp1sq * hl * (eigni + eignv)) )
      !
      ! Store the minimum timestep for this cell in the location
      ! of the first solution point for this cell
      !
      dtsp(l1:l2) = dtl
      !
      dtmin = min( dtmin , dtl )
      !
    end do
    !
  end do
  !
  do nf = nfbnd+1,nface
    !
    jl = cells_on_face(1,nf) ! left  cell on face
    jr = cells_on_face(2,nf) ! right cell on face
    !
    kbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    kbr = cells_on_face(6,nf) ! beginning flux point for right cell
    !
    l1 = pts_in_cell_ptr(jl)+1 ! beginning index for sol-pts in left  cell
    l2 = pts_in_cell_ptr(jl+1) ! ending    index for sol-pts in right cell
    !
    r1 = pts_in_cell_ptr(jr)+1 ! beginning index for sol-pts in left  cell
    r2 = pts_in_cell_ptr(jr+1) ! ending    index for sol-pts in right cell
    !
    npl = l2-l1+1 ! number of solution points in left  cell
    npr = r2-r1+1 ! number of solution points in right cell
    !
    wl = merge( 1 , 2 , npl == n_sp2dt ) ! check if jl is a tria or quad
    wr = merge( 1 , 2 , npr == n_sp2dt ) ! check if jr is a tria or quad
    !
    do k = 1,nsp
      !
      kl = kbl-1 + k ! index for left  cell
      kr = kbr+1 - k ! index for right cell
      !
      vl(1:nq) = usp2v_sp( faceusp(1:nq,k,1,nf) )
      vr(1:nq) = usp2v_sp( faceusp(1:nq,k,2,nf) )
      !
      mul = viscosity_cv_sp( faceusp(1:nq,k,1,nf) )
      mur = viscosity_cv_sp( faceusp(1:nq,k,2,nf) )
      !
      ds = rnorm2( geofa(1:nr,k,nf) )
      dsl = ds * cl(wl) ! face length/area for left  cell
      dsr = ds * cl(wr) ! face length/area for right cell
      !
      ! Normal velocity at this flux point
      !
     !vell = dot( vl(2:nq-1) , unit_vector(geofa(1:nr,k,nf)) )
     !velr = dot( vr(2:nq-1) , unit_vector(geofa(1:nr,k,nf)) )
      vell = selfdot( vl(2:nq-1) )
      velr = selfdot( vr(2:nq-1) )
      !
      ! Speed of sound at this face
      !
      al = max( gam*vl(nq)/vl(1) , rndoff )
      ar = max( gam*vr(nq)/vr(1) , rndoff )
      !
      ! Interpolate the jacobian to the flux point
      !
      jacl = dot( interp_s2f(1:npl,kl,wl) , jacobian(l1:l2) )
      jacr = dot( interp_s2f(1:npr,kr,wr) , jacobian(r1:r2) )
      !
      hl = half * dsl / jacl
      hr = half * dsr / jacr
      !
      eigni = max( abs(vell) + al , abs(velr) + ar )
      eignv = pp1sq * max( hl * mul , hr * mur )
      !
      ! Compute the timestep for this flux point
      !
      dtl = min( dtsp(l1) , CFL / (pp1sq * hl * (eigni + eignv) ) )
      dtr = min( dtsp(r1) , CFL / (pp1sq * hr * (eigni + eignv) ) )
      !
      ! Store the minimum timestep for all solution points in both cells
      !
      dtsp(l1:l2) = dtl
      dtsp(r1:r2) = dtr
      !
      dtmin = min( dtmin , dtl , dtr )
      !
    end do
    !
  end do
  !
  ! Get the minimum timestep across all processors
  !
  if (ncpu > 1) then
    dttemp = dtmin
    call mpi_allreduce(dttemp,dtmin,1,mpi_flttyp, &
                       mpi_min,mpi_comm_world,mpierr)
  end if
  !
  ! Limit the local timesteps if using a global timestepping
  !
  if (Timestep_Type == Global_Cell_Timestep) then
    !
    ! Set all time step sizes to the global minimum timestep
    !
    dtsp(:) = dtmin
    !
    ! Make sure the timestep size doesnt make the final time exceed Final_Time
    !
    if (Final_Time > zero) then
      if (time + dtmin > Final_Time + eps9) then
        dtsp(:) = Final_Time - time
      end if
    end if
    !
  end if
  !
end subroutine get_cell_timestep
!
!###############################################################################
!
subroutine get_point_timestep(dtmin)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts,n_sp1de
  use ovar,    only : gam,Final_Time,time,CFL,Timestep_Type
  use metvar,  only : jacobian,dx,dy
  use flowvar, only : usp,dtsp
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: dtmin
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: dtpt,gm1,dtc,eign
  real(wp) :: dxdr,dydr,dxds,dyds
  real(wp) :: rho,vx,vy,pres,aspd
  real(wp) :: alfa,beta,gama,dttemp
  !
continue
  !
  gm1 = gam - one
  !
  dtmin = huge(zero)
  !
  dtc = CFL * two / real( n_sp1de**2 , kind=wp )
  !
  do n = 1,n_solpts
    !
    dxdr = dx(n)%dr
    dydr = dy(n)%dr
    dxds = dx(n)%ds
    dyds = dy(n)%ds
    !
    rho  = usp(1,n)
    vx   = usp(2,n)/rho
    vy   = usp(3,n)/rho
    pres = gm1 * (usp(4,n) - half*rho*(vx*vx + vy*vy))
    !
    aspd = sqrt( max( gam*pres/rho , rndoff ) )
    !
    alfa = abs( vx*dydr - vy*dxdr )
    beta = abs( vx*dyds - vy*dxds )
    !
    gama = aspd*sqrt(dxdr**2 + dydr**2 + dxds**2 + dyds**2)
    !
    eign = alfa + beta + gama
    !
    dtpt = dtc * jacobian(n) / eign
    !
    dtmin = min( dtpt , dtmin )
    !
    dtsp(n) = dtpt
    !
  end do
  !
  ! Get the minimum timestep across all processors
  !
  if (ncpu > 1) then
    dttemp = dtmin
    call mpi_allreduce(dttemp,dtmin,1,mpi_flttyp, &
                       mpi_min,mpi_comm_world,mpierr)
  end if
  !
  ! Limit the local timesteps if using a global timestepping
  !
  if (Timestep_Type == Global_Point_Timestep) then
    !
    ! Set all time step sizes to the minimum timestep
    !
    dtsp(:) = dtmin
    !
    ! Make sure the timestep size doesnt make the final time exceed Final_Time
    !
    if (Final_Time > zero) then
      if (time + dtmin > Final_Time + eps9) then
        dtsp(:) = Final_Time - time
      end if
    end if
    !
  end if
  !
end subroutine get_point_timestep
!
!###############################################################################
!
subroutine rresi
  !
  !... calculates d(utd)/dt
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts,neqns
  use geovar,  only : ncell,pts_in_cell_ptr
  use ovar,    only : itestcase
  use metvar,  only : jacobian
  use flowvar, only : residual
  !
  !.. Local Scalars ..
  integer :: n,m
  !
continue
  !
  do n = 1,n_solpts
    residual(:,n) = -residual(:,n) / jacobian(n)
  end do
  !
  ! Zero out the residuals for all but the density
  ! equation if running the density transport test case
  !
  if (itestcase == Density_Transport) then
    do n = 1,n_solpts
      do m = 2,neqns
        residual(m,n) = residual(1,n)
      end do
    end do
  end if
  !
  call project_to_porder(ncell,pts_in_cell_ptr,residual)
  !
end subroutine rresi
!
!###############################################################################
!
subroutine project_to_porder(ncell,pts_in_cell_ptr,array,ptype)
  !
  !... calculates d(utd)/dt
  !
  !.. Use Statements ..
  use geovar,  only : n_sp2dt,n_sp2dq
  use geovar,  only : p_order,n_order
  use geovar,  only : projct_s2p
  !
  ! PTYPE : gives the location of the variables needing projection
  !         = 1 -> at solution points
  !         = 2 -> at Gauss    points
  !         = 3 -> at Lobatto  points
  !         = 4 -> at Radau    points
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: ncell
  integer,  dimension(:),      intent(in) :: pts_in_cell_ptr
  real(wp), dimension(:,:), intent(inout) :: array
  integer,  optional,          intent(in) :: ptype
  !
  !.. Local Scalars ..
  integer :: k,n,l,n1,n2,nsp,nw,ndim,ploc
  !
  !.. Local Arrays ..
  real(wp), dimension(1:n_sp2dq,1:size(array,dim=1)) :: tmp
  !
continue
  !
  if (p_order < n_order) then
    !
    ploc = 1
    if (present(ptype)) then
      if (ptype > 0 .and. ptype < 5) then
        ploc = ptype
      end if
    end if
    !
    ndim = size(array,dim=1)
    !
    do n = 1,ncell
      !
      n1 = pts_in_cell_ptr(n)+1
      n2 = pts_in_cell_ptr(n+1)
      nsp = n2-n1+1
      nw = merge(1,2,nsp==n_sp2dt)
      !
      do k = 1,nsp
        do l = 1,ndim
          tmp(k,l) = array(l,n1-1+k)
        end do
      end do
      !
      do k = 1,nsp
        do l = 1,ndim
          array(l,n1-1+k) = dot( projct_s2p(1:nsp,k,nw,ploc) , tmp(1:nsp,l) )
        end do
      end do
      !
    end do
    !
  end if
  !
end subroutine project_to_porder
!
!###############################################################################
!
subroutine rupdate_Classic_RK
  !
  ! First, store d(utd)/dt for runge kutta time
  ! stepping then update usp
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts,neqns,write_variable
  use ovar,    only : nst,num_stages,Filter_Option
  use flowvar, only : rssprk,residual,usp,uoldsp,dtsp
  !
  !.. Local Scalars ..
  integer   :: n,m,nt
  real(wp) :: temp
  !
continue
  !
  ! Store the residual for this Runge-Kutta step in rssprk
  !
 !call write_variable(101,residual)
  !
  rssprk(:,:,nst) = residual(:,:)
  !
  ! Update usp
  !
  do n = 1,n_solpts
    do m = 1,neqns
      temp = zero
      do nt = 1,nst
        temp = temp + alphark(nt,nst,num_stages)*rssprk(m,n,nt)
      end do
      usp(m,n) = uoldsp(m,n) + dtsp(n)*temp
     !usp(m,n) = uoldsp(m,n) &
     !         + dtsp(n)*sum( alphark(1:nst,nst,num_stages)*rssprk(m,n,1:nst) )
    end do
  end do
  !
  ! Apply the filter
  !
  if (Filter_Option > 0) call filter_usp
  !
end subroutine rupdate_Classic_RK
!
!###############################################################################
!
subroutine rupdate_CK4_RK
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use ovar,    only : nst,Filter_Option
  use flowvar, only : rssprk,residual,usp,dtsp
  !
  !.. Local Scalars ..
  integer :: n
  !
  ! Carpenter-Kennedy 4th-order/5-stage Runge-Kutta coefficients
  !
  real(wp), parameter, dimension(1:5) :: &
             alfa = (/                                   0.0_wp, &
                        -567301805773.0_wp/ 1357537059087.0_wp, &
                       -2404267990393.0_wp/ 2016746695238.0_wp, &
                       -3550918686646.0_wp/ 2091501179385.0_wp, &
                       -1275806237668.0_wp/  842570457699.0_wp  /)

  real(wp), parameter, dimension(1:5) :: &
             beta = (/  1432997174477.0_wp/ 9575080441755.0_wp, &
                        5161836677717.0_wp/13612068292357.0_wp, &
                        1720146321549.0_wp/ 2090206949498.0_wp, &
                        3134564353537.0_wp/ 4481467310338.0_wp, &
                        2277821191437.0_wp/14882151754819.0_wp  /)
  !
continue
  !
  do n = 1,n_solpts
    !
    rssprk(:,n,1) = alfa(nst)*rssprk(:,n,1) + dtsp(n)*residual(:,n)
    !
    usp(:,n) = usp(:,n) + beta(nst)*rssprk(:,n,1)
    !
  end do
  !
  ! Apply the filter
  !
  if (Filter_Option > 0) call filter_usp
  !
end subroutine rupdate_CK4_RK
!
!###############################################################################
!
subroutine rupdate_TVD2_RK
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use ovar,    only : nst,Filter_Option
  use flowvar, only : residual,usp,uoldsp,dtsp
  !
  !.. Local Scalars ..
  integer   :: n
  !
  ! TVD 2nd-order/2-stage Runge-Kutta Coefficients
  !
  real(wp), parameter, dimension(1:2) :: alfa = (/ zero,one  /)
  real(wp), parameter, dimension(1:2) :: beta = (/ one ,half /)
  !
continue
  !
  do n = 1,n_solpts
    !
    usp(:,n) = beta(nst) * (alfa(nst)*uoldsp(:,n) + usp(:,n) + dtsp(n)*residual(:,n))
    !
  end do
  !
  ! Apply the filter
  !
  if (Filter_Option > 0) call filter_usp
  !
end subroutine rupdate_TVD2_RK
!
!###############################################################################
!
subroutine rupdate_TVD3_RK
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use ovar,    only : nst,Filter_Option
  use flowvar, only : residual,usp,uoldsp,dtsp
  !
  !.. Local Scalars ..
  integer :: n
  !
  ! TVD 2nd-order/2-stage Runge-Kutta Coefficients
  !
  real(wp), parameter, dimension(1:3) :: alfa = (/ zero,three,half /)
  real(wp), parameter, dimension(1:3) :: beta = (/ one ,one4 ,two3 /)
  !
continue
  !
  do n = 1,n_solpts
    !
    usp(:,n) = beta(nst) * (alfa(nst)*uoldsp(:,n) + usp(:,n) + dtsp(n)*residual(:,n))
    !
  end do
  !
  ! Apply the filter
  !
  if (Filter_Option > 0) call filter_usp
  !
end subroutine rupdate_TVD3_RK
!
!###############################################################################
!
subroutine rupdate_TVD4_RK
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use ovar,    only : nst,Filter_Option
  use flowvar, only : rssprk,residual,usp,uoldsp,usp_sum,dtsp
  !
  !.. Local Scalars ..
  integer :: n
  !
  ! TVD 4th-order/5-stage Runge-Kutta Coefficients
  !
  real(wp), parameter, dimension(1:2,1:5) :: &
    beta = reshape( (/ 0.00000000000000_wp,0.39175222700395_wp, &
                       0.00000000000000_wp,0.36841059262959_wp, &
                       0.00000000000000_wp,0.25189177424738_wp, &
                       0.00000000000000_wp,0.54497475021237_wp, &
                       0.08460416338212_wp,0.22600748319395_wp  /), &
                    (/2,5/) )
  !
  real(wp), parameter, dimension(1:3,1:5) :: &
    alfa = reshape( (/ 1.00000000000000_wp,0.00000000000000_wp,0.0_wp, &
                       0.44437049406734_wp,0.55562950593266_wp,0.0_wp, &
                       0.62010185138540_wp,0.37989814861460_wp,0.0_wp, &
                       0.17807995410773_wp,0.82192004589227_wp,0.0_wp, &
                       0.00683325884039_wp,0.34833675773694_wp,1.0_wp  /), &
                    (/3,5/) )
  !
  real(wp), parameter, dimension(1:5) :: gama = (/ 0.00000000000000_wp, &
                                                    0.51723167208978_wp, &
                                                    0.12759831133288_wp, &
                                                    0.00000000000000_wp, &
                                                    0.00000000000000_wp  /)
  !
continue
  !
  if (nst == 1) usp_sum(:,:) = zero
  !
  do n = 1,n_solpts
    !
    rssprk(:,n,1) = beta(1,nst)*rssprk(:,n,1) + beta(2,nst)*residual(:,n)
    !
    usp(:,n) = alfa(1,nst)* uoldsp(:,n) + &
               alfa(2,nst)*    usp(:,n) + &
               alfa(3,nst)*usp_sum(:,n) + dtsp(n)*rssprk(:,n,1)
    !
    usp_sum(:,n) = usp_sum(:,n) + gama(nst)*usp(:,n)
    !
  end do
  !
  ! Apply the filter
  !
  if (Filter_Option > 0) call filter_usp
  !
end subroutine rupdate_TVD4_RK
!
!###############################################################################
!
subroutine filter_usp
  !
  !.. Use Statements ..
  use geovar,  only : neqns,n_sp2dt,n_sp2dq,ncell, &
                      pts_in_cell_ptr,filter
  use flowvar, only : usp
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,np,m,w1
  !
  !.. Local Arrays ..
  real(wp), dimension(1:n_sp2dq) :: var
  !
continue
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1
    n2 = pts_in_cell_ptr(n+1)
    np = n2-n1+1
    !
    w1 = merge(1,2,np==n_sp2dt)
    !
    do m = 1,neqns
      var(1:np) = usp(m,n1:n2)
      var(1:np) = matmul( filter(1:np,1:np,w1) , var(1:np) )
      usp(m,n1:n2) = var(1:np)
    end do
    !
  end do
  !
end subroutine filter_usp
!
!###############################################################################
!
subroutine compute_final_error()
  !
  ! Compute the final error. The error is
  ! stored in the arrays ultwo and vltwo as
  !
  ! ultwo(1)      = density
  ! ultwo(2:nq-1) = directional momentums
  ! ultwo(nq)     = energy
  ! ultwo(nq+1)   = entropy
  !
  ! vltwo(1)      = density
  ! vltwo(2:nq-1) = directional velocities
  ! vltwo(nq)     = pressure
  ! vltwo(nq+1)   = magnitude of velocity
  !
  !.. Use Statements ..
  use geovar,  only : nr,neqns,ncell,pts_in_cell_ptr,n_sp1de,n_sp2dt,n_sp2dq
  use geovar,  only : weights2d,numer2gauss,gauss_wts_2d
  use ovar,    only : vlinf,vlone,vltwo,vnloc,vkloc,total_volume, &
                      ulinf,ulone,ultwo,unloc,ukloc,vxref,vyref,aref
  use metvar,  only : jacobian
  use flowvar, only : usp,uoldsp
  use parallel_mod, only : cell_map,cell_map_ptr
  !
  !.. Local Scalars ..
  integer  :: n,n1,n2,np,k,l,w1,m,nc1,ierr
  real(wp) :: rho_ratio,pres_ratio,weight
  real(wp) :: delx,dely,xoffset,yoffset,sn,se
  !
  !.. Local Arrays ..
  real(wp), dimension(1:neqns,1:n_sp2dq) :: uold
  real(wp), dimension(1:neqns,1:n_sp2dq) :: vold
  real(wp), dimension(1:neqns,1:n_sp2dq) :: unew
  real(wp), dimension(1:neqns,1:n_sp2dq) :: vnew
  real(wp), dimension(1:neqns+1,1:2) :: ett
  real(wp), dimension(1:neqns+1,1:2,1:5) :: error
  real(wp), dimension(1:neqns+1,1:2,1:5,1:ncpu) :: all_error
  !
  ! error(1:neqns+1,:,:) = error value for each flow equation + extra variable
  ! error(:,1,:) = error for conserved variables [rho,rho(u,v,..),rhoe,entropy]
  ! error(:,2,:) = error for primitive variables [rho,u,v,..,pressure,mag of vel]
  ! error(:,:,1) = L1 norm of total error for all solution points
  ! error(:,:,2) = L2 norm of total error for all solution points
  ! error(:,:,3) = maximum value of error for all solution points
  ! error(:,:,4) = ID of grid cell containing maximum value of error
  ! error(:,:,5) = ID of solution point with maximum value of error
  !
continue
  !
  if (.not. allocated(unloc)) then
    allocate ( unloc(1:neqns+1) , source=0 , stat=ierr )
    call alloc_error("compute_final_error","unloc",1,ierr)
  end if
  if (.not. allocated(ukloc)) then
    allocate ( ukloc(1:neqns+1) , source=0 , stat=ierr )
    call alloc_error("compute_final_error","ukloc",1,ierr)
  end if
  if (.not. allocated(ulinf)) then
    allocate ( ulinf(1:neqns+1) , source=zero , stat=ierr )
    call alloc_error("compute_final_error","ulinf",1,ierr)
  end if
  if (.not. allocated(ulone)) then
    allocate ( ulone(1:neqns+1) , source=zero , stat=ierr )
    call alloc_error("compute_final_error","ulone",1,ierr)
  end if
  if (.not. allocated(ultwo)) then
    allocate ( ultwo(1:neqns+1) , source=zero , stat=ierr )
    call alloc_error("compute_final_error","ultwo",1,ierr)
  end if
  if (.not. allocated(vnloc)) then
    allocate ( vnloc(1:neqns+1) , source=0 , stat=ierr )
    call alloc_error("compute_final_error","vnloc",1,ierr)
  end if
  if (.not. allocated(vkloc)) then
    allocate ( vkloc(1:neqns+1) , source=0 , stat=ierr )
    call alloc_error("compute_final_error","vkloc",1,ierr)
  end if
  if (.not. allocated(vlinf)) then
    allocate ( vlinf(1:neqns+1) , source=zero , stat=ierr )
    call alloc_error("compute_final_error","vlinf",1,ierr)
  end if
  if (.not. allocated(vlone)) then
    allocate ( vlone(1:neqns+1) , source=zero , stat=ierr )
    call alloc_error("compute_final_error","vlone",1,ierr)
  end if
  if (.not. allocated(vltwo)) then
    allocate ( vltwo(1:neqns+1) , source=zero , stat=ierr )
    call alloc_error("compute_final_error","vltwo",1,ierr)
  end if
  !
  error(:,:,:) = zero  ! intialize everything to zero
  error(:,:,3) = -eps6 ! make this slightly negative in case all error is zero
  !
  nc1 = cell_map_ptr(mypnum+1)
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1
    n2 = pts_in_cell_ptr(n+1)
    np = n2-n1+1
    !
    w1 = merge(1,2,np==n_sp2dt)
    !
    ! Interpolate the new and old solution to the Gauss points
    ! to use Gaussian quadrature to integrate the change in the
    ! solution since the previous time step
    !
    do m = 1,neqns
      unew(m,1:np) = matmul( numer2gauss , usp(m,n1:n2) )
      uold(m,1:np) = matmul( numer2gauss , uoldsp(m,n1:n2) )
    end do
    !
    ! Get the primitive variables for the numerical and analytical solutions
    !
    vold = usp2v(uold)
    vnew = usp2v(unew)
    !
    ! Now loop through the solution points of this cell to compute
    ! the error between the simulation and analytical solutions.
    !
    do l = n1,n2
      !
      k = l-n1+1
      !
      ! Error for conserved variables at the current solution point
      !
      sn = entropy_pv_sp( vnew(:,k) )
      se = entropy_pv_sp( vold(:,k) )
      !
      ett(1:neqns,1) = abs( unew(:,k) - uold(:,k) )
      ett(neqns+1,1) = sn - se
      !
      ! Error for primitive variables at the current solution point
      !
      ett(1:neqns,2) = abs( vnew(:,k) - vold(:,k) )
      ett(neqns+1,2) = rnorm2( ett(2:neqns-1,2) )
      !
      ! Now add the error for all variables at the current solution
      ! point to the running integral sum for the L1 and L2 error norms
      !
      ! ##############################################################
      ! ##############################################################
      ! ##############################################################
      ! ####################  !!!! WARNING !!!!  #####################
      ! ##############################################################
      ! ##############################################################
      ! ##############################################################
      ! ###   GAUSS_WTS_2D NEEDS TO BE SET UP FOR TRIANGLE CELLS   ###
      ! ##############################################################
      ! ##############################################################
      ! ##############################################################
      ! ####################  !!!! WARNING !!!!  #####################
      ! ##############################################################
      ! ##############################################################
      ! ##############################################################
      !
     !weight = jacobian(l)*weights2d(k,w1)
      weight = jacobian(l)*gauss_wts_2d(k)
      !
      error(:,:,1) = error(:,:,1) + weight*ett(:,:)
      error(:,:,2) = error(:,:,2) + weight*ett(:,:)*ett(:,:)
      !
      ! If any of these errors are greater than the previous maximum
      ! values, save the error value and record this location
      !
      where (abs(ett(:,:)) > error(:,:,3))
        error(:,:,3) = abs(ett(:,:))
        error(:,:,4) = real(cell_map(nc1+n),kind=wp)
        error(:,:,5) = real(k,kind=wp)
      end where
      !
    end do
    !
  end do
  !
  ! Collect the error information from all processors
  !
  if (ncpu > 1) then
    call mpi_allgather(error,    size(error),mpi_flttyp, &
                       all_error,size(error),mpi_flttyp,mpi_comm_world,mpierr)
  else
    all_error(:,:,:,1) = error(:,:,:)
  end if
  !
  ! Find the max infinity norm for each variable
  !
  ulinf(:) = maxval(all_error(:,1,3,:),dim=2)
  vlinf(:) = maxval(all_error(:,2,3,:),dim=2)
  !
  ! Finish computing the L1 and L2 norms for each variable
  !
  ulone(:) =       sum(all_error(:,1,1,:),dim=2) / total_volume
  ultwo(:) = sqrt( sum(all_error(:,1,2,:),dim=2) / total_volume )
  vlone(:) =       sum(all_error(:,2,1,:),dim=2) / total_volume
  vltwo(:) = sqrt( sum(all_error(:,2,2,:),dim=2) / total_volume )
  !
  ! Get the location of the max infinity norm for each variable
  !
  do m = 1,neqns+1
    unloc(m) = nint( all_error(m,1,4,maxloc(all_error(m,1,3,:),dim=1)) )
    ukloc(m) = nint( all_error(m,1,5,maxloc(all_error(m,1,3,:),dim=1)) )
    vnloc(m) = nint( all_error(m,2,4,maxloc(all_error(m,2,3,:),dim=1)) )
    vkloc(m) = nint( all_error(m,2,5,maxloc(all_error(m,2,3,:),dim=1)) )
  end do
  !
end subroutine compute_final_error
!
!###############################################################################
!
subroutine compute_final_vortex_error()
  !
  ! Compute the final error and store it in the vortex_error.
  !
  ! vortex_error(1,:,:) = error value for density
  ! vortex_error(2,:,:) = error value for vx
  ! vortex_error(3,:,:) = error value for vy
  ! vortex_error(4,:,:) = error value for pressure
  ! vortex_error(5,:,:) = error value for entropy
  ! vortex_error(:,1,:) = error computed using gauss quadrature
  ! vortex_error(:,2,:) = error computed using lobatto quadrature
  ! vortex_error(:,:,1) = L1 norm of total error for all solution points
  ! vortex_error(:,:,2) = L2 norm of total error for all solution points
  ! vortex_error(:,:,3) = infinity norm of total error for all solution points
  !                      i.e., maximum value of error for all solution points
  ! vortex_error(:,:,4) = ID of grid cell containing maximum value of error
  ! vortex_error(:,:,5) = ID of solution point with maximum value of error
  !
  !.. Use Statements ..
  use geovar,  only : nr,neqns,ncell,pts_in_cell_ptr,n_sp1de,n_sp2dt,n_sp2dq
  use geovar,  only : xyz,xyz_gauss,xyz_lobatto
  use geovar,  only : numer2gauss,numer2lobatto
  use geovar,  only : weights2d,gauss_wts_2d,lobatto_wts_2d
  use ovar,    only : total_volume,pv_in,gam,itestcase,itcur,vortex_error, &
                      Period_Timesteps,Period_Time,Grid_Min_Y,Grid_Length_Y, &
                      Grid_Min_X,Grid_Length_X,vxref,vyref,aref, &
                      Density_Maximum_Residual,vortex_has_collapsed
  use ovar,    only : bound_error,error_bound_limit
  use metvar,  only : jacobian,volj
  use flowvar, only : usp
  use parallel_mod, only : cell_map,cell_map_ptr
  !
  !.. Local Scalars ..
  integer   :: n,n1,n2,np,k,l,w1,m,nc1
  real(wp) :: rho_ratio,pres_ratio,c1,wg,wr,wl
  real(wp) :: delx,dely,xoffset,yoffset
  real(wp) :: sn,se
  real(wp) :: my_volume,error_volume
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: cell_center
  real(wp), dimension(1:nr,1:n_sp2dq) :: xyz_g
  real(wp), dimension(1:nr,1:n_sp2dq) :: xyz_l
  !
  real(wp), dimension(1:neqns,1:n_sp2dq) :: vge
  real(wp), dimension(1:neqns,1:n_sp2dq) :: vle
  !
  real(wp), dimension(1:neqns,1:n_sp2dq) :: vgn
  real(wp), dimension(1:neqns,1:n_sp2dq) :: vln
  !
  real(wp), dimension(1:neqns+1,1:2) :: ett
  real(wp), dimension(1:neqns+1,1:2,1:5) :: error
  real(wp), dimension(1:neqns+1,1:2,1:5,1:ncpu) :: all_error
  !
continue
  !
  my_volume = zero
  error_volume = zero
  !
  ! If running one of the vortex transport problems, compute the change
  ! in position of the vortex since the start of the current period.
  !
  if (any(itestcase == Transport_Problems)) call vortex_position(delx,dely)
  !
  error(:,:,:) = zero  ! intialize everything to zero
  error(:,:,3) = -eps6 ! make this slightly negative in case all error is zero
  !
  c1 = deltafun(itestcase==Density_Transport)
  !
  nc1 = cell_map_ptr(mypnum+1)
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1
    n2 = pts_in_cell_ptr(n+1)
    np = n2-n1+1
    !
    w1 = merge(1,2,np==n_sp2dt)
    !
    ! Copy the numerical solution to a local array
    !
    vgn(:,1:np) = usp(:,n1:n2)
    vln(:,1:np) = usp(:,n1:n2)
    !
    ! Get the solution point coordinates for this cell
    ! at the Gauss and Lobatto quadrature points
    !
    xyz_g(:,1:np) = xyz_gauss(:,n1:n2)
    xyz_l(:,1:np) = xyz_lobatto(:,n1:n2)
    !
    if (bound_error) then
      cell_center(:) = sum(xyz_g(:,1:np),dim=2)/real(np,kind=wp)
      cell_center(1) = cell_center(1) - &
                       merge(delx,delx-Grid_Length_X,cell_center(1)-delx >= Grid_Min_X)
      cell_center(2) = cell_center(2) - &
                       merge(dely,dely-Grid_Length_Y,cell_center(2)-dely >= Grid_Min_Y)
    else
      cell_center(:) = zero
    end if
    !
    ! Get the exact solution for this cell for both Gauss and Lobatto quadratures
    !
    do k = 1,np ! Gauss
      !
      xoffset = merge(delx,delx-Grid_Length_X,xyz_g(1,k)-delx >= Grid_Min_X)
      yoffset = merge(dely,dely-Grid_Length_Y,xyz_g(2,k)-dely >= Grid_Min_Y)
      !
      vge(:,k) = vortex2D( xyz_g(1,k)-xoffset , xyz_g(2,k)-yoffset )
      !
    end do
    !
    do k = 1,np ! Lobatto
      !
      xoffset = merge(delx,delx-Grid_Length_X,xyz_l(1,k)-delx >= Grid_Min_X)
      yoffset = merge(dely,dely-Grid_Length_Y,xyz_l(2,k)-dely >= Grid_Min_Y)
      !
      vle(:,k) = vortex2D( xyz_l(1,k)-xoffset , xyz_l(2,k)-yoffset )
      !
    end do
    !
    ! Project the exact solutions to p_order
    !
    call project_to_porder(1,[0,np],vge(1:neqns,1:np),2)
    call project_to_porder(1,[0,np],vle(1:neqns,1:np),3)
    !
    ! Interpolate the numerical solution to the solution
    ! nodes for Guass and Lobatto quadratures
    !
    do m = 1,neqns
      vgn(m,1:np) = matmul( numer2gauss   , vgn(m,1:np) )
      vln(m,1:np) = matmul( numer2lobatto , vln(m,1:np) )
    end do
    !
    ! Get the primitive variables for the numerical and analytical solutions
    !
    vge = usp2v(vge)
    vle = usp2v(vle)
    !
    vgn = usp2v(vgn)
    vln = usp2v(vln)
    !
    if (all(abs(cell_center) <= error_bound_limit)) then
      !
      my_volume = my_volume + volj(n)
      !
      ! Now loop through the solution points of this cell to compute
      ! the error between the simulation and analytical solutions.
      !
      do l = n1,n2
        !
        k = l-n1+1
        !
        ! Gauss error
        !
       !rho_ratio  = vgn(    1,k)/pv_in(    1)
       !pres_ratio = vgn(neqns,k)/pv_in(neqns)
        sn = entropy_pv_sp( vgn(:,k) )
        se = entropy_pv_sp( vge(:,k) )
        !
        ett(1:neqns,1) = abs( vgn(:,k) - vge(:,k) )
       !ett(neqns+1,1) = log( (c1+pres_ratio) / (rho_ratio**gam) )
        ett(neqns+1,1) = sn - se
        !
        ! Lobatto error
        !
       !rho_ratio  = vln(    1,k)/pv_in(    1)
       !pres_ratio = vln(neqns,k)/pv_in(neqns)
        sn = entropy_pv_sp( vln(:,k) )
        se = entropy_pv_sp( vle(:,k) )
        !
        ett(1:neqns,2) = abs( vln(:,k) - vle(:,k) )
       !ett(neqns+1,2) = log( (c1+pres_ratio) / (rho_ratio**gam) )
        ett(neqns+1,2) = sn - se
        !
        ! Now add the error for all variables at the current solution
        ! point to the running integral sum for the L1 and L2 error norms
        !
        wg = jacobian(l)*gauss_wts_2d(k)
        wl = jacobian(l)*lobatto_wts_2d(k)
        !
        error(:,1,1) = error(:,1,1) + wg*ett(:,1)
        error(:,1,2) = error(:,1,2) + wg*ett(:,1)**2
        !
        error(:,2,1) = error(:,2,1) + wl*ett(:,2)
        error(:,2,2) = error(:,2,2) + wl*ett(:,2)**2
        !
        ! If any of these errors are greater than the previous maximum
        ! values, save the error value and record this location
        !
        where (abs(ett(:,:)) > error(:,:,3))
          error(:,:,3) = abs(ett(:,:))
          error(:,:,4) = real(cell_map(nc1+n),kind=wp)
          error(:,:,5) = real(k,kind=wp)
        end where
        !
      end do
      !
    else if (any(abs(cell_center) > [Grid_Length_X,Grid_Length_Y]/two)) then
      !
      write (iout,*) "shouldnt be here!!!"
      !
    end if
    !
  end do
  !
  ! Collect the error information from all processors
  !
  if (ncpu > 1) then
    call mpi_allgather(error,    size(error),mpi_flttyp, &
                       all_error,size(error),mpi_flttyp,mpi_comm_world,mpierr)
    call mpi_allreduce(my_volume,error_volume,1,mpi_flttyp, &
                       mpi_sum,mpi_comm_world,mpierr)
  else
    all_error(:,:,:,1) = error(:,:,:)
    error_volume = my_volume
  end if
  !
  ! Finish computing the global error across all processors
  !
 !error_volume = total_volume
  !
  do n = 1,2
    do m = 1,neqns+1
      vortex_error(m,n,1) =       sum(all_error(m,n,1,:)) / error_volume
      vortex_error(m,n,2) = sqrt( sum(all_error(m,n,2,:)) / error_volume )
      vortex_error(m,n,3) = maxval(all_error(m,n,3,:))
      vortex_error(m,n,4) = all_error(m,n,4,maxloc(all_error(m,n,3,:),dim=1))
      vortex_error(m,n,5) = all_error(m,n,5,maxloc(all_error(m,n,3,:),dim=1))
    end do
  end do
  !
  ! Check for Vortex Collapse
  !
  if (vortex_error(1,1,2) > zero) then
    if (log10(vortex_error(1,1,2)) > Density_Maximum_Residual) then
      Density_Maximum_Residual = log10(vortex_error(1,1,2))
    end if
    vortex_has_collapsed = Density_Maximum_Residual > -two
  end if
  !
end subroutine compute_final_vortex_error
!
!###############################################################################
!
subroutine compute_final_advec_error()
  !
  ! Compute the final error and store it in the advec_error.
  !
  ! advec_error(1,:,:) = error value for density
  ! advec_error(2,:,:) = error value for pressure
  ! advec_error(3,:,:) = error value for entropy
  ! advec_error(:,1,:) = error computed using gauss quadrature
  ! advec_error(:,2,:) = error computed using radau quadrature
  ! advec_error(:,3,:) = error computed using lobatto quadrature
  ! advec_error(:,:,1) = L1 norm of total error for all solution points
  ! advec_error(:,:,2) = L2 norm of total error for all solution points
  ! advec_error(:,:,3) = infinity norm of total error for all solution points
  !                      i.e., maximum value of error for all solution points
  ! advec_error(:,:,4) = ID of grid cell containing maximum value of error
  ! advec_error(:,:,5) = ID of solution point with maximum value of error
  !
  !.. Use Statements ..
  use geovar,  only : nr,neqns,ncell,pts_in_cell_ptr,n_sp1de,n_sp2dt,n_sp2dq
  use geovar,  only : xyz,xyz_radau,xyz_gauss,xyz_lobatto
  use geovar,  only : numer2radau,numer2gauss,numer2lobatto
  use geovar,  only : weights2d,radau_wts_2d,gauss_wts_2d,lobatto_wts_2d
  use ovar,    only : total_volume,pv_in,gam,itestcase,itcur,advec_error, &
                      Period_Timesteps,Grid_Min_Y,Grid_Length_Y, &
                      Grid_Min_X,Grid_Length_X,vxref,vyref,aref, &
                      Density_Maximum_Residual,vortex_has_collapsed
  use metvar,  only : jacobian
  use flowvar, only : usp
  use parallel_mod, only : cell_map,cell_map_ptr
  !
  !.. Local Scalars ..
  integer   :: n,n1,n2,np,k,l,w1,m,nc1
  real(wp) :: rho_ratio,pres_ratio,c1,wg,wr,wl
  real(wp) :: delx,dely,xoffset,yoffset
  real(wp) :: sn,se
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr,1:n_sp2dq) :: xyz_r
  real(wp), dimension(1:nr,1:n_sp2dq) :: xyz_g
  real(wp), dimension(1:nr,1:n_sp2dq) :: xyz_l
  !
  real(wp), dimension(1:neqns,1:n_sp2dq) :: vre
  real(wp), dimension(1:neqns,1:n_sp2dq) :: vge
  real(wp), dimension(1:neqns,1:n_sp2dq) :: vle
  !
  real(wp), dimension(1:neqns,1:n_sp2dq) :: vrn
  real(wp), dimension(1:neqns,1:n_sp2dq) :: vgn
  real(wp), dimension(1:neqns,1:n_sp2dq) :: vln
  !
  real(wp), dimension(1:3,1:3) :: ett
  real(wp), dimension(1:3,1:3,1:5) :: error
  real(wp), dimension(1:3,1:3,1:5,1:ncpu) :: all_error
  !
continue
  !
  ! If running one of the vortex transport problems, compute the change
  ! in position of the vortex since the start of the current period.
  !
  if (any(itestcase == Transport_Problems)) call vortex_position(delx,dely)
  !
  error(:,:,:) = zero  ! intialize everything to zero
  error(:,:,3) = -eps6 ! make this slightly negative in case all error is zero
  !
  c1 = deltafun(itestcase==Density_Transport)
  !
  nc1 = cell_map_ptr(mypnum+1)
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1
    n2 = pts_in_cell_ptr(n+1)
    np = n2-n1+1
    !
    w1 = merge(1,2,np==n_sp2dt)
    !
    ! Copy the numerical solution to a local array
    !
    vrn(:,1:np) = usp(:,n1:n2)
    vgn(:,1:np) = usp(:,n1:n2)
    vln(:,1:np) = usp(:,n1:n2)
    !
    ! Get the solution point coordinates for this cell
    ! using the three different quadrature points
    !
    xyz_r(:,1:np) = xyz_radau(:,n1:n2)
    xyz_g(:,1:np) = xyz_gauss(:,n1:n2)
    xyz_l(:,1:np) = xyz_lobatto(:,n1:n2)
    !
    ! Get the exact solution for this cell for each of the three quadratures
    !
    do k = 1,np ! radau
      !
      xoffset = merge(delx,delx-Grid_Length_X,xyz_r(1,k)-delx >= Grid_Min_X)
      yoffset = merge(dely,dely-Grid_Length_Y,xyz_r(2,k)-dely >= Grid_Min_Y)
      !
      vre(:,k) = vortex2D( xyz_r(1,k)-xoffset , xyz_r(2,k)-yoffset )
      !
    end do
    !
    do k = 1,np ! gauss
      !
      xoffset = merge(delx,delx-Grid_Length_X,xyz_g(1,k)-delx >= Grid_Min_X)
      yoffset = merge(dely,dely-Grid_Length_Y,xyz_g(2,k)-dely >= Grid_Min_Y)
      !
      vge(:,k) = vortex2D( xyz_g(1,k)-xoffset , xyz_g(2,k)-yoffset )
      !
    end do
    !
    do k = 1,np ! lobatto
      !
      xoffset = merge(delx,delx-Grid_Length_X,xyz_l(1,k)-delx >= Grid_Min_X)
      yoffset = merge(dely,dely-Grid_Length_Y,xyz_l(2,k)-dely >= Grid_Min_Y)
      !
      vle(:,k) = vortex2D( xyz_l(1,k)-xoffset , xyz_l(2,k)-yoffset )
      !
    end do
    !
    ! Project the exact solutions to p_order
    !
    call project_to_porder(1,[0,np],vge(1:neqns,1:np),2)
    call project_to_porder(1,[0,np],vle(1:neqns,1:np),3)
    call project_to_porder(1,[0,np],vre(1:neqns,1:np),4)
    !
    ! Interpolate the numerical solution to the solution
    ! nodes of the three different quadratures
    !
    do m = 1,neqns
      vrn(m,1:np) = matmul( numer2radau , vrn(m,1:np) )
      vgn(m,1:np) = matmul( numer2gauss , vgn(m,1:np) )
      vln(m,1:np) = matmul( numer2lobatto , vln(m,1:np) )
    end do
    !
    ! Get the primitive variables for the numerical and analytical solutions
    !
    vre = usp2v(vre)
    vge = usp2v(vge)
    vle = usp2v(vle)
    !
    vrn = usp2v(vrn)
    vgn = usp2v(vgn)
    vln = usp2v(vln)
    !
    ! Now loop through the solution points of this cell to compute
    ! the error between the simulation and analytical solutions.
    !
    do l = n1,n2
      !
      k = l-n1+1
      !
      ! Gauss error
      !
     !rho_ratio  = vgn(    1,k)/pv_in(    1)
     !pres_ratio = vgn(neqns,k)/pv_in(neqns)
      sn = entropy_pv_sp( vgn(:,k) )
      se = entropy_pv_sp( vge(:,k) )
      !
      ett(1,1) = abs( vgn(    1,k) - vge(    1,k) )
      ett(2,1) = abs( vgn(neqns,k) - vge(neqns,k) )
     !ett(3,1) = log( (c1+pres_ratio) / (rho_ratio**gam) )
      ett(3,1) = sn - se
      !
      ! Radau error
      !
     !rho_ratio  = vrn(    1,k)/pv_in(    1)
     !pres_ratio = vrn(neqns,k)/pv_in(neqns)
      sn = entropy_pv_sp( vrn(:,k) )
      se = entropy_pv_sp( vre(:,k) )
      !
      ett(1,2) = abs( vrn(    1,k) - vre(    1,k) )
      ett(2,2) = abs( vrn(neqns,k) - vre(neqns,k) )
     !ett(3,2) = log( (c1+pres_ratio) / (rho_ratio**gam) )
      ett(3,2) = sn - se
      !
      ! Lobatto error
      !
     !rho_ratio  = vln(    1,k)/pv_in(    1)
     !pres_ratio = vln(neqns,k)/pv_in(neqns)
      sn = entropy_pv_sp( vln(:,k) )
      se = entropy_pv_sp( vle(:,k) )
      !
      ett(1,3) = abs( vln(    1,k) - vle(    1,k) )
      ett(2,3) = abs( vln(neqns,k) - vle(neqns,k) )
     !ett(3,3) = log( (c1+pres_ratio) / (rho_ratio**gam) )
      ett(3,3) = sn - se
      !
      ! Now add the error for all variables at the current solution
      ! point to the running integral sum for the L1 and L2 error norms
      !
      wg = jacobian(l)*gauss_wts_2d(k)
      wr = jacobian(l)*radau_wts_2d(k)
      wl = jacobian(l)*lobatto_wts_2d(k)
      !
      error(:,1,1) = error(:,1,1) + wg*ett(:,1)
      error(:,1,2) = error(:,1,2) + wg*ett(:,1)**2
      !
      error(:,2,1) = error(:,2,1) + wr*ett(:,2)
      error(:,2,2) = error(:,2,2) + wr*ett(:,2)**2
      !
      error(:,3,1) = error(:,3,1) + wl*ett(:,3)
      error(:,3,2) = error(:,3,2) + wl*ett(:,3)**2
      !
      ! If any of these errors are greater than the previous maximum
      ! values, save the error value and record this location
      !
      where (abs(ett(:,:)) > error(:,:,3))
        error(:,:,3) = abs(ett(:,:))
        error(:,:,4) = real(cell_map(nc1+n),kind=wp)
        error(:,:,5) = real(k,kind=wp)
      end where
      !
    end do
    !
  end do
  !
  ! Collect the error information from all processors
  !
  if (ncpu > 1) then
    call mpi_allgather(error,    size(error),mpi_flttyp, &
                       all_error,size(error),mpi_flttyp,mpi_comm_world,mpierr)
  else
    all_error(:,:,:,1) = error(:,:,:)
  end if
  !
  ! Finish computing the global error across all processors
  !
  do n = 1,3
    do m = 1,3
      advec_error(m,n,1) =       sum(all_error(m,n,1,:)) / total_volume
      advec_error(m,n,2) = sqrt( sum(all_error(m,n,2,:)) / total_volume )
      advec_error(m,n,3) = maxval(all_error(m,n,3,:))
      advec_error(m,n,4) = all_error(m,n,4,maxloc(all_error(m,n,3,:),dim=1))
      advec_error(m,n,5) = all_error(m,n,5,maxloc(all_error(m,n,3,:),dim=1))
    end do
  end do
  !
  ! Check for Vortex Collapse
  !
  if (advec_error(1,1,2) > zero) then
    if (log10(advec_error(1,1,2)) > Density_Maximum_Residual) then
      Density_Maximum_Residual = log10(advec_error(1,1,2))
    end if
    vortex_has_collapsed = Density_Maximum_Residual > -two
  end if
  !
end subroutine compute_final_advec_error
!
!###############################################################################
!
subroutine test_error()
  !
  ! Compute the final error and store it in the final_error.
  !
  ! final_error(1,:,:) = error value for density
  ! final_error(2,:,:) = error value for pressure
  ! final_error(3,:,:) = error value for entropy
  ! final_error(:,1,:) = error computed using gauss quadrature
  ! final_error(:,2,:) = error computed using radau quadrature
  ! final_error(:,3,:) = error computed using lobatto quadrature
  ! final_error(:,:,1) = L1 norm of total error for all solution points
  ! final_error(:,:,2) = L2 norm of total error for all solution points
  ! final_error(:,:,3) = infinity norm of total error for all solution points
  !                      i.e., maximum value of error for all solution points
  ! final_error(:,:,4) = ID of grid cell containing maximum value of error
  ! final_error(:,:,5) = ID of solution point with maximum value of error
  !
  !.. Use Statements ..
  use geovar,  only : ncell,pts_in_cell_ptr
  use geovar,  only : radau_wts_2d,gauss_wts_2d,lobatto_wts_2d
  use ovar,    only : total_volume
  use metvar,  only : jacobian
  !
  !.. Local Scalars ..
  integer   :: n,n1,n2,np,k,l
  real(wp) :: ett,wg,wr,wl,j,e
  !
  !.. Local Arrays ..
  real(wp), dimension(1:3,1:2) :: error
  real(wp), dimension(1:3,1:4) :: final_error
  real(wp), dimension(1:3,1:2,1:ncpu) :: all_error
  !
continue
  !
  error(:,:) = zero
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1
    n2 = pts_in_cell_ptr(n+1)
    np = n2-n1+1
    !
    ! Now loop through the solution points of this cell to compute
    ! the error between the simulation and analytical solutions.
    !
    do l = n1,n2
      !
      k = l-n1+1
      !
      ett = rndoff
      !
     !e = scale(ett,-exponent(eps6))
     !j = scale(jacobian(l),-exponent(eps6))
     !wg = scale(gauss_wts_2d(k),-exponent(eps6))
     !wr = scale(radau_wts_2d(k),-exponent(eps6))
     !wl = scale(lobatto_wts_2d(k),-exponent(eps6))
      !
      ! Now add the error for all variables at the current solution
      ! point to the running integral sum for the L1 and L2 error norms
      !
      wg = jacobian(l)*gauss_wts_2d(k)
      wr = jacobian(l)*radau_wts_2d(k)
      wl = jacobian(l)*lobatto_wts_2d(k)
      !
      error(1,1) = error(1,1) + wg*ett
      error(1,2) = error(1,2) + wg*ett**2
      !
      error(2,1) = error(2,1) + wr*ett
      error(2,2) = error(2,2) + wr*ett**2
      !
      error(3,1) = error(3,1) + wl*ett
      error(3,2) = error(3,2) + wl*ett**2
      !
     !error(1,1) = error(1,1) + j*wg*e
     !error(1,2) = error(1,2) + j*wg*e**2
     !!
     !error(2,1) = error(2,1) + j*wr*e
     !error(2,2) = error(2,2) + j*wr*e**2
     !!
     !error(3,1) = error(3,1) + j*wl*e
     !error(3,2) = error(3,2) + j*wl*e**2
     !!
    end do
    !
  end do
  !
  ! Collect the error information from all processors
  !
  if (ncpu > 1) then
    call mpi_allgather(error,    size(error),mpi_flttyp, &
                       all_error,size(error),mpi_flttyp,mpi_comm_world,mpierr)
  else
    all_error(:,:,1) = error(:,:)
  end if
  !
  write (*,*)
  write (*,*) all_error(:,:,:)
  write (*,*)
  !
  ! Finish computing the global error across all processors
  !
  do n = 1,3
    final_error(n,1) = sum(all_error(n,1,:))
    final_error(n,2) = sum(all_error(n,2,:))
   !final_error(n,1) = scale( sum(all_error(n,1,:)) , exponent(eps6) )
   !final_error(n,2) = scale( sum(all_error(n,2,:)) , exponent(eps6) )
    final_error(n,3) = final_error(n,1) / total_volume
    final_error(n,4) = sqrt( final_error(n,2) / total_volume )
  end do
  !
  call mpi_barrier(mpi_comm_world,mpierr)
  !
  if (mypnum == 0) then
    write (*,1) final_error(:,:)
  end if
  !
  call mpi_barrier(mpi_comm_world,mpierr)
  !
  call finalize_mpi_environment
  !
  1 format (/,30x,"GAUSS",18x,"RADAU",17x,"LOBATTO",/, &
              " SUM OF ERROR     : ",3es23.15,/, &
              " SUM OF ERROR**2  : ",3es23.15,/, &
              " L1 NORM OF ERROR : ",3es23.15,/, &
              " L2 NORM OF ERROR : ",3es23.15,///)
  !
end subroutine test_error
!
!###############################################################################
!
subroutine rgeterr()
  !
  !.. Use Statements ..
  use geovar,  only : neqns,ncell,pts_in_cell_ptr,n_sp2dt,weights2d,xyz
  use geovar,  only : global_cell_offset,global_xyz
  use ovar,    only : rlinf,rlone,rltwo,nloc,kloc,total_volume,gam,cv_in, &
                      nentmax,kentmax,xentmax,yentmax,entmax,entL2,itestcase
  use ovar,    only : itcur,d_t,time,rl2mx
  use ovar,    only : max_wall_time,wall_time_left
  use ovar,    only : beg_wall_time,wall_time
  use ovar,    only : walltime_expiring,convergence_reached
  use ovar,    only : convergence_order_max_res,convergence_order_abs_res
  use metvar,  only : jacobian
  use flowvar, only : usp,uoldsp
  use parallel_mod, only : cell_map,cell_map_ptr
  !
  !.. Local Scalars ..
  integer   :: n,n1,n2,k,w1,nc1
  integer   :: hrs_left,min_left,sec_left
  real(wp) :: rho,pres,c1,sn,se
  real(wp) :: next_time_est
  real(wp) :: wall_time_last
  !
  !.. Local Arrays ..
  real(wp), dimension(1:neqns) :: vv
  real(wp), dimension(1:2) :: ett
  real(wp), dimension(1:2,1:7) :: error
  real(wp), dimension(1:2,1:7,1:ncpu) :: all_error
  !
  !.. Local Parameters ..
  real(wp), parameter :: max_sec = real( 1000*3600-1 , kind=wp )
  !
  ! error(:,1) = L1 norm of total error for all solution points
  ! error(:,2) = L2 norm of total error for all solution points
  ! error(:,3) = maximum value of error for all solution points
  ! error(:,4) = ID of grid cell containing maximum value of error
  ! error(:,5) = ID of solution point with maximum value of error
  ! error(:,6) = x-coordinate of solution point with maximum value of error
  ! error(:,7) = y-coordinate of solution point with maximum value of error
  ! error(1,:) = error values for density
  ! error(2,:) = error values for entropy
  !
continue
  !
  error(:,:) = zero  ! intialize everything to zero
  error(:,3) = -eps6 ! make this slightly negative in case all error is zero
  !
  c1 = deltafun(itestcase==Density_Transport)
  !
  nc1 = cell_map_ptr(mypnum+1)
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1
    n2 = pts_in_cell_ptr(n+1)
    !
    w1 = merge(1,2,n2-n1+1==n_sp2dt)
    !
    do k = n1,n2
      !
      ! Get the change in density for this solution point
      !
      if (itestcase == Laminar_BndLyr) then
        ett(1) = abs( usp(2,k) - uoldsp(2,k) )
      else
        ett(1) = abs( usp(1,k) - uoldsp(1,k) )
      end if
      !
      ! Compute the entropy at this solution point
      !
      rho = usp(1,k)
      pres = pressure_cv_sp( usp(:,k) )
      !
      if (any([rho,pres] <= zero)) write (iout,111) n,k,rho,pres
      111 format ("cell = ",i0,"; SP = ",i0,"; ",*(1x,es14.6))
      !
      sn = entropy_cv_sp( usp(:,k) )
      se = entropy_cv_sp( cv_in(:) )
      !
      ett(2) = sn - se
      !
      ! Get the norms of density and entropy
      !
      error(:,1) = error(:,1) + jacobian(k)*weights2d(k-n1+1,w1)*ett(:)
      error(:,2) = error(:,2) + jacobian(k)*weights2d(k-n1+1,w1)*ett(:)*ett(:)
      !
      where (abs(ett(:)) > error(:,3))
        error(:,3) = abs(ett)
        error(:,4) = real(cell_map(nc1+n),kind=wp)
        error(:,5) = real(k-n1+1,kind=wp)
        error(:,6) = xyz(1,k)
        error(:,7) = xyz(2,k)
      end where
      !
    end do
  end do
  !
  ! Collect the error information from all processors
  !
  if (ncpu > 1) then
    call mpi_allgather(error,    size(error),mpi_flttyp, &
                       all_error,size(error),mpi_flttyp,mpi_comm_world,mpierr)
  else
    all_error(:,:,1) = error(:,:)
  end if
  !
  rlinf = maxval(all_error(1,3,:))
  nloc  = nint( all_error(1,4,maxloc(all_error(1,3,:),dim=1)) )
  kloc  = nint( all_error(1,5,maxloc(all_error(1,3,:),dim=1)) )
  !
  entmax  = maxval(all_error(2,3,:))
  nentmax = nint( all_error(2,4,maxloc(all_error(2,3,:),dim=1)) )
  kentmax = nint( all_error(2,5,maxloc(all_error(2,3,:),dim=1)) )
  xentmax =       all_error(2,6,maxloc(all_error(2,3,:),dim=1))
  yentmax =       all_error(2,7,maxloc(all_error(2,3,:),dim=1))
  !
  rlone =       sum(all_error(1,1,:)) / total_volume
  rltwo = sqrt( sum(all_error(1,2,:)) / total_volume )
  !
  rl2mx = max( rltwo , rl2mx )
  !
  if (any(itestcase == Transport_Problems)) then
    entL2 = sqrt( sum(all_error(2,2,:)) / total_volume )
  else
    entL2 = log10(rl2mx) - log10(rltwo)
    convergence_reached = (entL2 <= convergence_order_max_res) .or. &
                          (log10(rltwo) <= convergence_order_abs_res)
  end if
  !
  ! Output the results for the current time step
  !
  if (mypnum == 0) then
    !
#ifdef PBS_ENV
    !
    wall_time_last = wall_time_left
    wall_time = mpi_wtime() - beg_wall_time
    wall_time_left = min( max_wall_time - wall_time , max_sec )
    !
    next_time_est = (wall_time_last - wall_time_left) * (one+eps1)
    walltime_expiring = (wall_time_left-next_time_est <= zero)
    !
    sec_left = abs(nint(wall_time_left))
    !
    hrs_left = sec_left/3600
    min_left = (sec_left - hrs_left*3600)/60
    sec_left = sec_left - hrs_left*3600 - min_left*60
    ! Make hrs_left negative if wall_time_left happens to be negative
    hrs_left = sign(hrs_left,int(wall_time_left))
    !
    write (iout,iter_out_fmt) itcur,d_t,time,rltwo,entL2,nloc, &
                              global_cell_offset(nloc)+kloc, &
                              global_xyz(1,global_cell_offset(nloc)+kloc), &
                              global_xyz(2,global_cell_offset(nloc)+kloc), &
                              hrs_left,min_left,sec_left
    !
#else
    !
    walltime_expiring = fals
    !
    write (iout,iter_out_fmt) itcur,d_t,time,rltwo,entL2,nloc, &
                              global_cell_offset(nloc)+kloc, &
                              global_xyz(1,global_cell_offset(nloc)+kloc), &
                              global_xyz(2,global_cell_offset(nloc)+kloc)
    !
    !
#endif
    !
    if (convergence_reached) then
      if (entL2 <= convergence_order_max_res) write (iout,3)
      if (log10(rltwo) <= convergence_order_abs_res) write (iout,4)
    end if
    !
    if (walltime_expiring) write (iout,2)
    !
  end if
  !
  ! Broadcast walltime_expiring to all processors
  !
  call mpi_bcast(walltime_expiring,1,mpi_logical,0,mpi_comm_world,mpierr)
  !

  1 format (i9,4(1x,es11.4),2i8,2es9.1,1x,i3,2(":",i2.2))
  2 format (/,"*******************************************************",/, &
              " THE SIMULATION WALL TIME REQUESTED IS EXPIRING SOON!",/, &
              " STOPPING EXECUTION AND DUMPING THE FINAL RESULTS!",/, &
              "*******************************************************",/)
  3 format (/,"*******************************************************",/, &
              " THE CONVERGENCE CRITERION FOR THE ORDERS OF",/, &
              " MAGNITUDE REDUCTION OF THE RESIDUAL ERROR",/, &
              " RELATIVE TO THE MAXIMUM RESIDUAL ERROR HAS BEEN MET!",/, &
              "*******************************************************",/)
  4 format (/,"*******************************************************",/, &
              " THE CONVERGENCE CRITERION FOR THE ORDER OF MAGNITUDE",/, &
              " REDUCTION OF THE ABSOLUTE RESIDUAL ERROR HAS BEEN MET!",/, &
              "*******************************************************",/)
  !
end subroutine rgeterr
!
!###############################################################################
!
subroutine init_generic_mod
  !
  !.. Use Statements ..
  use geovar, only : n_global_cell,n_global_solpts
  use ovar,   only : itcur,num_ts,rl2mx,itestcase
  !
  !.. Local Scalars ..
  integer :: final_iter_digits
  integer :: max_cell_digits
  integer :: max_solpt_digits
  character(len=100) :: char_num
  character(len=1000) :: iter_header
  !
continue
  !
  rl2mx = -huge(zero)
  !
  write (char_num,1) itcur + num_ts
  final_iter_digits = max(4,len_trim(adjustl(char_num)))
  !
  write (char_num,1) n_global_cell
  max_cell_digits = max(4,len_trim(adjustl(char_num)))
  !
  write (char_num,1) n_global_solpts
  max_solpt_digits = max(5,len_trim(adjustl(char_num)))
  !
  write (iter_out_fmt,100) final_iter_digits,max_cell_digits,max_solpt_digits
  !
  ! Add on the remaining simulation time to the output if
  ! running in the PBS environment
  !
#ifdef PBS_ENV
  iter_out_fmt = iter_out_fmt(1:len_trim(iter_out_fmt)) // ',1x,i3,2(":",i2.2)'
#endif
  !
  ! Put the last ')' onto iter_out_fmt
  !
  iter_out_fmt = iter_out_fmt(1:len_trim(iter_out_fmt)) // ")"
  !
  ! Now get the format for the header
  !
  if (any(itestcase == Transport_Problems)) then
    write (iter_header,200) final_iter_digits-3, &
                            max_cell_digits-2, &
                            max_solpt_digits-3
  else if (itestcase == Laminar_BndLyr) then
    write (iter_header,201) final_iter_digits-3, &
                            max_cell_digits-2, &
                            max_solpt_digits-3
  else
    write (iter_header,202) final_iter_digits-3, &
                            max_cell_digits-2, &
                            max_solpt_digits-3
  end if
  !
#ifdef PBS_ENV
  iter_header = iter_header(1:len_trim(iter_header)) // ",4x,'Time Left'"
#endif
  !
  ! Put the last ')' onto iter_header
  !
  iter_header = iter_header(1:len_trim(iter_header)) // ")"
  !
  ! Write the header for the terminal output
  !
  if (mypnum == 0) write (iout,iter_header)
  !
  ! Format Statements
  !
  1 format (i0)
  100 format ("(1x,i",i0,",4es12.4,2x,i",i0,",2x,i",i0,",2es9.1")
  200 format ("(/,",i0,"x,'iter',5x,'d_t',8x,'Time',6x,'L2-Density',2x,", &
              "'L2-Entropy',",i0,"x,'Cell',",i0,"x,'Solpt',5x,'x',8x,'y'")
  201 format ("(/,",i0,"x,'iter',5x,'d_t',8x,'Time',6x,'L2-X-Mmntm',2x,", &
              "'OrdersDrop',",i0,"x,'Cell',",i0,"x,'Solpt',5x,'x',8x,'y'")
  202 format ("(/,",i0,"x,'iter',5x,'d_t',8x,'Time',6x,'L2-Density',2x,", &
              "'OrdersDrop',",i0,"x,'Cell',",i0,"x,'Solpt',5x,'x',8x,'y'")
  !
end subroutine init_generic_mod
!
!###############################################################################
!
subroutine getmaxent
  !
  !... get max entropy
  !
  !.. Use Statements ..
  use geovar,  only : neqns,ncell,pts_in_cell_ptr,n_sp2dt,weights2d,xyz
  use ovar,    only : nentmax,kentmax,entmax,gam,cv_in, &
                      xentmax,yentmax,entL2,total_volume,itestcase
  use metvar,  only : jacobian
  use flowvar, only : usp
  use parallel_mod, only : cell_map,cell_map_ptr
  !
  !.. Local Scalars ..
  integer  :: n,n1,n2,k,w1,nc1
  real(wp) :: ett,c1,sn,se
  !
  !.. Local Arrays ..
  real(wp), dimension(1:6) :: error
  real(wp), dimension(1:6,1:ncpu) :: all_error
  !
continue
  !
  error(:) = zero  ! intialize everything to zero
  error(2) = -eps6 ! make this slightly negative in case all error is zero
  !
  c1 = deltafun(itestcase==Density_Transport)
  !
  nc1 = cell_map_ptr(mypnum+1)
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1
    n2 = pts_in_cell_ptr(n+1)
    !
    w1 = merge(1,2,n2-n1+1==n_sp2dt)
    !
    do k = n1,n2
      !
      sn = entropy_cv_sp( usp(:,k) )
      se = entropy_cv_sp( cv_in(:) )
      !
      ett = sn - se
      !
      ! L2 norm of entropy
      !
      error(1) = error(1) + jacobian(k)*weights2d(k-n1+1,w1)*ett*ett
      !
      if (abs(ett) > error(2)) then
        error(2) = abs(ett)
        error(3) = real(cell_map(nc1+n),kind=wp)
        error(4) = real(k-n1+1,kind=wp)
        error(5) = xyz(1,k)
        error(6) = xyz(2,k)
      end if
      !
    end do
    !
  end do
  !
  ! Collect the error information from all processors
  !
  if (ncpu > 1) then
    call mpi_allgather(error,    size(error),mpi_flttyp, &
                       all_error,size(error),mpi_flttyp,mpi_comm_world,mpierr)
  else
    all_error(:,1) = error(:)
  end if
  !
  entmax  = maxval(all_error(2,:))
  nentmax = nint( all_error(3,maxloc(all_error(2,:),dim=1)) )
  kentmax = nint( all_error(4,maxloc(all_error(2,:),dim=1)) )
  xentmax =       all_error(5,maxloc(all_error(2,:),dim=1))
  yentmax =       all_error(6,maxloc(all_error(2,:),dim=1))
  !
  entL2 = sqrt( sum(all_error(1,:)) / total_volume )
  !
end subroutine getmaxent
!
!###############################################################################
!
subroutine getbc
  !
  !...at inflow,
  !
  !...kbc = 0 for solid boundary;
  !...      1 subsonic inflow (extrapolate p from interior);
  !...      2 supersonic inflow (specify all 4 quantities) ;
  !...      3 subsonic outflow (specify pressure) ;
  !...      4 supersonic outflow (extrapolate all 4 quantities) ;
  !...      7 periodic bc ;
  !
  !... get boundary conditions, i.e., usp, faceusp, and facedusp
  !... just outside the domain
  !
  !.. Use Statements ..
  use geovar
  use ovar,         only : gam,cv_in,pv_in
  use ovar,         only : flux_divergence
  use metvar,       only : dx,dy,geofa
  use flowvar,      only : usp,faceusp,facedusp
  use parallel_mod, only : nonblocking_exchange
  !
  !.. Local Scalars ..
  integer   :: gc,g1,g2,gp,gw,gkb,gk,k,p1,p2,ip1,nf
  integer   :: hc,h1,h2,hp,hw,hkb,hk,m,n1,n2,ip2,kg
  real(wp) :: ai,ao,vni,vno,rp,rm,uc,ac,vti,vto
  real(wp) :: u,v,ut,un,gm1,mach,utc,unc,sc
  character(len=8) :: flow
  character(len=6) :: side
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: rnrm
  real(wp), dimension(1:nr) :: rtan
  real(wp), dimension(1:neqns) :: varh
  real(wp), dimension(1:n_sp2dq,1:neqns) :: tusp
  real(wp), dimension(1:n_sp2dq,1:4*n_fp1de,1:2) :: interp_mat
  !
continue
  !
  gm1 = gam - one
  !
  if (n_fp1de /= n_sp1de) then
    interp_mat(:,:,:) = interp_s2g(:,:,:)
  else
    interp_mat(:,:,:) = interp_s2f(:,:,:)
  end if
  !
  ! First exchange information between adjacent partitions
  !
  if (ncpu > 1 .and. all(flux_divergence/=[TestQuadrature,TestViscous])) then
    call nonblocking_exchange(usp)
  end if
  !
  do nf = 1,nfbnd
    !
    hc = cells_on_face(1,nf) ! left/host   cell on boundary face
    gc = cells_on_face(2,nf) ! right/ghost cell on boundary face
    !
    hkb = cells_on_face(5,nf) ! beginning flxpt on current face for host  cell
    gkb = cells_on_face(6,nf) ! beginning flxpt on current face for ghost cell
    !
   !hf1 = flxpts_of_cell_ptr(hc)+1 ! beginning flxpt index for host  cell
   !hf2 = flxpts_of_cell_ptr(hc+1) ! ending    flxpt index for host  cell
   !!
   !gf1 = flxpts_of_cell_ptr(gc)+1 ! beginning flxpt index for ghost cell
   !gf2 = flxpts_of_cell_ptr(gc+1) ! ending    flxpt index for ghost cell
    !
    h1 = pts_in_cell_ptr(hc)+1 ! beginning index for sol-pts in host cell
    h2 = pts_in_cell_ptr(hc+1) ! ending    index for sol-pts in host cell
    !
    g1 = pts_in_cell_ptr(gc)+1 ! beginning index for sol-pts in ghost cell
    g2 = pts_in_cell_ptr(gc+1) ! ending    index for sol-pts in ghost cell
    !
    hp = h2-h1+1
    gp = g2-g1+1
    !
    hw = merge(1,2,hp == n_sp2dt)
    gw = merge(1,2,gp == n_sp2dt)
    !
    ! Fill the ghost cell depending on boundary condition
    !
    if (bface(1,nf) == bc_freestream) then
      !
      ! #####################################################
      ! #####################################################
      ! ##  Free stream boundary condition
      ! #####################################################
      ! #####################################################
      !
      do k = g1,g2
        usp(:,k) = cv_in(:)
      end do
      !
      do k = 1,n_fp1de
        !
        faceusp(:,k,2,nf) = cv_in(:)
        !
        facedusp(:,:,k,2,nf) = zero
        !
      end do
      !
    else if (bface(1,nf) == bc_characteristic) then
      !
      ! #####################################################
      ! #####################################################
      ! ##  Characteristic boundary condition
      ! #####################################################
      ! #####################################################
      !
      n1 = nodes_of_cell_ptr(hc)+1
      !
      ! Store the transpose of usp for the host cell to reduce cache misses
      !
      do k = h1,h2
        tusp(k-h1+1,:) = usp(:,k)
      end do
      !
      do k = g1,g2
        usp(:,k) = cv_in(:)
      end do
      !
      ! Get the unit tangent vector for this face
      !
      ip1 = nodes_on_face(1,nf)
      ip2 = nodes_on_face(2,nf)
      !
      rtan(:) = -unit_vector( xyz_nodes(:,ip2)-xyz_nodes(:,ip1) )
      !
      do k = 1,n_fp1de
        !
        ! Inward pointing unit normal at this flux point
        !
        rnrm(:) = -unit_vector( geofa(:,k,nf) )
        !
        ! Unit tangent at this flux point computed from unit normal
        ! NOTE: Usually t = (tx,ty) = (ny,-nx) but we computed the inward
        !       pointing unit normal so here (tx,ty) = (-ny,nx).
        !       
       !rtan(1) = -rnrm(2)
       !rtan(2) =  rnrm(1)
       !write (101,2) nf,k,rnrm(:),rtan(:)
       !2 format ("nf=",i0,"; k=",i0,"; rnrm(:)=(",es13.6,",",es13.6,")", &
       !                             "; rtan(:)=(",es13.6,",",es13.6,")")
        !
        hk = hkb-1 + k ! index for current host  cell flux point
        gk = gkb+1 - k ! index for current ghost cell flux point
        kg = g1-1 + gk ! index for current ghost cell solution point
        !
        ! Interpolate the conservative variables from the host cell
        ! solution points to the current host cell flux point.
        !
        do m = 1,neqns
          varh(m) = dot(interp_mat(1:hp,hk,hw),tusp(1:hp,m))
        end do
        !
        ! Now convert the current host cell flux point
        ! conservative variables to primitive variables
        !
        varh(:) = usp2v_sp( varh(:) )
        !
        ! Compute the speed of sound in the host and ghost cells
        !
        ai = sqrt(max(gam*varh(neqns)/varh(1),rndoff))
        ao = sqrt(max(gam*pv_in(neqns)/pv_in(1),rndoff))
        !
        ! Compute the interior and free stream normal velocity
        !
        vni = dot(rnrm,varh(2:neqns-1))
        vno = dot(rnrm,pv_in(2:neqns-1))
        !
        ! Compute the interior and free stream tangent velocity
        !
        vti = dot(rtan,varh(2:neqns-1))
        vto = dot(rtan,pv_in(2:neqns-1))
        !
        ! Apply boundary conditions based on the
        ! normal Mach number to the boundary face
        !
        mach = vni/ai
        !
        if (mach >= one) then ! Supersonic inflow
          !
          ! Set ghost cell flux point to free stream conditions
          !
          usp(:,kg) = cv_in(:) ! sup-in
          !
          faceusp(:,k,2,nf) = usp(:,kg)
          !
          facedusp(:,:,k,2,nf) = zero
          !
        else if (mach <= -one) then ! Supersonic outflow
          !
          ! Extrapolate interior to ghost cell flux point
          !
          usp(:,kg) = vsp2u_sp( varh(:) ) ! sup-out
          !
          faceusp(:,k,2,nf) = usp(:,kg)
          !
          facedusp(:,:,k,2,nf) = zero
          !
        else ! Subsonic inflow/outflow using Riemann invariants
          !
          if (mach > zero) then ! this is an inflow boundary
            !
            rp = vno + two*ao/gm1
            rm = vni - two*ai/gm1
            unc = half*(rp+rm)
            utc = vto
            ac = one4*gm1*(rp-rm)
            sc = pv_in(neqns)/(pv_in(1)**gam)
            !
          else ! this is an outflow boundary
            !
            rp = -abs(vno) + two*ao/gm1
            rm = vni - two*ai/gm1
            unc = half*(rp+rm)
            utc = vti
            ac = one4*gm1*(rp-rm)
            sc = varh(neqns)/(varh(1)**gam)
            !
          end if
          !
          varh(1) = ((ac*ac)/(gam*sc))**(one/gm1)
          varh(neqns) = sc*(varh(1)**gam)
          varh(2:neqns-1) = unc*rnrm(:) + utc*rtan(:)
          !
          usp(:,kg) = vsp2u_sp( varh(:) )
          !
          faceusp(:,k,2,nf) = usp(:,kg)
          !
          facedusp(:,:,k,2,nf) = zero
          !
        end if
        !
      end do
      !
    else if (bface(1,nf) == bc_generic_inflow) then
      !
      ! #####################################################
      ! #####################################################
      ! ##  Generic inflow boundary condition
      ! ##    - Use subsonic or supersonic inflow boundary
      ! ##      conditions based on the Mach number of the
      ! ##      normal inflow velocity
      ! #####################################################
      ! #####################################################
      !
      ! Store the transpose of usp for the host cell to reduce cache misses
      !
      do k = h1,h2
        tusp(k-h1+1,:) = usp(:,k)
      end do
      !
      ! Now loop through the boundary face flux points to get the boundary
      ! condition solution at each flux point
      !
      do k = 1,n_fp1de
        !
        ! Unit normal and tangent at this flux point
        !
        rnrm(:) = unit_vector( geofa(:,k,nf) )
        !
        hk = hkb-1 + k ! index for current host  cell flux point
        gk = gkb+1 - k ! index for current ghost cell flux point
        kg = g1-1 + gk ! index for current ghost cell solution point
        !
        ! Interpolate the conservative variables for the host cell
        ! solution points to the current host cell flux point.
        !
        do m = 1,neqns
          varh(m) = dot(interp_mat(1:hp,hk,hw),tusp(1:hp,m))
        end do
        !
        ! Now convert the current host cell flux point
        ! conservative variables to primitive variables
        !
        varh(:) = usp2v_sp( varh(:) )
        !
        ! Compute the speed of sound at this flux point
        !
        ai = sqrt(max(gam*varh(neqns)/varh(1),rndoff))
        !
        ! Compute the Mach number of the normal velocity at this flux point
        !
        mach = dot( varh(2:neqns-1) , rnrm(1:nr) ) / ai
        !
        if (mach > -one .and. mach <= zero) then
          !
          ! Subsonic inflow
          !
          varh(1:neqns-1) = pv_in(1:neqns-1)
         !varh(neqns) = varh(neqns)
          !
        else if (mach <= -one) then
          !
          ! Supersonic inflow
          !
          varh(1:neqns) = pv_in(1:neqns)
          !
        else
          !
          write (iout,*) "Error with bc_generic_inflow boundary condition!"
          !
        end if
        !
        usp(:,kg) = vsp2u_sp( varh(:) )
        !
        faceusp(:,k,2,nf) = usp(:,kg)
        !
        facedusp(:,:,k,2,nf) = zero
        !
      end do
      !
    else if (bface(1,nf) == bc_generic_outflow) then
      !
      ! #####################################################
      ! #####################################################
      ! ##  Generic outflow boundary condition
      ! ##    - Use subsonic or supersonic outflow boundary
      ! ##      conditions based on the Mach number of the
      ! ##      normal outflow velocity
      ! #####################################################
      ! #####################################################
      !
      ! Store the transpose of usp for the host cell to reduce cache misses
      !
      do k = h1,h2
        tusp(k-h1+1,:) = usp(:,k)
      end do
      !
      ! Now loop through the boundary face flux points to get the boundary
      ! condition solution at each flux point
      !
      do k = 1,n_fp1de
        !
        ! Unit normal and tangent at this flux point
        !
        rnrm(:) = unit_vector( geofa(:,k,nf) )
        !
        hk = hkb-1 + k ! index for current host  cell flux point
        gk = gkb+1 - k ! index for current ghost cell flux point
        kg = g1-1 + gk ! index for current ghost cell solution point
        !
        ! Interpolate the conservative variables for the host cell
        ! solution points to the current host cell flux point.
        !
        do m = 1,neqns
          varh(m) = dot(interp_mat(1:hp,hk,hw),tusp(1:hp,m))
        end do
        !
        ! Now convert the current host cell flux point
        ! conservative variables to primitive variables
        !
        varh(:) = usp2v_sp( varh(:) )
        !
        ! Compute the speed of sound at this flux point
        !
        ai = sqrt(max(gam*varh(neqns)/varh(1),rndoff))
        !
        ! Compute the Mach number of the normal velocity at this flux point
        !
        mach = dot( varh(2:neqns-1) , rnrm(1:nr) ) / ai
        !
        if (mach < one .and. mach >= zero) then
          !
          ! Subsonic outflow
          !
         !varh(1:neqns-1) = varh(1:neqns-1)
          varh(neqns) = pv_in(neqns)
          !
        else if (mach >= one) then
          !
          ! Supersonic outflow
          !
         !varh(1:neqns) = varh(1:neqns)
          !
        else
          !
          write (iout,*) "Error with bc_generic_outflow boundary condition!"
          !
        end if
        !
        usp(:,kg) = vsp2u_sp( varh(:) )
        !
        faceusp(:,k,2,nf) = usp(:,kg)
        !
        facedusp(:,:,k,2,nf) = zero
        !
      end do
      !
    else if (bface(1,nf) == bc_sub_inflow) then
      !
      ! #####################################################
      ! #####################################################
      ! ##  Subsonic inflow boundary condition
      ! ##    - Fix density to free stream
      ! ##    - Fix velocities to free stream
      ! ##    - Extrapolate pressure
      ! #####################################################
      ! #####################################################
      !
      ! Store the transpose of usp for the host cell to reduce cache misses
      !
      do k = h1,h2
        tusp(k-h1+1,:) = usp(:,k)
      end do
      !
      ! Now loop through the boundary face flux points to get the boundary
      ! condition solution at each flux point
      !
      do k = 1,n_fp1de
        !
        hk = hkb-1 + k ! index for current host  cell flux point
        gk = gkb+1 - k ! index for current ghost cell flux point
        kg = g1-1 + gk ! index for current ghost cell solution point
        !
        ! Interpolate the conservative variables for the host cell
        ! solution points to the current host cell flux point.
        !
        do m = 1,neqns
          varh(m) = dot(interp_mat(1:hp,hk,hw),tusp(1:hp,m))
        end do
        !
        ! Now convert the current host cell flux point
        ! conservative variables to primitive variables
        !
        varh(:) = usp2v_sp( varh(:) )
        !
        varh(1:neqns-1) = pv_in(1:neqns-1)
       !varh(neqns) = varh(neqns)
        !
        usp(:,kg) = vsp2u_sp( varh(:) )
        !
        faceusp(:,k,2,nf) = usp(:,kg)
        !
        facedusp(:,:,k,2,nf) = zero
        !
      end do
      !
    else if (bface(1,nf) == bc_sub_outflow) then
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
      ! Store the transpose of usp for the host cell to reduce cache misses
      !
      do k = h1,h2
        tusp(k-h1+1,:) = usp(:,k)
      end do
      !
      ! Now loop through the boundary face flux points to get the boundary
      ! condition solution at each flux point
      !
      do k = 1,n_fp1de
        !
        hk = hkb-1 + k ! index for current host  cell flux point
        gk = gkb+1 - k ! index for current ghost cell flux point
        kg = g1-1 + gk ! index for current ghost cell solution point
        !
        ! Interpolate the conservative variables for the host cell
        ! solution points to the current host cell flux point.
        !
        do m = 1,neqns
          varh(m) = dot(interp_mat(1:hp,hk,hw),tusp(1:hp,m))
        end do
        !
        ! Now convert the current host cell flux point
        ! conservative variables to primitive variables
        !
        varh(:) = usp2v_sp( varh(:) )
        !
       !varh(1:neqns-1) = varh(1:neqns-1)
        varh(neqns) = pv_in(neqns)
        !
        usp(:,kg) = vsp2u_sp( varh(:) )
        !
        faceusp(:,k,2,nf) = usp(:,kg)
        !
        facedusp(:,:,k,2,nf) = zero
        !
      end do
      !
    else if (bface(1,nf) == bc_sup_inflow) then
      !
      ! #####################################################
      ! #####################################################
      ! ##  Supersonic inflow boundary condition
      ! ##    - Fix all variables to free stream
      ! #####################################################
      ! #####################################################
      !
      ! Now loop through the boundary face flux points to get the boundary
      ! condition solution at each flux point
      !
      do k = 1,n_fp1de
        !
        gk = gkb+1 - k ! index for current ghost cell flux point
        kg = g1-1 + gk ! index for current ghost cell solution point
        !
        usp(:,kg) = cv_in(:)
        !
        faceusp(:,k,2,nf) = usp(:,kg)
        !
        facedusp(:,:,k,2,nf) = zero
        !
      end do
      !
    else if (bface(1,nf) == bc_sup_outflow) then
      !
      ! #####################################################
      ! #####################################################
      ! ##  Supersonic outflow boundary condition
      ! ##    - Extrapolate all variables
      ! #####################################################
      ! #####################################################
      !
      ! Store the transpose of usp for the host cell to reduce cache misses
      !
      do k = h1,h2
        tusp(k-h1+1,:) = usp(:,k)
      end do
      !
      ! Now loop through the boundary face flux points to get the boundary
      ! condition solution at each flux point
      !
      do k = 1,n_fp1de
        !
        hk = hkb-1 + k ! index for current host  cell flux point
        gk = gkb+1 - k ! index for current ghost cell flux point
        kg = g1-1 + gk ! index for current ghost cell solution point
        !
        ! Interpolate the conservative variables for the host cell
        ! solution points to the current host cell flux point.
        !
        do m = 1,neqns
          varh(m) = dot(interp_mat(1:hp,hk,hw),tusp(1:hp,m))
        end do
        !
        ! Now convert the current host cell flux point
        ! conservative variables to primitive variables
        !
        usp(:,kg) = varh(:)
        !
        faceusp(:,k,2,nf) = usp(:,kg)
        !
        facedusp(:,:,k,2,nf) = zero
        !
      end do
      !
    else if ( any(bface(1,nf) == [bc_slip_wall,bc_symmetry]) ) then
      !
      ! #####################################################
      ! #####################################################
      ! ##  Slip wall / symmetry boundary condition
      ! #####################################################
      ! #####################################################
      !
      ! Store the transpose of usp for the host cell to reduce cache misses
      !
      do k = h1,h2
        tusp(k-h1+1,:) = usp(:,k)
      end do
      !
      ! Now loop through the boundary face flux points to get the boundary
      ! condition solution at each flux point
      !
      do k = 1,n_fp1de
        !
        ! Unit normal and tangent at this flux point
        !
        rnrm(:) =  unit_vector( geofa(:,k,nf) )
        !
        rtan(:) = [rnrm(2),-rnrm(1)]
        !
        hk = hkb-1 + k ! index for current host  cell flux point
        gk = gkb+1 - k ! index for current ghost cell flux point
        kg = g1-1 + gk ! index for current ghost cell solution point
        !
        ! Interpolate the conservative variables for the host cell
        ! solution points to the current host cell flux point.
        !
        do m = 1,neqns
          varh(m) = dot(interp_mat(1:hp,hk,hw),tusp(1:hp,m))
        end do
        !
        ! Now convert the current host cell flux point
        ! conservative variables to primitive variables
        !
        varh(:) = usp2v_sp( varh(:) )
        !
        u = varh(2) ! host cell x-velocity at flux point
        v = varh(3) ! host cell y-velocity at flux point
        !
        un = u*rnrm(1) + v*rnrm(2) ! host cell normal  velocity at flux point
        ut = u*rnrm(2) - v*rnrm(1) ! host cell tangent velocity at flux point
        !
       !varh(1) = varh(1)
        varh(2) = -un*rnrm(1) + ut*rnrm(2)
        varh(3) = -un*rtan(1) + ut*rtan(2)
       !varh(4) = varh(4)
        !
        usp(:,kg) = vsp2u_sp( varh(:) )
        !
        faceusp(:,k,2,nf) = usp(:,kg)
        !
        facedusp(:,:,k,2,nf) = zero
        !
      end do
      !
    else if (bface(1,nf) == bc_adiabatic_wall) then
      !
      ! #####################################################
      ! #####################################################
      ! ##  No-slip adiabatic wall
      ! #####################################################
      ! #####################################################
      !
      ! Store the transpose of usp for the host cell to reduce cache misses
      !
      do k = h1,h2
        tusp(k-h1+1,:) = usp(:,k)
      end do
      !
      ! Now loop through the boundary face flux points to get the boundary
      ! condition solution at each flux point
      !
      do k = 1,n_fp1de
        !
        ! Unit normal and tangent at this flux point
        !
        rnrm(:) =  unit_vector( geofa(:,k,nf) )
        !
        hk = hkb-1 + k ! index for current host  cell flux point
        gk = gkb+1 - k ! index for current ghost cell flux point
        kg = g1-1 + gk ! index for current ghost cell solution point
        !
        ! Interpolate the conservative variables for the host cell
        ! solution points to the current host cell flux point.
        !
        do m = 1,neqns
          varh(m) = dot(interp_mat(1:hp,hk,hw),tusp(1:hp,m))
        end do
        !
        ! Now convert the current host cell flux point
        ! conservative variables to primitive variables
        !
        varh(:) = usp2v_sp( varh(:) )
        !
       !varh(1) = varh(1)
        varh(2) = zero
        varh(3) = zero
       !varh(4) = varh(4)
        !
        usp(:,kg) = vsp2u_sp( varh(:) )
        !
        faceusp(:,k,2,nf) = usp(:,kg)
        !
        facedusp(:,:,k,2,nf) = zero
        !
      end do
      !
    else if (bface(1,nf) == bc_isothermal_wall) then
      !
      ! #####################################################
      ! #####################################################
      ! ##  No-slip isothermal wall
      ! #####################################################
      ! #####################################################
      !
      ! THIS NEEDS TO BE SET UP!!!
      !
    else if (bface(1,nf) == bc_periodic) then
      !
      ! #####################################################
      ! #####################################################
      ! ##  Periodic boundary condition
      ! #####################################################
      ! #####################################################
      !
      p1 = pts_of_periodic_ptr(nf)+1
      p2 = pts_of_periodic_ptr(nf+1)
      !
      usp(:,g1:g2) = usp(:,pts_of_periodic(p1:p2))
      do k = g1,g2
        tusp(k-g1+1,:) = usp(:,k)
      end do
      !
      do k = 1,n_fp1de
        !
        gk = gkb+1 - k ! index for current ghost cell flux point
        !
        ! Use the Lagrange polynomial to interpolate the conservative variables from
        ! the ghost cell solution points to the current ghost cell flux point.
        !
        do m = 1,neqns
          faceusp(m,k,2,nf) = dot(interp_mat(1:gp,gk,gw),tusp(1:gp,m))
        end do
        !
      end do
      !
   !else if (bface(1,nf) == bc_cpu_bnd) then
     !!
     !do k = 1,n_sp1de
     !  !
     !  gk = gkb+1 - k ! index for current ghost cell flux point
     !  !
     !  do m = 1,neqns
     !    faceusp(m,k,2,nf) = dot(interp_s2f(1:gp,gk,gw),usp(m,g1:g2))
     !  end do
     !  !
     !end do
     !!
    else if (bface(1,nf) /= bc_cpu_bnd) then
      !
      ! This shouldn't be here and is an unknown boundary condition
      !
      write (iout,1) bface(1,nf)
      call stop_recon(abort,"getbc")
      !
    end if
    !
  end do
  !
  1 format (/," ***************** ERROR! *****************",/, &
              "  UNKNOWN BOUNDARY CONDITION IN SUB GETBC!",/, &
              "  BOUNDARY CONDITION VALUE IS: ",i0,/, &
              " ***************** ERROR! *****************",/)
  !
end subroutine getbc
!
!###############################################################################
!
include "Functions/vortex.f90"
include "Functions/usp2v.f90"
include "Functions/vsp2u.f90"
include "Functions/viscosity.f90"
include "Functions/entropy.f90"
include "Functions/pressure.f90"
!
!###############################################################################
!
end module generic_mod
