module time_mod
  !
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq
  !
  implicit none
  !
  private
  !
  public :: get_dtsp
  public :: rresi
  public :: rupdate_Classic_RK
  public :: rupdate_CK4_RK
  public :: rupdate_TVD2_RK
  public :: rupdate_TVD3_RK
  public :: rupdate_TVD4_RK
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
  use order_mod, only : nsp => n_sp1de
  !
  use geovar, only : nr,nfbnd,nface
  use geovar, only : cells_on_face,pts_in_cell_ptr
  use geovar, only : cell_geom,cell_order
  !
  use ovar, only : gam,Timestep_Type,Final_Time,time,CFL
  !
  use metrics_mod, only : geofa,jacobian
  !
  use flowvar, only : dtsp,faceusp
  !
  use interpolation_mod, only : interp
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: dtmin
  !
  !.. Local Scalars ..
  integer  :: nf,k
  integer  :: kbl,l1,l2,npl,kl
  integer  :: kbr,r1,r2,npr,kr
  integer  :: host_cell,host_geom,host_order
  integer  :: left_cell,left_geom,left_order
  integer  :: right_cell,right_geom,right_order
  real(wp) :: ds,eigni,eignv,dttemp,pp1sq
  real(wp) :: dsl,vell,al,jacl,dtl,mul,hl
  real(wp) :: dsr,velr,ar,jacr,dtr,mur,hr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nq) :: vl,vr
  !
  !.. Local Pointers ..
  real(wp), pointer, contiguous :: host_interp(:,:)
  real(wp), pointer, contiguous :: left_interp(:,:)
  real(wp), pointer, contiguous :: right_interp(:,:)
  !
continue
  !
  ! Initialize the local pointers to disassociated
  !
  host_interp => null()
  left_interp => null()
  right_interp => null()
  !
  pp1sq = real( nsp**2 , kind=wp )
  !
  dtmin = huge(zero)
  dtsp(:) = dtmin
  !
  do nf = 1,nfbnd
    !
    host_cell = cells_on_face(1,nf) ! left cell on face
    !
    kbl = cells_on_face(5,nf) ! beginning flux point for left cell
    !
    l1 = pts_in_cell_ptr(host_cell)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(host_cell+1) ! ending    index for sol-pts in left cell
    !
    npl = l2-l1+1 ! number of solution points in left cell
    !
    host_geom = cell_geom(host_cell)
    host_order = cell_order(host_cell)
    !
    ! Assign local pointer as an alias to the
    ! interpolation matrix to simplify the code
    !
    host_interp => interp(host_geom,host_order)%toFace(host_order)%mat
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
      dsl = ds * cl(host_geom) ! face length/area for left cell
      !
      ! Normal velocity at this flux point
      !
     !vell = dot( vl(nmb:nme) , unit_vector(geofa(1:nr,k,nf)) )
      vell = selfdot( vl(nmb:nme) )
      !
      ! Speed of sound at this face
      !
      al = max( gam*vl(nee)/vl(nec) , rndoff )
      !
      ! Interpolate the jacobian to the flux point
      !
      jacl = dot( host_interp(1:npl,kl) , jacobian(l1:l2) )
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
    ! Disassociate host_interp before continuing to the next boundary face
    !
    if (associated(host_interp)) host_interp => null()
    !
  end do
  !
  do nf = nfbnd+1,nface
    !
    left_cell = cells_on_face(1,nf) ! left  cell on face
    right_cell = cells_on_face(2,nf) ! right cell on face
    !
    kbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    kbr = cells_on_face(6,nf) ! beginning flux point for right cell
    !
    l1 = pts_in_cell_ptr(left_cell)+1 ! beginning index for sol-pts in left  cell
    l2 = pts_in_cell_ptr(left_cell+1) ! ending    index for sol-pts in right cell
    !
    r1 = pts_in_cell_ptr(right_cell)+1 ! beginning index for sol-pts in left  cell
    r2 = pts_in_cell_ptr(right_cell+1) ! ending    index for sol-pts in right cell
    !
    npl = l2-l1+1 ! number of solution points in left  cell
    npr = r2-r1+1 ! number of solution points in right cell
    !
    left_geom = cell_geom(left_cell)
    left_order = cell_order(left_cell)
    !
    right_geom = cell_geom(right_cell)
    right_order = cell_order(right_cell)
    !
    ! Assign local pointer as an alias to the
    ! interpolation matrix to simplify the code
    !
    left_interp => interp(left_geom,left_order)%toFace(left_order)%mat
    right_interp => interp(right_geom,right_order)%toFace(right_order)%mat
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
      dsl = ds * cl(left_geom) ! face length/area for left  cell
      dsr = ds * cl(right_geom) ! face length/area for right cell
      !
      ! Normal velocity at this flux point
      !
     !vell = dot( vl(nmb:nme) , unit_vector(geofa(1:nr,k,nf)) )
     !velr = dot( vr(nmb:nme) , unit_vector(geofa(1:nr,k,nf)) )
      vell = selfdot( vl(nmb:nme) )
      velr = selfdot( vr(nmb:nme) )
      !
      ! Speed of sound at this face
      !
      al = max( gam*vl(nee)/vl(nec) , rndoff )
      ar = max( gam*vr(nee)/vr(nec) , rndoff )
      !
      ! Interpolate the jacobian to the flux point
      !
      jacl = dot( left_interp(1:npl,kl) , jacobian(l1:l2) )
      jacr = dot( right_interp(1:npr,kr) , jacobian(r1:r2) )
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
    ! Disassociate left_interp and right_interp before
    ! continuing to the next interior face
    !
    if (associated(left_interp)) left_interp => null()
    if (associated(right_interp)) right_interp => null()
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
  use order_mod,  only : n_order
  use geovar,  only : nr,n_solpts
  use ovar,    only : gam,Final_Time,time,CFL,Timestep_Type
  use metrics_mod,  only : jacobian,metrics
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
  real(wp), dimension(1:nr) :: vel
  !
continue
  !
  gm1 = gam - one
  !
  dtmin = huge(zero)
  !
  dtc = CFL * two / real( (n_order+1)**2 , kind=wp )
  !
  do n = 1,n_solpts
    !
    dxdr = metrics(1,1,n)
    dydr = metrics(2,1,n)
    dxds = metrics(1,2,n)
    dyds = metrics(2,2,n)
    !
    rho  = usp(nec,n)
    vel  = usp(nmb:nme,n)/rho
    pres = gm1 * (usp(nee,n) - half*rho*selfdot(vel))
    !
    aspd = sqrt( max( gam*pres/rho , rndoff ) )
    !
    alfa = abs( vel(1)*dydr - vel(2)*dxdr )
    beta = abs( vel(1)*dyds - vel(2)*dxds )
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
  use geovar,  only : n_solpts
  use geovar,  only : ncell,pts_in_cell_ptr
  use geovar,  only : cell_geom,cell_order
  use ovar,    only : itestcase
  use metrics_mod,  only : jacobian
  use flowvar, only : residual
  use projection_mod, only : project_to_porder
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
      do m = nec+1,nq
        residual(m,n) = residual(nec,n)
      end do
    end do
  end if
  !
  call project_to_porder(ncell,cell_geom,cell_order,pts_in_cell_ptr,residual)
  !
end subroutine rresi
!
!###############################################################################
!
subroutine rupdate_Classic_RK
  !
  ! First, store d(utd)/dt for runge kutta time
  ! stepping then update usp
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use ovar,    only : nst,num_stages,Filter_Option
  use flowvar, only : rssprk,residual,usp,uoldsp,dtsp
  use filter_mod, only : filter_usp
  !
  !.. Local Scalars ..
  integer   :: n,m,nt
  real(wp) :: temp
  !
continue
  !
  ! Store the residual for this Runge-Kutta step in rssprk
  !
  rssprk(:,:,nst) = residual(:,:)
  !
  ! Update usp
  !
  do n = 1,n_solpts
    do m = 1,nq
      temp = zero
      do nt = 1,nst
        temp = temp + alphark(nt,nst,num_stages)*rssprk(m,n,nt)
      end do
      usp(m,n) = uoldsp(m,n) + dtsp(n)*temp
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
  use filter_mod, only : filter_usp
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
  use filter_mod, only : filter_usp
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
  use filter_mod, only : filter_usp
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
  use filter_mod, only : filter_usp
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
include "Functions/usp2v.f90"
include "Functions/viscosity.f90"
!
!###############################################################################
!
end module time_mod
