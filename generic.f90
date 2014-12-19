module generic_mod
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
  public :: rgeterr
  public :: getmaxent
  public :: compute_final_error
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
  use order_mod,  only : n_sp2dt
  use geovar,  only : nr,nfbnd,nface,interp_s2f
  use geovar,  only : cells_on_face,pts_in_cell_ptr
  use ovar,    only : gam,Timestep_Type,Final_Time,time,CFL
  use metvar,  only : geofa,jacobian
  use flowvar, only : dtsp,faceusp
  use order_mod,  only : nsp => n_sp1de
  use geovar, only : cell_geom
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: dtmin
  !
  !.. Local Scalars ..
  integer  :: nf,k
  integer  :: kbl,l1,l2,npl,wlocl,kl
  integer  :: kbr,r1,r2,npr,wlocr,kr
  integer  :: host_cell,left_cell,right_cell
  integer  :: host_geom,left_geom,right_geom
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
    !
    wlocl = merge( 1 , 2 , npl == n_sp2dt ) ! check if host_cell is a tria or quad
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
      jacl = dot( interp_s2f(1:npl,kl,wlocl) , jacobian(l1:l2) )
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
    right_geom = cell_geom(right_cell)
    !
    wlocl = merge( 1 , 2 , npl == n_sp2dt ) ! check if left_cell is a tria or quad
    wlocr = merge( 1 , 2 , npr == n_sp2dt ) ! check if right_cell is a tria or quad
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
      jacl = dot( interp_s2f(1:npl,kl,wlocl) , jacobian(l1:l2) )
      jacr = dot( interp_s2f(1:npr,kr,wlocr) , jacobian(r1:r2) )
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
  use order_mod,  only : n_sp1de
  use geovar,  only : nr,n_solpts
  use ovar,    only : gam,Final_Time,time,CFL,Timestep_Type
  use metvar,  only : jacobian,metrics
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
  dtc = CFL * two / real( n_sp1de**2 , kind=wp )
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
      do m = nec+1,nq
        residual(m,n) = residual(nec,n)
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
  use order_mod,  only : n_sp2dt,n_sp2dq
  use order_mod,  only : p_order,n_order
  use geovar,  only : projct_s2p
  !
  ! PTYPE : gives the location of the variables needing projection
  !         = 1 -> at solution points
  !         = 2 -> at Gauss    points
  !         = 3 -> at Lobatto  points
  !         = 4 -> at solution points (used to be Radau)
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: ncell
  integer,  dimension(:),      intent(in) :: pts_in_cell_ptr
  real(wp), dimension(:,:), intent(inout) :: array
  integer,  optional,          intent(in) :: ptype
  !
  !.. Local Scalars ..
  integer :: k,n,l,n1,n2,nsp,wlocl,ndim,ploc
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
      wlocl = merge(1,2,nsp==n_sp2dt)
      !
      do k = 1,nsp
        do l = 1,ndim
          tmp(k,l) = array(l,n1-1+k)
        end do
      end do
      !
      do k = 1,nsp
        do l = 1,ndim
          array(l,n1-1+k) = dot( projct_s2p(1:nsp,k,wlocl,ploc) , tmp(1:nsp,l) )
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
  use geovar,  only : n_solpts,write_variable
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
  use geovar,     only : ncell,pts_in_cell_ptr
  use geovar,     only : cell_geom,cell_order
  use filter_mod, only : filter
  use flowvar,    only : usp
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,m
  !
  !.. Local Pointers ..
  real(wp), pointer, contiguous :: fmat(:,:)
  !
continue
  !
  ! Initialize fmat to disassociated
  !
  fmat => null()
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1
    n2 = pts_in_cell_ptr(n+1)
    !
    ! Get the filter matrix for the geometry and order of this cell
    !
    fmat => filter(cell_geom(n),cell_order(n))%mat
    !
    do m = 1,nq
      usp(m,n1:n2) = matmul( fmat , usp(m,n1:n2) )
    end do
    !
    ! Disassociate fmat before going to the next cell
    !
    if (associated(fmat)) fmat => null()
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
  ! ultwo(nmb:nme) = directional momentums
  ! ultwo(nq)     = energy
  ! ultwo(nq+1)   = entropy
  !
  ! vltwo(1)      = density
  ! vltwo(nmb:nme) = directional velocities
  ! vltwo(nq)     = pressure
  ! vltwo(nq+1)   = magnitude of velocity
  !
  !.. Use Statements ..
  use geovar,  only : ncell,pts_in_cell_ptr
  use geovar,  only : cell_geom,cell_order
  use ovar,    only : vlinf,vlone,vltwo,vnloc,vkloc,total_volume
  use ovar,    only : ulinf,ulone,ultwo,unloc,ukloc
  use metvar,  only : jacobian
  use flowvar, only : usp,uoldsp
  use parallel_mod, only : cell_map,cell_map_ptr
  use quadrature_mod, only : std_elem
  !
  !.. Local Scalars ..
  integer  :: n,n1,n2,k,m,nc1
  integer  :: this_geom,this_order
  real(wp) :: sn,se
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nq) :: uold
  real(wp), dimension(1:nq) :: unew
  real(wp), dimension(1:nq+1,1:2) :: ett
  real(wp), dimension(1:nq+1,1:2,1:5) :: error
  real(wp), dimension(1:nq+1,1:2,1:5,1:ncpu) :: all_error
  !
  !.. Local Pointers ..
  real(wp), pointer, contiguous :: weights(:)
  !
  ! error(1:nq+1,:,:) = error value for each flow equation + extra variable
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
  weights => null()
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
    !
    this_geom = cell_geom(n)
    this_order = cell_order(n)
    !
    weights => std_elem(this_geom,this_order)%wts
    !
    ! Now loop through the solution points of this cell to compute the
    ! error between the solutions of the current and previous time steps
    !
    do k = n1,n2
      !
      unew = usp(:,k)
      uold = uoldsp(:,k)
      !
      ! Error for conserved variables at the current solution point
      !
      sn = entropy_cv_sp( unew )
      se = entropy_cv_sp( uold )
      !
      ett(1:nq,1) = abs( unew - uold )
      ett(nq+1,1) = sn - se
      !
      ! Error for primitive variables at the current solution point
      !
      ett(1:nq,2) = abs( usp2v_sp(unew) - usp2v_sp(uold) )
      ett(nq+1,2) = rnorm2( ett(nmb:nme,2) )
      !
      ! Now add the error for all variables at the current solution
      ! point to the running integral sum for the L1 and L2 error norms
      !
      error(:,:,1) = error(:,:,1) + jacobian(k)*weights(k-n1+1)*ett(:,:)
      error(:,:,2) = error(:,:,2) + jacobian(k)*weights(k-n1+1)*ett(:,:)*ett(:,:)
      !
      ! If any of these errors are greater than the previous maximum
      ! values, save the error value and record this location
      !
      where (abs(ett(:,:)) > error(:,:,3))
        error(:,:,3) = abs(ett(:,:))
        error(:,:,4) = real(cell_map(nc1+n),kind=wp)
        error(:,:,5) = real(k-n1+1,kind=wp)
      end where
      !
    end do
    !
    if (associated(weights)) weights => null()
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
  do m = 1,nq+1
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
subroutine rgeterr()
  !
  !.. Use Statements ..
  use order_mod,  only : n_sp2dt
  use geovar,  only : ncell,pts_in_cell_ptr,xyz
  use geovar,  only : cell_geom,cell_order
  use geovar,  only : global_pinc_ptr,global_xyz
  use ovar,    only : rlinf,rlone,rltwo,nloc,kloc,total_volume,gam,cv_in
  use ovar,    only : nentmax,kentmax,xentmax,yentmax,entmax,entL2,itestcase
  use ovar,    only : itcur,d_t,time,rl2mx
  use ovar,    only : max_wall_time,wall_time_left
  use ovar,    only : beg_wall_time,wall_time
  use ovar,    only : walltime_expiring,convergence_reached
  use ovar,    only : convergence_order_max_res,convergence_order_abs_res
  use metvar,  only : jacobian
  use flowvar, only : usp,uoldsp
  use parallel_mod, only : cell_map,cell_map_ptr
  use quadrature_mod, only : std_elem
  !
  !.. Local Scalars ..
  integer  :: n,n1,n2,k,nc1
  integer  :: hrs_left,min_left,sec_left
  integer  :: this_geom,this_order
  real(wp) :: rho,pres,sn,se
  real(wp) :: next_time_est
  real(wp) :: wall_time_last
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nq) :: vv
  real(wp), dimension(1:2) :: ett
  real(wp), dimension(1:2,1:7) :: error
  real(wp), dimension(1:2,1:7,1:ncpu) :: all_error
  !
  !.. Local Pointers ..
  real(wp), pointer, contiguous :: weights(:)
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
  weights => null()
  !
  error(:,:) = zero  ! intialize everything to zero
  error(:,3) = -eps6 ! make this slightly negative in case all error is zero
  !
  nc1 = cell_map_ptr(mypnum+1)
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1
    n2 = pts_in_cell_ptr(n+1)
    !
    this_geom = cell_geom(n)
    this_order = cell_order(n)
    !
    weights => std_elem(this_geom,this_order)%wts
    !
    do k = n1,n2
      !
      ! Get the change in density for this solution point
      !
      if (itestcase == Laminar_BndLyr) then
        ett(1) = abs( usp(nmb,k) - uoldsp(nmb,k) )
      else
        ett(1) = abs( usp(nec,k) - uoldsp(nec,k) )
      end if
      !
      ! Compute the entropy at this solution point
      !
      rho = usp(nec,k)
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
      error(:,1) = error(:,1) + jacobian(k)*weights(k-n1+1)*ett(:)
      error(:,2) = error(:,2) + jacobian(k)*weights(k-n1+1)*ett(:)*ett(:)
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
    !
    if (associated(weights)) weights => null()
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
                              global_pinc_ptr(nloc)+kloc, &
                              global_xyz(1,global_pinc_ptr(nloc)+kloc), &
                              global_xyz(2,global_pinc_ptr(nloc)+kloc), &
                              hrs_left,min_left,sec_left
    !
#else
    !
    walltime_expiring = fals
    !
    write (iout,iter_out_fmt) itcur,d_t,time,rltwo,entL2,nloc, &
                              global_pinc_ptr(nloc)+kloc, &
                              global_xyz(1,global_pinc_ptr(nloc)+kloc), &
                              global_xyz(2,global_pinc_ptr(nloc)+kloc)
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
  use ovar,   only : ulinf,ulone,ultwo,unloc,ukloc
  use ovar,   only : vlinf,vlone,vltwo,vnloc,vkloc
  !
  !.. Local Scalars ..
  integer :: ierr
  integer :: final_iter_digits
  integer :: max_cell_digits
  integer :: max_solpt_digits
  character(len=100) :: char_num
  character(len=1000) :: iter_header
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_generic_mod"
  !
continue
  !
  ! Allocate the error arrays
  !
  if (.not. allocated(unloc)) then
    allocate ( unloc(1:nq+1) , source=0 , stat=ierr )
    call alloc_error(pname,"unloc",1,ierr)
  end if
  if (.not. allocated(ukloc)) then
    allocate ( ukloc(1:nq+1) , source=0 , stat=ierr )
    call alloc_error(pname,"ukloc",1,ierr)
  end if
  if (.not. allocated(ulinf)) then
    allocate ( ulinf(1:nq+1) , source=zero , stat=ierr )
    call alloc_error(pname,"ulinf",1,ierr)
  end if
  if (.not. allocated(ulone)) then
    allocate ( ulone(1:nq+1) , source=zero , stat=ierr )
    call alloc_error(pname,"ulone",1,ierr)
  end if
  if (.not. allocated(ultwo)) then
    allocate ( ultwo(1:nq+1) , source=zero , stat=ierr )
    call alloc_error(pname,"ultwo",1,ierr)
  end if
  if (.not. allocated(vnloc)) then
    allocate ( vnloc(1:nq+1) , source=0 , stat=ierr )
    call alloc_error(pname,"vnloc",1,ierr)
  end if
  if (.not. allocated(vkloc)) then
    allocate ( vkloc(1:nq+1) , source=0 , stat=ierr )
    call alloc_error(pname,"vkloc",1,ierr)
  end if
  if (.not. allocated(vlinf)) then
    allocate ( vlinf(1:nq+1) , source=zero , stat=ierr )
    call alloc_error(pname,"vlinf",1,ierr)
  end if
  if (.not. allocated(vlone)) then
    allocate ( vlone(1:nq+1) , source=zero , stat=ierr )
    call alloc_error(pname,"vlone",1,ierr)
  end if
  if (.not. allocated(vltwo)) then
    allocate ( vltwo(1:nq+1) , source=zero , stat=ierr )
    call alloc_error(pname,"vltwo",1,ierr)
  end if
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
  use order_mod,  only : n_sp2dt
  use geovar,  only : ncell,pts_in_cell_ptr,xyz
  use geovar,  only : cell_geom,cell_order
  use ovar,    only : nentmax,kentmax,entmax,gam,cv_in
  use ovar,    only : xentmax,yentmax,entL2,total_volume,itestcase
  use metvar,  only : jacobian
  use flowvar, only : usp
  use parallel_mod, only : cell_map,cell_map_ptr
  use quadrature_mod, only : std_elem
  !
  !.. Local Scalars ..
  integer  :: n,n1,n2,k,nc1
  integer  :: this_geom,this_order
  real(wp) :: ett,sn,se
  !
  !.. Local Arrays ..
  real(wp), dimension(1:6) :: error
  real(wp), dimension(1:6,1:ncpu) :: all_error
  !
  !.. Local Pointers ..
  real(wp), pointer, contiguous :: weights(:)
  !
continue
  !
  !
  !
  weights => null()
  !
  error(:) = zero  ! intialize everything to zero
  error(2) = -eps6 ! make this slightly negative in case all error is zero
  !
  nc1 = cell_map_ptr(mypnum+1)
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1
    n2 = pts_in_cell_ptr(n+1)
    !
    this_geom = cell_geom(n)
    this_order = cell_order(n)
    !
    weights => std_elem(this_geom,this_order)%wts
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
      error(1) = error(1) + jacobian(k)*weights(k-n1+1)*ett*ett
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
    if (associated(weights)) weights => null()
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
