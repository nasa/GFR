module time_mod
  !
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq,ntk,ntl
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
  public :: update_time_ave
  public :: initialize_time_ave
  !
  real(wp), parameter :: o = zero
  !
  ! Classic N stage Runge-Kutta coefficients
  !
  real(wp), parameter, dimension(4,4,4) :: &
  alphaRK = reshape( &
                        ! 1-Stage Coefficients
       (/ONE ,o,o,o,o, o  ,o,o, o  ,o,  o   ,o, o  , o  , o  , o  ,   &
                        ! 2-Stage Coefficients
         HALF,o,o,o,o,ONE ,o,o, o  ,o,  o   ,o, o  , o  , o  , o  ,   &
                        ! 3-Stage Coefficients
         ONE3,o,o,o,o,TWO3,o,o,ONE4,o,THREE4,o, o  , o  , o  , o  ,   &
                        ! 4-Stage Coefficients
         HALF,o,o,o,o,HALF,o,o, o  ,o, ONE  ,o,ONE6,ONE3,ONE3,ONE6/), &
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
    alphaRK_TVD = reshape( &
         (/ 1.00000000000000_wp,0.00000000000000_wp,0.0_wp, &
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
  real(wp), save, allocatable :: ave(:)
  !
  type :: time_ave_var_t
    integer, allocatable :: idx(:)
  end type time_ave_var_t
  !
  type(time_ave_var_t), save, allocatable :: time_ave_var(:)
  !
contains
!
!###############################################################################
!
subroutine get_dtsp(dtmin)
  !
  !.. Use Statements ..
  use ovar,    only : Timestep_Type,time,d_t
  use ovar,    only : itestcase,time_ref
  use ovar,    only : write_TaylorGreen_8s_solution
  use ovar,    only : write_TaylorGreen_full_solution
  use flowvar, only : dtsp
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: dtmin
  !
  !.. Local Scalars ..
  real(wp) :: tstar,new_tstar
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_dtsp"
  !
continue
  !
  if (Timestep_Type == Constant_Timestep) then
    !
    call get_current_cfl
    !
    call get_constant_timestep(dtmin)
    !
  else if (abs(Timestep_Type) == Cell_Timestep) then
    !
    call get_current_cfl
    !
    call get_cell_timestep(dtmin)
    !
  else if (abs(Timestep_Type) == Point_Timestep) then
    !
    call get_current_cfl
    !
    call get_point_timestep(dtmin)
    !
  end if
  !
 !write (iout,1) time,dtmin,time+dtmin
  !
  if (itestcase == Taylor_Green_Vortex) then
    !
    tstar = time*time_ref
    new_tstar = tstar + dtmin*time_ref
    !
    if (tstar < eight .and. new_tstar >= eight) then
      !
      dtmin = (eight - tstar) / time_ref
      dtsp(:) = dtmin
      write_TaylorGreen_8s_solution = true
      !
    else if (tstar < three .and. new_tstar >= three) then
      !
      dtmin = (three - tstar) / time_ref
      dtsp(:) = dtmin
      write_TaylorGreen_full_solution = true
      !
    else if (tstar < five .and. new_tstar >= five) then
      !
      dtmin = (five - tstar) / time_ref
      dtsp(:) = dtmin
      write_TaylorGreen_full_solution = true
      !
    else if (tstar < seven .and. new_tstar >= seven) then
      !
      dtmin = (seven - tstar) / time_ref
      dtsp(:) = dtmin
      write_TaylorGreen_full_solution = true
      !
    else if (tstar < nine .and. new_tstar >= nine) then
      !
      dtmin = (nine - tstar) / time_ref
      dtsp(:) = dtmin
      write_TaylorGreen_full_solution = true
      !
    else if (tstar < 11.0_wp .and. new_tstar >= 11.0_wp) then
      !
      dtmin = (11.0_wp - tstar) / time_ref
      dtsp(:) = dtmin
      write_TaylorGreen_full_solution = true
      !
    else if (tstar < 15.0_wp .and. new_tstar >= 15.0_wp) then
      !
      dtmin = (15.0_wp - tstar) / time_ref
      dtsp(:) = dtmin
      write_TaylorGreen_full_solution = true
      !
    end if
    !
  end if
  !
  time = time + dtmin
  d_t = dtmin
  !
  ! Format Statements
  !
  1 format (" Old Time = ",es14.6,"; Timestep = ",es14.6,"; New Time = ",es14.6)
  !
end subroutine get_dtsp
!
!###############################################################################
!
subroutine get_current_cfl
  !
  !.. Use Statements ..
  use ovar, only : cfl,cfl_beg,cfl_end,cfl_cycles,itcur
  use ovar, only : cfl_cycle_start,cfl_cycle_end
  !
  !.. Local Scalars ..
  real(wp) :: cycle_ratio
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_current_cfl"
  !
continue
  !
  if (itcur < cfl_cycle_end) then
    if (cfl_cycles < 2) then
      cycle_ratio = one
    else
      cycle_ratio = real(itcur-cfl_cycle_start-1,kind=wp) / &
                    real(cfl_cycles-1,kind=wp)
    end if
    !
    cfl = cfl_beg + (cfl_end - cfl_beg) * cycle_ratio
  else
    cfl = cfl_end
  end if
  !
end subroutine get_current_cfl
!
!###############################################################################
!
subroutine get_constant_timestep(dtmin)
  !
  !.. Use Statements ..
  use ovar,    only : constant_dt,time_ref,Final_Time,time,CFL
  use ovar,    only : this_is_final_timestep
  use flowvar, only : dtsp
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: dtmin
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_constant_timestep"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  dtmin = CFL * constant_dt
  !
  if (Final_Time > zero) then
    if (time + dtmin > Final_Time/time_ref + eps9) then
      dtmin = Final_Time/time_ref - time
      this_is_final_timestep = true
    end if
  end if
  !
  dtsp(:) = dtmin
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_constant_timestep
!
!###############################################################################
!
subroutine get_cell_timestep(dtmin)
  !
  !.. Use Statements ..
  use order_mod, only : n_order,maxSP,geom_solpts
  !
  use geovar, only : nr,nfbnd,nface
  use geovar, only : face,cell
  !
  use ovar, only : gam,Timestep_Type,time_ref,Final_Time,time,CFL
  use ovar, only : governing_equations
  use ovar, only : this_is_final_timestep
  !
  use metrics_mod, only : geofa,metrics_dt
  !
  use flowvar, only : dtsp,usp
  !
  use interpolation_mod, only : interp
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: dtmin
  !
  !.. Local Scalars ..
  integer  :: nf,k,m
  integer  :: l1,l2,npl,kl
  integer  :: r1,r2,npr,kr
  integer  :: face_geom,face_order
  integer  :: host_cell,host_geom,host_order
  integer  :: left_cell,left_geom,left_order
  integer  :: right_cell,right_geom,right_order
  real(wp) :: ds,eigni,eignv,pp1sq
  real(wp) :: dsl,vell,al,jacl,dtl,mul,hl
  real(wp) :: dsr,velr,ar,jacr,dtr,mur,hr
  real(wp) :: psql,psqr
  !
  !.. Local Arrays ..
  real(wp) :: vl(1:nq),tlusp(1:maxSP,1:nq)
  real(wp) :: vr(1:nq),trusp(1:maxSP,1:nq)
  !
  integer, parameter :: make_p_plus_one = 0
 !integer, parameter :: make_p_plus_one = 1
  !
  logical(lk), parameter :: use_sqrt = true
 !logical(lk), parameter :: use_sqrt = fals
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_cell_timestep"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  mul = zero
  mur = zero
  !
  dtmin = huge(zero)
  dtsp(:) = dtmin
  !
  do nf = 1,nfbnd
    !
    host_cell = face(nf)%left%cell ! left cell on face
    !
    host_geom  = cell(host_cell)%geom  ! host cell geometry
    host_order = cell(host_cell)%order ! host cell order
    !
    l1 = cell(host_cell)%beg_sp ! beginning index for sol-pts in left cell
    l2 = cell(host_cell)%end_sp ! ending    index for sol-pts in left cell
    npl = l2-l1+1               ! number of solution points in left cell
    !
    face_geom = face(nf)%geom   ! geometry of the current face
    face_order = face(nf)%order ! order    of the current face
    !
    psql = real( (host_order+make_p_plus_one)**2 , kind=wp )
    !
    tlusp(1:npl,1:nq) = transpose( usp(1:nq,l1:l2) )
    !
    !
    do k = 1,geom_solpts(face_geom,face_order)
      !
      kl = face(nf)%left%get_fp_idx(face_geom,face_order,k)
      !
      ! Interpolate the conservative variables to this flux point
      vl(1:nq) = interp(host_geom,host_order)% &
                           toFace(face_order)% &
                 dot_mat( kl , tlusp(1:npl,1:nq) )
     !do m = 1,nq
     !  vl(m) = interp(host_geom,host_order)% &
     !                    toFace(face_order)% &
     !          dot( kl , tlusp(1:npl,m) )
     !end do
      ! Convert to primitive variables
      vl(1:nq) = usp2v_sp( vl(1:nq) )
      !
      if (governing_equations == NavierStokes_Eqns) then
        mul = viscosity_pv_sp( vl(1:nq) )
      end if
      !
      ds = norm2( geofa(nf)%v(1:nr,k) )
      dsl = ds * cl(host_geom) ! face length/area for left cell
      !
      ! Normal velocity at this flux point
      !
     !vell = dot( vl(nmb:nme) , unit_vector(geofa(nf)%v(1:nr,k)) )
      vell = selfdot( vl(nmb:nme) )
      if (use_sqrt) vell = sqrt(vell)
      !
      ! Speed of sound at this face
      !
      al = max( gam*vl(nee)/vl(nec) , rndoff )
      if (use_sqrt) al = sqrt(al)
      !
      ! Interpolate the jacobian to the flux point
      !
      jacl = interp(host_geom,host_order)% &
                       toFace(face_order)% &
             dot( kl , metrics_dt(host_cell)%jac(:) )
      !
      hl = half * dsl / jacl
      !
      eigni = abs(vell) + al
      eignv = psql * hl * mul
      !
      ! Compute the timestep for this flux point
      !
      dtl = min( dtsp(l1) , CFL / (psql * hl * (eigni + eignv)) )
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
    left_cell  = face(nf)%left%cell  ! left  cell on face
    right_cell = face(nf)%right%cell ! right cell on face
    !
    left_geom  = cell(left_cell)%geom  ! left cell geometry
    left_order = cell(left_cell)%order ! left cell order
    !
    l1 = cell(left_cell)%beg_sp ! beginning index for sol-pts in left  cell
    l2 = cell(left_cell)%end_sp ! ending    index for sol-pts in left cell
    npl = l2-l1+1               ! number of solution points in left  cell
    !
    right_geom  = cell(right_cell)%geom  ! right cell geometry
    right_order = cell(right_cell)%order ! right cell order
    !
    r1 = cell(right_cell)%beg_sp ! beginning index for sol-pts in right cell
    r2 = cell(right_cell)%end_sp ! ending    index for sol-pts in right cell
    npr = r2-r1+1                ! number of solution points in right cell
    !
    face_geom = face(nf)%geom   ! geometry of the current face
    face_order = face(nf)%order ! order    of the current face
    !
    psql = real( (left_order +make_p_plus_one)**2 , kind=wp )
    psqr = real( (right_order+make_p_plus_one)**2 , kind=wp )
    !
    tlusp(1:npl,1:nq) = transpose( usp(1:nq,l1:l2) )
    trusp(1:npr,1:nq) = transpose( usp(1:nq,r1:r2) )
    !
    !
    do k = 1,geom_solpts(face_geom,face_order)
      !
      kl = face(nf)%left%get_fp_idx(face_geom,face_order,k)
      kr = face(nf)%right%get_fp_idx(face_geom,face_order,k)
      !
      ! Interpolate the conservative variables to this flux point
      vl(1:nq) = interp(left_geom,left_order)% &
                           toFace(face_order)% &
                 dot_mat( kl , tlusp(1:npl,1:nq) )
      vr(1:nq) = interp(right_geom,right_order)% &
                             toFace(face_order)% &
                 dot_mat( kr , trusp(1:npr,1:nq) )
     !do m = 1,nq
     !  vl(m) = interp(left_geom,left_order)% &
     !                    toFace(face_order)% &
     !          dot( kl , tlusp(1:npl,m) )
     !  vr(m) = interp(right_geom,right_order)% &
     !                      toFace(face_order)% &
     !          dot( kr , trusp(1:npr,m) )
     !end do
      ! Convert to primitive variables
      vl(1:nq) = usp2v_sp( vl(1:nq) )
      vr(1:nq) = usp2v_sp( vr(1:nq) )
      !
      if (governing_equations == NavierStokes_Eqns) then
        mul = viscosity_pv_sp( vl(1:nq) )
        mur = viscosity_pv_sp( vr(1:nq) )
      end if
      !
      ds = norm2( geofa(nf)%v(1:nr,k) )
      dsl = ds * cl(left_geom) ! face length/area for left  cell
      dsr = ds * cl(right_geom) ! face length/area for right cell
      !
      ! Normal velocity at this flux point
      !
     !vell = dot( vl(nmb:nme) , unit_vector(geofa(nf)%v(1:nr,k)) )
     !velr = dot( vr(nmb:nme) , unit_vector(geofa(nf)%v(1:nr,k)) )
      vell = selfdot( vl(nmb:nme) )
      if (use_sqrt) vell = sqrt(vell)
      velr = selfdot( vr(nmb:nme) )
      if (use_sqrt) velr = sqrt(velr)
      !
      ! Speed of sound at this face
      !
      al = max( gam*vl(nee)/vl(nec) , rndoff )
      if (use_sqrt) al = sqrt(al)
      ar = max( gam*vr(nee)/vr(nec) , rndoff )
      if (use_sqrt) ar = sqrt(ar)
      !
      ! Interpolate the jacobian to the flux point for the left cell
      !
      jacl = interp(left_geom,left_order)% &
                       toFace(face_order)% &
             dot( kl , metrics_dt(left_cell)%jac(:) )
      !
      ! Interpolate the jacobian to the flux point for the right cell
      !
      jacr = interp(right_geom,right_order)% &
                         toFace(face_order)% &
             dot( kr , metrics_dt(right_cell)%jac(:) )
      !
      hl = half * dsl / jacl
      hr = half * dsr / jacr
      !
      eigni = max( abs(vell) + al , abs(velr) + ar )
      eignv = max( psql * hl * mul , psqr * hr * mur )
      !
      ! Compute the timestep for this flux point
      !
      dtl = min( dtsp(l1) , CFL / (psql * hl * (eigni + eignv) ) )
      dtr = min( dtsp(r1) , CFL / (psqr * hr * (eigni + eignv) ) )
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
    call mpi_allreduce(MPI_IN_PLACE,dtmin,1_int_mpi,mpi_flttyp, &
                       MPI_MIN,MPI_COMM_WORLD,mpierr)
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
  end if
  !
  ! Make sure the timestep size doesnt make the final time exceed Final_Time
  !
  if (Final_Time > zero) then
    if (time + dtmin > Final_Time/time_ref + eps9) then
      dtsp(:) = Final_Time/time_ref - time
      this_is_final_timestep = true
    end if
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_cell_timestep
!
!###############################################################################
!
subroutine get_point_timestep(dtmin)
  !
  !.. Use Statements ..
  use geovar,  only : nr,ncell,cell
  use ovar,    only : gam,time_ref,Final_Time,time,CFL,Timestep_Type
  use ovar,    only : this_is_final_timestep
  use ovar,    only : governing_equations
  use metrics_mod,  only : metrics_dt
  use flowvar, only : usp,dtsp
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: dtmin
  !
  !.. Local Scalars ..
  integer  :: n,nc,k,n1,n2,np
  integer  :: this_geom,this_order
  real(wp) :: dtpt,dtc,eigni,eignv
  real(wp) :: vmag,aspd,psq,visc
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: vel
  real(wp), dimension(1:nq) :: pv
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_point_timestep"
  !
 !integer, parameter :: make_p_plus_one = 0
  integer, parameter :: make_p_plus_one = 1
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  dtmin = huge(zero)
  !
  visc = zero
  !
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
    n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
    np = n2-n1+1         ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    psq = real( (this_order+make_p_plus_one)**2 , kind=wp )
    !
    dtc = CFL * two / psq
    !
    do k = 1,np
      !
      n = k + n1-1
      !
      pv(1:nq) = usp2v_sp( usp(1:nq,n) )
      !
      if (governing_equations == NavierStokes_Eqns) then
        visc = viscosity_cv_sp( usp(1:nq,n) )
      end if
      !
      vel(1:nr) = matmul( metrics_dt(nc)%met(k,1:nr,1:nr) , pv(nmb:nme) )
      !
      vmag = norm2(vel)
      !
      aspd = sqrt( max( gam*pv(nee)/pv(nec) , rndoff ) )
      !
      eigni = vmag + aspd
      !
      eignv = psq * visc
      !
      dtpt = dtc / (eigni + eignv)
      !
      dtmin = min( dtpt , dtmin )
      !
      dtsp(n) = dtpt
      !
    end do
    !
  end do
  !
  ! Get the minimum timestep across all processors
  !
  if (ncpu > 1) then
    call mpi_allreduce(MPI_IN_PLACE,dtmin,1_int_mpi,mpi_flttyp, &
                       MPI_MIN,MPI_COMM_WORLD,mpierr)
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
      if (time + dtmin > Final_Time/time_ref + eps9) then
        dtsp(:) = Final_Time/time_ref - time
        this_is_final_timestep = true
      end if
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_point_timestep
!
!###############################################################################
!
subroutine old_get_point_timestep(dtmin)
  !
  !.. Use Statements ..
  use order_mod,  only : n_order
  use geovar,  only : nr,ncell,cell
  use ovar,    only : gam,time_ref,Final_Time,time,CFL,Timestep_Type
  use ovar,    only : this_is_final_timestep
  use metrics_mod,  only : metrics_dt
  use flowvar, only : usp,dtsp
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: dtmin
  !
  !.. Local Scalars ..
  integer  :: k,n,nc,n1,n2
  real(wp) :: dtpt,gm1,dtc,eign
  real(wp) :: dxdr,dydr,dxds,dyds
  real(wp) :: rho,vx,vy,pres,aspd
  real(wp) :: alfa,beta,gama
  !
  real(wp) :: vel(1:nr)
  real(wp) :: met(1:nr,1:nr)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "old_get_point_timestep"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  gm1 = gam - one
  !
  dtmin = huge(zero)
  !
  dtc = CFL * two / real( (n_order+1)**2 , kind=wp )
  !
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    !
    do n = n1,n2
      !
      k = n-n1+1
      !
      met(1:nr,1:nr) = invert_matrix(metrics_dt(nc)%met(k,1:nr,1:nr))
      dxdr = met(1,1)
      dydr = met(2,1)
      dxds = met(1,2)
      dyds = met(2,1)
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
      dtpt = dtc * metrics_dt(nc)%jac(k) / eign
      !
      dtmin = min( dtpt , dtmin )
      !
      dtsp(n) = dtpt
      !
    end do
    !
  end do
  !
  ! Get the minimum timestep across all processors
  !
  if (ncpu > 1) then
    call mpi_allreduce(MPI_IN_PLACE,dtmin,1_int_mpi,mpi_flttyp, &
                       MPI_MIN,MPI_COMM_WORLD,mpierr)
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
      if (time + dtmin > Final_Time/time_ref + eps9) then
        dtsp(:) = Final_Time/time_ref - time
        this_is_final_timestep = true
      end if
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine old_get_point_timestep
!
!###############################################################################
!
subroutine rresi
  !
  !... calculates d(utd)/dt
  !
  !.. Use Statements ..
  use order_mod, only : p_order,n_order
  use geovar,  only : n_solpts,ncell,cell
  use ovar,    only : itestcase
  use metrics_mod,  only : metrics_dt
  use flowvar, only : residual
  use projection_mod, only : project_to_porder
  !
  !.. Local Scalars ..
  integer  :: m,n,nc,n1,n2
  real(wp) :: jac_inv
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "rresi"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  do nc = 1,ncell
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    do n = n1,n2
      jac_inv = one / metrics_dt(nc)%jac(n-n1+1)
      do m = 1,nq
        residual(m,n) = -residual(m,n) * jac_inv
      end do
    end do
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
  if (p_order < n_order) then
    call project_to_porder(cell,residual)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine rresi
!
!###############################################################################
!
subroutine stabilize_usp(nst,Apply_Filter,Apply_Limiter)
  !
  !.. Use Statements ..
  use ovar,            only : Limiter_Option,Filter_Option
  use filter_mod,      only : filter_usp
  use module_limiters, only : resolution_indicator
  use module_limiters, only : limit_by_projecting
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "stabilize_usp"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Apply the filter
  !
  if (Filter_Option > 0) then
    if (Apply_Filter) then
#ifdef PROFILE_ON
      if (mypnum == glb_root) write (iout,1) nst
      1 format (8x,"Applying Standard Filter : RK Stage = ",i0)
#endif
      call filter_usp
    end if
  end if
  !
  ! Apply the limiter
  !
  if (Limiter_Option /= 0) then
    if (Apply_Limiter) then
#ifdef PROFILE_ON
      if (mypnum == glb_root) write (iout,2) nst
      2 format (8x,"Applying Limiter : RK Stage = ",i0)
#endif
      if (Limiter_Option == -1) then
        call resolution_indicator(apply_projection=true)
      else if (Limiter_Option >= 1) then
        call limit_by_projecting
      end if
    end if
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine stabilize_usp
!
!###############################################################################
!
subroutine rupdate_Classic_RK(nst,Apply_Filter,Apply_Limiter)
  !
  ! First, store d(utd)/dt for runge kutta time
  ! stepping then update usp
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use ovar,    only : num_rk_stages
  use ovar,    only : Limiter_Option
  use ovar,    only : Filter_Option
  use flowvar, only : rssprk,residual,usp,uoldsp,dtsp
  use filter_mod, only : filter_usp
  use module_limiters, only : resolution_indicator
  use module_limiters, only : limit_by_projecting
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
  !
  !.. Local Scalars ..
  integer   :: n,m,nt
  real(wp) :: temp
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "rupdate_Classic_RK"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
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
        temp = temp + alphark(nt,nst,num_rk_stages)*rssprk(m,n,nt)
      end do
      usp(m,n) = uoldsp(m,n) + dtsp(n)*temp
    end do
  end do
  !
  ! Apply filtering or limiting to stabilize the solution
  !
  call stabilize_usp(nst,Apply_Filter,Apply_Limiter)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine rupdate_Classic_RK
!
!###############################################################################
!
subroutine rupdate_CK4_RK(nst,Apply_Filter,Apply_Limiter)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use ovar,    only : Limiter_Option,Filter_Option
  use flowvar, only : rssprk,residual,usp,dtsp
  use filter_mod, only : filter_usp
  use module_limiters, only : resolution_indicator
  use module_limiters, only : limit_by_projecting
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
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
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "rupdate_CK4_RK"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  do n = 1,n_solpts
    !
    rssprk(:,n,1) = alfa(nst)*rssprk(:,n,1) + dtsp(n)*residual(:,n)
    !
    usp(:,n) = usp(:,n) + beta(nst)*rssprk(:,n,1)
    !
  end do
  !
  ! Apply filtering or limiting to stabilize the solution
  !
  call stabilize_usp(nst,Apply_Filter,Apply_Limiter)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine rupdate_CK4_RK
!
!###############################################################################
!
subroutine rupdate_TVD2_RK(nst,Apply_Filter,Apply_Limiter)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use ovar,    only : Limiter_Option,Filter_Option
  use flowvar, only : residual,usp,uoldsp,dtsp
  use filter_mod, only : filter_usp
  use module_limiters, only : resolution_indicator
  use module_limiters, only : limit_by_projecting
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
  !
  !.. Local Scalars ..
  integer   :: n
  !
  ! TVD 2nd-order/2-stage Runge-Kutta Coefficients
  !
  real(wp), parameter, dimension(1:2) :: alfa = (/ zero,one  /)
  real(wp), parameter, dimension(1:2) :: beta = (/ one ,half /)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "rupdate_TVD2_RK"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  do n = 1,n_solpts
    !
    usp(:,n) = beta(nst) * (alfa(nst)*uoldsp(:,n) + &
                            usp(:,n) + &
                            dtsp(n)*residual(:,n))
    !
  end do
  !
  ! Apply filtering or limiting to stabilize the solution
  !
  call stabilize_usp(nst,Apply_Filter,Apply_Limiter)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine rupdate_TVD2_RK
!
!###############################################################################
!
subroutine rupdate_TVD3_RK(nst,Apply_Filter,Apply_Limiter)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use ovar,    only : Limiter_Option,Filter_Option
  use flowvar, only : residual,usp,uoldsp,dtsp
  use filter_mod, only : filter_usp
  use module_limiters, only : resolution_indicator
  use module_limiters, only : limit_by_projecting
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
  !
  !.. Local Scalars ..
  integer :: n
  !
  ! TVD 2nd-order/2-stage Runge-Kutta Coefficients
  !
  real(wp), parameter, dimension(1:3) :: alfa = (/ zero,three,half /)
  real(wp), parameter, dimension(1:3) :: beta = (/ one ,one4 ,two3 /)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "rupdate_TVD3_RK"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  do n = 1,n_solpts
    !
    usp(:,n) = beta(nst) * (alfa(nst)*uoldsp(:,n) + &
                            usp(:,n) + &
                            dtsp(n)*residual(:,n))
    !
  end do
  !
  ! Apply filtering or limiting to stabilize the solution
  !
  call stabilize_usp(nst,Apply_Filter,Apply_Limiter)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine rupdate_TVD3_RK
!
!###############################################################################
!
subroutine rupdate_TVD4_RK(nst,Apply_Filter,Apply_Limiter)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use ovar,    only : Limiter_Option,Filter_Option
  use flowvar, only : rssprk,residual,usp,uoldsp,usp_sum,dtsp
  use filter_mod, only : filter_usp
  use module_limiters, only : resolution_indicator
  use module_limiters, only : limit_by_projecting
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
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
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "rupdate_TVD4_RK"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
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
  ! Apply filtering or limiting to stabilize the solution
  !
  call stabilize_usp(nst,Apply_Filter,Apply_Limiter)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine rupdate_TVD4_RK
!
!###############################################################################
!
subroutine initialize_time_ave()
  !
  !.. Use Statements ..
  use geovar,  only : nr,n_solpts
  use ovar,    only : time,ave_start_time
  use ovar,    only : time_ave_vel_is_axisymm
  use flowvar, only : time_ave_variables
  use flowvar, only : uavesp
  !
  !.. Local Scalars ..
  integer :: l,m,nv,nvar,ierr
  character(len=1)  :: ta_char,v1,v2,v3
  character(len=30) :: array_name
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_time_ave"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Allocate and create the time_ave_variables array if it hasnt already
  ! been allocated; otherwise, count the number of time-averaged variables it
  ! contains. This character array identifies which time-averaged variables
  ! are computed and written to the time-averaged solution/restart files.
  !
  if (.not. allocated(time_ave_variables)) then
    !
    v1 = "x"
    v2 = merge( "r" , "y" , time_ave_vel_is_axisymm )
    v3 = merge( "f" , "z" , time_ave_vel_is_axisymm )
    !
    ! Since the time_ave_variables array isnt allocated, this means that a
    ! time-averaged restart file wasnt read and the start time of the
    ! time-averaging needs to be initialized to the current time
    !
    ave_start_time = time
    !
   !nvar = merge(23,11,nr==3)
   !allocate ( time_ave_variables(1:nvar) , stat=ierr , errmsg=error_message )
   !call alloc_error(pname,"time_ave_variables",1,__LINE__,__FILE__, &
   !                 ierr,error_message)
    !
    if (nr == 3) then
      !
      ! First 5 are the conserved variables
      !
      ! USING F2003 AUTO-REALLOCATION
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             "d" , "dx" , "dy" , "dz" , "e" ]
      !
      ! Next 5 are the primitive variables
      ! (with temperature instead of density)
      !
      ! USING F2003 AUTO-REALLOCATION
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , &
                             "t" , "x" , "y" , "z" , "p" ]
      !
      ! Remaining are derived variables
      !
      ! ALL OF THE BELOW ASSIGNMENTS USE F2003 AUTO-REALLOCATION TO
      ! APPEND A NEW ARRAY ELEMENT ONTO THE time_ave_variables ARRAY
      ! FOR EACH ADDITIONAL TIME AVERAGED VARIABLE
      !
      if (time_ave_vel_is_axisymm) then
        time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                               time_ave_variables(:) , "d" // v2 ]
        time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                               time_ave_variables(:) , "d" // v3 ]
        time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                               time_ave_variables(:) , v2 ]
        time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                               time_ave_variables(:) , v3 ]
      end if
      !
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , v1 // v1 ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , v2 // v2 ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , v3 // v3 ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , v1 // v2 ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , v1 // v3 ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , v2 // v3 ]
      !
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , v1 // "t" ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , v2 // "t" ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , v3 // "t" ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "dt" ]
      !
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "d" // v1 // v1 ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "d" // v2 // v2 ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "d" // v3 // v3 ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "d" // v1 // v2 ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "d" // v1 // v3 ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "d" // v2 // v3 ]
      !
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "d" // v1 // "t" ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "d" // v2 // "t" ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "d" // v3 // "t" ]
      !
      nvar = size(time_ave_variables)
      !
    else
      !
      ! First 4 are the conserved variables
      !
      ! USING F2003 AUTO-REALLOCATION
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             "d" , "dx" , "dy" , "e" ]
      !
      ! Next 4 are the primitive variables
      ! (with temperature instead of density)
      !
      ! USING F2003 AUTO-REALLOCATION
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "t" , "x" , "y" , "p" ]
      !
      ! Remaining are derived variables
      !
      ! USING F2003 AUTO-REALLOCATION
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "xx" ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "yy" ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "xy" ]
      !
      nvar = size(time_ave_variables)
      !
    end if
    !
  else
    !
    nvar = size(time_ave_variables)
    !
  end if
  !
  if (allocated(time_ave_var)) then
    deallocate ( time_ave_var , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"time_ave_var",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  allocate ( time_ave_var(1:nvar) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"time_ave_var",1,__LINE__,__FILE__,ierr,error_message)
  !
  do m = 1,size(time_ave_variables)
    !
    time_ave_variables(m) = trim(adjustl(time_ave_variables(m)))
    !
    nv = len_trim(time_ave_variables(m))
    !
    allocate ( time_ave_var(m)%idx(1:nv) , source=0 , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) m
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    do l = 1,nv
      !
      ta_char = uppercase( time_ave_variables(m)(l:l) )
      !
      select case (ta_char)
        case ("D")
          time_ave_var(m)%idx(l) = nec
        case ("X")
          time_ave_var(m)%idx(l) = nmx
        case ("Y")
          time_ave_var(m)%idx(l) = nmy
        case ("Z")
          time_ave_var(m)%idx(l) = nmz
        case ("E")
          time_ave_var(m)%idx(l) = nee
        case ("P")
          time_ave_var(m)%idx(l) = nee
        case ("T")
          time_ave_var(m)%idx(l) = nq+1
        case ("R")
          time_ave_var(m)%idx(l) = nq+2
        case ("F")
          time_ave_var(m)%idx(l) = nq+3
        case ("K")
          time_ave_var(m)%idx(l) = ntk
        case ("W")
          time_ave_var(m)%idx(l) = ntl
        case default
          ! ??? Should probably produce an error message and abort
          time_ave_var(m)%idx(l) = nec
      end select
    end do
  end do
  !
  ! Allocate the local module array used for computing
  ! the current update to the time-averaged variables.
  !
  if (allocated(ave)) then
    deallocate ( ave , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ave",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( ave(1:nvar) , source=zero , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"ave",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Allocate and initialize the time-averaged variables
  ! solution array if it hasnt already been allocated.
  !
  if (.not. allocated(uavesp)) then
    !
    allocate ( uavesp(1:nvar,1:n_solpts) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"uavesp",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Initialize the time-averaged variables to the initial/restart
    ! solution since they werent read from a time-averaged restart file.
    !
    call update_time_ave
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("time_ave_var(",i0,")%idx")
  !
end subroutine initialize_time_ave
!
!###############################################################################
!
subroutine update_time_ave()
  !
  !.. Use Statements ..
  use geovar,  only : xyz,n_solpts
  use flowvar, only : usp,uavesp
  use ovar,    only : time,prev_ave_time
  use ovar,    only : ave_start_time,saved_ave_time
  !
  !.. Local Scalars ..
  integer  :: m,n,nvar
  integer  :: nc1,nc2,nm1,nm2
  real(wp) :: a,b,y,z,st,ct,theta
  !
  !.. Local Arrays ..
  real(wp) :: time_val(1:2)
  real(wp) :: ave_var(1:nq+3)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_time_ave"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (time == zero .or. time == ave_start_time) then
    time_val(1) = zero
    time_val(2) = one
  else
    time_val(1) = (prev_ave_time-ave_start_time) / (time-ave_start_time)
    time_val(2) = one - time_val(1) ! to try and ensure that sum(time_val) = 1
  end if
  !
  saved_ave_time = prev_ave_time
  !
  ! Broadcast the time_val array so hopefully all MPI processes
  ! have the exact same values for the two averaging coefficients
  ! (in case there might be roundoff differences between processes).
  !
  call mpi_bcast(time_val,size(time_val,kind=int_mpi),mpi_flttyp, &
                 glb_root,MPI_COMM_WORLD,mpierr)
  !
  a = time_val(1)
  b = time_val(2)
  !
  nvar = size(time_ave_var)
  !
  ! Adjust the coordinate and momentum indices for conversion to cylindrical
  ! coordinates depending whether this is a 2D or 3D simulation
  !
  if (size(xyz,dim=1) == 2) then
    nc1 = 1; nc2 = 2; nm1 = nmx; nm2 = nmy
  else
    nc1 = 2; nc2 = 3; nm1 = nmy; nm2 = nmz
  end if
  !
  do n = 1,n_solpts
    !
    y = xyz(nc1,n)
    z = xyz(nc2,n)
    !
    if (all([y,z] == zero)) then
      !
      st = zero
      ct = zero
      !
    else
      !
      theta = atan2(z,y)
      if (theta < zero) theta = theta + two*pi
      !
      st = sin(theta)
      ct = cos(theta)
      !
    end if
    !
    ! Primitive variables
    ave_var(1:nq) = usp2v_sp( usp(1:nq,n) )
    !
    ! Temperature
    ave_var(nq+1) = temperature_pv_sp( ave_var(1:nq) )
    !
    ! Radial    instantaneous velocity
    ave_var(nq+2) = st*ave_var(nm2) + ct*ave_var(nm1)
    ! Azimuthal instantaneous velocity
    ave_var(nq+3) = ct*ave_var(nm2) - st*ave_var(nm1)
    !
    ! First, update the time-averaged conservative variables
    !
    do m = 1,nq
      uavesp(m,n) = a*uavesp(m,n) + b*usp(m,n)
    end do
    !
    ! Next, update the time-averaged primitive variables
    ! (with temperature instead of density)
    !
    uavesp(nq+1,n) = a*uavesp(nq+1,n) + b*ave_var(nq+1)
    do m = 2,nq
      uavesp(nq+m,n) = a*uavesp(nq+m,n) + b*ave_var(m)
    end do
    !
    ! Update the remaining time-averaged variables
    !
    do m = 2*nq+1,nvar
      uavesp(m,n) = a*uavesp(m,n) + b*product( ave_var(time_ave_var(m)%idx) )
    end do
    !
  end do
  !
  ! Update prev_ave_time to the current time
  !
  prev_ave_time = time
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine update_time_ave
!
!###############################################################################
!
include "Functions/usp2v.f90"
include "Functions/viscosity.f90"
include "Functions/temperature.f90"
!
!###############################################################################
!
end module time_mod
