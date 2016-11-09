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
continue
  !
 !call get_current_cfl
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
  !
continue
  !
  if (itcur < cfl_cycles) then
    cfl = cfl_beg + (cfl_end - cfl_beg) * real(itcur-1,kind=wp) / &
                                          real(cfl_cycles-1,kind=wp)
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
  use ovar,    only : dt_in,Final_Time,time,CFL
  use ovar,    only : this_is_final_timestep
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
      this_is_final_timestep = true
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
  use order_mod, only : n_order,maxSP
  !
  use geovar, only : nr,nfbnd,nface
  use geovar, only : face,cell
  !
  use ovar, only : gam,Timestep_Type,Final_Time,time,CFL
  use ovar, only : governing_equations
  use ovar, only : this_is_final_timestep
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
  integer  :: l1,l2,npl,kl
  integer  :: r1,r2,npr,kr,face_order
  integer  :: host_cell,host_geom,host_order
  integer  :: left_cell,left_geom,left_order
  integer  :: right_cell,right_geom,right_order
  real(wp) :: ds,eigni,eignv,pp1sq
  real(wp) :: dsl,vell,al,jacl,dtl,mul,hl
  real(wp) :: dsr,velr,ar,jacr,dtr,mur,hr
  real(wp) :: psql,psqr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nq)    :: vl,vr
  !
  integer, parameter :: make_p_plus_one = 0
 !integer, parameter :: make_p_plus_one = 1
  !
  logical(lk), parameter :: use_sqrt = true
 !logical(lk), parameter :: use_sqrt = fals
  !
continue
  !
  ! ######################################
  ! ######################################
  ! ######################################
  ! ######################################
  ! ##########                  ##########
  ! ##########   FROM HPMUSIC   ##########
  ! ##########   ||||||||||||   ##########
  ! ##########   VVVVVVVVVVVV   ##########
  ! ##########                  ##########
  ! ######################################
  ! ######################################
  ! ######################################
  ! ######################################
  !
  !void CPR_NS_2D_Solver<T>::ComputeLocalTimeStep(Real cfl, Real *dt_cell){
  !
  !  for(int i=0;i<n_cell_;++i){
  !    Real factor=0.;
  !    ElemType type=this->grid_->elem_type[i];
  !    for(int side=0; side<NumSideFromType(type); ++side){
  !      int face = this->grid_->elem_to_face[i][side];
  !      FaceData *face_data=&this->grid_->face[face];
  !      int c1 = face_data->left_cell;
  !      int c2 = face_data->right_cell;
  !
  !      Real dx,dy;
  !      if(face_data->face_type==BOUNDARY){
  !        dx=x_fc_[face]-x_cc_[c1];
  !        dy=y_fc_[face]-y_cc_[c1];
  !      }else{
  !        dx=x_cc_[c2]-x_cc_[c1];
  !        dy=y_cc_[c2]-y_cc_[c1];
  !      }
  !
  !      Real u=0.,v=0., c=0.,mu_k=0.,nx=0.,ny=0.,area=0.,p=0.;
  !      T *q=&u_fp_[i][0];
  !      int nfp_f=n_lfp_f_[face];
  !      for(int j=0;j<nfp_f; ++j){
  !        int ir=j*nv_;
  !        int ix=ir+1;
  !        int iy=ir+2;
  !        int ie=ir+3;
  !        T rhoi = 1./q[ir];
  !        u = u+q[ix]*rhoi;
  !        v = v+q[iy]*rhoi;
  !        p = (gas_.gamma-1.)*
  !            (q[ie]-0.5*(q[ix]*q[ix]+q[iy]*q[iy])*rhoi);//? why not +=
  !        c = c+sqrt(gas_.gamma*p*rhoi);
  !        mu_k = mu_k+gas_.mu*rhoi;
  !        nx+=fm_fp_.nx_[face][j];
  !        ny+=fm_fp_.ny_[face][j];
  !        area+=fm_fp_.area_[face][j];
  !      }
  !      u /= nfp_f;
  !      v /= nfp_f;
  !      c /= nfp_f;
  !      mu_k /= nfp_f;
  !      nx/=  nfp_f;
  !      ny/=  nfp_f;
  !      area/= nfp_f;
  !
  !      factor += (fabs(u*nx+v*ny)+c+2.*mu_k/fabs(dx*nx+dy*ny))*area;
  !    }
  !
  !    dt_cell[i] = cfl*vol_[i]/factor;
  !    if(dt_cell[i]!=dt_cell[i]){
  !      stringstream ss;
  !      ss<<"Local time step is NaN at cell "<<i<<" !!!"<<endl;
  !      fatal_error(ss.str());
  !    }
  !  }//end of cell i
  !}
  !
  ! ######################################
  ! ######################################
  ! ######################################
  ! ######################################
  ! ##########                  ##########
  ! ##########   AAAAAAAAAAAA   ##########
  ! ##########   ||||||||||||   ##########
  ! ##########   FROM HPMUSIC   ##########
  ! ##########                  ##########
  ! ######################################
  ! ######################################
  ! ######################################
  ! ######################################
  !
  mul = zero
  mur = zero
  !
 !pp1sq = real( (n_order+1)**2 , kind=wp )
 !pp1sq = real( (n_order+1)**nr , kind=wp )
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
    face_order = face(nf)%order ! face solution order
    !
    psql = real( (host_order+make_p_plus_one)**2 , kind=wp )
    !
    do k = 1,size( face(nf)%left%flxpts )
      !
      kl = face(nf)%left%flxpts(k)
      !
      vl(1:nq) = usp2v_sp( faceusp(nf)%v(1:nq,k,1) )
      !
      if (governing_equations == NavierStokes_Eqns) then
        mul = viscosity_cv_sp( faceusp(nf)%v(1:nq,k,1) )
      end if
      !
      ds = rnorm2( geofa(nf)%v(1:nr,k) )
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
             dot( kl , jacobian(l1:l2) )
      !
      hl = half * dsl / jacl
      !
      eigni = abs(vell) + al
     !eignv = pp1sq * hl * mul
      eignv = psql * hl * mul
      !
      ! Compute the timestep for this flux point
      !
     !dtl = min( dtsp(l1) , CFL / (pp1sq * hl * (eigni + eignv)) )
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
    face_order = face(nf)%order ! face solution order
    !
    psql = real( (left_order +make_p_plus_one)**2 , kind=wp )
    psqr = real( (right_order+make_p_plus_one)**2 , kind=wp )
    !
    do k = 1,size( face(nf)%left%flxpts )
      !
      kl = face(nf)%left%flxpts(k)
      kr = face(nf)%right%flxpts(k)
      !
      vl(1:nq) = usp2v_sp( faceusp(nf)%v(1:nq,k,1) )
      vr(1:nq) = usp2v_sp( faceusp(nf)%v(1:nq,k,2) )
      !
      if (governing_equations == NavierStokes_Eqns) then
        mul = viscosity_cv_sp( faceusp(nf)%v(1:nq,k,1) )
        mur = viscosity_cv_sp( faceusp(nf)%v(1:nq,k,2) )
      end if
      !
      ds = rnorm2( geofa(nf)%v(1:nr,k) )
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
             dot( kl , jacobian(l1:l2) )
      !
      ! Interpolate the jacobian to the flux point for the right cell
      !
      jacr = interp(right_geom,right_order)% &
                         toFace(face_order)% &
             dot( kr , jacobian(r1:r2) )
      !
      hl = half * dsl / jacl
      hr = half * dsr / jacr
      !
      eigni = max( abs(vell) + al , abs(velr) + ar )
     !eignv = pp1sq * max( hl * mul , hr * mur )
      eignv = max( psql * hl * mul , psqr * hr * mur )
      !
      ! Compute the timestep for this flux point
      !
     !dtl = min( dtsp(l1) , CFL / (pp1sq * hl * (eigni + eignv) ) )
      dtl = min( dtsp(l1) , CFL / (psql * hl * (eigni + eignv) ) )
     !dtr = min( dtsp(r1) , CFL / (pp1sq * hr * (eigni + eignv) ) )
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
    ! Make sure the timestep size doesnt make the final time exceed Final_Time
    !
    if (Final_Time > zero) then
      if (time + dtmin > Final_Time + eps9) then
        dtsp(:) = Final_Time - time
        this_is_final_timestep = true
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
  use geovar,  only : nr,ncell,cell
  use ovar,    only : gam,Final_Time,time,CFL,Timestep_Type
  use ovar,    only : this_is_final_timestep
  use ovar,    only : governing_equations
  use metrics_mod,  only : jacobian,metrics
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
  integer, parameter :: make_p_plus_one = 0
 !integer, parameter :: make_p_plus_one = 1
  !
continue
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
      vel(1:nr) = matmul( metrics(1:nr,1:nr,n,2) , pv(nmb:nme) ) / jacobian(n)
      !
      vmag = rnorm2(vel)
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
      if (time + dtmin > Final_Time + eps9) then
        dtsp(:) = Final_Time - time
        this_is_final_timestep = true
      end if
    end if
    !
  end if
  !
end subroutine get_point_timestep
!
!###############################################################################
!
subroutine old_get_point_timestep(dtmin)
  !
  !.. Use Statements ..
  use order_mod,  only : n_order
  use geovar,  only : nr,n_solpts
  use ovar,    only : gam,Final_Time,time,CFL,Timestep_Type
  use ovar,    only : this_is_final_timestep
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
  real(wp) :: alfa,beta,gama
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
    dxdr = metrics(1,1,n,1)
    dydr = metrics(2,1,n,1)
    dxds = metrics(1,2,n,1)
    dyds = metrics(2,2,n,1)
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
      if (time + dtmin > Final_Time + eps9) then
        dtsp(:) = Final_Time - time
        this_is_final_timestep = true
      end if
    end if
    !
  end if
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
  use geovar,  only : n_solpts
  use geovar,  only : ncell,cell
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
  if (p_order < n_order) then
    call project_to_porder(cell,residual)
  end if
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
    usp(:,n) = beta(nst) * (alfa(nst)*uoldsp(:,n) + &
                            usp(:,n) + &
                            dtsp(n)*residual(:,n))
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
    usp(:,n) = beta(nst) * (alfa(nst)*uoldsp(:,n) + &
                            usp(:,n) + &
                            dtsp(n)*residual(:,n))
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
