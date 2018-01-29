module module_bc
  !
#include <mpi_defs.h>
  !
  ! Module to contain the procedures and variables related
  ! to the evaluating the boundary conditions for the flow.
  !
  !.. Use Statements ..
  use module_kind_types
  !
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq,ntk,ntl
  use ovar,    only : gam
  !
  implicit none
  !
  private
  !
  !
  !
  ! ##################################
  ! ###  Module Public Procedures  ###
  ! ##################################
  !
  public :: get_faceusp_bc
  public :: get_facedusp_bc
  !
  !
  !
  ! ###########################
  ! ###  Module Parameters  ###
  ! ###########################
  !
  !
  !
  ! ###################################
  ! ###  Module Allocatable Arrays  ###
  ! ###################################
  !
  !
  !
  ! ################################
  ! ###  Module Local Variables  ###
  ! ################################
  !
  !
contains
!
!###############################################################################
!
subroutine get_faceusp_bc(exit_error)
  !
  ! Routine to compute the conservative variables at the boundary flux points
  ! and store the results in the faceusp array. This should be called after
  ! the interp_cv_to_flxpts routine so the internal states at the boundary
  ! flux points should already be defined.  We use these internal states
  ! to compute the boundary states based on the face boundary condition.
  !
  !.. Use Statements ..
  use geovar,       only : bface,nfbnd
  use ovar,         only : output_bface_array
  use ovar,         only : bc_in,sub_inflow_method
  use metrics_mod,  only : geofa
  use flowvar,      only : faceusp,facexyz
  !
  use parallel_mod, only : wait_faceusp
  !
  !.. Optional Arguments ..
  integer, optional, intent(out) :: exit_error
  !
  !.. Local Scalars ..
  integer :: nf,pf,bc_idx
  integer :: newton_failures
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_faceusp_bc"
  !
 !logical(lk), parameter :: use_old_iso_wall = true
  logical(lk), parameter :: use_old_iso_wall = fals
  !
continue
  !
  newton_failures = 0
  if (present(exit_error)) exit_error = 0
  !
  if (output_bface_array) call output_bface
  !
  ! Loop over the boundary faces and compute the conservative variables
  ! at the boundary flux points based on the face boundary condition.
  !
  do nf = 1,nfbnd
    !
    bc_idx = bface(6,nf)
    !
    ! Compute the solution at the boundary flux points
    ! based on the specified boundary condition
    !
    select case (bface(1,nf))
      !
      case (bc_freestream)
        !
        ! Free stream boundary condition
        !
        call get_bc_freestream_solution(faceusp(nf)%v(:,:,2), &
                                        bc_in(0)%cv)
        !
      case (bc_fixed)
        !
        ! Free stream boundary condition
        !
        call get_bc_freestream_solution(faceusp(nf)%v(:,:,2), &
                                        bc_in(bc_idx)%cv)
        !
      case (bc_characteristic)
        !
        ! Characteristic boundary condition
        !
       !call get_bc_characteristic_solution(faceusp(nf)%v(:,:,1), &
       !                                    faceusp(nf)%v(:,:,2), &
       !                                    geofa(nf)%v(:,:), &
       !                                    bc_in(bc_idx)%cv, &
       !                                    bc_in(bc_idx)%pv)
        call get_bc_characteristic_solution(faceusp(nf)%v(:,:,1), &
                                            faceusp(nf)%v(:,:,2), &
                                            geofa(nf)%v(:,:), &
                                            bc_in(bc_idx)%pv, &
                                            bc_in(bc_idx)%relax)
        !
      case (bc_generic_inflow)
        !
        ! Generic inflow boundary condition
        !
        call get_bc_generic_inflow_solution(faceusp(nf)%v(:,:,1), &
                                            faceusp(nf)%v(:,:,2), &
                                            geofa(nf)%v(:,:), &
                                            bc_in(bc_idx)%cv, &
                                            bc_in(bc_idx)%pv)
        !
      case (bc_generic_outflow)
        !
        ! Generic outflow boundary condition
        !
        call get_bc_generic_outflow_solution(faceusp(nf)%v(:,:,1), &
                                             faceusp(nf)%v(:,:,2), &
                                             geofa(nf)%v(:,:), &
                                             bc_in(bc_idx)%pv)
        !
      case (bc_generic_freeflow)
        !
        ! Generic freeflow boundary condition
        !
        call get_bc_generic_freeflow_solution(faceusp(nf)%v(:,:,1), &
                                              faceusp(nf)%v(:,:,2), &
                                              geofa(nf)%v(:,:), &
                                              bc_in(bc_idx)%cv, &
                                              bc_in(bc_idx)%pv)
        !
      case (bc_sub_inflow)
        !
        ! Subsonic inflow boundary condition
        !
        if (any(sub_inflow_method == [2,3,4])) then
          call get_bc_sub_inflow_solution(faceusp(nf)%v(:,:,1), &
                                          faceusp(nf)%v(:,:,2), &
                                          bc_in(bc_idx)%cv, &
                                          bc_in(bc_idx)%pv)
        else
          call get_bc_sub_inflow_newton(faceusp(nf)%v(:,:,1), &
                                        faceusp(nf)%v(:,:,2), &
                                        geofa(nf)%v(:,:), &
                                        bc_in(bc_idx)%pv, &
                                        newton_failures, &
                                        nf)
        end if
        !
      case (bc_sub_outflow)
        !
        ! Subsonic outflow boundary condition
        !
        call get_bc_sub_outflow_solution(faceusp(nf)%v(:,:,1), &
                                         faceusp(nf)%v(:,:,2), &
                                         geofa(nf)%v(:,:), &
                                         bc_in(bc_idx)%pv)
        !
      case (bc_sup_inflow)
        !
        ! Supersonic inflow boundary condition
        !
        call get_bc_sup_inflow_solution(faceusp(nf)%v(:,:,2), &
                                        bc_in(bc_idx)%cv)
        !
      case (bc_sup_outflow)
        !
        ! Supersonic outflow boundary condition
        !
        call get_bc_sup_outflow_solution(faceusp(nf)%v(:,:,1), &
                                         faceusp(nf)%v(:,:,2))
        !
      case (bc_slip_wall,bc_symmetry)
        !
        ! Slip wall / symmetry boundary condition
        !
        call get_bc_slip_wall_solution(faceusp(nf)%v(:,:,1), &
                                       faceusp(nf)%v(:,:,2), &
                                       geofa(nf)%v(:,:))
        !
      case (bc_adiabatic_wall)
        !
        ! No-slip adiabatic wall
        !
        call get_bc_adiabatic_wall_solution(faceusp(nf)%v(:,:,1), &
                                            faceusp(nf)%v(:,:,2), &
                                            bc_in(bc_idx)%pv)
        !
      case (bc_isothermal_wall)
        !
        ! No-slip isothermal wall
        !
        if (use_old_iso_wall) then
          call get_bc_iso_wall_solution_old(faceusp(nf)%v(:,:,1), &
                                            faceusp(nf)%v(:,:,2), &
                                            bc_in(bc_idx)%wall_temperature, &
                                            geofa(nf)%v(:,:))
        else
          call get_bc_isothermal_wall_solution(faceusp(nf)%v(:,:,1), &
                                               faceusp(nf)%v(:,:,2), &
                                               bc_in(bc_idx)%wall_temperature, &
                                               bc_in(bc_idx)%pv)
        end if
        !
      case (bc_periodic)
        !
        ! Periodic boundary condition
        !
        pf = bface(5,nf) ! index of connecting periodic boundary face
        !
        call get_bc_periodic_solution(faceusp(pf)%v(:,:,1), &
                                      faceusp(nf)%v(:,:,2))
        !
      case (bc_mms_dirichlet)
        !
        ! Dirichlet boundary condition for MMS
        !
        call get_bc_mms_dirichlet_solution(faceusp(nf)%v(:,:,1), &
                                           faceusp(nf)%v(:,:,2), &
                                           facexyz(nf)%v(:,:), &
                                           geofa(nf)%v(:,:))
        !
      case (bc_cpu_bnd)
        !
        ! We dont need to do anything here for communication boundaries
        !
      case default
        !
        ! This shouldn't be here and is an unknown boundary condition
        !
        write (iout,1) bface(1,nf)
        call stop_gfr(abort,pname,__LINE__,__FILE__)
        !
    end select
    !
  end do
  !
  if (ncpu > 1) then
    call mpi_allreduce(MPI_IN_PLACE,newton_failures,1_int_mpi, &
                       mpi_inttyp,MPI_SUM,MPI_COMM_WORLD,mpierr)
  end if
  !
  if (newton_failures > 0) then
    if (present(exit_error)) then
      exit_error = newton_failures
    else
      write (error_message,2) newton_failures
      call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
    end if
  end if
  !
  ! Now that we have finished getting the boundary conditions local to
  ! this cell, make sure that the communication of faceusp between
  ! processors has completed before continuing onwards.
  !
  if (ncpu > 1) call wait_faceusp(faceusp)
  !
  call exchange_bc_conflicts_usp(faceusp)
  !
  ! Format Statements
  !
  1 format (/," ********************* ERROR! **********************",/, &
              "  UNKNOWN BOUNDARY CONDITION IN SUB GET_FACEUSP_BC!",/, &
              "  BOUNDARY CONDITION VALUE IS: ",i0,/, &
              " ********************* ERROR! **********************",/)
  2 format (" FAILURES WERE DETECTED USING NEWTON'S METHOD", &
            " TO COMPUTE THE SUBSONIC INFLOW TEMPERATURE!", &
            "       newton_failures = ",i0)
  !
end subroutine get_faceusp_bc
!
!###############################################################################
!
subroutine get_facedusp_bc
  !
  ! Routine to compute the derivatives of the conservative variables at the
  ! boundary flux points and store the results in the facedusp array. This
  ! should be called after the interp_dusp_to_flxpts routine so the internal
  ! gradients at the boundary flux points should already be defined.  We use
  ! these internal gradients to compute the boundary gradients based on the
  ! face boundary condition.
  !
  !.. Use Statements ..
  use geovar,       only : bface,nfbnd
  use flowvar,      only : faceusp,facedusp
  use metrics_mod,  only : geofa
  !
  use parallel_mod, only : wait_facedusp
  !
  !.. Local Scalars ..
  integer :: nf,pf
  !
  !.. Local Parameters ..
  integer, parameter :: local_bc_ignore(1:3) = [ bc_cpu_bnd, &
                                                 bc_unknown, &
                                                 not_a_bc ]
  character(len=*), parameter :: pname = "get_facedusp_bc"
  !
continue
  !
  ! Loop over the boundary faces and compute the conservative variables
  ! at the boundary flux points based on the face boundary condition.
  !
  do nf = 1,nfbnd
    !
    ! Compute the solution at the boundary flux points
    ! based on the specified boundary condition
    !
    if (any(bface(1,nf) == bc_walls)) then
      !
      ! Wall boundary condition
      !
      if (bface(1,nf) == bc_adiabatic_wall) then
        !
        ! No-slip adiabatic wall boundary conditions
        !
        call get_bc_adiabatic_wall_gradient(faceusp(nf)%v(:,:,1), &
                                            facedusp(nf)%v(:,:,:,1), &
                                            facedusp(nf)%v(:,:,:,2))
        !
      else if (any(bface(1,nf) == [bc_slip_wall,bc_symmetry])) then
        !
        ! Slip wall / symmetry boundary condition
        !
        call get_bc_slip_wall_gradient(facedusp(nf)%v(:,:,:,1), &
                                       facedusp(nf)%v(:,:,:,2), &
                                       geofa(nf)%v(:,:))
        !
      else
        !
        ! Non-adiabatic wall boundary condition
        !
        call get_bc_wall_gradient(faceusp(nf)%v(:,:,1), &
                                  facedusp(nf)%v(:,:,:,1), &
                                  facedusp(nf)%v(:,:,:,2))
        !
      end if
      !
    else if (bface(1,nf) == bc_periodic) then
      !
      ! Periodic boundary condition
      !
      pf = bface(5,nf) ! index of connecting periodic boundary face
      !
      ! We just want to copy the gradients from face pf to face nf which is
      ! exactly what the subroutine get_bc_flowthrough_gradient does, so
      ! just use that instead of making a completely separate subroutine.
      !
      call get_bc_periodic_gradient(facedusp(pf)%v(:,:,:,1), &
                                    facedusp(nf)%v(:,:,:,2))
      !
    else if (all(bface(1,nf) /= local_bc_ignore)) then
      !
      ! General flow-through boundary conditions
      !
      call get_bc_flowthrough_gradient(facedusp(nf)%v(:,:,:,1), &
                                       facedusp(nf)%v(:,:,:,2))
      !
   !else
   !  !
   !  ! This shouldn't be here and is an unknown boundary condition
   !  !
   !  write (iout,1) bface(1,nf)
   !  call stop_gfr(abort,pname,__LINE__,__FILE__)
   !  !
    end if
    !
  end do
  !
  ! Now that we have finished getting the boundary conditions local to
  ! this cell, make sure that the communication of facedusp between
  ! processors has completed before continuing onwards.
  !
  if (ncpu > 1) call wait_facedusp(facedusp)
  !
  call exchange_bc_conflicts_dusp(facedusp)
  !
  ! Format Statements
  !
  1 format (/," ********************** ERROR! **********************",/, &
              "  UNKNOWN BOUNDARY CONDITION IN SUB GET_FACEDUSP_BC!",/, &
              "  BOUNDARY CONDITION VALUE IS: ",i0,/, &
              " ********************** ERROR! **********************",/)
  !
end subroutine get_facedusp_bc
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
! Functions that compute the boundary condition
! solution flow variables in conservative form
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine get_bc_freestream_solution(gufp,cv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Free stream boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! GUFP  : Conservative variables at the ghost cell flux points
  ! CV_IN : Boundary condition flow variables in conservative form
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: cv_in(:)
  real(wp), intent(inout) :: gufp(:,:)
  !
  !.. Local Scalars ..
  integer :: k
  !
continue
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Set all variables to free stream
    !
    gufp(:,k) = cv_in(:)
    !
  end do
  !
end subroutine get_bc_freestream_solution
!
!###############################################################################
!
subroutine get_bc_characteristic_solution_old(hufp,gufp,fnrm,cv_in,pv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Characteristic boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FNRM  : Face outward normal vector at the flux points
  ! HUFP  : Conservative variables at the host  cell flux points
  ! GUFP  : Conservative variables at the ghost cell flux points
  ! CV_IN : Boundary condition flow variables in conservative form
  ! PV_IN : Boundary condition flow variables in primitive    form
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: fnrm(:,:)
  real(wp),    intent(in) :: hufp(:,:)
  real(wp), intent(inout) :: gufp(:,:)
  real(wp),    intent(in) :: cv_in(:)
  real(wp),    intent(in) :: pv_in(:)
  !
  !.. Local Scalars ..
  integer  :: k
  real(wp) :: rp,rm,utc,unc,ac,sc,rho,pres
  real(wp) :: mach_int,gm1,gm1_inv
  real(wp) :: aspd_int,vt_int,vn_int
  real(wp) :: aspd_ext,vt_ext,vn_ext
  !
  !.. Local Arrays ..
  real(wp) :: rnrm(1:size(fnrm,dim=1))
  real(wp) :: rtan(1:size(fnrm,dim=1))
  real(wp) :: hvfp(1:size(gufp,dim=1))
  !
continue
  !
  ! Precompute gamma-1 and its inverse
  !
  gm1 = gam - one
  gm1_inv = one / gm1
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Inward pointing unit normal at this flux point
    !
    rnrm(:) = -unit_vector( fnrm(:,k) )
    !
    ! Get the host cell primitive variables at this flux point
    !
    hvfp(:) = usp2v_sp( hufp(:,k) )
    !
    ! Compute the speed of sound in the host cell
    !
    aspd_int = sqrt( max( gam*hvfp(nee)/hvfp(nec) , rndoff ) )
    !
    ! Compute the interior normal velocity at this flux point
    !
    vn_int = dot( rnrm , hvfp(nmb:nme) )
    !
    ! Apply boundary conditions based on the
    ! normal Mach number to the boundary face
    !
    mach_int = vn_int/aspd_int
    !
    if (mach_int >= one) then ! Supersonic inflow
      !
      ! Set ghost cell flux point to free stream conditions
      !
      gufp(:,k) = cv_in(:)
      !
    else if (mach_int <= -one) then ! Supersonic outflow
      !
      ! Extrapolate interior to ghost cell flux point
      !
      gufp(:,k) = hufp(:,k)
      !
    else ! Subsonic inflow/outflow using Riemann invariants
      !
      ! Unit tangent at this flux point computed from unit normal
      ! NOTE: Usually t = (tx,ty) = (ny,-nx) but we computed the inward
      !       pointing unit normal, so here (tx,ty) = (-ny,nx).
      !
      rtan(:) = [-rnrm(2),rnrm(1)]
      !
      ! Compute the exterior speed of sound
      !
      aspd_ext = sqrt( max( gam*pv_in(nee)/pv_in(nec) , rndoff ) )
      !
      ! Compute the exterior normal velocity at this flux point
      !
      vn_ext = dot( rnrm , pv_in(nmb:nme) )
      !
     !if (mach_int > zero) then ! this is an inflow boundary
     !  !
     !  rp = vn_ext + two*aspd_ext*gm1_inv
     !  rm = vn_int - two*aspd_int*gm1_inv
     !  unc = half*(rp+rm)
     !  utc = dot( rtan , pv_in(nmb:nme) )
     !  ac = one4*gm1*(rp-rm)
     !  sc = pv_in(nee)/(pv_in(nec)**gam)
     !  !
     !else ! this is an outflow boundary
     !  !
     !  rp = -abs(vn_ext) + two*aspd_ext*gm1_inv
     !  rm =      vn_int  - two*aspd_int*gm1_inv
     !  unc = half*(rp+rm)
     !  utc = dot( rtan , hvfp(nmb:nme) )
     !  ac = one4*gm1*(rp-rm)
     !  sc = hvfp(nee)/(hvfp(nec)**gam)
     !  !
     !end if
      !
      ! Compute the variables dependent on the interior flow direction
      !
      if (mach_int > zero) then
        rp   = vn_ext
        utc  = dot( rtan , pv_in(nmb:nme) )
        rho  = pv_in(nec)
        pres = pv_in(nee)
      else
        rp   = -abs(vn_ext)
        utc  = dot( rtan , hvfp(nmb:nme) )
        rho  = hvfp(nec)
        pres = hvfp(nee)
      end if
      !
      ! Compute the Riemann invariants
      !
      rp  = rp     + two*aspd_ext*gm1_inv
      rm  = vn_int - two*aspd_int*gm1_inv
      unc = half*(rp+rm)
      ac  = one4*gm1*(rp-rm)
      sc  = pres/(rho**gam)
      !
      ! Compute the density and pressure from the Riemann invariants
      !
      rho = ( (ac*ac) / (gam*sc) )**(gm1_inv)
      pres = sc * (rho**gam)
      !
      ! Compute the conservative variables at the boundary flux points
      !
      gufp(nec,k) = rho
      gufp(nmb:nme,k) = rho * (unc*rnrm(:) + utc*rtan(:))
      gufp(nee,k) = pres*gm1_inv + half * selfdot( gufp(nmb:nme,k) ) / rho
      !
    end if
    !
  end do
  !
end subroutine get_bc_characteristic_solution_old
!
!###############################################################################
!
subroutine get_bc_characteristic_solution(hufp,gufp,fnrm,pv_in,relax)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Characteristic boundary condition
  ! #####################################################
  ! #####################################################
  !
  !.. Use Statements ..
  use ovar, only : relaxation_t
  use ovar, only : itcur,time
  !
  ! FNRM  : Face outward normal vector at the flux points
  ! HUFP  : Conservative variables at the host  cell flux points
  ! GUFP  : Conservative variables at the ghost cell flux points
  ! PV_IN : Boundary condition flow variables in primitive    form
  ! RELAX : Relaxation conditions for this boundary condition
  !
  !.. Formal Arguments ..
  real(wp),           intent(in) :: fnrm(:,:)
  real(wp),           intent(in) :: hufp(:,:)
  real(wp),           intent(in) :: pv_in(:)
  type(relaxation_t), intent(in) :: relax
  real(wp),        intent(inout) :: gufp(:,:)
  !
  !.. Local Scalars ..
  integer  :: k
  real(wp) :: rp,rm,unc,ac,sc,rho,pres,vel
  real(wp) :: mach_int,mach_rlx,gm1,gm1_inv
  real(wp) :: aspd_int,vt_int,vn_int
  real(wp) :: aspd_rlx,vt_rlx,vn_rlx
  real(wp) :: rlx,rlx_iter,rlx_time
  !
  !.. Local Arrays ..
  real(wp) :: rnrm(1:size(fnrm,dim=1))
  real(wp) :: rtan(1:size(fnrm,dim=1))
  real(wp) :: hvfp(1:size(gufp,dim=1))
  real(wp) :: pv_in_rlx(1:size(pv_in))
  !
  !.. Local Parameters ..
  logical(lk), parameter :: veln_is_int = true
 !logical(lk), parameter :: veln_is_int = fals
  !
continue
  !
  ! Precompute gamma-1 and its inverse
  !
  gm1 = gam - one
  gm1_inv = one / gm1
  !
  ! Determine the relaxation coefficient
  !
  rlx = one
  rlx_iter = one
  rlx_time = one
  if (relax%iter > 0) then
    rlx_iter = max( one , real(itcur,kind=wp)/real(relax%iter,kind=wp) )
  end if
  if (relax%time > zero) then
    rlx_time = max( one , time/relax%time )
  end if
  rlx = max( zero , min(relax%value,rlx_iter,rlx_time,one) )
  !
  ! Compute the boundary solution at each flux point
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Inward pointing unit normal at this flux point
    !
    rnrm(:) = unit_vector( fnrm(:,k) )
    !
    ! Get the host cell primitive variables at this flux point
    !
    hvfp(:) = usp2v_sp( hufp(:,k) )
    !
    ! Compute the relaxed inflow condition
    !
    pv_in_rlx(:) = (one-rlx)*hvfp(:) + rlx*pv_in(:)
    !
    ! Compute the speed of sound in the host cell
    !
    ! interior speed of sound
    aspd_int = sqrt( max( gam*hvfp(nee)/hvfp(nec) , rndoff ) )
    ! exterior speed of sound
    aspd_rlx = sqrt( max( gam*pv_in_rlx(nee)/pv_in_rlx(nec) , rndoff ) )
    !
    ! Compute the interior normal velocity at this flux point
    !
    vn_int = dot( rnrm , hvfp(nmb:nme) )  ! interior normal velocity
    vn_rlx = dot( rnrm , pv_in_rlx(nmb:nme) ) ! exterior normal velocity
    !
    ! Apply boundary conditions based on the
    ! normal Mach number to the boundary face
    !
    rp = vn_int + two*aspd_int*gm1_inv ! Outgoing Riemann invariant
    rm = vn_rlx - two*aspd_rlx*gm1_inv ! Incoming Reimann invariant
    !
    if (vn_rlx <= -aspd_rlx) then
      ! Supersonic inflow boundary => Set R+ so that 0.5*(R+ + R-) = vn_rlx
      rp = vn_rlx + two*aspd_rlx*gm1_inv
    end if
    !
    if (vn_int >= aspd_int) then
      ! Supersonic outflow boundary => Set R- so that 0.5*(R+ + R-) = vn_int
      rm = vn_int - two*aspd_int*gm1_inv
    end if
    !
    unc = half*(rp+rm)
    ac  = one4*gm1*(rp-rm)
    !
    vel = merge( vn_int , unc , veln_is_int )
    !
    if (vel > zero) then ! outflow boundary
      gufp(nmb:nme,k) = hvfp(nmb:nme) + (unc - vn_int)*rnrm(:)
      sc = hvfp(nee) / (hvfp(nec)**gam)
    else
      gufp(nmb:nme,k) = pv_in_rlx(nmb:nme) + (unc - vn_rlx)*rnrm(:)
      sc = pv_in_rlx(nee) / (pv_in_rlx(nec)**gam)
    end if
    !
    rho = ( ac * ac / (gam*sc) )**gm1_inv
    pres = sc * (rho**gam)
    !
    gufp(nec,k) = rho
    gufp(nee,k) = pres
    !
    ! Convert the primitive variables to conservative variables
    !
    gufp(1:nq,k) = vsp2u_sp( gufp(1:nq,k) )
    !
  end do
  !
end subroutine get_bc_characteristic_solution
!
!###############################################################################
!
subroutine get_bc_generic_inflow_solution(hufp,gufp,fnrm,cv_in,pv_in)
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
  ! FNRM  : Face outward normal vector at the flux points
  ! HUFP  : Conservative variables at the host  cell flux points
  ! GUFP  : Conservative variables at the ghost cell flux points
  ! CV_IN : Boundary condition flow variables in conservative form
  ! PV_IN : Boundary condition flow variables in primitive    form
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: fnrm(:,:)
  real(wp),    intent(in) :: hufp(:,:)
  real(wp),    intent(in) :: cv_in(:)
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: gufp(:,:)
  !
  !.. Local Scalars ..
  integer  :: nr,k
  real(wp) :: gm1_inv,ske_inf,mach_int,aspd_int
  !
  !.. Local Arrays ..
  real(wp) :: rnrm(1:size(fnrm,dim=1))
  real(wp) :: hvfp(1:size(gufp,dim=1))
  !
continue
  !
  nr = size(fnrm,dim=1)
  !
  ! Precompute the free stream specific kinetic energy
  ! and the inverse of gamma-1
  !
  ske_inf = half * pv_in(nec) * selfdot( pv_in(nmb:nme) )
  !
  gm1_inv = one / (gam - one)
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Unit normal and tangent at this flux point
    !
    rnrm(:) = unit_vector( fnrm(:,k) )
    !
    ! Get the host cell primitive variables at this flux point
    !
    hvfp(:) = usp2v_sp( hufp(:,k) )
    !
    ! Compute the speed of sound at this flux point
    !
    aspd_int = sqrt( max( gam*hvfp(nee)/hvfp(nec) , rndoff ) )
    !
    ! Compute the Mach number of the normal velocity at this flux point
    !
    mach_int = dot( hvfp(nmb:nme) , rnrm(1:nr) ) / aspd_int
    !
    if (mach_int > -one .and. mach_int <= zero) then
      !
      ! Subsonic inflow
      !
      ! Set all variables to free stream
      !
      gufp(:,k) = cv_in(:)
      !
      ! Compute the total energy using interior pressure
      ! and free stream specific kinetic energy
      !
      gufp(nee,k) = hvfp(nee) * gm1_inv + ske_inf
      !
    else if (mach_int <= -one) then
      !
      ! Supersonic inflow
      !
      ! Set all variables to free stream
      !
      gufp(:,k) = cv_in(:)
      !
    else
      !
      ! The normal velocity at this flux point seems to
      ! indicate this should be an outflow boundary.
      !
      write (iout,1) mach_int,aspd_int,hvfp(nee),hvfp(nec), &
                     hufp(nmx,k),hufp(nmy,k),hufp(nee,k), &
                     rnrm(1),rnrm(2)
      1 format ("Error with bc_generic_inflow boundary condition!",/, &
                "   mach_int = ",es14.6,/, &
                "   aspd_int = ",es14.6,/, &
                "   pres_int = ",es14.6,/, &
                "   dens_int = ",es14.6,/, &
                "   momx_int = ",es14.6,/, &
                "   momy_int = ",es14.6,/, &
                "   enrg_int = ",es14.6,/, &
                "         nx = ",es14.6,/, &
                "         ny = ",es14.6)
      !
    end if
    !
  end do
  !
end subroutine get_bc_generic_inflow_solution
!
!###############################################################################
!
subroutine get_bc_generic_outflow_solution(hufp,gufp,fnrm,pv_in)
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
  ! FNRM  : Face outward normal vector at the flux points
  ! HUFP  : Conservative variables at the host  cell flux points
  ! GUFP  : Conservative variables at the ghost cell flux points
  ! PV_IN : Boundary condition flow variables in primitive form
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: fnrm(:,:)
  real(wp),    intent(in) :: hufp(:,:)
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: gufp(:,:)
  !
  !.. Local Scalars ..
  integer  :: nr,k
  real(wp) :: gm1_inv,ske_int,mach_int,aspd_int
  !
  !.. Local Arrays ..
  real(wp) :: rnrm(1:size(fnrm,dim=1))
  real(wp) :: hvfp(1:size(gufp,dim=1))
  !
continue
  !
  nr = size(fnrm,dim=1)
  !
  ! Precompute the inverse of gamma-1
  !
  gm1_inv = one / (gam - one)
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Unit normal and tangent at this flux point
    !
    rnrm(:) = unit_vector( fnrm(:,k) )
    !
    ! Now convert the current host cell flux point
    ! conservative variables to primitive variables
    !
    hvfp(:) = usp2v_sp( hufp(:,k) )
    !
    ! Compute the speed of sound at this flux point
    !
    aspd_int = sqrt( max( gam*hvfp(nee)/hvfp(nec) , rndoff ) )
    !
    ! Compute the Mach number of the normal velocity at this flux point
    !
    mach_int = dot( hvfp(nmb:nme) , rnrm(1:nr) ) / aspd_int
    !
    if (mach_int < one .and. mach_int >= zero) then
      !
      ! Subsonic outflow
      !
      ! Extrapolate interior density and velocities to boundary
      !
      gufp(:,k) = hufp(:,k)
      !
      ! Compute the interior specific kinetic energy
      !
      ske_int = half * selfdot( hufp(nmb:nme,k) ) / hufp(nec,k)
      !
      ! Compute the total energy using free stream pressure
      ! and interior specific kinetic energy
      !
      gufp(nee,k) = pv_in(nee) * gm1_inv + ske_int
      !
    else if (mach_int >= one) then
      !
      ! Supersonic outflow
      !
      ! Extrapolate all interior variables to the boundary
      !
      gufp(1:nq,k) = hufp(1:nq,k)
      !
    else
      !
      ! The normal velocity at this flux point seems to
      ! indicate this should be an inflow boundary.
      !
      write (iout,*) "Error with bc_generic_outflow boundary condition!"
      !
    end if
    !
  end do
  !
end subroutine get_bc_generic_outflow_solution
!
!###############################################################################
!
subroutine get_bc_generic_freeflow_solution(hufp,gufp,fnrm,cv_in,pv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Generic freeflow boundary condition
  ! ##    - At each flux point, use the normal velocity
  ! ##      Mach number to determine if it is an inflow
  ! ##      or outflow boundary condition and whether
  ! ##      it is subsonic or supersonic.
  ! #####################################################
  ! #####################################################
  !
  ! FNRM  : Face outward normal vector at the flux points
  ! HUFP  : Conservative variables at the host  cell flux points
  ! GUFP  : Conservative variables at the ghost cell flux points
  ! CV_IN : Boundary condition flow variables in conservative form
  ! PV_IN : Boundary condition flow variables in primitive    form
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: fnrm(:,:)
  real(wp),    intent(in) :: hufp(:,:)
  real(wp),    intent(in) :: cv_in(:)
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: gufp(:,:)
  !
  !.. Local Scalars ..
  integer  :: nr,k
  real(wp) :: gm1_inv,ske_inf,ske_int,mach_int,aspd_int
  !
  !.. Local Arrays ..
  real(wp) :: rnrm(1:size(fnrm,dim=1))
  real(wp) :: hvfp(1:size(gufp,dim=1))
  !
continue
  !
  nr = size(fnrm,dim=1)
  !
  ! Precompute the free stream specific kinetic energy
  ! and the inverse of gamma-1
  !
  ske_inf = half * pv_in(nec) * selfdot( pv_in(nmb:nme) )
  !
  gm1_inv = one / (gam - one)
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Unit normal and tangent at this flux point
    !
    rnrm(:) = unit_vector( fnrm(:,k) )
    !
    ! Get the host cell primitive variables at this flux point
    !
    hvfp(:) = usp2v_sp( hufp(:,k) )
    !
    ! Compute the speed of sound at this flux point
    !
    aspd_int = sqrt( max( gam*hvfp(nee)/hvfp(nec) , rndoff ) )
    !
    ! Compute the Mach number of the normal velocity at this flux point
    !
    mach_int = dot( hvfp(nmb:nme) , rnrm(1:nr) ) / aspd_int
    !
    if (mach_int >= one) then
      !
      ! Supersonic outflow
      !
      ! Extrapolate all interior variables to the boundary
      !
      gufp(1:nq,k) = hufp(1:nq,k)
      !
    else if (mach_int <= -one) then
      !
      ! Supersonic inflow
      !
      ! Set all variables to free stream
      !
      gufp(:,k) = cv_in(:)
      !
    else if (mach_int > -one .and. mach_int <= zero) then
      !
      ! Subsonic inflow
      !
      ! Set all variables but pressure to free stream
      !
      gufp(:,k) = cv_in(:)
      !
      ! Compute the total energy using interior pressure
      ! and free stream specific kinetic energy
      !
      gufp(nee,k) = hvfp(nee) * gm1_inv + ske_inf
      !
    else if (mach_int < one .and. mach_int > zero) then
      !
      ! Subsonic outflow
      !
      ! Extrapolate interior density and velocities to boundary
      !
      gufp(:,k) = hufp(:,k)
      !
      ! Compute the interior specific kinetic energy
      !
      ske_int = half * selfdot( hufp(nmb:nme,k) ) / hufp(nec,k)
      !
      ! Compute the total energy using free stream pressure
      ! and interior specific kinetic energy
      !
      gufp(nee,k) = pv_in(nee) * gm1_inv + ske_int
      !
    end if
    !
  end do
  !
end subroutine get_bc_generic_freeflow_solution
!
!###############################################################################
!
subroutine get_bc_sub_inflow_newton(hufp,gufp,fnrm,pv_in,ierr,nf)
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
  ! FNRM  : Face outward normal vector at the flux points
  ! HUFP  : Conservative variables at the host  cell flux points
  ! GUFP  : Conservative variables at the ghost cell flux points
  ! IERR  : count the number of failed convergences from Newton iterations
  ! PV_IN : Boundary condition flow variables in primitive form
  !
  !.. Use Statements ..
  use ovar, only : rgasref_nd
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: fnrm(:,:)
  real(wp),    intent(in) :: hufp(:,:)
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: gufp(:,:)
  integer,  intent(inout) :: ierr
  integer,     intent(in) :: nf
  !
  !.. Local Scalars ..
  integer  :: k,n
  real(wp) :: gm1,gam_over_1mg,tan_aoa,dterr
  real(wp) :: temp_int,aspd_int,veln_int
  real(wp) :: pres_ext,temp_ext,aspd_ext,mach_ext
  real(wp) :: rminus,rafac,veln_ext,func,dfunc,dtemp
  real(wp) :: gam_over_gm1,aspd_ref,mach_ref,total_cond
  real(wp) :: temp_ref,ttot_ref,ptot_ref,aoa,rc,tc
  !
  logical(lk) :: convergence_reached
  !
  !.. Local Arrays ..
  real(wp) :: rnrm(1:size(fnrm,dim=1))
  real(wp) :: dir_cos(1:size(fnrm,dim=1))
  real(wp) :: hvfp(1:size(gufp,dim=1))
  real(wp) :: xyz_fp(1:size(fnrm,dim=1))
  !
  !.. Local Parameters ..
  integer, parameter :: newton_iter = 100
  character(len=*), parameter :: cvel(1:3) = ["x","y","z"]
  !
  real(wp) :: dtemp_save(1:newton_iter)
  !
continue
  !
  ! Precompute the free stream specific kinetic energy
  ! and the inverse of gamma-1
  !
  gm1 = gam - one
  gam_over_gm1 = gam / gm1
  gam_over_1mg = -gam_over_gm1
  !
  ! Compute the inflow conditions
  !
  aspd_ref = sqrt(gam*pv_in(nee)/pv_in(nec))
  mach_ref = norm2(pv_in(nmb:nme))/aspd_ref
  !
  total_cond = one + half*gm1*mach_ref*mach_ref
  !
  temp_ref = temperature_pv_sp( pv_in(:) )
  ttot_ref = temp_ref * total_cond
  ptot_ref = pv_in(nee) * total_cond**gam_over_gm1
  !
  dterr = eps3*ttot_ref
  !
  ! Precompute the direction cosines
  !
  aoa = atan2(pv_in(nmy),pv_in(nmx))
  tan_aoa = tan( aoa * degrees2radians )
  !
  dir_cos(:) = zero
  dir_cos(1) = one / sqrt(one + tan_aoa**2)
  dir_cos(2) = tan_aoa * dir_cos(1)
  !
  do k = 1,size(gufp,dim=2)
    !
    convergence_reached = fals
    !
    ! Unit normal at this flux point
    !
    rnrm(:) = -unit_vector( fnrm(:,k) )
    !
    ! Get the host cell primitive variables at this flux point
    !
    hvfp(:) = usp2v_sp( hufp(:,k) )
    !
    ! Normal velocity of the interior/host cell at this flux point
    !
    veln_int = dot( hvfp(nmb:nme) , rnrm(:) )
    !
    temp_int = temperature_pv_sp( hvfp(:) )
    !
    ! Compute the speed of sound at this flux point
    !
    aspd_int = sqrt( max( gam*rgasref_nd*temp_int , rndoff ) )
    !
    rminus = veln_int - (two*aspd_int/gm1)
    !
    temp_ext = temp_ref
    !
    aspd_ext = aspd_ref
    !
    rafac = one / dot( dir_cos(:) , rnrm(:) )
    !
    veln_ext = rafac * (rminus + two*aspd_ext/gm1)
    !
    func = ttot_ref - temp_ext - half*gm1*veln_ext*veln_ext/(gam*rgasref_nd)
    dfunc = -(one + veln_ext*rafac/aspd_ext)
    !
    newton_iteration: do n = 1,newton_iter
      !
      dtemp = func/dfunc
      dtemp_save(n) = dtemp
      !
      if (abs(dtemp) < dterr) then
        convergence_reached = true
        exit newton_iteration
      end if
      !
      temp_ext = temp_ext - dtemp
      aspd_ext = sqrt( max( gam*rgasref_nd*temp_ext , rndoff ) )
      veln_ext = rafac * (rminus + two*aspd_ext/gm1)
      !
      func = ttot_ref - temp_ext - half*gm1*veln_ext*veln_ext/(gam*rgasref_nd)
      dfunc = -(one + veln_ext*rafac/aspd_ext)
      !
    end do newton_iteration
    !
    mach_ext = veln_ext / aspd_ext
    !
    pres_ext = ptot_ref*(one + half*gm1*mach_ext*mach_ext)**gam_over_1mg
    !
    gufp(nec,k) = pres_ext / (rgasref_nd * temp_ext)
    gufp(nmb:nme,k) = veln_ext*dir_cos(:)
    gufp(nee,k) = pres_ext
    !
    ! Check to make sure the temperature converged
    !
    if (.not. convergence_reached) then
      dtemp = func/dfunc
      call get_xyz_for_face_point(nf,k,xyz_fp)
      write (iout,1) "mypnum   = ",mypnum
      write (iout,1) "nf       = ",nf
      write (iout,1) "k        = ",k
      write (iout,2) "X        = ",xyz_fp(1)
      if (size(xyz_fp) == 3) then
        rc = sqrt(xyz_fp(2)*xyz_fp(2) + xyz_fp(3)*xyz_fp(3))
        tc = atan2(xyz_fp(2),xyz_fp(3)) * radians2degrees
        write (iout,2) "X        = ",xyz_fp(1)
        write (iout,2) "Y        = ",xyz_fp(2)
        write (iout,2) "Z        = ",xyz_fp(3)
        write (iout,2) "R        = ",rc
        write (iout,2) "Theta    = ",tc,         " "
      else if (size(xyz_fp) == 2) then
        write (iout,2) "X        = ",xyz_fp(1)
        write (iout,2) "Y        = ",xyz_fp(2),  " "
      else if (size(xyz_fp) == 1) then
        write (iout,2) "X        = ",xyz_fp(1),  " "
      end if
      write (iout,2) "aoa      = ",aoa
      write (iout,2) "tan_aoa  = ",tan_aoa
      write (iout,2) "dir_cos1 = ",dir_cos(1)
      write (iout,2) "dir_cos2 = ",dir_cos(2), " "
      write (iout,2) "dens_int = ",hvfp(nec)
      do n = nmb,nme
        write (iout,3) cvel(n-nmb+1),"int",hvfp(n)
      end do
      write (iout,2) "pres_int = ",hvfp(nee)
      write (iout,2) "temp_int = ",temp_int
      write (iout,2) "aspd_int = ",aspd_int
      write (iout,2) "veln_int = ",veln_int,   " "
      write (iout,2) "dterr    = ",dterr
      write (iout,2) "func     = ",func
      write (iout,2) "dfunc    = ",dfunc
      write (iout,2) "dtemp    = ",dtemp,      " "
      write (iout,2) "dens_ext = ",gufp(nec,k)
      do n = nmb,nme
        write (iout,3) cvel(n-nmb+1),"ext",gufp(n,k)
      end do
      write (iout,2) "pres_ext = ",pres_ext
      write (iout,2) "temp_ext = ",temp_ext
      write (iout,2) "aspd_ext = ",aspd_ext
      write (iout,2) "veln_ext = ",veln_ext,   " "
      do n = 1,size(dtemp_save)
        if (n == 1) write (iout,4)
        write (iout,5) n,dtemp_save(n)
        if (n == size(dtemp_save)) write (iout,5) n+1,dtemp, " "
      end do
      1 format (1x,a,i0,:,/,a)
      2 format (1x,a,es14.6,:,/,a)
      3 format (1x,"vel",a1,"_",a3," = ",es14.6)
      4 format (/,50("-"),/, &
                  4x,"Values of dtemp for each Newton Iteration",/, &
                  50("-"))
      5 format (4x,i3,4x,es23.15,:,/,50("-"),a,/)
      ierr = ierr + 1
    end if
    !
    gufp(1:nq,k) = vsp2u_sp( gufp(1:nq,k) )
    !
  end do
  !
end subroutine get_bc_sub_inflow_newton
!
!###############################################################################
!
pure &
subroutine get_xyz_for_face_point(nf,kf,xyz_fp)
  !
  !.. Use Statements ..
  use order_mod,         only : maxSP
  use geovar,            only : cell,face,fp,xyz
  use interpolation_mod, only : interp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nf
  integer,     intent(in) :: kf
  real(wp), intent(inout) :: xyz_fp(:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,lf,nc,np,n0
  integer  :: this_geom,this_order
  integer  :: face_geom,face_order
  integer  :: fp_offset,rotation
  real(wp) :: a
  !
  !.. Local Arrays ..
  real(wp) :: txyz(1:maxSP,1:size(xyz_fp))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_xyz_for_face_point"
  !
continue
  !
  nc = face(nf)%left%cell
  !
  this_geom = cell(nc)%geom
  this_order = cell(nc)%order
  !
  n0 = cell(nc)%beg_sp - 1
  np = cell(nc)%end_sp - n0
  !
  do k = 1,np
    do l = 1,size(xyz_fp)
      txyz(k,l) = xyz(l,n0+k)
    end do
  end do
  !
  face_geom = face(nf)%geom
  face_order = face(nf)%order
  !
  fp_offset = face(nf)%left%fp_offset
  rotation  = face(nf)%left%rotation
  !
  lf = fp_offset + fp%geord(face_geom,face_order)% &
                                    rot(rotation)% &
                   pts(kf)
  !
  xyz_fp(:) = zero
  !
  do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
         interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
    !
    j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
    a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
    !
    do l = 1,size(xyz_fp)
      xyz_fp(l) = xyz_fp(l) + a * txyz(j,l)
    end do
    !
  end do
  !

end subroutine get_xyz_for_face_point
!
!###############################################################################
!
subroutine get_bc_sub_inflow_solution(hufp,gufp,cv_in,pv_in)
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
  ! HUFP  : Conservative variables at the host  cell flux points
  ! GUFP  : Conservative variables at the ghost cell flux points
  ! CV_IN : Boundary condition flow variables in conservative form
  ! PV_IN : Boundary condition flow variables in primitive    form
  !
  !.. Use Statements ..
 !use ovar, only : ptotref_nd,ttotref_nd
  use ovar, only : rgasref_nd
  use ovar, only : sub_inflow_method
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: hufp(:,:)
  real(wp),    intent(in) :: cv_in(:)
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: gufp(:,:)
  !
  !.. Local Scalars ..
  integer  :: k
  real(wp) :: pres_int,ske_inf,gm1_inv
  real(wp) :: gm1,rt,rhot,a,b,c,velx2
  real(wp) :: rm,rmsq,aa,asq,gam_over_gm1
  real(wp) :: p0_ext,t_ext,t0_ext
  real(wp) :: aspd,mach,total_cond
  !
continue
  !
  ! Precompute the free stream specific kinetic energy
  ! and the inverse of gamma-1
  !
  gm1 = gam - one
  gm1_inv = one / gm1
  gam_over_gm1 = gam*gm1_inv
  !
  t_ext = temperature_pv_sp(pv_in)
  aspd = sqrt(gam*pv_in(nee)/pv_in(nec))
  mach = norm2(pv_in(nmb:nme))/aspd
  !
  total_cond = one + half*gm1*mach*mach
  !
  t0_ext = t_ext*total_cond
  p0_ext = pv_in(nee) * total_cond**gam_over_gm1
  !
  if (sub_inflow_method == 2) then
    !
    ! Precompute the free stream specific kinetic energy
    ! and the inverse of gamma-1
    !
    ske_inf = half * pv_in(nec) * selfdot( pv_in(nmb:nme) )
    !
    do k = 1,size(gufp,dim=2)
      !
      ! Set all variables but pressure to free stream
      !
      gufp(:,k) = cv_in(:)
      !
      ! Get the interior pressure
      !
      pres_int = pressure_cv_sp( hufp(:,k) )
      !
      ! Compute the total energy using interior pressure
      ! and free stream specific kinetic energy
      !
      gufp(nee,k) = pres_int * gm1_inv + ske_inf
      !
    end do
    !
  else if (sub_inflow_method == 3) then
    !
   !rt = rgasref_nd*ttotref_nd
    rt = rgasref_nd*t0_ext
    !
   !rhot = ptotref_nd/rt
    rhot = p0_ext/rt
    a = rt*gm1_inv
    b = two*gam*a
    c = half*gm1/gam
    !
    do k = 1,size(gufp,dim=2)
      !
      velx2 = hufp(nmx,k)**2
      !
      gufp(nec,k) = rhot * ( (one + velx2/(b-velx2))**(-gm1_inv) )
      gufp(nmb:nme,k) = zero
      gufp(nmx,k) = gufp(nec,k)*hufp(nmx,k)
      gufp(nee,k) = gufp(nec,k) * (a + c*velx2)
      !
    end do
    !
  else if (sub_inflow_method == 4) then
    !
    do k = 1,size(gufp,dim=2)
      !
      ! Get the interior pressure
      !
      pres_int = pressure_cv_sp( hufp(:,k) )
      !
     !rmsq = ( (ptotref_nd/pres_int)**(gm1/gam) - one ) * two/gm1
      rmsq = ( (p0_ext/pres_int)**(gm1/gam) - one ) * two/gm1
      rm   = sqrt(rmsq)
      !
      asq = one / (one+half*gm1*rmsq)
      aa = sqrt(asq)
      !
      ! Compute the total energy using interior pressure
      ! and free stream specific kinetic energy
      !
      gufp(nec,k) = gam*pres_int/asq
      gufp(nmb:nme,k) = zero
      gufp(nmx,k) = aa*rm
      gufp(nee,k) = pres_int
      !
      gufp(1:nq,k) = vsp2u_sp( gufp(1:nq,k) )
      !
    end do
    !
  end if
  !
end subroutine get_bc_sub_inflow_solution
!
!###############################################################################
!
subroutine get_bc_sub_outflow_solution(hufp,gufp,fnrm,pv_in)
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
  ! FNRM  : Face outward normal vector at the flux points
  ! HUFP  : Conservative variables at the host  cell flux points
  ! GUFP  : Conservative variables at the ghost cell flux points
  ! PV_IN : Boundary condition flow variables in primitive form
  !
  !.. Use Statements ..
  use ovar, only : ptotref_nd
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: fnrm(:,:)
  real(wp),    intent(in) :: hufp(:,:)
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: gufp(:,:)
  !
  !.. Local Scalars ..
  integer  :: k
  real(wp) :: gm1_inv,ske_int,vn_int,aspd_int
  !
  !.. Local Arrays ..
  real(wp) :: rnrm(1:size(fnrm,dim=1))
  real(wp) :: hvfp(1:size(gufp,dim=1))
  !
  !.. Local Parameters ..
  logical(lk), parameter :: use_old_sub_outflow = fals
 !logical(lk), parameter :: use_old_sub_outflow = true
  !
continue
  !
  ! Precompute the inverse of gamma-1
  !
  gm1_inv = one / (gam - one)
  !
  if (use_old_sub_outflow) then
    !
    do k = 1,size(gufp,dim=2)
      !
      ! Extrapolate interior density and velocities to boundary
      !
      gufp(:,k) = hufp(:,k)
      !
      ! Compute the interior specific kinetic energy
      !
      ske_int = half * selfdot( hufp(nmb:nme,k) ) / hufp(nec,k)
      !
      ! Compute the total energy using free stream pressure
      ! and interior specific kinetic energy
      !
      gufp(nee,k) = pv_in(nee) * gm1_inv + ske_int
      !
    end do
    !
  else
    !
    do k = 1,size(gufp,dim=2)
      !
      ! Unit normal and tangent at this flux point
      !
      rnrm(:) = unit_vector( fnrm(:,k) )
      !
      ! Now convert the current host cell flux point
      ! conservative variables to primitive variables
      !
      hvfp(:) = usp2v_sp( hufp(:,k) )
      !
      vn_int = dot( rnrm(:) , hvfp(nmb:nme) )
      aspd_int = sqrt( max( gam*hvfp(nee)/hvfp(nec) , rndoff ) )
      !
      ! Extrapolate interior primitive variables to boundary
      !
      gufp(1:nq,k) = hvfp(1:nq)
      !
      ! If the normal velocity from the interior says this is an inflow
      ! boundary, set the ghost velocity to the magnitude of the interior
      ! velocity multiplied by the unit normal
      !
      if (vn_int < zero) then
        gufp(nmb:nme,k) = norm2( hvfp(nmb:nme) ) * rnrm(:)
      end if
      !
      ! Set the boundary pressure to the reference pressure if the normal
      ! velocity from the interior is subsonic, otherwise set the boundary
      ! pressure to the reference total pressure
      !
      if (abs(vn_int)/aspd_int >= one) then
        gufp(nee,k) = ptotref_nd
      else
        gufp(nee,k) = pv_in(nee)
      end if
      !
      ! Convert the primitive variables to conservative variables
      !
      gufp(1:nq,k) = vsp2u_sp( gufp(1:nq,k) )
      !
    end do
    !
  end if
  !
end subroutine get_bc_sub_outflow_solution
!
!###############################################################################
!
subroutine get_bc_sup_inflow_solution(gufp,cv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Supersonic inflow boundary condition
  ! ##    - Fix all variables to free stream
  ! #####################################################
  ! #####################################################
  !
  ! GUFP  : Conservative variables at the ghost cell flux points
  ! CV_IN : Boundary condition flow variables in conservative form
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: cv_in(:)
  real(wp), intent(inout) :: gufp(:,:)
  !
  !.. Local Scalars ..
  integer :: k
  !
continue
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Set all variables to free stream
    !
    gufp(:,k) = cv_in(:)
    !
  end do
  !
end subroutine get_bc_sup_inflow_solution
!
!###############################################################################
!
subroutine get_bc_sup_outflow_solution(hufp,gufp)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Supersonic outflow boundary condition
  ! ##    - Extrapolate all variables
  ! #####################################################
  ! #####################################################
  !
  ! HUFP : Conservative variables at the host  cell flux points
  ! GUFP : Conservative variables at the ghost cell flux points
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: hufp(:,:)
  real(wp), intent(inout) :: gufp(:,:)
  !
  !.. Local Scalars ..
  integer :: k
  !
continue
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Extrapolate all interior variables to the boundary
    !
    gufp(1:nq,k) = hufp(1:nq,k)
    !
  end do
  !
end subroutine get_bc_sup_outflow_solution
!
!###############################################################################
!
subroutine get_bc_slip_wall_solution(hufp,gufp,fnrm)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Slip wall / symmetry boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FNRM : Face outward normal vector at the flux points
  ! HUFP : Conservative variables at the host  cell flux points
  ! GUFP : Conservative variables at the ghost cell flux points
  !
  !.. Use Statements ..
  use ovar, only : walls_are_exact
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: fnrm(:,:)
  real(wp),    intent(in) :: hufp(:,:)
  real(wp), intent(inout) :: gufp(:,:)
  !
  !.. Local Scalars ..
  integer  :: k
  real(wp) :: vn,wall_solution_location
  !
  !.. Local Arrays ..
  real(wp) :: rnrm(1:size(fnrm,dim=1))
  real(wp) :: hvfp(1:size(gufp,dim=1))
  !
continue
  !
  wall_solution_location = merge(zero,one,walls_are_exact)
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Unit normal at this flux point
    !
    rnrm(:) = unit_vector( fnrm(:,k) )
    !
    ! Get the interior primitive variables
    !
    hvfp(:) = usp2v_sp( hufp(:,k) )
    !
    ! Compute the normal velocity at this boundary flux point
    !
    vn = dot( hvfp(nmb:nme) , rnrm )
    !
    ! For velocities, subtract (wall_solution_location+one) * (the wall
    ! normal velocity) from the interior velocity to get the wall
    ! boundary velocity.
    !
    hvfp(nmb:nme) = hvfp(nmb:nme) - (wall_solution_location+one)*vn*rnrm(:)
    !
    ! Convert the primitive variables to conservative.
    ! This ensures that the interior density and pressure are
    ! extrapolated to the boundary and the total energy at the
    ! boundary is consistent with these two states and the
    ! computed boundary velocity.
    !
    gufp(:,k) = vsp2u_sp( hvfp(:) )
    !
  end do
  !
end subroutine get_bc_slip_wall_solution
!
!###############################################################################
!
subroutine get_bc_adiabatic_wall_solution(hufp,gufp,pv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip adiabatic wall
  ! #####################################################
  ! #####################################################
  !
  ! HUFP : Conservative variables at the host  cell flux points
  ! GUFP : Conservative variables at the ghost cell flux points
  !
  !.. Use Statements ..
  use ovar, only : walls_are_exact
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: hufp(:,:)
  real(wp), intent(inout) :: gufp(:,:)
  real(wp),    intent(in) :: pv_in(:)
  !
  !.. Local Scalars ..
  integer  :: k
  real(wp) :: velocity_const,moving_wall_const
  !
  !.. Local Arrays ..
  real(wp) :: hvfp(1:size(gufp,dim=1))
  !
continue
  !
  velocity_const = merge(zero,-one,walls_are_exact)
  moving_wall_const = merge(one,two,walls_are_exact)
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Get the interior primitive variables
    !
    hvfp(:) = usp2v_sp( hufp(:,k) )
    !
    ! Apply the no-slip condition to the boundary velocities
    !
    hvfp(nmb:nme) =  hvfp(nmb:nme) * velocity_const + &
                    pv_in(nmb:nme) * moving_wall_const
    !
    ! Convert the primitive variables to conservative.
    ! This ensures that the interior density and pressure are
    ! extrapolated to the boundary and the total energy at the
    ! boundary is consistent with these two states and the
    ! computed boundary velocity.
    !
    gufp(:,k) = vsp2u_sp( hvfp(:) )
    !
  end do
  !
end subroutine get_bc_adiabatic_wall_solution
!
!###############################################################################
!
subroutine get_bc_isothermal_wall_solution(hufp,gufp,wall_temperature,pv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip isothermal wall
  ! #####################################################
  ! #####################################################
  !
  ! HUFP : Conservative variables at the host  cell flux points
  ! GUFP : Conservative variables at the ghost cell flux points
  ! WALL_TEMPERATURE :: isothermal wall temperature
  !
  !.. Use Statements ..
  use ovar, only : rgasref_nd
  use ovar, only : walls_are_exact
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: hufp(:,:)
  real(wp), intent(inout) :: gufp(:,:)
  real(wp),    intent(in) :: wall_temperature
  real(wp),    intent(in) :: pv_in(:)
  !
  !.. Local Scalars ..
  integer  :: k
  real(wp) :: velocity_const,moving_wall_const
  !
  !.. Local Arrays ..
  real(wp) :: hvfp(1:size(gufp,dim=1))
  !
continue
  !
  velocity_const = merge(zero,-one,walls_are_exact)
  moving_wall_const = merge(one,two,walls_are_exact)
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Get the interior primitive variables
    !
    hvfp(:) = usp2v_sp( hufp(:,k) )
    !
    ! Apply the no-slip condition to the boundary velocities
    !
    hvfp(nmb:nme) =  hvfp(nmb:nme) * velocity_const + &
                    pv_in(nmb:nme) * moving_wall_const
    !
    ! Keeping the wall pressure constant, compute the wall
    ! density using the isothermal wall temperature.
    !
    hvfp(nec) = hvfp(nee) / rgasref_nd / wall_temperature
    !
    ! Convert the primitive variables to conservative.
    ! This ensures that the total energy at the boundary is
    ! consistent with the extrapolated internal pressure, the
    ! isothermal wall temperature, and the computed boundary
    ! velocity.
    !
    gufp(:,k) = vsp2u_sp( hvfp(:) )
    !
  end do
  !
end subroutine get_bc_isothermal_wall_solution
!
!###############################################################################
!
subroutine get_bc_iso_wall_solution_old(hufp,gufp,wall_temperature,fnrm)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip isothermal wall
  ! #####################################################
  ! #####################################################
  !
  ! HUFP : Conservative variables at the host  cell flux points
  ! GUFP : Conservative variables at the ghost cell flux points
  ! WALL_TEMPERATURE :: isothermal wall temperature
  !
  !.. Use Statements ..
  use ovar, only : tref,rgasref_nd
  use ovar, only : walls_are_exact
  use ovar, only : aref
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: hufp(:,:)
  real(wp),    intent(in) :: wall_temperature
  real(wp),    intent(in) :: fnrm(:,:)
  real(wp), intent(inout) :: gufp(:,:)
  !
  !.. Local Scalars ..
  integer  :: k
  real(wp) :: wall_solution_location
  real(wp) :: uwall
  !
  !.. Local Arrays ..
  real(wp) :: hvfp(1:size(gufp,dim=1))
  real(wp) :: unrm(1:size(fnrm,dim=1))
  !
  !.. Local Parameters ..
 !logical(lk), parameter :: use_moving_wall = fals
  logical(lk), parameter :: use_moving_wall = true
  !
continue
  !
  wall_solution_location = merge(zero,one,walls_are_exact)
  uwall = 30.0_wp / aref
  !
  do k = 1,size(gufp,dim=2)
    !
    unrm(:) = unit_vector(fnrm(:,k))
    !
    ! Get the interior primitive variables
    !
    hvfp(:) = usp2v_sp( hufp(:,k) )
    !
    ! Apply the no-slip condition to the boundary velocities
    !
    hvfp(nmb:nme) = -hvfp(nmb:nme) * wall_solution_location
    if (use_moving_wall) then
      if ( unrm(2) > half ) hvfp(nmb) = uwall
    end if
    !
    ! Keeping the wall pressure constant, compute the wall
    ! density using the isothermal wall temperature.
    !
    hvfp(nec) = hvfp(nee) / rgasref_nd / wall_temperature
    !
    ! Convert the primitive variables to conservative.
    ! This ensures that the total energy at the boundary is
    ! consistent with the extrapolated internal pressure, the
    ! isothermal wall temperature, and the computed boundary
    ! velocity.
    !
    gufp(:,k) = vsp2u_sp( hvfp(:) )
    !
  end do
  !
end subroutine get_bc_iso_wall_solution_old
!
!###############################################################################
!
subroutine get_bc_periodic_solution(pufp,gufp)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Periodic boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! PUFP : Conservative variables at the periodic cell flux points
  ! GUFP : Conservative variables at the ghost    cell flux points
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: pufp(:,:)
  real(wp), intent(inout) :: gufp(:,:)
  !
  !.. Local Scalars ..
  integer :: k,nfp
  !
continue
  !
  nfp = size(gufp,dim=2)
  !
  ! The conservative variables at this boundary flux point correspond
  ! to the conservative variables at the interior boundary flux points
  ! within the matching periodic cell. Therefore just copy these
  ! values over.
  !
  if (nme-nmb+1 == 2) then
    !
    ! For some reason we need to reverse the copy order for 2D
    !    i.e., gufp(:,k) = pufp(:,size(gufp,dim=2)+1-k)
    ! NOTE: We also need to reverse the y velocity direction
    !
    do k = 1,nfp
      !
      gufp(:,k) = pufp(:,nfp-k+1)
     !gufp(nme,k) = -pufp(nme,k)
      !
    end do
    !
  else
    !
    do k = 1,nfp
      !
      gufp(:,k) = pufp(:,k)
      !
    end do
    !
  end if
  !
end subroutine get_bc_periodic_solution
!
!###############################################################################
!
subroutine get_bc_mms_dirichlet_solution(hufp,gufp,fxyz,fnrm)
  !
  ! #########################################################
  ! #########################################################
  ! ##  MMS Dirichlet (generic freeflow) Boundary Condition
  ! ##    - At each flux point, use the normal velocity
  ! ##      Mach number to determine if it is an inflow
  ! ##      or outflow boundary condition and whether
  ! ##      it is subsonic or supersonic.
  ! #########################################################
  ! #########################################################
  !
  ! GUFP : Conservative variables at the ghost cell flux points
  ! HUFP : Conservative variables at the host  cell flux points
  ! FXYZ : Coordinates of the flux points
  ! FNRM : Face outward normal vector at the flux points
  !
  !.. Use Statements ..
  use mms_mod, only : mms_solution_at_xyz
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: hufp(:,:)
  real(wp),    intent(in) :: fxyz(:,:)
  real(wp),    intent(in) :: fnrm(:,:)
  real(wp), intent(inout) :: gufp(:,:)
  !
  !.. Local Scalars ..
  integer  :: nr,k
  real(wp) :: gm1_inv,ske_inf,ske_int,mach_int,aspd_int
  !
  !.. Local Arrays ..
  real(wp) :: rnrm(1:size(fnrm,dim=1))
  real(wp) :: hvfp(1:size(gufp,dim=1))
  real(wp) :: mms_sol(1:size(gufp,dim=1))
  !
continue
  !
  nr = size(fnrm,dim=1)
  !
  gm1_inv = one / (gam - one)
  !
  do k = 1,size(gufp,dim=2)
    !
    ! Get the 'freestream' conditions from the MMS solution
    !
    mms_sol(1:nq) = mms_solution_at_xyz( fxyz(:,k) , &
                                         return_nondim=true , &
                                         return_cv=true )
    !
    ! Unit normal and tangent at this flux point
    !
    rnrm(:) = unit_vector( fnrm(:,k) )
    !
    ! Get the host cell primitive variables at this flux point
    !
    hvfp(:) = usp2v_sp( hufp(:,k) )
    !
    ! Compute the speed of sound at this flux point
    !
    aspd_int = sqrt( max( gam*hvfp(nee)/hvfp(nec) , rndoff ) )
    !
    ! Compute the Mach number of the normal velocity at this flux point
    !
    mach_int = dot( hvfp(nmb:nme) , rnrm(1:nr) ) / aspd_int
    !
    if (mach_int >= one) then
      !
      ! Supersonic outflow
      !
      ! Extrapolate all interior variables to the boundary
      !
      gufp(1:nq,k) = hufp(1:nq,k)
      !
    else if (mach_int <= -one) then
      !
      ! Supersonic inflow
      !
      ! Set all variables to free stream
      !
      gufp(:,k) = mms_sol(1:nq)
      !
    else if (mach_int > -one .and. mach_int <= zero) then
      !
      ! Subsonic inflow
      !
      ! Set all variables but pressure to free stream
      !
      gufp(:,k) = mms_sol(:)
      !
      ! Compute the total energy using interior pressure
      ! and free stream specific kinetic energy
      !
      ske_inf = half * selfdot( mms_sol(nmb:nme) ) / mms_sol(nec)
      !
      gufp(nee,k) = hvfp(nee) * gm1_inv + ske_inf
      !
    else if (mach_int < one .and. mach_int > zero) then
      !
      ! Subsonic outflow
      !
      ! Extrapolate interior density and velocities to boundary
      !
      gufp(:,k) = hufp(:,k)
      !
      ! Compute the interior specific kinetic energy
      !
      ske_int = half * selfdot( hufp(nmb:nme,k) ) / hufp(nec,k)
      !
      ! Compute the total energy using free stream pressure
      ! and interior specific kinetic energy
      !
      gufp(nee,k) = pressure_cv_sp( mms_sol(1:nq) ) * gm1_inv + ske_int
      !
    end if
    !
  end do
  !
end subroutine get_bc_mms_dirichlet_solution
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
! Functions that compute the gradients of the
! conservative variables at the boundary flux points.
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine get_bc_periodic_gradient(hdufp,gdufp)
  !
  ! #####################################################
  ! #####################################################
  ! ##  General flow-through boundary conditions
  ! #####################################################
  ! #####################################################
  !
  ! HDUFP : Gradients of the conservative variables at the host  cell flux pts
  ! GDUFP : Gradients of the conservative variables at the ghost cell flux pts
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: hdufp(:,:,:)
  real(wp), intent(inout) :: gdufp(:,:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,k,nfp
  !
continue
  !
  ! Simply copy the interior gradients to the boundary flux points
  !
  nfp = size(gdufp,dim=3)
  !
  if (size(gdufp,dim=1) == 2) then
    !
    ! For some reason we need to reverse the copy order for 2D
    !
    do k = 1,nfp
      do j = 1,size(gdufp,dim=2)
        do i = 1,size(gdufp,dim=1)
          gdufp(i,j,k) = hdufp(i,j,nfp-k+1)
        end do
      end do
    end do
    !
  else
    !
    do k = 1,nfp
      do j = 1,size(gdufp,dim=2)
        do i = 1,size(gdufp,dim=1)
          gdufp(i,j,k) = hdufp(i,j,k)
        end do
      end do
    end do
    !
  end if
  !
end subroutine get_bc_periodic_gradient
!
!###############################################################################
!
subroutine get_bc_flowthrough_gradient(hdufp,gdufp)
  !
  ! #####################################################
  ! #####################################################
  ! ##  General flow-through boundary conditions
  ! #####################################################
  ! #####################################################
  !
  ! HDUFP : Gradients of the conservative variables at the host  cell flux pts
  ! GDUFP : Gradients of the conservative variables at the ghost cell flux pts
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: hdufp(:,:,:)
  real(wp), intent(inout) :: gdufp(:,:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,k
  !
continue
  !
  ! Simply copy the interior gradients to the boundary flux points
  !
  do k = 1,size(gdufp,dim=3)
    do j = 1,size(gdufp,dim=2)
      do i = 1,size(gdufp,dim=1)
        gdufp(i,j,k) = hdufp(i,j,k)
      end do
    end do
  end do
  !
end subroutine get_bc_flowthrough_gradient
!
!###############################################################################
!
subroutine get_bc_adiabatic_wall_gradient(hufp,hdufp,gdufp)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip adiabatic wall
  ! #####################################################
  ! #####################################################
  !
  ! HUFP  : Conservative variables at the host cell flux points
  ! HDUFP : Gradients of the conservative variables at the host  cell flux pts
  ! GDUFP : Gradients of the conservative variables at the ghost cell flux pts
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: hufp(:,:)
  real(wp),    intent(in) :: hdufp(:,:,:)
  real(wp), intent(inout) :: gdufp(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: nr,k,m,l
  real(wp) :: rho_inv,vel
 !real(wp) :: rho_inv2,vel_over_rho
  !
continue
  !
  nr = size(gdufp,dim=1)
  !
  do k = 1,size(gdufp,dim=3)
    !
    ! Copy over the density gradients
    !
    gdufp(1:nr,nec,k) = hdufp(1:nr,nec,k)
    !
    ! Convert the momentum gradients to velocity gradients
    !
    rho_inv = one/hufp(nec,k)
    !
    do m = nmb,nme
      vel = rho_inv * hufp(m,k)
      do l = 1,nr
        gdufp(l,m,k) = hdufp(l,m,k) - vel*hdufp(l,nec,k)
        gdufp(l,m,k) = rho_inv * gdufp(l,m,k)
      end do
    end do
    !
    ! Another algorithm to compute the velocity gradients. Unsure if
    ! performing the above subtraction step first before multiplying
    ! through by rho_inv could cause any precision errors.
    !
   !rho_inv = one/hufp(nec,k)
   !rho_inv2 = rho_inv*rho_inv
   !!
   !do m = nmb,nme
   !  vel_over_rho = rho_inv2 * hufp(m,k)
   !  do l = 1,nr
   !    gdufp(l,m,k) = rho_inv*hdufp(l,m,k) - vel_over_rho*hdufp(l,nec,k)
   !  end do
   !end do
    !
    ! Finally, since this is an adiabatic wall boundary,
    ! set the temperature gradients to zero.
    !
    gdufp(1:nr,nq,k) = zero
    !
  end do
  !
end subroutine get_bc_adiabatic_wall_gradient
!
!###############################################################################
!
subroutine get_bc_wall_gradient(hufp,hdufp,gdufp)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Non-adiabatic wall
  ! #####################################################
  ! #####################################################
  !
  ! HUFP  : Conservative variables at the host cell flux points
  ! HDUFP : Gradients of the conservative variables at the host  cell flux pts
  ! GDUFP : Gradients of the conservative variables at the ghost cell flux pts
  !
  !.. Use Statements ..
  use ovar, only : rgasref_nd
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: hufp(:,:)
  real(wp),    intent(in) :: hdufp(:,:,:)
  real(wp), intent(inout) :: gdufp(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: np,nr,k,m,l
  real(wp) :: rho_inv,vel,pvar,cv_inv
  real(wp) :: a,b
 !real(wp) :: rho_inv2,pvar_over_rho
  !
  !.. Local Arrays ..
  real(wp) :: dtdx(1:size(gdufp,dim=1))
  !
continue
  !
  nr = size(gdufp,dim=1)
  np = size(gdufp,dim=3)
  !
  ! Precompute the inverse of specific heat at constant volume
  !
  cv_inv = (gam - one) / rgasref_nd
  !
  ! Copy over the interior density gradient to the boundary flux points
  !
  do k = 1,np
    do l = 1,nr
      gdufp(l,nec,k) = hdufp(l,nec,k)
    end do
  end do
  !
  ! We will also copy the other conservative variable gradients from
  ! the interior cell to the boundary flux points. However, since
  ! we treat the viscous fluxes for wall boundaries differently from
  ! the flow-through boundaries, we need to ensure that the wall
  ! gradients are in the right form for the wall viscous flux function.
  ! This means pre-converting the conservative variable gradients to
  ! primitive variable gradients.
  !
  do k = 1,np
    !
    ! Compute the inverse of density
    !
    rho_inv = one/hufp(nec,k)
    !
    ! Using the chain rule:
    !
    ! 1. Compute the velocity gradients (m = 2,nq-1)
    ! 2. Compute the total energy gradient and store in the
    !    location for the temperature gradient (m = nq)
    !
    do m = nmb,nq
      pvar = rho_inv * hufp(m,k)
      do l = 1,nr
        gdufp(l,m,k) = hdufp(l,m,k) - pvar*hdufp(l,nec,k)
        gdufp(l,m,k) = rho_inv * gdufp(l,m,k)
      end do
    end do
   !!
   !! Another algorithm to compute the velocity gradients. Unsure if
   !! performing the above subtraction step first before multiplying
   !! through by rho_inv could cause any precision errors.
   !!
   !rho_inv = one/hufp(nec,k)
   !rho_inv2 = rho_inv*rho_inv
   !!
   !do m = nmb,nq
   !  pvar_over_rho = rho_inv2 * hufp(m,k)
   !  do l = 1,nr
   !    gdufp(l,m,k) = rho_inv*hdufp(l,m,k) - pvar_over_rho*hdufp(l,nec,k)
   !  end do
   !end do
    !
    ! The location for the temperature gradient currently contains
    ! the total energy gradient from the previous loop.
    ! Now use the velocity gradients to subtract the gradient of
    ! kinetic-energy-per-unit-mass from the total energy gradient.
    !
   !dtdx(1:nr) = gdufp(1:nr,nq,k)
   !do m = nmb,nme
   !  vel = rho_inv * hufp(m,k)
   !  do l = 1,nr
   !    dtdx(l) = dtdx(l) - vel * gdufp(l,m,k)
   !  end do
   !end do
   !!
   !! Finish computing the temperature gradient by dividing through by
   !! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
   !!
   !gdufp(1:nr,nq,k) = dtdx(1:nr) * cv_inv
    !
    do l = 1,nr
      a = rho_inv * selfdot( hufp(nmb:nme,k) ) * hdufp(l,nec,k)
      b = dot( hufp(nmb:nme,k) , hdufp(l,nmb:nme,k) )
      gdufp(l,nq,k) = gdufp(l,nq,k) + rho_inv*rho_inv*(a-b)
      gdufp(l,nq,k) = gdufp(l,nq,k) * cv_inv
    end do
    !
    !
  end do
  !
end subroutine get_bc_wall_gradient
!
!###############################################################################
!
subroutine get_bc_slip_wall_gradient(hdufp,gdufp,fnrm)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Slip wall / symmetry boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FNRM  : Face outward normal vector at the flux points
  ! HDUFP : Gradients of the conservative variables at the host  cell flux pts
  ! GDUFP : Gradients of the conservative variables at the ghost cell flux pts
  !
  !.. Use Statements ..
  use ovar, only : walls_are_exact
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: fnrm(:,:)
  real(wp),    intent(in) :: hdufp(:,:,:)
  real(wp), intent(inout) :: gdufp(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: nr,k,m,l
  real(wp) :: wall_solution_location,reflect_const
  !
  !.. Local Arrays ..
  real(wp) :: rnrm(1:size(fnrm,dim=1))
  !
continue
  !
  wall_solution_location = merge(zero,one,walls_are_exact)
  reflect_const = wall_solution_location + one
  !
  nr = size(gdufp,dim=1)
  !
  do k = 1,size(gdufp,dim=3)
    !
    ! Unit normal at this flux point
    !
    rnrm(:) = unit_vector( fnrm(:,k) )
    !
    ! Set the boundary gradient vector to the reflection of the
    ! interior gradient vector about the boundary face unit normal
    !
    do m = 1,nq
      do l = 1,nr
        gdufp(l,m,k) = hdufp(l,m,k) - &
                       reflect_const*rnrm(l)*dot( hdufp(1:nr,m,k) , rnrm(1:nr) )
      end do
    end do
    !
  end do
  !
end subroutine get_bc_slip_wall_gradient
!
!###############################################################################
!
subroutine output_bface
  !
  !.. Use Statements ..
  use geovar, only : bface
  !
  !.. Local Scalars ..
  integer :: i,n,maxint,myio,ierr
  character(len=100) :: write_format
  character(len=100) :: file_format
  character(len=100) :: my_file
  !
  !.. Local Arrays ..
  integer :: intlen(1:size(bface,dim=1)+1)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "output_bface"
  !
continue
  !
  myio = 1456 + mypnum
  !
  ! Get the maximum length of the string containing the name of the BC types
  !
  intlen(1) = maxval( len_trim( bc_integer_to_string(bface(1,:)) ) )
  !
  do i = 1,size(intlen)-1
    n = 1
    maxint = maxval(bface(i,:))
    check_bface_size: do
      if (maxint/(10**n) == 0) then
        intlen(i+1) = n+1
        exit check_bface_size
      else
        n = n + 1
        cycle check_bface_size
      end if
    end do check_bface_size
  end do
  !
  ! Primary write of the format statement which has an unlimited
  ! '*' repetition multiplier on the last quantity to allow for
  ! an arbitrarily large number of values contained in intlen.
  ! The ":" within the part that is unlimitedly repeated serves
  ! to terminate the format statement at that point preventing
  ! it from going back and adding on an extra ',3x,i' at the end.
  !
  write (write_format,1) (intlen(i),i=1,size(intlen))
 !write (iout,'(a,a,a)') "write_format = '",trim(write_format),"'"
  !
  ! Finally, add on the parentheses to the format statement
  !
  write (write_format,2) trim(write_format)
 !write (iout,'(a,a,a)') "write_format = '",trim(write_format),"'"
  !
  ! Get the number of digits of the integer for the number of cpus
  !
  n = 1
  check_ncpu_size: do
    if (ncpu/(10**n) == 0) then
      write (file_format,3) n
      exit check_ncpu_size
    else
      n = n + 1
      cycle check_ncpu_size
    end if
  end do check_ncpu_size
  !
  ! Create the name of the output file
  !
  write (my_file,fmt=trim(file_format)) mypnum
  !
  ! Open the output file
  !
  open (myio,file=trim(my_file),status="replace",action="write", &
             form="formatted",position="rewind", &
             iostat=ierr,iomsg=error_message)
  call io_error(pname,trim(my_file),1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Output the bface array
  !
  do n = 1,size(bface,dim=2)
    write (myio,fmt=trim(write_format)) &
                               trim(bc_integer_to_string(bface(1,n))), &
                                   (bface(i,n),i=1,size(bface,dim=1))
  end do
  !
  ! Close the output file
  !
  close (myio,iostat=ierr,iomsg=error_message)
  call io_error(pname,trim(my_file),2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Format Statements
  !
  1 format ('a',i0,',i',i0,',": ",1x,i',i0,',1x,i',i0,',1x,i',i0, &
            ',5x,1x,i',i0,',1x,i',i0,',1x,i',i0,',1x,i',i0,',1x,i',i0, &
            100(',3x,i',i0,:))
           !*(',3x,i',i0,:))
  2 format ('(',a,')')
  3 format ('("bface.",i0.',i0,',".dat")')
  !
end subroutine output_bface
!
!###############################################################################
!
subroutine exchange_bc_conflicts_usp(faceusp,override)
  !
  !.. Use Statements ..
  use geovar,  only : bc_conflict
  use ovar,    only : loc_flux_pts
  use ovar,    only : use_bc_conflicts
  use flowvar, only : face_var
  !
  !.. Formal Arguments ..
  type(face_var), intent(inout) :: faceusp(:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: override
  !
  !.. Local Scalars ..
  integer :: i,n,nf,km,ks,nfm,nfs
  integer(int_mpi) :: itag,send_cpu,recv_cpu
  !
  !.. Local Arrays ..
  _MPI_STATUS_TYPE_ :: istat(1:_MPI_STATUS_SIZE_)
  real(r8) :: exch_send(1:nq)
  real(r8) :: exch_recv(1:nq)
  !
continue
  !
  ! We have finished computing the gradients on the non-communication
  ! boundary faces and the communication of all other gradients between
  ! processors of adjacent partitions has completed. Now we have to
  ! correct any flux points on boundary faces that have conflicting
  ! boundary conditions.
  !
  ! This only needs to be done if we are using flux points that include
  ! the end points of the face (e.g., Lobatto) so skip all of this and
  ! return if we are not using a quadrature with flux points at the
  ! end points UNLESS override is present. If override is present, only
  ! exit return_loop if it is true. If it is not true, we dont care
  ! that it is present and we only care if we are using Lobatto points.
  !
  return_loop: do i = 1,1
    !
    if (use_bc_conflicts) then
      !
      if (present(override)) then
        if (override) exit return_loop
      end if
      !
      if (any(loc_flux_pts == Gauss_Lobatto)) exit return_loop
      !
    end if
    !
    ! Return if we have reached this point
    !
    return
    !
  end do return_loop
  !
  ! Now go through and exchange faceusp at bc conflicts
  !
  do n = 1,size(bc_conflict)
    !
    ! Continue to the next bc_conflict if this processor
    ! has zero involvment with this conflict.
    !
    if (bc_conflict(n)%not_my_conflict) cycle
    !
    ! Overwrite the boundary gradient depending on if this conflict
    ! is local or needs to be sent to another processor via MPI.
    !
    if (bc_conflict(n)%is_local) then
      !
      ! This conflict is local to this processor so just
      ! overwrite the boundary gradient in the slave flux points
      ! with the boundary gradient in the master flux point.
      !
      km = bc_conflict(n)%master%flux_point
      nfm = bc_conflict(n)%master%face
      !
      do i = 1,size(bc_conflict(n)%slave)
        !
        ks = bc_conflict(n)%slave(i)%flux_point
        nfs = bc_conflict(n)%slave(i)%face
        !
        faceusp(nfs)%v(:,ks,2) = faceusp(nfm)%v(:,km,2)
        !
      end do
      !
    else
      !
      ! This conflict is on a partition boundary so we need to
      ! use MPI to exchange the boundary gradients.
      !
      send_cpu = bc_conflict(n)%master%cpu
      !
      if (mypnum == send_cpu) then
        !
        km = bc_conflict(n)%master%flux_point
        nfm = bc_conflict(n)%master%face
        !
        exch_send = real(faceusp(nfm)%v(:,km,2),kind=r8)
        !
        do i = 1,size(bc_conflict(n)%slave)
          !
          recv_cpu = bc_conflict(n)%slave(i)%cpu
          !
          ! Off-chance that we have multiple slaves on a conflict
          ! but some are local and some need MPI communication
          !
          if (mypnum == recv_cpu) then
            !
            ks = bc_conflict(n)%slave(i)%flux_point
            nfs = bc_conflict(n)%slave(i)%face
            !
            faceusp(nfs)%v(:,ks,2) = real(exch_send,kind=wp)
            !
          else
            !
            itag = int( (send_cpu+n+1)*i + (recv_cpu+i+1)*n , kind=int_mpi )
            call mpi_send(exch_send,size(exch_send,kind=int_mpi),mpi_flttyp, &
                          recv_cpu,itag,MPI_COMM_WORLD,mpierr)
            !
          end if
          !
        end do
        !
      else
        !
        do i = 1,size(bc_conflict(n)%slave)
          !
          recv_cpu = bc_conflict(n)%slave(i)%cpu
          !
          if (mypnum == recv_cpu) then
            !
            itag = int( (send_cpu+n+1)*i + (recv_cpu+i+1)*n , kind=int_mpi )
            call mpi_recv(exch_recv,size(exch_recv,kind=int_mpi),mpi_flttyp, &
                          send_cpu,itag,MPI_COMM_WORLD,istat(1),mpierr)
            !
            ks = bc_conflict(n)%slave(i)%flux_point
            nfs = bc_conflict(n)%slave(i)%face
            !
            faceusp(nfs)%v(:,ks,2) = real(exch_recv,kind=wp)
            !
          end if
          !
        end do
        !
      end if
      !
    end if
    !
  end do
  !
end subroutine exchange_bc_conflicts_usp
!
!###############################################################################
!
subroutine exchange_bc_conflicts_dusp(facedusp,override)
  !
  !.. Use Statements ..
  use geovar,  only : nr,bc_conflict
  use ovar,    only : loc_flux_pts
  use ovar,    only : use_bc_conflicts
  use flowvar, only : face_vec
  !
  !.. Formal Arguments ..
  type(face_vec), intent(inout) :: facedusp(:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: override
  !
  !.. Local Scalars ..
  integer :: i,n,nf,km,ks,nfm,nfs
  integer(int_mpi) :: itag,send_cpu,recv_cpu
  !
  !.. Local Arrays ..
  _MPI_STATUS_TYPE_ :: istat(1:_MPI_STATUS_SIZE_)
  real(r8) :: exch_send(1:nr,1:nq)
  real(r8) :: exch_recv(1:nr,1:nq)
  !
continue
  !
  ! We have finished computing the gradients on the non-communication
  ! boundary faces and the communication of all other gradients between
  ! processors of adjacent partitions has completed. Now we have to
  ! correct any flux points on boundary faces that have conflicting
  ! boundary conditions.
  !
  ! This only needs to be done if we are using flux points that include
  ! the end points of the face (e.g., Lobatto) so skip all of this and
  ! return if we are not using a quadrature with flux points at the
  ! end points UNLESS override is present. If override is present, only
  ! exit return_loop if it is true. If it is not true, we dont care
  ! that it is present and we only care if we are using Lobatto points.
  !
  return_loop: do i = 1,1
    !
    if (use_bc_conflicts) then
      !
      if (present(override)) then
        if (override) exit return_loop
      end if
      !
      if (any(loc_flux_pts == Gauss_Lobatto)) exit return_loop
      !
    end if
    !
    ! Return if we have reached this point
    !
    return
    !
  end do return_loop
  !
  ! Now go through and exchange facedusp at bc conflicts
  !
  do n = 1,size(bc_conflict)
    !
    ! Continue to the next bc_conflict if this processor
    ! has zero involvment with this conflict.
    !
    if (bc_conflict(n)%not_my_conflict) cycle
    !
    ! Overwrite the boundary gradient depending on if this conflict
    ! is local or needs to be sent to another processor via MPI.
    !
    if (bc_conflict(n)%is_local) then
      !
      ! This conflict is local to this processor so just
      ! overwrite the boundary gradient in the slave flux points
      ! with the boundary gradient in the master flux point.
      !
      km = bc_conflict(n)%master%flux_point
      nfm = bc_conflict(n)%master%face
      !
      do i = 1,size(bc_conflict(n)%slave)
        !
        ks = bc_conflict(n)%slave(i)%flux_point
        nfs = bc_conflict(n)%slave(i)%face
        !
        facedusp(nfs)%v(:,:,ks,2) = facedusp(nfm)%v(:,:,km,2)
        !
      end do
      !
    else
      !
      ! This conflict is on a partition boundary so we need to
      ! use MPI to exchange the boundary gradients.
      !
      send_cpu = bc_conflict(n)%master%cpu
      !
      if (mypnum == send_cpu) then
        !
        km = bc_conflict(n)%master%flux_point
        nfm = bc_conflict(n)%master%face
        !
        exch_send = real(facedusp(nfm)%v(:,:,km,2),kind=r8)
        !
        do i = 1,size(bc_conflict(n)%slave)
          !
          recv_cpu = bc_conflict(n)%slave(i)%cpu
          !
          ! Off-chance that we have multiple slaves on a conflict
          ! but some are local and some need MPI communication
          !
          if (mypnum == recv_cpu) then
            !
            ks = bc_conflict(n)%slave(i)%flux_point
            nfs = bc_conflict(n)%slave(i)%face
            !
            facedusp(nfs)%v(:,:,ks,2) = real(exch_send,kind=wp)
            !
          else
            !
            itag = int( (send_cpu+n+1)*i + (recv_cpu+i+1)*n , kind=int_mpi )
            call mpi_send(exch_send,size(exch_send,kind=int_mpi),mpi_flttyp, &
                          recv_cpu,itag,MPI_COMM_WORLD,mpierr)
            !
          end if
          !
        end do
        !
      else
        !
        do i = 1,size(bc_conflict(n)%slave)
          !
          recv_cpu = bc_conflict(n)%slave(i)%cpu
          !
          if (mypnum == recv_cpu) then
            !
            itag = int( (send_cpu+n+1)*i + (recv_cpu+i+1)*n , kind=int_mpi )
            call mpi_recv(exch_recv,size(exch_recv,kind=int_mpi),mpi_flttyp, &
                          send_cpu,itag,MPI_COMM_WORLD,istat(1),mpierr)
            !
            ks = bc_conflict(n)%slave(i)%flux_point
            nfs = bc_conflict(n)%slave(i)%face
            !
            facedusp(nfs)%v(:,:,ks,2) = real(exch_recv,kind=wp)
            !
          end if
          !
        end do
        !
      end if
      !
    end if
    !
  end do
  !
end subroutine exchange_bc_conflicts_dusp
!
!###############################################################################
!
include "Functions/temperature.f90"
include "Functions/pressure.f90"
include "Functions/usp2v.f90"
include "Functions/vsp2u.f90"
!
!###############################################################################
!
end module module_bc
