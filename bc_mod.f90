module module_bc
  !
#include <mpi_defs.h>
!#define PROFILE_ON
  !
  ! Module to contain the procedures and variables related
  ! to the evaluating the boundary conditions for the flow.
  !
  !.. Use Statements ..
  use module_kind_types, global_debug_timer => debug_timer
  !
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq,ntk,ntl
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
  use metrics_mod,  only : facenrm
  use flowvar,      only : faceusp_array
  use flowvar,      only : facexyz_array
  use geovar,       only : bc_grps
  use ovar,         only : bc_in,itcur
  use ovar,         only : output_bface_array
!#define DEBUG_BC_GROUPS
#ifdef DEBUG_BC_GROUPS
  use geovar,       only : n_wall_flxpts,n_bnd_flxpts,n_flxpts,n_solpts
#endif
  !
  !.. Optional Arguments ..
  integer, optional, intent(out) :: exit_error
  !
  !.. Local Scalars ..
  integer :: i,j,n,m,n1,n2
  integer :: bc_idx,bc_type
  integer :: newton_failures
  !
  !.. Local Saved Variables ..
  logical(lk), save :: first_output = true
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_faceusp_bc"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  newton_failures = 0
  if (present(exit_error)) exit_error = 0
  !
  if (first_output) then
    if (output_bface_array) call output_bface
  end if
  !
  ! Loop over the boundary face groups and compute the conservative variables
  ! at the boundary flux points based on the face boundary condition.
  !
#ifdef DEBUG_BC_GROUPS
  if (first_output) then
    write (mypnum+500,101) n_wall_flxpts,n_bnd_flxpts,n_flxpts,n_solpts
    101 format (/,"In get_faceusp_bc:",//, &
                  "n_wall_flxpts = ",i0,/, &
                  " n_bnd_flxpts = ",i0,/, &
                  "     n_flxpts = ",i0,/, &
                  "     n_solpts = ",i0,/)
    flush (mypnum+500)
  end if
#endif
  do i = 1,size(bc_grps)
    !
    n1 = bc_grps(i)%beg_fp
    n2 = bc_grps(i)%end_fp
    !
    bc_type = bc_grps(i)%bc_type
    bc_idx = bc_grps(i)%bc_idx
#ifdef DEBUG_BC_GROUPS
    if (first_output) then
      write (mypnum+500,100) i,n1,i,n2,i,bc_type, &
                             trim(adjustl(bc_integer_to_string(bc_type))), &
                             i,bc_idx, &
            bc_idx,bc_idx, &
            lbound(bc_in(bc_idx)%cv,dim=1),ubound(bc_in(bc_idx)%cv,dim=1), &
            bc_idx,bc_idx, &
            lbound(bc_in(bc_idx)%pv,dim=1),ubound(bc_in(bc_idx)%pv,dim=1), &
            lbound(facenrm,dim=1),ubound(facenrm,dim=1), &
            lbound(facenrm,dim=2),ubound(facenrm,dim=2), &
            lbound(faceusp_array,dim=1),ubound(faceusp_array,dim=1), &
            lbound(faceusp_array,dim=2),ubound(faceusp_array,dim=2), &
            lbound(faceusp_array,dim=3),ubound(faceusp_array,dim=3), &
            bc_idx,(bc_in(bc_idx)%cv(j),j=1,nq), &
            bc_idx,(bc_in(bc_idx)%pv(j),j=1,nq)
      flush (mypnum+500)
      100 format (/, &
        "bc_grps(",i0,")%beg_fp  = ",i0,/, &
        "bc_grps(",i0,")%end_fp  = ",i0,/, &
        "bc_grps(",i0,")%bc_type = ",i0," (",a,")",/, &
        "bc_grps(",i0,")%bc_idx  = ",i0,/, &
        "dim(bc_in(",i0,")%cv) = bc_in(",i0,")%cv(",i0,":",i0,")",/, &
        "dim(bc_in(",i0,")%pv) = bc_in(",i0,")%pv(",i0,":",i0,")",/, &
        "dim(facenrm) = facenrm(",i0,":",i0,",",i0,":",i0,")",/, &
        "dim(faceusp) = faceusp(",i0,":",i0,",",i0,":",i0,",",i0,":",i0,")",/, &
        "bc_in(",i0,")%cv = ",5es14.6,/, &
        "bc_in(",i0,")%pv = ",5es14.6)
    end if
#endif
    !
    ! Compute the solution at the boundary flux points
    ! based on the specified boundary condition
    !
    select case (bc_type)
      !
      case (bc_freestream)
        !
        ! Free stream boundary condition
        !
        call get_bc_freestream_solution(n1,n2,faceusp_array(:,:,2), &
                                        bc_in(0)%cv)
      case (bc_fixed)
        !
        ! Fixed boundary condition
        !
        call get_bc_freestream_solution(n1,n2,faceusp_array(:,:,2), &
                                        bc_in(bc_idx)%cv)
        !
      case (bc_characteristic)
        !
        ! Characteristic boundary condition
        !
        call get_bc_characteristic_solution(n1,n2,facenrm,faceusp_array, &
                                            bc_in(bc_idx)%cv, &
                                            bc_in(bc_idx)%relax)
        !
      case (bc_generic_inflow)
        !
        ! Generic inflow boundary condition
        !
        call get_bc_generic_inflow_solution(n1,n2,facenrm,faceusp_array, &
                                            bc_in(bc_idx)%pv)
        !
      case (bc_generic_outflow)
        !
        ! Generic outflow boundary condition
        !
        call get_bc_generic_outflow_solution(n1,n2,facenrm,faceusp_array, &
                                             bc_in(bc_idx)%pv)
        !
      case (bc_generic_freeflow)
        !
        ! Generic freeflow boundary condition
        !
        call get_bc_generic_freeflow_solution(n1,n2,facenrm,faceusp_array, &
                                              bc_in(bc_idx)%cv)
        !
      case (bc_sub_inflow)
        !
        ! Subsonic inflow boundary condition
        !
        call get_bc_sub_inflow_solution(n1,n2,facenrm,faceusp_array, &
                                        bc_in(bc_idx)%cv, &
                                        bc_in(bc_idx)%pv, &
                                        newton_failures)
        !
      case (bc_sub_outflow)
        !
        ! Subsonic outflow boundary condition
        !
        call get_bc_sub_outflow_solution(n1,n2,facenrm,faceusp_array, &
                                         bc_in(bc_idx)%pv)
        !
      case (bc_sup_inflow)
        !
        ! Supersonic inflow boundary condition
        !
        call get_bc_sup_inflow_solution(n1,n2,faceusp_array(:,:,2), &
                                        bc_in(bc_idx)%cv)
        !
      case (bc_sup_outflow)
        !
        ! Supersonic outflow boundary condition
        !
        call get_bc_sup_outflow_solution(n1,n2,faceusp_array)
        !
      case (bc_slip_wall,bc_symmetry)
        !
        ! Slip wall / symmetry boundary condition
        !
        call get_bc_slip_wall_solution(n1,n2,facenrm,faceusp_array)
        !
      case (bc_adiabatic_wall)
        !
        ! No-slip adiabatic wall
        !
        call get_bc_adiabatic_wall_solution(n1,n2,faceusp_array, &
                                            bc_in(bc_idx)%pv)
        !
      case (bc_isothermal_wall)
        !
        ! No-slip isothermal wall
        !
        call get_bc_isothermal_wall_solution(n1,n2,faceusp_array, &
                                             bc_in(bc_idx)%wall_temperature, &
                                             bc_in(bc_idx)%pv)
        !
      case (bc_periodic)
        !
        ! Periodic boundary condition
        !
        call get_bc_periodic_solution(n1,n2,faceusp_array, &
                                      bc_grps(i)%periodic_list)
        !
      case (bc_mms_dirichlet)
        !
        ! Dirichlet boundary condition for MMS
        !
        call get_bc_mms_dirichlet_solution(n1,n2,facenrm, &
                                           facexyz_array, &
                                           faceusp_array)
        !
      case (bc_cpu_bnd)
        !
        ! We dont need to do anything here for communication boundaries
        !
      case default
        !
        ! This shouldn't be here and is an unknown boundary condition
        !
        write (iout,1) bc_type
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
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
  end if
  !
  first_output = fals
  !
#ifdef DEBUG_BC_GROUPS
 !call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,"after debugging bc_groups")
#endif
  !
 !call exchange_bc_conflicts_usp(faceusp)
  !
  call debug_timer(leaving_procedure,pname)
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
  use geovar,       only : bc_grps
  use metrics_mod,  only : facenrm
  use flowvar,      only : faceusp_array
  use flowvar,      only : facedusp_array
  !
  !.. Local Scalars ..
  integer :: i,n1,n2,bc_type
  !
  !.. Local Parameters ..
  integer, parameter :: local_bc_ignore(1:3) = [ bc_cpu_bnd, &
                                                 bc_unknown, &
                                                 not_a_bc ]
  character(len=*), parameter :: pname = "get_facedusp_bc"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Loop over the boundary face groups and compute the conservative variables
  ! at the boundary flux points based on the face boundary condition.
  !
  do i = 1,size(bc_grps)
    !
    n1 = bc_grps(i)%beg_fp
    n2 = bc_grps(i)%end_fp
    !
    bc_type = bc_grps(i)%bc_type
    !
    ! Compute the solution at the boundary flux points
    ! based on the specified boundary condition
    !
    if (any(bc_type == bc_walls)) then
      !
      ! Wall boundary condition
      !
      if (bc_type == bc_adiabatic_wall) then
        !
        ! No-slip adiabatic wall boundary conditions
        !
        call get_bc_adiabatic_wall_gradient(n1,n2, &
                                            faceusp_array(:,:,1), &
                                            facedusp_array)
        !
      else if (any(bc_type == [bc_slip_wall,bc_symmetry])) then
        !
        ! Slip wall / symmetry boundary condition
        !
        call get_bc_slip_wall_gradient(n1,n2,facenrm, &
                                       faceusp_array(:,:,1), &
                                       facedusp_array)
        !
      else
        !
        ! Non-adiabatic wall boundary condition
        !
        call get_bc_wall_gradient(n1,n2, &
                                  faceusp_array(:,:,1), &
                                  facedusp_array)
        !
      end if
      !
    else if (bc_type == bc_periodic) then
      !
      ! Periodic boundary condition
      !
      call get_bc_periodic_gradient(n1,n2,facedusp_array, &
                                    bc_grps(i)%periodic_list)
      !
    else if (all(bc_type /= local_bc_ignore)) then
      !
      ! General flow-through boundary conditions
      !
      call get_bc_flowthrough_gradient(n1,n2,facedusp_array)
      !
    end if
    !
  end do
  !
 !call exchange_bc_conflicts_dusp(facedusp)
  !
  call debug_timer(leaving_procedure,pname)
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
pure subroutine get_bc_freestream_solution(n1,n2,face_cv,cv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Free stream boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV : Conservative variables at the ghost cell flux points
  ! CV_IN   : Boundary condition flow variables in conservative form
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv_in(:)
  real(wp), intent(inout) :: face_cv(:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    call get_bc_freestream_solution_3d(n1,n2,face_cv,cv_in)
  else if (nr == 2) then
    call get_bc_freestream_solution_2d(n1,n2,face_cv,cv_in)
  else
    call get_bc_freestream_solution_nd(n1,n2,face_cv,cv_in)
  end if
  !
end subroutine get_bc_freestream_solution
!
!###############################################################################
!
pure subroutine get_bc_freestream_solution_3d(n1,n2,face_cv,cv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Free stream boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV : Conservative variables at the ghost cell flux points
  ! CV_IN   : Boundary condition flow variables in conservative form
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv_in(:)
  real(wp), intent(inout) :: face_cv(:,:)
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  ! Set all variables to free stream
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,cv_in)
  do n = n1,n2
    face_cv(n,nec) = cv_in(nec)
    face_cv(n,nmx) = cv_in(nmx)
    face_cv(n,nmy) = cv_in(nmy)
    face_cv(n,nmz) = cv_in(nmz)
    face_cv(n,nee) = cv_in(nee)
  end do
  !
end subroutine get_bc_freestream_solution_3d
!
!###############################################################################
!
pure subroutine get_bc_freestream_solution_2d(n1,n2,face_cv,cv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Free stream boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV : Conservative variables at the ghost cell flux points
  ! CV_IN   : Boundary condition flow variables in conservative form
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv_in(:)
  real(wp), intent(inout) :: face_cv(:,:)
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  ! Set all variables to free stream
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,cv_in)
  do n = n1,n2
    face_cv(n,nec) = cv_in(nec)
    face_cv(n,nmx) = cv_in(nmx)
    face_cv(n,nmy) = cv_in(nmy)
    face_cv(n,nee) = cv_in(nee)
  end do
  !
end subroutine get_bc_freestream_solution_2d
!
!###############################################################################
!
pure subroutine get_bc_freestream_solution_nd(n1,n2,face_cv,cv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Free stream boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV : Conservative variables at the ghost cell flux points
  ! CV_IN   : Boundary condition flow variables in conservative form
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv_in(:)
  real(wp), intent(inout) :: face_cv(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n
  !
continue
  !
  ! Set all variables to free stream
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) &
  !$ACC&              PRESENT(face_cv,cv_in)
  do m = 1,size(cv_in)
    do n = n1,n2
      face_cv(n,m) = cv_in(m)
    end do
  end do
  !
end subroutine get_bc_freestream_solution_nd
!
!###############################################################################
!
pure subroutine get_bc_characteristic_solution(n1,n2,facenrm,face_cv, &
                                               cv_in,relax)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Characteristic boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! CV_IN   : Boundary condition flow variables in conservative form
  ! RELAX   : Boundary condition relaxation parameters
  !
  !.. Use Statements ..
  use ovar, only : relaxation_t
  !
  !.. Formal Arguments ..
  integer,            intent(in) :: n1,n2
  real(wp),           intent(in) :: facenrm(:,:)
  real(wp),           intent(in) :: cv_in(:)
  type(relaxation_t), intent(in) :: relax
  real(wp),        intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
  !.. Local Parameters ..
  logical(lk), parameter :: veln_is_int = true
 !logical(lk), parameter :: veln_is_int = fals
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    call get_bc_characteristic_solution_3d(n1,n2,facenrm,face_cv, &
                                           cv_in,relax,veln_is_int)
  else if (nr == 2) then
    call get_bc_characteristic_solution_2d(n1,n2,facenrm,face_cv, &
                                           cv_in,relax,veln_is_int)
  else
    call get_bc_characteristic_solution_nd(n1,n2,facenrm,face_cv, &
                                           cv_in,relax,veln_is_int)
  end if
  !
end subroutine get_bc_characteristic_solution
!
!###############################################################################
!
pure subroutine get_bc_characteristic_solution_3d(n1,n2,facenrm,face_cv, &
                                                  cv_in,relax,veln_is_int)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Characteristic boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! CV_IN   : Boundary condition flow variables in conservative form
  ! RELAX   : Boundary condition relaxation parameters
  !
  !.. Use Statements ..
  use ovar, only : relaxation_t
  use ovar, only : itcur,time
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,            intent(in) :: n1,n2
  real(wp),           intent(in) :: facenrm(:,:)
  real(wp),           intent(in) :: cv_in(:)
  type(relaxation_t), intent(in) :: relax
  real(wp),        intent(inout) :: face_cv(:,:,:)
  logical(lk),        intent(in) :: veln_is_int
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: rlx,rlx_iter,rlx_time,gm1,gm1_inv
  real(wp) :: rho,vx,vy,vz,vn,pres,sc,rp,rm,unc,ac,vel
  real(wp) :: rho_int,rhou_int,rhov_int,rhow_int,rhoe_int
  real(wp) :: rho_rlx,rhou_rlx,rhov_rlx,rhow_rlx,rhoe_rlx
  real(wp) :: rho_ext,rhou_ext,rhov_ext,rhow_ext,rhoe_ext
  real(wp) :: rhoq2_int,pres_int,aspd_int,vn_int
  real(wp) :: rhoq2_rlx,pres_rlx,aspd_rlx,vn_rlx
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
  rho_ext = cv_in(nec)
  rhou_ext = cv_in(nmx)
  rhov_ext = cv_in(nmy)
  rhow_ext = cv_in(nmz)
  rhoe_ext = cv_in(nee)
  !
  ! Compute the boundary solution at each flux point
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
  do n = n1,n2
    !
    ! Unit normal at this flux point
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    rnz = facenrm(3,n)
    ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
    unx = rnx/ds
    uny = rny/ds
    unz = rnz/ds
    !
    ! Get the host cell variables at this flux point
    !
    rho_int = face_cv(n,nec,1)
    rhou_int = face_cv(n,nmx,1)
    rhov_int = face_cv(n,nmy,1)
    rhow_int = face_cv(n,nmz,1)
    rhoe_int = face_cv(n,nee,1)
    rhoq2_int = rhou_int*rhou_int + rhov_int*rhov_int + rhow_int*rhow_int
    pres_int = gm1*(rhoe_int - half*rhoq2_int/rho_int)
    !
    ! Compute the relaxed inflow condition
    !
    rho_rlx  = (one-rlx)*rho_int  + rlx*rho_ext
    rhou_rlx = (one-rlx)*rhou_int + rlx*rhou_ext
    rhov_rlx = (one-rlx)*rhov_int + rlx*rhov_ext
    rhow_rlx = (one-rlx)*rhow_int + rlx*rhow_ext
    rhoe_rlx = (one-rlx)*rhoe_int + rlx*rhoe_ext
    rhoq2_rlx = rhou_rlx*rhou_rlx + rhov_rlx*rhov_rlx + rhow_rlx*rhow_rlx
    pres_rlx = gm1*(rhoe_rlx - half*rhoq2_rlx/rho_rlx)
    !
    ! Compute the speed of sound in the host cell
    !
    aspd_int = sqrt( max(gam*pres_int/rho_int,rndoff) ) ! Int. speed of sound
    aspd_rlx = sqrt( max(gam*pres_rlx/rho_rlx,rndoff) ) ! Ext. speed of sound
    !
    ! Compute the normal velocity at this flux point
    !
    ! Interior normal velocity
    vn_int = (rhou_int*unx + rhov_int*uny + rhow_int*unz) / rho_int
    ! Exterior normal velocity
    vn_rlx = (rhou_rlx*unx + rhov_rlx*uny + rhow_rlx*unz) / rho_rlx
    !
    ! Apply boundary conditions based on the
    ! normal Mach number to the boundary face
    !
    rp = vn_int + two*aspd_int*gm1_inv ! Outgoing Riemann invariant
    rm = vn_rlx - two*aspd_rlx*gm1_inv ! Incoming Reimann invariant
    !
    ! Supersonic inflow boundary
    !    Set R+ so that 0.5*(R+ + R-) = vn_rlx
    !
    if (vn_rlx <= -aspd_rlx) then
      rp = vn_rlx + two*aspd_rlx*gm1_inv
    end if
    !
    ! Supersonic outflow boundary
    !    Set R- so that 0.5*(R+ + R-) = vn_int
    !
    if (vn_int >= aspd_int) then
      rm = vn_int - two*aspd_int*gm1_inv
    end if
    !
    unc = half*(rp+rm)
    ac  = one4*gm1*(rp-rm)
    !
    vel = merge( vn_int , unc , veln_is_int )
    if (vel > zero) then ! outflow boundary
      vx = rhou_int/rho_int
      vy = rhov_int/rho_int
      vz = rhow_int/rho_int
      vn = vn_int
      sc = pres_int / (rho_int**gam)
     !gufp(nmb:nme,k) = hvfp(nmb:nme) + (unc - vn_int)*rnrm(:)
    else                 ! inflow boundary
      vx = rhou_rlx/rho_rlx
      vy = rhov_rlx/rho_rlx
      vz = rhow_rlx/rho_rlx
      vn = vn_rlx
      sc = pres_rlx / (rho_rlx**gam)
     !gufp(nmb:nme,k) = pv_in_rlx(nmb:nme) + (unc - vn_rlx)*rnrm(:)
    end if
    !
    rho = ( ac * ac / (gam*sc) )**gm1_inv
    pres = sc * (rho**gam)
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rho*(vx + (unc - vn)*unx)
    face_cv(n,nmy,2) = rho*(vy + (unc - vn)*uny)
    face_cv(n,nmz,2) = rho*(vz + (unc - vn)*unz)
    face_cv(n,nee,2) = pres*gm1_inv + half*rho*(vx*vx + vy*vy + vz*vz)
    !
  end do
  !
end subroutine get_bc_characteristic_solution_3d
!
!###############################################################################
!
pure subroutine get_bc_characteristic_solution_2d(n1,n2,facenrm,face_cv, &
                                                  cv_in,relax,veln_is_int)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Characteristic boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! CV_IN   : Boundary condition flow variables in conservative form
  ! RELAX   : Boundary condition relaxation parameters
  !
  !.. Use Statements ..
  use ovar, only : relaxation_t
  use ovar, only : itcur,time
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,            intent(in) :: n1,n2
  real(wp),           intent(in) :: facenrm(:,:)
  real(wp),           intent(in) :: cv_in(:)
  type(relaxation_t), intent(in) :: relax
  real(wp),        intent(inout) :: face_cv(:,:,:)
  logical(lk),        intent(in) :: veln_is_int
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,rny,unx,uny
  real(wp) :: rlx,rlx_iter,rlx_time,gm1,gm1_inv
  real(wp) :: rho,vx,vy,vn,pres,sc,rp,rm,unc,ac,vel
  real(wp) :: rho_int,rhou_int,rhov_int,rhoe_int
  real(wp) :: rho_rlx,rhou_rlx,rhov_rlx,rhoe_rlx
  real(wp) :: rho_ext,rhou_ext,rhov_ext,rhoe_ext
  real(wp) :: rhoq2_int,pres_int,aspd_int,vn_int
  real(wp) :: rhoq2_rlx,pres_rlx,aspd_rlx,vn_rlx
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
  rho_ext = cv_in(nec)
  rhou_ext = cv_in(nmx)
  rhov_ext = cv_in(nmy)
  rhoe_ext = cv_in(nee)
  !
  ! Compute the boundary solution at each flux point
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
  do n = n1,n2
    !
    ! Unit normal at this flux point
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds = sqrt( rnx*rnx + rny*rny )
    unx = rnx/ds
    uny = rny/ds
    !
    ! Get the host cell variables at this flux point
    !
    rho_int = face_cv(n,nec,1)
    rhou_int = face_cv(n,nmx,1)
    rhov_int = face_cv(n,nmy,1)
    rhoe_int = face_cv(n,nee,1)
    rhoq2_int = rhou_int*rhou_int + rhov_int*rhov_int
    pres_int = gm1*(rhoe_int - half*rhoq2_int/rho_int)
    !
    ! Compute the relaxed inflow condition
    !
    rho_rlx  = (one-rlx)*rho_int  + rlx*rho_ext
    rhou_rlx = (one-rlx)*rhou_int + rlx*rhou_ext
    rhov_rlx = (one-rlx)*rhov_int + rlx*rhov_ext
    rhoe_rlx = (one-rlx)*rhoe_int + rlx*rhoe_ext
    rhoq2_rlx = rhou_rlx*rhou_rlx + rhov_rlx*rhov_rlx
    pres_rlx = gm1*(rhoe_rlx - half*rhoq2_rlx/rho_rlx)
    !
    ! Compute the speed of sound in the host cell
    !
    aspd_int = sqrt( max(gam*pres_int/rho_int,rndoff) ) ! Int. speed of sound
    aspd_rlx = sqrt( max(gam*pres_rlx/rho_rlx,rndoff) ) ! Ext. speed of sound
    !
    ! Compute the normal velocity at this flux point
    !
    vn_int = (rhou_int*unx + rhov_int*uny) / rho_int ! Interior normal velocity
    vn_rlx = (rhou_rlx*unx + rhov_rlx*uny) / rho_rlx ! Exterior normal velocity
    !
    ! Apply boundary conditions based on the
    ! normal Mach number to the boundary face
    !
    rp = vn_int + two*aspd_int*gm1_inv ! Outgoing Riemann invariant
    rm = vn_rlx - two*aspd_rlx*gm1_inv ! Incoming Reimann invariant
    !
    ! Supersonic inflow boundary
    !    Set R+ so that 0.5*(R+ + R-) = vn_rlx
    !
    if (vn_rlx <= -aspd_rlx) then
      rp = vn_rlx + two*aspd_rlx*gm1_inv
    end if
    !
    ! Supersonic outflow boundary
    !    Set R- so that 0.5*(R+ + R-) = vn_int
    !
    if (vn_int >= aspd_int) then
      rm = vn_int - two*aspd_int*gm1_inv
    end if
    !
    unc = half*(rp+rm)
    ac  = one4*gm1*(rp-rm)
    !
    vel = merge( vn_int , unc , veln_is_int )
    if (vel > zero) then ! outflow boundary
      vx = rhou_int/rho_int
      vy = rhov_int/rho_int
      vn = vn_int
      sc = pres_int / (rho_int**gam)
     !gufp(nmb:nme,k) = hvfp(nmb:nme) + (unc - vn_int)*rnrm(:)
    else                 ! inflow boundary
      vx = rhou_rlx/rho_rlx
      vy = rhov_rlx/rho_rlx
      vn = vn_rlx
      sc = pres_rlx / (rho_rlx**gam)
     !gufp(nmb:nme,k) = pv_in_rlx(nmb:nme) + (unc - vn_rlx)*rnrm(:)
    end if
    !
    rho = ( ac * ac / (gam*sc) )**gm1_inv
    pres = sc * (rho**gam)
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rho*(vx + (unc - vn)*unx)
    face_cv(n,nmy,2) = rho*(vy + (unc - vn)*uny)
    face_cv(n,nee,2) = pres*gm1_inv + half*rho*(vx*vx + vy*vy)
    !
  end do
  !
end subroutine get_bc_characteristic_solution_2d
!
!###############################################################################
!
pure subroutine get_bc_characteristic_solution_nd(n1,n2,facenrm,face_cv, &
                                                  cv_in,relax,veln_is_int)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Characteristic boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! CV_IN   : Boundary condition flow variables in conservative form
  ! RELAX   : Boundary condition relaxation parameters
  !
  !.. Use Statements ..
  use ovar, only : relaxation_t
  use ovar, only : itcur,time
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,            intent(in) :: n1,n2
  real(wp),           intent(in) :: facenrm(:,:)
  real(wp),           intent(in) :: cv_in(:)
  type(relaxation_t), intent(in) :: relax
  real(wp),        intent(inout) :: face_cv(:,:,:)
  logical(lk),        intent(in) :: veln_is_int
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: rlx,rlx_iter,rlx_time,gm1,gm1_inv
  real(wp) :: rho,vx,vy,vz,vn,pres,sc,rp,rm,unc,ac,vel
  real(wp) :: rho_int,rhou_int,rhov_int,rhow_int,rhoe_int
  real(wp) :: rho_rlx,rhou_rlx,rhov_rlx,rhow_rlx,rhoe_rlx
  real(wp) :: rho_ext,rhou_ext,rhov_ext,rhow_ext,rhoe_ext
  real(wp) :: rhoq2_int,pres_int,aspd_int,vn_int
  real(wp) :: rhoq2_rlx,pres_rlx,aspd_rlx,vn_rlx
  !
continue
  !
  nr = size(facenrm,dim=1)
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
  rho_ext = cv_in(nec)
  if (nr == 3) then
    rhou_ext = cv_in(nmx)
    rhov_ext = cv_in(nmy)
    rhow_ext = cv_in(nmz)
  else if (nr == 2) then
    rhou_ext = cv_in(nmx)
    rhov_ext = cv_in(nmy)
  else
    rhou_ext = cv_in(nmx)
  end if
  rhoe_ext = cv_in(nee)
  !
  ! Compute the boundary solution at each flux point
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
    do n = n1,n2
      !
      ! Unit normal at this flux point
      !
      rnx = facenrm(1,n)
      rny = facenrm(2,n)
      rnz = facenrm(3,n)
      ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
      unx = rnx/ds
      uny = rny/ds
      unz = rnz/ds
      !
      ! Get the host cell variables at this flux point
      !
      rho_int = face_cv(n,nec,1)
      rhou_int = face_cv(n,nmx,1)
      rhov_int = face_cv(n,nmy,1)
      rhow_int = face_cv(n,nmz,1)
      rhoe_int = face_cv(n,nee,1)
      rhoq2_int = rhou_int*rhou_int + rhov_int*rhov_int + rhow_int*rhow_int
      pres_int = gm1*(rhoe_int - half*rhoq2_int/rho_int)
      !
      ! Compute the relaxed inflow condition
      !
      rho_rlx  = (one-rlx)*rho_int  + rlx*rho_ext
      rhou_rlx = (one-rlx)*rhou_int + rlx*rhou_ext
      rhov_rlx = (one-rlx)*rhov_int + rlx*rhov_ext
      rhow_rlx = (one-rlx)*rhow_int + rlx*rhow_ext
      rhoe_rlx = (one-rlx)*rhoe_int + rlx*rhoe_ext
      rhoq2_rlx = rhou_rlx*rhou_rlx + rhov_rlx*rhov_rlx + rhow_rlx*rhow_rlx
      pres_rlx = gm1*(rhoe_rlx - half*rhoq2_rlx/rho_rlx)
      !
      ! Compute the speed of sound in the host cell
      !
      aspd_int = sqrt( max(gam*pres_int/rho_int,rndoff) ) ! Int. speed of sound
      aspd_rlx = sqrt( max(gam*pres_rlx/rho_rlx,rndoff) ) ! Ext. speed of sound
      !
      ! Compute the normal velocity at this flux point
      !
      ! Interior normal velocity
      vn_int = (rhou_int*unx + rhov_int*uny + rhow_int*unz) / rho_int
      ! Exterior normal velocity
      vn_rlx = (rhou_rlx*unx + rhov_rlx*uny + rhow_rlx*unz) / rho_rlx
      !
      ! Apply boundary conditions based on the
      ! normal Mach number to the boundary face
      !
      rp = vn_int + two*aspd_int*gm1_inv ! Outgoing Riemann invariant
      rm = vn_rlx - two*aspd_rlx*gm1_inv ! Incoming Reimann invariant
      !
      ! Supersonic inflow boundary
      !    Set R+ so that 0.5*(R+ + R-) = vn_rlx
      !
      if (vn_rlx <= -aspd_rlx) then
        rp = vn_rlx + two*aspd_rlx*gm1_inv
      end if
      !
      ! Supersonic outflow boundary
      !    Set R- so that 0.5*(R+ + R-) = vn_int
      !
      if (vn_int >= aspd_int) then
        rm = vn_int - two*aspd_int*gm1_inv
      end if
      !
      unc = half*(rp+rm)
      ac  = one4*gm1*(rp-rm)
      !
      vel = merge( vn_int , unc , veln_is_int )
      if (vel > zero) then ! outflow boundary
        vx = rhou_int/rho_int
        vy = rhov_int/rho_int
        vz = rhow_int/rho_int
        vn = vn_int
        sc = pres_int / (rho_int**gam)
       !gufp(nmb:nme,k) = hvfp(nmb:nme) + (unc - vn_int)*rnrm(:)
      else                 ! inflow boundary
        vx = rhou_rlx/rho_rlx
        vy = rhov_rlx/rho_rlx
        vz = rhow_rlx/rho_rlx
        vn = vn_rlx
        sc = pres_rlx / (rho_rlx**gam)
       !gufp(nmb:nme,k) = pv_in_rlx(nmb:nme) + (unc - vn_rlx)*rnrm(:)
      end if
      !
      rho = ( ac * ac / (gam*sc) )**gm1_inv
      pres = sc * (rho**gam)
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*(vx + (unc - vn)*unx)
      face_cv(n,nmy,2) = rho*(vy + (unc - vn)*uny)
      face_cv(n,nmz,2) = rho*(vz + (unc - vn)*unz)
      face_cv(n,nee,2) = pres*gm1_inv + half*rho*(vx*vx + vy*vy + vz*vz)
      !
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
    do n = n1,n2
      !
      ! Unit normal at this flux point
      !
      rnx = facenrm(1,n)
      rny = facenrm(2,n)
      ds = sqrt( rnx*rnx + rny*rny )
      unx = rnx/ds
      uny = rny/ds
      !
      ! Get the host cell variables at this flux point
      !
      rho_int = face_cv(n,nec,1)
      rhou_int = face_cv(n,nmx,1)
      rhov_int = face_cv(n,nmy,1)
      rhoe_int = face_cv(n,nee,1)
      rhoq2_int = rhou_int*rhou_int + rhov_int*rhov_int
      pres_int = gm1*(rhoe_int - half*rhoq2_int/rho_int)
      !
      ! Compute the relaxed inflow condition
      !
      rho_rlx  = (one-rlx)*rho_int  + rlx*rho_ext
      rhou_rlx = (one-rlx)*rhou_int + rlx*rhou_ext
      rhov_rlx = (one-rlx)*rhov_int + rlx*rhov_ext
      rhoe_rlx = (one-rlx)*rhoe_int + rlx*rhoe_ext
      rhoq2_rlx = rhou_rlx*rhou_rlx + rhov_rlx*rhov_rlx
      pres_rlx = gm1*(rhoe_rlx - half*rhoq2_rlx/rho_rlx)
      !
      ! Compute the speed of sound in the host cell
      !
      aspd_int = sqrt( max(gam*pres_int/rho_int,rndoff) ) ! Int. speed of sound
      aspd_rlx = sqrt( max(gam*pres_rlx/rho_rlx,rndoff) ) ! Ext. speed of sound
      !
      ! Compute the normal velocity at this flux point
      !
      vn_int = (rhou_int*unx + rhov_int*uny) / rho_int ! Int. normal velocity
      vn_rlx = (rhou_rlx*unx + rhov_rlx*uny) / rho_rlx ! Ext. normal velocity
      !
      ! Apply boundary conditions based on the
      ! normal Mach number to the boundary face
      !
      rp = vn_int + two*aspd_int*gm1_inv ! Outgoing Riemann invariant
      rm = vn_rlx - two*aspd_rlx*gm1_inv ! Incoming Reimann invariant
      !
      ! Supersonic inflow boundary
      !    Set R+ so that 0.5*(R+ + R-) = vn_rlx
      !
      if (vn_rlx <= -aspd_rlx) then
        rp = vn_rlx + two*aspd_rlx*gm1_inv
      end if
      !
      ! Supersonic outflow boundary
      !    Set R- so that 0.5*(R+ + R-) = vn_int
      !
      if (vn_int >= aspd_int) then
        rm = vn_int - two*aspd_int*gm1_inv
      end if
      !
      unc = half*(rp+rm)
      ac  = one4*gm1*(rp-rm)
      !
      vel = merge( vn_int , unc , veln_is_int )
      if (vel > zero) then ! outflow boundary
        vx = rhou_int/rho_int
        vy = rhov_int/rho_int
        vn = vn_int
        sc = pres_int / (rho_int**gam)
       !gufp(nmb:nme,k) = hvfp(nmb:nme) + (unc - vn_int)*rnrm(:)
      else                 ! inflow boundary
        vx = rhou_rlx/rho_rlx
        vy = rhov_rlx/rho_rlx
        vn = vn_rlx
        sc = pres_rlx / (rho_rlx**gam)
       !gufp(nmb:nme,k) = pv_in_rlx(nmb:nme) + (unc - vn_rlx)*rnrm(:)
      end if
      !
      rho = ( ac * ac / (gam*sc) )**gm1_inv
      pres = sc * (rho**gam)
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*(vx + (unc - vn)*unx)
      face_cv(n,nmy,2) = rho*(vy + (unc - vn)*uny)
      face_cv(n,nee,2) = pres*gm1_inv + half*rho*(vx*vx + vy*vy)
      !
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
    do n = n1,n2
      !
      ! Unit normal at this flux point
      !
      unx = facenrm(1,n)/abs(facenrm(1,n))
      !
      ! Get the host cell variables at this flux point
      !
      rho_int = face_cv(n,nec,1)
      rhou_int = face_cv(n,nmx,1)
      rhoe_int = face_cv(n,nee,1)
      pres_int = gm1*(rhoe_int - half*rhou_int*rhou_int/rho_int)
      !
      ! Compute the relaxed inflow condition
      !
      rho_rlx  = (one-rlx)*rho_int  + rlx*rho_ext
      rhou_rlx = (one-rlx)*rhou_int + rlx*rhou_ext
      rhoe_rlx = (one-rlx)*rhoe_int + rlx*rhoe_ext
      pres_rlx = gm1*(rhoe_rlx - half*rhou_rlx*rhou_rlx/rho_rlx)
      !
      ! Compute the speed of sound in the host cell
      !
      aspd_int = sqrt( max(gam*pres_int/rho_int,rndoff) ) ! Int. speed of sound
      aspd_rlx = sqrt( max(gam*pres_rlx/rho_rlx,rndoff) ) ! Ext. speed of sound
      !
      ! Compute the normal velocity at this flux point
      !
      vn_int = rhou_int*unx/rho_int ! Interior normal velocity
      vn_rlx = rhou_rlx*unx/rho_rlx ! Exterior normal velocity
      !
      ! Apply boundary conditions based on the
      ! normal Mach number to the boundary face
      !
      rp = vn_int + two*aspd_int*gm1_inv ! Outgoing Riemann invariant
      rm = vn_rlx - two*aspd_rlx*gm1_inv ! Incoming Reimann invariant
      !
      ! Supersonic inflow boundary
      !    Set R+ so that 0.5*(R+ + R-) = vn_rlx
      !
      if (vn_rlx <= -aspd_rlx) then
        rp = vn_rlx + two*aspd_rlx*gm1_inv
      end if
      !
      ! Supersonic outflow boundary
      !    Set R- so that 0.5*(R+ + R-) = vn_int
      !
      if (vn_int >= aspd_int) then
        rm = vn_int - two*aspd_int*gm1_inv
      end if
      !
      unc = half*(rp+rm)
      ac  = one4*gm1*(rp-rm)
      !
      vel = merge( vn_int , unc , veln_is_int )
      if (vel > zero) then ! outflow boundary
        vx = rhou_int/rho_int
        vn = vn_int
        sc = pres_int / (rho_int**gam)
       !gufp(nmb:nme,k) = hvfp(nmb:nme) + (unc - vn_int)*rnrm(:)
      else                 ! inflow boundary
        vx = rhou_rlx/rho_rlx
        vn = vn_rlx
        sc = pres_rlx / (rho_rlx**gam)
       !gufp(nmb:nme,k) = pv_in_rlx(nmb:nme) + (unc - vn_rlx)*rnrm(:)
      end if
      !
      rho = ( ac * ac / (gam*sc) )**gm1_inv
      pres = sc * (rho**gam)
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*(vx + (unc - vn)*unx)
      face_cv(n,nee,2) = pres*gm1_inv + half*rho*vx*vx
      !
    end do
  end if
  !
end subroutine get_bc_characteristic_solution_nd
!
!###############################################################################
!
pure subroutine get_bc_generic_inflow_solution(n1,n2,facenrm,face_cv,pv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Supersonic inflow boundary condition
  ! ##    - Fix all variables to the inflow condition
  ! #####################################################
  ! #####################################################
  !
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    call get_bc_generic_inflow_solution_3d(n1,n2,facenrm,face_cv,pv_in)
  else if (nr == 2) then
    call get_bc_generic_inflow_solution_2d(n1,n2,facenrm,face_cv,pv_in)
  else
    call get_bc_generic_inflow_solution_nd(n1,n2,facenrm,face_cv,pv_in)
  end if
  !
end subroutine get_bc_generic_inflow_solution
!
!###############################################################################
!
pure subroutine get_bc_generic_inflow_solution_3d(n1,n2,facenrm,face_cv,pv_in)
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
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,gm1_inv
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: rho,rhou,rhov,rhow,rhoe,rhoq2,pres,mach
  real(wp) :: rho_ext,vx_ext,vy_ext,vz_ext,pres_ext,q2_ext
  real(wp) :: ske_ext,rhoe_ext,a2_ext,total_cond,p0_ext
  !
continue
  !
  ! Precompute the reference specific kinetic energy
  ! and the inverse of gamma-1
  !
  gm1 = gam - one
  gm1_inv = one/gm1
  !
  rho_ext = pv_in(nec)
  pres_ext = pv_in(nee)
  !
  vx_ext = pv_in(nmx)
  vy_ext = pv_in(nmy)
  vz_ext = pv_in(nmz)
  q2_ext = vx_ext*vx_ext + vy_ext*vy_ext + vz_ext*vz_ext
  !
  ske_ext = half*rho_ext*q2_ext
  rhoe_ext = pres_ext*gm1_inv + ske_ext
  !
  a2_ext = max( gam*pres_ext/rho_ext , rndoff )
  !
  total_cond = one + half*gm1*(q2_ext/a2_ext)
  !
  p0_ext = pres_ext * total_cond**(gam*gm1_inv)
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
  do n = n1,n2
    !
    ! Unit normal at this flux point
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    rnz = facenrm(3,n)
    ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
    unx = rnx/ds
    uny = rny/ds
    unz = rnz/ds
    !
    ! Compute the Mach number of the normal velocity at this flux point
    !
    rho = face_cv(n,nec,1)
    rhou = face_cv(n,nmx,1)
    rhov = face_cv(n,nmy,1)
    rhow = face_cv(n,nmz,1)
    rhoq2 = rhou*rhou + rhov*rhov + rhow*rhow
    pres = gm1*(face_cv(n,nee,1) - half*rhoq2/rho)
    !
    ! NOTE: The speed of sound (denominator) is gam*pres*rho instead of
    !       gam*pres/rho because the normal velocity (numerator) is rho*vn
    !       instead of just vn.
    !
    mach = (rhou*unx + rhov*uny + rhow*unz) / sqrt( max(gam*pres*rho,rndoff) )
    !
    ! Initialize the BC state to the reference conditions
    !
    rho = rho_ext
    rhou = rho_ext*vx_ext
    rhov = rho_ext*vy_ext
    rhow = rho_ext*vz_ext
    rhoe = rhoe_ext
    !
    ! If mach <= -one, do nothing because in that case we want to use the
    ! the reference pressure for the boundary, so leave rhoe set to the
    ! reference value.
    !
    ! NOTE: The following prescribed BCs are probably wrong if the normal
    !       velocity indicates that this is an outflow instead of the
    !       expected inflow. However, the precribed BCs should be correct
    !       if the normal velocity indicates this is a subsonic inflow BC.
    !
    if (mach > -one) then
      !
      ! Subsonic inflow
      !
      ! If the normal velocity from the interior says this is an outflow
      ! boundary, set the ghost velocity to the magnitude of the reference
      ! velocity multiplied by the unit normal
      ! NOTE: q2 will be the same as q2_ext, and thus ske
      !       will be the same as ske_ext
      !
      if (mach > zero) then
        rhou = rho_ext*sqrt(q2_ext)*unx
        rhov = rho_ext*sqrt(q2_ext)*uny
        rhow = rho_ext*sqrt(q2_ext)*unz
      end if
      !
      ! Set the boundary pressure to the interior pressure if the normal
      ! velocity from the interior is subsonic, otherwise set the boundary
      ! pressure to the reference total pressure
      !
      ! If the normal velocity from the interior is supersonic, set the
      ! boundary pressure to the reference total pressure. Otherwise, leave
      ! the boundary pressure set to the interior value.
      !
      if (mach >= one) then
        pres = p0_ext
      end if
      !
      rhoe = pres*gm1_inv + ske_ext
      !
    end if
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rhou
    face_cv(n,nmy,2) = rhov
    face_cv(n,nmz,2) = rhow
    face_cv(n,nee,2) = rhoe
    !
  end do
  !
end subroutine get_bc_generic_inflow_solution_3d
!
!###############################################################################
!
pure subroutine get_bc_generic_inflow_solution_2d(n1,n2,facenrm,face_cv,pv_in)
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
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,gm1_inv
  real(wp) :: ds,rnx,rny,unx,uny
  real(wp) :: rho,rhou,rhov,rhoe,rhoq2,pres,mach
  real(wp) :: rho_ext,vx_ext,vy_ext,pres_ext,q2_ext
  real(wp) :: ske_ext,rhoe_ext,a2_ext,total_cond,p0_ext
  !
continue
  !
  ! Precompute the reference specific kinetic energy
  ! and the inverse of gamma-1
  !
  gm1 = gam - one
  gm1_inv = one/gm1
  !
  rho_ext = pv_in(nec)
  pres_ext = pv_in(nee)
  !
  vx_ext = pv_in(nmx)
  vy_ext = pv_in(nmy)
  q2_ext = vx_ext*vx_ext + vy_ext*vy_ext
  !
  ske_ext = half*rho_ext*q2_ext
  rhoe_ext = pres_ext*gm1_inv + ske_ext
  !
  a2_ext = max( gam*pres_ext/rho_ext , rndoff )
  !
  total_cond = one + half*gm1*(q2_ext/a2_ext)
  !
  p0_ext = pres_ext * total_cond**(gam*gm1_inv)
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
  do n = n1,n2
    !
    ! Unit normal at this flux point
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds = sqrt( rnx*rnx + rny*rny )
    unx = rnx/ds
    uny = rny/ds
    !
    ! Compute the Mach number of the normal velocity at this flux point
    !
    rho = face_cv(n,nec,1)
    rhou = face_cv(n,nmx,1)
    rhov = face_cv(n,nmy,1)
    rhoq2 = rhou*rhou + rhov*rhov
    pres = gm1*(face_cv(n,nee,1) - half*rhoq2/rho)
    !
    ! NOTE: The speed of sound (denominator) is gam*pres*rho instead of
    !       gam*pres/rho because the normal velocity (numerator) is rho*vn
    !       instead of just vn.
    !
    mach = (rhou*unx + rhov*uny) / sqrt( max(gam*pres*rho,rndoff) )
    !
    ! Initialize the BC state to the reference conditions
    !
    rho = rho_ext
    rhou = rho_ext*vx_ext
    rhov = rho_ext*vy_ext
    rhoe = rhoe_ext
    !
    ! If mach <= -one, do nothing because in that case we want to use the
    ! the reference pressure for the boundary, so leave rhoe set to the
    ! reference value.
    !
    ! NOTE: The following prescribed BCs are probably wrong if the normal
    !       velocity indicates that this is an outflow instead of the
    !       expected inflow. However, the precribed BCs should be correct
    !       if the normal velocity indicates this is a subsonic inflow BC.
    !
    if (mach > -one) then
      !
      ! Subsonic inflow
      !
      ! If the normal velocity from the interior says this is an outflow
      ! boundary, set the ghost velocity to the magnitude of the reference
      ! velocity multiplied by the unit normal
      ! NOTE: q2 will be the same as q2_ext, and thus ske
      !       will be the same as ske_ext
      !
      if (mach > zero) then
        rhou = rho_ext*sqrt(q2_ext)*unx
        rhov = rho_ext*sqrt(q2_ext)*uny
      end if
      !
      ! Set the boundary pressure to the interior pressure if the normal
      ! velocity from the interior is subsonic, otherwise set the boundary
      ! pressure to the reference total pressure
      !
      ! If the normal velocity from the interior is supersonic, set the
      ! boundary pressure to the reference total pressure. Otherwise, leave
      ! the boundary pressure set to the interior value.
      !
      if (mach >= one) then
        pres = p0_ext
      end if
      !
      rhoe = pres*gm1_inv + ske_ext
      !
    end if
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rhou
    face_cv(n,nmy,2) = rhov
    face_cv(n,nee,2) = rhoe
    !
  end do
  !
end subroutine get_bc_generic_inflow_solution_2d
!
!###############################################################################
!
pure subroutine get_bc_generic_inflow_solution_nd(n1,n2,facenrm,face_cv,pv_in)
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
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: gm1,gm1_inv
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: rho,rhou,rhov,rhow,rhoe,rhoq2,pres,mach
  real(wp) :: rho_ext,vx_ext,vy_ext,vz_ext,pres_ext,q2_ext
  real(wp) :: ske_ext,rhoe_ext,a2_ext,total_cond,p0_ext
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  ! Precompute the reference specific kinetic energy
  ! and the inverse of gamma-1
  !
  gm1 = gam - one
  gm1_inv = one/gm1
  !
  rho_ext = pv_in(nec)
  pres_ext = pv_in(nee)
  !
  if (nr == 3) then
    vx_ext = pv_in(nmx)
    vy_ext = pv_in(nmy)
    vz_ext = pv_in(nmz)
    q2_ext = vx_ext*vx_ext + vy_ext*vy_ext + vz_ext*vz_ext
  else if (nr == 2) then
    vx_ext = pv_in(nmx)
    vy_ext = pv_in(nmy)
    q2_ext = vx_ext*vx_ext + vy_ext*vy_ext
  else
    vx_ext = pv_in(nmx)
    q2_ext = vx_ext*vx_ext
  end if
  !
  ske_ext = half*rho_ext*q2_ext
  rhoe_ext = pres_ext*gm1_inv + ske_ext
  !
  a2_ext = max( gam*pres_ext/rho_ext , rndoff )
  !
  total_cond = one + half*gm1*(q2_ext/a2_ext)
  !
  p0_ext = pres_ext * total_cond**(gam*gm1_inv)
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
    do n = n1,n2
      !
      ! Unit normal at this flux point
      !
      rnx = facenrm(1,n)
      rny = facenrm(2,n)
      rnz = facenrm(3,n)
      ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
      unx = rnx/ds
      uny = rny/ds
      unz = rnz/ds
      !
      ! Compute the Mach number of the normal velocity at this flux point
      !
      rho = face_cv(n,nec,1)
      rhou = face_cv(n,nmx,1)
      rhov = face_cv(n,nmy,1)
      rhow = face_cv(n,nmz,1)
      rhoq2 = rhou*rhou + rhov*rhov + rhow*rhow
      pres = gm1*(face_cv(n,nee,1) - half*rhoq2/rho)
      !
      ! NOTE: The speed of sound (denominator) is gam*pres*rho instead of
      !       gam*pres/rho because the normal velocity (numerator) is rho*vn
      !       instead of just vn.
      !
      mach = (rhou*unx + rhov*uny + rhow*unz) / sqrt( max(gam*pres*rho,rndoff) )
      !
      ! Initialize the BC state to the reference conditions
      !
      rho = rho_ext
      rhou = rho_ext*vx_ext
      rhov = rho_ext*vy_ext
      rhow = rho_ext*vz_ext
      rhoe = rhoe_ext
      !
      ! If mach <= -one, do nothing because in that case we want to use the
      ! the reference pressure for the boundary, so leave rhoe set to the
      ! reference value.
      !
      ! NOTE: The following prescribed BCs are probably wrong if the normal
      !       velocity indicates that this is an outflow instead of the
      !       expected inflow. However, the precribed BCs should be correct
      !       if the normal velocity indicates this is a subsonic inflow BC.
      !
      if (mach > -one) then
        !
        ! Subsonic inflow
        !
        ! If the normal velocity from the interior says this is an outflow
        ! boundary, set the ghost velocity to the magnitude of the reference
        ! velocity multiplied by the unit normal
        ! NOTE: q2 will be the same as q2_ext, and thus ske
        !       will be the same as ske_ext
        !
        if (mach > zero) then
          rhou = rho_ext*sqrt(q2_ext)*unx
          rhov = rho_ext*sqrt(q2_ext)*uny
          rhow = rho_ext*sqrt(q2_ext)*unz
        end if
        !
        ! Set the boundary pressure to the interior pressure if the normal
        ! velocity from the interior is subsonic, otherwise set the boundary
        ! pressure to the reference total pressure
        !
        ! If the normal velocity from the interior is supersonic, set the
        ! boundary pressure to the reference total pressure. Otherwise, leave
        ! the boundary pressure set to the interior value.
        !
        if (mach >= one) then
          pres = p0_ext
        end if
        !
        rhoe = pres*gm1_inv + ske_ext
        !
      end if
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rhou
      face_cv(n,nmy,2) = rhov
      face_cv(n,nmz,2) = rhow
      face_cv(n,nee,2) = rhoe
      !
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
    do n = n1,n2
      !
      ! Unit normal at this flux point
      !
      rnx = facenrm(1,n)
      rny = facenrm(2,n)
      ds = sqrt( rnx*rnx + rny*rny )
      unx = rnx/ds
      uny = rny/ds
      !
      ! Compute the Mach number of the normal velocity at this flux point
      !
      rho = face_cv(n,nec,1)
      rhou = face_cv(n,nmx,1)
      rhov = face_cv(n,nmy,1)
      rhoq2 = rhou*rhou + rhov*rhov
      pres = gm1*(face_cv(n,nee,1) - half*rhoq2/rho)
      !
      ! NOTE: The speed of sound (denominator) is gam*pres*rho instead of
      !       gam*pres/rho because the normal velocity (numerator) is rho*vn
      !       instead of just vn.
      !
      mach = (rhou*unx + rhov*uny) / sqrt( max(gam*pres*rho,rndoff) )
      !
      ! Initialize the BC state to the reference conditions
      !
      rho = rho_ext
      rhou = rho_ext*vx_ext
      rhov = rho_ext*vy_ext
      rhoe = rhoe_ext
      !
      ! If mach <= -one, do nothing because in that case we want to use the
      ! the reference pressure for the boundary, so leave rhoe set to the
      ! reference value.
      !
      ! NOTE: The following prescribed BCs are probably wrong if the normal
      !       velocity indicates that this is an outflow instead of the
      !       expected inflow. However, the precribed BCs should be correct
      !       if the normal velocity indicates this is a subsonic inflow BC.
      !
      if (mach > -one) then
        !
        ! Subsonic inflow
        !
        ! If the normal velocity from the interior says this is an outflow
        ! boundary, set the ghost velocity to the magnitude of the reference
        ! velocity multiplied by the unit normal
        ! NOTE: q2 will be the same as q2_ext, and thus ske
        !       will be the same as ske_ext
        !
        if (mach > zero) then
          rhou = rho_ext*sqrt(q2_ext)*unx
          rhov = rho_ext*sqrt(q2_ext)*uny
        end if
        !
        ! Set the boundary pressure to the interior pressure if the normal
        ! velocity from the interior is subsonic, otherwise set the boundary
        ! pressure to the reference total pressure
        !
        ! If the normal velocity from the interior is supersonic, set the
        ! boundary pressure to the reference total pressure. Otherwise, leave
        ! the boundary pressure set to the interior value.
        !
        if (mach >= one) then
          pres = p0_ext
        end if
        !
        rhoe = pres*gm1_inv + ske_ext
        !
      end if
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rhou
      face_cv(n,nmy,2) = rhov
      face_cv(n,nee,2) = rhoe
      !
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
    do n = n1,n2
      !
      ! Unit normal at this flux point
      !
      unx = facenrm(1,n)/abs(facenrm(1,n))
      !
      ! Compute the Mach number of the normal velocity at this flux point
      !
      rho = face_cv(n,nec,1)
      rhou = face_cv(n,nmx,1)
      pres = gm1*(face_cv(n,nee,1) - half*rhou*rhou/rho)
      !
      ! NOTE: The speed of sound (denominator) is gam*pres*rho instead of
      !       gam*pres/rho because the normal velocity (numerator) is rho*vn
      !       instead of just vn.
      !
      mach = rhou*unx / sqrt( max(gam*pres*rho,rndoff) )
      !
      ! Initialize the BC state to the reference conditions
      !
      rho = rho_ext
      rhou = rho_ext*vx_ext
      rhoe = rhoe_ext
      !
      ! If mach <= -one, do nothing because in that case we want to use the
      ! the reference pressure for the boundary, so leave rhoe set to the
      ! reference value.
      !
      ! NOTE: The following prescribed BCs are probably wrong if the normal
      !       velocity indicates that this is an outflow instead of the
      !       expected inflow. However, the precribed BCs should be correct
      !       if the normal velocity indicates this is a subsonic inflow BC.
      !
      if (mach > -one) then
        !
        ! Subsonic inflow
        !
        ! If the normal velocity from the interior says this is an outflow
        ! boundary, set the ghost velocity to the magnitude of the reference
        ! velocity multiplied by the unit normal
        ! NOTE: q2 will be the same as q2_ext, and thus ske
        !       will be the same as ske_ext
        !
        if (mach > zero) then
          rhou = abs(rhou)*unx
        end if
        !
        ! Set the boundary pressure to the interior pressure if the normal
        ! velocity from the interior is subsonic, otherwise set the boundary
        ! pressure to the reference total pressure
        !
        ! If the normal velocity from the interior is supersonic, set the
        ! boundary pressure to the reference total pressure. Otherwise, leave
        ! the boundary pressure set to the interior value.
        !
        if (mach >= one) then
          pres = p0_ext
        end if
        !
        rhoe = pres*gm1_inv + ske_ext
        !
      end if
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rhou
      face_cv(n,nee,2) = rhoe
      !
    end do
  end if
  !
end subroutine get_bc_generic_inflow_solution_nd
!
!###############################################################################
!
pure subroutine get_bc_generic_outflow_solution(n1,n2,facenrm,face_cv,pv_in)
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
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    call get_bc_generic_outflow_solution_3d(n1,n2,facenrm,face_cv,pv_in)
  else if (nr == 2) then
    call get_bc_generic_outflow_solution_2d(n1,n2,facenrm,face_cv,pv_in)
  else
    call get_bc_generic_outflow_solution_nd(n1,n2,facenrm,face_cv,pv_in)
  end if
  !
end subroutine get_bc_generic_outflow_solution
!
!###############################################################################
!
pure subroutine get_bc_generic_outflow_solution_3d(n1,n2,facenrm,face_cv,pv_in)
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
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,gm1_inv
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: rho,rhou,rhov,rhow,rhoq2,pres,mach
  real(wp) :: pres_ext,q2_ext,a2_ext,p0_ext,total_cond
  !
continue
  !
  ! Precompute gamma-1 and its inverse
  !
  gm1 = gam - one
  gm1_inv = one/gm1
  !
  pres_ext = pv_in(nee)
  a2_ext = max( gam*pres_ext/pv_in(nec) , rndoff )
  q2_ext = pv_in(nmx)**2 + pv_in(nmy)**2 + pv_in(nmz)**2
  !
  total_cond = one + half*gm1*(q2_ext/a2_ext)
  !
  p0_ext = pres_ext * total_cond**(gam*gm1_inv)
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
  do n = n1,n2
    !
    ! Unit normal at this flux point
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    rnz = facenrm(3,n)
    ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
    unx = rnx/ds
    uny = rny/ds
    unz = rnz/ds
    !
    ! Compute the Mach number of the normal velocity at this flux point
    !
    rho = face_cv(n,nec,1)
    rhou = face_cv(n,nmx,1)
    rhov = face_cv(n,nmy,1)
    rhow = face_cv(n,nmz,1)
    rhoq2 = rhou*rhou + rhov*rhov + rhow*rhow
    pres = gm1*(face_cv(n,nee,1) - half*rhoq2/rho)
    !
    ! NOTE: The speed of sound (denominator) is gam*pres*rho instead of
    !       gam*pres/rho because the normal velocity (numerator) is rho*vn
    !       instead of just vn.
    !
    mach = (rhou*unx + rhov*uny + rhow*unz) / sqrt( max(gam*pres*rho,rndoff) )
    !
    ! If mach >= one, do nothing because in that case we want to extrapolate
    ! the interior pressure to the boundary, so leave pres set to the interior
    ! value.
    !
    if (mach < one) then
      !
      ! Subsonic outflow
      !
      ! If the normal velocity from the interior says this is an inflow
      ! boundary, set the ghost velocity to the magnitude of the interior
      ! velocity multiplied by the unit normal
      ! NOTE: rhoq2 will remain the same
      !
      if (mach < zero) then
        rhou = sqrt(rhoq2)*unx
        rhov = sqrt(rhoq2)*uny
        rhow = sqrt(rhoq2)*unz
      end if
      !
      ! Set the boundary pressure to the reference pressure if the normal
      ! velocity from the interior is subsonic, otherwise set the boundary
      ! pressure to the reference total pressure
      !
      pres = pres_ext
      !
      if (mach <= -one) then
        pres = p0_ext
      end if
      !
    end if
    !
    ! Compute the boundary condition state
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rhou
    face_cv(n,nmy,2) = rhov
    face_cv(n,nmz,2) = rhow
    face_cv(n,nee,2) = pres*gm1_inv + half*rhoq2/rho
    !
  end do
  !
end subroutine get_bc_generic_outflow_solution_3d
!
!###############################################################################
!
pure subroutine get_bc_generic_outflow_solution_2d(n1,n2,facenrm,face_cv,pv_in)
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
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,gm1_inv
  real(wp) :: ds,rnx,rny,unx,uny
  real(wp) :: rho,rhou,rhov,rhoq2,pres,mach
  real(wp) :: pres_ext,q2_ext,a2_ext,p0_ext,total_cond
  !
continue
  !
  ! Precompute gamma-1 and its inverse
  !
  gm1 = gam - one
  gm1_inv = one/gm1
  !
  pres_ext = pv_in(nee)
  a2_ext = max( gam*pres_ext/pv_in(nec) , rndoff )
  q2_ext = pv_in(nmx)**2 + pv_in(nmy)**2
  !
  total_cond = one + half*gm1*(q2_ext/a2_ext)
  !
  p0_ext = pres_ext * total_cond**(gam*gm1_inv)
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
  do n = n1,n2
    !
    ! Unit normal at this flux point
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds = sqrt( rnx*rnx + rny*rny )
    unx = rnx/ds
    uny = rny/ds
    !
    ! Compute the Mach number of the normal velocity at this flux point
    !
    rho = face_cv(n,nec,1)
    rhou = face_cv(n,nmx,1)
    rhov = face_cv(n,nmy,1)
    rhoq2 = rhou*rhou + rhov*rhov
    pres = gm1*(face_cv(n,nee,1) - half*rhoq2/rho)
    !
    ! NOTE: The speed of sound (denominator) is gam*pres*rho instead of
    !       gam*pres/rho because the normal velocity (numerator) is rho*vn
    !       instead of just vn.
    !
    mach = (rhou*unx + rhov*uny) / sqrt( max(gam*pres*rho,rndoff) )
    !
    ! If mach >= one, do nothing because in that case we want to extrapolate
    ! the interior pressure to the boundary, so leave pres set to the interior
    ! value.
    !
    if (mach < one) then
      !
      ! Subsonic outflow
      !
      ! If the normal velocity from the interior says this is an inflow
      ! boundary, set the ghost velocity to the magnitude of the interior
      ! velocity multiplied by the unit normal
      ! NOTE: rhoq2 will remain the same
      !
      if (mach < zero) then
        rhou = sqrt(rhoq2)*unx
        rhov = sqrt(rhoq2)*uny
      end if
      !
      ! Set the boundary pressure to the reference pressure if the normal
      ! velocity from the interior is subsonic, otherwise set the boundary
      ! pressure to the reference total pressure
      !
      pres = pres_ext
      !
      if (mach <= -one) then
        pres = p0_ext
      end if
      !
    end if
    !
    ! Compute the boundary condition state
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rhou
    face_cv(n,nmy,2) = rhov
    face_cv(n,nee,2) = pres*gm1_inv + half*rhoq2/rho
    !
  end do
  !
end subroutine get_bc_generic_outflow_solution_2d
!
!###############################################################################
!
pure subroutine get_bc_generic_outflow_solution_nd(n1,n2,facenrm,face_cv,pv_in)
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
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: gm1,gm1_inv
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: rho,rhou,rhov,rhow,rhoq2,pres,mach
  real(wp) :: pres_ext,q2_ext,a2_ext,p0_ext,total_cond
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  ! Precompute gamma-1 and its inverse
  !
  gm1 = gam - one
  gm1_inv = one/gm1
  !
  pres_ext = pv_in(nee)
  a2_ext = max( gam*pres_ext/pv_in(nec) , rndoff )
  !
  if (nr == 3) then
    q2_ext = pv_in(nmx)**2 + pv_in(nmy)**2 + pv_in(nmz)**2
  else if (nr == 2) then
    q2_ext = pv_in(nmx)**2 + pv_in(nmy)**2
  else
    q2_ext = pv_in(nmx)**2
  end if
  !
  total_cond = one + half*gm1*(q2_ext/a2_ext)
  !
  p0_ext = pres_ext * total_cond**(gam*gm1_inv)
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
    do n = n1,n2
      !
      ! Unit normal at this flux point
      !
      rnx = facenrm(1,n)
      rny = facenrm(2,n)
      rnz = facenrm(3,n)
      ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
      unx = rnx/ds
      uny = rny/ds
      unz = rnz/ds
      !
      ! Compute the Mach number of the normal velocity at this flux point
      !
      rho = face_cv(n,nec,1)
      rhou = face_cv(n,nmx,1)
      rhov = face_cv(n,nmy,1)
      rhow = face_cv(n,nmz,1)
      rhoq2 = rhou*rhou + rhov*rhov + rhow*rhow
      pres = gm1*(face_cv(n,nee,1) - half*rhoq2/rho)
      !
      ! NOTE: The speed of sound (denominator) is gam*pres*rho instead of
      !       gam*pres/rho because the normal velocity (numerator) is rho*vn
      !       instead of just vn.
      !
      mach = (rhou*unx + rhov*uny + rhow*unz) / sqrt( max(gam*pres*rho,rndoff) )
      !
      ! If mach >= one, do nothing because in that case we want to extrapolate
      ! the interior pressure to the boundary, so leave pres set to the interior
      ! value.
      !
      if (mach < one) then
        !
        ! Subsonic outflow
        !
        ! If the normal velocity from the interior says this is an inflow
        ! boundary, set the ghost velocity to the magnitude of the interior
        ! velocity multiplied by the unit normal
        ! NOTE: rhoq2 will remain the same
        !
        if (mach < zero) then
          rhou = sqrt(rhoq2)*unx
          rhov = sqrt(rhoq2)*uny
          rhow = sqrt(rhoq2)*unz
        end if
        !
        ! Set the boundary pressure to the reference pressure if the normal
        ! velocity from the interior is subsonic, otherwise set the boundary
        ! pressure to the reference total pressure
        !
        pres = pres_ext
        !
        if (mach <= -one) then
          pres = p0_ext
        end if
        !
      end if
      !
      ! Compute the boundary condition state
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rhou
      face_cv(n,nmy,2) = rhov
      face_cv(n,nmz,2) = rhow
      face_cv(n,nee,2) = pres*gm1_inv + half*rhoq2/rho
      !
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
    do n = n1,n2
      !
      ! Unit normal at this flux point
      !
      rnx = facenrm(1,n)
      rny = facenrm(2,n)
      ds = sqrt( rnx*rnx + rny*rny )
      unx = rnx/ds
      uny = rny/ds
      !
      ! Compute the Mach number of the normal velocity at this flux point
      !
      rho = face_cv(n,nec,1)
      rhou = face_cv(n,nmx,1)
      rhov = face_cv(n,nmy,1)
      rhoq2 = rhou*rhou + rhov*rhov
      pres = gm1*(face_cv(n,nee,1) - half*rhoq2/rho)
      !
      ! NOTE: The speed of sound (denominator) is gam*pres*rho instead of
      !       gam*pres/rho because the normal velocity (numerator) is rho*vn
      !       instead of just vn.
      !
      mach = (rhou*unx + rhov*uny) / sqrt( max(gam*pres*rho,rndoff) )
      !
      ! If mach >= one, do nothing because in that case we want to extrapolate
      ! the interior pressure to the boundary, so leave pres set to the interior
      ! value.
      !
      if (mach < one) then
        !
        ! Subsonic outflow
        !
        ! If the normal velocity from the interior says this is an inflow
        ! boundary, set the ghost velocity to the magnitude of the interior
        ! velocity multiplied by the unit normal
        ! NOTE: rhoq2 will remain the same
        !
        if (mach < zero) then
          rhou = sqrt(rhoq2)*unx
          rhov = sqrt(rhoq2)*uny
        end if
        !
        ! Set the boundary pressure to the reference pressure if the normal
        ! velocity from the interior is subsonic, otherwise set the boundary
        ! pressure to the reference total pressure
        !
        pres = pres_ext
        !
        if (mach <= -one) then
          pres = p0_ext
        end if
        !
      end if
      !
      ! Compute the boundary condition state
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rhou
      face_cv(n,nmy,2) = rhov
      face_cv(n,nee,2) = pres*gm1_inv + half*rhoq2/rho
      !
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
    do n = n1,n2
      !
      ! Unit normal at this flux point
      !
      unx = facenrm(1,n)/abs(facenrm(1,n))
      !
      ! Compute the Mach number of the normal velocity at this flux point
      !
      rho = face_cv(n,nec,1)
      rhou = face_cv(n,nmx,1)
      pres = gm1*(face_cv(n,nee,1) - half*rhou*rhou/rho)
      !
      ! NOTE: The speed of sound (denominator) is gam*pres*rho instead of
      !       gam*pres/rho because the normal velocity (numerator) is rho*vn
      !       instead of just vn.
      !
      mach = rhou*unx / sqrt( max(gam*pres*rho,rndoff) )
      !
      ! If mach >= one, do nothing because in that case we want to extrapolate
      ! the interior pressure to the boundary, so leave pres set to the interior
      ! value.
      !
      if (mach < one) then
        !
        ! Subsonic outflow
        !
        ! If the normal velocity from the interior says this is an inflow
        ! boundary, set the ghost velocity to the magnitude of the interior
        ! velocity multiplied by the unit normal
        ! NOTE: rhoq2 will remain the same
        !
        if (mach < zero) then
          rhou = abs(rhou)*unx
        end if
        !
        ! Set the boundary pressure to the reference pressure if the normal
        ! velocity from the interior is subsonic, otherwise set the boundary
        ! pressure to the reference total pressure
        !
        pres = pres_ext
        !
        if (mach <= -one) then
          pres = p0_ext
        end if
        !
      end if
      !
      ! Compute the boundary condition state
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rhou
      face_cv(n,nee,2) = pres*gm1_inv + half*rhou*rhou/rho
      !
    end do
  end if
  !
end subroutine get_bc_generic_outflow_solution_nd
!
!###############################################################################
!
pure subroutine get_bc_generic_freeflow_solution(n1,n2,facenrm,face_cv,cv_in)
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
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! CV_IN   : Boundary condition flow variables in conservative form
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: cv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    call get_bc_generic_freeflow_solution_3d(n1,n2,facenrm,face_cv,cv_in)
  else if (nr == 2) then
    call get_bc_generic_freeflow_solution_2d(n1,n2,facenrm,face_cv,cv_in)
  else
    call get_bc_generic_freeflow_solution_nd(n1,n2,facenrm,face_cv,cv_in)
  end if
  !
end subroutine get_bc_generic_freeflow_solution
!
!###############################################################################
!
pure subroutine get_bc_generic_freeflow_solution_3d(n1,n2,facenrm,face_cv,cv_in)
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
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! CV_IN   : Boundary condition flow variables in conservative form
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: cv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,gm1_inv
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: rho,rhou,rhov,rhow,rhoe,rhoq2,pres,ske
  real(wp) :: rho_ext,rhou_ext,rhov_ext,rhow_ext,rhoe_ext
  real(wp) :: rhoq2_ext,pres_ext,ske_ext,mach
  !
continue
  !
  gm1 = gam - one
  gm1_inv = one/gm1
  !
  ! Precompute the free stream specific kinetic energy
  ! and the inverse of gamma-1
  !
  rho_ext = cv_in(nec)
  rhou_ext = cv_in(nmx)
  rhov_ext = cv_in(nmy)
  rhow_ext = cv_in(nmz)
  rhoq2_ext = rhou_ext**2 + rhov_ext**2 + rhow_ext**2
  rhoe_ext = cv_in(nee)
  !
  ske_ext = half * rhoq2_ext / rho_ext
  pres_ext = gm1*(rhoe_ext - ske_ext)
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
  do n = n1,n2
    !
    ! Unit normal at this flux point
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    rnz = facenrm(3,n)
    ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
    unx = rnx/ds
    uny = rny/ds
    unz = rnz/ds
    !
    ! Compute the Mach number of the normal velocity at this flux point
    !
    rho = face_cv(n,nec,1)
    rhou = face_cv(n,nmx,1)
    rhov = face_cv(n,nmy,1)
    rhow = face_cv(n,nmz,1)
    rhoe = face_cv(n,nee,1)
    rhoq2 = rhou*rhou + rhov*rhov + rhow*rhow
    ske = half*rhoq2/rho
    pres = gm1*(rhoe - ske)
    !
    ! NOTE: The speed of sound (denominator) is gam*pres*rho instead of
    !       gam*pres/rho because the normal velocity (numerator) is rho*vn
    !       instead of just vn.
    !
    mach = (rhou*unx + rhov*uny + rhow*unz) / sqrt( max(gam*pres*rho,rndoff) )
    !
    ! If mach >= one, do nothing because in that case we want to extrapolate
    ! everything from the interior and we already have those values stored
    ! in rho, rhou, rhov, rhow, and rhoe.
    !
    if (mach <= -one) then
      !
      rho = rho_ext
      rhou = rhou_ext
      rhov = rhov_ext
      rhow = rhow_ext
      rhoe = rhoe_ext
      !
    else if (mach > -one .and. mach <= zero) then
      !
      rho = rho_ext
      rhou = rhou_ext
      rhov = rhov_ext
      rhow = rhow_ext
      rhoe = pres*gm1_inv + ske_ext
      !
    else if (mach < one .and. mach > zero) then
      !
      rhoe = pres_ext*gm1_inv + ske
      !
    end if
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rhou
    face_cv(n,nmy,2) = rhov
    face_cv(n,nmz,2) = rhow
    face_cv(n,nee,2) = rhoe
    !
  end do
  !
end subroutine get_bc_generic_freeflow_solution_3d
!
!###############################################################################
!
pure subroutine get_bc_generic_freeflow_solution_2d(n1,n2,facenrm,face_cv,cv_in)
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
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! CV_IN   : Boundary condition flow variables in conservative form
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: cv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,gm1_inv
  real(wp) :: ds,rnx,rny,unx,uny
  real(wp) :: rho,rhou,rhov,rhoe,rhoq2,pres,ske
  real(wp) :: rho_ext,rhou_ext,rhov_ext,rhoe_ext
  real(wp) :: rhoq2_ext,pres_ext,ske_ext,mach
  !
continue
  !
  gm1 = gam - one
  gm1_inv = one/gm1
  !
  ! Precompute the free stream specific kinetic energy
  ! and the inverse of gamma-1
  !
  rho_ext = cv_in(nec)
  rhou_ext = cv_in(nmx)
  rhov_ext = cv_in(nmy)
  rhoq2_ext = rhou_ext**2 + rhov_ext**2
  rhoe_ext = cv_in(nee)
  !
  ske_ext = half * rhoq2_ext / rho_ext
  pres_ext = gm1*(rhoe_ext - ske_ext)
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
  do n = n1,n2
    !
    ! Unit normal at this flux point
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds = sqrt( rnx*rnx + rny*rny )
    unx = rnx/ds
    uny = rny/ds
    !
    ! Compute the Mach number of the normal velocity at this flux point
    !
    rho = face_cv(n,nec,1)
    rhou = face_cv(n,nmx,1)
    rhov = face_cv(n,nmy,1)
    rhoe = face_cv(n,nee,1)
    rhoq2 = rhou*rhou + rhov*rhov
    ske = half*rhoq2/rho
    pres = gm1*(rhoe - ske)
    !
    ! NOTE: The speed of sound (denominator) is gam*pres*rho instead of
    !       gam*pres/rho because the normal velocity (numerator) is rho*vn
    !       instead of just vn.
    !
    mach = (rhou*unx + rhov*uny) / sqrt( max(gam*pres*rho,rndoff) )
    !
    ! If mach >= one, do nothing because in that case we want to extrapolate
    ! everything from the interior and we already have those values stored
    ! in rho, rhou, rhov, rhow, and rhoe.
    !
    if (mach <= -one) then
      !
      rho = rho_ext
      rhou = rhou_ext
      rhov = rhov_ext
      rhoe = rhoe_ext
      !
    else if (mach > -one .and. mach <= zero) then
      !
      rho = rho_ext
      rhou = rhou_ext
      rhov = rhov_ext
      rhoe = pres*gm1_inv + ske_ext
      !
    else if (mach < one .and. mach > zero) then
      !
      rhoe = pres_ext*gm1_inv + ske
      !
    end if
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rhou
    face_cv(n,nmy,2) = rhov
    face_cv(n,nee,2) = rhoe
    !
  end do
  !
end subroutine get_bc_generic_freeflow_solution_2d
!
!###############################################################################
!
pure subroutine get_bc_generic_freeflow_solution_nd(n1,n2,facenrm,face_cv,cv_in)
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
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! CV_IN   : Boundary condition flow variables in conservative form
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: cv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: gm1,gm1_inv
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: rho,rhou,rhov,rhow,rhoe,rhoq2,pres,ske
  real(wp) :: rho_ext,rhou_ext,rhov_ext,rhow_ext,rhoe_ext
  real(wp) :: rhoq2_ext,pres_ext,ske_ext,mach
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  gm1 = gam - one
  gm1_inv = one/gm1
  !
  ! Precompute the free stream specific kinetic energy
  ! and the inverse of gamma-1
  !
  rho_ext = cv_in(nec)
  if (nr == 3) then
    rhou_ext = cv_in(nmx)
    rhov_ext = cv_in(nmy)
    rhow_ext = cv_in(nmz)
    rhoq2_ext = rhou_ext**2 + rhov_ext**2 + rhow_ext**2
  else if (nr == 2) then
    rhou_ext = cv_in(nmx)
    rhov_ext = cv_in(nmy)
    rhoq2_ext = rhou_ext**2 + rhov_ext**2
  else
    rhou_ext = cv_in(nmx)
    rhoq2_ext = rhou_ext**2
  end if
  rhoe_ext = cv_in(nee)
  !
  ske_ext = half * rhoq2_ext / rho_ext
  pres_ext = gm1*(rhoe_ext - ske_ext)
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
    do n = n1,n2
      !
      ! Unit normal at this flux point
      !
      rnx = facenrm(1,n)
      rny = facenrm(2,n)
      rnz = facenrm(3,n)
      ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
      unx = rnx/ds
      uny = rny/ds
      unz = rnz/ds
      !
      ! Compute the Mach number of the normal velocity at this flux point
      !
      rho = face_cv(n,nec,1)
      rhou = face_cv(n,nmx,1)
      rhov = face_cv(n,nmy,1)
      rhow = face_cv(n,nmz,1)
      rhoe = face_cv(n,nee,1)
      rhoq2 = rhou*rhou + rhov*rhov + rhow*rhow
      ske = half*rhoq2/rho
      pres = gm1*(rhoe - ske)
      !
      ! NOTE: The speed of sound (denominator) is gam*pres*rho instead of
      !       gam*pres/rho because the normal velocity (numerator) is rho*vn
      !       instead of just vn.
      !
      mach = (rhou*unx + rhov*uny + rhow*unz) / sqrt( max(gam*pres*rho,rndoff) )
      !
      ! If mach >= one, do nothing because in that case we want to extrapolate
      ! everything from the interior and we already have those values stored
      ! in rho, rhou, rhov, rhow, and rhoe.
      !
      if (mach <= -one) then
        !
        rho = rho_ext
        rhou = rhou_ext
        rhov = rhov_ext
        rhow = rhow_ext
        rhoe = rhoe_ext
        !
      else if (mach > -one .and. mach <= zero) then
        !
        rho = rho_ext
        rhou = rhou_ext
        rhov = rhov_ext
        rhow = rhow_ext
        rhoe = pres*gm1_inv + ske_ext
        !
      else if (mach < one .and. mach > zero) then
        !
        rhoe = pres_ext*gm1_inv + ske
        !
      end if
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rhou
      face_cv(n,nmy,2) = rhov
      face_cv(n,nmz,2) = rhow
      face_cv(n,nee,2) = rhoe
      !
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
    do n = n1,n2
      !
      ! Unit normal at this flux point
      !
      rnx = facenrm(1,n)
      rny = facenrm(2,n)
      ds = sqrt( rnx*rnx + rny*rny )
      unx = rnx/ds
      uny = rny/ds
      !
      ! Compute the Mach number of the normal velocity at this flux point
      !
      rho = face_cv(n,nec,1)
      rhou = face_cv(n,nmx,1)
      rhov = face_cv(n,nmy,1)
      rhoe = face_cv(n,nee,1)
      rhoq2 = rhou*rhou + rhov*rhov
      ske = half*rhoq2/rho
      pres = gm1*(rhoe - ske)
      !
      ! NOTE: The speed of sound (denominator) is gam*pres*rho instead of
      !       gam*pres/rho because the normal velocity (numerator) is rho*vn
      !       instead of just vn.
      !
      mach = (rhou*unx + rhov*uny) / sqrt( max(gam*pres*rho,rndoff) )
      !
      ! If mach >= one, do nothing because in that case we want to extrapolate
      ! everything from the interior and we already have those values stored
      ! in rho, rhou, rhov, rhow, and rhoe.
      !
      if (mach <= -one) then
        !
        rho = rho_ext
        rhou = rhou_ext
        rhov = rhov_ext
        rhoe = rhoe_ext
        !
      else if (mach > -one .and. mach <= zero) then
        !
        rho = rho_ext
        rhou = rhou_ext
        rhov = rhov_ext
        rhoe = pres*gm1_inv + ske_ext
        !
      else if (mach < one .and. mach > zero) then
        !
        rhoe = pres_ext*gm1_inv + ske
        !
      end if
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rhou
      face_cv(n,nmy,2) = rhov
      face_cv(n,nee,2) = rhoe
      !
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
    do n = n1,n2
      !
      ! Unit normal at this flux point
      !
      unx = facenrm(1,n)/abs(facenrm(1,n))
      !
      ! Compute the Mach number of the normal velocity at this flux point
      !
      rho = face_cv(n,nec,1)
      rhou = face_cv(n,nmx,1)
      rhoe = face_cv(n,nee,1)
      ske = half*rhou*rhou/rho
      pres = gm1*(rhoe - ske)
      !
      ! NOTE: The speed of sound (denominator) is gam*pres*rho instead of
      !       gam*pres/rho because the normal velocity (numerator) is rho*vn
      !       instead of just vn.
      !
      mach = rhou*unx / sqrt( max(gam*pres*rho,rndoff) )
      !
      ! If mach >= one, do nothing because in that case we want to extrapolate
      ! everything from the interior and we already have those values stored
      ! in rho, rhou, rhov, rhow, and rhoe.
      !
      if (mach <= -one) then
        !
        rho = rho_ext
        rhou = rhou_ext
        rhoe = rhoe_ext
        !
      else if (mach > -one .and. mach <= zero) then
        !
        rho = rho_ext
        rhou = rhou_ext
        rhoe = pres*gm1_inv + ske_ext
        !
      else if (mach < one .and. mach > zero) then
        !
        rhoe = pres_ext*gm1_inv + ske
        !
      end if
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rhou
      face_cv(n,nee,2) = rhoe
      !
    end do
  end if
  !
end subroutine get_bc_generic_freeflow_solution_nd
!
!###############################################################################
!
!pure &
subroutine get_bc_sub_inflow_solution(n1,n2,facenrm,face_cv,cv_in,pv_in, &
                                      newton_failures)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Subsonic inflow boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  ! CV_IN   : Boundary condition flow variables in conservative form
  ! NEWTON_FAILURES : Number of failed convergences from Newton iterations
  !
  !.. Use Statements ..
  use ovar,    only : sub_inflow_method
  use gpu_mod, only : tvsp2u_gpu
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: cv_in(:)
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  integer,  intent(inout) :: newton_failures
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  newton_failures = 0
  !
  nr = size(facenrm,dim=1)
  !
  if (sub_inflow_method == 1) then
    !
    if (nr == 3) then
      call get_bc_sub_inflow_newton_3d(n1,n2,facenrm,face_cv, &
                                       pv_in,newton_failures)
    else if (nr == 2) then
      call get_bc_sub_inflow_newton_2d(n1,n2,facenrm,face_cv, &
                                       pv_in,newton_failures)
    else
      call get_bc_sub_inflow_newton_nd(n1,n2,facenrm,face_cv, &
                                       pv_in,newton_failures)
    end if
    !
    ! The Newton subsonic inflow boundary condition returns the BC states
    ! in primitive variables. Convert these to conservative variables.
    !
    call tvsp2u_gpu(n1,n2,face_cv(:,:,2))
    !
  else if (sub_inflow_method == 3) then
    !
    if (nr == 3) then
      call get_bc_sub_inflow_solution_v3_3d(n1,n2,face_cv,pv_in)
    else if (nr == 2) then
      call get_bc_sub_inflow_solution_v3_2d(n1,n2,face_cv,pv_in)
    else
      call get_bc_sub_inflow_solution_v3_nd(n1,n2,face_cv,pv_in)
    end if
    !
  else if (sub_inflow_method == 4) then
    !
    if (nr == 3) then
      call get_bc_sub_inflow_solution_v4_3d(n1,n2,face_cv,pv_in)
    else if (nr == 2) then
      call get_bc_sub_inflow_solution_v4_2d(n1,n2,face_cv,pv_in)
    else
      call get_bc_sub_inflow_solution_v4_nd(n1,n2,face_cv,pv_in)
    end if
    !
  else
    !
    if (nr == 3) then
      call get_bc_sub_inflow_solution_v2_3d(n1,n2,face_cv,cv_in)
    else if (nr == 2) then
      call get_bc_sub_inflow_solution_v2_2d(n1,n2,face_cv,cv_in)
    else
      call get_bc_sub_inflow_solution_v2_nd(n1,n2,face_cv,cv_in)
    end if
    !
  end if
  !
end subroutine get_bc_sub_inflow_solution
!
!###############################################################################
!
!pure &
subroutine get_bc_sub_inflow_newton_3d(n1,n2,facenrm,face_cv, &
                                       pv_in,newton_failures)
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
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  ! NEWTON_FAILURES : Number of failed convergences from Newton iterations
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd
#ifndef DISABLE_FOR_GPU
  use geovar, only : face
#endif
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  integer,  intent(inout) :: newton_failures
  !
  !.. Local Scalars ..
  integer  :: k,n
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: dir_cosx,dir_cosy,dir_cosz
  real(wp) :: rho,rho_inv,vx,vy,vz,q2,pres
  real(wp) :: gm1,gam_over_1mg,tan_aoa,dterr
  real(wp) :: temp_int,aspd_int,veln_int
  real(wp) :: pres_ext,temp_ext,aspd_ext,mach_ext
  real(wp) :: rminus,rafac,veln_ext,func,dfunc,dtemp,dtemp0
  real(wp) :: gam_over_gm1,aspd_ref,mach_ref,total_cond
  real(wp) :: temp_ref,ttot_ref,ptot_ref,aoa
  !
  logical(lk) :: convergence_reached
  !
  !.. Local Parameters ..
  integer, parameter :: newton_iter = 100
  !
#ifndef DISABLE_FOR_GPU
  integer :: i,nf,io,ierr
  real(wp) :: rc,tc
  character(len=100) :: fname
  logical(lk) :: file_is_opened
  real(wp) :: dtemp_save(1:6,1:newton_iter)
  real(wp) :: xyz_fp(1:3)
#endif
  !
continue
  !
#ifndef DISABLE_FOR_GPU
  file_is_opened = fals
#endif
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
  !
  mach_ref = sqrt( pv_in(nmx)*pv_in(nmx) + &
                   pv_in(nmy)*pv_in(nmy) + &
                   pv_in(nmz)*pv_in(nmz) ) / aspd_ref
  !
  total_cond = one + half*gm1*mach_ref*mach_ref
  !
  temp_ref = pv_in(nee)/pv_in(nec)/rgasref_nd
  ttot_ref = temp_ref * total_cond
  ptot_ref = pv_in(nee) * total_cond**gam_over_gm1
  !
  dterr = eps3*ttot_ref
  !
  ! Precompute the direction cosines
  !
  aoa = atan2(pv_in(nmy),pv_in(nmx))
  tan_aoa = tan( aoa )
  !
  dir_cosx = one / sqrt(one + tan_aoa**2)
  dir_cosy = tan_aoa * dir_cosx
  dir_cosz = zero
  !
  newton_failures = 0
  !
  !$ACC PARALLEL LOOP GANG WORKER VECTOR &
  !$ACC&              PRESENT(facenrm,face_cv) &
  !$ACC&              REDUCTION(+:newton_failures)
  do n = n1,n2
    !
    ! Unit normal at this flux point
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    rnz = facenrm(3,n)
    ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
    unx = -rnx/ds
    uny = -rny/ds
    unz = -rnz/ds
    !
    ! Get the host cell primitive variables at this flux point
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    vz = face_cv(n,nmz,1)*rho_inv
    q2 = vx*vx + vy*vy + vz*vz
    pres = gm1*(face_cv(n,nee,1) - half*rho*q2)
    !
    temp_int = pres*rho_inv/rgasref_nd
    !
    ! Normal velocity of the interior/host cell at this flux point
    !
    veln_int = vx*unx + vy*uny + vz*unz
    !
    ! Compute the speed of sound at this flux point
    !
   !aspd_int = sqrt( max( gam*pres*rho_inv , rndoff ) )
    aspd_int = sqrt( max( gam*rgasref_nd*temp_int , rndoff ) )
    !
    rminus = veln_int - (two*aspd_int/gm1)
    !
    temp_ext = temp_ref
    !
    aspd_ext = aspd_ref
    !
    rafac = one / (dir_cosx*unx + dir_cosy*uny + dir_cosz*unz)
    !
    veln_ext = rafac * (rminus + two*aspd_ext/gm1)
    !
    func = ttot_ref - temp_ext - half*gm1*veln_ext*veln_ext/(gam*rgasref_nd)
    dfunc = -(one + veln_ext*rafac/aspd_ext)
    !
    convergence_reached = fals
    newton_iteration: do k = 1,newton_iter
      !
      dtemp = func/dfunc
      if (k == 1) then
        dtemp0 = dtemp
      end if
#ifndef DISABLE_FOR_GPU
      dtemp_save(1,k) = dtemp
      dtemp_save(2,k) = func
      dtemp_save(3,k) = dfunc
      dtemp_save(4,k) = temp_ext
      dtemp_save(5,k) = aspd_ext
      dtemp_save(6,k) = veln_ext
#endif
      !
      if (abs(dtemp) < dterr .or. &
          (abs(dtemp) < ten*dterr .and. abs(dtemp0)/abs(dtemp) < eps1) ) then
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
    face_cv(n,nec,2) = pres_ext / temp_ext / rgasref_nd
    face_cv(n,nmx,2) = veln_ext * dir_cosx
    face_cv(n,nmy,2) = veln_ext * dir_cosy
    face_cv(n,nmz,2) = veln_ext * dir_cosz
    face_cv(n,nee,2) = pres_ext
    !
    ! Check to make sure the temperature converged
    !
    if (.not. convergence_reached) then
      newton_failures = newton_failures + 1
#ifndef DISABLE_FOR_GPU
      dtemp = func/dfunc
      find_face_loop: do k = 1,size(face)
        if (n <= face(k)%end_fp) then
          nf = k
          exit find_face_loop
        end if
      end do find_face_loop
      k = n - face(nf)%beg_fp + 1
      if (.not. file_is_opened) then
        write (fname,5) mypnum
        open (newunit=io,file=trim(adjustl(fname)),status="replace", &
                         action="write",form="formatted",iostat=ierr)
        file_is_opened = true
      end if
      call get_xyz_for_face_point_3d(nf,k,xyz_fp)
      rc = sqrt(xyz_fp(2)*xyz_fp(2) + xyz_fp(3)*xyz_fp(3))
      tc = atan2(xyz_fp(2),xyz_fp(3)) * radians2degrees
      write (io,1) "mypnum   = ",mypnum
      write (io,1) "glb_fp   = ",n
      write (io,1) "nf       = ",nf
      write (io,1) "loc_fp   = ",k
      write (io,2) "X        = ",xyz_fp(1)
      write (io,2) "Y        = ",xyz_fp(2)
      write (io,2) "Z        = ",xyz_fp(3)
      write (io,2) "R        = ",rc
      write (io,2) "Theta    = ",tc,         " "
      write (io,2) "aoa      = ",aoa
      write (io,2) "tan_aoa  = ",tan_aoa
      write (io,2) "dir_cos1 = ",dir_cosx
      write (io,2) "dir_cos2 = ",dir_cosy
      write (io,2) "dir_cos3 = ",dir_cosz
      write (io,2) "rafac    = ",rafac,      " "
      write (io,2) "dens_int = ",rho
      write (io,2) "vx_int   = ",vx
      write (io,2) "vy_int   = ",vy
      write (io,2) "vz_int   = ",vz
      write (io,2) "pres_int = ",pres
      write (io,2) "temp_int = ",temp_int
     !write (io,2) "temp_int = ",pres/rho/rgasref_nd
      write (io,2) "aspd_int = ",aspd_int
      write (io,2) "veln_int = ",veln_int
      write (io,2) "rminus   = ",rminus,     " "
      write (io,2) "dterr    = ",dterr
      write (io,2) "func     = ",func
      write (io,2) "dfunc    = ",dfunc
      write (io,2) "dtemp    = ",dtemp,      " "
      write (io,2) "dens_ext = ",face_cv(n,nec,2)
      write (io,2) "vx_ext   = ",face_cv(n,nmx,2)
      write (io,2) "vy_ext   = ",face_cv(n,nmy,2)
      write (io,2) "vz_ext   = ",face_cv(n,nmz,2)
      write (io,2) "pres_ext = ",pres_ext
      write (io,2) "temp_ext = ",temp_ext
      write (io,2) "aspd_ext = ",aspd_ext
      write (io,2) "veln_ext = ",veln_ext,   " "
      do k = 1,size(dtemp_save,dim=2)
        if (k == 1) write (io,3)
        write (io,4) k,(dtemp_save(i,k), i=1,size(dtemp_save,dim=1))
        if (k == size(dtemp_save,dim=2)) then
          write (io,4) k+1,dtemp,func,dfunc,temp_ext,aspd_ext,veln_ext," "
        end if
      end do
      flush (io)
#endif
    end if
    !
  end do
  !
#ifndef DISABLE_FOR_GPU
  if (file_is_opened) then
    flush (io)
    close (io,iostat=ierr)
  end if
#endif
  !
  ! Format Statements
  !
  1 format (a,i0,:,/,a)
  2 format (a,es14.6,:,/,a)
  3 format (/,101("-"),/, &
              4x,"VALUES AFFECTED BY EACH NEWTON ITERATION",/, &
              101("-"),/, &
              1x,"Iter",7x,"dtemp",11x,"func",12x,"dfunc", &
              9x,"temp_ext",8x,"aspd_ext",8x,"veln_ext",/, &
              101("-"))
  4 format (2x,i3,6(2x,es14.6),:,/,101("-"),a,/)
  5 format ("newton_bc_errors.CPU-",i0,".log")
  !
end subroutine get_bc_sub_inflow_newton_3d
!
!###############################################################################
!
pure &
subroutine get_xyz_for_face_point_3d(nf,kf,xyz_fp)
  !
  !.. Use Statements ..
  use order_mod,         only : maxSP
  use geovar,            only : cell,face,fp,xyz
  use interpolation_mod, only : interp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nf
  integer,     intent(in) :: kf
  real(wp), intent(inout) :: xyz_fp(1:3)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,lf,nc,np,n0
  integer  :: this_geom,this_order
  integer  :: face_geom,face_order
  integer  :: fp_offset,rotation
  real(wp) :: a
  !
  !.. Local Arrays ..
  real(wp) :: txyz(1:maxSP,1:3)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_xyz_for_face_point_3d"
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
    txyz(k,1) = xyz(1,n0+k)
    txyz(k,2) = xyz(2,n0+k)
    txyz(k,3) = xyz(3,n0+k)
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
    xyz_fp(1) = xyz_fp(1) + a * txyz(j,1)
    xyz_fp(2) = xyz_fp(2) + a * txyz(j,2)
    xyz_fp(3) = xyz_fp(3) + a * txyz(j,3)
    !
  end do
  !
end subroutine get_xyz_for_face_point_3d
!
!###############################################################################
!
!pure &
subroutine get_bc_sub_inflow_newton_2d(n1,n2,facenrm,face_cv, &
                                       pv_in,newton_failures)
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
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  ! NEWTON_FAILURES : Number of failed convergences from Newton iterations
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd
#ifndef DISABLE_FOR_GPU
  use geovar, only : face
#endif
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  integer,  intent(inout) :: newton_failures
  !
  !.. Local Scalars ..
  integer  :: k,n
  real(wp) :: ds,rnx,rny,unx,uny
  real(wp) :: dir_cosx,dir_cosy
  real(wp) :: rho,rho_inv,vx,vy,q2,pres
  real(wp) :: gm1,gam_over_1mg,tan_aoa,dterr
  real(wp) :: aspd_int,veln_int
  real(wp) :: pres_ext,temp_ext,aspd_ext,mach_ext
  real(wp) :: rminus,rafac,veln_ext,func,dfunc,dtemp
  real(wp) :: gam_over_gm1,aspd_ref,mach_ref,total_cond
  real(wp) :: temp_ref,ttot_ref,ptot_ref,aoa
  !
  logical(lk) :: convergence_reached
  !
  !.. Local Parameters ..
  integer, parameter :: newton_iter = 100
  !
#ifndef DISABLE_FOR_GPU
  integer :: nf,io,ierr
  character(len=100) :: fname
  logical(lk) :: file_is_opened
  real(wp) :: dtemp_save(1:newton_iter)
#endif
  !
continue
  !
#ifndef DISABLE_FOR_GPU
  file_is_opened = fals
#endif
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
  !
  mach_ref = sqrt( pv_in(nmx)*pv_in(nmx) + pv_in(nmy)*pv_in(nmy) ) / aspd_ref
  !
  total_cond = one + half*gm1*mach_ref*mach_ref
  !
  temp_ref = pv_in(nee)/pv_in(nec)/rgasref_nd
  ttot_ref = temp_ref * total_cond
  ptot_ref = pv_in(nee) * total_cond**gam_over_gm1
  !
  dterr = eps3*ttot_ref
  !
  ! Precompute the direction cosines
  !
  aoa = atan2(pv_in(nmy),pv_in(nmx))
  tan_aoa = tan( aoa )
  !
  dir_cosx = one / sqrt(one + tan_aoa**2)
  dir_cosy = tan_aoa * dir_cosx
  !
  newton_failures = 0
  !
  !$ACC PARALLEL LOOP GANG WORKER VECTOR &
  !$ACC&              PRESENT(facenrm,face_cv) &
  !$ACC&              REDUCTION(+:newton_failures)
  do n = n1,n2
    !
    ! Unit normal at this flux point
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds = sqrt( rnx*rnx + rny*rny )
    unx = -rnx/ds
    uny = -rny/ds
    !
    ! Get the host cell primitive variables at this flux point
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    q2 = vx*vx + vy*vy
    pres = gm1*(face_cv(n,nee,1) - half*rho*q2)
    !
    ! Normal velocity of the interior/host cell at this flux point
    !
    veln_int = vx*unx + vy*uny
    !
    ! Compute the speed of sound at this flux point
    !
    aspd_int = sqrt( max( gam*pres/rho , rndoff ) )
    !
    rminus = veln_int - (two*aspd_int/gm1)
    !
    temp_ext = temp_ref
    !
    aspd_ext = aspd_ref
    !
    rafac = one / (dir_cosx*unx + dir_cosy*uny)
    !
    veln_ext = rafac * (rminus + two*aspd_ext/gm1)
    !
    func = ttot_ref - temp_ext - half*gm1*veln_ext*veln_ext/(gam*rgasref_nd)
    dfunc = -(one + veln_ext*rafac/aspd_ext)
    !
    convergence_reached = fals
    newton_iteration: do k = 1,newton_iter
      !
      dtemp = func/dfunc
#ifndef DISABLE_FOR_GPU
      dtemp_save(k) = dtemp
#endif
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
    face_cv(n,nec,2) = pres_ext / temp_ext / rgasref_nd
    face_cv(n,nmx,2) = veln_ext * dir_cosx
    face_cv(n,nmy,2) = veln_ext * dir_cosy
    face_cv(n,nee,2) = pres_ext
    !
    ! Check to make sure the temperature converged
    !
    if (.not. convergence_reached) then
      newton_failures = newton_failures + 1
#ifndef DISABLE_FOR_GPU
      dtemp = func/dfunc
      find_face_loop: do k = 1,size(face)
        if (n <= face(k)%end_fp) then
          nf = k
          exit find_face_loop
        end if
      end do find_face_loop
      k = n - face(nf)%beg_fp + 1
      if (.not. file_is_opened) then
        write (fname,5) mypnum
        open (newunit=io,file=trim(adjustl(fname)),status="replace", &
                         action="write",form="formatted",iostat=ierr)
        file_is_opened = true
      end if
      write (io,1) "mypnum   = ",mypnum
      write (io,1) "glb_fp   = ",n
      write (io,1) "nf       = ",nf
      write (io,1) "loc_fp   = ",k,          " "
      write (io,2) "aoa      = ",aoa
      write (io,2) "tan_aoa  = ",tan_aoa
      write (io,2) "dir_cos1 = ",dir_cosx
      write (io,2) "dir_cos2 = ",dir_cosy, " "
      write (io,2) "dens_int = ",rho
      write (io,2) "vx_int   = ",vx
      write (io,2) "vy_int   = ",vy
      write (io,2) "pres_int = ",pres
      write (io,2) "temp_int = ",pres/rho/rgasref_nd
      write (io,2) "aspd_int = ",aspd_int
      write (io,2) "veln_int = ",veln_int,   " "
      write (io,2) "dterr    = ",dterr
      write (io,2) "func     = ",func
      write (io,2) "dfunc    = ",dfunc
      write (io,2) "dtemp    = ",dtemp,      " "
      write (io,2) "dens_ext = ",face_cv(n,nec,2)
      write (io,2) "vx_ext   = ",face_cv(n,nmx,2)
      write (io,2) "vy_ext   = ",face_cv(n,nmy,2)
      write (io,2) "pres_ext = ",pres_ext
      write (io,2) "temp_ext = ",temp_ext
      write (io,2) "aspd_ext = ",aspd_ext
      write (io,2) "veln_ext = ",veln_ext,   " "
      do k = 1,size(dtemp_save)
        if (k == 1) write (io,3)
        write (io,4) k,dtemp_save(k)
        if (k == size(dtemp_save)) write (io,4) k+1,dtemp, " "
      end do
      flush (io)
#endif
    end if
    !
  end do
  !
#ifndef DISABLE_FOR_GPU
  if (file_is_opened) then
    flush (io)
    close (io,iostat=ierr)
  end if
#endif
  !
  ! Format Statements
  !
  1 format (a,i0,:,/,a)
  2 format (a,es14.6,:,/,a)
  3 format (/,50("-"),/, &
              4x,"Values of dtemp for each Newton Iteration",/, &
              50("-"))
  4 format (4x,i3,4x,es23.15,:,/,50("-"),a,/)
  5 format ("newton_bc_errors.CPU-",i0,".log")
  !
end subroutine get_bc_sub_inflow_newton_2d
!
!###############################################################################
!
!pure &
subroutine get_bc_sub_inflow_newton_nd(n1,n2,facenrm,face_cv, &
                                       pv_in,newton_failures)
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
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  ! NEWTON_FAILURES : Number of failed convergences from Newton iterations
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd
#ifndef DISABLE_FOR_GPU
  use geovar, only : face
#endif
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  integer,  intent(inout) :: newton_failures
  !
  !.. Local Scalars ..
  integer  :: k,n,nr
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: dir_cosx,dir_cosy,dir_cosz
  real(wp) :: rho,rho_inv,vx,vy,vz,q2,pres
  real(wp) :: gm1,gam_over_1mg,tan_aoa,dterr
  real(wp) :: aspd_int,veln_int
  real(wp) :: pres_ext,temp_ext,aspd_ext,mach_ext
  real(wp) :: rminus,rafac,veln_ext,func,dfunc,dtemp
  real(wp) :: gam_over_gm1,aspd_ref,mach_ref,total_cond
  real(wp) :: temp_ref,ttot_ref,ptot_ref,aoa
  !
  logical(lk) :: convergence_reached
  !
  !.. Local Parameters ..
  integer, parameter :: newton_iter = 100
  !
#ifndef DISABLE_FOR_GPU
  integer :: nf,io,ierr
  character(len=100) :: fname
  logical(lk) :: file_is_opened
  real(wp) :: dtemp_save(1:newton_iter)
#endif
  !
continue
  !
#ifndef DISABLE_FOR_GPU
  file_is_opened = fals
#endif
  !
  nr = size(facenrm,dim=1)
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
  !
  if (nr == 3) then
    mach_ref = sqrt( pv_in(nmx)*pv_in(nmx) + &
                     pv_in(nmy)*pv_in(nmy) + &
                     pv_in(nmz)*pv_in(nmz) ) / aspd_ref
  else if (nr == 2) then
    mach_ref = sqrt( pv_in(nmx)*pv_in(nmx) + &
                     pv_in(nmy)*pv_in(nmy) ) / aspd_ref
  else
    mach_ref = abs( pv_in(nmx) ) / aspd_ref
  end if
  !
  total_cond = one + half*gm1*mach_ref*mach_ref
  !
  temp_ref = pv_in(nee)/pv_in(nec)/rgasref_nd
  ttot_ref = temp_ref * total_cond
  ptot_ref = pv_in(nee) * total_cond**gam_over_gm1
  !
  dterr = eps3*ttot_ref
  !
  ! Precompute the direction cosines
  !
  aoa = atan2(pv_in(nmy),pv_in(nmx))
  tan_aoa = tan( aoa )
  !
  dir_cosx = one / sqrt(one + tan_aoa**2)
  dir_cosy = tan_aoa * dir_cosx
  dir_cosz = zero
  !
  newton_failures = 0
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP GANG WORKER VECTOR &
    !$ACC&              PRESENT(facenrm,face_cv) &
    !$ACC&              REDUCTION(+:newton_failures)
    do n = n1,n2
      !
      ! Unit normal at this flux point
      !
      rnx = facenrm(1,n)
      rny = facenrm(2,n)
      rnz = facenrm(3,n)
      ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
      unx = -rnx/ds
      uny = -rny/ds
      unz = -rnz/ds
      !
      ! Get the host cell primitive variables at this flux point
      !
      rho = face_cv(n,nec,1)
      rho_inv = one/rho
      vx = face_cv(n,nmx,1)*rho_inv
      vy = face_cv(n,nmy,1)*rho_inv
      vz = face_cv(n,nmz,1)*rho_inv
      q2 = vx*vx + vy*vy + vz*vz
      pres = gm1*(face_cv(n,nee,1) - half*rho*q2)
      !
      ! Normal velocity of the interior/host cell at this flux point
      !
      veln_int = vx*unx + vy*uny + vz*unz
      !
      ! Compute the speed of sound at this flux point
      !
      aspd_int = sqrt( max( gam*pres/rho , rndoff ) )
      !
      rminus = veln_int - (two*aspd_int/gm1)
      !
      temp_ext = temp_ref
      !
      aspd_ext = aspd_ref
      !
      rafac = one / (dir_cosx*unx + dir_cosy*uny + dir_cosz*unz)
      !
      veln_ext = rafac * (rminus + two*aspd_ext/gm1)
      !
      func = ttot_ref - temp_ext - half*gm1*veln_ext*veln_ext/(gam*rgasref_nd)
      dfunc = -(one + veln_ext*rafac/aspd_ext)
      !
      convergence_reached = fals
      newton_iteration_3d: do k = 1,newton_iter
        !
        dtemp = func/dfunc
#ifndef DISABLE_FOR_GPU
        dtemp_save(k) = dtemp
#endif
        !
        if (abs(dtemp) < dterr) then
          convergence_reached = true
          exit newton_iteration_3d
        end if
        !
        temp_ext = temp_ext - dtemp
        aspd_ext = sqrt( max( gam*rgasref_nd*temp_ext , rndoff ) )
        veln_ext = rafac * (rminus + two*aspd_ext/gm1)
        !
        func = ttot_ref - temp_ext - half*gm1*veln_ext*veln_ext/(gam*rgasref_nd)
        dfunc = -(one + veln_ext*rafac/aspd_ext)
        !
      end do newton_iteration_3d
      !
      mach_ext = veln_ext / aspd_ext
      !
      pres_ext = ptot_ref*(one + half*gm1*mach_ext*mach_ext)**gam_over_1mg
      !
      face_cv(n,nec,2) = pres_ext / temp_ext / rgasref_nd
      face_cv(n,nmx,2) = veln_ext * dir_cosx
      face_cv(n,nmy,2) = veln_ext * dir_cosy
      face_cv(n,nmz,2) = veln_ext * dir_cosz
      face_cv(n,nee,2) = pres_ext
      !
      ! Check to make sure the temperature converged
      !
      if (.not. convergence_reached) then
        newton_failures = newton_failures + 1
#ifndef DISABLE_FOR_GPU
        dtemp = func/dfunc
        find_face_loop_3d: do k = 1,size(face)
          if (n <= face(k)%end_fp) then
            nf = k
            exit find_face_loop_3d
          end if
        end do find_face_loop_3d
        k = n - face(nf)%beg_fp + 1
        if (.not. file_is_opened) then
          write (fname,5) mypnum
          open (newunit=io,file=trim(adjustl(fname)),status="replace", &
                           action="write",form="formatted",iostat=ierr)
          file_is_opened = true
        end if
        write (io,1) "mypnum   = ",mypnum
        write (io,1) "glb_fp   = ",n
        write (io,1) "nf       = ",nf
        write (io,1) "loc_fp   = ",k,          " "
        write (io,2) "aoa      = ",aoa
        write (io,2) "tan_aoa  = ",tan_aoa
        write (io,2) "dir_cos1 = ",dir_cosx
        write (io,2) "dir_cos2 = ",dir_cosy, " "
        write (io,2) "dens_int = ",rho
        write (io,2) "vx_int   = ",vx
        write (io,2) "vy_int   = ",vy
        write (io,2) "vz_int   = ",vz
        write (io,2) "pres_int = ",pres
        write (io,2) "temp_int = ",pres/rho/rgasref_nd
        write (io,2) "aspd_int = ",aspd_int
        write (io,2) "veln_int = ",veln_int,   " "
        write (io,2) "dterr    = ",dterr
        write (io,2) "func     = ",func
        write (io,2) "dfunc    = ",dfunc
        write (io,2) "dtemp    = ",dtemp,      " "
        write (io,2) "dens_ext = ",face_cv(n,nec,2)
        write (io,2) "vx_ext   = ",face_cv(n,nmx,2)
        write (io,2) "vy_ext   = ",face_cv(n,nmy,2)
        write (io,2) "vz_ext   = ",face_cv(n,nmz,2)
        write (io,2) "pres_ext = ",pres_ext
        write (io,2) "temp_ext = ",temp_ext
        write (io,2) "aspd_ext = ",aspd_ext
        write (io,2) "veln_ext = ",veln_ext,   " "
        do k = 1,size(dtemp_save)
          if (k == 1) write (io,3)
          write (io,4) k,dtemp_save(k)
          if (k == size(dtemp_save)) write (io,4) k+1,dtemp, " "
        end do
        flush (io)
#endif
      end if
      !
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP GANG WORKER VECTOR &
    !$ACC&              PRESENT(facenrm,face_cv) &
    !$ACC&              REDUCTION(+:newton_failures)
    do n = n1,n2
      !
      ! Unit normal at this flux point
      !
      rnx = facenrm(1,n)
      rny = facenrm(2,n)
      ds = sqrt( rnx*rnx + rny*rny )
      unx = -rnx/ds
      uny = -rny/ds
      !
      ! Get the host cell primitive variables at this flux point
      !
      rho = face_cv(n,nec,1)
      rho_inv = one/rho
      vx = face_cv(n,nmx,1)*rho_inv
      vy = face_cv(n,nmy,1)*rho_inv
      q2 = vx*vx + vy*vy
      pres = gm1*(face_cv(n,nee,1) - half*rho*q2)
      !
      ! Normal velocity of the interior/host cell at this flux point
      !
      veln_int = vx*unx + vy*uny
      !
      ! Compute the speed of sound at this flux point
      !
      aspd_int = sqrt( max( gam*pres/rho , rndoff ) )
      !
      rminus = veln_int - (two*aspd_int/gm1)
      !
      temp_ext = temp_ref
      !
      aspd_ext = aspd_ref
      !
      rafac = one / (dir_cosx*unx + dir_cosy*uny)
      !
      veln_ext = rafac * (rminus + two*aspd_ext/gm1)
      !
      func = ttot_ref - temp_ext - half*gm1*veln_ext*veln_ext/(gam*rgasref_nd)
      dfunc = -(one + veln_ext*rafac/aspd_ext)
      !
      convergence_reached = fals
      newton_iteration_2d: do k = 1,newton_iter
        !
        dtemp = func/dfunc
#ifndef DISABLE_FOR_GPU
        dtemp_save(k) = dtemp
#endif
        !
        if (abs(dtemp) < dterr) then
          convergence_reached = true
          exit newton_iteration_2d
        end if
        !
        temp_ext = temp_ext - dtemp
        aspd_ext = sqrt( max( gam*rgasref_nd*temp_ext , rndoff ) )
        veln_ext = rafac * (rminus + two*aspd_ext/gm1)
        !
        func = ttot_ref - temp_ext - half*gm1*veln_ext*veln_ext/(gam*rgasref_nd)
        dfunc = -(one + veln_ext*rafac/aspd_ext)
        !
      end do newton_iteration_2d
      !
      mach_ext = veln_ext / aspd_ext
      !
      pres_ext = ptot_ref*(one + half*gm1*mach_ext*mach_ext)**gam_over_1mg
      !
      face_cv(n,nec,2) = pres_ext / temp_ext / rgasref_nd
      face_cv(n,nmx,2) = veln_ext * dir_cosx
      face_cv(n,nmy,2) = veln_ext * dir_cosy
      face_cv(n,nee,2) = pres_ext
      !
      ! Check to make sure the temperature converged
      !
      if (.not. convergence_reached) then
        newton_failures = newton_failures + 1
#ifndef DISABLE_FOR_GPU
        dtemp = func/dfunc
        find_face_loop_2d: do k = 1,size(face)
          if (n <= face(k)%end_fp) then
            nf = k
            exit find_face_loop_2d
          end if
        end do find_face_loop_2d
        k = n - face(nf)%beg_fp + 1
        if (.not. file_is_opened) then
          write (fname,5) mypnum
          open (newunit=io,file=trim(adjustl(fname)),status="replace", &
                           action="write",form="formatted",iostat=ierr)
          file_is_opened = true
        end if
        write (io,1) "mypnum   = ",mypnum
        write (io,1) "glb_fp   = ",n
        write (io,1) "nf       = ",nf
        write (io,1) "loc_fp   = ",k,          " "
        write (io,2) "aoa      = ",aoa
        write (io,2) "tan_aoa  = ",tan_aoa
        write (io,2) "dir_cos1 = ",dir_cosx
        write (io,2) "dir_cos2 = ",dir_cosy, " "
        write (io,2) "dens_int = ",rho
        write (io,2) "vx_int   = ",vx
        write (io,2) "vy_int   = ",vy
        write (io,2) "pres_int = ",pres
        write (io,2) "temp_int = ",pres/rho/rgasref_nd
        write (io,2) "aspd_int = ",aspd_int
        write (io,2) "veln_int = ",veln_int,   " "
        write (io,2) "dterr    = ",dterr
        write (io,2) "func     = ",func
        write (io,2) "dfunc    = ",dfunc
        write (io,2) "dtemp    = ",dtemp,      " "
        write (io,2) "dens_ext = ",face_cv(n,nec,2)
        write (io,2) "vx_ext   = ",face_cv(n,nmx,2)
        write (io,2) "vy_ext   = ",face_cv(n,nmy,2)
        write (io,2) "pres_ext = ",pres_ext
        write (io,2) "temp_ext = ",temp_ext
        write (io,2) "aspd_ext = ",aspd_ext
        write (io,2) "veln_ext = ",veln_ext,   " "
        do k = 1,size(dtemp_save)
          if (k == 1) write (io,3)
          write (io,4) k,dtemp_save(k)
          if (k == size(dtemp_save)) write (io,4) k+1,dtemp, " "
        end do
        flush (io)
#endif
      end if
      !
    end do
  else
    !$ACC PARALLEL LOOP GANG WORKER VECTOR &
    !$ACC&              PRESENT(facenrm,face_cv) &
    !$ACC&              REDUCTION(+:newton_failures)
    do n = n1,n2
      !
      ! Unit normal at this flux point
      !
      rnx = facenrm(1,n)
      ds = sqrt( rnx*rnx )
      unx = -rnx/ds
      !
      ! Get the host cell primitive variables at this flux point
      !
      rho = face_cv(n,nec,1)
      rho_inv = one/rho
      vx = face_cv(n,nmx,1)*rho_inv
      q2 = vx*vx
      pres = gm1*(face_cv(n,nee,1) - half*rho*q2)
      !
      ! Normal velocity of the interior/host cell at this flux point
      !
      veln_int = vx*unx
      !
      ! Compute the speed of sound at this flux point
      !
      aspd_int = sqrt( max( gam*pres/rho , rndoff ) )
      !
      rminus = veln_int - (two*aspd_int/gm1)
      !
      temp_ext = temp_ref
      !
      aspd_ext = aspd_ref
      !
      rafac = one / (dir_cosx*unx)
      !
      veln_ext = rafac * (rminus + two*aspd_ext/gm1)
      !
      func = ttot_ref - temp_ext - half*gm1*veln_ext*veln_ext/(gam*rgasref_nd)
      dfunc = -(one + veln_ext*rafac/aspd_ext)
      !
      convergence_reached = fals
      newton_iteration_1d: do k = 1,newton_iter
        !
        dtemp = func/dfunc
#ifndef DISABLE_FOR_GPU
        dtemp_save(k) = dtemp
#endif
        !
        if (abs(dtemp) < dterr) then
          convergence_reached = true
          exit newton_iteration_1d
        end if
        !
        temp_ext = temp_ext - dtemp
        aspd_ext = sqrt( max( gam*rgasref_nd*temp_ext , rndoff ) )
        veln_ext = rafac * (rminus + two*aspd_ext/gm1)
        !
        func = ttot_ref - temp_ext - half*gm1*veln_ext*veln_ext/(gam*rgasref_nd)
        dfunc = -(one + veln_ext*rafac/aspd_ext)
        !
      end do newton_iteration_1d
      !
      mach_ext = veln_ext / aspd_ext
      !
      pres_ext = ptot_ref*(one + half*gm1*mach_ext*mach_ext)**gam_over_1mg
      !
      face_cv(n,nec,2) = pres_ext / temp_ext / rgasref_nd
      face_cv(n,nmx,2) = veln_ext * dir_cosx
      face_cv(n,nee,2) = pres_ext
      !
      ! Check to make sure the temperature converged
      !
      if (.not. convergence_reached) then
        newton_failures = newton_failures + 1
#ifndef DISABLE_FOR_GPU
        dtemp = func/dfunc
        find_face_loop_1d: do k = 1,size(face)
          if (n <= face(k)%end_fp) then
            nf = k
            exit find_face_loop_1d
          end if
        end do find_face_loop_1d
        k = n - face(nf)%beg_fp + 1
        if (.not. file_is_opened) then
          write (fname,5) mypnum
          open (newunit=io,file=trim(adjustl(fname)),status="replace", &
                           action="write",form="formatted",iostat=ierr)
          file_is_opened = true
        end if
        write (io,1) "mypnum   = ",mypnum
        write (io,1) "glb_fp   = ",n
        write (io,1) "nf       = ",nf
        write (io,1) "loc_fp   = ",k,          " "
        write (io,2) "aoa      = ",aoa
        write (io,2) "tan_aoa  = ",tan_aoa
        write (io,2) "dir_cos1 = ",dir_cosx
        write (io,2) "dir_cos2 = ",dir_cosy, " "
        write (io,2) "dens_int = ",rho
        write (io,2) "vx_int   = ",vx
        write (io,2) "pres_int = ",pres
        write (io,2) "temp_int = ",pres/rho/rgasref_nd
        write (io,2) "aspd_int = ",aspd_int
        write (io,2) "veln_int = ",veln_int,   " "
        write (io,2) "dterr    = ",dterr
        write (io,2) "func     = ",func
        write (io,2) "dfunc    = ",dfunc
        write (io,2) "dtemp    = ",dtemp,      " "
        write (io,2) "dens_ext = ",face_cv(n,nec,2)
        write (io,2) "vx_ext   = ",face_cv(n,nmx,2)
        write (io,2) "pres_ext = ",pres_ext
        write (io,2) "temp_ext = ",temp_ext
        write (io,2) "aspd_ext = ",aspd_ext
        write (io,2) "veln_ext = ",veln_ext,   " "
        do k = 1,size(dtemp_save)
          if (k == 1) write (io,3)
          write (io,4) k,dtemp_save(k)
          if (k == size(dtemp_save)) write (io,4) k+1,dtemp, " "
        end do
        flush (io)
#endif
      end if
      !
    end do
  end if
  !
#ifndef DISABLE_FOR_GPU
  if (file_is_opened) then
    flush (io)
    close (io,iostat=ierr)
  end if
#endif
  !
  ! Format Statements
  !
  1 format (a,i0,:,/,a)
  2 format (a,es14.6,:,/,a)
  3 format (/,50("-"),/, &
              4x,"Values of dtemp for each Newton Iteration",/, &
              50("-"))
  4 format (4x,i3,4x,es23.15,:,/,50("-"),a,/)
  5 format ("newton_bc_errors.CPU-",i0,".log")
  !
end subroutine get_bc_sub_inflow_newton_nd
!
!###############################################################################
!
pure subroutine get_bc_sub_inflow_solution_v2_3d(n1,n2,face_cv,cv_in)
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
  ! FACE_CV : Conservative variables at the flux points
  ! CV_IN   : Boundary condition flow variables in conservative form
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ske_ext,ske_int
  real(wp) :: rho,rhou,rhov,rhow,rhoe
  !
continue
  !
  ! Precompute the free stream specific kinetic energy
  !
  ske_ext = half * ( cv_in(nmx)*cv_in(nmx) + &
                     cv_in(nmy)*cv_in(nmy) + &
                     cv_in(nmz)*cv_in(nmz) ) / cv_in(nec)
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,cv_in)
  do n = n1,n2
    !
    rho = face_cv(n,nec,1)
    rhou = face_cv(n,nmx,1)
    rhov = face_cv(n,nmy,1)
    rhow = face_cv(n,nmz,1)
    rhoe = face_cv(n,nee,1)
    !
    ske_int = half * (rhou*rhou + rhov*rhov + rhow*rhow) / rho
    !
    face_cv(n,nec,2) = cv_in(nec)
    face_cv(n,nmx,2) = cv_in(nmx)
    face_cv(n,nmy,2) = cv_in(nmy)
    face_cv(n,nmz,2) = cv_in(nmz)
    face_cv(n,nee,2) = rhoe - ske_int + ske_ext
    !
  end do
  !
end subroutine get_bc_sub_inflow_solution_v2_3d
!
!###############################################################################
!
pure subroutine get_bc_sub_inflow_solution_v2_2d(n1,n2,face_cv,cv_in)
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
  ! FACE_CV : Conservative variables at the flux points
  ! CV_IN   : Boundary condition flow variables in conservative form
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ske_ext,ske_int
  real(wp) :: rho,rhou,rhov,rhoe
  !
continue
  !
  ! Precompute the free stream specific kinetic energy
  !
  ske_ext = half * ( cv_in(nmx)*cv_in(nmx) + &
                     cv_in(nmy)*cv_in(nmy) ) / cv_in(nec)
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,cv_in)
  do n = n1,n2
    !
    rho = face_cv(n,nec,1)
    rhou = face_cv(n,nmx,1)
    rhov = face_cv(n,nmy,1)
    rhoe = face_cv(n,nee,1)
    !
    ske_int = half * (rhou*rhou + rhov*rhov) / rho
    !
    face_cv(n,nec,2) = cv_in(nec)
    face_cv(n,nmx,2) = cv_in(nmx)
    face_cv(n,nmy,2) = cv_in(nmy)
    face_cv(n,nee,2) = rhoe - ske_int + ske_ext
    !
  end do
  !
end subroutine get_bc_sub_inflow_solution_v2_2d
!
!###############################################################################
!
pure subroutine get_bc_sub_inflow_solution_v2_nd(n1,n2,face_cv,cv_in)
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
  ! FACE_CV : Conservative variables at the flux points
  ! CV_IN   : Boundary condition flow variables in conservative form
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: ske_ext,ske_int
  real(wp) :: rho,rhou,rhov,rhow,rhoe
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    !
    ! Precompute the free stream specific kinetic energy
    !
    ske_ext = half * ( cv_in(nmx)*cv_in(nmx) + &
                       cv_in(nmy)*cv_in(nmy) + &
                       cv_in(nmz)*cv_in(nmz) ) / cv_in(nec)
    !
    !$ACC PARALLEL LOOP PRESENT(face_cv,cv_in)
    do n = n1,n2
      !
      rho = face_cv(n,nec,1)
      rhou = face_cv(n,nmx,1)
      rhov = face_cv(n,nmy,1)
      rhow = face_cv(n,nmz,1)
      rhoe = face_cv(n,nee,1)
      !
      ske_int = half * (rhou*rhou + rhov*rhov + rhow*rhow) / rho
      !
      face_cv(n,nec,2) = cv_in(nec)
      face_cv(n,nmx,2) = cv_in(nmx)
      face_cv(n,nmy,2) = cv_in(nmy)
      face_cv(n,nmz,2) = cv_in(nmz)
      face_cv(n,nee,2) = rhoe - ske_int + ske_ext
      !
    end do
  else if (nr == 2) then
    !
    ! Precompute the free stream specific kinetic energy
    !
    ske_ext = half * ( cv_in(nmx)*cv_in(nmx) + &
                       cv_in(nmy)*cv_in(nmy) ) / cv_in(nec)
    !
    !$ACC PARALLEL LOOP PRESENT(face_cv,cv_in)
    do n = n1,n2
      !
      rho = face_cv(n,nec,1)
      rhou = face_cv(n,nmx,1)
      rhov = face_cv(n,nmy,1)
      rhoe = face_cv(n,nee,1)
      !
      ske_int = half * (rhou*rhou + rhov*rhov) / rho
      !
      face_cv(n,nec,2) = cv_in(nec)
      face_cv(n,nmx,2) = cv_in(nmx)
      face_cv(n,nmy,2) = cv_in(nmy)
      face_cv(n,nee,2) = rhoe - ske_int + ske_ext
      !
    end do
  else
    !
    ! Precompute the free stream specific kinetic energy
    !
    ske_ext = half * cv_in(nmx)*cv_in(nmx) / cv_in(nec)
    !
    !$ACC PARALLEL LOOP PRESENT(face_cv,cv_in)
    do n = n1,n2
      !
      rho = face_cv(n,nec,1)
      rhou = face_cv(n,nmx,1)
      rhoe = face_cv(n,nee,1)
      !
      ske_int = half * rhou * rhou / rho
      !
      face_cv(n,nec,2) = cv_in(nec)
      face_cv(n,nmx,2) = cv_in(nmx)
      face_cv(n,nee,2) = rhoe - ske_int + ske_ext
      !
    end do
  end if
  !
end subroutine get_bc_sub_inflow_solution_v2_nd
!
!###############################################################################
!
pure subroutine get_bc_sub_inflow_solution_v3_3d(n1,n2,face_cv,pv_in)
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
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,gm1_inv,gam_over_gm1
  real(wp) :: rho,vx,a,b,c
  real(wp) :: p0_ext,t0_ext,rho0_ext
  real(wp) :: aspd,mach,total_cond
  !
continue
  !
  ! Precompute various BC quantities
  !
  gm1 = gam - one
  gm1_inv = one / gm1
  gam_over_gm1 = gam*gm1_inv
  !
  aspd = sqrt(gam*pv_in(nee)/pv_in(nec))
  !
  mach = sqrt( pv_in(nmx)*pv_in(nmx) + &
               pv_in(nmy)*pv_in(nmy) + &
               pv_in(nmz)*pv_in(nmz) ) / aspd
  !
  total_cond = one + half*gm1*mach*mach
  !
  t0_ext = total_cond * (pv_in(nee)/pv_in(nec)/rgasref_nd)
  p0_ext = pv_in(nee) * total_cond**gam_over_gm1
  rho0_ext = p0_ext/(rgasref_nd*t0_ext)
  !
  a = gm1_inv*rgasref_nd*t0_ext
  b = two*gam*a
  c = half*gm1/gam
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv)
  do n = n1,n2
    !
    vx = face_cv(n,nmx,1)/face_cv(n,nec,1)
    !
    rho = rho0_ext * ( (one + vx*vx/(b-vx*vx))**(-gm1_inv) )
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rho*vx
    face_cv(n,nmy,2) = zero
    face_cv(n,nmz,2) = zero
    face_cv(n,nee,2) = rho * (a + c*vx*vx)
    !
  end do
  !
end subroutine get_bc_sub_inflow_solution_v3_3d
!
!###############################################################################
!
pure subroutine get_bc_sub_inflow_solution_v3_2d(n1,n2,face_cv,pv_in)
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
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,gm1_inv,gam_over_gm1
  real(wp) :: rho,vx,a,b,c
  real(wp) :: p0_ext,t0_ext,rho0_ext
  real(wp) :: aspd,mach,total_cond
  !
continue
  !
  ! Precompute various BC quantities
  !
  gm1 = gam - one
  gm1_inv = one / gm1
  gam_over_gm1 = gam*gm1_inv
  !
  aspd = sqrt(gam*pv_in(nee)/pv_in(nec))
  !
  mach = sqrt( pv_in(nmx)*pv_in(nmx) + pv_in(nmy)*pv_in(nmy) ) / aspd
  !
  total_cond = one + half*gm1*mach*mach
  !
  t0_ext = total_cond * (pv_in(nee)/pv_in(nec)/rgasref_nd)
  p0_ext = pv_in(nee) * total_cond**gam_over_gm1
  rho0_ext = p0_ext/(rgasref_nd*t0_ext)
  !
  a = gm1_inv*rgasref_nd*t0_ext
  b = two*gam*a
  c = half*gm1/gam
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv)
  do n = n1,n2
    !
    vx = face_cv(n,nmx,1)/face_cv(n,nec,1)
    !
    rho = rho0_ext * ( (one + vx*vx/(b-vx*vx))**(-gm1_inv) )
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rho*vx
    face_cv(n,nmy,2) = zero
    face_cv(n,nee,2) = rho * (a + c*vx*vx)
    !
  end do
  !
end subroutine get_bc_sub_inflow_solution_v3_2d
!
!###############################################################################
!
pure subroutine get_bc_sub_inflow_solution_v3_nd(n1,n2,face_cv,pv_in)
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
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: gm1,gm1_inv,gam_over_gm1
  real(wp) :: rho,vx,a,b,c
  real(wp) :: p0_ext,t0_ext,rho0_ext
  real(wp) :: aspd,mach,total_cond
  !
continue
  !
  nr = nme-nmb+1
  !
  ! Precompute various BC quantities
  !
  gm1 = gam - one
  gm1_inv = one / gm1
  gam_over_gm1 = gam*gm1_inv
  !
  aspd = sqrt(gam*pv_in(nee)/pv_in(nec))
  !
  if (nr == 3) then
    mach = sqrt( pv_in(nmx)*pv_in(nmx) + &
                 pv_in(nmy)*pv_in(nmy) + &
                 pv_in(nmz)*pv_in(nmz) ) / aspd
  else if (nr == 2) then
    mach = sqrt( pv_in(nmx)*pv_in(nmx) + &
                 pv_in(nmy)*pv_in(nmy) ) / aspd
  else
    mach = abs( pv_in(nmx) ) / aspd
  end if
  !
  total_cond = one + half*gm1*mach*mach
  !
  t0_ext = total_cond * (pv_in(nee)/pv_in(nec)/rgasref_nd)
  p0_ext = pv_in(nee) * total_cond**gam_over_gm1
  rho0_ext = p0_ext/(rgasref_nd*t0_ext)
  !
  a = gm1_inv*rgasref_nd*t0_ext
  b = two*gam*a
  c = half*gm1/gam
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(face_cv)
    do n = n1,n2
      !
      vx = face_cv(n,nmx,1)/face_cv(n,nec,1)
      !
      rho = rho0_ext * ( (one + vx*vx/(b-vx*vx))**(-gm1_inv) )
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*vx
      face_cv(n,nmy,2) = zero
      face_cv(n,nmz,2) = zero
      face_cv(n,nee,2) = rho * (a + c*vx*vx)
      !
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(face_cv)
    do n = n1,n2
      !
      vx = face_cv(n,nmx,1)/face_cv(n,nec,1)
      !
      rho = rho0_ext * ( (one + vx*vx/(b-vx*vx))**(-gm1_inv) )
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*vx
      face_cv(n,nmy,2) = zero
      face_cv(n,nee,2) = rho * (a + c*vx*vx)
      !
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(face_cv)
    do n = n1,n2
      !
      vx = face_cv(n,nmx,1)/face_cv(n,nec,1)
      !
      rho = rho0_ext * ( (one + vx*vx/(b-vx*vx))**(-gm1_inv) )
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*vx
      face_cv(n,nee,2) = rho * (a + c*vx*vx)
      !
    end do
  end if
  !
end subroutine get_bc_sub_inflow_solution_v3_nd
!
!###############################################################################
!
pure subroutine get_bc_sub_inflow_solution_v4_3d(n1,n2,face_cv,pv_in)
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
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: rho,rhou,rhov,rhow,rhoq2,pres
  real(wp) :: rho_ext,q2_ext,rhoe_ext,rmsq,asq
  real(wp) :: gm1,gm1_inv,gm1og,aspd,mach,p0_ext
  !
continue
  !
  ! Precompute various BC quantities
  !
  gm1 = gam - one
  gm1_inv = one / gm1
  gm1og = gm1 / gam
  !
  aspd = sqrt(gam*pv_in(nee)/pv_in(nec))
  !
  mach = sqrt( pv_in(nmx)*pv_in(nmx) + &
               pv_in(nmy)*pv_in(nmy) + &
               pv_in(nmz)*pv_in(nmz) ) / aspd
  !
  p0_ext = pv_in(nee) * (one + half*gm1*mach*mach)**(gam*gm1_inv)
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv)
  do n = n1,n2
    !
    rho = face_cv(n,nec,1)
    rhou = face_cv(n,nmx,1)
    rhov = face_cv(n,nmy,1)
    rhow = face_cv(n,nmz,1)
    !
    rhoq2 = rhou*rhou + rhov*rhov + rhow*rhow
    !
    pres = gm1*(face_cv(n,nee,1) - half*rhoq2/rho)
    !
    rmsq = ( (p0_ext/pres)**gm1og - one ) * two * gm1_inv
    !
    asq = one / (one + half*gm1*rmsq)
    !
    rho_ext = gam*pres/asq
    q2_ext = rmsq * asq
    rhoe_ext = pres*gm1_inv + half*rho_ext*q2_ext
    !
    face_cv(n,nec,2) = rho_ext
    face_cv(n,nmx,2) = rho_ext*sqrt(q2_ext)
    face_cv(n,nmy,2) = zero
    face_cv(n,nmz,2) = zero
    face_cv(n,nee,2) = rhoe_ext
    !
  end do
  !
end subroutine get_bc_sub_inflow_solution_v4_3d
!
!###############################################################################
!
pure subroutine get_bc_sub_inflow_solution_v4_2d(n1,n2,face_cv,pv_in)
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
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: rho,rhou,rhov,rhoq2,pres
  real(wp) :: rho_ext,q2_ext,rhoe_ext,rmsq,asq
  real(wp) :: gm1,gm1_inv,gm1og,aspd,mach,p0_ext
  !
continue
  !
  ! Precompute various BC quantities
  !
  gm1 = gam - one
  gm1_inv = one / gm1
  gm1og = gm1 / gam
  !
  aspd = sqrt(gam*pv_in(nee)/pv_in(nec))
  !
  mach = sqrt( pv_in(nmx)*pv_in(nmx) + pv_in(nmy)*pv_in(nmy) ) / aspd
  !
  p0_ext = pv_in(nee) * (one + half*gm1*mach*mach)**(gam*gm1_inv)
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv)
  do n = n1,n2
    !
    rho = face_cv(n,nec,1)
    rhou = face_cv(n,nmx,1)
    rhov = face_cv(n,nmy,1)
    !
    rhoq2 = rhou*rhou + rhov*rhov
    !
    pres = gm1*(face_cv(n,nee,1) - half*rhoq2/rho)
    !
    rmsq = ( (p0_ext/pres)**gm1og - one ) * two * gm1_inv
    !
    asq = one / (one + half*gm1*rmsq)
    !
    rho_ext = gam*pres/asq
    q2_ext = rmsq * asq
    rhoe_ext = pres*gm1_inv + half*rho_ext*q2_ext
    !
    face_cv(n,nec,2) = rho_ext
    face_cv(n,nmx,2) = rho_ext*sqrt(q2_ext)
    face_cv(n,nmy,2) = zero
    face_cv(n,nee,2) = rhoe_ext
    !
  end do
  !
end subroutine get_bc_sub_inflow_solution_v4_2d
!
!###############################################################################
!
pure subroutine get_bc_sub_inflow_solution_v4_nd(n1,n2,face_cv,pv_in)
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
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: rho,rhou,rhov,rhow,rhoq2,pres
  real(wp) :: rho_ext,q2_ext,rhoe_ext,rmsq,asq
  real(wp) :: gm1,gm1_inv,gm1og,aspd,mach,p0_ext
  !
continue
  !
  nr = nme-nmb+1
  !
  ! Precompute various BC quantities
  !
  gm1 = gam - one
  gm1_inv = one / gm1
  gm1og = gm1 / gam
  !
  aspd = sqrt(gam*pv_in(nee)/pv_in(nec))
  !
  if (nr == 3) then
    mach = sqrt( pv_in(nmx)*pv_in(nmx) + &
                 pv_in(nmy)*pv_in(nmy) + &
                 pv_in(nmz)*pv_in(nmz) ) / aspd
  else if (nr == 2) then
    mach = sqrt( pv_in(nmx)*pv_in(nmx) + &
                 pv_in(nmy)*pv_in(nmy) ) / aspd
  else
    mach = abs( pv_in(nmx) ) / aspd
  end if
  !
  p0_ext = pv_in(nee) * (one + half*gm1*mach*mach)**(gam*gm1_inv)
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(face_cv)
    do n = n1,n2
      !
      rho = face_cv(n,nec,1)
      rhou = face_cv(n,nmx,1)
      rhov = face_cv(n,nmy,1)
      rhow = face_cv(n,nmz,1)
      !
      rhoq2 = rhou*rhou + rhov*rhov + rhow*rhow
      !
      pres = gm1*(face_cv(n,nee,1) - half*rhoq2/rho)
      !
      rmsq = ( (p0_ext/pres)**gm1og - one ) * two * gm1_inv
      !
      asq = one / (one + half*gm1*rmsq)
      !
      rho_ext = gam*pres/asq
      q2_ext = rmsq * asq
      rhoe_ext = pres*gm1_inv + half*rho_ext*q2_ext
      !
      face_cv(n,nec,2) = rho_ext
      face_cv(n,nmx,2) = rho_ext*sqrt(q2_ext)
      face_cv(n,nmy,2) = zero
      face_cv(n,nmz,2) = zero
      face_cv(n,nee,2) = rhoe_ext
      !
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(face_cv)
    do n = n1,n2
      !
      rho = face_cv(n,nec,1)
      rhou = face_cv(n,nmx,1)
      rhov = face_cv(n,nmy,1)
      !
      rhoq2 = rhou*rhou + rhov*rhov
      !
      pres = gm1*(face_cv(n,nee,1) - half*rhoq2/rho)
      !
      rmsq = ( (p0_ext/pres)**gm1og - one ) * two * gm1_inv
      !
      asq = one / (one + half*gm1*rmsq)
      !
      rho_ext = gam*pres/asq
      q2_ext = rmsq * asq
      rhoe_ext = pres*gm1_inv + half*rho_ext*q2_ext
      !
      face_cv(n,nec,2) = rho_ext
      face_cv(n,nmx,2) = rho_ext*sqrt(q2_ext)
      face_cv(n,nmy,2) = zero
      face_cv(n,nee,2) = rhoe_ext
      !
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(face_cv)
    do n = n1,n2
      !
      rho = face_cv(n,nec,1)
      rhou = face_cv(n,nmx,1)
      !
      pres = gm1*(face_cv(n,nee,1) - half*rhou*rhou/rho)
      !
      rmsq = ( (p0_ext/pres)**gm1og - one ) * two * gm1_inv
      !
      asq = one / (one + half*gm1*rmsq)
      !
      rho_ext = gam*pres/asq
      q2_ext = rmsq * asq
      rhoe_ext = pres*gm1_inv + half*rho_ext*q2_ext
      !
      face_cv(n,nec,2) = rho_ext
      face_cv(n,nmx,2) = rho_ext*sqrt(q2_ext)
      face_cv(n,nee,2) = rhoe_ext
      !
    end do
  end if
  !
end subroutine get_bc_sub_inflow_solution_v4_nd
!
!###############################################################################
!
pure subroutine get_bc_sub_outflow_solution(n1,n2,facenrm,face_cv,pv_in)
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
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
  !.. Local Parameters ..
  logical(lk), parameter :: use_old_sub_outflow = fals
 !logical(lk), parameter :: use_old_sub_outflow = true
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  if (use_old_sub_outflow) then
    if (nr == 3) then
      call get_bc_sub_outflow_solution_old_3d(n1,n2,face_cv,pv_in)
    else if (nr == 2) then
      call get_bc_sub_outflow_solution_old_2d(n1,n2,face_cv,pv_in)
    else
      call get_bc_sub_outflow_solution_old_nd(n1,n2,face_cv,pv_in)
    end if
  else
    if (nr == 3) then
      call get_bc_sub_outflow_solution_3d(n1,n2,facenrm,face_cv,pv_in)
    else if (nr == 2) then
      call get_bc_sub_outflow_solution_2d(n1,n2,facenrm,face_cv,pv_in)
    else
      call get_bc_sub_outflow_solution_nd(n1,n2,facenrm,face_cv,pv_in)
    end if
  end if
  !
end subroutine get_bc_sub_outflow_solution
!
!###############################################################################
!
pure subroutine get_bc_sub_outflow_solution_old_3d(n1,n2,face_cv,pv_in)
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
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: gm1_inv,rho,rhou,rhov,rhow,rhoq2
  !
continue
  !
  nr = nme-nmb+1
  !
  ! Precompute the inverse of gamma-1
  !
  gm1_inv = one / (gam - one)
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,pv_in)
  do n = n1,n2
    !
    rho = face_cv(n,nec,1)
    rhou = face_cv(n,nmx,1)
    rhov = face_cv(n,nmy,1)
    rhow = face_cv(n,nmz,1)
    !
    rhoq2 = rhou*rhou + rhov*rhov + rhow*rhow
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rhou
    face_cv(n,nmy,2) = rhov
    face_cv(n,nmz,2) = rhow
    face_cv(n,nee,2) = pv_in(nee)*gm1_inv + half*rhoq2/rho
    !
  end do
  !
end subroutine get_bc_sub_outflow_solution_old_3d
!
!###############################################################################
!
pure subroutine get_bc_sub_outflow_solution_old_2d(n1,n2,face_cv,pv_in)
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
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: gm1_inv,rho,rhou,rhov,rhoq2
  !
continue
  !
  nr = nme-nmb+1
  !
  ! Precompute the inverse of gamma-1
  !
  gm1_inv = one / (gam - one)
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,pv_in)
  do n = n1,n2
    !
    rho = face_cv(n,nec,1)
    rhou = face_cv(n,nmx,1)
    rhov = face_cv(n,nmy,1)
    !
    rhoq2 = rhou*rhou + rhov*rhov
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rhou
    face_cv(n,nmy,2) = rhov
    face_cv(n,nee,2) = pv_in(nee)*gm1_inv + half*rhoq2/rho
    !
  end do
  !
end subroutine get_bc_sub_outflow_solution_old_2d
!
!###############################################################################
!
pure subroutine get_bc_sub_outflow_solution_old_nd(n1,n2,face_cv,pv_in)
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
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: gm1_inv,rho,rhou,rhov,rhow,rhoq2
  !
continue
  !
  nr = nme-nmb+1
  !
  ! Precompute the inverse of gamma-1
  !
  gm1_inv = one / (gam - one)
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(face_cv,pv_in)
    do n = n1,n2
      !
      rho = face_cv(n,nec,1)
      rhou = face_cv(n,nmx,1)
      rhov = face_cv(n,nmy,1)
      rhow = face_cv(n,nmz,1)
      !
      rhoq2 = rhou*rhou + rhov*rhov + rhow*rhow
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rhou
      face_cv(n,nmy,2) = rhov
      face_cv(n,nmz,2) = rhow
      face_cv(n,nee,2) = pv_in(nee)*gm1_inv + half*rhoq2/rho
      !
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(face_cv,pv_in)
    do n = n1,n2
      !
      rho = face_cv(n,nec,1)
      rhou = face_cv(n,nmx,1)
      rhov = face_cv(n,nmy,1)
      !
      rhoq2 = rhou*rhou + rhov*rhov
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rhou
      face_cv(n,nmy,2) = rhov
      face_cv(n,nee,2) = pv_in(nee)*gm1_inv + half*rhoq2/rho
      !
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(face_cv,pv_in)
    do n = n1,n2
      !
      rho = face_cv(n,nec,1)
      rhou = face_cv(n,nmx,1)
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rhou
      face_cv(n,nee,2) = pv_in(nee)*gm1_inv + half*rhou*rhou/rho
      !
    end do
  end if
  !
end subroutine get_bc_sub_outflow_solution_old_nd
!
!###############################################################################
!
pure subroutine get_bc_sub_outflow_solution_3d(n1,n2,facenrm,face_cv,pv_in)
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
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  !
  !.. Use Statements ..
  use ovar, only : gam,ptotref_nd
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,gm1_inv,ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: rho,rho_inv,vx,vy,vz,vn,q2,pres,pres_ext,aspd
  !
continue
  !
  ! Precompute the inverse of gamma-1
  !
  gm1 = gam - one
  gm1_inv = one/gm1
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,pv_in)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    rnz = facenrm(3,n)
    ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
    unx = rnx/ds
    uny = rny/ds
    unz = rnz/ds
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    vz = face_cv(n,nmz,1)*rho_inv
    !
    vn = vx*unx + vy*uny + vz*unz
    q2 = vx*vx + vy*vy + vz*vz
    !
    pres = gm1*(face_cv(n,nee,1) - half*rho*q2)
    aspd = sqrt(max(gam*pres*rho_inv,rndoff))
    !
    ! If the normal velocity from the interior says this is an inflow
    ! boundary, set the ghost velocity to the magnitude of the interior
    ! velocity multiplied by the unit normal
    ! NOTE: q2 will remain the same
    !
    if (vn < zero) then
      vx = sqrt(q2)*unx
      vy = sqrt(q2)*uny
      vz = sqrt(q2)*unz
    end if
    !
    ! Set the boundary pressure to the reference pressure if the normal
    ! velocity from the interior is subsonic, otherwise set the boundary
    ! pressure to the reference total pressure
    !
    pres_ext = pv_in(nee)
    !
    if (abs(vn)/aspd >= one) then
      pres_ext = ptotref_nd
    end if
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rho*vx
    face_cv(n,nmy,2) = rho*vy
    face_cv(n,nmz,2) = rho*vz
    face_cv(n,nee,2) = pres_ext*gm1_inv + half*rho*q2
    !
  end do
  !
end subroutine get_bc_sub_outflow_solution_3d
!
!###############################################################################
!
pure subroutine get_bc_sub_outflow_solution_2d(n1,n2,facenrm,face_cv,pv_in)
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
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  !
  !.. Use Statements ..
  use ovar, only : gam,ptotref_nd
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,gm1_inv,ds,rnx,rny,unx,uny
  real(wp) :: rho,rho_inv,vx,vy,vn,q2,pres,pres_ext,aspd
  !
continue
  !
  ! Precompute the inverse of gamma-1
  !
  gm1 = gam - one
  gm1_inv = one/gm1
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,pv_in)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds = sqrt( rnx*rnx + rny*rny )
    unx = rnx/ds
    uny = rny/ds
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    !
    vn = vx*unx + vy*uny
    q2 = vx*vx + vy*vy
    !
    pres = gm1*(face_cv(n,nee,1) - half*rho*q2)
    aspd = sqrt(max(gam*pres*rho_inv,rndoff))
    !
    ! If the normal velocity from the interior says this is an inflow
    ! boundary, set the ghost velocity to the magnitude of the interior
    ! velocity multiplied by the unit normal
    ! NOTE: q2 will remain the same
    !
    if (vn < zero) then
      vx = sqrt(q2)*unx
      vy = sqrt(q2)*uny
    end if
    !
    ! Set the boundary pressure to the reference pressure if the normal
    ! velocity from the interior is subsonic, otherwise set the boundary
    ! pressure to the reference total pressure
    !
    pres_ext = pv_in(nee)
    !
    if (abs(vn)/aspd >= one) then
      pres_ext = ptotref_nd
    end if
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rho*vx
    face_cv(n,nmy,2) = rho*vy
    face_cv(n,nee,2) = pres_ext*gm1_inv + half*rho*q2
    !
  end do
  !
end subroutine get_bc_sub_outflow_solution_2d
!
!###############################################################################
!
pure subroutine get_bc_sub_outflow_solution_nd(n1,n2,facenrm,face_cv,pv_in)
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
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  ! PV_IN   : Boundary condition flow variables in primitive form
  !
  !.. Use Statements ..
  use ovar, only : gam,ptotref_nd
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_in(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: gm1,gm1_inv,ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: rho,rho_inv,vx,vy,vz,vn,q2,pres,pres_ext,aspd
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  ! Precompute the inverse of gamma-1
  !
  gm1 = gam - one
  gm1_inv = one/gm1
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,pv_in)
    do n = n1,n2
      !
      rnx = facenrm(1,n)
      rny = facenrm(2,n)
      rnz = facenrm(3,n)
      ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
      unx = rnx/ds
      uny = rny/ds
      unz = rnz/ds
      !
      rho = face_cv(n,nec,1)
      rho_inv = one/rho
      vx = face_cv(n,nmx,1)*rho_inv
      vy = face_cv(n,nmy,1)*rho_inv
      vz = face_cv(n,nmz,1)*rho_inv
      !
      vn = vx*unx + vy*uny + vz*unz
      q2 = vx*vx + vy*vy + vz*vz
      !
      pres = gm1*(face_cv(n,nee,1) - half*rho*q2)
      aspd = sqrt(max(gam*pres*rho_inv,rndoff))
      !
      ! If the normal velocity from the interior says this is an inflow
      ! boundary, set the ghost velocity to the magnitude of the interior
      ! velocity multiplied by the unit normal
      ! NOTE: q2 will remain the same
      !
      if (vn < zero) then
        vx = sqrt(q2)*unx
        vy = sqrt(q2)*uny
        vz = sqrt(q2)*unz
      end if
      !
      ! Set the boundary pressure to the reference pressure if the normal
      ! velocity from the interior is subsonic, otherwise set the boundary
      ! pressure to the reference total pressure
      !
      pres_ext = pv_in(nee)
      !
      if (abs(vn)/aspd >= one) then
        pres_ext = ptotref_nd
      end if
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*vx
      face_cv(n,nmy,2) = rho*vy
      face_cv(n,nmz,2) = rho*vz
      face_cv(n,nee,2) = pres_ext*gm1_inv + half*rho*q2
      !
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,pv_in)
    do n = n1,n2
      !
      rnx = facenrm(1,n)
      rny = facenrm(2,n)
      ds = sqrt( rnx*rnx + rny*rny )
      unx = rnx/ds
      uny = rny/ds
      !
      rho = face_cv(n,nec,1)
      rho_inv = one/rho
      vx = face_cv(n,nmx,1)*rho_inv
      vy = face_cv(n,nmy,1)*rho_inv
      !
      vn = vx*unx + vy*uny
      q2 = vx*vx + vy*vy
      !
      pres = gm1*(face_cv(n,nee,1) - half*rho*q2)
      aspd = sqrt(max(gam*pres*rho_inv,rndoff))
      !
      ! If the normal velocity from the interior says this is an inflow
      ! boundary, set the ghost velocity to the magnitude of the interior
      ! velocity multiplied by the unit normal
      ! NOTE: q2 will remain the same
      !
      if (vn < zero) then
        vx = sqrt(q2)*unx
        vy = sqrt(q2)*uny
      end if
      !
      ! Set the boundary pressure to the reference pressure if the normal
      ! velocity from the interior is subsonic, otherwise set the boundary
      ! pressure to the reference total pressure
      !
      pres_ext = pv_in(nee)
      !
      if (abs(vn)/aspd >= one) then
        pres_ext = ptotref_nd
      end if
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*vx
      face_cv(n,nmy,2) = rho*vy
      face_cv(n,nee,2) = pres_ext*gm1_inv + half*rho*q2
      !
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,pv_in)
    do n = n1,n2
      !
      unx = facenrm(1,n)/abs(facenrm(1,n))
      !
      rho = face_cv(n,nec,1)
      rho_inv = one/rho
      vx = face_cv(n,nmx,1)*rho_inv
      !
      vn = vx*unx
      q2 = vx*vx
      !
      pres = gm1*(face_cv(n,nee,1) - half*rho*q2)
      aspd = sqrt(max(gam*pres*rho_inv,rndoff))
      !
      ! If the normal velocity from the interior says this is an inflow
      ! boundary, set the ghost velocity to the magnitude of the interior
      ! velocity multiplied by the unit normal
      ! NOTE: q2 will remain the same
      !
      if (vn < zero) then
        vx = sqrt(q2)*unx
      end if
      !
      ! Set the boundary pressure to the reference pressure if the normal
      ! velocity from the interior is subsonic, otherwise set the boundary
      ! pressure to the reference total pressure
      !
      pres_ext = pv_in(nee)
      !
      if (abs(vn)/aspd >= one) then
        pres_ext = ptotref_nd
      end if
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*vx
      face_cv(n,nee,2) = pres_ext*gm1_inv + half*rho*q2
      !
    end do
  end if
  !
end subroutine get_bc_sub_outflow_solution_nd
!
!###############################################################################
!
pure subroutine get_bc_sup_inflow_solution(n1,n2,face_cv,cv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Supersonic inflow boundary condition
  ! ##    - Fix all variables to the inflow condition
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV : Conservative variables at the ghost cell flux points
  ! CV_IN   : Boundary condition flow variables in conservative form
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv_in(:)
  real(wp), intent(inout) :: face_cv(:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    call get_bc_sup_inflow_solution_3d(n1,n2,face_cv,cv_in)
  else if (nr == 2) then
    call get_bc_sup_inflow_solution_2d(n1,n2,face_cv,cv_in)
  else
    call get_bc_sup_inflow_solution_nd(n1,n2,face_cv,cv_in)
  end if
  !
end subroutine get_bc_sup_inflow_solution
!
!###############################################################################
!
pure subroutine get_bc_sup_inflow_solution_3d(n1,n2,face_cv,cv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Supersonic inflow boundary condition
  ! ##    - Fix all variables to the inflow condition
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV : Conservative variables at the ghost cell flux points
  ! CV_IN   : Boundary condition flow variables in conservative form
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv_in(:)
  real(wp), intent(inout) :: face_cv(:,:)
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  ! Set all variables to the inflow condition
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,cv_in)
  do n = n1,n2
    face_cv(n,nec) = cv_in(nec)
    face_cv(n,nmx) = cv_in(nmx)
    face_cv(n,nmy) = cv_in(nmy)
    face_cv(n,nmz) = cv_in(nmz)
    face_cv(n,nee) = cv_in(nee)
  end do
  !
end subroutine get_bc_sup_inflow_solution_3d
!
!###############################################################################
!
pure subroutine get_bc_sup_inflow_solution_2d(n1,n2,face_cv,cv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Supersonic inflow boundary condition
  ! ##    - Fix all variables to the inflow condition
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV : Conservative variables at the ghost cell flux points
  ! CV_IN   : Boundary condition flow variables in conservative form
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv_in(:)
  real(wp), intent(inout) :: face_cv(:,:)
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  ! Set all variables to the inflow condition
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,cv_in)
  do n = n1,n2
    face_cv(n,nec) = cv_in(nec)
    face_cv(n,nmx) = cv_in(nmx)
    face_cv(n,nmy) = cv_in(nmy)
    face_cv(n,nee) = cv_in(nee)
  end do
  !
end subroutine get_bc_sup_inflow_solution_2d
!
!###############################################################################
!
pure subroutine get_bc_sup_inflow_solution_nd(n1,n2,face_cv,cv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Supersonic inflow boundary condition
  ! ##    - Fix all variables to the inflow condition
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV : Conservative variables at the ghost cell flux points
  ! CV_IN   : Boundary condition flow variables in conservative form
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv_in(:)
  real(wp), intent(inout) :: face_cv(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n
  !
continue
  !
  ! Set all variables to the inflow condition
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) &
  !$ACC&              PRESENT(face_cv,cv_in)
  do m = 1,size(cv_in)
    do n = n1,n2
      face_cv(n,m) = cv_in(m)
    end do
  end do
  !
end subroutine get_bc_sup_inflow_solution_nd
!
!###############################################################################
!
pure subroutine get_bc_sup_outflow_solution(n1,n2,face_cv)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Supersonic outflow boundary condition
  ! ##    - Extrapolate all variables
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV : Conservative variables at the flux points
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    call get_bc_sup_outflow_solution_3d(n1,n2,face_cv)
  else if (nr == 2) then
    call get_bc_sup_outflow_solution_2d(n1,n2,face_cv)
  else
    call get_bc_sup_outflow_solution_nd(n1,n2,face_cv)
  end if
  !
end subroutine get_bc_sup_outflow_solution
!
!###############################################################################
!
pure subroutine get_bc_sup_outflow_solution_3d(n1,n2,face_cv)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Supersonic outflow boundary condition
  ! ##    - Extrapolate all variables
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV : Conservative variables at the flux points
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv)
  do n = n1,n2
    face_cv(n,nec,2) = face_cv(n,nec,1)
    face_cv(n,nmx,2) = face_cv(n,nmx,1)
    face_cv(n,nmy,2) = face_cv(n,nmy,1)
    face_cv(n,nmz,2) = face_cv(n,nmz,1)
    face_cv(n,nee,2) = face_cv(n,nee,1)
  end do
  !
end subroutine get_bc_sup_outflow_solution_3d
!
!###############################################################################
!
pure subroutine get_bc_sup_outflow_solution_2d(n1,n2,face_cv)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Supersonic outflow boundary condition
  ! ##    - Extrapolate all variables
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV : Conservative variables at the flux points
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv)
  do n = n1,n2
    face_cv(n,nec,2) = face_cv(n,nec,1)
    face_cv(n,nmx,2) = face_cv(n,nmx,1)
    face_cv(n,nmy,2) = face_cv(n,nmy,1)
    face_cv(n,nee,2) = face_cv(n,nee,1)
  end do
  !
end subroutine get_bc_sup_outflow_solution_2d
!
!###############################################################################
!
pure subroutine get_bc_sup_outflow_solution_nd(n1,n2,face_cv)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Supersonic outflow boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV : Conservative variables at the flux points
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(face_cv)
    do n = n1,n2
      face_cv(n,nec,2) = face_cv(n,nec,1)
      face_cv(n,nmx,2) = face_cv(n,nmx,1)
      face_cv(n,nmy,2) = face_cv(n,nmy,1)
      face_cv(n,nmz,2) = face_cv(n,nmz,1)
      face_cv(n,nee,2) = face_cv(n,nee,1)
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(face_cv)
    do n = n1,n2
      face_cv(n,nec,2) = face_cv(n,nec,1)
      face_cv(n,nmx,2) = face_cv(n,nmx,1)
      face_cv(n,nmy,2) = face_cv(n,nmy,1)
      face_cv(n,nee,2) = face_cv(n,nee,1)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(face_cv)
    do n = n1,n2
      face_cv(n,nec,2) = face_cv(n,nec,1)
      face_cv(n,nmx,2) = face_cv(n,nmx,1)
      face_cv(n,nee,2) = face_cv(n,nee,1)
    end do
  end if
  !
end subroutine get_bc_sup_outflow_solution_nd
!
!###############################################################################
!
pure subroutine get_bc_slip_wall_solution(n1,n2,facenrm,face_cv)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Slip wall / symmetry boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  if (nr == 3) then
    call get_bc_slip_wall_solution_3d(n1,n2,facenrm,face_cv)
  else if (nr == 2) then
    call get_bc_slip_wall_solution_2d(n1,n2,facenrm,face_cv)
  else
    call get_bc_slip_wall_solution_nd(n1,n2,facenrm,face_cv)
  end if
  !
end subroutine get_bc_slip_wall_solution
!
!###############################################################################
!
pure subroutine get_bc_slip_wall_solution_3d(n1,n2,facenrm,face_cv)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Slip wall / symmetry boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  !
  !.. Use Statements ..
  use ovar, only : gam,walls_are_exact
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: gm1,gm1_inv
  real(wp) :: rho,rho_inv,vx,vy,vz,vn,rhoe
  real(wp) :: velocity_const,energy_const
  !
continue
  !
  velocity_const = merge(one,two,walls_are_exact)
  energy_const = merge(-one,zero,walls_are_exact)
  !
  gm1 = gam - one
  gm1_inv = one/gm1
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    rnz = facenrm(3,n)
    ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
    unx = rnx/ds
    uny = rny/ds
    unz = rnz/ds
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    vz = face_cv(n,nmz,1)*rho_inv
    rhoe = face_cv(n,nee,1)
    !
    vn = vx*unx + vy*uny + vz*unz
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rho*(vx - velocity_const*vn*unx)
    face_cv(n,nmy,2) = rho*(vy - velocity_const*vn*uny)
    face_cv(n,nmz,2) = rho*(vz - velocity_const*vn*unz)
    face_cv(n,nee,2) = rhoe + energy_const*half*rho*vn*vn
    !
  end do
  !
end subroutine get_bc_slip_wall_solution_3d
!
!###############################################################################
!
pure subroutine get_bc_slip_wall_solution_2d(n1,n2,facenrm,face_cv)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Slip wall / symmetry boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  !
  !.. Use Statements ..
  use ovar, only : gam,walls_are_exact
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,rny,unx,uny
  real(wp) :: gm1,gm1_inv
  real(wp) :: rho,rho_inv,vx,vy,vn,rhoe
  real(wp) :: velocity_const,energy_const
  !
continue
  !
  velocity_const = merge(one,two,walls_are_exact)
  energy_const = merge(-one,zero,walls_are_exact)
  !
  gm1 = gam - one
  gm1_inv = one/gm1
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds = sqrt( rnx*rnx + rny*rny )
    unx = rnx/ds
    uny = rny/ds
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    rhoe = face_cv(n,nee,1)
    !
    vn = vx*unx + vy*uny
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rho*(vx - velocity_const*vn*unx)
    face_cv(n,nmy,2) = rho*(vy - velocity_const*vn*uny)
    face_cv(n,nee,2) = rhoe + energy_const*half*rho*vn*vn
    !
  end do
  !
end subroutine get_bc_slip_wall_solution_2d
!
!###############################################################################
!
pure subroutine get_bc_slip_wall_solution_nd(n1,n2,facenrm,face_cv)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Slip wall / symmetry boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conservative variables at the flux points
  !
  !.. Use Statements ..
  use ovar, only : gam,walls_are_exact
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: gm1,gm1_inv
  real(wp) :: rho,rho_inv,vx,vy,vz,vn,rhoe
  real(wp) :: velocity_const,energy_const
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  velocity_const = merge(one,two,walls_are_exact)
  energy_const = merge(-one,zero,walls_are_exact)
  !
  gm1 = gam - one
  gm1_inv = one/gm1
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
    do n = n1,n2
      !
      rnx = facenrm(1,n)
      rny = facenrm(2,n)
      rnz = facenrm(3,n)
      ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
      unx = rnx/ds
      uny = rny/ds
      unz = rnz/ds
      !
      rho = face_cv(n,nec,1)
      rho_inv = one/rho
      vx = face_cv(n,nmx,1)*rho_inv
      vy = face_cv(n,nmy,1)*rho_inv
      vz = face_cv(n,nmz,1)*rho_inv
      rhoe = face_cv(n,nee,1)
      !
      vn = vx*unx + vy*uny + vz*unz
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*(vx - velocity_const*vn*unx)
      face_cv(n,nmy,2) = rho*(vy - velocity_const*vn*uny)
      face_cv(n,nmz,2) = rho*(vz - velocity_const*vn*unz)
      face_cv(n,nee,2) = rhoe + energy_const*half*rho*vn*vn
      !
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
    do n = n1,n2
      !
      rnx = facenrm(1,n)
      rny = facenrm(2,n)
      ds = sqrt( rnx*rnx + rny*rny )
      unx = rnx/ds
      uny = rny/ds
      !
      rho = face_cv(n,nec,1)
      rho_inv = one/rho
      vx = face_cv(n,nmx,1)*rho_inv
      vy = face_cv(n,nmy,1)*rho_inv
      rhoe = face_cv(n,nee,1)
      !
      vn = vx*unx + vy*uny
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*(vx - velocity_const*vn*unx)
      face_cv(n,nmy,2) = rho*(vy - velocity_const*vn*uny)
      face_cv(n,nee,2) = rhoe + energy_const*half*rho*vn*vn
      !
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(face_cv)
    do n = n1,n2
      !
      rho = face_cv(n,nec,1)
      vx = face_cv(n,nmx,1)/rho
      rhoe = face_cv(n,nee,1)
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rho*(one - velocity_const)*vx
      face_cv(n,nee,2) = rhoe + energy_const*half*rho*vx*vx
      !
    end do
  end if
  !
end subroutine get_bc_slip_wall_solution_nd
!
!###############################################################################
!
pure subroutine get_bc_adiabatic_wall_solution(n1,n2,face_cv,pv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip adiabatic wall
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV : Conservative variables at the flux points
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: pv_in(:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    call get_bc_adiabatic_wall_solution_3d(n1,n2,face_cv,pv_in)
  else if (nr == 2) then
    call get_bc_adiabatic_wall_solution_2d(n1,n2,face_cv,pv_in)
  else
    call get_bc_adiabatic_wall_solution_nd(n1,n2,face_cv,pv_in)
  end if
  !
  ! The above subroutines will return the boundary condition states as
  ! primitive variables. These need to be converted to conservative variables.
  !
  call tvsp2u_gpu(n1,n2,face_cv(:,:,2))
  !
end subroutine get_bc_adiabatic_wall_solution
!
!###############################################################################
!
pure subroutine get_bc_adiabatic_wall_solution_3d(n1,n2,face_cv,pv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip adiabatic wall
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV : Conservative variables at the flux points
  !
  !.. Use Statements ..
  use ovar, only : gam,walls_are_exact
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: pv_in(:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: rho,rho_inv,vx,vy,vz,q2,pres
  real(wp) :: gm1,wall_vx,wall_vy,wall_vz
  real(wp) :: velocity_const
  real(wp) :: moving_wall_const
  !
continue
  !
  velocity_const = merge(zero,-one,walls_are_exact)
  moving_wall_const = merge(one,two,walls_are_exact)
  !
  wall_vx = pv_in(nmx)
  wall_vy = pv_in(nmy)
  wall_vz = pv_in(nmz)
  !
  gm1 = gam - one
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv)
  do n = n1,n2
    !
    rho  = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    vz = face_cv(n,nmz,1)*rho_inv
    !
    q2 = vx*vx + vy*vy + vz*vz
    !
    pres = gm1*(face_cv(n,nee,1) - half*rho*q2)
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = vx*velocity_const + wall_vx*moving_wall_const
    face_cv(n,nmy,2) = vy*velocity_const + wall_vy*moving_wall_const
    face_cv(n,nmz,2) = vz*velocity_const + wall_vz*moving_wall_const
    face_cv(n,nee,2) = pres
    !
  end do
  !
end subroutine get_bc_adiabatic_wall_solution_3d
!
!###############################################################################
!
pure subroutine get_bc_adiabatic_wall_solution_2d(n1,n2,face_cv,pv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip adiabatic wall
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV : Conservative variables at the flux points
  !
  !.. Use Statements ..
  use ovar, only : gam,walls_are_exact
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: pv_in(:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: rho,rho_inv,vx,vy,q2,pres
  real(wp) :: gm1,wall_vx,wall_vy
  real(wp) :: velocity_const
  real(wp) :: moving_wall_const
  !
continue
  !
  velocity_const = merge(zero,-one,walls_are_exact)
  moving_wall_const = merge(one,two,walls_are_exact)
  !
  wall_vx = pv_in(nmx)
  wall_vy = pv_in(nmy)
  !
  gm1 = gam - one
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv)
  do n = n1,n2
    !
    rho  = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    !
    q2 = vx*vx + vy*vy
    !
    pres = gm1*(face_cv(n,nee,1) - half*rho*q2)
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = vx*velocity_const + wall_vx*moving_wall_const
    face_cv(n,nmy,2) = vy*velocity_const + wall_vy*moving_wall_const
    face_cv(n,nee,2) = pres
    !
  end do
  !
end subroutine get_bc_adiabatic_wall_solution_2d
!
!###############################################################################
!
pure subroutine get_bc_adiabatic_wall_solution_nd(n1,n2,face_cv,pv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip adiabatic wall
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV : Conservative variables at the flux points
  !
  !.. Use Statements ..
  use ovar, only : gam,walls_are_exact
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: pv_in(:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: rho,rho_inv,vx,vy,vz,q2,pres
  real(wp) :: gm1,wall_vx,wall_vy,wall_vz
  real(wp) :: velocity_const
  real(wp) :: moving_wall_const
  !
continue
  !
  nr = nme-nmb+1
  velocity_const = merge(zero,-one,walls_are_exact)
  moving_wall_const = merge(one,two,walls_are_exact)
  !
  gm1 = gam - one
  !
  if (nr == 3) then
    !
    wall_vx = pv_in(nmx)
    wall_vy = pv_in(nmy)
    wall_vz = pv_in(nmz)
    !
    !$ACC PARALLEL LOOP PRESENT(face_cv)
    do n = n1,n2
      !
      rho = face_cv(n,nec,1)
      rho_inv = one/rho
      vx = face_cv(n,nmx,1)*rho_inv
      vy = face_cv(n,nmy,1)*rho_inv
      vz = face_cv(n,nmz,1)*rho_inv
      !
      q2 = vx*vx + vy*vy + vz*vz
      !
      pres = gm1*(face_cv(n,nee,1) - half*rho*q2)
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = vx*velocity_const + wall_vx*moving_wall_const
      face_cv(n,nmy,2) = vy*velocity_const + wall_vy*moving_wall_const
      face_cv(n,nmz,2) = vz*velocity_const + wall_vz*moving_wall_const
      face_cv(n,nee,2) = pres
      !
    end do
    !
  else if (nr == 2) then
    !
    wall_vx = pv_in(nmx)
    wall_vy = pv_in(nmy)
    !
    !$ACC PARALLEL LOOP PRESENT(face_cv)
    do n = n1,n2
      !
      rho = face_cv(n,nec,1)
      rho_inv = one/rho
      vx = face_cv(n,nmx,1)*rho_inv
      vy = face_cv(n,nmy,1)*rho_inv
      !
      q2 = vx*vx + vy*vy
      !
      pres = gm1*(face_cv(n,nee,1) - half*rho*q2)
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = vx*velocity_const + wall_vx*moving_wall_const
      face_cv(n,nmy,2) = vy*velocity_const + wall_vy*moving_wall_const
      face_cv(n,nee,2) = pres
      !
    end do
    !
  else
    !
    wall_vx = pv_in(nmx)
    !
    !$ACC PARALLEL LOOP PRESENT(face_cv)
    do n = n1,n2
      !
      rho = face_cv(n,nec,1)
      rho_inv = one/rho
      vx = face_cv(n,nmx,1)*rho_inv
      !
      pres = gm1*(face_cv(n,nee,1) - half*rho*vx*vx)
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = vx*velocity_const + wall_vx*moving_wall_const
      face_cv(n,nee,2) = pres
      !
    end do
    !
  end if
  !
end subroutine get_bc_adiabatic_wall_solution_nd
!
!###############################################################################
!
pure subroutine get_bc_isothermal_wall_solution(n1,n2,face_cv, &
                                                wall_temperature,pv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip isothermal wall
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV : Conservative variables at the flux points
  ! WALL_TEMPERATURE : isothermal wall temperature
  !
  !.. Use Statements ..
  use gpu_mod, only : tvsp2u_gpu
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: wall_temperature
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: pv_in(:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    call get_bc_isothermal_wall_solution_3d(n1,n2,face_cv, &
                                            wall_temperature,pv_in)
  else if (nr == 2) then
    call get_bc_isothermal_wall_solution_2d(n1,n2,face_cv, &
                                            wall_temperature,pv_in)
  else
    call get_bc_isothermal_wall_solution_nd(n1,n2,face_cv, &
                                            wall_temperature,pv_in)
  end if
  !
  ! The above subroutines will return the boundary condition states as
  ! primitive variables. These need to be converted to conservative variables.
  !
  call tvsp2u_gpu(n1,n2,face_cv(:,:,2))
  !
end subroutine get_bc_isothermal_wall_solution
!
!###############################################################################
!
pure subroutine get_bc_isothermal_wall_solution_3d(n1,n2,face_cv, &
                                                   wall_temperature,pv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip isothermal wall
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV : Conservative variables at the flux points
  ! WALL_TEMPERATURE : isothermal wall temperature
  !
  !.. Use Statements ..
  use ovar, only : gam,aref,rgasref_nd
  use ovar, only : walls_are_exact
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: wall_temperature
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: pv_in(:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,rtw_inv
  real(wp) :: rho,rho_inv,vx,vy,vz,q2,pres
  real(wp) :: wall_vx,wall_vy,wall_vz
  real(wp) :: velocity_const
  real(wp) :: moving_wall_const
  !
continue
  !
  velocity_const = merge(zero,-one,walls_are_exact)
  moving_wall_const = merge(one,two,walls_are_exact)
  !
  gm1 = gam - one
  rtw_inv = one / rgasref_nd / wall_temperature
  !
  wall_vx = pv_in(nmx)
  wall_vy = pv_in(nmy)
  wall_vz = pv_in(nmz)
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv)
  do n = n1,n2
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    vz = face_cv(n,nmz,1)*rho_inv
    !
    q2 = vx*vx + vy*vy + vz*vz
    !
    pres = gm1*(face_cv(n,nee,1) - half*rho*q2)
    !
    face_cv(n,nec,2) = pres * rtw_inv
    face_cv(n,nmx,2) = vx*velocity_const + wall_vx*moving_wall_const
    face_cv(n,nmy,2) = vy*velocity_const + wall_vy*moving_wall_const
    face_cv(n,nmz,2) = vz*velocity_const + wall_vz*moving_wall_const
    face_cv(n,nee,2) = pres
    !
  end do
  !
end subroutine get_bc_isothermal_wall_solution_3d
!
!###############################################################################
!
pure subroutine get_bc_isothermal_wall_solution_2d(n1,n2,face_cv, &
                                                   wall_temperature,pv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip isothermal wall
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV : Conservative variables at the flux points
  ! WALL_TEMPERATURE : isothermal wall temperature
  !
  !.. Use Statements ..
  use ovar, only : gam,aref,rgasref_nd
  use ovar, only : walls_are_exact
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: wall_temperature
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: pv_in(:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,rtw_inv
  real(wp) :: rho,rho_inv,vx,vy,q2,pres
  real(wp) :: wall_vx,wall_vy
  real(wp) :: velocity_const
  real(wp) :: moving_wall_const
  !
continue
  !
  velocity_const = merge(zero,-one,walls_are_exact)
  moving_wall_const = merge(one,two,walls_are_exact)
  !
  gm1 = gam - one
  rtw_inv = one / rgasref_nd / wall_temperature
  !
  wall_vx = pv_in(nmx)
  wall_vy = pv_in(nmy)
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv)
  do n = n1,n2
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    !
    q2 = vx*vx + vy*vy
    !
    pres = gm1*(face_cv(n,nee,1) - half*rho*q2)
    !
    face_cv(n,nec,2) = pres * rtw_inv
    face_cv(n,nmx,2) = vx*velocity_const + wall_vx*moving_wall_const
    face_cv(n,nmy,2) = vy*velocity_const + wall_vy*moving_wall_const
    face_cv(n,nee,2) = pres
    !
  end do
  !
end subroutine get_bc_isothermal_wall_solution_2d
!
!###############################################################################
!
pure subroutine get_bc_isothermal_wall_solution_nd(n1,n2,face_cv, &
                                                   wall_temperature,pv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip isothermal wall
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV : Conservative variables at the flux points
  ! WALL_TEMPERATURE : isothermal wall temperature
  !
  !.. Use Statements ..
  use ovar, only : gam,aref,rgasref_nd
  use ovar, only : walls_are_exact
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: wall_temperature
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: pv_in(:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: gm1,rtw_inv
  real(wp) :: rho,rho_inv,vx,vy,vz,q2,pres
  real(wp) :: wall_vx,wall_vy,wall_vz
  real(wp) :: velocity_const
  real(wp) :: moving_wall_const
  !
continue
  !
  nr = nme-nmb+1
  !
  velocity_const = merge(zero,-one,walls_are_exact)
  moving_wall_const = merge(one,two,walls_are_exact)
  !
  gm1 = gam - one
  rtw_inv = one / rgasref_nd / wall_temperature
  !
  if (nr == 3) then
    !
    wall_vx = pv_in(nmx)
    wall_vy = pv_in(nmy)
    wall_vz = pv_in(nmz)
    !
    !$ACC PARALLEL LOOP PRESENT(face_cv)
    do n = n1,n2
      !
      rho = face_cv(n,nec,1)
      rho_inv = one/rho
      vx = face_cv(n,nmx,1)*rho_inv
      vy = face_cv(n,nmy,1)*rho_inv
      vz = face_cv(n,nmz,1)*rho_inv
      !
      q2 = vx*vx + vy*vy + vz*vz
      !
      pres = gm1*(face_cv(n,nee,1) - half*rho*q2)
      !
      face_cv(n,nec,2) = pres * rtw_inv
      face_cv(n,nmx,2) = vx*velocity_const + wall_vx*moving_wall_const
      face_cv(n,nmy,2) = vy*velocity_const + wall_vy*moving_wall_const
      face_cv(n,nmz,2) = vz*velocity_const + wall_vz*moving_wall_const
      face_cv(n,nee,2) = pres
      !
    end do
    !
  else if (nr == 2) then
    !
    wall_vx = pv_in(nmx)
    wall_vy = pv_in(nmy)
    !
    !$ACC PARALLEL LOOP PRESENT(face_cv)
    do n = n1,n2
      !
      rho = face_cv(n,nec,1)
      rho_inv = one/rho
      vx = face_cv(n,nmx,1)*rho_inv
      vy = face_cv(n,nmy,1)*rho_inv
      !
      q2 = vx*vx + vy*vy
      !
      pres = gm1*(face_cv(n,nee,1) - half*rho*q2)
      !
      face_cv(n,nec,2) = pres * rtw_inv
      face_cv(n,nmx,2) = vx*velocity_const + wall_vx*moving_wall_const
      face_cv(n,nmy,2) = vy*velocity_const + wall_vy*moving_wall_const
      face_cv(n,nee,2) = pres
      !
    end do
    !
  else
    !
    wall_vx = pv_in(nmx)
    !
    !$ACC PARALLEL LOOP PRESENT(face_cv)
    do n = n1,n2
      !
      rho = face_cv(n,nec,1)
      rho_inv = one/rho
      vx = face_cv(n,nmx,1)*rho_inv
      !
      pres = gm1*(face_cv(n,nee,1) - half*rho*vx*vx)
      !
      face_cv(n,nec,2) = pres * rtw_inv
      face_cv(n,nmx,2) = vx*velocity_const + wall_vx*moving_wall_const
      face_cv(n,nee,2) = pres
      !
    end do
    !
  end if
  !
end subroutine get_bc_isothermal_wall_solution_nd
!
!###############################################################################
!

pure subroutine get_bc_periodic_solution(n1,n2,face_cv,periodic_list)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: periodic_list(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    call get_bc_periodic_solution_3d(n1,n2,face_cv,periodic_list)
  else if (nr == 2) then
    call get_bc_periodic_solution_2d(n1,n2,face_cv,periodic_list)
  else
    call get_bc_periodic_solution_nd(n1,n2,face_cv,periodic_list)
  end if
  !
end subroutine get_bc_periodic_solution
!
!###############################################################################
!
pure subroutine get_bc_periodic_solution_3d(n1,n2,face_cv,periodic_list)
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
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: periodic_list(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  ! NOTE: We may need to reverse the copy order, i.e.,
  !          gufp(:,k) = pufp(:,nfp+1-k)
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,periodic_list)
  do n = n1,n2
    face_cv(n,nec,2) =  face_cv(periodic_list(n-n1+1),nec,1)
    face_cv(n,nmx,2) =  face_cv(periodic_list(n-n1+1),nmx,1)
    face_cv(n,nmy,2) =  face_cv(periodic_list(n-n1+1),nmy,1)
#ifdef REVERSE_PERIODIC
    face_cv(n,nmz,2) = -face_cv(periodic_list(n-n1+1),nmz,1)
#else
    face_cv(n,nmz,2) =  face_cv(periodic_list(n-n1+1),nmz,1)
#endif
    face_cv(n,nee,2) =  face_cv(periodic_list(n-n1+1),nee,1)
  end do
  !
end subroutine get_bc_periodic_solution_3d
!
!###############################################################################
!
pure subroutine get_bc_periodic_solution_2d(n1,n2,face_cv,periodic_list)
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
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: periodic_list(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  ! NOTE: We may need to reverse the copy order, i.e.,
  !          gufp(:,k) = pufp(:,nfp+1-k)
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,periodic_list)
  do n = n1,n2
    face_cv(n,nec,2) =  face_cv(periodic_list(n-n1+1),nec,1)
    face_cv(n,nmx,2) =  face_cv(periodic_list(n-n1+1),nmx,1)
#ifdef REVERSE_PERIODIC
    face_cv(n,nmy,2) = -face_cv(periodic_list(n-n1+1),nmy,1)
#else
    face_cv(n,nmy,2) =  face_cv(periodic_list(n-n1+1),nmy,1)
#endif
    face_cv(n,nee,2) =  face_cv(periodic_list(n-n1+1),nee,1)
  end do
  !
end subroutine get_bc_periodic_solution_2d
!
!###############################################################################
!
pure subroutine get_bc_periodic_solution_nd(n1,n2,face_cv,periodic_list)
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
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: periodic_list(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: n,m
  !
continue
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) &
  !$ACC&              PRESENT(face_cv,periodic_list)
  do m = 1,size(face_cv,dim=2)
    do n = n1,n2
      face_cv(n,m,2) = face_cv(periodic_list(n-n1+1),m,1)
    end do
  end do
#ifdef REVERSE_PERIODIC
  !
  ! NOTE: We may need to reverse the copy order, i.e.,
  !          gufp(:,k) = pufp(:,nfp+1-k)
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,periodic_list)
  do n = n1,n2
    face_cv(n,nme,2) = -face_cv(periodic_list(n-n1+1),nme,1)
  end do
#endif
  !
end subroutine get_bc_periodic_solution_nd
!
!###############################################################################
!
pure subroutine get_bc_mms_dirichlet_solution(n1,n2,facenrm,facexyz,face_cv)
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
  ! FACENRM : Face outward normal vector at the flux points
  ! FACEXYZ : XYZ coordinates for all the flux points
  ! FACE_CV : Conservative variables at the flux points
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: facexyz(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    call get_bc_mms_dirichlet_solution_3d(n1,n2,facenrm,facexyz,face_cv)
  else if (nr == 2) then
    call get_bc_mms_dirichlet_solution_2d(n1,n2,facenrm,facexyz,face_cv)
  else
    call get_bc_mms_dirichlet_solution_nd(n1,n2,facenrm,facexyz,face_cv)
  end if
  !
end subroutine get_bc_mms_dirichlet_solution
!
!###############################################################################
!
pure subroutine get_bc_mms_dirichlet_solution_3d(n1,n2,facenrm,facexyz,face_cv)
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
  ! FACENRM : Face outward normal vector at the flux points
  ! FACEXYZ : XYZ coordinates for all the flux points
  ! FACE_CV : Conservative variables at the flux points
  !
  !.. Use Statements ..
  use ovar,    only : gam
  use mms_mod, only : mms_solution_at_xyz_gpu
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: facexyz(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,gm1_inv
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: rho,rhou,rhov,rhow,rhoe,rhoq2,ske,pres,mach
  real(wp) :: rho_mms,rhou_mms,rhov_mms,rhow_mms,rhoe_mms
  real(wp) :: rhoq2_mms,ske_mms,pres_mms
  !
continue
  !
  gm1 = gam - one
  gm1_inv = one/gm1
  !
  call mms_solution_at_xyz_gpu(n1,n2,facexyz,face_cv(:,:,2), &
                               return_nondim=true, &
                               return_cv=true)
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
  do n = n1,n2
    !
    ! Unit normal at this flux point
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    rnz = facenrm(3,n)
    ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
    unx = rnx/ds
    uny = rny/ds
    unz = rnz/ds
    !
    rho_mms = face_cv(n,nec,2)
    rhou_mms = face_cv(n,nmx,2)
    rhov_mms = face_cv(n,nmy,2)
    rhow_mms = face_cv(n,nmz,2)
    rhoe_mms = face_cv(n,nee,2)
    rhoq2_mms = rhou_mms*rhou_mms + rhov_mms*rhov_mms + rhow_mms*rhow_mms
    ske_mms = half*rhoq2_mms/rho_mms
    pres_mms = gm1*(rhoe_mms - ske_mms)
    !
    ! Compute the Mach number of the normal velocity at this flux point
    !
    rho = face_cv(n,nec,1)
    rhou = face_cv(n,nmx,1)
    rhov = face_cv(n,nmy,1)
    rhow = face_cv(n,nmz,1)
    rhoe = face_cv(n,nee,1)
    rhoq2 = rhou*rhou + rhov*rhov + rhow*rhow
    ske = half*rhoq2/rho
    pres = gm1*(rhoe - ske)
    !
    ! NOTE: The speed of sound (denominator) is gam*pres*rho instead of
    !       gam*pres/rho because the normal velocity (numerator) is rho*vn
    !       instead of just vn.
    !
    mach = (rhou*unx + rhov*uny + rhow*unz) / sqrt( max(gam*pres*rho,rndoff) )
    !
    ! If mach >= one, do nothing because in that case we want to extrapolate
    ! everything from the interior and we already have those values stored
    ! in rho, rhou, rhov, rhow, and rhoe.
    !
    if (mach <= -one) then
      !
      ! Supersonic inflow
      !
      ! Set all variables to free stream
      !
      rho = rho_mms
      rhou = rhou_mms
      rhov = rhov_mms
      rhow = rhow_mms
      rhoe = rhoe_mms
      !
    else if (mach > -one .and. mach <= zero) then
      !
      ! Subsonic inflow
      !
      ! Set all variables but pressure to free stream
      !
      rho = rho_mms
      rhou = rhou_mms
      rhov = rhov_mms
      rhow = rhow_mms
      rhoe = pres*gm1_inv + ske_mms
      !
    else if (mach < one .and. mach > zero) then
      !
      ! Subsonic outflow
      !
      ! Extrapolate interior density and velocities to boundary
      ! Compute the total energy using free stream pressure
      ! and interior specific kinetic energy
      !
      rhoe = pres_mms*gm1_inv + ske
      !
    end if
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rhou
    face_cv(n,nmy,2) = rhov
    face_cv(n,nmz,2) = rhow
    face_cv(n,nee,2) = rhoe
    !
  end do
  !
end subroutine get_bc_mms_dirichlet_solution_3d
!
!###############################################################################
!
pure subroutine get_bc_mms_dirichlet_solution_2d(n1,n2,facenrm,facexyz,face_cv)
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
  ! FACENRM : Face outward normal vector at the flux points
  ! FACEXYZ : XYZ coordinates for all the flux points
  ! FACE_CV : Conservative variables at the flux points
  !
  !.. Use Statements ..
  use ovar,    only : gam
  use mms_mod, only : mms_solution_at_xyz_gpu
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: facexyz(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,gm1_inv
  real(wp) :: ds,rnx,rny,unx,uny
  real(wp) :: rho,rhou,rhov,rhoe,rhoq2,ske,pres,mach
  real(wp) :: rho_mms,rhou_mms,rhov_mms,rhoe_mms
  real(wp) :: rhoq2_mms,ske_mms,pres_mms
  !
continue
  !
  gm1 = gam - one
  gm1_inv = one/gm1
  !
  call mms_solution_at_xyz_gpu(n1,n2,facexyz,face_cv(:,:,2), &
                               return_nondim=true, &
                               return_cv=true)
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
  do n = n1,n2
    !
    ! Unit normal at this flux point
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds = sqrt( rnx*rnx + rny*rny )
    unx = rnx/ds
    uny = rny/ds
    !
    rho_mms = face_cv(n,nec,2)
    rhou_mms = face_cv(n,nmx,2)
    rhov_mms = face_cv(n,nmy,2)
    rhoe_mms = face_cv(n,nee,2)
    rhoq2_mms = rhou_mms*rhou_mms + rhov_mms*rhov_mms
    ske_mms = half*rhoq2_mms/rho_mms
    pres_mms = gm1*(rhoe_mms - ske_mms)
    !
    ! Compute the Mach number of the normal velocity at this flux point
    !
    rho = face_cv(n,nec,1)
    rhou = face_cv(n,nmx,1)
    rhov = face_cv(n,nmy,1)
    rhoe = face_cv(n,nee,1)
    rhoq2 = rhou*rhou + rhov*rhov
    ske = half*rhoq2/rho
    pres = gm1*(rhoe - ske)
    !
    ! NOTE: The speed of sound (denominator) is gam*pres*rho instead of
    !       gam*pres/rho because the normal velocity (numerator) is rho*vn
    !       instead of just vn.
    !
    mach = (rhou*unx + rhov*uny) / sqrt( max(gam*pres*rho,rndoff) )
    !
    ! If mach >= one, do nothing because in that case we want to extrapolate
    ! everything from the interior and we already have those values stored
    ! in rho, rhou, rhov, rhow, and rhoe.
    !
    if (mach <= -one) then
      !
      ! Supersonic inflow
      !
      ! Set all variables to free stream
      !
      rho = rho_mms
      rhou = rhou_mms
      rhov = rhov_mms
      rhoe = rhoe_mms
      !
    else if (mach > -one .and. mach <= zero) then
      !
      ! Subsonic inflow
      !
      ! Set all variables but pressure to free stream
      !
      rho = rho_mms
      rhou = rhou_mms
      rhov = rhov_mms
      rhoe = pres*gm1_inv + ske_mms
      !
    else if (mach < one .and. mach > zero) then
      !
      ! Subsonic outflow
      !
      ! Extrapolate interior density and velocities to boundary
      ! Compute the total energy using free stream pressure
      ! and interior specific kinetic energy
      !
      rhoe = pres_mms*gm1_inv + ske
      !
    end if
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rhou
    face_cv(n,nmy,2) = rhov
    face_cv(n,nee,2) = rhoe
    !
  end do
  !
end subroutine get_bc_mms_dirichlet_solution_2d
!
!###############################################################################
!
pure subroutine get_bc_mms_dirichlet_solution_nd(n1,n2,facenrm,facexyz,face_cv)
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
  ! FACENRM : Face outward normal vector at the flux points
  ! FACEXYZ : XYZ coordinates for all the flux points
  ! FACE_CV : Conservative variables at the flux points
  !
  !.. Use Statements ..
  use ovar,    only : gam
  use mms_mod, only : mms_solution_at_xyz_gpu
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: facexyz(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: gm1,gm1_inv
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: rho,rhou,rhov,rhow,rhoe,rhoq2,ske,pres,mach
  real(wp) :: rho_mms,rhou_mms,rhov_mms,rhow_mms,rhoe_mms
  real(wp) :: rhoq2_mms,ske_mms,pres_mms
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  gm1 = gam - one
  gm1_inv = one/gm1
  !
  call mms_solution_at_xyz_gpu(n1,n2,facexyz,face_cv(:,:,2), &
                               return_nondim=true, &
                               return_cv=true)
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
    do n = n1,n2
      !
      ! Unit normal at this flux point
      !
      rnx = facenrm(1,n)
      rny = facenrm(2,n)
      rnz = facenrm(3,n)
      ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
      unx = rnx/ds
      uny = rny/ds
      unz = rnz/ds
      !
      rho_mms = face_cv(n,nec,2)
      rhou_mms = face_cv(n,nmx,2)
      rhov_mms = face_cv(n,nmy,2)
      rhow_mms = face_cv(n,nmz,2)
      rhoe_mms = face_cv(n,nee,2)
      rhoq2_mms = rhou_mms*rhou_mms + rhov_mms*rhov_mms + rhow_mms*rhow_mms
      ske_mms = half*rhoq2_mms/rho_mms
      pres_mms = gm1*(rhoe_mms - ske_mms)
      !
      ! Compute the Mach number of the normal velocity at this flux point
      !
      rho = face_cv(n,nec,1)
      rhou = face_cv(n,nmx,1)
      rhov = face_cv(n,nmy,1)
      rhow = face_cv(n,nmz,1)
      rhoe = face_cv(n,nee,1)
      rhoq2 = rhou*rhou + rhov*rhov + rhow*rhow
      ske = half*rhoq2/rho
      pres = gm1*(rhoe - ske)
      !
      ! NOTE: The speed of sound (denominator) is gam*pres*rho instead of
      !       gam*pres/rho because the normal velocity (numerator) is rho*vn
      !       instead of just vn.
      !
      mach = (rhou*unx + rhov*uny + rhow*unz) / sqrt( max(gam*pres*rho,rndoff) )
      !
      ! If mach >= one, do nothing because in that case we want to extrapolate
      ! everything from the interior and we already have those values stored
      ! in rho, rhou, rhov, rhow, and rhoe.
      !
      if (mach <= -one) then
        !
        ! Supersonic inflow
        !
        ! Set all variables to free stream
        !
        rho = rho_mms
        rhou = rhou_mms
        rhov = rhov_mms
        rhow = rhow_mms
        rhoe = rhoe_mms
        !
      else if (mach > -one .and. mach <= zero) then
        !
        ! Subsonic inflow
        !
        ! Set all variables but pressure to free stream
        !
        rho = rho_mms
        rhou = rhou_mms
        rhov = rhov_mms
        rhow = rhow_mms
        rhoe = pres*gm1_inv + ske_mms
        !
      else if (mach < one .and. mach > zero) then
        !
        ! Subsonic outflow
        !
        ! Extrapolate interior density and velocities to boundary
        ! Compute the total energy using free stream pressure
        ! and interior specific kinetic energy
        !
        rhoe = pres_mms*gm1_inv + ske
        !
      end if
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rhou
      face_cv(n,nmy,2) = rhov
      face_cv(n,nmz,2) = rhow
      face_cv(n,nee,2) = rhoe
      !
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
    do n = n1,n2
      !
      ! Unit normal at this flux point
      !
      rnx = facenrm(1,n)
      rny = facenrm(2,n)
      ds = sqrt( rnx*rnx + rny*rny )
      unx = rnx/ds
      uny = rny/ds
      !
      rho_mms = face_cv(n,nec,2)
      rhou_mms = face_cv(n,nmx,2)
      rhov_mms = face_cv(n,nmy,2)
      rhoe_mms = face_cv(n,nee,2)
      rhoq2_mms = rhou_mms*rhou_mms + rhov_mms*rhov_mms
      ske_mms = half*rhoq2_mms/rho_mms
      pres_mms = gm1*(rhoe_mms - ske_mms)
      !
      ! Compute the Mach number of the normal velocity at this flux point
      !
      rho = face_cv(n,nec,1)
      rhou = face_cv(n,nmx,1)
      rhov = face_cv(n,nmy,1)
      rhoe = face_cv(n,nee,1)
      rhoq2 = rhou*rhou + rhov*rhov
      ske = half*rhoq2/rho
      pres = gm1*(rhoe - ske)
      !
      ! NOTE: The speed of sound (denominator) is gam*pres*rho instead of
      !       gam*pres/rho because the normal velocity (numerator) is rho*vn
      !       instead of just vn.
      !
      mach = (rhou*unx + rhov*uny) / sqrt( max(gam*pres*rho,rndoff) )
      !
      ! If mach >= one, do nothing because in that case we want to extrapolate
      ! everything from the interior and we already have those values stored
      ! in rho, rhou, rhov, rhow, and rhoe.
      !
      if (mach <= -one) then
        !
        ! Supersonic inflow
        !
        ! Set all variables to free stream
        !
        rho = rho_mms
        rhou = rhou_mms
        rhov = rhov_mms
        rhoe = rhoe_mms
        !
      else if (mach > -one .and. mach <= zero) then
        !
        ! Subsonic inflow
        !
        ! Set all variables but pressure to free stream
        !
        rho = rho_mms
        rhou = rhou_mms
        rhov = rhov_mms
        rhoe = pres*gm1_inv + ske_mms
        !
      else if (mach < one .and. mach > zero) then
        !
        ! Subsonic outflow
        !
        ! Extrapolate interior density and velocities to boundary
        ! Compute the total energy using free stream pressure
        ! and interior specific kinetic energy
        !
        rhoe = pres_mms*gm1_inv + ske
        !
      end if
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rhou
      face_cv(n,nmy,2) = rhov
      face_cv(n,nee,2) = rhoe
      !
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
    do n = n1,n2
      !
      ! Unit normal at this flux point
      !
      unx = facenrm(1,n)/abs(facenrm(1,n))
      !
      rho_mms = face_cv(n,nec,2)
      rhou_mms = face_cv(n,nmx,2)
      rhoe_mms = face_cv(n,nee,2)
      ske_mms = half*rhou_mms*rhou_mms/rho_mms
      pres_mms = gm1*(rhoe_mms - ske_mms)
      !
      ! Compute the Mach number of the normal velocity at this flux point
      !
      rho = face_cv(n,nec,1)
      rhou = face_cv(n,nmx,1)
      rhoe = face_cv(n,nee,1)
      ske = half*rhou*rhou/rho
      pres = gm1*(rhoe - ske)
      !
      ! NOTE: The speed of sound (denominator) is gam*pres*rho instead of
      !       gam*pres/rho because the normal velocity (numerator) is rho*vn
      !       instead of just vn.
      !
      mach = rhou*unx / sqrt( max(gam*pres*rho,rndoff) )
      !
      ! If mach >= one, do nothing because in that case we want to extrapolate
      ! everything from the interior and we already have those values stored
      ! in rho, rhou, rhov, rhow, and rhoe.
      !
      if (mach <= -one) then
        !
        ! Supersonic inflow
        !
        ! Set all variables to free stream
        !
        rho = rho_mms
        rhou = rhou_mms
        rhoe = rhoe_mms
        !
      else if (mach > -one .and. mach <= zero) then
        !
        ! Subsonic inflow
        !
        ! Set all variables but pressure to free stream
        !
        rho = rho_mms
        rhou = rhou_mms
        rhoe = pres*gm1_inv + ske_mms
        !
      else if (mach < one .and. mach > zero) then
        !
        ! Subsonic outflow
        !
        ! Extrapolate interior density and velocities to boundary
        ! Compute the total energy using free stream pressure
        ! and interior specific kinetic energy
        !
        rhoe = pres_mms*gm1_inv + ske
        !
      end if
      !
      face_cv(n,nec,2) = rho
      face_cv(n,nmx,2) = rhou
      face_cv(n,nee,2) = rhoe
      !
    end do
  end if
  !
end subroutine get_bc_mms_dirichlet_solution_nd
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
pure subroutine get_bc_periodic_gradient(n1,n2,face_dcvdx,periodic_list)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: periodic_list(:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = size(face_dcvdx,dim=2)
  !
  if (nr == 3) then
    call get_bc_periodic_gradient_3d(n1,n2,face_dcvdx,periodic_list)
  else if (nr == 2) then
    call get_bc_periodic_gradient_2d(n1,n2,face_dcvdx,periodic_list)
  else
    call get_bc_periodic_gradient_nd(n1,n2,face_dcvdx,periodic_list)
  end if
  !
end subroutine get_bc_periodic_gradient
!
!###############################################################################
!
pure subroutine get_bc_periodic_gradient_3d(n1,n2,face_dcvdx,periodic_list)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: periodic_list(:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  ! Simply copy the interior gradients to the boundary flux points
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx,periodic_list)
  do n = n1,n2
    !
    face_dcvdx(n,1,nec,2) = face_dcvdx(periodic_list(n-n1+1),1,nec,1)
    face_dcvdx(n,2,nec,2) = face_dcvdx(periodic_list(n-n1+1),2,nec,1)
    face_dcvdx(n,3,nec,2) = face_dcvdx(periodic_list(n-n1+1),3,nec,1)
    !
    face_dcvdx(n,1,nmx,2) = face_dcvdx(periodic_list(n-n1+1),1,nmx,1)
    face_dcvdx(n,2,nmx,2) = face_dcvdx(periodic_list(n-n1+1),2,nmx,1)
    face_dcvdx(n,3,nmx,2) = face_dcvdx(periodic_list(n-n1+1),3,nmx,1)
    !
    face_dcvdx(n,1,nmy,2) = face_dcvdx(periodic_list(n-n1+1),1,nmy,1)
    face_dcvdx(n,2,nmy,2) = face_dcvdx(periodic_list(n-n1+1),2,nmy,1)
    face_dcvdx(n,3,nmy,2) = face_dcvdx(periodic_list(n-n1+1),3,nmy,1)
    !
    face_dcvdx(n,1,nmz,2) = face_dcvdx(periodic_list(n-n1+1),1,nmz,1)
    face_dcvdx(n,2,nmz,2) = face_dcvdx(periodic_list(n-n1+1),2,nmz,1)
    face_dcvdx(n,3,nmz,2) = face_dcvdx(periodic_list(n-n1+1),3,nmz,1)
    !
    face_dcvdx(n,1,nee,2) = face_dcvdx(periodic_list(n-n1+1),1,nee,1)
    face_dcvdx(n,2,nee,2) = face_dcvdx(periodic_list(n-n1+1),2,nee,1)
    face_dcvdx(n,3,nee,2) = face_dcvdx(periodic_list(n-n1+1),3,nee,1)
    !
  end do
  !
end subroutine get_bc_periodic_gradient_3d
!
!###############################################################################
!
pure subroutine get_bc_periodic_gradient_2d(n1,n2,face_dcvdx,periodic_list)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: periodic_list(:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  ! Simply copy the interior gradients to the boundary flux points
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx,periodic_list)
  do n = n1,n2
    !
    face_dcvdx(n,1,nec,2) = face_dcvdx(periodic_list(n-n1+1),1,nec,1)
    face_dcvdx(n,2,nec,2) = face_dcvdx(periodic_list(n-n1+1),2,nec,1)
    !
    face_dcvdx(n,1,nmx,2) = face_dcvdx(periodic_list(n-n1+1),1,nmx,1)
    face_dcvdx(n,2,nmx,2) = face_dcvdx(periodic_list(n-n1+1),2,nmx,1)
    !
    face_dcvdx(n,1,nmy,2) = face_dcvdx(periodic_list(n-n1+1),1,nmy,1)
    face_dcvdx(n,2,nmy,2) = face_dcvdx(periodic_list(n-n1+1),2,nmy,1)
    !
    face_dcvdx(n,1,nee,2) = face_dcvdx(periodic_list(n-n1+1),1,nee,1)
    face_dcvdx(n,2,nee,2) = face_dcvdx(periodic_list(n-n1+1),2,nee,1)
    !
  end do
  !
end subroutine get_bc_periodic_gradient_2d
!
!###############################################################################
!
pure subroutine get_bc_periodic_gradient_nd(n1,n2,face_dcvdx,periodic_list)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: periodic_list(:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: n,m,l
  !
continue
  !
  ! Simply copy the interior gradients to the boundary flux points
  !
  !$ACC PARALLEL LOOP COLLAPSE(3) &
  !$ACC&              PRESENT(face_dcvdx,periodic_list)
  do m = 1,size(face_dcvdx,dim=3)
    do l = 1,size(face_dcvdx,dim=2)
      do n = n1,n2
        face_dcvdx(n,l,m,2) = face_dcvdx(periodic_list(n-n1+1),l,m,1)
      end do
    end do
  end do
  !
end subroutine get_bc_periodic_gradient_nd
!
!###############################################################################
!
pure subroutine get_bc_flowthrough_gradient(n1,n2,face_dcvdx)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = size(face_dcvdx,dim=2)
  !
  if (nr == 3) then
    call get_bc_flowthrough_gradient_3d(n1,n2,face_dcvdx)
  else if (nr == 2) then
    call get_bc_flowthrough_gradient_2d(n1,n2,face_dcvdx)
  else
    call get_bc_flowthrough_gradient_nd(n1,n2,face_dcvdx)
  end if
  !
end subroutine get_bc_flowthrough_gradient
!
!###############################################################################
!
pure subroutine get_bc_flowthrough_gradient_3d(n1,n2,face_dcvdx)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  ! Simply copy the interior gradients to the boundary flux points
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
  do n = n1,n2
    !
    face_dcvdx(n,1,nec,2) = face_dcvdx(n,1,nec,1)
    face_dcvdx(n,2,nec,2) = face_dcvdx(n,2,nec,1)
    face_dcvdx(n,3,nec,2) = face_dcvdx(n,3,nec,1)
    !
    face_dcvdx(n,1,nmx,2) = face_dcvdx(n,1,nmx,1)
    face_dcvdx(n,2,nmx,2) = face_dcvdx(n,2,nmx,1)
    face_dcvdx(n,3,nmx,2) = face_dcvdx(n,3,nmx,1)
    !
    face_dcvdx(n,1,nmy,2) = face_dcvdx(n,1,nmy,1)
    face_dcvdx(n,2,nmy,2) = face_dcvdx(n,2,nmy,1)
    face_dcvdx(n,3,nmy,2) = face_dcvdx(n,3,nmy,1)
    !
    face_dcvdx(n,1,nmz,2) = face_dcvdx(n,1,nmz,1)
    face_dcvdx(n,2,nmz,2) = face_dcvdx(n,2,nmz,1)
    face_dcvdx(n,3,nmz,2) = face_dcvdx(n,3,nmz,1)
    !
    face_dcvdx(n,1,nee,2) = face_dcvdx(n,1,nee,1)
    face_dcvdx(n,2,nee,2) = face_dcvdx(n,2,nee,1)
    face_dcvdx(n,3,nee,2) = face_dcvdx(n,3,nee,1)
    !
  end do
  !
end subroutine get_bc_flowthrough_gradient_3d
!
!###############################################################################
!
pure subroutine get_bc_flowthrough_gradient_2d(n1,n2,face_dcvdx)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  ! Simply copy the interior gradients to the boundary flux points
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
  do n = n1,n2
    !
    face_dcvdx(n,1,nec,2) = face_dcvdx(n,1,nec,1)
    face_dcvdx(n,2,nec,2) = face_dcvdx(n,2,nec,1)
    !
    face_dcvdx(n,1,nmx,2) = face_dcvdx(n,1,nmx,1)
    face_dcvdx(n,2,nmx,2) = face_dcvdx(n,2,nmx,1)
    !
    face_dcvdx(n,1,nmy,2) = face_dcvdx(n,1,nmy,1)
    face_dcvdx(n,2,nmy,2) = face_dcvdx(n,2,nmy,1)
    !
    face_dcvdx(n,1,nee,2) = face_dcvdx(n,1,nee,1)
    face_dcvdx(n,2,nee,2) = face_dcvdx(n,2,nee,1)
    !
  end do
  !
end subroutine get_bc_flowthrough_gradient_2d
!
!###############################################################################
!
pure subroutine get_bc_flowthrough_gradient_nd(n1,n2,face_dcvdx)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: n,m,l
  !
continue
  !
  ! Simply copy the interior gradients to the boundary flux points
  !
  !$ACC PARALLEL LOOP COLLAPSE(3) &
  !$ACC&              PRESENT(face_dcvdx)
  do m = 1,size(face_dcvdx,dim=3)
    do l = 1,size(face_dcvdx,dim=2)
      do n = n1,n2
        face_dcvdx(n,l,m,2) = face_dcvdx(n,l,m,1)
      end do
    end do
  end do
  !
end subroutine get_bc_flowthrough_gradient_nd
!
!###############################################################################
!
pure subroutine get_bc_adiabatic_wall_gradient(n1,n2,face_cv,face_dcvdx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip adiabatic wall
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV    : Conservative variables at the host cell flux points
  ! FACE_DCVDX : Gradients of the conservative variables at the flux pts
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = size(face_dcvdx,dim=2)
  !
  if (nr == 3) then
    call get_bc_adiabatic_wall_gradient_3d(n1,n2,face_cv,face_dcvdx)
  else if (nr == 2) then
    call get_bc_adiabatic_wall_gradient_2d(n1,n2,face_cv,face_dcvdx)
  else
    call get_bc_adiabatic_wall_gradient_nd(n1,n2,face_cv,face_dcvdx)
  end if
  !
end subroutine get_bc_adiabatic_wall_gradient
!
!###############################################################################
!
pure subroutine get_bc_adiabatic_wall_gradient_3d(n1,n2,face_cv,face_dcvdx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip adiabatic wall
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV    : Conservative variables at the host cell flux points
  ! FACE_DCVDX : Gradients of the conservative variables at the flux pts
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: r_inv
  !
continue
  !
  ! Copy over the interior density gradient to the boundary flux points
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
  do n = n1,n2
    face_dcvdx(n,1,nec,2) = face_dcvdx(n,1,nec,1)
    face_dcvdx(n,2,nec,2) = face_dcvdx(n,2,nec,1)
    face_dcvdx(n,3,nec,2) = face_dcvdx(n,3,nec,1)
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
  ! Using the chain rule,convert the momentum gradients to velocity gradients
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx,face_cv)
  do n = n1,n2
    !
    r_inv = one/face_cv(n,nec)
    !
    face_dcvdx(n,1,nmx,2) = r_inv*face_dcvdx(n,1,nmx,1) &
                          - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nmx)*r_inv
    face_dcvdx(n,2,nmx,2) = r_inv*face_dcvdx(n,2,nmx,1) &
                          - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nmx)*r_inv
    face_dcvdx(n,3,nmx,2) = r_inv*face_dcvdx(n,3,nmx,1) &
                          - r_inv*face_dcvdx(n,3,nec,1)*face_cv(n,nmx)*r_inv
    !
    face_dcvdx(n,1,nmy,2) = r_inv*face_dcvdx(n,1,nmy,1) &
                          - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nmy)*r_inv
    face_dcvdx(n,2,nmy,2) = r_inv*face_dcvdx(n,2,nmy,1) &
                          - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nmy)*r_inv
    face_dcvdx(n,3,nmy,2) = r_inv*face_dcvdx(n,3,nmy,1) &
                          - r_inv*face_dcvdx(n,3,nec,1)*face_cv(n,nmy)*r_inv
    !
    face_dcvdx(n,1,nmz,2) = r_inv*face_dcvdx(n,1,nmz,1) &
                          - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nmz)*r_inv
    face_dcvdx(n,2,nmz,2) = r_inv*face_dcvdx(n,2,nmz,1) &
                          - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nmz)*r_inv
    face_dcvdx(n,3,nmz,2) = r_inv*face_dcvdx(n,3,nmz,1) &
                          - r_inv*face_dcvdx(n,3,nec,1)*face_cv(n,nmz)*r_inv
    !
  end do
  !
  ! Finally, since this is an adiabatic wall boundary,
  ! set the temperature gradients to zero.
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
  do n = n1,n2
    face_dcvdx(n,1,nee,2) = zero
    face_dcvdx(n,2,nee,2) = zero
    face_dcvdx(n,3,nee,2) = zero
  end do
  !
end subroutine get_bc_adiabatic_wall_gradient_3d
!
!###############################################################################
!
pure subroutine get_bc_adiabatic_wall_gradient_2d(n1,n2,face_cv,face_dcvdx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip adiabatic wall
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV    : Conservative variables at the host cell flux points
  ! FACE_DCVDX : Gradients of the conservative variables at the flux pts
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: r_inv
  !
continue
  !
  ! Copy over the interior density gradient to the boundary flux points
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
  do n = n1,n2
    face_dcvdx(n,1,nec,2) = face_dcvdx(n,1,nec,1)
    face_dcvdx(n,2,nec,2) = face_dcvdx(n,2,nec,1)
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
  ! Using the chain rule,convert the momentum gradients to velocity gradients
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx,face_cv)
  do n = n1,n2
    !
    r_inv = one/face_cv(n,nec)
    !
    face_dcvdx(n,1,nmx,2) = r_inv*face_dcvdx(n,1,nmx,1) &
                          - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nmx)*r_inv
    face_dcvdx(n,2,nmx,2) = r_inv*face_dcvdx(n,2,nmx,1) &
                          - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nmx)*r_inv
    !
    face_dcvdx(n,1,nmy,2) = r_inv*face_dcvdx(n,1,nmy,1) &
                          - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nmy)*r_inv
    face_dcvdx(n,2,nmy,2) = r_inv*face_dcvdx(n,2,nmy,1) &
                          - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nmy)*r_inv
    !
  end do
  !
  ! Finally, since this is an adiabatic wall boundary,
  ! set the temperature gradients to zero.
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
  do n = n1,n2
    face_dcvdx(n,1,nee,2) = zero
    face_dcvdx(n,2,nee,2) = zero
  end do
  !
end subroutine get_bc_adiabatic_wall_gradient_2d
!
!###############################################################################
!
pure subroutine get_bc_adiabatic_wall_gradient_nd(n1,n2,face_cv,face_dcvdx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip adiabatic wall
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV    : Conservative variables at the host cell flux points
  ! FACE_DCVDX : Gradients of the conservative variables at the flux pts
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: r_inv
  !
continue
  !
  nr = size(face_dcvdx,dim=2)
  !
  ! Copy over the interior density gradient to the boundary flux points
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
    do n = n1,n2
      face_dcvdx(n,1,nec,2) = face_dcvdx(n,1,nec,1)
      face_dcvdx(n,2,nec,2) = face_dcvdx(n,2,nec,1)
      face_dcvdx(n,3,nec,2) = face_dcvdx(n,3,nec,1)
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
    do n = n1,n2
      face_dcvdx(n,1,nec,2) = face_dcvdx(n,1,nec,1)
      face_dcvdx(n,2,nec,2) = face_dcvdx(n,2,nec,1)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
    do n = n1,n2
      face_dcvdx(n,1,nec,2) = face_dcvdx(n,1,nec,1)
    end do
  end if
  !
  ! We will also copy the other conservative variable gradients from
  ! the interior cell to the boundary flux points. However, since
  ! we treat the viscous fluxes for wall boundaries differently from
  ! the flow-through boundaries, we need to ensure that the wall
  ! gradients are in the right form for the wall viscous flux function.
  ! This means pre-converting the conservative variable gradients to
  ! primitive variable gradients.
  !
  ! Using the chain rule,convert the momentum gradients to velocity gradients
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(face_dcvdx,face_cv)
    do n = n1,n2
      !
      r_inv = one/face_cv(n,nec)
      !
      face_dcvdx(n,1,nmx,2) = r_inv*face_dcvdx(n,1,nmx,1) &
                            - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nmx)*r_inv
      face_dcvdx(n,2,nmx,2) = r_inv*face_dcvdx(n,2,nmx,1) &
                            - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nmx)*r_inv
      face_dcvdx(n,3,nmx,2) = r_inv*face_dcvdx(n,3,nmx,1) &
                            - r_inv*face_dcvdx(n,3,nec,1)*face_cv(n,nmx)*r_inv
      !
      face_dcvdx(n,1,nmy,2) = r_inv*face_dcvdx(n,1,nmy,1) &
                            - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nmy)*r_inv
      face_dcvdx(n,2,nmy,2) = r_inv*face_dcvdx(n,2,nmy,1) &
                            - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nmy)*r_inv
      face_dcvdx(n,3,nmy,2) = r_inv*face_dcvdx(n,3,nmy,1) &
                            - r_inv*face_dcvdx(n,3,nec,1)*face_cv(n,nmy)*r_inv
      !
      face_dcvdx(n,1,nmz,2) = r_inv*face_dcvdx(n,1,nmz,1) &
                            - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nmz)*r_inv
      face_dcvdx(n,2,nmz,2) = r_inv*face_dcvdx(n,2,nmz,1) &
                            - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nmz)*r_inv
      face_dcvdx(n,3,nmz,2) = r_inv*face_dcvdx(n,3,nmz,1) &
                            - r_inv*face_dcvdx(n,3,nec,1)*face_cv(n,nmz)*r_inv
      !
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(face_dcvdx,face_cv)
    do n = n1,n2
      !
      r_inv = one/face_cv(n,nec)
      !
      face_dcvdx(n,1,nmx,2) = r_inv*face_dcvdx(n,1,nmx,1) &
                            - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nmx)*r_inv
      face_dcvdx(n,2,nmx,2) = r_inv*face_dcvdx(n,2,nmx,1) &
                            - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nmx)*r_inv
      !
      face_dcvdx(n,1,nmy,2) = r_inv*face_dcvdx(n,1,nmy,1) &
                            - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nmy)*r_inv
      face_dcvdx(n,2,nmy,2) = r_inv*face_dcvdx(n,2,nmy,1) &
                            - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nmy)*r_inv
      !
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(face_dcvdx,face_cv)
    do n = n1,n2
      !
      r_inv = one/face_cv(n,nec)
      !
      face_dcvdx(n,1,nmx,2) = r_inv*face_dcvdx(n,1,nmx,1) &
                            - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nmx)*r_inv
      !
    end do
  end if
  !
  ! Finally, since this is an adiabatic wall boundary,
  ! set the temperature gradients to zero.
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
    do n = n1,n2
      face_dcvdx(n,1,nee,2) = zero
      face_dcvdx(n,2,nee,2) = zero
      face_dcvdx(n,3,nee,2) = zero
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
    do n = n1,n2
      face_dcvdx(n,1,nee,2) = zero
      face_dcvdx(n,2,nee,2) = zero
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
    do n = n1,n2
      face_dcvdx(n,1,nee,2) = zero
    end do
  end if
  !
end subroutine get_bc_adiabatic_wall_gradient_nd
!
!###############################################################################
!
pure subroutine get_bc_wall_gradient(n1,n2,face_cv,face_dcvdx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Non-adiabatic wall
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV    : Conservative variables at the host cell flux points
  ! FACE_DCVDX : Gradients of the conservative variables at the flux pts
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = size(face_dcvdx,dim=2)
  !
  if (nr == 3) then
    call get_bc_wall_gradient_3d(n1,n2,face_cv,face_dcvdx)
  else if (nr == 2) then
    call get_bc_wall_gradient_2d(n1,n2,face_cv,face_dcvdx)
  else
    call get_bc_wall_gradient_nd(n1,n2,face_cv,face_dcvdx)
  end if
  !
end subroutine get_bc_wall_gradient
!
!###############################################################################
!
pure subroutine get_bc_wall_gradient_3d(n1,n2,face_cv,face_dcvdx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Non-adiabatic wall
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV    : Conservative variables at the host cell flux points
  ! FACE_DCVDX : Gradients of the conservative variables at the flux pts
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: r_inv,cv_inv
  !
continue
  !
  ! Precompute the inverse of specific heat at constant volume
  !
  cv_inv = (gam - one) / rgasref_nd
  !
  ! Copy over the interior density gradient to the boundary flux points
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
  do n = n1,n2
    face_dcvdx(n,1,nec,2) = face_dcvdx(n,1,nec,1)
    face_dcvdx(n,2,nec,2) = face_dcvdx(n,2,nec,1)
    face_dcvdx(n,3,nec,2) = face_dcvdx(n,3,nec,1)
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
  ! Using the chain rule:
  !
  ! 1. Compute the velocity gradients (m = 2,nq-1)
  ! 2. Compute the total energy gradient and store in the
  !    location for the temperature gradient (m = nq)
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx,face_cv)
  do n = n1,n2
    !
    r_inv = one/face_cv(n,nec)
    !
    face_dcvdx(n,1,nmx,2) = r_inv*face_dcvdx(n,1,nmx,1) &
                          - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nmx)*r_inv
    face_dcvdx(n,2,nmx,2) = r_inv*face_dcvdx(n,2,nmx,1) &
                          - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nmx)*r_inv
    face_dcvdx(n,3,nmx,2) = r_inv*face_dcvdx(n,3,nmx,1) &
                          - r_inv*face_dcvdx(n,3,nec,1)*face_cv(n,nmx)*r_inv
    !
    face_dcvdx(n,1,nmy,2) = r_inv*face_dcvdx(n,1,nmy,1) &
                          - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nmy)*r_inv
    face_dcvdx(n,2,nmy,2) = r_inv*face_dcvdx(n,2,nmy,1) &
                          - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nmy)*r_inv
    face_dcvdx(n,3,nmy,2) = r_inv*face_dcvdx(n,3,nmy,1) &
                          - r_inv*face_dcvdx(n,3,nec,1)*face_cv(n,nmy)*r_inv
    !
    face_dcvdx(n,1,nmz,2) = r_inv*face_dcvdx(n,1,nmz,1) &
                          - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nmz)*r_inv
    face_dcvdx(n,2,nmz,2) = r_inv*face_dcvdx(n,2,nmz,1) &
                          - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nmz)*r_inv
    face_dcvdx(n,3,nmz,2) = r_inv*face_dcvdx(n,3,nmz,1) &
                          - r_inv*face_dcvdx(n,3,nec,1)*face_cv(n,nmz)*r_inv
    !
    face_dcvdx(n,1,nee,2) = r_inv*face_dcvdx(n,1,nee,1) &
                          - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nee)*r_inv
    face_dcvdx(n,2,nee,2) = r_inv*face_dcvdx(n,2,nee,1) &
                          - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nee)*r_inv
    face_dcvdx(n,3,nee,2) = r_inv*face_dcvdx(n,3,nee,1) &
                          - r_inv*face_dcvdx(n,3,nec,1)*face_cv(n,nee)*r_inv
    !
  end do
  !
  ! The location for the temperature gradient currently contains
  ! the total energy gradient from the previous loop.
  ! Now use the velocity gradients to subtract the gradient of
  ! kinetic-energy-per-unit-mass from the total energy gradient.
  !
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx,face_cv)
  do n = n1,n2
    !
    r_inv = one/face_cv(n,nec)
    !
    face_dcvdx(n,1,nee,2) = face_dcvdx(n,1,nee,2) &
                          - face_dcvdx(n,1,nmx,2)*face_cv(n,nmx)*r_inv &
                          - face_dcvdx(n,1,nmy,2)*face_cv(n,nmy)*r_inv &
                          - face_dcvdx(n,1,nmz,2)*face_cv(n,nmz)*r_inv
    !
    face_dcvdx(n,2,nee,2) = face_dcvdx(n,2,nee,2) &
                          - face_dcvdx(n,2,nmx,2)*face_cv(n,nmx)*r_inv &
                          - face_dcvdx(n,2,nmy,2)*face_cv(n,nmy)*r_inv &
                          - face_dcvdx(n,2,nmz,2)*face_cv(n,nmz)*r_inv
    !
    face_dcvdx(n,3,nee,2) = face_dcvdx(n,3,nee,2) &
                          - face_dcvdx(n,3,nmx,2)*face_cv(n,nmx)*r_inv &
                          - face_dcvdx(n,3,nmy,2)*face_cv(n,nmy)*r_inv &
                          - face_dcvdx(n,3,nmz,2)*face_cv(n,nmz)*r_inv
    !
  end do
  !
  ! Finish computing the temperature gradient by dividing through by
  ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
  do n = n1,n2
    face_dcvdx(n,1,nee,2) = face_dcvdx(n,1,nee,2) * cv_inv
    face_dcvdx(n,2,nee,2) = face_dcvdx(n,2,nee,2) * cv_inv
    face_dcvdx(n,3,nee,2) = face_dcvdx(n,3,nee,2) * cv_inv
  end do
  !
end subroutine get_bc_wall_gradient_3d
!
!###############################################################################
!
pure subroutine get_bc_wall_gradient_2d(n1,n2,face_cv,face_dcvdx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Non-adiabatic wall
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV    : Conservative variables at the host cell flux points
  ! FACE_DCVDX : Gradients of the conservative variables at the flux pts
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: r_inv,cv_inv
  !
continue
  !
  ! Precompute the inverse of specific heat at constant volume
  !
  cv_inv = (gam - one) / rgasref_nd
  !
  ! Copy over the interior density gradient to the boundary flux points
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
  do n = n1,n2
    face_dcvdx(n,1,nec,2) = face_dcvdx(n,1,nec,1)
    face_dcvdx(n,2,nec,2) = face_dcvdx(n,2,nec,1)
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
  ! Using the chain rule:
  !
  ! 1. Compute the velocity gradients (m = 2,nq-1)
  ! 2. Compute the total energy gradient and store in the
  !    location for the temperature gradient (m = nq)
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx,face_cv)
  do n = n1,n2
    !
    r_inv = one/face_cv(n,nec)
    !
    face_dcvdx(n,1,nmx,2) = r_inv*face_dcvdx(n,1,nmx,1) &
                          - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nmx)*r_inv
    face_dcvdx(n,2,nmx,2) = r_inv*face_dcvdx(n,2,nmx,1) &
                          - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nmx)*r_inv
    !
    face_dcvdx(n,1,nmy,2) = r_inv*face_dcvdx(n,1,nmy,1) &
                          - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nmy)*r_inv
    face_dcvdx(n,2,nmy,2) = r_inv*face_dcvdx(n,2,nmy,1) &
                          - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nmy)*r_inv
    !
    face_dcvdx(n,1,nee,2) = r_inv*face_dcvdx(n,1,nee,1) &
                          - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nee)*r_inv
    face_dcvdx(n,2,nee,2) = r_inv*face_dcvdx(n,2,nee,1) &
                          - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nee)*r_inv
    !
  end do
  !
  ! The location for the temperature gradient currently contains
  ! the total energy gradient from the previous loop.
  ! Now use the velocity gradients to subtract the gradient of
  ! kinetic-energy-per-unit-mass from the total energy gradient.
  !
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx,face_cv)
  do n = n1,n2
    !
    r_inv = one/face_cv(n,nec)
    !
    face_dcvdx(n,1,nee,2) = face_dcvdx(n,1,nee,2) &
                          - face_dcvdx(n,1,nmx,2)*face_cv(n,nmx)*r_inv &
                          - face_dcvdx(n,1,nmy,2)*face_cv(n,nmy)*r_inv
    !
    face_dcvdx(n,2,nee,2) = face_dcvdx(n,2,nee,2) &
                          - face_dcvdx(n,2,nmx,2)*face_cv(n,nmx)*r_inv &
                          - face_dcvdx(n,2,nmy,2)*face_cv(n,nmy)*r_inv
    !
  end do
  !
  ! Finish computing the temperature gradient by dividing through by
  ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
  do n = n1,n2
    face_dcvdx(n,1,nee,2) = face_dcvdx(n,1,nee,2) * cv_inv
    face_dcvdx(n,2,nee,2) = face_dcvdx(n,2,nee,2) * cv_inv
  end do
  !
end subroutine get_bc_wall_gradient_2d
!
!###############################################################################
!
pure subroutine get_bc_wall_gradient_nd(n1,n2,face_cv,face_dcvdx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Non-adiabatic wall
  ! #####################################################
  ! #####################################################
  !
  ! FACE_CV    : Conservative variables at the host cell flux points
  ! FACE_DCVDX : Gradients of the conservative variables at the flux pts
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: r_inv,cv_inv
  !
continue
  !
  nr = size(face_dcvdx,dim=2)
  !
  ! Precompute the inverse of specific heat at constant volume
  !
  cv_inv = (gam - one) / rgasref_nd
  !
  ! Copy over the interior density gradient to the boundary flux points
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
    do n = n1,n2
      face_dcvdx(n,1,nec,2) = face_dcvdx(n,1,nec,1)
      face_dcvdx(n,2,nec,2) = face_dcvdx(n,2,nec,1)
      face_dcvdx(n,3,nec,2) = face_dcvdx(n,3,nec,1)
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
    do n = n1,n2
      face_dcvdx(n,1,nec,2) = face_dcvdx(n,1,nec,1)
      face_dcvdx(n,2,nec,2) = face_dcvdx(n,2,nec,1)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
    do n = n1,n2
      face_dcvdx(n,1,nec,2) = face_dcvdx(n,1,nec,1)
    end do
  end if
  !
  ! We will also copy the other conservative variable gradients from
  ! the interior cell to the boundary flux points. However, since
  ! we treat the viscous fluxes for wall boundaries differently from
  ! the flow-through boundaries, we need to ensure that the wall
  ! gradients are in the right form for the wall viscous flux function.
  ! This means pre-converting the conservative variable gradients to
  ! primitive variable gradients.
  !
  ! Using the chain rule:
  !
  ! 1. Compute the velocity gradients (m = 2,nq-1)
  ! 2. Compute the total energy gradient and store in the
  !    location for the temperature gradient (m = nq)
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(face_dcvdx,face_cv)
    do n = n1,n2
      !
      r_inv = one/face_cv(n,nec)
      !
      face_dcvdx(n,1,nmx,2) = r_inv*face_dcvdx(n,1,nmx,1) &
                            - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nmx)*r_inv
      face_dcvdx(n,2,nmx,2) = r_inv*face_dcvdx(n,2,nmx,1) &
                            - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nmx)*r_inv
      face_dcvdx(n,3,nmx,2) = r_inv*face_dcvdx(n,3,nmx,1) &
                            - r_inv*face_dcvdx(n,3,nec,1)*face_cv(n,nmx)*r_inv
      !
      face_dcvdx(n,1,nmy,2) = r_inv*face_dcvdx(n,1,nmy,1) &
                            - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nmy)*r_inv
      face_dcvdx(n,2,nmy,2) = r_inv*face_dcvdx(n,2,nmy,1) &
                            - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nmy)*r_inv
      face_dcvdx(n,3,nmy,2) = r_inv*face_dcvdx(n,3,nmy,1) &
                            - r_inv*face_dcvdx(n,3,nec,1)*face_cv(n,nmy)*r_inv
      !
      face_dcvdx(n,1,nmz,2) = r_inv*face_dcvdx(n,1,nmz,1) &
                            - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nmz)*r_inv
      face_dcvdx(n,2,nmz,2) = r_inv*face_dcvdx(n,2,nmz,1) &
                            - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nmz)*r_inv
      face_dcvdx(n,3,nmz,2) = r_inv*face_dcvdx(n,3,nmz,1) &
                            - r_inv*face_dcvdx(n,3,nec,1)*face_cv(n,nmz)*r_inv
      !
      face_dcvdx(n,1,nee,2) = r_inv*face_dcvdx(n,1,nee,1) &
                            - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nee)*r_inv
      face_dcvdx(n,2,nee,2) = r_inv*face_dcvdx(n,2,nee,1) &
                            - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nee)*r_inv
      face_dcvdx(n,3,nee,2) = r_inv*face_dcvdx(n,3,nee,1) &
                            - r_inv*face_dcvdx(n,3,nec,1)*face_cv(n,nee)*r_inv
      !
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(face_dcvdx,face_cv)
    do n = n1,n2
      !
      r_inv = one/face_cv(n,nec)
      !
      face_dcvdx(n,1,nmx,2) = r_inv*face_dcvdx(n,1,nmx,1) &
                            - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nmx)*r_inv
      face_dcvdx(n,2,nmx,2) = r_inv*face_dcvdx(n,2,nmx,1) &
                            - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nmx)*r_inv
      !
      face_dcvdx(n,1,nmy,2) = r_inv*face_dcvdx(n,1,nmy,1) &
                            - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nmy)*r_inv
      face_dcvdx(n,2,nmy,2) = r_inv*face_dcvdx(n,2,nmy,1) &
                            - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nmy)*r_inv
      !
      face_dcvdx(n,1,nee,2) = r_inv*face_dcvdx(n,1,nee,1) &
                            - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nee)*r_inv
      face_dcvdx(n,2,nee,2) = r_inv*face_dcvdx(n,2,nee,1) &
                            - r_inv*face_dcvdx(n,2,nec,1)*face_cv(n,nee)*r_inv
      !
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(face_dcvdx,face_cv)
    do n = n1,n2
      !
      r_inv = one/face_cv(n,nec)
      !
      face_dcvdx(n,1,nmx,2) = r_inv*face_dcvdx(n,1,nmx,1) &
                            - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nmx)*r_inv
      !
      face_dcvdx(n,1,nee,2) = r_inv*face_dcvdx(n,1,nee,1) &
                            - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nee)*r_inv
      !
    end do
  end if
  !
  ! The location for the temperature gradient currently contains
  ! the total energy gradient from the previous loop.
  ! Now use the velocity gradients to subtract the gradient of
  ! kinetic-energy-per-unit-mass from the total energy gradient.
  !
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(face_dcvdx,face_cv)
    do n = n1,n2
      !
      r_inv = one/face_cv(n,nec)
      !
      face_dcvdx(n,1,nee,2) = face_dcvdx(n,1,nee,2) &
                            - face_dcvdx(n,1,nmx,2)*face_cv(n,nmx)*r_inv &
                            - face_dcvdx(n,1,nmy,2)*face_cv(n,nmy)*r_inv &
                            - face_dcvdx(n,1,nmz,2)*face_cv(n,nmz)*r_inv
      !
      face_dcvdx(n,2,nee,2) = face_dcvdx(n,2,nee,2) &
                            - face_dcvdx(n,2,nmx,2)*face_cv(n,nmx)*r_inv &
                            - face_dcvdx(n,2,nmy,2)*face_cv(n,nmy)*r_inv &
                            - face_dcvdx(n,2,nmz,2)*face_cv(n,nmz)*r_inv
      !
      face_dcvdx(n,3,nee,2) = face_dcvdx(n,3,nee,2) &
                            - face_dcvdx(n,3,nmx,2)*face_cv(n,nmx)*r_inv &
                            - face_dcvdx(n,3,nmy,2)*face_cv(n,nmy)*r_inv &
                            - face_dcvdx(n,3,nmz,2)*face_cv(n,nmz)*r_inv
      !
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(face_dcvdx,face_cv)
    do n = n1,n2
      !
      r_inv = one/face_cv(n,nec)
      !
      face_dcvdx(n,1,nee,2) = face_dcvdx(n,1,nee,2) &
                            - face_dcvdx(n,1,nmx,2)*face_cv(n,nmx)*r_inv &
                            - face_dcvdx(n,1,nmy,2)*face_cv(n,nmy)*r_inv
      !
      face_dcvdx(n,2,nee,2) = face_dcvdx(n,2,nee,2) &
                            - face_dcvdx(n,2,nmx,2)*face_cv(n,nmx)*r_inv &
                            - face_dcvdx(n,2,nmy,2)*face_cv(n,nmy)*r_inv
      !
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(face_dcvdx,face_cv)
    do n = n1,n2
      r_inv = one/face_cv(n,nec)
      face_dcvdx(n,1,nee,2) = face_dcvdx(n,1,nee,2) &
                            - face_dcvdx(n,1,nmx,2)*face_cv(n,nmx)*r_inv
    end do
  end if
  !
  ! Finish computing the temperature gradient by dividing through by
  ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
    do n = n1,n2
      face_dcvdx(n,1,nee,2) = face_dcvdx(n,1,nee,2) * cv_inv
      face_dcvdx(n,2,nee,2) = face_dcvdx(n,2,nee,2) * cv_inv
      face_dcvdx(n,3,nee,2) = face_dcvdx(n,3,nee,2) * cv_inv
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
    do n = n1,n2
      face_dcvdx(n,1,nee,2) = face_dcvdx(n,1,nee,2) * cv_inv
      face_dcvdx(n,2,nee,2) = face_dcvdx(n,2,nee,2) * cv_inv
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
    do n = n1,n2
      face_dcvdx(n,1,nee,2) = face_dcvdx(n,1,nee,2) * cv_inv
    end do
  end if
  !
end subroutine get_bc_wall_gradient_nd
!
!###############################################################################
!
pure subroutine get_bc_slip_wall_gradient(n1,n2,facenrm,face_cv,face_dcvdx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Slip wall / symmetry boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FACENRM    : Face outward normal vector at the flux points
  ! FACE_CV    : Conservative variables at the flux points
  ! FACE_DCVDX : Gradients of the conservative variables at the flux pts
  !
  !.. Use Statements ..
  use gpu_mod, only : grad_cv_to_grad_pv_gpu
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr,np
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: dpvdx(:,:,:)
  !
continue
  !
  np = size(face_dcvdx,dim=1)
  nr = size(face_dcvdx,dim=2)
  !
  allocate ( dpvdx(1:np,1:nr,1:nq) )
  !$ACC DATA CREATE(dpvdx)
  !
  ! Convert the interior host cell conservative gradients to primitive
  !
  call grad_cv_to_grad_pv_gpu(n1,n2,face_cv,face_dcvdx(:,:,:,1),dpvdx)
  !
  ! Get the primitive variable gradients at the ghost cell flux points
  !
  if (nr == 3) then
    call get_bc_slip_wall_gradient_3d(n1,n2,facenrm,dpvdx,face_dcvdx(:,:,:,2))
  else if (nr == 2) then
    call get_bc_slip_wall_gradient_2d(n1,n2,facenrm,dpvdx,face_dcvdx(:,:,:,2))
  else
    call get_bc_slip_wall_gradient_nd(n1,n2,facenrm,dpvdx,face_dcvdx(:,:,:,2))
  end if
  !
  !$ACC END DATA
  deallocate ( dpvdx )
  !
end subroutine get_bc_slip_wall_gradient
!
!###############################################################################
!
pure subroutine get_bc_slip_wall_gradient_3d(n1,n2,facenrm,dpvdx,face_dcvdx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Slip wall / symmetry boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FACENRM    : Face outward normal vector at the flux points
  ! DPVDX      : Gradients of the primitive variables at the host  cell flux pts
  ! FACE_DCVDX : Gradients of the primitive variables at the ghost cell flux pts
  !
  !.. Use Statements ..
  use ovar, only : walls_are_exact
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: dpvdx(:,:,:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: drdn,dudn,dvdn,dwdn,dtdn
  real(wp) :: reflect_const
  !
continue
  !
  reflect_const = merge(one,two,walls_are_exact)
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_dcvdx,dpvdx)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    rnz = facenrm(3,n)
    ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
    unx = rnx/ds
    uny = rny/ds
    unz = rnz/ds
    !
    drdn = dpvdx(n,1,nec)*unx + dpvdx(n,2,nec)*uny + dpvdx(n,3,nec)*unz
    dudn = dpvdx(n,1,nmx)*unx + dpvdx(n,2,nmx)*uny + dpvdx(n,3,nmx)*unz
    dvdn = dpvdx(n,1,nmy)*unx + dpvdx(n,2,nmy)*uny + dpvdx(n,3,nmy)*unz
    dwdn = dpvdx(n,1,nmz)*unx + dpvdx(n,2,nmz)*uny + dpvdx(n,3,nmz)*unz
    dtdn = dpvdx(n,1,nee)*unx + dpvdx(n,2,nee)*uny + dpvdx(n,3,nee)*unz
    !
    face_dcvdx(n,1,nec) = dpvdx(n,1,nec) - reflect_const*unx*drdn
    face_dcvdx(n,2,nec) = dpvdx(n,2,nec) - reflect_const*uny*drdn
    face_dcvdx(n,3,nec) = dpvdx(n,3,nec) - reflect_const*unz*drdn
    !
    face_dcvdx(n,1,nmx) = dpvdx(n,1,nmx) - reflect_const*unx*dudn
    face_dcvdx(n,2,nmx) = dpvdx(n,2,nmx) - reflect_const*uny*dudn
    face_dcvdx(n,3,nmx) = dpvdx(n,3,nmx) - reflect_const*unz*dudn
    !
    face_dcvdx(n,1,nmy) = dpvdx(n,1,nmy) - reflect_const*unx*dvdn
    face_dcvdx(n,2,nmy) = dpvdx(n,2,nmy) - reflect_const*uny*dvdn
    face_dcvdx(n,3,nmy) = dpvdx(n,3,nmy) - reflect_const*unz*dvdn
    !
    face_dcvdx(n,1,nmz) = dpvdx(n,1,nmz) - reflect_const*unx*dwdn
    face_dcvdx(n,2,nmz) = dpvdx(n,2,nmz) - reflect_const*uny*dwdn
    face_dcvdx(n,3,nmz) = dpvdx(n,3,nmz) - reflect_const*unz*dwdn
    !
    face_dcvdx(n,1,nee) = dpvdx(n,1,nee) - reflect_const*unx*dtdn
    face_dcvdx(n,2,nee) = dpvdx(n,2,nee) - reflect_const*uny*dtdn
    face_dcvdx(n,3,nee) = dpvdx(n,3,nee) - reflect_const*unz*dtdn
    !
  end do
  !
end subroutine get_bc_slip_wall_gradient_3d
!
!###############################################################################
!
pure subroutine get_bc_slip_wall_gradient_2d(n1,n2,facenrm,dpvdx,face_dcvdx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Slip wall / symmetry boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FACENRM    : Face outward normal vector at the flux points
  ! DPVDX      : Gradients of the primitive variables at the host  cell flux pts
  ! FACE_DCVDX : Gradients of the primitive variables at the ghost cell flux pts
  !
  !.. Use Statements ..
  use ovar, only : walls_are_exact
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: dpvdx(:,:,:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,rny,unx,uny
  real(wp) :: drdn,dudn,dvdn,dtdn
  real(wp) :: reflect_const
  !
continue
  !
  reflect_const = merge(one,two,walls_are_exact)
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_dcvdx,dpvdx)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds = sqrt( rnx*rnx + rny*rny )
    unx = rnx/ds
    uny = rny/ds
    !
    drdn = dpvdx(n,1,nec)*unx + dpvdx(n,2,nec)*uny
    dudn = dpvdx(n,1,nmx)*unx + dpvdx(n,2,nmx)*uny
    dvdn = dpvdx(n,1,nmy)*unx + dpvdx(n,2,nmy)*uny
    dtdn = dpvdx(n,1,nee)*unx + dpvdx(n,2,nee)*uny
    !
    face_dcvdx(n,1,nec) = dpvdx(n,1,nec) - reflect_const*unx*drdn
    face_dcvdx(n,2,nec) = dpvdx(n,2,nec) - reflect_const*uny*drdn
    !
    face_dcvdx(n,1,nmx) = dpvdx(n,1,nmx) - reflect_const*unx*dudn
    face_dcvdx(n,2,nmx) = dpvdx(n,2,nmx) - reflect_const*uny*dudn
    !
    face_dcvdx(n,1,nmy) = dpvdx(n,1,nmy) - reflect_const*unx*dvdn
    face_dcvdx(n,2,nmy) = dpvdx(n,2,nmy) - reflect_const*uny*dvdn
    !
    face_dcvdx(n,1,nee) = dpvdx(n,1,nee) - reflect_const*unx*dtdn
    face_dcvdx(n,2,nee) = dpvdx(n,2,nee) - reflect_const*uny*dtdn
    !
  end do
  !
end subroutine get_bc_slip_wall_gradient_2d
!
!###############################################################################
!
pure subroutine get_bc_slip_wall_gradient_nd(n1,n2,facenrm,dpvdx,face_dcvdx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Slip wall / symmetry boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! FACENRM    : Face outward normal vector at the flux points
  ! DPVDX      : Gradients of the primitive variables at the host  cell flux pts
  ! FACE_DCVDX : Gradients of the primitive variables at the ghost cell flux pts
  !
  !.. Use Statements ..
  use ovar, only : walls_are_exact
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: dpvdx(:,:,:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: drdn,dudn,dvdn,dwdn,dtdn
  real(wp) :: reflect_const
  !
continue
  !
  reflect_const = merge(one,two,walls_are_exact)
  !
  nr = size(face_dcvdx,dim=2)
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_dcvdx,dpvdx)
    do n = n1,n2
      !
      rnx = facenrm(1,n)
      rny = facenrm(2,n)
      rnz = facenrm(3,n)
      ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
      unx = rnx/ds
      uny = rny/ds
      unz = rnz/ds
      !
      drdn = dpvdx(n,1,nec)*unx + dpvdx(n,2,nec)*uny + dpvdx(n,3,nec)*unz
      dudn = dpvdx(n,1,nmx)*unx + dpvdx(n,2,nmx)*uny + dpvdx(n,3,nmx)*unz
      dvdn = dpvdx(n,1,nmy)*unx + dpvdx(n,2,nmy)*uny + dpvdx(n,3,nmy)*unz
      dwdn = dpvdx(n,1,nmz)*unx + dpvdx(n,2,nmz)*uny + dpvdx(n,3,nmz)*unz
      dtdn = dpvdx(n,1,nee)*unx + dpvdx(n,2,nee)*uny + dpvdx(n,3,nee)*unz
      !
      face_dcvdx(n,1,nec) = dpvdx(n,1,nec) - reflect_const*unx*drdn
      face_dcvdx(n,2,nec) = dpvdx(n,2,nec) - reflect_const*uny*drdn
      face_dcvdx(n,3,nec) = dpvdx(n,3,nec) - reflect_const*unz*drdn
      !
      face_dcvdx(n,1,nmx) = dpvdx(n,1,nmx) - reflect_const*unx*dudn
      face_dcvdx(n,2,nmx) = dpvdx(n,2,nmx) - reflect_const*uny*dudn
      face_dcvdx(n,3,nmx) = dpvdx(n,3,nmx) - reflect_const*unz*dudn
      !
      face_dcvdx(n,1,nmy) = dpvdx(n,1,nmy) - reflect_const*unx*dvdn
      face_dcvdx(n,2,nmy) = dpvdx(n,2,nmy) - reflect_const*uny*dvdn
      face_dcvdx(n,3,nmy) = dpvdx(n,3,nmy) - reflect_const*unz*dvdn
      !
      face_dcvdx(n,1,nmz) = dpvdx(n,1,nmz) - reflect_const*unx*dwdn
      face_dcvdx(n,2,nmz) = dpvdx(n,2,nmz) - reflect_const*uny*dwdn
      face_dcvdx(n,3,nmz) = dpvdx(n,3,nmz) - reflect_const*unz*dwdn
      !
      face_dcvdx(n,1,nee) = dpvdx(n,1,nee) - reflect_const*unx*dtdn
      face_dcvdx(n,2,nee) = dpvdx(n,2,nee) - reflect_const*uny*dtdn
      face_dcvdx(n,3,nee) = dpvdx(n,3,nee) - reflect_const*unz*dtdn
      !
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_dcvdx,dpvdx)
    do n = n1,n2
      !
      rnx = facenrm(1,n)
      rny = facenrm(2,n)
      ds = sqrt( rnx*rnx + rny*rny )
      unx = rnx/ds
      uny = rny/ds
      !
      drdn = dpvdx(n,1,nec)*unx + dpvdx(n,2,nec)*uny
      dudn = dpvdx(n,1,nmx)*unx + dpvdx(n,2,nmx)*uny
      dvdn = dpvdx(n,1,nmy)*unx + dpvdx(n,2,nmy)*uny
      dtdn = dpvdx(n,1,nee)*unx + dpvdx(n,2,nee)*uny
      !
      face_dcvdx(n,1,nec) = dpvdx(n,1,nec) - reflect_const*unx*drdn
      face_dcvdx(n,2,nec) = dpvdx(n,2,nec) - reflect_const*uny*drdn
      !
      face_dcvdx(n,1,nmx) = dpvdx(n,1,nmx) - reflect_const*unx*dudn
      face_dcvdx(n,2,nmx) = dpvdx(n,2,nmx) - reflect_const*uny*dudn
      !
      face_dcvdx(n,1,nmy) = dpvdx(n,1,nmy) - reflect_const*unx*dvdn
      face_dcvdx(n,2,nmy) = dpvdx(n,2,nmy) - reflect_const*uny*dvdn
      !
      face_dcvdx(n,1,nee) = dpvdx(n,1,nee) - reflect_const*unx*dtdn
      face_dcvdx(n,2,nee) = dpvdx(n,2,nee) - reflect_const*uny*dtdn
      !
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(face_dcvdx,dpvdx)
    do n = n1,n2
      face_dcvdx(n,1,nec) = (one - reflect_const)*dpvdx(n,1,nec)
      face_dcvdx(n,1,nmx) = (one - reflect_const)*dpvdx(n,1,nmx)
      face_dcvdx(n,1,nee) = (one - reflect_const)*dpvdx(n,1,nee)
    end do
  end if
  !
end subroutine get_bc_slip_wall_gradient_nd
!
!###############################################################################
!
subroutine output_bface
  !
  !.. Use Statements ..
  use geovar, only : bface,face
  !
  !.. Local Scalars ..
  integer :: i,n,maxint,myio,ierr
  character(len=100) :: write_format
  character(len=100) :: file_format
  character(len=100) :: my_file
  !
  !.. Local Arrays ..
  integer :: intlen(1:size(bface,dim=1)+3)
  integer :: temp_bface(1:size(bface,dim=1)+2,1:size(bface,dim=2))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "output_bface"
  !
continue
  !
  myio = 1456 + mypnum
  !
  do n = 1,size(bface,dim=2)
    do i = 1,size(bface,dim=1)
      temp_bface(i,n) = bface(i,n)
    end do
    temp_bface(size(bface,dim=1)+1,n) = face(n)%beg_fp
    temp_bface(size(bface,dim=1)+2,n) = face(n)%end_fp
  end do
  !
  ! Get the maximum length of the string containing the name of the BC types
  !
  intlen(1) = maxval( len_trim( bc_integer_to_string(temp_bface(1,:)) ) )
  !
  do i = 1,size(intlen)-1
    n = 1
    maxint = maxval(temp_bface(i,:))
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
  do n = 1,size(temp_bface,dim=2)
    write (myio,fmt=trim(write_format)) &
                               trim(bc_integer_to_string(temp_bface(1,n))), &
                                   (temp_bface(i,n),i=1,size(temp_bface,dim=1))
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
#include "Functions/debug_timer.f90"
!
!###############################################################################
!
end module module_bc
