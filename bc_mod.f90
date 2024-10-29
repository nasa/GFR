module bc_mod
  !
#include <mpi_defs.h>
!#define PROFILE_ON
  !
  ! Module to contain the procedures and variables related
  ! to the evaluating the boundary conditions for the flow.
  !
  !.. Use Statements ..
  use module_kind_types
  !
  use eqn_idx_mod, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq
  use eqn_idx_mod, only : ncs
  use timer_mod, only : solver_timer
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
  use metrics_mod,   only : facenrm
  use flowvar,       only : faceusp
  use flowvar,       only : facexyz
  use flowvar,       only : facegpsp
  use geovar,        only : bc_grps,n_bnd_flxpts
  use geovar,        only : output_bface
  use ovar,          only : bc_in,itcur,aref
  use ovar,          only : current_rk_stage
  use ovar,          only : moving_ref_frame
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  use synthetic_eddy_mod, only : sem
  !
  use bc_fixed_mod,            only : get_bc_fixed_solution
  use bc_riemann_mod,          only : get_bc_characteristic_solution
  use bc_gen_in_mod,           only : get_bc_generic_inflow_solution
  use bc_gen_out_mod,          only : get_bc_generic_outflow_solution
  use bc_gen_free_mod,         only : get_bc_generic_freeflow_solution
  use bc_sub_in_mod,           only : get_bc_sub_inflow_solution
  use bc_sub_out_mod,          only : get_bc_sub_outflow_solution
  use bc_sup_out_mod,          only : get_bc_sup_outflow_solution
  use bc_wall_mod,             only : get_bc_wall_solution
  use bc_periodic_mod,         only : get_bc_periodic_solution
  use bc_mms_mod,              only : get_bc_mms_dirichlet_solution
  !
  !.. Optional Arguments ..
  integer, optional, intent(out) :: exit_error
  !
  !.. Local Scalars ..
  integer :: grp_idx,n1,n2
  integer :: bc_idx,bc_type
  integer :: newton_failures
  !
  !.. Local Saved Scalars ..
  logical(lk), save :: first_output = true
  integer, save :: idx = -1
  !
  !.. Local Arrays ..
  real(wp) :: rotation_rate(1:3)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: pv_in(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_faceusp_bc"
  !
continue
  !
  call solver_timer(entering_procedure,pname,idx)
  !
  call debug_timer(entering_procedure,pname)
  !
  newton_failures = 0
  if (present(exit_error)) exit_error = 0
  !
  if (first_output) then
    call output_bface(strng="bc_")
  end if
  !
  ! Default the rotation rate to zero
  !
  rotation_rate(:) = zero
  !
  ! First, update the thermodynamic and transport gas properties for the
  ! interior flux points of the boundary faces
  !
  n1 = 1
  n2 = n_bnd_flxpts
  !
  call update_thermo_properties(n1,n2,faceusp(:,:,1),facegpsp(:,:,1), &
                                      loc="face",caller_name=pname)
  call update_transport_properties(n1,n2,faceusp(:,:,1),facegpsp(:,:,1))
  !
  ! Call the SEM turbulence inflow generator before looping
  ! through bc_grps to set the boundary conditions
  !
  if (current_rk_stage == 1) then
    call sem(faceusp,facegpsp)
  end if
  !
  ! Loop over the boundary face groups and compute the conservative variables
  ! at the boundary flux points based on the face boundary condition.
  !
  loop_over_boundary_groups: &
  do grp_idx = 1,size(bc_grps)
    !
    n1 = bc_grps(grp_idx)%beg_fp ! beginning face point for this BC group
    n2 = bc_grps(grp_idx)%end_fp ! ending    face point for this BC group
    !
    bc_type = bc_grps(grp_idx)%bc_type ! BC for this group
    bc_idx = bc_grps(grp_idx)%bc_idx   ! index within bc_in for this BC group
    !
!#define DEBUG_BC_GROUPS
#ifdef DEBUG_BC_GROUPS
    if (first_output) then
      call output_to_debug_bc_groups(grp_idx,bc_idx,bc_type,n1,n2, &
                                     facenrm,faceusp)
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
        call get_bc_fixed_solution(n1,n2,faceusp(:,:,2),bc_in(0)%cv)
        !
      case (bc_fixed)
        !
        ! Fixed boundary condition
        !
        call get_bc_fixed_solution(n1,n2,faceusp(:,:,2),bc_in(bc_idx)%cv)
        !
      case (bc_characteristic)
        !
        ! Characteristic boundary condition
        ! NOTE: The characteristic boundary condition is not really valid
        !       for gas mixtures because the Riemann invariants assume a
        !       frozen gamma.
        ! NOTE: We need to output a WARNING of some kind before the
        !       simulation starts to warn the user that the characteristic
        !       boundary condition is not really valid for gas mixtures.
        !
        call get_bc_characteristic_solution(n1,n2,facenrm, &
                                            faceusp,facegpsp, &
                                            bc_in(bc_idx)%pv, &
                                            bc_in(0)%pv, &
                                            bc_in(bc_idx)%relax)
        !
      case (bc_generic_inflow)
        !
        ! Generic inflow boundary condition
        !
        call get_bc_generic_inflow_solution(n1,n2,facenrm, &
                                            faceusp,facegpsp, &
                                            bc_in(bc_idx)%pv, &
                                            bc_in(bc_idx)%gp, &
                                            bc_in(0)%pv, &
                                            bc_in(0)%gp, &
                                            bc_in(bc_idx)%relax, &
                                            newton_failures)
        !
      case (bc_generic_outflow)
        !
        ! Generic outflow boundary condition
        !
        call get_bc_generic_outflow_solution(n1,n2,facenrm, &
                                             faceusp,facegpsp, &
                                             bc_in(bc_idx)%pv, &
                                             bc_in(bc_idx)%gp, &
                                             bc_in(bc_idx)%relax)
        !
      case (bc_generic_freeflow)
        !
        ! Generic freeflow boundary condition
        !
        call get_bc_generic_freeflow_solution(n1,n2,facenrm, &
                                              faceusp,facegpsp, &
                                              bc_in(bc_idx)%pv, &
                                              bc_in(bc_idx)%gp, &
                                              bc_in(0)%pv, &
                                              bc_in(0)%gp, &
                                              bc_in(bc_idx)%relax, &
                                              newton_failures)
        !
      case (bc_sub_inflow)
        !
        ! Subsonic inflow boundary condition
        !
        call get_bc_sub_inflow_solution(n1,n2,facenrm, &
                                        faceusp,facegpsp, &
                                        bc_in(bc_idx)%cv, &
                                        bc_in(bc_idx)%pv, &
                                        bc_in(bc_idx)%gp, &
                                        bc_in(0)%pv, &
                                        bc_in(0)%gp, &
                                        bc_in(bc_idx)%relax, &
                                        newton_failures)
        !
      case (bc_sub_outflow)
        !
        ! Subsonic outflow boundary condition
        !
        call get_bc_sub_outflow_solution(n1,n2,facenrm, &
                                         faceusp,facegpsp, &
                                         bc_in(bc_idx)%pv, &
                                         bc_in(bc_idx)%gp, &
                                         bc_in(bc_idx)%relax)
        !
      case (bc_sup_inflow)
        !
        ! Supersonic inflow boundary condition
        !
        call get_bc_fixed_solution(n1,n2,faceusp(:,:,2),bc_in(bc_idx)%cv)
        !
      case (bc_sup_outflow)
        !
        ! Supersonic outflow boundary condition
        !
        call get_bc_sup_outflow_solution(n1,n2,faceusp)
        !
      case (bc_slip_wall,bc_symmetry,bc_adiabatic_wall,bc_isothermal_wall)
        !
        ! All wall-type BCs
        !  - Slip wall
        !  - Symmetry BC
        !  - No-slip adiabatic  wall, with/without imposed wall velocity
        !  - No-slip isothermal wall, with/without imposed wall velocity
        !
        pv_in = bc_in(bc_idx)%pv(:) ! F2003 AUTO-REALLOC
        !
        ! If bc_idx==0, that means this wall-type boundary condition doesnt
        ! specify any flow conditions (e.g., imposed wall velocity) and just
        ! uses interior values and a zero-velocity wall for a no-slip wall.
        ! That means we need to zero out the velocity in pv_in when bc_idx==0,
        ! because otherwise the no-slip wall BCs will impose the free stream
        ! velocity at the wall, since that is what is saved in
        ! bc_in(0)%pv(nmb:nme).
        !
        if (bc_idx == 0) pv_in(nmb:nme) = zero
        !
        ! If this BC has a specified wall rotation rate, multiply it by the
        ! moving reference frame spinup ratio to make sure any BCs specified
        ! in the MRF use a consistent rotation rate with the MRF.
        !
        if (any(abs(bc_in(bc_idx)%rotation_rate(:)) >= rpm_to_mps/aref)) then
          rotation_rate(:) = bc_in(bc_idx)%rotation_rate(:) * &
                             moving_ref_frame%spinup%get_spinup_ratio()
        end if
        !
        call get_bc_wall_solution(n1,n2,bc_type,facexyz,facenrm, &
                                  faceusp,facegpsp, &
                                  bc_in(bc_idx)%relax,pv_in, &
                                  bc_in(bc_idx)%wall_temperature, &
                                  rotation_rate)
        !
      case (bc_periodic)
        !
        ! Periodic boundary condition
        !
        call get_bc_periodic_solution(n1,n2,faceusp, &
                                      bc_grps(grp_idx)%periodic_list, &
                                      bc_in(bc_idx)%pbc_rotation%mat)
        !
      case (bc_mms_dirichlet)
        !
        ! Dirichlet boundary condition for MMS
        !
        if (ncs > 0) then
          write (error_message,3) &
                     trim(adjustl(bc_integer_to_string(bc_characteristic)))
          call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
        end if
        !
        call get_bc_mms_dirichlet_solution(n1,n2,facenrm,facexyz,faceusp)
        !
      case (bc_cpu_bnd)
        !
        ! We dont need to do anything here for communication boundaries
        !
      case (bc_sem)
        !
        ! We already set the SEM boundary conditions when calling subroutine sem
        ! just before the current loop, so we dont need to do anything here.
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
  end do loop_over_boundary_groups
  !
  ! NOTE: We dont need to update the thermo/transport properties for the
  !       boundary states [i.e., face_gp(:,:,2)] because that is done just
  !       before computing the inviscid common fluxes for the boundary faces
  !       in subroutine get_invs_common_fluxes_bc, which is the first place
  !       these values would be used.
  !
  ! Sum the total number of Newton failures across all MPI processes
  !
  if (ncpu > 1) then
    call mpi_allreduce(MPI_IN_PLACE,newton_failures,1_int_mpi, &
                       mpi_inttyp,MPI_SUM,MPI_COMM_WORLD,mpierr)
  end if
  !
  ! Copy the number of Newton failures into exit error if the optional
  ! argument is present; otherwise, write the number of Newton failures
  ! to an error message and stop execution.
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
  ! Set first_output to false so that any BC debug output is
  ! only written the first time this subroutine is called.
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
  call solver_timer(leaving_procedure,pname,idx)
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
  3 format (a," BCs ARE NOT SUPPORTED FOR NON-CALORICALLY PERFECT GASES!!!")
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
  use ovar,         only : bc_in
  use metrics_mod,  only : facenrm
  use flowvar,      only : faceusp
  use flowvar,      only : facegpsp
  use flowvar,      only : facedusp
  use flowvar,      only : facegpsp
  !
  use bc_wall_mod,             only : get_bc_wall_gradient
  use bc_periodic_mod,         only : get_bc_periodic_gradient
  use bc_flowthrough_grad_mod, only : get_bc_flowthrough_gradient
  !
  !.. Local Scalars ..
  integer :: grp_idx,n1,n2,bc_type,bc_idx
  !
  !. Local Saved Scalars ..
  integer, save :: idx = -1
  !
  !.. Local Parameters ..
  integer, parameter :: local_bc_ignore(1:3) = [ bc_cpu_bnd, &
                                                 bc_unknown, &
                                                 not_a_bc ]
  character(len=*), parameter :: pname = "get_facedusp_bc"
  !
continue
  !
  call solver_timer(entering_procedure,pname,idx)
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Loop over the boundary face groups and compute the conservative variables
  ! at the boundary flux points based on the face boundary condition.
  !
  do grp_idx = 1,size(bc_grps)
    !
    n1 = bc_grps(grp_idx)%beg_fp ! beginning face point for this BC group
    n2 = bc_grps(grp_idx)%end_fp ! ending    face point for this BC group
    !
    bc_type = bc_grps(grp_idx)%bc_type ! BC for this group
    bc_idx = bc_grps(grp_idx)%bc_idx   ! index within bc_in for this BC group
    !
    ! Compute the solution at the boundary flux points
    ! based on the specified boundary condition
    !
    if (any(bc_type == bc_walls)) then
      !
      ! All wall-type boundary conditions
      !  - Slip wall
      !  - Symmetry BC
      !  - No-slip adiabatic  wall
      !  - No-slip isothermal wall
      !
      call get_bc_wall_gradient(n1,n2,bc_type,facenrm, &
                                faceusp(:,:,1), &
                                facegpsp(:,:,1), &
                                facedusp)
      !
    else if (bc_type == bc_periodic) then
      !
      ! Periodic boundary condition
      !
      call get_bc_periodic_gradient(n1,n2,facedusp, &
                                    bc_grps(grp_idx)%periodic_list, &
                                    bc_in(bc_idx)%pbc_rotation%mat)
      !
    else if (bc_type == bc_sem) then
      !
      ! For now, use the general flow-through boundary conditions for
      ! getting the gradient boundary conditions for SEM boundaries.
      !
      call get_bc_flowthrough_gradient(n1,n2,facedusp)
      !
    else if (all(bc_type /= local_bc_ignore)) then
      !
      ! General flow-through boundary conditions
      !
      call get_bc_flowthrough_gradient(n1,n2,facedusp)
      !
    end if
    !
  end do
  !
 !call exchange_bc_conflicts_dusp(facedusp)
  !
  call debug_timer(leaving_procedure,pname)
  !
  call solver_timer(leaving_procedure,pname,idx)
  !
end subroutine get_facedusp_bc
!
!###############################################################################
!
subroutine output_to_debug_bc_groups(grp_idx,bc_idx,bc_type,beg_fp,end_fp, &
                                     facenrm,faceusp)
  !
  !.. Use Statements ..
  use geovar, only : n_solpts,n_flxpts,n_wall_flxpts,n_bnd_flxpts
  use ovar,   only : bc_in
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: grp_idx
  integer,  intent(in) :: bc_idx
  integer,  intent(in) :: bc_type
  integer,  intent(in) :: beg_fp
  integer,  intent(in) :: end_fp
  real(wp), intent(in) :: facenrm(:,:)
  real(wp), intent(in) :: faceusp(:,:,:)
  !
  !.. Local Scalars ..
  integer :: m
  !
  !.. Local Saved Scalars ..
  logical(lk), save :: very_first_output = true
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "output_to_debug_bc_groups"
  !
continue
  !
  if (very_first_output) then
    write (mypnum+500,1) n_wall_flxpts,n_bnd_flxpts,n_flxpts,n_solpts
    flush (mypnum+500)
  end if
  !
  very_first_output = fals
  !
  write (mypnum+500,2) &
        grp_idx, beg_fp, &
        grp_idx, end_fp, &
        grp_idx, bc_type, trim(adjustl(bc_integer_to_string(bc_type))), &
        grp_idx, bc_idx, &
        bc_idx, bc_idx, lbound(bc_in(bc_idx)%cv,dim=1), &
                        ubound(bc_in(bc_idx)%cv,dim=1), &
        bc_idx, bc_idx, lbound(bc_in(bc_idx)%pv,dim=1), &
                        ubound(bc_in(bc_idx)%pv,dim=1), &
        lbound(facenrm,dim=1), ubound(facenrm,dim=1), &
        lbound(facenrm,dim=2), ubound(facenrm,dim=2), &
        lbound(faceusp,dim=1), ubound(faceusp,dim=1), &
        lbound(faceusp,dim=2), ubound(faceusp,dim=2), &
        lbound(faceusp,dim=3), ubound(faceusp,dim=3)
  !
  if (nq < 10) then
    write (mypnum+500,3) bc_idx, (m,bc_in(bc_idx)%cv(m), m=1,nq)
    write (mypnum+500,4) bc_idx, (m,bc_in(bc_idx)%pv(m), m=1,nq)
  else
    write (mypnum+500,5) bc_idx, (m,bc_in(bc_idx)%cv(m), m=1,nq)
    write (mypnum+500,6) bc_idx, (m,bc_in(bc_idx)%pv(m), m=1,nq)
  end if
  !
  flush (mypnum+500)
  !
  ! Format Statements
  !
  1 format (/,"In get_faceusp_bc:",//, &
              "n_wall_flxpts = ",i0,/, &
              " n_bnd_flxpts = ",i0,/, &
              "     n_flxpts = ",i0,/, &
              "     n_solpts = ",i0,/)
  !
  2 format (/, &
        "bc_grps(",i0,")%beg_fp  = ",i0,/, &
        "bc_grps(",i0,")%end_fp  = ",i0,/, &
        "bc_grps(",i0,")%bc_type = ",i0," (",a,")",/, &
        "bc_grps(",i0,")%bc_idx  = ",i0,/, &
        "dim(bc_in(",i0,")%cv) = bc_in(",i0,")%cv(",i0,":",i0,")",/, &
        "dim(bc_in(",i0,")%pv) = bc_in(",i0,")%pv(",i0,":",i0,")",/, &
        "dim(facenrm) = facenrm(",i0,":",i0,",",i0,":",i0,")",/, &
        "dim(faceusp) = faceusp(",i0,":",i0,",",i0,":",i0,",",i0,":",i0,")")
  !
  3 format ("bc_in(",i2,")%cv(",i1,") = ",es14.6,/, &
                   20(10x,"cv(",i1,") = ",es14.6,:,/))
  4 format ("bc_in(",i2,")%pv(",i1,") = ",es14.6,/, &
                   20(10x,"pv(",i1,") = ",es14.6,:,/))
  5 format ("bc_in(",i2,")%cv(",i2,") = ",es14.6,/, &
                   20(10x,"cv(",i2,") = ",es14.6,:,/))
  6 format ("bc_in(",i2,")%pv(",i2,") = ",es14.6,/, &
                   20(10x,"pv(",i2,") = ",es14.6,:,/))
  !
end subroutine output_to_debug_bc_groups
!
!###############################################################################
!
end module bc_mod
