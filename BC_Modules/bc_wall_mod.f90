module bc_wall_mod
  !
  !.. Use Statements
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nee,nq
  use eqn_idx, only : ncs,npres,ntemp,ngamm
  !
  implicit none
  !
  private
  !
  public :: get_bc_wall_solution
  public :: get_bc_wall_gradient
  !
contains
!
!###############################################################################
!
#ifdef DISABLE_FOR_GPU
pure &
#endif
subroutine get_bc_wall_solution(n1,n2,bc_type,facexyz,facenrm,face_cv,face_gp, &
                                relax,pv_in,wall_temperature,rotation_speed)
  !
  ! ####################################################################
  ! ####################################################################
  ! ##  This is for all wall-type boundary conditions
  ! ##   - Slip wall / symmetry BC
  ! ##   - No-slip adiabatic  wall, with/without imposed wall velocity
  ! ##   - No-slip isothermal wall, with/without imposed wall velocity
  ! ####################################################################
  ! ####################################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! BC_TYPE : Integer identifying this boundary condition type
  ! FACEXYZ : Coordinates for the boundary flux points
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conserved    variables at the flux points
  ! FACE_GP : Gas property variables at the flux points
  ! RELAX   : Relaxation parameters for this BC
  ! PV_IN   : Primitive    variables specified for this BC
  ! WALL_TEMPERATURE : Specified wall temperature for isothermal wall BC
  ! ROTATION_SPEED : RPM for rotating non-slip wall
  !
  !.. Use Statements ..
  use ovar,         only : relaxation_t
  use gpu_mod,      only : pv2cv_gpu
  use bc_relax_mod, only : get_relax_coefficient
  !
  !.. Formal Arguments ..
  integer,            intent(in) :: n1,n2
  integer,            intent(in) :: bc_type
  real(wp),           intent(in) :: facexyz(:,:)
  real(wp),           intent(in) :: facenrm(:,:)
  real(wp),        intent(inout) :: face_cv(:,:,:)
  real(wp),        intent(inout) :: face_gp(:,:,:)
  type(relaxation_t), intent(in) :: relax
  real(wp),           intent(in) :: pv_in(:)
  real(wp),           intent(in) :: wall_temperature
  real(wp),           intent(in) :: rotation_speed
  !
  !.. Local Scalars ..
  integer  :: nr
  real(wp) :: rlx
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  ! Compute the relaxation coefficient for this BC
  !
  call get_relax_coefficient(relax,rlx)
  !
  ! Get the velocity part of the boundary condition states
  !
  if (nr == 3) then
    call get_bc_wall_vel_solution_3d(n1,n2,bc_type,facexyz,facenrm, &
                                     face_cv,rlx,pv_in,rotation_speed)
  else if (nr == 2) then
    call get_bc_wall_vel_solution_2d(n1,n2,bc_type,facexyz,facenrm, &
                                     face_cv,rlx,pv_in,rotation_speed)
  else
    ! NOTE: This branch and the called subroutine should probably just
    !       be removed entirely because there isnt really a 1D wall
    !       boundary condition.
    call get_bc_wall_vel_solution_1d(n1,n2,bc_type,facenrm,face_cv,rlx,pv_in)
  end if
  !
  ! Get the gas properties (density, pressure, temperature) for this wall BC
  !
  call get_bc_wall_gp_solution(n1,n2,bc_type,face_cv,face_gp, &
                               rlx,wall_temperature)
  !
  ! The above subroutines will return the boundary condition states as
  ! primitive variables. These need to be converted to conservative variables.
  !
  call pv2cv_gpu(n1,n2,face_cv(:,:,2),face_gp(:,:,2))
  !
  ! NOTE: We dont need to update the thermo/transport properties for the
  !       boundary states [i.e., face_gp(:,:,2)] because that is done just
  !       before computing the inviscid common fluxes for the boundary faces
  !       in subroutine get_invs_common_fluxes_bc, which is the first place
  !       these values would be used.
  !
end subroutine get_bc_wall_solution
!
!###############################################################################
!
#ifdef DISABLE_FOR_GPU
pure &
#endif
subroutine get_bc_wall_gp_solution(n1,n2,bc_type,face_cv,face_gp, &
                                   rlx,wall_temperature)
  !
  ! ####################################################################
  ! ####################################################################
  ! ##  This is for all wall-type boundary conditions
  ! ##   - Slip wall / symmetry BC
  ! ##   - No-slip adiabatic  wall, with/without imposed wall velocity
  ! ##   - No-slip isothermal wall, with/without imposed wall velocity
  ! ####################################################################
  ! ####################################################################
  ! #####################################################
  ! #####################################################
  ! ##  No-slip adiabatic wall
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! BC_TYPE : Integer identifying this boundary condition type
  ! FACE_CV : Conserved    variables at the flux points
  !           NOTE: The wall states will be saved as primitive variables
  !                 i.e., face_cv(:,:,2) = primitive variables
  ! FACE_GP : Gas property variables at the flux points
  ! RLX     : BC relaxation coefficient
  ! WALL_TEMPERATURE : specified isothermal wall temperature
  !
  !.. Use Statements ..
  use ovar,       only : rgasref_nd
  use thermo_mod, only : compute_mixture_rgas
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: bc_type
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp), intent(inout) :: face_gp(:,:,:)
  real(wp),    intent(in) :: rlx
  real(wp),    intent(in) :: wall_temperature
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: rho,pres,rgas,temp
  !
continue
  !
  ! 1) Copy the interior mass fractions to the wall (if TPG)
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv)
    do m = 1,ncs
      do n = n1,n2
        face_cv(n,m,2) = face_cv(n,m,1)/face_cv(n,nec,1)
      end do
    end do
  end if
  !
  ! 2) Compute the mixture gas constant and temporarily save it in
  !    the gamma location of face_gp (mostly for TPG)
  !
  if (ncs > 0) then
    call compute_mixture_rgas(n1,n2,face_cv(:,1:ncs,2),face_gp(:,ngamm,2))
  else
    !$ACC PARALLEL LOOP PRESENT(face_gp)
    do n = n1,n2
      face_gp(n,ngamm,2) = rgasref_nd
    end do
  end if
  !
  ! 3) Copy the interior pressure to the wall.
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,face_gp)
  do n = n1,n2
    face_cv(n,nee,2)   = face_gp(n,npres,1)
    face_gp(n,npres,2) = face_gp(n,npres,1)
  end do
  !
  ! 4) Get the wall density depending on whether or not this is an
  !    isothermal wall BC. We also need to save the wall temperature to
  !    face_gp(:,ntemp,2) which will be used to convert the primitive
  !    variables at the wall to conserved variables.
  !
  if (bc_type == bc_isothermal_wall) then
    !
    ! For isothermal walls, use the ideal gas law to compute the wall density
    ! using the wall temperature, interior pressure, and interior gas constant
    ! (based on interior mass fraction mixture for TPG).
    !
    ! For the wall temperature, use the relaxation coefficient to ease into
    ! the specified wall temperature for this BC.
    !
    !$ACC PARALLEL LOOP PRESENT(face_cv,face_gp)
    do n = n1,n2
      !
      pres = face_cv(n,nee,2)
      rgas = face_gp(n,ngamm,2)
      temp = face_gp(n,ntemp,1)*(one-rlx) + wall_temperature*rlx
      !
      face_cv(n,nec,2) = pres / (rgas * temp)
      face_gp(n,ntemp,2) = temp
      !
    end do
    !
  else
    !
    ! For adiabatic and slip walls, just copy the internal density to the wall
    !
    !$ACC PARALLEL LOOP PRESENT(face_cv,face_gp)
    do n = n1,n2
      !
      rho  = face_cv(n,nec,1)
      pres = face_cv(n,nee,2)
      rgas = face_gp(n,ngamm,2)
      !
      face_cv(n,nec,2) = rho
      face_gp(n,ntemp,2) = pres / (rgas * rho)
      !
    end do
    !
  end if
  !
end subroutine get_bc_wall_gp_solution
!
!###############################################################################
!
#ifdef DISABLE_FOR_GPU
pure &
#endif
subroutine get_bc_wall_vel_solution_3d(n1,n2,bc_type,facexyz,facenrm, &
                                       face_cv,rlx,pv_in,rotation_speed)
  !
  ! ####################################################################
  ! ####################################################################
  ! ##  This is for all wall-type boundary conditions
  ! ##   - Slip wall / symmetry BC
  ! ##   - No-slip adiabatic  wall, with/without imposed wall velocity
  ! ##   - No-slip isothermal wall, with/without imposed wall velocity
  ! ####################################################################
  ! ####################################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! BC_TYPE : Integer identifying this boundary condition type
  ! FACEXYZ : Coordinates for the boundary flux points
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conserved    variables at the flux points
  !           NOTE: The wall states will be saved as primitive variables
  !                 i.e., face_cv(:,:,2) = primitive variables
  ! RLX     : BC relaxation coefficient
  ! PV_IN   : Primitive    variables specified for this BC
  ! ROTATION_SPEED : RPM for rotating non-slip wall
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: bc_type
  real(wp),    intent(in) :: facexyz(:,:)
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: rlx
  real(wp),    intent(in) :: pv_in(:)
  real(wp),    intent(in) :: rotation_speed
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz,y,z,speed
  real(wp) :: vx,vx_slip,vx_wall,vx_bc,vtx_bc
  real(wp) :: vy,vy_slip,vy_wall,vy_bc,vty_bc
  real(wp) :: vz,vz_slip,vz_wall,vz_bc,vtz_bc
  real(wp) :: rho,rho_inv,vn,vn_bc,vmag_bc,vmag_tan
  !
  logical(lk) :: use_relaxation
  logical(lk) :: no_wall_velocity
  !
  !.. Local Parameters ..
  ! RPM_CONVERSION : converts RPM to velocity magnitude when multiplied
  !                  with the distance r from the rotational axis.
  real(wp), parameter :: rpm_conversion = six*ten*two*pi
  !
continue
  !
  rotating_wall_check: &
  if (abs(rotation_speed) > eps1) then
    !
    speed = rotation_speed*rpm_conversion
    !
    !#######################################################################
    !#######################################################################
    ! NOTE: This assumes that the rotation axis is aligned with the x-axis
    !#######################################################################
    !#######################################################################
    !
    !$ACC PARALLEL LOOP PRESENT(facexyz,face_cv)
    do n = n1,n2
      !
      y = facexyz(2,n)
      z = facexyz(3,n)
      !
      face_cv(n,nmx,2) = zero
      face_cv(n,nmy,2) = -z * speed
      face_cv(n,nmz,2) =  y * speed
      !
    end do
    !
  else rotating_wall_check
    !
    vx_bc = pv_in(nmx)
    vy_bc = pv_in(nmy)
    vz_bc = pv_in(nmz)
    vmag_bc = sqrt( vx_bc*vx_bc + vy_bc*vy_bc + vz_bc*vz_bc )
    !
    use_relaxation = (abs(one-rlx) > eps6)
    no_wall_velocity = (vmag_bc < eps6)
    !
    fixed_wall_vel_if: &
    if (any(bc_type == [bc_slip_wall,bc_symmetry])) then
      !
      ! For a slip wall or symmetry BC, set the boundary velocity to 
      ! the projection of the interior velocity tangent to the boundary.
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
        !
        vn = vx*unx + vy*uny + vz*unz
        !
        vx_slip = vx - vn*unx
        vy_slip = vy - vn*uny
        vz_slip = vz - vn*unz
        !
        face_cv(n,nmx,2) = vx_slip
        face_cv(n,nmy,2) = vy_slip
        face_cv(n,nmz,2) = vz_slip
        !
      end do
      !
    else if (no_wall_velocity .and. .not.use_relaxation) then fixed_wall_vel_if
      !
      ! No-slip wall with no imposed wall velocity, and no relaxation
      !
      !$ACC PARALLEL LOOP PRESENT(face_cv)
      do n = n1,n2
        face_cv(n,nmx,2) = zero
        face_cv(n,nmy,2) = zero
        face_cv(n,nmz,2) = zero
      end do
      !
    else if (no_wall_velocity .and. use_relaxation) then fixed_wall_vel_if
      !
      ! No-slip wall with no imposed wall velocity, with a non-zero relaxation
      ! coefficient to ease between a slip wall and no-slip wall
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
        !
        vn = vx*unx + vy*uny + vz*unz
        !
        vx_slip = vx - vn*unx
        vy_slip = vy - vn*uny
        vz_slip = vz - vn*unz
        !
        face_cv(n,nmx,2) = vx_slip*(one-rlx) ! + zero*rlx
        face_cv(n,nmy,2) = vy_slip*(one-rlx) ! + zero*rlx
        face_cv(n,nmz,2) = vz_slip*(one-rlx) ! + zero*rlx
        !
      end do
      !
    else fixed_wall_vel_if
      !
      ! No-slip wall with non-zero imposed wall velocity,
      ! with or without relaxation
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
        ! Magnitude of the imposed wall velocity normal to the surface
        vn_bc = vx_bc*unx + vy_bc*uny + vz_bc*unz
        !
        ! Imposed wall velocity components tangent to the surface
        vtx_bc = vx_bc - vn_bc*unx
        vty_bc = vy_bc - vn_bc*uny
        vtz_bc = vz_bc - vn_bc*unz
        !
        ! Magnitude of the imposed wall velocity tangent to the surface
        vmag_tan = sqrt( vtx_bc*vtx_bc + vty_bc*vty_bc + vtz_bc*vtz_bc )
        !
        if (vmag_tan < eps6) then
          ! Not really sure what to do if the magnitude of the imposed wall
          ! velocity tangent to the surface is very close to zero. This would
          ! be an odd case with the boundary normal being in the same direction
          ! as the imposed wall velocity, which I dont think should happen.
          ! For now, just set the imposed wall velocity to zero.
          vx_wall = zero
          vy_wall = zero
          vz_wall = zero
        else
          vx_wall = vmag_bc * (vtx_bc/vmag_tan)
          vy_wall = vmag_bc * (vty_bc/vmag_tan)
          vz_wall = vmag_bc * (vtz_bc/vmag_tan)
        end if
        !
        ! Compute the slip wall boundary velocity
        !
        rho = face_cv(n,nec,1)
        rho_inv = one/rho
        vx = face_cv(n,nmx,1)*rho_inv
        vy = face_cv(n,nmy,1)*rho_inv
        vz = face_cv(n,nmz,1)*rho_inv
        !
        vn = vx*unx + vy*uny + vz*unz
        !
        vx_slip = vx - vn*unx
        vy_slip = vy - vn*uny
        vz_slip = vz - vn*unz
        !
        ! Set the wall boundary velocity to a relaxation combination of the
        ! slip wall boundary velocity and the imposed wall velocity
        !
        face_cv(n,nmx,2) = vx_slip*(one-rlx) + vx_wall*rlx
        face_cv(n,nmy,2) = vy_slip*(one-rlx) + vy_wall*rlx
        face_cv(n,nmz,2) = vz_slip*(one-rlx) + vz_wall*rlx
        !
      end do
      !
    end if fixed_wall_vel_if
    !
  end if rotating_wall_check
  !
end subroutine get_bc_wall_vel_solution_3d
!
!###############################################################################
!
#ifdef DISABLE_FOR_GPU
pure &
#endif
subroutine get_bc_wall_vel_solution_2d(n1,n2,bc_type,facexyz,facenrm, &
                                       face_cv,rlx,pv_in,rotation_speed)
  !
  ! ####################################################################
  ! ####################################################################
  ! ##  This is for all wall-type boundary conditions
  ! ##   - Slip wall / symmetry BC
  ! ##   - No-slip adiabatic  wall, with/without imposed wall velocity
  ! ##   - No-slip isothermal wall, with/without imposed wall velocity
  ! ####################################################################
  ! ####################################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! BC_TYPE : Integer identifying this boundary condition type
  ! FACEXYZ : Coordinates for the boundary flux points
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conserved    variables at the flux points
  !           NOTE: The wall states will be saved as primitive variables
  !                 i.e., face_cv(:,:,2) = primitive variables
  ! RLX     : BC relaxation coefficient
  ! PV_IN   : Primitive    variables specified for this BC
  ! ROTATION_SPEED : RPM for rotating non-slip wall
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: bc_type
  real(wp),    intent(in) :: facexyz(:,:)
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: rlx
  real(wp),    intent(in) :: pv_in(:)
  real(wp),    intent(in) :: rotation_speed
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,rny,unx,uny,x,y,speed
  real(wp) :: vx,vx_slip,vx_wall,vx_bc,vtx_bc
  real(wp) :: vy,vy_slip,vy_wall,vy_bc,vty_bc
  real(wp) :: rho,rho_inv,vn,vn_bc,vmag_bc,vmag_tan
  !
  logical(lk) :: use_relaxation
  logical(lk) :: no_wall_velocity
  !
  !.. Local Parameters ..
  ! RPM_CONVERSION : converts RPM to velocity magnitude when multiplied
  !                  with the distance r from the rotational axis.
  real(wp), parameter :: rpm_conversion = six*ten*two*pi
  !
continue
  !
  rotating_wall_check: &
  if (abs(rotation_speed) > eps1) then
    !
    speed = rotation_speed*rpm_conversion
    !
    !#######################################################################
    !#######################################################################
    ! NOTE: This assumes that the rotation axis is aligned with the x-axis
    !#######################################################################
    !#######################################################################
    !
    !$ACC PARALLEL LOOP PRESENT(facexyz,face_cv)
    do n = n1,n2
      !
      x = facexyz(1,n)
      y = facexyz(2,n)
      !
      face_cv(n,nmx,2) = -y * speed
      face_cv(n,nmy,2) =  x * speed
      !
    end do
    !
  else rotating_wall_check
    !
    vx_bc = pv_in(nmx)
    vy_bc = pv_in(nmy)
    vmag_bc = sqrt( vx_bc*vx_bc + vy_bc*vy_bc )
    !
    use_relaxation = (abs(one-rlx) > eps6)
    no_wall_velocity = (vmag_bc < eps6)
    !
    fixed_wall_vel_if: &
    if (any(bc_type == [bc_slip_wall,bc_symmetry])) then
      !
      ! For a slip wall or symmetry BC, set the boundary velocity to 
      ! the projection of the interior velocity tangent to the boundary.
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
        !
        vn = vx*unx + vy*uny
        !
        vx_slip = vx - vn*unx
        vy_slip = vy - vn*uny
        !
        face_cv(n,nmx,2) = vx_slip
        face_cv(n,nmy,2) = vy_slip
        !
      end do
      !
    else if (no_wall_velocity .and. .not.use_relaxation) then fixed_wall_vel_if
      !
      ! No-slip wall with no imposed wall velocity, and no relaxation
      !
      !$ACC PARALLEL LOOP PRESENT(face_cv)
      do n = n1,n2
        face_cv(n,nmx,2) = zero
        face_cv(n,nmy,2) = zero
      end do
      !
    else if (no_wall_velocity .and. use_relaxation) then fixed_wall_vel_if
      !
      ! No-slip wall with no imposed wall velocity, with a non-zero relaxation
      ! coefficient to ease between a slip wall and no-slip wall
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
        !
        vn = vx*unx + vy*uny
        !
        vx_slip = vx - vn*unx
        vy_slip = vy - vn*uny
        !
        face_cv(n,nmx,2) = vx_slip*(one-rlx) ! + zero*rlx
        face_cv(n,nmy,2) = vy_slip*(one-rlx) ! + zero*rlx
        !
      end do
      !
    else fixed_wall_vel_if
      !
      ! No-slip wall with non-zero imposed wall velocity,
      ! with or without relaxation
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
        ! Magnitude of the imposed wall velocity normal to the surface
        vn_bc = vx_bc*unx + vy_bc*uny
        !
        ! Imposed wall velocity components tangent to the surface
        vtx_bc = vx_bc - vn_bc*unx
        vty_bc = vy_bc - vn_bc*uny
        !
        ! Magnitude of the imposed wall velocity tangent to the surface
        vmag_tan = sqrt( vtx_bc*vtx_bc + vty_bc*vty_bc )
        !
        if (vmag_tan < eps6) then
          ! Not really sure what to do if the magnitude of the imposed wall
          ! velocity tangent to the surface is very close to zero. This would
          ! be an odd case with the boundary normal being in the same direction
          ! as the imposed wall velocity, which I dont think should happen.
          ! For now, just set the imposed wall velocity to zero.
          vx_wall = zero
          vy_wall = zero
        else
          vx_wall = vmag_bc * (vtx_bc/vmag_tan)
          vy_wall = vmag_bc * (vty_bc/vmag_tan)
        end if
        !
        ! Compute the slip wall boundary velocity
        !
        rho = face_cv(n,nec,1)
        rho_inv = one/rho
        vx = face_cv(n,nmx,1)*rho_inv
        vy = face_cv(n,nmy,1)*rho_inv
        !
        vn = vx*unx + vy*uny
        !
        vx_slip = vx - vn*unx
        vy_slip = vy - vn*uny
        !
        ! Set the wall boundary velocity to a relaxation combination of the
        ! slip wall boundary velocity and the imposed wall velocity
        !
        face_cv(n,nmx,2) = vx_slip*(one-rlx) + vx_wall*rlx
        face_cv(n,nmy,2) = vy_slip*(one-rlx) + vy_wall*rlx
        !
      end do
      !
    end if fixed_wall_vel_if
    !
  end if rotating_wall_check
  !
end subroutine get_bc_wall_vel_solution_2d
!
!###############################################################################
!
#ifdef DISABLE_FOR_GPU
pure &
#endif
subroutine get_bc_wall_vel_solution_1d(n1,n2,bc_type,facenrm, &
                                       face_cv,rlx,pv_in)
  !
  ! ####################################################################
  ! ####################################################################
  ! ##  This is for all wall-type boundary conditions
  ! ##   - Slip wall / symmetry BC
  ! ##   - No-slip adiabatic  wall, with/without imposed wall velocity
  ! ##   - No-slip isothermal wall, with/without imposed wall velocity
  ! ####################################################################
  ! ####################################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! BC_TYPE : Integer identifying this boundary condition type
  ! FACENRM : Face outward normal vector at the flux points
  ! FACE_CV : Conserved    variables at the flux points
  !           NOTE: The wall states will be saved as primitive variables
  !                 i.e., face_cv(:,:,2) = primitive variables
  ! RLX     : BC relaxation coefficient
  ! PV_IN   : Primitive    variables specified for this BC
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: bc_type
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: rlx
  real(wp),    intent(in) :: pv_in(:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,unx
  real(wp) :: vx,vx_slip,vx_wall,vx_bc,vtx_bc
  real(wp) :: vn,vn_bc,vmag_bc,vmag_tan
  !
  logical(lk) :: use_relaxation
  logical(lk) :: no_wall_velocity
  !
continue
  !
  vx_bc = pv_in(nmx)
  vmag_bc = sqrt( vx_bc*vx_bc )
  !
  use_relaxation = (abs(one-rlx) > eps6)
  no_wall_velocity = (vmag_bc < eps6)
  !
  fixed_wall_vel_if: &
  if (any(bc_type == [bc_slip_wall,bc_symmetry])) then
    !
    ! For a slip wall or symmetry BC, set the boundary velocity to 
    ! the projection of the interior velocity tangent to the boundary.
    !
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
    do n = n1,n2
      !
      rnx = facenrm(1,n)
      ds = sqrt( rnx*rnx )
      unx = rnx/ds
      !
      vx = face_cv(n,nmx,1)/face_cv(n,nec,1)
      !
      vn = vx*unx
      !
      vx_slip = vx - vn*unx
      !
      face_cv(n,nmx,2) = vx_slip
      !
    end do
    !
  else if (no_wall_velocity .and. .not.use_relaxation) then fixed_wall_vel_if
    !
    ! No-slip wall with no imposed wall velocity, and no relaxation
    !
    !$ACC PARALLEL LOOP PRESENT(face_cv)
    do n = n1,n2
      face_cv(n,nmx,2) = zero
    end do
    !
  else if (no_wall_velocity .and. use_relaxation) then fixed_wall_vel_if
    !
    ! No-slip wall with no imposed wall velocity, with a non-zero relaxation
    ! coefficient to ease between a slip wall and no-slip wall
    !
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
    do n = n1,n2
      !
      rnx = facenrm(1,n)
      ds = sqrt( rnx*rnx )
      unx = rnx/ds
      !
      vx = face_cv(n,nmx,1)/face_cv(n,nec,1)
      !
      vn = vx*unx
      !
      vx_slip = vx - vn*unx
      !
      face_cv(n,nmx,2) = vx_slip*(one-rlx) ! + zero*rlx
      !
    end do
    !
  else fixed_wall_vel_if
    !
    ! No-slip wall with non-zero imposed wall velocity,
    ! with or without relaxation
    !
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv)
    do n = n1,n2
      !
      rnx = facenrm(1,n)
      ds = sqrt( rnx*rnx )
      unx = rnx/ds
      !
      ! Magnitude of the imposed wall velocity normal to the surface
      vn_bc = vx_bc*unx
      !
      ! Imposed wall velocity components tangent to the surface
      vtx_bc = vx_bc - vn_bc*unx
      !
      ! Magnitude of the imposed wall velocity tangent to the surface
      vmag_tan = sqrt( vtx_bc*vtx_bc )
      !
      if (vmag_tan < eps6) then
        ! Not really sure what to do if the magnitude of the imposed wall
        ! velocity tangent to the surface is very close to zero. This would
        ! be an odd case with the boundary normal being in the same direction
        ! as the imposed wall velocity, which I dont think should happen.
        ! For now, just set the imposed wall velocity to zero.
        vx_wall = zero
      else
        vx_wall = vmag_bc * (vtx_bc/vmag_tan)
      end if
      !
      ! Compute the slip wall boundary velocity
      !
      vx = face_cv(n,nmx,1)/face_cv(n,nec,1)
      !
      vn = vx*unx
      !
      vx_slip = vx - vn*unx
      !
      ! Set the wall boundary velocity to a relaxation combination of the
      ! slip wall boundary velocity and the imposed wall velocity
      !
      face_cv(n,nmx,2) = vx_slip*(one-rlx) + vx_wall*rlx
      !
    end do
    !
  end if fixed_wall_vel_if
  !
end subroutine get_bc_wall_vel_solution_1d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine get_bc_wall_gradient(n1,n2,bc_type,facenrm, &
                                face_cv,face_gp,face_dcvdx)
  !
  ! ####################################################################
  ! ####################################################################
  ! ##  Gradients for all wall-type boundary conditions
  ! ##   - Slip wall / symmetry BC
  ! ##   - No-slip adiabatic  wall
  ! ##   - No-slip isothermal wall
  ! ####################################################################
  ! ####################################################################
  !
  ! N1         : Beginning face point for this boundary condition
  ! N2         : Ending    face point for this boundary condition
  ! BC_TYPE    : Integer identifying this boundary condition type
  ! FACENRM    : Face outward normal vector at the flux points
  ! FACE_CV    : Conserved    variables at the flux points from the interior
  ! FACE_GP    : Gas property variables at the flux points from the interior
  ! FACE_DCVDX : Gradients of the conserved variables at the flux points
  !              NOTE: The wall gradients will be saved as gradients of the
  !                    primitive variables, i.e.,
  !                    face_dcvdx(:,:,:,2) = gradients of primitive variables
  !
  !.. Use Statements
  use gpu_mod, only : grad_cv_to_grad_pv_gpu
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: bc_type
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:)
  real(wp),    intent(in) :: face_gp(:,:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: l,m,n,nr,np
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: dpvdx(:,:,:)
  !
  !.. Local Parameters ..
 !logical(lk), parameter :: use_simple_zero_temp_gradient = true
  logical(lk), parameter :: use_simple_zero_temp_gradient = fals
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
  call grad_cv_to_grad_pv_gpu(n1,n2,face_cv,face_gp,face_dcvdx(:,:,:,1),dpvdx)
  !
  ! Get the primitive variable gradients at the wall using the interior values
  ! based on the wall BC
  !
  if (any(bc_type == [bc_slip_wall,bc_symmetry])) then
    !
    ! Reflect the interior primitive variable gradients using the face normals
    ! to get the wall vaules
    !
    if (nr == 3) then
      call get_slip_wall_gradient_3d(n1,n2,facenrm,dpvdx,face_dcvdx(:,:,:,2))
    else
      call get_slip_wall_gradient_2d(n1,n2,facenrm,dpvdx,face_dcvdx(:,:,:,2))
    end if
    !
  else ! Solid/no-slip wall BC
    !
    ! First, set the mass fraction densities to zero at the wall (if TPG)
    !
    if (ncs > 0) then
      !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(face_dcvdx)
      do m = 1,ncs
        do l = 1,nr
          do n = n1,n2
            face_dcvdx(n,l,m,2) = zero
          end do
        end do
      end do
    end if
    !
    ! Next, copy over the interior primitive variable gradients
    ! over to the wall
    !
    !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(dpvdx,face_dcvdx)
    do m = nec,nq ! skip species densities
      do l = 1,nr
        do n = n1,n2
          face_dcvdx(n,l,m,2) = dpvdx(n,l,m)
        end do
      end do
    end do
    !
    ! For adiabatic walls, zero out the temperature gradient normal to the wall
    !
    if (bc_type == bc_adiabatic_wall) then
      !
      if (use_simple_zero_temp_gradient) then
        !
        !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_dcvdx)
        do l = 1,nr
          do n = n1,n2
            face_dcvdx(n,l,nee,2) = zero
          end do
        end do
        !
      else
        !
        call get_adiabatic_temp_gradient(n1,n2,facenrm,dpvdx(:,:,nee), &
                                         face_dcvdx(:,:,nee,2))
        !
      end if
      !
    end if
    !
  end if
  !
  !$ACC END DATA
  deallocate ( dpvdx )
  !
end subroutine get_bc_wall_gradient
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine get_adiabatic_temp_gradient(n1,n2,facenrm,dtdxi,dtdxw)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Temperature gradient for no-slip adiabatic wall
  ! #####################################################
  ! #####################################################
  !
  ! N1         : Beginning face point for this boundary condition
  ! N2         : Ending    face point for this boundary condition
  ! FACENRM    : Face outward normal vector at the flux points
  ! DTDXI      : Temperature gradient from the interior
  ! DTDXW      : Temperature gradient at the wall
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: dtdxi(:,:)
  real(wp), intent(inout) :: dtdxw(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: dtdn
  !
continue
  !
  nr = size(dtdxi,dim=2)
  !
  ! Compute the temperature gradient at the wall for an adiabatic no-slip wall
  ! An adiabatic wall is supposed to have a zero temperature gradient normal
  ! to the wall, so subtract out the interior normal temperature gradient at
  ! the wall from the interior temperature gradient, which will result in the
  ! temperature gradient tangent to the wall.
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(facenrm,dtdxi,dtdxw)
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
      dtdn = dtdxi(n,1)*unx + dtdxi(n,2)*uny + dtdxi(n,3)*unz
      !
      dtdxw(n,1) = dtdxi(n,1) - unx*dtdn
      dtdxw(n,2) = dtdxi(n,2) - uny*dtdn
      dtdxw(n,3) = dtdxi(n,3) - unz*dtdn
      !
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(facenrm,dtdxi,dtdxw)
    do n = n1,n2
      !
      rnx = facenrm(1,n)
      rny = facenrm(2,n)
      ds = sqrt( rnx*rnx + rny*rny )
      unx = rnx/ds
      uny = rny/ds
      !
      dtdn = dtdxi(n,1)*unx + dtdxi(n,2)*uny
      !
      dtdxw(n,1) = dtdxi(n,1) - unx*dtdn
      dtdxw(n,2) = dtdxi(n,2) - uny*dtdn
      !
    end do
  end if
  !
end subroutine get_adiabatic_temp_gradient
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine get_slip_wall_gradient_3d(n1,n2,facenrm,dpvdx,face_dcvdx)
  !
  ! ###########################################################
  ! ###########################################################
  ! ##  Slip wall / symmetry boundary condition for gradients
  ! ###########################################################
  ! ###########################################################
  !
  ! N1         : Beginning face point for this boundary condition
  ! N2         : Ending    face point for this boundary condition
  ! FACENRM    : Face outward normal vector at the flux points
  ! DPVDX      : Gradients of the primitive variables from the interior
  ! FACE_DCVDX : Gradients of the primitive variables at the boundary
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: dpvdx(:,:,:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: drdn,dudn,dvdn,dwdn,dtdn,dsdn
  !
continue
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
    face_dcvdx(n,1,nec) = dpvdx(n,1,nec) - unx*drdn
    face_dcvdx(n,2,nec) = dpvdx(n,2,nec) - uny*drdn
    face_dcvdx(n,3,nec) = dpvdx(n,3,nec) - unz*drdn
    !
    face_dcvdx(n,1,nmx) = dpvdx(n,1,nmx) - unx*dudn
    face_dcvdx(n,2,nmx) = dpvdx(n,2,nmx) - uny*dudn
    face_dcvdx(n,3,nmx) = dpvdx(n,3,nmx) - unz*dudn
    !
    face_dcvdx(n,1,nmy) = dpvdx(n,1,nmy) - unx*dvdn
    face_dcvdx(n,2,nmy) = dpvdx(n,2,nmy) - uny*dvdn
    face_dcvdx(n,3,nmy) = dpvdx(n,3,nmy) - unz*dvdn
    !
    face_dcvdx(n,1,nmz) = dpvdx(n,1,nmz) - unx*dwdn
    face_dcvdx(n,2,nmz) = dpvdx(n,2,nmz) - uny*dwdn
    face_dcvdx(n,3,nmz) = dpvdx(n,3,nmz) - unz*dwdn
    !
    face_dcvdx(n,1,nee) = dpvdx(n,1,nee) - unx*dtdn
    face_dcvdx(n,2,nee) = dpvdx(n,2,nee) - uny*dtdn
    face_dcvdx(n,3,nee) = dpvdx(n,3,nee) - unz*dtdn
    !
  end do
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(facenrm,face_dcvdx,dpvdx)
    do m = 1,ncs
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
        dsdn = dpvdx(n,1,m)*unx + dpvdx(n,2,m)*uny + dpvdx(n,3,m)*unz
        !
        face_dcvdx(n,1,m) = dpvdx(n,1,m) - unx*dsdn
        face_dcvdx(n,2,m) = dpvdx(n,2,m) - uny*dsdn
        face_dcvdx(n,3,m) = dpvdx(n,3,m) - unz*dsdn
        !
      end do
    end do
  end if
  !
end subroutine get_slip_wall_gradient_3d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine get_slip_wall_gradient_2d(n1,n2,facenrm,dpvdx,face_dcvdx)
  !
  ! ###########################################################
  ! ###########################################################
  ! ##  Slip wall / symmetry boundary condition for gradients
  ! ###########################################################
  ! ###########################################################
  !
  ! N1         : Beginning face point for this boundary condition
  ! N2         : Ending    face point for this boundary condition
  ! FACENRM    : Face outward normal vector at the flux points
  ! DPVDX      : Gradients of the primitive variables from the interior
  ! FACE_DCVDX : Gradients of the primitive variables at the boundary
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: dpvdx(:,:,:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: ds,rnx,rny,unx,uny
  real(wp) :: drdn,dudn,dvdn,dtdn,dsdn
  !
continue
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
    face_dcvdx(n,1,nec) = dpvdx(n,1,nec) - unx*drdn
    face_dcvdx(n,2,nec) = dpvdx(n,2,nec) - uny*drdn
    !
    face_dcvdx(n,1,nmx) = dpvdx(n,1,nmx) - unx*dudn
    face_dcvdx(n,2,nmx) = dpvdx(n,2,nmx) - uny*dudn
    !
    face_dcvdx(n,1,nmy) = dpvdx(n,1,nmy) - unx*dvdn
    face_dcvdx(n,2,nmy) = dpvdx(n,2,nmy) - uny*dvdn
    !
    face_dcvdx(n,1,nee) = dpvdx(n,1,nee) - unx*dtdn
    face_dcvdx(n,2,nee) = dpvdx(n,2,nee) - uny*dtdn
    !
  end do
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(facenrm,face_dcvdx,dpvdx)
    do m = 1,ncs
      do n = n1,n2
        !
        rnx = facenrm(1,n)
        rny = facenrm(2,n)
        ds = sqrt( rnx*rnx + rny*rny )
        unx = rnx/ds
        uny = rny/ds
        !
        dsdn = dpvdx(n,1,m)*unx + dpvdx(n,2,m)*uny
        !
        face_dcvdx(n,1,m) = dpvdx(n,1,m) - unx*dsdn
        face_dcvdx(n,2,m) = dpvdx(n,2,m) - uny*dsdn
        !
      end do
    end do
  end if
  !
end subroutine get_slip_wall_gradient_2d
!
!###############################################################################
!
end module bc_wall_mod
