module bc_isothermal_wall_mod
  !
  !.. Use Statements
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee
  use eqn_idx, only : ncs,npres,ntemp,ngamm
  !
  implicit none
  !
  private
  !
  public :: get_bc_isothermal_wall_solution
  public :: get_bc_wall_gradient
  !
contains
!
!###############################################################################
!
pure subroutine get_bc_isothermal_wall_solution(n1,n2,face_cv,face_gp, &
                                                wall_temperature,pv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip isothermal wall
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACE_CV : Conserved    variables at the flux points
  ! FACE_GP : Gas property variables at the flux points
  ! PV_IN   : Primitive    variables specified for this BC
  ! WALL_TEMPERATURE : isothermal wall temperature
  !
  !.. Use Statements ..
  use ovar, only : rgasref_nd
  !
  use gpu_mod,    only : pv2cv_gpu
  use thermo_mod, only : compute_mixture_rgas
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: wall_temperature
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp), intent(inout) :: face_gp(:,:,:)
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: pv_in(:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  ! Compute the mixture gas constant and store it in the gamma location
  ! of the gas properties for the exterior state
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv)
    do m = 1,ncs
      do n = n1,n2
        face_cv(n,m,2) = face_cv(n,m,1)/face_cv(n,nec,1)
      end do
    end do
    call compute_mixture_rgas(n1,n2,face_cv(:,1:ncs,2),face_gp(:,ngamm,2))
  else
    !$ACC PARALLEL LOOP PRESENT(face_gp)
    do n = n1,n2
      face_gp(n,ngamm,2) = rgasref_nd
    end do
  end if
  !
  ! Compute the isothermal wall boundary conditions
  !
  if (nr == 3) then
    call get_bc_isothermal_wall_solution_3d(n1,n2,face_cv,face_gp, &
                                            wall_temperature,pv_in)
  else if (nr == 2) then
    call get_bc_isothermal_wall_solution_2d(n1,n2,face_cv,face_gp, &
                                            wall_temperature,pv_in)
  else
    call get_bc_isothermal_wall_solution_1d(n1,n2,face_cv,face_gp, &
                                            wall_temperature,pv_in)
  end if
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,face_gp)
  do n = n1,n2
    face_gp(n,npres,2) = face_cv(n,nee,2)
    face_gp(n,ntemp,2) = wall_temperature
  end do
  !
  ! The above subroutines will return the boundary condition states as
  ! primitive variables. These need to be converted to conservative variables.
  !
  call pv2cv_gpu(n1,n2,face_cv(:,:,2),face_gp(:,:,2))
  !
end subroutine get_bc_isothermal_wall_solution
!
!###############################################################################
!
pure subroutine get_bc_isothermal_wall_solution_3d(n1,n2,face_cv,face_gp, &
                                                   wall_temperature,pv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip isothermal wall
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACE_CV : Conserved    variables at the flux points
  ! FACE_GP : Gas property variables at the flux points
  ! PV_IN   : Primitive    variables specified for this BC
  ! WALL_TEMPERATURE : isothermal wall temperature
  !
  !.. Use Statements ..
  use ovar, only : walls_are_exact
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: wall_temperature
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: pv_in(:)
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: rho,rho_inv,vx,vy,vz,pres,rgas
  real(wp) :: wall_vx,wall_vy,wall_vz
  real(wp) :: velocity_const
  real(wp) :: moving_wall_const
  !
continue
  !
  velocity_const = merge(zero,-one,walls_are_exact)
  moving_wall_const = merge(one,two,walls_are_exact)
  !
  wall_vx = zero
  wall_vy = zero
  wall_vz = zero
  if (present(pv_in)) then
    wall_vx = pv_in(nmx)*moving_wall_const
    wall_vy = pv_in(nmy)*moving_wall_const
    wall_vz = pv_in(nmz)*moving_wall_const
  end if
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,face_gp)
  do n = n1,n2
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    vz = face_cv(n,nmz,1)*rho_inv
    !
    pres = face_gp(n,npres,1)
    rgas = face_gp(n,ngamm,2)
    !
    face_cv(n,nec,2) = pres / (rgas * wall_temperature)
    face_cv(n,nmx,2) = vx*velocity_const + wall_vx
    face_cv(n,nmy,2) = vy*velocity_const + wall_vy
    face_cv(n,nmz,2) = vz*velocity_const + wall_vz
    face_cv(n,nee,2) = pres
    !
  end do
  !
end subroutine get_bc_isothermal_wall_solution_3d
!
!###############################################################################
!
pure subroutine get_bc_isothermal_wall_solution_2d(n1,n2,face_cv,face_gp, &
                                                   wall_temperature,pv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip isothermal wall
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACE_CV : Conserved    variables at the flux points
  ! FACE_GP : Gas property variables at the flux points
  ! PV_IN   : Primitive    variables specified for this BC
  ! WALL_TEMPERATURE : isothermal wall temperature
  !
  !.. Use Statements ..
  use ovar, only : walls_are_exact
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: wall_temperature
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: pv_in(:)
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: rho,rho_inv,vx,vy,pres,rgas
  real(wp) :: wall_vx,wall_vy
  real(wp) :: velocity_const
  real(wp) :: moving_wall_const
  !
continue
  !
  velocity_const = merge(zero,-one,walls_are_exact)
  moving_wall_const = merge(one,two,walls_are_exact)
  !
  wall_vx = zero
  wall_vy = zero
  if (present(pv_in)) then
    wall_vx = pv_in(nmx)*moving_wall_const
    wall_vy = pv_in(nmy)*moving_wall_const
  end if
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,face_gp)
  do n = n1,n2
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    !
    pres = face_gp(n,npres,1)
    rgas = face_gp(n,ngamm,2)
    !
    face_cv(n,nec,2) = pres / (rgas * wall_temperature)
    face_cv(n,nmx,2) = vx*velocity_const + wall_vx
    face_cv(n,nmy,2) = vy*velocity_const + wall_vy
    face_cv(n,nee,2) = pres
    !
  end do
  !
end subroutine get_bc_isothermal_wall_solution_2d
!
!###############################################################################
!
pure subroutine get_bc_isothermal_wall_solution_1d(n1,n2,face_cv,face_gp, &
                                                   wall_temperature,pv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  No-slip isothermal wall
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACE_CV : Conserved    variables at the flux points
  ! FACE_GP : Gas property variables at the flux points
  ! PV_IN   : Primitive    variables specified for this BC
  ! WALL_TEMPERATURE : isothermal wall temperature
  !
  !.. Use Statements ..
  use ovar, only : walls_are_exact
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: wall_temperature
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: pv_in(:)
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: pres,vx,wall_vx,rgas
  real(wp) :: velocity_const
  real(wp) :: moving_wall_const
  !
continue
  !
  velocity_const = merge(zero,-one,walls_are_exact)
  moving_wall_const = merge(one,two,walls_are_exact)
  !
  wall_vx = zero
  if (present(pv_in)) then
    wall_vx = pv_in(nmx)*moving_wall_const
  end if
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,face_gp)
  do n = n1,n2
    !
    vx = face_cv(n,nmx,1)/face_cv(n,nec,1)
    !
    pres = face_gp(n,npres,1)
    rgas = face_gp(n,ngamm,2)
    !
    face_cv(n,nec,2) = pres / (rgas * wall_temperature)
    face_cv(n,nmx,2) = vx*velocity_const + wall_vx
    face_cv(n,nee,2) = pres
    !
  end do
  !
end subroutine get_bc_isothermal_wall_solution_1d
!
!###############################################################################
!
pure &
subroutine get_bc_wall_gradient(n1,n2,face_cv,face_gp,face_dcvdx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Gradients for non-adiabatic, no-slip wall
  ! #####################################################
  ! #####################################################
  !
  ! N1         : Beginning face point for this boundary condition
  ! N2         : Ending    face point for this boundary condition
  ! FACE_CV    : Conserved    variables at the flux points
  ! FACE_GP    : Gas property variables at the flux points
  ! FACE_DCVDX : Gradients of the conserved variables at the flux points
  !
  !.. Use Statements
  use gpu_mod, only : grad_cv_to_grad_pv_gpu
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:)
  real(wp),    intent(in) :: face_gp(:,:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: l,m,n,nr
  !
  !.. Local Parameters ..
 !logical(lk), parameter :: use_bc_wall_grad_subs = true
  logical(lk), parameter :: use_bc_wall_grad_subs = fals
  !
continue
  !
  nr = size(face_dcvdx,dim=2)
  !
  if (use_bc_wall_grad_subs) then
    !
    if (nr == 3) then
      call get_bc_wall_gradient_3d(n1,n2,face_cv,face_dcvdx)
    else if (nr == 2) then
      call get_bc_wall_gradient_2d(n1,n2,face_cv,face_dcvdx)
    else
      call get_bc_wall_gradient_1d(n1,n2,face_cv,face_dcvdx)
    end if
    !
  else
    !
    ! This should do the same thing that the get_bc_wall_gradient_[123]d
    ! subroutines do, except that this will return a non-zero gradient
    ! for the species equations. We can just simply zero those out afterwards.
    !
    call grad_cv_to_grad_pv_gpu(n1,n2,face_cv,face_gp, &
                                face_dcvdx(:,:,:,1), &
                                face_dcvdx(:,:,:,2))
    !
    ! Set the mass fraction gradients to zero for wall boundaries
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
  ! ##  Gradients for non-adiabatic, no-slip wall
  ! #####################################################
  ! #####################################################
  !
  ! N1         : Beginning face point for this boundary condition
  ! N2         : Ending    face point for this boundary condition
  ! FACE_CV    : Conserved    variables at the flux points
  ! FACE_DCVDX : Gradients of the conserved variables at the flux points
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
  integer  :: m,n
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
  ! Using the product rule:
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
  ! Set the mass fraction gradients to zero for wall boundaries
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_dcvdx)
    do m = 1,ncs
      do n = n1,n2
        face_dcvdx(n,1,m,2) = zero
        face_dcvdx(n,2,m,2) = zero
        face_dcvdx(n,3,m,2) = zero
      end do
    end do
  end if
  !
end subroutine get_bc_wall_gradient_3d
!
!###############################################################################
!
pure subroutine get_bc_wall_gradient_2d(n1,n2,face_cv,face_dcvdx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Gradients for non-adiabatic, no-slip wall
  ! #####################################################
  ! #####################################################
  !
  ! N1         : Beginning face point for this boundary condition
  ! N2         : Ending    face point for this boundary condition
  ! FACE_CV    : Conserved    variables at the flux points
  ! FACE_DCVDX : Gradients of the conserved variables at the flux points
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
  integer  :: m,n
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
  ! Using the product rule:
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
  ! Set the mass fraction gradients to zero for wall boundaries
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_dcvdx)
    do m = 1,ncs
      do n = n1,n2
        face_dcvdx(n,1,m,2) = zero
        face_dcvdx(n,2,m,2) = zero
      end do
    end do
  end if
  !
end subroutine get_bc_wall_gradient_2d
!
!###############################################################################
!
pure subroutine get_bc_wall_gradient_1d(n1,n2,face_cv,face_dcvdx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Gradients for non-adiabatic, no-slip wall
  ! #####################################################
  ! #####################################################
  !
  ! N1         : Beginning face point for this boundary condition
  ! N2         : Ending    face point for this boundary condition
  ! FACE_CV    : Conserved    variables at the flux points
  ! FACE_DCVDX : Gradients of the conserved variables at the flux points
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
  integer  :: m,n
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
  ! Using the product rule:
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
    !
    face_dcvdx(n,1,nee,2) = r_inv*face_dcvdx(n,1,nee,1) &
                          - r_inv*face_dcvdx(n,1,nec,1)*face_cv(n,nee)*r_inv
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
    r_inv = one/face_cv(n,nec)
    face_dcvdx(n,1,nee,2) = face_dcvdx(n,1,nee,2) &
                          - face_dcvdx(n,1,nmx,2)*face_cv(n,nmx)*r_inv
  end do
  !
  ! Finish computing the temperature gradient by dividing through by
  ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
  do n = n1,n2
    face_dcvdx(n,1,nee,2) = face_dcvdx(n,1,nee,2) * cv_inv
  end do
  !
  ! Set the mass fraction gradients to zero for wall boundaries
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_dcvdx)
    do m = 1,ncs
      do n = n1,n2
        face_dcvdx(n,1,m,2) = zero
      end do
    end do
  end if
  !
end subroutine get_bc_wall_gradient_1d
!
!###############################################################################
!
end module bc_isothermal_wall_mod
