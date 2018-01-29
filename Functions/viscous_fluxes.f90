!
!###############################################################################
!###############################################################################
!###############################################################################
!
! Procedures for viscous fluxes
!
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function visc_flx_cv1(uu,dcvdx) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : cpref_nd,Pr
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nq),      intent(in) :: uu
  real(wp), dimension(1:nr,1:nq), intent(in) :: dcvdx
  !
  !.. Function Result ..
  real(wp), dimension(1:nr,1:nq) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,l
  real(wp) :: mvsc,mcdt
  real(wp) :: ave_ext_rate
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: vel
  real(wp), dimension(1:nr,1:nq) :: dpvdx
  !
continue
  !
  ! Initialize the viscous fluxes to zero.
  ! NOTE: This takes care of the continuity equation
  !
  return_value = zero
  !
  ! Get the gradients of the primitive variables using the
  ! gradients of the conservative variables and the product rule
  !
  dpvdx(1:nr,1:nq) = grad_cv_to_grad_pv_sp(uu,dcvdx)
  !
  ! Get the velocity vector
  !
  vel(:) = uu(nmb:nme)/uu(nec)
  !
  ! Compute the molecular viscosity and molecular conductivity
  !
  mvsc = viscosity_cv_sp(uu)
  mcdt = cpref_nd * mvsc / Pr
  !
  ! Compute the average extension rate
  !
  ave_ext_rate = zero
  do l = 1,nr
    ave_ext_rate = ave_ext_rate + dpvdx(l,nmb-1+l)
  end do
  ave_ext_rate = two3 * ave_ext_rate
  !
  ! Compute the viscous fluxes for the momentum equations which are just
  ! the components of the viscous stress tensor.
  ! NOTE: We are initially excluding the molecular viscosity term which
  !       will be multiplied to the whole tau array after these loops.
  !
  do i = 1,nr
    do j = 1,nr
        return_value(j,nmb-1+i) = dpvdx(i,nmb-1+j) + dpvdx(j,nmb-1+i)
    end do
  end do
  !
  ! Subtract the average extension rate from the normal stresses
  !
  do l = 1,nr
    return_value(l,nmb-1+l) = return_value(l,nmb-1+l) - ave_ext_rate
  end do
  !
  ! Multiply the viscous stresses by the molecular viscosity
  !
  return_value(1:nr,nmb:nme) = mvsc * return_value(1:nr,nmb:nme)
  !
  ! Compute the viscous fluxes for the energy equation
  !
  do l = 1,nr
    return_value(l,nq) = mcdt*dpvdx(l,nq) + dot( vel , return_value(l,nmb:nme) )
  end do
  !
end function visc_flx_cv1
!
!###############################################################################
!
pure function visc_flx_cv(uu,dcvdx) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : cpref_nd,Pr
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nq),      intent(in) :: uu
  real(wp), dimension(1:nr,1:nq), intent(in) :: dcvdx
  !
  !.. Function Result ..
  real(wp), dimension(1:nq,1:nr) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,l
  real(wp) :: mvsc,mcdt
  real(wp) :: ave_ext_rate
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: vel
  real(wp), dimension(1:nr,1:nq) :: dpvdx
  !
continue
  !
  ! Initialize the viscous fluxes to zero.
  ! NOTE: This takes care of the continuity equation
  !
  return_value = zero
  !
  ! Get the gradients of the primitive variables using the
  ! gradients of the conservative variables and the product rule
  !
  dpvdx(1:nr,1:nq) = grad_cv_to_grad_pv_sp(uu,dcvdx)
  !
  ! Get the velocity vector
  !
  vel(:) = uu(nmb:nme)/uu(nec)
  !
  ! Compute the molecular viscosity and molecular conductivity
  !
  mvsc = viscosity_cv_sp(uu)
  mcdt = cpref_nd * mvsc / Pr
  !
  ! Compute the average extension rate
  !
  ave_ext_rate = zero
  do l = 1,nr
    ave_ext_rate = ave_ext_rate + dpvdx(l,nmb-1+l)
  end do
  ave_ext_rate = two3 * ave_ext_rate
  !
  ! Compute the viscous fluxes for the momentum equations which are just
  ! the components of the viscous stress tensor.
  ! NOTE: We are initially excluding the molecular viscosity term which
  !       will be multiplied to the whole tau array after these loops.
  !
  do j = 1,nr
    do i = 1,nr
        return_value(nmb-1+i,j) = dpvdx(i,nmb-1+j) + dpvdx(j,nmb-1+i)
    end do
  end do
  !
  ! Subtract the average extension rate from the normal stresses
  !
  do l = 1,nr
    return_value(nmb-1+l,l) = return_value(nmb-1+l,l) - ave_ext_rate
  end do
  !
  ! Multiply the viscous stresses by the molecular viscosity
  !
  return_value(nmb:nme,1:nr) = mvsc * return_value(nmb:nme,1:nr)
  !
  ! Compute the viscous fluxes for the energy equation
  !
  do l = 1,nr
    return_value(nq,l) = mcdt*dpvdx(l,nq) + dot( vel , return_value(nmb:nme,l) )
  end do
  !
end function visc_flx_cv
!
!###############################################################################
!
pure function normal_visc_flx_cv(rnrm,uu,dcvdx) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nq),      intent(in) :: uu
  real(wp), dimension(1:nr),         intent(in) :: rnrm
  real(wp), dimension(1:nr,1:nq), intent(in) :: dcvdx
  !
  !.. Function Result ..
  real(wp), dimension(1:nq) :: return_value
  !
  !.. Local Scalars ..
  integer :: l
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)         :: unrm
  real(wp), dimension(1:nq,1:nr) :: vflx
  !
continue
  !
  unrm(:) = unit_vector( rnrm(:) )
  !
  vflx(1:nq,1:nr) = visc_flx_cv(uu,dcvdx)
  !
  return_value(1:nq) = unrm(1)*vflx(1:nq,1)
  do l = 2,nr
    return_value(1:nq) = return_value(1:nq) + unrm(l)*vflx(1:nq,l)
  end do
  !
end function normal_visc_flx_cv
!
!###############################################################################
!
pure function visc_wall_flx_cv(unrm,uwall,dpvdx) result(return_value)
  !
  ! UNRM  : Outward pointing unit normal to the wall
  !         boundary face for current flux point
  ! UWALL : Conservative variables exact solution at the wall
  ! DPVDX : Prescribed wall gradients of the primitive variables at the
  !         current flux point on the wall boundary face. Most gradients
  !         are extrapolated from the interior cell flux point on the
  !         boundary except for certain special cases, e.g., an adiabatic
  !         wall where dT/dx = dT/dy = 0.
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : cpref_nd,Pr
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nr),      intent(in) :: unrm
  real(wp), dimension(1:nq),      intent(in) :: uwall
  real(wp), dimension(1:nr,1:nq), intent(in) :: dpvdx
  !
  !.. Function Result ..
  real(wp), dimension(1:nq) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,l,m
  real(wp) :: mvsc,mcdt
  real(wp) :: ave_ext_rate
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)      :: vel
  real(wp), dimension(1:nr,1:nq) :: vflx
  !
continue
  !
  ! For slip wall: all gradients /= zero
  !                uwall /= zero
  !
  ! For adiabatic wall: Temperature gradient == zero
  !                     Velocity gradients /= zero
  !                     uwall(nmb:nme) [i.e., wall velocities] == zero
  !
  ! For isothermal wall: all gradients /= zero
  !                      uwall(nmb:nme) [i.e., wall velocities] == zero
  !
  ! Initialize the viscous fluxes to zero.
  ! NOTE: This takes care of the continuity equation
  !
  vflx = zero
  !
  ! Get the velocity vector
  !
  vel(:) = uwall(nmb:nme)/uwall(nec)
  !
  ! Compute the molecular viscosity and molecular conductivity
  !
  mvsc = viscosity_cv_sp(uwall)
  mcdt = cpref_nd * mvsc / Pr
  !
  ! Compute the average extension rate
  !
  ave_ext_rate = zero
  do l = 1,nr
    ave_ext_rate = ave_ext_rate + dpvdx(l,nmb-1+l)
  end do
  ave_ext_rate = two3 * ave_ext_rate
  !
  ! Compute the viscous fluxes for the momentum equations which are just
  ! the components of the viscous stress tensor.
  ! NOTE: We are initially excluding the molecular viscosity term which
  !       will be multiplied to the whole tau array after these loops.
  !
  do j = 1,nr
    do i = 1,nr
        vflx(i,nmb-1+j) = dpvdx(i,nmb-1+j) + dpvdx(j,nmb-1+i)
    end do
  end do
  !
  ! Subtract the average extension rate from the normal stresses
  !
  do l = 1,nr
    vflx(l,nmb-1+l) = vflx(l,nmb-1+l) - ave_ext_rate
  end do
  !
  ! Multiply the viscous stresses by the molecular viscosity
  !
  vflx(1:nr,nmb:nme) = mvsc * vflx(1:nr,nmb:nme)
  !
  ! Compute the viscous fluxes for the energy equation
  !
  do l = 1,nr
    vflx(l,nq) = mcdt*dpvdx(l,nq) + dot( vel , vflx(l,nmb:nme) )
  end do
  !
  ! Compute the normal viscous flux
  !
  do m = 1,nq
    return_value(m) = dot( vflx(1:nr,m) , unrm(1:nr) )
  end do
  !
end function visc_wall_flx_cv
!
!###############################################################################
!
pure function common_visc_flx_cv(unrm,ul,duldx,ur,durdx) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : visc_variables_to_ave
  use ovar,   only : LDG_tau, LDG_beta
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nr),         intent(in) :: unrm
  real(wp), dimension(1:nq),      intent(in) :: ul,ur
  real(wp), dimension(1:nr,1:nq), intent(in) :: duldx,durdx
  !
  !.. Function Result ..
  real(wp), dimension(1:nq) :: return_value
  !
  !.. Local Scalars ..
  integer  :: l,m
  real(wp) :: vjump,qjump,penalty
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nq)      :: ucom
  real(wp), dimension(1:nr,1:nq) :: qcom
  real(wp), dimension(1:nq,1:nr) :: lvflx
  real(wp), dimension(1:nq,1:nr) :: rvflx
  real(wp), dimension(1:nr,1:nq) :: vflux
  real(wp), dimension(1:nq,1:nr) :: fvcom
  !
continue
  !
  if (visc_variables_to_ave == visc_interface_fluxes) then
    !
    ! Get the left and right viscous fluxes
    !
    lvflx(1:nq,1:nr) = visc_flx_cv(ul,duldx)
    rvflx(1:nq,1:nr) = visc_flx_cv(ur,durdx)
    !
    ! Use the LDG method to compute the common viscous flux vector
    !
    ! Start by applying the averaging operator {{.}} to the viscous flux vectors
    !
    do m = 1,nq
      do l = 1,nr
        vflux(l,m) = half * (lvflx(m,l) + rvflx(m,l))
      end do
    end do
    !
    ! Apply the differencing (jump) operator [[.]] to the interface viscous
    ! flux vectors and then add it to the common viscous flux vector
    !
    do m = 1,nq
      !
      vjump = dot( unrm(1:nr) , lvflx(m,1:nr)-rvflx(m,1:nr) )
      !
      do l = 1,nr
        vflux(l,m) = vflux(l,m) + LDG_beta(l)*vjump
      end do
      !
    end do
    !
    ! Finally, add the penalty term, (simply the differencing (jump)
    ! operator [[.]] applied to the interface solutions) to the
    ! common viscous flux vector
    !
    do m = 1,nq
      penalty = LDG_tau * ( ur(m) - ul(m) )
      do l = 1,nr
        vflux(l,m) = vflux(l,m) + penalty * unrm(l)
      end do
    end do
    !
    ! Finally, dot the common viscous flux vector with the face unit normal
    !
    do m = 1,nq
      return_value(m) = dot( vflux(1:nr,m) , unrm(1:nr) )
    end do
    !
  else if (visc_variables_to_ave == visc_solution_and_gradients) then
    !
    ! Get the common solution at the interface using a simple average
    !
    do m = 1,nq
      ucom(m) = half * ( ul(m) + ur(m) )
    end do
    !
    ! Use the LDG method to compute the common gradient
    !
    ! Start by applying the averaging operator {{.}} to the gradients
    !
    do m = 1,nq
      do l = 1,nr
        qcom(l,m) = half * (duldx(l,m) + durdx(l,m))
      end do
    end do
    !
    ! Apply the differencing (jump) operator [[.]] to the
    ! gradients and then add it to the common gradient
    !
    do m = 1,nq
      !
      qjump = dot( unrm(1:nr) , duldx(1:nr,m)-durdx(1:nr,m) )
      !
      do l = 1,nr
        qcom(l,m) = qcom(l,m) + LDG_beta(l)*qjump
      end do
      !
    end do
    !
    ! Finally, add the penalty term, (simply the differencing (jump) operator
    ! [[.]] applied to the interface solutions) to the common gradient
    !
    do m = 1,nq
      penalty = LDG_tau * ( ur(m) - ul(m) )
      do l = 1,nr
        qcom(l,m) = qcom(l,m) + penalty * unrm(l)
      end do
    end do
    !
    ! Compute the common viscous flux vector using
    ! the common solution and common gradient
    !
    fvcom(1:nq,1:nr) = visc_flx_cv(ucom,qcom)
    !
    ! Finally, dot the common viscous flux vector with the face unit normal
    !
    do m = 1,nq
      return_value(m) = dot( fvcom(m,1:nr) , unrm(1:nr) )
    end do
    !
  end if
  !
end function common_visc_flx_cv
!
!###############################################################################
!
