module mrf_mod
  !
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq,ncs
  use timer_mod, only : solver_timer
  !
  implicit none
  !
  private
  !
  public :: is_in_mrf
  public :: is_not_in_mrf
  public :: mrf_flx_cv
  public :: mrf_wall_flx_cv
  public :: mrf_src_cv
  public :: switch_vel_ref_frame
  !
  type :: mrf_sign_t
    real(wp) :: add_sub
  end type mrf_sign_t
  !
  type(mrf_sign_t), public, parameter :: abs_to_mrf = mrf_sign_t(-one)
  type(mrf_sign_t), public, parameter :: mrf_to_abs = mrf_sign_t( one)
  !
  logical(lk) :: wall_flux_check_not_completed = true
  logical(lk) :: mrf_wall_fluxes_are_zero = fals
  !
  real(wp), parameter :: wall_default_sign = -one
  real(wp), parameter :: flx_default_sign = -one
  real(wp), parameter :: src_default_sign = one
  real(wp), parameter :: upwind_default_sign = -one
  !
  real(wp), parameter :: switch_sign = one
 !real(wp), parameter :: switch_sign = -one
 !real(wp), parameter :: switch_sign = zero
  !
  real(wp), parameter :: wall_sign =  switch_sign*wall_default_sign
 !real(wp), parameter :: wall_sign = -switch_sign*wall_default_sign
 !real(wp), parameter :: wall_sign = zero*wall_default_sign
  !
  real(wp), parameter :: flx_sign =  switch_sign*flx_default_sign
 !real(wp), parameter :: flx_sign = -switch_sign*flx_default_sign
 !real(wp), parameter :: flx_sign = zero*flx_default_sign
  !
  real(wp), parameter :: src_sign =  switch_sign*src_default_sign
 !real(wp), parameter :: src_sign = -switch_sign*src_default_sign
 !real(wp), parameter :: src_sign = zero*src_default_sign
  !
  real(wp), public, parameter :: upwind_sign =  switch_sign*upwind_default_sign
 !real(wp), public, parameter :: upwind_sign = -switch_sign*upwind_default_sign
 !real(wp), public, parameter :: upwind_sign = zero*upwind_default_sign
  !
  ! #########################################
  ! #####   MODULE GENERIC INTERFACES   #####
  ! #########################################
  !
  interface switch_vel_ref_frame
    module procedure switch_vel_ref_frame_fp, &
                     switch_vel_ref_frame_sp
  end interface switch_vel_ref_frame
  !
contains
!
!###############################################################################
!
pure function is_in_mrf() result(return_value)
  !.. Use Statements ..
  use ovar, only : moving_ref_frame
  !.. Function Result ..
  logical(lk) :: return_value
continue
  return_value = moving_ref_frame%is_active()
end function is_in_mrf
!
!###############################################################################
!
pure function is_not_in_mrf() result(return_value)
  !.. Use Statements ..
  use ovar, only : moving_ref_frame
  !.. Function Result ..
  logical(lk) :: return_value
continue
  return_value = moving_ref_frame%is_not_active()
end function is_not_in_mrf
!
!###############################################################################
!
subroutine mrf_flx_cv(cv,xyz,iflx)
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: cv(:,:)
  real(wp),    intent(in) :: xyz(:,:)
  ! Upon entering this subroutine, iflx should contain the inviscid fluxes.
  real(wp), intent(inout) :: iflx(:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
  !. Local Saved Scalars ..
  integer, save :: idx = -1
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "mrf_flx_cv"
  !
continue
  !
  ! If not in a moving reference frame, return without computing
  ! the flux terms needed for a rotating reference frame.
  !
  if ( is_not_in_mrf() ) return
  !
  call solver_timer(entering_procedure,pname,idx)
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = size(iflx,dim=2)
  !
  if (nr == 3) then
    call mrf_flx_cv_3d(cv,xyz,iflx)
  else if (nr == 2) then
   !call mrf_flx_cv_2d(cv,xyz,iflx)
  else
   !call mrf_flx_cv_1d(cv,xyz,iflx)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  call solver_timer(leaving_procedure,pname,idx)
  !
end subroutine mrf_flx_cv
!
!###############################################################################
!
subroutine mrf_flx_cv_3d(cv,xyz,iflx)
  !
  !.. Use Statements ..
  use ovar, only : moving_ref_frame
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: cv(:,:)
  real(wp),    intent(in) :: xyz(:,:)
  real(wp), intent(inout) :: iflx(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,m,np
  real(wp) :: x,x_r,omegax
  real(wp) :: y,y_r,omegay
  real(wp) :: z,z_r,omegaz
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "mrf_flx_cv_3d"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  np = size(iflx,dim=1)
  !
  call moving_ref_frame%get_rotation(x_r,y_r,z_r)
  !
!#ifdef DEBUG_ON
!!#define DEBUG_MRF_FREESTREAM_TCF
!#ifdef DEBUG_MRF_FREESTREAM_TCF
!block
!  use ovar, only : itcur,current_rk_stage
!  write (341,1) itcur,current_rk_stage
!  write (342,1) itcur,current_rk_stage
!  write (343,1) itcur,current_rk_stage
!  write (347,1) itcur,current_rk_stage
!  write (348,1) itcur,current_rk_stage
!  write (349,1) itcur,current_rk_stage
!  1 format (/,2x,"Iteration = ",i4,5x,"RK stage = ",i1,/,80("-"))
!end block
!#endif
!#endif
  !
  !$ACC PARALLEL LOOP PRESENT(cv,xyz,iflx)
  do n = 1,np
    !
    x = xyz(1,n)
    y = xyz(2,n)
    z = xyz(3,n)
    !
    omegax = y_r*z - z_r*y
    omegay = z_r*x - x_r*z
    omegaz = x_r*y - y_r*x
    !
!#ifdef DEBUG_MRF_FREESTREAM_TCF
!block
!    real(wp) :: fec_x,fmx_x,fmy_x,fmz_x,fee_x
!    real(wp) :: fec_y,fmx_y,fmy_y,fmz_y,fee_y
!    real(wp) :: fec_z,fmx_z,fmy_z,fmz_z,fee_z
!    real(wp) :: r,t
!    !
!    r = sqrt(y*y+z*z)
!    t = atan2(z,y)
!    t = (t + merge(two*pi,zero,t<zero))*radians2degrees
!    !
!    fec_x = flx_sign*cv(n,nec)*omegax
!    fmx_x = flx_sign*cv(n,nmx)*omegax
!    fmy_x = flx_sign*cv(n,nmy)*omegax
!    fmz_x = flx_sign*cv(n,nmz)*omegax
!    fee_x = flx_sign*cv(n,nee)*omegax
!    !
!    fec_y = flx_sign*cv(n,nec)*omegay
!    fmx_y = flx_sign*cv(n,nmx)*omegay
!    fmy_y = flx_sign*cv(n,nmy)*omegay
!    fmz_y = flx_sign*cv(n,nmz)*omegay
!    fee_y = flx_sign*cv(n,nee)*omegay
!    !
!    fec_z = flx_sign*cv(n,nec)*omegaz
!    fmx_z = flx_sign*cv(n,nmx)*omegaz
!    fmy_z = flx_sign*cv(n,nmy)*omegaz
!    fmz_z = flx_sign*cv(n,nmz)*omegaz
!    fee_z = flx_sign*cv(n,nee)*omegaz
!    !
!    write (341,2) n,x,r,t,fec_x,fmx_x,fmy_x,fmz_x,fee_x
!    write (342,2) n,x,r,t,fec_y,fmx_y,fmy_y,fmz_y,fee_y
!    write (343,2) n,x,r,t,fec_z,fmx_z,fmy_z,fmz_z,fee_z
!    !
!    write (347,2) n,x,r,t,(iflx(n,1,m), m=1,nq)
!    write (348,2) n,x,r,t,(iflx(n,2,m), m=1,nq)
!    write (349,2) n,x,r,t,(iflx(n,3,m), m=1,nq)
!    !
!    iflx(n,1,nec) = iflx(n,1,nec) + fec_x
!    iflx(n,1,nmx) = iflx(n,1,nmx) + fmx_x
!    iflx(n,1,nmy) = iflx(n,1,nmy) + fmy_x
!    iflx(n,1,nmz) = iflx(n,1,nmz) + fmz_x
!    iflx(n,1,nee) = iflx(n,1,nee) + fee_x
!    !
!    iflx(n,2,nec) = iflx(n,2,nec) + fec_y
!    iflx(n,2,nmx) = iflx(n,2,nmx) + fmx_y
!    iflx(n,2,nmy) = iflx(n,2,nmy) + fmy_y
!    iflx(n,2,nmz) = iflx(n,2,nmz) + fmz_y
!    iflx(n,2,nee) = iflx(n,2,nee) + fee_y
!    !
!    iflx(n,3,nec) = iflx(n,3,nec) + fec_z
!    iflx(n,3,nmx) = iflx(n,3,nmx) + fmx_z
!    iflx(n,3,nmy) = iflx(n,3,nmy) + fmy_z
!    iflx(n,3,nmz) = iflx(n,3,nmz) + fmz_z
!    iflx(n,3,nee) = iflx(n,3,nee) + fee_z
!    !
!    2 format (i8,3x,3(2x,f12.6),3x,5es20.10)
!end block
!#else
    iflx(n,1,nec) = iflx(n,1,nec) + flx_sign*cv(n,nec)*omegax
    iflx(n,1,nmx) = iflx(n,1,nmx) + flx_sign*cv(n,nmx)*omegax
    iflx(n,1,nmy) = iflx(n,1,nmy) + flx_sign*cv(n,nmy)*omegax
    iflx(n,1,nmz) = iflx(n,1,nmz) + flx_sign*cv(n,nmz)*omegax
    iflx(n,1,nee) = iflx(n,1,nee) + flx_sign*cv(n,nee)*omegax
    !
    iflx(n,2,nec) = iflx(n,2,nec) + flx_sign*cv(n,nec)*omegay
    iflx(n,2,nmx) = iflx(n,2,nmx) + flx_sign*cv(n,nmx)*omegay
    iflx(n,2,nmy) = iflx(n,2,nmy) + flx_sign*cv(n,nmy)*omegay
    iflx(n,2,nmz) = iflx(n,2,nmz) + flx_sign*cv(n,nmz)*omegay
    iflx(n,2,nee) = iflx(n,2,nee) + flx_sign*cv(n,nee)*omegay
    !
    iflx(n,3,nec) = iflx(n,3,nec) + flx_sign*cv(n,nec)*omegaz
    iflx(n,3,nmx) = iflx(n,3,nmx) + flx_sign*cv(n,nmx)*omegaz
    iflx(n,3,nmy) = iflx(n,3,nmy) + flx_sign*cv(n,nmy)*omegaz
    iflx(n,3,nmz) = iflx(n,3,nmz) + flx_sign*cv(n,nmz)*omegaz
    iflx(n,3,nee) = iflx(n,3,nee) + flx_sign*cv(n,nee)*omegaz
!#endif
    !
  end do
  !
  ! Species continuity equations
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(cv,xyz,iflx)
    do m = 1,ncs
      do n = 1,np
        iflx(n,1,m) = iflx(n,1,m) + flx_sign*cv(n,m)*(y_r*xyz(3,n)-z_r*xyz(2,n))
        iflx(n,2,m) = iflx(n,2,m) + flx_sign*cv(n,m)*(z_r*xyz(1,n)-x_r*xyz(3,n))
        iflx(n,3,m) = iflx(n,3,m) + flx_sign*cv(n,m)*(x_r*xyz(2,n)-y_r*xyz(1,n))
      end do
    end do
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine mrf_flx_cv_3d
!
!###############################################################################
!
subroutine mrf_wall_flx_cv(n1,n2,facenrm,facexyz,face_cv,faceflx)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: facexyz(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  ! Upon entering this subroutine, faceflx should
  ! contain the inviscid wall fluxes.
  real(wp), intent(inout) :: faceflx(:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
  !. Local Saved Scalars ..
  integer, save :: idx = -1
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "mrf_wall_flx_cv"
  !
continue
  !
  ! If not in a moving reference frame, return without computing
  ! the wall fluxes needed for a rotating reference frame.
  !
  if ( is_not_in_mrf() ) return
  !
  ! Check to see if the wall fluxes for a rotating reference frame
  ! will actually be needed or if they will always be zero.
  !
!!if (wall_flux_check_not_completed) then
!!  call check_if_wall_fluxes_are_needed(n1,n2,facenrm,facexyz)
!!end if
  !
  ! Return without computing the wall fluxes for a rotating reference frame
  ! if the wall fluxes were previously determined to all be zero.
  !
  if (mrf_wall_fluxes_are_zero) return
  !
  ! Continue on to compute the wall fluxes for a rotating reference frame
  !
  call solver_timer(entering_procedure,pname,idx)
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = size(facexyz,dim=1)
  !
  if (nr == 3) then
    call mrf_wall_flx_cv_3d(n1,n2,facenrm,facexyz,face_cv,faceflx)
  else if (nr == 2) then
   !call mrf_wall_flx_cv_2d(n1,n2,facenrm,facexyz,face_cv,faceflx)
  else
   !call mrf_wall_flx_cv_1d(n1,n2,facenrm,facexyz,face_cv,faceflx)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  call solver_timer(leaving_procedure,pname,idx)
  !
end subroutine mrf_wall_flx_cv
!
!###############################################################################
!
subroutine mrf_wall_flx_cv_3d(n1,n2,facenrm,facexyz,face_cv,faceflx)
  !
  !.. Use Statements ..
  use ovar, only : moving_ref_frame
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: facexyz(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  ! Upon entering this subroutine, faceflx should
  ! contain the inviscid wall fluxes.
  real(wp), intent(inout) :: faceflx(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,m
  real(wp) :: x_r,y_r,z_r,omega_nrm
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "mrf_wall_flx_cv_3d"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  call moving_ref_frame%get_rotation(x_r,y_r,z_r)
  !
!#ifdef DEBUG_MRF_FREESTREAM_TCF
!block
!  use ovar, only : itcur,current_rk_stage
!  write (351,1) itcur,current_rk_stage
!  write (357,1) itcur,current_rk_stage
!  1 format (/,2x,"Iteration = ",i4,5x,"RK stage = ",i1,/,80("-"))
!end block
!#endif
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,facexyz,face_cv,faceflx)
  do n = n1,n2
    !
    ! NOTE: We dont normalize facenrm to get the unit normal to compute
    !       omega_nrm, because we need the flux term we add to faceflx
    !       to include the area of the face point since this is ultimately
    !       computing an integral of the flux over the face.
    !
    omega_nrm = facenrm(1,n)*(y_r*facexyz(3,n) - z_r*facexyz(2,n)) + &
                facenrm(2,n)*(z_r*facexyz(1,n) - x_r*facexyz(3,n)) + &
                facenrm(3,n)*(x_r*facexyz(2,n) - y_r*facexyz(1,n))
    !
!#ifdef DEBUG_MRF_FREESTREAM_TCF
!block
!    real(wp) :: fec,fmx,fmy,fmz,fee
!    real(wp) :: x,y,z,r,t
!    !
!    x = facexyz(1,n)
!    y = facexyz(2,n)
!    z = facexyz(3,n)
!    r = sqrt(y*y+z*z)
!    t = atan2(z,y)
!    t = (t + merge(two*pi,zero,t<zero))*radians2degrees
!    !
!    fec = wall_sign*face_cv(n,nec,2)*omega_nrm
!    fmx = wall_sign*face_cv(n,nmx,2)*omega_nrm
!    fmy = wall_sign*face_cv(n,nmy,2)*omega_nrm
!    fmz = wall_sign*face_cv(n,nmz,2)*omega_nrm
!    fee = wall_sign*face_cv(n,nee,2)*omega_nrm
!    !
!    write (351,2) n,x,r,t,fec,fmx,fmy,fmz,fee
!    !
!    write (357,2) n,x,r,t,(faceflx(n,m), m=1,nq)
!    !
!    faceflx(n,nec) = faceflx(n,nec) + fec
!    faceflx(n,nmx) = faceflx(n,nmx) + fmx
!    faceflx(n,nmy) = faceflx(n,nmy) + fmy
!    faceflx(n,nmz) = faceflx(n,nmz) + fmz
!    faceflx(n,nee) = faceflx(n,nee) + fee
!    !
!    2 format (i8,3x,3(2x,f12.6),3x,5es20.10)
!end block
!#else
    faceflx(n,nec) = faceflx(n,nec) + wall_sign*face_cv(n,nec,2)*omega_nrm
    faceflx(n,nmx) = faceflx(n,nmx) + wall_sign*face_cv(n,nmx,2)*omega_nrm
    faceflx(n,nmy) = faceflx(n,nmy) + wall_sign*face_cv(n,nmy,2)*omega_nrm
    faceflx(n,nmz) = faceflx(n,nmz) + wall_sign*face_cv(n,nmz,2)*omega_nrm
    faceflx(n,nee) = faceflx(n,nee) + wall_sign*face_cv(n,nee,2)*omega_nrm
!#endif
    !
  end do
  !
  ! Species continuity equations
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(facenrm,facexyz,face_cv,faceflx)
    do m = 1,ncs
      do n = n1,n2
        !
        ! NOTE: We dont normalize facenrm to get the unit normal to compute
        !       omega_nrm, because we need the flux term we add to faceflx
        !       to include the area of the face point since this is ultimately
        !       computing an integral of the flux over the face.
        !
        omega_nrm = facenrm(1,n)*(y_r*facexyz(3,n) - z_r*facexyz(2,n)) + &
                    facenrm(2,n)*(z_r*facexyz(1,n) - x_r*facexyz(3,n)) + &
                    facenrm(3,n)*(x_r*facexyz(2,n) - y_r*facexyz(1,n))
        !
        faceflx(n,m) = faceflx(n,m) + wall_sign*face_cv(n,m,2)*omega_nrm
        !
      end do
    end do
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine mrf_wall_flx_cv_3d
!
!###############################################################################
!
subroutine check_if_wall_fluxes_are_needed(n1,n2,facenrm,facexyz)
  !
  !.. Use Statements ..
  use ovar, only : moving_ref_frame
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: facexyz(:,:)
  !
  !.. Local Scalars ..
  integer  :: nr,n
  real(wp) :: x_r,y_r,z_r,omega_nrm,max_omega_nrm
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_if_wall_fluxes_are_needed"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  wall_flux_check_not_completed = fals
  !
  ! Get the full speed (i.e., after any spinup) rotation rates
  !
  call moving_ref_frame%get_rotation(x_r,y_r,z_r,full_speed=true)
  !
  nr = size(facexyz,dim=1)
  !
  ! Find the maximum value of dot_product( omega , face normal )
  ! across all the face/flux points.
  !
  max_omega_nrm = zero
  if (nr == 3) then
    !
    !$ACC PARALLEL LOOP REDUCTION(max:max_omega_nrm) &
    !$ACC&              PRESENT(facenrm,facexyz)
    do n = n1,n2
      !
      omega_nrm = facenrm(1,n)*(y_r*facexyz(3,n) - z_r*facexyz(2,n)) + &
                  facenrm(2,n)*(z_r*facexyz(1,n) - x_r*facexyz(3,n)) + &
                  facenrm(3,n)*(x_r*facexyz(2,n) - y_r*facexyz(1,n))
      !
!#ifdef DEBUG_MRF_FREESTREAM_TCF
!block
!      integer  :: l
!      real(wp) :: x,y,z,r,t
!      x = facexyz(1,n)
!      y = facexyz(2,n)
!      z = facexyz(3,n)
!      r = sqrt(y*y+z*z)
!      t = (atan2(z,y)+merge(two*pi,zero,atan2(z,y)<zero))*radians2degrees
!      write (381,2) n,x,r,t,omega_nrm,(facenrm(l,n),l=1,3), &
!                    facenrm(2,n),facexyz(3,n), &
!                    facenrm(3,n),facexyz(2,n), &
!                    -facenrm(2,n)*x_r*facexyz(3,n), &
!                     facenrm(3,n)*x_r*facexyz(2,n)
!      2 format (i8,3x,3(1x,f12.6),3x,es20.10,3x,3es20.10,3(3x,2es20.10))
!end block
!#endif
      max_omega_nrm = max( max_omega_nrm , abs(omega_nrm) )
      !
    end do
    !
  else if (nr == 2) then
    !
  else
    !
  end if
!#ifdef DEBUG_MRF_FREESTREAM_TCF
!  write (iout,1) mypnum,"max_omega_nrm = ",max_omega_nrm
!  1 format ("CPU #",i5," ",a,g15.8)
!#endif
  !
  ! Set mrf_wall_fluxes_are_zero to true if the maximum omega_nrm was smaller
  ! than eps12. This will prevent the wasted computations of computing the
  ! MRF wall fluxes since they will always be zero.
  !
  if (abs(max_omega_nrm) < eps12) then
    mrf_wall_fluxes_are_zero = true
  end if
 !mrf_wall_fluxes_are_zero = true
! write (iout,1) mypnum,"mrf_wall_fluxes_are_zero = ",mrf_wall_fluxes_are_zero
  !
! call debug_mrf_walls(n1,n2,x_r,y_r,z_r,facenrm,facexyz)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine check_if_wall_fluxes_are_needed
!
!###############################################################################
!
subroutine debug_mrf_walls(n1,n2,x_r,y_r,z_r,facenrm,facexyz)
  !
  !.. Use Statements ..
  use geovar, only : nfbnd,face,bface
  use ovar,   only : aref,grid_scaling_factor
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: n1,n2
  real(wp), intent(in) :: x_r,y_r,z_r
  real(wp), intent(in) :: facenrm(:,:)
  real(wp), intent(in) :: facexyz(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nf,f1,f2
  real(wp) :: r,t,ds,omega_nrm
  real(wp) :: x,x_rpm,rnx,unx,omegax
  real(wp) :: y,y_rpm,rny,uny,omegay
  real(wp) :: z,z_rpm,rnz,unz,omegaz
  !
continue
  !
  x_rpm = x_r*aref/rpm_to_mps
  y_rpm = y_r*aref/rpm_to_mps
  z_rpm = z_r*aref/rpm_to_mps
  !
  write (1000+mypnum,1) n1,n2,x_r,x_rpm,y_r,y_rpm,z_r,z_rpm
  !
  bnd_face_loop: &
  do nf = 1,nfbnd
    !
    if (face(nf)%beg_fp > n2) exit bnd_face_loop
    !
    f1 = face(nf)%beg_fp
    f2 = face(nf)%end_fp
    !
    do n = f1,f2
      !
      x = facexyz(1,n)
      y = facexyz(2,n)
      z = facexyz(3,n)
      !
      rnx = facenrm(1,n)
      rny = facenrm(2,n)
      rnz = facenrm(3,n)
      ds = sqrt(rnx*rnx + rny*rny + rnz*rnz)
      unx = rnx/ds
      uny = rny/ds
      unz = rnz/ds
      !
      ! [x_r,y_r,z_r] are in rad/s, [x,y,z] are in meters,
      ! so [omegax,omegay,omegaz] are in m/s (and nondimensionalized by aref)
      !
      omegax = y_r*z - z_r*y
      omegay = z_r*x - x_r*z
      omegaz = x_r*y - y_r*x
      !
      ! [rnx,rny,rnz] are in m**2
      !
      omega_nrm = rnx*omegax + rny*omegay + rnz*omegaz
      !
      x = x/grid_scaling_factor
      y = y/grid_scaling_factor
      z = z/grid_scaling_factor
      !
      r = sqrt(y*y+z*z)
      t = (atan2(z,y)+merge(two*pi,zero,atan2(z,y)<zero))*radians2degrees
      !
      write (1000+mypnum,2) n,bface(1,nf),x,y,z,r,t,unx,uny,unz,rnx,rny,rnz, &
                            omegax,omegay,omegaz,omega_nrm
      !
    end do
    !
  end do bnd_face_loop
  !
  1 format ("[n1,n2] = [",i0,",",i0,"]",/, &
            "x_r = ",es20.10," ( ",f12.6," RPM )",/, &
            "y_r = ",es20.10," ( ",f12.6," RPM )",/, &
            "z_r = ",es20.10," ( ",f12.6," RPM )",/, &
            100("="))
  2 format (i8,1x,i3,3x,3(1x,f10.6),2x,2(1x,f10.6),3(3x,3es14.6),3x,es20.10)
  !
end subroutine debug_mrf_walls
!
!###############################################################################
!
subroutine mrf_src_cv()
  !
  !.. Use Statements ..
  use geovar,      only : nr
  use flowvar,     only : usp,residual
  use metrics_mod, only : jac
  !
  !. Local Saved Scalars ..
  integer, save :: idx = -1
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "mrf_src_cv"
  !
continue
  !
  ! If not in a moving reference frame, return without computing
  ! the source terms needed for a rotating reference frame.
  !
  if ( is_not_in_mrf() ) return
  !
  call solver_timer(entering_procedure,pname,idx)
  !
  call debug_timer(entering_procedure,pname)
  !
  if (nr == 3) then
    call mrf_src_cv_3d(jac,usp,residual)
  else if (nr == 2) then
   !call mrf_src_cv_2d(jac,usp,residual)
  else
   !call mrf_src_cv_1d(jac,usp,residual)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  call solver_timer(leaving_procedure,pname,idx)
  !
end subroutine mrf_src_cv
!
!###############################################################################
!
subroutine mrf_src_cv_3d(jac,cv,rhs)
  !
  !.. Use Statements ..
  use ovar, only : moving_ref_frame
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: jac(:)
  real(wp),    intent(in) :: cv(:,:)
  real(wp), intent(inout) :: rhs(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,np
  real(wp) :: x,x_r,omegax
  real(wp) :: y,y_r,omegay
  real(wp) :: z,z_r,omegaz
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "mrf_src_cv_3d"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  np = size(rhs,dim=1)
  !
  call moving_ref_frame%get_rotation(x_r,y_r,z_r)
  !
!#ifdef DEBUG_MRF_FREESTREAM_TCF
!block
!  use ovar, only : itcur,current_rk_stage
!  write (331,1) itcur,current_rk_stage
!  1 format (/,2x,"Iteration = ",i4,5x,"RK stage = ",i1,/,80("-"))
!end block
!#endif
  !
  ! NOTE: In the governing equations, the source terms have a leading negative
  !       sign when they are located on the RHS of the equation. However, the
  !       rhs/residual array is actually formed as if it is on the LHS of the
  !       equations. Just before updating each time step, the rhs/residual
  !       array is multiplied by -1 to move it to the RHS of the equation.
  !       Because of this, we need to remove the negative sign on these source
  !       terms since we are essentially moving them to the LHS to add to the
  !       rhs/residual array.
  !
  !$ACC PARALLEL LOOP PRESENT(cv,rhs)
  do n = 1,np
    !
!#ifdef DEBUG_MRF_FREESTREAM_TCF
!block
!    use geovar, only : xyz
!    real(wp) :: sec,smx,smy,smz,see
!    real(wp) :: x,y,z,r,t
!    !
!    x = xyz(1,n)
!    y = xyz(2,n)
!    z = xyz(3,n)
!    r = sqrt(y*y+z*z)
!    t = atan2(z,y)
!    t = (t + merge(two*pi,zero,t<zero))*radians2degrees
!    !
!    sec = zero
!    smx = src_sign*(y_r*cv(n,nmz) - z_r*cv(n,nmy))*jac(n)
!    smy = src_sign*(z_r*cv(n,nmx) - x_r*cv(n,nmz))*jac(n)
!    smz = src_sign*(x_r*cv(n,nmy) - y_r*cv(n,nmx))*jac(n)
!    see = zero
!    !
!    write (331,2) n,x,r,t,sec,smx,smy,smz,see
!    !
!    rhs(n,nmx) = rhs(n,nmx) + smx
!    rhs(n,nmy) = rhs(n,nmy) + smy
!    rhs(n,nmz) = rhs(n,nmz) + smz
!    !
!    2 format (i8,3x,3(2x,f12.6),3x,5es20.10)
!end block
!#else
    rhs(n,nmx) = rhs(n,nmx) + src_sign*(y_r*cv(n,nmz) - z_r*cv(n,nmy))*jac(n)
    rhs(n,nmy) = rhs(n,nmy) + src_sign*(z_r*cv(n,nmx) - x_r*cv(n,nmz))*jac(n)
    rhs(n,nmz) = rhs(n,nmz) + src_sign*(x_r*cv(n,nmy) - y_r*cv(n,nmx))*jac(n)
!#endif
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine mrf_src_cv_3d
!
!###############################################################################
!
pure subroutine switch_vel_ref_frame_fp(n1,n2,cv,fxyz,mrf_sign)
  !
  !.. Use Statements ..
  use ovar, only : moving_ref_frame
  !
  ! This subroutine switches the momentum terms in the cv array from the
  ! absolute/inertial reference frame to the rotating reference frame, or
  ! vice versa depending on the value of argument mrf_sign%add_sub.
  !
  !.. Formal Arguments ..
  integer,          intent(in) :: n1,n2
  real(wp),      intent(inout) :: cv(:,:,:)
  real(wp),         intent(in) :: fxyz(:,:)
  type(mrf_sign_t), intent(in) :: mrf_sign ! this is either +1 or -1
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: add_sub,x_r,y_r,z_r
  real(wp) :: omegax,omegay,omegaz
  !
continue
  !
  ! NOTE: The argument mrf_sign%add_sub is a coefficient that determines
  !       if we are either adding/subtracting the rotating reference frame
  !       angular velocity to/from the momentum value stored in cv.
  !         mrf_sign%add_sub = +1  ->  this is converting rotating to absolute
  !         mrf_sign%add_sub = -1  ->  this is converting absolute to rotating
  !
  nr = nme-nmb+1
  !
  add_sub = mrf_sign%add_sub
  !
  call moving_ref_frame%get_rotation(x_r,y_r,z_r)
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(cv,fxyz)
    do n = n1,n2
      !
      omegax = y_r*fxyz(3,n) - z_r*fxyz(2,n)
      omegay = z_r*fxyz(1,n) - x_r*fxyz(3,n)
      omegaz = x_r*fxyz(2,n) - y_r*fxyz(1,n)
      !
      ! Left side of face
      cv(n,nmx,1) = cv(n,nmx,1) + add_sub*cv(n,nec,1)*omegax
      cv(n,nmy,1) = cv(n,nmy,1) + add_sub*cv(n,nec,1)*omegay
      cv(n,nmz,1) = cv(n,nmz,1) + add_sub*cv(n,nec,1)*omegaz
      !
      ! Right side of face
      cv(n,nmx,2) = cv(n,nmx,2) + add_sub*cv(n,nec,2)*omegax
      cv(n,nmy,2) = cv(n,nmy,2) + add_sub*cv(n,nec,2)*omegay
      cv(n,nmz,2) = cv(n,nmz,2) + add_sub*cv(n,nec,2)*omegaz
      !
    end do
  else if (nr == 2) then
   !!$ACC PARALLEL LOOP PRESENT(cv,fxyz)
   !do n = n1,n2
   !  !
   !  omegax = y_r*fxyz(1,n) - x_r*fxyz(2,n)
   !  omegay = x_r*fxyz(2,n) - y_r*fxyz(1,n)
   !  !
   !  ! Left side of face
   !  cv(n,nmx,1) = cv(n,nmx,1) + add_sub*cv(n,nec,1)*omegax
   !  cv(n,nmy,1) = cv(n,nmy,1) + add_sub*cv(n,nec,1)*omegay
   !  !
   !  ! Right side of face
   !  cv(n,nmx,2) = cv(n,nmx,2) + add_sub*cv(n,nec,2)*omegax
   !  cv(n,nmy,2) = cv(n,nmy,2) + add_sub*cv(n,nec,2)*omegay
   !  !
   !end do
  else
   !!$ACC PARALLEL LOOP PRESENT(cv,fxyz)
   !do n = n1,n2
   !  !
   !  omegax = x_r*fxyz(1,n)
   !  !
   !  ! Left side of face
   !  cv(n,nmx,1) = cv(n,nmx,1) + add_sub*cv(n,nec,1)*omegax
   !  !
   !  ! Right side of face
   !  cv(n,nmx,2) = cv(n,nmx,2) + add_sub*cv(n,nec,2)*omegax
   !  !
   !end do
  end if
  !
end subroutine switch_vel_ref_frame_fp
!
!###############################################################################
!
pure subroutine switch_vel_ref_frame_sp(n1,n2,cv,xyz,mrf_sign)
  !
  ! This subroutine switches the momentum terms in the cv array from the
  ! absolute/inertial reference frame to the rotating reference frame, or
  ! vice versa depending on the value of argument mrf_sign%add_sub.
  !
  !.. Use Statements ..
  use ovar, only : moving_ref_frame
  !
  !.. Formal Arguments ..
  integer,          intent(in) :: n1,n2
  real(wp),      intent(inout) :: cv(:,:)
  real(wp),         intent(in) :: xyz(:,:)
  type(mrf_sign_t), intent(in) :: mrf_sign ! this is either +1 or -1
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: add_sub,x_r,y_r,z_r
  !
continue
  !
  ! NOTE: The argument mrf_sign%add_sub is a coefficient that determines
  !       if we are either adding/subtracting the rotating reference frame
  !       angular velocity to/from the momentum value stored in cv.
  !         mrf_sign%add_sub = +1  ->  this is converting rotating to absolute
  !         mrf_sign%add_sub = -1  ->  this is converting absolute to rotating
  !
  nr = nme-nmb+1
  !
  add_sub = mrf_sign%add_sub
  !
  call moving_ref_frame%get_rotation(x_r,y_r,z_r)
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(cv,xyz)
    do n = n1,n2
      !
      cv(n,nmx) = cv(n,nmx) + add_sub*cv(n,nec)*(y_r*xyz(3,n) - z_r*xyz(2,n))
      cv(n,nmy) = cv(n,nmy) + add_sub*cv(n,nec)*(z_r*xyz(1,n) - x_r*xyz(3,n))
      cv(n,nmz) = cv(n,nmz) + add_sub*cv(n,nec)*(x_r*xyz(2,n) - y_r*xyz(1,n))
      !
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(cv,xyz)
    do n = n1,n2
      !
     !cv(n,nmx) = cv(n,nmx) + add_sub*cv(n,nec)*(y_r*xyz(1,n) - x_r*xyz(2,n))
     !cv(n,nmy) = cv(n,nmy) + add_sub*cv(n,nec)*(x_r*xyz(2,n) - y_r*xyz(1,n))
      !
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(cv,xyz)
    do n = n1,n2
      !
     !cv(n,nmx) = cv(n,nmx) + add_sub*cv(n,nec)*x_r*xyz(1,n)
      !
    end do
  end if
  !
end subroutine switch_vel_ref_frame_sp
!
!###############################################################################
!
end module mrf_mod
!
!###############################################################################
!
