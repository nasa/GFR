module mms_mod
  !
  ! Module containing procedures for computing the
  ! Method of Manufactured Solution source terms (MMS)
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq,ntk,ntl
  !
  implicit none
  !
  private
  !
  !.. Public procedures
  public :: mms_src
  public :: initialize_mms
  public :: mms_solution_at_xyz
  !
  !
  !.. Local Module Variables ..
  real(wp), save :: dens0,densx,densy,densz,densxy,densxz,densyz
  real(wp), save :: pres0,presx,presy,presz,presxy,presxz,presyz
  real(wp), save :: uvel0,uvelx,uvely,uvelz,uvelxy,uvelxz,uvelyz
  real(wp), save :: vvel0,vvelx,vvely,vvelz,vvelxy,vvelxz,vvelyz
  real(wp), save :: wvel0,wvelx,wvely,wvelz,wvelxy,wvelxz,wvelyz
  real(wp), save :: adx,ady,adz,adxy,adxz,adyz
  real(wp), save :: apx,apy,apz,apxy,apxz,apyz
  real(wp), save :: aux,auy,auz,auxy,auxz,auyz
  real(wp), save :: avx,avy,avz,avxy,avxz,avyz
  real(wp), save :: awx,awy,awz,awxy,awxz,awyz
  !
  real(wp), save :: rmu = zero
  real(wp), save :: L = one
  real(wp), save :: R
  real(wp), save :: gm1
  !
contains
!
!###############################################################################
!
subroutine initialize_mms()
  !
  ! Routine to apply source terms for the method of manufactured solutions
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar, only : gam,rgasref,muref
  use ovar, only : governing_equations,mms_opt
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_mms"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  R = rgasref
  !
  gm1 = gam - one
  !
  if (governing_equations == NavierStokes_Eqns) then
    rmu = muref
  end if
  !
  if (mms_opt == 1) then
    !
    uvel0 = 70.0_wp
    uvelx = four
    uvely = -twelve
    uvelz = ten
    uvelxy = seven
    uvelxz = seven
    uvelyz = seven
    !
    vvel0 = 90.0_wp
    vvelx = -20.0_wp
    vvely = four
    vvelz = -17.5_wp
    vvelxy = -11.0_wp
    vvelxz = -11.0_wp
    vvelyz = -11.0_wp
    !
    wvel0 = 8.0e+1_wp
    wvelx = -eight
    wvely = -four
    wvelz = -3.75_wp
    wvelxy = -two
    wvelxz = -two
    wvelyz = -two
    !
    dens0 = one
    densx = one/ten
    densy = 0.15_wp
    densz = one8
    densxy = 0.08_wp
    densxz = 0.08_wp
    densyz = 0.08_wp
    !
    pres0 = 1.0e+5_wp
    presx = -0.3e+5_wp
    presy = 0.2e+5_wp
    presz = -0.5e+4_wp
    presxy = -0.25e+5_wp
    presxz = -0.25e+5_wp
    presyz = +0.25e+5_wp
    !
    adx = three4
    ady = one
    adz = seven/eight
    adxy = five/four
    adxz = five/four
    adyz = five/four
    !
    apx = one
    apy = five/four
    apz = nine/eight
    apxy = three4
    apxz = three4
    apyz = three4
    !
    aux = five/three
    auy = three/two
    auz = 1.58_wp
    auxy = three/five
    auxz = three/five
    auyz = three/five
    !
    avx = three/two
    avy = one
    avz = five/four
    avxy = nine/ten
    avxz = nine/ten
    avyz = nine/ten
    !
    awx = 1.58_wp
    awy = five/four
    awz = 1.415_wp
    awxy = three4
    awxz = three4
    awyz = three4
    !
  end if
  !
  if (mms_opt == 2 .or. mms_opt == 3 .or. mms_opt == 4) then
    !
    uvel0 = 800.0_wp
    uvelx = 50.0_wp
    uvely = -30.0_wp
    uvelz = ten
    uvelxy = seven
    uvelxz = seven
    uvelyz = seven
    !
    vvel0 = 800.0_wp
    vvelx = -75.0_wp
    vvely = 40.0_wp
    vvelz = -17.5_wp
    vvelxy = -11.0_wp
    vvelxz = -11.0_wp
    vvelyz = -11.0_wp
    !
    wvel0 = 800.0_wp
    wvelx = -12.5_wp
    wvely = five
    wvelz = -3.75_wp
    wvelxy = -two
    wvelxz = -two
    wvelyz = -two
    !
    dens0 = one
    densx = 0.15_wp
    densy = -one/ten
    densz = one8
    densxy = 0.08_wp
    densxz = 0.08_wp
    densyz = 0.08_wp
    !
    pres0 = 1.0e+5_wp
    presx = 0.2e+5_wp
    presy = 0.5e+5_wp
    presz = -0.5e+4_wp
    presxy = -0.25e+5_wp
    presxz = -0.25e+5_wp
    presyz = +0.25e+5_wp
    !
    adx = one
    ady = half
    adz = seven/eight
    adxy = five/four
    adxz = five/four
    adyz = five/four
    !
    apx = two
    apy = one
    apz = nine/eight
    apxy = three4
    apxz = three4
    apyz = three4
    !
    aux = three/two
    auy = three/five
    auz = 1.58_wp
    auxy = three/five
    auxz = three/five
    auyz = three/five
    !
    avx = half
    avy = two3
    avz = five/four
    avxy = nine/ten
    avxz = nine/ten
    avyz = nine/ten
    !
    awx = 1.58_wp
    awy = five/four
    awz = 1.415_wp
    awxy = three4
    awxz = three4
    awyz = three4
    !
    if (mms_opt == 3) then
      !
      uvely = zero
      uvelxy = zero
      !
      vvel0 = zero
      vvelx = zero
      vvely = zero
      vvelxy = zero
      !
      densy = zero
      densxy = zero
      !
      presy = zero
      presxy = zero
      !
      ady = zero
      adxy = zero
      !
      apy = zero
      apxy = zero
      !
      auy = zero
      auxy = zero
      !
      avx = zero
      avy = zero
      avxy = zero
      !
    end if
    !
    if (mms_opt == 4) then
      !
      uvel0 = zero
      uvelx = zero
      uvely = zero
      uvelxy = zero
      !
      vvelx = zero
      vvelxy = zero
      !
      densx = zero
      densxy = zero
      !
      presx = zero
      presxy = zero
      !
      adx = zero
      adxy = zero
      !
      apx = zero
      apxy = zero
      !
      aux = zero
      auy = zero
      auxy = zero
      !
      avx = zero
      avxy = zero
      !
    end if
    !
  end if
  !
  if (nr <= 2) then
    !
    uvelz = zero
    uvelxz = zero
    uvelyz = zero
    !
    vvelz = zero
    vvelxz = zero
    vvelyz = zero
    !
    wvel0 = zero
    wvelx = zero
    wvely = zero
    wvelz = zero
    wvelxy = zero
    wvelxz = zero
    wvelyz = zero
    !
    densz = zero
    densxz = zero
    densyz = zero
    !
    presz = zero
    presxz = zero
    presyz = zero
    !
    adz = zero
    adxz = zero
    adyz = zero
    !
    apz = zero
    apxz = zero
    apyz = zero
    !
    auz = zero
    auxz = zero
    auyz = zero
    !
    avz = zero
    avxz = zero
    avyz = zero
    !
    awx = zero
    awy = zero
    awz = zero
    awxy = zero
    awxz = zero
    awyz = zero
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine initialize_mms
!
!###############################################################################
!
pure function mms_solution_at_xyz(xyz,return_nondim,return_cv) &
                           result(return_value)
  !
  ! Function to evaluate the exact MMS solution at a set of coordinates
  !
  !.. Use Statements ..
  use ovar, only : rhoref,aref
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xyz
  !
  !.. Optional Arguments ..
  logical(lk), optional,      intent(in) :: return_nondim
  logical(lk), optional,      intent(in) :: return_cv
  !
  !.. Function Result ..
  real(wp), dimension(1:nq) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: x,y,z
  !
continue
  !
  x = xyz(1)
  y = zero
  z = zero
  !
  if (size(xyz) >= 2) y = xyz(2)
  if (size(xyz) >= 3) z = xyz(3)
  !
  return_value(nec) = dens0 + densx*sin((adx*pi*x)/L) &
                            + densy*cos((ady*pi*y)/L) &
                            + densz*sin((adz*pi*z)/L) &
                            + densxy*cos((adxy*pi*x*y)/L**2) &
                            + densxz*sin((adxz*pi*x*z)/L**2) &
                            + densyz*cos((adyz*pi*y*z)/L**2)
  !
  return_value(nee) = pres0 + presx*cos((apx*pi*x)/L) &
                            + presy*sin((apy*pi*y)/L) &
                            + presz*cos((apz*pi*z)/L) &
                            + presxy*sin((apxy*pi*x*y)/L**2) &
                            + presxz*cos((apxz*pi*x*z)/L**2) &
                            + presyz*cos((apyz*pi*y*z)/L**2)
  !
  return_value(nmx) = uvel0 + uvelx*sin((aux*pi*x)/L) &
                            + uvely*cos((auy*pi*y)/L) &
                            + uvelz*cos((auz*pi*z)/L) &
                            + uvelxy*cos((auxy*pi*x*y)/L**2) &
                            + uvelxz*sin((auxz*pi*x*z)/L**2) &
                            + uvelyz*sin((auyz*pi*y*z)/L**2)
  !
  if (size(xyz) >= 2) then
    return_value(nmy) = vvel0 + vvelx*cos((avx*pi*x)/L) &
                              + vvely*sin((avy*pi*y)/L) &
                              + vvelz*sin((avz*pi*z)/L) &
                              + vvelxy*cos((avxy*pi*x*y)/L**2) &
                              + vvelxz*cos((avxz*pi*x*z)/L**2) &
                              + vvelyz*sin((avyz*pi*y*z)/L**2)
  end if
  !
  if (size(xyz) >= 3) then
    return_value(nmz) = wvel0 + wvelx*sin((awx*pi*x)/L) &
                              + wvely*sin((awy*pi*y)/L) &
                              + wvelz*cos((awz*pi*z)/L) &
                              + wvelxy*sin((awxy*pi*x*y)/L**2) &
                              + wvelxz*sin((awxz*pi*x*z)/L**2) &
                              + wvelyz*cos((awyz*pi*y*z)/L**2)
  end if
  !
  if (present(return_nondim)) then
    if (return_nondim) then
      return_value(nec) = return_value(nec) / rhoref
      return_value(nmb:nme) = return_value(nmb:nme) / aref
      return_value(nee) = return_value(nee) / (rhoref * aref**2)
    end if
  end if
  !
  if (present(return_cv)) then
    if (return_cv) then
      return_value(1:nq) = vsp2u_sp( return_value(1:nq) )
    end if
  end if
  !
end function mms_solution_at_xyz
!
!###############################################################################
!
subroutine mms_src()
  !
  ! Routine to apply source terms for the method of manufactured solutions
  !
  !.. Use Statements ..
  use geovar, only : nr,xyz,ncell,cell
  use ovar, only : gam,rhoref,aref,tmin,tmax,Pr
  use flowvar, only : residual
  use metrics_mod, only : metrics_dt
  !
  !.. Local Scalars ..
  integer :: k,n,nc,n1,n2
  !
  real(wp) :: x,y,z
  real(wp) :: hc
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nq) :: mms_sol
  real(wp), dimension(1:nq+1) :: var_min,var_max
  !
  !.. Local Parameters ..
 !real(wp), parameter :: add_or_subtract =  one
  real(wp), parameter :: add_or_subtract = -one
  character(len=*), parameter :: pname = "mms_src"
 !integer, parameter :: jac_exp = 0
  integer, parameter :: jac_exp = 1
 !integer, parameter :: jac_exp = 2
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  var_min(:) = +huge(zero)
  var_max(:) = -huge(zero)
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
      mms_sol(1:nq) = mms_solution_at_xyz( xyz(1:nr,n) )
      !
      var_min(1:nq) = min(mms_sol(1:nq),var_min(1:nq))
      var_max(1:nq) = min(mms_sol(1:nq),var_max(1:nq))
      !
      var_min(nq+1) = min(mms_sol(nee)/(mms_sol(nec)*R),var_min(nq+1))
      var_max(nq+1) = max(mms_sol(nee)/(mms_sol(nec)*R),var_max(nq+1))
      !
      x = xyz(1,n)
      y = zero
      z = zero
      !
      if (nr >= 2) y = xyz(2,n)
      if (nr >= 3) z = xyz(3,n)
!
! #######################################################################
! ########################  Continuity equation  ########################
! #######################################################################
!
      hc = ((aux*Pi*uvelx*Cos((aux*Pi*x)/L))/L +                    &
        (auxz*Pi*uvelxz*z*Cos((auxz*Pi*x*z)/L**2))/L**2 -           &
        (auxy*Pi*uvelxy*y*Sin((auxy*Pi*x*y)/L**2))/L**2)*           &
      (dens0 + densy*Cos((ady*Pi*y)/L) +                            &
        densxy*Cos((adxy*Pi*x*y)/L**2) +                            &
        densyz*Cos((adyz*Pi*y*z)/L**2) + densx*Sin((adx*Pi*x)/L) +  &
        densz*Sin((adz*Pi*z)/L) + densxz*Sin((adxz*Pi*x*z)/L**2)) + &
     ((avy*Pi*vvely*Cos((avy*Pi*y)/L))/L +                          &
        (avyz*Pi*vvelyz*z*Cos((avyz*Pi*y*z)/L**2))/L**2 -           &
        (avxy*Pi*vvelxy*x*Sin((avxy*Pi*x*y)/L**2))/L**2)*           &
      (dens0 + densy*Cos((ady*Pi*y)/L) +                            &
        densxy*Cos((adxy*Pi*x*y)/L**2) +                            &
        densyz*Cos((adyz*Pi*y*z)/L**2) + densx*Sin((adx*Pi*x)/L) +  &
        densz*Sin((adz*Pi*z)/L) + densxz*Sin((adxz*Pi*x*z)/L**2)) + &
     (wvel0 + wvelz*Cos((awz*Pi*z)/L) +                             &
        wvelyz*Cos((awyz*Pi*y*z)/L**2) + wvelx*Sin((awx*Pi*x)/L) +  &
        wvely*Sin((awy*Pi*y)/L) + wvelxy*Sin((awxy*Pi*x*y)/L**2) +  &
        wvelxz*Sin((awxz*Pi*x*z)/L**2))*                            &
      ((adz*densz*Pi*Cos((adz*Pi*z)/L))/L +                         &
        (adxz*densxz*Pi*x*Cos((adxz*Pi*x*z)/L**2))/L**2 -           &
        (adyz*densyz*Pi*y*Sin((adyz*Pi*y*z)/L**2))/L**2) +          &
     ((adx*densx*Pi*Cos((adx*Pi*x)/L))/L +                          &
        (adxz*densxz*Pi*z*Cos((adxz*Pi*x*z)/L**2))/L**2 -           &
        (adxy*densxy*Pi*y*Sin((adxy*Pi*x*y)/L**2))/L**2)*           &
      (uvel0 + uvely*Cos((auy*Pi*y)/L) +                            &
        uvelxy*Cos((auxy*Pi*x*y)/L**2) + uvelz*Cos((auz*Pi*z)/L) +  &
        uvelx*Sin((aux*Pi*x)/L) + uvelxz*Sin((auxz*Pi*x*z)/L**2) +  &
        uvelyz*Sin((auyz*Pi*y*z)/L**2)) +                           &
     (-((ady*densy*Pi*Sin((ady*Pi*y)/L))/L) -                       &
        (adxy*densxy*Pi*x*Sin((adxy*Pi*x*y)/L**2))/L**2 -           &
        (adyz*densyz*Pi*z*Sin((adyz*Pi*y*z)/L**2))/L**2)*           &
      (vvel0 + vvelx*Cos((avx*Pi*x)/L) +                            &
        vvelxy*Cos((avxy*Pi*x*y)/L**2) +                            &
        vvelxz*Cos((avxz*Pi*x*z)/L**2) + vvely*Sin((avy*Pi*y)/L) +  &
        vvelz*Sin((avz*Pi*z)/L) + vvelyz*Sin((avyz*Pi*y*z)/L**2)) + &
     (dens0 + densy*Cos((ady*Pi*y)/L) +                             &
        densxy*Cos((adxy*Pi*x*y)/L**2) +                            &
        densyz*Cos((adyz*Pi*y*z)/L**2) + densx*Sin((adx*Pi*x)/L) +  &
        densz*Sin((adz*Pi*z)/L) + densxz*Sin((adxz*Pi*x*z)/L**2))*  &
      ((awxz*Pi*wvelxz*x*Cos((awxz*Pi*x*z)/L**2))/L**2 -            &
        (awz*Pi*wvelz*Sin((awz*Pi*z)/L))/L -                        &
        (awyz*Pi*wvelyz*y*Sin((awyz*Pi*y*z)/L**2))/L**2)
      !
      hc = metrics_dt(nc)%jac(k)**jac_exp * hc / (rhoref * aref)
      !
      residual(nec,n) = residual(nec,n) + add_or_subtract * hc
!
! #######################################################################
! ########################  X-momentum equation  ########################
! #######################################################################
!
      hc = (apxy*Pi*presxy*y*Cos((apxy*Pi*x*y)/L**2))/L**2 -        &
     (apx*Pi*presx*Sin((apx*Pi*x)/L))/L -                           &
     (apxz*Pi*presxz*z*Sin((apxz*Pi*x*z)/L**2))/L**2 +              &
     ((auxz*Pi*uvelxz*x*Cos((auxz*Pi*x*z)/L**2))/L**2 +             &
        (auyz*Pi*uvelyz*y*Cos((auyz*Pi*y*z)/L**2))/L**2 -           &
        (auz*Pi*uvelz*Sin((auz*Pi*z)/L))/L)*                        &
      (dens0 + densy*Cos((ady*Pi*y)/L) +                            &
        densxy*Cos((adxy*Pi*x*y)/L**2) +                            &
        densyz*Cos((adyz*Pi*y*z)/L**2) + densx*Sin((adx*Pi*x)/L) +  &
        densz*Sin((adz*Pi*z)/L) + densxz*Sin((adxz*Pi*x*z)/L**2))*  &
      (wvel0 + wvelz*Cos((awz*Pi*z)/L) +                            &
        wvelyz*Cos((awyz*Pi*y*z)/L**2) + wvelx*Sin((awx*Pi*x)/L) +  &
        wvely*Sin((awy*Pi*y)/L) + wvelxy*Sin((awxy*Pi*x*y)/L**2) +  &
        wvelxz*Sin((awxz*Pi*x*z)/L**2)) -                           &
     (two*rmu*((avxy**2*Pi**2*vvelxy*x*y*Cos((avxy*Pi*x*y)/L**2))/    &
           L**4 - (awxz*Pi*wvelxz*Cos((awxz*Pi*x*z)/L**2))/L**2 +   &
          (avxy*Pi*vvelxy*Sin((avxy*Pi*x*y)/L**2))/L**2 +           &
          two*(-((auxy**2*Pi**2*uvelxy*y**2*Cos((auxy*Pi*x*y)/L**2))/ &
                L**4) - (aux**2*Pi**2*uvelx*Sin((aux*Pi*x)/L))/     &
              L**2 - (auxz**2*Pi**2*uvelxz*z**2*                    &
                Sin((auxz*Pi*x*z)/L**2))/L**4) +                    &
          (awxz**2*Pi**2*wvelxz*x*z*Sin((awxz*Pi*x*z)/L**2))/L**4))/&
      three + two*((aux*Pi*uvelx*Cos((aux*Pi*x)/L))/L +                  &
        (auxz*Pi*uvelxz*z*Cos((auxz*Pi*x*z)/L**2))/L**2 -           &
        (auxy*Pi*uvelxy*y*Sin((auxy*Pi*x*y)/L**2))/L**2)*           &
      (dens0 + densy*Cos((ady*Pi*y)/L) +                            &
        densxy*Cos((adxy*Pi*x*y)/L**2) +                            &
        densyz*Cos((adyz*Pi*y*z)/L**2) + densx*Sin((adx*Pi*x)/L) +  &
        densz*Sin((adz*Pi*z)/L) + densxz*Sin((adxz*Pi*x*z)/L**2))*  &
      (uvel0 + uvely*Cos((auy*Pi*y)/L) +                            &
        uvelxy*Cos((auxy*Pi*x*y)/L**2) + uvelz*Cos((auz*Pi*z)/L) +  &
        uvelx*Sin((aux*Pi*x)/L) + uvelxz*Sin((auxz*Pi*x*z)/L**2) +  &
        uvelyz*Sin((auyz*Pi*y*z)/L**2)) +                           &
     ((avy*Pi*vvely*Cos((avy*Pi*y)/L))/L +                          &
        (avyz*Pi*vvelyz*z*Cos((avyz*Pi*y*z)/L**2))/L**2 -           &
        (avxy*Pi*vvelxy*x*Sin((avxy*Pi*x*y)/L**2))/L**2)*           &
      (dens0 + densy*Cos((ady*Pi*y)/L) +                            &
        densxy*Cos((adxy*Pi*x*y)/L**2) +                            &
        densyz*Cos((adyz*Pi*y*z)/L**2) + densx*Sin((adx*Pi*x)/L) +  &
        densz*Sin((adz*Pi*z)/L) + densxz*Sin((adxz*Pi*x*z)/L**2))*  &
      (uvel0 + uvely*Cos((auy*Pi*y)/L) +                            &
        uvelxy*Cos((auxy*Pi*x*y)/L**2) + uvelz*Cos((auz*Pi*z)/L) +  &
        uvelx*Sin((aux*Pi*x)/L) + uvelxz*Sin((auxz*Pi*x*z)/L**2) +  &
        uvelyz*Sin((auyz*Pi*y*z)/L**2)) +                           &
     (wvel0 + wvelz*Cos((awz*Pi*z)/L) +                             &
        wvelyz*Cos((awyz*Pi*y*z)/L**2) + wvelx*Sin((awx*Pi*x)/L) +  &
        wvely*Sin((awy*Pi*y)/L) + wvelxy*Sin((awxy*Pi*x*y)/L**2) +  &
        wvelxz*Sin((awxz*Pi*x*z)/L**2))*                            &
      ((adz*densz*Pi*Cos((adz*Pi*z)/L))/L +                         &
        (adxz*densxz*Pi*x*Cos((adxz*Pi*x*z)/L**2))/L**2 -           &
        (adyz*densyz*Pi*y*Sin((adyz*Pi*y*z)/L**2))/L**2)*           &
      (uvel0 + uvely*Cos((auy*Pi*y)/L) +                            &
        uvelxy*Cos((auxy*Pi*x*y)/L**2) + uvelz*Cos((auz*Pi*z)/L) +  &
        uvelx*Sin((aux*Pi*x)/L) + uvelxz*Sin((auxz*Pi*x*z)/L**2) +  &
        uvelyz*Sin((auyz*Pi*y*z)/L**2)) +                           &
     ((adx*densx*Pi*Cos((adx*Pi*x)/L))/L +                          &
        (adxz*densxz*Pi*z*Cos((adxz*Pi*x*z)/L**2))/L**2 -           &
        (adxy*densxy*Pi*y*Sin((adxy*Pi*x*y)/L**2))/L**2)*           &
      (uvel0 + uvely*Cos((auy*Pi*y)/L) +                            &
         uvelxy*Cos((auxy*Pi*x*y)/L**2) + uvelz*Cos((auz*Pi*z)/L) + &
         uvelx*Sin((aux*Pi*x)/L) + uvelxz*Sin((auxz*Pi*x*z)/L**2) + &
         uvelyz*Sin((auyz*Pi*y*z)/L**2))**2 -                       &
     rmu*(-((auz**2*Pi**2*uvelz*Cos((auz*Pi*z)/L))/L**2) +          &
        (awxz*Pi*wvelxz*Cos((awxz*Pi*x*z)/L**2))/L**2 -             &
        (auxz**2*Pi**2*uvelxz*x**2*Sin((auxz*Pi*x*z)/L**2))/L**4 -  &
        (awxz**2*Pi**2*wvelxz*x*z*Sin((awxz*Pi*x*z)/L**2))/L**4 -   &
        (auyz**2*Pi**2*uvelyz*y**2*Sin((auyz*Pi*y*z)/L**2))/L**4) - &
     rmu*(-((auy**2*Pi**2*uvely*Cos((auy*Pi*y)/L))/L**2) -          &
        (auxy**2*Pi**2*uvelxy*x**2*Cos((auxy*Pi*x*y)/L**2))/L**4 -  &
        (avxy**2*Pi**2*vvelxy*x*y*Cos((avxy*Pi*x*y)/L**2))/L**4 -   &
        (avxy*Pi*vvelxy*Sin((avxy*Pi*x*y)/L**2))/L**2 -             &
        (auyz**2*Pi**2*uvelyz*z**2*Sin((auyz*Pi*y*z)/L**2))/L**4) + &
     ((auyz*Pi*uvelyz*z*Cos((auyz*Pi*y*z)/L**2))/L**2 -             &
        (auy*Pi*uvely*Sin((auy*Pi*y)/L))/L -                        &
        (auxy*Pi*uvelxy*x*Sin((auxy*Pi*x*y)/L**2))/L**2)*           &
      (dens0 + densy*Cos((ady*Pi*y)/L) +                            &
        densxy*Cos((adxy*Pi*x*y)/L**2) +                            &
        densyz*Cos((adyz*Pi*y*z)/L**2) + densx*Sin((adx*Pi*x)/L) +  &
        densz*Sin((adz*Pi*z)/L) + densxz*Sin((adxz*Pi*x*z)/L**2))*  &
      (vvel0 + vvelx*Cos((avx*Pi*x)/L) +                            &
        vvelxy*Cos((avxy*Pi*x*y)/L**2) +                            &
        vvelxz*Cos((avxz*Pi*x*z)/L**2) + vvely*Sin((avy*Pi*y)/L) +  &
        vvelz*Sin((avz*Pi*z)/L) + vvelyz*Sin((avyz*Pi*y*z)/L**2)) + &
     (-((ady*densy*Pi*Sin((ady*Pi*y)/L))/L) -                       &
        (adxy*densxy*Pi*x*Sin((adxy*Pi*x*y)/L**2))/L**2 -           &
        (adyz*densyz*Pi*z*Sin((adyz*Pi*y*z)/L**2))/L**2)*           &
      (uvel0 + uvely*Cos((auy*Pi*y)/L) +                            &
        uvelxy*Cos((auxy*Pi*x*y)/L**2) + uvelz*Cos((auz*Pi*z)/L) +  &
        uvelx*Sin((aux*Pi*x)/L) + uvelxz*Sin((auxz*Pi*x*z)/L**2) +  &
        uvelyz*Sin((auyz*Pi*y*z)/L**2))*                            &
      (vvel0 + vvelx*Cos((avx*Pi*x)/L) +                            &
        vvelxy*Cos((avxy*Pi*x*y)/L**2) +                            &
        vvelxz*Cos((avxz*Pi*x*z)/L**2) + vvely*Sin((avy*Pi*y)/L) +  &
        vvelz*Sin((avz*Pi*z)/L) + vvelyz*Sin((avyz*Pi*y*z)/L**2)) + &
     (dens0 + densy*Cos((ady*Pi*y)/L) +                             &
        densxy*Cos((adxy*Pi*x*y)/L**2) +                            &
        densyz*Cos((adyz*Pi*y*z)/L**2) + densx*Sin((adx*Pi*x)/L) +  &
        densz*Sin((adz*Pi*z)/L) + densxz*Sin((adxz*Pi*x*z)/L**2))*  &
      (uvel0 + uvely*Cos((auy*Pi*y)/L) +                            &
        uvelxy*Cos((auxy*Pi*x*y)/L**2) + uvelz*Cos((auz*Pi*z)/L) +  &
        uvelx*Sin((aux*Pi*x)/L) + uvelxz*Sin((auxz*Pi*x*z)/L**2) +  &
        uvelyz*Sin((auyz*Pi*y*z)/L**2))*                            &
      ((awxz*Pi*wvelxz*x*Cos((awxz*Pi*x*z)/L**2))/L**2 -            &
        (awz*Pi*wvelz*Sin((awz*Pi*z)/L))/L -                        &
        (awyz*Pi*wvelyz*y*Sin((awyz*Pi*y*z)/L**2))/L**2)
      !
      hc = metrics_dt(nc)%jac(k)**jac_exp * hc / (rhoref * aref**2)
      !
      residual(nmx,n) = residual(nmx,n) + add_or_subtract * hc
!
! #######################################################################
! ########################  Y-momentum equation  ########################
! #######################################################################
!
      if (nr >= 2) then
        !
        hc = (apy*Pi*presy*Cos((apy*Pi*y)/L))/L +                     &
       (apxy*Pi*presxy*x*Cos((apxy*Pi*x*y)/L**2))/L**2 -              &
       rmu*(-((avx**2*Pi**2*vvelx*Cos((avx*Pi*x)/L))/L**2) -          &
          (auxy**2*Pi**2*uvelxy*x*y*Cos((auxy*Pi*x*y)/L**2))/L**4 -   &
          (avxy**2*Pi**2*vvelxy*y**2*Cos((avxy*Pi*x*y)/L**2))/L**4 -  &
          (avxz**2*Pi**2*vvelxz*z**2*Cos((avxz*Pi*x*z)/L**2))/L**4 -  &
          (auxy*Pi*uvelxy*Sin((auxy*Pi*x*y)/L**2))/L**2) +            &
       (dens0 + densy*Cos((ady*Pi*y)/L) +                             &
          densxy*Cos((adxy*Pi*x*y)/L**2) +                            &
          densyz*Cos((adyz*Pi*y*z)/L**2) + densx*Sin((adx*Pi*x)/L) +  &
          densz*Sin((adz*Pi*z)/L) + densxz*Sin((adxz*Pi*x*z)/L**2))*  &
        ((avz*Pi*vvelz*Cos((avz*Pi*z)/L))/L +                         &
          (avyz*Pi*vvelyz*y*Cos((avyz*Pi*y*z)/L**2))/L**2 -           &
          (avxz*Pi*vvelxz*x*Sin((avxz*Pi*x*z)/L**2))/L**2)*           &
        (wvel0 + wvelz*Cos((awz*Pi*z)/L) +                            &
          wvelyz*Cos((awyz*Pi*y*z)/L**2) + wvelx*Sin((awx*Pi*x)/L) +  &
          wvely*Sin((awy*Pi*y)/L) + wvelxy*Sin((awxy*Pi*x*y)/L**2) +  &
          wvelxz*Sin((awxz*Pi*x*z)/L**2)) -                           &
       (apyz*Pi*presyz*z*Sin((apyz*Pi*y*z)/L**2))/L**2 +              &
       (dens0 + densy*Cos((ady*Pi*y)/L) +                             &
          densxy*Cos((adxy*Pi*x*y)/L**2) +                            &
          densyz*Cos((adyz*Pi*y*z)/L**2) + densx*Sin((adx*Pi*x)/L) +  &
          densz*Sin((adz*Pi*z)/L) + densxz*Sin((adxz*Pi*x*z)/L**2))*  &
        (-((avx*Pi*vvelx*Sin((avx*Pi*x)/L))/L) -                      &
          (avxy*Pi*vvelxy*y*Sin((avxy*Pi*x*y)/L**2))/L**2 -           &
          (avxz*Pi*vvelxz*z*Sin((avxz*Pi*x*z)/L**2))/L**2)*           &
        (uvel0 + uvely*Cos((auy*Pi*y)/L) +                            &
          uvelxy*Cos((auxy*Pi*x*y)/L**2) + uvelz*Cos((auz*Pi*z)/L) +  &
          uvelx*Sin((aux*Pi*x)/L) + uvelxz*Sin((auxz*Pi*x*z)/L**2) +  &
          uvelyz*Sin((auyz*Pi*y*z)/L**2)) +                           &
       ((aux*Pi*uvelx*Cos((aux*Pi*x)/L))/L +                          &
          (auxz*Pi*uvelxz*z*Cos((auxz*Pi*x*z)/L**2))/L**2 -           &
          (auxy*Pi*uvelxy*y*Sin((auxy*Pi*x*y)/L**2))/L**2)*           &
        (dens0 + densy*Cos((ady*Pi*y)/L) +                            &
          densxy*Cos((adxy*Pi*x*y)/L**2) +                            &
          densyz*Cos((adyz*Pi*y*z)/L**2) + densx*Sin((adx*Pi*x)/L) +  &
          densz*Sin((adz*Pi*z)/L) + densxz*Sin((adxz*Pi*x*z)/L**2))*  &
        (vvel0 + vvelx*Cos((avx*Pi*x)/L) +                            &
          vvelxy*Cos((avxy*Pi*x*y)/L**2) +                            &
          vvelxz*Cos((avxz*Pi*x*z)/L**2) + vvely*Sin((avy*Pi*y)/L) +  &
          vvelz*Sin((avz*Pi*z)/L) + vvelyz*Sin((avyz*Pi*y*z)/L**2)) + &
       two*((avy*Pi*vvely*Cos((avy*Pi*y)/L))/L +                        &
          (avyz*Pi*vvelyz*z*Cos((avyz*Pi*y*z)/L**2))/L**2 -           &
          (avxy*Pi*vvelxy*x*Sin((avxy*Pi*x*y)/L**2))/L**2)*           &
        (dens0 + densy*Cos((ady*Pi*y)/L) +                            &
          densxy*Cos((adxy*Pi*x*y)/L**2) +                            &
          densyz*Cos((adyz*Pi*y*z)/L**2) + densx*Sin((adx*Pi*x)/L) +  &
          densz*Sin((adz*Pi*z)/L) + densxz*Sin((adxz*Pi*x*z)/L**2))*  &
        (vvel0 + vvelx*Cos((avx*Pi*x)/L) +                            &
          vvelxy*Cos((avxy*Pi*x*y)/L**2) +                            &
          vvelxz*Cos((avxz*Pi*x*z)/L**2) + vvely*Sin((avy*Pi*y)/L) +  &
          vvelz*Sin((avz*Pi*z)/L) + vvelyz*Sin((avyz*Pi*y*z)/L**2)) + &
       (wvel0 + wvelz*Cos((awz*Pi*z)/L) +                             &
          wvelyz*Cos((awyz*Pi*y*z)/L**2) + wvelx*Sin((awx*Pi*x)/L) +  &
          wvely*Sin((awy*Pi*y)/L) + wvelxy*Sin((awxy*Pi*x*y)/L**2) +  &
          wvelxz*Sin((awxz*Pi*x*z)/L**2))*                            &
        ((adz*densz*Pi*Cos((adz*Pi*z)/L))/L +                         &
          (adxz*densxz*Pi*x*Cos((adxz*Pi*x*z)/L**2))/L**2 -           &
          (adyz*densyz*Pi*y*Sin((adyz*Pi*y*z)/L**2))/L**2)*           &
        (vvel0 + vvelx*Cos((avx*Pi*x)/L) +                            &
          vvelxy*Cos((avxy*Pi*x*y)/L**2) +                            &
          vvelxz*Cos((avxz*Pi*x*z)/L**2) + vvely*Sin((avy*Pi*y)/L) +  &
          vvelz*Sin((avz*Pi*z)/L) + vvelyz*Sin((avyz*Pi*y*z)/L**2)) + &
       ((adx*densx*Pi*Cos((adx*Pi*x)/L))/L +                          &
          (adxz*densxz*Pi*z*Cos((adxz*Pi*x*z)/L**2))/L**2 -           &
          (adxy*densxy*Pi*y*Sin((adxy*Pi*x*y)/L**2))/L**2)*           &
        (uvel0 + uvely*Cos((auy*Pi*y)/L) +                            &
          uvelxy*Cos((auxy*Pi*x*y)/L**2) + uvelz*Cos((auz*Pi*z)/L) +  &
          uvelx*Sin((aux*Pi*x)/L) + uvelxz*Sin((auxz*Pi*x*z)/L**2) +  &
          uvelyz*Sin((auyz*Pi*y*z)/L**2))*                            &
        (vvel0 + vvelx*Cos((avx*Pi*x)/L) +                            &
          vvelxy*Cos((avxy*Pi*x*y)/L**2) +                            &
          vvelxz*Cos((avxz*Pi*x*z)/L**2) + vvely*Sin((avy*Pi*y)/L) +  &
          vvelz*Sin((avz*Pi*z)/L) + vvelyz*Sin((avyz*Pi*y*z)/L**2)) + &
       (-((ady*densy*Pi*Sin((ady*Pi*y)/L))/L) -                       &
          (adxy*densxy*Pi*x*Sin((adxy*Pi*x*y)/L**2))/L**2 -           &
          (adyz*densyz*Pi*z*Sin((adyz*Pi*y*z)/L**2))/L**2)*           &
        (vvel0 + vvelx*Cos((avx*Pi*x)/L) +                            &
           vvelxy*Cos((avxy*Pi*x*y)/L**2) +                           &
           vvelxz*Cos((avxz*Pi*x*z)/L**2) + vvely*Sin((avy*Pi*y)/L) + &
           vvelz*Sin((avz*Pi*z)/L) + vvelyz*Sin((avyz*Pi*y*z)/L**2))**2 &
         - rmu*(-((avxz**2*Pi**2*vvelxz*x**2*                       &
               Cos((avxz*Pi*x*z)/L**2))/L**4) -                       &
          (awyz**2*Pi**2*wvelyz*y*z*Cos((awyz*Pi*y*z)/L**2))/L**4 -   &
          (avz**2*Pi**2*vvelz*Sin((avz*Pi*z)/L))/L**2 -               &
          (avyz**2*Pi**2*vvelyz*y**2*Sin((avyz*Pi*y*z)/L**2))/L**4 -  &
          (awyz*Pi*wvelyz*Sin((awyz*Pi*y*z)/L**2))/L**2) -            &
       (two*rmu*((auxy**2*Pi**2*uvelxy*x*y*Cos((auxy*Pi*x*y)/L**2))/    &
             L**4 + (awyz**2*Pi**2*wvelyz*y*z*                        &
               Cos((awyz*Pi*y*z)/L**2))/L**4 +                        &
            (auxy*Pi*uvelxy*Sin((auxy*Pi*x*y)/L**2))/L**2 +           &
            two*(-((avxy**2*Pi**2*vvelxy*x**2*Cos((avxy*Pi*x*y)/L**2))/ &
                  L**4) - (avy**2*Pi**2*vvely*Sin((avy*Pi*y)/L))/     &
                L**2 - (avyz**2*Pi**2*vvelyz*z**2*                    &
                  Sin((avyz*Pi*y*z)/L**2))/L**4) +                    &
            (awyz*Pi*wvelyz*Sin((awyz*Pi*y*z)/L**2))/L**2))/three +      &
       (dens0 + densy*Cos((ady*Pi*y)/L) +                             &
          densxy*Cos((adxy*Pi*x*y)/L**2) +                            &
          densyz*Cos((adyz*Pi*y*z)/L**2) + densx*Sin((adx*Pi*x)/L) +  &
          densz*Sin((adz*Pi*z)/L) + densxz*Sin((adxz*Pi*x*z)/L**2))*  &
        (vvel0 + vvelx*Cos((avx*Pi*x)/L) +                            &
          vvelxy*Cos((avxy*Pi*x*y)/L**2) +                            &
          vvelxz*Cos((avxz*Pi*x*z)/L**2) + vvely*Sin((avy*Pi*y)/L) +  &
          vvelz*Sin((avz*Pi*z)/L) + vvelyz*Sin((avyz*Pi*y*z)/L**2))*  &
        ((awxz*Pi*wvelxz*x*Cos((awxz*Pi*x*z)/L**2))/L**2 -            &
          (awz*Pi*wvelz*Sin((awz*Pi*z)/L))/L -                        &
          (awyz*Pi*wvelyz*y*Sin((awyz*Pi*y*z)/L**2))/L**2)
        !
        hc = metrics_dt(nc)%jac(k)**jac_exp * hc / (rhoref * aref**2)
        !
        residual(nmy,n) = residual(nmy,n) + add_or_subtract * hc
        !
      end if
!
! #######################################################################
! ########################  Z-momentum equation  ########################
! #######################################################################
!
      if (nr >= 3) then
        !
        hc = -((apz*Pi*presz*Sin((apz*Pi*z)/L))/L) -                  &
       (apxz*Pi*presxz*x*Sin((apxz*Pi*x*z)/L**2))/L**2 +              &
       ((aux*Pi*uvelx*Cos((aux*Pi*x)/L))/L +                          &
          (auxz*Pi*uvelxz*z*Cos((auxz*Pi*x*z)/L**2))/L**2 -           &
          (auxy*Pi*uvelxy*y*Sin((auxy*Pi*x*y)/L**2))/L**2)*           &
        (dens0 + densy*Cos((ady*Pi*y)/L) +                            &
          densxy*Cos((adxy*Pi*x*y)/L**2) +                            &
          densyz*Cos((adyz*Pi*y*z)/L**2) + densx*Sin((adx*Pi*x)/L) +  &
          densz*Sin((adz*Pi*z)/L) + densxz*Sin((adxz*Pi*x*z)/L**2))*  &
        (wvel0 + wvelz*Cos((awz*Pi*z)/L) +                            &
          wvelyz*Cos((awyz*Pi*y*z)/L**2) + wvelx*Sin((awx*Pi*x)/L) +  &
          wvely*Sin((awy*Pi*y)/L) + wvelxy*Sin((awxy*Pi*x*y)/L**2) +  &
          wvelxz*Sin((awxz*Pi*x*z)/L**2)) +                           &
       ((avy*Pi*vvely*Cos((avy*Pi*y)/L))/L +                          &
          (avyz*Pi*vvelyz*z*Cos((avyz*Pi*y*z)/L**2))/L**2 -           &
          (avxy*Pi*vvelxy*x*Sin((avxy*Pi*x*y)/L**2))/L**2)*           &
        (dens0 + densy*Cos((ady*Pi*y)/L) +                            &
          densxy*Cos((adxy*Pi*x*y)/L**2) +                            &
          densyz*Cos((adyz*Pi*y*z)/L**2) + densx*Sin((adx*Pi*x)/L) +  &
          densz*Sin((adz*Pi*z)/L) + densxz*Sin((adxz*Pi*x*z)/L**2))*  &
        (wvel0 + wvelz*Cos((awz*Pi*z)/L) +                            &
          wvelyz*Cos((awyz*Pi*y*z)/L**2) + wvelx*Sin((awx*Pi*x)/L) +  &
          wvely*Sin((awy*Pi*y)/L) + wvelxy*Sin((awxy*Pi*x*y)/L**2) +  &
          wvelxz*Sin((awxz*Pi*x*z)/L**2)) -                           &
       rmu*((auxz*Pi*uvelxz*Cos((auxz*Pi*x*z)/L**2))/L**2 -           &
          (awx**2*Pi**2*wvelx*Sin((awx*Pi*x)/L))/L**2 -               &
          (awxy**2*Pi**2*wvelxy*y**2*Sin((awxy*Pi*x*y)/L**2))/L**4 -  &
          (auxz**2*Pi**2*uvelxz*x*z*Sin((auxz*Pi*x*z)/L**2))/L**4 -   &
          (awxz**2*Pi**2*wvelxz*z**2*Sin((awxz*Pi*x*z)/L**2))/L**4) + &
       (wvel0 + wvelz*Cos((awz*Pi*z)/L) +                             &
           wvelyz*Cos((awyz*Pi*y*z)/L**2) + wvelx*Sin((awx*Pi*x)/L) + &
           wvely*Sin((awy*Pi*y)/L) + wvelxy*Sin((awxy*Pi*x*y)/L**2) + &
           wvelxz*Sin((awxz*Pi*x*z)/L**2))**2*                        &
        ((adz*densz*Pi*Cos((adz*Pi*z)/L))/L +                         &
          (adxz*densxz*Pi*x*Cos((adxz*Pi*x*z)/L**2))/L**2 -           &
          (adyz*densyz*Pi*y*Sin((adyz*Pi*y*z)/L**2))/L**2) -          &
       (apyz*Pi*presyz*y*Sin((apyz*Pi*y*z)/L**2))/L**2 +              &
       ((awx*Pi*wvelx*Cos((awx*Pi*x)/L))/L +                          &
          (awxy*Pi*wvelxy*y*Cos((awxy*Pi*x*y)/L**2))/L**2 +           &
          (awxz*Pi*wvelxz*z*Cos((awxz*Pi*x*z)/L**2))/L**2)*           &
        (dens0 + densy*Cos((ady*Pi*y)/L) +                            &
          densxy*Cos((adxy*Pi*x*y)/L**2) +                            &
          densyz*Cos((adyz*Pi*y*z)/L**2) + densx*Sin((adx*Pi*x)/L) +  &
          densz*Sin((adz*Pi*z)/L) + densxz*Sin((adxz*Pi*x*z)/L**2))*  &
        (uvel0 + uvely*Cos((auy*Pi*y)/L) +                            &
          uvelxy*Cos((auxy*Pi*x*y)/L**2) + uvelz*Cos((auz*Pi*z)/L) +  &
          uvelx*Sin((aux*Pi*x)/L) + uvelxz*Sin((auxz*Pi*x*z)/L**2) +  &
          uvelyz*Sin((auyz*Pi*y*z)/L**2)) +                           &
       ((adx*densx*Pi*Cos((adx*Pi*x)/L))/L +                          &
          (adxz*densxz*Pi*z*Cos((adxz*Pi*x*z)/L**2))/L**2 -           &
          (adxy*densxy*Pi*y*Sin((adxy*Pi*x*y)/L**2))/L**2)*           &
        (wvel0 + wvelz*Cos((awz*Pi*z)/L) +                            &
          wvelyz*Cos((awyz*Pi*y*z)/L**2) + wvelx*Sin((awx*Pi*x)/L) +  &
          wvely*Sin((awy*Pi*y)/L) + wvelxy*Sin((awxy*Pi*x*y)/L**2) +  &
          wvelxz*Sin((awxz*Pi*x*z)/L**2))*                            &
        (uvel0 + uvely*Cos((auy*Pi*y)/L) +                            &
          uvelxy*Cos((auxy*Pi*x*y)/L**2) + uvelz*Cos((auz*Pi*z)/L) +  &
          uvelx*Sin((aux*Pi*x)/L) + uvelxz*Sin((auxz*Pi*x*z)/L**2) +  &
          uvelyz*Sin((auyz*Pi*y*z)/L**2)) +                           &
       (wvel0 + wvelz*Cos((awz*Pi*z)/L) +                             &
          wvelyz*Cos((awyz*Pi*y*z)/L**2) + wvelx*Sin((awx*Pi*x)/L) +  &
          wvely*Sin((awy*Pi*y)/L) + wvelxy*Sin((awxy*Pi*x*y)/L**2) +  &
          wvelxz*Sin((awxz*Pi*x*z)/L**2))*                            &
        (-((ady*densy*Pi*Sin((ady*Pi*y)/L))/L) -                      &
          (adxy*densxy*Pi*x*Sin((adxy*Pi*x*y)/L**2))/L**2 -           &
          (adyz*densyz*Pi*z*Sin((adyz*Pi*y*z)/L**2))/L**2)*           &
        (vvel0 + vvelx*Cos((avx*Pi*x)/L) +                            &
          vvelxy*Cos((avxy*Pi*x*y)/L**2) +                            &
          vvelxz*Cos((avxz*Pi*x*z)/L**2) + vvely*Sin((avy*Pi*y)/L) +  &
          vvelz*Sin((avz*Pi*z)/L) + vvelyz*Sin((avyz*Pi*y*z)/L**2)) - &
       rmu*((avyz*Pi*vvelyz*Cos((avyz*Pi*y*z)/L**2))/L**2 -           &
          (awyz**2*Pi**2*wvelyz*z**2*Cos((awyz*Pi*y*z)/L**2))/L**4 -  &
          (awy**2*Pi**2*wvely*Sin((awy*Pi*y)/L))/L**2 -               &
          (awxy**2*Pi**2*wvelxy*x**2*Sin((awxy*Pi*x*y)/L**2))/L**4 -  &
          (avyz**2*Pi**2*vvelyz*y*z*Sin((avyz*Pi*y*z)/L**2))/L**4) -  &
       (two*rmu*(-((auxz*Pi*uvelxz*Cos((auxz*Pi*x*z)/L**2))/L**2) -     &
            (avyz*Pi*vvelyz*Cos((avyz*Pi*y*z)/L**2))/L**2 +           &
            (auxz**2*Pi**2*uvelxz*x*z*Sin((auxz*Pi*x*z)/L**2))/L**4 + &
            two*(-((awz**2*Pi**2*wvelz*Cos((awz*Pi*z)/L))/L**2) -       &
               (awyz**2*Pi**2*wvelyz*y**2*Cos((awyz*Pi*y*z)/L**2))/   &
                L**4 - (awxz**2*Pi**2*wvelxz*x**2*                    &
                  Sin((awxz*Pi*x*z)/L**2))/L**4) +                    &
            (avyz**2*Pi**2*vvelyz*y*z*Sin((avyz*Pi*y*z)/L**2))/L**4))/&
        three + two*(dens0 + densy*Cos((ady*Pi*y)/L) +                     &
          densxy*Cos((adxy*Pi*x*y)/L**2) +                            &
          densyz*Cos((adyz*Pi*y*z)/L**2) + densx*Sin((adx*Pi*x)/L) +  &
          densz*Sin((adz*Pi*z)/L) + densxz*Sin((adxz*Pi*x*z)/L**2))*  &
        (wvel0 + wvelz*Cos((awz*Pi*z)/L) +                            &
          wvelyz*Cos((awyz*Pi*y*z)/L**2) + wvelx*Sin((awx*Pi*x)/L) +  &
          wvely*Sin((awy*Pi*y)/L) + wvelxy*Sin((awxy*Pi*x*y)/L**2) +  &
          wvelxz*Sin((awxz*Pi*x*z)/L**2))*                            &
        ((awxz*Pi*wvelxz*x*Cos((awxz*Pi*x*z)/L**2))/L**2 -            &
          (awz*Pi*wvelz*Sin((awz*Pi*z)/L))/L -                        &
          (awyz*Pi*wvelyz*y*Sin((awyz*Pi*y*z)/L**2))/L**2) +          &
       (dens0 + densy*Cos((ady*Pi*y)/L) +                             &
          densxy*Cos((adxy*Pi*x*y)/L**2) +                            &
          densyz*Cos((adyz*Pi*y*z)/L**2) + densx*Sin((adx*Pi*x)/L) +  &
          densz*Sin((adz*Pi*z)/L) + densxz*Sin((adxz*Pi*x*z)/L**2))*  &
        (vvel0 + vvelx*Cos((avx*Pi*x)/L) +                            &
          vvelxy*Cos((avxy*Pi*x*y)/L**2) +                            &
          vvelxz*Cos((avxz*Pi*x*z)/L**2) + vvely*Sin((avy*Pi*y)/L) +  &
          vvelz*Sin((avz*Pi*z)/L) + vvelyz*Sin((avyz*Pi*y*z)/L**2))*  &
        ((awy*Pi*wvely*Cos((awy*Pi*y)/L))/L +                         &
          (awxy*Pi*wvelxy*x*Cos((awxy*Pi*x*y)/L**2))/L**2 -           &
          (awyz*Pi*wvelyz*z*Sin((awyz*Pi*y*z)/L**2))/L**2)
        !
        hc = metrics_dt(nc)%jac(k)**jac_exp * hc / (rhoref * aref**2)
        !
        if (nr == 3) then
          residual(nmz,n) = residual(nmz,n) + add_or_subtract * hc
        end if
        !
      end if
!
! #########################################################################
! ########################  Total Energy equation  ########################
! #########################################################################
!
      hc = -(rmu*((awx*Pi*wvelx*Cos((awx*Pi*x)/L))/L +              &
          (awxy*Pi*wvelxy*y*Cos((awxy*Pi*x*y)/L**2))/L**2 +         &
          (awxz*Pi*wvelxz*z*Cos((awxz*Pi*x*z)/L**2))/L**2)*         &
        ((awx*Pi*wvelx*Cos((awx*Pi*x)/L))/L +                       &
          (awxy*Pi*wvelxy*y*Cos((awxy*Pi*x*y)/L**2))/L**2 +         &
          (auxz*Pi*uvelxz*x*Cos((auxz*Pi*x*z)/L**2))/L**2 +         &
          (awxz*Pi*wvelxz*z*Cos((awxz*Pi*x*z)/L**2))/L**2 +         &
          (auyz*Pi*uvelyz*y*Cos((auyz*Pi*y*z)/L**2))/L**2 -         &
          (auz*Pi*uvelz*Sin((auz*Pi*z)/L))/L)) -                    &
     rmu*((auxz*Pi*uvelxz*x*Cos((auxz*Pi*x*z)/L**2))/L**2 +         &
        (auyz*Pi*uvelyz*y*Cos((auyz*Pi*y*z)/L**2))/L**2 -           &
        (auz*Pi*uvelz*Sin((auz*Pi*z)/L))/L)*                        &
      ((awx*Pi*wvelx*Cos((awx*Pi*x)/L))/L +                         &
        (awxy*Pi*wvelxy*y*Cos((awxy*Pi*x*y)/L**2))/L**2 +           &
        (auxz*Pi*uvelxz*x*Cos((auxz*Pi*x*z)/L**2))/L**2 +           &
        (awxz*Pi*wvelxz*z*Cos((awxz*Pi*x*z)/L**2))/L**2 +           &
        (auyz*Pi*uvelyz*y*Cos((auyz*Pi*y*z)/L**2))/L**2 -           &
        (auz*Pi*uvelz*Sin((auz*Pi*z)/L))/L) -                       &
     (gam*R*rmu*((two*((adx*densx*Pi*Cos((adx*Pi*x)/L))/L +         &
                (adxz*densxz*Pi*z*Cos((adxz*Pi*x*z)/L**2))/L**2 -   &
                (adxy*densxy*Pi*y*Sin((adxy*Pi*x*y)/L**2))/L**2)**2*&
             (pres0 + presx*Cos((apx*Pi*x)/L) +                     &
               presz*Cos((apz*Pi*z)/L) +                            &
               presxz*Cos((apxz*Pi*x*z)/L**2) +                     &
               presyz*Cos((apyz*Pi*y*z)/L**2) +                     &
               presy*Sin((apy*Pi*y)/L) +                            &
               presxy*Sin((apxy*Pi*x*y)/L**2)))/                    &
           (R*(dens0 + densy*Cos((ady*Pi*y)/L) +                    &
                densxy*Cos((adxy*Pi*x*y)/L**2) +                    &
                densyz*Cos((adyz*Pi*y*z)/L**2) +                    &
                densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) + &
                densxz*Sin((adxz*Pi*x*z)/L**2))**3) +               &
          (-((apx**2*Pi**2*presx*Cos((apx*Pi*x)/L))/L**2) -         &
             (apxz**2*Pi**2*presxz*z**2*Cos((apxz*Pi*x*z)/L**2))/   &
              L**4 - (apxy**2*Pi**2*presxy*y**2*                    &
                Sin((apxy*Pi*x*y)/L**2))/L**4)/                     &
           (R*(dens0 + densy*Cos((ady*Pi*y)/L) +                    &
               densxy*Cos((adxy*Pi*x*y)/L**2) +                     &
               densyz*Cos((adyz*Pi*y*z)/L**2) +                     &
               densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) +  &
               densxz*Sin((adxz*Pi*x*z)/L**2))) -                   &
          ((pres0 + presx*Cos((apx*Pi*x)/L) +                       &
               presz*Cos((apz*Pi*z)/L) +                            &
               presxz*Cos((apxz*Pi*x*z)/L**2) +                     &
               presyz*Cos((apyz*Pi*y*z)/L**2) +                     &
               presy*Sin((apy*Pi*y)/L) +                            &
               presxy*Sin((apxy*Pi*x*y)/L**2))*                     &
             (-((adxy**2*densxy*Pi**2*y**2*Cos((adxy*Pi*x*y)/L**2))/&
                  L**4) -                                           &
               (adx**2*densx*Pi**2*Sin((adx*Pi*x)/L))/L**2 -        &
               (adxz**2*densxz*Pi**2*z**2*Sin((adxz*Pi*x*z)/L**2))/ &
                L**4))/                                             &
           (R*(dens0 + densy*Cos((ady*Pi*y)/L) +                    &
                densxy*Cos((adxy*Pi*x*y)/L**2) +                    &
                densyz*Cos((adyz*Pi*y*z)/L**2) +                    &
                densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) + &
                densxz*Sin((adxz*Pi*x*z)/L**2))**2) -               &
          (two*((adx*densx*Pi*Cos((adx*Pi*x)/L))/L +                &
               (adxz*densxz*Pi*z*Cos((adxz*Pi*x*z)/L**2))/L**2 -    &
               (adxy*densxy*Pi*y*Sin((adxy*Pi*x*y)/L**2))/L**2)*    &
             ((apxy*Pi*presxy*y*Cos((apxy*Pi*x*y)/L**2))/L**2 -     &
               (apx*Pi*presx*Sin((apx*Pi*x)/L))/L -                 &
               (apxz*Pi*presxz*z*Sin((apxz*Pi*x*z)/L**2))/L**2))/   &
           (R*(dens0 + densy*Cos((ady*Pi*y)/L) +                    &
                densxy*Cos((adxy*Pi*x*y)/L**2) +                    &
                densyz*Cos((adyz*Pi*y*z)/L**2) +                    &
                densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) + &
                densxz*Sin((adxz*Pi*x*z)/L**2))**2))) / (gm1*Pr)
      !--------------------------------------------------------------
      hc = hc -                                                     &
     rmu*((auyz*Pi*uvelyz*z*Cos((auyz*Pi*y*z)/L**2))/L**2 -         &
        (auy*Pi*uvely*Sin((auy*Pi*y)/L))/L -                        &
        (auxy*Pi*uvelxy*x*Sin((auxy*Pi*x*y)/L**2))/L**2)*           &
      ((auyz*Pi*uvelyz*z*Cos((auyz*Pi*y*z)/L**2))/L**2 -            &
        (avx*Pi*vvelx*Sin((avx*Pi*x)/L))/L -                        &
        (auy*Pi*uvely*Sin((auy*Pi*y)/L))/L -                        &
        (auxy*Pi*uvelxy*x*Sin((auxy*Pi*x*y)/L**2))/L**2 -           &
        (avxy*Pi*vvelxy*y*Sin((avxy*Pi*x*y)/L**2))/L**2 -           &
        (avxz*Pi*vvelxz*z*Sin((avxz*Pi*x*z)/L**2))/L**2) -          &
     rmu*(-((avx*Pi*vvelx*Sin((avx*Pi*x)/L))/L) -                   &
        (avxy*Pi*vvelxy*y*Sin((avxy*Pi*x*y)/L**2))/L**2 -           &
        (avxz*Pi*vvelxz*z*Sin((avxz*Pi*x*z)/L**2))/L**2)*           &
      ((auyz*Pi*uvelyz*z*Cos((auyz*Pi*y*z)/L**2))/L**2 -            &
        (avx*Pi*vvelx*Sin((avx*Pi*x)/L))/L -                        &
        (auy*Pi*uvely*Sin((auy*Pi*y)/L))/L -                        &
        (auxy*Pi*uvelxy*x*Sin((auxy*Pi*x*y)/L**2))/L**2 -           &
        (avxy*Pi*vvelxy*y*Sin((avxy*Pi*x*y)/L**2))/L**2 -           &
        (avxz*Pi*vvelxz*z*Sin((avxz*Pi*x*z)/L**2))/L**2) -          &
     rmu*(wvel0 + wvelz*Cos((awz*Pi*z)/L) +                         &
        wvelyz*Cos((awyz*Pi*y*z)/L**2) + wvelx*Sin((awx*Pi*x)/L) +  &
        wvely*Sin((awy*Pi*y)/L) + wvelxy*Sin((awxy*Pi*x*y)/L**2) +  &
        wvelxz*Sin((awxz*Pi*x*z)/L**2))*                            &
      ((auxz*Pi*uvelxz*Cos((auxz*Pi*x*z)/L**2))/L**2 -              &
        (awx**2*Pi**2*wvelx*Sin((awx*Pi*x)/L))/L**2 -               &
        (awxy**2*Pi**2*wvelxy*y**2*Sin((awxy*Pi*x*y)/L**2))/L**4 -  &
        (auxz**2*Pi**2*uvelxz*x*z*Sin((auxz*Pi*x*z)/L**2))/L**4 -   &
        (awxz**2*Pi**2*wvelxz*z**2*Sin((awxz*Pi*x*z)/L**2))/L**4) - &
     (gam*R*rmu*((-((apz**2*Pi**2*presz*Cos((apz*Pi*z)/L))/       &
                L**2) - (apxz**2*Pi**2*presxz*x**2*                 &
                Cos((apxz*Pi*x*z)/L**2))/L**4 -                     &
             (apyz**2*Pi**2*presyz*y**2*Cos((apyz*Pi*y*z)/L**2))/   &
              L**4)/                                                &
           (R*(dens0 + densy*Cos((ady*Pi*y)/L) +                    &
               densxy*Cos((adxy*Pi*x*y)/L**2) +                     &
               densyz*Cos((adyz*Pi*y*z)/L**2) +                     &
               densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) +  &
               densxz*Sin((adxz*Pi*x*z)/L**2))) -                   &
          ((pres0 + presx*Cos((apx*Pi*x)/L) +                       &
               presz*Cos((apz*Pi*z)/L) +                            &
               presxz*Cos((apxz*Pi*x*z)/L**2) +                     &
               presyz*Cos((apyz*Pi*y*z)/L**2) +                     &
               presy*Sin((apy*Pi*y)/L) +                            &
               presxy*Sin((apxy*Pi*x*y)/L**2))*                     &
             (-((adyz**2*densyz*Pi**2*y**2*Cos((adyz*Pi*y*z)/L**2))/&
                  L**4) -                                           &
               (adz**2*densz*Pi**2*Sin((adz*Pi*z)/L))/L**2 -        &
               (adxz**2*densxz*Pi**2*x**2*Sin((adxz*Pi*x*z)/L**2))/ &
                L**4))/                                             &
           (R*(dens0 + densy*Cos((ady*Pi*y)/L) +                    &
                densxy*Cos((adxy*Pi*x*y)/L**2) +                    &
                densyz*Cos((adyz*Pi*y*z)/L**2) +                    &
                densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) + &
                densxz*Sin((adxz*Pi*x*z)/L**2))**2) +               &
          (two*(pres0 + presx*Cos((apx*Pi*x)/L) +                     &
               presz*Cos((apz*Pi*z)/L) +                            &
               presxz*Cos((apxz*Pi*x*z)/L**2) +                     &
               presyz*Cos((apyz*Pi*y*z)/L**2) +                     &
               presy*Sin((apy*Pi*y)/L) +                            &
               presxy*Sin((apxy*Pi*x*y)/L**2))*                     &
             ((adz*densz*Pi*Cos((adz*Pi*z)/L))/L +                  &
                (adxz*densxz*Pi*x*Cos((adxz*Pi*x*z)/L**2))/L**2 -   &
                (adyz*densyz*Pi*y*Sin((adyz*Pi*y*z)/L**2))/L**2)**2)&
            /(R*(dens0 + densy*Cos((ady*Pi*y)/L) +                  &
                densxy*Cos((adxy*Pi*x*y)/L**2) +                    &
                densyz*Cos((adyz*Pi*y*z)/L**2) +                    &
                densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) + &
                densxz*Sin((adxz*Pi*x*z)/L**2))**3) -               &
          (two*((adz*densz*Pi*Cos((adz*Pi*z)/L))/L +                  &
               (adxz*densxz*Pi*x*Cos((adxz*Pi*x*z)/L**2))/L**2 -    &
               (adyz*densyz*Pi*y*Sin((adyz*Pi*y*z)/L**2))/L**2)*    &
             (-((apz*Pi*presz*Sin((apz*Pi*z)/L))/L) -               &
               (apxz*Pi*presxz*x*Sin((apxz*Pi*x*z)/L**2))/L**2 -    &
               (apyz*Pi*presyz*y*Sin((apyz*Pi*y*z)/L**2))/L**2))/   &
           (R*(dens0 + densy*Cos((ady*Pi*y)/L) +                    &
                densxy*Cos((adxy*Pi*x*y)/L**2) +                    &
                densyz*Cos((adyz*Pi*y*z)/L**2) +                    &
                densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) + &
                densxz*Sin((adxz*Pi*x*z)/L**2))**2))) / (gm1*Pr)
      !--------------------------------------------------------------
      hc = hc -                                                     &
     (gam*R*rmu*(-(((-((ady**2*densy*Pi**2*Cos((ady*Pi*y)/L))/    &
                    L**2) -                                         &
                 (adxy**2*densxy*Pi**2*x**2*                        &
                    Cos((adxy*Pi*x*y)/L**2))/L**4 -                 &
                 (adyz**2*densyz*Pi**2*z**2*                        &
                    Cos((adyz*Pi*y*z)/L**2))/L**4)*                 &
               (pres0 + presx*Cos((apx*Pi*x)/L) +                   &
                 presz*Cos((apz*Pi*z)/L) +                          &
                 presxz*Cos((apxz*Pi*x*z)/L**2) +                   &
                 presyz*Cos((apyz*Pi*y*z)/L**2) +                   &
                 presy*Sin((apy*Pi*y)/L) +                          &
                 presxy*Sin((apxy*Pi*x*y)/L**2)))/                  &
             (R*(dens0 + densy*Cos((ady*Pi*y)/L) +                  &
                  densxy*Cos((adxy*Pi*x*y)/L**2) +                  &
                  densyz*Cos((adyz*Pi*y*z)/L**2) +                  &
                  densx*Sin((adx*Pi*x)/L) +                         &
                  densz*Sin((adz*Pi*z)/L) +                         &
                  densxz*Sin((adxz*Pi*x*z)/L**2))**2)) +            &
          (-((apyz**2*Pi**2*presyz*z**2*Cos((apyz*Pi*y*z)/L**2))/   &
                L**4) - (apy**2*Pi**2*presy*Sin((apy*Pi*y)/L))/     &
              L**2 - (apxy**2*Pi**2*presxy*x**2*                    &
                Sin((apxy*Pi*x*y)/L**2))/L**4)/                     &
           (R*(dens0 + densy*Cos((ady*Pi*y)/L) +                    &
               densxy*Cos((adxy*Pi*x*y)/L**2) +                     &
               densyz*Cos((adyz*Pi*y*z)/L**2) +                     &
               densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) +  &
               densxz*Sin((adxz*Pi*x*z)/L**2))) +                   &
          (two*(pres0 + presx*Cos((apx*Pi*x)/L) +                     &
               presz*Cos((apz*Pi*z)/L) +                            &
               presxz*Cos((apxz*Pi*x*z)/L**2) +                     &
               presyz*Cos((apyz*Pi*y*z)/L**2) +                     &
               presy*Sin((apy*Pi*y)/L) +                            &
               presxy*Sin((apxy*Pi*x*y)/L**2))*                     &
             (-((ady*densy*Pi*Sin((ady*Pi*y)/L))/L) -               &
                (adxy*densxy*Pi*x*Sin((adxy*Pi*x*y)/L**2))/L**2 -   &
                (adyz*densyz*Pi*z*Sin((adyz*Pi*y*z)/L**2))/L**2)**2)&
            /(R*(dens0 + densy*Cos((ady*Pi*y)/L) +                  &
                densxy*Cos((adxy*Pi*x*y)/L**2) +                    &
                densyz*Cos((adyz*Pi*y*z)/L**2) +                    &
                densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) + &
                densxz*Sin((adxz*Pi*x*z)/L**2))**3) -               &
          (two*(-((ady*densy*Pi*Sin((ady*Pi*y)/L))/L) -               &
               (adxy*densxy*Pi*x*Sin((adxy*Pi*x*y)/L**2))/L**2 -    &
               (adyz*densyz*Pi*z*Sin((adyz*Pi*y*z)/L**2))/L**2)*    &
             ((apy*Pi*presy*Cos((apy*Pi*y)/L))/L +                  &
               (apxy*Pi*presxy*x*Cos((apxy*Pi*x*y)/L**2))/L**2 -    &
               (apyz*Pi*presyz*z*Sin((apyz*Pi*y*z)/L**2))/L**2))/   &
           (R*(dens0 + densy*Cos((ady*Pi*y)/L) +                    &
                densxy*Cos((adxy*Pi*x*y)/L**2) +                    &
                densyz*Cos((adyz*Pi*y*z)/L**2) +                    &
                densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) + &
                densxz*Sin((adxz*Pi*x*z)/L**2))**2))) / (gm1*Pr)
      !--------------------------------------------------------------
      hc = hc -                                                     &
     (two*rmu*((avxy**2*Pi**2*vvelxy*x*y*Cos((avxy*Pi*x*y)/L**2))/    &
           L**4 - (awxz*Pi*wvelxz*Cos((awxz*Pi*x*z)/L**2))/L**2 +   &
          (avxy*Pi*vvelxy*Sin((avxy*Pi*x*y)/L**2))/L**2 +           &
          two*(-((auxy**2*Pi**2*uvelxy*y**2*Cos((auxy*Pi*x*y)/L**2))/ &
                L**4) - (aux**2*Pi**2*uvelx*Sin((aux*Pi*x)/L))/     &
              L**2 - (auxz**2*Pi**2*uvelxz*z**2*                    &
                Sin((auxz*Pi*x*z)/L**2))/L**4) +                    &
          (awxz**2*Pi**2*wvelxz*x*z*Sin((awxz*Pi*x*z)/L**2))/L**4)* &
        (uvel0 + uvely*Cos((auy*Pi*y)/L) +                          &
          uvelxy*Cos((auxy*Pi*x*y)/L**2) +                          &
          uvelz*Cos((auz*Pi*z)/L) + uvelx*Sin((aux*Pi*x)/L) +       &
          uvelxz*Sin((auxz*Pi*x*z)/L**2) +                          &
          uvelyz*Sin((auyz*Pi*y*z)/L**2)))/three -                     &
     rmu*(uvel0 + uvely*Cos((auy*Pi*y)/L) +                         &
        uvelxy*Cos((auxy*Pi*x*y)/L**2) + uvelz*Cos((auz*Pi*z)/L) +  &
        uvelx*Sin((aux*Pi*x)/L) + uvelxz*Sin((auxz*Pi*x*z)/L**2) +  &
        uvelyz*Sin((auyz*Pi*y*z)/L**2))*                            &
      (-((auz**2*Pi**2*uvelz*Cos((auz*Pi*z)/L))/L**2) +             &
        (awxz*Pi*wvelxz*Cos((awxz*Pi*x*z)/L**2))/L**2 -             &
        (auxz**2*Pi**2*uvelxz*x**2*Sin((auxz*Pi*x*z)/L**2))/L**4 -  &
        (awxz**2*Pi**2*wvelxz*x*z*Sin((awxz*Pi*x*z)/L**2))/L**4 -   &
        (auyz**2*Pi**2*uvelyz*y**2*Sin((auyz*Pi*y*z)/L**2))/L**4) - &
     rmu*(uvel0 + uvely*Cos((auy*Pi*y)/L) +                         &
        uvelxy*Cos((auxy*Pi*x*y)/L**2) + uvelz*Cos((auz*Pi*z)/L) +  &
        uvelx*Sin((aux*Pi*x)/L) + uvelxz*Sin((auxz*Pi*x*z)/L**2) +  &
        uvelyz*Sin((auyz*Pi*y*z)/L**2))*                            &
      (-((auy**2*Pi**2*uvely*Cos((auy*Pi*y)/L))/L**2) -             &
        (auxy**2*Pi**2*uvelxy*x**2*Cos((auxy*Pi*x*y)/L**2))/L**4 -  &
        (avxy**2*Pi**2*vvelxy*x*y*Cos((avxy*Pi*x*y)/L**2))/L**4 -   &
        (avxy*Pi*vvelxy*Sin((avxy*Pi*x*y)/L**2))/L**2 -             &
        (auyz**2*Pi**2*uvelyz*z**2*Sin((auyz*Pi*y*z)/L**2))/L**4) - &
     rmu*(-((avx**2*Pi**2*vvelx*Cos((avx*Pi*x)/L))/L**2) -          &
        (auxy**2*Pi**2*uvelxy*x*y*Cos((auxy*Pi*x*y)/L**2))/L**4 -   &
        (avxy**2*Pi**2*vvelxy*y**2*Cos((avxy*Pi*x*y)/L**2))/L**4 -  &
        (avxz**2*Pi**2*vvelxz*z**2*Cos((avxz*Pi*x*z)/L**2))/L**4 -  &
        (auxy*Pi*uvelxy*Sin((auxy*Pi*x*y)/L**2))/L**2)*             &
      (vvel0 + vvelx*Cos((avx*Pi*x)/L) +                            &
        vvelxy*Cos((avxy*Pi*x*y)/L**2) +                            &
        vvelxz*Cos((avxz*Pi*x*z)/L**2) + vvely*Sin((avy*Pi*y)/L) +  &
        vvelz*Sin((avz*Pi*z)/L) + vvelyz*Sin((avyz*Pi*y*z)/L**2)) - &
     rmu*(wvel0 + wvelz*Cos((awz*Pi*z)/L) +                         &
        wvelyz*Cos((awyz*Pi*y*z)/L**2) + wvelx*Sin((awx*Pi*x)/L) +  &
        wvely*Sin((awy*Pi*y)/L) + wvelxy*Sin((awxy*Pi*x*y)/L**2) +  &
        wvelxz*Sin((awxz*Pi*x*z)/L**2))*                            &
      ((avyz*Pi*vvelyz*Cos((avyz*Pi*y*z)/L**2))/L**2 -              &
        (awyz**2*Pi**2*wvelyz*z**2*Cos((awyz*Pi*y*z)/L**2))/L**4 -  &
        (awy**2*Pi**2*wvely*Sin((awy*Pi*y)/L))/L**2 -               &
        (awxy**2*Pi**2*wvelxy*x**2*Sin((awxy*Pi*x*y)/L**2))/L**4 -  &
        (avyz**2*Pi**2*vvelyz*y*z*Sin((avyz*Pi*y*z)/L**2))/L**4) -  &
     (two*rmu*(wvel0 + wvelz*Cos((awz*Pi*z)/L) +                      &
          wvelyz*Cos((awyz*Pi*y*z)/L**2) +                          &
          wvelx*Sin((awx*Pi*x)/L) + wvely*Sin((awy*Pi*y)/L) +       &
          wvelxy*Sin((awxy*Pi*x*y)/L**2) +                          &
          wvelxz*Sin((awxz*Pi*x*z)/L**2))*                          &
        (-((auxz*Pi*uvelxz*Cos((auxz*Pi*x*z)/L**2))/L**2) -         &
          (avyz*Pi*vvelyz*Cos((avyz*Pi*y*z)/L**2))/L**2 +           &
          (auxz**2*Pi**2*uvelxz*x*z*Sin((auxz*Pi*x*z)/L**2))/L**4 + &
          two*(-((awz**2*Pi**2*wvelz*Cos((awz*Pi*z)/L))/L**2) -       &
             (awyz**2*Pi**2*wvelyz*y**2*Cos((awyz*Pi*y*z)/L**2))/   &
              L**4 - (awxz**2*Pi**2*wvelxz*x**2*                    &
                Sin((awxz*Pi*x*z)/L**2))/L**4) +                    &
          (avyz**2*Pi**2*vvelyz*y*z*Sin((avyz*Pi*y*z)/L**2))/L**4))/&
      three + (uvel0 + uvely*Cos((auy*Pi*y)/L) +                       &
        uvelxy*Cos((auxy*Pi*x*y)/L**2) + uvelz*Cos((auz*Pi*z)/L) +  &
        uvelx*Sin((aux*Pi*x)/L) + uvelxz*Sin((auxz*Pi*x*z)/L**2) +  &
        uvelyz*Sin((auyz*Pi*y*z)/L**2))*                            &
      ((apxy*Pi*presxy*y*Cos((apxy*Pi*x*y)/L**2))/L**2 -            &
        (apx*Pi*presx*Sin((apx*Pi*x)/L))/L -                        &
        (apxz*Pi*presxz*z*Sin((apxz*Pi*x*z)/L**2))/L**2 +           &
        (dens0 + densy*Cos((ady*Pi*y)/L) +                          &
           densxy*Cos((adxy*Pi*x*y)/L**2) +                         &
           densyz*Cos((adyz*Pi*y*z)/L**2) +                         &
           densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) +      &
           densxz*Sin((adxz*Pi*x*z)/L**2))*                         &
         (-((((adx*densx*Pi*Cos((adx*Pi*x)/L))/L +                  &
                  (adxz*densxz*Pi*z*Cos((adxz*Pi*x*z)/L**2))/L**2 - &
                  (adxy*densxy*Pi*y*Sin((adxy*Pi*x*y)/L**2))/L**2)* &
                (pres0 + presx*Cos((apx*Pi*x)/L) +                  &
                  presz*Cos((apz*Pi*z)/L) +                         &
                  presxz*Cos((apxz*Pi*x*z)/L**2) +                  &
                  presyz*Cos((apyz*Pi*y*z)/L**2) +                  &
                  presy*Sin((apy*Pi*y)/L) +                         &
                  presxy*Sin((apxy*Pi*x*y)/L**2)))/                 &
              (gm1*                                        &
                (dens0 + densy*Cos((ady*Pi*y)/L) +                  &
                   densxy*Cos((adxy*Pi*x*y)/L**2) +                 &
                   densyz*Cos((adyz*Pi*y*z)/L**2) +                 &
                   densx*Sin((adx*Pi*x)/L) +                        &
                   densz*Sin((adz*Pi*z)/L) +                        &
                   densxz*Sin((adxz*Pi*x*z)/L**2))**2)) +           &
           ((apxy*Pi*presxy*y*Cos((apxy*Pi*x*y)/L**2))/L**2 -       &
              (apx*Pi*presx*Sin((apx*Pi*x)/L))/L -                  &
              (apxz*Pi*presxz*z*Sin((apxz*Pi*x*z)/L**2))/L**2)/     &
            (gm1*                                          &
              (dens0 + densy*Cos((ady*Pi*y)/L) +                    &
                densxy*Cos((adxy*Pi*x*y)/L**2) +                    &
                densyz*Cos((adyz*Pi*y*z)/L**2) +                    &
                densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) + &
                densxz*Sin((adxz*Pi*x*z)/L**2))) +                  &
           (two*((awx*Pi*wvelx*Cos((awx*Pi*x)/L))/L +                 &
                 (awxy*Pi*wvelxy*y*Cos((awxy*Pi*x*y)/L**2))/L**2 +  &
                 (awxz*Pi*wvelxz*z*Cos((awxz*Pi*x*z)/L**2))/L**2)*  &
               (wvel0 + wvelz*Cos((awz*Pi*z)/L) +                   &
                 wvelyz*Cos((awyz*Pi*y*z)/L**2) +                   &
                 wvelx*Sin((awx*Pi*x)/L) +                          &
                 wvely*Sin((awy*Pi*y)/L) +                          &
                 wvelxy*Sin((awxy*Pi*x*y)/L**2) +                   &
                 wvelxz*Sin((awxz*Pi*x*z)/L**2)) +                  &
              two*((aux*Pi*uvelx*Cos((aux*Pi*x)/L))/L +               &
                 (auxz*Pi*uvelxz*z*Cos((auxz*Pi*x*z)/L**2))/L**2 -  &
                 (auxy*Pi*uvelxy*y*Sin((auxy*Pi*x*y)/L**2))/L**2)*  &
               (uvel0 + uvely*Cos((auy*Pi*y)/L) +                   &
                 uvelxy*Cos((auxy*Pi*x*y)/L**2) +                   &
                 uvelz*Cos((auz*Pi*z)/L) +                          &
                 uvelx*Sin((aux*Pi*x)/L) +                          &
                 uvelxz*Sin((auxz*Pi*x*z)/L**2) +                   &
                 uvelyz*Sin((auyz*Pi*y*z)/L**2)) +                  &
              two*(-((avx*Pi*vvelx*Sin((avx*Pi*x)/L))/L) -            &
                 (avxy*Pi*vvelxy*y*Sin((avxy*Pi*x*y)/L**2))/L**2 -  &
                 (avxz*Pi*vvelxz*z*Sin((avxz*Pi*x*z)/L**2))/L**2)*  &
               (vvel0 + vvelx*Cos((avx*Pi*x)/L) +                   &
                 vvelxy*Cos((avxy*Pi*x*y)/L**2) +                   &
                 vvelxz*Cos((avxz*Pi*x*z)/L**2) +                   &
                 vvely*Sin((avy*Pi*y)/L) +                          &
                 vvelz*Sin((avz*Pi*z)/L) +                          &
                 vvelyz*Sin((avyz*Pi*y*z)/L**2)))/two) +             &
        ((adx*densx*Pi*Cos((adx*Pi*x)/L))/L +                       &
           (adxz*densxz*Pi*z*Cos((adxz*Pi*x*z)/L**2))/L**2 -        &
           (adxy*densxy*Pi*y*Sin((adxy*Pi*x*y)/L**2))/L**2)*        &
         ((pres0 + presx*Cos((apx*Pi*x)/L) +                        &
              presz*Cos((apz*Pi*z)/L) +                             &
              presxz*Cos((apxz*Pi*x*z)/L**2) +                      &
              presyz*Cos((apyz*Pi*y*z)/L**2) +                      &
              presy*Sin((apy*Pi*y)/L) +                             &
              presxy*Sin((apxy*Pi*x*y)/L**2))/                      &
            (gm1*                                          &
              (dens0 + densy*Cos((ady*Pi*y)/L) +                    &
                densxy*Cos((adxy*Pi*x*y)/L**2) +                    &
                densyz*Cos((adyz*Pi*y*z)/L**2) +                    &
                densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) + &
                densxz*Sin((adxz*Pi*x*z)/L**2))) +                  &
           ((wvel0 + wvelz*Cos((awz*Pi*z)/L) +                      &
                 wvelyz*Cos((awyz*Pi*y*z)/L**2) +                   &
                 wvelx*Sin((awx*Pi*x)/L) +                          &
                 wvely*Sin((awy*Pi*y)/L) +                          &
                 wvelxy*Sin((awxy*Pi*x*y)/L**2) +                   &
                 wvelxz*Sin((awxz*Pi*x*z)/L**2))**2 +               &
              (uvel0 + uvely*Cos((auy*Pi*y)/L) +                    &
                 uvelxy*Cos((auxy*Pi*x*y)/L**2) +                   &
                 uvelz*Cos((auz*Pi*z)/L) +                          &
                 uvelx*Sin((aux*Pi*x)/L) +                          &
                 uvelxz*Sin((auxz*Pi*x*z)/L**2) +                   &
                 uvelyz*Sin((auyz*Pi*y*z)/L**2))**2 +               &
              (vvel0 + vvelx*Cos((avx*Pi*x)/L) +                    &
                 vvelxy*Cos((avxy*Pi*x*y)/L**2) +                   &
                 vvelxz*Cos((avxz*Pi*x*z)/L**2) +                   &
                 vvely*Sin((avy*Pi*y)/L) +                          &
                 vvelz*Sin((avz*Pi*z)/L) +                          &
                 vvelyz*Sin((avyz*Pi*y*z)/L**2))**2)/two))
      !--------------------------------------------------------------
      hc = hc +                                                     &
     ((aux*Pi*uvelx*Cos((aux*Pi*x)/L))/L +                          &
        (auxz*Pi*uvelxz*z*Cos((auxz*Pi*x*z)/L**2))/L**2 -           &
        (auxy*Pi*uvelxy*y*Sin((auxy*Pi*x*y)/L**2))/L**2)*           &
      (pres0 + presx*Cos((apx*Pi*x)/L) + presz*Cos((apz*Pi*z)/L) +  &
        presxz*Cos((apxz*Pi*x*z)/L**2) +                            &
        presyz*Cos((apyz*Pi*y*z)/L**2) + presy*Sin((apy*Pi*y)/L) +  &
        presxy*Sin((apxy*Pi*x*y)/L**2) +                            &
        (dens0 + densy*Cos((ady*Pi*y)/L) +                          &
           densxy*Cos((adxy*Pi*x*y)/L**2) +                         &
           densyz*Cos((adyz*Pi*y*z)/L**2) +                         &
           densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) +      &
           densxz*Sin((adxz*Pi*x*z)/L**2))*                         &
         ((pres0 + presx*Cos((apx*Pi*x)/L) +                        &
              presz*Cos((apz*Pi*z)/L) +                             &
              presxz*Cos((apxz*Pi*x*z)/L**2) +                      &
              presyz*Cos((apyz*Pi*y*z)/L**2) +                      &
              presy*Sin((apy*Pi*y)/L) +                             &
              presxy*Sin((apxy*Pi*x*y)/L**2))/                      &
            (gm1*                                          &
              (dens0 + densy*Cos((ady*Pi*y)/L) +                    &
                densxy*Cos((adxy*Pi*x*y)/L**2) +                    &
                densyz*Cos((adyz*Pi*y*z)/L**2) +                    &
                densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) + &
                densxz*Sin((adxz*Pi*x*z)/L**2))) +                  &
           ((wvel0 + wvelz*Cos((awz*Pi*z)/L) +                      &
                 wvelyz*Cos((awyz*Pi*y*z)/L**2) +                   &
                 wvelx*Sin((awx*Pi*x)/L) +                          &
                 wvely*Sin((awy*Pi*y)/L) +                          &
                 wvelxy*Sin((awxy*Pi*x*y)/L**2) +                   &
                 wvelxz*Sin((awxz*Pi*x*z)/L**2))**2 +               &
              (uvel0 + uvely*Cos((auy*Pi*y)/L) +                    &
                 uvelxy*Cos((auxy*Pi*x*y)/L**2) +                   &
                 uvelz*Cos((auz*Pi*z)/L) +                          &
                 uvelx*Sin((aux*Pi*x)/L) +                          &
                 uvelxz*Sin((auxz*Pi*x*z)/L**2) +                   &
                 uvelyz*Sin((auyz*Pi*y*z)/L**2))**2 +               &
              (vvel0 + vvelx*Cos((avx*Pi*x)/L) +                    &
                 vvelxy*Cos((avxy*Pi*x*y)/L**2) +                   &
                 vvelxz*Cos((avxz*Pi*x*z)/L**2) +                   &
                 vvely*Sin((avy*Pi*y)/L) +                          &
                 vvelz*Sin((avz*Pi*z)/L) +                          &
                 vvelyz*Sin((avyz*Pi*y*z)/L**2))**2)/two))
      !--------------------------------------------------------------
      hc = hc +                                                     &
     ((avy*Pi*vvely*Cos((avy*Pi*y)/L))/L +                          &
        (avyz*Pi*vvelyz*z*Cos((avyz*Pi*y*z)/L**2))/L**2 -           &
        (avxy*Pi*vvelxy*x*Sin((avxy*Pi*x*y)/L**2))/L**2)*           &
      (pres0 + presx*Cos((apx*Pi*x)/L) + presz*Cos((apz*Pi*z)/L) +  &
        presxz*Cos((apxz*Pi*x*z)/L**2) +                            &
        presyz*Cos((apyz*Pi*y*z)/L**2) + presy*Sin((apy*Pi*y)/L) +  &
        presxy*Sin((apxy*Pi*x*y)/L**2) +                            &
        (dens0 + densy*Cos((ady*Pi*y)/L) +                          &
           densxy*Cos((adxy*Pi*x*y)/L**2) +                         &
           densyz*Cos((adyz*Pi*y*z)/L**2) +                         &
           densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) +      &
           densxz*Sin((adxz*Pi*x*z)/L**2))*                         &
         ((pres0 + presx*Cos((apx*Pi*x)/L) +                        &
              presz*Cos((apz*Pi*z)/L) +                             &
              presxz*Cos((apxz*Pi*x*z)/L**2) +                      &
              presyz*Cos((apyz*Pi*y*z)/L**2) +                      &
              presy*Sin((apy*Pi*y)/L) +                             &
              presxy*Sin((apxy*Pi*x*y)/L**2))/                      &
            (gm1*                                          &
              (dens0 + densy*Cos((ady*Pi*y)/L) +                    &
                densxy*Cos((adxy*Pi*x*y)/L**2) +                    &
                densyz*Cos((adyz*Pi*y*z)/L**2) +                    &
                densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) + &
                densxz*Sin((adxz*Pi*x*z)/L**2))) +                  &
           ((wvel0 + wvelz*Cos((awz*Pi*z)/L) +                      &
                 wvelyz*Cos((awyz*Pi*y*z)/L**2) +                   &
                 wvelx*Sin((awx*Pi*x)/L) +                          &
                 wvely*Sin((awy*Pi*y)/L) +                          &
                 wvelxy*Sin((awxy*Pi*x*y)/L**2) +                   &
                 wvelxz*Sin((awxz*Pi*x*z)/L**2))**2 +               &
              (uvel0 + uvely*Cos((auy*Pi*y)/L) +                    &
                 uvelxy*Cos((auxy*Pi*x*y)/L**2) +                   &
                 uvelz*Cos((auz*Pi*z)/L) +                          &
                 uvelx*Sin((aux*Pi*x)/L) +                          &
                 uvelxz*Sin((auxz*Pi*x*z)/L**2) +                   &
                 uvelyz*Sin((auyz*Pi*y*z)/L**2))**2 +               &
              (vvel0 + vvelx*Cos((avx*Pi*x)/L) +                    &
                 vvelxy*Cos((avxy*Pi*x*y)/L**2) +                   &
                 vvelxz*Cos((avxz*Pi*x*z)/L**2) +                   &
                 vvely*Sin((avy*Pi*y)/L) +                          &
                 vvelz*Sin((avz*Pi*z)/L) +                          &
                 vvelyz*Sin((avyz*Pi*y*z)/L**2))**2)/two))
      !--------------------------------------------------------------
      hc = hc -                                                     &
     rmu*(vvel0 + vvelx*Cos((avx*Pi*x)/L) +                         &
        vvelxy*Cos((avxy*Pi*x*y)/L**2) +                            &
        vvelxz*Cos((avxz*Pi*x*z)/L**2) + vvely*Sin((avy*Pi*y)/L) +  &
        vvelz*Sin((avz*Pi*z)/L) + vvelyz*Sin((avyz*Pi*y*z)/L**2))*  &
      (-((avxz**2*Pi**2*vvelxz*x**2*Cos((avxz*Pi*x*z)/L**2))/       &
           L**4) - (awyz**2*Pi**2*wvelyz*y*z*                       &
           Cos((awyz*Pi*y*z)/L**2))/L**4 -                          &
        (avz**2*Pi**2*vvelz*Sin((avz*Pi*z)/L))/L**2 -               &
        (avyz**2*Pi**2*vvelyz*y**2*Sin((avyz*Pi*y*z)/L**2))/L**4 -  &
        (awyz*Pi*wvelyz*Sin((awyz*Pi*y*z)/L**2))/L**2) -            &
     (two*rmu*(vvel0 + vvelx*Cos((avx*Pi*x)/L) +                      &
          vvelxy*Cos((avxy*Pi*x*y)/L**2) +                          &
          vvelxz*Cos((avxz*Pi*x*z)/L**2) +                          &
          vvely*Sin((avy*Pi*y)/L) + vvelz*Sin((avz*Pi*z)/L) +       &
          vvelyz*Sin((avyz*Pi*y*z)/L**2))*                          &
        ((auxy**2*Pi**2*uvelxy*x*y*Cos((auxy*Pi*x*y)/L**2))/L**4 +  &
          (awyz**2*Pi**2*wvelyz*y*z*Cos((awyz*Pi*y*z)/L**2))/L**4 + &
          (auxy*Pi*uvelxy*Sin((auxy*Pi*x*y)/L**2))/L**2 +           &
          two*(-((avxy**2*Pi**2*vvelxy*x**2*Cos((avxy*Pi*x*y)/L**2))/ &
                L**4) - (avy**2*Pi**2*vvely*Sin((avy*Pi*y)/L))/     &
              L**2 - (avyz**2*Pi**2*vvelyz*z**2*                    &
                Sin((avyz*Pi*y*z)/L**2))/L**4) +                    &
          (awyz*Pi*wvelyz*Sin((awyz*Pi*y*z)/L**2))/L**2))/three +      &
     (pres0 + presx*Cos((apx*Pi*x)/L) + presz*Cos((apz*Pi*z)/L) +   &
        presxz*Cos((apxz*Pi*x*z)/L**2) +                            &
        presyz*Cos((apyz*Pi*y*z)/L**2) + presy*Sin((apy*Pi*y)/L) +  &
        presxy*Sin((apxy*Pi*x*y)/L**2) +                            &
        (dens0 + densy*Cos((ady*Pi*y)/L) +                          &
           densxy*Cos((adxy*Pi*x*y)/L**2) +                         &
           densyz*Cos((adyz*Pi*y*z)/L**2) +                         &
           densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) +      &
           densxz*Sin((adxz*Pi*x*z)/L**2))*                         &
         ((pres0 + presx*Cos((apx*Pi*x)/L) +                        &
              presz*Cos((apz*Pi*z)/L) +                             &
              presxz*Cos((apxz*Pi*x*z)/L**2) +                      &
              presyz*Cos((apyz*Pi*y*z)/L**2) +                      &
              presy*Sin((apy*Pi*y)/L) +                             &
              presxy*Sin((apxy*Pi*x*y)/L**2))/                      &
            (gm1*                                          &
              (dens0 + densy*Cos((ady*Pi*y)/L) +                    &
                densxy*Cos((adxy*Pi*x*y)/L**2) +                    &
                densyz*Cos((adyz*Pi*y*z)/L**2) +                    &
                densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) + &
                densxz*Sin((adxz*Pi*x*z)/L**2))) +                  &
           ((wvel0 + wvelz*Cos((awz*Pi*z)/L) +                      &
                 wvelyz*Cos((awyz*Pi*y*z)/L**2) +                   &
                 wvelx*Sin((awx*Pi*x)/L) +                          &
                 wvely*Sin((awy*Pi*y)/L) +                          &
                 wvelxy*Sin((awxy*Pi*x*y)/L**2) +                   &
                 wvelxz*Sin((awxz*Pi*x*z)/L**2))**2 +               &
              (uvel0 + uvely*Cos((auy*Pi*y)/L) +                    &
                 uvelxy*Cos((auxy*Pi*x*y)/L**2) +                   &
                 uvelz*Cos((auz*Pi*z)/L) +                          &
                 uvelx*Sin((aux*Pi*x)/L) +                          &
                 uvelxz*Sin((auxz*Pi*x*z)/L**2) +                   &
                 uvelyz*Sin((auyz*Pi*y*z)/L**2))**2 +               &
              (vvel0 + vvelx*Cos((avx*Pi*x)/L) +                    &
                 vvelxy*Cos((avxy*Pi*x*y)/L**2) +                   &
                 vvelxz*Cos((avxz*Pi*x*z)/L**2) +                   &
                 vvely*Sin((avy*Pi*y)/L) +                          &
                 vvelz*Sin((avz*Pi*z)/L) +                          &
                 vvelyz*Sin((avyz*Pi*y*z)/L**2))**2)/two))*          &
      ((awxz*Pi*wvelxz*x*Cos((awxz*Pi*x*z)/L**2))/L**2 -            &
        (awz*Pi*wvelz*Sin((awz*Pi*z)/L))/L -                        &
        (awyz*Pi*wvelyz*y*Sin((awyz*Pi*y*z)/L**2))/L**2) -          &
     (two*rmu*((aux*Pi*uvelx*Cos((aux*Pi*x)/L))/L +                   &
          (auxz*Pi*uvelxz*z*Cos((auxz*Pi*x*z)/L**2))/L**2 -         &
          (auxy*Pi*uvelxy*y*Sin((auxy*Pi*x*y)/L**2))/L**2)*         &
        (-((avy*Pi*vvely*Cos((avy*Pi*y)/L))/L) -                    &
          (awxz*Pi*wvelxz*x*Cos((awxz*Pi*x*z)/L**2))/L**2 -         &
          (avyz*Pi*vvelyz*z*Cos((avyz*Pi*y*z)/L**2))/L**2 +         &
          two*((aux*Pi*uvelx*Cos((aux*Pi*x)/L))/L +                   &
             (auxz*Pi*uvelxz*z*Cos((auxz*Pi*x*z)/L**2))/L**2 -      &
             (auxy*Pi*uvelxy*y*Sin((auxy*Pi*x*y)/L**2))/L**2) +     &
          (avxy*Pi*vvelxy*x*Sin((avxy*Pi*x*y)/L**2))/L**2 +         &
          (awz*Pi*wvelz*Sin((awz*Pi*z)/L))/L +                      &
          (awyz*Pi*wvelyz*y*Sin((awyz*Pi*y*z)/L**2))/L**2))/three
      !--------------------------------------------------------------
      hc = hc -                                                     &
     (two*rmu*((avy*Pi*vvely*Cos((avy*Pi*y)/L))/L +                   &
          (avyz*Pi*vvelyz*z*Cos((avyz*Pi*y*z)/L**2))/L**2 -         &
          (avxy*Pi*vvelxy*x*Sin((avxy*Pi*x*y)/L**2))/L**2)*         &
        (-((aux*Pi*uvelx*Cos((aux*Pi*x)/L))/L) -                    &
          (auxz*Pi*uvelxz*z*Cos((auxz*Pi*x*z)/L**2))/L**2 -         &
          (awxz*Pi*wvelxz*x*Cos((awxz*Pi*x*z)/L**2))/L**2 +         &
          (auxy*Pi*uvelxy*y*Sin((auxy*Pi*x*y)/L**2))/L**2 +         &
          two*((avy*Pi*vvely*Cos((avy*Pi*y)/L))/L +                   &
             (avyz*Pi*vvelyz*z*Cos((avyz*Pi*y*z)/L**2))/L**2 -      &
             (avxy*Pi*vvelxy*x*Sin((avxy*Pi*x*y)/L**2))/L**2) +     &
          (awz*Pi*wvelz*Sin((awz*Pi*z)/L))/L +                      &
          (awyz*Pi*wvelyz*y*Sin((awyz*Pi*y*z)/L**2))/L**2))/three -    &
     rmu*((avz*Pi*vvelz*Cos((avz*Pi*z)/L))/L +                      &
        (avyz*Pi*vvelyz*y*Cos((avyz*Pi*y*z)/L**2))/L**2 -           &
        (avxz*Pi*vvelxz*x*Sin((avxz*Pi*x*z)/L**2))/L**2)*           &
      ((awy*Pi*wvely*Cos((awy*Pi*y)/L))/L +                         &
        (awxy*Pi*wvelxy*x*Cos((awxy*Pi*x*y)/L**2))/L**2 +           &
        (avz*Pi*vvelz*Cos((avz*Pi*z)/L))/L +                        &
        (avyz*Pi*vvelyz*y*Cos((avyz*Pi*y*z)/L**2))/L**2 -           &
        (avxz*Pi*vvelxz*x*Sin((avxz*Pi*x*z)/L**2))/L**2 -           &
        (awyz*Pi*wvelyz*z*Sin((awyz*Pi*y*z)/L**2))/L**2) -          &
     rmu*((awy*Pi*wvely*Cos((awy*Pi*y)/L))/L +                      &
        (awxy*Pi*wvelxy*x*Cos((awxy*Pi*x*y)/L**2))/L**2 -           &
        (awyz*Pi*wvelyz*z*Sin((awyz*Pi*y*z)/L**2))/L**2)*           &
      ((awy*Pi*wvely*Cos((awy*Pi*y)/L))/L +                         &
        (awxy*Pi*wvelxy*x*Cos((awxy*Pi*x*y)/L**2))/L**2 +           &
        (avz*Pi*vvelz*Cos((avz*Pi*z)/L))/L +                        &
        (avyz*Pi*vvelyz*y*Cos((avyz*Pi*y*z)/L**2))/L**2 -           &
        (avxz*Pi*vvelxz*x*Sin((avxz*Pi*x*z)/L**2))/L**2 -           &
        (awyz*Pi*wvelyz*z*Sin((awyz*Pi*y*z)/L**2))/L**2) -          &
     (two*rmu*((awxz*Pi*wvelxz*x*Cos((awxz*Pi*x*z)/L**2))/L**2 -      &
          (awz*Pi*wvelz*Sin((awz*Pi*z)/L))/L -                      &
          (awyz*Pi*wvelyz*y*Sin((awyz*Pi*y*z)/L**2))/L**2)*         &
        (-((aux*Pi*uvelx*Cos((aux*Pi*x)/L))/L) -                    &
          (avy*Pi*vvely*Cos((avy*Pi*y)/L))/L -                      &
          (auxz*Pi*uvelxz*z*Cos((auxz*Pi*x*z)/L**2))/L**2 -         &
          (avyz*Pi*vvelyz*z*Cos((avyz*Pi*y*z)/L**2))/L**2 +         &
          (auxy*Pi*uvelxy*y*Sin((auxy*Pi*x*y)/L**2))/L**2 +         &
          (avxy*Pi*vvelxy*x*Sin((avxy*Pi*x*y)/L**2))/L**2 +         &
          two*((awxz*Pi*wvelxz*x*Cos((awxz*Pi*x*z)/L**2))/L**2 -      &
             (awz*Pi*wvelz*Sin((awz*Pi*z)/L))/L -                   &
             (awyz*Pi*wvelyz*y*Sin((awyz*Pi*y*z)/L**2))/L**2)))/three  &
      + (wvel0 + wvelz*Cos((awz*Pi*z)/L) +                          &
        wvelyz*Cos((awyz*Pi*y*z)/L**2) + wvelx*Sin((awx*Pi*x)/L) +  &
        wvely*Sin((awy*Pi*y)/L) + wvelxy*Sin((awxy*Pi*x*y)/L**2) +  &
        wvelxz*Sin((awxz*Pi*x*z)/L**2))*                            &
      (-((apz*Pi*presz*Sin((apz*Pi*z)/L))/L) -                      &
        (apxz*Pi*presxz*x*Sin((apxz*Pi*x*z)/L**2))/L**2 -           &
        (apyz*Pi*presyz*y*Sin((apyz*Pi*y*z)/L**2))/L**2 +           &
        ((adz*densz*Pi*Cos((adz*Pi*z)/L))/L +                       &
           (adxz*densxz*Pi*x*Cos((adxz*Pi*x*z)/L**2))/L**2 -        &
           (adyz*densyz*Pi*y*Sin((adyz*Pi*y*z)/L**2))/L**2)*        &
         ((pres0 + presx*Cos((apx*Pi*x)/L) +                        &
              presz*Cos((apz*Pi*z)/L) +                             &
              presxz*Cos((apxz*Pi*x*z)/L**2) +                      &
              presyz*Cos((apyz*Pi*y*z)/L**2) +                      &
              presy*Sin((apy*Pi*y)/L) +                             &
              presxy*Sin((apxy*Pi*x*y)/L**2))/                      &
            (gm1*                                          &
              (dens0 + densy*Cos((ady*Pi*y)/L) +                    &
                densxy*Cos((adxy*Pi*x*y)/L**2) +                    &
                densyz*Cos((adyz*Pi*y*z)/L**2) +                    &
                densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) + &
                densxz*Sin((adxz*Pi*x*z)/L**2))) +                  &
           ((wvel0 + wvelz*Cos((awz*Pi*z)/L) +                      &
                 wvelyz*Cos((awyz*Pi*y*z)/L**2) +                   &
                 wvelx*Sin((awx*Pi*x)/L) +                          &
                 wvely*Sin((awy*Pi*y)/L) +                          &
                 wvelxy*Sin((awxy*Pi*x*y)/L**2) +                   &
                 wvelxz*Sin((awxz*Pi*x*z)/L**2))**2 +               &
              (uvel0 + uvely*Cos((auy*Pi*y)/L) +                    &
                 uvelxy*Cos((auxy*Pi*x*y)/L**2) +                   &
                 uvelz*Cos((auz*Pi*z)/L) +                          &
                 uvelx*Sin((aux*Pi*x)/L) +                          &
                 uvelxz*Sin((auxz*Pi*x*z)/L**2) +                   &
                 uvelyz*Sin((auyz*Pi*y*z)/L**2))**2 +               &
              (vvel0 + vvelx*Cos((avx*Pi*x)/L) +                    &
                 vvelxy*Cos((avxy*Pi*x*y)/L**2) +                   &
                 vvelxz*Cos((avxz*Pi*x*z)/L**2) +                   &
                 vvely*Sin((avy*Pi*y)/L) +                          &
                 vvelz*Sin((avz*Pi*z)/L) +                          &
                 vvelyz*Sin((avyz*Pi*y*z)/L**2))**2)/two) +          &
        (dens0 + densy*Cos((ady*Pi*y)/L) +                          &
           densxy*Cos((adxy*Pi*x*y)/L**2) +                         &
           densyz*Cos((adyz*Pi*y*z)/L**2) +                         &
           densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) +      &
           densxz*Sin((adxz*Pi*x*z)/L**2))*                         &
         (-(((pres0 + presx*Cos((apx*Pi*x)/L) +                     &
                  presz*Cos((apz*Pi*z)/L) +                         &
                  presxz*Cos((apxz*Pi*x*z)/L**2) +                  &
                  presyz*Cos((apyz*Pi*y*z)/L**2) +                  &
                  presy*Sin((apy*Pi*y)/L) +                         &
                  presxy*Sin((apxy*Pi*x*y)/L**2))*                  &
                ((adz*densz*Pi*Cos((adz*Pi*z)/L))/L +               &
                  (adxz*densxz*Pi*x*Cos((adxz*Pi*x*z)/L**2))/L**2 - &
                  (adyz*densyz*Pi*y*Sin((adyz*Pi*y*z)/L**2))/L**2))/&
              (gm1*                                        &
                (dens0 + densy*Cos((ady*Pi*y)/L) +                  &
                   densxy*Cos((adxy*Pi*x*y)/L**2) +                 &
                   densyz*Cos((adyz*Pi*y*z)/L**2) +                 &
                   densx*Sin((adx*Pi*x)/L) +                        &
                   densz*Sin((adz*Pi*z)/L) +                        &
                   densxz*Sin((adxz*Pi*x*z)/L**2))**2)) +           &
           (-((apz*Pi*presz*Sin((apz*Pi*z)/L))/L) -                 &
              (apxz*Pi*presxz*x*Sin((apxz*Pi*x*z)/L**2))/L**2 -     &
              (apyz*Pi*presyz*y*Sin((apyz*Pi*y*z)/L**2))/L**2)/     &
            (gm1*                                          &
              (dens0 + densy*Cos((ady*Pi*y)/L) +                    &
                densxy*Cos((adxy*Pi*x*y)/L**2) +                    &
                densyz*Cos((adyz*Pi*y*z)/L**2) +                    &
                densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) + &
                densxz*Sin((adxz*Pi*x*z)/L**2))) +                  &
           (two*((auxz*Pi*uvelxz*x*Cos((auxz*Pi*x*z)/L**2))/L**2 +    &
                 (auyz*Pi*uvelyz*y*Cos((auyz*Pi*y*z)/L**2))/L**2 -  &
                 (auz*Pi*uvelz*Sin((auz*Pi*z)/L))/L)*               &
               (uvel0 + uvely*Cos((auy*Pi*y)/L) +                   &
                 uvelxy*Cos((auxy*Pi*x*y)/L**2) +                   &
                 uvelz*Cos((auz*Pi*z)/L) +                          &
                 uvelx*Sin((aux*Pi*x)/L) +                          &
                 uvelxz*Sin((auxz*Pi*x*z)/L**2) +                   &
                 uvelyz*Sin((auyz*Pi*y*z)/L**2)) +                  &
              two*((avz*Pi*vvelz*Cos((avz*Pi*z)/L))/L +               &
                 (avyz*Pi*vvelyz*y*Cos((avyz*Pi*y*z)/L**2))/L**2 -  &
                 (avxz*Pi*vvelxz*x*Sin((avxz*Pi*x*z)/L**2))/L**2)*  &
               (vvel0 + vvelx*Cos((avx*Pi*x)/L) +                   &
                 vvelxy*Cos((avxy*Pi*x*y)/L**2) +                   &
                 vvelxz*Cos((avxz*Pi*x*z)/L**2) +                   &
                 vvely*Sin((avy*Pi*y)/L) +                          &
                 vvelz*Sin((avz*Pi*z)/L) +                          &
                 vvelyz*Sin((avyz*Pi*y*z)/L**2)) +                  &
              two*(wvel0 + wvelz*Cos((awz*Pi*z)/L) +                  &
                 wvelyz*Cos((awyz*Pi*y*z)/L**2) +                   &
                 wvelx*Sin((awx*Pi*x)/L) +                          &
                 wvely*Sin((awy*Pi*y)/L) +                          &
                 wvelxy*Sin((awxy*Pi*x*y)/L**2) +                   &
                 wvelxz*Sin((awxz*Pi*x*z)/L**2))*                   &
               ((awxz*Pi*wvelxz*x*Cos((awxz*Pi*x*z)/L**2))/L**2 -   &
                 (awz*Pi*wvelz*Sin((awz*Pi*z)/L))/L -               &
                 (awyz*Pi*wvelyz*y*Sin((awyz*Pi*y*z)/L**2))/L**2))/ &
            two)) + (vvel0 + vvelx*Cos((avx*Pi*x)/L) +               &
        vvelxy*Cos((avxy*Pi*x*y)/L**2) +                            &
        vvelxz*Cos((avxz*Pi*x*z)/L**2) + vvely*Sin((avy*Pi*y)/L) +  &
        vvelz*Sin((avz*Pi*z)/L) + vvelyz*Sin((avyz*Pi*y*z)/L**2))*  &
      ((apy*Pi*presy*Cos((apy*Pi*y)/L))/L +                         &
        (apxy*Pi*presxy*x*Cos((apxy*Pi*x*y)/L**2))/L**2 -           &
        (apyz*Pi*presyz*z*Sin((apyz*Pi*y*z)/L**2))/L**2 +           &
        (-((ady*densy*Pi*Sin((ady*Pi*y)/L))/L) -                    &
           (adxy*densxy*Pi*x*Sin((adxy*Pi*x*y)/L**2))/L**2 -        &
           (adyz*densyz*Pi*z*Sin((adyz*Pi*y*z)/L**2))/L**2)*        &
         ((pres0 + presx*Cos((apx*Pi*x)/L) +                        &
              presz*Cos((apz*Pi*z)/L) +                             &
              presxz*Cos((apxz*Pi*x*z)/L**2) +                      &
              presyz*Cos((apyz*Pi*y*z)/L**2) +                      &
              presy*Sin((apy*Pi*y)/L) +                             &
              presxy*Sin((apxy*Pi*x*y)/L**2))/                      &
            (gm1*                                          &
              (dens0 + densy*Cos((ady*Pi*y)/L) +                    &
                densxy*Cos((adxy*Pi*x*y)/L**2) +                    &
                densyz*Cos((adyz*Pi*y*z)/L**2) +                    &
                densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) + &
                densxz*Sin((adxz*Pi*x*z)/L**2))) +                  &
           ((wvel0 + wvelz*Cos((awz*Pi*z)/L) +                      &
                 wvelyz*Cos((awyz*Pi*y*z)/L**2) +                   &
                 wvelx*Sin((awx*Pi*x)/L) +                          &
                 wvely*Sin((awy*Pi*y)/L) +                          &
                 wvelxy*Sin((awxy*Pi*x*y)/L**2) +                   &
                 wvelxz*Sin((awxz*Pi*x*z)/L**2))**2 +               &
              (uvel0 + uvely*Cos((auy*Pi*y)/L) +                    &
                 uvelxy*Cos((auxy*Pi*x*y)/L**2) +                   &
                 uvelz*Cos((auz*Pi*z)/L) +                          &
                 uvelx*Sin((aux*Pi*x)/L) +                          &
                 uvelxz*Sin((auxz*Pi*x*z)/L**2) +                   &
                 uvelyz*Sin((auyz*Pi*y*z)/L**2))**2 +               &
              (vvel0 + vvelx*Cos((avx*Pi*x)/L) +                    &
                 vvelxy*Cos((avxy*Pi*x*y)/L**2) +                   &
                 vvelxz*Cos((avxz*Pi*x*z)/L**2) +                   &
                 vvely*Sin((avy*Pi*y)/L) +                          &
                 vvelz*Sin((avz*Pi*z)/L) +                          &
                 vvelyz*Sin((avyz*Pi*y*z)/L**2))**2)/two) +          &
        (dens0 + densy*Cos((ady*Pi*y)/L) +                          &
           densxy*Cos((adxy*Pi*x*y)/L**2) +                         &
           densyz*Cos((adyz*Pi*y*z)/L**2) +                         &
           densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) +      &
           densxz*Sin((adxz*Pi*x*z)/L**2))*                         &
         (-(((pres0 + presx*Cos((apx*Pi*x)/L) +                     &
                  presz*Cos((apz*Pi*z)/L) +                         &
                  presxz*Cos((apxz*Pi*x*z)/L**2) +                  &
                  presyz*Cos((apyz*Pi*y*z)/L**2) +                  &
                  presy*Sin((apy*Pi*y)/L) +                         &
                  presxy*Sin((apxy*Pi*x*y)/L**2))*                  &
                (-((ady*densy*Pi*Sin((ady*Pi*y)/L))/L) -            &
                  (adxy*densxy*Pi*x*Sin((adxy*Pi*x*y)/L**2))/L**2 - &
                  (adyz*densyz*Pi*z*Sin((adyz*Pi*y*z)/L**2))/L**2))/&
              (gm1*                                        &
                (dens0 + densy*Cos((ady*Pi*y)/L) +                  &
                   densxy*Cos((adxy*Pi*x*y)/L**2) +                 &
                   densyz*Cos((adyz*Pi*y*z)/L**2) +                 &
                   densx*Sin((adx*Pi*x)/L) +                        &
                   densz*Sin((adz*Pi*z)/L) +                        &
                   densxz*Sin((adxz*Pi*x*z)/L**2))**2)) +           &
           ((apy*Pi*presy*Cos((apy*Pi*y)/L))/L +                    &
              (apxy*Pi*presxy*x*Cos((apxy*Pi*x*y)/L**2))/L**2 -     &
              (apyz*Pi*presyz*z*Sin((apyz*Pi*y*z)/L**2))/L**2)/     &
            (gm1*                                          &
              (dens0 + densy*Cos((ady*Pi*y)/L) +                    &
                densxy*Cos((adxy*Pi*x*y)/L**2) +                    &
                densyz*Cos((adyz*Pi*y*z)/L**2) +                    &
                densx*Sin((adx*Pi*x)/L) + densz*Sin((adz*Pi*z)/L) + &
                densxz*Sin((adxz*Pi*x*z)/L**2))) +                  &
           (two*((auyz*Pi*uvelyz*z*Cos((auyz*Pi*y*z)/L**2))/L**2 -    &
                 (auy*Pi*uvely*Sin((auy*Pi*y)/L))/L -               &
                 (auxy*Pi*uvelxy*x*Sin((auxy*Pi*x*y)/L**2))/L**2)*  &
               (uvel0 + uvely*Cos((auy*Pi*y)/L) +                   &
                 uvelxy*Cos((auxy*Pi*x*y)/L**2) +                   &
                 uvelz*Cos((auz*Pi*z)/L) +                          &
                 uvelx*Sin((aux*Pi*x)/L) +                          &
                 uvelxz*Sin((auxz*Pi*x*z)/L**2) +                   &
                 uvelyz*Sin((auyz*Pi*y*z)/L**2)) +                  &
              two*((avy*Pi*vvely*Cos((avy*Pi*y)/L))/L +               &
                 (avyz*Pi*vvelyz*z*Cos((avyz*Pi*y*z)/L**2))/L**2 -  &
                 (avxy*Pi*vvelxy*x*Sin((avxy*Pi*x*y)/L**2))/L**2)*  &
               (vvel0 + vvelx*Cos((avx*Pi*x)/L) +                   &
                 vvelxy*Cos((avxy*Pi*x*y)/L**2) +                   &
                 vvelxz*Cos((avxz*Pi*x*z)/L**2) +                   &
                 vvely*Sin((avy*Pi*y)/L) +                          &
                 vvelz*Sin((avz*Pi*z)/L) +                          &
                 vvelyz*Sin((avyz*Pi*y*z)/L**2)) +                  &
              two*(wvel0 + wvelz*Cos((awz*Pi*z)/L) +                  &
                 wvelyz*Cos((awyz*Pi*y*z)/L**2) +                   &
                 wvelx*Sin((awx*Pi*x)/L) +                          &
                 wvely*Sin((awy*Pi*y)/L) +                          &
                 wvelxy*Sin((awxy*Pi*x*y)/L**2) +                   &
                 wvelxz*Sin((awxz*Pi*x*z)/L**2))*                   &
               ((awy*Pi*wvely*Cos((awy*Pi*y)/L))/L +                &
                 (awxy*Pi*wvelxy*x*Cos((awxy*Pi*x*y)/L**2))/L**2 -  &
                 (awyz*Pi*wvelyz*z*Sin((awyz*Pi*y*z)/L**2))/L**2))/ &
            two))
      !
      hc = metrics_dt(nc)%jac(k)**jac_exp * hc / (rhoref * aref**3)
      !
      residual(nee,n) = residual(nee,n) + add_or_subtract * hc
      !
    end do
    !
  end do
  !
  ! Check for realizeability in the mms solution
  !
  if (var_min(nec) < small*rhoref) then
    write (iout,*)
    write (iout,*) "**************** ERROR ***************** "
    write (iout,*) "MMS SOLUTION DENSITY VALUE WAS TOO SMALL "
    write (iout,*) "ERROR OCCURRED ON PROCESSOR NO. ", mypnum
    write (iout,*) "**************** ERROR ***************** "
    write (iout,*)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  end if
  !
  if (var_min(nq+1) < tmin) then
    write (iout,*)
    write (iout,*) "****************** ERROR ******************* "
    write (iout,*) "MMS SOLUTION TEMPERATURE VALUE WAS TOO SMALL "
    write (iout,*) "ERROR OCCURRED ON PROCESSOR NO. ", mypnum
    write (iout,*) "****************** ERROR ******************* "
    write (iout,*)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  end if
  !
  if (var_max(nq+1) > tmax) then
    write (iout,*)
    write (iout,*) "****************** ERROR ******************* "
    write (iout,*) "MMS SOLUTION TEMPERATURE VALUE WAS TOO LARGE "
    write (iout,*) "ERROR OCCURRED ON PROCESSOR NO. ", mypnum
    write (iout,*) "****************** ERROR ******************* "
    write (iout,*)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Finished applying source terms for the method of manufactured solutions
  !
end subroutine mms_src
!
!###############################################################################
!
include "Functions/vsp2u.f90"
!
!###############################################################################
!
end module mms_mod
!
!###############################################################################
!
