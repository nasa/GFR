module bc_sub_in_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nee
  use eqn_idx, only : ncs,npres,ntemp,ncp,ngamm,nenth
  !
  implicit none
  !
  private
  !
  public :: get_bc_sub_inflow_solution
  !
contains
!
!###############################################################################
!
#ifdef DISABLE_FOR_GPU
pure &
#endif
subroutine get_bc_sub_inflow_solution(n1,n2,facenrm,face_cv,face_gp, &
                                      cv_in,pv_in,gp_in,pv_fs,gp_fs, &
                                      relax,newton_failures)
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
  use ovar, only : sub_inflow_method
  use ovar, only : relaxation_t
  !
  use bc_relax_mod, only : get_relaxed_total_bc
  !
  !.. Formal Arguments ..
  integer,            intent(in) :: n1,n2
  real(wp),           intent(in) :: facenrm(:,:)
  real(wp),           intent(in) :: cv_in(:)
  real(wp),           intent(in) :: pv_in(:)
  real(wp),           intent(in) :: gp_in(:)
  real(wp),           intent(in) :: pv_fs(:)
  real(wp),           intent(in) :: gp_fs(:)
  type(relaxation_t), intent(in) :: relax
  real(wp),        intent(inout) :: face_cv(:,:,:)
  real(wp),        intent(inout) :: face_gp(:,:,:)
  integer,         intent(inout) :: newton_failures
  !
  !.. Local Scalars ..
  integer  :: nr
  !
  !.. Local Arrays ..
  real(wp) :: pv_bc(1:size(pv_in))
  !
continue
  !
  newton_failures = 0
  !
  nr = size(facenrm,dim=1)
  !
  ! NEED TO GET PV_BC, WITH:
  !           PV_BC(1:NCS)   = MASS_FRACTIONS
  !           PV_BC(NEC)     = TOTAL DENSITY
  !           PV_BC(NMB:NME) = VELOCITIES
  !           PV_BC(NEE)     = TOTAL PRESSURE
  !
  if (ncs > 0 .or. sub_inflow_method /= 2) then
    call get_relaxed_total_bc(relax,pv_fs,gp_fs,pv_in,gp_in,pv_rlx=pv_bc)
  end if
  !
  if (ncs > 0) then
    !
    ! Compute the subsonic inflow boundary conditions
    ! for a thermally perfect gas mixture.
    !
    if (nr == 3) then
      call get_bc_sub_inflow_sol_tpg_3d(n1,n2,facenrm, &
                                        face_cv,face_gp,pv_bc, &
                                        newton_failures)
    else if (nr == 2) then
      call get_bc_sub_inflow_sol_tpg_2d(n1,n2,facenrm, &
                                        face_cv,face_gp,pv_bc, &
                                        newton_failures)
    else
      call get_bc_sub_inflow_sol_tpg_1d(n1,n2,facenrm, &
                                        face_cv,face_gp,pv_bc, &
                                        newton_failures)
    end if
    !
  else
    !
    ! Use the input parameter sub_inflow_method to determine which
    ! method to use to compute the subsonic inflow boundary conditions
    ! for a calorically perfect gas.
    !
    if (sub_inflow_method == 1) then
      !
      ! This method uses a Newton iteration on the total temperature at the
      ! inflow boundary based on the total temperature specified for this
      ! boundary condition.
      !
      if (nr == 3) then
        call get_bc_sub_inflow_newton_3d(n1,n2,facenrm,face_cv,face_gp, &
                                         pv_bc,newton_failures)
      else if (nr == 2) then
        call get_bc_sub_inflow_newton_2d(n1,n2,facenrm,face_cv,face_gp, &
                                         pv_bc,newton_failures)
      else
        call get_bc_sub_inflow_newton_1d(n1,n2,facenrm,face_cv,face_gp, &
                                         pv_bc,newton_failures)
      end if
      !
    else if (sub_inflow_method == 2) then
      !
      ! This method using the specified boundary condition density and
      ! velocities, and extrapolates interior pressure to the boundary.
      !
      if (nr == 3) then
        call get_bc_sub_inflow_sol_cpg_v2_3d(n1,n2,face_cv,cv_in)
      else if (nr == 2) then
        call get_bc_sub_inflow_sol_cpg_v2_2d(n1,n2,face_cv,cv_in)
      else
        call get_bc_sub_inflow_sol_cpg_v2_1d(n1,n2,face_cv,cv_in)
      end if
      !
    else
      !
      if (sub_inflow_method == 3) then
        !
        ! This method uses the method specified in NASA TM-2011-217181
        ! and was the old method used by VULCAN (it is now commented out).
        !   1) The interior pressure, total enthalpy, and the outward
        !      propagating Riemann invariant, R+, are extrapolated
        !      from the interior to the boundary.
        !   2) Uses the total pressure and total temperature specified
        !      for this boundary condition.
        !   3) The boundary speed of sound is found by solving a quadratic
        !      equation created by using some adiabatic and isentropic
        !      assumptions/relations with the conditions from 1) and 2).
        !   4) The boundary speed of sound is then used to compute the
        !      inflow velocity and Mach number.
        !   5) The inflow Mach number is used with the boundary total
        !      conditions to get the boundary static density, and total
        !      energy.
        !
        if (nr == 3) then
          call get_bc_sub_inflow_sol_cpg_v3_3d(n1,n2,facenrm, &
                                               face_cv,face_gp,pv_bc)
        else if (nr == 2) then
          call get_bc_sub_inflow_sol_cpg_v3_2d(n1,n2,facenrm, &
                                               face_cv,face_gp,pv_bc)
        else
          call get_bc_sub_inflow_sol_cpg_v3_1d(n1,n2,facenrm, &
                                               face_cv,face_gp,pv_bc)
        end if
        !
      else
        !
        ! This is the new method used by VULCAN in their subroutine sbibc.
        ! This method is similar to 'v3' above, except that it extrapolates
        ! the velocity normal to the boundary face from the interior to get
        ! the boundary Mach number. From there, the specified boundary total
        ! pressure, total temperature, and angle of attack can be used to
        ! back out the inflow conserved variables at the boundary.
        !
        if (nr == 3) then
          call get_bc_sub_inflow_sol_cgp_v4_3d(n1,n2,facenrm, &
                                               face_cv,face_gp,pv_bc)
        else if (nr == 2) then
          call get_bc_sub_inflow_sol_cgp_v4_2d(n1,n2,facenrm, &
                                               face_cv,face_gp,pv_bc)
        else
          call get_bc_sub_inflow_sol_cgp_v4_1d(n1,n2,facenrm, &
                                               face_cv,face_gp,pv_bc)
        end if
        !
      end if
      !
    end if
    !
  end if
  !
end subroutine get_bc_sub_inflow_solution
!
!###############################################################################
!
#define DTERR_USES_TEMP_EXT
#ifdef DISABLE_FOR_GPU
pure &
#endif
subroutine get_bc_sub_inflow_newton_3d(n1,n2,facenrm,face_cv,face_gp, &
                                       pv_bc,newton_failures)
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
  ! PV_BC   : Boundary condition flow variables in primitive form
  ! NEWTON_FAILURES : Number of failed convergences from Newton iterations
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd
  use ovar, only : tmin,tmax,tref
#ifndef DISABLE_FOR_GPU
  use ovar, only : rhoref,aref
  use geovar, only : face
#endif
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  integer,  intent(inout) :: newton_failures
  !
  !.. Local Scalars ..
  integer  :: k,n
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: vmag,dir_cosx,dir_cosy,dir_cosz
  real(wp) :: aoa,tan_aoa,boa,tan_boa
  real(wp) :: gm1,gam_over_1mg,one_over_gm1,cp_inv
  real(wp) :: rho_int,vx_int,vy_int,vz_int,q2_int
  real(wp) :: rho_ext,vx_ext,vy_ext,vz_ext,q2_ext
  real(wp) :: pres_int,temp_int,aspd_int,veln_int
  real(wp) :: pres_ext,temp_ext,aspd_ext,veln_ext
  real(wp) :: rho_inv,ttot,rminus,rafac,mach_ext
  real(wp) :: func,dfunc,dtemp,dtemp0,dterr,dterr0
  real(wp) :: rtot_bc,ptot_bc,ttot_bc,temp_bc,aspd_bc
  real(wp) :: tmin_nd,tmax_nd,conv_temp
  !
  logical(lk) :: convergence_reached
  !
  !.. Local Parameters ..
  integer, parameter :: newton_iter = 100
  !
#ifndef DISABLE_FOR_GPU
  integer :: i,nf,io,ierr,last_fp
  real(wp) :: rc,tc,rhoa2
  character(len=100) :: fname
  logical(lk) :: file_is_opened
  real(wp) :: dtemp_save(1:7,1:newton_iter)
  real(wp) :: xyz_fp(1:3)
#endif
  !
continue
  !
  tmin_nd = tmin/tref
  tmax_nd = tmax/tref
  !
#ifndef DISABLE_FOR_GPU
  file_is_opened = fals
  rhoa2 = rhoref*aref*aref
  last_fp = 0
#endif
  !
  ! Precompute some constants involving gamma and the reference gas constant
  !
  gm1 = gam - one
  one_over_gm1 = one / gm1
  gam_over_1mg = -gam * one_over_gm1
  cp_inv = gm1/(gam*rgasref_nd)
  !
  ! Get the directional cosines of the specified boundary condition velocity
  !
  vmag = sqrt( pv_bc(nmx)*pv_bc(nmx) + &
               pv_bc(nmy)*pv_bc(nmy) + &
               pv_bc(nmz)*pv_bc(nmz) )
  !
  dir_cosx = pv_bc(nmx) / vmag
  dir_cosy = pv_bc(nmy) / vmag
  dir_cosz = pv_bc(nmz) / vmag
  !
  aoa = atan2(dir_cosy,dir_cosx)
  tan_aoa = tan( aoa )
  !
  boa = atan2(dir_cosz,sqrt(dir_cosx*dir_cosx+dir_cosy*dir_cosy))
  tan_boa = tan( boa )
 !!
 !! Precompute the direction cosines
 !!
 !aoa = atan2(pv_bc(nmy),pv_bc(nmx))
 !tan_aoa = tan( aoa )
 !!
 !dir_cosx = one / sqrt(one + tan_aoa**2)
 !dir_cosy = tan_aoa * dir_cosx
 !dir_cosz = zero
  !
  ! Get the boundary total density and total pressure
  ! NOTE: These total conditions have been saved into static
  !       density and static pressure locations of pv_bc.
  !
  rtot_bc = pv_bc(nec)
  ptot_bc = pv_bc(nee)
  !
  ! Compute the boundary total temperature
  !
  ttot_bc = ptot_bc / (rtot_bc*rgasref_nd)
  !
  ! Compute the boundary static temperature by solving
  !     T0 = T * (1 + half*gm1*(mach**2))
  ! for T by replacing mach**2 with
  !     mach**2 = vmag**2 / aspd**2 = vmag**2 / (gam*rgasref_nd*T)
  ! and rearranging for T gives
  !     T = T0 - half*gm1*(vmag**2)/(gam*rgasref_nd)
  ! This can also be derived by using
  !     h = Cp*T = H0 - half*(vmag**2) = Cp*T0 - half*(vmag**2)
  ! Solving for T by dividing through by Cp gives
  !     T = T0 - half*(vmag**2)/Cp
  ! which is the same because Cp = gam*rgasref_nd/gm1.
  !
  temp_bc = ttot_bc - half*vmag*vmag*cp_inv
  !
  ! Compute the boundary speeed of sound
  !
  aspd_bc = sqrt(max(gam*temp_bc*rgasref_nd,rndoff))
  !
  !
#ifdef DTERR_USES_TEMP_EXT
  conv_temp = eps6
#else
  conv_temp = eps3*ttot_bc
#endif
  !
  newton_failures = 0
  !
  !$ACC PARALLEL LOOP GANG WORKER VECTOR &
  !$ACC&              PRESENT(facenrm,face_cv,face_gp) &
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
    rho_int = face_cv(n,nec,1)
    rho_inv = one/rho_int
    vx_int = face_cv(n,nmx,1)*rho_inv
    vy_int = face_cv(n,nmy,1)*rho_inv
    vz_int = face_cv(n,nmz,1)*rho_inv
    !
    pres_int = face_gp(n,npres,1)
    !
    temp_int = pres_int*rho_inv/rgasref_nd
    !
    ! Normal velocity of the interior/host cell at this flux point
    !
    veln_int = vx_int*unx + vy_int*uny + vz_int*unz
    !
    ! Compute the speed of sound at this flux point
    !
   !aspd_int = sqrt( max( gam*pres_int*rho_inv , rndoff ) )
    aspd_int = sqrt( max( gam*rgasref_nd*temp_int , rndoff ) )
    !
    rminus = veln_int - (two*aspd_int*one_over_gm1)
    !
    temp_ext = temp_bc
    !
    aspd_ext = aspd_bc
    !
    rafac = one / (dir_cosx*unx + dir_cosy*uny + dir_cosz*unz)
    !
    veln_ext = rafac * (rminus + two*aspd_ext*one_over_gm1)
    !
    ttot = temp_ext + half*veln_ext*veln_ext*cp_inv
    !
    func = ttot - ttot_bc
    dfunc = one + veln_ext*rafac/aspd_ext
    !
    convergence_reached = fals
    newton_iteration: do k = 1,newton_iter
      !
      dtemp = func/dfunc
      if (k == 1) then
        dtemp0 = dtemp
        dterr0 = get_dterr_func(dtemp0,temp_ext)
      end if
      !
      dterr = get_dterr_func(dtemp,temp_ext)
      !
#ifndef DISABLE_FOR_GPU
      dtemp_save(1,k) = dtemp
      dtemp_save(2,k) = dterr
      dtemp_save(3,k) = func
      dtemp_save(4,k) = dfunc
      dtemp_save(5,k) = temp_ext
      dtemp_save(6,k) = aspd_ext
      dtemp_save(7,k) = veln_ext
#endif
      !
      if (dterr < conv_temp) then
        convergence_reached = true
      end if
#ifndef DTERR_USES_TEMP_EXT
      if (dterr < ten*conv_temp .and. dterr0/dterr < eps1) then
        convergence_reached = true
      end if
#endif
      !
      if (convergence_reached) then
        exit newton_iteration
      end if
      !
      temp_ext = min( tmax_nd , max(tmin_nd,temp_ext-dtemp) )
      !
      aspd_ext = sqrt( max( gam*rgasref_nd*temp_ext , rndoff ) )
      veln_ext = rafac * (rminus + two*aspd_ext*one_over_gm1)
      !
      ttot = temp_ext + half*veln_ext*veln_ext*cp_inv
      !
      func = ttot - ttot_bc
      dfunc = one + veln_ext*rafac/aspd_ext
      !
     !func = ttot_bc - temp_ext - half*veln_ext*veln_ext*cp_inv
     !dfunc = -(one + veln_ext*rafac/aspd_ext)
      !
    end do newton_iteration
    !
    mach_ext = veln_ext / aspd_ext
    !
    pres_ext = ptot_bc*(one + half*gm1*mach_ext*mach_ext)**gam_over_1mg
    !
    rho_ext = pres_ext / temp_ext / rgasref_nd
    !
    ! Set the boundary velocity in the same direction
    ! as the specified boundary condition velocity
    !
    vx_ext = veln_ext*dir_cosx
    vy_ext = veln_ext*dir_cosy
    vz_ext = veln_ext*dir_cosz
    !
    ! Set the boundary velocity normal to the boundary face
    !
   !vx_ext = veln_ext*unx
   !vy_ext = veln_ext*uny
   !vz_ext = veln_ext*unz
    !
    face_cv(n,nec,2) = rho_ext
    face_cv(n,nmx,2) = rho_ext*vx_ext
    face_cv(n,nmy,2) = rho_ext*vy_ext
    face_cv(n,nmz,2) = rho_ext*vz_ext
    face_cv(n,nee,2) = pres_ext*one_over_gm1 + half*rho_ext*veln_ext*veln_ext
    !
    ! Check to make sure the temperature converged
    !
    if (.not. convergence_reached) then
      newton_failures = newton_failures + 1
#ifndef DISABLE_FOR_GPU
      dtemp = func/dfunc
      dterr = get_dterr_func(dtemp,temp_ext)
      find_face_loop: do k = 1,size(face)
        if (n <= face(k)%end_fp) then
          nf = k
          exit find_face_loop
        end if
      end do find_face_loop
      last_fp = face(nf)%end_fp
      k = n - face(nf)%beg_fp + 1
      if (.not. file_is_opened) then
        write (fname,5) mypnum
        open (newunit=io,file=trim(adjustl(fname)),status="replace", &
                         action="write",form="formatted",iostat=ierr)
        file_is_opened = true
      end if
      call get_xyz_for_bnd_face_point(nf,k,xyz_fp)
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
      write (io,2) "boa      = ",boa
      write (io,2) "tan_boa  = ",tan_boa
      write (io,2) "dir_cos1 = ",dir_cosx
      write (io,2) "dir_cos2 = ",dir_cosy
      write (io,2) "dir_cos3 = ",dir_cosz
      write (io,2) "rafac    = ",rafac,      " "
      write (io,6) "Non-dimensional values"
      write (io,2) "dens_int  = ",rho_int
      write (io,2) "vx_int    = ",vx_int
      write (io,2) "vy_int    = ",vy_int
      write (io,2) "vz_int    = ",vz_int
      write (io,2) "pres_int  = ",pres_int
      write (io,2) "temp_int  = ",temp_int
     !write (io,2) "temp_int  = ",pres_int/rho_int/rgasref_nd
      write (io,2) "aspd_int  = ",aspd_int
      write (io,2) "veln_int  = ",veln_int
      write (io,2) "rminus    = ",rminus,     " "
      write (io,2) "conv_temp = ",conv_temp
      write (io,2) "func      = ",func
      write (io,2) "dfunc     = ",dfunc
      write (io,2) "dtemp     = ",dtemp
      write (io,2) "dterr     = ",dterr,      " "
      write (io,2) "dens_ext  = ",rho_ext
      write (io,2) "vx_ext    = ",vx_ext
      write (io,2) "vy_ext    = ",vy_ext
      write (io,2) "vz_ext    = ",vz_ext
      write (io,2) "pres_ext  = ",pres_ext
      write (io,2) "temp_ext  = ",temp_ext
      write (io,2) "aspd_ext  = ",aspd_ext
      write (io,2) "veln_ext  = ",veln_ext
      write (io,2) "mach_ext  = ",mach_ext,   " "
      write (io,6) "Dimensional values"
      write (io,7) "dens_int  = ",rho_int*rhoref,           " kg/m^3"
      write (io,7) "vx_int    = ",vx_int*aref,              " m/s"
      write (io,7) "vy_int    = ",vy_int*aref,              " m/s"
      write (io,7) "vz_int    = ",vz_int*aref,              " m/s"
      write (io,7) "pres_int  = ",pres_int*rhoa2,           " Pa"
      write (io,7) "temp_int  = ",temp_int*tref,            " K"
      write (io,7) "aspd_int  = ",aspd_int*aref,            " m/s"
      write (io,7) "veln_int  = ",veln_int*aref,            " m/s"
      write (io,7) "rminus    = ",rminus*aref,              " m/s",     " "
      write (io,7) "conv_temp = ",conv_temp*tref,           " K"
      write (io,7) "func      = ",func*tref,                " K"
      write (io,7) "dtemp     = ",dtemp*tref,               " K"
      write (io,7) "dterr     = ",dterr*tref,               " K",       " "
      write (io,7) "dens_ext  = ",rho_ext*rhoref,           " kg/m^3"
      write (io,7) "vx_ext    = ",vx_ext*aref,              " m/s"
      write (io,7) "vy_ext    = ",vy_ext*aref,              " m/s"
      write (io,7) "vz_ext    = ",vz_ext*aref,              " m/s"
      write (io,7) "pres_ext  = ",pres_ext*rhoa2,           " Pa"
      write (io,7) "temp_ext  = ",temp_ext*tref,            " K"
      write (io,7) "aspd_ext  = ",aspd_ext*aref,            " m/s"
      write (io,7) "veln_ext  = ",veln_ext*aref,            " m/s",     " "
      do k = 1,size(dtemp_save,dim=2)
        if (k == 1) write (io,3)
        write (io,4) k,(dtemp_save(i,k), i=1,size(dtemp_save,dim=1))
        if (k == size(dtemp_save,dim=2)) then
          write (io,4) k+1,dtemp,dterr,func,dfunc,temp_ext,aspd_ext,veln_ext,"-"
        end if
      end do
      flush (io)
#endif
    end if
    !
#ifndef DISABLE_FOR_GPU
    ! If this is the last face point on this face, output the final
    ! BC state (primitive variables) for all the face points on this
    ! boundary face
    ! NOTE: last_fp will only have a non-zero value if convergence_reached
    !       was false for one of the face points of the current face.
    !       This means the following will only be written if the above
    !       Newton iteration failed to converge for any of the face points.
    !
    if (n == last_fp) then
      write (io,8) nf
      do k = face(nf)%beg_fp,face(nf)%end_fp
        !
        i = k-face(nf)%beg_fp+1
        !
        call get_xyz_for_bnd_face_point(nf,i,xyz_fp)
        !
        rho_int = face_cv(k,nec,1)
        vx_int = face_cv(k,nmx,1)/rho_int
        vy_int = face_cv(k,nmy,1)/rho_int
        vz_int = face_cv(k,nmz,1)/rho_int
        q2_int = vx_int*vx_int + vy_int*vy_int + vz_int*vz_int
        pres_int = gm1*(face_cv(k,nee,1) - half*rho_int*q2_int)
        !
        rho_ext = face_cv(k,nec,2)
        vx_ext = face_cv(k,nmx,2)/rho_ext
        vy_ext = face_cv(k,nmy,2)/rho_ext
        vz_ext = face_cv(k,nmz,2)/rho_ext
        q2_ext = vx_ext*vx_ext + vy_ext*vy_ext + vz_ext*vz_ext
        pres_ext = gm1*(face_cv(k,nee,2) - half*rho_ext*q2_ext)
        !
        write (io,9) k,xyz_fp(1),xyz_fp(2),xyz_fp(3), &
                 "rho",  rho_int,  rho_ext,  rho_int*rhoref, rho_ext*rhoref, &
                  "vx",   vx_int,   vx_ext,   vx_int*aref,    vx_ext*aref, &
                  "vy",   vy_int,   vy_ext,   vy_int*aref,    vy_ext*aref, &
                  "vz",   vz_int,   vz_ext,   vz_int*aref,    vz_ext*aref, &
                "pres", pres_int, pres_ext, pres_int*rhoa2, pres_ext*rhoa2
        !
      end do
      write (io,8)
      flush (io)
    end if
#endif
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
  2 format (a,g15.8,:,/,a)
  3 format (/,101("-"),/, &
              4x,"VALUES AFFECTED BY EACH NEWTON ITERATION",/, &
              101("-"),/, &
              1x,"Iter",7x,"dtemp",11x,"dterr",11x,"func",12x,"dfunc", &
              9x,"temp_ext",8x,"aspd_ext",8x,"veln_ext",/, &
              101("-"))
  4 format (2x,i3,7(2x,es14.6),:,/,100("-"),a,/)
  5 format ("newton_bc_errors.CPU-",i0,".log")
  6 format (4x,a,/)
  7 format (a,g15.8,a,:,/,a)
  8 format (/,101("="),/,:, &
              4x,"Final states for the face points on face ",i0,/, &
              101("="),/, &
              10x,"Non-Dim Int State",3x,"Non-Dim Ext State", &
              12x,"Dim Int State",5x,"Dim Ext State",/, &
              101("="))
  9 format (/,5x,"Flux point # ",i0, &
              t49,"[x,y,z] = [",es13.6,",",es13.6,",",es13.6,"]",/, &
              1x,a,t11,g15.8,t31,g15.8,t58,g15.8,t76,g15.8,/, &
            3(1x,a,t11,g15.8,t31,g15.8,t58,g15.8,t76,g15.8,/), &
              1x,a,t11,g15.8,t31,g15.8,t58,g15.8,t76,g15.8)
  !
contains
  !
  pure function get_dterr_func(a,b) result(return_value)
    real(wp), intent(in) :: a
    real(wp), intent(in) :: b
    real(wp) :: return_value
  continue
#ifdef DTERR_USES_TEMP_EXT
    return_value = abs(a)/b
#else
    return_value = abs(a)
#endif
  end function get_dterr_func
  !
end subroutine get_bc_sub_inflow_newton_3d
!
!###############################################################################
!
#ifdef DISABLE_FOR_GPU
pure &
#endif
subroutine get_bc_sub_inflow_newton_2d(n1,n2,facenrm,face_cv,face_gp, &
                                       pv_bc,newton_failures)
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
  ! PV_BC   : Boundary condition flow variables in primitive form
  ! NEWTON_FAILURES : Number of failed convergences from Newton iterations
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd
  use ovar, only : tmin,tmax,tref
#ifndef DISABLE_FOR_GPU
  use ovar, only : rhoref,aref
  use geovar, only : face
#endif
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  integer,  intent(inout) :: newton_failures
  !
  !.. Local Scalars ..
  integer  :: k,n
  real(wp) :: ds,rnx,rny,unx,uny
  real(wp) :: vmag,dir_cosx,dir_cosy,aoa,tan_aoa
  real(wp) :: gm1,gam_over_1mg,one_over_gm1,cp_inv
  real(wp) :: rho_int,vx_int,vy_int,q2_int
  real(wp) :: rho_ext,vx_ext,vy_ext,q2_ext
  real(wp) :: pres_int,temp_int,aspd_int,veln_int
  real(wp) :: pres_ext,temp_ext,aspd_ext,veln_ext
  real(wp) :: rho_inv,ttot,rminus,rafac,mach_ext
  real(wp) :: func,dfunc,dtemp,dtemp0,dterr,dterr0
  real(wp) :: rtot_bc,ptot_bc,ttot_bc,temp_bc,aspd_bc
  real(wp) :: tmin_nd,tmax_nd,conv_temp
  !
  logical(lk) :: convergence_reached
  !
  !.. Local Parameters ..
  integer, parameter :: newton_iter = 100
  !
#ifndef DISABLE_FOR_GPU
  integer :: i,nf,io,ierr,last_fp
  real(wp) :: rhoa2
  character(len=100) :: fname
  logical(lk) :: file_is_opened
  real(wp) :: dtemp_save(1:7,1:newton_iter)
  real(wp) :: xyz_fp(1:2)
#endif
  !
continue
  !
  tmin_nd = tmin/tref
  tmax_nd = tmax/tref
  !
#ifndef DISABLE_FOR_GPU
  file_is_opened = fals
  rhoa2 = rhoref*aref*aref
  last_fp = 0
#endif
  !
  ! Precompute some constants involving gamma and the reference gas constant
  !
  gm1 = gam - one
  one_over_gm1 = one / gm1
  gam_over_1mg = -gam * one_over_gm1
  cp_inv = gm1/(gam*rgasref_nd)
  !
  ! Get the directional cosines of the specified boundary condition velocity
  !
  vmag = sqrt( pv_bc(nmx)*pv_bc(nmx) + &
               pv_bc(nmy)*pv_bc(nmy) )
  !
  dir_cosx = pv_bc(nmx) / vmag
  dir_cosy = pv_bc(nmy) / vmag
  !
  aoa = atan2(dir_cosy,dir_cosx)
  tan_aoa = tan( aoa )
 !!
 !! Precompute the direction cosines
 !!
 !aoa = atan2(pv_bc(nmy),pv_bc(nmx))
 !tan_aoa = tan( aoa )
 !!
 !dir_cosx = one / sqrt(one + tan_aoa**2)
 !dir_cosy = tan_aoa * dir_cosx
  !
  ! Get the boundary total density and total pressure
  ! NOTE: These total conditions have been saved into static
  !       density and static pressure locations of pv_bc.
  !
  rtot_bc = pv_bc(nec)
  ptot_bc = pv_bc(nee)
  !
  ! Compute the boundary total temperature
  !
  ttot_bc = ptot_bc / (rtot_bc*rgasref_nd)
  !
  ! Compute the boundary static temperature by solving
  !     T0 = T * (1 + half*gm1*(mach**2))
  ! for T by replacing mach**2 with
  !     mach**2 = vmag**2 / aspd**2 = vmag**2 / (gam*rgasref_nd*T)
  ! and rearranging for T gives
  !     T = T0 - half*gm1*(vmag**2)/(gam*rgasref_nd)
  ! This can also be derived by using
  !     h = Cp*T = H0 - half*(vmag**2) = Cp*T0 - half*(vmag**2)
  ! Solving for T by dividing through by Cp gives
  !     T = T0 - half*(vmag**2)/Cp
  ! which is the same because Cp = gam*rgasref_nd/gm1.
  !
  temp_bc = ttot_bc - half*vmag*vmag*cp_inv
  !
  ! Compute the boundary speeed of sound
  !
  aspd_bc = sqrt(max(gam*temp_bc*rgasref_nd,rndoff))
  !
  !
#ifdef DTERR_USES_TEMP_EXT
  conv_temp = eps6
#else
  conv_temp = eps3*ttot_bc
#endif
  !
  newton_failures = 0
  !
  !$ACC PARALLEL LOOP GANG WORKER VECTOR &
  !$ACC&              PRESENT(facenrm,face_cv,face_gp) &
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
    rho_int = face_cv(n,nec,1)
    rho_inv = one/rho_int
    vx_int = face_cv(n,nmx,1)*rho_inv
    vy_int = face_cv(n,nmy,1)*rho_inv
    !
    pres_int = face_gp(n,npres,1)
    !
    temp_int = pres_int*rho_inv/rgasref_nd
    !
    ! Normal velocity of the interior/host cell at this flux point
    !
    veln_int = vx_int*unx + vy_int*uny
    !
    ! Compute the speed of sound at this flux point
    !
   !aspd_int = sqrt( max( gam*pres_int*rho_inv , rndoff ) )
    aspd_int = sqrt( max( gam*rgasref_nd*temp_int , rndoff ) )
    !
    rminus = veln_int - (two*aspd_int*one_over_gm1)
    !
    temp_ext = temp_bc
    !
    aspd_ext = aspd_bc
    !
    rafac = one / (dir_cosx*unx + dir_cosy*uny)
    !
    veln_ext = rafac * (rminus + two*aspd_ext*one_over_gm1)
    !
    ttot = temp_ext + half*veln_ext*veln_ext*cp_inv
    !
    func = ttot - ttot_bc
    dfunc = one + veln_ext*rafac/aspd_ext
    !
    convergence_reached = fals
    newton_iteration: do k = 1,newton_iter
      !
      dtemp = func/dfunc
      if (k == 1) then
        dtemp0 = dtemp
        dterr0 = get_dterr_func(dtemp0,temp_ext)
      end if
      !
      dterr = get_dterr_func(dtemp,temp_ext)
      !
#ifndef DISABLE_FOR_GPU
      dtemp_save(1,k) = dtemp
      dtemp_save(2,k) = dterr
      dtemp_save(3,k) = func
      dtemp_save(4,k) = dfunc
      dtemp_save(5,k) = temp_ext
      dtemp_save(6,k) = aspd_ext
      dtemp_save(7,k) = veln_ext
#endif
      !
      if (dterr < conv_temp) then
        convergence_reached = true
      end if
#ifndef DTERR_USES_TEMP_EXT
      if (dterr < ten*conv_temp .and. dterr0/dterr < eps1) then
        convergence_reached = true
      end if
#endif
      !
      if (convergence_reached) then
        exit newton_iteration
      end if
      !
      temp_ext = min( tmax_nd , max(tmin_nd,temp_ext-dtemp) )
      !
      aspd_ext = sqrt( max( gam*rgasref_nd*temp_ext , rndoff ) )
      veln_ext = rafac * (rminus + two*aspd_ext*one_over_gm1)
      !
      ttot = temp_ext + half*veln_ext*veln_ext*cp_inv
      !
      func = ttot - ttot_bc
      dfunc = one + veln_ext*rafac/aspd_ext
      !
     !func = ttot_bc - temp_ext - half*veln_ext*veln_ext*cp_inv
     !dfunc = -(one + veln_ext*rafac/aspd_ext)
      !
    end do newton_iteration
    !
    mach_ext = veln_ext / aspd_ext
    !
    pres_ext = ptot_bc*(one + half*gm1*mach_ext*mach_ext)**gam_over_1mg
    !
    rho_ext = pres_ext / temp_ext / rgasref_nd
    !
    ! Set the boundary velocity in the same direction
    ! as the specified boundary condition velocity
    !
    vx_ext = veln_ext*dir_cosx
    vy_ext = veln_ext*dir_cosy
    !
    ! Set the boundary velocity normal to the boundary face
    !
   !vx_ext = veln_ext*unx
   !vy_ext = veln_ext*uny
    !
    face_cv(n,nec,2) = rho_ext
    face_cv(n,nmx,2) = rho_ext*vx_ext
    face_cv(n,nmy,2) = rho_ext*vy_ext
    face_cv(n,nee,2) = pres_ext*one_over_gm1 + half*rho_ext*veln_ext*veln_ext
    !
    ! Check to make sure the temperature converged
    !
    if (.not. convergence_reached) then
      newton_failures = newton_failures + 1
#ifndef DISABLE_FOR_GPU
      dtemp = func/dfunc
      dterr = get_dterr_func(dtemp,temp_ext)
      find_face_loop: do k = 1,size(face)
        if (n <= face(k)%end_fp) then
          nf = k
          exit find_face_loop
        end if
      end do find_face_loop
      last_fp = face(nf)%end_fp
      k = n - face(nf)%beg_fp + 1
      if (.not. file_is_opened) then
        write (fname,5) mypnum
        open (newunit=io,file=trim(adjustl(fname)),status="replace", &
                         action="write",form="formatted",iostat=ierr)
        file_is_opened = true
      end if
      call get_xyz_for_bnd_face_point(nf,k,xyz_fp)
      write (io,1) "mypnum   = ",mypnum
      write (io,1) "glb_fp   = ",n
      write (io,1) "nf       = ",nf
      write (io,1) "loc_fp   = ",k
      write (io,2) "X        = ",xyz_fp(1)
      write (io,2) "Y        = ",xyz_fp(2),  " "
      write (io,2) "aoa      = ",aoa
      write (io,2) "tan_aoa  = ",tan_aoa
      write (io,2) "dir_cos1 = ",dir_cosx
      write (io,2) "dir_cos2 = ",dir_cosy
      write (io,2) "rafac    = ",rafac,      " "
      write (io,6) "Non-dimensional values"
      write (io,2) "dens_int  = ",rho_int
      write (io,2) "vx_int    = ",vx_int
      write (io,2) "vy_int    = ",vy_int
      write (io,2) "pres_int  = ",pres_int
      write (io,2) "temp_int  = ",temp_int
     !write (io,2) "temp_int  = ",pres_int/rho_int/rgasref_nd
      write (io,2) "aspd_int  = ",aspd_int
      write (io,2) "veln_int  = ",veln_int
      write (io,2) "rminus    = ",rminus,     " "
      write (io,2) "conv_temp = ",conv_temp
      write (io,2) "func      = ",func
      write (io,2) "dfunc     = ",dfunc
      write (io,2) "dtemp     = ",dtemp
      write (io,2) "dterr     = ",dterr,      " "
      write (io,2) "dens_ext  = ",rho_ext
      write (io,2) "vx_ext    = ",vx_ext
      write (io,2) "vy_ext    = ",vy_ext
      write (io,2) "pres_ext  = ",pres_ext
      write (io,2) "temp_ext  = ",temp_ext
      write (io,2) "aspd_ext  = ",aspd_ext
      write (io,2) "veln_ext  = ",veln_ext
      write (io,2) "mach_ext  = ",mach_ext,   " "
      write (io,6) "Dimensional values"
      write (io,7) "dens_int  = ",rho_int*rhoref,           " kg/m^3"
      write (io,7) "vx_int    = ",vx_int*aref,              " m/s"
      write (io,7) "vy_int    = ",vy_int*aref,              " m/s"
      write (io,7) "pres_int  = ",pres_int*rhoa2,           " Pa"
      write (io,7) "temp_int  = ",temp_int*tref,            " K"
      write (io,7) "aspd_int  = ",aspd_int*aref,            " m/s"
      write (io,7) "veln_int  = ",veln_int*aref,            " m/s"
      write (io,7) "rminus    = ",rminus*aref,              " m/s",     " "
      write (io,7) "conv_temp = ",conv_temp*tref,           " K"
      write (io,7) "func      = ",func*tref,                " K"
      write (io,7) "dtemp     = ",dtemp*tref,               " K"
      write (io,7) "dterr     = ",dterr*tref,               " K",       " "
      write (io,7) "dens_ext  = ",rho_ext*rhoref,           " kg/m^3"
      write (io,7) "vx_ext    = ",vx_ext*aref,              " m/s"
      write (io,7) "vy_ext    = ",vy_ext*aref,              " m/s"
      write (io,7) "pres_ext  = ",pres_ext*rhoa2,           " Pa"
      write (io,7) "temp_ext  = ",temp_ext*tref,            " K"
      write (io,7) "aspd_ext  = ",aspd_ext*aref,            " m/s"
      write (io,7) "veln_ext  = ",veln_ext*aref,            " m/s",     " "
      do k = 1,size(dtemp_save,dim=2)
        if (k == 1) write (io,3)
        write (io,4) k,(dtemp_save(i,k), i=1,size(dtemp_save,dim=1))
        if (k == size(dtemp_save,dim=2)) then
          write (io,4) k+1,dtemp,dterr,func,dfunc,temp_ext,aspd_ext,veln_ext,"-"
        end if
      end do
      flush (io)
#endif
    end if
    !
#ifndef DISABLE_FOR_GPU
    ! If this is the last face point on this face, output the final
    ! BC state (primitive variables) for all the face points on this
    ! boundary face
    ! NOTE: last_fp will only have a non-zero value if convergence_reached
    !       was false for one of the face points of the current face.
    !       This means the following will only be written if the above
    !       Newton iteration failed to converge for any of the face points.
    !
    if (n == last_fp) then
      write (io,8) nf
      do k = face(nf)%beg_fp,face(nf)%end_fp
        !
        i = k-face(nf)%beg_fp+1
        !
        call get_xyz_for_bnd_face_point(nf,i,xyz_fp)
        !
        rho_int = face_cv(k,nec,1)
        vx_int = face_cv(k,nmx,1)/rho_int
        vy_int = face_cv(k,nmy,1)/rho_int
        q2_int = vx_int*vx_int + vy_int*vy_int
        pres_int = gm1*(face_cv(k,nee,1) - half*rho_int*q2_int)
        !
        rho_ext = face_cv(k,nec,2)
        vx_ext = face_cv(k,nmx,2)/rho_ext
        vy_ext = face_cv(k,nmy,2)/rho_ext
        q2_ext = vx_ext*vx_ext + vy_ext*vy_ext
        pres_ext = gm1*(face_cv(k,nee,2) - half*rho_ext*q2_ext)
        !
        write (io,9) k,xyz_fp(1),xyz_fp(2), &
                 "rho",  rho_int,  rho_ext,  rho_int*rhoref, rho_ext*rhoref, &
                  "vx",   vx_int,   vx_ext,   vx_int*aref,    vx_ext*aref, &
                  "vy",   vy_int,   vy_ext,   vy_int*aref,    vy_ext*aref, &
                "pres", pres_int, pres_ext, pres_int*rhoa2, pres_ext*rhoa2
        !
      end do
      write (io,8)
      flush (io)
    end if
#endif
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
  2 format (a,g15.8,:,/,a)
  3 format (/,101("-"),/, &
              4x,"VALUES AFFECTED BY EACH NEWTON ITERATION",/, &
              101("-"),/, &
              1x,"Iter",7x,"dtemp",11x,"dterr",11x,"func",12x,"dfunc", &
              9x,"temp_ext",8x,"aspd_ext",8x,"veln_ext",/, &
              101("-"))
  4 format (2x,i3,7(2x,es14.6),:,/,100("-"),a,/)
  5 format ("newton_bc_errors.CPU-",i0,".log")
  6 format (4x,a,/)
  7 format (a,g15.8,a,:,/,a)
  8 format (/,101("="),/,:, &
              4x,"Final states for the face points on face ",i0,/, &
              101("="),/, &
              10x,"Non-Dim Int State",3x,"Non-Dim Ext State", &
              12x,"Dim Int State",5x,"Dim Ext State",/, &
              101("="))
  9 format (/,5x,"Flux point # ",i0, &
              t49,"[x,y] = [",es13.6,",",es13.6,"]",/, &
              1x,a,t11,g15.8,t31,g15.8,t58,g15.8,t76,g15.8,/, &
            2(1x,a,t11,g15.8,t31,g15.8,t58,g15.8,t76,g15.8,/), &
              1x,a,t11,g15.8,t31,g15.8,t58,g15.8,t76,g15.8)
  !
contains
  !
  pure function get_dterr_func(a,b) result(return_value)
    real(wp), intent(in) :: a
    real(wp), intent(in) :: b
    real(wp) :: return_value
  continue
#ifdef DTERR_USES_TEMP_EXT
    return_value = abs(a)/b
#else
    return_value = abs(a)
#endif
  end function get_dterr_func
  !
end subroutine get_bc_sub_inflow_newton_2d
!
!###############################################################################
!
#ifdef DISABLE_FOR_GPU
pure &
#endif
subroutine get_bc_sub_inflow_newton_1d(n1,n2,facenrm,face_cv,face_gp, &
                                       pv_bc,newton_failures)
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
  ! PV_BC   : Boundary condition flow variables in primitive form
  ! NEWTON_FAILURES : Number of failed convergences from Newton iterations
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd
  use ovar, only : tmin,tmax,tref
#ifndef DISABLE_FOR_GPU
  use ovar, only : rhoref,aref
  use geovar, only : face
#endif
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  integer,  intent(inout) :: newton_failures
  !
  !.. Local Scalars ..
  integer  :: k,n
  real(wp) :: ds,rnx,unx
  real(wp) :: vmag,dir_cosx
  real(wp) :: gm1,gam_over_1mg,one_over_gm1,cp_inv
  real(wp) :: rho_int,vx_int
  real(wp) :: rho_ext,vx_ext
  real(wp) :: pres_int,temp_int,aspd_int,veln_int
  real(wp) :: pres_ext,temp_ext,aspd_ext,veln_ext
  real(wp) :: rho_inv,ttot,rminus,rafac,mach_ext
  real(wp) :: func,dfunc,dtemp,dtemp0,dterr,dterr0
  real(wp) :: rtot_bc,ptot_bc,ttot_bc,temp_bc,aspd_bc
  real(wp) :: tmin_nd,tmax_nd,conv_temp
  !
  logical(lk) :: convergence_reached
  !
  !.. Local Parameters ..
  integer, parameter :: newton_iter = 100
  !
#ifndef DISABLE_FOR_GPU
  integer :: i,nf,io,ierr,last_fp
  real(wp) :: rhoa2
  character(len=100) :: fname
  logical(lk) :: file_is_opened
  real(wp) :: dtemp_save(1:7,1:newton_iter)
  real(wp) :: xyz_fp(1:1)
#endif
  !
continue
  !
  tmin_nd = tmin/tref
  tmax_nd = tmax/tref
  !
#ifndef DISABLE_FOR_GPU
  file_is_opened = fals
  rhoa2 = rhoref*aref*aref
  last_fp = 0
#endif
  !
  ! Precompute some constants involving gamma and the reference gas constant
  !
  gm1 = gam - one
  one_over_gm1 = one / gm1
  gam_over_1mg = -gam * one_over_gm1
  cp_inv = gm1/(gam*rgasref_nd)
  !
  ! Get the directional cosines of the specified boundary condition velocity
  !
  vmag = sqrt( pv_bc(nmx)*pv_bc(nmx) )
  !
  dir_cosx = pv_bc(nmx) / vmag
  !
  ! Get the boundary total density and total pressure
  ! NOTE: These total conditions have been saved into static
  !       density and static pressure locations of pv_bc.
  !
  rtot_bc = pv_bc(nec)
  ptot_bc = pv_bc(nee)
  !
  ! Compute the boundary total temperature
  !
  ttot_bc = ptot_bc / (rtot_bc*rgasref_nd)
  !
  ! Compute the boundary static temperature by solving
  !     T0 = T * (1 + half*gm1*(mach**2))
  ! for T by replacing mach**2 with
  !     mach**2 = vmag**2 / aspd**2 = vmag**2 / (gam*rgasref_nd*T)
  ! and rearranging for T gives
  !     T = T0 - half*gm1*(vmag**2)/(gam*rgasref_nd)
  ! This can also be derived by using
  !     h = Cp*T = H0 - half*(vmag**2) = Cp*T0 - half*(vmag**2)
  ! Solving for T by dividing through by Cp gives
  !     T = T0 - half*(vmag**2)/Cp
  ! which is the same because Cp = gam*rgasref_nd/gm1.
  !
  temp_bc = ttot_bc - half*vmag*vmag*cp_inv
  !
  ! Compute the boundary speeed of sound
  !
  aspd_bc = sqrt(max(gam*temp_bc*rgasref_nd,rndoff))
  !
  !
#ifdef DTERR_USES_TEMP_EXT
  conv_temp = eps6
#else
  conv_temp = eps3*ttot_bc
#endif
  !
  newton_failures = 0
  !
  !$ACC PARALLEL LOOP GANG WORKER VECTOR &
  !$ACC&              PRESENT(facenrm,face_cv,face_gp) &
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
    rho_int = face_cv(n,nec,1)
    rho_inv = one/rho_int
    vx_int = face_cv(n,nmx,1)*rho_inv
    !
    pres_int = face_gp(n,npres,1)
    !
    temp_int = pres_int*rho_inv/rgasref_nd
    !
    ! Normal velocity of the interior/host cell at this flux point
    !
    veln_int = vx_int*unx
    !
    ! Compute the speed of sound at this flux point
    !
   !aspd_int = sqrt( max( gam*pres_int*rho_inv , rndoff ) )
    aspd_int = sqrt( max( gam*rgasref_nd*temp_int , rndoff ) )
    !
    rminus = veln_int - (two*aspd_int*one_over_gm1)
    !
    temp_ext = temp_bc
    !
    aspd_ext = aspd_bc
    !
    rafac = one / (dir_cosx*unx)
    !
    veln_ext = rafac * (rminus + two*aspd_ext*one_over_gm1)
    !
    ttot = temp_ext + half*veln_ext*veln_ext*cp_inv
    !
    func = ttot - ttot_bc
    dfunc = one + veln_ext*rafac/aspd_ext
    !
    convergence_reached = fals
    newton_iteration: do k = 1,newton_iter
      !
      dtemp = func/dfunc
      if (k == 1) then
        dtemp0 = dtemp
        dterr0 = get_dterr_func(dtemp0,temp_ext)
      end if
      !
      dterr = get_dterr_func(dtemp,temp_ext)
      !
#ifndef DISABLE_FOR_GPU
      dtemp_save(1,k) = dtemp
      dtemp_save(2,k) = dterr
      dtemp_save(3,k) = func
      dtemp_save(4,k) = dfunc
      dtemp_save(5,k) = temp_ext
      dtemp_save(6,k) = aspd_ext
      dtemp_save(7,k) = veln_ext
#endif
      !
      if (dterr < conv_temp) then
        convergence_reached = true
      end if
#ifndef DTERR_USES_TEMP_EXT
      if (dterr < ten*conv_temp .and. dterr0/dterr < eps1) then
        convergence_reached = true
      end if
#endif
      !
      if (convergence_reached) then
        exit newton_iteration
      end if
      !
      temp_ext = min( tmax_nd , max(tmin_nd,temp_ext-dtemp) )
      !
      aspd_ext = sqrt( max( gam*rgasref_nd*temp_ext , rndoff ) )
      veln_ext = rafac * (rminus + two*aspd_ext*one_over_gm1)
      !
      ttot = temp_ext + half*veln_ext*veln_ext*cp_inv
      !
      func = ttot - ttot_bc
      dfunc = one + veln_ext*rafac/aspd_ext
      !
     !func = ttot_bc - temp_ext - half*veln_ext*veln_ext*cp_inv
     !dfunc = -(one + veln_ext*rafac/aspd_ext)
      !
    end do newton_iteration
    !
    mach_ext = veln_ext / aspd_ext
    !
    pres_ext = ptot_bc*(one + half*gm1*mach_ext*mach_ext)**gam_over_1mg
    !
    rho_ext = pres_ext / temp_ext / rgasref_nd
    !
    ! Set the boundary velocity in the same direction
    ! as the specified boundary condition velocity
    !
    vx_ext = veln_ext*dir_cosx
    !
    ! Set the boundary velocity normal to the boundary face
    !
   !vx_ext = veln_ext*unx
    !
    face_cv(n,nec,2) = rho_ext
    face_cv(n,nmx,2) = rho_ext*vx_ext
    face_cv(n,nee,2) = pres_ext*one_over_gm1 + half*rho_ext*veln_ext*veln_ext
    !
    ! Check to make sure the temperature converged
    !
    if (.not. convergence_reached) then
      newton_failures = newton_failures + 1
#ifndef DISABLE_FOR_GPU
      dtemp = func/dfunc
      dterr = get_dterr_func(dtemp,temp_ext)
      find_face_loop: do k = 1,size(face)
        if (n <= face(k)%end_fp) then
          nf = k
          exit find_face_loop
        end if
      end do find_face_loop
      last_fp = face(nf)%end_fp
      k = n - face(nf)%beg_fp + 1
      if (.not. file_is_opened) then
        write (fname,5) mypnum
        open (newunit=io,file=trim(adjustl(fname)),status="replace", &
                         action="write",form="formatted",iostat=ierr)
        file_is_opened = true
      end if
      call get_xyz_for_bnd_face_point(nf,k,xyz_fp)
      write (io,1) "mypnum   = ",mypnum
      write (io,1) "glb_fp   = ",n
      write (io,1) "nf       = ",nf
      write (io,1) "loc_fp   = ",k
      write (io,2) "X        = ",xyz_fp(1)
      write (io,2) "dir_cos1 = ",dir_cosx
      write (io,2) "rafac    = ",rafac,      " "
      write (io,6) "Non-dimensional values"
      write (io,2) "dens_int  = ",rho_int
      write (io,2) "vx_int    = ",vx_int
      write (io,2) "pres_int  = ",pres_int
      write (io,2) "temp_int  = ",temp_int
     !write (io,2) "temp_int  = ",pres_int/rho_int/rgasref_nd
      write (io,2) "aspd_int  = ",aspd_int
      write (io,2) "veln_int  = ",veln_int
      write (io,2) "rminus    = ",rminus,     " "
      write (io,2) "conv_temp = ",conv_temp
      write (io,2) "func      = ",func
      write (io,2) "dfunc     = ",dfunc
      write (io,2) "dtemp     = ",dtemp
      write (io,2) "dterr     = ",dterr,      " "
      write (io,2) "dens_ext  = ",rho_ext
      write (io,2) "vx_ext    = ",vx_ext
      write (io,2) "pres_ext  = ",pres_ext
      write (io,2) "temp_ext  = ",temp_ext
      write (io,2) "aspd_ext  = ",aspd_ext
      write (io,2) "veln_ext  = ",veln_ext
      write (io,2) "mach_ext  = ",mach_ext,   " "
      write (io,6) "Dimensional values"
      write (io,7) "dens_int  = ",rho_int*rhoref,           " kg/m^3"
      write (io,7) "vx_int    = ",vx_int*aref,              " m/s"
      write (io,7) "pres_int  = ",pres_int*rhoa2,           " Pa"
      write (io,7) "temp_int  = ",temp_int*tref,            " K"
      write (io,7) "aspd_int  = ",aspd_int*aref,            " m/s"
      write (io,7) "veln_int  = ",veln_int*aref,            " m/s"
      write (io,7) "rminus    = ",rminus*aref,              " m/s",     " "
      write (io,7) "conv_temp = ",conv_temp*tref,           " K"
      write (io,7) "func      = ",func*tref,                " K"
      write (io,7) "dtemp     = ",dtemp*tref,               " K"
      write (io,7) "dterr     = ",dterr*tref,               " K",       " "
      write (io,7) "dens_ext  = ",rho_ext*rhoref,           " kg/m^3"
      write (io,7) "vx_ext    = ",vx_ext*aref,              " m/s"
      write (io,7) "pres_ext  = ",pres_ext*rhoa2,           " Pa"
      write (io,7) "temp_ext  = ",temp_ext*tref,            " K"
      write (io,7) "aspd_ext  = ",aspd_ext*aref,            " m/s"
      write (io,7) "veln_ext  = ",veln_ext*aref,            " m/s",     " "
      do k = 1,size(dtemp_save,dim=2)
        if (k == 1) write (io,3)
        write (io,4) k,(dtemp_save(i,k), i=1,size(dtemp_save,dim=1))
        if (k == size(dtemp_save,dim=2)) then
          write (io,4) k+1,dtemp,dterr,func,dfunc,temp_ext,aspd_ext,veln_ext,"-"
        end if
      end do
      flush (io)
#endif
    end if
    !
#ifndef DISABLE_FOR_GPU
    ! If this is the last face point on this face, output the final
    ! BC state (primitive variables) for all the face points on this
    ! boundary face
    ! NOTE: last_fp will only have a non-zero value if convergence_reached
    !       was false for one of the face points of the current face.
    !       This means the following will only be written if the above
    !       Newton iteration failed to converge for any of the face points.
    !
    if (n == last_fp) then
      write (io,8) nf
      do k = face(nf)%beg_fp,face(nf)%end_fp
        !
        i = k-face(nf)%beg_fp+1
        !
        call get_xyz_for_bnd_face_point(nf,i,xyz_fp)
        !
        rho_int = face_cv(k,nec,1)
        vx_int = face_cv(k,nmx,1)/rho_int
        pres_int = gm1*(face_cv(k,nee,1) - half*rho_int*vx_int*vx_int)
        !
        rho_ext = face_cv(k,nec,2)
        vx_ext = face_cv(k,nmx,2)/rho_ext
        pres_ext = gm1*(face_cv(k,nee,2) - half*rho_ext*vx_ext*vx_ext)
        !
        write (io,9) k,xyz_fp(1), &
                 "rho",  rho_int,  rho_ext,  rho_int*rhoref, rho_ext*rhoref, &
                  "vx",   vx_int,   vx_ext,   vx_int*aref,    vx_ext*aref, &
                "pres", pres_int, pres_ext, pres_int*rhoa2, pres_ext*rhoa2
        !
      end do
      write (io,8)
      flush (io)
    end if
#endif
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
  2 format (a,g15.8,:,/,a)
  3 format (/,101("-"),/, &
              4x,"VALUES AFFECTED BY EACH NEWTON ITERATION",/, &
              101("-"),/, &
              1x,"Iter",7x,"dtemp",11x,"dterr",11x,"func",12x,"dfunc", &
              9x,"temp_ext",8x,"aspd_ext",8x,"veln_ext",/, &
              101("-"))
  4 format (2x,i3,7(2x,es14.6),:,/,100("-"),a,/)
  5 format ("newton_bc_errors.CPU-",i0,".log")
  6 format (4x,a,/)
  7 format (a,g15.8,a,:,/,a)
  8 format (/,101("="),/,:, &
              4x,"Final states for the face points on face ",i0,/, &
              101("="),/, &
              10x,"Non-Dim Int State",3x,"Non-Dim Ext State", &
              12x,"Dim Int State",5x,"Dim Ext State",/, &
              101("="))
  9 format (/,5x,"Flux point # ",i0, &
              t49,"[x] = [",es13.6,"]",/, &
              1x,a,t11,g15.8,t31,g15.8,t58,g15.8,t76,g15.8,/, &
              1x,a,t11,g15.8,t31,g15.8,t58,g15.8,t76,g15.8,/, &
              1x,a,t11,g15.8,t31,g15.8,t58,g15.8,t76,g15.8)
  !
contains
  !
  pure function get_dterr_func(a,b) result(return_value)
    real(wp), intent(in) :: a
    real(wp), intent(in) :: b
    real(wp) :: return_value
  continue
#ifdef DTERR_USES_TEMP_EXT
    return_value = abs(a)/b
#else
    return_value = abs(a)
#endif
  end function get_dterr_func
  !
end subroutine get_bc_sub_inflow_newton_1d
!
!###############################################################################
!
pure &
subroutine get_xyz_for_bnd_face_point(nf,kf,xyz_fp)
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
  character(len=*), parameter :: pname = "get_xyz_for_bnd_face_point"
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
  do l = 1,size(xyz_fp)
    do k = 1,np
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
  do l = 1,size(xyz_fp)
    xyz_fp(l) = zero
  end do
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
end subroutine get_xyz_for_bnd_face_point
!
!###############################################################################
!
pure subroutine get_bc_sub_inflow_sol_cpg_v2_3d(n1,n2,face_cv,cv_in)
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
end subroutine get_bc_sub_inflow_sol_cpg_v2_3d
!
!###############################################################################
!
pure subroutine get_bc_sub_inflow_sol_cpg_v2_2d(n1,n2,face_cv,cv_in)
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
end subroutine get_bc_sub_inflow_sol_cpg_v2_2d
!
!###############################################################################
!
pure subroutine get_bc_sub_inflow_sol_cpg_v2_1d(n1,n2,face_cv,cv_in)
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
  real(wp) :: rho,rhou,rhoe
  !
continue
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
  !
end subroutine get_bc_sub_inflow_sol_cpg_v2_1d
!
!###############################################################################
!
pure &
subroutine get_bc_sub_inflow_sol_cpg_v3_3d(n1,n2,facenrm,face_cv,face_gp,pv_bc)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: a,b,c,d,rplus,htot,aspd_p,aspd_m
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
 !real(wp) :: aoa,tan_aoa
  real(wp) :: vmag,dir_cosx,dir_cosy,dir_cosz
  real(wp) :: total_cond,rtot_bc,ptot_bc
  real(wp) :: gm1,gam_over_1mg
  real(wp) :: one_over_gm1,one_over_1mg
  real(wp) :: rho,vx,vy,vz,q2,vel,pres,aspd,mach
  real(wp) :: rho_inv,aspd_int,mach_int,vel_int
  !
continue
  !
  ! Precompute some constants using gamma
  !
  gm1 = gam - one
  one_over_gm1 = one / gm1
  one_over_1mg = -one_over_gm1
  gam_over_1mg = gam * one_over_1mg
  !
  a = one + two*one_over_gm1
  !
  ! Get the boundary total density and total pressure
  ! NOTE: These total conditions have been saved into static
  !       density and static pressure locations of pv_bc.
  !
  rtot_bc = pv_bc(nec)
  ptot_bc = pv_bc(nee)
  !
  ! Get the directional cosines of the specified boundary condition velocity
  !
  vmag = sqrt( pv_bc(nmx)*pv_bc(nmx) + &
               pv_bc(nmy)*pv_bc(nmy) + &
               pv_bc(nmz)*pv_bc(nmz) )
  !
  dir_cosx = pv_bc(nmx) / vmag
  dir_cosy = pv_bc(nmy) / vmag
  dir_cosz = pv_bc(nmz) / vmag
 !!
 !aoa = atan2(pv_bc(nmy),pv_bc(nmx))
 !tan_aoa = tan( aoa )
 !!
 !dir_cosx = one / sqrt(one + tan_aoa**2)
 !dir_cosy = tan_aoa * dir_cosx
 !dir_cosz = zero
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp)
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
    ! Get the host cell primitive variables at this flux point
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    vz = face_cv(n,nmz,1)*rho_inv
    q2 = vx*vx + vy*vy + vz*vz
    !
    pres = face_gp(n,npres,1)
    !
    ! Magnitude of velocity of the interior/host cell at this flux point
    !
   !vel_int = sqrt( q2 )
    vel_int = vx*unx + vy*uny + vz*unz
    !
    ! Compute the speed of sound at this flux point
    !
    aspd_int = sqrt( max( gam*pres*rho_inv , rndoff ) )
    mach_int = vel_int / aspd_int
    !
    rplus = -vel_int - two*aspd_int*one_over_gm1
    htot  = aspd_int*aspd_int*one_over_gm1 + half*q2
    !
    b = two*rplus
    c = half*gm1*(rplus*rplus - two*htot)
    !
    d = sqrt( max( b*b-four*a*c , rndoff ) )
    !
    if (d >= zero) then
      aspd_p =  half * ( d - b ) / a
      aspd_m = -half * ( d + b ) / a
      aspd = max( aspd_p , aspd_m )
      if (aspd < zero) then
        aspd = aspd_int
      end if
    else
      aspd_p = -999.0_wp
      aspd_m = -999.0_wp
      aspd = aspd_int
    end if
    !
    vel = -rplus - two*aspd*one_over_gm1
    mach = abs(vel)/aspd
    !
    total_cond = one + half*gm1*mach*mach
    !
    ! Compute the density and pressure using the inflow total density,
    ! total pressure, and interior normal Mach number
    !
    rho  = rtot_bc * total_cond**one_over_1mg
    pres = ptot_bc * total_cond**gam_over_1mg
    !
    ! Set the boundary velocity in the same direction
    ! as the specified boundary condition velocity
    !
    vx = abs(vel)*dir_cosx
    vy = abs(vel)*dir_cosy
    vz = abs(vel)*dir_cosz
    !
    ! Set the boundary velocity normal to the boundary face
    !
   !vx = vel*unx
   !vy = vel*uny
   !vz = vel*unz
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rho*vx
    face_cv(n,nmy,2) = rho*vy
    face_cv(n,nmz,2) = rho*vz
    face_cv(n,nee,2) = pres*one_over_gm1 + half*rho*vel*vel
    !
  end do
  !
end subroutine get_bc_sub_inflow_sol_cpg_v3_3d
!
!###############################################################################
!
pure &
subroutine get_bc_sub_inflow_sol_cpg_v3_2d(n1,n2,facenrm,face_cv,face_gp,pv_bc)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: a,b,c,d,rplus,htot,aspd_p,aspd_m
  real(wp) :: ds,rnx,rny,unx,uny
 !real(wp) :: aoa,tan_aoa
  real(wp) :: vmag,dir_cosx,dir_cosy
  real(wp) :: total_cond,rtot_bc,ptot_bc
  real(wp) :: gm1,gam_over_1mg
  real(wp) :: one_over_gm1,one_over_1mg
  real(wp) :: rho,vx,vy,q2,vel,pres,aspd,mach
  real(wp) :: rho_inv,aspd_int,mach_int,vel_int
  !
continue
  !
  ! Precompute some constants using gamma
  !
  gm1 = gam - one
  one_over_gm1 = one / gm1
  one_over_1mg = -one_over_gm1
  gam_over_1mg = gam * one_over_1mg
  !
  a = one + two*one_over_gm1
  !
  ! Get the boundary total density and total pressure
  ! NOTE: These total conditions have been saved into static
  !       density and static pressure locations of pv_bc.
  !
  rtot_bc = pv_bc(nec)
  ptot_bc = pv_bc(nee)
  !
  ! Get the directional cosines of the specified boundary condition velocity
  !
  vmag = sqrt( pv_bc(nmx)*pv_bc(nmx) + &
               pv_bc(nmy)*pv_bc(nmy) )
  !
  dir_cosx = pv_bc(nmx) / vmag
  dir_cosy = pv_bc(nmy) / vmag
 !!
 !aoa = atan2(pv_bc(nmy),pv_bc(nmx))
 !tan_aoa = tan( aoa )
 !!
 !dir_cosx = one / sqrt(one + tan_aoa**2)
 !dir_cosy = tan_aoa * dir_cosx
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp)
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
    ! Get the host cell primitive variables at this flux point
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    q2 = vx*vx + vy*vy
    !
    pres = face_gp(n,npres,1)
    !
    ! Magnitude of velocity of the interior/host cell at this flux point
    !
   !vel_int = sqrt( q2 )
    vel_int = vx*unx + vy*uny
    !
    ! Compute the speed of sound at this flux point
    !
    aspd_int = sqrt( max( gam*pres*rho_inv , rndoff ) )
    mach_int = vel_int / aspd_int
    !
    rplus = -vel_int - two*aspd_int*one_over_gm1
    htot  = aspd_int*aspd_int*one_over_gm1 + half*q2
    !
    b = two*rplus
    c = half*gm1*(rplus*rplus - two*htot)
    !
    d = sqrt( max( b*b-four*a*c , rndoff ) )
    !
    if (d >= zero) then
      aspd_p =  half * ( d - b ) / a
      aspd_m = -half * ( d + b ) / a
      aspd = max( aspd_p , aspd_m )
      if (aspd < zero) then
        aspd = aspd_int
      end if
    else
      aspd_p = -999.0_wp
      aspd_m = -999.0_wp
      aspd = aspd_int
    end if
    !
    vel = -rplus - two*aspd*one_over_gm1
    mach = abs(vel)/aspd
    !
    total_cond = one + half*gm1*mach*mach
    !
    ! Compute the density and pressure using the inflow total density,
    ! total pressure, and interior normal Mach number
    !
    rho  = rtot_bc * total_cond**one_over_1mg
    pres = ptot_bc * total_cond**gam_over_1mg
    !
    ! Set the boundary velocity in the same direction
    ! as the specified boundary condition velocity
    !
    vx = abs(vel)*dir_cosx
    vy = abs(vel)*dir_cosy
    !
    ! Set the boundary velocity normal to the boundary face
    !
   !vx = vel*unx
   !vy = vel*uny
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rho*vx
    face_cv(n,nmy,2) = rho*vy
    face_cv(n,nee,2) = pres*one_over_gm1 + half*rho*vel*vel
    !
  end do
  !
end subroutine get_bc_sub_inflow_sol_cpg_v3_2d
!
!###############################################################################
!
pure &
subroutine get_bc_sub_inflow_sol_cpg_v3_1d(n1,n2,facenrm,face_cv,face_gp,pv_bc)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: a,b,c,d,rplus,htot,aspd_p,aspd_m
  real(wp) :: ds,rnx,unx
 !real(wp) :: aoa,tan_aoa
  real(wp) :: vmag,dir_cosx
  real(wp) :: total_cond,rtot_bc,ptot_bc
  real(wp) :: gm1,gam_over_1mg
  real(wp) :: one_over_gm1,one_over_1mg
  real(wp) :: rho,vx,q2,vel,pres,aspd,mach
  real(wp) :: rho_inv,aspd_int,mach_int,vel_int
  !
continue
  !
  ! Precompute some constants using gamma
  !
  gm1 = gam - one
  one_over_gm1 = one / gm1
  one_over_1mg = -one_over_gm1
  gam_over_1mg = gam * one_over_1mg
  !
  a = one + two*one_over_gm1
  !
  ! Get the boundary total density and total pressure
  ! NOTE: These total conditions have been saved into static
  !       density and static pressure locations of pv_bc.
  !
  rtot_bc = pv_bc(nec)
  ptot_bc = pv_bc(nee)
  !
  ! Get the directional cosines of the specified boundary condition velocity
  !
  vmag = sqrt( pv_bc(nmx)*pv_bc(nmx) )
  !
  dir_cosx = pv_bc(nmx) / vmag
 !!
 !aoa = atan2(pv_bc(nmy),pv_bc(nmx))
 !tan_aoa = tan( aoa )
 !!
 !dir_cosx = one / sqrt(one + tan_aoa**2)
 !dir_cosy = tan_aoa * dir_cosx
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp)
  do n = n1,n2
    !
    ! Unit normal at this flux point
    !
    rnx = facenrm(1,n)
    ds = sqrt( rnx*rnx )
    unx = rnx/ds
    !
    ! Get the host cell primitive variables at this flux point
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    q2 = vx*vx
    !
    pres = face_gp(n,npres,1)
    !
    ! Magnitude of velocity of the interior/host cell at this flux point
    !
   !vel_int = sqrt( q2 )
    vel_int = vx*unx
    !
    ! Compute the speed of sound at this flux point
    !
    aspd_int = sqrt( max( gam*pres*rho_inv , rndoff ) )
    mach_int = vel_int / aspd_int
    !
    rplus = -vel_int - two*aspd_int*one_over_gm1
    htot  = aspd_int*aspd_int*one_over_gm1 + half*q2
    !
    b = two*rplus
    c = half*gm1*(rplus*rplus - two*htot)
    !
    d = sqrt( max( b*b-four*a*c , rndoff ) )
    !
    if (d >= zero) then
      aspd_p =  half * ( d - b ) / a
      aspd_m = -half * ( d + b ) / a
      aspd = max( aspd_p , aspd_m )
      if (aspd < zero) then
        aspd = aspd_int
      end if
    else
      aspd_p = -999.0_wp
      aspd_m = -999.0_wp
      aspd = aspd_int
    end if
    !
    vel = -rplus - two*aspd*one_over_gm1
    mach = abs(vel)/aspd
    !
    total_cond = one + half*gm1*mach*mach
    !
    ! Compute the density and pressure using the inflow total density,
    ! total pressure, and interior normal Mach number
    !
    rho  = rtot_bc * total_cond**one_over_1mg
    pres = ptot_bc * total_cond**gam_over_1mg
    !
    ! Set the boundary velocity in the same direction
    ! as the specified boundary condition velocity
    !
    vx = abs(vel)*dir_cosx
    !
    ! Set the boundary velocity normal to the boundary face
    !
   !vx = vel*unx
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rho*vx
    face_cv(n,nee,2) = pres*one_over_gm1 + half*rho*vel*vel
    !
  end do
  !
end subroutine get_bc_sub_inflow_sol_cpg_v3_1d
!
!###############################################################################
!
pure &
subroutine get_bc_sub_inflow_sol_cgp_v4_3d(n1,n2,facenrm,face_cv,face_gp,pv_bc)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
 !real(wp) :: aoa,tan_aoa
  real(wp) :: vmag,dir_cosx,dir_cosy,dir_cosz
  real(wp) :: total_cond,rtot_bc,ptot_bc
  real(wp) :: rho,vx,vy,vz,vel,pres
  real(wp) :: rho_inv,aspd_int,mach_int,veln_int
  real(wp) :: gm1,one_over_gm1,one_over_1mg,gam_over_1mg
  !
continue
  !
  gm1 = gam - one
  one_over_gm1 = one / gm1
  one_over_1mg = -one_over_gm1
  gam_over_1mg = gam * one_over_1mg
  !
  ! Get the boundary total density and total pressure
  ! NOTE: These total conditions have been saved into static
  !       density and static pressure locations of pv_bc.
  !
  rtot_bc = pv_bc(nec)
  ptot_bc = pv_bc(nee)
  !
  ! Get the directional cosines of the specified boundary condition velocity
  !
  vmag = sqrt( pv_bc(nmx)*pv_bc(nmx) + &
               pv_bc(nmy)*pv_bc(nmy) + &
               pv_bc(nmz)*pv_bc(nmz) )
  !
  dir_cosx = pv_bc(nmx) / vmag
  dir_cosy = pv_bc(nmy) / vmag
  dir_cosz = pv_bc(nmz) / vmag
 !!
 !aoa = atan2(pv_bc(nmy),pv_bc(nmx))
 !tan_aoa = tan( aoa )
 !!
 !dir_cosx = one / sqrt(one + tan_aoa**2)
 !dir_cosy = tan_aoa * dir_cosx
 !dir_cosz = zero
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp)
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
    ! Get the host cell primitive variables at this flux point
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    vz = face_cv(n,nmz,1)*rho_inv
    !
    pres = face_gp(n,npres,1)
    !
    ! Normal velocity of the interior/host cell at this flux point
    !
    veln_int = vx*unx + vy*uny + vz*unz
    !
    ! Compute the speed of sound at this flux point
    !
    aspd_int = sqrt( max( gam*pres*rho_inv , rndoff ) )
    mach_int = veln_int / aspd_int
    !
    total_cond = one + half*gm1*mach_int*mach_int
    !
    ! Compute the density and pressure using the inflow total density,
    ! total pressure, and interior normal Mach number
    !
    rho  = rtot_bc * total_cond**one_over_1mg
    pres = ptot_bc * total_cond**gam_over_1mg
    !
    vel = mach_int * sqrt( max( gam*pres/rho , rndoff ) )
    !
    ! Set the boundary velocity in the same direction
    ! as the specified boundary condition velocity
    !
    vx = abs(vel)*dir_cosx
    vy = abs(vel)*dir_cosy
    vz = abs(vel)*dir_cosz
    !
    ! Set the boundary velocity normal to the boundary face
    !
   !vx = vel*unx
   !vy = vel*uny
   !vz = vel*unz
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rho*vx
    face_cv(n,nmy,2) = rho*vy
    face_cv(n,nmz,2) = rho*vz
    face_cv(n,nee,2) = pres*one_over_gm1 + half*rho*vel*vel
    !
  end do
  !
end subroutine get_bc_sub_inflow_sol_cgp_v4_3d
!
!###############################################################################
!
pure &
subroutine get_bc_sub_inflow_sol_cgp_v4_2d(n1,n2,facenrm,face_cv,face_gp,pv_bc)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,rny,unx,uny
 !real(wp) :: aoa,tan_aoa
  real(wp) :: vmag,dir_cosx,dir_cosy
  real(wp) :: total_cond,rtot_bc,ptot_bc
  real(wp) :: rho,vx,vy,vel,pres
  real(wp) :: rho_inv,aspd_int,mach_int,veln_int
  real(wp) :: gm1,one_over_gm1,one_over_1mg,gam_over_1mg
  !
continue
  !
  gm1 = gam - one
  one_over_gm1 = one / gm1
  one_over_1mg = -one_over_gm1
  gam_over_1mg = gam * one_over_1mg
  !
  ! Get the boundary total density and total pressure
  ! NOTE: These total conditions have been saved into static
  !       density and static pressure locations of pv_bc.
  !
  rtot_bc = pv_bc(nec)
  ptot_bc = pv_bc(nee)
  !
  ! Get the directional cosines of the specified boundary condition velocity
  !
  vmag = sqrt( pv_bc(nmx)*pv_bc(nmx) + &
               pv_bc(nmy)*pv_bc(nmy) )
  !
  dir_cosx = pv_bc(nmx) / vmag
  dir_cosy = pv_bc(nmy) / vmag
 !!
 !aoa = atan2(pv_bc(nmy),pv_bc(nmx))
 !tan_aoa = tan( aoa )
 !!
 !dir_cosx = one / sqrt(one + tan_aoa**2)
 !dir_cosy = tan_aoa * dir_cosx
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp)
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
    ! Get the host cell primitive variables at this flux point
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    !
    pres = face_gp(n,npres,1)
    !
    ! Normal velocity of the interior/host cell at this flux point
    !
    veln_int = vx*unx + vy*uny
    !
    ! Compute the speed of sound at this flux point
    !
    aspd_int = sqrt( max( gam*pres*rho_inv , rndoff ) )
    mach_int = veln_int / aspd_int
    !
    total_cond = one + half*gm1*mach_int*mach_int
    !
    ! Compute the density and pressure using the inflow total density,
    ! total pressure, and interior normal Mach number
    !
    rho  = rtot_bc * total_cond**one_over_1mg
    pres = ptot_bc * total_cond**gam_over_1mg
    !
    vel = mach_int * sqrt( max( gam*pres/rho , rndoff ) )
    !
    ! Set the boundary velocity in the same direction
    ! as the specified boundary condition velocity
    !
    vx = abs(vel)*dir_cosx
    vy = abs(vel)*dir_cosy
    !
    ! Set the boundary velocity normal to the boundary face
    !
   !vx = vel*unx
   !vy = vel*uny
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rho*vx
    face_cv(n,nmy,2) = rho*vy
    face_cv(n,nee,2) = pres*one_over_gm1 + half*rho*vel*vel
    !
  end do
  !
end subroutine get_bc_sub_inflow_sol_cgp_v4_2d
!
!###############################################################################
!
pure &
subroutine get_bc_sub_inflow_sol_cgp_v4_1d(n1,n2,facenrm,face_cv,face_gp,pv_bc)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: pv_bc(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,unx
 !real(wp) :: aoa,tan_aoa
  real(wp) :: vmag,dir_cosx
  real(wp) :: total_cond,rtot_bc,ptot_bc
  real(wp) :: rho,vx,vel,pres
  real(wp) :: rho_inv,aspd_int,mach_int,veln_int
  real(wp) :: gm1,one_over_gm1,one_over_1mg,gam_over_1mg
  !
continue
  !
  gm1 = gam - one
  one_over_gm1 = one / gm1
  one_over_1mg = -one_over_gm1
  gam_over_1mg = gam * one_over_1mg
  !
  ! Get the boundary total density and total pressure
  ! NOTE: These total conditions have been saved into static
  !       density and static pressure locations of pv_bc.
  !
  rtot_bc = pv_bc(nec)
  ptot_bc = pv_bc(nee)
  !
  ! Get the directional cosines of the specified boundary condition velocity
  !
  vmag = sqrt( pv_bc(nmx)*pv_bc(nmx) )
  !
  dir_cosx = pv_bc(nmx) / vmag
 !!
 !aoa = atan2(pv_bc(nmy),pv_bc(nmx))
 !tan_aoa = tan( aoa )
 !!
 !dir_cosx = one / sqrt(one + tan_aoa**2)
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp)
  do n = n1,n2
    !
    ! Unit normal at this flux point
    !
    rnx = facenrm(1,n)
    ds = sqrt( rnx*rnx )
    unx = rnx/ds
    !
    ! Get the host cell primitive variables at this flux point
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    vx = face_cv(n,nmx,1)*rho_inv
    !
    pres = face_gp(n,npres,1)
    !
    ! Normal velocity of the interior/host cell at this flux point
    !
    veln_int = vx*unx
    !
    ! Compute the speed of sound at this flux point
    !
    aspd_int = sqrt( max( gam*pres*rho_inv , rndoff ) )
    mach_int = veln_int / aspd_int
    !
    total_cond = one + half*gm1*mach_int*mach_int
    !
    ! Compute the density and pressure using the inflow total density,
    ! total pressure, and interior normal Mach number
    !
    rho  = rtot_bc * total_cond**one_over_1mg
    pres = ptot_bc * total_cond**gam_over_1mg
    !
    vel = mach_int * sqrt( max( gam*pres/rho , rndoff ) )
    !
    ! Set the boundary velocity in the same direction
    ! as the specified boundary condition velocity
    !
    vx = abs(vel)*dir_cosx
    !
    ! Set the boundary velocity normal to the boundary face
    !
   !vx = vel*unx
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rho*vx
    face_cv(n,nee,2) = pres*one_over_gm1 + half*rho*vel*vel
    !
  end do
  !
end subroutine get_bc_sub_inflow_sol_cgp_v4_1d
!
!###############################################################################
!
pure &
subroutine get_bc_sub_inflow_sol_tpg_3d(n1,n2,facenrm,face_cv,face_gp, &
                                              pv_bc,newton_failures)
  !
  ! Routine to set subsonic inflow boundary conditions on any subset of a
  ! block face
  !
  ! face_cv     : conserved variable array to be set
  ! face_gp    : gas variable array to be set
  ! facenrm  : metrics in the b.c. coordinate system
  ! pv_bc   : array containing primitive variables for this boundary condition
  !           NOTE: For this subroutine, the pv_bc array should contain:
  !                 pv_bc(1:ncs)   = relaxed species mass fractions
  !                 pv_bc(nec)     = relaxed total density
  !                 pv_bc(nmb:nme) = specified inflow velocities
  !                 pv_bc(nee)     = relaxed total pressure
  !
  !.. Use Statements ..
  use ovar, only : tref,tmin,tmax
  !
  use thermo_mod, only : compute_mixture_gamma
  use thermo_mod, only : compute_mixture_rgas
  use thermo_mod, only : compute_mixture_enthalpy
  use thermo_mod, only : compute_mixture_cp
  use thermo_mod, only : compute_cp_over_temp_integral
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: pv_bc(:)
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp), intent(inout) :: face_gp(:,:,:)
  integer,  intent(inout) :: newton_failures
  !
  !.. Local Scalars ..
  integer  :: k,m,n,ndnc
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: vmag,dir_cosx,dir_cosy,dir_cosz
  real(wp) :: rho,rho_inv,vx,vy,vz,vn,vel,pres
  real(wp) :: gamm,htot,mach,temp,enth,cp,rmix
  real(wp) :: rmix_bc,rtot_bc,ptot_bc,ttot_bc,htot_bc
  real(wp) :: dtemp,dtmx,dterr,func,dfunc,tmin_nd,tmax_nd
  !
  logical(lk) :: convergence_reached
  !
  !.. Local Parameters ..
  integer,  parameter :: newton_iter = 100
  real(wp), parameter :: conv_temp = eps6
  !
continue
  !
  tmin_nd = tmin / tref
  tmax_nd = tmax / tref
  !
  newton_failures = 0
  !
  ! Get the inflow total pressure and total density
  !
  rtot_bc = pv_bc(nec)
  ptot_bc = pv_bc(nee)
  !
  ! Get the directional cosines of the specified boundary condition velocity
  !
  vmag = sqrt( pv_bc(nmx)*pv_bc(nmx) + &
               pv_bc(nmy)*pv_bc(nmy) + &
               pv_bc(nmz)*pv_bc(nmz) )
  !
  dir_cosx = pv_bc(nmx) / vmag
  dir_cosy = pv_bc(nmy) / vmag
  dir_cosz = pv_bc(nmz) / vmag
  !
  ! Stagnation pressure and temperature specified, Mach number extrapolated
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp)
  do n = n1,n2
    !
    ! Face normal at this boundary point
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
    !
    ! Extract velocity from the interior face point
    !
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    vz = face_cv(n,nmz,1)*rho_inv
    !
    vn = vx*unx + vy*uny + vz*unz
    !
    ! Extract gas properties from the interior face point
    !
    pres = face_gp(n,npres,1)
    temp = face_gp(n,ntemp,1)
    gamm = face_gp(n,ngamm,1)
    !
    mach = abs(vn) / sqrt(max(gamm*pres*rho_inv,rndoff))
    !
    ! Use the total energy array for temporary storage of Mach number
    !
    face_cv(n,nee,2) = mach
    face_gp(n,ntemp,2) = temp
    !
  end do
  !
  ! Compute the mixture gas constant at the inflow plane
  !
  call compute_mixture_rgas(pv_bc(1:ncs),rmix_bc)
  !
  ! Compute the total temperature (dimensional) at the inflow plane
  !
  ttot_bc = ptot_bc / (rtot_bc*rmix_bc)
  !
  ! Compute the total enthalpy at the inflow plane
  !
  call compute_mixture_enthalpy(pv_bc(1:ncs),ttot_bc,htot_bc)
  !
  ! Copy the specified species mass fractions into the boundary state
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv,pv_bc)
    do m = 1,ncs
      do n = n1,n2
        face_cv(n,m,2) = pv_bc(m)
      end do
    end do
  end if
  !
  ! Compute static temperature using a Newton-Raphson method
  !
  convergence_reached = fals
  newton_iteration: &
  do k = 1,newton_iter
    !
    call compute_mixture_enthalpy(n1,n2,face_cv(:,1:ncs,2), &
                                        face_gp(:,ntemp,2), &
                                        face_gp(:,nenth,2))
    !
    call compute_mixture_cp(n1,n2,face_cv(:,1:ncs,2), &
                                  face_gp(:,ntemp,2), &
                                  face_gp(:,ncp,2))
    !
    !$ACC PARALLEL LOOP PRESENT(face_gp)
    do n = n1,n2
      face_gp(n,ngamm,2) = rmix_bc
    end do
    !
    call compute_mixture_gamma(n1,n2,face_gp(:,ncp,2),face_gp(:,ngamm,2))
    !
    ndnc = 0
    dtmx = zero
    !
    !$ACC PARALLEL LOOP PRESENT(face_cv,face_gp) &
    !$ACC&              REDUCTION(+:ndnc) &
    !$ACC&              REDUCTION(max:dtmx)
    do n = n1,n2
      !
      mach = face_cv(n,nee,2)
      !
      temp = face_gp(n,ntemp,2)
      gamm = face_gp(n,ngamm,2)
      enth = face_gp(n,nenth,2)
      cp   = face_gp(n,ncp,2)
      !
      htot = enth + half * gamm*rmix_bc*temp * mach*mach
      !
      func = htot - htot_bc
      dfunc = cp + half*gamm*rmix_bc * mach*mach
      !
      dtemp = func/dfunc
     !func = (htot - htot_bc) / (cp + half*gamm*rmix_bc*mach**2)
      !
      face_gp(n,ntemp,2) = min( tmax_nd , max(tmin_nd,temp-dtemp) )
      !
      dterr = abs(dtemp)/face_gp(n,ntemp,2)
      !
      if (dterr > dtmx)  dtmx = dterr
      if (dterr > conv_temp) ndnc = ndnc + 1
      !
    end do
    !
    if (dtmx <= conv_temp) then
      convergence_reached = true
      exit newton_iteration
    end if
    !
  end do newton_iteration
  !
  if (.not. convergence_reached) then
    newton_failures = ndnc
  end if
  !
  ! Temporarily store the boundary condition total temperature in
  ! the enthalpy location of the boundary state gas properties
  !
  !$ACC PARALLEL LOOP PRESENT(face_gp)
  do n = n1,n2
    face_gp(n,nenth,2) = ttot_bc
  end do
  !
  ! Compute the integral of Cp/T for each species and sum
  ! NOTE: face_gp(nenth) contains nondimensional BC total temperature
  !       face_gp(ntemp) contains nondimensional static temperature
  !       face_gp(npres) contains nondimensional integral of Cp/T
  !
  call compute_cp_over_temp_integral(n1,n2,face_cv(:,1:ncs,2), &
                                           face_gp(:,npres,2), &
                                           face_gp(:,ntemp,2), &
                                           face_gp(:,nenth,2))
  !
  ! Compute static pressure
  !
  !$ACC PARALLEL LOOP PRESENT(face_gp)
  do n = n1,n2
    face_gp(n,npres,2) = ptot_bc / exp(face_gp(n,npres,2)/rmix_bc)
  end do
  !
  ! Compute thermodynamic properties given pressure, temperature, and
  ! mass fractions
  !
  ! Compute gas constant of the mixture
  ! NOTE: Rmix is stored in face_gp(:,ngamm,2) temporarily
  !
  call compute_mixture_rgas(n1,n2,face_cv(:,1:ncs,2),face_gp(:,ngamm,2))
  !
  ! Use the ideal gas law to compute the updated density now that we
  ! have the updated mixture gas constant, temperature, and pressure.
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,face_gp)
  do n = n1,n2
    rmix = face_gp(n,ngamm,2)
    temp = face_gp(n,ntemp,2)
    pres = face_gp(n,npres,2)
    face_cv(n,nec,2) = pres / (rmix*temp)
  end do
  !
  ! Compute the mixture static enthalpy
  !
  call compute_mixture_enthalpy(n1,n2,face_cv(:,1:ncs,2), &
                                      face_gp(:,ntemp,2), &
                                      face_gp(:,nenth,2))
  !
  ! Compute the mixture specific heat at constant pressure, Cp
  !
  call compute_mixture_cp(n1,n2,face_cv(:,1:ncs,2), &
                                face_gp(:,ntemp,2), &
                                face_gp(:,ncp,2))
  !
  ! Compute the ratio of the specific heats for the gas mixture
  ! NOTE: Upon entry to compute_mixture_gamma:
  !           face_gp(:,ngamm,2) = Rgas_mix (the mixture gas constant)
  !       After returning from compute_mixture_gamma:
  !           face_gp(:,ngamm,2) = gamma_mix (the mixture gamma)
  !
  call compute_mixture_gamma(n1,n2,face_gp(:,ncp,2),face_gp(:,ngamm,2))
  !
  ! Assign boundary state density, momentum, and provisional total energy
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,face_gp)
  do n = n1,n2
    !
    ! Extract interior mach number which was temporarily
    ! stored in the total energy location of face_cv
    !
    mach = face_cv(n,nee,2)
    !
    ! Get the other flow conditions needed to compute the
    ! conserved variables at the boundary
    !
    rho = face_cv(n,nec,2)
    !
    pres = face_gp(n,npres,2)
    temp = face_gp(n,ntemp,2)
    gamm = face_gp(n,ngamm,2)
    enth = face_gp(n,nenth,2)
    !
    vel = mach * sqrt(max(gamm*pres/rho,rndoff))
    !
    ! Compute the conserved variables for this boundary point
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rho*vel*dir_cosx
    face_cv(n,nmy,2) = rho*vel*dir_cosy
    face_cv(n,nmz,2) = rho*vel*dir_cosz
   !face_cv(n,nee,2) = rho*enth - pres + half*rho*vel*vel
    face_cv(n,nee,2) = rho*htot_bc - rho*rmix_bc*face_gp(n,ntemp,2)
    !
  end do
  !
  ! Convert mass fractions to densities
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv)
    do m = 1,ncs
      do n = n1,n2
        face_cv(n,m,2) = face_cv(n,m,2) * face_cv(n,nec,2)
      end do
    end do
  end if
  !
  ! Finished setting subsonic inflow boundary conditions
  !
end subroutine get_bc_sub_inflow_sol_tpg_3d
!
!###############################################################################
!
pure &
subroutine get_bc_sub_inflow_sol_tpg_2d(n1,n2,facenrm,face_cv,face_gp, &
                                              pv_bc,newton_failures)
  !
  ! Routine to set subsonic inflow boundary conditions on any subset of a
  ! block face
  !
  ! face_cv     : conserved variable array to be set
  ! face_gp    : gas variable array to be set
  ! facenrm  : metrics in the b.c. coordinate system
  ! pv_bc   : array containing primitive variables for this boundary condition
  !           NOTE: For this subroutine, the pv_bc array should contain:
  !                 pv_bc(1:ncs)   = relaxed species mass fractions
  !                 pv_bc(nec)     = relaxed total density
  !                 pv_bc(nmb:nme) = specified inflow velocities
  !                 pv_bc(nee)     = relaxed total pressure
  !
  !.. Use Statements ..
  use ovar, only : tref,tmin,tmax
  !
  use thermo_mod, only : compute_mixture_gamma
  use thermo_mod, only : compute_mixture_rgas
  use thermo_mod, only : compute_mixture_enthalpy
  use thermo_mod, only : compute_mixture_cp
  use thermo_mod, only : compute_cp_over_temp_integral
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: pv_bc(:)
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp), intent(inout) :: face_gp(:,:,:)
  integer,  intent(inout) :: newton_failures
  !
  !.. Local Scalars ..
  integer  :: k,m,n,ndnc
  real(wp) :: ds,rnx,rny,unx,uny
  real(wp) :: vmag,dir_cosx,dir_cosy
  real(wp) :: rho,rho_inv,vx,vy,vn,vel,pres
  real(wp) :: gamm,htot,mach,temp,enth,cp,rmix
  real(wp) :: rmix_bc,rtot_bc,ptot_bc,ttot_bc,htot_bc
  real(wp) :: dtemp,dtmx,dterr,func,dfunc,tmin_nd,tmax_nd
  !
  logical(lk) :: convergence_reached
  !
  !.. Local Parameters ..
  integer,  parameter :: newton_iter = 100
  real(wp), parameter :: conv_temp = eps6
  !
continue
  !
  tmin_nd = tmin / tref
  tmax_nd = tmax / tref
  !
  newton_failures = 0
  !
  ! Get the inflow total pressure and total density
  !
  rtot_bc = pv_bc(nec)
  ptot_bc = pv_bc(nee)
  !
  ! Get the directional cosines of the specified boundary condition velocity
  !
  vmag = sqrt( pv_bc(nmx)*pv_bc(nmx) + &
               pv_bc(nmy)*pv_bc(nmy) )
  !
  dir_cosx = pv_bc(nmx) / vmag
  dir_cosy = pv_bc(nmy) / vmag
  !
  ! Stagnation pressure and temperature specified, Mach number extrapolated
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp)
  do n = n1,n2
    !
    ! Face normal at this boundary point
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds = sqrt( rnx*rnx + rny*rny )
    unx = rnx/ds
    uny = rny/ds
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    !
    ! Extract velocity from the interior face point
    !
    vx = face_cv(n,nmx,1)*rho_inv
    vy = face_cv(n,nmy,1)*rho_inv
    !
    vn = vx*unx + vy*uny
    !
    ! Extract gas properties from the interior face point
    !
    pres = face_gp(n,npres,1)
    temp = face_gp(n,ntemp,1)
    gamm = face_gp(n,ngamm,1)
    !
    mach = abs(vn) / sqrt(max(gamm*pres*rho_inv,rndoff))
    !
    ! Use the total energy array for temporary storage of Mach number
    !
    face_cv(n,nee,2) = mach
    face_gp(n,ntemp,2) = temp
    !
  end do
  !
  ! Compute the mixture gas constant at the inflow plane
  !
  call compute_mixture_rgas(pv_bc(1:ncs),rmix_bc)
  !
  ! Compute the total temperature (dimensional) at the inflow plane
  !
  ttot_bc = ptot_bc / (rtot_bc*rmix_bc)
  !
  ! Compute the total enthalpy at the inflow plane
  !
  call compute_mixture_enthalpy(pv_bc(1:ncs),ttot_bc,htot_bc)
  !
  ! Copy the specified species mass fractions into the boundary state
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv,pv_bc)
    do m = 1,ncs
      do n = n1,n2
        face_cv(n,m,2) = pv_bc(m)
      end do
    end do
  end if
  !
  ! Compute static temperature using a Newton-Raphson method
  !
  convergence_reached = fals
  newton_iteration: &
  do k = 1,newton_iter
    !
    call compute_mixture_enthalpy(n1,n2,face_cv(:,1:ncs,2), &
                                        face_gp(:,ntemp,2), &
                                        face_gp(:,nenth,2))
    !
    call compute_mixture_cp(n1,n2,face_cv(:,1:ncs,2), &
                                  face_gp(:,ntemp,2), &
                                  face_gp(:,ncp,2))
    !
    !$ACC PARALLEL LOOP PRESENT(face_gp)
    do n = n1,n2
      face_gp(n,ngamm,2) = rmix_bc
    end do
    !
    call compute_mixture_gamma(n1,n2,face_gp(:,ncp,2),face_gp(:,ngamm,2))
    !
    ndnc = 0
    dtmx = zero
    !
    !$ACC PARALLEL LOOP PRESENT(face_cv,face_gp) &
    !$ACC&              REDUCTION(+:ndnc) &
    !$ACC&              REDUCTION(max:dtmx)
    do n = n1,n2
      !
      mach = face_cv(n,nee,2)
      !
      temp = face_gp(n,ntemp,2)
      gamm = face_gp(n,ngamm,2)
      enth = face_gp(n,nenth,2)
      cp   = face_gp(n,ncp,2)
      !
      htot = enth + half * gamm*rmix_bc*temp * mach*mach
      !
      func = htot - htot_bc
      dfunc = cp + half*gamm*rmix_bc * mach*mach
      !
      dtemp = func/dfunc
     !func = (htot - htot_bc) / (cp + half*gamm*rmix_bc*mach**2)
      !
      face_gp(n,ntemp,2) = min( tmax_nd , max(tmin_nd,temp-dtemp) )
      !
      dterr = abs(dtemp)/face_gp(n,ntemp,2)
      !
      if (dterr > dtmx)  dtmx = dterr
      if (dterr > conv_temp) ndnc = ndnc + 1
      !
    end do
    !
    if (dtmx <= conv_temp) then
      convergence_reached = true
      exit newton_iteration
    end if
    !
  end do newton_iteration
  !
  if (.not. convergence_reached) then
    newton_failures = ndnc
  end if
  !
  ! Temporarily store the boundary condition total temperature in
  ! the enthalpy location of the boundary state gas properties
  !
  !$ACC PARALLEL LOOP PRESENT(face_gp)
  do n = n1,n2
    face_gp(n,nenth,2) = ttot_bc
  end do
  !
  ! Compute the integral of Cp/T for each species and sum
  ! NOTE: face_gp(nenth) contains nondimensional BC total temperature
  !       face_gp(ntemp) contains nondimensional static temperature
  !       face_gp(npres) contains nondimensional integral of Cp/T
  !
  call compute_cp_over_temp_integral(n1,n2,face_cv(:,1:ncs,2), &
                                           face_gp(:,npres,2), &
                                           face_gp(:,ntemp,2), &
                                           face_gp(:,nenth,2))
  !
  ! Compute static pressure
  !
  !$ACC PARALLEL LOOP PRESENT(face_gp)
  do n = n1,n2
    face_gp(n,npres,2) = ptot_bc / exp(face_gp(n,npres,2)/rmix_bc)
  end do
  !
  ! Compute thermodynamic properties given pressure, temperature, and
  ! mass fractions
  !
  ! Compute gas constant of the mixture
  ! NOTE: Rmix is stored in face_gp(:,ngamm,2) temporarily
  !
  call compute_mixture_rgas(n1,n2,face_cv(:,1:ncs,2),face_gp(:,ngamm,2))
  !
  ! Use the ideal gas law to compute the updated density now that we
  ! have the updated mixture gas constant, temperature, and pressure.
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,face_gp)
  do n = n1,n2
    rmix = face_gp(n,ngamm,2)
    temp = face_gp(n,ntemp,2)
    pres = face_gp(n,npres,2)
    face_cv(n,nec,2) = pres / (rmix*temp)
  end do
  !
  ! Compute the mixture static enthalpy
  !
  call compute_mixture_enthalpy(n1,n2,face_cv(:,1:ncs,2), &
                                      face_gp(:,ntemp,2), &
                                      face_gp(:,nenth,2))
  !
  ! Compute the mixture specific heat at constant pressure, Cp
  !
  call compute_mixture_cp(n1,n2,face_cv(:,1:ncs,2), &
                                face_gp(:,ntemp,2), &
                                face_gp(:,ncp,2))
  !
  ! Compute the ratio of the specific heats for the gas mixture
  ! NOTE: Upon entry to compute_mixture_gamma:
  !           face_gp(:,ngamm,2) = Rgas_mix (the mixture gas constant)
  !       After returning from compute_mixture_gamma:
  !           face_gp(:,ngamm,2) = gamma_mix (the mixture gamma)
  !
  call compute_mixture_gamma(n1,n2,face_gp(:,ncp,2),face_gp(:,ngamm,2))
  !
  ! Assign boundary state density, momentum, and provisional total energy
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,face_gp)
  do n = n1,n2
    !
    ! Extract interior mach number which was temporarily
    ! stored in the total energy location of face_cv
    !
    mach = face_cv(n,nee,2)
    !
    ! Get the other flow conditions needed to compute the
    ! conserved variables at the boundary
    !
    rho = face_cv(n,nec,2)
    !
    pres = face_gp(n,npres,2)
    temp = face_gp(n,ntemp,2)
    gamm = face_gp(n,ngamm,2)
    enth = face_gp(n,nenth,2)
    !
    vel = mach * sqrt(max(gamm*pres/rho,rndoff))
    !
    ! Compute the conserved variables for this boundary point
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rho*vel*dir_cosx
    face_cv(n,nmy,2) = rho*vel*dir_cosy
   !face_cv(n,nee,2) = rho*enth - pres + half*rho*vel*vel
    face_cv(n,nee,2) = rho*htot_bc - rho*rmix_bc*face_gp(n,ntemp,2)
    !
  end do
  !
  ! Convert mass fractions to densities
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv)
    do m = 1,ncs
      do n = n1,n2
        face_cv(n,m,2) = face_cv(n,m,2) * face_cv(n,nec,2)
      end do
    end do
  end if
  !
  ! Finished setting subsonic inflow boundary conditions
  !
end subroutine get_bc_sub_inflow_sol_tpg_2d
!
!###############################################################################
!
pure &
subroutine get_bc_sub_inflow_sol_tpg_1d(n1,n2,facenrm,face_cv,face_gp, &
                                              pv_bc,newton_failures)
  !
  ! Routine to set subsonic inflow boundary conditions on any subset of a
  ! block face
  !
  ! face_cv     : conserved variable array to be set
  ! face_gp    : gas variable array to be set
  ! facenrm  : metrics in the b.c. coordinate system
  ! pv_bc   : array containing primitive variables for this boundary condition
  !           NOTE: For this subroutine, the pv_bc array should contain:
  !                 pv_bc(1:ncs)   = relaxed species mass fractions
  !                 pv_bc(nec)     = relaxed total density
  !                 pv_bc(nmb:nme) = specified inflow velocities
  !                 pv_bc(nee)     = relaxed total pressure
  !
  !.. Use Statements ..
  use ovar, only : tref,tmin,tmax
  !
  use thermo_mod, only : compute_mixture_gamma
  use thermo_mod, only : compute_mixture_rgas
  use thermo_mod, only : compute_mixture_enthalpy
  use thermo_mod, only : compute_mixture_cp
  use thermo_mod, only : compute_cp_over_temp_integral
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: pv_bc(:)
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp), intent(inout) :: face_gp(:,:,:)
  integer,  intent(inout) :: newton_failures
  !
  !.. Local Scalars ..
  integer  :: k,m,n,ndnc
  real(wp) :: ds,rnx,unx
  real(wp) :: vmag,dir_cosx
  real(wp) :: rho,rho_inv,vx,vn,vel,pres
  real(wp) :: gamm,htot,mach,temp,enth,cp,rmix
  real(wp) :: rmix_bc,rtot_bc,ptot_bc,ttot_bc,htot_bc
  real(wp) :: dtemp,dtmx,dterr,func,dfunc,tmin_nd,tmax_nd
  !
  logical(lk) :: convergence_reached
  !
  !.. Local Parameters ..
  integer,  parameter :: newton_iter = 100
  real(wp), parameter :: conv_temp = eps6
  !
continue
  !
  tmin_nd = tmin / tref
  tmax_nd = tmax / tref
  !
  newton_failures = 0
  !
  ! Get the inflow total pressure and total density
  !
  rtot_bc = pv_bc(nec)
  ptot_bc = pv_bc(nee)
  !
  ! Get the directional cosines of the specified boundary condition velocity
  !
  vmag = sqrt( pv_bc(nmx)*pv_bc(nmx) )
  !
  dir_cosx = pv_bc(nmx) / vmag
  !
  ! Stagnation pressure and temperature specified, Mach number extrapolated
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp)
  do n = n1,n2
    !
    ! Face normal at this boundary point
    !
    rnx = facenrm(1,n)
    ds = sqrt( rnx*rnx )
    unx = rnx/ds
    !
    rho = face_cv(n,nec,1)
    rho_inv = one/rho
    !
    ! Extract velocity from the interior face point
    !
    vx = face_cv(n,nmx,1)*rho_inv
    !
    vn = vx*unx
    !
    ! Extract gas properties from the interior face point
    !
    pres = face_gp(n,npres,1)
    temp = face_gp(n,ntemp,1)
    gamm = face_gp(n,ngamm,1)
    !
    mach = abs(vn) / sqrt(max(gamm*pres*rho_inv,rndoff))
    !
    ! Use the total energy array for temporary storage of Mach number
    !
    face_cv(n,nee,2) = mach
    face_gp(n,ntemp,2) = temp
    !
  end do
  !
  ! Compute the mixture gas constant at the inflow plane
  !
  call compute_mixture_rgas(pv_bc(1:ncs),rmix_bc)
  !
  ! Compute the total temperature (dimensional) at the inflow plane
  !
  ttot_bc = ptot_bc / (rtot_bc*rmix_bc)
  !
  ! Compute the total enthalpy at the inflow plane
  !
  call compute_mixture_enthalpy(pv_bc(1:ncs),ttot_bc,htot_bc)
  !
  ! Copy the specified species mass fractions into the boundary state
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv,pv_bc)
    do m = 1,ncs
      do n = n1,n2
        face_cv(n,m,2) = pv_bc(m)
      end do
    end do
  end if
  !
  ! Compute static temperature using a Newton-Raphson method
  !
  convergence_reached = fals
  newton_iteration: &
  do k = 1,newton_iter
    !
    call compute_mixture_enthalpy(n1,n2,face_cv(:,1:ncs,2), &
                                        face_gp(:,ntemp,2), &
                                        face_gp(:,nenth,2))
    !
    call compute_mixture_cp(n1,n2,face_cv(:,1:ncs,2), &
                                  face_gp(:,ntemp,2), &
                                  face_gp(:,ncp,2))
    !
    !$ACC PARALLEL LOOP PRESENT(face_gp)
    do n = n1,n2
      face_gp(n,ngamm,2) = rmix_bc
    end do
    !
    call compute_mixture_gamma(n1,n2,face_gp(:,ncp,2),face_gp(:,ngamm,2))
    !
    ndnc = 0
    dtmx = zero
    !
    !$ACC PARALLEL LOOP PRESENT(face_cv,face_gp) &
    !$ACC&              REDUCTION(+:ndnc) &
    !$ACC&              REDUCTION(max:dtmx)
    do n = n1,n2
      !
      mach = face_cv(n,nee,2)
      !
      temp = face_gp(n,ntemp,2)
      gamm = face_gp(n,ngamm,2)
      enth = face_gp(n,nenth,2)
      cp   = face_gp(n,ncp,2)
      !
      htot = enth + half * gamm*rmix_bc*temp * mach*mach
      !
      func = htot - htot_bc
      dfunc = cp + half*gamm*rmix_bc * mach*mach
      !
      dtemp = func/dfunc
     !func = (htot - htot_bc) / (cp + half*gamm*rmix_bc*mach**2)
      !
      face_gp(n,ntemp,2) = min( tmax_nd , max(tmin_nd,temp-dtemp) )
      !
      dterr = abs(dtemp)/face_gp(n,ntemp,2)
      !
      if (dterr > dtmx)  dtmx = dterr
      if (dterr > conv_temp) ndnc = ndnc + 1
      !
    end do
    !
    if (dtmx <= conv_temp) then
      convergence_reached = true
      exit newton_iteration
    end if
    !
  end do newton_iteration
  !
  if (.not. convergence_reached) then
    newton_failures = ndnc
  end if
  !
  ! Temporarily store the boundary condition total temperature in
  ! the enthalpy location of the boundary state gas properties
  !
  !$ACC PARALLEL LOOP PRESENT(face_gp)
  do n = n1,n2
    face_gp(n,nenth,2) = ttot_bc
  end do
  !
  ! Compute the integral of Cp/T for each species and sum
  ! NOTE: face_gp(nenth) contains nondimensional BC total temperature
  !       face_gp(ntemp) contains nondimensional static temperature
  !       face_gp(npres) contains nondimensional integral of Cp/T
  !
  call compute_cp_over_temp_integral(n1,n2,face_cv(:,1:ncs,2), &
                                           face_gp(:,npres,2), &
                                           face_gp(:,ntemp,2), &
                                           face_gp(:,nenth,2))
  !
  ! Compute static pressure
  !
  !$ACC PARALLEL LOOP PRESENT(face_gp)
  do n = n1,n2
    face_gp(n,npres,2) = ptot_bc / exp(face_gp(n,npres,2)/rmix_bc)
  end do
  !
  ! Compute thermodynamic properties given pressure, temperature, and
  ! mass fractions
  !
  ! Compute gas constant of the mixture
  ! NOTE: Rmix is stored in face_gp(:,ngamm,2) temporarily
  !
  call compute_mixture_rgas(n1,n2,face_cv(:,1:ncs,2),face_gp(:,ngamm,2))
  !
  ! Use the ideal gas law to compute the updated density now that we
  ! have the updated mixture gas constant, temperature, and pressure.
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,face_gp)
  do n = n1,n2
    rmix = face_gp(n,ngamm,2)
    temp = face_gp(n,ntemp,2)
    pres = face_gp(n,npres,2)
    face_cv(n,nec,2) = pres / (rmix*temp)
  end do
  !
  ! Compute the mixture static enthalpy
  !
  call compute_mixture_enthalpy(n1,n2,face_cv(:,1:ncs,2), &
                                      face_gp(:,ntemp,2), &
                                      face_gp(:,nenth,2))
  !
  ! Compute the mixture specific heat at constant pressure, Cp
  !
  call compute_mixture_cp(n1,n2,face_cv(:,1:ncs,2), &
                                face_gp(:,ntemp,2), &
                                face_gp(:,ncp,2))
  !
  ! Compute the ratio of the specific heats for the gas mixture
  ! NOTE: Upon entry to compute_mixture_gamma:
  !           face_gp(:,ngamm,2) = Rgas_mix (the mixture gas constant)
  !       After returning from compute_mixture_gamma:
  !           face_gp(:,ngamm,2) = gamma_mix (the mixture gamma)
  !
  call compute_mixture_gamma(n1,n2,face_gp(:,ncp,2),face_gp(:,ngamm,2))
  !
  ! Assign boundary state density, momentum, and provisional total energy
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,face_gp)
  do n = n1,n2
    !
    ! Extract interior mach number which was temporarily
    ! stored in the total energy location of face_cv
    !
    mach = face_cv(n,nee,2)
    !
    ! Get the other flow conditions needed to compute the
    ! conserved variables at the boundary
    !
    rho = face_cv(n,nec,2)
    !
    pres = face_gp(n,npres,2)
    temp = face_gp(n,ntemp,2)
    gamm = face_gp(n,ngamm,2)
    enth = face_gp(n,nenth,2)
    !
    vel = mach * sqrt(max(gamm*pres/rho,rndoff))
    !
    ! Compute the conserved variables for this boundary point
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rho*vel*dir_cosx
   !face_cv(n,nee,2) = rho*enth - pres + half*rho*vel*vel
    face_cv(n,nee,2) = rho*htot_bc - rho*rmix_bc*face_gp(n,ntemp,2)
    !
  end do
  !
  ! Convert mass fractions to densities
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv)
    do m = 1,ncs
      do n = n1,n2
        face_cv(n,m,2) = face_cv(n,m,2) * face_cv(n,nec,2)
      end do
    end do
  end if
  !
  ! Finished setting subsonic inflow boundary conditions
  !
end subroutine get_bc_sub_inflow_sol_tpg_1d
!
!###############################################################################
!
end module bc_sub_in_mod
