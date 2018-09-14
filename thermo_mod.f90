module thermo_mod
  !
  !.. Use Statements ..
  use module_kind_types
  !
  use eqn_idx, only : nec,nmb,nme,nmx,nmy,nmz,nee,nq
  use eqn_idx, only : ncs
  use eqn_idx, only : npres,ntemp,ncp,ngamm,nenth,ng
  !
  implicit none
  !
  private
  !
  public :: init_thermo_mod
  public :: update_thermo_properties
  public :: compute_species_enthalpy
  public :: compute_species_cp
  public :: compute_mixture_dcpdt
  !
  !
  real(wp), parameter :: runivl = 8314.5_wp
  !
  integer, protected :: ncpc = 5
  !
  !
  !
  real(wp), public, protected, save, allocatable :: wtmol(:)
  !
  character(len=16), public, protected, save, allocatable :: species_name(:)
  !
  real(wp), public, protected, save, allocatable :: rspec(:)
  !
  real(wp), protected, save, allocatable :: cp_coef(:,:,:)
  real(wp), protected, save, allocatable :: enth_coef(:,:,:)
  !
  real(wp), protected, save :: temperature_switch
  !
  real(wp), protected, save, allocatable :: mass_frac_ref(:)
  !
  real(wp), protected, save, allocatable :: t_limit(:,:)
  !
  integer, parameter :: Use_VULCAN_Coefficients = 1
  integer, parameter :: Use_NCC_Coefficients = 2
  !
  integer, parameter :: Thermo_Coefficients = Use_NCC_Coefficients
  !
  !
  ! #########################################
  ! #####   MODULE GENERIC INTERFACES   #####
  ! #########################################
  !
  interface update_thermo_properties
    module procedure update_thermo_properties_sp, &
                     update_thermo_properties_fp, &
                     update_thermo_properties_1p
  end interface update_thermo_properties
  !
contains
!
!###############################################################################
!
subroutine init_thermo_mod()
  !
  !.. Use Statements ..
  use ovar, only : itestcase
  !
  !.. Local Scalars ..
  integer  :: i,ierr
  real(wp) :: mf_sum
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_thermo_mod"
  !
continue
  !
  allocate ( species_name(1:ncs) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"species_name",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( mass_frac_ref(1:ncs) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"mass_frac_ref",1,__LINE__,__FILE__,ierr,error_message)
  !
  if (itestcase == CH4_3D_Diffusion_Flame) then
    !
    species_name(1) = "CO2"
    species_name(2) = "CO"
    species_name(3) = "H2O"
    species_name(4) = "CH4"
    species_name(5) = "N2"
    species_name(6) = "O2"
    !
    mass_frac_ref(1) = zero   ! CO2
    mass_frac_ref(2) = zero   ! CO
    mass_frac_ref(3) = zero   ! H2O
    mass_frac_ref(4) = 0.1_wp ! CH4
    mass_frac_ref(5) = 0.7_wp ! N2
    mass_frac_ref(6) = 0.2_wp ! O2
    !
  else if (itestcase == H2_2D_Diffusion_Flame) then
    !
    species_name(1) = "N2"
    species_name(2) = "O"
    species_name(3) = "OH"
    species_name(4) = "O2"
    species_name(5) = "H"
    species_name(6) = "H2"
    species_name(7) = "H2O"
    !
    mass_frac_ref(1) = 0.7686_wp ! N2
    mass_frac_ref(2) = zero      ! O
    mass_frac_ref(3) = zero      ! OH
    mass_frac_ref(4) = 0.2314_wp ! O2
    mass_frac_ref(5) = zero      ! H
    mass_frac_ref(6) = zero      ! H2
    mass_frac_ref(7) = zero      ! H2O
    !
  else if (itestcase == Reacting_Backstep) then
    !
    species_name(1) = "N2"
    species_name(2) = "O2"
    !
    mass_frac_ref(1) = 0.7686 ! N2
    mass_frac_ref(2) = 0.2314 ! O2
    !
  end if
  !
  ! Make sure the reference mass fractions sum to 1
  !
  mf_sum = sum( mass_frac_ref(1:ncs) )
  do i = 1,ncs
    mass_frac_ref(i) = mass_frac_ref(i) / mf_sum
  end do
  !
  ! Allocate the remaining arrays for species thermodynamic properties
  !
  allocate ( wtmol(1:ncs) , source=zero , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"wtmol",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( t_limit(1:2,1:ncs) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"t_limit",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( rspec(1:ncs) , source=zero , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"rspec",1,__LINE__,__FILE__,ierr,error_message)
  !
  ncpc = 5
  !
  allocate ( cp_coef(1:ncpc,1:2,1:ncs) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cp_coef",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( enth_coef(1:ncpc+1,1:2,1:ncs) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"enth_coef",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Fill the arrays with the species thermodynamic properties
  !
  if (Thermo_Coefficients == Use_NCC_Coefficients) then
    call hardwired_NCC_coefficients
  else if (Thermo_Coefficients == Use_Vulcan_Coefficients) then
    call hardwired_VULCAN_coefficients
  end if
  !
end subroutine init_thermo_mod
!
!###############################################################################
!
pure subroutine update_thermo_properties_1p(cv_1p,gp_1p)
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: cv_1p(:)
  real(wp), intent(inout) :: gp_1p(:)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: cv(:,:)
  real(wp), allocatable :: gp(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_thermo_properties_1p"
  !
continue
  !
  allocate ( cv(1:1,1:nq) , source=zero )
  allocate ( gp(1:1,1:ng) , source=zero )
  !
  cv(1,1:nq) = cv_1p(1:nq)
  gp(1,1:nq) = gp_1p(1:ng)
  !$ACC DATA CREATE(cv,gp)
  !
  call update_thermo_properties_sp(1,1,cv,gp)
  !
  gp_1p(1:ng) = gp(1,1:ng)
  !
  !$ACC END DATA
  deallocate ( gp )
  deallocate ( cv )
  !
end subroutine update_thermo_properties_1p
!
!###############################################################################
!
pure subroutine update_thermo_properties_fp(n1,n2,cv,gp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp), intent(inout) :: cv(:,:,:)
  real(wp), intent(inout) :: gp(:,:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_thermo_properties_fp"
  !
continue
  !
  call update_thermo_properties_sp(n1,n2,cv(:,:,1),gp(:,:,1))
  if (all([size(cv,dim=3),size(gp,dim=3)] > 1)) then
    call update_thermo_properties_sp(n1,n2,cv(:,:,2),gp(:,:,2))
  end if
  !
end subroutine update_thermo_properties_fp
!
!###############################################################################
!
pure subroutine update_thermo_properties_sp(n1,n2,cv,gp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp), intent(inout) :: cv(:,:)
  real(wp), intent(inout) :: gp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,n,num_newton_errors
  real(wp) :: rho,max_newton_error
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_thermo_properties_sp"
  !
continue
  !
  ! 1. Update the mixture density by summing the species densities
  !    and convert the species densities to mass fractions
  !
  if (ncs > 0) then
    !
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(cv)
    do i = 1,ncs
      do n = n1,n2
        cv(n,i) = max(cv(n,i),zero)
      end do
    end do
    !
    !$ACC PARALLEL LOOP PRESENT(cv)
    do n = n1,n2
      cv(n,nec) = zero
    end do
    !
    !$ACC PARALLEL LOOP PRESENT(cv)
    do n = n1,n2
      rho = zero
      !$ACC LOOP SEQ REDUCTION(+:rho)
      do i = 1,ncs
        rho = rho + cv(n,i)
      end do
      cv(n,nec) = rho
    end do
    !
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(cv)
    do i = 1,ncs
      do n = n1,n2
        cv(n,i) = cv(n,i)/cv(n,nec)
      end do
    end do
    !
    ! 2. Compute the internal energy. Temporarily store the internal energy
    !    in the pressure location of gp.
    !
    call compute_internal_energy(n1,n2,cv,gp(:,npres))
    !
    ! 3. Compute the mixture gas constant. Temporarily store the mixture
    !    gas constant in the gamma location of gp.
    !
    call compute_mixture_rgas(n1,n2,cv(:,1:ncs),gp(:,ngamm))
    !
    ! 4. Compute static temperature and gas variables using Newton iteration
    !    NOTE: Currently, the gamma and pressure locations of gp respectively
    !          contain the mixture gas constant and the internal energy.
    !
    call newton_solve_for_temperature(n1,n2,cv(:,1:ncs),gp(:,ngamm), &
                                            gp(:,npres),gp(:,ntemp), &
                                            gp(:,nenth),gp(:,ncp), &
                                            num_newton_errors, &
                                            max_newton_error, &
                                            dim_t_upon_entry=true, &
                                            nondim_t_upon_return=true)
    !
    ! 5. Compute the static pressure
    !    NOTE: The gamma location of gp currently
    !          contains the mixture gas constant.
    !
    !$ACC PARALLEL LOOP PRESENT(gp,cv)
    do n = n1,n2
      gp(n,npres) = cv(n,nec) * gp(n,ngamm) * gp(n,ntemp)
    end do
    !
    ! 6. Compute the ratio of the specific heats of the mixture.
    !    NOTE: gp(:,ngamm) contains the mixture gas constant going into
    !          compute_mixture_gamma, but returns containing the ratio of
    !          the specific heats of the mixture.
    !
    call compute_mixture_gamma(n1,n2,gp(:,ncp),gp(:,ngamm))
   !!
   !! 7. Nondimensionalize temperature, pressure, enthalpy, cp
   !!
   !!$ACC PARALLEL LOOP PRESENT(gp)
   !do n = n1,n2
   !  gp(n,npres) = gp(n,npres) * rhoa2_inv
   !  gp(n,ntemp) = gp(n,ntemp) * tref_inv
   !  gp(n,nenth) = gp(n,nenth) * aref2_inv
   !  gp(n,ncp)   = gp(n,ncp)   * cpref_inv
   !end do
    !
    ! 9. Convert species mass fractions back to species densities
    !
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(cv)
    do i = 1,ncs
      do n = n1,n2
        cv(n,i) = cv(n,i)*cv(n,nec)
      end do
    end do
    !
  else
    !
    ! Compute the pressure and temperature of a calorically perfect gas
    ! 
    call update_thermo_cpg(n1,n2,cv,gp)
    !
  end if
  !
end subroutine update_thermo_properties_sp
!
!###############################################################################
!
pure &
subroutine update_thermo_cpg(n1,n2,cv,gp)
  !
  ! N1 : beginning solution point
  ! N1 : ending    solution point
  ! CV : conserved variables at points
  ! GP : gas properties variables at points
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: cv(:,:)
  real(wp), intent(inout) :: gp(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: gm1,rgasref_nd_inv
  real(wp) :: rho,rhou,rhov,rhow,rhoe
  real(wp) :: rho_inv,rhoq2,pres,temp
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_thermo_cpg"
  !
continue
  !
  nr = nme-nmb+1
  !
  gm1 = gam - one
  rgasref_nd_inv = one/rgasref_nd
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(cv,gp)
    do n = n1,n2
      !
      rho = cv(n,nec)
      rho_inv = one/rho
      !
      rhou = cv(n,nmx)
      rhov = cv(n,nmy)
      rhow = cv(n,nmz)
      rhoe = cv(n,nee)
      !
      rhoq2 = rhou*rhou + rhov*rhov + rhow*rhow
      !
      pres = gm1*(rhoe - half*rhoq2*rho_inv)
      temp = pres * rgasref_nd_inv * rho_inv
      !
      gp(n,npres) = pres
      gp(n,ntemp) = temp
      !
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(cv,gp)
    do n = n1,n2
      !
      rho = cv(n,nec)
      rho_inv = one/rho
      !
      rhou = cv(n,nmx)
      rhov = cv(n,nmy)
      rhoe = cv(n,nee)
      !
      rhoq2 = rhou*rhou + rhov*rhov
      !
      pres = gm1*(rhoe - half*rhoq2*rho_inv)
      temp = pres * rgasref_nd_inv * rho_inv
      !
      gp(n,npres) = pres
      gp(n,ntemp) = temp
      !
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(cv,gp)
    do n = n1,n2
      !
      rho = cv(n,nec)
      rho_inv = one/rho
      !
      rhou = cv(n,nmx)
      rhoe = cv(n,nee)
      !
      pres = gm1*(rhoe - half*rhou*rhou*rho_inv)
      temp = pres * rgasref_nd_inv * rho_inv
      !
      gp(n,npres) = pres
      gp(n,ntemp) = temp
      !
    end do
  end if
  !
end subroutine update_thermo_cpg
!
!###############################################################################
!
pure &
subroutine compute_mixture_gamma(n1,n2,cp,gam)
  !
  ! N1  : beginning solution point
  ! N1  : ending    solution point
  ! CP  : mixture specific heat at constant pressure
  ! GAM : in  - mixture gas constant
  !       out - ratio of the specific heats for the mixture
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: cp(:)
  real(wp), intent(inout) :: gam(:)
  !
  !.. Local Scalars ..
  integer :: n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_mixture_gamma"
  !
continue
  !
  !$ACC PARALLEL LOOP PRESENT(cp,gam)
  do n = n1,n2
    gam(n) = cp(n) / (cp(n)-gam(n))
  end do
  !
end subroutine compute_mixture_gamma
!
!###############################################################################
!
pure &
subroutine compute_mixture_rgas(n1,n2,mass_frac,rmix)
  !
  !.. Use Statements ..
  use ovar, only : aref,tref
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: mass_frac(:,:)
  real(wp), intent(inout) :: rmix(:)
  !
  !.. Local Scalars ..
  integer  :: i,n
  real(wp) :: rsum,tda2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_mixture_rgas"
  !
continue
  !
  tda2 = tref/aref/aref
  !
  !$ACC PARALLEL LOOP PRESENT(mass_frac,rmix,rspec)
  do n = n1,n2
    rsum = zero
    !$ACC LOOP SEQ REDUCTION(+:rsum)
    do i = 1,ncs
      rsum = rsum + rspec(i)*mass_frac(n,i)
    end do
    rmix(n) = rsum*tda2
  end do
  !
end subroutine compute_mixture_rgas
!
!###############################################################################
!
pure &
subroutine newton_solve_for_temperature(n1,n2,mass_frac,rmix,eint,t,h,cp, &
                                        ndnc,dtmx, &
                                        dim_t_upon_entry,nondim_t_upon_return)
  !
  !.. Use Statements ..
  use ovar, only : tref,tmin,tmax
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: mass_frac(:,:)
  real(wp),    intent(in) :: rmix(:)
  real(wp),    intent(in) :: eint(:)
  real(wp), intent(inout) :: t(:)
  real(wp), intent(inout) :: h(:)
  real(wp), intent(inout) :: cp(:)
  integer,  intent(inout) :: ndnc
  real(wp), intent(inout) :: dtmx
  !
  !.. Optional Arguments
  logical(lk), optional, intent(in) :: dim_t_upon_entry
  logical(lk), optional, intent(in) :: nondim_t_upon_return
  !
  !.. Local Scalars ..
  integer  :: k,n
  real(wp) :: func,dfunc,dtemp,dterr
  real(wp) :: tmin_nd,tmax_nd
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "newton_solve_for_temperature"
  !
  integer, parameter :: newton_iter = 20
  !
 !real(wp), parameter :: conv_temp = eps2  ! NCC code
  real(wp), parameter :: conv_temp = eps6
 !real(wp), parameter :: conv_temp = eps8  ! VULCAN if ref. Mach > 0.01
 !real(wp), parameter :: conv_temp = eps12 ! VULCAN if ref. Mach <= 0.01
  !
continue
  !
  tmin_nd = tmin / tref
  tmax_nd = tmax / tref
  !
  ! Compute initial static enthalpy and specific heat of the mixture
  !
  call compute_mixture_enthalpy(n1,n2,mass_frac,t,h)
  !
  call compute_mixture_cp(n1,n2,mass_frac,t,cp)
  !
  ! Compute the temperature using Newton solve
  !
  newton_iteration : &
  do k = 1,newton_iter
    !
    ndnc = 0
    dtmx = zero
    !
    do n = n1,n2
      !
      func = h(n) - rmix(n)*t(n) - eint(n)
      dfunc = cp(n) - rmix(n)
      !
      dtemp = func/dfunc
      !
      t(n) = min( tmax_nd , max(tmin_nd,t(n)-dtemp) )
      !
      dterr = abs(dtemp) / t(n)
      !
      if (dterr > dtmx) dtmx = dterr
      !
      if (dterr > conv_temp) ndnc = ndnc + 1
      !
    end do
    !
    ! Update the static enthalpy and specific heat of the mixture
    !
    call compute_mixture_enthalpy(n1,n2,mass_frac,t,h)
    !
    call compute_mixture_cp(n1,n2,mass_frac,t,cp)
    !
    ! Check for convergence
    !
    if (dtmx <= conv_temp) then
      exit newton_iteration
    end if
    !
   !if (k == newton_iter) then
   !  write (iout,100) ndnc,dtmx
   !end if
    !
  end do newton_iteration
  !
  ! Format Statements
  !
  100 format (/, &
   1x,64("*"),/, &
   1x,5("*"),7x,"NEWTON ITERATIONTO SOLVE FOR THE STATIC",t61,5("*"),/, &
   1x,5("*"),4x,"TEMPERATURE OF A GAS MIXTURE DID NOT CONVERGE!",t61,5("*"),/, &
   1x,5("*"),7x,"NUMBER OF UNCONVERGED POINTS = ",i0,t61,5("*"),/, &
   1x,5("*"),7x,"MAXIMUM ERROR IN TEMPERATURE = ",es13.6,t61,5("*"),/)
  !
end subroutine newton_solve_for_temperature
!
!###############################################################################
!
pure &
subroutine compute_internal_energy(n1,n2,cv,e)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: cv(:,:)
  real(wp), intent(inout) :: e(:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: rho_inv,velx,vely,velz,q2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_internal_energy"
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    !
    !$ACC PARALLEL LOOP PRESENT(cv,e)
    do n = n1,n2
      !
      rho_inv = one/cv(n,nec)
      !
      velx = cv(n,nmx)*rho_inv
      vely = cv(n,nmy)*rho_inv
      velz = cv(n,nmz)*rho_inv
      !
      q2 = velx*velx + vely*vely + velz*velz
      !
      e(n) = cv(n,nee)*rho_inv - half*q2
      !
    end do
    !
  else if (nr == 2) then
    !
    !$ACC PARALLEL LOOP PRESENT(cv,e)
    do n = n1,n2
      !
      rho_inv = one/cv(n,nec)
      !
      velx = cv(n,nmx)*rho_inv
      vely = cv(n,nmy)*rho_inv
      !
      q2 = velx*velx + vely*vely
      !
      e(n) = cv(n,nee)*rho_inv - half*q2
      !
    end do
    !
  else
    !
    !$ACC PARALLEL LOOP PRESENT(cv,e)
    do n = n1,n2
      !
      rho_inv = one/cv(n,nec)
      !
      velx = cv(n,nmx)*rho_inv
      !
      e(n) = cv(n,nee)*rho_inv - half*velx*velx
      !
    end do
    !
  end if
  !
end subroutine compute_internal_energy
!
!###############################################################################
!
pure &
subroutine compute_mixture_enthalpy(n1,n2,mass_frac,t,h)
  !
  !.. Use Statements ..
  use ovar, only : aref,tref
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: mass_frac(:,:)
  real(wp),    intent(in) :: t(:)
  real(wp), intent(inout) :: h(:)
  !
  !.. Local Scalars ..
  integer  :: i,n,ni
  real(wp) :: hi,ts,aref2_inv
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: intrvl(:,:)
  !
continue
  !
  aref2_inv = one/aref/aref
  !
  allocate ( intrvl(n1:n2,1:ncs) , source=1 )
  !$ACC DATA CREATE(intrvl)
  !
  !$ACC PARALLEL LOOP PRESENT(h)
  do n = n1,n2
    h(n) = zero
  end do
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(t,intrvl)
  do i = 1,ncs
    do n = n1,n2
      if (t(n) > temperature_switch) then
        intrvl(n,i) = 2
      end if
    end do
  end do
  !
  do i = 1,ncs
    !$ACC PARALLEL LOOP PRESENT(mass_frac,intrvl,t,enth_coef,h)
    do n = n1,n2
      !
      ni = intrvl(n,i)
      ts = t(n)*tref
      !
      hi = enth_coef(6,ni,i) + &
           ts*(enth_coef(1,ni,i) + &
               ts*(enth_coef(2,ni,i) + &
                   ts*(enth_coef(3,ni,i) + &
                       ts*(enth_coef(4,ni,i) + &
                           ts*(enth_coef(5,ni,i))))))
      !
      h(n) = h(n) + mass_frac(n,i)*hi*aref2_inv
      !
    end do
  end do
  !
  !$ACC END DATA
  deallocate ( intrvl )
  !
end subroutine compute_mixture_enthalpy
!
!###############################################################################
!
pure &
subroutine compute_species_enthalpy(n1,n2,t,hsi)
  !
  !.. Use Statements ..
  use ovar, only : aref,tref
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: t(:)
  real(wp), intent(inout) :: hsi(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,n,ni
  real(wp) :: ts,aref2_inv
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: intrvl(:,:)
  !
continue
  !
  aref2_inv = one/aref/aref
  !
  allocate ( intrvl(n1:n2,1:ncs) , source=1 )
  !$ACC DATA CREATE(intrvl)
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(t,intrvl,temperature_switch)
  do i = 1,ncs
    do n = n1,n2
      if (t(n) > temperature_switch) then
        intrvl(n,i) = 2
      end if
    end do
  end do
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(intrvl,t,enth_coef,hsi)
  do i = 1,ncs
    do n = n1,n2
      !
      ni = intrvl(n,i)
      ts = t(n)*tref
      !
      hsi(n,i) = enth_coef(6,ni,i) + &
                 ts*(enth_coef(1,ni,i) + &
                     ts*(enth_coef(2,ni,i) + &
                         ts*(enth_coef(3,ni,i) + &
                             ts*(enth_coef(4,ni,i) + &
                                 ts*(enth_coef(5,ni,i))))))
      !
    end do
  end do
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(hsi)
  do i = 1,ncs
    do n = n1,n2
      hsi(n,i) = hsi(n,i)*aref2_inv
    end do
  end do
  !
  !$ACC END DATA
  deallocate ( intrvl )
  !
end subroutine compute_species_enthalpy
!
!###############################################################################
!
pure &
subroutine compute_mixture_cp(n1,n2,mass_frac,t,cp)
  !
  !.. Use Statements ..
  use ovar, only : aref,tref
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: mass_frac(:,:)
  real(wp),    intent(in) :: t(:)
  real(wp), intent(inout) :: cp(:)
  !
  !.. Local Scalars ..
  integer  :: i,n,ni
  real(wp) :: cpi,ts,tda2
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: intrvl(:,:)
  !
continue
  !
  tda2 = tref/aref/aref
  !
  allocate ( intrvl(n1:n2,1:ncs) , source=1 )
  !$ACC DATA CREATE(intrvl)
  !
  !$ACC PARALLEL LOOP PRESENT(cp)
  do n = n1,n2
    cp(n) = zero
  end do
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(t,intrvl,temperature_switch)
  do i = 1,ncs
    do n = n1,n2
      if (t(n) > temperature_switch) then
        intrvl(n,i) = 2
      end if
    end do
  end do
  !
  do i = 1,ncs
    !$ACC PARALLEL LOOP PRESENT(mass_frac,intrvl,t,cp_coef,cp)
    do n = n1,n2
      !
      ni = intrvl(n,i)
      ts = t(n)*tref
      !
      cpi = cp_coef(1,ni,i) + &
            ts*(cp_coef(2,ni,i) + &
                ts*(cp_coef(3,ni,i) + &
                    ts*(cp_coef(4,ni,i) + &
                        ts*(cp_coef(5,ni,i)))))
      !
      cp(n) = cp(n) + mass_frac(n,i)*cpi*tda2
      !
    end do
  end do
  !
  !$ACC END DATA
  deallocate ( intrvl )
  !
end subroutine compute_mixture_cp
!
!###############################################################################
!
pure &
subroutine compute_species_cp(n1,n2,t,cpi)
  !
  !.. Use Statements ..
  use ovar, only : aref,tref
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: t(:)
  real(wp), intent(inout) :: cpi(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,n,ni
  real(wp) :: ts,tda2
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: intrvl(:,:)
  !
continue
  !
  tda2 = tref/aref/aref
  !
  allocate ( intrvl(n1:n2,1:ncs) , source=1 )
  !$ACC DATA CREATE(intrvl)
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(t,intrvl,temperature_switch)
  do i = 1,ncs
    do n = n1,n2
      if (t(n) > temperature_switch) then
        intrvl(n,i) = 2
      end if
    end do
  end do
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(intrvl,t,cp_coef,cpi)
  do i = 1,ncs
    do n = n1,n2
      !
      ni = intrvl(n,i)
      ts = t(n)*tref
      !
      cpi(n,i) = cp_coef(1,ni,i) + &
                 ts*(cp_coef(2,ni,i) + &
                     ts*(cp_coef(3,ni,i) + &
                         ts*(cp_coef(4,ni,i) + &
                             ts*(cp_coef(5,ni,i)))))
      !
    end do
  end do
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(cpi)
  do i = 1,ncs
    do n = n1,n2
      cpi(n,i) = cpi(n,i)*tda2
    end do
  end do
  !
  !$ACC END DATA
  deallocate ( intrvl )
  !
end subroutine compute_species_cp
!
!###############################################################################
!
pure &
subroutine compute_mixture_dcpdt(n1,n2,mass_frac,t,dcpdt)
  !
  ! Computes the following sum for points n1 to n2:
  !
  !       NCS
  !     ------
  !     \                         d [Cp(i)]
  !      \     MassFraction(i) * -----------
  !      /                           dT
  !     /
  !     ------
  !       i=1
  !
  !.. Use Statements ..
  use ovar, only : aref,tref
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: mass_frac(:,:)
  real(wp),    intent(in) :: t(:)
  real(wp), intent(inout) :: dcpdt(:)
  !
  !.. Local Scalars ..
  integer  :: i,n,ni
  real(wp) :: dcpidt,ts,aref2_inv
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: intrvl(:,:)
  !
continue
  !
  aref2_inv = one/aref/aref
  !
  allocate ( intrvl(n1:n2,1:ncs) , source=1 )
  !$ACC DATA CREATE(intrvl)
  !
  !$ACC PARALLEL LOOP PRESENT(dcpdt)
  do n = n1,n2
    dcpdt(n) = zero
  end do
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(t,intrvl,temperature_switch)
  do i = 1,ncs
    do n = n1,n2
      if (t(n) > temperature_switch) then
        intrvl(n,i) = 2
      end if
    end do
  end do
  !
  do i = 1,ncs
    !$ACC PARALLEL LOOP PRESENT(mass_frac,intrvl,t,cp_coef,dcpdt)
    do n = n1,n2
      !
      ni = intrvl(n,i)
      ts = t(n)*tref
      !
      dcpidt = cp_coef(2,ni,i) + &
               ts*(cp_coef(3,ni,i)*two + &
                   ts*(cp_coef(4,ni,i)*three + &
                       ts*(cp_coef(5,ni,i)*four)))
      !
      dcpdt(n) = dcpdt(n) + mass_frac(n,i)*dcpidt*aref2_inv
      !
    end do
  end do
  !
  !$ACC END DATA
  deallocate ( intrvl )
  !
end subroutine compute_mixture_dcpdt
!
!###############################################################################
!
subroutine hardwired_NCC_coefficients
  !
  !.. Use Statements ..
  use ovar, only : tref
  !
  real(wp), parameter :: NCC_interval_switch = 1000.0_wp
  !
  real(wp), parameter :: CO2_wtmol = 44.00995_wp
  !
  real(wp), parameter :: CO2_cp(1:6,1:2) = reshape( [ 0.02275724e+02_wp, &
                                                      0.09922072e-01_wp, &
                                                     -0.10409113e-04_wp, &
                                                      0.06866686e-07_wp, &
                                                     -0.02117280e-10_wp, &
                                                     -0.04837314e+06_wp, &
                                                      0.04453623e+02_wp, &
                                                      0.03140168e-01_wp, &
                                                     -0.12784105e-05_wp, &
                                                      0.02393996e-08_wp, &
                                                     -0.16690333e-13_wp, &
                                                     -0.04896696e+06_wp] , &
                                                    [6,2] )
  !
  !
  real(wp), parameter :: CO_wtmol = 28.01055_wp
  !
  real(wp), parameter :: CO_cp(1:6,1:2) = reshape( [ 0.03262451e+02_wp, &
                                                     0.15119409e-02_wp, &
                                                    -0.03881755e-04_wp, &
                                                     0.05581944e-07_wp, &
                                                    -0.02474951e-10_wp, &
                                                    -0.14310539e+05_wp, &
                                                     0.03025078e+02_wp, &
                                                     0.14426885e-02_wp, &
                                                    -0.05630827e-05_wp, &
                                                     0.10185813e-09_wp, &
                                                    -0.06910951e-13_wp, &
                                                    -0.14268350e+05_wp] , &
                                                   [6,2] )
  !
  !
  real(wp), parameter :: H2O_wtmol = 18.01534_wp
  !
  real(wp), parameter :: H2O_cp(1:6,1:2) = reshape( [ 0.03386842e+02_wp, &
                                                      0.03474982e-01_wp, &
                                                     -0.06354696e-04_wp, &
                                                      0.06968581e-07_wp, &
                                                     -0.02506588e-10_wp, &
                                                     -0.03020811e+06_wp, &
                                                      0.02672145e+02_wp, &
                                                      0.03056293e-01_wp, &
                                                     -0.08730260e-05_wp, &
                                                      0.12009964e-09_wp, &
                                                     -0.06391618e-13_wp, &
                                                     -0.02989921e+06_wp] , &
                                                    [6,2] )
  !
  !
  real(wp), parameter :: CH4_wtmol = 16.04303_wp
  !
  real(wp), parameter :: CH4_cp(1:6,1:2) = reshape( [ 0.07787415e+01_wp, &
                                                      0.01747668e+00_wp, &
                                                     -0.02783409e-03_wp, &
                                                      0.03049708e-06_wp, &
                                                     -0.12239307e-10_wp, &
                                                     -0.09825229e+05_wp, &
                                                      0.01683478e+02_wp, &
                                                      0.10237236e-01_wp, &
                                                     -0.03875128e-04_wp, &
                                                      0.06785585e-08_wp, &
                                                     -0.04503423e-12_wp, &
                                                     -0.10080787e+05_wp] , &
                                                    [6,2] )
  !
  !
  real(wp), parameter :: N2_wtmol = 28.01340_wp
  !
  real(wp), parameter :: N2_cp(1:6,1:2) = reshape( [ 0.03298677e+02_wp, &
                                                     0.14082404e-02_wp, &
                                                    -0.03963222e-04_wp, &
                                                     0.05641515e-07_wp, &
                                                    -0.02444854e-10_wp, &
                                                    -0.10208999e+04_wp, &
                                                     0.02926640e+02_wp, &
                                                     0.14879768e-02_wp, &
                                                    -0.05684760e-05_wp, &
                                                     0.10097038e-09_wp, &
                                                    -0.06753351e-13_wp, &
                                                    -0.09227977e+04_wp] , &
                                                   [6,2] )
  !
  !
  real(wp), parameter :: O2_wtmol = 31.99880_wp
  !
  real(wp), parameter :: O2_cp(1:6,1:2) = reshape( [ 0.03212936e+02_wp, &
                                                     0.11274864e-02_wp, &
                                                    -0.05756150e-05_wp, &
                                                     0.13138773e-08_wp, &
                                                    -0.08768554e-11_wp, &
                                                    -0.10052490e+04_wp, &
                                                     0.03697578e+02_wp, &
                                                     0.06135197e-02_wp, &
                                                    -0.12588420e-06_wp, &
                                                     0.01775281e-09_wp, &
                                                    -0.11364354e-14_wp, &
                                                    -0.12339301e+04_wp] , &
                                                   [6,2] )
  !
  real(wp), parameter :: H_wtmol = 1.0079400_wp
  !
  real(wp), parameter :: H_cp(1:6,1:2) = reshape( [ 0.02500000E+02_wp, &
                                                    zero, &
                                                    zero, &
                                                    zero, &
                                                    zero, &
                                                    0.02547162E+06_wp, &
                                                    0.02500000E+02_wp, &
                                                    zero, &
                                                    zero, &
                                                    zero, &
                                                    zero, &
                                                    0.02547162E+06_wp] , &
                                                  [6,2] )
  !
  real(wp), parameter :: H2_wtmol = 2.0158800_wp
  !
  real(wp), parameter :: H2_cp(1:6,1:2) = reshape( [ 0.02991423E+02_wp, &
                                                     0.07000644E-02_wp, &
                                                    -0.05633828E-06_wp, &
                                                    -0.09231578E-10_wp, &
                                                     0.15827519E-14_wp, &
                                                    -0.08350340E+04_wp, &
                                                     0.03298124E+02_wp, &
                                                     0.08249441E-02_wp, &
                                                    -0.08143015E-05_wp, &
                                                    -0.09475434E-09_wp, &
                                                     0.04134872E-11_wp, &
                                                    -0.10125209E+04_wp] , &
                                                   [6,2] )
  !
  real(wp), parameter :: O_wtmol = 15.9994000_wp
  !
  real(wp), parameter :: O_cp(1:6,1:2) = reshape( [ 0.02542059E+02_wp, &
                                                   -0.02755061E-03_wp, &
                                                   -0.03102803E-07_wp, &
                                                    0.04551067E-10_wp, &
                                                   -0.04368051E-14_wp, &
                                                    0.02923080E+06_wp, &
                                                    0.02946428E+02_wp, &
                                                   -0.16381665E-02_wp, &
                                                    0.02421031E-04_wp, &
                                                   -0.16028431E-08_wp, &
                                                    0.03890696E-11_wp, &
                                                    0.02914764E+06_wp] , &
                                                  [6,2] )
  !
  real(wp), parameter :: OH_wtmol = 17.0073400_wp
  !
  real(wp), parameter :: OH_cp(1:6,1:2) = reshape( [ 0.02882730E+02_wp, &
                                                     0.10139743E-02_wp, &
                                                    -0.02276877E-05_wp, &
                                                     0.02174683E-09_wp, &
                                                    -0.05126305E-14_wp, &
                                                     0.03886888E+05_wp, &
                                                     0.03637266E+02_wp, &
                                                     0.01850910E-02_wp, &
                                                    -0.16761646E-05_wp, &
                                                     0.02387202E-07_wp, &
                                                    -0.08431442E-11_wp, &
                                                     0.03606781E+05_wp] , &
                                                   [6,2] )
  !
  integer  :: i,j,n,ncoef
  real(wp) :: denom
  !
  character(len=len(species_name)) :: species
  !
  real(wp), allocatable :: ncc_cp(:,:,:)
  !
continue
  !
  allocate ( ncc_cp(1:6,1:2,1:ncs) , source=zero )
  !
  do i = 1,ncs
    !
    species = uppercase( trim(adjustl(species_name(i))) )
    !
    select case (trim(adjustl(species)))
      !
      case ("CO2")
        !
        wtmol(i) = CO2_wtmol
        ncc_cp(1:6,1:2,i) = CO2_cp(1:6,1:2)
        !
      case ("CO")
        !
        wtmol(i) = CO_wtmol
        ncc_cp(1:6,1:2,i) = CO_cp(1:6,1:2)
        !
      case ("H2O")
        !
        wtmol(i) = H2O_wtmol
        ncc_cp(1:6,1:2,i) = H2O_cp(1:6,1:2)
        !
      case ("CH4")
        !
        wtmol(i) = CH4_wtmol
        ncc_cp(1:6,1:2,i) = CH4_cp(1:6,1:2)
        !
      case ("N2")
        !
        wtmol(i) = N2_wtmol
        ncc_cp(1:6,1:2,i) = N2_cp(1:6,1:2)
        !
      case ("O2")
        !
        wtmol(i) = O2_wtmol
        ncc_cp(1:6,1:2,i) = O2_cp(1:6,1:2)
        !
      case ("H")
        !
        wtmol(i) = H_wtmol
        ncc_cp(1:6,1:2,i) = H_cp(1:6,1:2)
        !
      case ("H2")
        !
        wtmol(i) = H2_wtmol
        ncc_cp(1:6,1:2,i) = H2_cp(1:6,1:2)
        !
      case ("O")
        !
        wtmol(i) = O_wtmol
        ncc_cp(1:6,1:2,i) = O_cp(1:6,1:2)
        !
      case ("OH")
        !
        wtmol(i) = OH_wtmol
        ncc_cp(1:6,1:2,i) = OH_cp(1:6,1:2)
        !
    end select
    !
  end do
  !
  do n = 1,ncs
    rspec(n) = runivl / wtmol(n)
  end do
  !
  do n = 1,ncs
    do i = 1,size(ncc_cp,dim=2)
      do j = 1,ncpc
        cp_coef(j,i,n) = ncc_cp(j,i,n) * rspec(n)
      end do
    end do
  end do
  !
  do n = 1,ncs
    do i = 1,size(ncc_cp,dim=2)
      do j = 1,ncpc
        ncoef = j + min(6-(ncpc+1),0)
        denom = max(real(ncoef,kind=wp),one)
        enth_coef(j,i,n) = ncc_cp(j,i,n) * rspec(n) / denom
      end do
      enth_coef(ncpc+1,i,n) = ncc_cp(ncpc+1,i,n) * rspec(n)
    end do
  end do
  !
  temperature_switch = NCC_interval_switch/tref
  !
  deallocate ( ncc_cp )
  !
end subroutine hardwired_NCC_coefficients
!
!###############################################################################
!
subroutine hardwired_VULCAN_coefficients
  !
  !.. Use Statements ..
  use ovar, only : tref
  !
  real(wp), parameter :: vul_interval_switch = 1000.0_wp
  !
  real(wp), parameter :: CO2_wtmol = 44.0095000_wp
  !
  real(wp), parameter :: CO2_cp(1:7,1:2) = reshape( [ 2.25257700E+00_wp, &
                                                      9.82353535E-03_wp, &
                                                     -9.40516149E-06_wp, &
                                                      4.99311980E-09_wp, &
                                                     -1.13098807E-12_wp, &
                                                     -4.83625753E+04_wp, &
                                                      1.03255562E+01_wp, &
                                                      4.69229277E+00_wp, &
                                                      2.64831933E-03_wp, &
                                                     -9.50767176E-07_wp, &
                                                      1.51842546E-10_wp, &
                                                     -8.60487968E-15_wp, &
                                                     -4.90469718E+04_wp, &
                                                     -2.24622115E+00_wp] , &
                                                    [7,2] )
  !
  ! Species : Carbon Monoxide                Gurvich,1979 pt1 p25 pt2 p29.
  ! Transport data fit from 300 K to 1500 K  Svehla,1994.
  !
  real(wp), parameter :: CO_wtmol = 28.0101000_wp
  !
  real(wp), parameter :: CO_cp(1:7,1:2) = reshape( [ 3.64231417E+00_wp, &
                                                    -1.13846700E-03_wp, &
                                                     2.47202434E-06_wp, &
                                                    -7.01205135E-10_wp, &
                                                    -2.88361728E-13_wp, &
                                                    -1.43499254E+04_wp, &
                                                     3.25676435E+00_wp, &
                                                     3.05476345E+00_wp, &
                                                     1.33701825E-03_wp, &
                                                    -4.78351638E-07_wp, &
                                                     7.74904656E-11_wp, &
                                                    -4.61587625E-15_wp, &
                                                    -1.42680818E+04_wp, &
                                                     5.98462203E+00_wp] , &
                                                    [7,2] )
  !
  ! Species : Steam            Cox,1989. Woolley,1987. TRC(10/88) tuv25.
  ! Transport data fit from 300 K to 1500 K  Sengers & Watson,1986, Svehla,1994.
  !
  real(wp), parameter :: H2O_wtmol = 18.0152800_wp
  !
  real(wp), parameter :: H2O_cp(1:7,1:2) = reshape( [ 4.21338451E+00_wp, &
                                                     -2.13843352E-03_wp, &
                                                      6.78402160E-06_wp, &
                                                     -5.78650372E-09_wp, &
                                                      1.89595906E-12_wp, &
                                                     -3.02953868E+04_wp, &
                                                     -9.11950745E-01_wp, &
                                                      2.63316752E+00_wp, &
                                                      3.05564198E-03_wp, &
                                                     -8.19557831E-07_wp, &
                                                      1.04135761E-10_wp, &
                                                     -4.95949013E-15_wp, &
                                                     -2.98701572E+04_wp, &
                                                      7.12319891E+00_wp] , &
                                                     [7,2] )
  !
  ! Species : Methane                        Gurvich,1991 pt1 p44 pt2 p36.
  ! Transport data fit from 300 K to 1500 K  Boushehri,1987, Svehla,1994.
  !
  real(wp), parameter :: CH4_wtmol = 16.0424600_wp
  !
  real(wp), parameter :: CH4_cp(1:7,1:2) = reshape( [ 5.02395784E+00_wp, &
                                                     -1.25330563E-02_wp, &
                                                      4.59743916E-05_wp, &
                                                     -4.49448239E-08_wp, &
                                                      1.53601897E-11_wp, &
                                                     -1.02377171E+04_wp, &
                                                     -4.14928349E+00_wp, &
                                                      1.84141025E+00_wp, &
                                                      9.74004923E-03_wp, &
                                                     -3.18244614E-06_wp, &
                                                      5.11502852E-10_wp, &
                                                     -2.98571724E-14_wp, &
                                                     -1.00921820E+04_wp, &
                                                      8.83569259E+00_wp] , &
                                                     [7,2] )
  !
  ! Species : Nitrogen               Ref-Elm. Gurvich,1978 pt1 p280 pt2 p207.
  ! Transport data fit from 300 K to 1500 K  Boushehri,1987, Svehla,1994.
  !
  real(wp), parameter :: N2_wtmol = 28.0134000_wp
  !
  real(wp), parameter :: N2_cp(1:7,1:2) = reshape( [ 3.57565174E+00_wp, &
                                                    -5.08294729E-04_wp, &
                                                     5.68705824E-07_wp, &
                                                     1.24596694E-09_wp, &
                                                    -9.53471933E-13_wp, &
                                                    -1.05052493E+03_wp, &
                                                     2.78974534E+00_wp, &
                                                     2.95737286E+00_wp, &
                                                     1.38466908E-03_wp, &
                                                    -4.86418174E-07_wp, &
                                                     7.74752085E-11_wp, &
                                                    -4.54114158E-15_wp, &
                                                    -9.24683183E+02_wp, &
                                                     5.84752726E+00_wp] , &
                                                    [7,2] )
  !
  ! Species : Oxygen                   Ref-Elm. Gurvich,1989 pt1 p94 pt2 p9.
  ! Transport data fit from 300 K to 1500 K  Boushehri,1987, Svehla,1994.
  !
  real(wp), parameter :: O2_wtmol = 31.9988000_wp
  !
  real(wp), parameter :: O2_cp(1:7,1:2) = reshape( [ 3.80103438E+00_wp, &
                                                    -3.14189698E-03_wp, &
                                                     1.02483417E-05_wp, &
                                                    -1.01443245E-08_wp, &
                                                     3.43356459E-12_wp, &
                                                    -1.06574832E+03_wp, &
                                                     3.58099673E+00_wp, &
                                                     3.69191906E+00_wp, &
                                                     6.07103434E-04_wp, &
                                                    -1.17310026E-07_wp, &
                                                     1.60002783E-11_wp, &
                                                    -9.93511564E-16_wp, &
                                                    -1.22908552E+03_wp, &
                                                     3.24042872E+00_wp] , &
                                                    [7,2] )
  !
  real(wp), parameter :: H_wtmol = 1.0079400_wp
  !
  real(wp), parameter :: H_cp(1:7,1:2) = reshape( [ 2.50000000E+00_wp, &
                                                    zero, &
                                                    zero, &
                                                    zero, &
                                                    zero, &
                                                    2.54737080E+04_wp, &
                                                   -4.46682853E-01_wp, &
                                                    2.50000000E+00_wp, &
                                                    zero, &
                                                    zero, &
                                                    zero, &
                                                    zero, &
                                                    2.54737080E+04_wp, &
                                                   -4.46682853E-01_wp] , &
                                                  [7,2] )
  !
  real(wp), parameter :: H2_wtmol = 2.0158800_wp
  !
  real(wp), parameter :: H2_cp(1:7,1:2) = reshape( [ 2.50571195E+00_wp, &
                                                     6.60361016E-03_wp, &
                                                    -1.56447491E-05_wp, &
                                                     1.58868648E-08_wp, &
                                                    -5.72989941E-12_wp, &
                                                    -9.31057594E+02_wp, &
                                                     3.81454289E-02_wp, &
                                                     2.86238922E+00_wp, &
                                                     9.31641036E-04_wp, &
                                                    -1.96008073E-07_wp, &
                                                     2.48300032E-11_wp, &
                                                    -1.31375840E-15_wp, &
                                                    -7.81539065E+02_wp, &
                                                    -6.22896707E-01_wp] , &
                                                   [7,2] )
  !
  real(wp), parameter :: O_wtmol = 15.9994000_wp
  !
  real(wp), parameter :: O_cp(1:7,1:2) = reshape( [ 3.06486352E+00_wp, &
                                                   -2.40286593E-03_wp, &
                                                    4.16607889E-06_wp, &
                                                   -3.29904637E-09_wp, &
                                                    9.86514943E-13_wp, &
                                                    2.91309512E+04_wp, &
                                                    2.46709645E+00_wp, &
                                                    2.53870743E+00_wp, &
                                                   -2.12005320E-05_wp, &
                                                   -6.91650640E-09_wp, &
                                                    5.46894449E-12_wp, &
                                                   -5.14276968E-16_wp, &
                                                    2.92285500E+04_wp, &
                                                    4.95173847E+00_wp] , &
                                                  [7,2] )
  !
  real(wp), parameter :: OH_wtmol = 17.0073400_wp
  !
  real(wp), parameter :: OH_cp(1:7,1:2) = reshape( [ 3.93471682E+00_wp, &
                                                    -1.92496370E-03_wp, &
                                                     3.31005131E-06_wp, &
                                                    -2.43377158E-09_wp, &
                                                     8.07066518E-13_wp, &
                                                     3.37111948E+03_wp, &
                                                     1.26959021E-01_wp, &
                                                     2.83337597E+00_wp, &
                                                     1.11922731E-03_wp, &
                                                    -3.00339776E-07_wp, &
                                                     4.33461483E-11_wp, &
                                                    -2.51028340E-15_wp, &
                                                     3.69646446E+03_wp, &
                                                     5.87244494E+00_wp] , &
                                                   [7,2] )
  !
  integer  :: i,j,n,ncoef
  real(wp) :: denom
  !
  character(len=len(species_name)) :: species
  !
  real(wp), allocatable :: vul_cp(:,:,:)
  !
continue
  !
  allocate ( vul_cp(1:7,1:2,1:ncs) , source=zero )
  !
  do i = 1,ncs
    !
    species = uppercase( trim(adjustl(species_name(i))) )
    !
    select case (trim(adjustl(species)))
      !
      case ("CO2")
        !
        wtmol(i) = CO2_wtmol
        vul_cp(1:7,1:2,i) = CO2_cp(1:7,1:2)
        !
      case ("CO")
        !
        wtmol(i) = CO_wtmol
        vul_cp(1:7,1:2,i) = CO_cp(1:7,1:2)
        !
      case ("H2O")
        !
        wtmol(i) = H2O_wtmol
        vul_cp(1:7,1:2,i) = H2O_cp(1:7,1:2)
        !
      case ("CH4")
        !
        wtmol(i) = CH4_wtmol
        vul_cp(1:7,1:2,i) = CH4_cp(1:7,1:2)
        !
      case ("N2")
        !
        wtmol(i) = N2_wtmol
        vul_cp(1:7,1:2,i) = N2_cp(1:7,1:2)
        !
      case ("O2")
        !
        wtmol(i) = O2_wtmol
        vul_cp(1:7,1:2,i) = O2_cp(1:7,1:2)
        !
      case ("H")
        !
        wtmol(i) = H_wtmol
        vul_cp(1:7,1:2,i) = H_cp(1:7,1:2)
        !
      case ("H2")
        !
        wtmol(i) = H2_wtmol
        vul_cp(1:7,1:2,i) = H2_cp(1:7,1:2)
        !
      case ("O")
        !
        wtmol(i) = O_wtmol
        vul_cp(1:7,1:2,i) = O_cp(1:7,1:2)
        !
      case ("OH")
        !
        wtmol(i) = OH_wtmol
        vul_cp(1:7,1:2,i) = OH_cp(1:7,1:2)
        !
    end select
    !
  end do
  !
  do n = 1,ncs
    rspec(n) = runivl / wtmol(n)
  end do
  !
  do n = 1,ncs
    do i = 1,size(vul_cp,dim=2)
      do j = 1,ncpc
        cp_coef(j,i,n) = vul_cp(j,i,n) * rspec(n)
      end do
    end do
  end do
  !
  do n = 1,ncs
    do i = 1,size(vul_cp,dim=2)
      do j = 1,ncpc
        ncoef = j + min(6-(ncpc+1),0)
        denom = max(real(ncoef,kind=wp),one)
        enth_coef(j,i,n) = vul_cp(j,i,n) * rspec(n) / denom
      end do
      enth_coef(ncpc+1,i,n) = vul_cp(ncpc+1,i,n) * rspec(n)
    end do
  end do
  !
  temperature_switch = vul_interval_switch/tref
  !
  deallocate ( vul_cp )
  !
end subroutine hardwired_VULCAN_coefficients
!
!###############################################################################
!
end module thermo_mod
