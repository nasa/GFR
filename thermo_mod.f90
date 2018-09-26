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
  public :: compute_mixture_enthalpy
  public :: compute_species_enthalpy
  public :: compute_mixture_cp
  public :: compute_species_cp
  public :: compute_mixture_dcpdt
  public :: compute_mixture_gamma
  public :: compute_mixture_rgas
  public :: compute_cp_over_temp_integral
  !
  !
  real(wp), parameter :: runivl = 8314.5_wp
  !
  real(wp), parameter :: cp_air_wtmol = 28.9651159_wp
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
  interface compute_mixture_enthalpy
    module procedure compute_mixture_enthalpy_sp, &
                     compute_mixture_enthalpy_1p
  end interface compute_mixture_enthalpy
  !
  interface compute_species_enthalpy
    module procedure compute_species_enthalpy_sp, &
                     compute_species_enthalpy_1p
  end interface compute_species_enthalpy
  !
  interface compute_mixture_cp
    module procedure compute_mixture_cp_sp, &
                     compute_mixture_cp_1p
  end interface compute_mixture_cp
  !
  interface compute_species_cp
    module procedure compute_species_cp_sp, &
                     compute_species_cp_1p
  end interface compute_species_cp
  !
  interface compute_mixture_gamma
    module procedure compute_mixture_gamma_sp, &
                     compute_mixture_gamma_1p, &
                     compute_mixture_gamma_from_temp_1p
  end interface compute_mixture_gamma
  !
  interface compute_mixture_rgas
    module procedure compute_mixture_rgas_sp, &
                     compute_mixture_rgas_1p
  end interface compute_mixture_rgas
  !
  interface compute_cp_over_temp_integral
    module procedure compute_cp_over_temp_integral_sp, &
                     compute_cp_over_temp_integral_1p
  end interface compute_cp_over_temp_integral
  !
contains
!
!###############################################################################
!
subroutine init_thermo_mod(cv,pv,gp)
  !
  !.. Use Statements ..
  use ovar, only : itestcase
  use ovar, only : rgasref,gam,rhoref,pref,tref,cpref
  use ovar, only : aref,uref,machref,velref,velref_nd
  use ovar, only : alpha_aoaref,beta_aoaref
  use ovar, only : pref_nd,rgasref_nd,cpref_nd
  use ovar, only : bc_input
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: cv(:)
  real(wp), intent(inout) :: pv(:)
  real(wp), intent(inout) :: gp(:)
  !
  !.. Local Scalars ..
  integer  :: ics,ierr,m,n,nr
  real(wp) :: mf_sum,href,href_nd,q2,temp,pres_check
#ifdef ALLOCATABLE_MASS_FRAC
  !
  logical(lk) :: error_found
  !
  character(len=50) :: array_name
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: tmp_mass_frac(:)
#endif
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_thermo_mod"
  !
continue
  !
  nr = nme-nmb+1
  !
  if (ncs > 0) then
    !
    allocate ( species_name(1:ncs) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"species_name",1,__LINE__,__FILE__, &
                     ierr,error_message)
    !
    allocate ( mass_frac_ref(1:ncs) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"mass_frac_ref",1,__LINE__,__FILE__, &
                     ierr,error_message)
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
    do m = 1,ncs
      mass_frac_ref(m) = mass_frac_ref(m) / mf_sum
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
    do m = 1,ncs
      pv(m) = mass_frac_ref(m)
    end do
    !
    call compute_mixture_rgas(pv(1:ncs),rgasref,dimensional=true)
    !
  else
    !
    gam = 1.4_wp
    rgasref = runivl / cp_air_wtmol
    !
  end if
  !
  ! Set up the base reference conditions based
  ! on the input flow reference conditions
  !
  ! Check the reference conditions for consistency
  !
  if (rhoref < zero) then
    rhoref = pref/(rgasref*tref)
  else if (pref < zero) then
    pref = rhoref*rgasref*tref
  else if (tref < zero) then
    tref = pref/(rgasref*rhoref)
  else
    pres_check = rhoref*rgasref*tref
    if (abs(pres_check-pref) > eps6) then
     !write (iout,102) pref,pres_check
      write (error_message,100) pref,pres_check
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
  end if
  !
  if (ncs > 0) then
    !
    temperature_switch = temperature_switch/tref
    !
    temp = tref
    !
    call compute_mixture_cp(pv(1:ncs),temp,cpref,dimensional=true)
    !
    gam = rgasref
    call compute_mixture_gamma(cpref,gam)
    !
    call compute_mixture_enthalpy(pv(1:ncs),temp,href,dimensional=true)
    !
  else
    !
    cpref = rgasref * gam/(gam-one)
    href = cpref*tref
    !
  end if
  !
  ! Compute the velocity reference conditions
  !
  aref = sqrt( max( gam*rgasref*tref , rndoff ) )
  !
  uref = machref * aref
  velref(1) = uref * cos(beta_aoaref  * degrees2radians) * &
                     cos(alpha_aoaref * degrees2radians)
  velref(2) = uref * cos(beta_aoaref  * degrees2radians) * &
                     sin(alpha_aoaref * degrees2radians)
  velref(3) = uref * sin(beta_aoaref  * degrees2radians)
  !
  ! Get the nondimensional reference conditions
  !
  ! Nondimensional reference velocities
  velref_nd = velref / aref
  ! Nondimensional reference pressure
  pref_nd = pref / (rhoref*aref*aref)
  ! Nondimensional reference gas constant
  rgasref_nd = rgasref * tref/(aref*aref)
  ! Nondimensional reference specific heat at constant pressure
  cpref_nd = cpref * tref/(aref*aref)
  ! Nondimensional reference static enthalpy
  href_nd = href / (aref*aref)
  !
  ! Set the cv, pv, and gp arrays
  !
  pv(nec)     = rhoref / rhoref
  pv(nmb:nme) = velref_nd(1:nr)
  pv(nee)     = pref_nd
  !
  gp(npres) = pref_nd
  gp(ntemp) = tref / tref
  gp(ngamm) = gam
  !
  q2 = dot_product( pv(nmb:nme) , pv(nmb:nme) )
  !
  cv(nec)     = pv(nec)
  cv(nmb:nme) = pv(nec)*pv(nmb:nme)
  cv(nee)     = pv(nec)*href_nd - pref_nd + half*pv(nec)*q2
  !
  if (ncs > 0) then
    !
    do m = 1,ncs
      cv(m) = pv(m)*pv(nec)
    end do
    !
    gp(nenth) = href_nd
    gp(ncp)   = cpref_nd
    !
  end if
  !
  ! Make sure that the mass_frac component of bc_input makes sense
  !
  if (ncs > 0) then
    !
#ifdef ALLOCATABLE_MASS_FRAC
    error_found = fals
    do n = 1,size(bc_input)
      if (.not. allocated(bc_input(n)%mass_frac)) then
        error_found = true
        if (mypnum == glb_root) write (iout,101) n
      else
        ics = size(bc_input(n)%mass_frac)
        if (ics < ncs) then
          error_found = true
          if (mypnum == glb_root) write (iout,102) n,n,ics,ncs
        else if (ics > ncs) then
          if (any(bc_input(n)%mass_frac(ncs+1:ics) > zero)) then
            error_found = true
            if (mypnum == glb_root) write (iout,103) n,ics,ncs
          else
            tmp_mass_frac = bc_input(n)%mass_frac(1:ncs) ! F2003 AUTO-REALLOC
            bc_input(n)%mass_frac = tmp_mass_frac(1:ncs) ! F2003 AUTO-REALLOC
          end if
       !else if (ics == ncs) then ! this is ideal, do nothing
        end if
      end if
    end do
    flush (iout)
    !
    if (error_found) then
      write (error_message,104)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
#endif
    !
    ! Make sure the values of the mass fractions make sense
    !
    do n = 1,size(bc_input)
      !
      ! Make sure all mass fractions are greater than or equal to zero
      !
      do m = 1,ncs
        bc_input(n)%mass_frac(m) = max( bc_input(n)%mass_frac(m) , zero )
      end do
      !
      ! Adjust the mass fractions if they do not sum to one
      !
      if (abs(sum(bc_input(n)%mass_frac(:))-one) > eps6) then
        if (sum(bc_input(n)%mass_frac(:)) < eps1) then
          ! Set the mass fractions to the reference mass fractions
          ! if the sum is close to zero (i.e., < 0.1)
          do m = 1,ncs
            bc_input(n)%mass_frac(m) = mass_frac_ref(m)
          end do
        else
          ! Make sure the mass fractions sum to one if
          ! the sum is not close to zero (i.e., >= 0.1)
          mf_sum = sum(bc_input(n)%mass_frac(:))
          do m = 1,ncs
            bc_input(n)%mass_frac(m) = bc_input(n)%mass_frac(m) / mf_sum
          end do
        end if
      end if
      !
    end do
    !
  else
    !
    ! Deallocate the mass_frac arrays if this is a calorically perfect gas
    !
    do n = 1,size(bc_input)
#ifdef ALLOCATABLE_MASS_FRAC
      if (allocated(bc_input(n)%mass_frac)) then
        deallocate ( bc_input(n)%mass_frac , stat=ierr , errmsg=error_message )
        write (array_name,1) n
        call alloc_error(pname,array_name,2,__LINE__,__FILE__, &
                         ierr,error_message)
      end if
#else
      bc_input(n)%mass_frac(:) = zero
#endif
    end do
    !
  end if
  !
  ! Format Statements
  !
  1 format ("bc_input(",i0,")%mass_frac")
  !
  100 format ( &
     " Free-stream conditions in the input file are not consistent!",'/&', &
     " Reference pressure given    : ",es14.6,'/&', &
     " Reference pressure computed : ",es14.6,'/&', &
     " Please correct the input file or set one of density, pressure,", &
     " and temperature to a negative value so it is explicitly computed.")
  101 format (/,4x,"ERROR: bc_input(",i0,")%mass_frac is not allocated",/, &
                4x,"        when it needs to be!!!",/)
  102 format (/,4x,"ERROR: bc_input(",i0,")%mass_frac has been allocated",/, &
                4x,"       to a size smaller than the number of chemical",/, &
                4x,"       species present for the current simulation!!!",/, &
                4x,"          size(bc_input(",i0,")%mass_frac = ",i0,/, &
                4x,"               number of chemical species = ",i0,/)
  103 format (/,4x,"ERROR: bc_input(",i0,")%mass_frac has been allocated",/, &
                4x,"       to a size larger than the number of chemical",/, &
                4x,"       species present for the current simulation,",/, &
                4x,"       and the mass fraction is non-zero for one of",/, &
                4x,"       the species larger than the number of chemical",/, &
                4x,"       species.",/, &
                4x,"          size(bc_input(",i0,")%mass_frac = ",i0,/, &
                4x,"               number of chemical species = ",i0,/)
  104 format ("ERROR FOUND WITH THE MASS_FRAC COMPONENT OF BC_INPUT!!!")
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
  integer  :: m,n,num_newton_errors
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
    do m = 1,ncs
      do n = n1,n2
        cv(n,m) = max(cv(n,m),zero)
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
      do m = 1,ncs
        rho = rho + cv(n,m)
      end do
      cv(n,nec) = rho
    end do
    !
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(cv)
    do m = 1,ncs
      do n = n1,n2
        cv(n,m) = cv(n,m)/cv(n,nec)
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
    do m = 1,ncs
      do n = n1,n2
        cv(n,m) = cv(n,m)*cv(n,nec)
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
  ! N2 : ending    solution point
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
subroutine compute_mixture_gamma_sp(n1,n2,cp,gam)
  !
  ! N1  : beginning solution point
  ! N2  : ending    solution point
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
  character(len=*), parameter :: pname = "compute_mixture_gamma_sp"
  !
continue
  !
  !$ACC PARALLEL LOOP PRESENT(cp,gam)
  do n = n1,n2
    gam(n) = cp(n) / (cp(n)-gam(n))
  end do
  !
end subroutine compute_mixture_gamma_sp
!
!###############################################################################
!
pure &
subroutine compute_mixture_gamma_1p(cp,gam)
  !
  ! CP  : mixture specific heat at constant pressure
  ! GAM : in  - mixture gas constant
  !       out - ratio of the specific heats for the mixture
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: cp
  real(wp), intent(inout) :: gam
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_mixture_gamma_1p"
  !
continue
  !
  gam = cp / (cp-gam)
  !
end subroutine compute_mixture_gamma_1p
!
!###############################################################################
!
pure &
subroutine compute_mixture_gamma_from_temp_1p(mass_frac,temp,gam,dimensional)
  !
  ! CP  : mixture specific heat at constant pressure
  ! GAM : in  - mixture gas constant
  !       out - ratio of the specific heats for the mixture
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: mass_frac(:)
  real(wp),    intent(in) :: temp
  real(wp), intent(inout) :: gam
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: dimensional
  !
  !.. Local Scalars ..
  real(wp) :: cp,rgas
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_mixture_gamma_1p"
  !
continue
  !
  call compute_mixture_rgas(mass_frac,rgas,dimensional)
  !
  call compute_mixture_cp(mass_frac,temp,cp,dimensional)
  !
  gam = cp / (cp-rgas)
  !
end subroutine compute_mixture_gamma_from_temp_1p
!
!###############################################################################
!
pure &
subroutine compute_mixture_rgas_sp(n1,n2,mass_frac,rmix)
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
  integer  :: m,n
  real(wp) :: rsum,tda2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_mixture_rgas_sp"
  !
continue
  !
  tda2 = tref/aref/aref
  !
  !$ACC PARALLEL LOOP PRESENT(mass_frac,rmix,rspec)
  do n = n1,n2
    rsum = zero
    !$ACC LOOP SEQ REDUCTION(+:rsum)
    do m = 1,ncs
      rsum = rsum + rspec(m)*mass_frac(n,m)
    end do
    rmix(n) = rsum*tda2
  end do
  !
end subroutine compute_mixture_rgas_sp
!
!###############################################################################
!
pure &
subroutine compute_mixture_rgas_1p(mass_frac,rmix,dimensional)
  !
  !.. Use Statements ..
  use ovar, only : aref,tref
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: mass_frac(:)
  real(wp), intent(inout) :: rmix
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: dimensional
  !
  !.. Local Scalars ..
  logical(lk) :: nondim
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_mixture_rgas_1p"
  !
continue
  !
  nondim = true
  if (present(dimensional)) then
    nondim = (.not. dimensional)
  end if
  !
  rmix = dot_product( rspec(1:ncs) , mass_frac(1:ncs) )
  !
  if (nondim) then
    rmix = rmix * tref/aref/aref
  end if
  !
end subroutine compute_mixture_rgas_1p
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
   1x,5("*"),7x,"NEWTON ITERATION TO SOLVE FOR THE STATIC",t61,5("*"),/, &
   1x,5("*"),4x,"TEMPERATURE OF A GAS MIXTURE DID NOT CONVERGE!",t61,5("*"),/, &
   1x,5("*"),7x,"NUMBER OF UNCONVERGED POINTS = ",i0,t61,5("*"),/, &
   1x,5("*"),7x,"MAXIMUM ERROR IN TEMPERATURE = ",es13.6,t61,5("*"),/, &
   1x,64("*"),/)
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
subroutine compute_mixture_enthalpy_sp(n1,n2,mass_frac,t,h)
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
  integer  :: m,n,ni
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
  do m = 1,ncs
    do n = n1,n2
      if (t(n) > temperature_switch) then
        intrvl(n,m) = 2
      end if
    end do
  end do
  !
  do m = 1,ncs
    !$ACC PARALLEL LOOP PRESENT(mass_frac,intrvl,t,enth_coef,h)
    do n = n1,n2
      !
      ni = intrvl(n,m)
      ts = t(n)*tref
      !
      hi = enth_coef(6,ni,m) + &
           ts*(enth_coef(1,ni,m) + &
               ts*(enth_coef(2,ni,m) + &
                   ts*(enth_coef(3,ni,m) + &
                       ts*(enth_coef(4,ni,m) + &
                           ts*(enth_coef(5,ni,m))))))
      !
      h(n) = h(n) + mass_frac(n,m)*hi*aref2_inv
      !
    end do
  end do
  !
  !$ACC END DATA
  deallocate ( intrvl )
  !
end subroutine compute_mixture_enthalpy_sp
!
!###############################################################################
!
pure &
subroutine compute_mixture_enthalpy_1p(mass_frac,t,h,dimensional)
  !
  !.. Use Statements ..
  use ovar, only : aref,tref
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: mass_frac(:)
  real(wp),    intent(in) :: t
  real(wp), intent(inout) :: h
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: dimensional
  !
  !.. Local Scalars ..
  logical(lk) :: nondim
  integer     :: m,ni
  real(wp)    :: hi,ts,aref2_inv,tdim,tswd
  !
continue
  !
  nondim = true
  if (present(dimensional)) then
    nondim = (.not. dimensional)
  end if
  !
  if (nondim) then
    aref2_inv = one/aref/aref
    tdim = tref
    tswd = one
  else
    aref2_inv = one
    tdim = one
    tswd = tref
  end if
  !
  aref2_inv = one/aref/aref
  !
  h = zero
  !
  do m = 1,ncs
    !
    ni = merge( 1 , 2 , t <= temperature_switch*tswd )
    ts = t*tdim
    !
    hi = enth_coef(6,ni,m) + &
         ts*(enth_coef(1,ni,m) + &
             ts*(enth_coef(2,ni,m) + &
                 ts*(enth_coef(3,ni,m) + &
                     ts*(enth_coef(4,ni,m) + &
                         ts*(enth_coef(5,ni,m))))))
    !
    h = h + mass_frac(m)*hi*aref2_inv
    !
  end do
  !
end subroutine compute_mixture_enthalpy_1p
!
!###############################################################################
!
pure &
subroutine compute_species_enthalpy_sp(n1,n2,t,hsi)
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
  integer  :: m,n,ni
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
  do m = 1,ncs
    do n = n1,n2
      if (t(n) > temperature_switch) then
        intrvl(n,m) = 2
      end if
    end do
  end do
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(intrvl,t,enth_coef,hsi)
  do m = 1,ncs
    do n = n1,n2
      !
      ni = intrvl(n,m)
      ts = t(n)*tref
      !
      hsi(n,m) = enth_coef(6,ni,m) + &
                 ts*(enth_coef(1,ni,m) + &
                     ts*(enth_coef(2,ni,m) + &
                         ts*(enth_coef(3,ni,m) + &
                             ts*(enth_coef(4,ni,m) + &
                                 ts*(enth_coef(5,ni,m))))))
      !
    end do
  end do
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(hsi)
  do m = 1,ncs
    do n = n1,n2
      hsi(n,m) = hsi(n,m)*aref2_inv
    end do
  end do
  !
  !$ACC END DATA
  deallocate ( intrvl )
  !
end subroutine compute_species_enthalpy_sp
!
!###############################################################################
!
pure &
subroutine compute_species_enthalpy_1p(t,hsi)
  !
  !.. Use Statements ..
  use ovar, only : aref,tref
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: t
  real(wp), intent(inout) :: hsi(:)
  !
  !.. Local Scalars ..
  integer  :: m,ni
  real(wp) :: ts,aref2_inv
  !
continue
  !
  aref2_inv = one/aref/aref
  !
  do m = 1,ncs
    !
    ni = merge( 1 , 2 , t <= temperature_switch )
    ts = t*tref
    !
    hsi(m) = enth_coef(6,ni,m) + &
             ts*(enth_coef(1,ni,m) + &
                 ts*(enth_coef(2,ni,m) + &
                     ts*(enth_coef(3,ni,m) + &
                         ts*(enth_coef(4,ni,m) + &
                             ts*(enth_coef(5,ni,m))))))
    !
  end do
  !
  do m = 1,ncs
    hsi(m) = hsi(m)*aref2_inv
  end do
  !
end subroutine compute_species_enthalpy_1p
!
!###############################################################################
!
pure &
subroutine compute_mixture_cp_sp(n1,n2,mass_frac,t,cp)
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
  integer  :: m,n,ni
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
  do m = 1,ncs
    do n = n1,n2
      if (t(n) > temperature_switch) then
        intrvl(n,m) = 2
      end if
    end do
  end do
  !
  do m = 1,ncs
    !$ACC PARALLEL LOOP PRESENT(mass_frac,intrvl,t,cp_coef,cp)
    do n = n1,n2
      !
      ni = intrvl(n,m)
      ts = t(n)*tref
      !
      cpi = cp_coef(1,ni,m) + &
            ts*(cp_coef(2,ni,m) + &
                ts*(cp_coef(3,ni,m) + &
                    ts*(cp_coef(4,ni,m) + &
                        ts*(cp_coef(5,ni,m)))))
      !
      cp(n) = cp(n) + mass_frac(n,m)*cpi*tda2
      !
    end do
  end do
  !
  !$ACC END DATA
  deallocate ( intrvl )
  !
end subroutine compute_mixture_cp_sp
!
!###############################################################################
!
pure &
subroutine compute_mixture_cp_1p(mass_frac,t,cp,dimensional)
  !
  !.. Use Statements ..
  use ovar, only : aref,tref
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: mass_frac(:)
  real(wp),    intent(in) :: t
  real(wp), intent(inout) :: cp
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: dimensional
  !
  !.. Local Scalars ..
  logical(lk) :: nondim
  integer     :: m,ni
  real(wp)    :: cpi,ts,tda2,tdim,tswd
  !
continue
  !
  nondim = true
  if (present(dimensional)) then
    nondim = (.not. dimensional)
  end if
  !
  if (nondim) then
    tda2 = tref/aref/aref
    tdim = tref
    tswd = one
  else
    tda2 = one
    tdim = one
    tswd = tref
  end if
  !
  cp = zero
  !
  do m = 1,ncs
    !
    ni = merge( 1 , 2 , t <= temperature_switch*tswd )
    ts = t*tdim
    !
    cpi = cp_coef(1,ni,m) + &
          ts*(cp_coef(2,ni,m) + &
              ts*(cp_coef(3,ni,m) + &
                  ts*(cp_coef(4,ni,m) + &
                      ts*(cp_coef(5,ni,m)))))
    !
    cp = cp + mass_frac(m)*cpi*tda2
    !
  end do
  !
end subroutine compute_mixture_cp_1p
!
!###############################################################################
!
pure &
subroutine compute_species_cp_sp(n1,n2,t,cpi)
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
  integer  :: m,n,ni
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
  do m = 1,ncs
    do n = n1,n2
      if (t(n) > temperature_switch) then
        intrvl(n,m) = 2
      end if
    end do
  end do
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(intrvl,t,cp_coef,cpi)
  do m = 1,ncs
    do n = n1,n2
      !
      ni = intrvl(n,m)
      ts = t(n)*tref
      !
      cpi(n,m) = cp_coef(1,ni,m) + &
                 ts*(cp_coef(2,ni,m) + &
                     ts*(cp_coef(3,ni,m) + &
                         ts*(cp_coef(4,ni,m) + &
                             ts*(cp_coef(5,ni,m)))))
      !
    end do
  end do
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(cpi)
  do m = 1,ncs
    do n = n1,n2
      cpi(n,m) = cpi(n,m)*tda2
    end do
  end do
  !
  !$ACC END DATA
  deallocate ( intrvl )
  !
end subroutine compute_species_cp_sp
!
!###############################################################################
!
pure &
subroutine compute_species_cp_1p(t,cpi)
  !
  !.. Use Statements ..
  use ovar, only : aref,tref
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: t
  real(wp), intent(inout) :: cpi(:)
  !
  !.. Local Scalars ..
  integer  :: m,ni
  real(wp) :: ts,tda2
  !
continue
  !
  tda2 = tref/aref/aref
  !
  do m = 1,ncs
    !
    ni = merge( 1 , 2 , t <= temperature_switch )
    ts = t*tref
    !
    cpi(m) = cp_coef(1,ni,m) + &
             ts*(cp_coef(2,ni,m) + &
                 ts*(cp_coef(3,ni,m) + &
                     ts*(cp_coef(4,ni,m) + &
                         ts*(cp_coef(5,ni,m)))))
    !
  end do
  !
  do m = 1,ncs
    cpi(m) = cpi(m)*tda2
  end do
  !
end subroutine compute_species_cp_1p
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
  !     \                         d [Cp(m)]
  !      \     MassFraction(m) * -----------
  !      /                           dT
  !     /
  !     ------
  !       m=1
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
  integer  :: m,n,ni
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
  do m = 1,ncs
    do n = n1,n2
      if (t(n) > temperature_switch) then
        intrvl(n,m) = 2
      end if
    end do
  end do
  !
  do m = 1,ncs
    !$ACC PARALLEL LOOP PRESENT(mass_frac,intrvl,t,cp_coef,dcpdt)
    do n = n1,n2
      !
      ni = intrvl(n,m)
      ts = t(n)*tref
      !
      dcpidt = cp_coef(2,ni,m) + &
               ts*(cp_coef(3,ni,m)*two + &
                   ts*(cp_coef(4,ni,m)*three + &
                       ts*(cp_coef(5,ni,m)*four)))
      !
      dcpdt(n) = dcpdt(n) + mass_frac(n,m)*dcpidt*aref2_inv
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
pure &
subroutine compute_cp_over_temp_integral_sp(n1,n2,mass_frac, &
                                            return_value,temp,ttot)
  !
  ! Routine to compute the integral required for calculating total pressure
  ! for a mixture of gases in thermal equilibrium using the Gordon and
  ! McBride curve fits
  !
  ! NPTS : no. of points at which to calculate
  ! M  : species no.
  ! MASS_FRAC  : species mass fraction
  ! TTOT : total temperature
  ! TEMP : static temperature
  ! RETURN_VALUE : Final intergral of Cp/T summed over all species
  !
  !.. Use Statements ..
  use ovar, only : aref,tref
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: mass_frac(:,:)
  real(wp),    intent(in) :: ttot(:)
  real(wp),    intent(in) :: temp(:)
  real(wp), intent(inout) :: return_value(:)
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: tbl,tel,ccl,cel,cbl,al,bl,cl,dl,el
  real(wp) :: tbh,teh,cch,ceh,cbh,ah,bh,ch,dh,eh
  real(wp) :: val,tda2
  !
  !.. Local Parameters ..
  real(wp), parameter :: three2 = three/two
  real(wp), parameter :: four3  = four/three
  real(wp), parameter :: five4  = five/four
  !
continue
  !
  tda2 = tref/aref/aref
  !
  !$ACC PARALLEL LOOP PRESENT(return_value)
  do n = n1,n2
    return_value(n) = zero
  end do
  !
!$define LOOP_PTS_THEN_SPECIES
#ifdef LOOP_PTS_THEN_SPECIES
  !$ACC PARALLEL LOOP PRESENT(ttot,temp,mass_frac,return_value,enth_coef)
  do n = n1,n2
    !
    val = zero
    !
    if (ttot(n) <= temperature_switch) then
      !
      ! Static and total temperature are both smaller than temperature switch
      ! NOTE: We need to dimensionalize the temperatures because
      !       the curve fits require dimensional temperature.
      !
      tbl = temp(n)*tref
      tel = ttot(n)*tref
      !
      !$ACC LOOP SEQ
      do m = 1,ncs
        !
        ccl = log(tel/tbl)*enth_coef(1,1,m)
        !
        cel = tel*(enth_coef(2,1,m)*two + &
                  tel*(enth_coef(3,1,m)*three2 + &
                      tel*(enth_coef(4,1,m)*four3 + &
                          tel*(enth_coef(5,1,m)*five4))))
        !
        cbl = tbl*(enth_coef(2,1,m)*two + &
                  tbl*(enth_coef(3,1,m)*three2 + &
                      tbl*(enth_coef(4,1,m)*four3 + &
                          tbl*(enth_coef(5,1,m)*five4))))
        !
        val = val + mass_frac(n,m)*(ccl + cel - cbl)
        !
      end do
      !
    else if (temp(n) > temperature_switch) then
      !
      ! Static and total temperature are both larger than temperature switch
      ! NOTE: We need to dimensionalize the temperatures because
      !       the curve fits require dimensional temperature.
      !
      tbh = temp(n)*tref
      teh = ttot(n)*tref
      !
      !$ACC LOOP SEQ
      do m = 1,ncs
        !
        cch = log(teh/tbh)*enth_coef(1,2,m)
        !
        ceh = teh*(enth_coef(2,2,m)*two + &
                  teh*(enth_coef(3,2,m)*three2 + &
                      teh*(enth_coef(4,2,m)*four3 + &
                          teh*(enth_coef(5,2,m)*five4))))
        !
        cbh = tbh*(enth_coef(2,2,m)*two + &
                  tbh*(enth_coef(3,2,m)*three2 + &
                      tbh*(enth_coef(4,2,m)*four3 + &
                          tbh*(enth_coef(5,2,m)*five4))))
        !
        val = val + mass_frac(n,m)*(cch + ceh - cbh)
        !
      end do
      !
    else
      !
      ! Static temperature is smaller than temperature switch and
      ! total temperature is larger than temperature switch
      ! NOTE: We need to dimensionalize the temperatures because
      !       the curve fits require dimensional temperature.
      !
      tbl = temp(n)*tref
      tel = temperature_switch*tref
      !
      tbh = temperature_switch*tref
      teh = ttot(n)*tref
      !
      !$ACC LOOP SEQ
      do m = 1,ncs
        !
        ccl = log(tel/tbl)*enth_coef(1,1,m)
        !
        cel = tel*(enth_coef(2,1,m)*two + &
                  tel*(enth_coef(3,1,m)*three2 + &
                      tel*(enth_coef(4,1,m)*four3 + &
                          tel*(enth_coef(5,1,m)*five4))))
        !
        cbl = tbl*(enth_coef(2,1,m)*two + &
                  tbl*(enth_coef(3,1,m)*three2 + &
                      tbl*(enth_coef(4,1,m)*four3 + &
                          tbl*(enth_coef(5,1,m)*five4))))
        !
        cch = log(teh/tbh)*enth_coef(1,2,m)
        !
        ceh = teh*(enth_coef(2,2,m)*two + &
                  teh*(enth_coef(3,2,m)*three2 + &
                      teh*(enth_coef(4,2,m)*four3 + &
                          teh*(enth_coef(5,2,m)*five4))))
        !
        cbh = tbh*(enth_coef(2,2,m)*two + &
                  tbh*(enth_coef(3,2,m)*three2 + &
                      tbh*(enth_coef(4,2,m)*four3 + &
                          tbh*(enth_coef(5,2,m)*five4))))
        !
        val = val + mass_frac(n,m)*( (ccl+cch) + (cel+ceh) - (cbl+cbh) )
        !
      end do
      !
    end if
    !
    return_value(n) = val*tda2
    !
  end do
#else
  do m = 1,ncs
    !
    al = enth_coef(1,1,m)
    bl = enth_coef(2,1,m)*two
    cl = enth_coef(3,1,m)*three2
    dl = enth_coef(4,1,m)*four3
    el = enth_coef(5,1,m)*five4
    !
    ah = enth_coef(1,2,m)
    bh = enth_coef(2,2,m)*two
    ch = enth_coef(3,2,m)*three2
    dh = enth_coef(4,2,m)*four3
    eh = enth_coef(5,2,m)*five4
    !
    !$ACC PARALLEL LOOP PRESENT(ttot,temp,mass_frac,return_value)
    do n = n1,n2
      !
      if (ttot(n) <= temperature_switch) then
        !
        ! Static and total temperature are both smaller than temperature switch
        ! NOTE: We need to dimensionalize the temperatures because
        !       the curve fits require dimensional temperature.
        !
        tbl = temp(n)*tref
        tel = ttot(n)*tref
        !
        ccl = log(tel/tbl)*al
        !
        cel = tel*(bl + tel*(cl + tel*(dl + tel*el)))
        cbl = tbl*(bl + tbl*(cl + tbl*(dl + tbl*el)))
        !
        return_value(n) = return_value(n) + mass_frac(n,m)*(ccl + cel - cbl)
        !
      else if (temp(n) > temperature_switch) then
        !
        ! Static and total temperature are both larger than temperature switch
        ! NOTE: We need to dimensionalize the temperatures because
        !       the curve fits require dimensional temperature.
        !
        tbh = temp(n)*tref
        teh = ttot(n)*tref
        !
        cch = log(teh/tbh)*ah
        !
        ceh = teh*(bh + teh*(ch + teh*(dh + teh*eh)))
        cbh = tbh*(bh + tbh*(ch + tbh*(dh + tbh*eh)))
        !
        return_value(n) = return_value(n) + mass_frac(n,m)*(cch + ceh - cbh)
        !
      else
        !
        ! Static temperature is smaller than temperature switch and
        ! total temperature is larger than temperature switch
        ! NOTE: We need to dimensionalize the temperatures because
        !       the curve fits require dimensional temperature.
        !
        tbl = temp(n)*tref
        tel = temperature_switch*tref
        !
        ccl = log(tel/tbl)*al
        !
        cel = tel*(bl + tel*(cl + tel*(dl + tel*el)))
        cbl = tbl*(bl + tbl*(cl + tbl*(dl + tbl*el)))
        !
        tbh = temperature_switch*tref
        teh = ttot(n)*tref
        !
        cch = log(teh/tbh)*ah
        !
        ceh = teh*(bh + teh*(ch + teh*(dh + teh*eh)))
        cbh = tbh*(bh + tbh*(ch + tbh*(dh + tbh*eh)))
        !
        return_value(n) = return_value(n) + mass_frac(n,m)*( (ccl+cel-cbl) + &
                                                             (cch+ceh-cbh) )
        !
      end if
      !
    end do
    !
  end do
  !
  ! Non-dimensionalize return_value
  !
  !$ACC PARALLEL LOOP PRESENT(return_value)
  do n = n1,n2
    return_value(n) = return_value(n) * tda2
  end do
#endif
  !
  ! Finished computing the integral required for calculating total pressure
  ! for a mixture of gases in thermal equilibrium using the Gordon and
  ! McBride curve fits
  !
end subroutine compute_cp_over_temp_integral_sp
!
!###############################################################################
!
pure &
subroutine compute_cp_over_temp_integral_1p(mass_frac,return_value,temp,ttot)
  !
  ! Routine to compute the integral required for calculating total pressure
  ! for a mixture of gases in thermal equilibrium using the Gordon and
  ! McBride curve fits
  !
  ! NPTS : no. of points at which to calculate
  ! M  : species no.
  ! MASS_FRAC  : species mass fraction
  ! TTOT : total temperature
  ! TEMP : static temperature
  ! RETURN_VALUE : Final intergral of Cp/T summed over all species
  !
  !.. Use Statements ..
  use ovar, only : aref,tref
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: mass_frac(:)
  real(wp),    intent(in) :: ttot
  real(wp),    intent(in) :: temp
  real(wp), intent(inout) :: return_value
  !
  !.. Local Scalars ..
  integer  :: m
  real(wp) :: tbl,tel,ccl,cel,cbl,al,bl,cl,dl,el
  real(wp) :: tbh,teh,cch,ceh,cbh,ah,bh,ch,dh,eh
  real(wp) :: tda2
  !
  !.. Local Parameters ..
  real(wp), parameter :: three2 = three/two
  real(wp), parameter :: four3  = four/three
  real(wp), parameter :: five4  = five/four
  !
continue
  !
  tda2 = tref/aref/aref
  !
  return_value = zero
  !
  do m = 1,ncs
    !
    al = enth_coef(1,1,m)
    bl = enth_coef(2,1,m)*two
    cl = enth_coef(3,1,m)*three2
    dl = enth_coef(4,1,m)*four3
    el = enth_coef(5,1,m)*five4
    !
    ah = enth_coef(1,2,m)
    bh = enth_coef(2,2,m)*two
    ch = enth_coef(3,2,m)*three2
    dh = enth_coef(4,2,m)*four3
    eh = enth_coef(5,2,m)*five4
    !
    if (ttot <= temperature_switch) then
      !
      ! Static and total temperature are both smaller than temperature switch
      ! NOTE: We need to dimensionalize the temperatures because
      !       the curve fits require dimensional temperature.
      !
      tbl = temp*tref
      tel = ttot*tref
      !
      ccl = log(tel/tbl)*al
      !
      cel = tel*(bl + tel*(cl + tel*(dl + tel*el)))
      cbl = tbl*(bl + tbl*(cl + tbl*(dl + tbl*el)))
      !
      return_value = return_value + mass_frac(m)*(ccl + cel - cbl)
      !
    else if (temp > temperature_switch) then
      !
      ! Static and total temperature are both larger than temperature switch
      ! NOTE: We need to dimensionalize the temperatures because
      !       the curve fits require dimensional temperature.
      !
      tbh = temp*tref
      teh = ttot*tref
      !
      cch = log(teh/tbh)*ah
      !
      ceh = teh*(bh + teh*(ch + teh*(dh + teh*eh)))
      cbh = tbh*(bh + tbh*(ch + tbh*(dh + tbh*eh)))
      !
      return_value = return_value + mass_frac(m)*(cch + ceh - cbh)
      !
    else
      !
      ! Static temperature is smaller than temperature switch and
      ! total temperature is larger than temperature switch
      ! NOTE: We need to dimensionalize the temperatures because
      !       the curve fits require dimensional temperature.
      !
      tbl = temp*tref
      tel = temperature_switch*tref
      !
      ccl = log(tel/tbl)*al
      !
      cel = tel*(bl + tel*(cl + tel*(dl + tel*el)))
      cbl = tbl*(bl + tbl*(cl + tbl*(dl + tbl*el)))
      !
      tbh = temperature_switch*tref
      teh = ttot*tref
      !
      cch = log(teh/tbh)*ah
      !
      ceh = teh*(bh + teh*(ch + teh*(dh + teh*eh)))
      cbh = tbh*(bh + tbh*(ch + tbh*(dh + tbh*eh)))
      !
      return_value = return_value + mass_frac(m)*((ccl+cel-cbl) + (cch+ceh-cbh))
      !
    end if
    !
  end do
  !
  ! Non-dimensionalize return_value
  !
  return_value = return_value * tda2
  !
  ! Finished computing the integral required for calculating total pressure
  ! for a mixture of gases in thermal equilibrium using the Gordon and
  ! McBride curve fits
  !
end subroutine compute_cp_over_temp_integral_1p
!
!###############################################################################
!
subroutine hardwired_NCC_coefficients
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
  integer  :: i,j,m,ncoef
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
  do m = 1,ncs
    !
    species = uppercase( trim(adjustl(species_name(m))) )
    !
    select case (trim(adjustl(species)))
      !
      case ("CO2")
        !
        wtmol(m) = CO2_wtmol
        ncc_cp(1:6,1:2,m) = CO2_cp(1:6,1:2)
        !
      case ("CO")
        !
        wtmol(m) = CO_wtmol
        ncc_cp(1:6,1:2,m) = CO_cp(1:6,1:2)
        !
      case ("H2O")
        !
        wtmol(m) = H2O_wtmol
        ncc_cp(1:6,1:2,m) = H2O_cp(1:6,1:2)
        !
      case ("CH4")
        !
        wtmol(m) = CH4_wtmol
        ncc_cp(1:6,1:2,m) = CH4_cp(1:6,1:2)
        !
      case ("N2")
        !
        wtmol(m) = N2_wtmol
        ncc_cp(1:6,1:2,m) = N2_cp(1:6,1:2)
        !
      case ("O2")
        !
        wtmol(m) = O2_wtmol
        ncc_cp(1:6,1:2,m) = O2_cp(1:6,1:2)
        !
      case ("H")
        !
        wtmol(m) = H_wtmol
        ncc_cp(1:6,1:2,m) = H_cp(1:6,1:2)
        !
      case ("H2")
        !
        wtmol(m) = H2_wtmol
        ncc_cp(1:6,1:2,m) = H2_cp(1:6,1:2)
        !
      case ("O")
        !
        wtmol(m) = O_wtmol
        ncc_cp(1:6,1:2,m) = O_cp(1:6,1:2)
        !
      case ("OH")
        !
        wtmol(m) = OH_wtmol
        ncc_cp(1:6,1:2,m) = OH_cp(1:6,1:2)
        !
    end select
    !
  end do
  !
  do m = 1,ncs
    rspec(m) = runivl / wtmol(m)
  end do
  !
  do m = 1,ncs
    do i = 1,size(ncc_cp,dim=2)
      do j = 1,ncpc
        cp_coef(j,i,m) = ncc_cp(j,i,m) * rspec(m)
      end do
    end do
  end do
  !
  do m = 1,ncs
    do i = 1,size(ncc_cp,dim=2)
      do j = 1,ncpc
        ncoef = j + min(6-(ncpc+1),0)
        denom = max(real(ncoef,kind=wp),one)
        enth_coef(j,i,m) = ncc_cp(j,i,m) * rspec(m) / denom
      end do
      enth_coef(ncpc+1,i,m) = ncc_cp(ncpc+1,i,m) * rspec(m)
    end do
  end do
  !
  temperature_switch = NCC_interval_switch
  !
  deallocate ( ncc_cp )
  !
end subroutine hardwired_NCC_coefficients
!
!###############################################################################
!
subroutine hardwired_VULCAN_coefficients
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
  integer  :: i,j,m,ncoef
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
  do m = 1,ncs
    !
    species = uppercase( trim(adjustl(species_name(m))) )
    !
    select case (trim(adjustl(species)))
      !
      case ("CO2")
        !
        wtmol(m) = CO2_wtmol
        vul_cp(1:7,1:2,m) = CO2_cp(1:7,1:2)
        !
      case ("CO")
        !
        wtmol(m) = CO_wtmol
        vul_cp(1:7,1:2,m) = CO_cp(1:7,1:2)
        !
      case ("H2O")
        !
        wtmol(m) = H2O_wtmol
        vul_cp(1:7,1:2,m) = H2O_cp(1:7,1:2)
        !
      case ("CH4")
        !
        wtmol(m) = CH4_wtmol
        vul_cp(1:7,1:2,m) = CH4_cp(1:7,1:2)
        !
      case ("N2")
        !
        wtmol(m) = N2_wtmol
        vul_cp(1:7,1:2,m) = N2_cp(1:7,1:2)
        !
      case ("O2")
        !
        wtmol(m) = O2_wtmol
        vul_cp(1:7,1:2,m) = O2_cp(1:7,1:2)
        !
      case ("H")
        !
        wtmol(m) = H_wtmol
        vul_cp(1:7,1:2,m) = H_cp(1:7,1:2)
        !
      case ("H2")
        !
        wtmol(m) = H2_wtmol
        vul_cp(1:7,1:2,m) = H2_cp(1:7,1:2)
        !
      case ("O")
        !
        wtmol(m) = O_wtmol
        vul_cp(1:7,1:2,m) = O_cp(1:7,1:2)
        !
      case ("OH")
        !
        wtmol(m) = OH_wtmol
        vul_cp(1:7,1:2,m) = OH_cp(1:7,1:2)
        !
    end select
    !
  end do
  !
  do m = 1,ncs
    rspec(m) = runivl / wtmol(m)
  end do
  !
  do m = 1,ncs
    do i = 1,size(vul_cp,dim=2)
      do j = 1,ncpc
        cp_coef(j,i,m) = vul_cp(j,i,m) * rspec(m)
      end do
    end do
  end do
  !
  do m = 1,ncs
    do i = 1,size(vul_cp,dim=2)
      do j = 1,ncpc
        ncoef = j + min(6-(ncpc+1),0)
        denom = max(real(ncoef,kind=wp),one)
        enth_coef(j,i,m) = vul_cp(j,i,m) * rspec(m) / denom
      end do
      enth_coef(ncpc+1,i,m) = vul_cp(ncpc+1,i,m) * rspec(m)
    end do
  end do
  !
  temperature_switch = vul_interval_switch
  !
  deallocate ( vul_cp )
  !
end subroutine hardwired_VULCAN_coefficients
!
!###############################################################################
!
end module thermo_mod
