module thermo_mod
  !
  !.. Use Statements ..
  use module_kind_types
  !
  use eqn_idx, only : nec,nmb,nme,nmx,nmy,nmz,nee,nq
  use eqn_idx, only : ncs,ncr
  use eqn_idx, only : npres,ntemp,ncp,ngamm,nenth,ng
  !
  implicit none
  !
  private
  !
  public :: init_thermo_mod
  public :: update_thermo_properties
  public :: update_mixture_density
  public :: compute_mixture_enthalpy
  public :: compute_species_enthalpy
  public :: compute_species_energy
  public :: compute_mixture_cp
  public :: compute_species_cp
  public :: compute_mixture_dcpdt
  public :: compute_species_dcpdt
  public :: compute_mixture_gamma
  public :: compute_mixture_rgas
  public :: compute_cp_over_temp_integral
  public :: compute_gibbs_free_energy
  public :: compute_dgbdt
  !
  !
  ! RUNIVL : Universal gas constant
  !          NOTE: This is normally given in J/(K*moles) [i.e., (kg*m^2)/(s^2)],
  !                but we are using mJoules/(K*mol) [i.e., (g*m^2)/(s^2)].
  !                This is because this is only really used with the molecular
  !                weights of gas species which are given in [g/mole].
  !
  real(wp), public, parameter :: runivl = 8314.5_wp ! mJoules/(K*moles)
  !
  ! CP_AIR_WTMOL : molecular weight for calorically perfect air (g/mole)
  real(wp), parameter :: cp_air_wtmol = 28.9651159_wp
  !
  ! RGAS_CP : Specific gas constant for calorically perfect air using
  !           cp_air_wtmol for the molecular weight of air.
  !           rgas_cpg = 287.05219 J/(kg*K)
  real(wp), public, parameter :: rgas_cpg = runivl/cp_air_wtmol
  !
  integer, protected :: ncpc = 5
  !
  !
  !
  real(wp), public, protected, save, allocatable :: wtmol(:)
  !
  character(len=16), public, protected, save, allocatable :: species_name(:)
  !
  integer, public, protected, save, allocatable :: inert_species_map(:)
  !
  real(wp), public, protected, save, allocatable :: rspec(:)
  !
  real(wp), protected, save, allocatable :: cp_coef(:,:,:)
  real(wp), protected, save, allocatable :: enth_coef(:,:,:)
  real(wp), protected, save, allocatable :: gbfe_coef(:,:,:)
  real(wp), protected, save, allocatable :: dgbfe_coef(:,:,:)
  !
  real(wp), public, protected, save :: temperature_high
  real(wp), public, protected, save :: temperature_low
  !
  real(wp), protected, save, allocatable :: mass_frac_ref(:)
  !
  integer, parameter :: Use_VULCAN_Coefficients = 1
  integer, parameter :: Use_NCC_Coefficients = 2
  !
 !integer, parameter :: Thermo_Coefficients = Use_NCC_Coefficients
  integer, parameter :: Thermo_Coefficients = Use_Vulcan_Coefficients
  !
  logical(lk), save :: file_is_opened = fals
  integer, save :: last_itcur_output = 0
  integer, save :: last_rk_stage_output = 0
  integer, save :: io = 3682
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
  interface compute_species_energy
    module procedure compute_species_energy_sp
  end interface compute_species_energy
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
  use ovar, only : bc_input,disable_chemical_reactions
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: cv(:)
  real(wp), intent(inout) :: pv(:)
  real(wp), intent(inout) :: gp(:)
  !
  !.. Local Scalars ..
  integer  :: ics,jcs,ierr,m,n,nr
  real(wp) :: mf_sum,href,href_nd,q2,temp,pres_check
  real(wp) :: tdim,aref2_inv
  !
  character(len=10) :: this_species,input_species
  !
#ifdef ALLOCATABLE_MASS_FRAC
  logical(lk) :: error_found
  !
  character(len=50) :: array_name
#endif
  !
  !.. Local Allocatable Arrays ..
  real(wp),         allocatable :: tmp_mass_frac(:)
  character(len=:), allocatable :: species_list(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_thermo_mod"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = nme-nmb+1
  !
  if (ncs > 0) then
    !
    if (.not. allocated(species_name)) then
      allocate ( species_name(1:ncs) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"species_name",1,__LINE__,__FILE__, &
                       ierr,error_message)
    end if
    !
    if (.not. allocated(mass_frac_ref)) then
      allocate ( mass_frac_ref(1:ncs) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"mass_frac_ref",1,__LINE__,__FILE__, &
                       ierr,error_message)
    end if
    !
    if (abs(itestcase) == CH4_3D_Diffusion_Flame .or. &
        abs(itestcase) == CH4_Flame_Speed_Test) then
      !
      if (disable_chemical_reactions) then
        !
        species_name(1) = "CH4"
        species_name(2) = "N2"
        species_name(3) = "O2"
        !
        if (abs(itestcase) == CH4_Flame_Speed_Test) then
          mass_frac_ref(1) = 0.055_wp  ! CH4
          mass_frac_ref(2) = 0.724_wp  ! N2
          mass_frac_ref(3) = 0.221_wp  ! O2
        else
          mass_frac_ref(1) = zero      ! CH4
          mass_frac_ref(2) = 0.7686_wp ! N2
          mass_frac_ref(3) = 0.2314_wp ! O2
        end if
        !
        inert_species_map = [1,2,3] ! F2003 AUTO REALLOCATION
        !
      else
        !
        species_name(1) = "CO2"
        species_name(2) = "CO"
       !species_name(3) = "H2"
        species_name(3) = "H2O"
        species_name(4) = "CH4"
        species_name(5) = "N2"
        species_name(6) = "O2"
        !
        mass_frac_ref(1) = zero      ! CO2
        mass_frac_ref(2) = zero      ! CO
        mass_frac_ref(3) = zero      ! H2O
        if (abs(itestcase) == CH4_Flame_Speed_Test) then
          mass_frac_ref(4) = 0.055_wp  ! CH4
          mass_frac_ref(5) = 0.724_wp  ! N2
          mass_frac_ref(6) = 0.221_wp  ! O2
        else
          mass_frac_ref(4) = zero      ! CH4
          mass_frac_ref(5) = 0.7686_wp ! N2
          mass_frac_ref(6) = 0.2314_wp ! O2
        end if
        !
        ! Create the map from inert species indices to
        ! reacting species indices. This is needed if a
        ! previous simulation was run with reactions turned off, 
        ! and the new simulation is restarting with reactions
        ! turned on.
        !
        inert_species_map = [4,5,6] ! F2003 AUTO REALLOCATION
        !
      end if
      !
    else if (itestcase == H2_2D_Diffusion_Flame) then
      !
      if (disable_chemical_reactions) then
        !
        species_name(1) = "H2"
        species_name(2) = "O2"
        species_name(3) = "N2"
        !
        mass_frac_ref(1) = zero      ! H2
        mass_frac_ref(2) = 0.2314_wp ! O2
        mass_frac_ref(3) = 0.7686_wp ! N2
        !
        inert_species_map = [1,2,3] ! F2003 AUTO REALLOCATION
        !
      else
        !
        species_name(1) = "H2"
        species_name(2) = "O2"
        species_name(3) = "H2O"
        species_name(4) = "OH"
        species_name(5) = "H"
        species_name(6) = "O"
        species_name(7) = "N2"
        !
        mass_frac_ref(1) = zero      ! H2
        mass_frac_ref(2) = 0.2314_wp ! O2
        mass_frac_ref(3) = zero      ! H2O
        mass_frac_ref(4) = zero      ! OH
        mass_frac_ref(5) = zero      ! H
        mass_frac_ref(6) = zero      ! O
        mass_frac_ref(7) = 0.7686_wp ! N2
        !
        ! Create the map from inert species indices to
        ! reacting species indices. This is needed if a
        ! previous simulation was run with reactions turned off, 
        ! and the new simulation is restarting with reactions
        ! turned on.
        !
        inert_species_map = [1,2,7] ! F2003 AUTO REALLOCATION
        !
      end if
      !
    else if (itestcase == Compression_Corner) then
      !
      if (disable_chemical_reactions) then
        !
        species_name(1) = "O2"
        species_name(2) = "H2"
        species_name(3) = "N2"
        !
        mass_frac_ref(1) = 0.24_wp  ! O2
        mass_frac_ref(2) = 0.005_wp ! H2
        mass_frac_ref(3) = 0.755_wp ! N2
        !
        inert_species_map = [1,2,3] ! F2003 AUTO REALLOCATION
        !
      else
        !
        species_name(1) = "O2"
        species_name(2) = "H2"
        species_name(3) = "H2O"
        species_name(4) = "OH"
        species_name(5) = "O"
        species_name(6) = "H"
        species_name(7) = "N2"
        !
        mass_frac_ref(1) = 0.24_wp  ! O2
        mass_frac_ref(2) = 0.005_wp ! H2
        mass_frac_ref(3) = zero     ! H2O
        mass_frac_ref(4) = zero     ! OH
        mass_frac_ref(5) = zero     ! O
        mass_frac_ref(6) = zero     ! H
        mass_frac_ref(7) = 0.755_wp ! N2
        !
        ! Create the map from inert species indices to
        ! reacting species indices. This is needed if a
        ! previous simulation was run with reactions turned off, 
        ! and the new simulation is restarting with reactions
        ! turned on.
        !
        inert_species_map = [1,2,7] ! F2003 AUTO REALLOCATION
        !
      end if
      !
    else if (itestcase == Reacting_Backstep) then
      !
      species_name(1) = "N2"
      species_name(2) = "O2"
      !
      mass_frac_ref(1) = 0.7686 ! N2
      mass_frac_ref(2) = 0.2314 ! O2
      !
      inert_species_map = [1,2] ! F2003 AUTO REALLOCATION
      !
    else if (itestcase == N2_O2_Diffusion) then
      !
      species_name(1) = "N2"
      species_name(2) = "O2"
      !
      mass_frac_ref(1) = 1.0_wp ! N2
      mass_frac_ref(2) = 0.0_wp ! O2
      !
      inert_species_map = [1,2] ! F2003 AUTO REALLOCATION
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
    if (.not. allocated(wtmol)) then
      allocate ( wtmol(1:ncs) , source=zero , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"wtmol",1,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    if (.not. allocated(rspec)) then
      allocate ( rspec(1:ncs) , source=zero , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"rspec",1,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    ncpc = 5
    !
    if (.not. allocated(cp_coef)) then
      allocate ( cp_coef(1:ncpc,0:2,1:ncs) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"cp_coef",1,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    if (.not. allocated(enth_coef)) then
      allocate ( enth_coef(1:ncpc+1,0:2,1:ncs) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"enth_coef",1,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    if (.not. allocated(gbfe_coef)) then
      allocate ( gbfe_coef(1:ncpc+2,0:2,1:ncs) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"gbfe_coef",1,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    if (.not. allocated(dgbfe_coef)) then
      allocate ( dgbfe_coef(1:ncpc+1,0:2,1:ncs) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"dgbfe_coef",1,__LINE__,__FILE__, &
                       ierr,error_message)
    end if
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
    ! NOTE: The rspec array has not been non-dimensionalized yet, so the
    !       following call to compute_mixture_rgas will return a dimensional
    !       value for rgasref.
    !
    call compute_mixture_rgas(pv(1:ncs),rgasref)
    !
  else
    !
    if (gam <= zero) then
      gam = 1.4_wp
    end if
    !
    if (rgasref <= zero) then
      rgasref = rgas_cpg
    end if
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
    temp = tref
    !
    call compute_mixture_cpref(pv(1:ncs),temp,cpref)
    !
    gam = rgasref
    call compute_mixture_gamma(cpref,gam)
    !
    call compute_mixture_href(pv(1:ncs),temp,href)
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
    ! Non-dimensionalize temperature_high, temperature_low, rspec,
    !                    cp_coef, enth_coef, gbfe_coef, and dgbfe_coef
    !
    temperature_low  = temperature_low /tref
    temperature_high = temperature_high/tref
    !
    aref2_inv = one/aref/aref
    !
    ! Non-dimensionalize rspec
    !
    do ics = 1,ncs
      rspec(ics) = rspec(ics) * tref*aref2_inv
    end do
    !
    ! Non-dimensionalize cp_coef
    !
    do ics = 1,ncs
      do n = lbound(cp_coef,dim=2),ubound(cp_coef,dim=2)
        tdim = one
        do m = 1,ncpc
          tdim = tdim*tref
          cp_coef(m,n,ics) = cp_coef(m,n,ics) * tdim * aref2_inv
        end do
      end do
    end do
    !
    ! Non-dimensionalize enth_coef
    !
    do ics = 1,ncs
      do n = lbound(enth_coef,dim=2),ubound(enth_coef,dim=2)
        !
        enth_coef(ncpc+1,n,ics) = enth_coef(ncpc+1,n,ics) * aref2_inv
        !
        tdim = one
        do m = 1,ncpc
          tdim = tdim*tref
          enth_coef(m,n,ics) = enth_coef(m,n,ics) * tdim * aref2_inv
        end do
        !
      end do
    end do
    !
    ! Non-dimensionalize gbfe_coef
    !
    if (ncpc == 5) then
      !
      do ics = 1,ncs
        do n = lbound(gbfe_coef,dim=2),ubound(gbfe_coef,dim=2)
          !
          ! NO NONDIMENSIONALIZATION NEEDED FOR gbfe_coef(1,:,:)
          !
          tdim = one
          do m = 2,ncpc
            tdim = tdim*tref
            gbfe_coef(m,n,ics) = gbfe_coef(m,n,ics) * tdim
          end do
          !
          gbfe_coef(6,n,ics) = gbfe_coef(6,n,ics) / tref
          !
          gbfe_coef(7,n,ics) = gbfe_coef(7,n,ics) + &
                               gbfe_coef(1,n,ics)*(log(tref)-one)
          !
        end do
      end do
      !
    end if
    !
    ! Non-dimensionalize dgbfe_coef
    !
    if (ncpc == 5) then
      !
      do ics = 1,ncs
        do n = lbound(dgbfe_coef,dim=2),ubound(dgbfe_coef,dim=2)
          !
          ! NO NONDIMENSIONALIZATION NEEDED FOR dgbfe_coef(2,:,:)
          !
          tdim = one/(tref*tref)
          !
          dgbfe_coef(6,n,ics) = dgbfe_coef(6,n,ics) * tdim
          !
          do m = 1,ncpc
            tdim = tdim*tref
           !if (m == 2) cycle ! ensure there is no round off error
            dgbfe_coef(m,n,ics) = dgbfe_coef(m,n,ics) * tdim
          end do
          !
        ! dgbfe_coef(1,n,ics) = dgbfe_coef(1,n,ics) / tref
        !!dgbfe_coef(2,n,ics) = dgbfe_coef(2,n,ics)
        ! dgbfe_coef(3,n,ics) = dgbfe_coef(3,n,ics) * tref
        ! dgbfe_coef(4,n,ics) = dgbfe_coef(4,n,ics) * tref * tref
        ! dgbfe_coef(5,n,ics) = dgbfe_coef(5,n,ics) * tref * tref * tref
        ! dgbfe_coef(6,n,ics) = dgbfe_coef(6,n,ics) / tref / tref
          !
        end do
      end do
      !
    end if
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
    ! Reorder the mass fractions information in bc_input so that it
    ! matches the order of the chemical species defined above.
    !
    do n = 1,size(bc_input)
      !
      ! Save the mass fractions and accompanying species names that
      ! were input for this boundary condition into temporary arrays.
      !
      tmp_mass_frac = bc_input(n)%mass_frac   ! F2003 AUTO-REALLOCATION
      species_list = bc_input(n)%species_name ! F2003 AUTO-REALLOCATION
      !
      ! Loop over the chemical species defined above for this simulation.
      !
      do ics = 1,ncs
        !
        this_species = uppercase( trim(adjustl(species_name(ics))) )
        !
        ! Default the current species mass fraction to zero for this
        ! boundary condition in case no mass fraction was provided
        ! in the input file.
        !
        bc_input(n)%mass_frac(ics) = zero
        bc_input(n)%species_name(ics) = this_species
        !
        ! Loop through the species given for this boundary condition in
        ! the input file to find one that matches the current species.
        ! NOTE: If no matches are found, we've already set the mass
        !       fraction for this species to zero.
        !
        input_loop: &
        do jcs = 1,size(tmp_mass_frac)
          !
          ! Extract the name of the species that accompanies the
          ! input mass fraction.
          !
          input_species = uppercase( trim(adjustl(species_list(jcs))) )
          !
          if (input_species == this_species) then
            !
            ! We found a match, so assign the input mass fraction back into
            ! the mass_frac component of bc_input for this boundary condition.
            !
            bc_input(n)%mass_frac(ics) = tmp_mass_frac(jcs)
            exit input_loop
            !
          end if
          !
        end do input_loop
        !
      end do
      !
      ! Zero out the mass fractions and delete the species names from
      ! bc_input(n)%mass_frac(ncs+1:) and bc_input(n)%species_name(ncs+1:).
      !
      do ics = ncs+1,size(bc_input(n)%mass_frac)
        bc_input(n)%mass_frac(ics) = zero
      end do
      !
      do ics = ncs+1,size(bc_input(n)%species_name)
        bc_input(n)%species_name(ics) = ""
      end do
      !
    end do
    !
    ! Make sure the values of the mass fractions make sense
    !
    do n = 1,size(bc_input)
      !
      if ( bc_input(n)%has_default_flow_conditions() ) cycle
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
  call debug_timer(leaving_procedure,pname)
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
#ifndef DISABLE_PURE
pure &
#endif
subroutine update_thermo_properties_1p(cv_1p,gp_1p,caller_name)
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: cv_1p(:)
  real(wp), intent(inout) :: gp_1p(:)
  !
  !.. Optional Arguments ..
  character(len=*), optional, intent(in) :: caller_name
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
  gp(1,1:ng) = gp_1p(1:ng)
  !$ACC DATA CREATE(cv,gp)
  !
  call update_thermo_properties_sp(1,1,cv,gp,loc="node",caller_name=caller_name)
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
#ifndef DISABLE_PURE
pure &
#endif
subroutine update_thermo_properties_fp(n1,n2,cv,gp,caller_name)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp), intent(inout) :: cv(:,:,:)
  real(wp), intent(inout) :: gp(:,:,:)
  !
  !.. Optional Arguments ..
  character(len=*), optional, intent(in) :: caller_name
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_thermo_properties_fp"
  !
continue
  !
  call update_thermo_properties_sp(n1,n2,cv(:,:,1),gp(:,:,1), &
                                         loc="face", &
                                         caller_name=caller_name)
  if (all([size(cv,dim=3),size(gp,dim=3)] > 1)) then
    call update_thermo_properties_sp(n1,n2,cv(:,:,2),gp(:,:,2), &
                                           loc="face", &
                                           caller_name=caller_name)
  end if
  !
end subroutine update_thermo_properties_fp
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine update_thermo_properties_sp(n1,n2,cv,gp,loc,caller_name)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp), intent(inout) :: cv(:,:)
  real(wp), intent(inout) :: gp(:,:)
  !
  !.. Optional Arguments ..
  character(len=*), optional, intent(in) :: loc
  character(len=*), optional, intent(in) :: caller_name
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
  call pure_debug_timer(entering_procedure,pname)
  !
  ! 1. Update the mixture density by summing the species densities
  !    and convert the species densities to mass fractions
  !
  if (ncs > 0) then
    !
    call pure_debug_timer(start_timer,"update_mixture_density")
    call update_mixture_density(n1,n2,cv)
    call pure_debug_timer( stop_timer,"update_mixture_density")
    !
    ! Convert the species densities to mass fractions
    !
    call pure_debug_timer(start_timer,"Convert species densities to SMF")
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(cv)
    do m = 1,ncs
      do n = n1,n2
        cv(n,m) = cv(n,m)/cv(n,nec)
      end do
    end do
    call pure_debug_timer( stop_timer,"Convert species densities to SMF")
    !
    ! 2. Compute the internal energy. Temporarily store the internal energy
    !    in the pressure location of gp.
    !
   !call compute_internal_energy(n1,n2,cv,gp(:,npres))
    !
    ! 3. Compute the mixture gas constant. Temporarily store the mixture
    !    gas constant in the gamma location of gp.
    !
    call pure_debug_timer(start_timer,"compute_mixture_rgas")
    call compute_mixture_rgas(n1,n2,cv(:,1:ncs),gp(:,ngamm))
    call pure_debug_timer( stop_timer,"compute_mixture_rgas")
    !
    ! 4. Compute static temperature and gas variables using Newton iteration
    !    NOTE: Currently, the gamma and pressure locations of gp respectively
    !          contain the mixture gas constant and the internal energy.
    !
#ifndef DISABLE_PURE
    call newton_solve_for_temperature(n1,n2,cv,gp(:,ngamm), &
                                            gp(:,npres),gp(:,ntemp), &
                                            gp(:,nenth),gp(:,ncp), &
                                            num_newton_errors, &
                                            max_newton_error, &
                                            loc=loc, &
                                            caller_name=caller_name, &
                                            dim_t_upon_entry=true, &
                                            nondim_t_upon_return=true)
#else
    call impure_newton_solve_for_temperature(n1,n2,cv,gp(:,ngamm), &
                                                   gp(:,npres),gp(:,ntemp), &
                                                   gp(:,nenth),gp(:,ncp), &
                                                   num_newton_errors, &
                                                   max_newton_error, &
                                                   loc=loc, &
                                                   caller_name=caller_name, &
                                                   dim_t_upon_entry=true, &
                                                   nondim_t_upon_return=true)
#endif
    !
    ! 5. Compute the static pressure
    !    NOTE: The gamma location of gp currently
    !          contains the mixture gas constant.
    !
    call pure_debug_timer(start_timer,"Compute static pressure")
    !$ACC PARALLEL LOOP PRESENT(gp,cv)
    do n = n1,n2
      gp(n,npres) = cv(n,nec) * gp(n,ngamm) * gp(n,ntemp)
    end do
    call pure_debug_timer( stop_timer,"Compute static pressure")
    !
    ! 6. Compute the ratio of the specific heats of the mixture.
    !    NOTE: gp(:,ngamm) contains the mixture gas constant going into
    !          compute_mixture_gamma, but returns containing the ratio of
    !          the specific heats of the mixture.
    !
    call pure_debug_timer(start_timer,"compute_mixture_gamma")
    call compute_mixture_gamma(n1,n2,gp(:,ncp),gp(:,ngamm))
    call pure_debug_timer( stop_timer,"compute_mixture_gamma")
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
    call pure_debug_timer(start_timer,"Convert SMF to species densities")
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(cv)
    do m = 1,ncs
      do n = n1,n2
        cv(n,m) = cv(n,m)*cv(n,nec)
      end do
    end do
    call pure_debug_timer( stop_timer,"Convert SMF to species densities")
    !
  else
    !
    ! Compute the pressure and temperature of a calorically perfect gas
    !
    call pure_debug_timer(start_timer,"update_thermo_cpg")
    call update_thermo_cpg(n1,n2,cv,gp)
    call pure_debug_timer( stop_timer,"update_thermo_cpg")
    !
  end if
  !
  call pure_debug_timer(leaving_procedure,pname)
  !
end subroutine update_thermo_properties_sp
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
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
      gp(n,ngamm) = gam
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
      gp(n,ngamm) = gam
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
      gp(n,ngamm) = gam
      !
    end do
  end if
  !
end subroutine update_thermo_cpg
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine update_mixture_density(n1,n2,cv)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp), intent(inout) :: cv(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: rho
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_mixture_density"
  !
continue
  !
  ! Update the mixture density by summing the species densities
  !
  if (ncs > 0) then
    !
    ! Make sure the species densities are >= zero
    !
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(cv)
    do m = 1,ncs
      do n = n1,n2
        cv(n,m) = max(cv(n,m),zero)
      end do
    end do
    !
    ! Sum the species densities to get the mixture density
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
  end if
  !
end subroutine update_mixture_density
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
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
#ifndef DISABLE_PURE
pure &
#endif
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
#ifndef DISABLE_PURE
pure &
#endif
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
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_mixture_rgas_sp(n1,n2,mass_frac,rmix,dimensional)
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
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: dimensional
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: rsum,rmix_dim
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_mixture_rgas_sp"
  !
continue
  !
  !$ACC PARALLEL LOOP PRESENT(mass_frac,rmix,rspec)
  do n = n1,n2
    rsum = zero
    !$ACC LOOP SEQ REDUCTION(+:rsum)
    do m = 1,ncs
      rsum = rsum + rspec(m)*mass_frac(n,m)
    end do
    rmix(n) = rsum
  end do
  !
  if (present(dimensional)) then
    if (dimensional) then
      rmix_dim = aref*aref/tref
      !$ACC PARALLEL LOOP PRESENT(rmix)
      do n = n1,n2
        rmix(n) = rmix(n) * rmix_dim
      end do
    end if
  end if
  !
end subroutine compute_mixture_rgas_sp
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
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
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_mixture_rgas_1p"
  !
continue
  !
  rmix = dot_product( rspec(1:ncs) , mass_frac(1:ncs) )
  !
  if (present(dimensional)) then
    if (dimensional) then
      rmix = rmix * aref*aref/tref
    end if
  end if
  !
end subroutine compute_mixture_rgas_1p
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine newton_solve_for_temperature(n1,n2,cv,rmix,p,t,h,cp, &
                                        ndnc,dtmx,loc,caller_name, &
                                        dim_t_upon_entry,nondim_t_upon_return)
  !
  !.. Use Statements ..
  use ovar, only : tref,tmin,tmax
  use ovar, only : newton_rlx,conv_temp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: cv(:,:) ! mass fractions and conserved variables
  real(wp),    intent(in) :: rmix(:)
  real(wp),    intent(in) :: p(:)
  real(wp), intent(inout) :: t(:)
  real(wp), intent(inout) :: h(:)
  real(wp), intent(inout) :: cp(:)
  integer,  intent(inout) :: ndnc
  real(wp), intent(inout) :: dtmx
  !
  !.. Optional Arguments
  character(len=*), optional, intent(in) :: loc
  character(len=*), optional, intent(in) :: caller_name
  logical(lk),      optional, intent(in) :: dim_t_upon_entry
  logical(lk),      optional, intent(in) :: nondim_t_upon_return
  !
  !.. Local Scalars ..
  integer  :: k,n
  real(wp) :: func,dfunc,dtemp,dterr
  real(wp) :: tmin_nd,tmax_nd
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: eint(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "newton_solve_for_temperature"
  !
  integer, parameter :: newton_iter = 20
  !
continue
  !
  call pure_debug_timer(entering_procedure,pname)
  !
  tmin_nd = tmin / tref
  tmax_nd = tmax / tref
  !
  ! Reset the temperature to the reference temperature wherever
  ! the temperature is less-than tmin_nd or greater-than tmax_nd.
  !
  where ( t(:) < tmin_nd ) t(:) = one
  where ( t(:) > tmax_nd ) t(:) = one
  !
  allocate ( eint(1:size(t)) , source=zero )
  !
  call compute_internal_energy(n1,n2,cv,eint)
  !
  ! Compute initial static enthalpy and specific heat of the mixture
  !
  call compute_mixture_enthalpy(n1,n2,cv(:,1:ncs),t,h)
  !
  call compute_mixture_cp(n1,n2,cv(:,1:ncs),t,cp)
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
      if (dterr > conv_temp) then
        ndnc = ndnc + 1
      end if
      !
    end do
    !
    ! Update the static enthalpy and specific heat of the mixture
    !
    call compute_mixture_enthalpy(n1,n2,cv(:,1:ncs),t,h)
    !
    call compute_mixture_cp(n1,n2,cv(:,1:ncs),t,cp)
    !
    ! Add some special logic to prevent false positives for failing
    ! to converge when the "not converged" temperature is still 
    ! safely within the min/max temperature bounds.
    ! NOTE: This is known to happen when the temperature is very close
    !       to the temperature that seperates the curve fit intervals,
    !       e.g., 1000 K if the curve fit intervals are 300 K to 1000 K
    !       and 1000 K to 5000 K.
    !
    if (k == newton_iter) then
      ndnc = 0
      do n = n1,n2
        if (any( [ abs(t(n)-tmin_nd) , abs(t(n)-tmax_nd) ] <= eps9 )) then
          ndnc = ndnc + 1
        end if
      end do
      if (ndnc == 0) then
        dtmx = half*conv_temp
      end if
    end if
    !
    ! Check for convergence
    !
    if (dtmx <= conv_temp) then
      exit newton_iteration
    end if
    !
  end do newton_iteration
  !
  if (allocated(eint)) then
    deallocate ( eint )
  end if
  !
  call pure_debug_timer(leaving_procedure,pname)
  !
end subroutine newton_solve_for_temperature
!
!###############################################################################
!
subroutine impure_newton_solve_for_temperature(n1,n2,cv,rmix,p,t,h,cp, &
                                               ndnc,dtmx,loc,caller_name, &
                                               dim_t_upon_entry, &
                                               nondim_t_upon_return)
  !
  !.. Use Statements ..
  use ovar, only : tref,tmin,tmax
  use ovar, only : newton_rlx,conv_temp
  use ovar, only : traceback_unconverged_newton_temps
  use ovar, only : ignore_thermo_newton_warnings
#ifdef SPECIAL_FOR_INTEL
  use ifcore, only : tracebackqq
#endif
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: cv(:,:) ! mass fractions and conserved variables
  real(wp),    intent(in) :: rmix(:)
  real(wp),    intent(in) :: p(:)
  real(wp), intent(inout) :: t(:)
  real(wp), intent(inout) :: h(:)
  real(wp), intent(inout) :: cp(:)
  integer,  intent(inout) :: ndnc
  real(wp), intent(inout) :: dtmx
  !
  !.. Optional Arguments
  character(len=*), optional, intent(in) :: loc
  character(len=*), optional, intent(in) :: caller_name
  logical(lk),      optional, intent(in) :: dim_t_upon_entry
  logical(lk),      optional, intent(in) :: nondim_t_upon_return
  !
  !.. Local Scalars ..
  integer  :: k,n
  real(wp) :: func,dfunc,dtemp,dterr
  real(wp) :: tmin_nd,tmax_nd
  !
  character(len=100) :: message
  !
  !.. Local Saved Scalars ..
  logical(lk), save :: output_thermo_coef = true
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: eint(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "impure_newton_solve_for_temperature"
  !
  integer, parameter :: newton_iter = 20
  !
 !real(wp), parameter :: conv_temp = eps2  ! NCC code
 !real(wp), parameter :: conv_temp = eps6
 !real(wp), parameter :: conv_temp = eps8  ! VULCAN if ref. Mach > 0.01
 !real(wp), parameter :: conv_temp = eps12 ! VULCAN if ref. Mach <= 0.01
  !
#ifdef DEBUG_ON
!#define DEBUG_TEMPERATURE_NEWTON_SOLVE
#endif
#ifdef DEBUG_TEMPERATURE_NEWTON_SOLVE
  real(wp),    allocatable :: dtemp_save(:,:,:)
  logical(lk), allocatable :: did_not_converge(:)
  integer,     allocatable :: ndnc_save(:)
  real(wp),    allocatable :: dtmx_save(:)
#endif
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  tmin_nd = tmin / tref
  tmax_nd = tmax / tref
  !
  ! Reset the temperature to the reference temperature wherever
  ! the temperature is less-than tmin_nd or greater-than tmax_nd.
  !
  where ( t(:) < tmin_nd ) t(:) = one
  where ( t(:) > tmax_nd ) t(:) = one
  !
  allocate ( eint(1:size(t)) , source=zero )
  !
  call debug_timer(start_timer,"compute_internal_energy")
  call compute_internal_energy(n1,n2,cv,eint)
  call debug_timer(stop_timer,"compute_internal_energy")
  !
  ! Compute initial static enthalpy and specific heat of the mixture
  !
  call debug_timer(start_timer,"compute_mixture_enthalpy")
  call compute_mixture_enthalpy(n1,n2,cv(:,1:ncs),t,h)
  call debug_timer(stop_timer,"compute_mixture_enthalpy")
  !
  call debug_timer(start_timer,"compute_mixture_cp")
  call compute_mixture_cp(n1,n2,cv(:,1:ncs),t,cp)
  call debug_timer(stop_timer,"compute_mixture_cp")
  !
#ifdef DEBUG_TEMPERATURE_NEWTON_SOLVE
  allocate ( did_not_converge(n1:n2) , source=fals )
  allocate ( dtemp_save(1:8,n1:n2,1:newton_iter) )
  allocate ( ndnc_save(1:newton_iter) , source=0 )
  allocate ( dtmx_save(1:newton_iter) , source=zero )
#endif
  !
  ! Compute the temperature using Newton solve
  !
  newton_iteration : &
  do k = 1,newton_iter
    !
    ndnc = 0
    dtmx = zero
    !
    write (message,'(a,i0)') "thermo Newton iteration # ",k
    call debug_timer(start_timer,message)
    do n = n1,n2
      !
      func = h(n) - rmix(n)*t(n) - eint(n)
      dfunc = cp(n) - rmix(n)
      !
      dtemp = func/dfunc
      !
#ifdef DEBUG_TEMPERATURE_NEWTON_SOLVE
      dtemp_save(1,n,k) = h(n)
      dtemp_save(2,n,k) = cp(n)
      dtemp_save(3,n,k) = func
      dtemp_save(4,n,k) = dfunc
      dtemp_save(5,n,k) = dtemp
      dtemp_save(6,n,k) = t(n)
#endif
      !
      t(n) = min( tmax_nd , max(tmin_nd,t(n)-dtemp) )
      !
      dterr = abs(dtemp) / t(n)
      !
#ifdef DEBUG_TEMPERATURE_NEWTON_SOLVE
      dtemp_save(7,n,k) = t(n)
      dtemp_save(8,n,k) = dterr
      if (any( [ abs(t(n)-tmin_nd) , abs(t(n)-tmax_nd) ] <= eps9 )) then
        ! Set dtemp_save(7,n,k) to the negative value of the updated temperature
        ! if the updated temperature was bounded by tmin_dt or tmax_dt.
        dtemp_save(7,n,k) = -t(n)
      end if
#endif
      !
      if (dterr > dtmx) dtmx = dterr
      !
      if (dterr > conv_temp) then
        ndnc = ndnc + 1
#ifdef DEBUG_TEMPERATURE_NEWTON_SOLVE
        did_not_converge(n) = true
      else
        did_not_converge(n) = fals
#endif
      end if
      !
    end do
    call debug_timer(stop_timer,message)
    !
    ! Update the static enthalpy and specific heat of the mixture
    !
    write (message,'(a,i0)') &
                "compute_mixture_enthalpy Newton iteration # ",k
    call debug_timer(start_timer,message)
    call compute_mixture_enthalpy(n1,n2,cv(:,1:ncs),t,h)
    call debug_timer(stop_timer,message)
    !
    write (message,'(a,i0)') &
                "compute_mixture_cp Newton iteration # ",k
    call debug_timer(start_timer,message)
    call compute_mixture_cp(n1,n2,cv(:,1:ncs),t,cp)
    call debug_timer(stop_timer,message)
    !
    ! Add some special logic to prevent false positives for failing
    ! to converge when the "not converged" temperature is still 
    ! safely within the min/max temperature bounds.
    ! NOTE: This is known to happen when the temperature is very close
    !       to the temperature that seperates the curve fit intervals,
    !       e.g., 1000 K if the curve fit intervals are 300 K to 1000 K
    !       and 1000 K to 5000 K.
    !
    if (k == newton_iter) then
      ndnc = 0
      do n = n1,n2
        if (any( [ abs(t(n)-tmin_nd) , abs(t(n)-tmax_nd) ] <= eps9 )) then
          ndnc = ndnc + 1
#ifdef DEBUG_TEMPERATURE_NEWTON_SOLVE
          did_not_converge(n) = true
        else
          did_not_converge(n) = fals
#endif
        end if
      end do
      if (ndnc == 0) then
        dtmx = half*conv_temp
      end if
    end if
    !
#ifdef DEBUG_TEMPERATURE_NEWTON_SOLVE
    ndnc_save(k) = ndnc
    dtmx_save(k) = dtmx
#endif
    !
    ! Check for convergence
    !
    if (dtmx <= conv_temp) then
      exit newton_iteration
    end if
    !
    if (dtmx <= conv_temp .or. k == newton_iter) then
#ifdef DEBUG_TEMPERATURE_NEWTON_SOLVE
    block
      use geovar, only : nr,xyz
      use flowvar, only : facexyz
      use ovar, only : itcur,current_rk_stage,aref,rhoref
      character(len=50) :: fname
      integer :: i,j,m,n,ierr
      real(wp) :: vx,vy,vz,q2,r,theta
      character(len=5), parameter :: cxyz(1:5) = ["X    ", &
                                                  "Y    ", &
                                                  "Z    ", &
                                                  "R    ", &
                                                  "THETA"]
      if (.not. file_is_opened) then
        write (fname,6) mypnum
        open (newunit=io,file=trim(adjustl(fname)),status="replace", &
                         action="write",form="formatted",iostat=ierr)
        file_is_opened = true
      end if
      do n = n1,n2
        if (did_not_converge(n)) then
          if (last_itcur_output /= itcur) then
            write (io,10) itcur
            last_itcur_output = itcur
          end if
          if (last_rk_stage_output /= current_rk_stage) then
            write (io,11) current_rk_stage
            last_rk_stage_output = current_rk_stage
          end if
          write (io,1) "n",n
          write (io,1) "n1",n1
          write (io,1) "n2",n2
          if (present(loc)) then
            if (trim(adjustl(loc)) == "cell") then
              if (n >= lbound(xyz,dim=2) .and. n <= ubound(xyz,dim=2)) then
                do j = 1,nr
                  write (io,7) "S",cxyz(j),xyz(j,n)
                end do
                if (nr >= 2) then
                  r = norm2(xyz(2:,n))
                  write (io,7) "S",cxyz(4),r
                  if (nr == 3) then
                    theta = atan2(xyz(3,n),xyz(2,n))
                    if (theta < zero) theta = theta + two*pi
                    theta = theta*radians2degrees
                    write (io,7) "S",cxyz(5),theta
                  end if
                end if
              end if
            else if (trim(adjustl(loc)) == "face") then
              if (allocated(facexyz)) then
                if (n >= lbound(facexyz,dim=2) .and. &
                    n <= ubound(facexyz,dim=2)) then
                  do j = 1,nr
                    write (io,7) "F",cxyz(j),facexyz(j,n)
                  end do
                  if (nr >= 2) then
                    r = norm2(facexyz(2:,n))
                    write (io,7) "F",cxyz(4),r
                    if (nr == 3) then
                      theta = atan2(facexyz(3,n),facexyz(2,n))
                      if (theta < zero) theta = theta + two*pi
                      theta = theta*radians2degrees
                      write (io,7) "F",cxyz(5),theta
                    end if
                  end if
                end if
              end if
            end if
          end if
          !
          vy = zero ; vz = zero
          vx = cv(n,nmx)/cv(n,nec)
          if (nr >= 2) vy = cv(n,nmy)/cv(n,nec)
          if (nr == 3) vz = cv(n,nmz)/cv(n,nec)
          q2 = vx*vx + vy*vy + vz*vz
          !
          write (io,2) "rmix",rmix(n),rmix(n)*aref*aref/tref
          write (io,2) "eint",eint(n),eint(n)*rhoref*aref*aref
          write (io,2) "rho",cv(n,nec),cv(n,nec)*rhoref
          write (io,2) "vx",vx,vx*aref
          if (nr >= 2) &
          write (io,2) "vy",vy,vy*aref
          if (nr == 3) &
          write (io,2) "vz",vz,vz*aref
          write (io,2) "q^2",q2,q2*aref*aref
          write (io,2) "rhoE",cv(n,nee),cv(n,nee)*rhoref*aref*aref
          write (io,2) "pres",p(n),p(n)*rhoref*aref*aref
          write (io,2) "t_beg",dtemp_save(6,n,1), &
                               dtemp_save(6,n,1)*tref
          write (io,2) "h_beg",dtemp_save(1,n,1), &
                               dtemp_save(1,n,1)*aref*aref
          write (io,2) "cp_beg",dtemp_save(2,n,1), &
                                dtemp_save(2,n,1)*aref*aref/tref
          write (io,2) "t_end",t(n),t(n)*tref
          write (io,2) "h_end",h(n),h(n)*aref*aref
          write (io,2) "cp_end",cp(n),cp(n)*aref*aref/tref
          do m = 1,ncs
            write (io,3) trim(adjustl(species_name(m))),cv(n,m)
          end do
          if (output_thermo_coef) then
            do m = 1,ncs
              !
             !if (cv(n,m) <= eps6) cycle
              !
              write (io,20) trim(adjustl(species_name(m))),"enthalpy"
              write (io,21) 0, (enth_coef(6,i,m), i=lbound(enth_coef,dim=2), &
                                                    ubound(enth_coef,dim=2))
              do j = 1,ncpc
                write (io,21) j, (enth_coef(j,i,m), i=lbound(enth_coef,dim=2), &
                                                      ubound(enth_coef,dim=2))
              end do
              !
              write (io,20) trim(adjustl(species_name(m))),"cp"
              do j = 1,ncpc
                write (io,21) j-1, (cp_coef(j,i,m), i=lbound(cp_coef,dim=2), &
                                                      ubound(cp_coef,dim=2))
              end do
              !
              write (io,20) trim(adjustl(species_name(m))),"d(cp)/dT"
              do j = 2,ncpc
                write (io,21) j-2, &
                  (real(j-1,kind=wp)*cp_coef(j,i,m), i=lbound(cp_coef,dim=2), &
                                                       ubound(cp_coef,dim=2))
              end do
              !
            end do
            output_thermo_coef = fals
          end if
          do j = 1,size(dtemp_save,dim=3)
            if (j == 1) write (io,4)
            if (j == size(dtemp_save,dim=3)) then
              write (io,5) j,(dtemp_save(i,n,j), i=1,size(dtemp_save,dim=1)),"-"
            else
              write (io,5) j,(dtemp_save(i,n,j), i=1,size(dtemp_save,dim=1))
            end if
          end do
          flush (io)
  1 format (4x,a,t20," = ",i0)
  2 format (4x,a,t20," = ",es23.15," ; ",es23.15)
  3 format (4x,a," mass frac",t20," = ",es23.15)
  4 format (/,101("-"),/, &
              4x,"VALUES AFFECTED BY EACH NEWTON ITERATION",/, &
              101("-"),/, &
              1x,"Iter",8x,"enth",13x,"cp",12x,"func",12x,"dfunc", &
              11x,"dtemp",9x,"temp_beg",8x,"temp_end",8x,"temp_err",/, &
              101("-"))
  5 format (2x,i3,8(2x,es14.6),:,/,100("-"),a,/)
  6 format ("newton_temp_errors.CPU-",i0,".log")
  7 format (4x,"(",a1,"P) ",a5,t20," = ",es23.15)
  !
  10 format (/,2(101("#"),/), 20x,"TIME STEP # ",i0, /,2(101("#"),/))
  11 format (/,  101("="),/ , 10x, "RK STAGE # ",i0, /,  101("="),/ )
  20 format (/,4x,a," curve fit coefficients")
  21 format (6x,"T^",i0," coefficient : ",5(es23.15,:,2x))

        end if
      end do
    end block
#endif
      if (.not. ignore_thermo_newton_warnings) then
        if (present(caller_name)) then
          write (iout,101) ndnc,dtmx,trim(adjustl(caller_name))
#ifdef SPECIAL_FOR_INTEL
          if (traceback_unconverged_newton_temps) then
            call tracebackqq(user_exit_code=-1)
          end if
#endif
        else
          write (iout,100) ndnc,dtmx
#ifdef SPECIAL_FOR_INTEL
          if (traceback_unconverged_newton_temps) then
            call tracebackqq(user_exit_code=-1)
          end if
#endif
        end if
      end if
      exit newton_iteration
    end if
    !
  end do newton_iteration
  !
  if (allocated(eint)) then
    deallocate ( eint )
  end if
  !
  call debug_timer(leaving_procedure,pname)
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
  101 format (/, &
   1x,64("*"),/, &
   1x,5("*"),7x,"NEWTON ITERATION TO SOLVE FOR THE STATIC",t61,5("*"),/, &
   1x,5("*"),4x,"TEMPERATURE OF A GAS MIXTURE DID NOT CONVERGE!",t61,5("*"),/, &
   1x,5("*"),7x,"NUMBER OF UNCONVERGED POINTS = ",i0,t61,5("*"),/, &
   1x,5("*"),7x,"MAXIMUM ERROR IN TEMPERATURE = ",es13.6,t61,5("*"),/, &
   1x,5("*"),7x,"CALLED FROM PROCEDURE = '",a,"'",/, &
   1x,64("*"),/)
   !
end subroutine impure_newton_solve_for_temperature
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
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
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_mixture_enthalpy_sp(n1,n2,mass_frac,t,h,dimensional)
  !
  !.. Use Statements ..
  use ovar, only : aref
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: mass_frac(:,:)
  real(wp),    intent(in) :: t(:)
  real(wp), intent(inout) :: h(:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: dimensional
  !
  !.. Local Scalars ..
  integer  :: m,n,ni
  real(wp) :: hi,ts,h_dim
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: intrvl(:)
  !
continue
  !
  allocate ( intrvl(n1:n2) , source=1 )
  !$ACC DATA CREATE(intrvl)
  !
  !$ACC PARALLEL LOOP PRESENT(h)
  do n = n1,n2
    h(n) = zero
  end do
  !
  !$ACC PARALLEL LOOP PRESENT(t,intrvl)
  do n = n1,n2
    if (t(n) > temperature_high) then
      intrvl(n) = 2
    else if (t(n) < temperature_low) then
      intrvl(n) = 0
    end if
  end do
  !
  do m = 1,ncs
    !$ACC PARALLEL LOOP PRESENT(mass_frac,intrvl,t,enth_coef,h)
    do n = n1,n2
      !
      ni = intrvl(n)
      ts = t(n)
      !
      hi = enth_coef(6,ni,m) + &
           ts*(enth_coef(1,ni,m) + &
               ts*(enth_coef(2,ni,m) + &
                   ts*(enth_coef(3,ni,m) + &
                       ts*(enth_coef(4,ni,m) + &
                           ts*(enth_coef(5,ni,m))))))
      !
      h(n) = h(n) + mass_frac(n,m)*hi
      !
    end do
  end do
  !
  !$ACC END DATA
  deallocate ( intrvl )
  !
  if (present(dimensional)) then
    if (dimensional) then
      h_dim = aref*aref
      !$ACC PARALLEL LOOP PRESENT(h)
      do n = n1,n2
        h(n) = h(n)*h_dim
      end do
    end if
  end if
  !
end subroutine compute_mixture_enthalpy_sp
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_mixture_enthalpy_1p(mass_frac,t,h,dimensional)
  !
  !.. Use Statements ..
  use ovar, only : aref
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
  integer  :: m,ni
  real(wp) :: hi
  !
continue
  !
  h = zero
  !
  ni = 1
  if (t > temperature_high) then
    ni = 2
  else if (t < temperature_low) then
    ni = 0
  end if
  !
  do m = 1,ncs
    !
    hi = enth_coef(6,ni,m) + &
         t*(enth_coef(1,ni,m) + &
             t*(enth_coef(2,ni,m) + &
                 t*(enth_coef(3,ni,m) + &
                     t*(enth_coef(4,ni,m) + &
                         t*(enth_coef(5,ni,m))))))
    !
    h = h + mass_frac(m)*hi
    !
  end do
  !
  if (present(dimensional)) then
    if (dimensional) then
      h = h*aref*aref
    end if
  end if
  !
end subroutine compute_mixture_enthalpy_1p
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_species_enthalpy_sp(n1,n2,t,hsi,dimensional)
  !
  !.. Use Statements ..
  use ovar, only : aref
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: t(:)
  real(wp), intent(inout) :: hsi(:,:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: dimensional
  !
  !.. Local Scalars ..
  integer  :: m,n,ni
  real(wp) :: ts,h_dim
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: intrvl(:)
  !
continue
  !
  allocate ( intrvl(n1:n2) , source=1 )
  !$ACC DATA CREATE(intrvl)
  !
  !$ACC PARALLEL LOOP PRESENT(t,intrvl)
  do n = n1,n2
    if (t(n) > temperature_high) then
      intrvl(n) = 2
    else if (t(n) < temperature_low) then
      intrvl(n) = 0
    end if
  end do
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(intrvl,t,enth_coef,hsi)
  do m = 1,ncs
    do n = n1,n2
      !
      ni = intrvl(n)
      ts = t(n)
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
  !$ACC END DATA
  deallocate ( intrvl )
  !
  if (present(dimensional)) then
    if (dimensional) then
      h_dim = aref*aref
      !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(hsi)
      do m = 1,ncs
        do n = n1,n2
          hsi(n,m) = hsi(n,m)*h_dim
        end do
      end do
    end if
  end if
  !
end subroutine compute_species_enthalpy_sp
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_species_enthalpy_1p(t,hsi,dimensional)
  !
  !.. Use Statements ..
  use ovar, only : aref
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: t
  real(wp), intent(inout) :: hsi(:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: dimensional
  !
  !.. Local Scalars ..
  integer  :: m,ni
  real(wp) :: h_dim
  !
continue
  !
  ni = 1
  if (t > temperature_high) then
    ni = 2
  else if (t < temperature_low) then
    ni = 0
  end if
  !
  do m = 1,ncs
    !
    hsi(m) = enth_coef(6,ni,m) + &
             t*(enth_coef(1,ni,m) + &
                 t*(enth_coef(2,ni,m) + &
                     t*(enth_coef(3,ni,m) + &
                         t*(enth_coef(4,ni,m) + &
                             t*(enth_coef(5,ni,m))))))
    !
  end do
  !
  if (present(dimensional)) then
    if (dimensional) then
      h_dim = aref*aref
      do m = 1,ncs
        hsi(m) = hsi(m)*h_dim
      end do
    end if
  end if
  !
end subroutine compute_species_enthalpy_1p
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_mixture_href(mass_frac,t,h)
  !
  ! NOTE: This version is different from the others in
  !       that the input temperature (dummy argument 't'), the
  !       module variable temperature_high, and the module array
  !       enth_coef are assumed have dimensional values.
  !
  !       The other subroutines that compute either species or
  !       mixture enthalpy assume that the input temperature,
  !       temperature_high, and enth_coef are all non-dimensional.
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: mass_frac(:)
  real(wp),    intent(in) :: t
  real(wp), intent(inout) :: h
  !
  !.. Local Scalars ..
  integer  :: m,ni
  real(wp) :: hi
  !
continue
  !
  h = zero
  !
  ni = 1
  if (t > temperature_high) then
    ni = 2
  else if (t < temperature_low) then
    ni = 0
  end if
  !
  do m = 1,ncs
    !
    hi = enth_coef(6,ni,m) + &
         t*(enth_coef(1,ni,m) + &
             t*(enth_coef(2,ni,m) + &
                 t*(enth_coef(3,ni,m) + &
                     t*(enth_coef(4,ni,m) + &
                         t*(enth_coef(5,ni,m))))))
    !
    h = h + mass_frac(m)*hi
    !
  end do
  !
end subroutine compute_mixture_href
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_species_energy_sp(n1,n2,t,esi,dimensional)
  !
  !.. Use Statements ..
  use ovar, only : aref
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: t(:)
  real(wp), intent(inout) :: esi(:,:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: dimensional
  !
  !.. Local Scalars ..
  integer  :: m,n,ni
  real(wp) :: ts,e_dim
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: intrvl(:)
  !
continue
  !
  allocate ( intrvl(n1:n2) , source=1 )
  !$ACC DATA CREATE(intrvl)
  !
  !$ACC PARALLEL LOOP PRESENT(t,intrvl)
  do n = n1,n2
    if (t(n) > temperature_high) then
      intrvl(n) = 2
    else if (t(n) < temperature_low) then
      intrvl(n) = 0
    end if
  end do
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(intrvl,t,enth_coef,esi,rspec)
  do m = 1,ncs
    do n = n1,n2
      !
      ni = intrvl(n)
      ts = t(n)
      !
      esi(n,m) = enth_coef(6,ni,m) + &
                 ts*(enth_coef(1,ni,m) + &
                     ts*(enth_coef(2,ni,m) + &
                         ts*(enth_coef(3,ni,m) + &
                             ts*(enth_coef(4,ni,m) + &
                                 ts*(enth_coef(5,ni,m))))))
      !
      ! esi currently contains the species enthalpy, now subtract
      ! from that [species gas constant]*T to get the species energy.
      !
      esi(n,m) = esi(n,m) - rspec(m)*ts
      !
    end do
  end do
  !
  !$ACC END DATA
  deallocate ( intrvl )
  !
  if (present(dimensional)) then
    if (dimensional) then
      e_dim = aref*aref
      !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(esi)
      do m = 1,ncs
        do n = n1,n2
          esi(n,m) = esi(n,m)*e_dim
        end do
      end do
    end if
  end if
  !
end subroutine compute_species_energy_sp
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_mixture_cp_sp(n1,n2,mass_frac,t,cp,dimensional)
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
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: dimensional
  !
  !.. Local Scalars ..
  integer  :: m,n,ni
  real(wp) :: cpi,ts,cp_dim
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: intrvl(:)
  !
continue
  !
  allocate ( intrvl(n1:n2) , source=1 )
  !$ACC DATA CREATE(intrvl)
  !
  !$ACC PARALLEL LOOP PRESENT(cp)
  do n = n1,n2
    cp(n) = zero
  end do
  !
  !$ACC PARALLEL LOOP PRESENT(t,intrvl)
  do n = n1,n2
    if (t(n) > temperature_high) then
      intrvl(n) = 2
    else if (t(n) < temperature_low) then
      intrvl(n) = 0
    end if
  end do
  !
  do m = 1,ncs
    !$ACC PARALLEL LOOP PRESENT(mass_frac,intrvl,t,cp_coef,cp)
    do n = n1,n2
      !
      ni = intrvl(n)
      ts = t(n)
      !
      cpi = cp_coef(1,ni,m) + &
            ts*(cp_coef(2,ni,m) + &
                ts*(cp_coef(3,ni,m) + &
                    ts*(cp_coef(4,ni,m) + &
                        ts*(cp_coef(5,ni,m)))))
      !
      cp(n) = cp(n) + mass_frac(n,m)*cpi
      !
    end do
  end do
  !
  !$ACC END DATA
  deallocate ( intrvl )
  !
  if (present(dimensional)) then
    if (dimensional) then
      cp_dim = aref*aref/tref
      !$ACC PARALLEL LOOP PRESENT(cp)
      do n = n1,n2
        cp(n) = cp(n)*cp_dim
      end do
    end if
  end if
  !
end subroutine compute_mixture_cp_sp
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
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
  integer  :: m,ni
  real(wp) :: cpi
  !
continue
  !
  cp = zero
  !
  ni = 1
  if (t > temperature_high) then
    ni = 2
  else if (t < temperature_low) then
    ni = 0
  end if
  !
  do m = 1,ncs
    !
    cpi = cp_coef(1,ni,m) + &
          t*(cp_coef(2,ni,m) + &
              t*(cp_coef(3,ni,m) + &
                  t*(cp_coef(4,ni,m) + &
                      t*(cp_coef(5,ni,m)))))
    !
    cp = cp + mass_frac(m)*cpi
    !
  end do
  !
  if (present(dimensional)) then
    if (dimensional) then
      cp = cp * aref*aref/tref
    end if
  end if
  !
end subroutine compute_mixture_cp_1p
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_species_cp_sp(n1,n2,t,cpi,dimensional)
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
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: dimensional
  !
  !.. Local Scalars ..
  integer  :: m,n,ni
  real(wp) :: ts,cp_dim
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: intrvl(:)
  !
continue
  !
  allocate ( intrvl(n1:n2) , source=1 )
  !$ACC DATA CREATE(intrvl)
  !
  !$ACC PARALLEL LOOP PRESENT(t,intrvl)
  do n = n1,n2
    if (t(n) > temperature_high) then
      intrvl(n) = 2
    else if (t(n) < temperature_low) then
      intrvl(n) = 0
    end if
  end do
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(intrvl,t,cp_coef,cpi)
  do m = 1,ncs
    do n = n1,n2
      !
      ni = intrvl(n)
      ts = t(n)
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
  if (present(dimensional)) then
    if (dimensional) then
      cp_dim = aref*aref/tref
      !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(cpi)
      do m = 1,ncs
        do n = n1,n2
          cpi(n,m) = cpi(n,m)*cp_dim
        end do
      end do
    end if
  end if
  !
  !$ACC END DATA
  deallocate ( intrvl )
  !
end subroutine compute_species_cp_sp
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_species_cp_1p(t,cpi,dimensional)
  !
  !.. Use Statements ..
  use ovar, only : aref,tref
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: t
  real(wp), intent(inout) :: cpi(:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: dimensional
  !
  !.. Local Scalars ..
  integer  :: m,ni
  real(wp) :: cp_dim
  !
continue
  !
  ni = 1
  if (t > temperature_high) then
    ni = 2
  else if (t < temperature_low) then
    ni = 0
  end if
  !
  do m = 1,ncs
    !
    cpi(m) = cp_coef(1,ni,m) + &
             t*(cp_coef(2,ni,m) + &
                 t*(cp_coef(3,ni,m) + &
                     t*(cp_coef(4,ni,m) + &
                         t*(cp_coef(5,ni,m)))))
    !
  end do
  !
  !
  if (present(dimensional)) then
    if (dimensional) then
      cp_dim = aref*aref/tref
      do m = 1,ncs
        cpi(m) = cpi(m)*cp_dim
      end do
    end if
  end if
  !
end subroutine compute_species_cp_1p
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_mixture_cpref(mass_frac,t,cp)
  !
  ! NOTE: This version is different from the others in
  !       that the input temperature (dummy argument 't'), the
  !       module variable temperature_high, and the module array
  !       cp_coef are assumed have dimensional values.
  !
  !       The other subroutines that compute either species or
  !       mixture Cp assume that the input temperature,
  !       temperature_high, and cp_coef are all non-dimensional.
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: mass_frac(:)
  real(wp),    intent(in) :: t
  real(wp), intent(inout) :: cp
  !
  !.. Local Scalars ..
  integer  :: m,ni
  real(wp) :: cpi
  !
continue
  !
  cp = zero
  !
  ni = 1
  if (t > temperature_high) then
    ni = 2
  else if (t < temperature_low) then
    ni = 0
  end if
  !
  do m = 1,ncs
    !
    cpi = cp_coef(1,ni,m) + &
          t*(cp_coef(2,ni,m) + &
              t*(cp_coef(3,ni,m) + &
                  t*(cp_coef(4,ni,m) + &
                      t*(cp_coef(5,ni,m)))))
    !
    cp = cp + mass_frac(m)*cpi
    !
  end do
  !
end subroutine compute_mixture_cpref
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
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
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: mass_frac(:,:)
  real(wp),    intent(in) :: t(:)
  real(wp), intent(inout) :: dcpdt(:)
  !
  !.. Local Scalars ..
  integer  :: m,n,ni
  real(wp) :: dcpidt,ts
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: intrvl(:)
  !
continue
  !
  allocate ( intrvl(n1:n2) , source=1 )
  !$ACC DATA CREATE(intrvl)
  !
  !$ACC PARALLEL LOOP PRESENT(dcpdt)
  do n = n1,n2
    dcpdt(n) = zero
  end do
  !
  !$ACC PARALLEL LOOP PRESENT(t,intrvl)
  do n = n1,n2
    if (t(n) > temperature_high) then
      intrvl(n) = 2
    else if (t(n) < temperature_low) then
      intrvl(n) = 0
    end if
  end do
  !
  do m = 1,ncs
    !$ACC PARALLEL LOOP PRESENT(mass_frac,intrvl,t,cp_coef,dcpdt)
    do n = n1,n2
      !
      ni = intrvl(n)
      ts = t(n)
      !
      dcpidt = cp_coef(2,ni,m) + &
               ts*(cp_coef(3,ni,m)*two + &
                   ts*(cp_coef(4,ni,m)*three + &
                       ts*(cp_coef(5,ni,m)*four)))
      !
      dcpdt(n) = dcpdt(n) + mass_frac(n,m)*dcpidt
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
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_species_dcpdt(n1,n2,t,dcpidt)
  !
  ! Computes the partial derivative of the species Cp w.r.t. temperature:
  !
  !             d [Cp(m)]
  !   dcpidt = -----------
  !                dT
  !
  ! for points n1 to n2.
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: t(:)
  real(wp), intent(inout) :: dcpidt(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,ni
  real(wp) :: ts
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: intrvl(:)
  !
continue
  !
  allocate ( intrvl(n1:n2) , source=1 )
  !$ACC DATA CREATE(intrvl)
  !
  !$ACC PARALLEL LOOP PRESENT(t,intrvl)
  do n = n1,n2
    if (t(n) > temperature_high) then
      intrvl(n) = 2
    else if (t(n) < temperature_low) then
      intrvl(n) = 0
    end if
  end do
  !
  do m = 1,ncs
    !$ACC PARALLEL LOOP PRESENT(intrvl,t,cp_coef,dcpidt)
    do n = n1,n2
      !
      ni = intrvl(n)
      ts = t(n)
      !
      dcpidt(n,m) = cp_coef(2,ni,m) + &
                    ts*(cp_coef(3,ni,m)*two + &
                        ts*(cp_coef(4,ni,m)*three + &
                            ts*(cp_coef(5,ni,m)*four)))
      !
    end do
  end do
  !
  !$ACC END DATA
  deallocate ( intrvl )
  !
end subroutine compute_species_dcpdt
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_cp_over_temp_integral_sp(n1,n2,mass_frac, &
                                            return_value,temp,ttot)
  !
  ! Routine to compute the integral required for calculating total pressure
  ! for a mixture of gases in thermal equilibrium using the Gordon and
  ! McBride curve fits
  !
  ! N1           : beginning point
  ! N2           : ending    point
  ! MASS_FRAC    : species mass fraction
  ! TTOT         : total temperature
  ! TEMP         : static temperature
  ! RETURN_VALUE : Final intergral of Cp/T summed over all species
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
  real(wp) :: tb0,te0,cc0,ce0,cb0,a0,b0,c0,d0,e0,cf0
  real(wp) :: tb1,te1,cc1,ce1,cb1,a1,b1,c1,d1,e1,cf1
  real(wp) :: tb2,te2,cc2,ce2,cb2,a2,b2,c2,d2,e2,cf2
  real(wp) :: ts,tt
  !
  !.. Local Parameters ..
  real(wp), parameter :: three2 = three/two
  real(wp), parameter :: four3  = four/three
  real(wp), parameter :: five4  = five/four
  !
continue
  !
  !$ACC PARALLEL LOOP PRESENT(return_value)
  do n = n1,n2
    return_value(n) = zero
  end do
  !
  ! 6 POSSIBLE SCENARIOS
  !
  !   ------------------------------------------------
  !   | SCENARIO |   NI=0    |   NI=1    |   NI=2    |
  !   |----------|-----------|-----------|-----------|
  !   |     1    | temp,ttot |           |           |
  !   |     2    |   temp    |   ttot    |           |
  !   |     3    |   temp    |           |   ttot    |
  !   |     4    |           | temp,ttot |           |
  !   |     5    |           |   temp    |   ttot    |
  !   |     6    |           |           | temp,ttot |
  !   ------------------------------------------------
  !
  do m = 1,ncs
    !
    a0 = enth_coef(1,0,m)
    b0 = enth_coef(2,0,m)*two
    c0 = enth_coef(3,0,m)*three2
    d0 = enth_coef(4,0,m)*four3
    e0 = enth_coef(5,0,m)*five4
    !
    a1 = enth_coef(1,1,m)
    b1 = enth_coef(2,1,m)*two
    c1 = enth_coef(3,1,m)*three2
    d1 = enth_coef(4,1,m)*four3
    e1 = enth_coef(5,1,m)*five4
    !
    a2 = enth_coef(1,2,m)
    b2 = enth_coef(2,2,m)*two
    c2 = enth_coef(3,2,m)*three2
    d2 = enth_coef(4,2,m)*four3
    e2 = enth_coef(5,2,m)*five4
    !
    !$ACC PARALLEL LOOP PRESENT(ttot,temp,mass_frac,return_value)
    do n = n1,n2
      !
      ts = temp(n)
      tt = ttot(n)
      !
      cf0 = zero
      cf1 = zero
      cf2 = zero
      !
      ! Compute the contribution, if any, from the constant Cp curvefit interval
      !
      if (ts < temperature_low) then
        !
        tb0 = ts
        te0 = min(tt,temperature_low)
        !
        cc0 = log(te0/tb0)*a0
        !
        ce0 = te0*(b0 + te0*(c0 + te0*(d0 + te0*e0)))
        cb0 = tb0*(b0 + tb0*(c0 + tb0*(d0 + tb0*e0)))
        !
        cf0 = mass_frac(n,m)*(cc0 + ce0 - cb0)
        !
      end if
      !
      ! Compute the contribution, if any, from the
      ! low temperature Cp curvefit interval
      !
      if (ts <= temperature_high .and. tt >= temperature_low) then
        !
        tb1 = max(ts,temperature_low)
        te1 = min(tt,temperature_high)
        !
        cc1 = log(te1/tb1)*a1
        !
        ce1 = te1*(b1 + te1*(c1 + te1*(d1 + te1*e1)))
        cb1 = tb1*(b1 + tb1*(c1 + tb1*(d1 + tb1*e1)))
        !
        cf1 = mass_frac(n,m)*(cc1 + ce1 - cb1)
        !
      end if
      !
      ! Compute the contribution, if any, from the
      ! high temperature Cp curvefit interval
      !
      if (tt > temperature_high) then
        !
        tb2 = max(ts,temperature_high)
        te2 = tt
        !
        cc2 = log(te2/tb2)*a2
        !
        ce2 = te2*(b2 + te2*(c2 + te2*(d2 + te2*e2)))
        cb2 = tb2*(b2 + tb2*(c2 + tb2*(d2 + tb2*e2)))
        !
        cf2 = mass_frac(n,m)*(cc2 + ce2 - cb2)
        !
      end if
      !
      return_value(n) = return_value(n) + cf0 + cf1 + cf2
      !
    end do
    !
  end do
  !
end subroutine compute_cp_over_temp_integral_sp
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_cp_over_temp_integral_1p(mass_frac,return_value,temp,ttot)
  !
  ! Routine to compute the integral required for calculating total pressure
  ! for a mixture of gases in thermal equilibrium using the Gordon and
  ! McBride curve fits
  !
  ! MASS_FRAC    : species mass fraction
  ! TTOT         : total temperature
  ! TEMP         : static temperature
  ! RETURN_VALUE : Final intergral of Cp/T summed over all species
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: mass_frac(:)
  real(wp),    intent(in) :: ttot
  real(wp),    intent(in) :: temp
  real(wp), intent(inout) :: return_value
  !
  !.. Local Scalars ..
  integer  :: m
  real(wp) :: tb0,te0,cc0,ce0,cb0,a0,b0,c0,d0,e0,cf0
  real(wp) :: tb1,te1,cc1,ce1,cb1,a1,b1,c1,d1,e1,cf1
  real(wp) :: tb2,te2,cc2,ce2,cb2,a2,b2,c2,d2,e2,cf2
  real(wp) :: ts,tt
  !
  !.. Local Parameters ..
  real(wp), parameter :: three2 = three/two
  real(wp), parameter :: four3  = four/three
  real(wp), parameter :: five4  = five/four
  !
continue
  !
  return_value = zero
  !
  ! 6 POSSIBLE SCENARIOS
  !
  !   ------------------------------------------------
  !   | SCENARIO |   NI=0    |   NI=1    |   NI=2    |
  !   |----------|-----------|-----------|-----------|
  !   |     1    | temp,ttot |           |           |
  !   |     2    |   temp    |   ttot    |           |
  !   |     3    |   temp    |           |   ttot    |
  !   |     4    |           | temp,ttot |           |
  !   |     5    |           |   temp    |   ttot    |
  !   |     6    |           |           | temp,ttot |
  !   ------------------------------------------------
  !
  do m = 1,ncs
    !
    a0 = enth_coef(1,0,m)
    b0 = enth_coef(2,0,m)*two
    c0 = enth_coef(3,0,m)*three2
    d0 = enth_coef(4,0,m)*four3
    e0 = enth_coef(5,0,m)*five4
    !
    a1 = enth_coef(1,1,m)
    b1 = enth_coef(2,1,m)*two
    c1 = enth_coef(3,1,m)*three2
    d1 = enth_coef(4,1,m)*four3
    e1 = enth_coef(5,1,m)*five4
    !
    a2 = enth_coef(1,2,m)
    b2 = enth_coef(2,2,m)*two
    c2 = enth_coef(3,2,m)*three2
    d2 = enth_coef(4,2,m)*four3
    e2 = enth_coef(5,2,m)*five4
    !
    cf0 = zero
    cf1 = zero
    cf2 = zero
    !
    ts = temp
    tt = ttot
    !
    ! Compute the contribution, if any, from the constant Cp curvefit interval
    !
    if (ts < temperature_low) then
      !
      tb0 = ts
      te0 = min(tt,temperature_low)
      !
      cc0 = log(te0/tb0)*a0
      !
      ce0 = te0*(b0 + te0*(c0 + te0*(d0 + te0*e0)))
      cb0 = tb0*(b0 + tb0*(c0 + tb0*(d0 + tb0*e0)))
      !
      cf0 = mass_frac(m)*(cc0 + ce0 - cb0)
      !
    end if
    !
    ! Compute the contribution, if any, from the
    ! low temperature Cp curvefit interval
    !
    if (ts <= temperature_high .and. tt >= temperature_low) then
      !
      tb1 = max(ts,temperature_low)
      te1 = min(tt,temperature_high)
      !
      cc1 = log(te1/tb1)*a1
      !
      ce1 = te1*(b1 + te1*(c1 + te1*(d1 + te1*e1)))
      cb1 = tb1*(b1 + tb1*(c1 + tb1*(d1 + tb1*e1)))
      !
      cf1 = mass_frac(m)*(cc1 + ce1 - cb1)
      !
    end if
    !
    ! Compute the contribution, if any, from the
    ! high temperature Cp curvefit interval
    !
    if (tt > temperature_high) then
      !
      tb2 = max(ts,temperature_high)
      te2 = tt
      !
      cc2 = log(te2/tb2)*a2
      !
      ce2 = te2*(b2 + te2*(c2 + te2*(d2 + te2*e2)))
      cb2 = tb2*(b2 + tb2*(c2 + tb2*(d2 + tb2*e2)))
      !
      cf2 = mass_frac(m)*(cc2 + ce2 - cb2)
      !
    end if
    !
    return_value = return_value + cf0 + cf1 + cf2
    !
  end do
  !
end subroutine compute_cp_over_temp_integral_1p
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_gibbs_free_energy(n1,n2,t,gbfe)
  !
  ! Routine to calculate the Gibbs free energy of a specie
  !
  ! N1   : starting point
  ! N2   : ending   point
  ! T    : nondimensional static temperature
  ! GBFE : Gibbs free energy / (R*T)
  !
  !
  ! NOTE: The form is that of an assumed equilibrium model:
  !
  !       kb(ir,T) = kf(ir,T) / ke(ir,T)
  !
  !       and
  !
  !       ke(ir,T) = (Ro*T)**(-DSCSUM)*exp(-deltaG(ir)/(Ro*T))
  !
  !       where:
  !
  !           -deltaG(ir)/(Ro*T) = sum of the Gibbs free energy times
  !                                the difference of the stoichiometic
  !                                coeffifients
  !
  !                                 ns
  !           -deltaG(ir)/(Ro*T) =  sum[-g(is)/(Ro*T)*(delsco(is,ir))]
  !                                 is=1
  !
  !    IF USING A DIMENSIONAL TEMPERATURE (T)
  !           -g(is)/(Ro*T) = -F(is)*T^-1 +
  !                            H(is)      +
  !                            B(is)*T    +
  !                            C(is)*T^2  +
  !                            D(is)*T^3  +
  !                            E(is)*T^4  +
  !                            A(is)*(ln(T)-1)
  !
  !    IF USING A NONDIMENSIONAL TEMPERATURE (T_nd)
  !           -g(is)/(Ro*T) = -F_nd(is)*T_nd^-1 +
  !                            H_nd(is)         +
  !                            B_nd(is)*T_nd    +
  !                            C_nd(is)*T_nd^2  +
  !                            D_nd(is)*T_nd^3  +
  !                            E_nd(is)*T_nd^4  +
  !                            A_nd(is)*ln(T_nd)
  !            where
  !               F_nd(is) = F(is)*tref^-1
  !               H_nd(is) = H(is) + A(is)*(ln(tref)-1)
  !               B_nd(is) = B(is)*tref
  !               C_nd(is) = C(is)*tref^2
  !               D_nd(is) = D(is)*tref^3
  !               E_nd(is) = E(is)*tref^4
  !               A_nd(is) = A(is)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: t(:)
  real(wp), intent(inout) :: gbfe(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,ni
  real(wp) :: ts,tlog,tinv
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: intrvl(:)
  !
continue
  !
  ! Compute the temperature interval for each point
  !
  allocate ( intrvl(n1:n2) , source=1 )
  !
  do n = n1,n2
    if (t(n) > temperature_high) then
      intrvl(n) = 2
    else if (t(n) < temperature_low) then
      intrvl(n) = 0
    end if
  end do
  !
  ! Compute the Gibbs free energy / (Runiv*T)
  !
  do m = 1,ncs
    do n = n1,n2
      !
      ni = intrvl(n)
      ts = t(n)
      tlog = log(ts)
      tinv = -one/ts
      !
      !   Compute the (-Gibbs free energy / (R*T)) using the NASA SP-3001 fit
      !
      gbfe(n,m) = gbfe_coef(1,ni,m)*tlog &
                + gbfe_coef(6,ni,m)*tinv &
                + gbfe_coef(7,ni,m) &
                + ts*(gbfe_coef(2,ni,m) + &
                      ts*(gbfe_coef(3,ni,m) + &
                          ts*(gbfe_coef(4,ni,m) + &
                              ts*(gbfe_coef(5,ni,m)))))
      !
    end do
  end do
  !
  deallocate ( intrvl )
  !
end subroutine compute_gibbs_free_energy
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_dgbdt(n1,n2,t,delsco,dgbdt,dimensional)
  !
  ! Routine to calculate the derivtive of the Gibbs free energy
  ! of a specie with respect to temperature
  !
  ! N1     : starting point
  ! N2     : ending   point
  ! DELSCO : stoichiometric coefficient of the forward reaction -
  !          stoichiometric coefficient of the backward reaction
  ! T      : nondimensional static temperature
  ! DGBDT  : derivative of Gibbs free energy / (R*T) w.r.t. T
  !
  ! NOTE: The form is that of an assumed equilibrium model:
  !
  !       kb(ir,T) = kf(ir,T) / ke(ir,T)
  !
  !       and
  !
  !       ke(ir,T) = (Ro*T)**(-DSCSUM)*exp(-deltaG(ir)/(Ro*T))
  !
  !       where:
  !
  !        d/dT[-deltaG(ir)/(Ro*T)] = sum of the Gibbs free energy times
  !                                   the difference of the stoichiometic
  !                                   coeffifients
  !
  !                                         ns
  !        d/dT[-deltaG(ir)/(Ro*T)] =  d/dT{sum[-g(is)/(Ro*T)*(delsco(is,ir))]}
  !                                         is=1
  !
  !    IF USING A DIMENSIONAL TEMPERATURE (T)
  !        d/dT[-g(is)/(Ro*T)] = F(is)*T^-2 +
  !                              A(is)*T^-1 +
  !                              B(is)      +
  !                              C(is)*T    +
  !                              D(is)*T^2  +
  !                              E(is)*T^3
  !
  !    IF USING A NONDIMENSIONAL TEMPERATURE (T_nd)
  !        d/dT[-g(is)/(Ro*T)] = F_nd(is)*T_nd^-2 +
  !                              A_nd(is)*T_nd^-1 +
  !                              B_nd(is)         +
  !                              C_nd(is)*T_nd    +
  !                              D_nd(is)*T_nd^2  +
  !                              E_nd(is)*T_nd^3
  !            where
  !               F_nd(is) = F(is)*tref^-2
  !               A_nd(is) = A(is)*tref^-1
  !               B_nd(is) = B(is)
  !               C_nd(is) = C(is)*tref
  !               D_nd(is) = D(is)*tref^2
  !               E_nd(is) = E(is)*tref^3
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: t(:)
  real(wp),    intent(in) :: delsco(:,:)
  real(wp), intent(inout) :: dgbdt(:,:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: dimensional
  !
  !.. Local Scalars ..
  integer  :: ir,m,n,ni
  real(wp) :: dgdt,ts,tinv,delta_stoich_coef
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: intrvl(:)
  !
continue
  !
  ! Initialize the derivative array
  !
  do ir = 1,ncr
    do n = n1,n2
      dgbdt(n,ir) = zero
    end do
  end do
  !
  ! Compute the temperature interval for each point
  !
  allocate ( intrvl(n1:n2) , source=1 )
  !
  do n = n1,n2
    if (t(n) > temperature_high) then
      intrvl(n) = 2
    else if (t(n) < temperature_low) then
      intrvl(n) = 0
    end if
  end do
  !
  ! Compute the derivative of [Gibbs free energy / (Runiv*T)] wrt T
  ! for each chemical reaction.
  !
  do ir = 1,ncr
    !
    species_loop: &
    do m = 1,ncs
      !
      delta_stoich_coef = delsco(m,ir)
      !
      ! If (forward stoich coef - backward stoich coef) is zero, there
      ! is no point in evaluating dgdt just to multiply it by zero and
      ! add it to dgbdt, so continue on to the next species.
      !
      if (abs(delta_stoich_coef) < rndoff) cycle species_loop
      !
      do n = n1,n2
        !
        ni = intrvl(n)
        ts = t(n)
        tinv = one/ts
        !
        !       Compute d/dT(-Gibbs free energy / (R*T)) from
        !       a curve fit based on NASA SP-3001
        !
        dgdt = dgbfe_coef(2,ni,m) + &               ! B(ni)
               tinv*(dgbfe_coef(1,ni,m) + &         ! A(ni)*T^-1
                     tinv*dgbfe_coef(6,ni,m)) + &   ! F(ni)*T^-2
               ts*(dgbfe_coef(3,ni,m) + &           ! C(ni)*T
                   ts*(dgbfe_coef(4,ni,m) + &       ! D(ni)*T^2
                       ts*(dgbfe_coef(5,ni,m))))    ! E(ni)*T^3
        !
        dgbdt(n,ir) = dgbdt(n,ir) + delta_stoich_coef*dgdt
        !
      end do
      !
    end do species_loop
    !
  end do
  !
  deallocate ( intrvl )
  !
end subroutine compute_dgbdt
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_cp_low_temperature(cpi)
  !
  ! Routine to calculate the specific heat at constant pressure
  ! of a species at the low temperature limit in order to form
  ! the special low temperature coefficients.
  !
  ! NOTE: At the point this subroutine should be called, the temperature_low
  !       variable should still be dimensional and the cp_coef array should
  !       still contain coefficients for curve fits that require dimensional
  !       temperature as input.
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: cpi(:)
  !
  !.. Local Scalars ..
  integer  :: m,ni
  real(wp) :: t
  !
continue
  !
  ni = 1
  !
  t = temperature_low
  !
  do m = 1,ncs
    !
    cpi(m) = cp_coef(1,ni,m) + &
             t*(cp_coef(2,ni,m) + &
                 t*(cp_coef(3,ni,m) + &
                     t*(cp_coef(4,ni,m) + &
                         t*(cp_coef(5,ni,m)))))
    !
  end do
  !
end subroutine compute_cp_low_temperature
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_enthalpy_low_temperature(hsi)
  !
  ! Routine to calculate the enthalpy of a species at the low temperature
  ! limit in order to form the special low temperature coefficients.
  !
  ! NOTE: At the point this subroutine should be called, the temperature_low
  !       variable should still be dimensional and the enth_coef array should
  !       still contain coefficients for curve fits that require dimensional
  !       temperature as input.
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: hsi(:)
  !
  !.. Local Scalars ..
  integer  :: m,ni
  real(wp) :: t
  !
continue
  !
  ni = 1
  !
  t = temperature_low
  !
  do m = 1,ncs
    !
    hsi(m) = enth_coef(6,ni,m) + &
             t*(enth_coef(1,ni,m) + &
                 t*(enth_coef(2,ni,m) + &
                     t*(enth_coef(3,ni,m) + &
                         t*(enth_coef(4,ni,m) + &
                             t*(enth_coef(5,ni,m))))))
    !
  end do
  !
end subroutine compute_enthalpy_low_temperature
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_gbfe_low_temperature(gbfe)
  !
  ! Routine to calculate the Gibbs free energy of a
  ! species at the low temperature limit in order to
  ! form the special low temperature coefficients.
  !
  ! NOTE: At the point this subroutine should be called, the temperature_low
  !       variable should still be dimensional and the gbfe_coef array should
  !       still contain coefficients for curve fits that require dimensional
  !       temperature as input.
  !
  ! GBFE : Gibbs free energy / (R*T)
  !
  !           -g(is)/(Ro*T) = -F(is)*T^-1 +
  !                            H(is)      +
  !                            B(is)*T    +
  !                            C(is)*T^2  +
  !                            D(is)*T^3  +
  !                            E(is)*T^4  +
  !                            A(is)*(ln(T)-1)
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: gbfe(:)
  !
  !.. Local Scalars ..
  integer  :: m,ni
  real(wp) :: ts,tlm1,tinv
  !
continue
  !
  ! Set the temperature interval to 1 since we want to evaluate the low
  ! temperature curve fit at the low temperature limit
  !
  ni = 1
  !
  ! Compute the (-Gibbs free energy / (R*T)) for
  ! each species using the NASA SP-3001 fit
  !
  ts = temperature_low
  tlm1 = log(ts) - one
  tinv = -one/ts
  !
  do m = 1,ncs
      !
      gbfe(m) = gbfe_coef(1,ni,m)*tlm1 &
              + gbfe_coef(6,ni,m)*tinv &
              + gbfe_coef(7,ni,m) &
              + ts*(gbfe_coef(2,ni,m) + &
                    ts*(gbfe_coef(3,ni,m) + &
                        ts*(gbfe_coef(4,ni,m) + &
                            ts*(gbfe_coef(5,ni,m)))))
      !
  end do
  !
end subroutine compute_gbfe_low_temperature
!
!###############################################################################
!
subroutine hardwired_NCC_coefficients
  !
  real(wp), parameter :: NCC_min_temperature = 300.0_wp
  real(wp), parameter :: NCC_interval_switch = 1000.0_wp
  !
  real(wp), parameter :: CO2_wtmol = 44.00995_wp
  !
  real(wp), parameter :: CO2_cp(1:7,1:2) = reshape( [ 0.02275724E+02_wp, &
                                                      0.09922072E-01_wp, &
                                                     -0.10409113E-04_wp, &
                                                      0.06866686E-07_wp, &
                                                     -0.02117280E-10_wp, &
                                                     -0.04837314E+06_wp, &
                                                      0.10188488E+02_wp, &
                                                      0.04453623E+02_wp, &
                                                      0.03140168E-01_wp, &
                                                     -0.12784105E-05_wp, &
                                                      0.02393996E-08_wp, &
                                                     -0.16690333E-13_wp, &
                                                     -0.04896696E+06_wp, &
                                                     -0.09553959E+01_wp] , &
                                                    [7,2] )
  !
  !
  real(wp), parameter :: CO_wtmol = 28.01055_wp
  !
  real(wp), parameter :: CO_cp(1:7,1:2) = reshape( [ 0.03262451E+02_wp, &
                                                     0.15119409E-02_wp, &
                                                    -0.03881755E-04_wp, &
                                                     0.05581944E-07_wp, &
                                                    -0.02474951E-10_wp, &
                                                    -0.14310539E+05_wp, &
                                                     0.04848897E+02_wp, &
                                                     0.03025078E+02_wp, &
                                                     0.14426885E-02_wp, &
                                                    -0.05630827E-05_wp, &
                                                     0.10185813E-09_wp, &
                                                    -0.06910951E-13_wp, &
                                                    -0.14268350E+05_wp, &
                                                     0.06108217E+02_wp] , &
                                                   [7,2] )
  !
  !
  real(wp), parameter :: H2O_wtmol = 18.01534_wp
  !
  real(wp), parameter :: H2O_cp(1:7,1:2) = reshape( [ 0.03386842E+02_wp, &
                                                      0.03474982E-01_wp, &
                                                     -0.06354696E-04_wp, &
                                                      0.06968581E-07_wp, &
                                                     -0.02506588E-10_wp, &
                                                     -0.03020811E+06_wp, &
                                                      0.02590232E+02_wp, &
                                                      0.02672145E+02_wp, &
                                                      0.03056293E-01_wp, &
                                                     -0.08730260E-05_wp, &
                                                      0.12009964E-09_wp, &
                                                     -0.06391618E-13_wp, &
                                                     -0.02989921E+06_wp, &
                                                      0.06862817E+02_wp] , &
                                                    [7,2] )
  !
  !
  real(wp), parameter :: CH4_wtmol = 16.04303_wp
  !
  real(wp), parameter :: CH4_cp(1:7,1:2) = reshape( [ 0.07787415E+01_wp, &
                                                      0.01747668E+00_wp, &
                                                     -0.02783409E-03_wp, &
                                                      0.03049708E-06_wp, &
                                                     -0.12239307E-10_wp, &
                                                     -0.09825229E+05_wp, &
                                                      0.13722195E+02_wp, &
                                                      0.01683478E+02_wp, &
                                                      0.10237236E-01_wp, &
                                                     -0.03875128E-04_wp, &
                                                      0.06785585E-08_wp, &
                                                     -0.04503423E-12_wp, &
                                                     -0.10080787E+05_wp, &
                                                      0.09623395E+02_wp] , &
                                                    [7,2] )
  !
  !
  real(wp), parameter :: N2_wtmol = 28.01340_wp
  !
  real(wp), parameter :: N2_cp(1:7,1:2) = reshape( [ 0.03298677E+02_wp, &
                                                     0.14082404E-02_wp, &
                                                    -0.03963222E-04_wp, &
                                                     0.05641515E-07_wp, &
                                                    -0.02444854E-10_wp, &
                                                    -0.10208999E+04_wp, &
                                                     0.03950372E+02_wp, &
                                                     0.02926640E+02_wp, &
                                                     0.14879768E-02_wp, &
                                                    -0.05684760E-05_wp, &
                                                     0.10097038E-09_wp, &
                                                    -0.06753351E-13_wp, &
                                                    -0.09227977E+04_wp, &
                                                     0.05980528E+02_wp] , &
                                                   [7,2] )
  !
  !
  real(wp), parameter :: O2_wtmol = 31.99880_wp
  !
  real(wp), parameter :: O2_cp(1:7,1:2) = reshape( [ 0.03212936E+02_wp, &
                                                     0.11274864E-02_wp, &
                                                    -0.05756150E-05_wp, &
                                                     0.13138773E-08_wp, &
                                                    -0.08768554E-11_wp, &
                                                    -0.10052490E+04_wp, &
                                                     0.06034737E+02_wp, &
                                                     0.03697578E+02_wp, &
                                                     0.06135197E-02_wp, &
                                                    -0.12588420E-06_wp, &
                                                     0.01775281E-09_wp, &
                                                    -0.11364354E-14_wp, &
                                                    -0.12339301E+04_wp, &
                                                     0.03189165E+02_wp] , &
                                                   [7,2] )
  !
  real(wp), parameter :: H_wtmol = 1.0079400_wp
  !
  real(wp), parameter :: H_cp(1:7,1:2) = reshape( [ 0.02500000E+02_wp, &
                                                    zero, &
                                                    zero, &
                                                    zero, &
                                                    zero, &
                                                    0.02547162E+06_wp, &
                                                   -0.04601176E+01_wp, &
                                                    0.02500000E+02_wp, &
                                                    zero, &
                                                    zero, &
                                                    zero, &
                                                    zero, &
                                                    0.02547162E+06_wp, &
                                                   -0.04601176E+01_wp] , &
                                                  [7,2] )
  !
  real(wp), parameter :: H2_wtmol = 2.0158800_wp
  !
  real(wp), parameter :: H2_cp(1:7,1:2) = reshape( [ 0.03298124E+02_wp, &
                                                     0.08249441E-02_wp, &
                                                    -0.08143015E-05_wp, &
                                                    -0.09475434E-09_wp, &
                                                     0.04134872E-11_wp, &
                                                    -0.10125209E+04_wp, &
                                                    -0.03294094E+02_wp, &
                                                     0.02991423E+02_wp, &
                                                     0.07000644E-02_wp, &
                                                    -0.05633828E-06_wp, &
                                                    -0.09231578E-10_wp, &
                                                     0.15827519E-14_wp, &
                                                    -0.08350340E+04_wp, &
                                                    -0.13551101E+01_wp] , &
                                                   [7,2] )
  !
  real(wp), parameter :: O_wtmol = 15.9994000_wp
  !
  real(wp), parameter :: O_cp(1:7,1:2) = reshape( [ 0.02946428E+02_wp, &
                                                   -0.16381665E-02_wp, &
                                                    0.02421031E-04_wp, &
                                                   -0.16028431E-08_wp, &
                                                    0.03890696E-11_wp, &
                                                    0.02914764E+06_wp, &
                                                    0.02963995E+02_wp, &
                                                    0.02542059E+02_wp, &
                                                   -0.02755061E-03_wp, &
                                                   -0.03102803E-07_wp, &
                                                    0.04551067E-10_wp, &
                                                   -0.04368051E-14_wp, &
                                                    0.02923080E+06_wp, &
                                                    0.04920308E+02_wp] , &
                                                  [7,2] )
  !
  real(wp), parameter :: OH_wtmol = 17.0073400_wp
  !
  real(wp), parameter :: OH_cp(1:7,1:2) = reshape( [ 0.03637266E+02_wp, &
                                                     0.01850910E-02_wp, &
                                                    -0.16761646E-05_wp, &
                                                     0.02387202E-07_wp, &
                                                    -0.08431442E-11_wp, &
                                                     0.03606781E+05_wp, &
                                                     0.13588605E+01_wp, &
                                                     0.02882730E+02_wp, &
                                                     0.10139743E-02_wp, &
                                                    -0.02276877E-05_wp, &
                                                     0.02174683E-09_wp, &
                                                    -0.05126305E-14_wp, &
                                                     0.03886888E+05_wp, &
                                                     0.05595712E+02_wp] , &
                                                   [7,2] )
  !
  !.. Local Scalars ..
  integer  :: i,j,m,ni
  real(wp) :: t,tlog
  !
  character(len=len(species_name)) :: species
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: ncc_cp(:,:,:)
  real(wp), allocatable :: tmp_coef(:)
  real(wp), allocatable :: cpi(:)
  real(wp), allocatable :: hsi(:)
  real(wp), allocatable :: gbfe(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "hardwired_NCC_coefficients"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  allocate ( ncc_cp(1:7,1:2,1:ncs) , source=zero )
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
        ncc_cp(1:7,1:2,m) = CO2_cp(1:7,1:2)
        !
      case ("CO")
        !
        wtmol(m) = CO_wtmol
        ncc_cp(1:7,1:2,m) = CO_cp(1:7,1:2)
        !
      case ("H2O")
        !
        wtmol(m) = H2O_wtmol
        ncc_cp(1:7,1:2,m) = H2O_cp(1:7,1:2)
        !
      case ("CH4")
        !
        wtmol(m) = CH4_wtmol
        ncc_cp(1:7,1:2,m) = CH4_cp(1:7,1:2)
        !
      case ("N2")
        !
        wtmol(m) = N2_wtmol
        ncc_cp(1:7,1:2,m) = N2_cp(1:7,1:2)
        !
      case ("O2")
        !
        wtmol(m) = O2_wtmol
        ncc_cp(1:7,1:2,m) = O2_cp(1:7,1:2)
        !
      case ("H")
        !
        wtmol(m) = H_wtmol
        ncc_cp(1:7,1:2,m) = H_cp(1:7,1:2)
        !
      case ("H2")
        !
        wtmol(m) = H2_wtmol
        ncc_cp(1:7,1:2,m) = H2_cp(1:7,1:2)
        !
      case ("O")
        !
        wtmol(m) = O_wtmol
        ncc_cp(1:7,1:2,m) = O_cp(1:7,1:2)
        !
      case ("OH")
        !
        wtmol(m) = OH_wtmol
        ncc_cp(1:7,1:2,m) = OH_cp(1:7,1:2)
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
  if (ncpc == 5) then
    tmp_coef = one / [one,two,three,four,five,one] ! F2003 AUTO-REALLOC
  end if
  !
  do m = 1,ncs
    do i = 1,size(ncc_cp,dim=2)
      do j = 1,ncpc+1
        enth_coef(j,i,m) = ncc_cp(j,i,m) * rspec(m) * tmp_coef(j)
      end do
    end do
  end do
  !
  if (ncpc == 5) then
    tmp_coef = one / [one,two,six,twelve,two*ten,one,one] ! F2003 AUTO-REALLOC
  end if
  !
  do m = 1,ncs
    do i = 1,size(ncc_cp,dim=2)
      do j = 1,ncpc+2
        gbfe_coef(j,i,m) = ncc_cp(j,i,m) * tmp_coef(j)
      end do
    end do
  end do
  !
  if (ncpc == 5) then
    tmp_coef = [one,one,two,three,four,one] ! F2003 AUTO-REALLOC
  end if
  !
  do m = 1,ncs
    do i = 1,size(ncc_cp,dim=2)
      do j = 1,ncpc+1
        dgbfe_coef(j,i,m) = gbfe_coef(j,i,m) * tmp_coef(j)
      end do
    end do
  end do
  !
  temperature_low  = NCC_min_temperature
  temperature_high = NCC_interval_switch
  !
  ! Build special low temperature coefficients
  ! for each of the thermodynamic curve fits
  !
  allocate ( cpi(1:ncs) , source=zero )
  allocate ( hsi(1:ncs) , source=zero )
  allocate ( gbfe(1:ncs) , source=zero )
  !
  call compute_cp_low_temperature(cpi)
  call compute_enthalpy_low_temperature(hsi)
  call compute_gbfe_low_temperature(gbfe)
  !
  ni = 0
  !
  t = temperature_low
  tlog = log(t)
  !
  do m = 1,ncs
    !
    ! Null the coefficients for Cp(T)
    !
    do j = 1,size(cp_coef,dim=1)
      cp_coef(j,ni,m) = zero
    end do
    !
    ! Null coefficients for Hs(T)
    !
    do j = 1,size(enth_coef,dim=1)
      enth_coef(j,ni,m) = zero
    end do
    !
    ! Null coefficients for Gb(T)
    !
    do j = 1,size(gbfe_coef,dim=1)
      gbfe_coef(j,ni,m) = zero
    end do
    !
    ! Set coefficients for d/dT[Gb(T)] = coefficients from
    ! low temperature interval of d/dT[Gb(T)]
    !
    do j = 1,size(dgbfe_coef,dim=1)
      dgbfe_coef(j,ni,m) = dgbfe_coef(j,ni+1,m)
    end do
    !
    ! Compute special low temperature coefficients for Cp(T) , Hs(T) , and Gb(T)
    !
    j = 1
    !
    cp_coef(j,ni,m) = cpi(m)
    enth_coef(j,ni,m) = cpi(m)
    gbfe_coef(j,ni,m) = cpi(m)/ rspec(m)
    !
    j = size(enth_coef,dim=1)
    enth_coef(j,ni,m) = hsi(m)- cpi(m)*t
    !
    j = size(gbfe_coef,dim=1)-1
    gbfe_coef(j,ni,m) = (hsi(m)-cpi(m)*t) / rspec(m)
    !
    j = size(gbfe_coef,dim=1)
    gbfe_coef(j,ni,m) = gbfe(m) + (hsi(m)-cpi(m)*t*tlog) / (rspec(m)*t)
    !
  end do
  !
  if (allocated(gbfe)) deallocate ( gbfe )
  if (allocated(hsi)) deallocate ( hsi )
  if (allocated(cpi)) deallocate ( cpi )
  if (allocated(ncc_cp)) deallocate ( ncc_cp )
  if (allocated(tmp_coef)) deallocate ( tmp_coef )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine hardwired_NCC_coefficients
!
!###############################################################################
!
subroutine hardwired_VULCAN_coefficients
  !
  real(wp), parameter :: vul_min_temperature = 300.0_wp
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
  !.. Local Scalars ..
  integer  :: i,j,m,ni
  real(wp) :: t,tlog
  !
  character(len=len(species_name)) :: species
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: vul_cp(:,:,:)
  real(wp), allocatable :: tmp_coef(:)
  real(wp), allocatable :: cpi(:)
  real(wp), allocatable :: hsi(:)
  real(wp), allocatable :: gbfe(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "hardwired_VULCAN_coefficients"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
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
  if (ncpc == 5) then
    tmp_coef = one / [one,two,three,four,five,one] ! F2003 AUTO-REALLOC
  end if
  !
  do m = 1,ncs
    do i = 1,size(vul_cp,dim=2)
      do j = 1,ncpc+1
        enth_coef(j,i,m) = vul_cp(j,i,m) * rspec(m) * tmp_coef(j)
      end do
    end do
  end do
  !
  if (ncpc == 5) then
    tmp_coef = one / [one,two,six,twelve,two*ten,one,one] ! F2003 AUTO-REALLOC
  end if
  !
  do m = 1,ncs
    do i = 1,size(vul_cp,dim=2)
      do j = 1,ncpc+2
        gbfe_coef(j,i,m) = vul_cp(j,i,m) * tmp_coef(j)
      end do
    end do
  end do
  !
  if (ncpc == 5) then
    tmp_coef = [one,one,two,three,four,one] ! F2003 AUTO-REALLOC
  end if
  !
  do m = 1,ncs
    do i = 1,size(vul_cp,dim=2)
      do j = 1,ncpc+1
        dgbfe_coef(j,i,m) = gbfe_coef(j,i,m) * tmp_coef(j)
      end do
    end do
  end do
  !
  temperature_low  = vul_min_temperature
  temperature_high = vul_interval_switch
  !
  ! Build special low temperature coefficients
  ! for each of the thermodynamic curve fits
  !
  allocate ( cpi(1:ncs) , source=zero )
  allocate ( hsi(1:ncs) , source=zero )
  allocate ( gbfe(1:ncs) , source=zero )
  !
  call compute_cp_low_temperature(cpi)
  call compute_enthalpy_low_temperature(hsi)
  call compute_gbfe_low_temperature(gbfe)
  !
  ni = 0
  !
  t = temperature_low
  tlog = log(t)
  !
  do m = 1,ncs
    !
    ! Null the coefficients for Cp(T)
    !
    do j = 1,size(cp_coef,dim=1)
      cp_coef(j,ni,m) = zero
    end do
    !
    ! Null coefficients for Hs(T)
    !
    do j = 1,size(enth_coef,dim=1)
      enth_coef(j,ni,m) = zero
    end do
    !
    ! Null coefficients for Gb(T)
    !
    do j = 1,size(gbfe_coef,dim=1)
      gbfe_coef(j,ni,m) = zero
    end do
    !
    ! Set coefficients for d/dT[Gb(T)] = coefficients from
    ! low temperature interval of d/dT[Gb(T)]
    !
    do j = 1,size(dgbfe_coef,dim=1)
      dgbfe_coef(j,ni,m) = dgbfe_coef(j,ni+1,m)
    end do
    !
    ! Compute special low temperature coefficients for Cp(T) , Hs(T) , and Gb(T)
    !
    j = 1
    !
    cp_coef(j,ni,m) = cpi(m)
    enth_coef(j,ni,m) = cpi(m)
    gbfe_coef(j,ni,m) = cpi(m)/ rspec(m)
    !
    j = size(enth_coef,dim=1)
    enth_coef(j,ni,m) = hsi(m)- cpi(m)*t
    !
    j = size(gbfe_coef,dim=1)-1
    gbfe_coef(j,ni,m) = (hsi(m)-cpi(m)*t) / rspec(m)
    !
    j = size(gbfe_coef,dim=1)
    gbfe_coef(j,ni,m) = gbfe(m) + (hsi(m)-cpi(m)*t*tlog) / (rspec(m)*t)
    !
  end do
  !
#ifdef DEBUG_ON
#ifdef DEBUG_THERMO_ON
  if (mypnum == glb_root) then
    do m = 1,ncs
      !
      write (610,50) trim(adjustl(species_name(m)))
      !
      write (610,53) "cpi ",cpi(m)
      write (610,53) "hsi ",hsi(m)
      write (610,53) "gpfe",gbfe(m)
      !
      write (610,52) "cp_coef"
      do ni = lbound(cp_coef,dim=2),ubound(cp_coef,dim=2)
        write (610,51) (cp_coef(j,ni,m), j=1,size(cp_coef,dim=1))
      end do
      !
      write (610,52) "enth_coef"
      do ni = lbound(enth_coef,dim=2),ubound(enth_coef,dim=2)
        write (610,51) (enth_coef(j,ni,m), j=1,size(enth_coef,dim=1))
      end do
      !
      write (610,52) "gbfe_coef"
      do ni = lbound(gbfe_coef,dim=2),ubound(gbfe_coef,dim=2)
        write (610,51) (gbfe_coef(j,ni,m), j=1,size(gbfe_coef,dim=1))
      end do
      !
      write (610,52) "dgbfe_coef"
      do ni = lbound(dgbfe_coef,dim=2),ubound(dgbfe_coef,dim=2)
        write (610,51) (dgbfe_coef(j,ni,m), j=1,size(dgbfe_coef,dim=1))
      end do
      !
    end do
    flush (610)
  end if
  50 format (//,6x,a,/)
  51 format (4x,10es14.6)
  52 format (/,a,:," = ",es14.6)
  53 format (a," = ",es14.6)
#endif
#endif
  !
  if (allocated(gbfe)) deallocate ( gbfe )
  if (allocated(hsi)) deallocate ( hsi )
  if (allocated(cpi)) deallocate ( cpi )
  if (allocated(vul_cp)) deallocate ( vul_cp )
  if (allocated(tmp_coef)) deallocate ( tmp_coef )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine hardwired_VULCAN_coefficients
!
!###############################################################################
!
end module thermo_mod
