module reactions_mod
  !
  !.. Use Statements ..
  use module_kind_types
  !
  use eqn_idx_mod, only : nec,nee,nq
  use eqn_idx_mod, only : ncs,ntemp,ngamm,ncp,ng
  use eqn_idx_mod, only : ncr
  !
  use thermo_mod, only : runivl_millijoules => runivl
  !
  implicit none
  !
  private
  !
  public :: init_reactions_mod
  public :: compute_wdot
  !
  !
  logical(lk), save, allocatable :: nonreacting_species(:)
  !
  integer, save, allocatable :: reactant_species_list(:)
  integer, save, allocatable :: product_species_list(:)
  !
  ! Reaction arrays for CH4
  !
  real(wp), save, allocatable :: fwd_arrh_coef(:)
  real(wp), save, allocatable :: fwd_actv_temp(:)
  real(wp), save, allocatable :: fwd_exp(:,:)
  !
  real(wp), save, allocatable :: bwd_arrh_coef(:)
  real(wp), save, allocatable :: bwd_actv_temp(:)
  real(wp), save, allocatable :: bwd_exp(:,:)
  !
  real(wp), save, allocatable :: fwd_stoich_coef(:,:)
  real(wp), save, allocatable :: bwd_stoich_coef(:,:)
  real(wp), save, allocatable :: stoich_coef(:,:)
  !
  ! Reaction arrays for H2
  !
  real(wp), save, allocatable :: scoerp(:,:,:) ! scoerp(1:2,1:ncs+1,1:ncr)
  real(wp), save, allocatable :: frcoef(:,:)   ! frcoef(1:3,1:ncr)
  real(wp), save, allocatable :: cbeffc(:,:)   ! cbeffc(1:ncs,1:ncr)
  !
  ! N_REAC_SUB_ITER: # of sub-iterations for computing
  !                  source terms from reaction kinetics
  !
  integer, save :: n_reac_sub_iter = 10
  !
  !
  !
  integer, save :: nCO2  = 0
  integer, save :: nCO   = 0
  integer, save :: nH2O  = 0
  integer, save :: nCH4  = 0
  integer, save :: nC2H4 = 0
  integer, save :: nN2   = 0
  integer, save :: nO2   = 0
  !
  integer, save :: nH2   = 0
  integer, save :: nOH   = 0
  integer, save :: nO    = 0
  integer, save :: nH    = 0
  !
  ! nHP : Hydrodgen product from CH4 reaction.
  !
  !       If using CH4 reaction mechanism, the CH4 reaction equation is either
  !
  !          2*CH4 + 1*O2 => 2*CO + 4*H2
  !                      or
  !          2*CH4 + 3*O2 => 2*CO + 4*H2O
  !
  !       If using the first  equation, nHP = nH2
  !       If using the second equation, nHP = nH2O
  !
  integer, save :: nHP  = 0
  !
  !
  !
  ! ###########################
  ! ###  Module Parameters  ###
  ! ###########################
  !
  real(wp), parameter :: dens_cgs_to_si = 1000.0_wp
  real(wp), parameter :: dens_si_to_cgs = one/dens_cgs_to_si
  !
  real(wp), parameter :: joules_per_calorie = 4.184_wp
  real(wp), parameter :: calories_per_joule = one/joules_per_calorie
  !
  real(wp), parameter :: runivl_calories = runivl_millijoules &
                                         * calories_per_joule &
                                         / 1000.0_wp ! convert mJ to J
  !
contains
!
!###############################################################################
!
subroutine init_reactions_mod()
  !
  !.. Use Statements ..
  use ovar, only : itestcase
  use ovar, only : chemistry_sub_iterations
  !
  use thermo_mod, only : species_name
  !
  !.. Local Scalars ..
  integer :: m,ierr
  !
  character(len=10) :: species
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_reactions_mod"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (ncr > 0) then
    !
    if (chemistry_sub_iterations > 0) then
      n_reac_sub_iter = chemistry_sub_iterations
    end if
    !
    do m = 1,ncs
      !
      species = uppercase( trim(adjustl(species_name(m))) )
      !
      select case (trim(adjustl(species)))
        !
        case ("CO2")
          nCO2 = m
        case ("CO")
          nCO  = m
        case ("H2O")
          nH2O = m
        case ("CH4")
          nCH4 = m
        case ("C2H4")
          nC2H4 = m
        case ("N2")
          nN2  = m
        case ("O2")
          nO2  = m
        case ("O")
          nO   = m
        case ("H2")
          nH2  = m
        case ("H")
          nH   = m
        case ("OH")
          nOH  = m
          !
      end select
      !
    end do
    !
    if (.not. allocated(nonreacting_species)) then
      allocate ( nonreacting_species(1:ncs) , source=true , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"nonreacting_species",1,__LINE__,__FILE__, &
                       ierr,error_message)
    end if
    !
    if (abs(itestcase) == CH4_3D_Diffusion_Flame .or. &
        abs(itestcase) == CH4_Flame_Speed_Test) then
      !
      if (nH2O /= 0 .and. nH2 == 0) then
        nHP = nH2O
      else if (nH2 /= 0 .and. nH2O == 0) then
        nHP = nH2
      end if
      !
      !----------------------------------------
      ! Species:  CO2  CO  H2O  CH4  N2  O2
      !----------------------------------------
      ! Reaction 1:  2CH4 + 3O2 => 2CO  + 4H2O
      ! Reaction 2:  2CO2       => 2CO  +  O2
      ! Reaction 3:  2CO  +  O2 => 2CO2
      !----------------------------------------
      !
      if (itestcase > 0) then
        call initialize_CH4_6x3_reactions
        call check_for_bad_number(__LINE__,__FILE__)
      else
        call initialize_CH4_6x3_Vulcan_reactions
        call check_for_bad_number(__LINE__,__FILE__)
      end if
      !
      ! Create the lists of reactant and product species
      !
      reactant_species_list = [nCH4,nO2,nN2] ! F2003 AUTO REALLOCATION
      product_species_list  = [nCO2,nCO,nHP] ! F2003 AUTO REALLOCATION
      !
    else if ( itestcase == H2_2D_Diffusion_Flame .or. &
              itestcase == Compression_Corner ) then
      !
      !-----------------------------------------------
      ! Species:  O2  H2  H2O  OH  O  H  N2
      !-----------------------------------------------
      ! Reaction 1:   O2 + H2      <=>  2OH
      ! Reaction 2:   O2 + H       <=>   OH  +  O
      ! Reaction 3:   H2 + OH      <=>   H2O +  H
      ! Reaction 4:   H2 + O       <=>   OH  +  H
      ! Reaction 5:  2OH           <=>   H2O +  O
      ! Reaction 6:   OH + H  + M  <=>   H2O +  M
      ! Reaction 7:   H  + H  + M  <=>   H2  +  M
      !-----------------------------------------------
      !
      call initialize_H2_7x7_reactions
      !
      ! Create the lists of reactant and product species
      !
      reactant_species_list = [nO2,nH2,nN2] ! F2003 AUTO REALLOCATION
      product_species_list  = [nOH,nH2O,nO,nH] ! F2003 AUTO REALLOCATION
      !
    else if ( itestcase == UIUC_Scramjet) then
      !
      !-----------------------------------------------
      ! Species:  O2  H2  H2O  C2H4  N2 CO CO2
      !-----------------------------------------------
      ! Reaction 1:   C2H4 + 02    <=>  2CO + 2H2
      ! Reaction 2:   2CO + O2     <=>  2CO2
      ! Reaction 3:   2H2 + O2     <=>  2H2O
      !-----------------------------------------------
      !
      call initialize_Ethylene_7x3_reactions
      !
      ! Create the lists of reactant and product species
      !
      reactant_species_list = [nO2,nC2H4,nN2] ! F2003 AUTO REALLOCATION
      product_species_list  = [nCO,nCO2,nH2O,nH2] ! F2003 AUTO REALLOCATION
      !
    else
      !
      write (error_message,1)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      !
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  1 format ("ERROR! Test case is not valid with chemical reactions! ERROR!")
  !
end subroutine init_reactions_mod
!
!###############################################################################
!
subroutine initialize_CH4_6x3_reactions
  !
  !.. Use Statements ..
  use ovar, only : rhoref,aref,tref
  use thermo_mod, only : wtmol
  use thermo_mod, only : species_name
  !
  !.. Local Scalars ..
  integer  :: i,j,ir,ierr
  real(wp) :: fwd_mwp,bwd_mwp,c1,c2,c3,val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_CH4_6x3_reactions"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  c1 = one
 !c1 = half
  c2 = one/c1
  c3 = c1
  !
  ! Allocate all the arrays needed for CH4 reactions
  !
  if (.not. allocated(fwd_arrh_coef)) then
    allocate ( fwd_arrh_coef(1:ncr) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"fwd_arrh_coef",1,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (.not. allocated(bwd_arrh_coef)) then
    allocate ( bwd_arrh_coef(1:ncr) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"bwd_arrh_coef",1,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  !
  !
  if (.not. allocated(fwd_actv_temp)) then
    allocate ( fwd_actv_temp(1:ncr) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"fwd_actv_temp",1,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (.not. allocated(bwd_actv_temp)) then
    allocate ( bwd_actv_temp(1:ncr) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"bwd_actv_temp",1,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  !
  !
  if (.not. allocated(fwd_exp)) then
    allocate ( fwd_exp(1:ncs,1:ncr) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"fwd_exp",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (.not. allocated(bwd_exp)) then
    allocate ( bwd_exp(1:ncs,1:ncr) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"bwd_exp",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  !
  !
  if (.not. allocated(fwd_stoich_coef)) then
    allocate ( fwd_stoich_coef(1:ncs,1:ncr) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"fwd_stoich_coef",1,__LINE__,__FILE__, &
                   ierr,error_message)
  end if
  !
  if (.not. allocated(bwd_stoich_coef)) then
    allocate ( bwd_stoich_coef(1:ncs,1:ncr) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"bwd_stoich_coef",1,__LINE__,__FILE__, &
                   ierr,error_message)
  end if
  !
  if (.not. allocated(stoich_coef)) then
    allocate ( stoich_coef(1:ncs,1:ncr) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"stoich_coef",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  !##############################################
  !
  ! Species:  CO2  CO  H2O  CH4  N2  O2
  !
  !##############################################
  !
  ! Reaction 1:  2CH4 + 3O2 => 2CO + 4H2O
  !
  ir = 1
  !
  fwd_arrh_coef(ir) = 1.0e+15_wp * c2 ! cgs units
  fwd_actv_temp(ir) = 35000.0_wp * c3 ! [calories/mole]
  !
  bwd_arrh_coef(ir) = zero
  bwd_actv_temp(ir) = zero
  !
  fwd_exp(nCH4,ir) = 0.9_wp ! CH4
  fwd_exp(nO2 ,ir) = 1.1_wp ! O2
  !
  if (nH2O == nHP) then
    fwd_stoich_coef(nCH4,ir) = two   * c1 ! CH4
    fwd_stoich_coef(nO2 ,ir) = three * c1 ! O2
    !
    bwd_stoich_coef(nCO ,ir) = two   * c1 ! CO
    bwd_stoich_coef(nH2O,ir) = four  * c1 ! H2O
  else if (nH2 == nHP) then
    fwd_stoich_coef(nCH4,ir) = two   * c1 ! CH4
    fwd_stoich_coef(nO2 ,ir) = one   * c1 ! O2
    !
    bwd_stoich_coef(nCO ,ir) = two   * c1 ! CO
    bwd_stoich_coef(nH2 ,ir) = four  * c1 ! H2
  end if
  !
  !##############################################
  !
  ! Reaction 2:  2CO2 => 2CO + O2
  !
  ir = 2
  !
  fwd_arrh_coef(ir) = 4.0552e+10_wp * c2 ! cgs units
  fwd_actv_temp(ir) = 77194.0_wp    * c3 ! [calories/mole]
  !
  bwd_arrh_coef(ir) = zero
  bwd_actv_temp(ir) = zero
  !
  fwd_exp(nCO2,ir) = one ! CO2
  !
  fwd_stoich_coef(nCO2,ir) = two * c1 ! CO2
  !
  bwd_stoich_coef(nCO ,ir) = two * c1 ! CO
  bwd_stoich_coef(nO2 ,ir) = one * c1 ! O2
  !
  !##############################################
  !
  ! Reaction 3:  2CO + O2 => 2CO2
  !
  ir = 3
  !
  fwd_arrh_coef(ir) = 1.0e+9_wp  * c2 ! cgs units
  fwd_actv_temp(ir) = 12000.0_wp * c3 ! [calories/mole]
  !
  bwd_arrh_coef(ir) = zero
  bwd_actv_temp(ir) = zero
  !
  fwd_exp(nCO ,ir) = one  ! CO
  fwd_exp(nO2 ,ir) = half ! O2
  !
  fwd_stoich_coef(nCO ,ir) = two * c1 ! CO
  fwd_stoich_coef(nO2 ,ir) = one * c1 ! O2
  !
  bwd_stoich_coef(nCO2,ir) = two * c1 ! CO2
  !
  call check_for_bad_number(__LINE__,__FILE__)
  !
  !################################################
  !
  ! Pre-multiply the forwards and backwards Arrhenius coefficients by
  !     product( (dens_si_to_cgs*rhoref/wtmol(1:ncs))**fwd_exp(1:ncs,j) ) / aref
  !     product( (dens_si_to_cgs*rhoref/wtmol(1:ncs))**bwd_exp(1:ncs,j) ) / aref
  !
  do j = 1,ncr
    fwd_mwp = one
    bwd_mwp = one
    do i = 1,ncs
      val = dens_si_to_cgs*rhoref/wtmol(i)
      fwd_mwp = fwd_mwp * (val**fwd_exp(i,j))
      bwd_mwp = bwd_mwp * (val**bwd_exp(i,j))
    end do
    fwd_arrh_coef(j) = fwd_arrh_coef(j)*fwd_mwp / aref
    bwd_arrh_coef(j) = bwd_arrh_coef(j)*bwd_mwp / aref
  end do
  !
  call check_for_bad_number(__LINE__,__FILE__)
  !
  ! Pre-multiply the stoichiometric coefficients by
  !  the molar mass * [density conversion from kg/m^3 to g/cm^3] * rhoref
  !
  do j = 1,ncr
    do i = 1,ncs
      val = wtmol(i)*dens_cgs_to_si/rhoref
      fwd_stoich_coef(i,j) = fwd_stoich_coef(i,j) * val
      bwd_stoich_coef(i,j) = bwd_stoich_coef(i,j) * val
    end do
  end do
  !
  ! Nondimensionalize the activation temperatures
  ! NOTE: Currently the activation temperatures are actually the
  !       activation energies in calories/mole, so divide them by
  !       the universal gas constant in calories/(mole*Kelvin) to
  !       get the activation temperatures.
  !
  do j = 1,ncr
    fwd_actv_temp(j) = fwd_actv_temp(j)/tref/runivl_calories
    bwd_actv_temp(j) = bwd_actv_temp(j)/tref/runivl_calories
  end do
  !
  ! Compute the net stoichiometric coefficient of each species for each reaction
  !
  do j = 1,ncr
    do i = 1,ncs
      stoich_coef(i,j) = bwd_stoich_coef(i,j) - fwd_stoich_coef(i,j)
    end do
  end do
  !
  !##############################################
  !
  ! Identify any species that are non-reacting. All species are initially
  ! marked as non-reacting and are marked reacting if either the forward
  ! or backward stoichiometric coefficents for the species are non-zero
  ! for any of the chemical reactions.
  !
  do i = 1,ncs
    nonreacting_species(i) = all(abs(fwd_stoich_coef(i,1:ncr)) < eps6) .and. &
                             all(abs(bwd_stoich_coef(i,1:ncr)) < eps6)
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine initialize_CH4_6x3_reactions
!
!###############################################################################
!
subroutine initialize_CH4_6x3_Vulcan_reactions
  !
  !.. Use Statements ..
  use thermo_mod, only : species_name
  !
  !.. Local Scalars ..
  integer  :: m,ir,ierr
  real(wp) :: c1,c2,c3
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_CH4_6x3_Vulcan_reactions"
  !
continue
  !
  c1 = one
 !c1 = half
  c2 = one/c1
  c3 = c1
 !c3 = one
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Allocate all the arrays needed for H2 reactions
  ! NOTE: cbeffc is initialized to 1 instead of 0
  !
  if (.not. allocated(frcoef)) then
    allocate ( frcoef(1:3,1:ncr) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"frcoef",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (.not. allocated(scoerp)) then
    allocate ( scoerp(1:2,1:ncs+1,1:ncr) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"scoerp",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (.not. allocated(fwd_exp)) then
    allocate ( fwd_exp(1:ncs,1:ncr) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"fwd_exp",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (.not. allocated(bwd_exp)) then
    allocate ( bwd_exp(1:ncs,1:ncr) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"bwd_exp",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (.not. allocated(cbeffc)) then
    allocate ( cbeffc(1:ncs,1:ncr) , source=one , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cbeffc",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  !##############################################
  !
  ! Species:  CH4 CO2 CO H2O O2 N2
  !
  !##############################################
  !
  ! Reaction 1:  2CH4 + 3O2 => 2CO + 4H2O
  !
  ir = 1
  !
  if (nH2O == nHP) then
    scoerp(1,nCH4,ir) = two    * c1
    scoerp(1,nO2 ,ir) = three  * c1
    !
    scoerp(2,nCO ,ir) = two    * c1
    scoerp(2,nH2O,ir) = four   * c1
  else if (nH2 == nHP) then
    scoerp(1,nCH4,ir) = two    * c1
    scoerp(1,nO2 ,ir) = one    * c1
    !
    scoerp(2,nCO ,ir) = two    * c1
    scoerp(2,nH2 ,ir) = four   * c1
  end if
  !
  fwd_exp(nCH4,ir) = 0.9_wp
  fwd_exp(nO2 ,ir) = 1.1_wp
  !
  frcoef(1,ir) = 1.0e15_wp  * c2 ! cgs units
  frcoef(2,ir) = zero
  frcoef(3,ir) = 35000.0_wp * c3 ! [calories/mole]
  !
  !##############################################
  !
  ! Reaction 2:  2CO2 => 2CO + O2
  !
  ir = 2
  !
  scoerp(1,nCO2,ir) = two * c1
  !
  scoerp(2,nCO ,ir) = two * c1
  scoerp(2,nO2 ,ir) = one * c1
  !
  fwd_exp(nCO2,ir) = one
  !
  frcoef(1,ir) = 4.0552e10_wp * c2 ! cgs units
  frcoef(2,ir) = zero
  frcoef(3,ir) = 77194.0_wp   * c3 ! [calories/mole]
  !
  !##############################################
  !
  ! Reaction 3:  2CO + O2 => 2CO2
  !
  ir = 3
  !
  scoerp(1,nCO ,ir) = two * c1
  scoerp(1,nO2 ,ir) = one * c1
  !
  scoerp(2,nCO2,ir) = two * c1
  !
  fwd_exp(nCO ,ir) = one
  fwd_exp(nO2 ,ir) = half
  !
  frcoef(1,ir) = 1.0e9_wp   * c2 ! cgs units
  frcoef(2,ir) = zero
  frcoef(3,ir) = 12000.0_wp * c3 ! [calories/mole]
  !
  !##############################################
  !
  do ir = 1,ncr
    frcoef(1,ir) = log( frcoef(1,ir) )
  end do
  !
  ! Currently the activation temperature are actually the activation energies
  ! in calories/mole, so divide them by the universal gas constant in
  ! calories/(mole*Kelvin) to get the activation temperatures.
  !
  do ir = 1,ncr
    frcoef(3,ir) = frcoef(3,ir)/runivl_calories
  end do
  !
  !##############################################
  !
  ! Identify any species that are non-reacting. All species are initially
  ! marked as non-reacting and are marked reacting if either the forward
  ! or backward stoichiometric coefficents for the species are non-zero
  ! for any of the chemical reactions.
  !
  do m = 1,ncs
    nonreacting_species(m) = all(abs(scoerp(1:2,m,1:ncr)) < eps6)
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine initialize_CH4_6x3_Vulcan_reactions
!
!###############################################################################
!
subroutine initialize_H2_7x7_reactions
  !
  !.. Use Statements ..
  use thermo_mod, only : species_name
  !
  !.. Local Scalars ..
  integer :: m,ir,ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_H2_7x7_reactions"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Allocate all the arrays needed for H2 reactions
  ! NOTE: cbeffc is initialized to 1 instead of 0
  !
  if (.not. allocated(frcoef)) then
    allocate ( frcoef(1:3,1:ncr) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"frcoef",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (.not. allocated(scoerp)) then
    allocate ( scoerp(1:2,1:ncs+1,1:ncr) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"scoerp",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (.not. allocated(fwd_exp)) then
    allocate ( fwd_exp(1:ncs,1:ncr) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"fwd_exp",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (.not. allocated(bwd_exp)) then
    allocate ( bwd_exp(1:ncs,1:ncr) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"bwd_exp",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (.not. allocated(cbeffc)) then
    allocate ( cbeffc(1:ncs,1:ncr) , source=one , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cbeffc",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  !##############################################
  !
  ! Species:  H2  O2  H2O  OH  H  O  N2
  !
  !##############################################
  !
  ! Reaction 1:  H2 + O2 <=> 2OH
  !
  ir = 1
  !
  scoerp(1,nH2,ir) = one
  scoerp(1,nO2,ir) = one
  !
  scoerp(2,nOH,ir) = two
  !
  frcoef(1,ir) = 1.70e13_wp
  frcoef(2,ir) = zero
  frcoef(3,ir) = 24154.6_wp
  !
  !##############################################
  !
  ! Reaction 2:  H + O2 <=> OH + O
  !
  ir = 2
  !
  scoerp(1,nH ,ir) = one
  scoerp(1,nO2,ir) = one
  !
  scoerp(2,nOH,ir) = one
  scoerp(2,nO ,ir) = one
  !
  frcoef(1,ir) = 1.20e17_wp
  frcoef(2,ir) = -0.91_wp
  frcoef(3,ir) = 8309.7_wp
  !
  !##############################################
  !
  ! Reaction 3:  OH + H2 <=> H2O + H
  !
  ir = 3
  !
  scoerp(1,nOH ,ir) = one
  scoerp(1,nH2 ,ir) = one
  !
  scoerp(2,nH2O,ir) = one
  scoerp(2,nH  ,ir) = one
  !
  frcoef(1,ir) = 2.20e13_wp
  frcoef(2,ir) = zero
  frcoef(3,ir) = 2591.6_wp
  !
  !##############################################
  !
  ! Reaction 4:  O + H2 <=> OH  + H
  !
  ir = 4
  !
  scoerp(1,nO ,ir) = one
  scoerp(1,nH2,ir) = one
  !
  scoerp(2,nOH,ir) = one
  scoerp(2,nH ,ir) = one
  !
  frcoef(1,ir) = 5.06e4_wp
  frcoef(2,ir) = 2.67_wp
  frcoef(3,ir) = 3165.3_wp
  !
  !##############################################
  !
  ! Reaction 5:  2OH <=> H2O + O
  !
  ir = 5
  !
  scoerp(1,nOH ,ir) = two
  !
  scoerp(2,nH2O,ir) = one
  scoerp(2,nO  ,ir) = one
  !
  frcoef(1,ir) = 6.30e12_wp
  frcoef(2,ir) = zero
  frcoef(3,ir) = 548.5_wp
  !
  !##############################################
  !
  ! Reaction 6:  H + OH + M <=> H2O + M
  !
  ir = 6
  !
  scoerp(1,nH   ,ir) = one
  scoerp(1,nOH  ,ir) = one
  scoerp(1,ncs+1,ir) = one
  !
  scoerp(2,nH2O ,ir) = one
  scoerp(2,ncs+1,ir) = one
  !
  cbeffc(nH2 ,ir) = half*five
  cbeffc(nH2O,ir) = sixteen
  !
  frcoef(1,ir) = 2.21e22_wp
  frcoef(2,ir) = -two
  frcoef(3,ir) = zero
  !
  !##############################################
  !
  ! Reaction 7:  H + H + M <=> H2 + M
  !
  ir = 7
  !
  scoerp(1,nH   ,ir) = two
  scoerp(1,ncs+1,ir) = one
  !
  scoerp(2,nH2  ,ir) = one
  scoerp(2,ncs+1,ir) = one
  !
  cbeffc(nH2 ,ir) = half*five
  cbeffc(nH2O,ir) = sixteen
  !
  frcoef(1,ir) = 7.30e17_wp
  frcoef(2,ir) = -one
  frcoef(3,ir) = zero
  !
  !##############################################
  !
  do ir = 1,ncr
    frcoef(1,ir) = log( frcoef(1,ir) )
  end do
  !
  !##############################################
  !
  ! Identify any species that are non-reacting. All species are initially
  ! marked as non-reacting and are marked reacting if either the forward
  ! or backward stoichiometric coefficents for the species are non-zero
  ! for any of the chemical reactions.
  !
  do m = 1,ncs
    nonreacting_species(m) = all(abs(scoerp(1:2,m,1:ncr)) < eps6)
  end do
  !
  ! The forwards and backwards exponents are
  ! the same as the stoichiometric coefficients.
  !
  fwd_exp(1:ncs,1:ncr) = scoerp(1,1:ncs,1:ncr)
  bwd_exp(1:ncs,1:ncr) = scoerp(2,1:ncs,1:ncr)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine initialize_H2_7x7_reactions
!
!###############################################################################
!
subroutine initialize_Ethylene_7x3_reactions
  !
  !.. Use Statements ..
  use thermo_mod, only : species_name
  !
  !.. Local Scalars ..
  integer :: m,ir,ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_Ethylene_7x3_reactions"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Allocate all the arrays needed for H2 reactions
  ! NOTE: cbeffc is initialized to 1 instead of 0
  !
  if (.not. allocated(frcoef)) then
    allocate ( frcoef(1:3,1:ncr) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"frcoef",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (.not. allocated(scoerp)) then
    allocate ( scoerp(1:2,1:ncs+1,1:ncr) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"scoerp",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (.not. allocated(fwd_exp)) then
    allocate ( fwd_exp(1:ncs,1:ncr) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"fwd_exp",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (.not. allocated(bwd_exp)) then
    allocate ( bwd_exp(1:ncs,1:ncr) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"bwd_exp",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (.not. allocated(cbeffc)) then
    allocate ( cbeffc(1:ncs,1:ncr) , source=one , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cbeffc",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  !##############################################
  !
  !     Species:  O2  H2  H2O  C2H4  N2 CO CO2
  !
  !##############################################
  !
  ! Reaction 1:  C2H4 + O2 <=> 2CO + 2H2
  !
  ir = 1
  !
  scoerp(1,nC2H4,ir) = one
  scoerp(1,nO2,ir)   = one
  !
  scoerp(2,nCO,ir)   = two
  scoerp(2,nH2O,ir)  = two
  !
  frcoef(1,ir) = 1.04e13_wp
  frcoef(2,ir) = zero
  frcoef(3,ir) = 21511.0_wp
  !
  !##############################################
  !
  ! Reaction 2:  2CO + O2 <=> 2CO2
  !
  ir = 2
  !
  scoerp(1,nCO,ir)  = two
  scoerp(1,nO2,ir)  = one
  !
  scoerp(2,nCO2,ir) = two
  !
  frcoef(1,ir) = 7.21e7_wp
  frcoef(2,ir) = 1.92_wp
  frcoef(3,ir) = 16478.0_wp
  !
  !##############################################
  !
  ! Reaction 3:  OH + H2 <=> H2O + H
  ! Reaction 3:  2H2 + O2 <=> 2H2O
  !
  ir = 3
  !
  scoerp(1,nH2 ,ir) = two
  scoerp(1,nO2 ,ir) = one
  !
  scoerp(2,nH2O,ir) = two
  !
  frcoef(1,ir) =  9.47e16_wp
  frcoef(2,ir) = -0.9_wp
  frcoef(3,ir) =  1.00_wp
  !
  !##############################################
  !
  do ir = 1,ncr
    frcoef(1,ir) = log( frcoef(1,ir) )
  end do
  !
  !##############################################
  !
  ! Identify any species that are non-reacting. All species are initially
  ! marked as non-reacting and are marked reacting if either the forward
  ! or backward stoichiometric coefficents for the species are non-zero
  ! for any of the chemical reactions.
  !
  do m = 1,ncs
    nonreacting_species(m) = all(abs(scoerp(1:2,m,1:ncr)) < eps6)
  end do
  !
  ! The forwards and backwards exponents are
  ! the same as the stoichiometric coefficients.
  !
  fwd_exp(1:ncs,1:ncr) = scoerp(1,1:ncs,1:ncr)
  bwd_exp(1:ncs,1:ncr) = scoerp(2,1:ncs,1:ncr)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine initialize_Ethylene_7x3_reactions
!
!###############################################################################
!
subroutine compute_wdot(iter)
  !
  !.. Use Statements ..
  use ovar,        only : itestcase,itcur
  use ovar,        only : semi_implicit_chemistry
  use ovar,        only : ignition_start_timestep
  use ovar,        only : ignition_duration
  use ovar,        only : ignition_source_magnitude
  use geomvar_mod, only : n_solpts
  use flowvar_mod, only : usp,gpsp,residual
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iter
  !
  !.. Local Scalars ..
  integer :: n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_wdot"
  !
 !logical(lk), parameter :: use_hardcoded_H2_wdot = true
  logical(lk), parameter :: use_hardcoded_H2_wdot = fals
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (ncr > 0) then
    !
    n1 = 1
    n2 = n_solpts
    !
    if (abs(itestcase) == CH4_3D_Diffusion_Flame .or. &
        abs(itestcase) == CH4_Flame_Speed_Test) then
      if (itestcase > 0) then
        call compute_CH4_hardcoded_wdot(n1,n2,usp,gpsp,residual,iter)
      else
        call compute_H2_wdot(n1,n2,usp,gpsp,residual,iter)
      end if
    else if (itestcase == H2_2D_Diffusion_Flame .or. &
             itestcase == Compression_Corner) then
      if (semi_implicit_chemistry == -1) then
        call compute_H2_hardcoded_wdot(n1,n2,usp,gpsp,residual)
      else
        call compute_H2_wdot(n1,n2,usp,gpsp,residual,iter)
      end if
    else if (itestcase == UIUC_Scramjet) then
        call compute_H2_wdot(n1,n2,usp,gpsp,residual,iter)
    else
      call compute_generic_wdot(n1,n2,usp,gpsp,residual)
    end if
    !
    ! Compute the ignition source term for the energy equation
    !
    if (ignition_duration > 0) then
      if (itcur >= ignition_start_timestep) then
        if (ignition_source_magnitude > zero) then
          call compute_ignition_source_ncc(n1,n2,usp,gpsp,residual)
        else
          call compute_ignition_source(n1,n2,usp,gpsp,residual)
        end if
      end if
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_wdot
!
!###############################################################################
!
subroutine compute_ignition_source_ncc(n1,n2,cv,gp,residual)
  !
  !======================================================================
  !
  !   Calculate the species source term (omega) for each species
  !   For chemical integration, KIVA scheme is used.
  !
  !         Omega = (New_con - Old_con)/time-step
  !
  !      where New_con is the updated mass fraction.
  !      where Old_con is the original mass fraction.
  !======================================================================
  !
  use ovar,        only : uverr,cv_ref,itcur,rhoref,aref,tref
  use ovar,        only : num_rk_stages
  use ovar,        only : ignition_prod2reac_L2_ratio_goal
  use ovar,        only : ignition_max_temperature
  use ovar,        only : ignition_spark_cap
  use ovar,        only : ignition_source_magnitude
  use ovar,        only : ignition_source_radius
  use ovar,        only : ignition_source_location
  use ovar,        only : ignition_duration
  use geomvar_mod, only : xyz
  use flowvar_mod, only : dtsp
  use metrics_mod, only : jac
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv(:,:)
  real(wp),    intent(in) :: gp(:,:)
  real(wp), intent(inout) :: residual(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: r,x,y,z,spark
  real(wp) :: p2r,energy_source
  real(wp) :: products_L2_sum
  real(wp) :: reactants_L2_sum
  !
  !.. Local Saved Scalars ..
  logical(lk), save :: ignition_started = fals
  logical(lk), save :: ignition_reached_max_temp = fals
  logical(lk), save :: ignition_output_max_temp = fals
  integer,     save :: ignition_iter_max_temp = 0
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_ignition_source_ncc"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Perform some initializations if this is
  ! the first time calling this procedure.
  !
  if (.not. ignition_started) then
    !
    ! If ignition_source_radius is zero, change
    ! ignition_prod2reac_L2_ratio_goal to a negative
    ! value to turn disable ignition.
    !
    if (ignition_source_radius == zero) then
      ignition_prod2reac_L2_ratio_goal = -one
      ignition_duration = 0
      if (mypnum == glb_root) then
        write (iout,3)
      end if
      call debug_timer(leaving_procedure,pname)
      return
    end if
    !
    ! Initialize the ignition input values
    !
    ignition_started = true
    ignition_reached_max_temp = fals
    !
    ignition_source_magnitude = ignition_source_magnitude / tref / &
                                real(num_rk_stages,kind=wp)
    ignition_max_temperature = ignition_max_temperature / tref
    ignition_spark_cap = ignition_spark_cap / (rhoref*aref*aref)
    !
    if (mypnum == glb_root) then
      write (iout,1) itcur
    end if
    !
  end if
 !!
 !! Check if the ratio of [the sum of the product species L2 norms] to
 !! [the sum of the reactant species L2 norms] has met the input goal.
 !!
 !if (allocated(uverr)) then
 !  !
 !  ! Compute the ratio of [the sum of the product species L2 norms] to
 !  ! [the sum of the reactant species L2 norms]
 !  !
 !  products_L2_sum  = sum( uverr( product_species_list(:),2,2) )
 !  reactants_L2_sum = sum( uverr(reactant_species_list(:),2,2) )
 !  !
 !  if (reactants_L2_sum /= zero) then
 !    p2r = products_L2_sum / reactants_L2_sum
 !  else
 !    p2r = zero
 !  end if
 ! !p2r = sum( uverr( product_species_list(:),2,2) ) / &
 ! !      sum( uverr(reactant_species_list(:),2,2) )
 !  !
 !  ! If the above ratio is greater than the value of the input goal, we
 !  ! will assume that ignition was successful and reactions are occurring.
 !  ! Set the value of the input goal to a negative number so that ignition
 !  ! source terms will no longer be computed.
 !  !
 !  if (p2r > ignition_prod2reac_L2_ratio_goal) then
 !    ignition_prod2reac_L2_ratio_goal = -one
 !    if (mypnum == glb_root) then
 !      write (iout,2) itcur
 !    end if
 !    call debug_timer(leaving_procedure,pname)
 !    return
 !  end if
 !  !
 !end if
  !
  x = zero
  y = zero
  z = zero
  !
  do n = n1,n2
    !
    x = xyz(1,n) - ignition_source_location(1)
    if (size(xyz,dim=1) >= 2) y = xyz(2,n) - ignition_source_location(2)
    if (size(xyz,dim=1) == 3) z = xyz(3,n) - ignition_source_location(3)
    !
    r = sqrt( x*x + y*y + z*z )
    !
    energy_source = zero
    !
    if (gp(n,ntemp) >= ignition_max_temperature) then
      !
      ! Set the energy source term to zero if the current temperature
      ! is greater than or equal to the ignition max temperature. Also,
      ! mark that the ignition max temperature has been reached to
      ! begin the countdown of the ignition duration.
      !
      ignition_reached_max_temp = true
      !
      energy_source = zero
      !
    else
      !
      ! Compute the value of the energy source term, with the method used
      ! depending on the value of the input ignition_source_radius.
      !
      if (ignition_source_radius < zero) then
        !
        ! If ignition_source_radius is less than zero, use a Gaussian function
        ! based on radial distance from the input ignition location as the
        ! magnitude of the energy source term. The input value of
        ! ignition_source_radius should be approximately the radial distance
        ! from the ignition location at which the Gaussian energy source term
        ! reaches machine precision.
        !
        if (r <= abs(ignition_source_radius)) then
          !
          r = r/abs(ignition_source_radius)
          !
          spark = cv(n,nec) * gp(n,ncp) * ignition_source_magnitude
          !
          energy_source = min(spark,ignition_spark_cap)
          !
          energy_source = energy_source*exp(-36.0_wp*r*r)
          !
        end if
        !
        ! If the ignition_source_radius is greater than zero, every point that
        ! is within the ignition_source_radius from the input ignition location
        ! will get the full value of the energy source term based on the input
        ! value of ignition_source_magnitude. The energy source term for all
        ! points outside the ignition_source_radius from the input ignition
        ! location will be an exponential decay based on radius.
        !
      else if (ignition_source_radius > zero) then
        !
        spark = cv(n,nec) * gp(n,ncp) * ignition_source_magnitude
        !
        energy_source = min(spark,ignition_spark_cap)
        !
        if (r >= ignition_source_radius) then
          !
          r = (r-ignition_source_radius) / ignition_source_radius
          !
          energy_source = energy_source*exp(-36.0_wp*r*r*r)
          !
        end if
        !
      end if
      !
    end if
    !
    ! Add the ignition source term to the RHS of the energy equation
    ! NOTE: We actually need to subtract it from the RHS because the
    !       residual array will be multiplied by -1 at the end of the
    !       current time stage/step before using it to update the
    !       solution variables.
    !
    residual(n,nee) = residual(n,nee) - energy_source*jac(n)/dtsp(n)
    !
  end do
  !
  if (ncpu > 1) then
    n = merge(1,0,ignition_reached_max_temp)
    call mpi_allreduce(MPI_IN_PLACE,n,1_int_mpi,mpi_inttyp, &
                       MPI_MAX,MPI_COMM_WORLD,mpierr)
    if (n == 1) then
      ignition_reached_max_temp = true
    end if
  end if
  !
  if (ignition_reached_max_temp) then
    if (.not. ignition_output_max_temp) then
      ignition_iter_max_temp = itcur
      ignition_output_max_temp = true
      if (mypnum == glb_root) then
        write (iout,4) itcur
      end if
    end if
    if (itcur == ignition_iter_max_temp+ignition_duration) then
      ignition_duration = 0
    end if
    if (ignition_duration == 0) then
      if (mypnum == glb_root) then
        write (iout,2) itcur
      end if
    end if
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  1 format (//,2x,80("#"),//, &
               6x,"STARTING IGNITION SOURCE TERM AT TIMESTEP # ",i0,//, &
               2x,80("#"),//)
  2 format (//,2x,80("#"),//, &
               6x,"TURNING OFF IGNITION SOURCE TERM AT TIMESTEP # ",i0,//, &
               2x,80("#"),//)
  3 format (//,2x,80("#"),//, &
               6x,"DISABLING THE IGNITION SOURCE TERM BECAUSE",/, &
               6x,"THE IGNITION SOURCE TERM RADIUS IS ZERO!!",//, &
               2x,80("#"),//)
  4 format (//,2x,80("#"),//, &
               6x,"REACHED IGNITION MAX TEMPERATURE AT TIMESTEP # ",i0,//, &
               2x,80("#"),//)
  !
end subroutine compute_ignition_source_ncc
!
!###############################################################################
!
subroutine compute_ignition_source(n1,n2,cv,gp,residual)
  !
  !======================================================================
  !
  !   Calculate the species source term (omega) for each species
  !   For chemical integration, KIVA scheme is used.
  !
  !         Omega = (New_con - Old_con)/time-step
  !
  !      where New_con is the updated mass fraction.
  !      where Old_con is the original mass fraction.
  !======================================================================
  !
  use ovar,        only : uverr,cv_ref,itcur,rhoref,aref,tref
  use ovar,        only : ignition_prod2reac_L2_ratio_goal
  use ovar,        only : ignition_max_temperature
  use ovar,        only : ignition_spark_cap
  use ovar,        only : ignition_source_magnitude
  use ovar,        only : ignition_source_radius
  use ovar,        only : ignition_source_location
  use geomvar_mod, only : xyz
  use flowvar_mod, only : dtsp
  use metrics_mod, only : jac
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv(:,:)
  real(wp),    intent(in) :: gp(:,:)
  real(wp), intent(inout) :: residual(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: r,x,y,z,spark
  real(wp) :: p2r,energy_source
  real(wp) :: products_L2_sum
  real(wp) :: reactants_L2_sum
  !
  !.. Local Saved Scalars ..
  real(wp),    save :: ign_mag = zero
  logical(lk), save :: ignition_started = fals
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_ignition_source"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Perform some initializations if this is
  ! the first time calling this procedure.
  !
  if (.not. ignition_started) then
    !
    ! If ignition_source_radius is zero, change
    ! ignition_prod2reac_L2_ratio_goal to a negative
    ! value to turn disable ignition.
    !
    if (ignition_source_radius == zero) then
      ignition_prod2reac_L2_ratio_goal = -one
      if (mypnum == glb_root) then
        write (iout,3)
      end if
      call debug_timer(leaving_procedure,pname)
      return
    end if
    !
    ! Initialize the ignition input values
    !
    ignition_started = true
    !
    ignition_max_temperature = ignition_max_temperature / tref
    ignition_spark_cap = ignition_spark_cap / (rhoref*aref*aref)
    !
    if (mypnum == glb_root) then
      write (iout,1) itcur
    end if
    !
  end if
  !
  ! Check if the ratio of [the sum of the product species L2 norms] to
  ! [the sum of the reactant species L2 norms] has met the input goal.
  !
  !
  if (.not. ignition_started) then
    ignition_started = true
    if (mypnum == glb_root) then
      write (iout,1) itcur
    end if
  end if
  !
  if (allocated(uverr)) then
    !
    ! Compute the ratio of [the sum of the product species L2 norms] to
    ! [the sum of the reactant species L2 norms]
    !
    products_L2_sum  = sum( uverr( product_species_list(:),2,2) )
    reactants_L2_sum = sum( uverr(reactant_species_list(:),2,2) )
    !
    if (reactants_L2_sum /= zero) then
      p2r = products_L2_sum / reactants_L2_sum
    else
      p2r = zero
    end if
   !p2r = sum( uverr( product_species_list(:),2,2) ) / &
   !      sum( uverr(reactant_species_list(:),2,2) )
    !
    ! If the above ratio is greater than the value of the input goal, we
    ! will assume that ignition was successful and reactions are occurring.
    ! Set the value of the input goal to a negative number so that ignition
    ! source terms will no longer be computed.
    !
    if (p2r > ignition_prod2reac_L2_ratio_goal) then
      ignition_prod2reac_L2_ratio_goal = -one
      if (mypnum == glb_root) then
        write (iout,2) itcur
      end if
      call debug_timer(leaving_procedure,pname)
      return
    end if
    !
  end if
  !
  ign_mag = ign_mag + abs(ignition_source_magnitude) * cv_ref(nee)
  !
  spark = min(ign_mag,ignition_spark_cap)
  !
  x = zero
  y = zero
  z = zero
  !
  do n = n1,n2
    !
    x = xyz(1,n) - ignition_source_location(1)
    if (size(xyz,dim=1) >= 2) y = xyz(2,n) - ignition_source_location(2)
    if (size(xyz,dim=1) == 3) z = xyz(3,n) - ignition_source_location(3)
    !
    r = sqrt( x*x + y*y + z*z )
    !
    energy_source = zero
    !
    if (ignition_source_radius < zero) then
      !
      ! If ignition_source_radius is less than zero, use a Gaussian function
      ! based on radial distance from the input ignition location as the
      ! magnitude of the energy source term. The input value of
      ! ignition_source_radius should be approximately the radial distance from
      ! the ignition location at which the Gaussian energy source term reaches
      ! machine precision.
      !
      if (r <= abs(ignition_source_radius)) then
        !
        r = r/abs(ignition_source_radius)
        !
        energy_source = spark*exp(-36.0_wp*r*r)
        !
        ! If the temperature at this point is greater than the input maximum
        ! ignition temperature, change ignition_prod2reac_L2_ratio_goal to
        ! a negative value to turn off ignition.
        !
        if (gp(n,ntemp) >= ignition_max_temperature) then
          ignition_prod2reac_L2_ratio_goal = -one
        end if
        !
      end if
      !
      ! If the ignition_source_radius is greater than zero, every point that is
      ! within the ignition_source_radius from the input ignition location will
      ! get the full value of the energy source term based on the input value of
      ! ignition_source_magnitude. The energy source term for all points outside
      ! the ignition_source_radius from the input ignition location will be an
      ! exponential decay based on radius.
      !
    else if (ignition_source_radius > zero) then
      !
      energy_source = spark
      !
      if (r >= ignition_source_radius) then
        !
        r = (r-ignition_source_radius)/ignition_source_radius
        !
        energy_source = energy_source*exp(-36.0_wp*r*r*r)
        !
      end if
      !
      ! If the temperature at this point is greater than the input maximum
      ! ignition temperature, change ignition_prod2reac_L2_ratio_goal to
      ! a negative value to turn off ignition.
      !
      if (gp(n,ntemp) >= ignition_max_temperature) then
        ignition_prod2reac_L2_ratio_goal = -one
      end if
      !
    end if
    !
    ! Add the ignition source term to the RHS of the energy equation
    ! NOTE: We actually need to subtract it from the RHS because the
    !       residual array will be multiplied by -1 at the end of the
    !       current time stage/step before using it to update the
    !       solution variables.
    !
    residual(n,nee) = residual(n,nee) - energy_source*jac(n)/dtsp(n)
    !
  end do
  !
  if (ncpu > 1) then
    call mpi_allreduce(MPI_IN_PLACE,ignition_prod2reac_L2_ratio_goal, &
                       1_int_mpi,mpi_flttyp,MPI_MIN,MPI_COMM_WORLD,mpierr)
  end if
  !
  if (ignition_prod2reac_L2_ratio_goal < zero) then
    if (mypnum == glb_root) then
      write (iout,2) itcur
    end if
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  1 format (//,2x,80("#"),//, &
               6x,"STARTING IGNITION SOURCE TERM AT TIMESTEP # ",i0,//, &
               2x,80("#"),//)
  2 format (//,2x,80("#"),//, &
               6x,"TURNING OFF IGNITION SOURCE TERM AT TIMESTEP # ",i0,//, &
               2x,80("#"),//)
  3 format (//,2x,80("#"),//, &
               6x,"DISABLING THE IGNITION SOURCE TERM BECAUSE",/, &
               6x,"THE IGNITION SOURCE TERM RADIUS IS ZERO!!",//, &
               2x,80("#"),//)
  !
end subroutine compute_ignition_source
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_CH4_hardcoded_wdot(n1,n2,cv,gp,residual,iter)
  !
  !======================================================================
  !
  !   Calculate the species source term (omega) for each species
  !   For chemical integration, KIVA scheme is used.
  !
  !         Omega = (New_con - Old_con)/time-step
  !
  !      where New_con is the updated mass fraction.
  !      where Old_con is the original mass fraction.
  !======================================================================
  !
  use ovar,        only : rhoref,aref
  use flowvar_mod, only : dtsp
  use metrics_mod, only : jac
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: iter
  real(wp),    intent(in) :: cv(:,:)
  real(wp),    intent(in) :: gp(:,:)
  real(wp), intent(inout) :: residual(:,:)
  !
  !.. Local Scalars ..
  integer  :: k,m,n,ir
  real(wp) :: dt_inv,dt_chem
  real(wp) :: temp,arrh,tact,qf,q1,q2,qref,domega
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: rhoi(:,:)
  real(wp), allocatable :: rho_ratio(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_CH4_hardcoded_wdot"
  !
  logical(lk), parameter :: nondimensional_wdot = true
 !logical(lk), parameter :: nondimensional_wdot = fals
  !
continue
  !
  call pure_debug_timer(entering_procedure,pname)
  !
  !
  allocate ( rhoi(n1:n2,1:ncs) )
  allocate ( rho_ratio(n1:n2) )
  !
  ! Initialize the species densities array
  !
  do m = 1,ncs
    do n = n1,n2
      rhoi(n,m) = cv(n,m)
    end do
  end do
  !
  chem_sub_iter: &
  do k = 1,n_reac_sub_iter
    !
    !======================================================================
    !
    !  Use the method defined in Appendix I, page 128 of Loa Alomas Report
    !  LA-11560-MS, "KIVA-II: A Computer Program for Chemically Reacting
    !  Flows with Sprays", Amsden, O'Rourke, Bultler, May 1996,
    !  and defined in the source code for subroutine chem() of kiva3re.f
    !
    !  Input:
    !       chem_model
    !       chem_turb_model,  ! added 07/11/2016
    !       number_of_reaction
    !       ncs
    !       dt_chem       - time step for chemistry integration {sec}
    !       mass_frac     - array of mass fractions {non-dimensional}
    !                       before integration over dt_chem
    !
    !  Output:
    !       mass_frac     - array of mass fractions {non-dimensional}
    !                       after integration over dt_chem
    !
    !  Implemented into NCC by: Jeffrey Moder
    !                           Aug 2, 2004
    !
    !  Updated to prevent overflows and to account for possibility
    !  of negative mass densities spd(isp) < 0 by:
    !                           Jeffrey Moder
    !                           Sept 29, 2004
    !
    !======================================================================
    !
    ! Use KIVA II method of updating species mass density using contributions
    ! from one reaction at a time, with each reaction progress rate defined
    ! in terms of what KIVA calls a "reference species".
    !
    do n = n1,n2
      !
      dt_chem = dtsp(n) / real(n_reac_sub_iter,kind=wp)
      !
      temp = gp(n,ntemp)
      !
      !=========================================================================
      !=========================================================================
      !=========================================================================
      !
      !    REACTION 1: CH4 + 1.5 O2 => CO + 2 H2O
      !
      !-------------------------------------------------------------------------
      !
      ir = 1
      !
      ! Arrhenius coefficient pre-multiplied by bunch of stuff
      arrh = fwd_arrh_coef(ir)
      ! Activation temperature, nondimensionalized by tref
      tact = fwd_actv_temp(ir)
      !
      ! Compute the forward kinetic reaction rate
      !
      qf = arrh*exp(-tact/temp) * (rhoi(n,nCH4)**fwd_exp(nCH4,ir)) &
                                * (rhoi(n,nO2 )**fwd_exp(nO2 ,ir))
      !
      ! If ( rhoi(nCH4) = 0 ) and/or ( rhoi(nO2) = 0 ), then qf = 0,
      ! so q1 = q2 = qref = 0, and then domega = 0.
      !
      q1 = qf*fwd_stoich_coef(nCH4,ir)/max(rhoi(n,nCH4),rndoff)
      q2 = qf*fwd_stoich_coef(nO2 ,ir)/max(rhoi(n,nO2 ),rndoff)
      !
      qref = max( q1 , q2 )
      !
      domega = dt_chem * qf / (one + dt_chem*qref)
      !
      rhoi(n,nCO ) = max( zero , rhoi(n,nCO ) + stoich_coef(nCO ,ir)*domega )
      rhoi(n,nHP ) = max( zero , rhoi(n,nHP ) + stoich_coef(nHP ,ir)*domega )
      rhoi(n,nCH4) = max( zero , rhoi(n,nCH4) + stoich_coef(nCH4,ir)*domega )
      rhoi(n,nO2 ) = max( zero , rhoi(n,nO2 ) + stoich_coef(nO2 ,ir)*domega )
      !
      !=========================================================================
      !=========================================================================
      !=========================================================================
      !
      !    REACTION 2: CO2 => CO + 0.5 O2
      !
      !-------------------------------------------------------------------------
      !
      ir = 2
      !
      ! Arrhenius coefficient pre-multiplied by bunch of stuff
      arrh = fwd_arrh_coef(ir)
      ! Activation temperature, nondimensionalized by tref
      tact = fwd_actv_temp(ir)
      !
      ! Compute the forward kinetic reaction rate
      !
      qf = arrh*exp(-tact/temp) * (rhoi(n,nCO2)**fwd_exp(nCO2,ir))
      !
      ! If rhoi(nCO2) = 0, qf = 0, so qref = 0, and then domega = 0
      !
      qref = qf*fwd_stoich_coef(nCO2,ir)/max(rhoi(n,nCO2),rndoff)
      !
      domega = dt_chem * qf / (one + dt_chem*qref)
      !
      rhoi(n,nCO2) = max( zero , rhoi(n,nCO2) + stoich_coef(nCO2,ir)*domega )
      rhoi(n,nCO ) = max( zero , rhoi(n,nCO ) + stoich_coef(nCO ,ir)*domega )
      rhoi(n,nO2 ) = max( zero , rhoi(n,nO2 ) + stoich_coef(nO2 ,ir)*domega )
      !
      !=========================================================================
      !=========================================================================
      !=========================================================================
      !
      !    REACTION 3: CO + 0.5 O2 => CO2
      !
      !-------------------------------------------------------------------------
      !
      ir = 3
      !
      ! Arrhenius coefficient pre-multiplied by bunch of stuff
      arrh = fwd_arrh_coef(ir)
      ! Activation temperature, nondimensionalized by tref
      tact = fwd_actv_temp(ir)
      !
      ! Compute the forward kinetic reaction rate
      !
      qf = arrh*exp(-tact/temp) * (rhoi(n,nCO)**fwd_exp(nCO,ir)) &
                                * (rhoi(n,nO2)**fwd_exp(nO2,ir))
      !
      ! If ( rhoi(nCO) = 0 ) and/or ( rhoi(nO2) = 0 ), then qf = 0,
      ! so q1 = q2 = qref = 0, and then domega = 0.
      !
      q1 = qf*fwd_stoich_coef(nCO,ir)/max(rhoi(n,nCO),rndoff)
      q2 = qf*fwd_stoich_coef(nO2,ir)/max(rhoi(n,nO2),rndoff)
      !
      qref = max( q1 , q2 )
      !
      domega = dt_chem * qf / (one + dt_chem*qref)
      !
      rhoi(n,nCO2) = max( zero , rhoi(n,nCO2) + stoich_coef(nCO2,ir)*domega )
      rhoi(n,nCO ) = max( zero , rhoi(n,nCO ) + stoich_coef(nCO ,ir)*domega )
      rhoi(n,nO2 ) = max( zero , rhoi(n,nO2 ) + stoich_coef(nO2 ,ir)*domega )
      !
      !=========================================================================
      !=========================================================================
      !=========================================================================
      !
    end do
    !
    ! Compute the new mixture density by summing the species densities
    !
    do n = n1,n2
      rho_ratio(n) = rhoi(n,1)
    end do
    !
    do m = 2,ncs
      do n = n1,n2
        rho_ratio(n) = rho_ratio(n) + rhoi(n,m)
      end do
    end do
    !
    ! Compute the ratio of old mixture density to new mixture density
    !
    do n = n1,n2
      rho_ratio(n) = cv(n,nec)/rho_ratio(n)
    end do
    !
    ! Keeping the mixture density the same as upon entering this routine,
    !
    ! Update the species densities by multiplying the new species densities
    ! by the ratio of the old mixture density to the new mixture density.
    ! NOTE: This basically holds the mixture density constant to the old
    !       mixture density value, and updates the species densities based
    !       on the changes to the species mass fractions based on the
    !       chemical reaction source terms.
    !
    do m = 1,ncs
      do n = n1,n2
        rhoi(n,m) = rhoi(n,m) * rho_ratio(n)
      end do
    end do
    !
  end do chem_sub_iter
  !
  if (allocated(rho_ratio)) deallocate ( rho_ratio )
  !
  ! Add the change in species densities to the RHS residual
  !
  do m = 1,ncs
    !
    ! Leave the residual unchanged if this is a nonreacting species
    !
    if (nonreacting_species(m)) cycle
    !
    do n = n1,n2
      residual(n,m) = residual(n,m) - jac(n)*(rhoi(n,m) - cv(n,m))/dtsp(n)
    end do
    !
  end do
  !
  if (allocated(rhoi)) deallocate ( rhoi )
  !
  call pure_debug_timer(leaving_procedure,pname)
  !
end subroutine compute_CH4_hardcoded_wdot
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_H2_hardcoded_wdot(n1,n2,cv,gp,residual)
  !
  ! Routine to compute the chemical source terms
  !
  ! Subroutine COMPUTE_H2_HARDCODED_WDOT computes:
  !  2) The species concentrations
  !  2) The forward and backward reaction rate coefficients
  !  3) The species production rate
  !
  ! N1     : beginning point
  ! N2     : ending    point
  ! CV     : conserved variables
  ! T      : static temperature
  ! SPECON : species concentrations
  ! FRRATE : forward reaction rate coefficient
  ! BRRATE : backward reaction rate coefficient
  ! WDOT   : species source terms
  !
  !.. Use Statements ..
  use ovar,        only : rhoref,aref,tref
  use metrics_mod, only : jac
  use thermo_mod,  only : runivl,wtmol
  use thermo_mod,  only : compute_gibbs_free_energy
  !
  ! .. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv(:,:)
  real(wp),    intent(in) :: gp(:,:)
  real(wp), intent(inout) :: residual(:,:)
  !
  ! .. Local Scalars ..
  integer  :: ir,m,n,np
  real(wp) :: val,log_a,expnt,tpwr,rulog,tactiv
  real(wp) :: temp,logt,tbe
  real(wp) :: qf,qb,dgibrt,prdctf,prdctb
  real(wp) :: scO2,scH2,scOH,scH2O,scO,scH,scN2
  real(wp) :: wdotO2,wdotH2,wdotOH,wdotH2O,wdotO,wdotH
  real(wp) :: prd1,prd2,prd3,prd4,prd5,prd6,prd7
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: t(:)
  real(wp), allocatable :: tlog(:)
  real(wp), allocatable :: gibrt(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_H2_hardcoded_wdot"
  !
continue
  !
  call pure_debug_timer(entering_procedure,pname)
  !
  np = size(cv,dim=1)
  !
  rulog = log(runivl)
  !
  val = dens_cgs_to_si  / (rhoref*aref)
  !
  ! Compute the Gibbs free energy / (Runiv*T)
  !
  allocate ( gibrt(1:np,1:ncs) , source=zero )
  !
  call compute_gibbs_free_energy(n1,n2,gp(:,ntemp),gibrt)
  !
  ! Calculate the forwards and backward rates
  !
  do n = n1,n2
    !
    ! Compute the forward reaction rate coefficient
    !
    temp = gp(n,ntemp)*tref
    logt = log(temp)
    !
    ! Compute the dimensional species concentrations in moles/cm^3
    !
    scO2  = cv(n,nO2 ) * dens_si_to_cgs * rhoref / wtmol(nO2 )
    scH2  = cv(n,nH2 ) * dens_si_to_cgs * rhoref / wtmol(nH2 )
    scOH  = cv(n,nOH ) * dens_si_to_cgs * rhoref / wtmol(nOH )
    scH2O = cv(n,nH2O) * dens_si_to_cgs * rhoref / wtmol(nH2O)
    scO   = cv(n,nO  ) * dens_si_to_cgs * rhoref / wtmol(nO  )
    scH   = cv(n,nH  ) * dens_si_to_cgs * rhoref / wtmol(nH  )
    scN2  = cv(n,nN2 ) * dens_si_to_cgs * rhoref / wtmol(nN2 )
    !
    ! Compute the third body efficiencies
    !
    tbe = scO2 + half*five*scH2 + scOH + sixteen*scH2O + scO + scH + scN2
    !
    !=========================================================================
    !=========================================================================
    !=========================================================================
    !
    !    REACTION 1:  O2 + H2 => 2OH
    !
    !-------------------------------------------------------------------------
    !
    ir = 1
    !
    log_a = frcoef(1,ir)
    tpwr = frcoef(2,ir)
    tactiv = frcoef(3,ir)
    !
    expnt = log_a + tpwr*logt - tactiv/temp
    qf = exp( min(exp_lim,expnt) )
    !
    ! Compute the backwards reaction rate coefficient
    ! using an equilibrium assumption
    !
    dgibrt = (gibrt(n,nO2) + gibrt(n,nH2)) - two*gibrt(n,nOH)
    qb = exp( min(exp_lim,expnt+dgibrt) )
    !
    ! Calculate the species source terms
    !
    prdctf = qf * scO2 * scH2
    prdctb = qb * scOH * scOH
    !
    ! Compute net effect of forward and backward reactions
    !
    prd1 = prdctf - prdctb
    !
    !=========================================================================
    !=========================================================================
    !=========================================================================
    !
    !    REACTION 2:  O2 + H <=> OH + O
    !
    !-------------------------------------------------------------------------
    !
    ir = 2
    !
    log_a = frcoef(1,ir)
    tpwr = frcoef(2,ir)
    tactiv = frcoef(3,ir)
    !
    expnt = log_a + tpwr*logt - tactiv/temp
    qf = exp( min(exp_lim,expnt) )
    !
    ! Compute the backwards reaction rate coefficient
    ! using an equilibrium assumption
    !
    dgibrt = (gibrt(n,nO2) + gibrt(n,nH)) - (gibrt(n,nOH) + gibrt(n,nO))
    qb = exp( min(exp_lim,expnt+dgibrt) )
    !
    ! Calculate the species source terms
    !
    prdctf = qf * scO2 * scH
    prdctb = qb * scOH * scO
    !
    prd2 = prdctf - prdctb
    !
    !=========================================================================
    !=========================================================================
    !=========================================================================
    !
    !    REACTION 3:  H2 + OH <=> H2O + H
    !
    !-------------------------------------------------------------------------
    !
    ir = 3
    !
    log_a = frcoef(1,ir)
    tpwr = frcoef(2,ir)
    tactiv = frcoef(3,ir)
    !
    expnt = log_a + tpwr*logt - tactiv/temp
    qf = exp( min(exp_lim,expnt) )
    !
    ! Compute the backwards reaction rate coefficient
    ! using an equilibrium assumption
    !
    dgibrt = (gibrt(n,nH2) + gibrt(n,nOH)) - (gibrt(n,nH2O) + gibrt(n,nH))
    qb = exp( min(exp_lim,expnt+dgibrt) )
    !
    ! Calculate the species source terms
    !
    prdctf = qf * scH2  * scOH
    prdctb = qb * scH2O * scH
    !
    ! Compute net effect of forward and backward reactions
    !
    prd3 = prdctf - prdctb
    !
    !=========================================================================
    !=========================================================================
    !=========================================================================
    !
    !    REACTION 4:  H2 + O <=> OH  + H
    !
    !-------------------------------------------------------------------------
    !
    ir = 4
    !
    log_a = frcoef(1,ir)
    tpwr = frcoef(2,ir)
    tactiv = frcoef(3,ir)
    !
    expnt = log_a + tpwr*logt - tactiv/temp
    qf = exp( min(exp_lim,expnt) )
    !
    ! Compute the backwards reaction rate coefficient
    ! using an equilibrium assumption
    !
    dgibrt = (gibrt(n,nH2) + gibrt(n,nO)) - (gibrt(n,nOH) + gibrt(n,nH))
    qb = exp( min(exp_lim,expnt+dgibrt) )
    !
    ! Calculate the species source terms
    !
    prdctf = qf * scH2 * scO
    prdctb = qb * scOH * scH
    !
    ! Compute net effect of forward and backward reactions
    !
    prd4 = prdctf - prdctb
    !
    !=========================================================================
    !=========================================================================
    !=========================================================================
    !
    !    REACTION 5:  2OH <=> H2O + O
    !
    !-------------------------------------------------------------------------
    !
    ir = 5
    !
    log_a = frcoef(1,ir)
    tpwr = frcoef(2,ir)
    tactiv = frcoef(3,ir)
    !
    expnt = log_a + tpwr*logt - tactiv/temp
    qf = exp( min(exp_lim,expnt) )
    !
    ! Compute the backwards reaction rate coefficient
    ! using an equilibrium assumption
    !
    dgibrt = two*gibrt(n,nOH) - (gibrt(n,nH2O) + gibrt(n,nO))
    qb = exp( min(exp_lim,expnt+dgibrt) )
    !
    ! Calculate the species source terms
    !
    prdctf = qf * scOH  * scOH
    prdctb = qb * scH2O * scO
    !
    ! Compute net effect of forward and backward reactions
    !
    prd5 = prdctf - prdctb
    !
    !=========================================================================
    !=========================================================================
    !=========================================================================
    !
    !    REACTION 6:  OH + H + M <=> H2O + M
    !
    !-------------------------------------------------------------------------
    !
    ir = 6
    !
    log_a = frcoef(1,ir)
    tpwr = frcoef(2,ir)
    tactiv = frcoef(3,ir)
    !
    expnt = log_a + tpwr*logt - tactiv/temp
    qf = exp( min(exp_lim,expnt) )
    !
    ! Compute the backwards reaction rate coefficient
    ! using an equilibrium assumption
    !
    dgibrt = (gibrt(n,nOH) + gibrt(n,nH)) - gibrt(n,nH2O)
    qb = exp( min(exp_lim,expnt + dgibrt - logt - rulog) )
    !
    ! Calculate the species source terms
    !
    prdctf = qf * scOH * scH
    prdctb = qb * scH2O
    !
    ! Compute net effect of forward and backward reactions
    !
    prd6 = (prdctf - prdctb) * tbe
    !
    !=========================================================================
    !=========================================================================
    !=========================================================================
    !
    !    REACTION 7:  H + H + M <=> H2 + M
    !
    !-------------------------------------------------------------------------
    !
    ir = 7
    !
    log_a = frcoef(1,ir)
    tpwr = frcoef(2,ir)
    tactiv = frcoef(3,ir)
    !
    expnt = log_a + tpwr*logt - tactiv/temp
    qf = exp( min(exp_lim,expnt) )
    !
    ! Compute the backwards reaction rate coefficient
    ! using an equilibrium assumption
    !
    dgibrt = two*gibrt(n,nH) - gibrt(n,nH2)
    qb = exp( min(exp_lim,expnt + dgibrt - logt - rulog) )
    !
    ! Calculate the species source terms
    !
    prdctf = qf * scH * scH
    prdctb = qb * scH2
    !
    ! Compute net effect of forward and backward reactions
    !
    prd7 = (prdctf - prdctb) * tbe
    !
    !=========================================================================
    !=========================================================================
    !=========================================================================
    !
    ! Compute net effect of forward and backward reactions
    !
    wdotO2  =    -prd1 - prd2
    wdotH2  =    -prd1        - prd3 - prd4                   +     prd7
    wdotOH  = two*prd1 + prd2 - prd3 + prd4 - two*prd5 - prd6
    wdotH2O =                   prd3 +            prd5 + prd6
    wdotO   =            prd2 -        prd4 +     prd5
    wdotH   =          - prd2 + prd3 + prd4            - prd6 - two*prd7
    !
    ! Convert species production rate to kg/(m**3-sec), non-dimensionalize,
    ! and add the species source terms to the residual
    !
    residual(n,nO2 ) = residual(n,nO2 ) - jac(n)*wdotO2 *wtmol(nO2 )*val
    residual(n,nH2 ) = residual(n,nH2 ) - jac(n)*wdotH2 *wtmol(nH2 )*val
    residual(n,nOH ) = residual(n,nOH ) - jac(n)*wdotOH *wtmol(nOH )*val
    residual(n,nH2O) = residual(n,nH2O) - jac(n)*wdotH2O*wtmol(nH2O)*val
    residual(n,nO  ) = residual(n,nO  ) - jac(n)*wdotO  *wtmol(nO  )*val
    residual(n,nH  ) = residual(n,nH  ) - jac(n)*wdotH  *wtmol(nH  )*val
    !
  end do
  !
  if (allocated(gibrt)) deallocate ( gibrt )
  !
  call pure_debug_timer(leaving_procedure,pname)
  !
end subroutine compute_H2_hardcoded_wdot
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_H2_wdot(n1,n2,cv,gp,residual,iter)
  !
  ! N1       : beginning point
  ! N2       : ending    point
  ! CV       : conserved variables
  ! GP       : gas properties
  ! RESIDUAL : RHS residual
  !
  !.. Use Statements ..
  use ovar,        only : semi_implicit_chemistry
  use metrics_mod, only : jac
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: iter
  real(wp),    intent(in) :: cv(:,:)
  real(wp), intent(inout) :: gp(:,:)
  real(wp), intent(inout) :: residual(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,np
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: wdot(:,:)
  real(wp), allocatable :: wdotsv(:,:)
  real(wp), allocatable :: frrate(:,:)
  real(wp), allocatable :: brrate(:,:)
  real(wp), allocatable :: specon(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_H2_wdot"
  !
continue
  !
  call pure_debug_timer(entering_procedure,pname)
  !
  np = size(cv,dim=1)
  !
  allocate ( wdot(1:np,1:ncs) , source=zero )
  allocate ( frrate(1:np,1:ncs) , source=zero )
  allocate ( brrate(1:np,1:ncs) , source=zero )
  allocate ( specon(1:np,1:ncs) , source=zero )
  !
  ! Calculate the species source terms
  !
  call compute_species_source_terms(n1,n2,cv,gp(:,ntemp),specon, &
                                    frrate,brrate,wdot)
  !
  if (semi_implicit_chemistry > 0) then
    !
    if (semi_implicit_chemistry == 2) then
      wdotsv = wdot(1:np,1:ncs) ! F2003 AUTO-REALLOCATION
    end if
    !
    ! Update the source residual
    !
    do m = 1,ncs
      do n = n1,n2
        wdot(n,m) = residual(n,m) - jac(n)*wdot(n,m)
      end do
    end do
    !
    ! Implicit source terms for explicit solver
    !
    if (semi_implicit_chemistry == 1) then
      !
      ! Update the residual using a semi-implicit formulation for the
      ! species source terms using an analytical Jacobian for the
      ! species equations.
      !
      call semi_implicit_exact_jac(n1,n2,cv,gp,frrate,brrate,specon,wdot)
      !
    else if (semi_implicit_chemistry == 2) then
      !
      ! Update the residual using a semi-implicit formulation for the
      ! species source terms using a numerical approximation of the
      ! Jacobian for the species equations.
      !
      call semi_implicit_numerical_jac(n1,n2,cv,gp,frrate,brrate, &
                                       specon,wdot,wdotsv)
      !
    end if
    !
    if (allocated(wdotsv)) deallocate ( wdotsv )
    !
    ! Update residual
    !
    do m = 1,ncs
      !
      ! Leave the residual unchanged if this is a nonreacting species
      !
      if (nonreacting_species(m)) cycle
      !
      do n = n1,n2
        residual(n,m) = wdot(n,m)
      end do
      !
    end do
    !
  else
    !
    do m = 1,ncs
      !
      ! Leave the residual unchanged if this is a nonreacting species
      !
      if (nonreacting_species(m)) cycle
      !
      do n = n1,n2
        residual(n,m) = residual(n,m) - jac(n)*wdot(n,m)
      end do
      !
    end do
    !
  end if
  !
  if (allocated(specon)) deallocate ( specon )
  if (allocated(brrate)) deallocate ( brrate )
  if (allocated(frrate)) deallocate ( frrate )
  if (allocated(wdot))  deallocate ( wdot )
  !
  call pure_debug_timer(leaving_procedure,pname)
  !
end subroutine compute_H2_wdot
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_species_source_terms(n1,n2,cv,t,specon,frrate,brrate,wdot)
  !
  ! Routine to compute the chemical source terms
  !
  ! Subroutine COMPUTE_SPECIES_SOURCE_TERMS computes:
  !  2) The species concentrations
  !  2) The forward and backward reaction rate coefficients
  !  3) The species production rate
  !
  ! N1     : beginning point
  ! N2     : ending    point
  ! CV     : conserved variables
  ! T      : static temperature
  ! SPECON : species concentrations
  ! FRRATE : forward reaction rate coefficient
  ! BRRATE : backward reaction rate coefficient
  ! WDOT   : species source terms
  !
  !.. Use Statements ..
  use ovar,       only : rhoref
  use thermo_mod, only : wtmol
  !
  ! .. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv(:,:)
  real(wp),    intent(in) :: t(:)
  real(wp), intent(inout) :: frrate(:,:)
  real(wp), intent(inout) :: brrate(:,:)
  real(wp), intent(inout) :: wdot(:,:)
  real(wp), intent(inout) :: specon(:,:)
  !
  ! .. Local Scalars ..
  integer  :: n,m
  real(wp) :: val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_species_source_terms"
  !
continue
  !
  call pure_debug_timer(entering_procedure,pname)
  !
  ! Compute the species concentrations in (moles/cm^3)
  !
  do m = 1,ncs
    val = dens_si_to_cgs * rhoref / wtmol(m)
    do n = n1,n2
      specon(n,m) = cv(n,m) * val
    end do
  end do
  !
  ! Calculate the forwards and backward rates
  !
  call compute_reaction_rates(n1,n2,t,frrate,brrate)
  !
  ! Calculate the species source terms
  !
  call compute_species_production_rates(n1,n2,specon,frrate,brrate,wdot)
  !
  call pure_debug_timer(leaving_procedure,pname)
  !
end subroutine compute_species_source_terms
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_reaction_rates(n1,n2,t,frrate,brrate)
  !
  ! Routine to calculate the forward and backward reaction rate
  ! coefficients assuming no effect of temperature fluctuations
  !
  ! N1     : beginning point
  ! N2     : ending    point
  ! T      : static temperature
  ! FRRATE : forward reaction rate coefficient
  ! BRRATE : backward reaction rate coefficient
  !
  !.. Use Statements ..
  use ovar, only : itestcase
  !
  ! .. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: t(:)
  real(wp), intent(inout) :: brrate(:,:)
  real(wp), intent(inout) :: frrate(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_reaction_rates"
  !
continue
  !
  call pure_debug_timer(entering_procedure,pname)
  !
  ! Compute the forward reaction rate coefficient
  !
  call compute_arrhenius_rate(n1,n2,t,frrate)
  !
  ! Compute the backwards reaction rate coefficient
  ! using an equilibrium assumption
  !
  if (abs(itestcase) /= CH4_3D_Diffusion_Flame .and. &
      abs(itestcase) /= CH4_Flame_Speed_Test) then
    call compute_equilibrium_rate(n1,n2,t,frrate,brrate)
  end if
  !
  call pure_debug_timer(leaving_procedure,pname)
  !
end subroutine compute_reaction_rates
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_arrhenius_rate(n1,n2,t,frrate)
  !
  ! Routine to calculate the reaction rate coefficient
  ! based on an Arrhenius model.
  !
  ! N1     : beginning point
  ! N2     : ending    point
  ! T      : static temperature
  ! FRRATE : forward reaction rate coefficient
  !
  ! NOTE: The form is that of an Arrehnius model:
  !
  !        k(ir,T) = A(ir)*T^B(ir)*e^(-Ta(ir)/T)
  !
  !        where:
  !
  !        A(ir)  = frcoef(1,ir) = fwd_arrh_coef(ir)
  !        B(ir)  = frcoef(2,ir) = fwd_temp_exp(ir)
  !        Ta(ir) = frcoef(3,ir) = fwd_actv_temp(ir)
  !
  !.. Use Statements ..
  use ovar, only : tref
  !
  ! .. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: t(:)
  real(wp), intent(inout) :: frrate(:,:)
  !
  ! .. Local Scalars ..
  integer  :: icr,n
  real(wp) :: log_a,expnt,tpwr,tactiv,tlog,ts
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_arrhenius_rate"
  !
continue
  !
  call pure_debug_timer(entering_procedure,pname)
  !
  ! Compute the reaction rate coefficient
  !
  do icr = 1,ncr
    !
    log_a = frcoef(1,icr)
    tpwr = frcoef(2,icr)
    tactiv = frcoef(3,icr)
    !
    do n = n1,n2
      !
      ts = t(n)*tref
      tlog = log(ts)
      !
      expnt = log_a + tpwr*tlog - tactiv/ts
      frrate(n,icr) = exp( min(exp_lim,expnt) )
      !
    end do
    !
  end do
  !
  call pure_debug_timer(leaving_procedure,pname)
  !
end subroutine compute_arrhenius_rate
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_equilibrium_rate(n1,n2,t,frrate,brrate)
  !
  ! Routine to calculate the backward reaction rate coefficient
  ! based on knowing the forward reaction rate and using an
  ! equilibrium assumption.
  !
  ! N1     : beginning point
  ! N2     : ending    point
  ! T      : static temperature
  ! FRRATE : forward reaction rate coefficient
  ! BRRATE : backward reaction rate coefficient
  !
  !
  !  NOTE: The form is that of an assumed equilibrium model:
  !
  !        kb(ir,T) = kf(ir,T) / ke(ir,T)
  !
  !        and
  !
  !        ke(ir,T) = (Ro*T)**(-DSCSUM)*exp(-deltaG(ir)/(Ro*T))
  !
  !        where:
  !
  !            -deltaG(ir)/(Ro*T) = sum of the Gibbs free energy times
  !                                 the difference of the stoichiometic
  !                                 coeffifients
  !
  !                                  ns
  !            -deltaG(ir)/(Ro*T) =  sum[-g(is)/(Ro*T)*(delsco(is,ir))]
  !                                  is=1
  !
  !           -g(is)/(Ro*T) = -F(is)*T^-1 +
  !                            H(is)      +
  !                            B(is)*T    +
  !                            C(is)*T^2  +
  !                            D(is)*T^3  +
  !                            E(is)*T^4  +
  !                            A(is)*(ln(T)-1)
  !
  !
  !.. Use Statements ..
  use ovar,       only : tref
  use thermo_mod, only : runivl
  use thermo_mod, only : compute_gibbs_free_energy
  !
  ! .. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: t(:)
  real(wp),    intent(in) :: frrate(:,:)
  real(wp), intent(inout) :: brrate(:,:)
  !
  ! .. Local Scalars ..
  integer  :: icr,m,n,np
  real(wp) :: val,log_a,expnt,tpwr,rulog,tactiv,dscsum_val
  real(wp) :: ts,tlog
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: dgibrt(:)
  real(wp), allocatable :: gibrt(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_equilibrium_rate"
  !
continue
  !
  call pure_debug_timer(entering_procedure,pname)
  !
  rulog = log(runivl)
  !
  np = size(t)
  !
  allocate ( gibrt(1:np,1:ncs) , source=zero )
  allocate ( dgibrt(n1:n2) , source=zero )
  !
  ! Compute the Gibbs free energy / (Runiv*T)
  !
  call compute_gibbs_free_energy(n1,n2,t,gibrt)
  !
  ! Compute backward reaction rate coefficient
  !
  do icr = 1,ncr
    !
    ! Zero dgibrt
    !
    do n = n1,n2
      dgibrt(n) = zero
    end do
    !
    dscsum_val = zero
    !
    do m = 1,ncs
      val = scoerp(2,m,icr) - scoerp(1,m,icr)
      if (abs(val) > eps17) then
        dscsum_val = dscsum_val + val
        do n = n1,n2
          dgibrt(n) = dgibrt(n) + gibrt(n,m)*val
        end do
      end if
    end do
    !
    ! Extract Arrhenius reaction rate parameters
    !
    log_a = frcoef(1,icr)
    tpwr = frcoef(2,icr)
    tactiv = frcoef(3,icr)
    !
    val = dscsum_val
    do n = n1,n2
      !
      ts = t(n)*tref
      tlog = log(ts)
      !
      expnt = dgibrt(n) - val*(tlog+rulog)
      expnt = log_a + tpwr*tlog - tactiv/ts - expnt
      brrate(n,icr) = exp( min(exp_lim,expnt) )
      !
    end do
    !
  end do
  !
  if (allocated(gibrt)) deallocate ( gibrt )
  if (allocated(dgibrt)) deallocate ( dgibrt )
  !
  call pure_debug_timer(leaving_procedure,pname)
  !
end subroutine compute_equilibrium_rate
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_species_production_rates(n1,n2,specon,frrate,brrate,wdot)
  !
  ! Routine to calculate the species production term based on
  ! the law of mass action.
  !                        .
  !  Note: At the end the production term is converted from
  !        [g-mole / (cm**3-sec)] into [kg / (m**3-sec)],
  !        and then non-dimensionalized.
  !
  ! This routine is organized as follows:
  !
  !      1. Correct the concentrations for catalytic bodies if required.
  !      2. Calculate the species production due to the forward reaction
  !      3. Calculate the species production due to the backward reaction
  !      4. Calculate total production rate of the specie
  !
  ! N1     : beginning point
  ! N2     : ending    point
  ! SPECON : species concentrations
  ! FRRATE : forward reaction rate coefficient
  ! BRRATE : backward reaction rate coefficient
  ! WDOT   : species source terms
  !
  !
  !.. Use Statements ..
  use ovar,       only : rhoref,aref
  use thermo_mod, only : wtmol
  !
  ! .. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: brrate(:,:)
  real(wp),    intent(in) :: frrate(:,:)
  real(wp),    intent(in) :: specon(:,:)
  real(wp), intent(inout) :: wdot(:,:)
  !
  ! .. Local Scalars ..
  integer  :: icatb,icatf,icr,m,istoic,n
  real(wp) :: val,rstoic,remain,spcnlg
  !
  !.. Local Allocatable Arrays ..
  ! PRDCTF : production rate of species (forward reaction) [mol/cm^3/s]
  real(wp), allocatable :: prdctf(:)
  ! PRDCTB : production rate of species (backward reaction) [mol/cm^3/s]
  real(wp), allocatable :: prdctb(:)
  ! SPCNCB : species concentration of catalytic bodies [mol/cm^3]
  real(wp), allocatable :: spcncb(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_species_production_rates"
  !
continue
  !
  call pure_debug_timer(entering_procedure,pname)
  !
  allocate ( prdctb(n1:n2) , source=zero )
  allocate ( prdctf(n1:n2) , source=zero )
  !
  ! Compute species production rates in the forward and backward
  ! direction including the effect of catalytic body reactions
  !
  do icr = 1,ncr
    !
    ! Initialize forward and backward prdct arrays
    !
    do n = n1,n2
      prdctf(n) = frrate(n,icr)
    end do
    !
    do n = n1,n2
      prdctb(n) = brrate(n,icr)
    end do
    !
    do m = 1,ncs
      !
      ! Forward reaction
      !
      rstoic = fwd_exp(m,icr)
      istoic = nint(rstoic)
      remain = rstoic - anint(rstoic)
      !
      if (abs(remain) > rndoff) then
        !
        if (abs(rstoic) > rndoff) then
          do n = n1,n2
            if (specon(n,m) > rndoff) then
              spcnlg = log(specon(n,m))
              prdctf(n) = prdctf(n) * exp(rstoic*spcnlg)
            else
              prdctf(n) = zero
            end if
          end do
        end if
        !
      else if (istoic /= 0) then
        !
        do n = n1,n2
          prdctf(n) = prdctf(n) * specon(n,m)**istoic
        end do
        !
      end if
      !
      ! Backward reaction
      !
      rstoic = bwd_exp(m,icr)
      istoic = nint(rstoic)
      remain = rstoic - anint(rstoic)
      !
      if (abs(remain) > rndoff) then
        !
        if (abs(rstoic) > rndoff) then
          do n = n1,n2
            if (specon(n,m) > rndoff) then
              spcnlg = log(specon(n,m))
              prdctb(n) = prdctb(n) * exp(rstoic*spcnlg)
            else
              prdctb(n) = zero
            end if
          end do
        end if
        !
      else if (istoic /= 0) then
        !
        do n = n1,n2
          prdctb(n) = prdctb(n) * specon(n,m)**istoic
        end do
        !
      end if
      !
    end do
    !
    ! Add effects of catalytic body on forward and backward reactions
    !
    icatf = nint(scoerp(1,ncs+1,icr))
    icatb = nint(scoerp(2,ncs+1,icr))
    !
    if (icatf==1 .and. icatb==1) then
      !
      allocate ( spcncb(n1:n2) , source=zero )
      !
      ! Compute the catalytic body concentrations
      !
      val = cbeffc(1,icr)
      do n = n1,n2
        spcncb(n) = val*specon(n,1)
      end do
      !
      do m = 2,ncs
        val = cbeffc(m,icr)
        do n = n1,n2
          spcncb(n) = spcncb(n) + val*specon(n,m)
        end do
      end do
      !
      ! Forward and backward reactions production rates
      !
      do n = n1,n2
        val = spcncb(n)
        prdctf(n) = prdctf(n) * val
        prdctb(n) = prdctb(n) * val
      end do
      !
      if (allocated(spcncb)) deallocate ( spcncb )
      !
    end if
    !
    ! Compute net effect of forward and backward reactions
    !
    do m = 1,ncs
      val = scoerp(2,m,icr) - scoerp(1,m,icr)
      if (abs(val) > eps17) then
        do n = n1,n2
          wdot(n,m) = wdot(n,m) + val*(prdctf(n)-prdctb(n))
        end do
      end if
    end do
    !
  end do
  !
  ! Convert species production rate to kg/(m**3-sec) and non-dimensionalize
  !
  do m = 1,ncs
    !
    ! Set the species production rate to zero if this is a nonreacting species
    !
    if (nonreacting_species(m)) then
      !
      do n = n1,n2
        wdot(n,m) = zero
      end do
      !
    else
      !
      val = wtmol(m) * dens_cgs_to_si  / (rhoref*aref)
      do n = n1,n2
        wdot(n,m) = val * wdot(n,m)
      end do
      !
    end if
    !
  end do
  !
  if (allocated(prdctb)) deallocate ( prdctb )
  if (allocated(prdctf)) deallocate ( prdctf )
  !
  call pure_debug_timer(leaving_procedure,pname)
  !
end subroutine compute_species_production_rates
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine compute_generic_wdot(n1,n2,usp,gpsp,residual)
  !
  !.. Use Statements ..
  use flowvar_mod, only : dtsp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: usp(:,:)
  real(wp),    intent(in) :: gpsp(:,:)
  real(wp), intent(inout) :: residual(:,:)
  !
  !.. Local Scalars ..
  integer :: n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_generic_wdot"
  !
continue
  !
  call pure_debug_timer(entering_procedure,pname)
  !
  call pure_debug_timer(leaving_procedure,pname)
  !
end subroutine compute_generic_wdot
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine semi_implicit_exact_jac(n1,n2,cv,gp,frrate,brrate,specon,wdot)
  !
  ! Routine to perform the point implicit treatment for the source terms
  ! of the species transport equations. There is no turbulence/chemistry
  ! interaction modeling currently being taken into account here, i.e.
  ! this is based on the "laminar" W term.
  !
  ! Note: The order in which events occur is as follows:
  !
  ! 1) The Jacobian of the source term of the species transport
  !    equations is computed.
  !
  ! 2) The partial pivoting  of a Gaussian elimination to solve the
  ! equation,              (GAUSEL)
  !
  ! A * x = B
  ! where :
  !
  ! A - (I-Jacobian*dt) of the source terms of the
  ! species continuity eq.s
  ! B - Untreated Residual of the species continuity eq.s
  ! x - Point implicit treated residual of the
  ! species continuity eq.s
  !
  ! 1
  ! 3) Solve for x = A * B    (GESOLV)
  !
  ! 1) Invert the following using Householders method
  !
  ! A * x = B
  ! where :
  !
  ! A - (I-Jacobian*dt) of the source terms of the
  ! species continuity eq.s
  ! B - Untreated Residual of the species continuity eq.s
  ! x - Point implicit treated residual of the
  ! species continuity eq.s
  !
  ! ############################################################################
  ! ############################################################################
  !
  !   FROM CHEMICAL_JACOBIAN (ORIGINALLY CJACOB)
  !
  ! ----------------------------------------------------------------------------
  !
  ! Routine to calculate the Jacobian of the species production term based
  ! on the law of mass action. There is no turbulence/chemistry interaction
  ! being modeled here.
  !
  ! This is based on the "laminar" W term.
  !
  ! The species concentrations and reaction rate coefficients
  ! are provided to this routine from which it calculates
  ! the terms of the chemical Jacobian.
  !
  ! This routine is organized as follows:
  ! 1. Calculate the species concentration
  ! 2. Correct the concentrations for catalytic bodies if required.
  ! 3. Calculate the elements of the Jacobian of species production
  ! 4. Calculate the matrix (I - dt*Del(W)/Del(Q))
  !
  ! ############################################################################
  ! ############################################################################
  !
  ! N1     : beginning point
  ! N2     : ending    point
  ! CV     : conserved variables
  ! GP     : gas properties
  ! FRRATE : forward reaction rate coefficient
  ! BRRATE : backward reaction rate coefficient
  ! SPECON : species concentrations
  ! WDOT   : residual of the species continuity eq.s
  !
  !.. Use Statements ..
  use ovar,        only : rhoref,aref,tref
  use flowvar_mod, only : dtsp
  !
  use thermo_mod,  only : wtmol
  use thermo_mod,  only : compute_dgbdt
  use thermo_mod,  only : compute_species_energy
  !
  ! Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv(:,:)
  real(wp),    intent(in) :: gp(:,:)
  real(wp),    intent(in) :: brrate(:,:)
  real(wp),    intent(in) :: frrate(:,:)
  real(wp),    intent(in) :: specon(:,:)
  real(wp), intent(inout) :: wdot(:,:)
  !
  ! Local Scalars ..
  integer  :: icatb,icatf,icr,l,istoic,n,np,m
  real(wp) :: const,sumb,sumf,rstoic,remain,spcnlg
  real(wp) :: term,termb,termf,termfb
  real(wp) :: val,dscsum_val,tpwr,tactiv,tinv,dimval
  !
  !
  !.. Local Allocatable Arrays ..
  ! CJAC   : Jacobian of the chemical source terms
  real(wp), allocatable :: cjac(:,:,:)                  ! cjac(1:np,1:ncs,1:ncs)
  ! PRDCTF : production rate of species by the forward reaction
  real(wp), allocatable :: prdctf(:)                    ! prdctf(1:np)
  ! PRDCTB : production rate of species by the backward reaction
  real(wp), allocatable :: prdctb(:)                    ! prdctb(1:np)
  ! SPCNCB : species concentrations to add effect of catalytic bodies
  real(wp), allocatable :: spcncb(:)                    ! spcncb(1:np)
  ! SPECEN : storage for the species internal energys
  real(wp), allocatable :: specen(:,:)                  ! specen(1:np,1:ncs)
  ! DDGDT  : derivative of delta Gibbs free energy w.r.t temperature
  real(wp), allocatable :: ddgdt(:,:)                   ! ddgdt(1:np,1:ncr)
  ! DELSCO : stoichiometric coefficient of the forward  reaction -
  !          stoichiometric coefficient of the backward reaction
  real(wp), allocatable :: delsco(:,:)                  ! delsco(1:ncs,1:ncr)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "semi_implicit_exact_jac"
  !
continue
  !
  call pure_debug_timer(entering_procedure,pname)
  !
  np = size(wdot,dim=1)
  !
  dimval = tref/rhoref
  !
  allocate ( cjac(1:np,1:ncs,1:ncs) , source=zero )
  !
  ! Compute the Jacobian of the species production rate including
  ! including the effect of catalytic body reactions.
  !
  ! Compute species energies
  !
  allocate ( specen(1:np,1:ncs) , source=zero )
  !
  call compute_species_energy(n1,n2,gp(:,ntemp),specen)
  !
  ! Compute derivative of [(-delta G)/(Ru T)]
  !
  allocate ( ddgdt(1:np,1:ncr) , source=zero )
  !
  delsco = scoerp(2,1:ncs,1:ncr) - scoerp(1,1:ncs,1:ncr) ! F2003 AUTO-REALLOC
  !
  call compute_dgbdt(n1,n2,gp(:,ntemp),delsco,ddgdt)
  !
  if (allocated(delsco)) deallocate ( delsco )
  !
  ! Compute the derivative of the species production rates
  ! in the forward and backward direction, including the effect
  ! of catalytic body reactions, w.r.t. species
  !
  allocate ( prdctf(n1:n2) , source=zero )
  allocate ( prdctb(n1:n2) , source=zero )
  allocate ( spcncb(n1:n2) , source=zero )
  !
  reactions_loop: &
  do icr = 1,ncr
    !
    dscsum_val = sum( scoerp(2,1:ncs,icr) - scoerp(1,1:ncs,icr) )
    !
    tpwr = frcoef(2,icr)
    tactiv = frcoef(3,icr)
    !
    outer_species_loop: &
    do m = 1,ncs
      !
     !const = dens_si_to_cgs / wtmol(m)
      const = one / (dens_cgs_to_si * wtmol(m))
      !
      ! Initialize forward and backward prdct arrays
      !
      do n = n1,n2
        prdctf(n) = frrate(n,icr)
        prdctb(n) = brrate(n,icr)
      end do
      !
      do l = 1,ncs
        !
        ! Forward reaction
        !
        if (abs(fwd_exp(l,icr)) > rndoff) then
          !
          if (l == m) then
            rstoic = fwd_exp(l,icr) - one
          else
            rstoic = fwd_exp(l,icr)
          end if
          !
          istoic = nint(rstoic)
          remain = rstoic - anint(rstoic)
          !
          if (abs(remain) > rndoff) then
            if (abs(rstoic) > rndoff) then
              do n = n1,n2
                if (specon(n,l) > rndoff) then
                  spcnlg = log(specon(n,l))
                  prdctf(n) = prdctf(n) * exp(rstoic*spcnlg)
                else
                  prdctf(n) = zero
                end if
              end do
            end if
          else if (istoic /= 0) then
            do n = n1,n2
              prdctf(n) = prdctf(n) * specon(n,l)**istoic
            end do
          end if
          !
        end if
        !
        ! Backward reaction
        !
        if (abs(bwd_exp(l,icr)) > rndoff) then
          !
          if (l == m) then
            rstoic = bwd_exp(l,icr) - one
          else
            rstoic = bwd_exp(l,icr)
          end if
          !
          istoic = nint(rstoic)
          remain = rstoic - anint(rstoic)
          !
          if (abs(remain) > rndoff) then
            if (abs(rstoic) > rndoff) then
              do n = n1,n2
                if (specon(n,l) > rndoff) then
                  spcnlg = log(specon(n,l))
                  prdctb(n) = prdctb(n) * exp(rstoic*spcnlg)
                else
                  prdctb(n) = zero
                end if
              end do
            end if
          else if (istoic /= 0) then
            do n = n1,n2
              prdctb(n) = prdctb(n) * specon(n,l)**istoic
            end do
          end if
          !
        end if
        !
      end do
      !
      ! Add effects of catalytic body on forward and backward reactions
      !
      icatf = nint(scoerp(1,ncs+1,icr))
      icatb = nint(scoerp(2,ncs+1,icr))
      !
      third_body_reaction_check: &
      if (icatf==1 .and. icatb==1) then ! this is a third body reaction
        !
        do n = n1,n2
          spcncb(n) = specon(n,1) * cbeffc(1,icr)
        end do
        !
        do l = 2,ncs
          do n = n1,n2
            spcncb(n) = spcncb(n) + specon(n,l)*cbeffc(l,icr)
          end do
        end do
        !
        ! Derivative of forward reactions production rates
        ! w.r.t. species for catalytic body reactions
        !
        rstoic = fwd_exp(m,icr)
        if (abs(rstoic) > rndoff) then
          term = one
        else
          term = zero
        end if
        !
        do n = n1,n2
          !
          tinv = one / (gp(n,ntemp)*tref)
          !
          termfb = dimval * gp(n,ngamm) * specen(n,m) / cv(n,nec) / gp(n,ncp)
          !
          termf = (tactiv*tinv + tpwr) * tinv
          !
          termf = (term*specon(n,m) + (one-term))*termf
          !
          sumf = (cbeffc(m,icr)*(term*specon(n,m) + (one-term)) &
               + spcncb(n)*rstoic)*const &
               - spcncb(n)*termf*termfb
          !
          prdctf(n) = prdctf(n) * sumf
          !
        end do
        !
        ! Derivative of backward reactions production rates
        ! w.r.t. species for catalytic body reactions
        !
        rstoic = bwd_exp(m,icr)
        if (abs(rstoic) > rndoff) then
          term = one
        else
          term = zero
        end if
        !
        ! Equilibrium backward rate model
        !
        do n = n1,n2
          !
          tinv = one / (gp(n,ntemp)*tref)
          !
          termfb = dimval * gp(n,ngamm) * specen(n,m) / cv(n,nec) / gp(n,ncp)
          !
          termf = (tactiv*tinv + tpwr) * tinv
          !
          termb = termf - ddgdt(n,icr) + dscsum_val*tinv
          !
          termb = (term*specon(n,m) + (one-term))*termb
          !
          sumb = (cbeffc(m,icr)*(term*specon(n,m) + (one-term)) &
               + spcncb(n)*rstoic)*const &
               - spcncb(n)*termb*termfb
          !
          prdctb(n) = prdctb(n) * sumb
          !
        end do
        !
      else third_body_reaction_check ! this is NOT a third body reaction
        !
        ! Derivative of forward reactions production rates
        ! w.r.t. species for non-catalytic body reactions
        !
        rstoic = fwd_exp(m,icr)
        if (abs(rstoic) > rndoff) then
          term = one
        else
          term = zero
        end if
        !
        do n = n1,n2
          !
          tinv = one / (gp(n,ntemp)*tref)
          !
          termfb = dimval * gp(n,ngamm) * specen(n,m) / cv(n,nec) / gp(n,ncp)
          !
          termf = (tactiv*tinv + tpwr) * tinv
          !
          termf = (term*specon(n,m) + (one-term))*termf
          !
          sumf = rstoic*const - termf*termfb
          !
          prdctf(n) = prdctf(n) * sumf
          !
        end do
        !
        ! Derivative of backward reactions production rates
        ! w.r.t. species for non-catalytic body reactions
        !
        rstoic = bwd_exp(m,icr)
        if (abs(rstoic) > rndoff) then
          term = one
        else
          term = zero
        end if
        !
        ! Equilibrium backward rate model
        !
        do n = n1,n2
          !
          tinv = one / (gp(n,ntemp)*tref)
          !
          termfb = dimval * gp(n,ngamm) * specen(n,m) / cv(n,nec) / gp(n,ncp)
          !
          termf = (tactiv*tinv + tpwr) * tinv
          !
          termb = termf - ddgdt(n,icr) + dscsum_val*tinv
          !
          termb = (term*specon(n,m) + (one-term))*termb
          !
          sumb = rstoic*const - termb*termfb
          !
          prdctb(n) = prdctb(n) * sumb
          !
        end do
        !
      end if third_body_reaction_check
      !
      ! Compute net contribution of forward and backward reaction derivatives
      ! if both forward and backward reactions are present
      !
      do l = 1,ncs
        !
        val = scoerp(2,l,icr) - scoerp(1,l,icr)
        !
        if (abs(val) > rndoff) then
          const = val * dens_cgs_to_si * wtmol(l) / aref
          do n = n1,n2
            cjac(n,l,m) = cjac(n,l,m) + const*(prdctf(n)-prdctb(n))
          end do
        end if
        !
      end do
      !
    end do outer_species_loop
    !
  end do reactions_loop
  !
  if (allocated(spcncb)) deallocate ( spcncb )
  if (allocated(prdctb)) deallocate ( prdctb )
  if (allocated(prdctf)) deallocate ( prdctf )
  if (allocated(ddgdt )) deallocate ( ddgdt  )
  if (allocated(specen)) deallocate ( specen )
  !
  ! Multiply the chemical Jacobian matrix by the time step and
  ! add 1 to the diagonal terms.
  !
  do m = 1,ncs
    do l = 1,ncs
      if (m == l) then
        do n = n1,n2
         cjac(n,l,m) = one - dtsp(n)*cjac(n,l,m)
        end do
      else
        do n = n1,n2
         cjac(n,l,m) = -dtsp(n)*cjac(n,l,m)
        end do
      end if
    end do
  end do
  !
  ! Solve for the point implicit residual using Householders method
  !
  call hholdr(n1,n2,cjac,wdot)
  !
  if (allocated(cjac)) deallocate ( cjac )
  !
  call pure_debug_timer(leaving_procedure,pname)
  !
end subroutine semi_implicit_exact_jac
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine semi_implicit_numerical_jac(n1,n2,cv,gp,frrate,brrate, &
                                       specon,wdot,wdotsv)
  !
  ! Routine to perform the point implicit treatment for the source terms
  ! of the species transport equations.
  !
  ! The chemical Jacobian is evaluated numerically.
  !
  ! Note: The order in which events occur is as follows:
  !
  ! NCELLS : number of cell centers
  ! QC     : species mass fraction and total density
  ! QGC    : gas variables at cell centers
  ! SPECON : species concentrations
  ! DGIBRT : delta Gibbs over Runiv2*T
  ! FRRATE : forward reaction rate coefficient
  ! BRRATE : backward reaction rate coefficient
  ! WDOT   : residual of the species continuity eq.s
  !
  !.. Use Statements ..
  use ovar,        only : rhoref,tref
  use flowvar_mod, only : dtsp
  !
  use thermo_mod,  only : wtmol
  use thermo_mod,  only : compute_species_energy
  !
  ! Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv(:,:)
  real(wp),    intent(in) :: gp(:,:)
  real(wp),    intent(in) :: specon(:,:)
  real(wp), intent(inout) :: frrate(:,:)
  real(wp), intent(inout) :: brrate(:,:)
  real(wp), intent(inout) :: wdot(:,:)
  real(wp),    intent(in) :: wdotsv(:,:)
  !
  ! Local Scalars ..
  integer  :: i,m,n,np
  real(wp) :: absol,const,relat
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: perturb(:)
  real(wp), allocatable :: tempnd(:)
  real(wp), allocatable :: wdot2(:,:)
  real(wp), allocatable :: wdot3(:,:)
  real(wp), allocatable :: specen(:,:)
  real(wp), allocatable :: cjac(:,:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "semi_implicit_numerical_jac"
  !
continue
  !
  call pure_debug_timer(entering_procedure,pname)
  !
  ! Compute small numbers
  !
  absol = sqrt(rndoff)
  relat = sqrt(rndoff)
  !
  np = size(wdot,dim=1)
  !
  allocate ( cjac(1:np,1:ncs,1:ncs) , source=zero )
  allocate ( wdot3(n1:n2,1:ncs) , source=zero )
  allocate ( wdot2(n1:n2,1:ncs) , source=zero )
  allocate ( specen(1:np,1:ncs) , source=zero )
  allocate ( perturb(n1:n2) , source=zero )
  !
  !
  ! Compute the temperature independent contribution to the chemical Jacobian
  !
  !
  do m = 1,ncs
    !
    ! Modify species density for numerical derivative
    !
    do n = n1,n2
      ! [kg/m^3]
      perturb(n) = cv(n,m)*rhoref*relat + absol
      ! [1]
      perturb(n) = perturb(n) / rhoref
    end do
    !
    ! Initialize wdot2 array
    !
    do i = 1,ncs
      do n = n1,n2
        wdot2(n,i) = zero
      end do
    end do
    !
    ! Compute the source terms for the modified composition
    !
    ! Modify composition for numerical derivative (specen used for storage)
    !
    do i = 1,ncs
      do n = n1,n2
        specen(n,i) = specon(n,i)
      end do
    end do
    const = rhoref * dens_si_to_cgs / wtmol(m)
    do n = n1,n2
      specen(n,m) = specon(n,m) + perturb(n)*const
    end do
    !
    call compute_species_production_rates(n1,n2,specen,frrate,brrate,wdot2)
    !
    ! The production terms, wdotsv and wdot2: [(kg/m^3 s)(rlref/(rhoref*aref))]
    ! Also, dtsp: [1] ==> [(s) * (aref/rlref)]
    !
    do i = 1,ncs
      do n = n1,n2
        cjac(n,i,m) = -dtsp(n)*(wdot2(n,i) - wdotsv(n,i))/perturb(n)
      end do
    end do
    !
    do n = n1,n2
      cjac(n,m,m) = cjac(n,m,m) + one
    end do
    !
  end do
  !
  !
  ! Compute the temperature dependent contribution to the chemical Jacobian
  !
  !
  ! Initialize wdot2 array
  !
  do m = 1,ncs
    do n = n1,n2
      wdot2(n,m) = zero
    end do
  end do
  !
  ! Modify temperature for numerical derivative
  !
  allocate ( tempnd(1:np) , source=zero )
  !
  do n = n1,n2
    ! [K]
    perturb(n) = gp(n,ntemp)*tref*relat + absol
    ! [1]
    perturb(n) = perturb(n) / tref
    tempnd(n) = gp(n,ntemp) + perturb(n)
  end do
  !
  ! Compute the reaction rate coefficients based on the mean flow
  !
  call compute_reaction_rates(n1,n2,tempnd,frrate,brrate)
  !
  if (allocated(tempnd)) deallocate ( tempnd )
  !
  ! Compute the source terms for the modified temperature
  !
  call compute_species_production_rates(n1,n2,specon,frrate,brrate,wdot2)
  !
  ! Compute wdot2 = (rlref/aref) * {del wdot [kg/(m^3 s] / del T [K)]}
  !
  do m = 1,ncs
    do n = n1,n2
      wdot2(n,m) = (wdot2(n,m)-wdotsv(n,m)) / perturb(n)
    end do
  end do
  !
  if (allocated(perturb)) deallocate ( perturb )
  !
  ! Compute species energies
  !
  call compute_species_energy(n1,n2,gp(:,ntemp),specen)
  !
  !
  !
  do m = 1,ncs
    do n = n1,n2
      wdot3(n,m) = -gp(n,ngamm)*specen(n,m) / (cv(n,nec)*gp(n,ncp))
    end do
  end do
  !
  if (allocated(specen)) deallocate ( specen )
  !
  ! Add temperature partials to the chemical Jacobian
  ! delt * (del wdot_i / del T) * (del T / del rho_o)]
  !
  do m = 1,ncs
    do i = 1,ncs
      do n = n1,n2
        cjac(n,i,m) = cjac(n,i,m) - dtsp(n)*wdot2(n,i)*wdot3(n,m)
      end do
    end do
  end do
  !
  if (allocated(wdot2)) deallocate ( wdot2 )
  if (allocated(wdot3)) deallocate ( wdot3 )
  !
  ! Solve for the point implicit residual using Householders method
  !
  call hholdr(n1,n2,cjac,wdot)
  !
  if (allocated(cjac)) deallocate ( cjac )
  !
  call pure_debug_timer(leaving_procedure,pname)
  !
  ! Finished performing the point implicit treatment for the source terms
  ! of the species transport equations. There is no turbulence/chemistry
  ! interaction modeling currently being taken into account here.
  !
end subroutine semi_implicit_numerical_jac
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine hholdr(n1,n2,a,b)
!subroutine hholdr(a,mmax,m,mf,n,nmax,b,d,ifac,tmp,iblog)
  !
  ! This routine solves Ax=B using the householder method. This routine
  ! should be used if the size of the system is larger than five
  !
  !   This subroutine from HH modified may 1, 1988 in several ways
  !
  !   a. change implementation of householder transformations to use
  !      description given in ortegas new book,  Introduction to
  !      Parallel and Vector Solutions of Linear System which I felt
  !      was slightly more efficient.
  !
  !   b. modified implementation of TMP(M,3) array. now this program will
  !      invoke use of TMP(M,3) if IBLOG (a new parameter) is set true by
  !      the calling program and ignore IB otherwise. if IBLOG has
  !      been set true, the subroutine replaces the contents of the
  !      a array so that the particular system for which TMP(M,3) = 0 will
  !      be set to the identity matrix. this also means that the
  !      output vector B will not change at that point. this setup
  !      represents overhead which in many cases is unnecessary. it
  !      is included in case there are systems in the collection of
  !      M systems which contain garbage which might cause an abort
  !      due to a singular system or overflow from undefined elements.
  !
  !   c. most importantly, I unrolled the outer loops in order to
  !      "enri.h" the inner loops. this decreases the number of
  !      loads and stores, and provides more opportunity for the
  !      compiler to overlap +,* and load/store operations. the
  !      overall effect was to increase the mflops by 20 to 70 %
  !      with the increase the greatest when IBLOG = .false. and
  !      N getting large.
  !
  !   d. the unrolling in HHLE5 is much more efficient than here and
  !      should be used if, in fact, N <= 5 .
  !
  ! description of arguments
  !
  !    A   input array dimensioned as A(MMAX,NMAX,NMAX) which contains
  !        the coefficeints for the M independant N by N systems of
  !        equations to be solved.
  !
  ! MMAX   first dimension of multi-dimensional arrays A, B, and D.
  !
  !    M   actual number of systems to be solved --  M <= MMAX
  !
  !   MF   integer identifying first system (amongst the MMAX systems)
  !        to be solved. the last system to be solved is then
  !        MF + M - 1. vector lengths are M.
  !
  !    N  actual size of system to be solved
  !
  ! NMAX  actual dimension used in the calling programs for arrays A,
  !       B, and D.
  !
  !    B  array dimensioned as B(MMAX,NMAX) used for input and output
  !
  !       input - contains the right-hand-sided elements for the systems
  !               to be solved.
  !
  !       output - contains the solution.
  !                (when iblog = .true. , the B array is unchanged from
  !                 its input value at those elements corresponding to
  !                 TMP(M,3) = 0 since the identity matrix is being solve there
  !
  !    D  array dimensioned D(MMAX,NMAX) used internally by the subroutine.
  !       this array needs to be preserved if the IFAC=1 option is to be
  !       used.
  !
  ! IFAC  input integer
  !
  !       = 0  decompose matrices and solve systems
  !       = 1  do back solve only using information in d and b
  !            preserved from previous call.
  !
  ! TMP(M,1->2) elements of a internally used array dimensioned
  !             as TMP(MMAX,3)
  !
  ! TMP(M,3)    array identifying valid systems to be solved.
  !
  !     TMP(M,3) = 1. system should be solved
  !     TMP(M,3) = 0. ignore this system (attempt to solve it might
  !                   cause an abort - the subroutine will change the
  !                   coefficients in a to the identity matrix and then
  !                   solve it along with the rest of the valid systems.
  !
  ! IBLOG logical variable input from caller
  !        IBLOG = .true. if ib information is to be used
  !                .false. if subroutine is to ignore ib information.
  !                         this is equivalent to setting all IB to 1
  !                         and iblog= .true. but much more efficient.
  !
!subroutine hholdr(a,mmax,m,mf,n,nmax,b,d,ifac,tmp,iblog)
 !!
 !!.. Formal Arguments ..
 !logical(lk), intent(in) :: iblog
 !integer,     intent(in) :: ifac
 !integer,     intent(in) :: m
 !integer,     intent(in) :: mf
 !integer,     intent(in) :: mmax
 !integer,     intent(in) :: n
 !integer,     intent(in) :: nmax
 !real(wp), intent(inout) :: tmp(mmax,3)
 !real(wp), intent(inout) :: b(mmax,nmax)
 !real(wp), intent(inout) :: d(mmax,nmax)
 !real(wp), intent(inout) :: a(mmax,nmax,nmax)
  !
 !sub  hholdr(a         ,mmax  ,m     ,mf,n  ,nmax,b  ,d         ,tmp       )
 !call hholdr(work(iwk0),ncells,ncells,1 ,ncs,ncs ,res,work(iwk1),work(iwk2))
 !call hholdr(cjac      ,ncells,ncells,1 ,ncs,ncs ,res,work(iwk0),work(iwk1))
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: b(:,:)
  real(wp), intent(inout) :: a(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,m,n
  integer  :: kp1,ml,k2,ldo
  integer  :: li,lil,lj,llp1,lm1,lp1,lrem,lvar
  real(wp) :: tm,tm2,tm3,tmp4
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: d(:,:)
  real(wp), allocatable :: tmp(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "hholdr"
  !
continue
  !
  call pure_debug_timer(entering_procedure,pname)
  !
  lvar = size(b,dim=2)
  !
  allocate ( d(n1:n2,1:lvar) , source=zero )
  allocate ( tmp(n1:n2,1:3) , source=zero )
  !
  lp1 = lvar + 1
  lm1 = lvar - 1
  !
  !
  do k = 1,lm1
    !
    do n = n1,n2
      tmp(n,1) = zero
    end do
    !
    do i = k,lvar
      do n = n1,n2
        tmp(n,1) = tmp(n,1) + a(n,i,k)*a(n,i,k)
      end do
    end do
    !
    do n = n1,n2
      !
      tmp4 = -sign(sqrt(tmp(n,1)),a(n,k,k))
      tmp(n,2) = one / (tmp4*(tmp4-a(n,k,k)))
      !
      ! compute u1; store it in a(k,k)
      ! u2,u3  etc. are just remainder of kth col of a
      !
      a(n,k,k) = a(n,k,k) - tmp4
      !
      ! store new a(k,k) into d(k) for back substitution
      !
      d(n,k) = tmp4
      !
    end do
    !
    kp1 = k + 1
    !
    ! modify k+1, k+2, ..,lvar col of a  using u (i.e.,the kth col of a)
    !
    do j = kp1,lvar
      !
      ldo = lvar - k + 1
      lrem = mod(ldo,5)
      k2 = lrem + k - 1
      !
      do n = n1,n2
        tmp(n,1) = zero
      end do
      !
      if (lrem /= 0) then
        do i = k,k2
          do n = n1,n2
            tmp(n,1) = tmp(n,1) + a(n,i,k)*a(n,i,j)
          end do
        end do
      end if
      !
      do i = k2+1,lvar,5
        do n = n1,n2
          tmp(n,1) = tmp(n,1) + a(n,i  ,k)*a(n,i  ,j) &
                                  + a(n,i+1,k)*a(n,i+1,j) &
                                  + a(n,i+2,k)*a(n,i+2,j) &
                                  + a(n,i+3,k)*a(n,i+3,j) &
                                  + a(n,i+4,k)*a(n,i+4,j)
        end do
      end do
      !
      i = k
      !
      do
        !
        lil = lp1 - i
        li = lil - 1
        !
        if (lil == 0) exit
        !
        if (lil == 6) then
          li = 3
        end if
        !
        if (li > 4) then
          li = 4
        end if
        !
        if (li == 1) then
          !
          do n = n1,n2
            tm = tmp(n,1)
            tm = tm * tmp(n,2)
            a(n,i,j) = a(n,i,j) - a(n,i,k)*tm
            a(n,i+1,j) = a(n,i+1,j) - a(n,i+1,k)*tm
          end do
          !
        else if (li == 2) then
          !
          do n = n1,n2
            tm = tmp(n,1)
            tm = tm * tmp(n,2)
            a(n,i,j) = a(n,i,j) - a(n,i,k)*tm
            a(n,i+1,j) = a(n,i+1,j) - a(n,i+1,k)*tm
            a(n,i+2,j) = a(n,i+2,j) - a(n,i+2,k)*tm
          end do
          !
        else if (li == 3) then
          !
          do n = n1,n2
            tm = tmp(n,1)
            tm = tm * tmp(n,2)
            a(n,i,j) = a(n,i,j) - a(n,i,k)*tm
            a(n,i+1,j) = a(n,i+1,j) - a(n,i+1,k)*tm
            a(n,i+2,j) = a(n,i+2,j) - a(n,i+2,k)*tm
            a(n,i+3,j) = a(n,i+3,j) - a(n,i+3,k)*tm
          end do
          !
        else if (li == 4) then
          !
          do n = n1,n2
            tm = tmp(n,1)
            tm = tm * tmp(n,2)
            a(n,i,j) = a(n,i,j) - a(n,i,k)*tm
            a(n,i+1,j) = a(n,i+1,j) - a(n,i+1,k)*tm
            a(n,i+2,j) = a(n,i+2,j) - a(n,i+2,k)*tm
            a(n,i+3,j) = a(n,i+3,j) - a(n,i+3,k)*tm
            a(n,i+4,j) = a(n,i+4,j) - a(n,i+4,k)*tm
          end do
          !
        end if
        !
        i = i + li + 1
        !
      end do
      !
    end do
    !
    ! do same for rhs
    !
    do n = n1,n2
      tmp(n,1) = zero
    end do
    !
    if (lrem /= 0) then
      do i = k,k2
        do n = n1,n2
          tmp(n,1) = tmp(n,1) + a(n,i,k)*b(n,i)
        end do
      end do
    end if
    !
    do i = k2+1,lvar,5
      do n = n1,n2
        tmp(n,1) = tmp(n,1) + a(n,i  ,k)*b(n,i  ) &
                                + a(n,i+1,k)*b(n,i+1) &
                                + a(n,i+2,k)*b(n,i+2) &
                                + a(n,i+3,k)*b(n,i+3) &
                                + a(n,i+4,k)*b(n,i+4)
      end do
    end do
    !
    i = k
    !
    do
      !
      lil = lp1 - i
      li = lil - 1
      !
      if (lil == 0) exit
      !
      if (lil == 6) then
        li = 3
      end if
      !
      if (li > 4) then
        li = 4
      end if
      !
      if (li == 1) then
        !
        do n = n1,n2
          tm = tmp(n,1)
          tm = tm * tmp(n,2)
          b(n,i) = b(n,i) - a(n,i,k)*tm
          b(n,i+1) = b(n,i+1) - a(n,i+1,k)*tm
        end do
        !
      else if (li == 2) then
        !
        do n = n1,n2
          tm = tmp(n,1)
          tm = tm * tmp(n,2)
          b(n,i) = b(n,i) - a(n,i,k)*tm
          b(n,i+1) = b(n,i+1) - a(n,i+1,k)*tm
          b(n,i+2) = b(n,i+2) - a(n,i+2,k)*tm
        end do
        !
      else if (li == 3) then
        !
        do n = n1,n2
          tm = tmp(n,1)
          tm = tm * tmp(n,2)
          b(n,i) = b(n,i) - a(n,i,k)*tm
          b(n,i+1) = b(n,i+1) - a(n,i+1,k)*tm
          b(n,i+2) = b(n,i+2) - a(n,i+2,k)*tm
          b(n,i+3) = b(n,i+3) - a(n,i+3,k)*tm
        end do
        !
      else if (li == 4) then
        !
        do n = n1,n2
          tm = tmp(n,1)
          tm = tm * tmp(n,2)
          b(n,i) = b(n,i) - a(n,i,k)*tm
          b(n,i+1) = b(n,i+1) - a(n,i+1,k)*tm
          b(n,i+2) = b(n,i+2) - a(n,i+2,k)*tm
          b(n,i+3) = b(n,i+3) - a(n,i+3,k)*tm
          b(n,i+4) = b(n,i+4) - a(n,i+4,k)*tm
        end do
        !
      end if
      !
      i = i + li + 1
      !
    end do
    !
  end do
  !
  do n = n1,n2
    d(n,lvar) = a(n,lvar,lvar)
  end do
  !
  if (allocated(tmp)) deallocate ( tmp )
  !
  do j = 2,lvar
    !
    lj = lvar - j + 2
    !
    do n = n1,n2
      b(n,lj) = b(n,lj) / d(n,lj)
    end do
    !
    lil = lvar - j + 1
    llp1 = lil + 1
    !
    do
      !
      li = lil
      !
      if (li == 0) exit
      !
      if (li > 4) then
        li = 4
      end if
      !
      if (li == 1) then
        !
        do n = n1,n2
          b(n,lil) = b(n,lil) - a(n,lil,llp1)*b(n,llp1)
        end do
        !
      else if (li == 2) then
        !
        do n = n1,n2
          b(n,lil-1) = b(n,lil-1) - b(n,llp1)*a(n,lil-1,llp1)
          b(n,lil  ) = b(n,lil  ) - b(n,llp1)*a(n,lil  ,llp1)
        end do
        !
      else if (li == 3) then
        !
        do n = n1,n2
          b(n,lil-2) = b(n,lil-2) - b(n,llp1)*a(n,lil-2,llp1)
          b(n,lil-1) = b(n,lil-1) - b(n,llp1)*a(n,lil-1,llp1)
          b(n,lil  ) = b(n,lil  ) - b(n,llp1)*a(n,lil  ,llp1)
        end do
        !
      else if (li == 4) then
        !
        do n = n1,n2
          b(n,lil-3) = b(n,lil-3) - b(n,llp1)*a(n,lil-3,llp1)
          b(n,lil-2) = b(n,lil-2) - b(n,llp1)*a(n,lil-2,llp1)
          b(n,lil-1) = b(n,lil-1) - b(n,llp1)*a(n,lil-1,llp1)
          b(n,lil  ) = b(n,lil  ) - b(n,llp1)*a(n,lil  ,llp1)
        end do
        !
      end if
      !
      lil = lil - li
      !
    end do
    !
  end do
  !
  do n = n1,n2
    b(n,1) = b(n,1) / d(n,1)
  end do
  !
  if (allocated(d)) deallocate ( d )
  !
  call pure_debug_timer(leaving_procedure,pname)
  !
  ! Finished solving Ax=B using the householder method
  !
end subroutine hholdr
!
!###############################################################################
!
#include "Functions/check_for_bad_number.f90"
!
!###############################################################################
!
end module reactions_mod
