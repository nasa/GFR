module reactions_mod
  !
  !.. Use Statements ..
  use module_kind_types
  !
  use eqn_idx, only : nec,nq
  use eqn_idx, only : ncs,ntemp,ng
  use eqn_idx, only : ncr
  !
  implicit none
  !
  private
  !
  public :: init_reactions_mod
  public :: compute_wdot
  !
  !
  !
  real(wp), save, allocatable :: wdot(:,:)
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
  ! N_REAC_SUB_ITER: # of sub-iterations for computing
  !                  source terms from reaction kinetics
  !
  integer, save :: n_reac_sub_iter = 10
  !
  !
  !
  integer, save :: nCO2 = 0
  integer, save :: nCO  = 0
  integer, save :: nH2O = 0
  integer, save :: nCH4 = 0
  integer, save :: nN2  = 0
  integer, save :: nO2  = 0
  !
  integer, save :: nH2  = 0
  integer, save :: nOH  = 0
  integer, save :: nO   = 0
  integer, save :: nH   = 0
  !
contains
!
!###############################################################################
!
subroutine init_reactions_mod()
  !
  !.. Use Statements ..
  use ovar, only : tref
  use ovar, only : itestcase
  use ovar, only : chemistry_sub_iterations
  !
  use thermo_mod, only : wtmol
  use thermo_mod, only : species_name
  !
  !.. Local Scalars ..
  integer  :: i,j,m,ir,ierr
  real(wp) :: fwd_mwp,bwd_mwp
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
  !
  !
  allocate ( fwd_arrh_coef(1:ncr) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fwd_arrh_coef",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( bwd_arrh_coef(1:ncr) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bwd_arrh_coef",1,__LINE__,__FILE__,ierr,error_message)
  !
  !
  !
  allocate ( fwd_actv_temp(1:ncr) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fwd_actv_temp",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( bwd_actv_temp(1:ncr) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bwd_actv_temp",1,__LINE__,__FILE__,ierr,error_message)
  !
  !
  !
  allocate ( fwd_exp(1:ncs,1:ncr) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fwd_exp",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( bwd_exp(1:ncs,1:ncr) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bwd_exp",1,__LINE__,__FILE__,ierr,error_message)
  !
  !
  !
  allocate ( fwd_stoich_coef(1:ncs,1:ncr) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fwd_stoich_coef",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  allocate ( bwd_stoich_coef(1:ncs,1:ncr) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bwd_stoich_coef",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  allocate ( stoich_coef(1:ncs,1:ncr) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"stoich_coef",1,__LINE__,__FILE__,ierr,error_message)
  !
  !
  !
  if (itestcase == CH4_3D_Diffusion_Flame) then
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
    fwd_arrh_coef(ir) = 1.0e+15_wp
    fwd_actv_temp(ir) = 35000.0_wp
    !
    bwd_arrh_coef(ir) = zero
    bwd_actv_temp(ir) = zero
    !
    fwd_exp(nCH4,ir) = 0.9_wp ! CH4
    fwd_exp(nO2 ,ir) = 1.1_wp ! O2
    !
    fwd_stoich_coef(nCH4,ir) = two   ! CH4
    fwd_stoich_coef(nO2 ,ir) = three ! O2
    !
    bwd_stoich_coef(nCO ,ir) = two   ! CO
    bwd_stoich_coef(nH2O,ir) = four  ! H2O
    !
    !##############################################
    !
    ! Reaction 2:  2CO2 => 2CO + O2
    !
    ir = 2
    !
    fwd_arrh_coef(ir) = 4.0552e+10_wp
    fwd_actv_temp(ir) = 77194.0_wp
    !
    bwd_arrh_coef(ir) = zero
    bwd_actv_temp(ir) = zero
    !
    fwd_exp(nCO2,ir) = one ! CO2
    !
    fwd_stoich_coef(nCO2,ir) = two   ! CO2
    !
    bwd_stoich_coef(nCO ,ir) = two   ! CO
    bwd_stoich_coef(nO2 ,ir) = one   ! O2
    !
    !##############################################
    !
    ! Reaction 3:  2CO + O2 => 2CO2
    !
    ir = 3
    !
    fwd_arrh_coef(ir) = 1.0e+9_wp
    fwd_actv_temp(ir) = 12000.0_wp
    !
    bwd_arrh_coef(ir) = zero
    bwd_actv_temp(ir) = zero
    !
    fwd_exp(nCO ,ir) = one  ! CO
    fwd_exp(nO2 ,ir) = half ! O2
    !
    fwd_stoich_coef(nCO ,ir) = two   ! CO
    fwd_stoich_coef(nO2 ,ir) = one   ! O2
    !
    bwd_stoich_coef(nCO2,ir) = two   ! CO2
    !
  else if (itestcase == H2_2D_Diffusion_Flame) then
    !
    !
  else
    !
    write (error_message,1)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    !
  end if
  !
  !################################################
  !
  ! Pre-divide the forwards and backwards Arrhenius coefficients by
  !      product( wtmol(1:ncs)**fwd_exp(1:ncs,j) )
  !      product( wtmol(1:ncs)**bwd_exp(1:ncs,j) )
  !
  do j = 1,ncr
    fwd_mwp = one
    bwd_mwp = one
    do i = 1,ncs
      fwd_mwp = fwd_mwp * (wtmol(i)**fwd_exp(i,j))
      bwd_mwp = bwd_mwp * (wtmol(i)**bwd_exp(i,j))
    end do
    fwd_arrh_coef(j) = fwd_arrh_coef(j)/fwd_mwp
    bwd_arrh_coef(j) = bwd_arrh_coef(j)/bwd_mwp
  end do
  !
  ! Pre-multiply the stoichiometric coefficients by the molar mass
  !
  do j = 1,ncr
    do i = 1,ncs
      fwd_stoich_coef(i,j) = fwd_stoich_coef(i,j)*wtmol(i)
      bwd_stoich_coef(i,j) = bwd_stoich_coef(i,j)*wtmol(i)
    end do
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
  ! Nondimensionalize the activation temperatures
  !
  do j = 1,ncr
    fwd_actv_temp(j) = fwd_actv_temp(j)/tref
    bwd_actv_temp(j) = bwd_actv_temp(j)/tref
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  1 format ("ERROR! Test case is not valid with chemical reactions! ERROR!")
  !
end subroutine init_reactions_mod
!
!###############################################################################
!
subroutine compute_wdot(n1,n2,usp,gpsp,residual)
  !
  !.. Use Statements ..
  use ovar, only : itestcase
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
  character(len=*), parameter :: pname = "compute_wdot"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !'
  if (itestcase == CH4_3D_Diffusion_Flame) then
    call compute_CH4_hardcoded_wdot(n1,n2,usp,gpsp,residual)
  else if (itestcase == H2_2D_Diffusion_Flame) then
    call compute_H2_hardcoded_wdot(n1,n2,usp,gpsp,residual)
  else
    call compute_generic_wdot(n1,n2,usp,gpsp,residual)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_wdot
!
!###############################################################################
!
subroutine compute_CH4_hardcoded_wdot(n1,n2,usp,gpsp,residual)
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
  use flowvar,    only : dtsp
  use thermo_mod, only : species_name
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: usp(:,:)
  real(wp),    intent(in) :: gpsp(:,:)
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
  real(wp), parameter :: dens_cgs_to_si = 1000.0_wp
  real(wp), parameter :: dens_si_to_cgs = one/dens_cgs_to_si
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  !
  !
  allocate ( rhoi(n1:n2,1:ncs) )
  allocate ( rho_ratio(n1:n2) )
  !
  ! Get the species densities in g/cm^3
  !
  do m = 1,ncs
    do n = n1,n2
      rhoi(n,m) = usp(n,m)*dens_si_to_cgs
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
      temp = gpsp(n,ntemp)
      !
      !=========================================================================
      !=========================================================================
      !=========================================================================
      !
      !    REACTION 1: 2CH4 + 3O2 => 2CO + 4H2O
      !
      !-------------------------------------------------------------------------
      !
      ir = 1
      !
      ! Arrhenius coefficient pre-multiplied by 1/product(wtmol(:)**fwd_exp(:))
      arrh = fwd_arrh_coef(ir)
      ! Activation temperature, nondimensionalized by tref
      tact = fwd_actv_temp(ir)
      !
      ! Compute the forward kinetic reaction rate
      !
      qf = arrh*exp(-tact/temp) &
         * (rhoi(n,nCH4)**fwd_exp(nCH4,ir)) &
         * (rhoi(n,nO2 )**fwd_exp(nO2 ,ir))
      !
      !  fwd_stoich_coef(nCH4,1) = 2 pre-multiplied by wtmol(nCH4)
      !  fwd_stoich_coef(nO2 ,1) = 3 pre-multiplied by wtmol(nO2 )
      !
      !  bwd_stoich_coef(nCO ,1) = 2 pre-multiplied by wtmol(nCO )
      !  bwd_stoich_coef(nH2O,1) = 4 pre-multiplied by wtmol(nH2O)
      !
      !  stoich_coef(nCH4,1) = -2 pre-multiplied by wtmol(nCH4)
      !  stoich_coef(nO2 ,1) = -3 pre-multiplied by wtmol(nO2 )
      !  stoich_coef(nCO ,1) =  2 pre-multiplied by wtmol(nCO )
      !  stoich_coef(nH2O,1) =  4 pre-multiplied by wtmol(nH2O)
      !
      if (qf > zero) then
        !
        ! If qf > 0, then rhoi(nCH4) > 0 and rhoi(nO2) > 0
        !
        q1 = qf*fwd_stoich_coef(nCH4,ir)/rhoi(n,nCH4)
        q2 = qf*fwd_stoich_coef(nO2 ,ir)/rhoi(n,nO2 )
        !
        qref = max( q1 , q2 )
        !
        domega = dt_chem * qf / (one + dt_chem*qref)
        !
        rhoi(n,nCO ) = max( zero , rhoi(n,nCO ) + stoich_coef(nCO ,ir)*domega )
        rhoi(n,nH2O) = max( zero , rhoi(n,nH2O) + stoich_coef(nH2O,ir)*domega )
        rhoi(n,nCH4) = max( zero , rhoi(n,nCH4) + stoich_coef(nCH4,ir)*domega )
        rhoi(n,nO2 ) = max( zero , rhoi(n,nO2 ) + stoich_coef(nO2 ,ir)*domega )
        !
      end if
      !
      !=========================================================================
      !=========================================================================
      !=========================================================================
      !
      !    REACTION 2: 2CO2 => 2CO + O2
      !
      !-------------------------------------------------------------------------
      !
      ir = 2
      !
      ! Arrhenius coefficient pre-multiplied by 1/product(wtmol(:)**fwd_exp(:))
      arrh = fwd_arrh_coef(ir)
      ! Activation temperature, nondimensionalized by tref
      tact = fwd_actv_temp(ir)
      !
      ! Compute the forward kinetic reaction rate
      !
      qf = arrh*exp(-tact/temp) &
         * (rhoi(n,nCO2)**fwd_exp(nCO2,ir))
      !
      !  fwd_stoich_coef(nCO2,2) = 2 pre-multiplied by wtmol(nCO2)
      !
      !  bwd_stoich_coef(nCO ,2) = 2 pre-multiplied by wtmol(nCO )
      !  bwd_stoich_coef(nO2 ,2) = 1 pre-multiplied by wtmol(nO2 )
      !
      !  stoich_coef(nCO2,2) = -2 pre-multiplied by wtmol(nCO2)
      !  stoich_coef(nCO ,2) =  2 pre-multiplied by wtmol(nCO )
      !  stoich_coef(nO2 ,2) =  1 pre-multiplied by wtmol(nO2 )
      !
      if (qf > zero) then
        !
        ! If qf > 0, then rhoi(nCO2) > 0
        !
        qref = qf*fwd_stoich_coef(nCO2,ir)/rhoi(n,nCO2)
        !
        domega = dt_chem * qf / (one + dt_chem*qref)
        !
        rhoi(n,nCO2) = max( zero , rhoi(n,nCO2) + stoich_coef(nCO2,ir)*domega )
        rhoi(n,nCO ) = max( zero , rhoi(n,nCO ) + stoich_coef(nCO ,ir)*domega )
        rhoi(n,nO2 ) = max( zero , rhoi(n,nO2 ) + stoich_coef(nO2 ,ir)*domega )
        !
      end if
      !
      !=========================================================================
      !=========================================================================
      !=========================================================================
      !
      !    REACTION 3: 2CO + O2 => 2CO2
      !
      !-------------------------------------------------------------------------
      !
      ir = 3
      !
      ! Arrhenius coefficient pre-multiplied by 1/product(wtmol(:)**fwd_exp(:))
      arrh = fwd_arrh_coef(ir)
      ! Activation temperature, nondimensionalized by tref
      tact = fwd_actv_temp(ir)
      !
      ! Compute the forward kinetic reaction rate
      !
      qf = arrh*exp(-tact/temp) &
         * (rhoi(n,nCO)**fwd_exp(nCO,ir)) &
         * (rhoi(n,nO2)**fwd_exp(nO2,ir))
      !
      !  fwd_stoich_coef(nCO ,3) = 2 pre-multiplied by wtmol(nCO )
      !  fwd_stoich_coef(nO2 ,3) = 1 pre-multiplied by wtmol(nO2 )
      !
      !  bwd_stoich_coef(nCO2,3) = 2 pre-multiplied by wtmol(nCO2)
      !
      !  stoich_coef(nCO2,3) =  2 pre-multiplied by wtmol(nCO2)
      !  stoich_coef(nCO ,3) = -2 pre-multiplied by wtmol(nCO )
      !  stoich_coef(nO2 ,3) = -1 pre-multiplied by wtmol(nO2 )
      !
      if (qf > zero) then
        !
        ! If qf > 0, then rhoi(nCO) > 0 and rhoi(nO2) > 0
        !
        q1 = qf*fwd_stoich_coef(nCO,ir)/rhoi(n,nCO)
        q2 = qf*fwd_stoich_coef(nO2,ir)/rhoi(n,nO2)
        !
        qref = max( q1 , q2 )
        !
        domega = dt_chem * qf / (one + dt_chem*qref)
        !
        rhoi(n,nCO2) = max( zero , rhoi(n,nCO2) + stoich_coef(nCO2,ir)*domega )
        rhoi(n,nCO ) = max( zero , rhoi(n,nCO ) + stoich_coef(nCO ,ir)*domega )
        rhoi(n,nO2 ) = max( zero , rhoi(n,nO2 ) + stoich_coef(nO2 ,ir)*domega )
        !
      end if
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
      rho_ratio(n) = zero
    end do
    !
    do m = 1,ncs
      do n = n1,n2
        rho_ratio(n) = rho_ratio(n) + rhoi(n,m)
      end do
    end do
    !
    ! Compute the ratio of old mixture density to new mixture density
    !
    do n = n1,n2
      rho_ratio(n) = usp(n,nec)*dens_si_to_cgs/rho_ratio(n)
    end do
    !
    ! Initialize mass density {g/cm^3} for each species
    !
    ! Get the species densities in g/cm^3 and prevent negative density
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
  ! Compute the updated species densities (in kg/m^3)
  ! resulting from the reaction sub-iterations
  !
  do m = 1,ncs
    do n = n1,n2
      rhoi(n,m) = rhoi(n,m)*dens_cgs_to_si * rho_ratio(n)
    end do
  end do
  !
  ! Add the change in species densities to the RHS residual
  !
  do m = 1,ncs
    if (species_name(m) /= 'N2') then
      do n = n1,n2
       !dt_inv = one/dtsp(n)
       !omega_dot(n,m) = dt_inv*(rhoi(n,m)*rho_ratio(n) - usp(n,m))/usp(n,nec)
        residual(n,m) = residual(n,m) + (rhoi(n,m) - usp(n,m))
      end do
    end if
  end do
  !
  !
  !
  deallocate ( rho_ratio )
  deallocate ( rhoi )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_CH4_hardcoded_wdot
!
!###############################################################################
!
subroutine compute_H2_hardcoded_wdot(n1,n2,usp,gpsp,residual)
  !
  !.. Use Statements ..
  use flowvar, only : dtsp
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
  character(len=*), parameter :: pname = "compute_H2_hardcoded_wdot"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_H2_hardcoded_wdot
!
!###############################################################################
!
subroutine compute_generic_wdot(n1,n2,usp,gpsp,residual)
  !
  !.. Use Statements ..
  use flowvar, only : dtsp
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
  call debug_timer(entering_procedure,pname)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_generic_wdot
!
!###############################################################################
!
end module reactions_mod
