module eqn_idx_mod
  !
  !.. Use Statements ..
  use module_kind_types
  !
  implicit none
  !
  private
  !
  ! Public Procedures
  !
  public :: create_equation_indices
  public :: get_species_indices
  public :: set_reference_coordinate_directions
  !
  ! NCS : number of chemical species
  ! NCR : number of chemical reactions
  !
  integer, public, protected, save :: ncs
  !$ACC DECLARE CREATE(ncs)
  integer, public, protected, save :: ncr
  !$ACC DECLARE CREATE(ncr)
  !
  ! NCB : index for first species continuity equation
  ! NCE : index for last  species continuity equation
  !
  integer, public, protected, save :: ncb
  integer, public, protected, save :: nce
  !$ACC DECLARE CREATE(ncb,nce)
  !
  ! NEC : index for continuity equation
  !
  integer, public, protected, save :: nec
  !$ACC DECLARE CREATE(nec)
  !
  ! NMX : index for x-momentum equation
  ! NMY : index for y-momentum equation
  ! NMZ : index for z-momentum equation
  !
  integer, public, protected, save :: nmx
  integer, public, protected, save :: nmy
  integer, public, protected, save :: nmz
  !$ACC DECLARE CREATE(nmx,nmy,nmz)
  !
  ! NMB : index for first momentum equation
  ! NME : index for last  momentum equation
  !
  integer, public, protected, save :: nmb
  integer, public, protected, save :: nme
  !$ACC DECLARE CREATE(nmb,nme)
  !
  ! NEE : index for energy equation
  !
  integer, public, protected, save :: nee
  !$ACC DECLARE CREATE(nee)
  !
  ! NCD_Axial : index for the axial coordinate direction
  ! NCD_Vert  : index for the vertical coordinate direction
  !             i.e., the coordinate direction at theta = 0 degrees
  ! NCD_Span  : index for the spanwise coordinate direction
  !             i.e., the coordinate direction at theta = 90 degrees
  !
  integer, public, protected, save :: NCD_Axial
  integer, public, protected, save :: NCD_Vert
  integer, public, protected, save :: NCD_Span
  !
  ! NM_Axial : index for the momentum eq in the axial coordinate direction
  ! NM_Vert  : index for the momentum eq in the vertical coordinate direction
  ! NM_Span  : index for the momentum eq in the spanwise coordinate direction
  !
  integer, public, protected, save :: NM_Axial
  integer, public, protected, save :: NM_Vert
  integer, public, protected, save :: NM_Span
  !
  ! NTK : index for turbulence kinetic energy equation
  ! NTL : index for turbulence length scale equation
  !
  integer, public, protected, save :: ntk
  integer, public, protected, save :: ntl
  !$ACC DECLARE CREATE(ntk,ntl)
  !
  ! NQ : total number of equations
  !
  integer, public, protected, save :: nq = 4
  !$ACC DECLARE CREATE(nq)
  !
  ! NG : total number of gas properties
  !       1) pressure
  !       2) temperature
  !       3) gamma
  !       4) viscosity
  !       5) thermal conductivity
  !       6) enthalpy
  !       7) Cp
  !
  integer, public, protected, save :: ng = 7
  !$ACC DECLARE CREATE(ng)
  !
  ! NPRES : index for pressure
  ! NTEMP : index for temperature
  ! NGAMM : index for gamma
  ! NVISC : index for molecular viscosity
  ! NCOND : index for thermal conductivity
  ! NENTH : index for enthalpy
  ! NCP   : index for specific heat at constant pressure
  !
#define NPRES_AND_NTEMP_ARE_CONSTANTS
#ifdef NPRES_AND_NTEMP_ARE_CONSTANTS
  integer, public, parameter :: npres = 1
  integer, public, parameter :: ntemp = 2
  integer, public, parameter :: ngamm = 3
#else
  integer, public, protected, save :: npres = 1
  !$ACC DECLARE CREATE(npres)
  integer, public, protected, save :: ntemp = 2
  !$ACC DECLARE CREATE(ntemp)
  integer, public, protected, save :: ngamm = 3
  !$ACC DECLARE CREATE(ngamm)
#endif
  integer, public, protected, save :: nvisc = 4
  !$ACC DECLARE CREATE(nvisc)
  integer, public, protected, save :: ncond = 5
  !$ACC DECLARE CREATE(ncond)
  integer, public, protected, save :: nenth = 6
  !$ACC DECLARE CREATE(nenth)
  integer, public, protected, save :: ncp   = 7
  !$ACC DECLARE CREATE(ncp)
  !
  !############################################################################
  ! EQUATION INDICES FOR EACH OF THE POSSIBLE CHEMICAL SPECIES
  !----------------------------------------------------------------------------
  integer, public, protected, save :: nCO2  = 0
  !$ACC DECLARE CREATE(nCO2)
  integer, public, protected, save :: nCO   = 0
  !$ACC DECLARE CREATE(nCO)
  integer, public, protected, save :: nH2O  = 0
  !$ACC DECLARE CREATE(nH2O)
  integer, public, protected, save :: nCH4  = 0
  !$ACC DECLARE CREATE(nCH4)
  integer, public, protected, save :: nC2H4 = 0
  !$ACC DECLARE CREATE(nC2H4)
  integer, public, protected, save :: nN2   = 0
  !$ACC DECLARE CREATE(nN2)
  integer, public, protected, save :: nO2   = 0
  !$ACC DECLARE CREATE(nO2)
  integer, public, protected, save :: nH2   = 0
  !$ACC DECLARE CREATE(nH2)
  integer, public, protected, save :: nOH   = 0
  !$ACC DECLARE CREATE(nOH)
  integer, public, protected, save :: nO    = 0
  !$ACC DECLARE CREATE(nO)
  integer, public, protected, save :: nH    = 0
  !$ACC DECLARE CREATE(nH)
  !############################################################################
  !
  ! NGOV : Equation control parameters
  !
  !        NOTE: THIS IS COPIED FROM VULCAN MOST IS NOT USED
  !
  !        NGOV(1,1) = starting index of continuity equations
  !        NGOV(2,1) = ending   index of continuity equations
  !                    and the sum of all densities preceeding
  !        NGOV(3,1) = number of continuity equations
  !                  = 1 when thermally and calorically perfect
  !                  = no. of chemical species+1 when thermally perfect
  !                  = no. of chemical species+1 when thermally imperfect
  !        NGOV(4,1) = starting index of continuity equation source terms
  !        NGOV(5,1) = ending   index of continuity equation source terms
  !        NGOV(6,1) = no. of continuity equation source terms
  !                  = 0 for thermally and calorically perfect gas
  !        NGOV(7,1) = equation RHS solution procedure
  !                    0 = explicit
  !                    1 = point implicit
  !        NGOV(1,2) = starting index of momentum equations
  !        NGOV(2,2) = ending   index of momentum equations
  !        NGOV(3,2) = number of momentum equations
  !        NGOV(6,2) = no. of momentum equation source terms
  !                  = 0
  !        NGOV(7,2) = equation RHS solution procedure
  !                    0 = explicit
  !                    1 = point implicit
  !        NGOV(1,3) = starting index of energy equations
  !        NGOV(2,3) = ending   index of energy equations
  !        NGOV(3,3) = number of energy equations
  !                  = 1 when thermally and calorically perfect
  !                  = 1 when thermally perfect calorically imperfect
  !                  = 2 when thermally imperfect
  !        NGOV(6,3) = no. of energy equation source terms
  !                  = 0 when thermally and calorically perfect
  !                  = 0 when thermally and calorically imperfect
  !                  = N when thermally imperfect
  !        NGOV(1,4) = starting index of turbulence equations
  !        NGOV(2,4) = ending   index of turbulence equations
  !        NGOV(3,4) = number of turbulence equations
  !        NGOV(4,4) = starting index of turbulence equation source terms
  !        NGOV(5,4) = ending   index of turbulence equation source terms
  !        NGOV(6,4) = no. of turbulence equation source terms
  !        NGOV(7,4) = equation RHS solution procedure
  !                    0 = explicit
  !        NGOV(1,5) = starting index of "G" equations
  !        NGOV(2,5) = ending   index of "G" equations
  !        NGOV(3,5) = number of "G" equations
  !        NGOV(4,5) = starting index of "G" equation source terms
  !        NGOV(5,5) = ending   index of "G" equation source terms
  !        NGOV(6,5) = no. of "G" equation source terms
  !        NGOV(7,5) = equation RHS solution procedure
  !                    0 = explicit
  !                    1 = point implicit
  !
  integer, public, protected, save, dimension(3,4) :: ngov
  !
contains
!
!###############################################################################
!
subroutine create_equation_indices(num_dim,itestcase)
  !
  !.. Use Statements ..
  use ovar, only : governing_equations
  use ovar, only : turbulence_model
  use ovar, only : disable_chemical_reactions
  !
  !.. Formal Arguments ..
  integer, intent(in) :: num_dim
  integer, intent(in) :: itestcase
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_equation_indices"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ncs = 0
  ncr = 0
  if (itestcase == Reacting_Backstep) then
    ncs = 2
  else if (itestcase == H2_2D_Diffusion_Flame) then
    if (disable_chemical_reactions) then
      ncs = 3
      ncr = 0
    else
      ncs = 7
      ncr = 7
    end if
  else if (itestcase == Compression_Corner) then
    if (disable_chemical_reactions) then
      ncs = 3
      ncr = 0
    else
      ncs = 7
      ncr = 7
    end if
  else if (abs(itestcase) == CH4_3D_Diffusion_Flame .or. &
           abs(itestcase) == CH4_Flame_Speed_Test) then
    if (disable_chemical_reactions) then
      ncs = 3
      ncr = 0
    else
      ncs = 6
      ncr = 3
    end if
  else if (itestcase == UIUC_Scramjet) then
    if (disable_chemical_reactions) then
      ncs = 3
      ncr = 0
    else
      ncs = 7
      ncr = 3
    end if
  else if (itestcase == N2_O2_Diffusion) then
    ncs = 2
    ncr = 0
  end if
  !
  if (disable_chemical_reactions) then
    ncr = 0
  end if
  !
  ! Continuity equations
  !
  ngov(1,1) = 1
  ngov(2,1) = ngov(1,1) + ncs
  ngov(3,1) = ngov(2,1)
  !
  ! Momentum equations
  !
  ngov(3,2) = num_dim
  ngov(1,2) = ngov(2,1) + 1
  ngov(2,2) = ngov(1,2) + ngov(3,2) - 1
  !
  ! Energy equations
  !
  ngov(3,3) = 1
  ngov(1,3) = ngov(2,2) + 1
  ngov(2,3) = ngov(1,3) + ngov(3,3) - 1
  !
  ! Turbulence equations
  !
  ! Reset the turbulence model to laminar if flow is inviscid
  !
  if (governing_equations == Euler_Eqns) turbulence_model = laminar_flow
  !
  ! Get the number of turbulence model equations
  !
  if (turbulence_model == laminar_flow) then
    ngov(3,4) = 0
  else if (abs(turbulence_model) == k_omega) then
    ngov(3,4) = 2
  end if
  !
  ngov(1,4) = ngov(2,3) + 1
  ngov(2,4) = ngov(1,4) + ngov(3,4) - 1
  !
  ! Set the total number of equations
  !
  nq = ngov(3,1) + ngov(3,2) + ngov(3,3) + ngov(3,4)
  !
  ! Set the indices that will be used through the code
  !
  ncb = ngov(1,1)     ! index for first species continuity equation
  nce = ngov(2,1) - 1 ! index for last  species continuity equation
  !
  nec = ngov(2,1)
  nmx = ngov(1,2)
  nmy = merge( ngov(1,2) + 1 , 0 , ngov(3,2) >= 2 )
  nmz = merge( ngov(1,2) + 2 , 0 , ngov(3,2) == 3 )
  nee = ngov(2,3)
  !
  nmb = nmx
  nme = nec + ngov(3,2)
  !
  ntk = ngov(1,4)
  ntl = ngov(2,4)
  !
  ! Update the indices for the thermodynamic and transport gas properties
  !
#ifndef NPRES_AND_NTEMP_ARE_CONSTANTS
  npres = 1
  ntemp = 2
  ngamm = 3
  !$ACC UPDATE DEVICE(npres,ntemp,ngamm)
#endif
  ng = 3
  !
  if (ncs > 0) then
    !
    ! For a gas mixture, we also need to store enthalpy and Cp
    !
    nenth = ng + 1
    ncp   = ng + 2
    !
  else
    !
    ! For a calorically perfect gas, we dont need to really store gamma,
    ! enthalpy, or Cp. However, it would be convenient to have one auxiliary
    ! spot within the gas properties arrays to store one of these variables as
    ! needed. This would make it easier to create simple generic loops that work
    ! for both gas mixtures and a calorically perfect gas. To utilize this
    ! "auxiliary" variable, increase the number of gas properties variables by
    ! one and have ngamm, nenth, and ncp all point to the last variable.
    !
    nenth = ngamm
    ncp   = ngamm
    !
  end if
  !
  ng = ncp
  !
  ! Add indices for transport properties if necessary
  !
  if (governing_equations == NavierStokes_Eqns) then
    nvisc = ng + 1
    ncond = ng + 2
    ng = ncond
  else
    nvisc = 0
    ncond = 0
  end if
  !
  !$ACC UPDATE DEVICE(nq,nec,nmx,nmy,nmz,nee,nmb,nme,ntk,ntl)
  !$ACC UPDATE DEVICE(ncs,ncb,nce,ncr)
  !$ACC UPDATE DEVICE(ng,nvisc,ncond,ngamm,nenth,ncp)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine create_equation_indices
!
!###############################################################################
!
subroutine get_species_indices(species_name)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: species_name(:)
  !
  !.. Local Scalars ..
  integer :: m
  character(len=10) :: species
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_species_indices"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  do m = 1,size(species_name)
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
  !$ACC UPDATE DEVICE(nCO2,nCO,nH2O,nCH4,nC2H4,nN2,nO2,nH2,nOH,nO,nH)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_species_indices
!
!###############################################################################
!
subroutine set_reference_coordinate_directions()
  !
  !.. Use Statements ..
  use ovar, only : alpha_aoaref,beta_aoaref
  !
continue
  !
  ! Default to X being the axial coordinate direction
  !
  NCD_Axial = 1                       ! set X as axial    coordinate direction
  NCD_Vert  = merge( 2 , 0 , nmy > 0) ! set Y as vertical coordinate direction
  NCD_Span  = merge( 3 , 0 , nmz > 0) ! set Z as spanwise coordinate direction
  !
  NM_Axial = nmx ! set X-momentum as axial    momentum equation
  NM_Vert  = nmy ! set Y-momentum as vertical momentum equation
  NM_Span  = nmz ! set Z-momentum as spanwise momentum equation
  !
  ! Change the reference coordinate directions for
  ! specific reference angle of attack conditions
  !
  if (nmz == 0 .and. nmy > 0) then ! this is a 2D flow
    !
    if (abs(alpha_aoaref-nine*ten) < eps6) then
      !
      ! If alpha_aoaref is 90 degrees, set Y as the axial coordinate direction
      !
      NCD_Axial = 2  ! set Y as axial    coordinate direction
      NCD_Vert  = 1  ! set X as vertical coordinate direction
      !
      NM_Axial = nmy ! set Y-momentum as axial    momentum equation
      NM_Vert  = nmx ! set X-momentum as vertical momentum equation
      !
    end if
    !
  else if (nmz > 0) then ! this is a 3D flow
    !
    if (all(abs([beta_aoaref-nine*ten,alpha_aoaref]) < eps6)) then
      !
      ! If beta_aoaref is 90 degrees and alpha_aoaref is 0 degrees,
      ! set Z as the axial coordinate direction
      !
      NCD_Axial = 3  ! set Z as axial    coordinate direction
      NCD_Vert  = 1  ! set X as vertical coordinate direction
      NCD_Span  = 2  ! set Y as spanwise coordinate direction
      !
      NM_Axial = nmz ! set Z-momentum as axial    momentum equation
      NM_Vert  = nmx ! set X-momentum as vertical momentum equation
      NM_Span  = nmy ! set Y-momentum as spanwise momentum equation
      !
    else if (all(abs([alpha_aoaref-nine*ten,beta_aoaref]) < eps6)) then
      !
      ! If alpha_aoaref is 90 degrees and beta_aoaref is 0 degrees,
      ! set Y as the axial coordinate direction
      !
      NCD_Axial = 2  ! set Y as axial    coordinate direction
      NCD_Vert  = 3  ! set Z as vertical coordinate direction
      NCD_Span  = 1  ! set X as spanwise coordinate direction
      !
      NM_Axial = nmy ! set Y-momentum as axial    momentum equation
      NM_Vert  = nmz ! set Z-momentum as vertical momentum equation
      NM_Span  = nmx ! set X-momentum as spanwise momentum equation
      !
    end if
    !
  end if
  !
end subroutine set_reference_coordinate_directions
!
!###############################################################################
!
end module eqn_idx_mod
!
!###############################################################################
!
