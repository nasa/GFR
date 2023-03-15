module eqn_idx
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
  public :: set_auxiliary_variable
  public :: aux_var_is_not_gamma
  public :: aux_var_is_not_enthalpy
  public :: aux_var_is_not_cp
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
  !       3) viscosity
  !       4) thermal conductivity
  !       5) gamma
  !       6) enthalpy
  !       7) Cp
  !
  integer, public, protected, save :: ng = 7
  !$ACC DECLARE CREATE(ng)
  !
  ! NPRES : index for pressure
  ! NTEMP : index for temperature
  ! NVISC : index for molecular viscosity
  ! NCOND : index for thermal conductivity
  ! NGAMM : index for gamma
  ! NENTH : index for enthalpy
  ! NCP   : index for specific heat at constant pressure
  !
#define NPRES_AND_NTEMP_ARE_CONSTANTS
#ifdef NPRES_AND_NTEMP_ARE_CONSTANTS
  integer, public, parameter :: npres = 1
  integer, public, parameter :: ntemp = 2
#else
  integer, public, protected, save :: npres = 1
  !$ACC DECLARE CREATE(npres)
  integer, public, protected, save :: ntemp = 2
  !$ACC DECLARE CREATE(ntemp)
#endif
  integer, public, protected, save :: nvisc = 3
  !$ACC DECLARE CREATE(nvisc)
  integer, public, protected, save :: ncond = 4
  !$ACC DECLARE CREATE(ncond)
  integer, public, protected, save :: ngamm = 5
  !$ACC DECLARE CREATE(ngamm)
  integer, public, protected, save :: nenth = 6
  !$ACC DECLARE CREATE(nenth)
  integer, public, protected, save :: ncp   = 7
  !$ACC DECLARE CREATE(ncp)
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
  integer, private, save :: current_aux_variable = 0
  !
  type :: auxvar_t
    character(len=10) :: name
  end type auxvar_t
  !
  type(auxvar_t), public, parameter :: aux_gamma = auxvar_t("gamma")
  type(auxvar_t), public, parameter :: aux_enthalpy = auxvar_t("enthalpy")
  type(auxvar_t), public, parameter :: aux_cp = auxvar_t("Cp")
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
  !$ACC UPDATE DEVICE(npres,ntemp)
#endif
  ngamm = 3
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
subroutine set_auxiliary_variable(auxvar)
  !
  !.. Formal Arguments ..
  type(auxvar_t), intent(in) :: auxvar
  !
continue
  !
  ! NOTE: For now we are keeping current_aux_variable to always be 0 so that
  !       the code will never think that the desired auxiliary variable is
  !       currently in the auxiliary variable slot. This is because this logic
  !       works only if we have a single gas properties array, which is NOT the
  !       case. We have gpsp, facegpsp, facegpsp_dt, edgegpsp_dt, nodegpsp_dt,
  !       and possibly others. Therefore, we really need to have some way to
  !       track the current auxiliary variable for each separate array. Again,
  !       for now we are keeping current_aux_variable to always be 0 so that
  !       the code will never think that the desired auxiliary variable is
  !       currently in the auxiliary variable slot.
  !
  current_aux_variable = 0
 !if (ncs > 0) then
 !  current_aux_variable = 0
 !else
 !  if (auxvar%name = "gamma") then
 !    current_aux_variable = 1
 !  else if (auxvar%name = "enthalpy") then
 !    current_aux_variable = 2
 !  else if (auxvar%name = "Cp") then
 !    current_aux_variable = 3
 !  end if
 !end if
  !
end subroutine set_auxiliary_variable
!
!###############################################################################
!
pure function aux_var_is_not_gamma() result(return_value)
  logical(ldk) :: return_value
continue
  if (ncs == 0) then
    return_value = (current_aux_variable /= 1)
  else
    return_value = .true.
  end if
end function aux_var_is_not_gamma
!
!###############################################################################
!
pure function aux_var_is_not_enthalpy() result(return_value)
  logical(ldk) :: return_value
continue
  if (ncs == 0) then
    return_value = (current_aux_variable /= 2)
  else
    return_value = .true.
  end if
end function aux_var_is_not_enthalpy
!
!###############################################################################
!
pure function aux_var_is_not_cp() result(return_value)
  logical(ldk) :: return_value
continue
  if (ncs == 0) then
    return_value = (current_aux_variable /= 3)
  else
    return_value = .true.
  end if
end function aux_var_is_not_cp
!
!###############################################################################
!
end module eqn_idx
!
!###############################################################################
!
