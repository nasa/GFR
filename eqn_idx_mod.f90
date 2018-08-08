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
  integer, public, protected, save :: npres = 1
  !$ACC DECLARE CREATE(npres)
  integer, public, protected, save :: ntemp = 2
  !$ACC DECLARE CREATE(ntemp)
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
  if (num_dim == 1) then
    !
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__, &
                    "1D is currently not supported")
    !
  end if
  !
  ncs = 0
  ncr = 0
  if (itestcase == Reacting_Backstep) then
    ncs = 2
  else if (itestcase == H2_2D_Diffusion_Flame) then
    ncs = 7
    ncr = 7
  else if (itestcase == CH4_3D_Diffusion_Flame) then
    ncs = 6
    ncr = 3
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
  !$ACC UPDATE DEVICE(nq,ncs,ncb,nce,nec,nmx,nmy,nmz,nee,nmb,nme,ntk,ntl)
  !
end subroutine create_equation_indices
!
!###############################################################################
!
end module eqn_idx
!
!###############################################################################
!
