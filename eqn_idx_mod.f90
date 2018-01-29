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
  ! NEC : index for continuity equation
  !
  integer, public, protected, save :: nec
  !
  ! NMX : index for x-momentum equation
  ! NMY : index for y-momentum equation
  ! NMZ : index for z-momentum equation
  !
  integer, public, protected, save :: nmx
  integer, public, protected, save :: nmy
  integer, public, protected, save :: nmz
  !
  ! NMB : index for first momentum equation
  ! NME : index for last  momentum equation
  !
  integer, public, protected, save :: nmb
  integer, public, protected, save :: nme
  !
  ! NEC : index for energy equation
  !
  integer, public, protected, save :: nee
  !
  ! NTK : index for turbulence kinetic energy equation
  ! NTL : index for turbulence length scale equation
  !
  integer, public, protected, save :: ntk
  integer, public, protected, save :: ntl
  !
  ! NQ : total number of equations
  !
  integer, public, protected, save :: nq = 4
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
subroutine create_equation_indices(num_dim)
  !
  !.. Use Statements ..
  use ovar, only : governing_equations
  use ovar, only : turbulence_model
  !
  !.. Formal Arguments ..
  integer, intent(in) :: num_dim
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
  ! Continuity equations
  !
  ngov(1:3,1) = 1
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
  nec = ngov(2,1)
  nmx = ngov(1,2)
  nmy = ngov(1,2) + 1
  nmz = merge( ngov(1,2) + 2 , 0 , num_dim == 3 )
  nee = ngov(2,3)
  !
  nmb = nmx
  nme = nec + ngov(3,2)
  !
  ntk = ngov(1,4)
  ntl = ngov(2,4)
  !
end subroutine create_equation_indices
!
!###############################################################################
!
end module eqn_idx
!
!###############################################################################
!
