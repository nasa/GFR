!###############################################################################
!> Module for the bc_input data structure in version 0001 namelist.
!###############################################################################
module bc_input_0001_mod
  !
  use module_kind_types, only : wp,lk,fals,true
  use module_kind_types, only : zero,one,eps6
  use module_kind_types, only : bc_unknown
  !
  use ovar, only : relaxation_t
  use ovar, only : sem_bc_t
  use ovar, only : periodic_bc_t
  use ovar, only : ovar_bc_input_t => bc_input_t
  !
  implicit none
  !
  private
  !
  !
  !
  ! ##########################
  ! ###  Module Functions  ###
  ! ##########################
  !
  private :: is_not_default_bc_input
  private :: has_default_flow_conditions
  private :: bc_input_t_equals
  private :: bc_input_t_not_equal
  !
  public :: bc_input_v0001_to_ovar
  public :: bc_input_ovar_to_v0001
  !
  !
  ! ##############################
  ! ###  Module Derived Types  ###
  ! ##############################
  !
  !  bc_input_t :
  !> Derived type used in the input namelist to specify
  !> the flow parameters for each boundary condition.
  type, public :: bc_input_t
    !
    character(len=32) :: name = ""
    !
    integer :: set_bc_value = bc_unknown
    !
    real(wp) :: t_static   = -one
    real(wp) :: p_static   = -one
    real(wp) :: rho_static = -one
    real(wp) :: mach       = -one
    real(wp) :: vx         = zero
    real(wp) :: vy         = zero
    real(wp) :: vz         = zero
    real(wp) :: t_total    = -one
    real(wp) :: p_total    = -one
    real(wp) :: rho_total  = -one
    real(wp) :: alpha_aoa  = zero
    real(wp) :: beta_aoa   = zero
    real(wp) :: wall_temp  = -one
    !
    ! NOTE: WE NEED TO MAKE SURE WE ADD LOGIC THAT SETS ROTATION_RATE(1)
    !       TO THE ROTATION FOR 2D. TECHNICALLY, ROTATION_RATE(3) SHOULD
    !       SPECIFY THE ROTATION FOR 2D AND YOU USE THE 3D CROSS PRODUCT
    !       (OMEGA x XYZ) WITH THE Z COORDINATE SET TO ZERO TO GET THE
    !       X AND Y ANGULAR VELOCITIES, BUT WE DONT WANT TO BE USING A
    !       'Z' VALUE FOR 2D.
    real(wp) :: rotation_rate(1:3) = [zero,zero,zero]
    !
    logical :: is_specified_in_mrf = fals
    !
#ifdef ALLOCATABLE_MASS_FRAC
    real(wp), allocatable :: mass_frac(:)
    character(len=10), allocatable :: species_name(:)
#else
    real(wp) :: mass_frac(7) = zero
    character(len=10) :: species_name(7) = ""
#endif
    !
    type(relaxation_t) :: relax
    !
    type(sem_bc_t) :: sem
    !
    type(periodic_bc_t) :: periodic
    !
  contains
    !
    procedure, private, pass(a)    :: equals    => bc_input_t_equals
    procedure, private, pass(a)    :: not_equal => bc_input_t_not_equal
    procedure :: is_not_default => is_not_default_bc_input
    procedure :: has_default_flow_conditions
    generic :: operator(==)  => equals
    generic :: operator(/=)  => not_equal
    !
  end type bc_input_t
  !
  !
  ! ################################
  ! ###  Module Local Variables  ###
  ! ################################
  !
  ! INPUT BOUNDARY CONDITIONS
  type(bc_input_t), public, allocatable :: bc_input(:)
  !
  !
  ! #####################################
  ! ###  Module Function Definitions  ###
  ! #####################################
  !
contains
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!> Check if a bc_input derived type is not set to default values.
!###############################################################################
pure function is_not_default_bc_input(this) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(bc_input_t), intent(in) :: this
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
  !.. Local Allocatable Scalars ..
  type(bc_input_t), allocatable :: bc_input_default
  !
continue
  !
  ! Create a variable of type bc_input_t with default values
  !
  allocate ( bc_input_default )
  !
  ! Compare the input argument to the default values
  !
  return_value = ( this /= bc_input_default )
  !
end function is_not_default_bc_input
!
!###############################################################################
!> Check if the flow condition components of a bc_input derived type
!> have default values (within some tolerance).
!###############################################################################
pure function has_default_flow_conditions(this) result(return_value)
  !
  ! This function is slightly different than is_not_default_bc_input. This
  ! function is used to check only the flow condition components of the
  ! boundary condition against the default values. The function
  ! is_not_default_bc_input also checks bc_input%name and other non flow
  ! condition components against the default values.
  !
  !.. Passed Argument from Invoking Object ..
  class(bc_input_t), intent(in) :: this
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
  !.. Local Scalars ..
  integer :: n
  !
  !.. Local Allocatable Scalars ..
  type(bc_input_t), allocatable :: dflt
  !
continue
  !
  ! Create a variable of type bc_input_t with default values
  !
  allocate ( dflt )
  !
  ! Initialize the function result to false
  !
  return_value = fals
  !
  ! NOTE: Dont check this%name, this%set_bc_value, or this%periodic because
  !       those three components dont have anything to do with the flow
  !       conditions for this boundary condition
  !
  ! If any of the following conditionals are true, then the subcomponents
  ! of the input argument are not equal to the default values, so immediately
  ! return with a function result of false indicating they are not equal.
  !
  if (abs(this%t_static   - dflt%t_static  ) > eps6) return
  if (abs(this%p_static   - dflt%p_static  ) > eps6) return
  if (abs(this%rho_static - dflt%rho_static) > eps6) return
  if (abs(this%mach       - dflt%mach      ) > eps6) return
  if (abs(this%vx         - dflt%vx        ) > eps6) return
  if (abs(this%vy         - dflt%vy        ) > eps6) return
  if (abs(this%vz         - dflt%vz        ) > eps6) return
  if (abs(this%t_total    - dflt%t_total   ) > eps6) return
  if (abs(this%p_total    - dflt%p_total   ) > eps6) return
  if (abs(this%rho_total  - dflt%rho_total ) > eps6) return
  if (abs(this%alpha_aoa  - dflt%alpha_aoa ) > eps6) return
  if (abs(this%beta_aoa   - dflt%beta_aoa  ) > eps6) return
  if (abs(this%wall_temp  - dflt%wall_temp ) > eps6) return
  !
  do n = 1,size(this%rotation_rate)
    if (abs(this%rotation_rate(n) - dflt%rotation_rate(n)) > eps6) return
  end do
  !
  if (this%is_specified_in_mrf .neqv. dflt%is_specified_in_mrf) return
  !
  if ( this%relax%is_not_default() ) return
  !
  if ( this%sem%is_not_default() ) return
  !
  ! CHECK ON this%periodic%is_not_default() IS NOT DONE BECAUSE WHETHER OR
  ! NOT THE PERIODIC COMPONENT HAS DEFAULT CONDITIONS HAS NOTHING TO DO WITH
  ! WHETHER OR NOT THIS BC_INPUT ELEMENT HAS DEFAULT FLOW CONDITIONS WHICH
  ! IS WHAT THIS FUNCTION IS CHECKING.
  !
#ifdef ALLOCATABLE_MASS_FRAC
  !
  if (allocated(this%mass_frac   ) .neqv. allocated(dflt%mass_frac   )) return
  if (allocated(this%species_name) .neqv. allocated(dflt%species_name)) return
  !
  if (allocated(this%mass_frac)) then
    if (size(this%mass_frac) /= size(dflt%mass_frac)) return
    do n = 1,size(this%mass_frac)
      if (abs(this%mass_frac(n) - dflt%mass_frac(n)) > eps6) return
    end do
  end if
  !
  if (allocated(this%species_name)) then
    if (size(this%species_name) /= size(dflt%species_name)) return
    do n = 1,size(this%species_name)
      if (trim(adjustl(this%species_name(n))) /= &
          trim(adjustl(dflt%species_name(n)))) return
    end do
  end if
  !
#else
  !
  do n = 1,size(this%mass_frac)
    if (abs(this%mass_frac(n) - dflt%mass_frac(n)) > eps6) return
  end do
  !
  do n = 1,size(this%species_name)
    if (trim(adjustl(this%species_name(n))) /= &
        trim(adjustl(dflt%species_name(n)))) return
  end do
  !
#endif
  !
  return_value = true
  !
end function has_default_flow_conditions
!
!###############################################################################
!> Check if two bc_input derived type variables are equal.
!###############################################################################
pure function bc_input_t_equals(a,b) result(return_value)
  !
  !.. Formal Arguments ..
  class(bc_input_t), intent(in) :: a
  class(bc_input_t), intent(in) :: b
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  ! Initialize the function result to false
  !
  return_value = fals
  !
  ! If any of the following conditionals are true, then the subcomponents
  ! of the input arguments are not equal, so immediately return with a
  ! function result of false indicating they are not equal.
  !
  if (trim(adjustl(a%name)) /= trim(adjustl(b%name))) return
  !
  if (a%set_bc_value /= b%set_bc_value) return
  !
  if (abs(a%t_static   - b%t_static  ) > eps6) return
  if (abs(a%p_static   - b%p_static  ) > eps6) return
  if (abs(a%rho_static - b%rho_static) > eps6) return
  if (abs(a%mach       - b%mach      ) > eps6) return
  if (abs(a%vx         - b%vx        ) > eps6) return
  if (abs(a%vy         - b%vy        ) > eps6) return
  if (abs(a%vz         - b%vz        ) > eps6) return
  if (abs(a%t_total    - b%t_total   ) > eps6) return
  if (abs(a%p_total    - b%p_total   ) > eps6) return
  if (abs(a%rho_total  - b%rho_total ) > eps6) return
  if (abs(a%alpha_aoa  - b%alpha_aoa ) > eps6) return
  if (abs(a%beta_aoa   - b%beta_aoa  ) > eps6) return
  if (abs(a%wall_temp  - b%wall_temp ) > eps6) return
  !
  do n = 1,size(a%rotation_rate)
    if (abs(a%rotation_rate(n) - b%rotation_rate(n)) > eps6) return
  end do
  !
  if (a%is_specified_in_mrf .neqv. b%is_specified_in_mrf) return
  !
#ifdef ALLOCATABLE_MASS_FRAC
  !
  if (allocated(a%mass_frac   ) .neqv. allocated(b%mass_frac   )) return
  if (allocated(a%species_name) .neqv. allocated(b%species_name)) return
  !
  if (allocated(a%mass_frac)) then
    if (size(a%mass_frac) /= size(b%mass_frac)) return
    do n = 1,size(a%mass_frac)
      if (abs(a%mass_frac(n) - b%mass_frac(n)) > eps6) return
    end do
  end if
  !
  if (allocated(a%species_name)) then
    if (size(a%species_name) /= size(b%species_name)) return
    do n = 1,size(a%species_name)
      if (trim(adjustl(a%species_name(n))) /= &
          trim(adjustl(b%species_name(n)))) return
    end do
  end if
  !
#else
  !
  do n = 1,size(a%mass_frac)
    if (abs(a%mass_frac(n) - b%mass_frac(n)) > eps6) return
  end do
  !
  do n = 1,size(a%species_name)
    if (trim(adjustl(a%species_name(n))) /= &
        trim(adjustl(b%species_name(n)))) return
  end do
  !
#endif
  !
  if (a%relax    /= b%relax   ) return
  if (a%sem      /= b%sem     ) return
  if (a%periodic /= b%periodic) return
  !
  ! If we have reached this point, all subcomponents of the input arguments
  ! are equal, so return true to indicate this is the case.
  !
  return_value = true
  !
end function bc_input_t_equals
!
!###############################################################################
!> Check if two bc_input derived type variables are not equal.
!###############################################################################
pure function bc_input_t_not_equal(a,b) result(return_value)
  !
  !.. Formal Arguments ..
  class(bc_input_t), intent(in) :: a
  class(bc_input_t), intent(in) :: b
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
continue
  !
  return_value = .not. (a == b)
  !
end function bc_input_t_not_equal
!
!###############################################################################
!> Copy bc_input_t from ovar to local
!>
!> Converts the new bc_input structure from ovar
!> into the old version 1 namelist structure.
!###############################################################################
pure &
function bc_input_ovar_to_v0001(ovar_bc) result(return_value)
  !
  !.. Formal Arguments ..
  type(ovar_bc_input_t), allocatable, intent(in) :: ovar_bc(:)
  !
  !.. Function Result ..
  type(bc_input_t), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: n,nbc
  !
continue
  !
  nbc = 0
  if (allocated(ovar_bc)) nbc = size(ovar_bc)
  !
  allocate ( return_value(1:nbc) )
  !
  do n = 1,nbc
    !
    return_value(n)%name                = ovar_bc(n)%name
    !
    return_value(n)%set_bc_value        = ovar_bc(n)%set_bc_value
    !
    return_value(n)%t_static            = ovar_bc(n)%temperature
    return_value(n)%p_static            = ovar_bc(n)%pressure
    return_value(n)%rho_static          = ovar_bc(n)%density
    return_value(n)%mach                = ovar_bc(n)%mach
    return_value(n)%vx                  = ovar_bc(n)%velocityx
    return_value(n)%vy                  = ovar_bc(n)%velocityy
    return_value(n)%vz                  = ovar_bc(n)%velocityz
    return_value(n)%t_total             = ovar_bc(n)%total_temperature
    return_value(n)%p_total             = ovar_bc(n)%total_pressure
    return_value(n)%rho_total           = ovar_bc(n)%total_density
    return_value(n)%alpha_aoa           = ovar_bc(n)%aoa_alpha
    return_value(n)%beta_aoa            = ovar_bc(n)%aoa_beta
    return_value(n)%wall_temp           = ovar_bc(n)%wall_temp
    !
    return_value(n)%rotation_rate       = ovar_bc(n)%rotation_rate
    !
    return_value(n)%is_specified_in_mrf = ovar_bc(n)%is_specified_in_mrf
    !
    return_value(n)%mass_frac           = ovar_bc(n)%mass_frac
    return_value(n)%species_name        = ovar_bc(n)%species_name
    !
    return_value(n)%relax               = ovar_bc(n)%relax
    !
    return_value(n)%sem                 = ovar_bc(n)%sem
    !
    return_value(n)%periodic            = ovar_bc(n)%periodic
    !
  end do
  !
end function bc_input_ovar_to_v0001
!
!###############################################################################
!> Copy bc_input_t from local to ovar.
!>
!> Converts the old version 1 namelist structure
!> into the new bc_input structure from ovar.
!###############################################################################
pure &
function bc_input_v0001_to_ovar(bc_v1) result(return_value)
  !
  !.. Formal Arguments ..
  type(bc_input_t), allocatable, intent(in) :: bc_v1(:)
  !
  !.. Function Result ..
  type(ovar_bc_input_t), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: n,nbc
  !
continue
  !
  nbc = 0
  if (allocated(bc_v1)) nbc = size(bc_v1)
  !
  allocate ( return_value(1:nbc) )
  !
  do n = 1,nbc
    !
    return_value(n)%name                = bc_v1(n)%name
    !
    return_value(n)%set_bc_value        = bc_v1(n)%set_bc_value
    !
    return_value(n)%temperature         = bc_v1(n)%t_static
    return_value(n)%pressure            = bc_v1(n)%p_static
    return_value(n)%density             = bc_v1(n)%rho_static
    return_value(n)%mach                = bc_v1(n)%mach
    return_value(n)%velocityx           = bc_v1(n)%vx
    return_value(n)%velocityy           = bc_v1(n)%vy
    return_value(n)%velocityz           = bc_v1(n)%vz
    return_value(n)%total_temperature   = bc_v1(n)%t_total
    return_value(n)%total_pressure      = bc_v1(n)%p_total
    return_value(n)%total_density       = bc_v1(n)%rho_total
    return_value(n)%aoa_alpha           = bc_v1(n)%alpha_aoa
    return_value(n)%aoa_beta            = bc_v1(n)%beta_aoa
    return_value(n)%wall_temp           = bc_v1(n)%wall_temp
    !
    return_value(n)%rotation_rate       = bc_v1(n)%rotation_rate
    !
    return_value(n)%is_specified_in_mrf = bc_v1(n)%is_specified_in_mrf
    !
    return_value(n)%mass_frac           = bc_v1(n)%mass_frac
    return_value(n)%species_name        = bc_v1(n)%species_name
    !
    return_value(n)%relax               = bc_v1(n)%relax
    !
    return_value(n)%sem                 = bc_v1(n)%sem
    !
    return_value(n)%periodic            = bc_v1(n)%periodic
    !
  end do
  !
end function bc_input_v0001_to_ovar
!
!###############################################################################
!
end module bc_input_0001_mod
!
!###############################################################################
!
