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
  private :: bc_input_t_copyin
  private :: bc_input_t_copyout
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
    procedure,          pass(this) :: copyout => bc_input_t_copyout
    procedure, private, pass(this) :: copyin  => bc_input_t_copyin
    procedure, private, pass(a)    :: equals    => bc_input_t_equals
    procedure, private, pass(a)    :: not_equal => bc_input_t_not_equal
    procedure :: is_not_default => is_not_default_bc_input
    procedure :: has_default_flow_conditions
    generic :: assignment(=) => copyin
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
!> This has been made a generic type-bound assignment which can be used as:
!>    bc_input = ovar_bc_input(:)
!###############################################################################
elemental &
subroutine bc_input_t_copyin(this,from)
  !
  !.. Formal Arguments ..
  class(     bc_input_t), intent(inout) :: this
  type(ovar_bc_input_t),  intent(in   ) :: from
  !
continue
  !
  this%name                = from%name
  !
  this%set_bc_value        = from%set_bc_value
  !
  this%t_static            = from%temperature
  this%p_static            = from%pressure
  this%rho_static          = from%density
  this%mach                = from%mach
  this%vx                  = from%velocityx
  this%vy                  = from%velocityy
  this%vz                  = from%velocityz
  this%t_total             = from%total_temperature
  this%p_total             = from%total_pressure
  this%rho_total           = from%total_density
  this%alpha_aoa           = from%aoa_alpha
  this%beta_aoa            = from%aoa_beta
  this%wall_temp           = from%wall_temp
  !
  this%rotation_rate       = from%rotation_rate
  !
  this%is_specified_in_mrf = from%is_specified_in_mrf
  !
  this%mass_frac           = from%mass_frac
  this%species_name        = from%species_name
  !
  this%relax               = from%relax
  !
  this%sem                 = from%sem
  !
  this%periodic            = from%periodic
  !
end subroutine bc_input_t_copyin
!
!###############################################################################
!> Copy bc_input_t from local to ovar.
!>
!> Converts the old version 1 namelist structure
!> into the new bc_input structure from ovar.
!> This cant be made a generic type-bound assignment like copyin because
!> an assignment must be made to the bound derived type. Therefore, this
!> has been bound as a simple function that can be used as:
!>    ovar_bc_input = bc_input(:)%copyout
!###############################################################################
elemental &
function bc_input_t_copyout(this) result(return_value)
  !
  !.. Formal Arguments ..
  class(bc_input_t), intent(in) :: this
  !
  !.. Function Result ..
  type(ovar_bc_input_t) :: return_value
  !
continue
  !
  return_value%name                = this%name
  !
  return_value%set_bc_value        = this%set_bc_value
  !
  return_value%temperature         = this%t_static
  return_value%pressure            = this%p_static
  return_value%density             = this%rho_static
  return_value%mach                = this%mach
  return_value%velocityx           = this%vx
  return_value%velocityy           = this%vy
  return_value%velocityz           = this%vz
  return_value%total_temperature   = this%t_total
  return_value%total_pressure      = this%p_total
  return_value%total_density       = this%rho_total
  return_value%aoa_alpha           = this%alpha_aoa
  return_value%aoa_beta            = this%beta_aoa
  return_value%wall_temp           = this%wall_temp
  !
  return_value%rotation_rate       = this%rotation_rate
  !
  return_value%is_specified_in_mrf = this%is_specified_in_mrf
  !
  return_value%mass_frac           = this%mass_frac
  return_value%species_name        = this%species_name
  !
  return_value%relax               = this%relax
  !
  return_value%sem                 = this%sem
  !
  return_value%periodic            = this%periodic
  !
end function bc_input_t_copyout
!
!###############################################################################
!
end module bc_input_0001_mod
!
!###############################################################################
!
