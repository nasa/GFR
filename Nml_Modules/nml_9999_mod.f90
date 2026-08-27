!###############################################################################
!> Module for version 9999 of namelist io.
!> This is an alias for the most recent namelist module.
!###############################################################################
module nml_9999_mod
  !
  !.. Use Statements ..
  use nml_base_mod, only : namelist_procedures_t
  !
  implicit none
  !
  private
  !
  public :: set_latest_nml_version
  !
  !> Extension of the abstract namelist_procedures_t to house the
  !> procedures for processing the latest version of the namelist.
  type, public, extends(namelist_procedures_t) :: nml_9999_procedures_t
  contains
    !> This is not needed for the 9999/latest version,
    !> so it is set to an empty subroutine.
    !> NOTE: Trying set it to a subroutine that does something similar to
    !>       initialize_namelist_versions would result in a circular
    !>       build dependency.
    procedure, nopass :: init => dummy_init
  end type nml_9999_procedures_t
  !
  !> Set of procedure pointers that are assigned to the respective set
  !> of procedures for processing the latest version of the namelist.
  type(nml_9999_procedures_t), public, protected, save :: nml_ver_9999
  !
contains
!
!###############################################################################
!> Subroutine to set the nml_ver_9999 procedure pointers, but only if the
!> namelist version passed in to the argument nml_ver matches the parameter
!> latest_namelist_version from the module nml_base_mod.
!###############################################################################
subroutine set_latest_nml_version(nml_ver,nml_procedures)
  !
  !.. Use Statements ..
  use nml_base_mod, only : latest_namelist_version
  !
  !.. Formal Arguments ..
  integer,                      intent(in) :: nml_ver
  class(namelist_procedures_t), intent(in) :: nml_procedures
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "set_latest_nml_version"
  !
continue
  !
  ! This assumes that each version of the namelist has initialized its own
  ! namelist_procedures_t variable to assign all of these procedure pointers
  ! to the respective subroutines for that namelist version.
  ! NOTE: This assumption is easily made valid if this subroutine is called
  !       immediately after a namelist version has initialized its own
  !       procedure pointers. For examples, check out these subroutines:
  !           init_nml0001_procedures     in   nml_0001_mod
  !           init_nml0002_procedures     in   nml_0002_mod
  !
  if (nml_ver == latest_namelist_version) then
    !
    nml_ver_9999%alloc_data => nml_procedures%alloc_data
    nml_ver_9999%resiz_data => nml_procedures%resiz_data
    nml_ver_9999%  set_data => nml_procedures%  set_data
    nml_ver_9999% read_data => nml_procedures% read_data
    nml_ver_9999%write_data => nml_procedures%write_data
    nml_ver_9999% read_file => nml_procedures% read_file
    nml_ver_9999%write_file => nml_procedures%write_file
    !
  end if
  !
end subroutine set_latest_nml_version
!
!###############################################################################
!
subroutine dummy_init()
  !
  ! Blank subroutine to satisfy the deferred procedure init when
  ! extending derived type namelist_procedures_t
  !
continue
  !
end subroutine dummy_init
!
!###############################################################################
!
end module nml_9999_mod
