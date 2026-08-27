!###############################################################################
!> Module to house all the different namelist versions for easy access
!> elsewhere without having to put use statements for each version.
!###############################################################################
module nml_mod
  !
  !.. Use Statements ..
  use nml_0001_mod, only : nml_ver_0001
  use nml_0002_mod, only : nml_ver_0002
  !> Do not remove nml_ver_9999! This is included here so that this is
  !> the only module that needs to be used to access any version of the
  !> namelist including this 'alias' for accessing the latest version.
  use nml_9999_mod, only : nml_ver_9999 ! Do not remove!
  !
  implicit none
  !
  public
  !
contains
!
!###############################################################################
!> A single subroutine to initialize all namelist_procedures_t procedure
!> pointers for each namelist version.
!>
!> Any time a new namelist version is implemented, a call to the type-bound
!> init subroutine for its nml_ver_#### type (i.e., call nml_ver_####%init)
!> should be added to this subroutine to initialize the procedure pointers
!> for the new namelist version.
!###############################################################################
subroutine initialize_namelist_versions()
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_namelist_versions"
  !
continue
  !
  ! Initialize the procedure pointers for namelist version 0001
  call nml_ver_0001%init
  !
  ! Initialize the procedure pointers for namelist version 0002
  call nml_ver_0002%init
  !
  ! We do not need to initialize the procedure pointers for namelist
  ! version 9999 (latest version) because the above initializations will
  ! each call subroutine set_latest_nml_version which will initialize
  ! the procedure pointers in nml_ver_9999 if that namelist version
  ! matches the value of latest_namelist_version in nml_base_mod.
  !
end subroutine initialize_namelist_versions
!
!###############################################################################
!
end module nml_mod
!
!###############################################################################
!
