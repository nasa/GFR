!###############################################################################
!> Module to house all the different namelist versions for easy access
!> elsewhere without having to put use statements for each version.
!###############################################################################
module nml_mod
  !
  !.. Use Statements ..
  use nml_0001_mod, only : nml_ver_0001
  use nml_0002_mod, only : nml_ver_0002
  use nml_9999_mod, only : nml_ver_9999
  !
  implicit none
  !
  public
  !
contains
!
!###############################################################################
!
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
  ! each call subroutine set_current_nml_version which will initialize
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
