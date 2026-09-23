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
!> Convert all blocks of namelist input file to latest format.
!> NOTE: This should only be called by the global root
!>        MPI process if running in parallel.
!###############################################################################
subroutine convert_nml_file(fname,nml_ver_prev,mode)
  !
  !.. Use Statements ..
  use module_kind_types, only : debug_timer,start_timer
  use module_kind_types, only : io_error,error_message
  use module_kind_types, only : lk,fals,true
  !
  use nml_base_mod, only : new_nml_file1
  use nml_base_mod, only : new_nml_file2
  use nml_base_mod, only : new_nml_file3
  use nml_base_mod, only : new_nml_file4
  use nml_base_mod, only : nml_recl
  use nml_base_mod, only : merge_nml_files
  use nml_base_mod, only : namelist_procedures_t
  use nml_base_mod, only : set_NMLaction
  use nml_base_mod, only : NMLaction_SetDummy0,NMLaction_SetDummy1
  use nml_base_mod, only : NMLaction_SetValue,NMLaction_WritePretty
  use nml_base_mod, only : set_NML_mode
  use nml_base_mod, only : NMLmode_InputConverter
  use nml_base_mod, only : NMLmode_TrimDefaults
  !
  !.. Formal Arguments ..
  !> Name of the input file to convert.
  character(len=*),             intent(in) :: fname
  !> Namelist procedures for processing the input file
  !>  Note: The namelist version should have already been determined before
  !>        calling this subroutine since this argument needs to contain the
  !>        correct procedures to process the namelist(s) found in fname.
  class(namelist_procedures_t), intent(in) :: nml_ver_prev
  !> Conversion mode:
  !>   'short'  Produce a minimal file with just the user's inputs
  !>            and a few other things. Default.
  !>   'full'   Produce the full input file including every option.
  character(len=*), optional,   intent(in) :: mode
  !
  !.. Local Scalars ..
  integer :: action, io, io1, ierr, n
  integer :: int_mode
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: new_nml_file
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "convert_nml_file"
  !
  logical(lk), parameter :: debug_merge_files = fals
 !logical(lk), parameter :: debug_merge_files = true
  !
continue
  !
  ! NOTE: This should only be called by the global root
  !       MPI process if running in parallel.
  !
  ! Set the NML mode so the input variables are processed with with minimal
  ! modifications, and the only modifications that are kept when writing out the
  ! out the converted input file are those that are needed to correct for
  ! invalid input.
  !
  call set_NML_mode(NMLmode_InputConverter)
  !
  ! Open the input file
  !
  call debug_timer(start_timer,"actually reading the input file")
  open (newunit=io,file=trim(fname),status="old", &
                   action="read",form="formatted", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,trim(fname),1,__LINE__,__FILE__,ierr,error_message)
  !
  int_mode = 0  ! set default mode
  if ( present(mode) ) then
    if ( mode == 'short' ) int_mode = 0  ! short mode
    if ( mode == 'full'  ) int_mode = 1  ! full  mode
  end if
  !
  if ( int_mode == 0 ) then  ! short mode
    !
    ! Process the input file twice using different dummy default values.
    ! If we get the same input value, then user explicitly specified it.
    !
    do n = 1, 2
      ! Set namelist output file and action.
      if ( n == 1 ) then
        new_nml_file = trim(new_nml_file1)  ! F2003 realloc
        action = NMLaction_SetDummy0
      else
        new_nml_file = trim(new_nml_file2)  ! F2003 realloc
        action = NMLaction_SetDummy1
      end if
      call set_NMLaction(action)
      !
      ! Allocate namelist variables to default sizes.
      call nml_ver_9999%alloc_data()  ! alloc new arrays
      call nml_ver_prev%alloc_data()  ! alloc old arrays
      !
      ! Set namelist values.
      call nml_ver_9999%set_data(action)  ! set new arrays
      call nml_ver_prev%set_data(action)  ! set old arrays
      !
      ! Read namelist values.
      call nml_ver_prev%read_data(io)  ! read old vars
      !
      ! Resize namelist allocatables to their actual content.
      call nml_ver_9999%resiz_data()  ! resize new arrays
      !
      ! Write namelist values in new format.
     !call nml_ver_prev%write_file("old."//new_nml_file,NMLaction_WritePretty)
      call nml_ver_9999%write_file(new_nml_file,NMLaction_WritePretty)
     !open (newunit=io1,file=new_nml_file,status='replace', &
     !                  action='write',form='formatted', &
     !                  position="rewind",recl=nml_recl, &
     !                  iostat=ierr,iomsg=error_message)
     !call io_error(pname,new_nml_file,1,__LINE__,__FILE__, &
     !              ierr,error_message)
     !! write new vars
     !call nml_ver_9999%write_data(NMLaction_WritePretty,io1,banner=fals)
     !close (unit=io1,iostat=ierr,iomsg=error_message)
     !call io_error(pname,new_nml_file,2,__LINE__,__FILE__, &
     !              ierr,error_message)
    end do ! n
    !
    ! Merge the shared values of the two new namelist files into one file.
    ! This file will have the minimal content that the user specified.
    !
    call merge_nml_files(save_temp_files=debug_merge_files)
    !
    ! Finally, set the NML mode to trim the defaults out of the merged file.
    !
    call set_NML_mode(NMLmode_TrimDefaults)
    !
    ! Set the NML Action to set the default values
    !
    call set_NMLaction(NMLaction_SetValue)
    !
    ! Read the merged namelist file
    !
    call nml_ver_9999%read_file(new_nml_file3)
    !
    ! Write out the merged namelist file but with NML mode set to TrimDefaults,
    ! this will trim out any variables set to their default values 
    !
    call nml_ver_9999%write_file(new_nml_file4,NMLaction_WritePretty)
   !open (newunit=io1,file=new_nml_file4,status='replace', &
   !                  action='write',form='formatted', &
   !                  position="rewind",recl=nml_recl, &
   !                  iostat=ierr,iomsg=error_message)
   !call io_error(pname,new_nml_file4,1,__LINE__,__FILE__, &
   !              ierr,error_message)
   !! write new vars
   !call nml_ver_9999%write_data(NMLaction_WritePretty,io1,banner=fals)
   !close (unit=io1,iostat=ierr,iomsg=error_message)
   !call io_error(pname,new_nml_file4,2,__LINE__,__FILE__, &
   !              ierr,error_message)
    !
  else if ( int_mode == 1 ) then  ! full mode
    !
    ! Process the input file once using actual default values.
    !
    action = NMLaction_SetValue
    call set_NMLaction(action)
    !
    ! Allocate namelist variables to default sizes.
    call nml_ver_9999%alloc_data()  ! alloc new arrays
    call nml_ver_prev%alloc_data()  ! alloc old arrays
    !
    ! Set namelist values.
    call nml_ver_9999%set_data(action)  ! set new arrays
    call nml_ver_prev%set_data(action)  ! set old arrays
    !
    ! Read namelist values.
    call nml_ver_prev%read_data(io)  ! read old vars
    !
    ! Resize namelist allocatables to their actual content.
    call nml_ver_9999%resiz_data()  ! resize new arrays
    !
    ! Write namelist values in new format.
    call nml_ver_9999%write_file(new_nml_file3,NMLaction_WritePretty)
   !new_nml_file = new_nml_file3
   !open (newunit=io1,file=new_nml_file,status='unknown', &
   !                  action='write',form='formatted', &
   !                  position="rewind",recl=nml_recl, &
   !                  iostat=ierr,iomsg=error_message)
   !call io_error(pname,new_nml_file,1,__LINE__,__FILE__, &
   !              ierr,error_message)
   !! write new vars
   !call nml_ver_9999%write_data(NMLaction_WritePretty,io1,banner=fals)
   !close (unit=io1,iostat=ierr,iomsg=error_message)
   !call io_error(pname,new_nml_file,2,__LINE__,__FILE__, &
   !              ierr,error_message)
    !
  end if
  !
  close (unit=io,iostat=ierr,iomsg=error_message)
  call io_error(pname,trim(fname),2,__LINE__,__FILE__,ierr,error_message)
  !
end subroutine convert_nml_file
!
!###############################################################################
!
end module nml_mod
!
!###############################################################################
!
