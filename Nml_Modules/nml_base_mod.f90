!###############################################################################
!> Base module for all other namlist input modules.
!###############################################################################
module nml_base_mod
  !
  !.. Use Statements ..
  use module_kind_types
  !
  implicit none
  !
  private
  !
  ! ###############################
  ! #####  MODULE PROCEDURES  #####
  ! ###############################
  !
  public  :: get_namelist_version_from_file  ! fun(fname)
  public  :: get_namelist_version_from_io    ! fun(io)
  public  :: io_has_nml_block      ! sub(io,block_name)
  public  :: set_ior               ! sub(io)
  public  :: set_iow               ! sub(io)
  public  :: set_NMLaction         ! sub(action)
  public  :: set_NML_textcase      ! sub(mode)
  public  :: set_value             ! sub(nam,var,val)    generic interface
  public  :: iow_value             ! sub(nam,var)        generic interface
  public  :: iow_header            ! sub(text,lvar)
  public  :: iow_footer            ! sub(lvar)
  public  :: merge_nml_files       ! sub()
  !
  private :: set_value_i           ! sub(nam,var,val,act)
  private :: set_value_r           ! sub(nam,var,val,act)
  private :: set_value_c           ! sub(nam,var,val)      char(*)
  private :: set_value_s           ! sub(nam,var,val,len)  char(:)
  private :: set_value_l           ! sub(nam,var,val,act)
  private :: set_alloc_char        ! sub(nam,var,val,len)
  !
  private :: iow_value_i           ! sub(nam,var)
  private :: iow_value_r           ! sub(nam,var)
  private :: iow_value_c           ! sub(nam,var)
  private :: iow_value_l           ! sub(nam,var)
  private :: iow_value_lk          ! sub(nam,var)
  private :: iow_value_iarr        ! sub(nam,var)
  private :: iow_value_rarr        ! sub(nam,var)
  private :: iow_value_carr        ! sub(nam,var)
  !
  ! ###############################
  ! #####  MODULE PARAMETERS  #####
  ! ###############################
  !
  ! Latest namelist format.
  integer, public, parameter :: latest_namelist_version = 0001
  !
  ! Record length for temporary namelist files during old to new conversion.
  integer, public, parameter :: nml_recl = 5000
  !
  ! File names for temporary namelist files during old to new conversion.
  character(len=18), public, parameter :: new_nml_file1 = 'new_nml_dummy0.nml'
  character(len=18), public, parameter :: new_nml_file2 = 'new_nml_dummy1.nml'
  character(len=17), public, parameter :: new_nml_file3 = 'new_nml_input.nml'
  !
  ! Dummy default values for 2-pass namelist reads.
  integer,          public, parameter :: idum0 =  -9990  , idum1 =  -9991
  real,             public, parameter :: rdum0 =  -9990. , rdum1 =  -9991.
  character(len=5), public, parameter :: cdum0 = '-9990' , cdum1 = '-9991'
  logical,          public, parameter :: ldum0 =   fals  , ldum1 =   true
  !
  ! ###############################
  ! #####   MODULE VARIABLES  #####
  ! ###############################
  !
  ! Permissible actions for NML_action.
  integer, public, parameter :: NMLaction_SetDummy0   = 0 !< Set to dummy0.
  integer, public, parameter :: NMLaction_SetDummy1   = 1 !< Set to dummy1.
  integer, public, parameter :: NMLaction_SetValue    = 2 !< Set to value.
  integer, public, parameter :: NMLaction_Read        = 3 !< Standard nml read.
  integer, public, parameter :: NMLaction_Write       = 4 !< Standard nml write.
  integer, public, parameter :: NMLaction_WritePretty = 5 !< Pretty nml write.
  !> Controls the function of the set_value* subroutines.
  integer, public, protected :: NMLaction = NMLaction_SetValue
  !
  !> Minimum column width for namelist varables.
  integer, private :: NML_varname_width = 40
  !
  !> Text case for namelist blocks and variables using NMLaction_WritePretty:
  !> (0) titlecase, (1) lowercase, (2) uppercase.
  integer, private :: NML_textcase  = 0 
  !
  !> Write output iow when a variable is set to a value.
  logical(lk), public :: setverbose = fals
  !
  !> Unit for namelist read during NMLaction_Read.
  integer, public, protected :: ior = 0
  !
  !> Unit for namelsit write during NMLaction_Write or NMLaction_WritePretty.
  integer, public, protected :: iow = 0
  !
  !> Namelist block header (ie, '&input').
  character(len=:), allocatable, public, protected :: header
  !
  !> Namelist buffer.
  character(len=:), allocatable, public, protected :: buffer
  !
  ! ####################################
  ! #####   MODULE DERIVED TYPES   #####
  ! ####################################
  !
  !> Derived type that contains procedure pointers for the base set of
  !> procedures needed to process a version of the input namelist(s).
  !> It also contains the subroutine set_nml_procedures which assigns
  !> the pointers to the respective procedures for a namelist version
  !> which are sent in as arguments.
  type, abstract, public :: namelist_procedures_t
    integer :: version = 0 ! namelist version # for these procedures
    procedure(alloc_namelist_data), pointer, nopass :: alloc_data => null()
    procedure(resiz_namelist_data), pointer, nopass :: resiz_data => null()
    procedure(  set_namelist_data), pointer, nopass ::   set_data => null()
    procedure( read_namelist_data), pointer, nopass ::  read_data => null()
    procedure(write_namelist_data), pointer, nopass :: write_data => null()
    procedure( read_namelist_file), pointer, nopass ::  read_file => null()
    procedure(write_namelist_file), pointer, nopass :: write_file => null()
  contains
    procedure :: set_nml_procedures
    procedure(blank_sub), deferred, nopass :: init
  end type namelist_procedures_t
  !
  ! #########################################
  ! #####   MODULE GENERIC INTERFACES   #####
  ! #########################################
  !
  !> Generic interface to set a variable.
  interface set_value                   ! sub(nam, var, val)
    module procedure set_value_i,    &  ! sub(nam,ivar,ival)     integer
                     set_value_r,    &  ! sub(nam,rvar,rval)     real
                     set_value_c,    &  ! sub(nam,cvar,cval)     char
                     set_value_s,    &  ! sub(nam,cvar,cval,len) char(:)
                     set_value_l,    &  ! sub(nam,lvar,lval)     logical
                     set_value_lk,   &  ! sub(lnam,var,lval)     logical kind
                     set_value_iarr, &  ! sub(nam,ivar,ival)     integer array
                     set_value_rarr     ! sub(nam,rvar,rval)     real    array
  end interface set_value
  !
  !> Generic interface to Write a variable to unit iow.
  interface iow_value                   ! sub(nam, var)
    module procedure iow_value_i   , &  ! sub(nam,ivar) integer
                     iow_value_r   , &  ! sub(nam,rvar) real
                     iow_value_c   , &  ! sub(nam,cvar) char
                     iow_value_l   , &  ! sub(nam,lvar) logical
                     iow_value_lk  , &  ! sub(nam,lvar) logical kind
                     iow_value_iarr, &  ! sub(nam,lvar) integer array
                     iow_value_rarr, &  ! sub(nam,lvar) real    array
                     iow_value_carr     ! sub(nam,lvar) char    array
  end interface iow_value
  !
  !> Abstract interfaces for the procedure pointer components of the
  !> derived type namelist_procedures_t.
  abstract interface
    !
    !> Abstract interface for a subroutine to allocate all namelist
    !> variables and arrays before actually trying to read the namelist.
    subroutine alloc_namelist_data(opt_n_bc_input,opt_char_len, &
                                   opt_n_species,opt_n_point_probes, &
                                   opt_n_rake_probes,opt_n_cell_probes)
      !.. Optional Arguments ..
      integer, optional, intent(in) :: opt_n_bc_input
      integer, optional, intent(in) :: opt_char_len
      integer, optional, intent(in) :: opt_n_species
      integer, optional, intent(in) :: opt_n_point_probes
      integer, optional, intent(in) :: opt_n_rake_probes
      integer, optional, intent(in) :: opt_n_cell_probes
    end subroutine alloc_namelist_data
    !
    !> Abstract interface for a subroutine to resize all allocatable namelist
    !> variables and arrays to the minimum size needed based on the information
    !> read in from the namelist.
    subroutine resiz_namelist_data()
    end subroutine resiz_namelist_data
    !
    !> Subroutine to set the namelist variables based on the input action
    subroutine set_namelist_data(action,io)
      !.. Formal Arguments ..
      integer,           intent(in) :: action  !< An NMLaction.
      integer, optional, intent(in) :: io      !< Unit number to echo values.
    end subroutine set_namelist_data
    !
    !> Subroutine to read in the namelist data from an open file
    subroutine read_namelist_data(io)
      !.. Formal Arguments ..
      integer, intent(in) :: io      !< Unit number for input.
    end subroutine read_namelist_data
    !
    !> Subroutine to write the namelist data to an open file
    subroutine write_namelist_data(action,io,banner)
      import :: lk
      !.. Formal Arguments ..
      integer,     intent(in) :: action  !< An NMLaction.
      integer,     intent(in) :: io      !< Unit number for output.
      logical(lk), intent(in) :: banner  !< Flag to write banner.
    end subroutine write_namelist_data
    !
    !> Subroutine to perform all steps of reading data from a namelist file:
    !>   1. Allocate namelist variables,
    !>   2. initialize to default values,
    !>   3. open the namelist file,
    !>   4. read the namelist data from the file including any conversions if
    !>      the version of the namelist file is older than the latest version,
    !>   5. close the namelist file, and
    !>   6. resize all allocatable namelist variables.
    subroutine read_namelist_file(input_fname)
      !.. Formal Arguments ..
      character(len=*), intent(in) :: input_fname
    end subroutine read_namelist_file
    !
    !> Subroutine to perform all steps of writing the namelist data to a file:
    !>   1. Allocate namelist variables,
    !>   2. initialize to default values,
    !>   3. open the namelist file,
    !>   4. write the namelist data from the file including any conversions if
    !>      the version of the namelist file is older than the latest version,
    !>   5. close the namelist file.
    subroutine write_namelist_file(output_fname,action)
      !.. Formal Arguments ..
      character(len=*), intent(in) :: output_fname
      integer,          intent(in) :: action        !< An NMLaction.
    end subroutine write_namelist_file
    !
    subroutine blank_sub()
    end subroutine blank_sub
    !
  end interface
  !
contains
!
!###############################################################################
!> Get the namelist version number from named file.
!###############################################################################
function get_namelist_version_from_file(fname) result(namelist_version)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: fname  !< name of input file
  integer                      :: namelist_version
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_namelist_version_from_file"
  !
  !.. Local Scalars ..
  integer :: io, ierr
  !
continue
  !
  ! Open the input file
  !
  call debug_timer(start_timer,"actually reading the input file")
  open(newunit=io,file=trim(fname),status="old",action="read", &
                     form="formatted",iostat=ierr,iomsg=error_message)
  call io_error(pname,trim(fname),1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Get the namelist version number from the file.
  !
  namelist_version = get_namelist_version_from_io(io)
  !
  ! Close the input file.
  !
  close(unit=io,iostat=ierr,iomsg=error_message)
  call io_error(pname,trim(fname),2,__LINE__,__FILE__,ierr,error_message)
  !
end function get_namelist_version_from_file
  !
!###############################################################################
!> Get the namelist version number from file unit.
!###############################################################################
function get_namelist_version_from_io(io) result(namelist_version)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: io  !< unit number of input file
  integer             :: namelist_version
  !
  !.. Local Variables ..
  !
  namelist / nml_version / namelist_version
  !
continue
  !
  namelist_version = latest_namelist_version  ! default
  !
  if ( io_has_nml_block(io,'&input') ) then
    namelist_version = 1
  else if ( io_has_nml_block(io,'&nml_version') ) then
    read(io,nml=nml_version)
  end if
  !
end function get_namelist_version_from_io
!
!###############################################################################
!> Test if unit io contains the requested namelist input block
!> and rewind the file for subsequent use.
!###############################################################################
function io_has_nml_block(io,block_name) result(return_value)
  !
  !.. Formal Arguments ..
  integer,          intent(in) :: io          !< unit number of input file
  character(len=*), intent(in) :: block_name  !< name of namelist block,
                                              !< (ie, '&input')
  logical :: return_value
  !
  !.. Local Scalars ..
  logical(lk) :: match_some_block, found_some_block
  logical(lk) :: match_this_block, found_this_block
  character(len=1000) text
  !
continue
  !
  return_value = fals  ! namelist block not found
  !
  ! If no specific block name given, match some block instead.
  ! Otherwise, match this block.
  match_some_block = ( block_name == '' .or. block_name == '&' )
  match_this_block = .not. match_some_block
  !
  rewind(io)
  !
  file_read_loop: &
  do
    read(io,'(a)',end=100) text
    text = uppercase(trim(adjustl(text)))
    found_some_block = ( text(1:1) == '&' )
    found_this_block = ( text == uppercase(trim(adjustl(block_name))) )
    if ( (match_some_block .and. found_some_block) .or. &
         (match_this_block .and. found_this_block) ) then
      return_value = true  ! namelist block found
      exit file_read_loop
    end if
  end do file_read_loop
  !
100 continue
  rewind(io)
  !
end function io_has_nml_block
!
!###############################################################################
!> Subroutine that assigns the procedure pointers of a namelist_procedures_t
!> entity to the respective procedures passed in as arguments.
!###############################################################################
subroutine set_nml_procedures(this, &
                              alloc_data_proc,resiz_data_proc,set_data_proc, &
                               read_data_proc,write_data_proc, &
                               read_file_proc,write_file_proc)
  !
  !.. Formal Arguments ..
  class(namelist_procedures_t), intent(inout) :: this
  !
  procedure(alloc_namelist_data) :: alloc_data_proc
  procedure(resiz_namelist_data) :: resiz_data_proc
  procedure(  set_namelist_data) ::   set_data_proc
  procedure( read_namelist_data) ::  read_data_proc
  procedure(write_namelist_data) :: write_data_proc
  procedure( read_namelist_file) ::  read_file_proc
  procedure(write_namelist_file) :: write_file_proc
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "set_nml_procedures"
  !
continue
  !
  this%alloc_data => alloc_data_proc
  this%resiz_data => resiz_data_proc
  this%  set_data =>   set_data_proc
  this% read_data =>  read_data_proc
  this%write_data => write_data_proc
  this% read_file =>  read_file_proc
  this%write_file => write_file_proc
  !
end subroutine set_nml_procedures
!
!###############################################################################
!> Sets the module variable NMLaction.
!###############################################################################
subroutine set_NMLaction(val)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "set_NMLaction"
  !
continue
  !
  if ( any( val == [ NMLaction_SetDummy0, NMLaction_SetDummy1, &
                     NMLaction_SetValue,  NMLaction_Read,      &
                     NMLaction_Write,     NMLaction_WritePretty ] ) ) then
    NMLaction = val
  else
    if (mypnum == glb_root) write(iout,101) val
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  end if
  !
  101 format(/,"CODING ERROR: Invalid argument, val = ",i0)
end subroutine set_NMLaction
!
!###############################################################################
!> Sets the io for namelist read. Assumes the unit is already open.
!###############################################################################
subroutine set_ior(io)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: io
  !
continue
  !
  ior = io
  !
end subroutine set_ior
!
!###############################################################################
!> Sets the io for namelist write. Assumes the unit is already open.
!###############################################################################
subroutine set_iow(io)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: io
  !
continue
  !
  iow = io
  header = ''
  buffer = ''
  !
end subroutine set_iow
!
!###############################################################################
!> Set the module variable NML_textcase.
!###############################################################################
subroutine set_NML_textcase(val)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "set_NML_textcase"
  ! 
continue
  !
  select case (val)
    case (0) ; NML_textcase = val  ! use titlecase
    case (1) ; NML_textcase = val  ! use lowercase
    case (2) ; NML_textcase = val  ! use uppercase
    case default
      if (mypnum == glb_root) write(iout,101) val
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  end select
  !
  101 format(/,"CODING ERROR: Invalid argument, val = ",i0)
end subroutine set_NML_textcase
!
!###############################################################################
!###############################################################################
! Specific routines to set a value.
!###############################################################################
!###############################################################################
!
!###############################################################################
!> Set an integer value.  Behavior depends on module variable `NMLaction`.
!###############################################################################
subroutine set_value_i(nam,var,val,act)
  !
  !.. Formal Arguments ..
  character(len=*),  intent(in)    :: nam  !< name
  integer,           intent(inout) :: var  !< variable
  integer,           intent(in)    :: val  !< value
  integer, optional, intent(in)    :: act  !< NMLaction
  !
  integer :: loc_action
  !
continue
  !
  loc_action = NMLaction
  if ( present(act) ) loc_action = act
  !
  select case (loc_action)
    case (NMLaction_SetDummy0)  ; var = idum0
    case (NMLaction_SetDummy1)  ; var = idum1
    case (NMLaction_SetValue )  ; var = val
    case default                ; ! do nothing
  end select
  !
  if ( setverbose .or. loc_action == NMLaction_WritePretty ) then
    call iow_value(nam,var)
  end if
  !
end subroutine set_value_i
!
!###############################################################################
!> Set a real value.  Behavior depends on module variable `NMLaction`.
!###############################################################################
subroutine set_value_r(nam,var,val,act)
  !
  !.. Formal Arguments ..
  character(len=*),  intent(in)    :: nam  !< name
  real(wp),          intent(inout) :: var  !< variable
  real(wp),          intent(in)    :: val  !< value
  integer, optional, intent(in)    :: act  !< NMLaction
  !
  integer :: loc_action
  !
continue
  !
  loc_action = NMLaction
  if ( present(act) ) loc_action = act
  !
  select case (loc_action)
    case (NMLaction_SetDummy0)  ; var = rdum0
    case (NMLaction_SetDummy1)  ; var = rdum1
    case (NMLaction_SetValue )  ; var = val
    case default                ; ! do nothing
  end select
  !
  if ( setverbose .or. loc_action == NMLaction_WritePretty ) then
    call iow_value(nam,var)
  end if
  !
end subroutine set_value_r
!
!###############################################################################
!> Set a character value.  Behavior depends on module variable `NMLaction`.
!###############################################################################
subroutine set_value_c(nam,var,val)
  !
  !.. Formal Arguments ..
  character(len=*),  intent(in)    :: nam  !< name
  character(len=*),  intent(inout) :: var  !< variable
  character(len=*),  intent(in)    :: val  !< value
! integer, optional, intent(in)    :: act  !< NMLaction
    ! Disabled optional action, because of generic interface could not
    ! distinguish set_value_c from set_value_s.
  !
  integer :: loc_action
  !
continue
  !
  loc_action = NMLaction
! if ( present(act) ) loc_action = act
  !
  select case (loc_action)
    case (NMLaction_SetDummy0)  ; var = cdum0
    case (NMLaction_SetDummy1)  ; var = cdum1
    case (NMLaction_SetValue )  ; var = val
    case default                ; ! do nothing
  end select
  !
  if ( setverbose .or. loc_action == NMLaction_WritePretty ) then
    call iow_value(nam,var)
  end if
  !
end subroutine set_value_c
!
!###############################################################################
!> Set a allocatable character string value and length.
!> Behavior depends on module variable `NMLaction`.
!###############################################################################
subroutine set_value_s(nam,var,val,char_len)
  !
  !.. Formal Arguments ..
  character(len=*),              intent(in)    :: nam       !< name
  character(len=:), allocatable, intent(inout) :: var       !< variable
  character(len=*),              intent(in)    :: val       !< value
  integer,                       intent(in)    :: char_len  !< char length
! integer,             optional, intent(in)    :: act       !< NMLaction
    ! Disabled optional action, because of generic interface could not
    ! distinguish set_value_c from set_value_s.
  !
  integer :: loc_action
  !
continue
  !
  loc_action = NMLaction
! if ( present(act) ) loc_action = act
  !
  select case (loc_action)
    case (NMLaction_SetDummy0)  ; call set_alloc_char(var,cdum0,char_len)
    case (NMLaction_SetDummy1)  ; call set_alloc_char(var,cdum1,char_len)
    case (NMLaction_SetValue )  ; call set_alloc_char(var,val  ,char_len)
    case default                ; ! do nothing
  end select
  !
  if ( setverbose .or. loc_action == NMLaction_WritePretty ) then
    call iow_value(nam,var)
  end if
  !
end subroutine set_value_s
!
!###############################################################################
!> Set a logical value.  Behavior depends on module variable `NMLaction`.
!###############################################################################
subroutine set_value_l(nam,var,val,act)
  !
  !.. Formal Arguments ..
  character(len=*),  intent(in)    :: nam  !< name
  logical,           intent(inout) :: var  !< variable
  logical,           intent(in)    :: val  !< value
  integer, optional, intent(in)    :: act  !< NMLaction
  !
  integer :: loc_action
  !
continue
  !
  loc_action = NMLaction
  if ( present(act) ) loc_action = act
  !
  select case (loc_action)
    case (NMLaction_SetDummy0)  ; var = ldum0
    case (NMLaction_SetDummy1)  ; var = ldum1
    case (NMLaction_SetValue )  ; var = val
    case default                ; ! do nothing
  end select
  !
  if ( setverbose .or. NMLaction == NMLaction_WritePretty ) then
    call iow_value(nam,var)
  end if
  !
end subroutine set_value_l
!
!###############################################################################
!> Set a logical value.  Behavior depends on module variable `NMLaction`.
!###############################################################################
subroutine set_value_lk(nam,var,val,act)
  !
  !.. Formal Arguments ..
  character(len=*),  intent(in)    :: nam  !< name
  logical(lk),       intent(inout) :: var  !< variable
  logical(lk),       intent(in)    :: val  !< value
  integer, optional, intent(in)    :: act  !< NMLaction
  !
  integer :: loc_action
  !
continue
  !
  loc_action = NMLaction
  if ( present(act) ) loc_action = act
  !
  select case (loc_action)
    case (NMLaction_SetDummy0)  ; var = ldum0
    case (NMLaction_SetDummy1)  ; var = ldum1
    case (NMLaction_SetValue )  ; var = val
    case default                ; ! do nothing
  end select
  !
  if ( setverbose .or. loc_action == NMLaction_WritePretty ) then
    call iow_value(nam,var)
  end if
  !
end subroutine set_value_lk
!
!###############################################################################
!> Set an integer value array.  Behavior depends on module variable `NMLaction`.
!###############################################################################
subroutine set_value_iarr(nam,var,val,act)
  !
  !.. Formal Arguments ..
  character(len=*),  intent(in)    :: nam     !< name
  integer,           intent(inout) :: var(:)  !< variable
  integer,           intent(in)    :: val(:)  !< value
  integer, optional, intent(in)    :: act     !< NMLaction
  !
  integer :: loc_action
  !
continue
  !
  loc_action = NMLaction
  if ( present(act) ) loc_action = act
  !
  select case (loc_action)
    case (NMLaction_SetDummy0)  ; var = idum0
    case (NMLaction_SetDummy1)  ; var = idum1
    case (NMLaction_SetValue )  ; var = val
    case default                ; ! do nothing
  end select
  !
  if ( setverbose .or. loc_action == NMLaction_WritePretty ) then
    call iow_value(nam,var)
  end if
  !
end subroutine set_value_iarr
!
!###############################################################################
!> Set a real value array.  Behavior depends on module variable `NMLaction`.
!###############################################################################
subroutine set_value_rarr(nam,var,val,act)
  !
  !.. Formal Arguments ..
  character(len=*),  intent(in)    :: nam     !< name
  real(wp),          intent(inout) :: var(:)  !< variable
  real(wp),          intent(in)    :: val(:)  !< value
  integer, optional, intent(in)    :: act     !< NMLaction
  !
  integer :: loc_action
  !
continue
  !
  loc_action = NMLaction
  if ( present(act) ) loc_action = act
  !
  select case (loc_action)
    case (NMLaction_SetDummy0)  ; var = rdum0
    case (NMLaction_SetDummy1)  ; var = rdum1
    case (NMLaction_SetValue )  ; var = val
    case default                ; ! do nothing
  end select
  !
  if ( setverbose .or. loc_action == NMLaction_WritePretty ) then
    call iow_value(nam,var)
  end if
  !
end subroutine set_value_rarr
!
!###############################################################################
!> Set allocatable character string to a given value and length.
!> Behavior does NOT depend on module variable `NMLaction`.
!###############################################################################
subroutine set_alloc_char(var,val,char_len)
  !
  !.. Formal Arguments ..
  !> Character variable to set.
  character(len=:), allocatable, intent(inout) :: var
  !> Character value to assign.
  character(len=*),              intent(in)    :: val
  !> Character length to assign. Must be non-negative.
  integer,                       intent(in)    :: char_len
  !
  !.. Local Variables ..
  integer :: default_len
  integer :: ierr
  integer :: n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "set_alloc_char"
  !
continue
  default_len = char_len
  !
  ! Deallocate the existing character string.
  !
  if ( allocated(var) ) then
    deallocate( var , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"var",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Allocate the character string to the desired length.
  !
  allocate( character(len=default_len) :: var , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,'var',1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Init to blanks.
  !
  do n = 1, default_len
    var(n:n) = ' '
  end do
  !
  ! Fill with value.
  !
  if ( default_len >= len(val) ) then
    ! Preserve any traililng blanks if default_len is
    ! longer than the length of the input val.
    var(1:len(val)) = val
  else
    ! Increase the length of var to the length of the input val
    ! so the entire contents of val can be copied into var.
    var = val  ! F2003 AUTO-REALLOC
  end if
  !
end subroutine set_alloc_char
!
!###############################################################################
!###############################################################################
! Specific routines to write a value.
!###############################################################################
!###############################################################################
!
!###############################################################################
!> Write namelist header to the unit given by module variable `iow`.
!>
!> The default option writes the header to a file buffer which is only
!> output once a request is made to write a namelist variable.
!> This prevents empty namelists.
!>
!> The force_write option will force the header to be written immediately.
!###############################################################################
subroutine iow_header(nam,force_write)
  !
  !.. Formal Arguments ..
  !> Namelist header (ie, '&input').
  character(len=*),  intent(in) :: nam
  !> Option that causes the header to be written immediately.
  logical, optional, intent(in) :: force_write
  !
  !.. Local Scalars ..
  logical :: opt_force_write
  !
continue
  !
  if (mypnum == glb_root) then
    !
    opt_force_write = fals
    if ( present(force_write) ) opt_force_write = force_write
    !
    ! Save the header so we know what block we are in.
    select case (NML_textcase)
      case (0) ; header =           trim(nam)
      case (1) ; header = lowercase(trim(nam))
      case (2) ; header = uppercase(trim(nam))
    end select
    !
    if ( opt_force_write ) then
      ! Force write the header now, so the buffer is empty.
      buffer = ''
      write(iow,'(a)') header
    else
      ! Write the header to the buffer. Buffer is output during iow_name.
      buffer = header
    end if
    !
  end if
  !
end subroutine iow_header
!
!###############################################################################
!> Write namelist footer to the unit given by module variable `iow`.
!>
!> The default option only writes the namelist terminator if the namelist
!> block has written some variables. This prevents empty namelist blocks.
!>
!> The force_write option will always write the namelist terminator.
!###############################################################################
subroutine iow_footer(force_write)
  !
  !.. Formal Arguments ..
  !> Option that causes the footer to be written immediately.
  logical, optional, intent(in) :: force_write
  !
  !.. Local Scalars ..
  logical :: opt_force_write, header_was_used
  !
continue
  !
  if (mypnum == glb_root) then
    !
    opt_force_write = fals
    if ( present(force_write) ) opt_force_write = force_write
    !
    header_was_used = ( buffer /= header)
    if ( opt_force_write .or. header_was_used ) then
      write(iow,'(a)') '/'
    end if
    !
  end if
  !
end subroutine iow_footer
!
!###############################################################################
!> Write variable name to the unit given by module variable `iow`.
!###############################################################################
subroutine iow_name(nam)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: nam
  !
  !.. Local Scalars ..
  character(len=100) :: text
  character(len=20) :: fmt
  !
continue
  !
  if (mypnum == glb_root) then
    !
    ! Flush the buffer first.
    !
    if ( buffer /= '' ) then
      write(iow,'(a)') buffer
      buffer = ''
    end if
    !
    ! Create the output format.
    !
    write(fmt,1) max(len(nam),NML_varname_width)
  1 format('(2x,a',i0,',a)')
    !
    ! Copy variable name into long text so it can be left-justified on output.
    text = nam
    if ( nam == '' ) then
      write(iow,fmt,advance='no') text,'   '
    else
      write(iow,fmt,advance='no') text,' = '
    end if
    !
  end if
  !
end subroutine iow_name
!
!###############################################################################
!> Write variable to unit iow. Scalar integer.
!###############################################################################
subroutine iow_value_i(nam,var)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: nam
  integer         , intent(in) :: var
  !
continue
  !
  if (mypnum == glb_root) then
    call iow_name(nam)
    write(iow,*) var
  end if
  !
end subroutine iow_value_i
!
!###############################################################################
!> Write variable to unit iow. Scalar real.
!###############################################################################
subroutine iow_value_r(nam,var)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: nam
  real(wp)        , intent(in) :: var
  !
continue
  !
  if (mypnum == glb_root) then
    call iow_name(nam)
    write(iow,*) var
  end if
  !
end subroutine iow_value_r
!
!###############################################################################
!> Write variable to unit iow. Scalar character.
!###############################################################################
subroutine iow_value_c(nam,var)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: nam
  character(len=*), intent(in) :: var
  !
continue
  !
  if (mypnum == glb_root) then
    call iow_name(nam)
    write(iow,*) "'" // trim(var) // "'"
  end if
  !
end subroutine iow_value_c
!
!###############################################################################
!> Write variable to unit iow. Scalar logical.
!###############################################################################
subroutine iow_value_l(nam,var)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: nam
  logical,          intent(in) :: var
  !
continue
  !
  if (mypnum == glb_root) then
    call iow_name(nam)
    write(iow,*) var
  end if
  !
end subroutine iow_value_l
!
!###############################################################################
!> Write variable to unit iow. Scalar logical kind.
!###############################################################################
subroutine iow_value_lk(nam,var)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: nam
  logical(lk),      intent(in) :: var
  !
continue
  !
  if (mypnum == glb_root) then
    call iow_name(nam)
    write(iow,*) var
  end if
  !
end subroutine iow_value_lk
!
!###############################################################################
!> Write variable to unit iow. Integer array.
!###############################################################################
subroutine iow_value_iarr(nam,var)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: nam
  integer         , intent(in) :: var(:)
  !
  !.. Local Scalars ..
  integer                      :: n
  character(len=1000)          :: text
  !
continue
  !
  if (mypnum == glb_root) then
    do n = 1, size(var)
      write(text,1) trim(nam),n
      call iow_value_i(trim(text),var(n))
    end do
  end if
  !
1 format(a,'(',i0,')')
  !
end subroutine iow_value_iarr
!
!###############################################################################
!> Write variable to unit iow. Real array.
!###############################################################################
subroutine iow_value_rarr(nam,var)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: nam
  real(wp)        , intent(in) :: var(:)
  !
  !.. Local Scalars ..
  integer                      :: n
  character(len=1000)          :: text
  !
continue
  !
  if (mypnum == glb_root) then
    do n = 1, size(var)
      write(text,1) trim(nam),n
      call iow_value_r(trim(text),var(n))
    end do
  end if
  !
1 format(a,'(',i0,')')
  !
end subroutine iow_value_rarr
!
!###############################################################################
!> Write variable to unit iow. Character array.
!###############################################################################
subroutine iow_value_carr(nam,var)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: nam
  character(len=*), intent(in) :: var(:)
  !
  !.. Local Scalars ..
  integer                      :: n
  character(len=1000)          :: text
  !
continue
  !
  if (mypnum == glb_root) then
    do n = 1, size(var)
      write(text,1) trim(nam),n
      call iow_value_c(trim(text),var(n))
    end do
  end if
  !
1 format(a,'(',i0,')')
  !
end subroutine iow_value_carr
!
!###############################################################################
!> Merge identical values of two namelist input files into a new input file.
!>
!> Assume the two namelist input files were written using NMLaction_WritePretty.
!> - Block begin starts at column 1 with '&'.
!> - Block end   starts at column 1 with '/'.
!> - Each name/value pair        is on a separate line.
!> - Each array element          is on a separate line.
!> - Each derived type component is on a separate line.
!> - No comments.
!> - No blank lines.
!###############################################################################
subroutine merge_nml_files()
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "merge_nml_files"
  !
  !.. Local Scalars ..
  integer :: io1,io2,io3,ioerr1,ioerr2,ierr
  integer :: n1,n2
  integer :: block_used
  integer :: line_num
  character(len=nml_recl), allocatable :: line1, block1, vname1
  character(len=nml_recl), allocatable :: line2, block2, vname2
  !
continue
  !
  if (mypnum == glb_root) then
    open(newunit=io1,file=new_nml_file1,form='formatted',      &
                     status='old',action='read',recl=nml_recl, &
                     iostat=ierr,iomsg=error_message)
    call io_error(pname,new_nml_file1,1,__LINE__,__FILE__,ierr,error_message)
    open(newunit=io2,file=new_nml_file2,form='formatted',      &
                     status='old',action='read',recl=nml_recl, &
                     iostat=ierr,iomsg=error_message)
    call io_error(pname,new_nml_file2,1,__LINE__,__FILE__,ierr,error_message)
    open(newunit=io3,file=new_nml_file3,form='formatted',           &
                     status='unknown',action='write',               &
                     iostat=ierr,iomsg=error_message)
    call io_error(pname,new_nml_file3,1,__LINE__,__FILE__,ierr,error_message)
    !
    call set_NMLaction(NMLaction_Write)
!   call set_NMLaction(NMLaction_WritePretty)
    !
    ! Adding nml_version block should not be necessary. NML writer does that.
!   write(io3,101) '&nml_version'
!   write(io3,102) '  namelist_version = ',latest_namelist_version
!   write(io3,101) '/'
    !
    line1  = ''
    line2  = ''
    block1 = ''
    block2 = ''
    vname1 = ''
    vname2 = ''
    block_used = 0
    !
    line_num = 0
    line_loop: &
    do
      ! Count lines in file.
      line_num = line_num + 1
      !
      ! Read next line of file.
      read(io1,101,iostat=ioerr1) line1
      read(io2,101,iostat=ioerr2) line2
      !
      ! Check for end of file.
      if ( ioerr1 < 0 .and. ioerr2 < 0 ) then
        ! Both files ended.
        exit line_loop
      else if ( ioerr1 < 0 .or. ioerr2 < 0 ) then
        ! One file ended.
        write(error_message,111) trim(new_nml_file1), trim(new_nml_file2)
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
      !
      ! Check for read error.
      if ( ioerr1 > 0 ) then
        write(error_message,121) line_num, trim(new_nml_file1)
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
      if ( ioerr2 > 0 ) then
        write(error_message,121) line_num, trim(new_nml_file2)
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
      !
      ! Check for namelist block begin.
      if ( line1(1:1) == '&' ) then
        block1 = line1
        vname1 = ''
        block_used = 0
      end if
      if ( line2(1:1) == '&' ) then
        block2 = line2
        vname2 = ''
      end if
      !
      ! Check for namelist block end.
      if ( line1(1:1) == '/' ) then
        if ( block_used == 1 ) write(io3,101) trim(line1)
        block_used = 0
        block1 = ''
        vname1 = ''
      end if
      if ( line2(1:1) == '/' ) then
        block2 = ''
        vname2 = ''
      end if
      !
      ! Stop if block names become out of sync.
      if ( block1 /= block2 ) then
        write(error_message,131) line_num, new_nml_file1, &
                                 line_num, new_nml_file2
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
      !
      ! Skip to next line if namelist block begin.
      if ( line1(1:1) == '&' ) cycle line_loop
      !
      ! Skip to next line if namelist block end.
      if ( line1(1:1) == '/' ) cycle line_loop
      !
      ! Check for namelist variable name by searching for =.
      n1 = index('=',line1) ; if ( n1 > 0 ) vname1 = line1(1:n1-1)
      n2 = index('=',line2) ; if ( n2 > 0 ) vname2 = line2(1:n2-1)
      !
      ! Stop if variable names become out of sync.
      if ( vname1 /= vname2 ) then
        write(error_message,141) line_num, new_nml_file1, &
                                 line_num, new_nml_file2
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
      !
      ! If lines match, write to output file.
      if ( line1 == line2 ) then
        ! Write namelist block name on first use.
        if ( block_used == 0 ) write(io3,101) trim(block1)
        block_used = 1
        ! Write namelist variable and value.
        write(io3,101) trim(line1)
      end if
      !
    end do line_loop
    !
    close(io1,status='delete',iostat=ierr,iomsg=error_message)
    call io_error(pname,new_nml_file1,2,__LINE__,__FILE__,ierr,error_message)
    close(io2,status='delete',iostat=ierr,iomsg=error_message)
    call io_error(pname,new_nml_file2,2,__LINE__,__FILE__,ierr,error_message)
    close(io3,status='keep'  ,iostat=ierr,iomsg=error_message)
    call io_error(pname,new_nml_file3,2,__LINE__,__FILE__,ierr,error_message)
    !
  end if
  !
  ! Format statements ($=new_line)
  !
101 format(999a)
102 format(a,i0)
111 format("Namelist files of different lengths:","$", &
           "  file ",a,                           "$", &
           "  file ",a                                 )
121 format("Namelist read error:",                "$", &
           "  line ",i0,", file ",a                    )
131 format("Namelist block names not in sync:",   "$", &
           "  line ",i0,", file ",a,              "$", &
           "  line ",i0,", file ",a                    )
141 format("Namelist variable names not in sync:","$", &
           "  line ",i0,", file ",a,              "$", &
           "  line ",i0,", file ",a                    )
  !
end subroutine merge_nml_files
!
!###############################################################################
!
end module nml_base_mod
