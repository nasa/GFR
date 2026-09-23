!###############################################################################
!> Module for version 0002 of namelist io.
!###############################################################################
module nml_0002_mod
  !
  !.. Use Statements ..
  use module_kind_types
  !
  use nml_base_mod, only : io_has_nml_block
  use nml_base_mod, only : set_ior, ior
  use nml_base_mod, only : set_iow, iow
  use nml_base_mod, only : set_NMLaction
  use nml_base_mod, only : iow_header
  use nml_base_mod, only : iow_footer
  use nml_base_mod, only : set => set_value
  use nml_base_mod, only : NMLaction_SetDummy0
  use nml_base_mod, only : NMLaction_SetDummy1
  use nml_base_mod, only : NMLaction_SetValue
  use nml_base_mod, only : NMLaction_Read
  use nml_base_mod, only : NMLaction_Write
  use nml_base_mod, only : NMLaction_WritePretty
  use nml_base_mod, only : NMLaction
  use nml_base_mod, only : namelist_procedures_t
  !
  implicit none
  !
  private ! Make everything private by default
  !
  private ::   alloc_nml_data
  private ::   resiz_nml_data
  private ::     set_nml_data
  private ::    read_nml_data
  private ::   write_nml_data
  private ::   write_nml_file
  !
  private :: ionml
  private :: ionml_grid
  private :: ionml_interpolation
  private :: ionml_numerics
  private :: ionml_timestep
  private :: ionml_initialization
  private :: ionml_sem
  private :: ionml_output
  private :: ionml_output_variables
  private :: ionml_physics
  private :: ionml_chemistry
  private :: ionml_limiter
  private :: ionml_reference
  private :: ionml_debug
  private :: ionml_testcase_vortex
  private :: ionml_testcase_tcf
  private :: ionml_testcase_sod
  private :: ionml_testcase_pulse
  private :: ionml_testcase_channel
  private :: ionml_bc
  private :: ionml_time_ave
  private :: ionml_restart_interupt
  private :: ionml_unknown
  !
  !
  ! ###########################
  ! ###  MODULE PARAMETERS  ###
  ! ###########################
  !
  integer, public, parameter :: module_namelist_version = 0002
  !
  !
  ! ##############################
  ! ###  MODULE DERIVED TYPES  ###
  ! ##############################
  !
  type, public, extends(namelist_procedures_t) :: nml_0002_procedures_t
  contains
    procedure, nopass :: init => init_nml0002_procedures
   !procedure :: convert => convert_nml0002_file
  end type nml_0002_procedures_t
  !
  !
  ! #################################
  ! ###  MODULE PUBLIC VARIABLES  ###
  ! #################################
  !
  type(nml_0002_procedures_t), public, protected, save :: nml_ver_0002
  !
contains
!
!###############################################################################
!
subroutine init_nml0002_procedures()
  !
  !.. Use Statements ..
  use nml_9999_mod, only : set_latest_nml_version
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_nml0002_procedures"
  !
continue
  !
  nml_ver_0002%version = module_namelist_version
  !
  nml_ver_0002%alloc_data => alloc_nml_data
  nml_ver_0002%resiz_data => resiz_nml_data
  nml_ver_0002%  set_data =>   set_nml_data
  nml_ver_0002% read_data =>  read_nml_data
  nml_ver_0002%write_data => write_nml_data
  nml_ver_0002% read_file =>  read_nml_file
  nml_ver_0002%write_file => write_nml_file
  !
  ! Send this namelist version into the subroutine set_latest_nml_version
  ! to set up the nml_ver_9999 procedure pointers.
  ! NOTE: This will only set the nml_ver_9999 procedure pointers
  !       to the respective procedures from this module if the
  !       namelist version of this module matches the value of
  !       latest_namelist_version from the module nml_base_mod.
  !
  call set_latest_nml_version(module_namelist_version,nml_ver_0002)
  !
end subroutine init_nml0002_procedures
!
!###############################################################################
!> Allocate namelist variables that are array of character.
!###############################################################################
subroutine alloc_nml_data( opt_n_bc_input,opt_char_len, &
                           opt_n_species,opt_n_point_probes, &
                           opt_n_rake_probes,opt_n_cell_probes)
  !
  !.. Use Statements ..
  use ovar, only : bc_input!,bc_input_t
  use ovar, only : point_probes
  use ovar, only : rake_probes
  use ovar, only : cell_probes
  use ovar, only : grid_file => gridfile
  use ovar, only : plot3d_cutfile
  use ovar, only : plot3d_bcsfile
  use ovar, only : interpto_grid_file => interpto_gridfile
  use ovar, only : interpto_plot3d_cutfile
  use ovar, only : interpto_plot3d_bcsfile
  use ovar, only : restart_file
  use ovar, only : solution_file
  use ovar, only : output_dir
  use ovar, only : lustre_stripe_size
  use ovar, only : time_ave_file
  use ovar, only : restart_list_file
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: opt_n_bc_input
  integer, optional, intent(in) :: opt_char_len
  integer, optional, intent(in) :: opt_n_species
  integer, optional, intent(in) :: opt_n_point_probes
  integer, optional, intent(in) :: opt_n_rake_probes
  integer, optional, intent(in) :: opt_n_cell_probes
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "alloc_nml_data"
  integer, parameter :: max_default_len = 150
  integer, parameter :: max_n_species = 7
  integer, parameter :: max_n_bc_input = 20
  integer, parameter :: max_n_point_probes = 200
  integer, parameter :: max_n_rake_probes = 10
  integer, parameter :: max_n_cell_probes = 10
  !
  !.. Local Scalars ..
  integer :: n,ierr
  integer :: default_len,n_species,n_bc_input
  integer :: n_point_probes,n_rake_probes,n_cell_probes
  character(len=50) :: array_name
  !
continue
  !
  n_bc_input = max_n_bc_input
  n_species = max_n_species
  n_point_probes = max_n_point_probes
  n_rake_probes = max_n_rake_probes
  n_cell_probes = max_n_cell_probes
  default_len = max_default_len
  !
  if (present(opt_n_bc_input))     n_bc_input     = opt_n_bc_input
  if (present(opt_n_species))      n_species      = opt_n_species
  if (present(opt_n_point_probes)) n_point_probes = opt_n_point_probes
  if (present(opt_n_rake_probes))  n_rake_probes  = opt_n_rake_probes
  if (present(opt_n_cell_probes))  n_cell_probes  = opt_n_cell_probes
  if (present(opt_char_len))       default_len    = opt_char_len
  !
  call debug_timer(entering_procedure,pname)
  !
  ! -----------------------------------------
  ! Make sure everything is deallocated first
  ! -----------------------------------------
  !
  ! Next line deallocates bc_input and component arrays.
  !   bc_input(:)
  !   bc_input(:)%mass_frac(:)
  !   bc_input(:)%species_name(:)
  if (allocated(bc_input)) then
    deallocate ( bc_input , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"bc_input",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(point_probes)) then
    deallocate ( point_probes , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"point_probes",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(rake_probes)) then
    deallocate ( rake_probes , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"rake_probes",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(cell_probes)) then
    deallocate ( cell_probes , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cell_probes",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  ! ------------------------------------
  ! Allocate bc_input array: bc_input(:)
  ! ------------------------------------
  !
  allocate ( bc_input(1:n_bc_input) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bc_input",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Allocate bc_input component arrays:
  !   bcinput(:)%mass_frac(:)
  !   bcinput(:)%species_name(:)
  !
#ifdef ALLOCATABLE_MASS_FRAC
  do n = 1,size(bc_input)
    allocate ( bc_input(n)%mass_frac(1:n_species) , source=zero , &
               stat=ierr , errmsg=error_message )
    write (array_name,2) n,"mass_frac"
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    allocate ( bc_input(n)%species_name(1:n_species) , source="" , &
               stat=ierr , errmsg=error_message )
    write (array_name,2) n,"species_name"
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
  end do
#endif
  !
  ! ---------------
  ! Allocate probes
  ! ---------------
  !
  ! Allocate point_probes before reading the namelist
  !
  allocate ( point_probes(1:n_point_probes) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"point_probes",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Allocate rake_probes before reading the namelist
  !
  allocate ( rake_probes(1:n_rake_probes) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"rake_probes",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Allocate cell_probes before reading the namelist
  !
  allocate ( cell_probes(1:n_cell_probes) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_probes",1,__LINE__,__FILE__,ierr,error_message)
  !
  !
  ! ------------------------------------------------------------------------
  ! Set allocatable character variables to their default values and lengths.
  ! ------------------------------------------------------------------------
  !
  call set('Grid_File              ', &
            grid_file               , ''        , default_len )
  call set('Plot3d_Cutfile         ', &
            plot3d_cutfile          , ''        , default_len )
  call set('Plot3d_BCSfile         ', &
            plot3d_bcsfile          , ''        , default_len )
  call set('Interpto_Grid_file     ', &
            interpto_grid_file      , ''        , default_len )
  call set('Interpto_Plot3d_Cutfile', &
            interpto_plot3d_cutfile , ''        , default_len )
  call set('Interpto_Plot3d_BCSfile', &
            interpto_plot3d_bcsfile , ''        , default_len )
  call set('Restart_File           ', &
            restart_file            , 'out.rst' , default_len )
  call set('Solution_File          ', &
            solution_file           , ''        , default_len )
  call set('Output_Dir             ', &
            output_dir              , '.'       , default_len )
  call set('Lustre_Stripe_Size     ', &
            lustre_stripe_size      , '4m'      , default_len )
  call set('Time_Ave_File          ', &
            time_ave_file           , 'out.ave' , default_len )
  call set('Restart_List_File      ', &
            restart_list_file       , ''        , default_len )
  !
  call debug_timer(leaving_procedure,pname)
  !
  2 format ("bc_input(",i0,")%",a)
  !
end subroutine alloc_nml_data
!
!###############################################################################
!> Resize allocatable namelist variables to minimal size for content.
!###############################################################################
subroutine resiz_nml_data()
  !
  !.. Use Statements ..
  use ovar, only : bc_input_t
  use ovar, only : bc_input
  use ovar, only : point_probe_t
  use ovar, only : rake_probe_t
  use ovar, only : point_probes
  use ovar, only : rake_probes
  use ovar, only : cell_probes
  use ovar, only : gridfile
  use ovar, only : plot3d_cutfile
  use ovar, only : plot3d_bcsfile
  use ovar, only : interpto_gridfile
  use ovar, only : interpto_plot3d_cutfile
  use ovar, only : interpto_plot3d_bcsfile
  use ovar, only : restart_file
  use ovar, only : solution_file
  use ovar, only : output_dir
  use ovar, only : lustre_stripe_size
  use ovar, only : time_ave_file
  use ovar, only : restart_list_file
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "resiz_nml_data"
  !
  !.. Local Scalars ..
  integer :: i,n,ierr
  !
  !.. Local Arrays ..
  type(bc_input_t),    allocatable :: tmp_bc_input(:)
  type(point_probe_t), allocatable :: tmp_point_probes(:)
  type(rake_probe_t),  allocatable :: tmp_rake_probes(:)
  type(point_probe_t), allocatable :: tmp_cell_probes(:)
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! ============================================================
  ! Reallocate bc_input to the actual size needed for the inputs
  ! ============================================================
  !
  ! Count how many of bc_input(:) are not default values.
  n = 0
  do i = 1,size(bc_input)
    if (bc_input(i)%is_not_default()) n = n + 1
  end do
  !
  ! Allocate a new bc_input(:) to hold this minimal content.
  allocate ( tmp_bc_input(1:n) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_bc_input",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Copy the non-default bc_input(:) into the new array.
  ! Note the non-default bc_input(:) might be scattered within bc_input(:).
  n = 0
  do i = 1,size(bc_input)
    if (bc_input(i)%is_not_default()) then
      n = n + 1
      tmp_bc_input(n) = bc_input(i)
    end if
  end do
  !
  ! Replace bc_input(:) with the new array.
  call move_alloc( from=tmp_bc_input , to=bc_input )
  !
  !
  ! ----------------------
  ! Trim character strings
  ! ----------------------
  !
  do n = 1,size(bc_input)
    bc_input(n)%sem%profile_file = trim(adjustl(bc_input(n)%sem%profile_file))
  end do
  !
  !#############################################################################
  !
  ! Reallocate point_probes to the actual size neede for the inputs
  !
  n = 0
  do i = 1, size(point_probes)
    if ( point_probes(i)%is_not_default() ) n = n + 1
  end do
  !
  allocate( tmp_point_probes(1:n) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_point_probes",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  n = 0
  do i = 1, size(point_probes)
    if (point_probes(i)%is_not_default()) then
      n = n + 1
      tmp_point_probes(n) = point_probes(i)
    end if
  end do
  !
  call move_alloc( from=tmp_point_probes , to=point_probes )
  !
  !#############################################################################
  !
  ! Reallocate rake_probes to the actual size neede for the inputs
  !
  n = 0
  do i = 1, size(rake_probes)
    if ( rake_probes(i)%is_not_default() ) n = n + 1
  end do
  !
  allocate( tmp_rake_probes(1:n) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_rake_probes",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  n = 0
  do i = 1, size(rake_probes)
    if (rake_probes(i)%is_not_default()) then
      n = n + 1
      tmp_rake_probes(n) = rake_probes(i)
    end if
  end do
  !
  call move_alloc( from=tmp_rake_probes , to=rake_probes )
  !
  !#############################################################################
  !
  ! Reallocate cell_probes to the actual size neede for the inputs
  !
  n = 0
  do i = 1, size(cell_probes)
    if ( cell_probes(i)%is_not_default() ) n = n + 1
  end do
  !
  allocate( tmp_cell_probes(1:n) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_cell_probes",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  n = 0
  do i = 1, size(cell_probes)
    if (cell_probes(i)%is_not_default()) then
      n = n + 1
      tmp_cell_probes(n) = cell_probes(i)
    end if
  end do
  !
  call move_alloc( from=tmp_cell_probes , to=cell_probes )
  !
  !#############################################################################
  !
  ! Following are F2003 AUTO REALLOCATIONS
  !
  gridfile                = trim(adjustl( gridfile                ))
  plot3d_cutfile          = trim(adjustl( plot3d_cutfile          ))
  plot3d_bcsfile          = trim(adjustl( plot3d_bcsfile          ))
  interpto_gridfile       = trim(adjustl( interpto_gridfile       ))
  interpto_plot3d_cutfile = trim(adjustl( interpto_plot3d_cutfile ))
  interpto_plot3d_bcsfile = trim(adjustl( interpto_plot3d_bcsfile ))
  restart_file            = trim(adjustl( restart_file            ))
  solution_file           = trim(adjustl( solution_file           ))
  output_dir              = trim(adjustl( output_dir              ))
  lustre_stripe_size      = trim(adjustl( lustre_stripe_size      ))
  time_ave_file           = trim(adjustl( time_ave_file           ))
  restart_list_file       = trim(adjustl( restart_list_file       ))
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine resiz_nml_data
!
!###############################################################################
!> Set variables for all blocks of namelist input.
!> Optionally echo values to unit number.
!###############################################################################
subroutine set_nml_data(action,io)
  !
  !.. Use Statements ..
  use nml_base_mod, only : setverbose
  !
  !.. Formal Arguments ..
  integer,           intent(in) :: action  !< An NMLaction.
  integer, optional, intent(in) :: io      !< Unit number to echo values.
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "set_nml_data"
  !
continue
  !
  setverbose = fals
  if ( present(io) ) then
    call set_iow(io)
    setverbose = true
  end if
  !
  select case (action)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue )
      call set_NMLaction(action)
      call ionml()
    case default
      write (error_message,101) action
      ! This is called from convert_nml_file within an 'if (mypnum == glb_root)'
      ! block, so this needs to use abort instead of stop_mpi.
      call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
  end select
  !
  setverbose = fals
  !
  101 format("Coding Error:  Invalue arg, action = ",i0)
  !
end subroutine set_nml_data
!
!###############################################################################
!> Read variables for all blocks of namelist input from unit number.
!###############################################################################
subroutine read_nml_data(io)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: io      !< Unit number for input.
  !
continue
  !
  rewind(io)
  call set_ior(io)
  call set_NMLaction(NMLaction_Read)
  call ionml()
  !
end subroutine read_nml_data
!
!###############################################################################
!> Write variables for all blocks of namelist input to unit number.
!###############################################################################
subroutine write_nml_data(action,io,banner)
  !
  !.. Use Modules ..
  use ovar, only : namelist_version
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: action  !< An NMLaction.
  integer,     intent(in) :: io      !< Unit number for output.
  logical(lk), intent(in) :: banner  !< Flag to write banner.
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_nml_data"
  !
  !.. Local Scalars ..
  integer :: ovar_namelist_version
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (mypnum == glb_root) then
    select case (action)
      case (NMLaction_Write, NMLaction_WritePretty)
        ! Set namelist_version to reflect this writer.
        ovar_namelist_version = namelist_version
        namelist_version = module_namelist_version
        !
        if (banner) write (io,104)
        call set_iow(io)
        call set_NMLaction(action)
        call ionml()
        if (banner) write (io,105)
        !
        ! Restore namelist_version.
        namelist_version = ovar_namelist_version
      case default
        write(error_message,106) action
        ! This is within an 'if (mypnum == glb_root)' block,
        ! so it needs to use abort instead of stop_mpi.
        call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
    end select
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  104 format (/,"#######################################",/, &
                "#####      INPUT PARAMETERS       #####",/, &
                "#######################################",/)
  105 format (/,"#######################################",/, &
                "#####   END OF INPUT PARAMETERS   #####",/, &
                "#######################################",/)
  106 format("Coding Error:  Invalid arg, action = ",i0)
  !
end subroutine write_nml_data
!
!###############################################################################
!> Initialize and read all blocks of namelist input from unit number.
!###############################################################################
subroutine read_nml_file(input_fname)
  !
  !.. Use Statements ..
  !
! Note:
! - nml_9999_* references are placeholders for the future.
! - Use of nml_9999_mod creates a circular makefile reference
!   because it is currently an alias to nml_0002_mod routines.
! - When nml_9999_mod eventually aliases nml_0003_mod, then
!   calls to nml_9999 routines will be needed here.
!
! use nml_9999_mod, only : alloc_nml9999_data, set_nml9999_data
! use nml_9999_mod, only : resiz_nml_9999_data
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: input_fname
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_nml_file"
  !
  !.. Local Scalars ..
  integer :: inpt, ierr
  !
continue
  !
  ! Allocate namelist variables to default sizes.
  !
! call alloc_nml9999_data()  ! alloc new arrays
  call alloc_nml_data()  ! alloc old arrays
  !
  ! Set namelist variables to default values.
  !
! call set_nml9999_data(NMLaction_SetValue)  ! set new arrays
  call set_nml_data(NMLaction_SetValue)  ! set old arrays
  !
  ! Open the input file.
  !
  call debug_timer(start_timer,"actually reading the input file")
  open (newunit=inpt,file=input_fname,status="old",action="read", &
                     form="formatted",iostat=ierr,iomsg=error_message)
  call io_error(pname,input_fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Read the input file namelist.
  !
  call read_nml_data(inpt)
  !
  ! Close the input file now that we are finished reading it.
  !
  close (inpt,iostat=ierr,iomsg=error_message)
  call io_error(pname,input_fname,2,__LINE__,__FILE__,ierr,error_message)
  call debug_timer(stop_timer,"actually reading the input file")
  !
  ! Resize all the allocatable variables within the input
  ! namelist to sizes more appropriate to their input values
  !
! call resiz_nml9999_data()  ! resize new arrays, don't care about old arrays
  call resiz_nml_data()  ! resize new arrays, don't care about old arrays
  !
end subroutine read_nml_file
!
!###############################################################################
!> Write all blocks of namelist input file.
!###############################################################################
subroutine write_nml_file(output_fname,action,requires_allocations,print_banner)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: output_fname
  integer,          intent(in) :: action        !< An NMLaction.
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: requires_allocations
  logical(lk), optional, intent(in) :: print_banner
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_nml_file"
  character(len=:), allocatable :: delim
  !
  !.. Local Scalars ..
  integer     :: io, ierr
  logical(lk) :: banner
  logical(lk) :: allocations_required
  !
continue
  !
  banner = fals
  if (present(print_banner)) banner = print_banner
  !
  allocations_required = fals
  if (present(requires_allocations)) allocations_required = requires_allocations
  !
!-DEBUG_BEG:
!-DEBUG_CMT: Initializing the namelist should be done outside this routine.
  ! Allocate all the variables within the input namelist
  ! to their default sizes and values
  !
  ! NOTE: Changed the default number of the different opt_n_* arguments to
  !       at least 2, because otherwise if any of the resulting arrays are
  !       allocated to only 1 element, the compiler (at least Intel does this)  
  !       will write out the array without array indices when directly writing
  !       the namelist variables using its native format. For example, it
  !       writes 
  !          BC_INPUT%... = ...,
  !       instead of
  !          BC_INPUT(1)%... = ...,
  !       which could be confusing to someone just looking at the contents of
  !       the written namelist file, or worse trying to use the contents of
  !       the file to create a new input file without knowing they need to
  !       add the array indices.
  !
  if (allocations_required) then
    call alloc_nml_data(opt_n_bc_input=2, &
                        opt_n_species=3,  &
                        opt_char_len=0, &
                        opt_n_point_probes=2, &
                        opt_n_rake_probes=2, &
                        opt_n_cell_probes=2)
  end if
!-DEBUG_END:
  !
  if (mypnum == glb_root) then
    !
    select case (action)
!     case (NMLaction_Read)        ; delim = 'NONE'
      case (NMLaction_Write)       ; delim = 'APOSTROPHE'
      !> delim is set to 'none' for WritePretty because the subroutine that
      !> writes a character string to the file will add its own apostrophes,
      !> so this prevents it from writing double apostrophes before and after
      !> the character string.
      case (NMLaction_WritePretty) ; delim = 'NONE'
      case default
        write(error_message,101) action
        ! This is within an 'if (mypnum == glb_root)' block,
        ! so it needs to use abort instead of stop_mpi.
        call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
    end select
    !
    open(newunit=io,file=trim(output_fname),status="replace", &
                    action="write", position="rewind", &
                    delim=delim, recl=5000, &
                    iostat=ierr,iomsg=error_message)
    call io_error(pname,trim(adjustl(output_fname)),1,__LINE__,__FILE__, &
                  ierr,error_message)
    !
    call write_nml_data(action,io,banner)
    !
    close(io,iostat=ierr,iomsg=error_message)
    call io_error(pname,trim(output_fname),2,__LINE__,__FILE__, &
                  ierr,error_message)
    !
    write(iout,201) trim(output_fname)
    !
  end if
  !
  101 format("Coding Error:  Invalid arg, action = ",i0)
  201 format(/,2x,"Namelist input file written to:  ",A)
  !
end subroutine write_nml_file
!
!###############################################################################
!> Workhorse routine to set/read/write variables for all namelist blocks.
!>
!> - Requires that allocatable namelist variables are already allocated and
!>   values are set for NMLaction, ior, iow, setverbose.
!> - Do not call directly, but rather through the set/read/write subroutines
!>   so that these prerequisites are met.
!###############################################################################
subroutine ionml()
  !
continue
  !
  call ionml_version          ()
  call ionml_grid             ()
  call ionml_interpolation    ()
  call ionml_reference        ()
  call ionml_physics          ()
  call ionml_chemistry        ()
  call ionml_sem              ()
  call ionml_initialization   ()
  call ionml_timestep         ()
  call ionml_numerics         ()
  call ionml_limiter          ()
  call ionml_time_ave         ()
  call ionml_output           ()
  call ionml_output_variables ()
  call ionml_utility          ()
  call ionml_testcase_vortex  ()
  call ionml_testcase_tcf     ()
  call ionml_testcase_sod     ()
  call ionml_testcase_pulse   ()
  call ionml_testcase_channel ()
  call ionml_debug            ()
  call ionml_bc               ()
  call ionml_removed          ()
  !
end subroutine ionml
!
!###############################################################################
!{{{ subroutine ionml_version
!> I/O for namelist block 'nml_version'.
!###############################################################################
subroutine ionml_version()
  !
  !.. Use Statements ..
  use ovar, only : namelist_version
  !
  namelist / nml_version / namelist_version
  !
continue
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&nml_version')
      call set('Namelist_Version',namelist_version, 2 )
      call iow_footer()
    case ( NMLaction_Read )
      if ( io_has_nml_block(ior,'&nml_version') ) then
        read(ior,nml=nml_version)
      end if
    case ( NMLaction_Write )
      write(iow,nml=nml_version)
  end select
  !
end subroutine ionml_version
!}}}
!###############################################################################
!{{{ subroutine ionml_grid
!> I/O for namelist block 'nml_grid'.
!###############################################################################
subroutine ionml_grid()
  !
  !.. Use Statements ..
  use ovar, only : grid_format
  use ovar, only : grid_file                => gridfile
  use ovar, only : plot3d_cutfile
  use ovar, only : plot3d_bcsfile
  use ovar, only : grid_scaling_factor
  use ovar, only : internal_grid_dof
  use ovar, only : internal_grid_order
  use ovar, only : plot3d_stride_i          => plot3d_i_stride
  use ovar, only : plot3d_stride_j          => plot3d_j_stride
  use ovar, only : plot3d_stride_k          => plot3d_k_stride
  use ovar, only : plot3d_stride_all        => plot3d_all_strides       
  use ovar, only : plot3d_agglomerate_order
  use ovar, only : correct_cgns_grid_points
  use ovar, only : force_cgns_linear_grid_points
  use ovar, only : metis_option, use_parmetis
  use ovar, only : geom_wts!,geom_wts_t
  !
  !.. Local Scalars ..
  integer :: n
  character(len=:), allocatable :: text
  !
  ! GRID INFORMATION
  namelist / nml_grid / grid_scaling_factor, grid_format, grid_file, &
                        plot3d_cutfile, plot3d_bcsfile
  ! PLOT3D COARSENING OPTIONS
  namelist / nml_grid / plot3d_stride_i, plot3d_stride_j, plot3d_stride_k, &
                        plot3d_stride_all, plot3d_agglomerate_order
  ! CGNS INPUT OPTIONS
  namelist / nml_grid / correct_cgns_grid_points, force_cgns_linear_grid_points
  ! NUMBER OF CELLS OR DOF FOR INTERNALLY GENERATED GRIDS
  namelist / nml_grid / internal_grid_dof
  ! ORDER OF CURVED ELEMENTS FOR INTERNALLY GENERATED GRIDS
  namelist / nml_grid / internal_grid_order
  ! PARALLEL OPTIONS
  namelist / nml_grid / metis_option, use_parmetis
  ! GEOMETRY WEIGHTS FOR METIS PARTITIONING
  namelist / nml_grid / geom_wts
  !
continue
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&nml_grid')
      call set('Grid_Scaling_Factor          ',                           &
                grid_scaling_factor           , 1.0_wp                    )
      call set('Grid_Format                  ',                           &
                grid_format                   , Internally_Generated_Grid )
      call set('Grid_File                    ',                           &
                grid_file                     , ''                        )
      call set('Plot3d_Cutfile               ',                           &
                plot3d_cutfile                , ''                        )
      call set('Plot3d_BCSfile               ',                           &
                plot3d_bcsfile                , ''                        )
      call set('Plot3d_Stride_I              ',                           &
                plot3d_stride_i               , -1                        )
      call set('Plot3d_Stride_J              ',                           &
                plot3d_stride_j               , -1                        )
      call set('Plot3d_Stride_K              ',                           &
                plot3d_stride_k               , -1                        )
      call set('Plot3d_Stride_All            ',                           &
                plot3d_stride_all             , -1                        )
      call set('Plot3d_Agglomerate_Order     ',                           &
                plot3d_agglomerate_order      , -1                        )
      call set('Correct_CGNS_Grid_Points     ',                           &
                correct_cgns_grid_points      , fals                      )
      call set('Force_CGNS_Linear_Grid_Points',                           &
                force_cgns_linear_grid_points , fals                      )
      call set('Internal_Grid_DoF            ',                           &
                internal_grid_dof             , 0                         )
      call set('Internal_Grid_Order          ',                           &
                internal_grid_order           , -1                        )
      !> METIS_OPTION = 3  ! k-way part, min comm volume
      call set('Metis_Option                 ',                           &
                metis_option                  , 3                         )
      call set('Use_Parmetis                 ',                           &
                use_parmetis                  , 0                         )
      do n = 1, size(geom_wts)
        ! Must set text large enough to write into.
        text = '                               '
        write(text,101) 'Geom_Wts(',n,')'
        text = trim(text)
        call set( text //   '%Tria         ',geom_wts(n)%tria       , 1 )
        call set( text //   '%Quad         ',geom_wts(n)%quad       , 1 )
        call set( text //   '%Tetr         ',geom_wts(n)%tetr       , 1 )
        call set( text //   '%Pyra         ',geom_wts(n)%pyra       , 1 )
        call set( text //   '%Pris         ',geom_wts(n)%pris       , 1 )
        call set( text //   '%Hexa         ',geom_wts(n)%hexa       , 1 )
        call set( text //   '%Face_Tria    ',geom_wts(n)%face_tria  , 1 )
        call set( text //   '%Face_Quad    ',geom_wts(n)%face_quad  , 1 )
      end do
      call iow_footer()
    case ( NMLaction_Read )
      if ( io_has_nml_block(ior,'&nml_grid') ) then
        read(ior,nml=nml_grid)
      end if
    case ( NMLaction_Write )
      write(iow,nml=nml_grid)
  end select
  !
  ! Format Statements
  !
  101 format(a,i0,a)
  !
end subroutine ionml_grid
!}}}
!###############################################################################
!{{{ subroutine ionml_interpolation
!> I/O for namelist block 'nml_interpolation'.
!###############################################################################
subroutine ionml_interpolation()
  !
  !.. Use Statements ..
  use ovar, only : interpto_grid_format
  use ovar, only : interpto_grid_file         => interpto_gridfile
  use ovar, only : interpto_plot3d_cutfile
  use ovar, only : interpto_plot3d_bcsfile
  use ovar, only : interpto_grid_scaling_factor
  use ovar, only : interpto_plot3d_stride_i   => interpto_plot3d_i_stride
  use ovar, only : interpto_plot3d_stride_j   => interpto_plot3d_j_stride
  use ovar, only : interpto_plot3d_stride_k   => interpto_plot3d_k_stride
  use ovar, only : interpto_plot3d_stride_all => &
                   interpto_plot3d_all_strides
  use ovar, only : interpto_plot3d_agglomerate_order
  use ovar, only : interpto_solution_order
  !
  ! INTERPOLATION UTILITY GRID INFORMATION
  namelist / nml_interpolation / &
                     interpto_grid_scaling_factor, &
                     interpto_grid_format, interpto_grid_file, &
                     interpto_plot3d_cutfile, interpto_plot3d_bcsfile
  ! INTERPOLATION UTILITY PLOT3D COARSENING OPTIONS
  namelist / nml_interpolation / &
                     interpto_plot3d_stride_i, interpto_plot3d_stride_j, &
                     interpto_plot3d_stride_k, interpto_plot3d_stride_all, &
                     interpto_plot3d_agglomerate_order
  ! INTERPOLATION UTILITY SOLUTION ORDER
  namelist / nml_interpolation / &
                     interpto_solution_order
  !
continue
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&nml_interpolation')
      call set('Interpto_Grid_Scaling_Factor     ',                &
                interpto_grid_scaling_factor      , 1.0_wp         )
      call set('Interpto_Grid_Format             ',                &
                interpto_grid_format              , Unknown_Format )
      call set('Interpto_Grid_File               ',                &
                interpto_grid_file                , ''             )
      call set('Interpto_Plot3d_Cutfile          ',                &
                interpto_plot3d_cutfile           , ''             )
      call set('Interpto_Plot3d_BCSfile          ',                &
                interpto_plot3d_bcsfile           , ''             )
      call set('Interpto_Plot3d_Stride_I         ',                &
                interpto_plot3d_stride_i          ,-1              )
      call set('Interpto_Plot3d_Stride_J         ',                &
                interpto_plot3d_stride_j          ,-1              )
      call set('Interpto_Plot3d_Stride_K         ',                &
                interpto_plot3d_stride_k          ,-1              )
      call set('Interpto_Plot3d_Stride_All       ',                &
                interpto_plot3d_stride_all        ,-1              )
      call set('Interpto_Plot3d_Agglomerate_Order',                &
                interpto_plot3d_agglomerate_order ,-1              )
      call set('Interpto_Solution_Order          ',                &
                interpto_solution_order           ,-3              )
      call iow_footer()
    case ( NMLaction_Read )
      if ( io_has_nml_block(ior,'&nml_interpolation') ) then
        read(ior,nml=nml_interpolation)
      end if
    case ( NMLaction_Write )
      write(iow,nml=nml_interpolation)
  end select
  !
end subroutine ionml_interpolation
!}}}
!###############################################################################
!{{{ subroutine ionml_numerics
!> I/O for namelist block 'nml_numerics'.
!###############################################################################
subroutine ionml_numerics()
  !
  !.. Use Statements ..
  use ovar, only : solution_order
  use ovar, only : correction_function
  use ovar, only : loc_solution_pts
  use ovar, only : loc_flux_pts
  use ovar, only : loc_triangle_pts
  use ovar, only : loc_tetrahedra_pts
  use ovar, only : invs_flux_method
  use ovar, only : visc_flux_method
  use ovar, only : visc_variables_to_ave
  use ovar, only : projection_order
  use ovar, only : quadrature_order
  use ovar, only : loc_quadrature_pts
  use ovar, only : combine_common_flux_projections
  use ovar, only : error_order
  use ovar, only : loc_error_pts
  use ovar, only : all_orders
  use ovar, only : all_points
  !
  ! FORCE GLOBAL ORDER
  namelist / nml_numerics / all_orders
  ! BASE FR OPTIONS
  namelist / nml_numerics / solution_order
  ! OVER-INTEGRATION OPTIONS
  namelist / nml_numerics / projection_order, quadrature_order
  ! ERROR-INTEGRATION OPTIONS
  namelist / nml_numerics / error_order
  !
  ! FORCE GLOBAL SOLUTION POINTS
  namelist / nml_numerics / all_points
  ! LOCATION OF SOLUTION, FLUX, AND TRIANGULAR POINTS
  namelist / nml_numerics / loc_solution_pts, loc_flux_pts
  ! OVER-INTEGRATION OPTIONS
  namelist / nml_numerics / loc_quadrature_pts
  ! ERROR-INTEGRATION OPTIONS
  namelist / nml_numerics / loc_error_pts
  ! LOCATION OF TRIANGULAR AND TETRAHEDRAL POINTS
  namelist / nml_numerics / loc_triangle_pts, loc_tetrahedra_pts
  !
  ! COMMON FLUX OPTIONS
  namelist / nml_numerics / &
                     invs_flux_method, visc_flux_method, visc_variables_to_ave
  ! OVER-INTEGRATION OPTIONS
  namelist / nml_numerics / combine_common_flux_projections
  ! BASE FR OPTIONS
  namelist / nml_numerics / correction_function
  !
continue
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&nml_numerics')
      call set('All_Orders                     ',                         &
                all_orders                      ,-3                       )
      call set('Solution_Order                 ',                         &
                solution_order                  ,-3                       )
      call set('Projection_Order               ',                         &
                projection_order                ,-3                       )
      call set('Quadrature_Order               ',                         &
                quadrature_order                ,-3                       )
      call set('Error_Order                    ',                         &
                error_order                     ,-3                       )
      !
      call set('All_Points                     ',                         &
                all_points                      ,-1                       )
      call set('Loc_Solution_Pts               ',                         &
                loc_solution_pts                , Legendre_Gauss          )
      call set('Loc_Flux_Pts                   ',                         &
                loc_flux_pts                    , Legendre_Gauss          )
      call set('Loc_Quadrature_Pts             ',                         &
                loc_quadrature_pts              , Legendre_Gauss          )
      call set('Loc_Error_Pts                  ',                         &
                loc_error_pts                   , Legendre_Gauss          )
      call set('Loc_Triangle_Pts               ',                         &
                loc_triangle_pts                , WilliamsShunn_TriPoints )
      call set('Loc_Tetrahedra_Pts             ',                         &
                loc_tetrahedra_pts              , ShunnHam_TetPoints      )
      !
      call set('Invs_Flux_Method               ',                         &
                invs_flux_method                , invs_flux_ROE           )
      call set('Visc_Flux_Method               ',                         &
                visc_flux_method                , visc_flux_BR2           )
      call set('Visc_Variables_To_Ave          ',                         &
                visc_variables_to_ave           , visc_interface_fluxes   )
      call set('Combine_Common_Flux_Projections',                         &
                combine_common_flux_projections , 0                       )
      call set('Correction_Function            ',                         &
                correction_function             , Correction_gDG          )
      !
      call iow_footer()
    case ( NMLaction_Read )
      if ( io_has_nml_block(ior,'&nml_numerics') ) then
        read(ior,nml=nml_numerics)
      end if
    case ( NMLaction_Write )
      write(iow,nml=nml_numerics)
  end select
  !
end subroutine ionml_numerics
!}}}
!###############################################################################
!{{{ subroutine ionml_timestep
!> I/O for namelist block 'nml_timestep'.
!###############################################################################
subroutine ionml_timestep()
  !
  !.. Use Statements ..
  use ovar, only : Runge_Kutta_Scheme
  use ovar, only : runge_kutta_stages        =>  num_rk_stages
  use ovar, only : num_timesteps_min         => minimum_timesteps
  use ovar, only : num_timesteps
  use ovar, only : Final_Time
  use ovar, only : Timestep_Type
  use ovar, only : timestep                  => constant_dt
  use ovar, only : cfl
  use ovar, only : cfl_ramping_beg           => cfl_beg
  use ovar, only : cfl_ramping_end           => cfl_end
  use ovar, only : cfl_ramping_iterations    => cfl_cycles
  use ovar, only : cfl_ramping_start_iter    => cfl_start_iter
  use ovar, only : convergence_order_abs_res
  use ovar, only : convergence_order_max_res
  !
  ! RUNGE-KUTTA OPTIONS
  namelist / nml_timestep / Runge_Kutta_Scheme, runge_kutta_stages
  ! COMMON FLUX OPTIONS
  ! TIME STEPPING INFORMATION
  namelist / nml_timestep / num_timesteps_min, num_timesteps, &
                            Final_Time, Timestep_Type, timestep
  ! CFL CONTROL VARIABLES
  namelist / nml_timestep / cfl, cfl_ramping_beg, cfl_ramping_end, &
                            cfl_ramping_iterations, cfl_ramping_start_iter
  ! CONVERGENCE CRITERIA
  namelist / nml_timestep / convergence_order_abs_res, &
                            convergence_order_max_res
  !
continue
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&nml_timestep')
      call set('Runge_Kutta_Scheme       ',                         &
                Runge_Kutta_Scheme        , OERK_11Stage_3rdOrder   )
      call set('Runge_Kutta_Stages       ',                         &
                runge_kutta_stages        , 3                       )
      call set('Num_Timesteps_Min        ',                         &
                num_timesteps_min         , 100                     )
      call set('Num_Timesteps            ',                         &
                num_timesteps             , 25000                   )
      call set('Final_Time               ',                         &
                Final_Time                ,-10.0_wp                 )
      call set('Timestep_Type            ',                         &
                Timestep_Type             , Global_Cell_Timestep    )
      call set('Timestep                 ',                         &
                timestep                  , 1.0e-6_wp               )
      call set('CFL                      ',                         &
                cfl                       , 1.0_wp                  )
      call set('CFL_Ramping_Beg          ',                         &
                cfl_ramping_beg           , 0.1_wp                  )
      call set('CFL_Ramping_End          ',                         &
                cfl_ramping_end           , 1.0_wp                  )
      call set('CFL_Ramping_Iterations   ',                         &
                cfl_ramping_iterations    , 1                       )
      call set('CFL_Ramping_Start_Iter   ',                         &
                cfl_ramping_start_iter    , 0                       )
      call set('Convergence_Order_Abs_Res',                         &
                convergence_order_abs_res ,-16._wp                  )
      call set('Convergence_Order_Max_Res',                         &
                convergence_order_max_res ,-16._wp                  )
      call iow_footer()
    case ( NMLaction_Read )
      if ( io_has_nml_block(ior,'&nml_timestep') ) then
        read(ior,nml=nml_timestep)
      end if
    case ( NMLaction_Write )
      write(iow,nml=nml_timestep)
  end select
  !
end subroutine ionml_timestep
!}}}
!###############################################################################
!{{{ subroutine ionml_initialization
!> I/O for namelist block 'nml_initialization'.
!###############################################################################
subroutine ionml_initialization()
  !
  !.. Use Statements ..
  use ovar, only : initialization => itestcase
  use ovar, only : jet_init_func
  use ovar, only : load_restart_file
  use ovar, only : restart_file
  use ovar, only : solution_format
  use ovar, only : solution_file
  !
  ! INITIALIZATION / TEST CASE
  namelist / nml_initialization / &
             initialization     , &
             jet_init_func
  ! RESTART OPTIONS
  namelist / nml_initialization / &
             load_restart_file  , &
             restart_file
  ! OPTIONS TO USE SOLUTION FILE FOR INITIAL CONDITIONS
  namelist / nml_initialization / &
             solution_file      , &
             solution_format
  !
continue
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&nml_initialization')
      call set('Initialization                   ',                       &
                initialization                    , Generic_MeanFlow      )
      call set('Jet_Init_Func                    ',                       &
                jet_init_func                     , 2                     )
      call set('Load_Restart_File                ',                       &
                load_restart_file                 , fals                  )
      call set('Restart_File                     ',                       &
                restart_file                      , 'out.rst'             )
      call set('Solution_File                    ',                       &
                solution_file                     , ''                    )
      call set('Solution_Format                  ',                       &
                solution_format                   , Unknown_Format        )
      call iow_footer()
    case ( NMLaction_Read )
      if ( io_has_nml_block(ior,'&nml_initialization') ) then
        read(ior,nml=nml_initialization)
      end if
    case ( NMLaction_Write )
      write(iow,nml=nml_initialization)
  end select
  !
end subroutine ionml_initialization
!}}}
!###############################################################################
!{{{ subroutine ionml_sem
!> I/O for namelist block 'nml_sem'.
!###############################################################################
subroutine ionml_sem()
  !
  !.. Use Statements ..
  use ovar, only : &
    use_profile_for_initialization => use_SEM_profile_for_initialization , &
    use_random_seed_from_restart   => SEM_use_random_seed_from_restart   , &
    use_spherical_eddy_volume      => SEM_eddy_volume_is_spherical       , &
    include_thermo_fluctuations    => SEM_include_thermo_fluctuations    , &
    box_num_eddies_multiplier      => SEM_box_num_eddies_multiplier      , &
    use_restart_eddies             => SEM_use_restart_eddies
  !
  ! OPTIONS FOR SEM
  namelist / nml_sem / use_profile_for_initialization , &
                       use_random_seed_from_restart   , &
                       use_spherical_eddy_volume      , &
                       include_thermo_fluctuations    , &
                       box_num_eddies_multiplier      , &
                       use_restart_eddies
  !
continue
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&nml_sem')
      call set('Use_Profile_For_Initialization',                            &
                use_profile_for_initialization , fals                       )
      call set('Use_Random_Seed_From_Restart  ',                            &
                use_random_seed_from_restart   , fals                       )
      call set('Use_Spherical_Eddy_Volume     ',                            &
                use_spherical_eddy_volume      , fals                       )
      call set('Include_Thermo_Fluctuations   ',                            &
                include_thermo_fluctuations    , fals                       )
      call set('Box_Num_Eddies_Multiplier     ',                            &
                box_num_eddies_multiplier      , one                        )
      call set('Use_Restart_Eddies            ',                            &
                use_restart_eddies             , true                       )
      call iow_footer()
    case ( NMLaction_Read )
      if ( io_has_nml_block(ior,'&nml_sem') ) then
        read(ior,nml=nml_sem)
      end if
    case ( NMLaction_Write )
      write(iow,nml=nml_sem)
  end select
  !
end subroutine ionml_sem
!}}}
!###############################################################################
!{{{ subroutine ionml_output
!> I/O for namelist block 'nml_output'.
!###############################################################################
subroutine ionml_output()
  !
  !.. Use Statements ..
  use ovar, only : &
    time_scaling_factor                , &
    restart_interval                   , &
    restart_frequency                  , &
    restart_output_disabled            => completely_disable_restart_output, &
    display_interval                   => iter_out_interval , &
    residual_interval                  => results_interval  , &
    output_summary_file                , &
    output_dir                         , &
    output_interval                    , &
    continuous_output                  , &
    loc_output_pts                     , &
    output_order                       , &
    dmat_grad_temp_order               , &
    output_frequency                   , &
    interpolate_before_output_variable , &
    cgns_high_order_output             , &
    cgns_mixed_output                  , &
    cgns_output_restart_time_ave_data  , &
    cgns_output_disabled               => completely_disable_cgns , &
    probes_collection_interval         , &
    probes_write_interval              , &
    point_probes                       , &
    rake_probes                        , &
    cell_probes                        , &
    profile_io_all                     , &
    profile_io_restart                 , &
    profile_io_cgns                    , &
    profile_io_cgns_subs               , &
    lustre_stripe_count                , &
    lustre_stripe_size
  !
  integer :: n
  real(wp), parameter, dimension(3) :: probes_xyz_default = &
                                       [-999._wp,-999._wp,-999._wp]
  character(len=:), allocatable :: text
  !
  ! TIME SCALING FOR OUTPUT
  namelist / nml_output / time_scaling_factor
  ! STATISTICAL OUTPUT OPTIONS
  namelist / nml_output / display_interval, residual_interval
  ! RESTART OPTIONS
  namelist / nml_output / restart_interval                  , &
                          restart_frequency                 , &
                          restart_output_disabled
  ! STATISTICAL OUTPUT OPTIONS
  namelist / nml_output / output_summary_file
  ! SOLUTION OUTPUT OPTIONS
  namelist / nml_output / output_dir, output_interval, output_frequency, &
                          output_order, dmat_grad_temp_order, &
                          loc_output_pts, continuous_output
  ! CGNS OUTPUT OPTIONS
  namelist / nml_output / interpolate_before_output_variable, &
                          cgns_high_order_output, cgns_mixed_output, &
                          cgns_output_restart_time_ave_data, &
                          cgns_output_disabled
  ! PROBES
  namelist / nml_output / probes_collection_interval, probes_write_interval, &
                          point_probes, rake_probes, cell_probes
  ! PROFILE PARALLEL IO
  namelist / nml_output / profile_io_all, profile_io_restart, &
                          profile_io_cgns, profile_io_cgns_subs
  ! LUSTRE FILE SYSTEM STRIPING PARAMETERS
  namelist / nml_output / lustre_stripe_count, lustre_stripe_size
  !
continue
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&nml_output')
      call set('Time_Scaling_Factor               ',                       &
                time_scaling_factor                , 1.0_wp                )
      call set('Display_Interval                  ',                       &
                display_interval                   , 1                     )
      call set('Residual_Interval                 ',                       &
                residual_interval                  , 0                     )
      call set('Restart_Interval                  ',                       &
                restart_interval                   , 0                     )
      call set('Restart_Frequency                 ',                       &
                restart_frequency                  ,-one                   )
      call set('Restart_Output_Disabled           ',                       &
                restart_output_disabled            , fals                  )
      call set('Output_Summary_File               ',                       &
                output_summary_file                , fals                  )
      call set('Output_Dir                        ',                       &
                output_dir                         , '.'                   )
      call set('Output_Interval                   ',                       &
                output_interval                    ,-5                     )
      call set('Output_Frequency                  ',                       &
                output_frequency                   ,-1.0_wp                )
      call set('Output_Order                      ',                       &
                output_order                       ,-3                     )
      call set('Dmat_Grad_Temp_Order              ',                       &
                dmat_grad_temp_order               ,-1                     )
      call set('Loc_Output_Pts                    ',                       &
                loc_output_pts                     , Equi_Distant_Points   )
      call set('Continuous_Output                 ',                       &
                continuous_output                  , fals                  )
      call set('Interpolate_Before_Output_Variable',                       &
                interpolate_before_output_variable , true                  )
      call set('CGNS_High_Order_Output            ',                       &
                cgns_high_order_output             , fals                  )
      call set('CGNS_Mixed_Output                 ',                       &
                cgns_mixed_output                  , true                  )
      call set('CGNS_Output_Restart_Time_Ave_Data ',                       &
                cgns_output_restart_time_ave_data  , fals                  )
      call set('CGNS_Output_Disabled              ',                       &
                cgns_output_disabled               , fals                  )
      !
      call set('Probes_Collection_Interval        ',                       &
                probes_collection_interval         , 0                     )
      call set('Probes_Write_Interval             ',                       &
                probes_write_interval              , 100                   )
      !
      do n = 1, size(point_probes)
        text = '                               '
        write(text,101) 'Point_Probes(',n,')'
        text = trim(text)   ! F2003 realloc
        call set( text //   '%XYZ                 ',                       &
                  point_probes(n)%xyz              , probes_xyz_default    )
      end do
      !
      do n = 1, size(rake_probes)
        text = '                               '
        write(text,101) 'Rake_Probes(',n,')'
        text = trim(text)   ! F2003 realloc
        call set( text //   '%Num_Probes          ',                       &
                  rake_probes(n)%num_probes        , 0                     )
        call set( text //   '%XYZ_Beg             ',                       &
                  rake_probes(n)%xyz_beg           , probes_xyz_default    )
        call set( text //   '%XYZ_End             ',                       &
                  rake_probes(n)%xyz_end           , probes_xyz_default    )
      end do
      !
      do n = 1, size(cell_probes)
        text = '                               '
        write(text,101) 'Cell_Probes(',n,')'
        text = trim(text)   ! F2003 realloc
        call set( text //   '%XYZ                 ',                       &
                  cell_probes(n)%xyz               , probes_xyz_default    )
      end do
      !
      call set('Profile_IO_All                    ',                       &
                profile_io_all                     , fals                  )
      call set('Profile_IO_Restart                ',                       &
                profile_io_restart                 , fals                  )
      call set('Profile_IO_CGNS                   ',                       &
                profile_io_cgns                    , fals                  )
      call set('Profile_IO_CGNS_Subs              ',                       &
                profile_io_cgns_subs               , fals                  )
      call set('Lustre_Stripe_Count               ',                       &
                lustre_stripe_count                , 64                    )
      call set('Lustre_Stripe_Size                ',                       &
                lustre_stripe_size                 , '4m'                  )
      call iow_footer()
    case ( NMLaction_Read )
      if ( io_has_nml_block(ior,'&nml_output') ) then
        read(ior,nml=nml_output)
      end if
    case ( NMLaction_Write )
      write(iow,nml=nml_output)
  end select
  !
  ! Format Statements
  !
  101 format(a,i0,a)
  !
end subroutine ionml_output
!}}}
!###############################################################################
!{{{ subroutine ionml_output_variables
!> I/O for namelist block 'nml_output_variables'.
!###############################################################################
subroutine ionml_output_variables()
  !
  !.. Use Statements ..
  use ovar,   only : &
      conserved_variables            => output_conserved_variables           , &
      momentum                       => output_momentum                      , &
      vorticity                      => output_vorticity                     , &
      cylindrical_coordinates        => output_cylindrical_coordinates       , &
      cylindrical_velocity           => output_cylindrical_velocity          , &
      relative_velocity              => output_relative_velocity             , &
      relative_conditions            => output_relative_conditions           , &
      thermo_variables               => output_thermo_variables              , &
      temperature                    => output_temperature                   , &
      gradients                      => output_gradients                     , &
      density_gradient               => output_density_gradient              , &
      density_gradient_magnitude     => output_density_gradient_magnitude    , &
      pressure_gradient              => output_pressure_gradient             , &
      pressure_gradient_magnitude    => output_pressure_gradient_magnitude   , &
      pressure_dmat_gradient         => output_pressure_dmat_gradient        , &
      temperature_gradient           => output_temperature_gradient          , &
      temperature_gradient_magnitude => output_temperature_gradient_magnitude, &
      temperature_dmat_gradient      => output_temperature_dmat_gradient     , &
      schlieren                      => output_schlieren                     , &
      modified_schlieren             => output_modified_schlieren            , &
      transport_variables            => output_transport_variables           , &
      stagnation_variables           => output_stagnation_variables          , &
      stagnation_pressure            => output_stagnation_pressure           , &
      stagnation_temperature         => output_stagnation_temperature        , &
      marker_variables               => output_marker_variables              , &
      qcriterion                     => output_qcriterion                    , &
      grid_details                   => output_grid_details                  , &
      local_timestep                 => output_local_timestep                , &
      species_density                => output_species_density               , &
      species_mass_fractions         => output_species_mass_fractions        , &
      species_mole_fractions         => output_species_mole_fractions        , &
      extra_species_variables        => output_extra_species_variables       , &
      debug_species_variables        => output_debug_species_variables       , &
      grdnt_species_variables        => output_grdnt_species_variables       , &
      local_filter_dissipation       => Track_Limiter_Filter_Dissipation     , &
      cellave_density                => output_cellave_density               , &
      exact_density                  => output_exact_density                 , &
      exact_pressure                 => output_exact_pressure

  !
  ! SELECTION OF OUTPUT VARIABLES IN THE CGNS SOLUTION FILES
  namelist / nml_output_variables / &
    conserved_variables           , &
    momentum                      , &
    vorticity                     , &
    cylindrical_coordinates       , &
    cylindrical_velocity          , &
    relative_velocity             , &
    relative_conditions           , &
    thermo_variables              , &
    temperature                   , &
    gradients                     , &
    density_gradient              , &
    density_gradient_magnitude    , &
    pressure_gradient             , &
    pressure_gradient_magnitude   , &
    pressure_dmat_gradient        , &
    temperature_gradient          , &
    temperature_gradient_magnitude, &
    temperature_dmat_gradient     , &
    schlieren                     , &
    modified_schlieren            , &
    transport_variables           , &
    stagnation_variables          , &
    stagnation_pressure           , &
    stagnation_temperature        , &
    marker_variables              , &
    qcriterion                    , &
    grid_details                  , &
    local_timestep                , &
    species_density               , &
    species_mass_fractions        , &
    species_mole_fractions        , &
    extra_species_variables       , &
    debug_species_variables       , &
    grdnt_species_variables       , &
    local_filter_dissipation      , &
    cellave_density               , &
    exact_density                 , &
    exact_pressure
  !
continue
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&nml_output_variables')
      call set('Conserved_Variables           ',      &
                conserved_variables            , fals )
      call set('Momentum                      ',      &
                momentum                       , fals )
      call set('Vorticity                     ',      &
                vorticity                      , fals )
      call set('Cylindrical_Coordinates       ',      &
                cylindrical_coordinates        , fals )
      call set('Cylindrical_Velocity          ',      &
                cylindrical_velocity           , fals )
      call set('Relative_Velocity             ',      &
                relative_velocity              , fals )
      call set('Relative_Conditions           ',      &
                relative_conditions            , fals )
      call set('Thermo_Variables              ',      &
                thermo_variables               , fals )
      call set('Temperature                   ',      &
                temperature                    , fals )
      call set('Gradients                     ',      &
                gradients                      , fals )
      call set('Density_Gradient              ',      &
                density_gradient               , fals )
      call set('Density_Gradient_Magnitude    ',      &
                density_gradient_magnitude     , fals )
      call set('Pressure_Gradient             ',      &
                pressure_gradient              , fals )
      call set('Pressure_Gradient_Magnitude   ',      &
                pressure_gradient_magnitude    , fals )
      call set('Pressure_Dmat_Gradient        ',      &
                pressure_dmat_gradient         , fals )
      call set('Temperature_Gradient          ',      &
                temperature_gradient           , fals )
      call set('Temperature_Gradient_Magnitude',      &
                temperature_gradient_magnitude , fals )
      call set('Temperature_Dmat_Gradient     ',      &
                temperature_dmat_gradient      , fals )
      call set('Schlieren                     ',      &
                schlieren                      , fals )
      call set('Modified_Schlieren            ',      &
                modified_schlieren             , fals )
      call set('Transport_Variables           ',      &
                transport_variables            , fals )
      call set('Stagnation_Variables          ',      &
                stagnation_variables           , fals )
      call set('Stagnation_Pressure           ',      &
                stagnation_pressure            , fals )
      call set('Stagnation_Temperature        ',      &
                stagnation_temperature         , fals )
      call set('Marker_Variables              ',      &
                marker_variables               , fals )
      call set('Qcriterion                    ',      &
                qcriterion                     , fals )
      call set('Grid_Details                  ',      &
                grid_details                   , fals )
      call set('Local_Timestep                ',      &
                local_timestep                 , fals )
      call set('Species_Density               ',      &
                species_density                , fals )
      call set('Species_Mass_Fractions        ',      &
                species_mass_fractions         , true )
      call set('Species_Mole_Fractions        ',      &
                species_mole_fractions         , fals )
      call set('Extra_Species_Variables       ',      &
                extra_species_variables        , fals )
      call set('Debug_Species_Variables       ',      &
                debug_species_variables        , fals )
      call set('Grdnt_Species_Variables       ',      &
                grdnt_species_variables        , fals )
      call set('Local_Filter_Dissipation      ',      &
                local_filter_dissipation       , fals )
      call set('CellAve_Density               ',      &
                cellave_density                , fals )
      call set('Exact_Density                 ',      &
                exact_density                  , fals )
      call set('Exact_Pressure                ',      &
                exact_pressure                 , fals )
      call iow_footer()
    case ( NMLaction_Read )    
      if ( io_has_nml_block(ior,'&nml_output_variables') ) then
        read(ior,nml=nml_output_variables)
      end if
    case ( NMLaction_Write )   
      write(iow,nml=nml_output_variables)
  end select
  !
end subroutine ionml_output_variables
!}}}
!###############################################################################
!{{{ subroutine ionml_physics
!> I/O for namelist block 'nml_physics'.
!###############################################################################
subroutine ionml_physics()
  !
  !.. Use Statements ..
  use ovar, only : governing_equations, moving_ref_frame
  !
  ! GOVERNING EQUATIONS
  namelist / nml_physics / governing_equations
  ! MOVING REFERENCE FRAME
  namelist / nml_physics / moving_ref_frame
  !
continue
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&nml_physics')
      call set('Governing_Equations          ',                             &
                governing_equations           , NavierStokes_Eqns           )
      call set('Moving_Ref_Frame%Rotation    ',                             &
                moving_ref_frame%rotation     , [zero,zero,zero]            )
      call set('Moving_Ref_Frame%Translation ',                             &
                moving_ref_frame%translation  , [zero,zero,zero]            )
      call set('Moving_Ref_Frame%Spinup%Iter ',                             &
                moving_ref_frame%spinup%iter  , 0                           )
      call set('Moving_Ref_Frame%Spinup%Tstar',                             &
                moving_ref_frame%spinup%tstar , zero                        )
      call iow_footer()
    case ( NMLaction_Read )    
      if ( io_has_nml_block(ior,'&nml_physics') ) then
        read(ior,nml=nml_physics)
      end if
    case ( NMLaction_Write )
      write(iow,nml=nml_physics)
  end select
  !
end subroutine ionml_physics
!}}}
!###############################################################################
!{{{ subroutine ionml_chemistry
!> I/O for namelist block 'nml_chemistry
!###############################################################################
subroutine ionml_chemistry()
  !
  !.. Use Statements ..
  use ovar,    only :        disable_chemical_reactions         , &
                             semi_implicit_chemistry            , &
                             ignition_method                    , &
                             ignition_start_timestep            , &
                             ignition_duration                  , &
                             ignition_prod2reac_L2_ratio_goal   , &
                             ignition_max_temperature           , &
                             ignition_spark_cap                 , &
                             ignition_source_magnitude          , &
                             ignition_source_radius             , &
                             ignition_source_location
  !
  ! CHEMICAL REACTION SOURCE TERMS
  namelist / nml_chemistry / disable_chemical_reactions         , &
                             semi_implicit_chemistry
                            !reactions_start_iter
  !
  ! IGNITION CONTROL SOURCE TERMS
  namelist / nml_chemistry / ignition_method                    , &
                             ignition_start_timestep            , &
                             ignition_duration                  , &
                             ignition_prod2reac_L2_ratio_goal   , &
                             ignition_max_temperature           , &
                             ignition_spark_cap                 , &
                             ignition_source_magnitude          , &
                             ignition_source_radius             , &
                             ignition_source_location
  !
continue
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&nml_chemistry')
      call set('Disable_Chemical_Reactions        ',                      &
                disable_chemical_reactions         , fals                 )
      call set('Semi_Implicit_Chemistry           ',                      &
                semi_implicit_chemistry            , 0                    )
      call set('Ignition_Method                   ',                      &
                Ignition_Method                    , IGN_OFF              )
      call set('Ignition_Start_Timestep           ',                      &
                ignition_start_timestep            , 0                    )
      call set('Ignition_Duration                 ',                      &
                ignition_duration                  , 1                    )
      call set('Ignition_Prod2Reac_L2_Ratio_Goal  ',                      &
                ignition_prod2reac_L2_ratio_goal   ,-1.0_wp               )
      call set('Ignition_Max_Temperature          ',                      &
                ignition_max_temperature           , 1600.0_wp            )
      call set('Ignition_Spark_Cap                ',                      &
                ignition_spark_cap                 , 1.0e+12_wp           )
      call set('Ignition_Source_Magnitude         ',                      &
                ignition_source_magnitude          , 20.0_wp              )
      call set('Ignition_Source_Radius            ',                      &
                ignition_source_radius             , zero                 )
      call set('Ignition_Source_Location          ',                      &
                ignition_source_location           , [zero,zero,zero]     )
      call iow_footer()
    case ( NMLaction_Read )
      if ( io_has_nml_block(ior,'&nml_chemistry') ) then
        read(ior,nml=nml_chemistry)
      end if
    case ( NMLaction_Write )  
      write(iow,nml=nml_chemistry)
  end select
  !
end subroutine ionml_chemistry
!}}}
!###############################################################################
!{{{ subroutine ionml_limiter
!> I/O for namelist block 'nml_limiter'.
!###############################################################################
subroutine ionml_limiter()
  !
  !.. Use Statements ..
  use ovar, only : &
      Global_Filter_Interval           , &
      Global_Filter_Option             => Filter_Option         , &
      Global_Filter_Order              => Filter_Order          , &
      Global_Filter_Alpha              => Filter_Alpha          , &
      Global_Filter_Etac               => Filter_etac           , &
      Sponge_Filter_Interval           , &
      Sponge_Filter_Option             , &
      Sponge_Filter_Order              , &
      Sponge_Filter_Alpha              , &
      Sponge_Filter_Etac               , &
      Sponge_Filter_X_Start            => Sponge_Beyond_X , &
      Local_Limiter_Interval           => Limiter_Interval      , &
      Local_Limiter_Option             => Limiter_Option        , &
      Local_Filter_Option              => Limiter_Filter_Option , &
      Local_Filter_Order               => Limiter_Filter_Order  , &
      Local_Filter_Alpha               => Limiter_Filter_Alpha  , &
      Local_Filter_Etac                => Limiter_Filter_Etac   , &
      Marker_Stored_Interval           , &
      Marker_Option                    , &
      Marker_Variable                  , &
      Marker_Min                       => Limiter_Marker_Min    , &
      Marker_Max                       => Limiter_Marker_Max    , &
      Marker_dpsi                      => Limiter_Marker_dpsi   , &
      thermofix                        , &
      Positivity                       , &
      Entropy_Limiter                  , &
      Marker_MinMax_Default
  !
  ! FILTER OPTIONS
  namelist / nml_limiter / Global_Filter_Interval
  namelist / nml_limiter / Global_Filter_Option , &
                           Global_Filter_Order  , &
                           Global_Filter_Alpha  , &
                           Global_Filter_etac
  ! SPONGE FILTER OPTIONS
  namelist / nml_limiter / Sponge_Filter_Interval
  namelist / nml_limiter / Sponge_Filter_Option  , &
                           Sponge_Filter_Order   , &
                           Sponge_Filter_Alpha   , &
                           Sponge_Filter_etac    , &
                           Sponge_Filter_x_start
  ! LOCAL LIMITER OPTIONS
  namelist / nml_limiter / Local_Limiter_Interval
  namelist / nml_limiter / Local_Limiter_Option
  ! LOCAL FILTER OPTIONS
  namelist / nml_limiter / &
                           Local_Filter_Option , &
                           Local_Filter_Order  , &
                           Local_Filter_Alpha  , &
                           Local_Filter_etac
  ! MARKER OPTIONS
  namelist / nml_limiter / Marker_Stored_Interval
  namelist / nml_limiter / Marker_Option
  namelist / nml_limiter / Marker_Variable
  ! LIMITER MARKER SENSITIVITY OPTIONS
  namelist / nml_limiter / Marker_Min  , &
                           Marker_Max  , &
                           Marker_dpsi
  ! CHECK/FIX BAD THERMO PROPERTIES (PRES,TEMP,DENS)
  namelist / nml_limiter / thermofix
  ! POSITIVITY PRESERVING LIMITER OPTIONS AND SETTINGS
  namelist / nml_limiter / positivity
  ! ENTROPY LIMITER OPTIONS AND SETTINGS
  namelist / nml_limiter / entropy_limiter
  !
continue
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&nml_limiter')
!-rm  Changed:  Filtering_Interval -> Global_Filter_Interval
!-rm  call set('Filtering_Interval              ',                          &
!-rm            Filtering_Interval               , 0                        )
      call set('Global_Filter_Interval          ',                          &
                Global_Filter_Interval           , 0                        )
      call set('Global_Filter_Option            ',                          &
                Global_Filter_Option             , 0                        )
      call set('Global_Filter_Order             ',                          &
                Global_Filter_Order              , 16                       )
      call set('Global_Filter_Alpha             ',                          &
                Global_Filter_Alpha              , 36.0_wp                  )
      call set('Global_Filter_Etac              ',                          &
                Global_Filter_Etac               , 0.0_wp                   )
!-rm  Changed:  Use_Sponge_Filter -> Sponge_Filter_Option
!-rm  call set('Use_Sponge_Filter               ',                          &
!-rm            use_sponge_filter                , fals                     )
      call set('Sponge_Filter_Interval          ',                          &
                Sponge_Filter_Interval           , 0                        )
      call set('Sponge_Filter_Option            ',                          &
                Sponge_Filter_Option             , 0                        )
      call set('Sponge_Filter_Order             ',                          &
                Sponge_Filter_Order              , 8                        )
      call set('Sponge_Filter_Alpha             ',                          &
                Sponge_Filter_Alpha              , 36.0_wp                  )
      call set('Sponge_Filter_Etac              ',                          &
                Sponge_Filter_Etac               , 0.0_wp                   )
      call set('Sponge_Filter_X_Start           ',                          &
                Sponge_Filter_x_start            , huge(zero)               )
      call set('Local_Limiter_Interval          ',                          &
                Local_Limiter_Interval           , 0                        )
!-rm  Changed:  Limiting_Interval -> Local_Limiter_Interval
!-rm  call set('Limiting_Interval               ',                          &
!-rm            Limiting_Interval                , 0                        )
      call set('Local_Limiter_Option            ',                          &
                Local_Limiter_Option             , 0                        )
!-rm  Changed:  Limiter_Use_Dynamic_Filter -> Local_Limiter_Option
!-rm  call set('Limiter_Use_Dynamic_Filter      ',                          &
!-rm            Limiter_Use_Dynamic_Filter       , fals                     )
      call set('Local_Filter_Option             ',                          &
                Local_Filter_Option              , 0                        )
      call set('Local_Filter_Order              ',                          &
                Local_Filter_Order               , 16                       )
      call set('Local_Filter_Alpha              ',                          &
                Local_Filter_Alpha               , 36.0_wp                  )
      call set('Local_Filter_Etac               ',                          &
                Local_Filter_Etac                , 0.0_wp                   )
      call set('Marker_Stored_Interval          ',                          &
                Marker_Stored_Interval           , 1                        )
      call set('Marker_Option                   ',                          &
                Marker_Option                    , 0                        )
      call set('Marker_Variable                 ',                          &
                Marker_Variable                  , 0                        )
      call set('Marker_Min                      ',                          &
                Marker_Min                       , Marker_MinMax_Default    )
      call set('Marker_Max                      ',                          &
                Marker_Max                       , Marker_MinMax_Default    )
      call set('Marker_Dpsi                     ',                          &
                Marker_dpsi                      ,-1.0_wp                   )
      !
      call set('Thermofix%Method                ',                          &
                thermofix%method                 , TFM_OFF                  )
      call set('Thermofix%Iterations            ',                          &
                thermofix%iterations             , 10                       )
      call set('Thermofix%Check_Cells           ',                          &
                thermofix%check_cells            , true                     )
      call set('Thermofix%Check_Faces           ',                          &
                thermofix%check_faces            , true                     )
      call set('Thermofix%Check_Nodes           ',                          &
                thermofix%check_nodes            , true                     )
      call set('Thermofix%Verbose               ',                          &
                thermofix%verbose                , true                     )
      call set('Thermofix%Scale_Factor          ',                          &
                thermofix%scale_factor           , 0.5_wp                   )
      call set('Thermofix%Min_Pressure          ',                          &
                thermofix%min_pressure           , 10000._wp                )
      call set('Thermofix%Min_Temperature       ',                          &
                thermofix%min_temperature        ,   100._wp                )
      call set('Thermofix%Min_Density           ',                          &
                thermofix%min_density            ,  0.001_wp                )
      !
      call set('Positivity%Option               ',                          &
                Positivity%Option                , POSLIM_OPT__OFF          )
      call set('Positivity%Density_Limit        ',                          &
                Positivity%density_Limit         , eps6                     )
      call set('Positivity%Pressure_Limit       ',                          &
                Positivity%pressure_Limit        , eps6                     )
      !
      call set('Entropy_Limiter%Option          ',                          &
                Entropy_Limiter%Option           , ENTLIM_OPT__OFF          )
      call set('Entropy_Limiter%Interpolant     ',                          &
                Entropy_Limiter%Interpolant      , ENTLIM_FUN__LinearLvIhme )
      call set('Entropy_Limiter%Tolerance       ',                          &
                Entropy_Limiter%Tolerance        , eps3                     )
      call set('Entropy_Limiter%Check_Faces     ',                          &
                Entropy_Limiter%Check_Faces      , fals                     )
      call set('Entropy_Limiter%Check_Nodes     ',                          &
                Entropy_Limiter%Check_Nodes      , fals                     )
      call set('Entropy_Limiter%Check_Global    ',                          &
                Entropy_Limiter%Check_Global     , fals                     )
      !
      call iow_footer()
    case ( NMLaction_Read )
      if ( io_has_nml_block(ior,'&nml_limiter') ) then
        read(ior,nml=nml_limiter)
      end if
    case ( NMLaction_Write )
      write(iow,nml=nml_limiter)
  end select
  !
end subroutine ionml_limiter
!}}}
!###############################################################################
!{{{ subroutine ionml_reference
!> I/O for namelist block 'nml_reference'.
!###############################################################################
subroutine ionml_reference()
  !
  !.. Use Statements ..
  use ovar,    only : mach            => machref     , &
                      reynolds        => reyref      , &
                      gamma           => gam         , &
                      density         => rhoref      , &
                      temperature     => tref        , &
                      pressure        => pref        , &
                      gas_constant    => rgasref     , &
                      aoa_alpha       => aoa_alpharef, &
                      aoa_beta        => aoa_betaref , &
                      prandtl         => Pr          , &
                      schmidt         => schmidt_lam , &
                      suth_muref, suth_Tref, suth_Sref
  !
  ! NONDIMENSIONAL REFERENCE CONDITIONS
  namelist / nml_reference / mach, reynolds, gamma
  ! REFERENCE CONDITIONS
  namelist / nml_reference / gas_constant
  namelist / nml_reference / density, temperature, pressure
  ! REFERENCE ANGLES OF ATTACK
  namelist / nml_reference / aoa_alpha, aoa_beta
  ! LAMINAR PRANDTL NUMBER
  namelist / nml_reference / prandtl
  ! LAMINAR SCHMIDT NUMBER
  namelist / nml_reference / schmidt
  ! SUTHERLANDS LAW COEFFICIENTS
  namelist / nml_reference / suth_muref, suth_Tref, suth_Sref
  !
continue
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&nml_reference')
      call set('Mach           ',mach           , 0.4_wp        )
      call set('Reynolds       ',reynolds       ,-1.0_wp        )
      call set('Gamma          ',gamma          , 1.4_wp        )
      call set('Gas_Constant   ',gas_constant   ,-1.0_wp        )
      call set('Density        ',density        ,-1.0_wp        )
      call set('Temperature    ',temperature    , 300.0_wp      )
      call set('Pressure       ',pressure       , 100000.0_wp   )
      call set('AoA_Alpha      ',aoa_alpha      , 0.0_wp        )
      call set('AoA_Beta       ',aoa_beta       , 0.0_wp        )
      call set('Prandtl        ',prandtl        , 0.72_wp       )
      call set('Schmidt        ',schmidt        , 0.72_wp       )
      call set('Suth_Muref     ',suth_muref     , 1.716e-5_wp   )
      call set('Suth_Tref      ',suth_Tref      , 273.0_wp      )
      call set('Suth_Sref      ',suth_Sref      , 110.4_wp      )
      call iow_footer()
    case ( NMLaction_Read )
      if ( io_has_nml_block(ior,'&nml_reference') ) then
        read(ior,nml=nml_reference)
      end if
    case ( NMLaction_Write )
      write(iow,nml=nml_reference)
  end select
  !
end subroutine ionml_reference
!}}}
!###############################################################################
!{{{ subroutine ionml_debug
!> I/O for namelist block 'nml_debug'.
!###############################################################################
subroutine ionml_debug()
  !
  !.. Use Statements ..
  use ovar, only : &
      output_limiter_statistics     , &
      output_filter_diagonals       , &
      output_plot3d_to_tecplot      , &
      output_bface_array            , &
      check_flux_point_coordinates  , &
      force_nd_subroutines          , &
      newton_convergence_goal       => conv_temp  , &
      chemistry_sub_iterations      , &
      ignore_thermo_newton_warnings
  !
  ! OUTPUT THE STATISTICS OF THE LIMITER
  namelist / nml_debug / output_limiter_statistics
  ! OPTION TO OUTPUT FILTER DIAGONALS TO FILE
  namelist / nml_debug / output_filter_diagonals
  ! OPTION TO OUTPUT PLOT3D GRID TO TECPLOT IN UNSTRUCTURED FORMAT
  namelist / nml_debug / output_plot3d_to_tecplot
  ! OPTION TO OUTPUT CONTENTS OF BFACE ARRAY
  namelist / nml_debug / output_bface_array
  ! DUMP FLUX INFORMATION FOR DEBUGGING
  namelist / nml_debug / check_flux_point_coordinates
  ! USE GENERIC N-DIMENSIONAL NON-OPTIMIZED SUBROUTINES
  namelist / nml_debug / force_nd_subroutines
  ! CHEMISTRY
  namelist / nml_debug / newton_convergence_goal, chemistry_sub_iterations, &
                         ignore_thermo_newton_warnings
  !
continue
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&nml_debug')
      call set('Output_Limiter_Statistics    ',                             &
                output_limiter_statistics     , fals                        )
      call set('Output_Filter_Diagonals      ',                             &
                output_filter_diagonals       , fals                        )
      call set('Output_Plot3d_To_Tecplot     ',                             &
                output_plot3d_to_tecplot      , fals                        )
      call set('Output_Bface_Array           ',                             &
                output_bface_array            , fals                        )
      call set('Check_Flux_Point_Coordinates ',                             &
                check_flux_point_coordinates  , fals                        )
      call set('Force_Nd_Subroutines         ',                             &
                force_nd_subroutines          , fals                        )
      call set('Newton_Convergence_Goal      ',                             &
                newton_convergence_goal       , eps12                       )
      call set('Chemistry_Sub_Iterations     ',                             &
                chemistry_sub_iterations      ,-1                           )
      call set('Ignore_Thermo_Newton_Warnings',                             &
                ignore_thermo_newton_warnings , fals                        )
      call iow_footer()
    case ( NMLaction_Read )
      if ( io_has_nml_block(ior,'&nml_debug') ) then
        read(ior,nml=nml_debug)
      end if
    case ( NMLaction_Write )
      write(iow,nml=nml_debug)
  end select
  !
end subroutine ionml_debug
!}}}
!###############################################################################
!{{{ subroutine ionml_testcase_vortex
!> I/O for namelist block 'nml_testcase_vortex'.
!###############################################################################
subroutine ionml_testcase_vortex()
  !
  !.. Use Statements ..
  use ovar, only : &
             vortex_bigR                                        , &
             vortex_bigGam                                      , &
             Grid_Lref            => VortexGrid_Lref            , &
             Grid_Random_Perturb  => VortexGrid_random_perturb  , &
             Grid_Perturb_Percent => VortexGrid_perturb_percent , &
             Grid_Hex_to_Pyramids => VortexGrid_Hex_to_Pyramids , &
             Grid_is_3D           => VortexGrid_is_3D
  !
  ! ISENTROPIC VORTEX CONSTANTS
  namelist / nml_testcase_vortex /  &
             vortex_bigR          , &
             vortex_bigGam        , &
             Grid_Lref
  ! TAYLOR-GREEN VORTEX CONSTANTS
  namelist / nml_testcase_vortex /  &
             Grid_Random_Perturb  , &
             Grid_Perturb_Percent , &
             Grid_Hex_to_Pyramids , &
             Grid_is_3D
  !
continue
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&nml_testcase_vortex')
      call set('Vortex_BigR               ',vortex_bigR         , 1.0_wp )
      call set('Vortex_BigGam             ',vortex_bigGam       , 5.0_wp )
      call set('Grid_Lref                 ',Grid_Lref           ,10.0_wp )
      call set('Grid_Random_Perturb       ',Grid_Random_Perturb , fals   )
      call set('Grid_Perturb_Percent      ',Grid_Perturb_Percent, 30     )
      call set('Grid_Hex_To_Pyramids      ',Grid_Hex_to_Pyramids, fals   )
      call set('Grid_Is_3D                ',Grid_is_3D          , fals   )
      call iow_footer()
    case ( NMLaction_Read )
      if ( io_has_nml_block(ior,'&nml_testcase_vortex') ) then
        read(ior,nml=nml_testcase_vortex)
      end if
    case ( NMLaction_Write )
      write(iow,nml=nml_testcase_vortex)
  end select
  !
end subroutine ionml_testcase_vortex
!}}}
!###############################################################################
!{{{ subroutine ionml_testcase_tcf
!> I/O for namelist block 'nml_testcase_tcf'.
!###############################################################################
subroutine ionml_testcase_tcf()
  !
  !.. Use Statements ..
  use ovar, only : &
             Inner_Radius       => TCF_Inner_Radius       , &
             Inner_Wall_Spacing => TCF_Inner_Wall_Spacing , &
             Outer_Radius       => TCF_Outer_Radius       , &
             Outer_Wall_Spacing => TCF_Outer_Wall_Spacing , &
             Aspect_Ratio       => TCF_Aspect_Ratio       , &
             Axial_Cells        => TCF_Axial_Cells        , &
             Radial_Cells       => TCF_Radial_Cells       , &
             Azimuthal_Cells    => TCF_Azimuthal_Cells    , &
             Theta_Shift        => TCF_Theta_Shift        , &
             Periodic_Sectors   => TCF_Periodic_Sectors
  !
  ! TAYLOR-COUETTE FLOW CONSTANTS
  namelist / nml_testcase_tcf / &
             Inner_Radius, Inner_Wall_Spacing, &
             Outer_Radius, Outer_Wall_Spacing, &
             Aspect_Ratio, Axial_Cells, &
             Radial_Cells, Azimuthal_Cells, &
             Theta_Shift,  Periodic_Sectors
  !
continue
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&nml_testcase_tcf')
      call set('Inner_Radius      ',Inner_Radius      , half       )
      call set('Inner_Wall_Spacing',Inner_Wall_Spacing, 0.002_wp   )
      call set('Outer_Radius      ',Outer_Radius      , one        )
      call set('Outer_Wall_Spacing',Outer_Wall_Spacing, 0.0025_wp  )
      call set('Aspect_Ratio      ',Aspect_Ratio      , two*pi     )
      call set('Axial_Cells       ',Axial_Cells       , 24         )
      call set('Radial_Cells      ',Radial_Cells      , 20         )
      call set('Azimuthal_Cells   ',Azimuthal_Cells   , 48         )
      call set('Theta_Shift       ',Theta_Shift       , zero       )
      call set('Periodic_Sectors  ',Periodic_Sectors  , 1          )
      call iow_footer()
    case ( NMLaction_Read )
      if ( io_has_nml_block(ior,'&nml_testcase_tcf') ) then
        read(ior,nml=nml_testcase_tcf)
      end if
    case ( NMLaction_Write )
      write(iow,nml=nml_testcase_tcf)
  end select
  !
end subroutine ionml_testcase_tcf
!}}}
!###############################################################################
!{{{ subroutine ionml_testcase_sod
!> I/O for namelist block 'nml_testcase_sod'.
!###############################################################################
subroutine ionml_testcase_sod()
  !
  !.. Use Statements ..
  use ovar, only : right_vector => sod_test_right_vector
  !
  ! SOD SHOCK TUBE TEST PARAMETERS
  namelist / nml_testcase_sod / right_vector
  !
continue
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&nml_testcase_sod')
      call set('Right_Vector',right_vector, [one,zero,zero]  )
      call iow_footer()
    case ( NMLaction_Read )
      if ( io_has_nml_block(ior,'&nml_testcase_sod') ) then
        read(ior,nml=nml_testcase_sod)
      end if
    case ( NMLaction_Write )
      write(iow,nml=nml_testcase_sod)
  end select
  !
end subroutine ionml_testcase_sod
!}}}
!###############################################################################
!{{{ subroutine ionml_testcase_pulse
!> I/O for namelist block 'nml_testcase_pulse'.
!###############################################################################
subroutine ionml_testcase_pulse()
  !
  !.. Use Statements ..
  use ovar,    only : form     => pulse_test_form      , &
                      form_dim => pulse_test_form_dim  , &
                      length   => pulse_test_length    , &
                      center   => pulse_test_center    , &
                      width    => pulse_test_width     , &
                      gap      => pulse_test_gap       , &
                      repeat   => pulse_test_repeat
  !
  ! DENSITY PULSE CONVECTION TEST PARAMETERS
  namelist / nml_testcase_pulse / form, form_dim, length, center, width, &
                                  gap, repeat
  !
continue
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&nml_testcase_pulse')
      call set('Form        ',form        , 0                )
      call set('Form_Dim    ',form_dim    , 1                )
      call set('Length      ',length      , 1.0_wp           )
      call set('Center      ',center      , [half,zero,zero] )
      call set('Width       ',width       , 0.2_wp           )
      call set('Gap         ',gap         , 0.2_wp           )
      call set('Repeat      ',repeat      , 0                )
      call iow_footer()
    case ( NMLaction_Read )
      if ( io_has_nml_block(ior,'&nml_testcase_pulse') ) then
        read(ior,nml=nml_testcase_pulse)
      end if
    case ( NMLaction_Write )
      write(iow,nml=nml_testcase_pulse)
  end select
  !
end subroutine ionml_testcase_pulse
!}}}
!###############################################################################
!{{{ subroutine ionml_testcase_channel
!> I/O for namelist block 'nml_testcase_channel'.
!###############################################################################
subroutine ionml_testcase_channel()
  !
  !.. Use Statements ..
  use ovar, only : annulus_angle       => pipe_channel_angle
  use ovar, only : body_force          => channel_body_force          , &
                   nondim_output_coord => channel_nondim_output_coord , &
                   grid_ave_xplus      => channel_grid_ave_xplus      , &
                   grid_ave_yplus      => channel_grid_ave_yplus      , &
                   grid_ave_zplus      => channel_grid_ave_zplus      , &
                   grid_yplus_wall     => channel_grid_yplus_wall     , &
                   grid_yplus_of_one   => channel_yplus_of_one
  !
  ! CHANNEL FLOW PARAMETERS
  namelist / nml_testcase_channel / &
    annulus_angle, &
    grid_ave_xplus, &
    grid_ave_yplus, &
    grid_ave_zplus, &
    grid_yplus_wall, &
    grid_yplus_of_one, &
    body_force, &
    nondim_output_coord
  !
continue
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&nml_testcase_channel')
      call set('Annulus_Angle      ',annulus_angle      , five              )
      call set('Grid_Ave_Xplus     ',grid_ave_xplus     ,-30.0_wp           )
      call set('Grid_Ave_Yplus     ',grid_ave_yplus     ,-15.0_wp           )
      call set('Grid_Ave_Zplus     ',grid_ave_zplus     ,-15.0_wp           )
      call set('Grid_Yplus_Wall    ',grid_yplus_wall    ,  2.0_wp           )
      call set('Grid_Yplus_Of_One  ',grid_yplus_of_one  , 8.42134e-6_wp     )
      call set('Body_Force         ',body_force         , 0.0_wp            )
      call set('Nondim_Output_Coord',nondim_output_coord, fals              )
      call iow_footer()
    case ( NMLaction_Read )
      if ( io_has_nml_block(ior,'&nml_testcase_channel') ) then
        read(ior,nml=nml_testcase_channel)
      end if
    case ( NMLaction_Write )
      write(iow,nml=nml_testcase_channel)
  end select
  !
end subroutine ionml_testcase_channel
!}}}
!###############################################################################
!{{{ subroutine ionml_bc
!> I/O for namelist block 'nml_bc'.
!###############################################################################
subroutine ionml_bc()
  !
  !.. Use Statements ..
  use ovar, only : walls_are_exact, sub_inflow_method, use_bc_conflicts, &
                   bc_input
  !
  !.. Local Scalars ..
  integer :: i, n
  character(len=:), allocatable :: text, text2
  !
  ! OPTIONS FOR MODIFYING BC ALGORITHMS
  namelist / nml_bc / walls_are_exact, sub_inflow_method, use_bc_conflicts
  ! INPUT BOUNDARY CONDITIONS
  namelist / nml_bc / bc_input
  !
continue
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&nml_bc')
      call set('Walls_Are_Exact             ', &
                walls_are_exact              , true )
      call set('Sub_Inflow_Method           ', &
                sub_inflow_method            , 4    )
      call set('Use_BC_Conflicts            ', &
                use_bc_conflicts             , fals )
      loop_bc_input : &
      do n = 1, size(bc_input)
      ! Must set text large enough to write into.
      text = '                               '
      write(text,101) 'BC_Input(',n,')'
      text = trim(text)   ! F2003 realloc
101   format(a,i0,a)
      call set( text //   '%Name                   ', &
                bc_Input(n)%name                    , '' )
      call set( text //   '%Set_BC_Value           ', &
                bc_input(n)%set_bc_value            , bc_unknown )
      call set( text //   '%Temperature            ', &
                bc_input(n)%Temperature             ,-one  )
      call set( text //   '%Pressure               ', &
                bc_input(n)%Pressure                ,-one  )
      call set( text //   '%Density                ', &
                bc_input(n)%Density                 ,-one  )
      call set( text //   '%Mach                   ', &
                bc_input(n)%mach                    ,-one  )
      call set( text //   '%VelocityX              ', &
                bc_input(n)%VelocityX               , zero )
      call set( text //   '%VelocityY              ', &
                bc_input(n)%VelocityY               , zero )
      call set( text //   '%VelocityZ              ', &
                bc_input(n)%VelocityZ               , zero )
      call set( text //   '%Total_Temperature      ', &
                bc_input(n)%Total_Temperature       ,-one  )
      call set( text //   '%Total_Pressure         ', &
                bc_input(n)%Total_Pressure          ,-one  )
      call set( text //   '%Total_Density          ', &
                bc_input(n)%Total_Density           ,-one  )
      call set( text //   '%AoA_Alpha              ', &
                bc_input(n)%aoa_alpha               , zero )
      call set( text //   '%AoA_Beta               ', &
                bc_input(n)%aoa_beta                , zero )
      call set( text //   '%Wall_Temp              ', &
                bc_input(n)%wall_temp               ,-one  )
      call set( text //   '%Rotation_Rate          ', &
                bc_input(n)%rotation_rate           , [zero,zero,zero] )
      ! Next item must have default logical type.
      call set( text //   '%Is_Specified_In_MRF'    , &
                bc_input(n)%is_specified_in_mrf     , .false. )
      !
      loop_mass_frac : &
      do i = 1, size(bc_input(n)%mass_frac)
        text2 = '                                '
        write(text2,101) text // '%Mass_Frac(',i,')'
        text2 = trim(text2)   ! F2003 realloc
        call set( text2 //                   '       ', &
                  bc_input(n)%mass_frac(i)            , zero  )
      end do loop_mass_frac
      !
      loop_species_name : &
      do i = 1, size(bc_input(n)%species_name)
        text2 = '                                '
        write(text2,101) text // '%Species_Name(',i,')'
        text2 = trim(text2)   ! F2003 realloc
        call set( text2 //                      '    ', &
                  bc_input(n)%species_name(i)         , '' )
      end do loop_species_name
      !
      call set( text //   '%Relax%Value            ', &
                bc_input(n)%relax%value             , one )
      call set( text //   '%Relax%Time             ', &
                bc_input(n)%relax%time              , zero )
      call set( text //   '%Relax%Iter             ', &
                bc_input(n)%relax%iter              , 0 )
      !
      call set( text //   '%SEM%Is_Cylindrical     ', &
                bc_input(n)%sem%is_cylindrical      , fals )
      call set( text //   '%SEM%Length_Scale       ', &
                bc_input(n)%sem%length_scale        , zero )
      call set( text //   '%SEM%Convective_Velocity', &
                bc_input(n)%sem%convective_velocity , zero)
      call set( text //   '%SEM%Profile_Height     ', &
                bc_input(n)%sem%profile_height      , zero )
      call set( text //   '%SEM%Profile_File       ', &
                bc_input(n)%sem%profile_file        , '' )
      call set( text //   '%SEM%Sub_In_BC_Grp      ', &
                bc_input(n)%sem%sub_in_bc_grp       , 0 )
      call set( text //   '%SEM%Sub_In_BC_Idx      ', &
                bc_input(n)%sem%sub_in_bc_idx       , 0 )
      !
      call set( text //   '%Periodic%Matching_BC_Name', &
                bc_input(n)%Periodic%Matching_BC_Name , '' )
      call set( text //   '%Periodic%Matching_BC_Idx ', &
                bc_input(n)%Periodic%Matching_BC_Idx  , 0 )
      !
      loop_rotation_center : &
      do i = 1, size(bc_input(n)%Periodic%RotationCenter)
        text2 = '          '
        write(text2,101) '(',i,')'
        text2 = text // '%Periodic%RotationCenter' // text2
        text2 = trim(text2)   ! F2003 realloc
        call set( text2 //                         ' ', &
                  bc_input(n)%Periodic%RotationCenter(i), zero  )
      end do loop_rotation_center
      !
      loop_rotation_angle : &
      do i = 1, size(bc_input(n)%Periodic%RotationAngle)
        text2 = '          '
        write(text2,101) '(',i,')'
        text2 = text // '%Periodic%RotationAngle' // text2
        text2 = trim(text2)   ! F2003 realloc
        call set( text2 //                        '  ', &
                  bc_input(n)%Periodic%RotationAngle(i), zero  )
      end do loop_rotation_angle
      !
      loop_rotation_axis : &
      do i = 1, size(bc_input(n)%Periodic%RotationAxis)
        text2 = '          '      
        write(text2,101) '(',i,')'
        text2 = text // '%Periodic%RotationAxis' // text2
        text2 = trim(text2)   ! F2003 realloc
        call set( text2 //                       '   ', &
                  bc_input(n)%Periodic%RotationAxis(i), zero  )
      end do loop_rotation_axis
      !
      loop_translation : &
      do i = 1, size(bc_input(n)%Periodic%Translation)
        text2 = '          '
        write(text2,101) '(',i,')'
        text2 = text // '%Periodic%Translation' // text2
        text2 = trim(text2)   ! F2003 realloc
        call set( text2 //                      '    ', &
                  bc_input(n)%Periodic%Translation(i), zero  )
      end do loop_translation
      !
      end do loop_bc_input  ! n = 1, size(bc_input)
      call iow_footer()
    case ( NMLaction_Read )
      if ( io_has_nml_block(ior,'&nml_bc') ) then
        read(ior,nml=nml_bc)
      end if
    case ( NMLaction_Write )
      write(iow,nml=nml_bc)
  end select
  !
end subroutine ionml_bc
!}}}
!###############################################################################
!{{{ subroutine ionml_time_ave
!> I/O for namelist block 'nml_time_ave'.
!###############################################################################
subroutine ionml_time_ave()
  !
  !.. Use Statements ..
  use ovar,    only : compute_time_averaging => output_time_averaging, &
                      num_saved_time_ave_restarts                    , &
                      time_ave_file                                  , &
                      time_ave_vel_is_axisymm                        , &
                      start_averaging_timestep                       , &
                      start_averaging_tstar
  !
  !.. Local Scalars ..
  integer :: n
  character(len=:), allocatable :: text
  !
  ! TIME AVERAGING
  namelist / nml_time_ave / compute_time_averaging
  namelist / nml_time_ave / num_saved_time_ave_restarts
  namelist / nml_time_ave / time_ave_file, &
                            time_ave_vel_is_axisymm, &
                            start_averaging_timestep, start_averaging_tstar
  !
continue
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&nml_time_ave')
      call set('Compute_Time_Averaging            ',                        &
                compute_time_averaging             , fals                   )
      call set('Num_Saved_Time_Ave_Restarts       ',                        &
                num_saved_time_ave_restarts        , 2                      )
      call set('Time_Ave_File                     ',                        &
                time_Ave_File                      , 'out.ave'              )
      call set('Time_Ave_Vel_Is_Axisymm           ',                        &
                time_ave_vel_is_axisymm            , fals                   )
      call set('Start_Averaging_Timestep          ',                        &
                start_averaging_timestep           ,-1                      )
      call set('Start_Averaging_Tstar             ',                        &
                start_averaging_tstar              ,-one                    )
      call iow_footer()
    case ( NMLaction_Read )
      if ( io_has_nml_block(ior,'&nml_time_ave') ) then
        read(ior,nml=nml_time_ave)
      end if
    case ( NMLaction_Write )
      write(iow,nml=nml_time_ave)
  end select
  !
end subroutine ionml_time_ave
!}}}
!###############################################################################
!{{{ subroutine ionml_utility
!> I/O for namelist block 'nml_utility'.
!###############################################################################
subroutine ionml_utility()
  !
  !.. Use Statements ..
  use ovar, only : convert_restart_to_cgns
  use ovar, only : delete_converted_restart_files
  use ovar, only : restart_list_file
  use ovar, only : spatial_ave_restart_channel => channel_spatial_ave_restart
  use ovar, only : spatial_ave_restart_x_and_z_planes
  use ovar, only : spatial_ave_restart_z_planes
  use ovar, only : modify_time_ave_restart_history
  ! UTILITY TO PROCESS EXISTING RESTART FILES
  namelist / nml_utility / &
    convert_restart_to_cgns, &
    delete_converted_restart_files, &
    restart_list_file, &
    spatial_ave_restart_channel, &
    spatial_ave_restart_x_and_z_planes, &
    spatial_ave_restart_z_planes, &
    modify_time_ave_restart_history
  !
continue
  !
    select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&nml_utility')
      call set('Convert_Restart_To_CGNS           ',                        &
                convert_restart_to_cgns            , fals                   )
      call set('Delete_Converted_Restart_Files    ',                        &
                delete_converted_restart_files     , fals                   )
      call set('Restart_List_File                 ',                        &
                restart_list_file                  , ''                     )
      call set('Spatial_Ave_Restart_Channel       ',                        &
                spatial_ave_restart_channel        , fals                   )
      call set('Spatial_Ave_Restart_X_And_Z_Planes',                        &
                spatial_ave_restart_x_and_z_planes , fals                   )
      call set('Spatial_Ave_Restart_Z_Planes      ',                        &
                spatial_ave_restart_z_planes       , fals                   )
      call set('Modify_Time_Ave_Restart_History   ',                        &
                modify_time_ave_restart_history    , fals                   )
      call iow_footer()
    case ( NMLaction_Read )
      if ( io_has_nml_block(ior,'&nml_utility') ) then
        read(ior,nml=nml_utility)
      end if
    case ( NMLaction_Write )
      write(iow,nml=nml_utility)
  end select
!
end subroutine ionml_utility
!}}}
!###############################################################################
!{{{ subroutine ionml_restart_interupt
!> I/O for namelist block 'nml_restart_interupt'.
!###############################################################################
subroutine ionml_restart_interupt()
  !
  !.. Use Statements ..
  use ovar, only : restart_file
  !
  ! RESTART INTERUPTION NAMELIST
  !
  namelist / nml_restart_interupt / restart_file
  !
continue
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
!---  Init of restart_file is done in nml_initialzation.
!     call iow_header('&nml_restart_interupt')
!     call set('restart_file',restart_file, 'out.rst')
!     call iow_footer()
    case ( NMLaction_Read )
      if ( io_has_nml_block(ior,'&nml_restart_interupt') ) then
        read(ior,nml=nml_restart_interupt)
      end if
    case ( NMLaction_Write )
      write(iow,nml=nml_restart_interupt)
  end select
  !
end subroutine ionml_restart_interupt
!}}}
!###############################################################################
!{{{ subroutine ionml_unknown
!> Namelist inputs that have no home.
!###############################################################################
subroutine ionml_unknown()
  !
  !.. Use Statements ..
  !
continue
  !
! select case (NMLaction)
! end select
  !
end subroutine ionml_unknown
!}}}
!###############################################################################
!{{{ subroutine ionml_removed
!> Namelist inputs that have been removed as inputs.
!###############################################################################
subroutine ionml_removed()
  !
  !.. Use Statements ..
  !
  use ovar, only : loc_init_pts
  use ovar, only : CFL_reduction_for_NaN_restart
  use ovar, only : ideal_total_temperature
  use ovar, only : limit_total_temperature
  use ovar, only : Pyramid_Filter_Option
  use ovar, only : Limiter_Marker_c1
  use ovar, only : Limiter_Marker_c2
  use ovar, only : flux_divergence
  use ovar, only : restart_output_uses_host_roots
  use ovar, only : use_async_for_error
  use ovar, only : output_eperk_cfl
  use ovar, only : output_bnd_profiles
  use ovar, only : parametric_agglomeration
  use ovar, only : interpto_all_orders
  use ovar, only : interpto_projection_order
  use ovar, only : interpto_quadrature_order
  use ovar, only : interpto_output_order
  use ovar, only : interpto_error_order
  use ovar, only : Prt
  use ovar, only : schmidt_turb
  use ovar, only : tinref
  use ovar, only : tvrref
  use ovar, only : turbulence_model
  use ovar, only : traceback_unconverged_newton_temps
  use ovar, only : use_unrolled_dot_products
  use ovar, only : use_improved_SEM_normalization
  use ovar, only : average_sem_inflow_plane
  use ovar, only : time_average_restart_files
  use ovar, only : postproc_debug_grid
  use ovar, only : dump_memory_usage
  use ovar, only : dump_fluxes
  use ovar, only : disable_header_rewrites
  !
  ! These inputs are no longer used by the main code, so define them here.
  ! Must save them for init/read/write/convert passes to work properly if
  ! nml_removed is included in the namelist output.
  real(wp), save :: newton_rlx                 = 1.0_wp
  integer,  save :: VortexGrid_DOF             = 48
  integer,  save :: internal_grid_dof          = 0
  integer,  save :: flame_speed_test_dof       = 0
  integer,  save :: sod_test_dof               = 0
  integer,  save :: pulse_test_dof             = 0
  integer,  save :: shu_osher_test_dof         = 0
  !
  ! PLOT3D COARSENING OPTIONS
  namelist / nml_removed / parametric_agglomeration
  ! INTERNAL GRID DOF, now uses internal_grid_dof
  namelist / nml_removed / flame_speed_test_dof
  namelist / nml_removed / sod_test_dof
  namelist / nml_removed / pulse_test_dof
  namelist / nml_removed / shu_osher_test_dof
  namelist / nml_removed / VortexGrid_DOF
  ! INTERPOLATION UTILITY SOLUTION ORDER
  namelist / nml_removed / interpto_all_orders
  namelist / nml_removed / interpto_projection_order
  namelist / nml_removed / interpto_quadrature_order
  namelist / nml_removed / interpto_output_order
  namelist / nml_removed / interpto_error_order
  ! REFERENCE TURBULENCE CONDITIONS
  namelist / nml_removed / Prt
  namelist / nml_removed / schmidt_turb
  namelist / nml_removed / tinref, tvrref
  ! GOVERNING EQUATIONS
  namelist / nml_removed / turbulence_model
  ! CHEMISTRY
  namelist / nml_removed / traceback_unconverged_newton_temps
  ! CONTROL DOT_PRODUCT PROCEDURES BOUND TO MATRIX DERIVED-TYPE
  namelist / nml_removed / use_unrolled_dot_products
  ! SEM OPTIONS
  namelist / nml_removed / use_improved_SEM_normalization
  namelist / nml_removed / average_sem_inflow_plane
  ! 
  ! LOCATION OF INITIALIZATION POINTS
  namelist / nml_removed / loc_init_pts
  namelist / nml_removed / CFL_reduction_for_NaN_restart
  ! RESTART OPTIONS
  namelist / nml_removed / restart_output_uses_host_roots
  ! BASE FR OPTIONS
  namelist / nml_removed / flux_divergence
  ! REALIZABLE TOTAL TEMPERATURE LIMITER
  namelist / nml_removed / ideal_total_temperature, limit_total_temperature
  ! FILTER OPTIONS
  namelist / nml_removed / Pyramid_Filter_Option
  namelist / nml_removed / Limiter_Marker_c1
  namelist / nml_removed / Limiter_Marker_c2
  ! DISPLAY OPTIONS
  namelist / nml_removed / output_eperk_cfl
  ! OPTIONS TO OUTPUT BOUNDARY PROFILE
  namelist / nml_removed / output_bnd_profiles
  namelist / nml_removed / disable_header_rewrites
  ! OPENACC INPUT OPTIONS
  namelist / nml_removed / use_async_for_error
  ! UTILITY TO TIME AVERAGE EXISTING RESTART FILES
  namelist / nml_removed / time_average_restart_files
  ! ACTIVATE SPECIAL GRID FOR DEBUGGING THE POSTPROCESSOR
  namelist / nml_removed / postproc_debug_grid
  ! MEMORY USAGE
  namelist / nml_removed / dump_memory_usage
  ! MISC
  namelist / nml_removed / dump_fluxes
  namelist / nml_removed / newton_rlx
  !
  !
continue
  !
  ! Enabling the following line will omit output of this namelist block.
  ! The output is only intended for developers to validate old file conversion.
  ! Users should not include removed inputs in a new namelist file.
  ! Developers may include removed inputs in a new namelist file to 
  !   activate certain debugging options.
  ! 
! if ( NMLaction == NMLaction_Write .or. &
!      NMLaction == NMLaction_WritePretty ) return
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&nml_removed')
      ! Grid
      call set('Parametric_Agglomeration         ',                         &
                parametric_agglomeration          ,-1                       )
      call set('Flame_Speed_Test_DoF             ',                         &
                flame_speed_test_dof              , 0                       )
      call set('Sod_Test_DoF                     ',                         &
                sod_test_dof                      , 0                       )
      call set('Pulse_Test_DoF                   ',                         &
                pulse_test_dof                    , 0                       )
      call set('Shu_Osher_Test_DoF               ',                         &
                shu_osher_test_dof                , 0                       )
      call set('VortexGrid_DoF                   ',                         &
                VortexGrid_DOF                    , 48                      )
      ! Interpolation
      call set('Interpto_All_Orders              ',                         &
                interpto_all_orders               ,-3                       )
      call set('Interpto_Projection_Order        ',                         &
                interpto_projection_order         ,-3                       )
      call set('Interpto_Quadrature_Order        ',                         &
                interpto_quadrature_order         ,-3                       )
      call set('Interpto_Output_Order            ',                         &
                interpto_output_order             ,-3                       )
      call set('Interpto_Error_Order             ',                         &
                interpto_error_order              ,-3                       )
      ! Reference
      call set('Prt                              ',                         &
                Prt                               , 0.72_wp                 )
      call set('Schmidt_Turb                     ',                         &
                schmidt_turb                      , 0.72_wp                 )
      call set('Tinref                           ',                         &
                Tinref                            , 0.0_wp                  )
      call set('Tvrref                           ',                         &
                Tvrref                            , 1.0e-6_wp               )
      ! Physics
      call set('Turbulence_Model                 ',                         &
                turbulence_model                  , laminar_flow            )
      ! Chemistry
      call set('Traceback_Unconverged_Newton_Temps',                        &
                traceback_unconverged_newton_temps , fals                   )
      call set('Use_Unrolled_Dot_Products        ',                         &
                use_unrolled_dot_products         , true                    )
      ! SEM
      call set('Use_Improved_SEM_Normalization   ',                         &
                use_improved_SEM_normalization    , true                    )
      call set('Average_SEM_Inflow_Plane         ',                         &
                average_sem_inflow_plane          , fals                    )
      ! Initialization
      call set('Loc_Init_Pts                     ',                         &
                loc_init_pts                      , 0                       )
      call set('CFL_Reduction_For_NaN_Restart    ',                         &
                CFL_reduction_for_NaN_restart     , 2.0_wp                  ) 
      call set('Restart_Output_Uses_Host_Roots   ',                         &
                restart_output_uses_host_roots    , true                    ) 
      ! Numerics
      call set('Flux_Divergence                  ',                         &
                flux_divergence                   , LagrangePolynomial      )
      ! Limiters
      call set('Ideal_Total_Temperature          ',                         &
                ideal_total_temperature           ,-1.0_wp                  )
      call set('Limit_Total_Temperature          ',                         &
                limit_total_temperature           , 0.2_wp                  )
      call set('Pyramid_Filter_Option            ',                         &
                Pyramid_Filter_Option             , 1                       )
      call set('Limiter_Marker_C1                ',                         &
                Limiter_Marker_c1                 ,-1.0_wp                  )
      call set('Limiter_Marker_C2                ',                         &
                Limiter_Marker_c2                 ,-1.0_wp                  )
      ! Output
      call set('Output_EPERK_CFL                 ',                         &
                output_eperk_cfl                  , true                    )
      call set('Output_Bnd_Profiles              ',                         &
                output_bnd_profiles               , fals                    )
      call set('Disable_Header_Rewrites          ',                         &
                disable_header_rewrites           , fals                    )
      call set('Use_Async_For_Error              ',                         &
                use_async_for_error               , fals                    )
      ! Utility
      call set('Time_Average_Restart_Files       ',                         &
                time_average_restart_files        , fals                    )
      ! Debugging
      call set('Postproc_Debug_Grid              ',                       &
                postproc_debug_grid               , fals                  )
      call set('Dump_Memory_Usage                ',                         &
                dump_memory_usage                 , fals                    )
      call set('Dump_Fluxes                      ',                         &
                dump_fluxes                       , fals                    )
      ! Misc
      call set('Newton_Rlx                       ',                         &
                newton_rlx                        , 1.0_wp                  )
      call iow_footer()
    case ( NMLaction_Read )
      if ( io_has_nml_block(ior,'&nml_removed') ) then
        read(ior,nml=nml_removed)
      end if
    case ( NMLaction_Write )
      write(iow,nml=nml_removed)
  end select
  !
end subroutine ionml_removed
!}}}
!###############################################################################
!
end module nml_0002_mod
