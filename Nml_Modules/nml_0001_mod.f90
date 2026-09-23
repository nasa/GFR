!###############################################################################
!> Module for version 0001 of namelist io.
!###############################################################################
module nml_0001_mod
  !
  !.. Use Statements ..
  use module_kind_types
  !
  use bc_input_0001_mod, only : bc_input,bc_input_t
  use bc_input_0001_mod, only : bc_input_v0001_to_ovar
  use bc_input_0001_mod, only : bc_input_ovar_to_v0001
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
  !{{{ ovar use statements
  use ovar, only : grid_format
  use ovar, only : gridfile
  use ovar, only : plot3d_all_strides
  use ovar, only : plot3d_i_stride
  use ovar, only : plot3d_j_stride
  use ovar, only : plot3d_k_stride
  use ovar, only : plot3d_agglomerate_order
  use ovar, only : parametric_agglomeration
  use ovar, only : plot3d_cutfile
  use ovar, only : plot3d_bcsfile
  use ovar, only : grid_scaling_factor
  use ovar, only : internal_grid_dof
  use ovar, only : internal_grid_order
  use ovar, only : interpto_grid_format
  use ovar, only : interpto_gridfile
  use ovar, only : interpto_plot3d_all_strides
  use ovar, only : interpto_plot3d_i_stride
  use ovar, only : interpto_plot3d_j_stride
  use ovar, only : interpto_plot3d_k_stride
  use ovar, only : interpto_plot3d_agglomerate_order
  use ovar, only : interpto_plot3d_cutfile
  use ovar, only : interpto_plot3d_bcsfile
  use ovar, only : interpto_grid_scaling_factor
  use ovar, only : solution_format
  use ovar, only : solution_file
  use ovar, only : loc_solution_pts
  use ovar, only : loc_flux_pts
  use ovar, only : loc_triangle_pts
  use ovar, only : loc_tetrahedra_pts
  use ovar, only : loc_quadrature_pts
  use ovar, only : combine_common_flux_projections
  use ovar, only : flux_divergence
  use ovar, only : correction_function
  use ovar, only : Runge_Kutta_Scheme
  use ovar, only : num_rk_stages
  use ovar, only : invs_flux_method
  use ovar, only : visc_flux_method
  use ovar, only : visc_variables_to_ave
  use ovar, only : walls_are_exact
  use ovar, only : num_timesteps
  use ovar, only : Final_Time
  use ovar, only : constant_dt
  use ovar, only : cfl
  use ovar, only : cfl_beg
  use ovar, only : cfl_end
  use ovar, only : cfl_cycles
  use ovar, only : cfl_start_iter
  use ovar, only : Timestep_Type, minimum_timesteps
  use ovar, only : all_points, metis_option, use_parmetis
  use ovar, only : itestcase
  use ovar, only : jet_init_func
  use ovar, only : governing_equations
  use ovar, only : moving_ref_frame
  use ovar, only : disable_chemical_reactions, semi_implicit_chemistry
  use ovar, only : traceback_unconverged_newton_temps
  use ovar, only : ignition_method
  use ovar, only : ignition_start_timestep
  use ovar, only : ignition_duration
  use ovar, only : ignition_prod2reac_L2_ratio_goal
  use ovar, only : ignition_source_magnitude
  use ovar, only : ignition_max_temperature
  use ovar, only : ignition_spark_cap
  use ovar, only : ignition_source_radius
  use ovar, only : ignition_source_location
  use ovar, only : load_restart_file, restart_interval, restart_file
  use ovar, only : restart_output_uses_host_roots
  use ovar, only : CFL_reduction_for_NaN_restart
  use ovar, only : num_saved_time_ave_restarts
  use ovar, only : restart_frequency, completely_disable_restart_output
  use ovar, only : iter_out_interval
  use ovar, only : results_interval
  use ovar, only : output_eperk_cfl, disable_header_rewrites
  use ovar, only : output_summary_file
  use ovar, only : output_dir, output_interval, output_frequency
  use ovar, only : continuous_output
  use ovar, only : postproc_debug_grid
  use ovar, only : output_bnd_profiles,sub_inflow_method
  use ovar, only : machref
  use ovar, only : reyref
  use ovar, only : gam
  use ovar, only : Pr
  use ovar, only : Prt
  use ovar, only : schmidt_lam
  use ovar, only : schmidt_turb
  use ovar, only : rhoref
  use ovar, only : tref
  use ovar, only : pref
  use ovar, only : rgasref
  use ovar, only : turbulence_model
  use ovar, only : tinref
  use ovar, only : tvrref
  use ovar, only : alpha_aoaref => aoa_alpharef
  use ovar, only : beta_aoaref  => aoa_betaref
  use ovar, only : profile_io_all,profile_io_restart
  use ovar, only : profile_io_cgns,profile_io_cgns_subs
  use ovar, only : suth_muref, suth_Tref, suth_Sref
  use ovar, only : Filter_Option
  use ovar, only : Limiter_Filter_Option
  use ovar, only : Filter_etac
  use ovar, only : Limiter_Filter_etac
  use ovar, only : Filter_Order
  use ovar, only : Limiter_Filter_Order
  use ovar, only : Filter_Alpha
  use ovar, only : Limiter_Filter_Alpha
  use ovar, only : Limiter_Use_Dynamic_Filter
  use ovar, only : Track_Limiter_Filter_Dissipation
  use ovar, only : Pyramid_Filter_Option
  use ovar, only : Sponge_Filter_Option
  use ovar, only : Sponge_Filter_etac
  use ovar, only : Sponge_Filter_Order
  use ovar, only : Sponge_Filter_Alpha
  use ovar, only : Sponge_Beyond_X
  use ovar, only : Limiter_Option
  use ovar, only : Marker_Option
  use ovar, only : Marker_Variable
  use ovar, only : Limiter_Marker_Min, Limiter_Marker_Max
  use ovar, only : Limiter_Marker_dpsi
  use ovar, only : Limiter_Marker_c1, Limiter_Marker_c2
  use ovar, only : Filtering_Interval       => Global_Filter_Interval
  use ovar, only :                             Sponge_Filter_Interval
  use ovar, only : Limiting_Interval        => Limiter_Interval
  use ovar, only : Limiting_Interval_Output => Marker_Stored_Interval
  use ovar, only : thermofix
  use ovar, only : loc_output_pts, loc_init_pts, loc_error_pts
  use ovar, only : vortex_bigR, vortex_bigGam
  use ovar, only : output_plot3d_to_tecplot, output_bface_array
  use ovar, only : interpolate_before_output_variable
  use ovar, only : cgns_high_order_output,cgns_mixed_output
  use ovar, only : cgns_output_restart_time_ave_data
  use ovar, only : correct_cgns_grid_points
  use ovar, only : force_cgns_linear_grid_points
  use ovar, only : convergence_order_abs_res, convergence_order_max_res
  use ovar, only : dump_memory_usage,use_bc_conflicts
  use ovar, only : dump_fluxes,check_flux_point_coordinates
  use ovar, only : use_unrolled_dot_products
  use ovar, only : VortexGrid_random_perturb
  use ovar, only : VortexGrid_perturb_percent,VortexGrid_Lref
  use ovar, only : VortexGrid_Hex_to_Pyramids,VortexGrid_is_3D
  use ovar, only : TCF_Inner_Radius,TCF_Inner_Wall_Spacing
  use ovar, only : TCF_Outer_Radius,TCF_Outer_Wall_Spacing
  use ovar, only : TCF_Aspect_Ratio,TCF_Axial_Cells
  use ovar, only : TCF_Radial_Cells,TCF_Azimuthal_Cells
  use ovar, only : TCF_Theta_Shift, TCF_Periodic_Sectors
  use ovar, only : sod_test_right_vector
  use ovar, only : pulse_test_form, pulse_test_form_dim
  use ovar, only : pulse_test_length, pulse_test_center
  use ovar, only : pulse_test_width, pulse_test_gap, pulse_test_repeat
  use ovar, only : channel_body_force
  use ovar, only : channel_spatial_ave_restart
  use ovar, only : channel_nondim_output_coord
  use ovar, only : channel_yplus_of_one
  use ovar, only : channel_grid_ave_xplus, channel_grid_ave_zplus
  use ovar, only : channel_grid_ave_yplus, channel_grid_yplus_wall
  use ovar, only : pipe_channel_angle
  use ovar, only : convert_restart_to_cgns, completely_disable_cgns
  use ovar, only : delete_converted_restart_files
  use ovar, only : lustre_stripe_count, lustre_stripe_size
  use ovar, only : ovar_bc_input => bc_input
 !use ovar, only : ovar_bc_input_t => bc_input_t
  use ovar, only : geom_wts!,geom_wts_t
  use ovar, only : probes_collection_interval, probes_write_interval
  use ovar, only : point_probes, rake_probes, cell_probes
  use ovar, only : output_time_averaging
  use ovar, only : time_ave_file
  use ovar, only : start_averaging_timestep, start_averaging_tstar
  use ovar, only : time_ave_vel_is_axisymm
  use ovar, only : time_scaling_factor
  use ovar, only : spatial_ave_restart_x_and_z_planes
  use ovar, only : spatial_ave_restart_z_planes
  use ovar, only : time_average_restart_files, restart_list_file
  use ovar, only : modify_time_ave_restart_history
  use ovar, only : use_async_for_error
  use ovar, only : ideal_total_temperature,limit_total_temperature
  use ovar, only : output_limiter_statistics
  use ovar, only : output_filter_diagonals
  use ovar, only : output_vorticity
  use ovar, only : output_cylindrical_coordinates
  use ovar, only : output_cylindrical_velocity
  use ovar, only : output_relative_velocity
  use ovar, only : output_relative_conditions
  use ovar, only : output_qcriterion
  use ovar, only : output_conserved_variables
  use ovar, only : output_momentum
  use ovar, only : output_temperature
  use ovar, only : output_thermo_variables
  use ovar, only : output_gradients
  use ovar, only : output_density_gradient
  use ovar, only : output_density_gradient_magnitude
  use ovar, only : output_pressure_gradient
  use ovar, only : output_pressure_gradient_magnitude
  use ovar, only : output_pressure_dmat_gradient
  use ovar, only : output_temperature_gradient
  use ovar, only : output_temperature_gradient_magnitude
  use ovar, only : output_temperature_dmat_gradient
  use ovar, only : output_schlieren
  use ovar, only : output_modified_schlieren
  use ovar, only : output_transport_variables
  use ovar, only : output_stagnation_temperature
  use ovar, only : output_stagnation_pressure
  use ovar, only : output_stagnation_variables
  use ovar, only : output_marker_variables
  use ovar, only : output_grid_details
  use ovar, only : output_local_timestep
  use ovar, only : output_species_density
  use ovar, only : output_species_mass_fractions
  use ovar, only : output_species_mole_fractions
  use ovar, only : output_extra_species_variables
  use ovar, only : output_debug_species_variables
  use ovar, only : output_grdnt_species_variables
  use ovar, only : output_cellave_density
  use ovar, only : output_exact_density
  use ovar, only : output_exact_pressure
  use ovar, only : use_improved_SEM_normalization
  use ovar, only : use_SEM_profile_for_initialization
  use ovar, only : SEM_eddy_volume_is_spherical
  use ovar, only : SEM_box_num_eddies_multiplier
  use ovar, only : SEM_use_random_seed_from_restart
  use ovar, only : SEM_include_thermo_fluctuations
  use ovar, only : SEM_use_restart_eddies
  use ovar, only : average_sem_inflow_plane
  use ovar, only : force_nd_subroutines
  use ovar, only : conv_temp
  use ovar, only : chemistry_sub_iterations
  use ovar, only : ignore_thermo_newton_warnings
  use ovar, only : all_orders
  use ovar, only : solution_order
  use ovar, only : projection_order
  use ovar, only : quadrature_order
  use ovar, only : output_order
  use ovar, only : error_order
  use ovar, only : interpto_all_orders
  use ovar, only : interpto_solution_order
  use ovar, only : interpto_projection_order
  use ovar, only : interpto_quadrature_order
  use ovar, only : interpto_output_order
  use ovar, only : interpto_error_order
  use ovar, only : dmat_grad_temp_order
  !}}}
  !
  implicit none
  !
  private
  !
  !
  ! ###########################
  ! ###  MODULE PROCEDURES  ###
  ! ###########################
  !
  private :: alloc_nml_data
  private :: resiz_nml_data
  private ::   set_nml_data
  private ::  read_nml_data
  private :: write_nml_data
  private ::  read_nml_file
  private :: write_nml_file
  !
  private :: ionml
  private :: ionml_input
  !
  !
  ! ###########################
  ! ###  MODULE PARAMETERS  ###
  ! ###########################
  !
  integer, public, parameter :: module_namelist_version = 0001
  !
  !
  ! ##############################
  ! ###  MODULE DERIVED TYPES  ###
  ! ##############################
  !
  type, public, extends(namelist_procedures_t) :: nml_0001_procedures_t
  contains
    procedure, nopass :: init => init_nml0001_procedures
   !procedure :: convert => convert_nml0001_file
  end type nml_0001_procedures_t
  !
  !
  ! #################################
  ! ###  MODULE PUBLIC VARIABLES  ###
  ! #################################
  !
  type(nml_0001_procedures_t), public, protected, save :: nml_ver_0001
  !
  !
  ! ################################
  ! ###  MODULE LOCAL VARIABLES  ###
  ! ################################
  !
  ! Obsolete inputs defined locally in order to include in old namelist read.
  ! Explicitly save these, so their values will be retained.
  integer,     save :: Resolution_Marker_Variable = 0
  integer,     save :: Jump_Marker_Variable       = 0
  integer,     save :: MLP_Marker_Variable        = 0
  real(wp),    save :: newton_rlx                 = one
  logical(lk), save :: Use_Sponge_Filter          = fals
  integer,     save :: VortexGrid_DOF             = 48
  integer,     save :: flame_speed_test_DOF       = 0
  integer,     save :: sod_test_DOF               = 0
  integer,     save :: pulse_test_DOF             = 0
  integer,     save :: shu_osher_test_DOF         = 0
  !
contains
!
!###############################################################################
!
subroutine init_nml0001_procedures()
  !
  !.. Use Statements ..
  use nml_9999_mod, only : set_latest_nml_version
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_nml0001_procedures"
  !
continue
  !
  nml_ver_0001%version = module_namelist_version
  !
  nml_ver_0001%alloc_data => alloc_nml_data
  nml_ver_0001%resiz_data => resiz_nml_data
  nml_ver_0001%  set_data =>   set_nml_data
  nml_ver_0001% read_data =>  read_nml_data
  nml_ver_0001%write_data => write_nml_data
  nml_ver_0001% read_file =>  read_nml_file
  nml_ver_0001%write_file => write_nml_file
  !
  ! Send this namelist version into the subroutine set_latest_nml_version
  ! to set up the nml_ver_9999 procedure pointers.
  ! NOTE: This will only set the nml_ver_9999 procedure pointers
  !       to the respective procedures from this module if the
  !       namelist version of this module matches the value of
  !       latest_namelist_version from the module nml_base_mod.
  !
  call set_latest_nml_version(module_namelist_version,nml_ver_0001)
  !
end subroutine init_nml0001_procedures
!
!###############################################################################
!> Allocate namelist variables that are allocatable arrays or character scalars.
!###############################################################################
subroutine alloc_nml_data(opt_n_bc_input,opt_char_len, &
                          opt_n_species,opt_n_point_probes, &
                          opt_n_rake_probes,opt_n_cell_probes)
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: opt_n_bc_input
  integer, optional, intent(in) :: opt_char_len
  integer, optional, intent(in) :: opt_n_species
  integer, optional, intent(in) :: opt_n_point_probes
  integer, optional, intent(in) :: opt_n_rake_probes
  integer, optional, intent(in) :: opt_n_cell_probes
  !
  !.. Local Scalars ..
  integer :: n,ierr
  integer :: default_len,n_species,n_bc_input
  integer :: n_point_probes,n_rake_probes,n_cell_probes
  character(len=50) :: array_name
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
continue
  !
  n_bc_input = max_n_bc_input
  n_species = max_n_species
  n_point_probes = max_n_point_probes
  n_rake_probes = max_n_rake_probes
  n_cell_probes = max_n_cell_probes
  default_len = max_default_len
  !
  if (present(opt_n_bc_input)) n_bc_input = opt_n_bc_input
  if (present(opt_n_species)) n_species = opt_n_species
  if (present(opt_n_point_probes)) n_point_probes = opt_n_point_probes
  if (present(opt_n_rake_probes)) n_rake_probes = opt_n_rake_probes
  if (present(opt_n_cell_probes)) n_cell_probes = opt_n_cell_probes
  if (present(opt_char_len)) default_len = opt_char_len
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
  ! Copy bc_input into ovar_bc_input
  !
  ovar_bc_input = bc_input_v0001_to_ovar( bc_input ) ! F2003 AUTO-REALLOC
  !
  ! ---------------------------------------------------
  ! Allocate point_probes, rake_probes, and cell_probes
  ! ---------------------------------------------------
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
  ! ------------------------------------------------------------------------
  ! Set allocatable character variables to their default values and lengths.
  ! ------------------------------------------------------------------------
  !
  call set('Gridfile               ',                         &
            gridfile                , ''        , default_len )
  call set('Plot3d_Cutfile         ',                         &
            plot3d_cutfile          , ''        , default_len )
  call set('Plot3d_BCSfile         ',                         &
            plot3d_bcsfile          , ''        , default_len )
  call set('Interpto_Gridfile      ',                         &
            interpto_gridfile       , ''        , default_len )
  call set('Interpto_Plot3d_Cutfile',                         &
            interpto_plot3d_cutfile , ''        , default_len )
  call set('Interpto_Plot3d_BCSfile',                         &
            interpto_plot3d_bcsfile , ''        , default_len )
  call set('Restart_File           ',                         &
            restart_file            , 'out.rst' , default_len )
  call set('Solution_File          ',                         &
            solution_file           , ''        , default_len )
  call set('Output_Dir             ',                         &
            output_dir              , '.'       , default_len )
  call set('Lustre_Stripe_Size     ',                         &
            lustre_stripe_size      , '4m'      , default_len )
  call set('Time_Ave_File          ',                         &
            time_ave_file           , 'out.ave' , default_len )
  call set('Restart_List_File      ',                         &
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
  use ovar, only : point_probe_t, rake_probe_t
  !
  !.. Local Scalars ..
  integer :: i,n,ierr
  !
  !.. Local Allocatable Arrays ..
  type(bc_input_t),    allocatable :: tmp_bc_input(:)
  type(point_probe_t), allocatable :: tmp_point_probes(:)
  type(rake_probe_t),  allocatable :: tmp_rake_probes(:)
  type(point_probe_t), allocatable :: tmp_cell_probes(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "resiz_nml_data"
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
  ! ======================
  ! Trim character strings
  ! ======================
  !
  do n = 1,size(bc_input)
    bc_input(n)%sem%profile_file = trim(adjustl(bc_input(n)%sem%profile_file))
  end do
  !
  ! Copy bc_input into ovar_bc_input
  !
  ovar_bc_input = bc_input_v0001_to_ovar( bc_input ) ! F2003 AUTO-REALLOC
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
  integer, intent(in) :: io  !< Unit number for input.
  !
continue
  !
  rewind (io)
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
  !.. Local Scalars ..
  integer :: ovar_namelist_version
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_nml_data"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (mypnum == glb_root) then
    select case (action)
      case (NMLaction_Write, NMLaction_WritePretty)
        !
        ! Set namelist_version to reflect this module
        ovar_namelist_version = namelist_version
        namelist_version = module_namelist_version
        !
        if (banner) write (io,104)
        call set_iow(io)
        call set_NMLaction(action)
        call ionml()
        if (banner) write (io,105)
        flush (io)
        !
        ! Restore namelist_version
        namelist_version = ovar_namelist_version
        !
      case default
        !
        write (error_message,106) action
        ! This is within an 'if (mypnum == glb_root)' block,
        ! so it needs to use abort instead of stop_mpi.
        call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
        !
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
  use nml_9999_mod, only : nml_ver_9999
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: input_fname
  !
  !.. Local Scalars ..
  integer :: inpt,ierr
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: fname
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_nml_file"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  fname = trim(adjustl(input_fname)) ! F2003 AUTO-REALLOC
  !
  ! Allocate all the variables within the input namelist
  ! to their default sizes and values
  !
  call nml_ver_9999%alloc_data()  ! alloc new arrays
  call nml_ver_0001%alloc_data()  ! alloc old arrays
  !< Same as call alloc_nml_data
  !
  ! Set namelist variables to default values.
  !
  call nml_ver_9999%set_data(NMLaction_SetValue)  ! set new arrays
  call nml_ver_0001%set_data(NMLaction_SetValue)  ! set old arrays
  !< Same as call set_nml_data
  !
  ! Open the input file
  !
  call debug_timer(start_timer,"actually reading the input file")
  open (newunit=inpt,file=fname,status="old",action="read",form="formatted", &
                     iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Read in the input file using the input namelist
  !
  call nml_ver_0001%read_data(inpt)
  !< Same as call read_nml_data
  !
  ! Close the input file now that we are finished reading it
  !
  close (inpt,iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
  call debug_timer(stop_timer,"actually reading the input file")
  !
  ! Resize all the allocatable variables within the input
  ! namelist to sizes more appropriate to their input values
  !
  call nml_ver_9999%resiz_data() ! resize new arrays, dont care about old arrays
  !
  call debug_timer(leaving_procedure,pname)
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
  !.. Local Scalars ..
  integer     :: io,ierr
  logical(lk) :: banner
  logical(lk) :: allocations_required
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: delim
  character(len=:), allocatable :: fname
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_nml_file"
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
  ! NOTE: THIS SHOULD NOT BE USED FOR WRITING OUT THE FINAL PROCESSED
  !       NAMELIST TO STDOUT BEFORE STARTING THE SOLVER. Use write_nml_data
  !       to do that with iout as the unit argument.
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
                        opt_n_species=3, &
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
     !case (NMLaction_Read)        ; delim = 'NONE'       ! F2003 AUTO-REALLOC
      case (NMLaction_Write)       ; delim = 'APOSTROPHE' ! F2003 AUTO-REALLOC
      !> delim is set to 'none' for WritePretty because the subroutine that
      !> writes a character string to the file will add its own apostrophe,
      !> so this prevents it from writing double apostrophes before and after
      !> the character string.
      case (NMLaction_WritePretty) ; delim = 'NONE'       ! F2003 AUTO-REALLOC
      case default
        write (error_message,101) action
        ! This is within an 'if (mypnum == glb_root)' block,
        ! so it needs to use abort instead of stop_mpi.
        call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
    end select
    !
    fname = trim(adjustl(output_fname)) ! F2003 AUTO-REALLOC
    !
    open (newunit=io,file=fname,status="replace", &
                     action="write",position="rewind", &
                     delim=delim,recl=5000, &
                     iostat=ierr,iomsg=error_message)
    call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
    !
    call write_nml_data(action,io,banner)
    !
    close (io,iostat=ierr,iomsg=error_message)
    call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
    !
    write (iout,201) fname
    !
  end if
  !
  101 format("Coding Error:  Invalid arg, action = ",i0)
  201 format(/,"Default namelist input file written to:  ",A)
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
  call ionml_input()
  !
end subroutine ionml
!
!###############################################################################
!> I/O for namelist block 'input'.
!###############################################################################
subroutine ionml_input()
  !
  !.. Use Statements ..
  use ovar, only : Marker_MinMax_Default, Marker_Coeff_Default
  !
  ! ########################
  ! ###  Input Namelist  ###
  ! ########################
  !
  !{{{ input namelist definition
  ! GRID INFORMATION
  namelist / input / grid_format, gridfile, plot3d_cutfile, plot3d_bcsfile, &
                     grid_scaling_factor, internal_grid_dof, internal_grid_order
  ! PLOT3D COARSENING OPTIONS
  namelist / input / plot3d_i_stride, plot3d_j_stride, plot3d_k_stride, &
                     plot3d_all_strides, plot3d_agglomerate_order, &
                     parametric_agglomeration
  ! INTERPOLATION UTILITY GRID INFORMATION
  namelist / input / interpto_grid_format, interpto_gridfile, &
                     interpto_plot3d_cutfile, interpto_plot3d_bcsfile, &
                     interpto_grid_scaling_factor
  ! INTERPOLATION UTILITY PLOT3D COARSENING OPTIONS
  namelist / input / interpto_plot3d_i_stride, interpto_plot3d_j_stride, &
                     interpto_plot3d_k_stride, interpto_plot3d_all_strides, &
                     interpto_plot3d_agglomerate_order
  ! INTERPOLATION UTILITY SOLUTION ORDER
  namelist / input / interpto_solution_order, interpto_all_orders, &
                     interpto_projection_order, interpto_quadrature_order, &
                     interpto_output_order, interpto_error_order
  ! BASE FR OPTIONS
  namelist / input / solution_order, flux_divergence, correction_function
  ! LOCATION OF SOLUTION, FLUX, AND TRIANGULAR POINTS
  namelist / input / loc_solution_pts, loc_flux_pts, loc_triangle_pts, &
                     loc_tetrahedra_pts
  ! LOCATION OF INITIALIZATION POINTS
  namelist / input / loc_init_pts
  ! RUNGE-KUTTA OPTIONS
  namelist / input / Runge_Kutta_Scheme, num_rk_stages
  ! COMMON FLUX OPTIONS
  namelist / input / invs_flux_method, visc_flux_method, visc_variables_to_ave
  ! TIME STEPPING INFORMATION
  namelist / input / num_timesteps, Final_Time, constant_dt, &
                     Timestep_Type, minimum_timesteps
  ! CFL CONTROL VARIABLES
  namelist / input / cfl, cfl_beg, cfl_end, cfl_cycles, cfl_start_iter
  ! OVER-INTEGRATION OPTIONS
  namelist / input / projection_order, quadrature_order, loc_quadrature_pts, &
                     combine_common_flux_projections
  ! ERROR-INTEGRATION OPTIONS
  namelist / input / error_order, loc_error_pts
  ! FORCE GLOBAL ORDER AND SOLUTION POINTS
  namelist / input / all_orders, all_points
  ! PARALLEL OPTIONS
  namelist / input / metis_option, use_parmetis
  ! INITIALIZATION / TEST CASE
  namelist / input / itestcase, jet_init_func
  ! RESTART OPTIONS
  namelist / input / load_restart_file, restart_interval, &
                     restart_file, &
                     restart_output_uses_host_roots, &
                     CFL_reduction_for_NaN_restart, &
                     num_saved_time_ave_restarts, &
                     restart_frequency, &
                     completely_disable_restart_output
  ! OPTIONS TO USE SOLUTION FILE FOR INITIAL CONDITIONS
  namelist / input / solution_format, solution_file
  ! OPTIONS FOR SEM
  namelist / input / use_improved_SEM_normalization, &
                     use_SEM_profile_for_initialization, &
                     SEM_eddy_volume_is_spherical, &
                     SEM_box_num_eddies_multiplier, &
                     SEM_use_random_seed_from_restart, &
                     SEM_include_thermo_fluctuations, &
                     SEM_use_restart_eddies
  ! FLAG TO ACTIVATE A TEST OF THE SEM CODE
  namelist / input / average_sem_inflow_plane
  ! STATISTICAL OUTPUT OPTIONS
  namelist / input / iter_out_interval, results_interval, output_eperk_cfl, &
                     disable_header_rewrites, output_summary_file
  ! SOLUTION OUTPUT OPTIONS
  namelist / input / output_dir, output_interval, continuous_output, &
                     output_order, loc_output_pts, output_frequency
  ! SELECTION OF OUTPUT VARIABLES IN THE CGNS SOLUTION FILES
  namelist / input / output_vorticity, output_conserved_variables, &
                     output_momentum, output_cylindrical_coordinates, &
                     output_cylindrical_velocity, &
                     output_relative_velocity, output_relative_conditions, &
                     output_temperature, output_thermo_variables, &
                     output_gradients, output_transport_variables, &
                     output_density_gradient, &
                     output_density_gradient_magnitude, &
                     output_pressure_gradient, &
                     output_pressure_gradient_magnitude, &
                     output_pressure_dmat_gradient, &
                     output_temperature_gradient, &
                     output_temperature_gradient_magnitude, &
                     output_temperature_dmat_gradient, &
                     output_schlieren, output_modified_schlieren, &
                     output_stagnation_temperature, &
                     output_stagnation_pressure, &
                     output_stagnation_variables, output_marker_variables, &
                     output_grid_details, &
                     output_local_timestep, output_qcriterion, &
                     output_species_density, &
                     output_species_mass_fractions, &
                     output_species_mole_fractions, &
                     output_extra_species_variables, &
                     output_debug_species_variables, &
                     output_grdnt_species_variables, &
                     output_cellave_density, &
                     output_exact_density, &
                     output_exact_pressure
  ! GOVERNING EQUATIONS
  namelist / input / governing_equations, turbulence_model
  ! MOVING REFERENCE FRAME
  namelist / input / moving_ref_frame
  ! CHEMICAL REACTION SOURCE TERMS
  namelist / input / disable_chemical_reactions, semi_implicit_chemistry, &
                     traceback_unconverged_newton_temps
  !
  ! IGNITION CONTROL SOURCE TERMS
  namelist / input / ignition_method, &
                     ignition_start_timestep, &
                     ignition_duration, &
                     ignition_prod2reac_L2_ratio_goal, &
                     ignition_max_temperature, &
                     ignition_spark_cap, &
                     ignition_source_magnitude, &
                     ignition_source_radius, &
                     ignition_source_location
  ! NONDIMENSIONAL REFERENCE CONDITIONS
  namelist / input / machref, reyref, gam
  ! REALIZABLE TOTAL TEMPERATURE LIMITER
  namelist / input / ideal_total_temperature, limit_total_temperature
  ! REFERENCE CONDITIONS
  namelist / input / rhoref, tref, pref, rgasref
  ! REFERENCE TURBULENCE CONDITIONS
  namelist / input / tinref, tvrref
  ! REFERENCE ANGLES OF ATTACK
  namelist / input / alpha_aoaref, beta_aoaref
  ! LAMINAR AND TURBULENT PRANDTL NUMBERS
  namelist / input / Pr, Prt
  ! LAMINAR AND TURBULENT SCHMIDT NUMBERS
  namelist / input / schmidt_lam, schmidt_turb
  ! SUTHERLANDS LAW COEFFICIENTS
  namelist / input / suth_muref, suth_Tref, suth_Sref
  ! FILTER OPTIONS
  namelist / input / Filter_Option , &
                     Filter_etac   , &
                     Filter_Order  , &
                     Filter_Alpha  , &
                     Pyramid_Filter_Option
  ! LIMITER FILTER OPTIONS
  namelist / input / Limiter_Filter_Option      , &
                     Limiter_Filter_etac        , &
                     Limiter_Filter_Order       , &
                     Limiter_Filter_Alpha       , &
                     Limiter_Use_Dynamic_Filter , &
                     Track_Limiter_Filter_Dissipation
  ! SPONGE FILTER OPTIONS
  namelist / input / use_sponge_filter, &
                     Sponge_Filter_etac, &
                     Sponge_Filter_Order, &
                     Sponge_Filter_Alpha, &
                     Sponge_Beyond_X
  ! LIMITER OPTIONS
  namelist / input / Limiter_Option, Resolution_Marker_Variable, &
                     Jump_Marker_Variable, MLP_Marker_Variable
  ! LIMITER MARKER SENSITIVITY OPTIONS
  namelist / input / Limiter_Marker_Min, Limiter_Marker_Max, &
                     Limiter_Marker_c1, Limiter_Marker_c2, &
                     Limiter_Marker_dpsi
  ! CONTROL OVER FREQUENCY OF FILTERING/LIMITING
  namelist / input / Filtering_Interval, Sponge_Filter_Interval, &
                     Limiting_Interval, Limiting_Interval_Output
  ! CHECK/FIX BAD THERMO PROPERTIES (PRES,TEMP,DENS)
  namelist / input / thermofix
  ! OUTPUT THE STATISTICS OF THE LIMITER
  namelist / input / output_limiter_statistics
  ! OPTION TO OUTPUT FILTER DIAGONALS TO FILE
  namelist / input / output_filter_diagonals
  ! ISENTROPIC VORTEX CONSTANTS
  namelist / input / vortex_bigR, vortex_bigGam, VortexGrid_Lref
  ! TAYLOR-GREEN VORTEX CONSTANTS
  namelist / input / VortexGrid_DOF, VortexGrid_random_perturb, &
                     VortexGrid_perturb_percent, VortexGrid_Hex_to_Pyramids, &
                     VortexGrid_is_3D
  ! TAYLOR-COUETTE FLOW CONSTANTS
  namelist / input / TCF_Inner_Radius, TCF_Inner_Wall_Spacing, &
                     TCF_Outer_Radius, TCF_Outer_Wall_Spacing, &
                     TCF_Aspect_Ratio, TCF_Axial_Cells, &
                     TCF_Radial_Cells, TCF_Azimuthal_Cells, &
                     TCF_Theta_Shift,  TCF_Periodic_Sectors
  ! FLAME SPEED TEST PARAMETERS
  namelist / input / flame_speed_test_dof
  ! SOD SHOCK TUBE TEST PARAMETERS
  namelist / input / sod_test_dof, sod_test_right_vector
  ! DENSITY PULSE CONVECTION TEST PARAMETERS
  namelist / input / pulse_test_form, pulse_test_form_dim, &
                     pulse_test_dof, pulse_test_length, pulse_test_center, &
                     pulse_test_width, pulse_test_gap, pulse_test_repeat
  ! SHU-OSHER SHOCK TUBE TEST PARAMETERS
  namelist / input / shu_osher_test_dof
  ! CHANNEL FLOW PARAMETERS
  namelist / input / channel_body_force, &
                     channel_spatial_ave_restart, channel_nondim_output_coord, &
                     channel_grid_ave_xplus, channel_grid_ave_zplus, &
                     channel_grid_ave_yplus, channel_grid_yplus_wall, &
                     channel_yplus_of_one
  ! PIPE CHANNEL FLOW PARAMETERS
  namelist / input / pipe_channel_angle
  ! OPTIONS FOR MODIFYING BC ALGORITHMS
  namelist / input / walls_are_exact, sub_inflow_method, use_bc_conflicts
  ! CGNS OUTPUT OPTIONS
  namelist / input / interpolate_before_output_variable, &
                     cgns_high_order_output, cgns_mixed_output, &
                     cgns_output_restart_time_ave_data, &
                     convert_restart_to_cgns, completely_disable_cgns, &
                     delete_converted_restart_files
  ! CGNS INPUT OPTIONS
  namelist / input / correct_cgns_grid_points, force_cgns_linear_grid_points
  ! CONVERGENCE CRITERIA
  namelist / input / convergence_order_abs_res, convergence_order_max_res
  ! PROFILE PARALLEL IO
  namelist / input / profile_io_all, profile_io_restart, &
                     profile_io_cgns, profile_io_cgns_subs
  ! OPTIONS TO OUTPUT BOUNDARY PROFILE
  namelist / input / output_bnd_profiles
  ! CONTROL DOT_PRODUCT PROCEDURES BOUND TO MATRIX DERIVED-TYPE
  namelist / input / use_unrolled_dot_products
  ! LUSTRE FILE SYSTEM STRIPING PARAMETERS
  namelist / input / lustre_stripe_count, lustre_stripe_size
  ! OPENACC INPUT OPTIONS
  namelist / input / use_async_for_error
  ! POINT PROBE(S)
  namelist / input / probes_collection_interval, probes_write_interval, &
                     point_probes, rake_probes, cell_probes
  ! TIME AVERAGING
  namelist / input / output_time_averaging, time_ave_file, &
                     time_scaling_factor, time_ave_vel_is_axisymm, &
                     start_averaging_timestep, start_averaging_tstar
  ! SPATIAL AVERAGE RESTART FILES
  namelist / input / spatial_ave_restart_x_and_z_planes, &
                     spatial_ave_restart_z_planes
  ! UTILITY TO TIME AVERAGE EXISTING RESTART FILES
  namelist / input / time_average_restart_files, restart_list_file, &
                     modify_time_ave_restart_history
  ! MEMORY USAGE
  namelist / input / dump_memory_usage
  ! DUMP FLUX INFORMATION FOR DEBUGGING
  namelist / input / dump_fluxes, check_flux_point_coordinates
  ! ACTIVATE SPECIAL GRID FOR DEBUGGING THE POSTPROCESSOR
  namelist / input / postproc_debug_grid
  ! OPTION TO OUTPUT PLOT3D GRID TO TECPLOT IN UNSTRUCTURED FORMAT
  namelist / input / output_plot3d_to_tecplot
  ! OPTION TO OUTPUT CONTENTS OF BFACE ARRAY
  namelist / input / output_bface_array
  !
  namelist / input / force_nd_subroutines, dmat_grad_temp_order, &
                     newton_rlx, conv_temp, &
                     chemistry_sub_iterations, ignore_thermo_newton_warnings
  ! INPUT BOUNDARY CONDITIONS
  namelist / input / bc_input
  !
  ! GEOMETRY WEIGHTS FOR METIS PARTITIONING
  namelist / input / geom_wts
  !}}}
  !
  !.. Local Scalars ..
  integer  :: i,n,ierr
  integer  :: n_order
  integer  :: action_NewDefault, action_SetLimiter
  integer  :: error1, error2
  real(wp) :: l10p
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: text, text2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "ionml_input"
  real(wp), parameter :: probes_xyz_default(1:3) = [-999._wp,-999._wp,-999._wp]
  !
continue
  !
  ! Load local variables.
  !
  Use_Sponge_Filter = ( Sponge_Filter_Option == 1 ) ! logical to integer
  !
  ! Convert the bc_input structure from current to this version 1 reader.
  !   New: bc_input (ovar_mod)  =>  Old: bc_input (ver.1)
  !
  bc_input = bc_input_ovar_to_v0001( ovar_bc_input ) ! F2003 AUTO-REALLOC
  !
  ! Convert the new variables to the old.
  !   New: Limiter_Filter_Option  => Old: Limiter_Filter_Option
  !   New: Limiter_Option            Old: Limiter_Option
  !                                  Old: Limiter_Use_Dynamic_Filter
  !   New: Marker_Option          => Old: Resolution_Marker_Variable
  !   New: Marker_Variable           Old:       Jump_Marker_Variable
  !                                  Old:        MLP_Marker_Variable
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue )
      ! Conversion not needed, values will be set below.
      error1 = 0
    case default
      call convert_new_to_old(error1)
  end select
  !
  ! Some variables cannot have dummy values, so we impose NMLaction_SetValue.
  ! The reasons differ, so we define separate constants here:
  !   NewDefault  Default value changed.
  !   SetLimiter  Certain limiter/marker options need explicit values
  !               for subroutine convert_old_to_new() to work properly,
  !               because dummy values don't translate well.
  !
  select case (NMLaction)
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1 )
      action_NewDefault = NMLaction_SetValue
      action_SetLimiter = NMLaction_SetValue
    case default
      action_NewDefault = NMLaction
      action_SetLimiter = NMLaction
  end select
  !
  ! Perform an action on all of the &input namelist variables
  ! depending on the current value of NMLaction
  !
  select case (NMLaction)
    !
    ! If NMLaction is WritePretty, write out the namelist contents to a
    ! "pretty"-formatted file; otherwise, set the value for each variable
    ! to the default or dummy values.
    !
    case ( NMLaction_SetDummy0, NMLaction_SetDummy1, NMLaction_SetValue, &
           NMLaction_WritePretty )
      call iow_header('&input')
      !
!{{{  ! GRID INFORMATION
      call set('Grid_Format                      ',                          &
                grid_format                       ,Internally_Generated_Grid )
      call set('Gridfile                         ',                          &
                gridfile                          , ''                       )
      call set('Plot3d_Cutfile                   ',                          &
                plot3d_cutfile                    , ''                       )
      call set('Plot3d_BCSfile                   ',                          &
                plot3d_bcsfile                    , ''                       )
      call set('Grid_Scaling_Factor              ',                          &
                grid_scaling_factor               , 1.0_wp                   )
      call set('Internal_Grid_DoF                ',                          &
                internal_grid_dof                 , 0                        )
      call set('Internal_Grid_Order              ',                          &
                internal_grid_order               ,-1                        )
      !}}}
!{{{  ! PLOT3D COARSENING OPTIONS
      call set('Plot3d_I_Stride                  ',                          &
                plot3d_i_stride                   ,-1                        )
      call set('Plot3d_J_Stride                  ',                          &
                plot3d_j_stride                   ,-1                        )
      call set('Plot3d_K_Stride                  ',                          &
                plot3d_k_stride                   ,-1                        )
      call set('Plot3d_All_Strides               ',                          &
                plot3d_all_strides                ,-1                        )
      call set('Plot3d_Agglomerate_Order         ',                          &
                plot3d_agglomerate_order          ,-1                        )
      call set('Parametric_Agglomeration         ',                          &
                parametric_agglomeration          , PmetAgglom__OFF          )
      !}}}
!{{{  ! INTERPOLATION UTILITY GRID INFORMATION
      call set('Interpto_Grid_Format             ',                &
                interpto_grid_format              , Unknown_Format )
      call set('Interpto_Gridfile                ',                &
                interpto_gridfile                 , ''             )
      call set('Interpto_Plot3d_Cutfile          ',                &
                interpto_plot3d_cutfile           , ''             )
      call set('Interpto_Plot3d_BCSfile          ',                &
                interpto_plot3d_bcsfile           , ''             )
      call set('Interpto_Grid_Scaling_Factor     ',                &
                interpto_grid_scaling_factor      , 1.0_wp         )
      !}}}
!{{{  ! INTERPOLATION UTILITY PLOT3D COARSENING OPTIONS
      call set('Interpto_Plot3d_I_Stride         ',                &
                interpto_plot3d_i_stride          ,-1              )
      call set('Interpto_Plot3d_J_Stride         ',                &
                interpto_plot3d_j_stride          ,-1              )
      call set('Interpto_Plot3d_K_Stride         ',                &
                interpto_plot3d_k_stride          ,-1              )
      call set('Interpto_Plot3d_All_Strides      ',                &
                interpto_plot3d_all_strides       ,-1              )
      call set('Interpto_Plot3d_Agglomerate_Order',                &
                interpto_plot3d_agglomerate_order ,-1              )
      !}}}
!{{{  ! INTERPOLATION UTILITY SOLUTION ORDER
      call set('Interpto_Solution_Order          ',                &
                interpto_solution_order           , 3              )
      call set('Interpto_All_Orders              ',                &
                interpto_all_orders               ,-3              )
      call set('Interpto_Projection_Order        ',                &
                interpto_projection_order         ,-3              )
      call set('Interpto_Quadrature_Order        ',                &
                interpto_quadrature_order         ,-3              )
      call set('Interpto_Output_Order            ',                &
                interpto_output_order             ,-3              )
      call set('Interpto_Error_Order             ',                &
                interpto_error_order              ,-3              )
      !}}}
!{{{  ! BASE FR OPTIONS
      call set('Solution_Order                   ',                         &
                solution_order                    , 3                      ,&
                action_SetLimiter)
      call set('Flux_Divergence                  ',                         &
                flux_divergence                   , LagrangePolynomial      )
      call set('Correction_Function              ',                         &
                correction_function               , Correction_gDG          )
      !}}}
!{{{  ! LOCATION OF SOLUTION, FLUX, AND TRIANGULAR POINTS
      call set('Loc_Solution_Pts                 ',                         &
                loc_solution_pts                  , Legendre_Gauss          )
      call set('Loc_Flux_Pts                     ',                         &
                loc_flux_pts                      , Legendre_Gauss          )
      call set('Loc_Triangle_Pts                 ',                         &
                loc_triangle_pts                  , WilliamsShunn_TriPoints )
      call set('Loc_Tetrahedra_Pts               ',                         &
                loc_tetrahedra_pts                , ShunnHam_TetPoints      )
      !}}}
!{{{  ! LOCATION OF INITIALIZATION POINTS
      call set('Loc_Init_Pts                     ',                         &
                loc_init_pts                      , 0                       )
      !}}}
!{{{  ! RUNGE-KUTTA OPTIONS
      call set('Runge_Kutta_Scheme               ',                         &
                Runge_Kutta_Scheme                , OERK_11Stage_3rdOrder   )
              ! action_NewDefault, was TVD3_RK, now OERK_11Stage_3rdOrder
      call set('Num_RK_Stages                    ',                         &
                num_rk_stages                     , 3                       )
      !}}}
!{{{  ! COMMON FLUX OPTIONS
      call set('Invs_Flux_Method                 ',                         &
                invs_flux_method                  , invs_flux_ROE           )
      call set('Visc_Flux_Method                 ',                         &
                visc_flux_method                  , visc_flux_BR2           )
      call set('Visc_Variables_To_Ave            ',                         &
                visc_variables_to_ave             , visc_interface_fluxes   )
      !}}}
!{{{  ! TIME STEPPING INFORMATION
      call set('Num_Timesteps                    ',                         &
                num_timesteps                     , 25000                   )
      call set('Final_Time                       ',                         &
                Final_Time                        ,-10.0_wp                 )
      call set('Constant_Dt                      ',                         &
                constant_dt                       , 1.0e-6_wp               )
      call set('Timestep_Type                    ',                         &
                Timestep_Type                     , Global_Cell_Timestep    )
      call set('Minimum_Timesteps                ',                         &
                minimum_timesteps                 , 100                     )
      !}}}
!{{{  ! CFL CONTROL VARIABLES
      !   cfl     has new ver2 default.
      !   cfl_end has new ver2 default, so also force:  cfl_beg, cfl_cycles
      call set('CFL                              ',                         &
                cfl                               , 1.0_wp                  )
              ! action_NewDefault, was 0.1, now 1.0
      call set('CFL_Beg                          ',                         &
                cfl_beg                           , 0.1_wp                  )
      call set('CFL_End                          ',                         &
                cfl_end                           , 1.0_wp                  )
              ! action_NewDefault, was 0.1, now 1.0
      call set('CFL_Cycles                       ',                         &
                cfl_cycles                        , 1                       )
      call set('CFL_Start_Iter                   ',                         &
                cfl_start_iter                    , 0                       )
      !}}}
!{{{  ! OVER-INTEGRATION OPTIONS
      call set('Projection_Order                 ',                         &
                projection_order                  ,-3                       )
      call set('Quadrature_Order                 ',                         &
                quadrature_order                  ,-3                       )
      call set('Loc_Quadrature_Pts               ',                         &
                loc_quadrature_pts                , Legendre_Gauss          )
      call set('Combine_Common_Flux_Projections  ',                         &
                combine_common_flux_projections   , 0                       )
      !}}}
!{{{  ! ERROR-INTEGRATION OPTIONS
      call set('Error_Order                      ',                         &
                error_order                       ,-3                       )
      call set('Loc_Error_Pts                    ',                         &
                loc_error_pts                     , Legendre_Gauss          )
      !}}}
!{{{  ! FORCE GLOBAL ORDER AND SOLUTION POINTS
      call set('All_Orders                       ',                         &
                all_orders                        ,-3                      ,&
                action_SetLimiter)
      call set('All_Points                       ',                         &
                all_points                        ,-1                       )
      !}}}
!{{{  ! PARALLEL OPTIONS
!     !> METIS_OPTION = 3  ! k-way part, min comm volume
      call set('Metis_Option                     ',                         &
                metis_option                      , 3                       )
      call set('Use_Parmetis                     ',                         &
                use_parmetis                      , 0                       )
      !}}}
!{{{  ! INITIALIZATION / TEST CASE
      call set('Itestcase                        ',                         &
                itestcase                         , Generic_MeanFlow        )
              ! action_NewDefault, was Taylor_Green_Vortex, now Generic_Meanflow
      call set('Jet_Init_Func                    ',                         &
                jet_init_func                     , 2                       )
      !}}}
!{{{  ! RESTART OPTIONS
      call set('Load_Restart_File                ',                         &
                load_restart_file                 , fals                    )
      call set('Restart_Interval                 ',                         &
                restart_interval                  , 0                       )
      call set('Restart_File                     ',                         &
                restart_file                      , 'out.rst'               )
      call set('Restart_Output_Uses_Host_Roots   ',                         &
                restart_output_uses_host_roots    , true                    ) 
      call set('CFL_Reduction_For_NaN_Restart    ',                         &
                CFL_reduction_for_NaN_restart     , 2.0_wp                  ) 
      call set('Num_Saved_Time_Ave_Restarts      ',                         &
                num_saved_time_ave_restarts       , 2                       )
      call set('Restart_Frequency                ',                         &
                restart_frequency                 ,-one                     )
      call set('Completely_Disable_Restart_Output',                         &
                completely_disable_restart_output , fals                    )
      !}}}
!{{{  ! OPTIONS TO USE SOLUTION FILE FOR INITIAL CONDITIONS
      call set('Solution_Format                   ',                        &
                solution_format                    , Unknown_Format         )
      call set('Solution_File                     ',                        &
                solution_file                      , ''                     )
      !}}}
!{{{  ! OPTIONS FOR SEM
      call set('Use_Improved_SEM_Normalization    ',                        &
                use_improved_sem_normalization     , true                   )
              ! action_NewDefault, was false, now true
      call set('Use_SEM_Profile_For_Initialization',                        &
                use_sem_profile_for_initialization , fals                   )
      call set('SEM_Eddy_Volume_Is_Spherical      ',                        &
                sem_eddy_volume_is_spherical       , fals                   )
      call set('SEM_Box_Num_Eddies_Multiplier     ',                        &
                sem_box_num_eddies_multiplier      , one                    )
      call set('SEM_Use_Random_Seed_From_Restart  ',                        &
                sem_use_random_seed_from_restart   , fals                   )
      call set('SEM_Include_Thermo_Fluctuations   ',                        &
                sem_include_thermo_fluctuations    , fals                   )
      call set('SEM_Use_Restart_Eddies            ',                        &
                sem_use_restart_eddies             , true                   )
      !}}}
!{{{  ! FLAG TO ACTIVATE A TEST OF THE SEM CODE
      call set('Average_SEM_Inflow_Plane          ',                        &
                average_sem_inflow_plane           , fals                   )
      !}}}
!{{{  ! STATISTICAL OUTPUT OPTIONS
      call set('Iter_Out_Interval                 ',                        &
                iter_out_interval                  , 1                      )
      call set('Results_Interval                  ',                        &
                results_interval                   , 0                      )
      call set('Output_EPERK_CFL                  ',                        &
                output_eperk_cfl                   , true                   )
      call set('Disable_Header_Rewrites           ',                        &
                disable_header_rewrites            , fals                   )
      call set('Output_Summary_File               ',                        &
                output_summary_file                , fals                   )
      !}}}
!{{{  ! SOLUTION OUTPUT OPTIONS
      call set('Output_Dir                        ',                        &
                output_dir                         , '.'                    )
      call set('Output_Interval                   ',                        &
                output_interval                    ,-5                      )
      call set('Continuous_Output                 ',                        &
                continuous_output                  , fals                   )
      call set('Output_Order                      ',                        &
                output_order                       ,-3                      )
      call set('Loc_Output_Pts                    ',                        &
                loc_output_pts                     , Equi_Distant_Points    )
      call set('Output_Frequency                  ',                        &
                output_frequency                   ,-1.0_wp                 )
      !}}}
!{{{  ! SELECTION OF OUTPUT VARIABLES IN THE CGNS SOLUTION FILES
      call set('Output_Vorticity                     ',          &
                output_vorticity                      , fals     )
      call set('Output_Conserved_Variables           ',          &
                output_conserved_variables            , fals     )
      call set('Output_Momentum                      ',          &
                output_momentum                       , fals     )
      call set('Output_Cylindrical_Coordinates       ',          &
                output_cylindrical_coordinates        , fals     )
      call set('Output_Cylindrical_Velocity          ',          &
                output_cylindrical_velocity           , fals     )
      call set('Output_Relative_Velocity             ',          &
                output_relative_velocity              , fals     )
      call set('Output_Relative_Conditions           ',          &
                output_relative_conditions            , fals     )
      call set('Output_Temperature                   ',          &
                output_temperature                    , fals     )
      call set('Output_Thermo_Variables              ',          &
                output_thermo_variables               , fals     )
      call set('Output_Gradients                     ',          &
                output_gradients                      , fals     )
      call set('Output_Density_Gradient              ',          &
                output_density_gradient               , fals     )
      call set('Output_Density_Gradient_Magnitude    ',          &
                output_density_gradient_magnitude     , fals     )
      call set('Output_Pressure_Gradient             ',          &
                output_pressure_gradient              , fals     )
      call set('Output_Pressure_Gradient_Magnitude   ',          &
                output_pressure_gradient_magnitude    , fals     )
      call set('Output_Pressure_Dmat_Gradient        ',          &
                output_pressure_dmat_gradient         , fals     )
      call set('Output_Temperature_Gradient          ',          &
                output_temperature_gradient           , fals     )
      call set('Output_Temperature_Gradient_Magnitude',          &
                output_temperature_gradient_magnitude , fals     )
      call set('Output_Temperature_Dmat_Gradient     ',          &
                output_temperature_dmat_gradient      , fals     )
      call set('Output_Schlieren                     ',          &
                output_schlieren                      , fals     )
      call set('Output_Modified_Schlieren            ',          &
                output_modified_schlieren             , fals     )
      call set('Output_Transport_Variables           ',          &
                output_transport_variables            , fals     )
      call set('Output_Stagnation_Temperature        ',          &
                output_stagnation_temperature         , fals     )
      call set('Output_Stagnation_Pressure           ',          &
                output_stagnation_pressure            , fals     )
      call set('Output_Stagnation_Variables          ',          &
                output_stagnation_variables           , fals     )
      call set('Output_Marker_Variables              ',          &
                output_marker_variables               , fals     )
      call set('Output_Qcriterion                    ',          &
                output_qcriterion                     , fals     )
      call set('Output_Grid_Details                  ',          &
                output_grid_details                   , fals     )
      call set('Output_Local_Timestep                ',          &
                output_local_timestep                 , fals     )
      call set('Output_Species_Density               ',          &
                output_species_density                , fals     )
      call set('Output_Species_Mass_Fractions        ',          &
                output_species_mass_fractions         , true     )
              ! action_NewDefault, was false, now true
      call set('Output_Species_Mole_Fractions        ',          &
                output_species_mole_fractions         , fals     )
      call set('Output_Extra_Species_Variables       ',          &
                output_extra_species_variables        , fals     )
      call set('Output_Debug_Species_Variables       ',          &
                output_debug_species_variables        , fals     )
      call set('Output_Grdnt_Species_Variables       ',          &
                output_grdnt_species_variables        , fals     )
      call set('Output_CellAve_Density               ',          &
                output_cellave_density                , fals     )
      call set('Output_Exact_Density                 ',          &
                output_exact_density                  , fals     )
      call set('Output_Exact_Pressure                ',          &
                output_exact_pressure                 , fals     )
      !}}}
!{{{  ! GOVERNING EQUATIONS
      call set('Governing_Equations               ',                        &
                governing_equations                , NavierStokes_Eqns      )
      call set('Turbulence_Model                  ',                        &
                turbulence_model                   , laminar_flow           )
      !}}}
!{{{  ! MOVING REFERENCE FRAME
      call set('Moving_Ref_Frame%Rotation         ',                        &
                moving_ref_frame%rotation          , [zero,zero,zero]       )
      call set('Moving_Ref_Frame%Translation      ',                        &
                moving_ref_frame%translation       , [zero,zero,zero]       )
      call set('Moving_Ref_Frame%Spinup%Iter      ',                        &
                moving_ref_frame%spinup%iter       , 0                      )
      call set('Moving_Ref_Frame%Spinup%Tstar     ',                        &
                moving_ref_frame%spinup%tstar      , zero                   )
      !}}}
!{{{  ! CHEMICAL REACTION SOURCE TERMS
      call set('Disable_Chemical_Reactions        ',                        &
                disable_chemical_reactions         , fals                   )
      call set('Semi_Implicit_Chemistry           ',                        &
                semi_implicit_chemistry            , 0                      )
      call set('Traceback_Unconverged_Newton_Temps',                        &
                traceback_unconverged_newton_temps , fals                   )
      !}}}
!{{{  ! IGNITION CONTROL SOURCE TERMS
      call set('Ignition_Method                   ',                        &
                ignition_method                    , IGN_OFF                )
      call set('Ignition_Start_Timestep           ',                        &
                ignition_start_timestep            , 0                      )
      call set('Ignition_Duration                 ',                        &
                ignition_duration                  , 1                      )
      call set('Ignition_Prod2Reac_L2_Ratio_Goal  ',                        &
                ignition_prod2reac_L2_ratio_goal   ,-1.0_wp                 )
      call set('Ignition_Max_Temperature          ',                        &
                ignition_max_temperature           , 1600.0_wp              )
      call set('Ignition_Spark_Cap                ',                        &
                ignition_spark_cap                 , 1.0e+12_wp             )
      call set('Ignition_Source_Magnitude         ',                        &
                ignition_source_magnitude          , 20.0_wp                )
      call set('Ignition_Source_Radius            ',                        &
                ignition_source_radius             , zero                   )
      call set('Ignition_Source_Location          ',                        &
                ignition_source_location           , [zero,zero,zero]       )
      !}}}
!{{{  ! NONDIMENSIONAL REFERENCE CONDITIONS
      call set('Machref        ',machref           , 0.4_wp                 )
      call set('Reyref         ',reyref            ,-1.0_wp                 )
      call set('Gam            ',gam               , 1.4_wp                 )
      !}}}
!{{{  ! REALIZABLE TOTAL TEMPERATURE LIMITER
      call set('Ideal_Total_Temperature           ',                        &
                ideal_total_temperature            ,-1.0_wp                 )
      call set('Limit_Total_Temperature           ',                        &
                limit_total_temperature            , 0.2_wp                 )
      !}}}
!{{{  ! REFERENCE CONDITIONS
      call set('Rhoref         ',rhoref            ,-1.0_wp                 )
      call set('Tref           ',tref              , 300.0_wp               )
      call set('Pref           ',pref              , 100000.0_wp            )
      call set('Rgasref        ',rgasref           ,-1.0_wp                 )
      !}}}
!{{{  ! REFERENCE TURBULENCE CONDITIONS
      call set('Tinref         ',tinref            , 0.0_wp                 )
      call set('Tvrref         ',tvrref            , 1.0e-6_wp              )
      !}}}
!{{{  ! REFERENCE ANGLES OF ATTACK
      call set('Alpha_AoAref   ',alpha_aoaref      , 0.0_wp                 )
      call set('Beta_AoAref    ',beta_aoaref       , 0.0_wp                 )
      ! }}}
 !{{{  ! LAMINAR AND TURBULENT PRANDTL NUMBERS
      call set('Pr             ',pr                , 0.72_wp                )
      call set('Prt            ',prt               , 0.72_wp                )
      ! }}}
 !{{{  ! LAMINAR AND TURBULENT SCHMIDT NUMBERS
      call set('Schmidt_Lam    ',schmidt_lam       , 0.72_wp                )
      call set('Schmidt_Turb   ',schmidt_turb      , 0.72_wp                )
      ! }}}
 !{{{  ! SUTHERLANDS LAW COEFFICIENTS
      call set('Suth_Muref     ',suth_muref        , 1.716e-5_wp            )
      call set('Suth_Tref      ',suth_Tref         , 273.0_wp               )
      call set('Suth_Sref      ',suth_Sref         , 110.4_wp               )
      !}}}
!{{{  ! FILTER OPTIONS
      call set('Filter_Option                     ',                        &
                Filter_Option                      , 0                      )
      call set('Filter_Etac                       ',                        &
                Filter_Etac                        , 0.0_wp                 )
      call set('Filter_Order                      ',                        &
                Filter_Order                       , 16                     )
      call set('Filter_Alpha                      ',                        &
                Filter_Alpha                       , 36.0_wp                )
      call set('Pyramid_Filter_Option             ',                        &
                Pyramid_Filter_Option              , 1                      )
      !}}}
!{{{  ! LIMITER FILTER OPTIONS
      call set('Limiter_Filter_Option             ',                        &
                Limiter_Filter_Option              , 0                     ,&
                action_SetLimiter)
      call set('Limiter_Filter_Etac               ',                        &
                Limiter_Filter_Etac                , 0.0_wp                 )
      call set('Limiter_Filter_Order              ',                        &
                Limiter_Filter_Order               , 16                     )
      call set('Limiter_Filter_Alpha              ',                        &
                Limiter_Filter_Alpha               , 36.0_wp                )
      call set('Limiter_Use_Dynamic_Filter        ',                        &
                Limiter_Use_Dynamic_Filter         , fals                  ,&
                action_SetLimiter)
      call set('Track_Limiter_Filter_Dissipation  ',                        &
                Track_Limiter_Filter_Dissipation   , fals                   )
      !}}}
!{{{  ! SPONGE FILTER OPTIONS
      call set('Use_Sponge_Filter                 ',                        &
                Use_Sponge_Filter                  , fals                   )
      call set('Sponge_Filter_Etac                ',                        &
                Sponge_Filter_Etac                 , 0.0_wp                 )
      call set('Sponge_Filter_Order               ',                        &
                Sponge_Filter_Order                , 8                      )
      call set('Sponge_Filter_Alpha               ',                        &
                Sponge_Filter_Alpha                , 36.0_wp                )
      call set('Sponge_Beyond_X                   ',                        &
                Sponge_Beyond_X                    , huge(zero)             )
      !}}}
!{{{  ! LIMITER OPTIONS
      call set('Limiter_Option                    ',                        &
                Limiter_Option                     , 0                     ,&
                action_SetLimiter)
      call set('Resolution_Marker_Variable        ',                        &
                Resolution_Marker_Variable         , 0                     ,&
                action_SetLimiter)
      call set('Jump_Marker_Variable              ',                        &
                Jump_Marker_Variable               , 0                     ,&
                action_SetLimiter)
      call set('MLP_Marker_Variable               ',                        &
                MLP_Marker_Variable                , 0                     ,&
                action_SetLimiter)
      !}}}
!{{{  ! LIMITER MARKER SENSITIVITY OPTIONS
      call set('Limiter_Marker_Min                ',                        &
                Limiter_Marker_Min                 , Marker_MinMax_Default ,&
                action_SetLimiter)
      call set('Limiter_Marker_Max                ',                        &
                Limiter_Marker_Max                 , Marker_MinMax_Default ,&
                action_SetLimiter)
      call set('Limiter_Marker_C1                 ',                        &
                Limiter_Marker_c1                  , Marker_Coeff_Default  ,&
                action_SetLimiter)
      call set('Limiter_Marker_C2                 ',                        &
                Limiter_Marker_c2                  , Marker_Coeff_Default  ,&
                action_SetLimiter)
      call set('Limiter_Marker_Dpsi               ',                        &
                Limiter_Marker_dpsi                , Marker_Coeff_Default  ,&
                action_SetLimiter)
      !}}}
!{{{  ! CONTROL OVER FREQUENCY OF FILTERING/LIMITING
      call set('Filtering_Interval                ',                        &
                Filtering_Interval                 , 0                      )
      call set('Sponge_Filter_Interval            ',                        &
                Sponge_Filter_Interval             , 0                      )
      call set('Limiting_Interval                 ',                        &
                Limiting_Interval                  , 0                      )
      call set('Limiting_Interval_Output          ',                        &
                Limiting_Interval_Output           , 1                      )
      !}}}
!{{{  ! CHECK/FIX BAD THERMO PROPERTIES (PRES,TEMP,DENS)
      call set('Thermofix%Method                  ',                        &
                thermofix%method                   , TFM_OFF                )
      call set('Thermofix%Iterations              ',                        &
                thermofix%iterations               , 10                     )
      call set('Thermofix%Check_Cells             ',                        &
                thermofix%check_cells              , true                   )
      call set('Thermofix%Check_Faces             ',                        &
                thermofix%check_faces              , true                   )
      call set('Thermofix%Check_Nodes             ',                        &
                thermofix%check_nodes              , true                   )
      call set('Thermofix%Verbose                 ',                        &
                thermofix%verbose                  , true                   )
      call set('Thermofix%Scale_Factor            ',                        &
                thermofix%scale_factor             , 0.5_wp                 )
      call set('Thermofix%Min_Pressure            ',                        &
                thermofix%min_pressure             , 10000._wp              )
      call set('Thermofix%Min_Temperature         ',                        &
                thermofix%min_temperature          ,   100._wp              )
      call set('Thermofix%Min_Density             ',                        &
                thermofix%min_density              ,  0.001_wp              )
      !}}}
!{{{  ! OUTPUT THE STATISTICS OF THE LIMITER
      call set('Output_Limiter_Statistics         ',                        &
                output_limiter_statistics          , fals                   )
      !}}}
!{{{  ! OPTION TO OUTPUT FILTER DIAGONALS TO FILE
      call set('Output_Filter_Diagonals           ',                        &
                output_filter_diagonals            , fals                   )
      !}}}
!{{{  ! ISENTROPIC VORTEX CONSTANTS
      call set('Vortex_BigR               ',vortex_bigR               , 1.0_wp )
      call set('Vortex_BigGam             ',vortex_bigGam             , 5.0_wp )
      call set('VortexGrid_Lref           ',VortexGrid_Lref           ,10.0_wp )
      !}}}
!{{{  ! TAYLOR-GREEN VORTEX CONSTANTS
      call set('VortexGrid_DoF            ',VortexGrid_DOF            , 48     )
      call set('VortexGrid_Random_Perturb ',VortexGrid_random_perturb , fals   )
      call set('VortexGrid_Perturb_Percent',VortexGrid_perturb_percent, 30     )
      call set('VortexGrid_Hex_To_Pyramids',VortexGrid_Hex_to_Pyramids, fals   )
      call set('VortexGrid_Is_3D          ',VortexGrid_is_3D          , fals   )
      !}}}
!{{{  ! TAYLOR-COUETTE FLOW CONSTANTS
      call set('TCF_Inner_Radius          ',TCF_Inner_Radius      , half       )
      call set('TCF_Inner_Wall_Spacing    ',TCF_Inner_Wall_Spacing, 0.002_wp   )
      call set('TCF_Outer_Radius          ',TCF_Outer_Radius      , one        )
      call set('TCF_Outer_Wall_Spacing    ',TCF_Outer_Wall_Spacing, 0.0025_wp  )
      call set('TCF_Aspect_Ratio          ',TCF_Aspect_Ratio      , two*pi     )
      call set('TCF_Axial_Cells           ',TCF_Axial_Cells       , 24         )
      call set('TCF_Radial_Cells          ',TCF_Radial_Cells      , 20         )
      call set('TCF_Azimuthal_Cells       ',TCF_Azimuthal_Cells   , 48         )
      call set('TCF_Theta_Shift           ',TCF_Theta_Shift       , zero       )
      call set('TCF_Periodic_Sectors      ',TCF_Periodic_Sectors  , 1          )
      !}}}
!{{{  ! FLAME SPEED TEST PARAMETERS
      call set('Flame_Speed_Test_DoF      ',flame_speed_test_dof  , 0 )
      !}}}
!{{{  ! SOD SHOCK TUBE TEST PARAMETERS
      call set('Sod_Test_DoF              ',sod_test_dof          , 0 )
      call set('Sod_Test_Right_Vector     ',                          &
                sod_test_right_vector      , [one,zero,zero]          )
      !}}}
!{{{  ! DENSITY PULSE CONVECTION TEST PARAMETERS
      call set('Pulse_Test_Form    ',pulse_test_form    , 0                 )
      call set('Pulse_Test_Form_Dim',pulse_test_form_dim, 1                 )
      call set('Pulse_Test_DoF     ',pulse_test_dof     , 0                 )
      call set('Pulse_Test_Length  ',pulse_test_length  , 1.0_wp            )
      call set('Pulse_Test_Center  ',pulse_test_center  , [half,zero,zero]  )
      call set('Pulse_Test_Width   ',pulse_test_width   , 0.2_wp            )
      call set('Pulse_Test_Gap     ',pulse_test_gap     , 0.2_wp            )
      call set('Pulse_Test_Repeat  ',pulse_test_repeat  , 0                 )
      !}}}
!{{{  ! SHU-OSHER SHOCK TUBE TEST PARAMETERS
      call set('Shu_Osher_Test_DoF ',shu_osher_test_dof , 0 )
      !}}}
!{{{  ! CHANNEL FLOW PARAMETERS
      call set('Channel_Body_Force               ',                         &
                channel_body_force                , 0.0_wp                  )
      call set('Channel_Spatial_Ave_Restart      ',                         &
                channel_spatial_ave_restart       , fals                    )
      call set('Channel_Nondim_Output_Coord      ',                         &
                channel_nondim_output_coord       , fals                    )
      call set('Channel_Grid_Ave_Xplus           ',                         &
                channel_grid_ave_xplus            ,-30.0_wp                 )
      call set('Channel_Grid_Ave_Zplus           ',                         &
                channel_grid_ave_zplus            ,-15.0_wp                 )
      call set('Channel_Grid_Ave_Yplus           ',                         &
                channel_grid_ave_yplus            ,-15.0_wp                 )
      call set('Channel_Yplus_Of_One             ',                         &
                channel_yplus_of_one              , 8.42134e-6_wp           )
      call set('Channel_Grid_Yplus_Wall          ',                         &
                channel_grid_yplus_wall           ,  2.0_wp                 )
      !}}}
!{{{  ! PIPE CHANNEL FLOW PARAMETERS
      call set('Pipe_Channel_Angle               ',                         &
                pipe_channel_angle                , five                    )
      !}}}
!{{{  ! OPTIONS FOR MODIFYING BC ALGORITHMS
      call set('Walls_Are_Exact                   ',                        &
                walls_are_exact                    , true                   )
      call set('Sub_Inflow_Method                 ',                        &
                sub_inflow_method                  , 4                      )
      call set('Use_BC_Conflicts                  ',                        &
                use_bc_conflicts                   , fals                   )
      !}}}
!{{{  ! CGNS OUTPUT OPTIONS
      call set('Interpolate_Before_Output_Variable',                        &
                interpolate_before_output_variable , true                   )
      call set('CGNS_High_Order_Output            ',                        &
                cgns_high_order_output             , fals                   )
      call set('CGNS_Mixed_Output                 ',                        &
                cgns_mixed_output                  , true                   )
      call set('CGNS_Output_Restart_Time_Ave_Data ',                        &
                cgns_output_restart_time_ave_data  , true                   )
      call set('Convert_Restart_To_CGNS           ',                        &
                convert_restart_to_cgns            , fals                   )
      call set('Completely_Disable_CGNS           ',                        &
                completely_disable_cgns            , fals                   )
      call set('Delete_Converted_Restart_Files    ',                        &
                delete_converted_restart_files     , fals                   )
      !}}}
!{{{  ! CGNS INPUT OPTIONS
      call set('Correct_CGNS_Grid_Points          ',                        &
                correct_cgns_grid_points           , fals                   )
      call set('Force_CGNS_Linear_Grid_Points     ',                        &
                force_cgns_linear_grid_points      , fals                   )
      !}}}
!{{{  ! CONVERGENCE CRITERIA
      call set('Convergence_Order_Abs_Res         ',                        &
                convergence_order_abs_res          ,-16._wp                 )
      call set('Convergence_Order_Max_Res         ',                        &
                convergence_order_max_res          ,-16._wp                 )
      !}}}
!{{{  ! PROFILE PARALLEL IO
      call set('Profile_IO_All                    ',                        &
                profile_io_all                     , fals                   )
      call set('Profile_IO_Restart                ',                        &
                profile_io_restart                 , fals                   )
      call set('Profile_IO_CGNS                   ',                        &
                profile_io_cgns                    , fals                   )
      call set('Profile_IO_CGNS_Subs              ',                        &
                profile_io_cgns_subs               , fals                   )
      !}}}
!{{{  ! OPTIONS TO OUTPUT BOUNDARY PROFILE
      call set('Output_Bnd_Profiles               ',                        &
                output_bnd_profiles                , fals                   )
      !}}}
!{{{  ! CONTROL DOT_PRODUCT PROCEDURES BOUND TO MATRIX DERIVED-TYPE
      call set('Use_Unrolled_Dot_Products         ',                        &
                use_unrolled_dot_products          , true                   )
      !}}}
!{{{  ! LUSTRE FILE SYSTEM STRIPING PARAMETERS
      call set('Lustre_Stripe_Count               ',                        &
                lustre_stripe_count                , 64                     )
      call set('Lustre_Stripe_Size                ',                        &
                lustre_stripe_size                 , '4m'                   )
      !}}}
!{{{  ! OPENACC INPUT OPTIONS
      call set('Use_Async_For_Error               ',                        &
                use_async_for_error                , fals                   )
      !}}}
!{{{  ! PROBES
      call set('Probes_Collection_Interval        ',                        &
                Probes_Collection_Interval         , 0                      )
      call set('Probes_Write_Interval             ',                        &
                Probes_Write_Interval              , 100                    )
      !
      do n = 1, size(point_probes)
        text = '                               '
        write (text,101) 'Point_Probes(',n,')'
        text = trim(text)   ! F2003 realloc
        call set( text //       '%XYZ             ',                        &
                  point_probes(n)%xyz              , probes_xyz_default     )
      end do
      !
      do n = 1, size(rake_probes)
        text = '                               '
        write (text,101) 'Rake_Probes(',n,')'
        text = trim(text)   ! F2003 realloc
        call set( text //      '%Num_Probes       ',                        &
                  rake_probes(n)%num_probes        , 0                      )
        call set( text //      '%XYZ_Beg          ',                        &
                  rake_probes(n)%xyz_beg           , probes_xyz_default     )
        call set( text //      '%XYZ_End          ',                        &
                  rake_probes(n)%xyz_end           , probes_xyz_default     )
      end do
      !
      do n = 1, size(cell_probes)
        text = '                               '
        write (text,101) 'Cell_Probes(',n,')'
        text = trim(text)   ! F2003 realloc
        call set( text //      '%XYZ              ',                        &
                  cell_probes(n)%xyz               , probes_xyz_default     )
      end do
      !}}}
!{{{  ! TIME AVERAGING
      call set('Output_Time_Averaging             ',                        &
                output_time_averaging              , fals                   )
      call set('Time_Ave_File                     ',                        &
                time_ave_file                      , 'out.ave'              )
      call set('Time_Scaling_Factor               ',                        &
                time_scaling_factor                , 1.0_wp                 )
      call set('Time_Ave_Vel_Is_Axisymm           ',                        &
                time_ave_vel_is_axisymm            , fals                   )
      call set('Start_Averaging_Timestep          ',                        &
                start_averaging_timestep           ,-1                      )
      call set('Start_Averaging_Tstar             ',                        &
                start_averaging_tstar              ,-one                    )
      !}}}
!{{{  ! SPATIAL AVERAGE RESTART FILES
      call set('Spatial_Ave_Restart_X_And_Z_Planes',                        &
                spatial_ave_restart_x_and_z_planes , fals                   )
      call set('Spatial_Ave_Restart_Z_Planes      ',                        &
                spatial_ave_restart_z_planes       , fals                   )
      !}}}
!{{{  ! UTILITY TO TIME AVERAGE EXISTING RESTART FILES
      call set('Time_Average_Restart_Files        ',                        &
                time_average_restart_files         , fals                   )
      call set('Restart_List_File                 ',                        &
                restart_list_file                  , ''                     )
      call set('Modify_Time_Ave_Restart_History   ',                        &
                modify_time_ave_restart_history    , fals                   )
      !}}}
!{{{  ! MEMORY USAGE
      call set('Dump_Memory_Usage                 ',                        &
                dump_memory_usage                  , fals                   )
      !}}}
!{{{  ! DUMP FLUX INFORMATION FOR DEBUGGING
      call set('Dump_Fluxes                       ',                        &
                dump_fluxes                        , fals                   )
      call set('Check_Flux_Point_Coordinates      ',                        &
                check_flux_point_coordinates       , fals                   )
      !}}}
!{{{  ! ACTIVATE SPECIAL GRID FOR DEBUGGING THE POSTPROCESSOR
      call set('Postproc_Debug_Grid               ',                        &
                postproc_debug_grid                , fals                   )
      !}}}
!{{{  ! OPTION TO OUTPUT PLOT3D GRID TO TECPLOT IN UNSTRUCTURED FORMAT
      call set('Output_Plot3d_To_Tecplot          ',                        &
                output_plot3d_to_tecplot           , fals                   )
      !}}}
!{{{  ! OPTION TO OUTPUT CONTENTS OF BFACE ARRAY
      call set('Output_Bface_Array                ',                        &
                output_bface_array                 , fals                   )
      !}}}
!{{{  ! MISC
      call set('Force_Nd_Subroutines              ',                        &
                force_nd_subroutines               , fals                   )
      call set('Dmat_Grad_Temp_Order              ',                        &
                dmat_grad_temp_order               ,-1                      )
      call set('Newton_Rlx                        ',                        &
                newton_rlx                         , 1.0_wp                 )
      call set('Conv_Temp                         ',                        &
                conv_temp                          , eps12                  )
      call set('Chemistry_Sub_Iterations          ',                        &
                chemistry_sub_iterations           ,-1                      )
      call set('Ignore_Thermo_Newton_Warnings     ',                        &
                ignore_thermo_newton_warnings      , fals                   )
      !}}}
!{{{  ! INPUT BOUNDARY CONDITIONS
      !   bc_input(#) derived-type array
      !
      loop_bc_input : &
      do n = 1, size(bc_input)
        ! Must set text large enough to write into.
        text = '                               '
        write(text,101) 'BC_Input(',n,')'
        text = trim(text)   ! F2003 realloc
        call set( text //   '%Name                   ',                     &
                  bc_input(n)%name                    , ''                  )
        call set( text //   '%Set_BC_Value           ',                     &
                  bc_input(n)%set_bc_value            , bc_unknown          )
        call set( text //   '%T_Static               ',                     &
                  bc_input(n)%t_static                ,-one                 )
        call set( text //   '%P_Static               ',                     &
                  bc_input(n)%p_static                ,-one                 )
        call set( text //   '%Rho_Static             ',                     &
                  bc_input(n)%rho_static              ,-one                 )
        call set( text //   '%Mach                   ',                     &
                  bc_input(n)%mach                    ,-one                 )
        call set( text //   '%Vx                     ',                     &
                  bc_input(n)%vx                      , zero                )
        call set( text //   '%Vy                     ',                     &
                  bc_input(n)%vy                      , zero                )
        call set( text //   '%Vz                     ',                     &
                  bc_input(n)%vz                      , zero                )
        call set( text //   '%T_Total                ',                     &
                  bc_input(n)%t_total                 ,-one                 )
        call set( text //   '%P_Total                ',                     &
                  bc_input(n)%p_total                 ,-one                 )
        call set( text //   '%Rho_Total              ',                     &
                  bc_input(n)%rho_total               ,-one                 )
        call set( text //   '%Alpha_AoA              ',                     &
                  bc_input(n)%alpha_aoa               , zero                )
        call set( text //   '%Beta_AoA               ',                     &
                  bc_input(n)%beta_aoa                , zero                )
        call set( text //   '%Wall_Temp              ',                     &
                  bc_input(n)%wall_temp               ,-one                 )
        ! Was this once rotation_speed ?
        call set( text //   '%Rotation_Rate          ',                     &
                  bc_input(n)%rotation_rate           , [zero,zero,zero]    )
        ! Next item must have default logical type.
        call set( text //   '%Is_Specified_In_MRF'    ,                     &
                  bc_input(n)%is_specified_in_mrf     , .false.             )
        !
        loop_mass_frac : &
        do i = 1, size(bc_input(n)%mass_frac)
          text2 = '                                '
          write(text2,101) text // '%Mass_Frac(',i,')'
          text2 = trim(text2)   ! F2003 realloc
          call set( text2 //                   '      ',                    &
                    bc_input(n)%mass_frac(i)           , zero               )
        end do loop_mass_frac
        !
        loop_species_name : &
        do i = 1, size(bc_input(n)%species_name)
          text2 = '                                '
          write(text2,101) text // '%Species_Name(',i,')'
          text2 = trim(text2)   ! F2003 realloc
          call set( text2 //                      '   ',                    &
                    bc_input(n)%species_name(i)        , ''                 )
        end do loop_species_name
        !
        call set( text //   '%Relax%Value             ',                    &
                  bc_input(n)%relax%value              , one                )
        call set( text //   '%Relax%Time              ',                    &
                  bc_input(n)%relax%time               , zero               )
        call set( text //   '%Relax%Iter              ',                    &
                  bc_input(n)%relax%iter               , 0                  )
        !
        call set( text //   '%SEM%Is_Cylindrical      ',                    &
                  bc_input(n)%sem%is_cylindrical       , fals               )
        call set( text //   '%SEM%Length_Scale        ',                    &
                  bc_input(n)%sem%length_scale         , zero               )
        call set( text //   '%SEM%Convective_Velocity ',                    &
                  bc_input(n)%sem%convective_velocity  , zero               )
        call set( text //   '%SEM%Profile_Height      ',                    &
                  bc_input(n)%sem%profile_height       , zero               )
        call set( text //   '%SEM%Profile_File        ',                    &
                  bc_input(n)%sem%profile_file         , ''                 )
        call set( text //   '%SEM%Sub_In_BC_Grp       ',                    &
                  bc_input(n)%sem%sub_in_bc_grp        , 0                  )
        call set( text //   '%SEM%Sub_In_BC_Idx       ',                    &
                  bc_input(n)%sem%sub_in_bc_idx        , 0                  )
        !
        call set( text //   '%Periodic%Matching_BC_Name',                   &
                  bc_input(n)%periodic%matching_bc_name , ''                )
        call set( text //   '%Periodic%Matching_BC_Idx ',                   &
                  bc_input(n)%periodic%matching_bc_idx  , 0                 )
        !
        loop_rotation_center : &
        do i = 1, size(bc_input(n)%Periodic%RotationCenter)
          text2 = '          '
          write(text2,101) '(',i,')'
          text2 = text // '%Periodic%RotationCenter' // text2
          text2 = trim(text2)   ! F2003 realloc
          call set( text2 //                           ' ',                 &
                    bc_input(n)%Periodic%RotationCenter(i), zero            )
        end do loop_rotation_center
        !
        loop_rotation_angle : &
        do i = 1, size(bc_input(n)%Periodic%RotationAngle)
          text2 = '          '
          write(text2,101) '(',i,')'
          text2 = text // '%Periodic%RotationAngle' // text2
          text2 = trim(text2)   ! F2003 realloc
          call set( text2 //                         '  ' ,                 &
                    bc_input(n)%Periodic%RotationAngle(i) , zero            )
        end do loop_rotation_angle
        !
        loop_rotation_axis : &
        do i = 1, size(bc_input(n)%Periodic%RotationAxis)
          text2 = '          '      
          write(text2,101) '(',i,')'
          text2 = text // '%Periodic%RotationAxis' // text2
          text2 = trim(text2)   ! F2003 realloc
          call set( text2 //                       '   '  ,                 &
                    bc_input(n)%Periodic%RotationAxis(i)  , zero            )
        end do loop_rotation_axis
        !
        loop_translation : &
        do i = 1, size(bc_input(n)%Periodic%Translation)
          text2 = '          '
          write(text2,101) '(',i,')'
          text2 = text // '%Periodic%Translation' // text2
          text2 = trim(text2)   ! F2003 realloc
          call set( text2 //                     '    '   ,                 &
                    bc_input(n)%Periodic%Translation(i)   , zero            )
        end do loop_translation
        !
      end do loop_bc_input  ! n = 1, size(bc_input)
      !}}}
!{{{  ! GEOMETRY WEIGHTS FOR METIS PARTITIONING
      !   geom_wts(*) derived-type array
      !
      do n = 1, size(geom_wts)
        ! Must set text large enough to write into.
        text = '                               '
        write(text,101) 'Geom_Wts(',n,')'
        text = trim(text)
        call set( text //   '%Tria          ',geom_wts(n)%tria       , 1 )
        call set( text //   '%Quad          ',geom_wts(n)%quad       , 1 )
        call set( text //   '%Tetr          ',geom_wts(n)%tetr       , 1 )
        call set( text //   '%Pyra          ',geom_wts(n)%pyra       , 1 )
        call set( text //   '%Pris          ',geom_wts(n)%pris       , 1 )
        call set( text //   '%Hexa          ',geom_wts(n)%hexa       , 1 )
        call set( text //   '%Face_Tria     ',geom_wts(n)%face_tria  , 1 )
        call set( text //   '%Face_Quad     ',geom_wts(n)%face_quad  , 1 )
      end do
      !}}}
      call iow_footer()
    case ( NMLaction_Read )
      !
      ! Read in the namelist from the connected input file
      !
      if ( io_has_nml_block(ior,'&input') ) then
        read (ior,nml=input)
      end if
    case ( NMLaction_Write )
      !
      ! Have the compiler write out the namelist in the compilers default format
      !
      write (iow,nml=input,delim='APOSTROPHE')
  end select
  !
  ! Process the version 0001 namelist inputs to get the equivalent
  ! settings using the current GFR variables
  !
  call process_nml_data(ierr)
  !
  ! Convert the old variables to the new.
  !   Old: Limiter_Filter_Option          =>  New: Limiter_Filter_Option
  !   Old: Limiter_Option                     New: Limiter_Option
  !   Old: Limiter_Use_Dynamic_Filter
  !   Old: Resolution_Marker_Variable     =>  New: Marker_Option
  !   Old:      Jump_Marker_Variable          New: Marker_Variable
  !   Old:       MLP_Marker_Variable
  call convert_old_to_new(error2)
  !
  101 format(a,i0,a)
  !
end subroutine ionml_input
!
!  Note: subroutine ionml_input "contains":
!    subroutine convert_old_to_new
!    subroutine convert_new_to_old
!  in order to share scope with the following variables:
!    Limiter_Filter_Option
!    Limiter_Option
!    Limiter_Use_Dynamic_Filter
!    Resolution_Marker_Variable
!    Jump_Marker_Variable
!    MLP_Marker_Variable
!    Marker_Option
!    Marker_Variable
!
!###############################################################################
!> Convert:
!>   Old_Limiter_Filter_Option          =>  New_Limiter_Filter_Option
!>   Old_Limiter_Option                     New_Limiter_Option
!>   Old_Limiter_Use_Dynamic_Filter
!>   Old_Resolution_Marker_Variable     =>  New_Marker_Option
!>         Old_Jump_Marker_Variable         New_Marker_Variable
!>          Old_MLP_Marker_Variable
!> May encounter dummy values as old file is read twice with different defaults.
!>
!> The following old->new conversions should work.
!> Old variables               ->  New variables
!> LFO  LO  DF  Res  Jmp  MLP      LFO   LO  MO   MV
!> ---  -2  --  ---  0:5  ---        1 9002   2  0:5
!> ---  -1  --  0:5  ---  ---        1 9002   1  1:6
!> ---   0  --  0:5  ---  ---        1    0   1  1:6
!> ---   3  --  ---  ---  0:5        1    3   3  0:5
!>   0   1  --  0:5  ---  ---        1 9001   1  1:6
!>   0   2  --  ---  0:5  ---        1 9001   2  0:5
!> 1:2   1   0  0:5  ---  ---      1:2    1   1  1:6
!> 1:2   1   1  0:5  ---  ---      1:2    2   1  1:6
!> 1:2   2   0  ---  0:5  ---      1:2    1   2  0:5
!> 1:2   2   1  ---  0:5  ---      1:2    2   2  0:5
!
!###############################################################################
subroutine process_nml_data(error)
  !
  !.. Use Statements ..
  use ovar, only : process_marker_coefficients
  use nml_base_mod, only : idum0,idum1
  !
  !.. Formal Arguments ..
  integer, intent(out) :: error
  !
  ! Old variables
  integer :: Old_Limiter_Filter_Option
  integer :: Old_Limiter_Option
  logical :: Old_Limiter_Use_Dynamic_Filter
  integer :: Old_Resolution_Marker_Variable
  integer :: Old_Jump_Marker_Variable
  integer :: Old_MLP_Marker_Variable
  !
  ! New variables
  integer :: New_Limiter_Filter_Option
  integer :: New_Limiter_Option
  integer :: New_Marker_Option
  integer :: New_Marker_Variable
  !
continue
  !
  error = 0  ! ok
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! Import old variables
  Old_Limiter_Filter_Option      = Limiter_Filter_Option
  Old_Limiter_Option             = Limiter_Option
  Old_Limiter_Use_Dynamic_Filter = Limiter_Use_Dynamic_Filter
  Old_Resolution_Marker_Variable = Resolution_Marker_Variable
  Old_Jump_Marker_Variable       = Jump_Marker_Variable
  Old_MLP_Marker_Variable        = MLP_Marker_Variable
  !
  ! Convert:  Old_Limiter_Filter_Option -> New_Limiter_Filter_Option
  ! Note:     New_Limiter_Filter_Option = 0  not allowed, use 1 instead
  !
  New_Limiter_Filter_Option = Old_Limiter_Filter_Option
  if ( Old_Limiter_Filter_Option == 0 ) New_Limiter_Filter_Option = 1
  !
  ! Convert:  Old_Limiter_Filter_Option + Old_Limiter_Option,
  !           Old_Resolution_Marker_Variable, Old_Jump_Marker_Variable,
  !           Old_MLP_Marker_Variable
  !      ->   New_Limiter_Option, New_Marker_Option, New_Marker_Variable
  !
  select case (Old_Limiter_Option)
    case ( idum0, idum1)  ! None
      !
      ! Limiter_Option is not present, but has default=0 (None/Off).
      ! Treat marker as default=resolution, or not present?
      !
      New_Limiter_Option = Old_Limiter_Option  ! dummy value
!     New_Marker_Option  = Old_Limiter_Option  ! dummy value
      New_Marker_Option  = Marker_Option__Resolution
      !
    case ( 0)  ! None
      New_Limiter_Option = Limiter_Option__Off
      New_Marker_Option  = Marker_Option__Resolution
      !
    case (-2)  ! Jump Marker + Scaled Projection
      New_Limiter_Option = Limiter_Option__ProjectionScaled
      New_Marker_Option  = Marker_Option__Jump
      !
    case (-1)  ! Res. Marker + Scaled Projection
      New_Limiter_Option = Limiter_Option__ProjectionScaled
      New_Marker_Option  = Marker_Option__Resolution
      !
    case ( 1)  ! Res. Marker + Limiter_Filter_Option
      if ( Old_Limiter_Filter_Option == 0 ) then
        New_Limiter_Option = Limiter_Option__ProjectionRecursive
      else if ( Old_Limiter_Use_Dynamic_Filter ) then
        New_Limiter_Option = Limiter_Option__FilterScaled
      else
        New_Limiter_Option = Limiter_Option__Filter
      end if
      New_Marker_Option = Marker_Option__Resolution
      !
    case ( 2)  ! Jump Marker + Limiter_Filter_Option
      if ( Old_Limiter_Filter_Option == 0 ) then
        New_Limiter_Option = Limiter_Option__ProjectionRecursive
      else if ( Old_Limiter_Use_Dynamic_Filter ) then
        New_Limiter_Option = Limiter_Option__FilterScaled
      else
        New_Limiter_Option = Limiter_Option__Filter
      end if
      New_Marker_Option = Marker_Option__Jump
      !
    case ( 3)  ! MLP  Marker + MLP Limiter
      New_Limiter_Option = Limiter_Option__MLP
      New_Marker_Option  = Marker_Option__MLP
      !
    case ( 4)  ! Ducros/divergence Marker + Limiter_Filter_Option
      if ( Old_Limiter_Filter_Option == 0 ) then
        New_Limiter_Option = Limiter_Option__ProjectionRecursive
      else if ( Old_Limiter_Use_Dynamic_Filter ) then
        New_Limiter_Option = Limiter_Option__FilterScaled
      else
        New_Limiter_Option = Limiter_Option__Filter
      end if
      New_Marker_Option = Marker_Option__Ducros
      !
    case default
      ! Invalid value for Old_Limiter_Option.
      ! Unable to convert from old input.  Set values to off.
      New_Limiter_Option = Limiter_Option__Off
      New_Marker_Option  = Marker_Option__Off
      error = 1  ! fail
      !
  end select
  !
  ! Set:  New_Marker_Variable, adjust from old value to new.
  !
  select case (New_Marker_Option)
    case (idum0,idum1)
      New_Marker_Variable = New_Marker_Option  ! idum0,idum1
    case (Marker_Option__Off)
      New_Marker_Variable = Marker_Variable__Off
    case (Marker_Option__Resolution)
      if ( Old_Resolution_Marker_Variable == 0 ) then
        New_Marker_Variable = Marker_Variable__all
      else
        New_Marker_Variable = Old_Resolution_Marker_Variable
      end if
    case (Marker_Option__Jump)
      if ( Old_Jump_Marker_Variable == 0 ) then
        New_Marker_Variable = Marker_Variable__npres
      else
        New_Marker_Variable = Old_Jump_Marker_Variable
      end if
    case (Marker_Option__MLP)
      New_Marker_Variable = Old_MLP_Marker_Variable
    case (Marker_Option__Ducros)
      ! no variable option for Ducros/divergence marker
    case default
  end select
  !
  ! Export new variables
  Limiter_Filter_Option = New_Limiter_Filter_Option
  Limiter_Option        = New_Limiter_Option
  Marker_Option         = New_Marker_Option
  Marker_Variable       = New_Marker_Variable
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! Process the marker inputs to get the final marker
  ! coefficients that will be used to scale the raw marker
  ! values to the scaled_marker on the interval [0,1].
  !
  call process_marker_coefficients
  !
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! Get the sponge filter option
  !
  Sponge_Filter_Option = merge(1,0,Use_Sponge_Filter)  ! logical to int
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! Convert Dof inputs for specific testcases to generic internal_grid_dof
  !
  if (any(abs(itestcase) == Transport_Problems) .or. &
              itestcase  == Taylor_Green_Vortex) then
    internal_grid_dof = VortexGrid_DOF
  else
    select case (itestcase)
      case (CH4_Flame_Speed_Test) ; internal_grid_dof = flame_speed_test_dof
      case (Sod_Shock_Tube      ) ; internal_grid_dof = sod_test_dof
      case (Density_Pulse_Test  ) ; internal_grid_dof = pulse_test_dof
      case (Shu_Osher_Test      ) ; internal_grid_dof = shu_osher_test_dof
      case default
        ! use internal_grid_dof
    end select
  end if
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! Convert the bc_input structure from this version 1 reader to current.
  !   Old: bc_input (ver.1)  =>  New: bc_input (ovar_mod)
  !
  ovar_bc_input = bc_input_v0001_to_ovar( bc_input ) ! F2003 AUTO-REALLOC
  !
end subroutine process_nml_data
!
!###############################################################################
!> This subroutine preprocesses the marker coefficients before they are
!> actually processed for the solver when the NML action is to set dummy values,
!> because the dummy values will cause problems with the logic that computes
!> the final marker coefficients used by the solver.
!###############################################################################
subroutine preprocess_marker_coeff(p)
  !
  !.. Use Statements ..
  use ovar, only : Limiter_Marker_Min
  use ovar, only : Limiter_Marker_Max
  use ovar, only : Limiter_Marker_dpsi
  use ovar, only : Limiter_Marker_c1
  use ovar, only : Limiter_Marker_c2
  use ovar, only : Marker_MinMax_Default
  use ovar, only : Marker_Coeff_Default
  !
  use nml_base_mod, only : idum0,idum1
  use nml_base_mod, only : NMLaction_LatestDummy
  use nml_base_mod, only : NMLmode
  use nml_base_mod, only : NMLmode_InputConverter
  !
  !.. Formal Arguments ..
  integer,  allocatable, intent(out) :: p
  !
  !.. Local Scalars ..
  integer :: idummy_val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "preprocess_marker_coeff"
  !
continue
  !
  ! Only preprocess the marker coefficients if the current NML mode is running
  ! the input namelist converter. Otherwise return without doing anything
  !
  if (NMLmode /= NMLmode_InputConverter) return
  !
  select case (NMLaction)
    case (NMLaction_SetDummy0) ; idummy_val = idum0
    case (NMLaction_SetDummy1) ; idummy_val = idum1
    case (NMLaction_Read)
      if (NMLaction_LatestDummy == NMLaction_SetDummy1) then
        idummy_val = idum1
      else
        idummy_val = idum0
      end if
    case default
      ! Return without doing anything for all other NML actions:
      !   - Dont care about Write or WritePretty here.
      !   - I dont think BlockCheck will get here.
      !   - The SetValue action will be setting default values for the
      !     coefficients, so we can just return without doing anything
      !     since the defaults are compatible with the logic that
      !     processes the marker coefficients.
      return
  end select
  !
  ! Reset any of the marker coefficients that are set to the current dummy
  ! value back to their default value.
  !
  if (nint(Limiter_Marker_Min ) == idummy_val) then
    Limiter_Marker_Min = Marker_MinMax_Default
  end if
  if (nint(Limiter_Marker_Max ) == idummy_val) then
    Limiter_Marker_Max = Marker_MinMax_Default
  end if
  if (nint(Limiter_Marker_dpsi) == idummy_val) then
    Limiter_Marker_dpsi = Marker_Coeff_Default
  end if
  if (nint(Limiter_Marker_c1  ) == idummy_val) then
    Limiter_Marker_c1 = Marker_Coeff_Default
  end if
  if (nint(Limiter_Marker_c2  ) == idummy_val) then
    Limiter_Marker_c2 = Marker_Coeff_Default
  end if
  !
  ! Dont forget about the solution order. These will also cause problems when
  ! processing the marker coefficients if they are both set to the dummy value.
  ! If these are both the dummy value, allocate the scalar subroutine argument
  ! to the value 3 which is the default value for solution_order. This
  ! allocatable will be passed in as an optional argument to subroutine
  ! process_marker_coefficients, which will only be used in replace of
  ! solution_order and all_orders if it is present and allocated. It will not
  ! actually modify the value of either of those input variables.
  !
  if (all([solution_order,all_orders] == idummy_val)) then
    p = 3 ! F2003 AUTO-REALLOC
  end if
  !
end subroutine preprocess_marker_coeff
!
!###############################################################################
!> This subroutine postprocesses the marker coefficients after they are
!> actually processed for the solver when the NML action is to set dummy values.
!> This undoes any preprocessing done before the coefficients were processed
!> for the solver.
!###############################################################################
subroutine postprocess_marker_coeff()
  !
  !.. Use Statements ..
  use ovar, only : Limiter_Marker_Min
  use ovar, only : Limiter_Marker_Max
  use ovar, only : Limiter_Marker_dpsi
  use ovar, only : Limiter_Marker_c1
  use ovar, only : Limiter_Marker_c2
  use ovar, only : Marker_MinMax_Default
  use ovar, only : Marker_Coeff_Default
  !
  use nml_base_mod, only : rdum0,rdum1
  use nml_base_mod, only : NMLaction_LatestDummy
  use nml_base_mod, only : NMLmode
  use nml_base_mod, only : NMLmode_InputConverter
  !
  !.. Local Scalars ..
  real(wp) :: rdummy_val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "postprocess_marker_coeff"
  !
continue
  !
  ! Only postprocess the marker coefficients if the current NML mode is running
  ! the input namelist converter. Otherwise return without doing anything
  !
  if (NMLmode /= NMLmode_InputConverter) return
  !
  select case (NMLaction)
    case (NMLaction_SetDummy0) ; rdummy_val = rdum0
    case (NMLaction_SetDummy1) ; rdummy_val = rdum1
    case (NMLaction_Read)
      if (NMLaction_LatestDummy == NMLaction_SetDummy1) then
        rdummy_val = rdum1
      else
        rdummy_val = rdum0
      end if
    case default
      ! Return without doing anything for all other NML actions:
      !   - Dont care about Write or WritePretty here.
      !   - I dont think BlockCheck will get here.
      !   - The SetValue action will be setting default values for the
      !     coefficients, so we can just return without doing anything
      !     since the defaults are compatible with the logic that
      !     processes the marker coefficients.
      return
  end select
  !
  if (abs(Limiter_Marker_Min - Marker_MinMax_Default) < eps9) then
    Limiter_Marker_Min = rdummy_val
  end if
  if (abs(Limiter_Marker_Max - Marker_MinMax_Default) < eps9) then
    Limiter_Marker_Max = rdummy_val
  end if
  if (abs(Limiter_Marker_dpsi - Marker_Coeff_Default) < eps9) then
    Limiter_Marker_dpsi = rdummy_val
  end if
  if (abs(Limiter_Marker_c1 - Marker_Coeff_Default) < eps9) then
    Limiter_Marker_c1   = rdummy_val
  end if
  if (abs(Limiter_Marker_c2 - Marker_Coeff_Default) < eps9) then
    Limiter_Marker_c2   = rdummy_val
  end if
  !
end subroutine postprocess_marker_coeff
!
!###############################################################################
!
!  Note: subroutine ionml_input "contains":
!    subroutine convert_old_to_new
!    subroutine convert_new_to_old
!  in order to share scope with the following variables:
!    Limiter_Filter_Option
!    Limiter_Option
!    Limiter_Use_Dynamic_Filter
!    Resolution_Marker_Variable
!    Jump_Marker_Variable
!    MLP_Marker_Variable
!    Marker_Option
!    Marker_Variable
!
!###############################################################################
!> Convert:
!>   Old_Limiter_Filter_Option          =>  New_Limiter_Filter_Option
!>   Old_Limiter_Option                     New_Limiter_Option
!>   Old_Limiter_Use_Dynamic_Filter
!>   Old_Resolution_Marker_Variable     =>  New_Marker_Option
!>         Old_Jump_Marker_Variable         New_Marker_Variable
!>          Old_MLP_Marker_Variable
!> May encounter dummy values as old file is read twice with different defaults.
!>
!> The following old->new conversions should work.
!> Old variables               ->  New variables
!> LFO  LO  DF  Res  Jmp  MLP      LFO   LO  MO   MV
!> ---  -2  --  ---  0:5  ---        1 9002   2  0:5
!> ---  -1  --  0:5  ---  ---        1 9002   1  1:6
!> ---   0  --  0:5  ---  ---        1    0   1  1:6
!> ---   3  --  ---  ---  0:5        1    3   3  0:5
!>   0   1  --  0:5  ---  ---        1 9001   1  1:6
!>   0   2  --  ---  0:5  ---        1 9001   2  0:5
!> 1:2   1   0  0:5  ---  ---      1:2    1   1  1:6
!> 1:2   1   1  0:5  ---  ---      1:2    2   1  1:6
!> 1:2   2   0  ---  0:5  ---      1:2    1   2  0:5
!> 1:2   2   1  ---  0:5  ---      1:2    2   2  0:5
!
!###############################################################################
subroutine convert_old_to_new(error)
  !
  !.. Use Statements ..
  use nml_base_mod, only : idum0,idum1
  !
  !.. Formal Arguments ..
  integer, intent(out) :: error
  !
  ! Old variables
  integer :: Old_Limiter_Filter_Option
  integer :: Old_Limiter_Option
  logical :: Old_Limiter_Use_Dynamic_Filter
  integer :: Old_Resolution_Marker_Variable
  integer :: Old_Jump_Marker_Variable
  integer :: Old_MLP_Marker_Variable
  !
  ! New variables
  integer :: New_Limiter_Filter_Option
  integer :: New_Limiter_Option
  integer :: New_Marker_Option
  integer :: New_Marker_Variable
  !
continue
  !
  error = 0  ! ok
  !
  ! Import old variables
  Old_Limiter_Filter_Option      = Limiter_Filter_Option
  Old_Limiter_Option             = Limiter_Option
  Old_Limiter_Use_Dynamic_Filter = Limiter_Use_Dynamic_Filter
  Old_Resolution_Marker_Variable = Resolution_Marker_Variable
  Old_Jump_Marker_Variable       = Jump_Marker_Variable
  Old_MLP_Marker_Variable        = MLP_Marker_Variable
  !
  ! Convert:  Old_Limiter_Filter_Option -> New_Limiter_Filter_Option
  ! Note:     New_Limiter_Filter_Option = 0  not allowed, use 1 instead
  !
  New_Limiter_Filter_Option = Old_Limiter_Filter_Option
  if ( Old_Limiter_Filter_Option == 0 ) New_Limiter_Filter_Option = 1
  !
  ! Convert:  Old_Limiter_Filter_Option + Old_Limiter_Option,
  !           Old_Resolution_Marker_Variable, Old_Jump_Marker_Variable,
  !           Old_MLP_Marker_Variable
  !      ->   New_Limiter_Option, New_Marker_Option, New_Marker_Variable
  !
  select case (Old_Limiter_Option)
    case ( idum0, idum1)  ! None
      !
      ! Limiter_Option is not present, but has default=0 (None/Off).
      ! Treat marker as default=resolution, or not present?
      !
      New_Limiter_Option = Old_Limiter_Option  ! dummy value
!     New_Marker_Option  = Old_Limiter_Option  ! dummy value
      New_Marker_Option  = Marker_Option__Resolution
      !
    case ( 0)  ! None
      New_Limiter_Option = Limiter_Option__Off
      New_Marker_Option  = Marker_Option__Resolution
      !
    case (-2)  ! Jump Marker + Scaled Projection
      New_Limiter_Option = Limiter_Option__ProjectionScaled
      New_Marker_Option  = Marker_Option__Jump
      !
    case (-1)  ! Res. Marker + Scaled Projection
      New_Limiter_Option = Limiter_Option__ProjectionScaled
      New_Marker_Option  = Marker_Option__Resolution
      !
    case ( 1)  ! Res. Marker + Limiter_Filter_Option
      if ( Old_Limiter_Filter_Option == 0 ) then
        New_Limiter_Option = Limiter_Option__ProjectionRecursive
      else if ( Old_Limiter_Use_Dynamic_Filter ) then
        New_Limiter_Option = Limiter_Option__FilterScaled
      else
        New_Limiter_Option = Limiter_Option__Filter
      end if
      New_Marker_Option = Marker_Option__Resolution
      !
    case ( 2)  ! Jump Marker + Limiter_Filter_Option
      if ( Old_Limiter_Filter_Option == 0 ) then
        New_Limiter_Option = Limiter_Option__ProjectionRecursive
      else if ( Old_Limiter_Use_Dynamic_Filter ) then
        New_Limiter_Option = Limiter_Option__FilterScaled
      else
        New_Limiter_Option = Limiter_Option__Filter
      end if
      New_Marker_Option = Marker_Option__Jump
      !
    case ( 3)  ! MLP  Marker + MLP Limiter
      New_Limiter_Option = Limiter_Option__MLP
      New_Marker_Option  = Marker_Option__MLP
      !
    case ( 4)  ! Ducros/divergence Marker + Limiter_Filter_Option
      if ( Old_Limiter_Filter_Option == 0 ) then
        New_Limiter_Option = Limiter_Option__ProjectionRecursive
      else if ( Old_Limiter_Use_Dynamic_Filter ) then
        New_Limiter_Option = Limiter_Option__FilterScaled
      else
        New_Limiter_Option = Limiter_Option__Filter
      end if
      New_Marker_Option = Marker_Option__Ducros
      !
    case default
      ! Invalid value for Old_Limiter_Option.
      ! Unable to convert from old input.  Set values to off.
      New_Limiter_Option = Limiter_Option__Off
      New_Marker_Option  = Marker_Option__Off
      error = 1  ! fail
      !
  end select
  !
  ! Set:  New_Marker_Variable, adjust from old value to new.
  !
  select case (New_Marker_Option)
    case (idum0,idum1)
      New_Marker_Variable = New_Marker_Option  ! idum0,idum1
    case (Marker_Option__Off)
      New_Marker_Variable = Marker_Variable__Off
    case (Marker_Option__Resolution)
      if ( Old_Resolution_Marker_Variable == 0 ) then
        New_Marker_Variable = Marker_Variable__all
      else
        New_Marker_Variable = Old_Resolution_Marker_Variable
      end if
    case (Marker_Option__Jump)
      if ( Old_Jump_Marker_Variable == 0 ) then
        New_Marker_Variable = Marker_Variable__npres
      else
        New_Marker_Variable = Old_Jump_Marker_Variable
      end if
    case (Marker_Option__MLP)
      New_Marker_Variable = Old_MLP_Marker_Variable
    case (Marker_Option__Ducros)
      ! no variable option for Ducros/divergence marker
    case default
  end select
  !
  ! Export new variables
  Limiter_Filter_Option = New_Limiter_Filter_Option
  Limiter_Option        = New_Limiter_Option
  Marker_Option         = New_Marker_Option
  Marker_Variable       = New_Marker_Variable
  !
  ! Convert Dof inputs for specific testcases to generic internal_grid_dof
  !
  if (any(abs(itestcase) == Transport_Problems) .or. &
              itestcase  == Taylor_Green_Vortex) then
    internal_grid_dof = VortexGrid_DOF
  else
    select case (itestcase)
      case (CH4_Flame_Speed_Test) ; internal_grid_dof = flame_speed_test_dof
      case (Sod_Shock_Tube      ) ; internal_grid_dof = sod_test_dof
      case (Density_Pulse_Test  ) ; internal_grid_dof = pulse_test_dof
      case (Shu_Osher_Test      ) ; internal_grid_dof = shu_osher_test_dof
      case default
        ! use internal_grid_dof
    end select
  end if
  !
end subroutine convert_old_to_new
!==========================================================================
! Convert:
!   New_Limiter_Filter_Option  => Old_Limiter_Filter_Option
!   New_Limiter_Option            Old_Limiter_Option
!                                 Old_Limiter_Use_Dynamic_Filter
!   New_Marker_Option          => Old_Resolution_Marker_Variable
!   New_Marker_Variable                 Old_Jump_Marker_Variable
!                                        Old_MLP_Marker_Variable
!
! The following new->old conversion should work, unless marked "not supported".
! Some start/end with invalid new/old modes.
! New Variables       ->  Old variables
! LFO   LO  MO    MV      LFO   LO  DF  Res  Jmp  MLP
! --- 9001   1   1:6        0    1   0  0:5  ---  ---
! --- 9001   2   0:5        0    2   0  ---  0:5  ---
! --- 9002   1   1:6        1   -1   0  0:5  ---  ---
! --- 9002   2   0:5        1   -2   0  ---  0:5  ---
! 0:2    0   1   1:6      0:2    0   0  0:5  ---  ---
! 0:2    0   2  -1:5      0:2    0   0  ---  0:5  ---  Invalid Old.
! 0:2    0   3   0:5      0:2    0   0  ---  ---  0:5  Invalid Old.
! 0:2    0   4* -2,0      0:2    0   0  ---  --- -2,0  Not supported.
! 0:2    3   3   0:5      0:2    3   0  ---  ---  0:5
!   0*   1   1   1:6        1    1   0  0:5  ---  ---  Invalid New.
!   0*   1   2   0:5        1    2   0  ---  0:5  ---  Invalid New.
!   0*   2   1   1:6        1    1   1  0:5  ---  ---  Invalid New.
!   0*   2   2   0:5        1    2   1  ---  0:5  ---  Invalid New.
! 1:2    1   1   1:6      1:2    1   0  0:5  ---  ---
! 1:2    1   2   0:5      1:2    2   0  ---  0:5  ---
! 1:2    2   1   1:6      1:2    1   1  0:5  ---  ---
! 1:2    2   2   0:5      1:2    2   1  ---  0:5  ---
!
!==========================================================================
subroutine convert_new_to_old(error)
  !
  !.. Use Statements ..
  use nml_base_mod, only : idum0,idum1
  !
  !.. Formal Arguments ..
  integer, intent(out) :: error
  !
  !.. Local Scalars ..
  integer :: Default_Marker_Option, tmp_Marker_Option
  integer :: Default_Marker_Variable, tmp_Marker_Variable
  integer :: Old_Marker_Variable
  !
  ! Old variables
  integer :: Old_Limiter_Filter_Option
  integer :: Old_Limiter_Option
  logical :: Old_Limiter_Use_Dynamic_Filter
  integer :: Old_Resolution_Marker_Variable
  integer :: Old_Jump_Marker_Variable
  integer :: Old_MLP_Marker_Variable
  !
  ! New variables
  integer :: New_Limiter_Filter_Option
  integer :: New_Limiter_Option
  integer :: New_Marker_Option
  integer :: New_Marker_Variable
  !
continue
  !
  error = 0  ! ok
  !
  ! Import new variables
  New_Limiter_Filter_Option = Limiter_Filter_Option
  New_Limiter_Option        = Limiter_Option
  New_Marker_Option         = Marker_Option
  New_Marker_Variable       = Marker_Variable
  !
  ! Determine Default_Marker_Option from Limiter_Option.
  !
  select case (New_Limiter_Option)
    case (Limiter_Option__MLP)
      Default_Marker_Option = Marker_Option__MLP
    case default
      Default_Marker_Option = Marker_Option__Resolution
  end select
  !
  ! If Marker_Option==Off, use Default_Marker_Option instead.
  !
  tmp_Marker_Option = New_Marker_Option
  if ( New_Marker_Option == Marker_Option__Off ) then
    tmp_Marker_Option = Default_Marker_Option
  end if
  !
  ! Determine Default_Marker_Variable from Marker_Option.
  !
  select case (tmp_Marker_Option)
    case (Marker_Option__Resolution)
      Default_Marker_Variable = Marker_Variable__nec
    case (Marker_Option__Jump)
      Default_Marker_Variable = Marker_Variable__npres
    case (Marker_Option__MLP)
      Default_Marker_Variable = Marker_Variable__nec
    case (Marker_Option__Ducros)
      Default_Marker_Variable = Marker_Variable__Ducros
    case default
      Default_Marker_Variable = Marker_Variable__nec
  end select
  !
  ! If Marker_Variable==Off, use Default_Marker_Variable instead.
  !
  tmp_Marker_Variable = New_Marker_Variable
  if ( New_Marker_Variable == Marker_Variable__Off ) then
    tmp_Marker_Variable = Default_Marker_Variable
  end if
  !
  ! -----------------------------------
  ! Set Old_Resolution_Marker_Variable,
  !     Old_Jump_Marker_Variable,
  !     Old_MLP_Marker_Variable
  ! -----------------------------------
  !
  ! Determine the old value of the Marker_Variable.
  select case (tmp_Marker_Variable)
    case (idum0                 ) ; Old_Marker_Variable = idum0
    case (idum1                 ) ; Old_Marker_Variable = idum1
    case (Marker_Variable__off  ) ! Off or default
      select case (tmp_Marker_Option)
        case (Marker_Option__Resolution) ; Old_Marker_Variable =  1
        case (Marker_Option__Jump      ) ; Old_Marker_Variable =  0
        case (Marker_Option__MLP       ) ; Old_Marker_Variable =  1
        case default
          ! Cannot convert Marker_Variable to old code, so set to default value.
          error = 1  ! fail
      end select
    case (Marker_Variable__npres) ; Old_Marker_Variable = 0
    case (Marker_Variable__all  ) ; Old_Marker_Variable = 0
    case (Marker_Variable__nec  ) ; Old_Marker_Variable = 1
    case (Marker_Variable__nmx  ) ; Old_Marker_Variable = 2
    case (Marker_Variable__nmy  ) ; Old_Marker_Variable = 3
    case (Marker_Variable__nmz  ) ; Old_Marker_Variable = 4
    case (Marker_Variable__nee  ) ; Old_Marker_Variable = 5
    case default                  ; Old_Marker_Variable = 0
      ! Usupported Marker_Variable.
      ! Cannot convert Marker_Variable to old input, so set to default value.
      error = 2  ! fail
  end select
  !
  ! Set default for Old_*_Marker_Variable.
  Old_Resolution_Marker_Variable = 0  ! default, meaning all conserved
        Old_Jump_Marker_Variable = 0  ! default, meaning pressure
         Old_MLP_Marker_Variable = 0  ! default, meaning density
  !
  ! Determine which old Marker_Variable to assign the value to.
  ! NOTE: No compatibility check between Marker_Variable and Marker_Option.
  !
  ! If Marker_Option has dummy value, set all *_Marker_Variable to same dummy.
  ! Otherwise set individual *_Marker_Variable to current value.
  !
  if ( any(New_Marker_Option == [idum0,idum1]) ) then
    ! Set to the same dummy value idum0,idum1.
    Old_Resolution_Marker_Variable = New_Marker_Option
          Old_Jump_Marker_Variable = New_Marker_Option
           Old_MLP_Marker_Variable = New_Marker_Option
    !
  else if ( tmp_Marker_Option == Marker_Option__Resolution ) then
    Old_Resolution_Marker_Variable = Old_Marker_Variable
    !
  else if ( tmp_Marker_Option == Marker_Option__Jump ) then
    Old_Jump_Marker_Variable       = Old_Marker_Variable
    !
  else if ( tmp_Marker_Option == Marker_Option__MLP ) then
    Old_MLP_Marker_Variable        = Old_Marker_Variable
    !
  else
    ! Unsupported Marker_Option.
    ! Cannot convert Marker_Option to old input.
    error = 3  ! fail
  end if
  !
  ! Set default for Old_Limiter_Filter_Option.  Adjusted below.
  ! For most cases, except limit by projection noted below.
  !
  Old_Limiter_Filter_Option = New_Limiter_Filter_Option
  !
  ! Set default for Old_Limiter_Use_Dynamic_Filter.  Adjusted below.
  !
  Old_Limiter_Use_Dynamic_Filter = .false.
  !
  !
  ! Set:    Old_Limiter_Option
  ! Adjust: Old_Limiter_Filter_Option, Old_Limiter_Use_Dynamic_Filter.
  ! These are the known combinations that can be back-ported.
  !
  Old_Limiter_Option = 0  ! default
  if      ( New_Limiter_Option == idum0 ) then
    Old_Limiter_Option = idum0  ! dummy value
    !
  else if ( New_Limiter_Option == idum1 ) then
    Old_Limiter_Option = idum1  ! dummy value
    !
  else if ( New_Limiter_Option == Limiter_Option__Off ) then
    Old_Limiter_Option =  0  ! Off
    !
  else if ( New_Limiter_Option == Limiter_Option__ProjectionScaled &
    .and.   tmp_Marker_Option  ==  Marker_Option__Jump ) then
    Old_Limiter_Option = -2 ! Jump Marker + Scaled Projection
    Old_Limiter_Filter_Option = 0
    !
  else if ( New_Limiter_Option == Limiter_Option__ProjectionScaled &
    .and.   tmp_Marker_Option  ==  Marker_Option__Resolution ) then
    Old_Limiter_Option = -1 ! Resolution Marker + Scaled Projection
    Old_Limiter_Filter_Option = 0
    !
  else if ( New_Limiter_Option == Limiter_Option__ProjectionRecursive &
    .and.   tmp_Marker_Option  ==  Marker_Option__Resolution ) then
    Old_Limiter_Option =  1 ! Resolution Marker + Recursive Projection
    Old_Limiter_Filter_Option = 0
    !
  else if ( New_Limiter_Option == Limiter_Option__ProjectionRecursive &
    .and.   tmp_Marker_Option  ==  Marker_Option__Jump ) then
    Old_Limiter_Option =  2 ! Jump Marker + Recursive Projection
    Old_Limiter_Filter_Option = 0
    !
  else if ( New_Limiter_Option == Limiter_Option__Filter .and. &
            tmp_Marker_Option  ==  Marker_Option__Resolution ) then
    Old_Limiter_Option =  1 ! Resolution Marker + Filter
    Old_Limiter_Use_Dynamic_Filter = fals
    !
  else if ( New_Limiter_Option == Limiter_Option__FilterScaled .and. &
            tmp_Marker_Option  ==  Marker_Option__Resolution ) then
    Old_Limiter_Option =  1 ! Resolution Marker + Filter
    Old_Limiter_Use_Dynamic_Filter = true
    !
  else if ( New_Limiter_Option == Limiter_Option__Filter .and. &
            tmp_Marker_Option  ==  Marker_Option__Jump ) then
    Old_Limiter_Option =  2 ! Jump Marker + Filter
    Old_Limiter_Use_Dynamic_Filter = fals
    !
  else if ( New_Limiter_Option == Limiter_Option__FilterScaled .and. &
            tmp_Marker_Option  ==  Marker_Option__Jump ) then
    Old_Limiter_Option =  2 ! Jump Marker + Scaled Filter
    Old_Limiter_Use_Dynamic_Filter = true
    !
  else if ( New_Limiter_Option == Limiter_Option__MLP .and. &
            tmp_Marker_Option  ==  Marker_Option__MLP ) then
    Old_Limiter_Option =  3 ! MLP Marker + MLP Limiter
    !
  else
    ! This combination of options is not supported in the old input format.
    ! Cannot convert Limiter_Option to old input.
    error = 4  ! fail
  end if
  !
  ! Export old variables
  Limiter_Filter_Option      = Old_Limiter_Filter_Option
  Limiter_Option             = Old_Limiter_Option
  Limiter_Use_Dynamic_Filter = Old_Limiter_Use_Dynamic_Filter
  Resolution_Marker_Variable = Old_Resolution_Marker_Variable
  Jump_Marker_Variable       = Old_Jump_Marker_Variable
  MLP_Marker_Variable        = Old_MLP_Marker_Variable
  !
end subroutine convert_new_to_old
!
!###############################################################################
!
!
!###############################################################################
!
end module nml_0001_mod
!
!###############################################################################
!
