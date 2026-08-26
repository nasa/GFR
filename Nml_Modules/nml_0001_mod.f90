!###############################################################################
!> Module for version 0001 of namelist io.
!###############################################################################
module nml_0001_mod
  !
  !.. Use Statements ..
  use module_kind_types
  !
  use bc_input_0001_mod, only : bc_input,bc_input_t
  !
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
  use ovar, only : channel_body_force,channel_init_file
  use ovar, only : channel_spatial_ave_restart
  use ovar, only : channel_nondim_output_coord
  use ovar, only : channel_yplus_of_one, channel_rake_save_interval
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
  !
  ! ###########################
  ! ###  MODULE PARAMETERS  ###
  ! ###########################
  !
  integer, public, parameter :: namelist_version = 0001
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
  type(nml_0001_procedures_t), public, save :: nml_ver_0001
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
  ! ########################
  ! ###  Input Namelist  ###
  ! ########################
  !
  private :: input
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
  namelist / input / channel_body_force, channel_init_file, &
                     channel_spatial_ave_restart, channel_nondim_output_coord, &
                     channel_grid_ave_xplus, channel_grid_ave_zplus, &
                     channel_grid_ave_yplus, channel_grid_yplus_wall, &
                     channel_yplus_of_one, channel_rake_save_interval
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
contains
!
!###############################################################################
!
subroutine init_nml0001_procedures()
  !
  !.. Use Statements ..
  use nml_9999_mod, only : set_current_nml_version
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_nml0001_procedures"
  !
continue
  !
  nml_ver_0001%alloc_data => alloc_nml_data
  nml_ver_0001%resiz_data => resiz_nml_data
  nml_ver_0001%  set_data =>   set_nml_data
  nml_ver_0001% read_data =>  read_nml_data
  nml_ver_0001%write_data => write_nml_data
  nml_ver_0001% read_file =>  read_nml_file
  nml_ver_0001%write_file => write_nml_file
  !
  call set_current_nml_version(namelist_version,nml_ver_0001)
  !
end subroutine init_nml0001_procedures
!
!###############################################################################
!> Allocate namelist variables that are array of character.
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
  character(len=20) :: blanks
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
  ! Make sure everything is deallocated first
  !
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
  if (allocated(gridfile)) then
    deallocate ( gridfile , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"gridfile",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(plot3d_cutfile)) then
    deallocate ( plot3d_cutfile , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"plot3d_cutfile",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(plot3d_bcsfile)) then
    deallocate ( plot3d_bcsfile , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"plot3d_bcsfile",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(interpto_gridfile)) then
    deallocate ( interpto_gridfile , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"interpto_gridfile",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(interpto_plot3d_cutfile)) then
    deallocate ( interpto_plot3d_cutfile , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"interpto_plot3d_cutfile",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(interpto_plot3d_bcsfile)) then
    deallocate ( interpto_plot3d_bcsfile , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"interpto_plot3d_bcsfile",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(restart_file)) then
    deallocate ( restart_file , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"restart_file",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(solution_file)) then
    deallocate ( solution_file , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"solution_file",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(output_dir)) then
    deallocate ( output_dir , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"output_dir",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(channel_init_file)) then
    deallocate ( channel_init_file , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"channel_init_file",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(lustre_stripe_size)) then
    deallocate ( lustre_stripe_size , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"lustre_stripe_size",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(time_ave_file)) then
    deallocate ( time_ave_file , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"time_ave_file",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(restart_list_file)) then
    deallocate ( restart_list_file , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"restart_list_file",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  ! Allocate bc_input before reading the namelist
  !
  allocate ( bc_input(1:n_bc_input) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bc_input",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Allocate bc_input component arrays before reading the namelist:
  !   bcinput(:)%mass_frac(:), bcinput(:)%species_name(:)
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
  ovar_bc_input = bc_input(:)%copyout() ! F2003 AUTO-REALLOC
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
  ! Allocate all the character input variables
  !
  allocate ( character(len=default_len) :: gridfile , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"gridfile",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( character(len=default_len) :: plot3d_cutfile , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"plot3d_cutfile",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  allocate ( character(len=default_len) :: plot3d_bcsfile , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"plot3d_bcsfile",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  allocate ( character(len=default_len) :: interpto_gridfile , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"interpto_gridfile",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  allocate ( character(len=default_len) :: interpto_plot3d_cutfile , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"interpto_plot3d_cutfile",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  allocate ( character(len=default_len) :: interpto_plot3d_bcsfile , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"interpto_plot3d_bcsfile",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  allocate ( character(len=default_len) :: restart_file , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"restart_file",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( character(len=default_len) :: solution_file , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"solution_file",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( character(len=default_len) :: output_dir , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"output_dir",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( character(len=default_len) :: channel_init_file , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"channel_init_file",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  allocate ( character(len=default_len) :: lustre_stripe_size , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"lustre_stripe_size",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  allocate ( character(len=default_len) :: time_ave_file , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"time_ave_file",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( character(len=default_len) :: restart_list_file , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"restart_list_file",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  ! Set character input variables to empty string of specified length.
  !
  if (default_len > 0) then
    !
    write (blanks,1) default_len
    !
    write (gridfile                , blanks)
    write (plot3d_cutfile          , blanks)
    write (plot3d_bcsfile          , blanks)
    write (interpto_gridfile       , blanks)
    write (interpto_plot3d_cutfile , blanks)
    write (interpto_plot3d_bcsfile , blanks)
    write (restart_file            , blanks)
    write (solution_file           , blanks)
    write (output_dir              , blanks)
    write (channel_init_file       , blanks)
    write (lustre_stripe_size      , blanks)
    write (time_ave_file           , blanks)
    write (restart_list_file       , blanks)
    !
  end if
  !
  ! Set all the character input variables to their default values.
  ! For long  strings, preserve trailing blanks.
  ! For short strings, resize to minimal length for content.
  !
  if (len(restart_file) >= 7) then
    restart_file(1:7) = "out.rst"
  else
    restart_file = "out.rst" ! F2003 AUTO-REALLOC
  end if
  !
  if (len(output_dir) >= 1) then
    output_dir(1:1) = "."
  else
    output_dir = "." ! F2003 AUTO-REALLOC
  end if
  !
  if (len(lustre_stripe_size) >= 2) then
    lustre_stripe_size(1:2) = "4m"
  else
    lustre_stripe_size = "4m" ! F2003 AUTO-REALLOC
  end if
  !
  if (len(time_ave_file) >= 7) then
    time_ave_file(1:7) = "out.ave"
  else
    time_ave_file = "out.ave" ! F2003 AUTO-REALLOC
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  1 format ('(',i0,'x)')
  2 format ("bc_input(",i0,")%",a)
  !
end subroutine alloc_nml_data
!
!###############################################################################
!> Resize namelist allocatables to minimally necessary size for content.
!###############################################################################
subroutine resiz_nml_data
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
  ! ------------------------------------------------------------
  ! Reallocate bc_input to the actual size needed for the inputs
  ! ------------------------------------------------------------
  !
  ! Count how many of bc_input(:) are not default values.
  n = 0
  do i = 1,size(bc_input)
    if (bc_input(i)%is_not_default()) then
      n = n + 1
    end if
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
  ! Copy bc_input into ovar_bc_input
  !
  ovar_bc_input = bc_input(:)%copyout() ! F2003 AUTO-REALLOC
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
  channel_init_file       = trim(adjustl( channel_init_file       ))
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
    setverbose = true
  end if
  !
  setverbose = fals
  !
end subroutine set_nml_data
!
!###############################################################################
!> Read the namelist input data from the specified input unit.
!###############################################################################
subroutine read_nml_data(io)
  !
  !.. Arguments ..
  integer, intent(in) :: io  !< File unit number.
  !
  !.. Local Scalars ..
  integer :: ierr
  !
continue
  !
  read (io,nml=input)
  !
  ! Process the version 0001 namelist inputs to get the equivalent
  ! settings using the current GFR variables
  !
  call process_nml_data(ierr)
  !
end subroutine read_nml_data
!
!###############################################################################
!> Write the namelist input data to the specified output unit.
!###############################################################################
subroutine write_nml_data(action,io,banner)
  !
  !.. Arguments ..
  integer,     intent(in) :: action  !< An NMLaction
  integer,     intent(in) :: io      !< File unit number.
  logical(lk), intent(in) :: banner  !< Flag to write banner.
  !
continue
  !
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  if (mypnum == glb_root) then
    if (banner) write (io,104)
    write (io,nml=input,delim='apostrophe')
    if (banner) write (io,105)
    flush (io)
  end if
  !
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  104 format (/,"#######################################",/, &
                "#####      INPUT PARAMETERS       #####",/, &
                "#######################################",/)
  105 format (/,"#######################################",/, &
                "#####   END OF INPUT PARAMETERS   #####",/, &
                "#######################################",/)
  !
end subroutine write_nml_data
!
!###############################################################################
!
subroutine read_nml_file(input_fname)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: input_fname
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_nml_file"
  !
  !.. Local Scalars ..
  integer :: inpt,ierr
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Allocate all the variables within the input namelist
  ! to their default sizes and values
  !
  call alloc_nml_data
  !
  ! Open the input file
  !
  call debug_timer(start_timer,"actually reading the input file")
  open (newunit=inpt,file=input_fname,status="old",action="read", &
                     form="formatted",iostat=ierr,iomsg=error_message)
  call io_error(pname,input_fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Read in the input file using the input namelist
  !
  call read_nml_data(inpt)
  !
  ! Close the input file now that we are finished reading it
  !
  close (inpt,iostat=ierr,iomsg=error_message)
  call io_error(pname,input_fname,2,__LINE__,__FILE__,ierr,error_message)
  call debug_timer(stop_timer,"actually reading the input file")
  !
  ! Resize all the allocatable variables within the input
  ! namelist to sizes more appropriate to their input values
  !
  call resiz_nml_data
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_nml_file
!
!###############################################################################
!
subroutine write_nml_file(output_fname,action)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: output_fname
  integer,          intent(in) :: action        !< An NMLaction
  !
  !.. Local Scalars ..
  integer ionml,ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_nml_file"
  !
continue
  !
  ! Allocate all the variables within the input namelist
  ! to their default sizes and values
  !
  call alloc_nml_data(opt_n_bc_input=4, &
                      opt_n_species=3, &
                      opt_char_len=0, &
                      opt_n_point_probes=2, &
                      opt_n_rake_probes=1, &
                      opt_n_cell_probes=1)
  !
  if (mypnum == glb_root) then
    open (newunit=ionml,file=output_fname,status="new",action="write", &
                        delim='apostrophe', &
                        position="rewind",iostat=ierr,iomsg=error_message)
    call io_error(pname,trim(adjustl(output_fname)),1,__LINE__,__FILE__, &
                  ierr,error_message)
    call write_nml_data(0,ionml,banner=fals)
    write (iout,101) trim(output_fname)
    close(ionml)
  end if
  !
  101 format(/,"Default namelist input file written to:  ",A)
  !
end subroutine write_nml_file
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
  use ovar, only : marker_c1
  use ovar, only : marker_c2
  use ovar, only : marker_dpsi
  use ovar, only : Default_Marker_c1
  use ovar, only : Default_Marker_c1_JI
  use ovar, only : Default_Marker_c1_RI
  use ovar, only : Default_Marker_c2
  use ovar, only : Default_Marker_c2_JI
  use ovar, only : Default_Marker_c2_RI
  !
  !.. Formal Arguments ..
  integer, intent(out) :: error
  !
  !.. Local Scalars ..
  integer  :: n_order
  real(wp) :: l10p,max_marker_value
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
  ! Complete the marker coefficients that will be used to scale the raw
  ! marker values to the scaled_marker on the interval [0,1]
  !
  n_order = abs(solution_order)
  if ( all_orders >= 0 ) n_order = all_orders
  !
  if ( Marker_Option == Marker_Option__Jump ) then  ! Filter with jump marker.
    Default_Marker_c1 = Default_Marker_c1_JI
    Default_Marker_c2 = Default_Marker_c2_JI
    max_marker_value = log10(two)
  else ! resolution and Ducros/divergence markers
    Default_Marker_c1 = Default_Marker_c1_RI
    Default_Marker_c2 = Default_Marker_c2_RI
    max_marker_value = zero
  end if
  marker_c1 = Default_Marker_c1
  marker_c2 = Default_Marker_c2
  !
  ! DEBUG: Verbosely write out input marker coefficients before processing them
  call write_marker_inputs(1,n_order)
  !
  if (Limiter_Marker_c1   > zero) marker_c1   = Limiter_Marker_c1
  if (Limiter_Marker_c2   > zero) marker_c2   = Limiter_Marker_c2
  if (Limiter_Marker_dpsi > zero) marker_dpsi = Limiter_Marker_dpsi
  !
  ! Check if only one of Limiter_Marker_Min and Limiter_Marker_Max was
  ! set, and if so, compute the other using marker_dpsi
  !
  if (Limiter_Marker_Max <= max_marker_value .and. &
      Limiter_Marker_Min >= Limiter_Marker_Max) then
    Limiter_Marker_Min = Limiter_Marker_Max - two*marker_dpsi
  else if (Limiter_Marker_Min <= max_marker_value .and. &
           Limiter_Marker_Max >  max_marker_value) then
    Limiter_Marker_Max = min(max_marker_value, &
                             Limiter_Marker_Min + two*marker_dpsi)
  end if
  !
  l10p = log10(real(max(1,n_order),kind=wp))
  !
  ! If Limiter_Marker_Min and Limiter_Marker_Max are both less than the
  ! maximum possible value, recompute c1, c2, and dpsi.
  !
  if (all([Limiter_Marker_Min,Limiter_Marker_Max] <= max_marker_value)) then
    !
    marker_dpsi = half*abs(Limiter_Marker_Min-Limiter_Marker_Max)
    !
    ! Confirmed via Mathematica that the following options:
    !   1) holding marker_c1 constant and using it plus Limiter_Marker_min
    !      and marker_dpsi to compute a new marker_c2
    !   2) holding marker_c2 constant and using it plus Limiter_Marker_min
    !      and marker_dpsi to compute a new marker_c1
    ! result in the SAME psi0, which makes them equivalent.
    !
    if (n_order <= 1) then
      ! We have to compute marker_c1 for P<=1 because l10p will be zero
      marker_c1 = -one*(Limiter_Marker_Min + marker_dpsi + marker_c2*l10p)
    else
      if (marker_c1 /= Default_Marker_c1) then
        marker_c2 = -one*(Limiter_Marker_Min + marker_dpsi + marker_c1) / l10p
      else
        marker_c1 = -one*(Limiter_Marker_Min + marker_dpsi + marker_c2*l10p)
      end if
    end if
    !
  end if
  !
  ! Copy the computed values back into the input parameters.
  !
  Limiter_Marker_c1   = marker_c1
  Limiter_Marker_c2   = marker_c2
  Limiter_Marker_dpsi = marker_dpsi
  Limiter_Marker_Min  = -one*(marker_c1 + marker_c2*l10p) - marker_dpsi
  Limiter_Marker_Max  = -one*(marker_c1 + marker_c2*l10p) + marker_dpsi
  !
  ! DEBUG: Verbosely write out input marker coefficients after processing them
  call write_marker_inputs(2,n_order)
  !
  ! DEBUG: Write out the marker coefficients after processing them
  call write_marker_inputs(3,n_order)
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
  ovar_bc_input = bc_input(:)%copyout() ! F2003 AUTO-REALLOC
  !
end subroutine process_nml_data
!
!###############################################################################
!
subroutine write_marker_inputs(flag,p)
  !
  !.. Use Statements ..
  use ovar, only : marker_c1
  use ovar, only : marker_c2
  use ovar, only : marker_dpsi
  !
  !.. Formal Arguments ..
  integer, intent(in) :: flag
  integer, intent(in) :: p
  !
  !.. Local Scalars ..
  real(wp) :: psi0
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: text
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_marker_inputs"
  !
continue
#ifdef DUMP_MARKER_INPUTS
  !
  psi0 = -one*(marker_c1 + marker_c2*log10(real(max(1,p),kind=wp)))
  !
#if DUMP_MARKER_INPUTS == 1
  if (mypnum == glb_root) then
    if (flag == 3) then
      write (iout,1) marker_c1,marker_c2,marker_dpsi, &
                     psi0-marker_dpsi,psi0,psi0+marker_dpsi
      flush (iout)
    end if
  end if
#endif
#if DUMP_MARKER_INPUTS == 2
  if (mypnum == glb_root) then
    if (any(flag == [1,2])) then
      if (flag == 1) then
        text = "BEFORE" ! F2003 AUTO-REALLOC
      else
        text = "AFTER"  ! F2003 AUTO-REALLOC
      end if
      write (iout,2) text, &
                     Limiter_Marker_c1, &
                     Limiter_Marker_c2, &
                     Limiter_Marker_dpsi, &
                     Limiter_Marker_Min, &
                     Limiter_Marker_Max, &
                     marker_c1, &
                     marker_c2, &
                     marker_dpsi, &
                     psi0, &
                     psi0-marker_dpsi, &
                     psi0+marker_dpsi
      flush (iout)
    end if
  end if
#endif
#ifdef DUMP_MARKER_INPUTS_THEN_STOP
  flush (iout)
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,"Dumping marker parameters!")
#endif
  !
  1 format (/,2x,"marker_c1   = ",f11.6,/, &
              2x,"marker_c2   = ",f11.6,/, &
              2x,"marker_dpsi = ",f11.6,//, &
              2x,"psi0 - dpsi = ",f11.6,/, &
              2x,"psi0        = ",f11.6,/, &
              2x,"psi0 + dpsi = ",f11.6,/)
  !
  2 format (//,2x,a," modifying marker parameters:",//, &
               4x,"Limiter_Marker_c1   = ",f10.6,/, &
               4x,"Limiter_Marker_c2   = ",f10.6,/, &
               4x,"Limiter_Marker_dpsi = ",f10.6,/, &
               4x,"Limiter_Marker_Min  = ",f10.6,/, &
               4x,"Limiter_Marker_Max  = ",f10.6,//, &
               4x,"        marker_c1   = ",f10.6,/, &
               4x,"        marker_c2   = ",f10.6,/, &
               4x,"        marker_dpsi = ",f10.6,/, &
               4x,"        marker_psi0 = ",f10.6,/, &
               4x,"        marker_min  = ",f10.6,/, &
               4x,"        marker_max  = ",f10.6,//)
  !
#endif
end subroutine write_marker_inputs
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
  !.. Local Scalars ..
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
