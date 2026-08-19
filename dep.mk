$(OBJDIR)/main.o : main.f90 \
         $(OBJDIR)/kind_types_mod.o \
         $(OBJDIR)/eqn_idx_mod.o \
         $(OBJDIR)/pbs_mod.o \
         $(OBJDIR)/order_mod.o \
         $(OBJDIR)/geomvar_mod.o \
         $(OBJDIR)/ovar_mod.o \
         $(OBJDIR)/flowvar_mod.o \
         $(OBJDIR)/element_type_mod.o \
         $(OBJDIR)/quadrature_mod.o \
         $(OBJDIR)/vandermonde_mod.o \
         $(OBJDIR)/derivatives_mod.o \
         $(OBJDIR)/correction_mod.o \
         $(OBJDIR)/interpolation_mod.o \
         $(OBJDIR)/projection_mod.o \
         $(OBJDIR)/filter_mod.o \
         $(OBJDIR)/io_mod.o \
         $(OBJDIR)/periodic_mod.o \
         $(OBJDIR)/write_cgns_mod.o \
         $(OBJDIR)/metrics_mod.o \
         $(OBJDIR)/mappings_mod.o \
         $(OBJDIR)/connectivity_mod.o \
         $(OBJDIR)/initialization_mod.o \
         $(OBJDIR)/generic_mod.o \
         $(OBJDIR)/bc_mod.o \
         $(OBJDIR)/bnd_profiles_mod.o \
         $(OBJDIR)/flux_mod.o \
         $(OBJDIR)/parallel_mod.o \
         $(OBJDIR)/dt_mod.o \
         $(OBJDIR)/erk_mod.o \
         $(OBJDIR)/rhs_mod.o \
         $(OBJDIR)/oerk_mod.o \
         $(OBJDIR)/eperk_mod.o \
         $(OBJDIR)/channel_mod.o \
         $(OBJDIR)/limiters_mod.o \
         $(OBJDIR)/postproc_mod.o \
         $(OBJDIR)/probes_mod.o \
         $(OBJDIR)/restart_mod.o \
         $(OBJDIR)/plot3d_mod.o \
         $(OBJDIR)/gmsh_mod.o \
         $(OBJDIR)/cgns_mod.o \
         $(OBJDIR)/averaging_mod.o \
         $(OBJDIR)/reactions_mod.o \
         $(OBJDIR)/mrf_mod.o \
         $(OBJDIR)/spatial_ave_mod.o \
         $(OBJDIR)/synthetic_eddy_mod.o \
         $(OBJDIR)/cpu_info_mod.o \
         $(OBJDIR)/interpsol_mod.o \
         $(OBJDIR)/timer_mod.o \
         $(OBJDIR)/interpsol_query_mod.o \
         $(FUNDIR)/debug_timer.f90

$(OBJDIR)/averaging_mod.o : averaging_mod.f90 \
                  $(OBJDIR)/kind_types_mod.o \
                  $(OBJDIR)/eqn_idx_mod.o \
                  $(OBJDIR)/geomvar_mod.o \
                  $(OBJDIR)/ovar_mod.o \
                  $(OBJDIR)/flowvar_mod.o \
                  $(OBJDIR)/write_cgns_mod.o \
                  $(OBJDIR)/restart_mod.o \
                  $(OBJDIR)/mrf_mod.o

$(OBJDIR)/bc_mod.o : bc_mod.f90 \
           $(OBJDIR)/kind_types_mod.o \
           $(OBJDIR)/eqn_idx_mod.o \
           $(OBJDIR)/geomvar_mod.o \
           $(OBJDIR)/ovar_mod.o \
           $(OBJDIR)/flowvar_mod.o \
           $(OBJDIR)/metrics_mod.o \
           $(OBJDIR)/thermo_mod.o \
           $(OBJDIR)/transport_mod.o \
           $(OBJDIR)/synthetic_eddy_mod.o \
           $(OBJDIR)/timer_mod.o \
           $(OBJDIR)/bc_wall_mod.o \
           $(OBJDIR)/bc_fixed_mod.o \
           $(OBJDIR)/bc_flowthrough_grad_mod.o \
           $(OBJDIR)/bc_gen_free_mod.o \
           $(OBJDIR)/bc_gen_in_mod.o \
           $(OBJDIR)/bc_gen_out_mod.o \
           $(OBJDIR)/bc_periodic_mod.o \
           $(OBJDIR)/bc_riemann_mod.o \
           $(OBJDIR)/bc_sub_in_mod.o \
           $(OBJDIR)/bc_sub_out_mod.o \
           $(OBJDIR)/bc_sup_out_mod.o \
                     mpi_defs.h

$(OBJDIR)/bc_flow_stats_mod.o : bc_flow_stats_mod.f90 \
                      $(OBJDIR)/kind_types_mod.o \
                      $(OBJDIR)/eqn_idx_mod.o \
                      $(OBJDIR)/geomvar_mod.o \
                      $(OBJDIR)/ovar_mod.o \
                      $(OBJDIR)/flowvar_mod.o \
                      $(OBJDIR)/metrics_mod.o \
                      $(OBJDIR)/element_type_mod.o \
                      $(FUNDIR)/debug_timer.f90

$(OBJDIR)/bnd_profiles_mod.o : bnd_profiles_mod.f90 \
                     $(OBJDIR)/kind_types_mod.o \
                     $(OBJDIR)/eqn_idx_mod.o \
                     $(OBJDIR)/order_mod.o \
                     $(OBJDIR)/geomvar_mod.o \
                     $(OBJDIR)/ovar_mod.o \
                     $(OBJDIR)/flowvar_mod.o \
                     $(OBJDIR)/element_type_mod.o \
                     $(OBJDIR)/interpolation_mod.o \
                     $(OBJDIR)/gpu_mod.o \
                     $(OBJDIR)/thermo_mod.o \
                     $(OBJDIR)/transport_mod.o \
                     $(FUNDIR)/bc_string_to_integer.f90

$(OBJDIR)/cell_ave_mod.o : cell_ave_mod.f90 \
                 $(OBJDIR)/kind_types_mod.o \
                 $(OBJDIR)/eqn_idx_mod.o \
                 $(OBJDIR)/flowvar_mod.o \
                 $(OBJDIR)/geomvar_mod.o \
                 $(OBJDIR)/ovar_mod.o \
                 $(OBJDIR)/thermo_mod.o \
                 $(OBJDIR)/vandermonde_mod.o

$(OBJDIR)/cgns_mod.o : cgns_mod.f90 \
             $(OBJDIR)/kind_types_mod.o \
             $(OBJDIR)/cgnstypes_mod.o \
             $(OBJDIR)/geomvar_mod.o \
             $(OBJDIR)/ovar_mod.o \
             $(OBJDIR)/order_mod.o \
             $(OBJDIR)/connectivity_mod.o \
             $(FUNDIR)/bc_string_to_integer.f90 \
             $(FUNDIR)/smc000_interior_wall_radius.f90 \
             $(FUNDIR)/last.f90 \
             $(FUNDIR)/debug_timer.f90
$(OBJDIR)/cgns_mod.o : OTHERFLAGS += $(CGNSMOD)

$(OBJDIR)/cgnstypes_mod.o : cgnstypes_mod.f90 \
                  $(OBJDIR)/kind_types_mod.o
$(OBJDIR)/cgnstypes_mod.o : OTHERFLAGS += $(CGNSMOD)
#$(OBJDIR)/cgnstypes_mod.o : OTHERFPP += $(CGNSFPP)

$(OBJDIR)/channel_mod.o : channel_mod.f90 \
                $(OBJDIR)/kind_types_mod.o \
                $(OBJDIR)/eqn_idx_mod.o \
                $(OBJDIR)/geomvar_mod.o \
                $(OBJDIR)/ovar_mod.o \
                $(OBJDIR)/flowvar_mod.o \
                $(OBJDIR)/metrics_mod.o

$(OBJDIR)/connectivity_mod.o : connectivity_mod.f90 \
                 $(OBJDIR)/kind_types_mod.o \
                 $(OBJDIR)/geomvar_mod.o \
                 $(OBJDIR)/ovar_mod.o \
                 $(OBJDIR)/order_mod.o \
                 $(FUNDIR)/last.f90 \
                           mpi_defs.h

$(OBJDIR)/correction_mod.o : correction_mod.f90 \
                   $(OBJDIR)/kind_types_mod.o \
                   $(OBJDIR)/matrix_types_mod.o \
                   $(OBJDIR)/order_mod.o \
                   $(OBJDIR)/ovar_mod.o \
                   $(OBJDIR)/element_type_mod.o \
                   $(OBJDIR)/quadrature_mod.o \
                   $(OBJDIR)/vandermonde_mod.o \
                   $(OBJDIR)/poly_mod.o \
                   $(FUNDIR)/solpts2order.f90

$(OBJDIR)/cpu_info_mod.o : cpu_info_mod.f90 \
                 $(OBJDIR)/kind_types_mod.o \
                 $(OBJDIR)/ovar_mod.o

$(OBJDIR)/derivatives_mod.o : derivatives_mod.f90 \
                    $(OBJDIR)/kind_types_mod.o \
                    $(OBJDIR)/matrix_types_mod.o \
                    $(OBJDIR)/order_mod.o \
                    $(OBJDIR)/ovar_mod.o \
                    $(OBJDIR)/element_type_mod.o \
                    $(OBJDIR)/poly_mod.o \
                    $(OBJDIR)/vandermonde_mod.o \
                    $(FUNDIR)/solpts2order.f90 \
                    $(FUNDIR)/check_for_bad_number.f90 \
                    $(FUNDIR)/debug_timer.f90

$(OBJDIR)/dt_mod.o : dt_mod.f90 \
           $(OBJDIR)/kind_types_mod.o \
           $(OBJDIR)/eqn_idx_mod.o \
           $(OBJDIR)/geomvar_mod.o \
           $(OBJDIR)/ovar_mod.o \
           $(OBJDIR)/flowvar_mod.o \
           $(OBJDIR)/metrics_mod.o \
           $(OBJDIR)/mrf_mod.o \
           $(OBJDIR)/timer_mod.o \
           $(FUNDIR)/debug_timer.f90

$(OBJDIR)/element_type_mod.o : element_type_mod.f90 \
                     $(OBJDIR)/kind_types_mod.o

$(OBJDIR)/eperk_mod.o : eperk_mod.f90 \
              $(OBJDIR)/kind_types_mod.o \
              $(OBJDIR)/eqn_idx_mod.o \
              $(OBJDIR)/geomvar_mod.o \
              $(OBJDIR)/ovar_mod.o \
              $(OBJDIR)/flowvar_mod.o \
              $(OBJDIR)/dt_mod.o \
              $(OBJDIR)/timer_mod.o \
              $(OBJDIR)/rhs_mod.o

$(OBJDIR)/eqn_idx_mod.o : eqn_idx_mod.f90 \
                $(OBJDIR)/kind_types_mod.o \
                $(OBJDIR)/ovar_mod.o

$(OBJDIR)/erk_mod.o : erk_mod.f90 \
            $(OBJDIR)/kind_types_mod.o \
            $(OBJDIR)/eqn_idx_mod.o \
            $(OBJDIR)/geomvar_mod.o \
            $(OBJDIR)/ovar_mod.o \
            $(OBJDIR)/flowvar_mod.o \
            $(OBJDIR)/rhs_mod.o \
            $(OBJDIR)/timer_mod.o \
            $(FUNDIR)/debug_timer.f90

$(OBJDIR)/filter_mod.o : filter_mod.f90 \
               $(OBJDIR)/kind_types_mod.o \
               $(OBJDIR)/eqn_idx_mod.o \
               $(OBJDIR)/matrix_types_mod.o \
               $(OBJDIR)/ovar_mod.o \
               $(OBJDIR)/order_mod.o \
               $(OBJDIR)/geomvar_mod.o \
               $(OBJDIR)/flowvar_mod.o \
               $(OBJDIR)/vandermonde_mod.o \
               $(OBJDIR)/thermo_mod.o \
               $(OBJDIR)/transport_mod.o \
               $(FUNDIR)/solpts2order.f90 \
               $(FUNDIR)/debug_timer.f90
$(OBJDIR)/filter_mod.o : OTHERFPP += $(MKLFPP)
$(OBJDIR)/filter_mod.o : OTHERFLAGS += $(MKLFLAG)

$(OBJDIR)/flowvar_mod.o : flowvar_mod.f90 \
                $(OBJDIR)/kind_types_mod.o \
                $(OBJDIR)/eqn_idx_mod.o \
                $(OBJDIR)/geomvar_mod.o \
                $(OBJDIR)/ovar_mod.o \
                $(OBJDIR)/order_mod.o

$(OBJDIR)/flux_mod.o : flux_mod.f90 \
         $(OBJDIR)/kind_types_mod.o \
         $(OBJDIR)/eqn_idx_mod.o \
         $(OBJDIR)/geomvar_mod.o \
         $(OBJDIR)/order_mod.o \
         $(OBJDIR)/ovar_mod.o \
         $(OBJDIR)/flowvar_mod.o \
         $(OBJDIR)/metrics_mod.o \
         $(OBJDIR)/mrf_mod.o \
         $(OBJDIR)/interpolation_mod.o \
         $(OBJDIR)/correction_mod.o \
         $(OBJDIR)/derivatives_mod.o \
         $(OBJDIR)/invs_flux_mod.o \
         $(OBJDIR)/visc_flux_mod.o \
         $(OBJDIR)/upwind_flux_mod.o \
         $(OBJDIR)/projection_mod.o \
         $(OBJDIR)/matrix_types_mod.o \
         $(OBJDIR)/thermo_mod.o \
         $(OBJDIR)/transport_mod.o \
         $(OBJDIR)/thermofix_mod.o \
         $(OBJDIR)/timer_mod.o \
         $(FUNDIR)/debug_timer.f90 \
         $(FUNDIR)/check_for_bad_number.f90

$(OBJDIR)/generic_mod.o : generic_mod.f90 \
            $(OBJDIR)/kind_types_mod.o \
            $(OBJDIR)/pbs_mod.o \
            $(OBJDIR)/eqn_idx_mod.o \
            $(OBJDIR)/geomvar_mod.o \
            $(OBJDIR)/ovar_mod.o \
            $(OBJDIR)/order_mod.o \
            $(OBJDIR)/metrics_mod.o \
            $(OBJDIR)/flowvar_mod.o \
            $(OBJDIR)/element_type_mod.o \
            $(OBJDIR)/gpu_mod.o \
            $(OBJDIR)/thx_mod.o \
            $(OBJDIR)/bc_flow_stats_mod.o \
            $(OBJDIR)/timer_mod.o \
            $(OBJDIR)/reactions_mod.o \
            $(FUNDIR)/debug_timer.f90

$(OBJDIR)/geomvar_mod.o : geomvar_mod.f90 \
               $(OBJDIR)/kind_types_mod.o \
               $(OBJDIR)/ovar_mod.o

$(OBJDIR)/gmsh_mod.o : gmsh_mod.f90 \
             $(OBJDIR)/kind_types_mod.o \
             $(OBJDIR)/geomvar_mod.o \
             $(OBJDIR)/ovar_mod.o \
             $(OBJDIR)/order_mod.o \
             $(FUNDIR)/bc_string_to_integer.f90 \
             $(FUNDIR)/last.f90

$(OBJDIR)/gpu_mod.o : gpu_mod.f90 \
            $(OBJDIR)/kind_types_mod.o \
            $(OBJDIR)/eqn_idx_mod.o \
            $(OBJDIR)/ovar_mod.o \
            $(OBJDIR)/thermo_mod.o

$(OBJDIR)/ignition_mod.o : ignition_mod.f90 \
                 $(OBJDIR)/kind_types_mod.o \
                 $(OBJDIR)/eqn_idx_mod.o \
                 $(OBJDIR)/geomvar_mod.o \
                 $(OBJDIR)/ovar_mod.o \
                 $(OBJDIR)/flowvar_mod.o \
                 $(OBJDIR)/metrics_mod.o

$(OBJDIR)/initialization_mod.o : initialization_mod.f90 \
               $(OBJDIR)/kind_types_mod.o \
               $(OBJDIR)/matrix_types_mod.o \
               $(OBJDIR)/eqn_idx_mod.o \
               $(OBJDIR)/geomvar_mod.o \
               $(OBJDIR)/ovar_mod.o \
               $(OBJDIR)/order_mod.o \
               $(OBJDIR)/flowvar_mod.o \
               $(OBJDIR)/channel_mod.o \
               $(OBJDIR)/metrics_mod.o \
               $(OBJDIR)/element_type_mod.o \
               $(OBJDIR)/restart_mod.o \
               $(OBJDIR)/projection_mod.o \
               $(OBJDIR)/interpolation_mod.o \
               $(OBJDIR)/derivatives_mod.o \
               $(OBJDIR)/thermo_mod.o \
               $(OBJDIR)/transport_mod.o \
               $(OBJDIR)/reactions_mod.o \
               $(OBJDIR)/synthetic_eddy_mod.o \
               $(OBJDIR)/periodic_mod.o \
               $(OBJDIR)/gpu_mod.o \
               $(OBJDIR)/shu_osher_mod.o \
               $(OBJDIR)/advection_js_mod.o \
               $(OBJDIR)/sod_shock_tube_mod.o \
               $(OBJDIR)/taylor_couette_mod.o \
               $(FUNDIR)/taylor_green.f90 \
               $(FUNDIR)/smc000_interior_wall_radius.f90 \
               $(FUNDIR)/arn2_interior_wall_radius.f90 \
               $(FUNDIR)/ssjet_interior_wall_radius.f90 \
               $(FUNDIR)/sod.f90 \
               $(FUNDIR)/pulse.f90 \
               $(FUNDIR)/uiuc_combustor.f90 \
               $(FUNDIR)/check_for_bad_number.f90

$(OBJDIR)/internal_grid_mod.o : internal_grid_mod.f90 \
                      $(OBJDIR)/kind_types_mod.o \
                      $(OBJDIR)/order_mod.o \
                      $(OBJDIR)/geomvar_mod.o \
                      $(OBJDIR)/ovar_mod.o \
                      $(OBJDIR)/cgns_mod.o \
                      $(OBJDIR)/poly_mod.o \
                      $(OBJDIR)/advection_js_mod.o \
                      $(OBJDIR)/sod_shock_tube_mod.o
$(OBJDIR)/internal_grid_mod.o : OTHERFLAGS += $(CGNSMOD)

$(OBJDIR)/interpolation_mod.o : interpolation_mod.f90 \
                      $(OBJDIR)/kind_types_mod.o \
                      $(OBJDIR)/eqn_idx_mod.o \
                      $(OBJDIR)/matrix_types_mod.o \
                      $(OBJDIR)/geomvar_mod.o \
                      $(OBJDIR)/ovar_mod.o \
                      $(OBJDIR)/order_mod.o \
                      $(OBJDIR)/element_type_mod.o \
                      $(OBJDIR)/quadrature_mod.o \
                      $(OBJDIR)/vandermonde_mod.o \
                      $(OBJDIR)/poly_mod.o \
                      $(OBJDIR)/flowvar_mod.o \
                      $(FUNDIR)/solpts2order.f90 \
                      $(FUNDIR)/cell_solpts.f90

$(OBJDIR)/invs_flux_mod.o : invs_flux_mod.f90 \
                  $(OBJDIR)/kind_types_mod.o \
                  $(OBJDIR)/eqn_idx_mod.o \
                  $(OBJDIR)/geomvar_mod.o \
                  $(OBJDIR)/ovar_mod.o \
                  $(OBJDIR)/timer_mod.o \
                  $(OBJDIR)/mrf_mod.o

$(OBJDIR)/io_mod.o : io_mod.f90 \
       $(OBJDIR)/kind_types_mod.o \
       $(OBJDIR)/eqn_idx_mod.o \
       $(OBJDIR)/order_mod.o \
       $(OBJDIR)/geomvar_mod.o \
       $(OBJDIR)/ovar_mod.o \
       $(OBJDIR)/nml_base_mod.o \
       $(OBJDIR)/nml_0001_mod.o \
       $(OBJDIR)/flowvar_mod.o \
       $(OBJDIR)/gmsh_mod.o \
       $(OBJDIR)/plot3d_mod.o \
       $(OBJDIR)/cgns_mod.o \
       $(OBJDIR)/internal_grid_mod.o \
       $(OBJDIR)/parallel_mod.o \
       $(OBJDIR)/quadrature_mod.o \
       $(OBJDIR)/interpolation_mod.o \
       $(OBJDIR)/initialization_mod.o \
       $(OBJDIR)/oerk_mod.o \
       $(OBJDIR)/eperk_mod.o \
       $(OBJDIR)/restart_mod.o \
       $(OBJDIR)/thermo_mod.o \
       $(OBJDIR)/synthetic_eddy_mod.o \
       $(OBJDIR)/periodic_mod.o \
       $(FUNDIR)/last.f90 \
       $(FUNDIR)/bc_priority.f90 \
       $(FUNDIR)/cell_vol.f90 \
       $(FUNDIR)/bc_string_to_integer.f90 \
       $(FUNDIR)/check_for_bad_number.f90

$(OBJDIR)/kind_types_mod.o : kind_types_mod.f90 \
               $(FUNDIR)/all_are_equal.f90 \
               $(FUNDIR)/bc_integer_to_string.f90 \
               $(FUNDIR)/checks_for_bad_numbers.f90 \
               $(FUNDIR)/chop.f90 \
               $(FUNDIR)/cross_product.f90 \
               $(FUNDIR)/deltafun.f90 \
               $(FUNDIR)/face_normal.f90 \
               $(FUNDIR)/init_array.f90 \
               $(FUNDIR)/spread.f90 \
               $(FUNDIR)/intseq.f90 \
               $(FUNDIR)/jacobians.f90 \
               $(FUNDIR)/lowercase.f90 \
               $(FUNDIR)/matrix_ops.f90 \
               $(FUNDIR)/norms.f90 \
               $(FUNDIR)/random.f90 \
               $(FUNDIR)/reallocate.f90 \
               $(FUNDIR)/reciprocal.f90 \
               $(FUNDIR)/reformat_string.f90 \
               $(FUNDIR)/reverse_metrics.f90 \
               $(FUNDIR)/uppercase.f90 \
               $(FUNDIR)/debug_timer.f90 \
                         mpi_defs.h
$(OBJDIR)/kind_types_mod.o : OTHERFPP += $(MKLFPP)
$(OBJDIR)/kind_types_mod.o : OTHERFLAGS += $(MKLFLAG)

$(OBJDIR)/limiters_mod.o : limiters_mod.f90 \
             $(OBJDIR)/kind_types_mod.o \
             $(OBJDIR)/eqn_idx_mod.o \
             $(OBJDIR)/geomvar_mod.o \
             $(OBJDIR)/ovar_mod.o \
             $(OBJDIR)/order_mod.o \
             $(OBJDIR)/flowvar_mod.o \
             $(OBJDIR)/metrics_mod.o \
             $(OBJDIR)/vandermonde_mod.o \
             $(OBJDIR)/projection_mod.o \
             $(OBJDIR)/filter_mod.o \
             $(OBJDIR)/element_type_mod.o \
             $(OBJDIR)/thermo_mod.o \
             $(OBJDIR)/cell_ave_mod.o \
             $(OBJDIR)/flux_mod.o \
             $(FUNDIR)/uiuc_combustor.f90 \
             $(FUNDIR)/arn2_interior_wall_radius.f90 \
             $(FUNDIR)/debug_timer.f90

$(OBJDIR)/matrix_types_mod.o : matrix_types_mod.f90 \
                     $(OBJDIR)/kind_types_mod.o \
                     $(OBJDIR)/ovar_mod.o
$(OBJDIR)/matrix_types_mod.o : OTHERFPP += $(MKLFPP)
$(OBJDIR)/matrix_types_mod.o : OTHERFLAGS += $(MKLFLAG)

$(OBJDIR)/mappings_mod.o : mappings_mod.f90 \
             $(OBJDIR)/kind_types_mod.o \
             $(OBJDIR)/eqn_idx_mod.o \
             $(OBJDIR)/geomvar_mod.o \
             $(OBJDIR)/ovar_mod.o \
             $(OBJDIR)/order_mod.o \
             $(OBJDIR)/flowvar_mod.o \
             $(OBJDIR)/element_type_mod.o \
             $(OBJDIR)/matrix_types_mod.o \
             $(OBJDIR)/poly_mod.o \
             $(OBJDIR)/quadrature_mod.o \
             $(OBJDIR)/vandermonde_mod.o \
             $(FUNDIR)/solpts2order.f90

$(OBJDIR)/metis_mod.o : metis_mod.f90
$(OBJDIR)/metis_mod.o : OTHERFPP += $(METISFPP)

$(OBJDIR)/metrics_mod.o : metrics_mod.f90 \
            $(OBJDIR)/kind_types_mod.o \
            $(OBJDIR)/geomvar_mod.o \
            $(OBJDIR)/order_mod.o \
            $(OBJDIR)/ovar_mod.o \
            $(OBJDIR)/element_type_mod.o \
            $(OBJDIR)/interpolation_mod.o \
            $(OBJDIR)/derivatives_mod.o

$(OBJDIR)/mlp_mod.o : mlp_mod.f90 \
            $(OBJDIR)/kind_types_mod.o \
            $(OBJDIR)/cell_ave_mod.o \
            $(OBJDIR)/eqn_idx_mod.o \
            $(OBJDIR)/filter_mod.o \
            $(OBJDIR)/flowvar_mod.o \
            $(OBJDIR)/geomvar_mod.o \
            $(OBJDIR)/interpolation_mod.o \
            $(OBJDIR)/order_mod.o \
            $(OBJDIR)/ovar_mod.o \
            $(OBJDIR)/parallel_mod.o

$(OBJDIR)/mrf_mod.o : mrf_mod.f90 \
            $(OBJDIR)/kind_types_mod.o \
            $(OBJDIR)/eqn_idx_mod.o \
            $(OBJDIR)/geomvar_mod.o \
            $(OBJDIR)/ovar_mod.o \
            $(OBJDIR)/flowvar_mod.o \
            $(OBJDIR)/timer_mod.o \
            $(OBJDIR)/metrics_mod.o

$(OBJDIR)/mycpu.o : $(FUNDIR)/mycpu.c

$(OBJDIR)/oerk_mod.o : oerk_mod.f90 \
             $(OBJDIR)/kind_types_mod.o \
             $(OBJDIR)/eqn_idx_mod.o \
             $(OBJDIR)/geomvar_mod.o \
             $(OBJDIR)/ovar_mod.o \
             $(OBJDIR)/flowvar_mod.o \
             $(OBJDIR)/timer_mod.o \
             $(OBJDIR)/rhs_mod.o

$(OBJDIR)/openmp_mod.o : openmp_mod.f90

$(OBJDIR)/order_mod.o : order_mod.f90 \
              $(OBJDIR)/kind_types_mod.o \
              $(OBJDIR)/geomvar_mod.o \
              $(FUNDIR)/cell_solpts.f90

$(OBJDIR)/ovar_mod.o : ovar_mod.f90 \
             $(OBJDIR)/kind_types_mod.o

$(OBJDIR)/nml_0001_mod.o : nml_0001_mod.f90 \
             $(OBJDIR)/kind_types_mod.o \
             $(OBJDIR)/ovar_mod.o

$(OBJDIR)/parallel_mod.o : parallel_mod.f90 \
             $(OBJDIR)/kind_types_mod.o \
             $(OBJDIR)/eqn_idx_mod.o \
             $(OBJDIR)/metis_mod.o \
             $(OBJDIR)/parmetis_mod.o \
             $(OBJDIR)/geomvar_mod.o \
             $(OBJDIR)/ovar_mod.o \
             $(OBJDIR)/order_mod.o \
             $(OBJDIR)/flowvar_mod.o \
             $(OBJDIR)/connectivity_mod.o \
             $(OBJDIR)/gmsh_mod.o \
             $(OBJDIR)/cgns_mod.o \
             $(OBJDIR)/plot3d_mod.o \
             $(OBJDIR)/timer_mod.o \
             $(FUNDIR)/cell_solpts.f90 \
             $(FUNDIR)/last.f90 \
             $(FUNDIR)/debug_timer.f90 \
                       mpi_defs.h
$(OBJDIR)/parallel_mod.o : OTHERFPP += $(PARMETISFPP)

$(OBJDIR)/parmetis_mod.o : parmetis_mod.f90 \
                 $(OBJDIR)/kind_types_mod.o \
                 $(OBJDIR)/ovar_mod.o \
                 $(OBJDIR)/metis_mod.o
$(OBJDIR)/parmetis_mod.o : OTHERFPP += $(PARMETISFPP)

$(OBJDIR)/pbs_mod.o : pbs_mod.f90 \
            $(OBJDIR)/kind_types_mod.o

$(OBJDIR)/pbs_time_left.o : $(FUNDIR)/pbs_time_left.c

$(OBJDIR)/periodic_mod.o : periodic_mod.f90 \
                 $(OBJDIR)/kind_types_mod.o \
                 $(OBJDIR)/geomvar_mod.o \
                 $(OBJDIR)/ovar_mod.o \
                 $(OBJDIR)/cgnstypes_mod.o \
                 $(OBJDIR)/connectivity_mod.o
$(OBJDIR)/periodic_mod.o : OTHERFLAGS += $(CGNSMOD)

$(OBJDIR)/plot3d_mod.o : plot3d_mod.f90 \
               $(OBJDIR)/kind_types_mod.o \
               $(OBJDIR)/geomvar_mod.o \
               $(OBJDIR)/ovar_mod.o \
               $(OBJDIR)/order_mod.o \
               $(OBJDIR)/plot3dio_mod.o \
               $(FUNDIR)/bc_string_to_integer.f90 \
               $(FUNDIR)/last.f90
#              $(FUNDIR)/check_for_bad_number.f90 \

$(OBJDIR)/plot3dio_mod.o : plot3dio_mod.f90 \
                 $(OBJDIR)/kind_types_mod.o

$(OBJDIR)/poly_mod.o : poly_mod.f90 \
         $(OBJDIR)/kind_types_mod.o \
         $(FUNDIR)/check_for_bad_number.f90

$(OBJDIR)/postproc_mod.o : postproc_mod.f90 \
                 $(OBJDIR)/kind_types_mod.o \
                 $(OBJDIR)/cgnstypes_mod.o \
                 $(OBJDIR)/eqn_idx_mod.o \
                 $(OBJDIR)/geomvar_mod.o \
                 $(OBJDIR)/ovar_mod.o \
                 $(OBJDIR)/flowvar_mod.o \
                 $(OBJDIR)/bc_mod.o \
                 $(OBJDIR)/order_mod.o \
                 $(OBJDIR)/interpolation_mod.o \
                 $(OBJDIR)/parallel_mod.o \
                 $(OBJDIR)/connectivity_mod.o \
                 $(OBJDIR)/thermofix_mod.o \
                 $(OBJDIR)/thermo_mod.o \
                 $(OBJDIR)/transport_mod.o \
                 $(OBJDIR)/gpu_mod.o \
                 $(OBJDIR)/metrics_mod.o \
                 $(OBJDIR)/limiters_mod.o \
                 $(OBJDIR)/restart_mod.o \
                 $(FUNDIR)/cgns_write_field_data.f90 \
                 $(FUNDIR)/cgns_write_field_header.f90 \
                 $(FUNDIR)/last.f90 \
                 $(FUNDIR)/debug_timer.f90 \
                 $(FUNDIR)/cell_solpts.f90
$(OBJDIR)/postproc_mod.o : OTHERFLAGS += $(CGNSMOD)

$(OBJDIR)/probes_mod.o : probes_mod.f90 \
               $(OBJDIR)/kind_types_mod.o \
               $(OBJDIR)/geomvar_mod.o \
               $(OBJDIR)/ovar_mod.o \
               $(OBJDIR)/io_mod.o

$(OBJDIR)/prism_mod.o : prism_mod.f90 \
          $(OBJDIR)/kind_types_mod.o \
          $(OBJDIR)/poly_mod.o \
          $(OBJDIR)/triangle_mod.o

$(OBJDIR)/projection_mod.o : projection_mod.f90 \
                   $(OBJDIR)/kind_types_mod.o \
                   $(OBJDIR)/eqn_idx_mod.o \
                   $(OBJDIR)/matrix_types_mod.o \
                   $(OBJDIR)/geomvar_mod.o \
                   $(OBJDIR)/order_mod.o \
                   $(OBJDIR)/ovar_mod.o \
                   $(OBJDIR)/element_type_mod.o \
                   $(OBJDIR)/quadrature_mod.o \
                   $(OBJDIR)/vandermonde_mod.o \
                   $(OBJDIR)/thermo_mod.o \
                   $(OBJDIR)/transport_mod.o \
                   $(FUNDIR)/debug_timer.f90
$(OBJDIR)/projection_mod.o : OTHERFPP += $(MKLFPP)
$(OBJDIR)/projection_mod.o : OTHERFLAGS += $(MKLFLAG)

$(OBJDIR)/pyramid_mod.o : pyramid_mod.f90 \
            $(OBJDIR)/kind_types_mod.o \
            $(OBJDIR)/poly_mod.o

$(OBJDIR)/quadrature_mod.o : quadrature_mod.f90 \
                   $(OBJDIR)/kind_types_mod.o \
                   $(OBJDIR)/geomvar_mod.o \
                   $(OBJDIR)/order_mod.o \
                   $(OBJDIR)/ovar_mod.o \
                   $(OBJDIR)/element_type_mod.o \
                   $(OBJDIR)/poly_mod.o \
                   $(OBJDIR)/vandermonde_mod.o \
                   $(OBJDIR)/triangle_mod.o \
                   $(OBJDIR)/tetrahedra_mod.o \
                   $(OBJDIR)/pyramid_mod.o \
                   $(OBJDIR)/prism_mod.o \
                   $(FUNDIR)/solpts2order.f90 \
                   $(FUNDIR)/check_for_bad_number.f90 \
                   $(FUNDIR)/debug_timer.f90

$(OBJDIR)/random_mod.o : random_mod.f90 \
               $(OBJDIR)/kind_types_mod.o \
               $(OBJDIR)/ovar_mod.o

$(OBJDIR)/reactions_mod.o : reactions_mod.f90 \
                  $(OBJDIR)/kind_types_mod.o \
                  $(OBJDIR)/eqn_idx_mod.o \
                  $(OBJDIR)/ovar_mod.o \
                  $(OBJDIR)/flowvar_mod.o \
                  $(OBJDIR)/metrics_mod.o \
                  $(OBJDIR)/thermo_mod.o \
                  $(OBJDIR)/ignition_mod.o \
                  $(FUNDIR)/check_for_bad_number.f90

$(OBJDIR)/restart_mod.o : restart_mod.f90 \
                $(OBJDIR)/kind_types_mod.o \
                $(OBJDIR)/eqn_idx_mod.o \
                $(OBJDIR)/geomvar_mod.o \
                $(OBJDIR)/ovar_mod.o \
                $(OBJDIR)/flowvar_mod.o \
                $(OBJDIR)/metrics_mod.o \
                $(OBJDIR)/random_mod.o \
                $(OBJDIR)/thermo_mod.o \
                $(OBJDIR)/synthetic_eddy_mod.o \
                $(OBJDIR)/order_mod.o \
                $(OBJDIR)/interpolation_mod.o \
                $(OBJDIR)/projection_mod.o \
                $(OBJDIR)/matrix_types_mod.o \
                $(FUNDIR)/cell_solpts.f90 \
                $(FUNDIR)/sort.f90 \
                $(FUNDIR)/debug_timer.f90 \
                          mpi_defs.h

$(OBJDIR)/rhs_mod.o : rhs_mod.f90 \
            $(OBJDIR)/kind_types_mod.o \
            $(OBJDIR)/eqn_idx_mod.o \
            $(OBJDIR)/geomvar_mod.o \
            $(OBJDIR)/ovar_mod.o \
            $(OBJDIR)/order_mod.o \
            $(OBJDIR)/flowvar_mod.o \
            $(OBJDIR)/bc_mod.o \
            $(OBJDIR)/flux_mod.o \
            $(OBJDIR)/metrics_mod.o \
            $(OBJDIR)/filter_mod.o \
            $(OBJDIR)/limiters_mod.o \
            $(OBJDIR)/mlp_mod.o \
            $(OBJDIR)/parallel_mod.o \
            $(OBJDIR)/projection_mod.o \
            $(OBJDIR)/thermofix_mod.o \
            $(OBJDIR)/thermo_mod.o \
            $(OBJDIR)/transport_mod.o \
            $(OBJDIR)/timer_mod.o \
            $(FUNDIR)/debug_timer.f90

$(OBJDIR)/shu_osher_mod.o : shu_osher_mod.f90 \
                $(OBJDIR)/kind_types_mod.o \
                $(OBJDIR)/ovar_mod.o \
                $(OBJDIR)/thermo_mod.o

$(OBJDIR)/spatial_ave_mod.o : spatial_ave_mod.f90 \
                    $(OBJDIR)/kind_types_mod.o \
                    $(OBJDIR)/cgnstypes_mod.o \
                    $(OBJDIR)/order_mod.o \
                    $(OBJDIR)/geomvar_mod.o \
                    $(OBJDIR)/ovar_mod.o \
                    $(OBJDIR)/flowvar_mod.o \
                    $(OBJDIR)/element_type_mod.o \
                    $(OBJDIR)/metrics_mod.o \
                    $(OBJDIR)/derivatives_mod.o \
                    $(OBJDIR)/transport_mod.o \
                    $(OBJDIR)/synthetic_eddy_mod.o \
                    $(OBJDIR)/io_mod.o \
                    $(OBJDIR)/write_cgns_mod.o
$(OBJDIR)/spatial_ave_mod.o : OTHERFLAGS += $(CGNSMOD)

$(OBJDIR)/taylor_couette_mod.o : taylor_couette_mod.f90 \
                       $(OBJDIR)/kind_types_mod.o \
                       $(OBJDIR)/eqn_idx_mod.o \
                       $(OBJDIR)/geomvar_mod.o \
                       $(OBJDIR)/ovar_mod.o \
                       $(OBJDIR)/thermo_mod.o

$(OBJDIR)/tetrahedra_mod.o : tetrahedra_mod.f90 \
               $(OBJDIR)/kind_types_mod.o \
               $(OBJDIR)/poly_mod.o

$(OBJDIR)/thermofix_mod.o : thermofix_mod.f90 \
               $(OBJDIR)/kind_types_mod.o \
               $(OBJDIR)/cell_ave_mod.o \
               $(OBJDIR)/eqn_idx_mod.o \
               $(OBJDIR)/filter_mod.o \
               $(OBJDIR)/flowvar_mod.o \
               $(OBJDIR)/geomvar_mod.o \
               $(OBJDIR)/interpolation_mod.o \
               $(OBJDIR)/mlp_mod.o \
               $(OBJDIR)/order_mod.o \
               $(OBJDIR)/ovar_mod.o \
               $(OBJDIR)/thermo_mod.o \
               $(OBJDIR)/transport_mod.o \
               $(FUNDIR)/debug_timer.f90

$(OBJDIR)/thermo_mod.o : thermo_mod.f90 \
               $(OBJDIR)/kind_types_mod.o \
               $(OBJDIR)/eqn_idx_mod.o \
               $(OBJDIR)/ovar_mod.o \
               $(OBJDIR)/flowvar_mod.o

$(OBJDIR)/timer_mod.o : timer_mod.f90 \
              $(OBJDIR)/kind_types_mod.o \
              $(OBJDIR)/ovar_mod.o

$(OBJDIR)/transport_mod.o : transport_mod.f90 \
                  $(OBJDIR)/kind_types_mod.o \
                  $(OBJDIR)/eqn_idx_mod.o \
                  $(OBJDIR)/ovar_mod.o \
                  $(OBJDIR)/thermo_mod.o

$(OBJDIR)/triangle_mod.o : triangle_mod.f90 \
             $(OBJDIR)/kind_types_mod.o \
             $(OBJDIR)/poly_mod.o

$(OBJDIR)/synthetic_eddy_mod.o : synthetic_eddy_mod.f90 \
                       $(OBJDIR)/kind_types_mod.o \
                       $(OBJDIR)/random_mod.o \
                       $(OBJDIR)/eqn_idx_mod.o \
                       $(OBJDIR)/geomvar_mod.o \
                       $(OBJDIR)/ovar_mod.o \
                       $(OBJDIR)/flowvar_mod.o \
                       $(OBJDIR)/metrics_mod.o \
                       $(OBJDIR)/thermo_mod.o \
                       $(OBJDIR)/transport_mod.o \
                       $(FUNDIR)/check_for_bad_number.f90

$(OBJDIR)/tecio_mod.o : tecio_mod.f90

$(OBJDIR)/thx_mod.o : thx_mod.f90 \
            $(OBJDIR)/kind_types_mod.o \
            $(OBJDIR)/eqn_idx_mod.o \
            $(OBJDIR)/geomvar_mod.o \
            $(OBJDIR)/ovar_mod.o \
            $(OBJDIR)/flowvar_mod.o \
            $(OBJDIR)/metrics_mod.o \
            $(OBJDIR)/element_type_mod.o \
            $(FUNDIR)/debug_timer.f90

$(OBJDIR)/upwind_flux_mod.o : upwind_flux_mod.f90 \
         $(OBJDIR)/kind_types_mod.o \
         $(OBJDIR)/eqn_idx_mod.o \
         $(OBJDIR)/geomvar_mod.o \
         $(OBJDIR)/ovar_mod.o \
         $(OBJDIR)/timer_mod.o \
         $(OBJDIR)/mrf_mod.o

$(OBJDIR)/vandermonde_mod.o : vandermonde_mod.f90 \
                    $(OBJDIR)/kind_types_mod.o \
                    $(OBJDIR)/matrix_types_mod.o \
                    $(OBJDIR)/order_mod.o \
                    $(OBJDIR)/element_type_mod.o \
                    $(OBJDIR)/poly_mod.o \
                    $(FUNDIR)/solpts2order.f90 \
                    $(FUNDIR)/cell_solpts.f90 \
                    $(FUNDIR)/debug_timer.f90

$(OBJDIR)/visc_flux_mod.o : visc_flux_mod.f90 \
         $(OBJDIR)/kind_types_mod.o \
         $(OBJDIR)/eqn_idx_mod.o \
         $(OBJDIR)/geomvar_mod.o \
         $(OBJDIR)/ovar_mod.o \
         $(OBJDIR)/gpu_mod.o \
         $(OBJDIR)/thermo_mod.o \
         $(OBJDIR)/timer_mod.o \
         $(OBJDIR)/transport_mod.o

$(OBJDIR)/write_cgns_mod.o : write_cgns_mod.f90 \
                   $(OBJDIR)/kind_types_mod.o \
                   $(OBJDIR)/cgnstypes_mod.o \
                   $(OBJDIR)/eqn_idx_mod.o \
                   $(OBJDIR)/order_mod.o \
                   $(OBJDIR)/geomvar_mod.o \
                   $(OBJDIR)/ovar_mod.o \
                   $(OBJDIR)/metrics_mod.o \
                   $(OBJDIR)/flowvar_mod.o \
                   $(OBJDIR)/interpolation_mod.o \
                   $(OBJDIR)/derivatives_mod.o \
                   $(OBJDIR)/projection_mod.o \
                   $(OBJDIR)/limiters_mod.o \
                   $(OBJDIR)/mlp_mod.o \
                   $(OBJDIR)/advection_js_mod.o \
                   $(OBJDIR)/sod_shock_tube_mod.o \
                   $(OBJDIR)/thermofix_mod.o \
                   $(OBJDIR)/thermo_mod.o \
                   $(OBJDIR)/transport_mod.o \
                   $(OBJDIR)/flux_mod.o \
                   $(OBJDIR)/gpu_mod.o \
                   $(OBJDIR)/restart_mod.o \
                   $(OBJDIR)/taylor_couette_mod.o \
                   $(OBJDIR)/mrf_mod.o \
                   $(FUNDIR)/cgns_write_field_data.f90 \
                   $(FUNDIR)/cgns_write_field_header.f90 \
                             mpi_defs.h
$(OBJDIR)/write_cgns_mod.o : OTHERFLAGS += $(CGNSMOD)

# Source files in BC_Modules
$(OBJDIR)/bc_wall_mod.o : $(BCDIR)/bc_wall_mod.f90 \
                         $(OBJDIR)/kind_types_mod.o \
                         $(OBJDIR)/eqn_idx_mod.o \
                         $(OBJDIR)/ovar_mod.o \
                         $(OBJDIR)/gpu_mod.o \
                         $(OBJDIR)/thermo_mod.o \
                         $(OBJDIR)/bc_relax_mod.o

$(OBJDIR)/bc_fixed_mod.o : $(BCDIR)/bc_fixed_mod.f90 \
                 $(OBJDIR)/kind_types_mod.o \
                 $(OBJDIR)/eqn_idx_mod.o

$(OBJDIR)/bc_flowthrough_grad_mod.o : $(BCDIR)/bc_flowthrough_grad_mod.f90 \
                            $(OBJDIR)/kind_types_mod.o \
                            $(OBJDIR)/eqn_idx_mod.o

$(OBJDIR)/bc_gen_free_mod.o : $(BCDIR)/bc_gen_free_mod.f90 \
                    $(OBJDIR)/kind_types_mod.o \
                    $(OBJDIR)/eqn_idx_mod.o \
                    $(OBJDIR)/ovar_mod.o \
                    $(OBJDIR)/thermo_mod.o \
                    $(OBJDIR)/bc_relax_mod.o

$(OBJDIR)/bc_gen_in_mod.o : $(BCDIR)/bc_gen_in_mod.f90 \
                  $(OBJDIR)/kind_types_mod.o \
                  $(OBJDIR)/eqn_idx_mod.o \
                  $(OBJDIR)/ovar_mod.o \
                  $(OBJDIR)/thermo_mod.o \
                  $(OBJDIR)/bc_relax_mod.o

$(OBJDIR)/bc_gen_out_mod.o : $(BCDIR)/bc_gen_out_mod.f90 \
                   $(OBJDIR)/kind_types_mod.o \
                   $(OBJDIR)/eqn_idx_mod.o \
                   $(OBJDIR)/ovar_mod.o \
                   $(OBJDIR)/bc_relax_mod.o

$(OBJDIR)/bc_periodic_mod.o : $(BCDIR)/bc_periodic_mod.f90 \
                    $(OBJDIR)/kind_types_mod.o \
                    $(OBJDIR)/eqn_idx_mod.o

$(OBJDIR)/bc_relax_mod.o : $(BCDIR)/bc_relax_mod.f90 \
                 $(OBJDIR)/kind_types_mod.o \
                 $(OBJDIR)/eqn_idx_mod.o \
                 $(OBJDIR)/ovar_mod.o \
                 $(OBJDIR)/thermo_mod.o \
                 $(OBJDIR)/transport_mod.o

$(OBJDIR)/bc_riemann_mod.o : $(BCDIR)/bc_riemann_mod.f90 \
                   $(OBJDIR)/kind_types_mod.o \
                   $(OBJDIR)/eqn_idx_mod.o \
                   $(OBJDIR)/ovar_mod.o \
                   $(OBJDIR)/bc_relax_mod.o

$(OBJDIR)/bc_sub_in_mod.o : $(BCDIR)/bc_sub_in_mod.f90 \
                  $(OBJDIR)/kind_types_mod.o \
                  $(OBJDIR)/eqn_idx_mod.o \
                  $(OBJDIR)/order_mod.o \
                  $(OBJDIR)/geomvar_mod.o \
                  $(OBJDIR)/ovar_mod.o \
                  $(OBJDIR)/interpolation_mod.o \
                  $(OBJDIR)/thermo_mod.o \
                  $(OBJDIR)/bc_relax_mod.o

$(OBJDIR)/bc_sub_out_mod.o : $(BCDIR)/bc_sub_out_mod.f90 \
                   $(OBJDIR)/kind_types_mod.o \
                   $(OBJDIR)/eqn_idx_mod.o \
                   $(OBJDIR)/ovar_mod.o \
                   $(OBJDIR)/bc_relax_mod.o

$(OBJDIR)/bc_sup_out_mod.o : $(BCDIR)/bc_sup_out_mod.f90 \
                   $(OBJDIR)/kind_types_mod.o \
                   $(OBJDIR)/eqn_idx_mod.o

# Source files in Nml_Modules
$(OBJDIR)/nml_base_mod.o : $(NMLDIR)/nml_base_mod.f90 \
                 $(OBJDIR)/kind_types_mod.o

# Source files in SolInterp_Modules
$(OBJDIR)/interpsol_mod.o : $(INTDIR)/interpsol_mod.f90 \
                  $(OBJDIR)/kind_types_mod.o \
                  $(OBJDIR)/eqn_idx_mod.o \
                  $(OBJDIR)/geomvar_mod.o \
                  $(OBJDIR)/ovar_mod.o \
                  $(OBJDIR)/order_mod.o \
                  $(OBJDIR)/flowvar_mod.o \
                  $(OBJDIR)/flux_mod.o \
                  $(OBJDIR)/element_type_mod.o \
                  $(OBJDIR)/vandermonde_mod.o \
                  $(OBJDIR)/derivatives_mod.o \
                  $(OBJDIR)/correction_mod.o \
                  $(OBJDIR)/interpolation_mod.o \
                  $(OBJDIR)/projection_mod.o \
                  $(OBJDIR)/filter_mod.o \
                  $(OBJDIR)/interpsolvar_mod.o \
                  $(OBJDIR)/metrics_mod.o \
                  $(OBJDIR)/octedat_facexyz_mod.o \
                  $(OBJDIR)/octstat_mod.o \
                  $(OBJDIR)/octtree_mod.o \
                  $(OBJDIR)/interpsol_query_mod.o \
                  $(FUNDIR)/debug_timer.f90

$(OBJDIR)/interpsolvar_mod.o : $(INTDIR)/interpsolvar_mod.f90 \
                     $(OBJDIR)/kind_types_mod.o \
                     $(OBJDIR)/geomvar_mod.o \
                     $(OBJDIR)/flowvar_mod.o \
                     $(OBJDIR)/matrix_types_mod.o \
                     $(OBJDIR)/element_type_mod.o \
                     $(OBJDIR)/vandermonde_mod.o \
                     $(OBJDIR)/derivatives_mod.o \
                     $(OBJDIR)/correction_mod.o \
                     $(OBJDIR)/interpolation_mod.o \
                     $(OBJDIR)/projection_mod.o \
                     $(OBJDIR)/filter_mod.o

$(OBJDIR)/interpsol_query_mod.o : $(INTDIR)/interpsol_query_mod.f90 \
                        $(OBJDIR)/kind_types_mod.o \
                        $(OBJDIR)/eqn_idx_mod.o \
                        $(OBJDIR)/geomvar_mod.o \
                        $(OBJDIR)/flowvar_mod.o \
                        $(OBJDIR)/interpsolvar_mod.o

$(OBJDIR)/octedat_facexyz_mod.o : $(INTDIR)/octedat_facexyz_mod.f90 \
                        $(OBJDIR)/kind_types_mod.o \
                        $(OBJDIR)/geomvar_mod.o \
                        $(OBJDIR)/flowvar_mod.o \
                        $(OBJDIR)/octedat_mod.o \
                        $(OBJDIR)/octstat_mod.o

$(OBJDIR)/octedat_mod.o : $(INTDIR)/octedat_mod.f90 \
                $(OBJDIR)/kind_types_mod.o \
                $(OBJDIR)/octstat_mod.o

$(OBJDIR)/octnode_mod.o : $(INTDIR)/octnode_mod.f90 \
                $(OBJDIR)/kind_types_mod.o \
                $(OBJDIR)/geomvar_mod.o \
                $(OBJDIR)/octedat_mod.o \
                $(OBJDIR)/octroot_mod.o \
                $(OBJDIR)/octstat_mod.o

$(OBJDIR)/octroot_mod.o : $(INTDIR)/octroot_mod.f90 \
                $(OBJDIR)/octstat_mod.o

$(OBJDIR)/octstat_mod.o : $(INTDIR)/octstat_mod.f90

$(OBJDIR)/octtree_mod.o : $(INTDIR)/octtree_mod.f90 \
                $(OBJDIR)/kind_types_mod.o \
                $(OBJDIR)/geomvar_mod.o \
                $(OBJDIR)/octedat_mod.o \
                $(OBJDIR)/octnode_mod.o \
                $(OBJDIR)/octroot_mod.o \
                $(OBJDIR)/octstat_mod.o \
                $(FUNDIR)/debug_timer.f90

# Source files in Testcase_Modules
$(OBJDIR)/advection_js_mod.o : $(TSTDIR)/advection_js_mod.f90 \
                     $(OBJDIR)/kind_types_mod.o \
                     $(OBJDIR)/geomvar_mod.o \
                     $(OBJDIR)/eqn_idx_mod.o \
                     $(OBJDIR)/ovar_mod.o \
                     $(OBJDIR)/thermo_mod.o

$(OBJDIR)/sod_shock_tube_mod.o : $(TSTDIR)/sod_shock_tube_mod.f90 \
                       $(OBJDIR)/kind_types_mod.o \
                       $(OBJDIR)/geomvar_mod.o \
                       $(OBJDIR)/eqn_idx_mod.o \
                       $(OBJDIR)/ovar_mod.o \
                       $(OBJDIR)/thermo_mod.o

####################################################################
# ADD ANY OBJECT FILES TO THE VARIABLE CUSTOM_DEBUG_FILES THAT YOU #
# WANT TO SEPARATELY COMPILE WITH DEBUG FLAGS (I.E., FULL BOUNDS   #
# CHECKING, NO OPTIMIZATIONS, FLOATING-POINT EXCEPTIONS, ETC.)     #
# WHEN COMPILING ANY OF THE OPTIMIZED BUILD TARGETS. ALL OTHER     #
# OBJECT FILES WILL BE COMPILED AS USUAL WITH OPTIMIZATIONS.       #
####################################################################
#CUSTOM_DEBUG_FILES := $(OBJDIR)/restart_mod.o \
#                      $(OBJDIR)/averaging_mod.o

$(CUSTOM_DEBUG_FILES) : FCFLAG = $(FCBASE) $(DBGFLAG) $(ADDFLAGS)
