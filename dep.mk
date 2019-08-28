$(OBJDIR)/main.o : main.f90 \
         $(OBJDIR)/kind_types.o \
         $(OBJDIR)/eqn_idx_mod.o \
         $(OBJDIR)/pbs_mod.o \
         $(OBJDIR)/order_mod.o \
         $(OBJDIR)/geovar_mod.o \
         $(OBJDIR)/ovar_mod.o \
         $(OBJDIR)/flowvar_mod.o \
         $(OBJDIR)/quadrature_mod.o \
         $(OBJDIR)/vandermonde_mod.o \
         $(OBJDIR)/derivatives_mod.o \
         $(OBJDIR)/correction_mod.o \
         $(OBJDIR)/interpolation_mod.o \
         $(OBJDIR)/projection_mod.o \
         $(OBJDIR)/filter_mod.o \
         $(OBJDIR)/io.o \
         $(OBJDIR)/metrics.o \
         $(OBJDIR)/mappings.o \
         $(OBJDIR)/connectivity.o \
         $(OBJDIR)/initialize.o \
         $(OBJDIR)/generic.o \
         $(OBJDIR)/bc_mod.o \
         $(OBJDIR)/bnd_profiles_mod.o \
         $(OBJDIR)/flux.o \
         $(OBJDIR)/parallel.o \
         $(OBJDIR)/time_mod.o \
         $(OBJDIR)/mms_mod.o \
         $(OBJDIR)/channel_mod.o \
         $(OBJDIR)/limiters.o \
         $(OBJDIR)/postproc_mod.o \
         $(OBJDIR)/restart_mod.o \
         $(OBJDIR)/plot3d_mod.o \
         $(OBJDIR)/gmsh_mod.o \
         $(OBJDIR)/cgns_mod.o \
         $(OBJDIR)/averaging_mod.o \
         $(OBJDIR)/reactions_mod.o \
         $(OBJDIR)/turb_inflow_mod.o \
         $(OBJDIR)/cpu_info_mod.o \
         $(OBJDIR)/interpsol_mod.o \
         $(OBJDIR)/mydebug_mod.o

$(OBJDIR)/averaging_mod.o : averaging_mod.f90 \
                  $(OBJDIR)/kind_types.o \
                  $(OBJDIR)/eqn_idx_mod.o \
                  $(OBJDIR)/geovar_mod.o \
                  $(OBJDIR)/ovar_mod.o \
                  $(OBJDIR)/flowvar_mod.o \
                  $(OBJDIR)/io.o \
                  $(OBJDIR)/restart_mod.o \
                  $(OBJDIR)/metrics.o

$(OBJDIR)/bc_mod.o : bc_mod.f90 \
           $(OBJDIR)/kind_types.o \
           $(OBJDIR)/eqn_idx_mod.o \
           $(OBJDIR)/geovar_mod.o \
           $(OBJDIR)/ovar_mod.o \
           $(OBJDIR)/flowvar_mod.o \
           $(OBJDIR)/metrics.o \
           $(OBJDIR)/thermo_mod.o \
           $(OBJDIR)/transport_mod.o \
           $(OBJDIR)/turb_inflow_mod.o \
           $(OBJDIR)/bc_adiabatic_wall_mod.o \
           $(OBJDIR)/bc_fixed_mod.o \
           $(OBJDIR)/bc_flowthrough_grad_mod.o \
           $(OBJDIR)/bc_gen_free_mod.o \
           $(OBJDIR)/bc_gen_in_mod.o \
           $(OBJDIR)/bc_gen_out_mod.o \
           $(OBJDIR)/bc_isothermal_wall_mod.o \
           $(OBJDIR)/bc_mms_mod.o \
           $(OBJDIR)/bc_periodic_mod.o \
           $(OBJDIR)/bc_riemann_mod.o \
           $(OBJDIR)/bc_slip_wall_mod.o \
           $(OBJDIR)/bc_sub_in_mod.o \
           $(OBJDIR)/bc_sub_out_mod.o \
           $(OBJDIR)/bc_sup_out_mod.o \
                     mpi_defs.h

$(OBJDIR)/bnd_profiles_mod.o : bnd_profiles_mod.f90 \
                     $(OBJDIR)/kind_types.o \
                     $(OBJDIR)/eqn_idx_mod.o \
                     $(OBJDIR)/order_mod.o \
                     $(OBJDIR)/geovar_mod.o \
                     $(OBJDIR)/ovar_mod.o \
                     $(OBJDIR)/flowvar_mod.o \
                     $(OBJDIR)/quadrature_mod.o \
                     $(OBJDIR)/interpolation_mod.o \
                     $(OBJDIR)/gpu_mod.o \
                     $(OBJDIR)/thermo_mod.o \
                     $(OBJDIR)/transport_mod.o \
                     $(FUNDIR)/bc_string_to_integer.f90

$(OBJDIR)/cgns_mod.o : cgns_mod.f90 \
             $(OBJDIR)/kind_types.o \
             $(OBJDIR)/cgnstypes_mod.o \
             $(OBJDIR)/geovar_mod.o \
             $(OBJDIR)/ovar_mod.o \
             $(OBJDIR)/order_mod.o \
             $(OBJDIR)/connectivity.o \
             $(FUNDIR)/bc_string_to_integer.f90 \
             $(FUNDIR)/smc000_interior_wall_radius.f90 \
             $(FUNDIR)/last.f90 \
             $(FUNDIR)/debug_timer.f90

$(OBJDIR)/cgnstypes_mod.o : cgnstypes_mod.f90 \
                  $(OBJDIR)/kind_types.o

$(OBJDIR)/channel_mod.o : channel_mod.f90 \
                $(OBJDIR)/kind_types.o \
                $(OBJDIR)/eqn_idx_mod.o \
                $(OBJDIR)/order_mod.o \
                $(OBJDIR)/geovar_mod.o \
                $(OBJDIR)/ovar_mod.o \
                $(OBJDIR)/flowvar_mod.o \
                $(OBJDIR)/quadrature_mod.o \
                $(OBJDIR)/metrics.o \
                $(OBJDIR)/derivatives_mod.o \
                $(OBJDIR)/transport_mod.o \
                $(OBJDIR)/turb_inflow_mod.o

$(OBJDIR)/connectivity.o : connectivity.f90 \
                 $(OBJDIR)/kind_types.o \
                 $(OBJDIR)/geovar_mod.o \
                 $(OBJDIR)/ovar_mod.o \
                 $(OBJDIR)/order_mod.o \
                 $(FUNDIR)/last.f90

$(OBJDIR)/correction_mod.o : correction_mod.f90 \
                   $(OBJDIR)/kind_types.o \
                   $(OBJDIR)/typesgen_mod.o \
                   $(OBJDIR)/order_mod.o \
                   $(OBJDIR)/ovar_mod.o \
                   $(OBJDIR)/quadrature_mod.o \
                   $(OBJDIR)/vandermonde_mod.o \
                   $(OBJDIR)/poly.o \
                   $(FUNDIR)/np2n.f90

$(OBJDIR)/cpu_info_mod.o : cpu_info_mod.f90 \
                 $(OBJDIR)/kind_types.o \
                 $(OBJDIR)/ovar_mod.o

$(OBJDIR)/derivatives_mod.o : derivatives_mod.f90 \
                    $(OBJDIR)/kind_types.o \
                    $(OBJDIR)/typesgen_mod.o \
                    $(OBJDIR)/order_mod.o \
                    $(OBJDIR)/ovar_mod.o \
                    $(OBJDIR)/quadrature_mod.o \
                    $(OBJDIR)/poly.o \
                    $(OBJDIR)/vandermonde_mod.o \
                    $(FUNDIR)/check_for_bad_number.f90 \
                    $(FUNDIR)/np2n.f90 \
                    $(FUNDIR)/debug_timer.f90

$(OBJDIR)/eqn_idx_mod.o : eqn_idx_mod.f90 \
                $(OBJDIR)/kind_types.o \
                $(OBJDIR)/ovar_mod.o

$(OBJDIR)/filter_mod.o : filter_mod.f90 \
               $(OBJDIR)/kind_types.o \
               $(OBJDIR)/eqn_idx_mod.o \
               $(OBJDIR)/typesgen_mod.o \
               $(OBJDIR)/ovar_mod.o \
               $(OBJDIR)/order_mod.o \
               $(OBJDIR)/geovar_mod.o \
               $(OBJDIR)/flowvar_mod.o \
               $(OBJDIR)/quadrature_mod.o \
               $(OBJDIR)/vandermonde_mod.o \
               $(OBJDIR)/thermo_mod.o \
               $(FUNDIR)/solpts2order.f90 \
               $(FUNDIR)/debug_timer.f90

$(OBJDIR)/flowvar_mod.o : flowvar_mod.f90 \
                $(OBJDIR)/kind_types.o \
                $(OBJDIR)/eqn_idx_mod.o \
                $(OBJDIR)/geovar_mod.o \
                $(OBJDIR)/ovar_mod.o \
                $(OBJDIR)/order_mod.o

$(OBJDIR)/flux.o : flux.f90 \
         $(OBJDIR)/kind_types.o \
         $(OBJDIR)/eqn_idx_mod.o \
         $(OBJDIR)/geovar_mod.o \
         $(OBJDIR)/order_mod.o \
         $(OBJDIR)/ovar_mod.o \
         $(OBJDIR)/flowvar_mod.o \
         $(OBJDIR)/metrics.o \
         $(OBJDIR)/interpolation_mod.o \
         $(OBJDIR)/correction_mod.o \
         $(OBJDIR)/derivatives_mod.o \
         $(OBJDIR)/invs_flux_mod.o \
         $(OBJDIR)/visc_flux_mod.o \
         $(OBJDIR)/upwind_flux_mod.o \
         $(OBJDIR)/projection_mod.o \
         $(OBJDIR)/typesgen_mod.o \
         $(OBJDIR)/thermo_mod.o \
         $(OBJDIR)/transport_mod.o \
         $(FUNDIR)/debug_timer.f90 \
         $(FUNDIR)/check_for_bad_number.f90

$(OBJDIR)/generic.o : generic.f90 \
            $(OBJDIR)/kind_types.o \
            $(OBJDIR)/pbs_mod.o \
            $(OBJDIR)/eqn_idx_mod.o \
            $(OBJDIR)/geovar_mod.o \
            $(OBJDIR)/ovar_mod.o \
            $(OBJDIR)/order_mod.o \
            $(OBJDIR)/metrics.o \
            $(OBJDIR)/flowvar_mod.o \
            $(OBJDIR)/quadrature_mod.o \
            $(OBJDIR)/gpu_mod.o \
            $(OBJDIR)/mms_mod.o \
            $(FUNDIR)/debug_timer.f90

$(OBJDIR)/geovar_mod.o : geovar_mod.f90 \
               $(OBJDIR)/kind_types.o

$(OBJDIR)/gmsh_mod.o : gmsh_mod.f90 \
             $(OBJDIR)/kind_types.o \
             $(OBJDIR)/geovar_mod.o \
             $(OBJDIR)/ovar_mod.o \
             $(OBJDIR)/order_mod.o \
             $(FUNDIR)/bc_string_to_integer.f90 \
             $(FUNDIR)/last.f90

$(OBJDIR)/gpu_mod.o : gpu_mod.f90 \
            $(OBJDIR)/kind_types.o \
            $(OBJDIR)/eqn_idx_mod.o \
            $(OBJDIR)/ovar_mod.o \
            $(OBJDIR)/thermo_mod.o

$(OBJDIR)/initialize.o : initialize.f90 \
               $(OBJDIR)/kind_types.o \
               $(OBJDIR)/typesgen_mod.o \
               $(OBJDIR)/eqn_idx_mod.o \
               $(OBJDIR)/geovar_mod.o \
               $(OBJDIR)/ovar_mod.o \
               $(OBJDIR)/order_mod.o \
               $(OBJDIR)/flowvar_mod.o \
               $(OBJDIR)/mms_mod.o \
               $(OBJDIR)/channel_mod.o \
               $(OBJDIR)/metrics.o \
               $(OBJDIR)/quadrature_mod.o \
               $(OBJDIR)/restart_mod.o \
               $(OBJDIR)/projection_mod.o \
               $(OBJDIR)/interpolation_mod.o \
               $(OBJDIR)/derivatives_mod.o \
               $(OBJDIR)/thermo_mod.o \
               $(OBJDIR)/transport_mod.o \
               $(OBJDIR)/reactions_mod.o \
               $(OBJDIR)/turb_inflow_mod.o \
               $(FUNDIR)/vortex.f90 \
               $(FUNDIR)/taylor_green.f90 \
               $(FUNDIR)/smc000_interior_wall_radius.f90 \
               $(FUNDIR)/arn2_interior_wall_radius.f90 \
               $(FUNDIR)/check_for_bad_number.f90

$(OBJDIR)/interpolation_mod.o : interpolation_mod.f90 \
                      $(OBJDIR)/kind_types.o \
                      $(OBJDIR)/typesgen_mod.o \
                      $(OBJDIR)/geovar_mod.o \
                      $(OBJDIR)/ovar_mod.o \
                      $(OBJDIR)/order_mod.o \
                      $(OBJDIR)/quadrature_mod.o \
                      $(OBJDIR)/poly.o \
                      $(OBJDIR)/flowvar_mod.o \
                      $(FUNDIR)/np2n.f90 \
                      $(FUNDIR)/cell_solpts.f90

$(OBJDIR)/interpsol_mod.o : interpsol_mod.f90 \
                  $(OBJDIR)/kind_types.o \
                  $(OBJDIR)/eqn_idx_mod.o \
                  $(OBJDIR)/geovar_mod.o \
                  $(OBJDIR)/ovar_mod.o \
                  $(OBJDIR)/order_mod.o \
                  $(OBJDIR)/flowvar_mod.o \
                  $(OBJDIR)/flux.o \
                  $(OBJDIR)/interpolation_mod.o \
                  $(OBJDIR)/interpsolvar_mod.o \
                  $(OBJDIR)/metrics.o \
                  $(OBJDIR)/octedat_facexyz_mod.o \
                  $(OBJDIR)/octstat_mod.o \
                  $(OBJDIR)/octtree_mod.o \
                  $(OBJDIR)/mydebug_mod.o \
                  $(FUNDIR)/debug_timer.f90

$(OBJDIR)/interpsolvar_mod.o : interpsolvar_mod.f90 \
                     $(OBJDIR)/kind_types.o \
                     $(OBJDIR)/geovar_mod.o \
                     $(OBJDIR)/flowvar_mod.o \
                     $(OBJDIR)/typesgen_mod.o \
                     $(OBJDIR)/quadrature_mod.o \
                     $(OBJDIR)/vandermonde_mod.o \
                     $(OBJDIR)/derivatives_mod.o \
                     $(OBJDIR)/correction_mod.o \
                     $(OBJDIR)/interpolation_mod.o \
                     $(OBJDIR)/projection_mod.o \
                     $(OBJDIR)/filter_mod.o

$(OBJDIR)/invs_flux_mod.o : invs_flux_mod.f90 \
         $(OBJDIR)/kind_types.o \
         $(OBJDIR)/eqn_idx_mod.o \
         $(OBJDIR)/geovar_mod.o \
         $(OBJDIR)/ovar_mod.o

$(OBJDIR)/io.o : io.f90 \
       $(OBJDIR)/kind_types.o \
       $(OBJDIR)/cgnstypes_mod.o \
       $(OBJDIR)/eqn_idx_mod.o \
       $(OBJDIR)/order_mod.o \
       $(OBJDIR)/geovar_mod.o \
       $(OBJDIR)/ovar_mod.o \
       $(OBJDIR)/metrics.o \
       $(OBJDIR)/flowvar_mod.o \
       $(OBJDIR)/gmsh_mod.o \
       $(OBJDIR)/plot3d_mod.o \
       $(OBJDIR)/cgns_mod.o \
       $(OBJDIR)/parallel.o \
       $(OBJDIR)/quadrature_mod.o \
       $(OBJDIR)/interpolation_mod.o \
       $(OBJDIR)/initialize.o \
       $(OBJDIR)/restart_mod.o \
       $(OBJDIR)/limiters.o \
       $(OBJDIR)/thermo_mod.o \
       $(OBJDIR)/transport_mod.o \
       $(OBJDIR)/turb_inflow_mod.o \
       $(OBJDIR)/flux.o \
       $(OBJDIR)/gpu_mod.o \
       $(FUNDIR)/last.f90 \
       $(FUNDIR)/bc_priority.f90 \
       $(FUNDIR)/cell_vol.f90 \
       $(FUNDIR)/bc_string_to_integer.f90 \
       $(FUNDIR)/check_for_bad_number.f90 \
       $(FUNDIR)/cgns_write_field_data.f90 \
       $(FUNDIR)/cgns_write_field_header.f90 \
                 mpi_defs.h

$(OBJDIR)/kind_types.o : kind_types.f90 \
               $(FUNDIR)/all_are_equal.f90 \
               $(FUNDIR)/bc_integer_to_string.f90 \
               $(FUNDIR)/big_to_little_endian.f90 \
               $(FUNDIR)/checks_for_bad_numbers.f90 \
               $(FUNDIR)/chop.f90 \
               $(FUNDIR)/cross_product.f90 \
               $(FUNDIR)/deltafun.f90 \
               $(FUNDIR)/face_normal.f90 \
               $(FUNDIR)/findloc.f90 \
               $(FUNDIR)/ghost_nodes.f90 \
               $(FUNDIR)/init_array.f90 \
               $(FUNDIR)/intseq.f90 \
               $(FUNDIR)/jacobians.f90 \
               $(FUNDIR)/matrix_ops.f90 \
               $(FUNDIR)/node_center.f90 \
               $(FUNDIR)/norms.f90 \
               $(FUNDIR)/outerprod.f90 \
               $(FUNDIR)/random.f90 \
               $(FUNDIR)/reallocate.f90 \
               $(FUNDIR)/reciprocal.f90 \
               $(FUNDIR)/reflect.f90 \
               $(FUNDIR)/reformat_string.f90 \
               $(FUNDIR)/reverse_metrics.f90 \
               $(FUNDIR)/rotate.f90 \
               $(FUNDIR)/rotate_2d.f90 \
               $(FUNDIR)/uppercase.f90 \
               $(FUNDIR)/debug_timer.f90 \
                         mpi_defs.h

$(OBJDIR)/limiters.o : limiters.f90 \
             $(OBJDIR)/kind_types.o \
             $(OBJDIR)/eqn_idx_mod.o \
             $(OBJDIR)/geovar_mod.o \
             $(OBJDIR)/ovar_mod.o \
             $(OBJDIR)/order_mod.o \
             $(OBJDIR)/flowvar_mod.o \
             $(OBJDIR)/metrics.o \
             $(OBJDIR)/vandermonde_mod.o \
             $(OBJDIR)/projection_mod.o \
             $(OBJDIR)/filter_mod.o \
             $(OBJDIR)/quadrature_mod.o \
             $(FUNDIR)/debug_timer.f90

$(OBJDIR)/mappings.o : mappings.f90 \
             $(OBJDIR)/kind_types.o \
             $(OBJDIR)/eqn_idx_mod.o \
             $(OBJDIR)/geovar_mod.o \
             $(OBJDIR)/ovar_mod.o \
             $(OBJDIR)/order_mod.o \
             $(OBJDIR)/flowvar_mod.o \
             $(OBJDIR)/quadrature_mod.o \
             $(OBJDIR)/typesgen_mod.o \
             $(OBJDIR)/poly.o

$(OBJDIR)/metis.o : metis.f90

$(OBJDIR)/metrics.o : metrics.f90 \
            $(OBJDIR)/kind_types.o \
            $(OBJDIR)/geovar_mod.o \
            $(OBJDIR)/order_mod.o \
            $(OBJDIR)/ovar_mod.o \
            $(OBJDIR)/quadrature_mod.o \
            $(OBJDIR)/interpolation_mod.o \
            $(OBJDIR)/derivatives_mod.o

$(OBJDIR)/mms_mod.o : mms_mod.f90 \
            $(OBJDIR)/kind_types.o \
            $(OBJDIR)/eqn_idx_mod.o \
            $(OBJDIR)/geovar_mod.o \
            $(OBJDIR)/ovar_mod.o \
            $(OBJDIR)/flowvar_mod.o \
            $(OBJDIR)/metrics.o \
            $(OBJDIR)/gpu_mod.o

$(OBJDIR)/mycpu.o : $(FUNDIR)/mycpu.c

$(OBJDIR)/mydebug_mod.o : mydebug_mod.f90 \
                $(OBJDIR)/kind_types.o \
                $(OBJDIR)/eqn_idx_mod.o \
                $(OBJDIR)/geovar_mod.o \
                $(OBJDIR)/flowvar_mod.o \
                $(OBJDIR)/interpsolvar_mod.o

$(OBJDIR)/octedat_facexyz_mod.o : octedat_facexyz_mod.f90 \
                        $(OBJDIR)/kind_types.o \
                        $(OBJDIR)/geovar_mod.o \
                        $(OBJDIR)/flowvar_mod.o \
                        $(OBJDIR)/octedat_mod.o \
                        $(OBJDIR)/octstat_mod.o

$(OBJDIR)/octedat_mod.o : octedat_mod.f90 \
                $(OBJDIR)/kind_types.o \
                $(OBJDIR)/octstat_mod.o

$(OBJDIR)/octnode_mod.o : octnode_mod.f90 \
                $(OBJDIR)/kind_types.o \
                $(OBJDIR)/geovar_mod.o \
                $(OBJDIR)/octedat_mod.o \
                $(OBJDIR)/octroot_mod.o \
                $(OBJDIR)/octstat_mod.o

$(OBJDIR)/octroot_mod.o : octroot_mod.f90 \
                $(OBJDIR)/octstat_mod.o

$(OBJDIR)/octstat_mod.o : octstat_mod.f90

$(OBJDIR)/octtree_mod.o : octtree_mod.f90 \
                $(OBJDIR)/kind_types.o \
                $(OBJDIR)/geovar_mod.o \
                $(OBJDIR)/octedat_mod.o \
                $(OBJDIR)/octnode_mod.o \
                $(OBJDIR)/octroot_mod.o \
                $(OBJDIR)/octstat_mod.o \
                $(FUNDIR)/debug_timer.f90

$(OBJDIR)/order_mod.o : order_mod.f90 \
              $(OBJDIR)/kind_types.o \
              $(FUNDIR)/cell_solpts.f90

$(OBJDIR)/ovar_mod.o : ovar_mod.f90 \
             $(OBJDIR)/kind_types.o \
             $(OBJDIR)/cgnstypes_mod.o

$(OBJDIR)/parallel.o : parallel.f90 \
             $(OBJDIR)/kind_types.o \
             $(OBJDIR)/eqn_idx_mod.o \
             $(OBJDIR)/metis.o \
             $(OBJDIR)/geovar_mod.o \
             $(OBJDIR)/ovar_mod.o \
             $(OBJDIR)/order_mod.o \
             $(OBJDIR)/flowvar_mod.o \
             $(OBJDIR)/connectivity.o \
             $(OBJDIR)/gmsh_mod.o \
             $(OBJDIR)/cgns_mod.o \
             $(OBJDIR)/plot3d_mod.o \
             $(FUNDIR)/cell_solpts.f90 \
             $(FUNDIR)/last.f90 \
             $(FUNDIR)/debug_timer.f90 \
                       mpi_defs.h

$(OBJDIR)/pbs_mod.o : pbs_mod.f90 \
            $(OBJDIR)/kind_types.o

$(OBJDIR)/pbs_time_left.o : $(FUNDIR)/pbs_time_left.c

$(OBJDIR)/plot3d_mod.o : plot3d_mod.f90 \
               $(OBJDIR)/kind_types.o \
               $(OBJDIR)/geovar_mod.o \
               $(OBJDIR)/ovar_mod.o \
               $(OBJDIR)/order_mod.o \
               $(OBJDIR)/plot3dio_mod.o \
               $(FUNDIR)/bc_string_to_integer.f90 \
               $(FUNDIR)/last.f90
#              $(FUNDIR)/check_for_bad_number.f90 \

$(OBJDIR)/plot3dio_mod.o : plot3dio_mod.f90 \
                 $(OBJDIR)/kind_types.o

$(OBJDIR)/poly.o : poly.f90 \
         $(OBJDIR)/kind_types.o \
         $(FUNDIR)/check_for_bad_number.f90

$(OBJDIR)/postproc_mod.o : postproc_mod.f90 \
                 $(OBJDIR)/kind_types.o \
                 $(OBJDIR)/cgnstypes_mod.o \
                 $(OBJDIR)/eqn_idx_mod.o \
                 $(OBJDIR)/geovar_mod.o \
                 $(OBJDIR)/ovar_mod.o \
                 $(OBJDIR)/flowvar_mod.o \
                 $(OBJDIR)/bc_mod.o \
                 $(OBJDIR)/order_mod.o \
                 $(OBJDIR)/interpolation_mod.o \
                 $(OBJDIR)/parallel.o \
                 $(OBJDIR)/connectivity.o \
                 $(OBJDIR)/thermo_mod.o \
                 $(OBJDIR)/transport_mod.o \
                 $(OBJDIR)/gpu_mod.o \
                 $(OBJDIR)/metrics.o \
                 $(OBJDIR)/limiters.o \
                 $(FUNDIR)/cgns_write_field_data.f90 \
                 $(FUNDIR)/cgns_write_field_header.f90 \
                 $(FUNDIR)/last.f90 \
                 $(FUNDIR)/cell_solpts.f90


$(OBJDIR)/projection_mod.o : projection_mod.f90 \
                   $(OBJDIR)/kind_types.o \
                   $(OBJDIR)/typesgen_mod.o \
                   $(OBJDIR)/geovar_mod.o \
                   $(OBJDIR)/order_mod.o \
                   $(OBJDIR)/ovar_mod.o \
                   $(OBJDIR)/quadrature_mod.o \
                   $(OBJDIR)/vandermonde_mod.o \
                   $(FUNDIR)/debug_timer.f90

$(OBJDIR)/quadrature_mod.o : quadrature_mod.f90 \
                   $(OBJDIR)/kind_types.o \
                   $(OBJDIR)/geovar_mod.o \
                   $(OBJDIR)/order_mod.o \
                   $(OBJDIR)/ovar_mod.o \
                   $(OBJDIR)/poly.o \
                   $(OBJDIR)/triangle.o \
                   $(FUNDIR)/np2n.f90 \
                   $(FUNDIR)/check_for_bad_number.f90 \
                   $(FUNDIR)/debug_timer.f90

$(OBJDIR)/random_mod.o : random_mod.f90 \
             $(OBJDIR)/kind_types.o

$(OBJDIR)/reactions_mod.o : reactions_mod.f90 \
                  $(OBJDIR)/kind_types.o \
                  $(OBJDIR)/eqn_idx_mod.o \
                  $(OBJDIR)/ovar_mod.o \
                  $(OBJDIR)/flowvar_mod.o \
                  $(OBJDIR)/metrics.o \
                  $(OBJDIR)/thermo_mod.o \
                  $(FUNDIR)/check_for_bad_number.f90

$(OBJDIR)/restart_mod.o : restart_mod.f90 \
                $(OBJDIR)/kind_types.o \
                $(OBJDIR)/eqn_idx_mod.o \
                $(OBJDIR)/geovar_mod.o \
                $(OBJDIR)/ovar_mod.o \
                $(OBJDIR)/flowvar_mod.o \
                $(OBJDIR)/random_mod.o \
                $(OBJDIR)/parallel.o \
                $(OBJDIR)/thermo_mod.o \
                $(OBJDIR)/turb_inflow_mod.o \
                $(FUNDIR)/cell_solpts.f90 \
                $(FUNDIR)/sort.f90 \
                          mpi_defs.h

$(OBJDIR)/thermo_mod.o : thermo_mod.f90 \
               $(OBJDIR)/kind_types.o \
               $(OBJDIR)/eqn_idx_mod.o \
               $(OBJDIR)/ovar_mod.o \
               $(OBJDIR)/flowvar_mod.o

$(OBJDIR)/time_mod.o : time_mod.f90 \
             $(OBJDIR)/kind_types.o \
             $(OBJDIR)/eqn_idx_mod.o \
             $(OBJDIR)/geovar_mod.o \
             $(OBJDIR)/ovar_mod.o \
             $(OBJDIR)/order_mod.o \
             $(OBJDIR)/flowvar_mod.o \
             $(OBJDIR)/metrics.o \
             $(OBJDIR)/filter_mod.o \
             $(OBJDIR)/limiters.o \
             $(OBJDIR)/projection_mod.o \
             $(OBJDIR)/thermo_mod.o \
             $(OBJDIR)/transport_mod.o \
             $(FUNDIR)/debug_timer.f90

$(OBJDIR)/transport_mod.o : transport_mod.f90 \
                  $(OBJDIR)/kind_types.o \
                  $(OBJDIR)/eqn_idx_mod.o \
                  $(OBJDIR)/ovar_mod.o \
                  $(OBJDIR)/thermo_mod.o

$(OBJDIR)/triangle.o : triangle.f90 \
             $(OBJDIR)/kind_types.o \
             $(OBJDIR)/poly.o

$(OBJDIR)/turb_inflow_mod.o : turb_inflow_mod.f90 \
                    $(OBJDIR)/kind_types.o \
                    $(OBJDIR)/random_mod.o \
                    $(OBJDIR)/eqn_idx_mod.o \
                    $(OBJDIR)/geovar_mod.o \
                    $(OBJDIR)/ovar_mod.o \
                    $(OBJDIR)/flowvar_mod.o \
                    $(OBJDIR)/metrics.o \
                    $(OBJDIR)/thermo_mod.o \
                    $(OBJDIR)/transport_mod.o

$(OBJDIR)/typesgen_mod.o : typesgen_mod.f90 \
                 $(OBJDIR)/kind_types.o \
                 $(OBJDIR)/ovar_mod.o

$(OBJDIR)/upwind_flux_mod.o : upwind_flux_mod.f90 \
         $(OBJDIR)/kind_types.o \
         $(OBJDIR)/eqn_idx_mod.o \
         $(OBJDIR)/geovar_mod.o \
         $(OBJDIR)/ovar_mod.o

$(OBJDIR)/vandermonde_mod.o : vandermonde_mod.f90 \
                    $(OBJDIR)/kind_types.o \
                    $(OBJDIR)/typesgen_mod.o \
                    $(OBJDIR)/order_mod.o \
                    $(OBJDIR)/quadrature_mod.o \
                    $(OBJDIR)/poly.o \
                    $(FUNDIR)/np2n.f90 \
                    $(FUNDIR)/cell_solpts.f90 \
                    $(FUNDIR)/debug_timer.f90

$(OBJDIR)/visc_flux_mod.o : visc_flux_mod.f90 \
         $(OBJDIR)/kind_types.o \
         $(OBJDIR)/eqn_idx_mod.o \
         $(OBJDIR)/geovar_mod.o \
         $(OBJDIR)/ovar_mod.o \
         $(OBJDIR)/gpu_mod.o \
         $(OBJDIR)/thermo_mod.o \
         $(OBJDIR)/transport_mod.o

$(OBJDIR)/bc_adiabatic_wall_mod.o : $(BCDIR)/bc_adiabatic_wall_mod.f90 \
                          $(OBJDIR)/kind_types.o \
                          $(OBJDIR)/eqn_idx_mod.o \
                          $(OBJDIR)/ovar_mod.o \
                          $(OBJDIR)/gpu_mod.o \
                          $(OBJDIR)/thermo_mod.o \
                          $(OBJDIR)/bc_relax_mod.o

$(OBJDIR)/bc_fixed_mod.o : $(BCDIR)/bc_fixed_mod.f90 \
                 $(OBJDIR)/kind_types.o \
                 $(OBJDIR)/eqn_idx_mod.o

$(OBJDIR)/bc_flowthrough_grad_mod.o : $(BCDIR)/bc_flowthrough_grad_mod.f90 \
                            $(OBJDIR)/kind_types.o \
                            $(OBJDIR)/eqn_idx_mod.o

$(OBJDIR)/bc_gen_free_mod.o : $(BCDIR)/bc_gen_free_mod.f90 \
                    $(OBJDIR)/kind_types.o \
                    $(OBJDIR)/eqn_idx_mod.o \
                    $(OBJDIR)/ovar_mod.o \
                    $(OBJDIR)/thermo_mod.o \
                    $(OBJDIR)/bc_relax_mod.o

$(OBJDIR)/bc_gen_in_mod.o : $(BCDIR)/bc_gen_in_mod.f90 \
                  $(OBJDIR)/kind_types.o \
                  $(OBJDIR)/eqn_idx_mod.o \
                  $(OBJDIR)/ovar_mod.o \
                  $(OBJDIR)/thermo_mod.o \
                  $(OBJDIR)/bc_relax_mod.o

$(OBJDIR)/bc_gen_out_mod.o : $(BCDIR)/bc_gen_out_mod.f90 \
                   $(OBJDIR)/kind_types.o \
                   $(OBJDIR)/eqn_idx_mod.o \
                   $(OBJDIR)/ovar_mod.o \
                   $(OBJDIR)/bc_relax_mod.o

$(OBJDIR)/bc_isothermal_wall_mod.o : $(BCDIR)/bc_isothermal_wall_mod.f90 \
                           $(OBJDIR)/kind_types.o \
                           $(OBJDIR)/eqn_idx_mod.o \
                           $(OBJDIR)/ovar_mod.o \
                           $(OBJDIR)/gpu_mod.o \
                           $(OBJDIR)/thermo_mod.o

$(OBJDIR)/bc_mms_mod.o : $(BCDIR)/bc_mms_mod.f90 \
               $(OBJDIR)/kind_types.o \
               $(OBJDIR)/eqn_idx_mod.o \
               $(OBJDIR)/ovar_mod.o \
               $(OBJDIR)/mms_mod.o

$(OBJDIR)/bc_periodic_mod.o : $(BCDIR)/bc_periodic_mod.f90 \
                    $(OBJDIR)/kind_types.o \
                    $(OBJDIR)/eqn_idx_mod.o

$(OBJDIR)/bc_relax_mod.o : $(BCDIR)/bc_relax_mod.f90 \
                 $(OBJDIR)/kind_types.o \
                 $(OBJDIR)/eqn_idx_mod.o \
                 $(OBJDIR)/ovar_mod.o \
                 $(OBJDIR)/thermo_mod.o \
                 $(OBJDIR)/transport_mod.o

$(OBJDIR)/bc_riemann_mod.o : $(BCDIR)/bc_riemann_mod.f90 \
                   $(OBJDIR)/kind_types.o \
                   $(OBJDIR)/eqn_idx_mod.o \
                   $(OBJDIR)/ovar_mod.o \
                   $(OBJDIR)/bc_relax_mod.o

$(OBJDIR)/bc_slip_wall_mod.o : $(BCDIR)/bc_slip_wall_mod.f90 \
                     $(OBJDIR)/kind_types.o \
                     $(OBJDIR)/eqn_idx_mod.o \
                     $(OBJDIR)/ovar_mod.o \
                     $(OBJDIR)/gpu_mod.o

$(OBJDIR)/bc_sub_in_mod.o : $(BCDIR)/bc_sub_in_mod.f90 \
                  $(OBJDIR)/kind_types.o \
                  $(OBJDIR)/eqn_idx_mod.o \
                  $(OBJDIR)/order_mod.o \
                  $(OBJDIR)/geovar_mod.o \
                  $(OBJDIR)/ovar_mod.o \
                  $(OBJDIR)/interpolation_mod.o \
                  $(OBJDIR)/thermo_mod.o \
                  $(OBJDIR)/bc_relax_mod.o

$(OBJDIR)/bc_sub_out_mod.o : $(BCDIR)/bc_sub_out_mod.f90 \
                   $(OBJDIR)/kind_types.o \
                   $(OBJDIR)/eqn_idx_mod.o \
                   $(OBJDIR)/ovar_mod.o \
                   $(OBJDIR)/bc_relax_mod.o

$(OBJDIR)/bc_sup_out_mod.o : $(BCDIR)/bc_sup_out_mod.f90 \
                   $(OBJDIR)/kind_types.o \
                   $(OBJDIR)/eqn_idx_mod.o

$(OBJDIR)/tecio_mod.o : tecio_mod.f90

