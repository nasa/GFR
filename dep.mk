$(OBJDIR)/main.o : main.f90 \
         $(OBJDIR)/kind_types.o \
         $(OBJDIR)/pbs_mod.o \
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
         $(OBJDIR)/cpu_info_mod.o

$(OBJDIR)/averaging_mod.o : averaging_mod.f90 \
                  $(OBJDIR)/kind_types.o \
                  $(OBJDIR)/geovar_mod.o \
                  $(OBJDIR)/ovar_mod.o \
                  $(OBJDIR)/flowvar_mod.o \
                  $(OBJDIR)/io.o \
                  $(OBJDIR)/time_mod.o \
                  $(OBJDIR)/restart_mod.o \
                  $(OBJDIR)/metrics.o

$(OBJDIR)/bc_mod.o : bc_mod.f90 \
           $(OBJDIR)/kind_types.o \
           $(OBJDIR)/eqn_idx_mod.o \
           $(OBJDIR)/geovar_mod.o \
           $(OBJDIR)/ovar_mod.o \
           $(OBJDIR)/metrics.o \
           $(OBJDIR)/flowvar_mod.o \
           $(OBJDIR)/interpolation_mod.o \
           $(OBJDIR)/parallel.o \
           $(OBJDIR)/mms_mod.o \
           $(FUNDIR)/pressure.f90 \
           $(FUNDIR)/temperature.f90 \
           $(FUNDIR)/usp2v.f90 \
           $(FUNDIR)/vsp2u.f90 \
                     mpi_defs.h

$(OBJDIR)/bnd_profiles_mod.o : bnd_profiles_mod.f90 \
                     $(OBJDIR)/kind_types.o \
                     $(OBJDIR)/eqn_idx_mod.o \
                     $(OBJDIR)/order_mod.o \
                     $(OBJDIR)/geovar_mod.o \
                     $(OBJDIR)/ovar_mod.o \
                     $(OBJDIR)/flowvar_mod.o \
                     $(OBJDIR)/metrics.o \
                     $(OBJDIR)/interpolation_mod.o \
                     $(FUNDIR)/grad_cv_to_grad_pv.f90 \
                     $(FUNDIR)/normal_streamwise_velocity_gradient.f90 \
                     $(FUNDIR)/bc_string_to_integer.f90 \
                     $(FUNDIR)/temperature.f90 \
                     $(FUNDIR)/viscosity.f90 \
                     $(FUNDIR)/usp2v.f90

$(OBJDIR)/cgns_mod.o : cgns_mod.f90 \
             $(OBJDIR)/kind_types.o \
             $(OBJDIR)/cgnstypes_mod.o \
             $(OBJDIR)/geovar_mod.o \
             $(OBJDIR)/ovar_mod.o \
             $(OBJDIR)/order_mod.o \
             $(OBJDIR)/connectivity.o \
             $(FUNDIR)/bc_string_to_integer.f90 \
             $(FUNDIR)/smc000_interior_wall_radius.f90 \
             $(FUNDIR)/last.f90

$(OBJDIR)/cgnstypes_mod.o : cgnstypes_mod.f90 \
                  $(OBJDIR)/kind_types.o

$(OBJDIR)/channel_mod.o : channel_mod.f90 \
                $(OBJDIR)/kind_types.o \
                $(OBJDIR)/eqn_idx_mod.o \
                $(OBJDIR)/geovar_mod.o \
                $(OBJDIR)/ovar_mod.o \
                $(OBJDIR)/flowvar_mod.o \
                $(OBJDIR)/metrics.o

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
                   $(OBJDIR)/geovar_mod.o \
                   $(OBJDIR)/quadrature_mod.o \
                   $(OBJDIR)/vandermonde_mod.o \
                   $(OBJDIR)/poly.o \
                   $(FUNDIR)/np2n.f90

$(OBJDIR)/cpu_info_mod.o : cpu_info_mod.f90 \
                 $(OBJDIR)/kind_types.o \
                 $(OBJDIR)/ovar_mod.o \
                 $(OBJDIR)/mycpu.o

$(OBJDIR)/derivatives_mod.o : derivatives_mod.f90 \
                    $(OBJDIR)/kind_types.o \
                    $(OBJDIR)/typesgen_mod.o \
                    $(OBJDIR)/geovar_mod.o \
                    $(OBJDIR)/order_mod.o \
                    $(OBJDIR)/ovar_mod.o \
                    $(OBJDIR)/quadrature_mod.o \
                    $(OBJDIR)/poly.o \
                    $(OBJDIR)/vandermonde_mod.o \
                    $(FUNDIR)/np2n.f90

$(OBJDIR)/eqn_idx_mod.o : eqn_idx_mod.f90 \
                $(OBJDIR)/kind_types.o \
                $(OBJDIR)/ovar_mod.o

$(OBJDIR)/filter_mod.o : filter_mod.f90 \
               $(OBJDIR)/kind_types.o \
               $(OBJDIR)/typesgen_mod.o \
               $(OBJDIR)/ovar_mod.o \
               $(OBJDIR)/order_mod.o \
               $(OBJDIR)/geovar_mod.o \
               $(OBJDIR)/flowvar_mod.o \
               $(OBJDIR)/quadrature_mod.o \
               $(OBJDIR)/vandermonde_mod.o \
               $(FUNDIR)/solpts2order.f90

$(OBJDIR)/flowvar_mod.o : flowvar_mod.f90 \
                $(OBJDIR)/kind_types.o \
                $(OBJDIR)/eqn_idx_mod.o \
                $(OBJDIR)/geovar_mod.o \
                $(OBJDIR)/ovar_mod.o \
                $(OBJDIR)/order_mod.o


$(OBJDIR)/flux.o : flux.f90 \
         $(OBJDIR)/kind_types.o \
         $(OBJDIR)/typesgen_mod.o \
         $(OBJDIR)/eqn_idx_mod.o \
         $(OBJDIR)/geovar_mod.o \
         $(OBJDIR)/order_mod.o \
         $(OBJDIR)/ovar_mod.o \
         $(OBJDIR)/flowvar_mod.o \
         $(OBJDIR)/bc_mod.o \
         $(OBJDIR)/metrics.o \
         $(OBJDIR)/parallel.o \
         $(OBJDIR)/quadrature_mod.o \
         $(OBJDIR)/interpolation_mod.o \
         $(OBJDIR)/correction_mod.o \
         $(OBJDIR)/derivatives_mod.o \
         $(OBJDIR)/projection_mod.o \
         $(FUNDIR)/upwind_fluxes.f90 \
         $(FUNDIR)/inviscid_fluxes.f90 \
         $(FUNDIR)/viscous_fluxes.f90 \
         $(FUNDIR)/grad_cv_to_grad_pv.f90 \
         $(FUNDIR)/viscosity.f90

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
            $(OBJDIR)/interpolation_mod.o \
            $(OBJDIR)/parallel.o \
            $(OBJDIR)/mms_mod.o \
            $(FUNDIR)/entropy.f90 \
            $(FUNDIR)/pressure.f90 \
            $(FUNDIR)/vortex.f90 \
            $(FUNDIR)/grad_cv_to_grad_pv.f90 \
            $(FUNDIR)/usp2v.f90

$(OBJDIR)/geovar_mod.o : geovar_mod.f90 \
               $(OBJDIR)/kind_types.o

$(OBJDIR)/gmsh_mod.o : gmsh_mod.f90 \
             $(OBJDIR)/kind_types.o \
             $(OBJDIR)/geovar_mod.o \
             $(OBJDIR)/ovar_mod.o \
             $(OBJDIR)/order_mod.o \
             $(FUNDIR)/bc_string_to_integer.f90 \
             $(FUNDIR)/last.f90

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
               $(OBJDIR)/restart_mod.o \
               $(OBJDIR)/projection_mod.o \
               $(OBJDIR)/interpolation_mod.o \
               $(OBJDIR)/derivatives_mod.o \
               $(FUNDIR)/vortex.f90 \
               $(FUNDIR)/taylor_green.f90 \
               $(FUNDIR)/smc000_interior_wall_radius.f90 \
               $(FUNDIR)/arn2_interior_wall_radius.f90 \
               $(FUNDIR)/usp2v.f90 \
               $(FUNDIR)/viscosity.f90
              #$(OBJDIR)/time_mod.o \

$(OBJDIR)/interpolation_mod.o : interpolation_mod.f90 \
                      $(OBJDIR)/kind_types.o \
                      $(OBJDIR)/typesgen_mod.o \
                      $(OBJDIR)/geovar_mod.o \
		      $(OBJDIR)/order_mod.o \
                      $(OBJDIR)/quadrature_mod.o \
                      $(OBJDIR)/poly.o \
                      $(OBJDIR)/flowvar_mod.o \
                      $(FUNDIR)/np2n.f90

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
       $(OBJDIR)/derivatives_mod.o \
       $(OBJDIR)/initialize.o \
       $(OBJDIR)/channel_mod.o \
       $(OBJDIR)/restart_mod.o \
       $(OBJDIR)/limiters.o \
       $(FUNDIR)/last.f90 \
       $(FUNDIR)/bc_priority.f90 \
       $(FUNDIR)/entropy.f90 \
       $(FUNDIR)/pressure.f90 \
       $(FUNDIR)/viscosity.f90 \
       $(FUNDIR)/temperature.f90 \
       $(FUNDIR)/mach_number.f90 \
       $(FUNDIR)/cell_vol.f90 \
       $(FUNDIR)/grad_cv_to_grad_pv.f90 \
       $(FUNDIR)/bc_string_to_integer.f90 \
       $(FUNDIR)/cgns_write_field_parallel.f90 \
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
                         mpi_defs.h

$(OBJDIR)/limiters.o : limiters.f90 \
             $(OBJDIR)/kind_types.o \
             $(OBJDIR)/eqn_idx_mod.o \
             $(OBJDIR)/geovar_mod.o \
             $(OBJDIR)/order_mod.o \
             $(OBJDIR)/flowvar_mod.o \
             $(OBJDIR)/vandermonde_mod.o \
             $(OBJDIR)/projection_mod.o \
             $(OBJDIR)/filter_mod.o \
             $(OBJDIR)/parallel.o \
             $(FUNDIR)/pressure.f90

$(OBJDIR)/mappings.o : mappings.f90 \
             $(OBJDIR)/kind_types.o \
             $(OBJDIR)/geovar_mod.o \
             $(OBJDIR)/order_mod.o \
             $(OBJDIR)/quadrature_mod.o \
             $(OBJDIR)/parallel.o \
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
            $(FUNDIR)/vsp2u.f90

$(OBJDIR)/mycpu.o : $(FUNDIR)/mycpu.c

$(OBJDIR)/order_mod.o : order_mod.f90 \
              $(OBJDIR)/kind_types.o \
              $(FUNDIR)/cell_solpts.f90

$(OBJDIR)/ovar_mod.o : ovar_mod.f90 \
             $(OBJDIR)/kind_types.o \
             $(OBJDIR)/cgnstypes_mod.o \
             $(OBJDIR)/order_mod.o

$(OBJDIR)/parallel.o : parallel.f90 \
             $(OBJDIR)/kind_types.o \
             $(OBJDIR)/eqn_idx_mod.o \
             $(OBJDIR)/metis.o \
             $(OBJDIR)/geovar_mod.o \
             $(OBJDIR)/order_mod.o \
             $(OBJDIR)/flowvar_mod.o \
             $(OBJDIR)/interpolation_mod.o \
             $(OBJDIR)/connectivity.o \
             $(OBJDIR)/gmsh_mod.o \
             $(OBJDIR)/cgns_mod.o \
             $(OBJDIR)/plot3d_mod.o \
             $(FUNDIR)/cell_solpts.f90 \
             $(FUNDIR)/last.f90 \
                       mpi_defs.h

$(OBJDIR)/pbs_mod.o : pbs_mod.f90 \
            $(OBJDIR)/kind_types.o

$(OBJDIR)/plot3d_mod.o : plot3d_mod.f90 \
               $(OBJDIR)/kind_types.o \
               $(OBJDIR)/geovar_mod.o \
               $(OBJDIR)/ovar_mod.o \
               $(OBJDIR)/order_mod.o \
               $(FUNDIR)/bc_string_to_integer.f90 \
               $(FUNDIR)/last.f90

$(OBJDIR)/poly.o : poly.f90 \
         $(OBJDIR)/kind_types.o

$(OBJDIR)/postproc_mod.o : postproc_mod.f90 \
                 $(OBJDIR)/kind_types.o \
                 $(OBJDIR)/cgnstypes_mod.o \
                 $(OBJDIR)/eqn_idx_mod.o \
                 $(OBJDIR)/geovar_mod.o \
                 $(OBJDIR)/flowvar_mod.o \
                 $(OBJDIR)/flux.o \
                 $(OBJDIR)/bc_mod.o \
                 $(OBJDIR)/order_mod.o \
                 $(OBJDIR)/interpolation_mod.o \
                 $(OBJDIR)/parallel.o \
                 $(OBJDIR)/connectivity.o \
                 $(FUNDIR)/last.f90 \
                 $(FUNDIR)/pressure.f90 \
                 $(FUNDIR)/temperature.f90 \
                 $(FUNDIR)/mach_number.f90 \
                 $(FUNDIR)/entropy.f90 \
                 $(FUNDIR)/grad_cv_to_grad_pv.f90 \
                 $(FUNDIR)/cell_solpts.f90


$(OBJDIR)/projection_mod.o : projection_mod.f90 \
                   $(OBJDIR)/kind_types.o \
                   $(OBJDIR)/typesgen_mod.o \
                   $(OBJDIR)/geovar_mod.o \
                   $(OBJDIR)/order_mod.o \
                   $(OBJDIR)/ovar_mod.o \
                   $(OBJDIR)/metrics.o \
                   $(OBJDIR)/poly.o \
                   $(OBJDIR)/quadrature_mod.o \
                   $(OBJDIR)/vandermonde_mod.o \
                   $(OBJDIR)/interpolation_mod.o

$(OBJDIR)/quadrature_mod.o : quadrature_mod.f90 \
                   $(OBJDIR)/kind_types.o \
                   $(OBJDIR)/geovar_mod.o \
                   $(OBJDIR)/order_mod.o \
                   $(OBJDIR)/ovar_mod.o \
                   $(OBJDIR)/poly.o \
                   $(OBJDIR)/triangle.o \
                   $(FUNDIR)/np2n.f90

$(OBJDIR)/restart_mod.o : restart_mod.f90 \
                $(OBJDIR)/kind_types.o \
                $(OBJDIR)/eqn_idx_mod.o \
                $(OBJDIR)/geovar_mod.o \
                $(OBJDIR)/ovar_mod.o \
                $(OBJDIR)/flowvar_mod.o \
                $(OBJDIR)/parallel.o \
                $(FUNDIR)/cell_solpts.f90 \
                $(FUNDIR)/sort.f90 \
                          mpi_defs.h

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
             $(OBJDIR)/interpolation_mod.o \
             $(FUNDIR)/usp2v.f90 \
             $(FUNDIR)/temperature.f90 \
             $(FUNDIR)/viscosity.f90

$(OBJDIR)/triangle.o : triangle.f90 \
             $(OBJDIR)/kind_types.o \
             $(OBJDIR)/poly.o

$(OBJDIR)/typesgen_mod.o : typesgen_mod.f90 \
                 $(OBJDIR)/kind_types.o \
                 $(OBJDIR)/ovar_mod.o

$(OBJDIR)/vandermonde_mod.o : vandermonde_mod.f90 \
                    $(OBJDIR)/kind_types.o \
                    $(OBJDIR)/typesgen_mod.o \
                    $(OBJDIR)/order_mod.o \
                    $(OBJDIR)/geovar_mod.o \
                    $(OBJDIR)/quadrature_mod.o \
                    $(OBJDIR)/poly.o \
                    $(FUNDIR)/np2n.f90

