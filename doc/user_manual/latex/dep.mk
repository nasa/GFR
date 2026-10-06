FILES = $(MainFile) \
        preamble.tex \
        latex/gfr.sty \
        latex/turbulence.sty \
        overview/overview.tex \
        intro/introduction.tex \
        geometry_mesh/geometry_mesh.tex \
        reference_frame/reference_frame.tex \
        reference_conditions/reference_conditions.tex \
        physical_models/physical_models.tex $(PhysModelsFiles) \
        BCs/boundary_conditions.tex $(BCsFiles) \
        time_stepping/time_stepping.tex $(TimeSteppingFiles) \
        numerics/numerics.tex $(NumericsFiles) \
        files/files.tex $(FilesFiles) \
        parallel/parallel.tex \
        postprocessing/postprocessing.tex \
        nml/nml.tex $(NMLFiles) \
        test_cases/test_cases.tex $(TestCasesFiles) \
        limitations/limitations.tex \
        appendix/appendix.tex \
        tikz/EnergySpectra.pdf \
        tikz/couette_schematic.pdf \
        tikz/Jump_Marker.pdf \
        tikz/sem_partition_boundary.pdf \
        tikz/tikz_frames.pdf

BCsFiles := BCs/name_match.tex \
            BCs/farfield/farfield.tex \
            BCs/farfield/characteristic.tex \
            BCs/farfield/fixed.tex \
            BCs/farfield/freestream.tex \
            BCs/farfield/general.tex \
            BCs/inflow/inflow.tex \
            BCs/inflow/general.tex \
            BCs/inflow/sem.tex \
            BCs/inflow/subsonic.tex \
            BCs/inflow/supersonic.tex \
            BCs/numerical/numerical.tex \
            BCs/numerical/periodic/periodic.tex \
            BCs/numerical/periodic/common.tex \
            BCs/numerical/periodic/examples.tex \
            BCs/numerical/periodic/rotated.tex \
            BCs/numerical/periodic/translated.tex \
            BCs/numerical/symmetry.tex \
            BCs/outflow/outflow.tex \
            BCs/outflow/general.tex \
            BCs/outflow/subsonic.tex \
            BCs/outflow/supersonic.tex \
            BCs/walls/walls.tex \
            BCs/walls/adiabatic.tex \
            BCs/walls/isothermal.tex \
            BCs/walls/slip.tex

FilesFiles := files/grid.tex \
              files/grid_cgns.tex \
              files/grid_gambit.tex \
              files/grid_gmsh.tex \
              files/grid_plot3d.tex \
              files/interpolation_grid.tex \
              files/interrupt.tex \
              files/namelist.tex \
              files/probes.tex \
              files/residual.tex \
              files/restart.tex \
              files/solution_input.tex \
              files/solution_output.tex \
              files/stdout.tex \
              files/time_ave.tex

NMLFiles := nml/changes.tex \
            nml/nml_bc.tex \
            nml/nml_chemistry.tex \
            nml/nml_debug.tex \
            nml/nml_grid.tex \
            nml/nml_initialization.tex \
            nml/nml_interpolation.tex \
            nml/nml_limiter.tex \
            nml/nml_mms.tex \
            nml/nml_numerics.tex \
            nml/nml_output.tex \
            nml/nml_output_variables.tex \
            nml/nml_physics.tex \
            nml/nml_reference.tex \
            nml/nml_removed.tex \
            nml/nml_sem.tex \
            nml/nml_testcase_channel.tex \
            nml/nml_testcase_pulse.tex \
            nml/nml_testcase_sod.tex \
            nml/nml_testcase_tcf.tex \
            nml/nml_testcase_vortex.tex \
            nml/nml_time_ave.tex \
            nml/nml_timestep.tex \
            nml/nml_utility.tex \
            nml/nml_version.tex

NumericsFiles := numerics/aliasing.tex \
                 numerics/filtering.tex \
                 numerics/limiting.tex \
                 numerics/markers/markers.tex \
                 numerics/markers/ducros.tex \
                 numerics/markers/jump.tex \
                 numerics/markers/mlp.tex \
                 numerics/markers/resolution.tex \
                 numerics/markers/scaling.tex

PhysModelsFiles := sem/faq.tex \
                   sem/sem.tex \
                   sem/theory.tex \
                   sem/usage.tex
#sem/implementation.tex \

TestCasesFiles := test_cases/bump.tex \
                  test_cases/couette_inp.tex \
                  test_cases/couette.tex \
                  test_cases/density_pulse.tex \
                  test_cases/isentropic_vortex.tex \
                  test_cases/les_over_cylinder.tex \
                  test_cases/shu_osher.tex \
                  test_cases/sod.tex \
                  test_cases/taylor_couette.tex \
                  test_cases/taylor_green_vortex.tex

TimeSteppingFiles := time_stepping/algorithms/algorithms.tex \
                     time_stepping/algorithms/adaptive.tex \
                     time_stepping/algorithms/cfl_based.tex \
                     time_stepping/algorithms/constant_dt.tex \
                     time_stepping/controlling_dt.tex \
                     time_stepping/methods/methods.tex \
                     time_stepping/methods/recommended.tex \
                     time_stepping/methods/other.tex \
                     time_stepping/methods/bad.tex
