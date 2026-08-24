
ifndef VERBOSE
  QUIET := @
endif

include system_variables.mk

OBJDIR = objects/$(BUILDDIR)
MODDIR = modules/$(BUILDDIR)
FUNDIR = Functions
BCDIR  = BC_Modules
INTDIR = SolInterp_Modules
NMLDIR = Nml_Modules
TSTDIR = Testcase_Modules

vpath %.o $(OBJDIR)
vpath %.f90 ./
vpath %.mod $(MODDIR)

ifeq "$(FC)" "ftn"
  COMPILER = $(shell $(FC) -craype-verbose 2>/dev/null | awk '{print $$1}')
endif

MOD = $(MODFLAG) $(MODDIR)
FSFXFLG = .f90
CSFXFLG = .c

ifneq "$(FC)" "ftn"
  HDF5LIB  := -lhdf5 -ldl
endif
CGNSLIB  := -lcgns
ifeq "$(COMPFAMILY)" "aocc"
  CGNSMOD := -I $(CGNSDIR)/include
endif
METISLIB := -lmetis
ifeq "$(PARMETISFPP)" "-DENABLE_PARMETIS"
  METISLIB := -lparmetis $(METISLIB)
endif
ifdef TECDIR
  TECLIB  := -lteciompi -lstdc++
endif

ifeq "$(COMPFAMILY)" "intel"
  ifdef MKLDIR
    ifdef ENABLE_MKL_BLAS
      MKLFPP += -DUSE_INTEL_MKL_BLAS
      MKLLIB += -lmkl_blas95_lp64
    endif
    ifdef ENABLE_MKL_LAPACK
      MKLFPP += -DUSE_INTEL_MKL_LAPACK
      MKLLIB += -lmkl_lapack95_lp64
    endif
    ifdef MKLFPP
      MKLFLAG := -qmkl=sequential
    endif
  endif
endif

ifdef ENABLEOMP
  OMPADD = $(OMPFLAG) -DENABLE_OPENMP
endif

ifdef HPCFPP
  ifeq "$(FC)" "ftn"
    HPCFPP += -DATHENA_PBS
  endif
endif

FPP = $(strip $(FPPFLAG) $(HPCFPP) $(MPIFPP) $(OTHERFPP))

ifeq "$(COMPFAMILY)" "gnu"
  INCDIR = $(sort $(CGNSDIR))
else
 #INCDIR = $(sort $(CGNSDIR))
endif
LIBDIR = $(sort $(HDF5DIR) $(CGNSDIR) $(METISDIR) $(PARMETISDIR) $(TECDIR))

FINC = $(strip $(foreach i,$(INCDIR),-I$(i)/include) -I./ $(OTHERINCDIR))
FLIB = $(strip $(foreach i,$(LIBDIR),-L$(i)/lib) $(OTHERLIBDIR))
FLIB += $(strip $(MKLLIB) $(METISLIB) $(CGNSLIB) $(HDF5LIB) $(TECLIB) $(OTHERLIB))

FCBASE = $(strip $(REALLOCFLAG) $(TRACEFLAG))

ADDFLAGS = $(strip $(OTHERFLAGS) $(MOD) $(FINC) $(OMPADD) $(FCADD) $(CMPFLAGS))

COBJS := $(OBJDIR)/mycpu.o \
         $(OBJDIR)/pbs_time_left.o

OBJS := $(OBJDIR)/kind_types_mod.o \
        $(OBJDIR)/main.o \
        $(OBJDIR)/random_mod.o \
        $(OBJDIR)/openmp_mod.o \
        $(OBJDIR)/order_mod.o \
        $(OBJDIR)/element_type_mod.o \
        $(OBJDIR)/quadrature_mod.o \
        $(OBJDIR)/geomvar_mod.o \
        $(OBJDIR)/flowvar_mod.o \
        $(OBJDIR)/ovar_mod.o \
        $(OBJDIR)/bc_input_0001_mod.o \
        $(OBJDIR)/nml_base_mod.o \
        $(OBJDIR)/nml_0001_mod.o \
        $(OBJDIR)/gmsh_mod.o \
        $(OBJDIR)/plot3d_mod.o \
        $(OBJDIR)/cgnstypes_mod.o \
        $(OBJDIR)/cgns_mod.o \
        $(OBJDIR)/internal_grid_mod.o \
        $(OBJDIR)/io_mod.o \
        $(OBJDIR)/write_cgns_mod.o \
        $(OBJDIR)/metrics_mod.o \
        $(OBJDIR)/poly_mod.o \
        $(OBJDIR)/initialization_mod.o \
        $(OBJDIR)/generic_mod.o \
        $(OBJDIR)/bc_mod.o \
        $(OBJDIR)/bc_flow_stats_mod.o \
        $(OBJDIR)/bnd_profiles_mod.o \
        $(OBJDIR)/triangle_mod.o \
        $(OBJDIR)/tetrahedra_mod.o \
        $(OBJDIR)/pyramid_mod.o \
        $(OBJDIR)/prism_mod.o \
        $(OBJDIR)/flux_mod.o \
        $(OBJDIR)/filter_mod.o \
        $(OBJDIR)/dt_mod.o \
        $(OBJDIR)/erk_mod.o \
        $(OBJDIR)/rhs_mod.o \
        $(OBJDIR)/oerk_mod.o \
        $(OBJDIR)/eperk_mod.o \
        $(OBJDIR)/matrix_types_mod.o \
        $(OBJDIR)/eqn_idx_mod.o \
        $(OBJDIR)/correction_mod.o \
        $(OBJDIR)/derivatives_mod.o \
        $(OBJDIR)/interpolation_mod.o \
        $(OBJDIR)/cell_ave_mod.o \
        $(OBJDIR)/projection_mod.o \
        $(OBJDIR)/vandermonde_mod.o \
        $(OBJDIR)/periodic_mod.o \
        $(OBJDIR)/metis_mod.o \
        $(OBJDIR)/parmetis_mod.o \
        $(OBJDIR)/parallel_mod.o \
        $(OBJDIR)/mappings_mod.o \
        $(OBJDIR)/connectivity_mod.o \
        $(OBJDIR)/limiters_mod.o \
        $(OBJDIR)/mlp_mod.o \
        $(OBJDIR)/channel_mod.o \
        $(OBJDIR)/restart_mod.o \
        $(OBJDIR)/postproc_mod.o \
        $(OBJDIR)/probes_mod.o \
        $(OBJDIR)/pbs_mod.o \
        $(OBJDIR)/cpu_info_mod.o \
        $(OBJDIR)/averaging_mod.o \
        $(OBJDIR)/gpu_mod.o \
        $(OBJDIR)/invs_flux_mod.o \
        $(OBJDIR)/visc_flux_mod.o \
        $(OBJDIR)/upwind_flux_mod.o \
        $(OBJDIR)/plot3dio_mod.o \
        $(OBJDIR)/octstat_mod.o \
        $(OBJDIR)/octroot_mod.o \
        $(OBJDIR)/octedat_mod.o \
        $(OBJDIR)/octnode_mod.o \
        $(OBJDIR)/octtree_mod.o \
        $(OBJDIR)/octedat_facexyz_mod.o \
        $(OBJDIR)/interpsolvar_mod.o \
        $(OBJDIR)/interpsol_mod.o \
        $(OBJDIR)/thermofix_mod.o \
        $(OBJDIR)/thermo_mod.o \
        $(OBJDIR)/transport_mod.o \
        $(OBJDIR)/reactions_mod.o \
        $(OBJDIR)/ignition_mod.o \
        $(OBJDIR)/mrf_mod.o \
        $(OBJDIR)/timer_mod.o \
        $(OBJDIR)/spatial_ave_mod.o \
        $(OBJDIR)/synthetic_eddy_mod.o \
        $(OBJDIR)/thx_mod.o \
        $(OBJDIR)/taylor_couette_mod.o \
        $(OBJDIR)/advection_js_mod.o \
        $(OBJDIR)/sod_shock_tube_mod.o \
        $(OBJDIR)/shu_osher_mod.o \
        $(OBJDIR)/interpsol_query_mod.o \
        $(OBJDIR)/bc_wall_mod.o \
        $(OBJDIR)/bc_fixed_mod.o \
        $(OBJDIR)/bc_flowthrough_grad_mod.o \
        $(OBJDIR)/bc_gen_free_mod.o \
        $(OBJDIR)/bc_gen_in_mod.o \
        $(OBJDIR)/bc_gen_out_mod.o \
        $(OBJDIR)/bc_periodic_mod.o \
        $(OBJDIR)/bc_relax_mod.o \
        $(OBJDIR)/bc_riemann_mod.o \
        $(OBJDIR)/bc_sub_in_mod.o \
        $(OBJDIR)/bc_sub_out_mod.o \
        $(OBJDIR)/bc_sup_out_mod.o \
        $(OBJDIR)/tecio_mod.o


HighlightName := | /bin/grep -hi --color '[^ ]\+\.\(f[0-9]*\|c\)\>'

.PHONY: opt
opt: OPTIM = $(HOSTOPT)
opt: optimize

.PHONY: aitken
aitken: OPTIM = $(AITKENOPT)
aitken: optimize

.PHONY: athena
athena: OPTIM = $(ATHENAOPT)
athena: optimize

.PHONY: gpu
gpu: OPTIM = $(AITKENOPT)
gpu: optimize

.PHONY: optimize
optimize: FCFLAG = $(strip $(FCBASE) $(OPTIM) $(ADDFLAGS))
optimize: CFLAG = -O3
optimize: all

.PHONY: debug
debug: FCFLAG = $(strip $(FCBASE) $(DBGFLAG) $(ADDFLAGS))
debug: CFLAG = -O0 -g
debug: all

.PHONY: all
all: Command-Beg := echo
all: Command-End  = $(HighlightName); $(FC) $(FPP) $(FCFLAG) -c -o $@ $<
all: Command-Cend  = $(HighlightName); $(CC) $(CFLAG) -c -o $@ $<
all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS) $(COBJS)
	$(FC) $(FPP) $(FCFLAG) $(MKLFLAG) -o $(EXECUTABLE) $(OBJS) $(COBJS) $(EXTOBJS) $(FLIB)
	if [ -d $(EXEDIR) ] ; then $(CP) $(EXECUTABLE) $(EXEDIR)/$(EXECUTABLE) ; fi

%.o:%$(FSFXFLG)
	$(QUIET) echo $(FC) $(FPP) $(FCFLAG) -c -o $@ $< $(HighlightName);\
	$(FC) $(FPP) $(FCFLAG) -c -o $@ $<

%.o:%$(CSFXFLG)
	$(QUIET) echo $(CC) $(CFLAG) -c -o $@ $< $(HighlightName);\
	$(CC) $(CFLAG) -c -o $@ $<

$(COBJS):
	$(QUIET) $(Command-Beg) $(CC) $(CFLAG) -c -o $@ $< $(Command-Cend)

$(OBJS):
	$(QUIET) $(Command-Beg) $(FC) $(FPP) $(FCFLAG) -c -o $@ $< $(Command-End)

.PHONY: clean
clean:
	$(RM) $(OBJDIR)/*.o $(MODDIR)/*.mod $(EXECUTABLE)
	$(QUIET) if [ -f $(EXEDIR)/$(EXECUTABLE) ] ; then \
		   $(RM) -v $(EXEDIR)/$(EXECUTABLE) ; fi

.PHONY: dump-variables
dump-variables:
	$(QUIET) echo "FC = $(FC)"
	$(QUIET) echo "FPP = $(FPP)"
	$(QUIET) echo "FCBASE = $(FCBASE)"
	$(QUIET) echo "ADDFLAGS = $(ADDFLAGS)"
	$(QUIET) echo "DBGFLAG = $(DBGFLAG)"
	$(QUIET) echo "HOSTOPT = $(HOSTOPT)"
	$(QUIET) echo "AITKENOPT = $(AITKENOPT)"
	$(QUIET) echo "ATHENAOPT = $(ATHENAOPT)"
	$(QUIET) echo "OBJDIR = $(OBJDIR)"
	$(QUIET) echo "MODDIR = $(MODDIR)"
	$(QUIET) echo "EXECUTABLE = $(EXECUTABLE)"
	$(QUIET) echo "COMPFAMILY = '$(COMPFAMILY)'"
	$(QUIET) echo "COMPILER = '$(COMPILER)'"
	$(QUIET) echo "FINC = '$(FINC)'"
	$(QUIET) echo "LIBDIR = '$(LIBDIR)'"
	$(QUIET) echo "FLIB = '$(FLIB)'"

################################################################################
# Unique definitions for each compiler
################################################################################
ifeq "$(COMPFAMILY)" "intel"
  # Intel classic compiler (ifort) flags
  #------------------------------------------------------------
  FPPFLAG     := -fpp
  MODFLAG     := -module
  OMPFLAG     := -qopenmp
  REALLOCFLAG := -assume realloc_lhs -assume fpe_summary
  TRACEFLAG   := -traceback -g
  CMPFLAGS    := -DSPECIAL_FOR_INTEL
  #------------------------------------------------------------
  DBGFLAG := -check all -debug all -ftrapuv \
	     -fpe0 -fpe-all=0 -fp-stack-check \
	     -diag-disable=10182,10448
  #------------------------------------------------------------
  OPTFLAGS  := -O3 -ip -fma
  AITKENOPT := $(OPTFLAGS) -march=core-avx2
  ATHENAOPT := $(AITKENOPT) -diag-disable=10121,10448 \
	       -mavx512f -mavx512dq -mavx512bw -mavx512vl
  HOSTOPT   := $(OPTFLAGS) -xHost
#===============================================================================
else ifeq "$(COMPFAMILY)" "oneapi"
  # Intel OneAPI compiler (ifx) flags
  #------------------------------------------------------------
  FPPFLAG     := -fpp
  MODFLAG     := -module
  OMPFLAG     := -qopenmp
  REALLOCFLAG := -assume realloc_lhs -assume fpe_summary
  TRACEFLAG   := -traceback -g
  CMPFLAGS    := -DSPECIAL_FOR_INTEL
  #------------------------------------------------------------
  DBGFLAG := -check all,nouninit -debug all -ftrapuv \
	     -fpe0 -fpe-all=0 -diag-disable=10182
  #------------------------------------------------------------
  OPTFLAGS  := -O3 -fma
  AITKENOPT := $(OPTFLAGS) -march=core-avx2
  ATHENAOPT := $(AITKENOPT) -diag-disable=10121 \
	       -mavx512f -mavx512dq -mavx512bw -mavx512vl
  HOSTOPT   := $(OPTFLAGS) -xHost
#===============================================================================
else ifeq "$(COMPFAMILY)" "aocc"
  # AOCC compiler (flang) flags
  #------------------------------------------------------------
  FPPFLAG     := -Mpreprocess
  MODFLAG     := -module
  OMPFLAG     := -mp
  REALLOCFLAG := -Mallocatable=03
  TRACEFLAG   := -g
 #TRACEFLAG   := -fbacktrace -g
  CMPFLAGS    := -DSPECIAL_FOR_AOCC -DDISABLE_DTIO -Mnobackslash
## NOTE: The -Mnobackslash flag is equal to pgi/nvidia -Mbackslash
##       -Mnobackslash causes \ in strings to be treated like a regular character
##       -Mbackslash   causes \ in strings to be treated like a C-style escape character
###------------------------------------------------------------
  DBGFLAG := -Wall -Wextra -Wuninitialized
 #DBGFLAG := -fcheck=all -Wall -Wextra -Wuninitialized
###------------------------------------------------------------
  OPTFLAGS  :=-O3 -favoid-fpe-causing-opt
  BROOPT    := broadwell
  SKYOPT    := skylake
  CASOPT    := cascadelake
  ROMOPT    := znver2
  MILOPT    := znver3
  AITKENOPT := $(OPTFLAGS) -march=$(ROMOPT)
  ATHENAOPT := $(OPTFLAGS) -march=znver5
  HOSTOPT   := $(OPTFLAGS) -march=native -mtune=znver2
#===============================================================================
else ifeq "$(COMPFAMILY)" "cray"
  # Cray compiler (crayftn) flags
  #------------------------------------------------------------
  FPPFLAG     := -eT
  MODFLAG     := -ef -J
  OMPFLAG     := -fopenmp
  REALLOCFLAG := -ew
  TRACEFLAG   := -g
  CMPFLAGS    := -DSPECIAL_FOR_CRAY -DDISABLE_DTIO -hnoacc -xacc -M 7405,7206,7409,7418,296
###------------------------------------------------------------
  DBGFLAG := -O0 -m3 -R abcdps -Ktrap=fp
 #DBGFLAG += -fsanitize=address
###------------------------------------------------------------
  OPTMESGS  := -hlist=dmioe
  OPTFLAGS  := -O3 -Oipa5 -hfma #$(OPTMESGS)
  BROOPT    := broadwell
  SKYOPT    := x86-skylake
  CASOPT    := x86-cascadelake
  ROMOPT    := x86-rome
  MILOPT    := x86-milan
  AITKENOPT := $(OPTFLAGS) -hcpu=$(ROMOPT)
  ATHENAOPT := $(OPTFLAGS) -hcpu=x86-turin
  HOSTOPT   := $(OPTFLAGS)
#===============================================================================
else ifeq "$(COMPFAMILY)" "gnu"
  # GNU compiler (gfortran) flags
  #------------------------------------------------------------
  FPPFLAG     := -cpp
  MODFLAG     := -J
  OMPFLAG     := -fopenmp
  REALLOCFLAG := -frealloc-lhs
  TRACEFLAG   := -fbacktrace -g
  CMPFLAGS    := -DSPECIAL_FOR_GNU -DDISABLE_DTIO
  #------------------------------------------------------------
  DBGFLAG := -fcheck=all -Wall -Wextra -Wuninitialized
  #------------------------------------------------------------
  OPTFLAGS  := -O3
  AITKENOPT := $(OPTFLAGS) -march=znver2
  ATHENAOPT := $(OPTFLAGS) -march=znver5
  HOSTOPT   := $(OPTFLAGS) -march=native -mtune=znver2
#===============================================================================
else ifeq "$(COMPFAMILY)" "nvidia"
  # NVIDIA compiler flags
  #------------------------------------------------------------
  FPPFLAG     := -Mpreprocess
  MODFLAG     := -module
  OMPFLAG     := -mp
  REALLOCFLAG := -Mallocatable=03
  TRACEFLAG   := -traceback -g
  CMPFLAGS    := -DSPECIAL_FOR_NVIDIA -DDISABLE_DTIO -Mbackslash
## NOTE: The -Mbackslash flag is equal to AOCC -Mnobackslash
##       -Mbackslash   causes \ in strings to be treated like a regular character
##       -Mnobackslash causes \ in strings to be treated like a C-style escape character
###------------------------------------------------------------
  DBGFLAG := -Mbounds -Mchkptr -Mchkstk -Minform=inform -Ktrap=fp
###------------------------------------------------------------
  OPTFLAGS  :=-O3 -Mfma
  BROOPT    := broadwell
  SKYOPT    := skylake
  CASOPT    := cascadelake
  ROMOPT    := znver2
  MILOPT    := znver3
  AITKENOPT := $(OPTFLAGS) -tp=$(ROMOPT)
  ATHENAOPT := $(OPTFLAGS) -tp=znver4
  HOSTOPT   := $(OPTFLAGS) -tp=native
endif
#===============================================================================

##############################################################
###  Variable definitions for some common system commands  ###
##############################################################

RM := /bin/rm -f
CP := /bin/cp -fp

################################################################################

include dep.mk
