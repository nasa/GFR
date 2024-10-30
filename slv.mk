
ifndef VERBOSE
  QUIET := @
endif

include system_variables.mk

OBJDIR = objects/$(BUILDDIR)
MODDIR = modules/$(BUILDDIR)
FUNDIR = Functions
BCDIR = BC_Modules

vpath %.o $(OBJDIR)
vpath %.f90 ./
vpath %.mod $(MODDIR)

MOD = $(MODFLAG) $(MODDIR)
FSFXFLG = .f90

CSFXFLG = .c
#CFLAG = -O3
#CFLAG += -g

#HDF5LIB  := -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -ldl -lz
#HDF5LIB  := -lhdf5 -ldl -lz
HDF5LIB  := -lhdf5 -ldl
CGNSLIB  := -lcgns
METISLIB := -lmetis
ifdef TECDIR
  TECLIB  := -lteciompi -lstdc++
endif

ifdef SZIPDIR
  SZIPLIB := -lsz
endif

ifdef MKLFPP
  MKLFLAG := -qmkl=sequential
  MKLINC  := -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
  MKLLIB  := -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
endif
BLASLIB := $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a
LAPACKLIB := $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a

ifdef ENABLEOMP
  OMPADD = $(OMPFLAG) -DENABLE_OPENMP
endif

ifdef PBSFPP
 #PBSDIR := -L/u/scicon/tools/lib 
 #PBSLIB := -lpbs_time_left
 #ifneq "$(COMPFAMILY)" "intel"
 #  INTELVER := 2016.2.181
 #  INTELDIR := /nasa/intel/Compiler/$(INTELVER)
 # #PBSDIR += -L$(INTELDIR)/compilers_and_libraries/linux/lib/intel64_lin
 # #PBSLIB += -lirc
 #endif
  TARGET = $(AUTOTARGETARCH)
endif


TARGET ?= $(HOSTARCH)

FPP = $(FPPFLAG) $(METISFPP) $(PBSFPP) $(CGNSFPP) $(MPIFPP) $(MKLFPP) $(OTHERFPP)

ifeq "$(COMPFAMILY)" "gcc"
  INCDIR = $(sort $(CGNSDIR))
else
 #INCDIR = $(sort $(CGNSDIR))
endif
LIBDIR = $(sort $(HDF5DIR) $(SZIPDIR) $(CGNSDIR) $(METISDIR) $(TECDIR))

FINC = $(foreach i,$(INCDIR),-I$(i)/include) -I./ $(MKLINC) $(OTHERINCDIR)
FLIB = $(foreach i,$(LIBDIR),-L$(i)/lib) $(PBSDIR) $(OTHERLIBDIR)
FLIB += $(MKLLIB) $(METISLIB) $(CGNSLIB) $(HDF5LIB) $(SZIPLIB) $(PBSLIB) $(TECLIB) $(OTHERLIB)
ifeq "$(COMPFAMILY)" "intel"
  FLIB += $(LAPACKLIB) $(BLASLIB)
endif

FCBASE = $(REALLOCFLAG)
#FCBASE += $(TRACEFLAG)
#FCBASE += -g

ADDFLAGS = $(MKLFLAG) $(MOD) $(FINC) $(OMPADD) $(FCADD) $(CMPFLAGS)

COBJS := $(OBJDIR)/mycpu.o \
         $(OBJDIR)/pbs_time_left.o

OBJS := $(OBJDIR)/kind_types.o \
        $(OBJDIR)/main.o \
        $(OBJDIR)/random_mod.o \
        $(OBJDIR)/order_mod.o \
        $(OBJDIR)/element_type_mod.o \
        $(OBJDIR)/quadrature_mod.o \
        $(OBJDIR)/geomvar_mod.o \
        $(OBJDIR)/flowvar_mod.o \
        $(OBJDIR)/ovar_mod.o \
        $(OBJDIR)/gmsh_mod.o \
        $(OBJDIR)/plot3d_mod.o \
        $(OBJDIR)/cgnstypes_mod.o \
        $(OBJDIR)/cgns_mod.o \
        $(OBJDIR)/internal_grid_mod.o \
        $(OBJDIR)/io.o \
        $(OBJDIR)/write_cgns_mod.o \
        $(OBJDIR)/metrics.o \
        $(OBJDIR)/poly.o \
        $(OBJDIR)/initialize.o \
        $(OBJDIR)/generic.o \
        $(OBJDIR)/bc_mod.o \
        $(OBJDIR)/bc_flow_stats_mod.o \
        $(OBJDIR)/bnd_profiles_mod.o \
        $(OBJDIR)/triangle.o \
        $(OBJDIR)/tetrahedra.o \
        $(OBJDIR)/pyramid.o \
        $(OBJDIR)/prism.o \
        $(OBJDIR)/flux.o \
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
        $(OBJDIR)/metis.o \
        $(OBJDIR)/parallel.o \
        $(OBJDIR)/mappings.o \
        $(OBJDIR)/connectivity.o \
        $(OBJDIR)/limiters_mod.o \
        $(OBJDIR)/mlp_mod.o \
        $(OBJDIR)/mms_mod.o \
        $(OBJDIR)/channel_mod.o \
        $(OBJDIR)/restart_mod.o \
        $(OBJDIR)/postproc_mod.o \
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
        $(OBJDIR)/mrf_mod.o \
        $(OBJDIR)/timer_mod.o \
        $(OBJDIR)/turb_inflow_mod.o \
        $(OBJDIR)/thx_mod.o \
        $(OBJDIR)/taylor_couette_mod.o \
        $(OBJDIR)/shu_osher_mod.o \
        $(OBJDIR)/mydebug_mod.o \
        $(OBJDIR)/bc_wall_mod.o \
        $(OBJDIR)/bc_fixed_mod.o \
        $(OBJDIR)/bc_flowthrough_grad_mod.o \
        $(OBJDIR)/bc_gen_free_mod.o \
        $(OBJDIR)/bc_gen_in_mod.o \
        $(OBJDIR)/bc_gen_out_mod.o \
        $(OBJDIR)/bc_mms_mod.o \
        $(OBJDIR)/bc_periodic_mod.o \
        $(OBJDIR)/bc_relax_mod.o \
        $(OBJDIR)/bc_riemann_mod.o \
        $(OBJDIR)/bc_sub_in_mod.o \
        $(OBJDIR)/bc_sub_out_mod.o \
        $(OBJDIR)/bc_sup_out_mod.o \
        $(OBJDIR)/tecio_mod.o


HighlightName := | /bin/grep -hi --color '[^ ]\+\.\(f[0-9]*\|c\)\>'

.PHONY: opt
opt: OPTIM = $(OPTFLAG) $(TARGET)
opt: optimize

.PHONY: opt-host
opt-host: OPTIM = $(OPTFLAG) $(HOSTARCH)
opt-host: optimize

.PHONY: rome
rome: OPTIM = $(OPTFLAG) $(ROMOPT)
rome: optimize

.PHONY: skylake
skylake: OPTIM = $(OPTFLAG) $(SKYOPT)
skylake: optimize

.PHONY: broadwell
broadwell: OPTIM = $(OPTFLAG) $(BROOPT)
broadwell: optimize

.PHONY: haswell
haswell: OPTIM = $(OPTFLAG) $(HASOPT)
haswell: optimize

.PHONY: ivybridge
ivybridge: OPTIM = $(OPTFLAG) $(IVYOPT)
ivybridge: optimize

.PHONY: sandybridge
sandybridge: OPTIM = $(OPTFLAG) $(SANOPT)
sandybridge: optimize

.PHONY: westmere
ifeq "$(COMPFAMILY)" "pgi"
#westmere: ADDFLAGS += -pg
endif
westmere: OPTIM = $(OPTFLAG) $(WESOPT)
#westmere: OPTIM = -O3 $(OPTFLAG) $(ARCHFLAG)$(WESARCH)
westmere: optimize

.PHONY: merope
merope: OPTIM = $(OPTFLAG) $(WESOPT)
#merope: OPTIM = $(OPTFLAG) $(ARCHFLAG)$(WESARCH)
merope: optimize

.PHONY: gpu-mm
gpu-mm: OPTIM = $(OPTFLAG) $(ARCHFLAG)$(WESARCH)
gpu-mm: optimize

.PHONY: gpu
gpu: OPTIM = $(OPTFLAG) $(ARCHFLAG)$(WESARCH)
gpu: optimize

.PHONY: optimize
optimize: FCFLAG = $(FCBASE) $(OPTIM) $(ADDFLAGS)
optimize: CFLAG = -O3
optimize: all

.PHONY: debug
debug: FPP += -DDEBUG_ON
#debug: FPP += -DDEBUG_ON -DDISABLE_QP
debug: FCFLAG = $(FCBASE) $(DBGFLAG) $(ADDFLAGS)
debug: CFLAG = -O0 -g
debug: all

.PHONY: deplist
deplist: FCBASE += -gen-dep=deplist.mk
deplist: debug

.PHONY: all
all: Command-Beg := echo
all: Command-End  = $(HighlightName); $(FC) $(FPP) $(FCFLAG) -c -o $@ $<
all: Command-Cend  = $(HighlightName); $(CC) $(CFLAG) -c -o $@ $<
all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS) $(COBJS)
	$(FC) $(FPP) $(FCFLAG) -o $(EXECUTABLE) $(OBJS) $(COBJS) $(EXTOBJS) $(FLIB)
	if [ -d $(EXEDIR) ] ; then /bin/cp -f $(EXECUTABLE) $(EXEDIR)/$(EXECUTABLE) ; fi

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

.PHONY: link-to-haswell
link-to-haswell:
	if [ -d $(EXEDIR) ] ; then /bin/ln -sf gfr_haswell $(EXEDIR)/gfr_broadwell ; fi

.PHONY: clean
clean:
	/bin/rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod $(EXECUTABLE)

.PHONY: dump-variables
dump-variables:
	$(QUIET) echo "FC = $(FC)"
	$(QUIET) echo "FPP = $(FPP)"
	$(QUIET) echo "FCBASE = $(FCBASE)"
	$(QUIET) echo "ADDFLAGS = $(ADDFLAGS)"
	$(QUIET) echo "FCFLAG = $(FCFLAG)"
	$(QUIET) echo "OBJDIR = $(OBJDIR)"
	$(QUIET) echo "MODDIR = $(MODDIR)"
	$(QUIET) echo "EXECUTABLE = $(EXECUTABLE)"

# Unique definitions for each compiler
ifeq "$(COMPFAMILY)" "intel"
  # Intel flags
 #FPSTACK := -fp-stack-check
  FPE     := -fpe0 -fpe-all=0 -fp-stack-check
 #FPEALL  := -fpe-all=0
  MODFLAG := -module
  OMPFLAG := -qopenmp
  ARCHFLAG := -x
  ARCHFLAGROM := -march=
  WESARCH := SSE4.2
  SANARCH := AVX
  IVYARCH := CORE-AVX-I
  HASARCH := CORE-AVX2
  BROARCH := CORE-AVX2
 #SKYARCH := CORE-AVX512
  SKYARCH := SKYLAKE
 #SKYARCH := SKYLAKE-AVX512
 #SKYOPT  := $(ARCHFLAG)$(SKYARCH) $(ALIGN64) $(FPE)
 #SKYOPT  := $(ARCHFLAG)$(SKYARCH) $(ALIGN64) $(FPE) -qopt-zmm-usage=high
  ROMARCH := core-avx2 -fma
  FPPFLAG := -fpp
  IPFLAG  := -ip
  IPOFLAG := -ipo
  ALIGN32 := -align array32byte
  ALIGN64 := -align array64byte
  OPTFLAG := -traceback -g -O3 $(IPFLAG)
  DBGFLAG := -check all -debug all -ftrapuv $(FPE) -g -traceback
  TRACEFLAG := -traceback
  REALLOCFLAG := -assume realloc_lhs
 #CMPFLAGS := -DSPECIAL_FOR_INTEL
  CMPFLAGS := -DSPECIAL_FOR_INTEL -DDISABLE_QP -qmkl=sequential -DOCTID
 #CPMFLAGS += -guide -parallel -qopt-report=5 -qopt-report-file=stderr
  AUTOTARGETARCH := $(ARCHFLAG)$(SANARCH) -ax$(SKYARCH),$(HASARCH),$(IVYARCH) $(ALIGN32)
  HOSTARCH := $(ARCHFLAG)Host
else ifeq "$(COMPFAMILY)" "gcc"
  # GCC flags
  MODFLAG := -J
  OMPFLAG := -fopenmp
  ARCHFLAG := -march=
  ARCHFLAGROM := -march=
  ######################################################
  #### USE THESE ARCH FLAGS FOR GCC 4.9.4 OR NEWER
  WESARCH := westmere
  SANARCH := sandybridge
  IVYARCH := ivybridge
  HASARCH := haswell
  BROARCH := broadwell
  SKYARCH := skylake-avx512 # ONLY WORKS FOR 6.1.0 or newer
  ROMARCH := znver2
  #### USE THESE ARCH FLAGS FOR GCC 4.8.5 OR OLDER
 #WESARCH := corei7
 #SANARCH := corei7-avx
 #IVYARCH := corei7-avx-i
 #HASARCH := corei7-avx2
 #BROARCH := corei7-avx2
 #SKYARCH := corei7-avx2
  #######################################################
  FPPFLAG := -cpp
  OPTFLAG := -O3
  DBGFLAG := -fcheck=all -Wall -Wextra -Wuninitialized
  TRACEFLAG := -fbacktrace
  REALLOCFLAG := -frealloc-lhs
  CMPFLAGS := -DDISABLE_DTIO -DSPECIAL_FOR_GCC -DDISABLE_QP -DOCTID
  AUTOTARGETARCH := $(ARCHFLAG)$(SANARCH)
  AUTOTARGETARCH += -mtune=$(IVYARCH) -mtune=$(HASARCH) -mtune=$(BROARCH) -mtune=$(SKYARCH) -mtune=$(ROMARCH)
  HOSTARCH := $(ARCHFLAG)native
else ifeq "$(COMPFAMILY)" "pgi"
  # PGI flags
  MODFLAG := -module
  OMPFLAG := -mp
  ARCHFLAG := -tp=
  ARCHFLAGROM := -tp=
  WESARCH := nehalem
  SANARCH := sandybridge
  IVYARCH := sandybridge
  HASARCH := haswell
  BROARCH := haswell
  SKYARCH := haswell
  ROMARCH := zen
  FPPFLAG := -Mpreprocess
  OPTFLAG := -O3
  DBGFLAG := -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -Minform=inform -Ktrap=fp
  TRACEFLAG := #-traceback
  REALLOCFLAG := -Mallocatable=03
  CMPFLAGS := -DSPECIAL_FOR_PGI -DDISABLE_QP -DDISABLE_DTIO -Mbackslash
  AUTOTARGETARCH := -tp=$(SANARCH),$(HASARCH)
  HOSTARCH :=
endif

WESOPT  := $(ARCHFLAG)$(WESARCH) $(ALIGN32)
SANOPT  := $(ARCHFLAG)$(SANARCH) $(ALIGN32)
IVYOPT  := $(ARCHFLAG)$(IVYARCH) $(ALIGN32)
HASOPT  := $(ARCHFLAG)$(HASARCH) $(ALIGN32)
BROOPT  := $(ARCHFLAG)$(WESARCH) -ax$(BROARCH) -fma $(ALIGN32)
SKYOPT  := $(ARCHFLAG)$(SKYARCH) $(ALIGN64) -qopt-zmm-usage=high
ROMOPT  := $(ARCHFLAGROM)$(ROMARCH)

include dep.mk
