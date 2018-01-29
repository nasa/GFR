
ifndef VERBOSE
  QUIET := @
endif

include system_variables.mk

OBJDIR = objects/$(BUILDDIR)
MODDIR = modules/$(BUILDDIR)
FUNDIR = Functions

vpath %.o $(OBJDIR)
vpath %.f90 ./
vpath %.mod $(MODDIR)

MOD = $(MODFLAG) $(MODDIR)
FSFXFLG = .f90

CSFXFLG = .c
#CFLAG = -O3
#CFLAG += -g

HDF5LIB  := -lhdf5 -ldl -lz
CGNSLIB  := -lcgns
METISLIB := -lmetis

ifdef SZIPDIR
  SZIPLIB := -lsz
endif

ifdef MKLFPP
  MKLFLAG := -mkl=sequential
  MKLINC  := -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
  MKLLIB  := -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
endif
BLASLIB := $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a
LAPACKLIB := $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a

#ifneq "$(COMPFAMILY)" "pgi"
  ifdef PBSFPP
    PBSDIR := -L/u/scicon/tools/lib 
    PBSLIB := -lpbs_time_left
    ifneq "$(COMPFAMILY)" "intel"
      INTELVER := 2016.2.181
      INTELDIR := /nasa/intel/Compiler/$(INTELVER)
      PBSDIR += -L$(INTELDIR)/compilers_and_libraries/linux/lib/intel64_lin
      PBSLIB += -lirc
    endif
    TARGET = $(AUTOTARGETARCH)
  endif
#endif

TARGET ?= $(HOSTARCH)

FPP = $(FPPFLAG) $(METISFPP) $(PBSFPP) $(CGNSFPP) $(MPIFPP) $(MKLFPP) $(OTHERFPP)

#INCDIR = $(sort $(CGNSDIR))
LIBDIR = $(sort $(HDF5DIR) $(SZIPDIR) $(CGNSDIR) $(METISDIR))

FINC = $(foreach i,$(INCDIR),-I$(i)/include) -I./ $(MKLINC) $(OTHERINCDIR)
FLIB = $(foreach i,$(LIBDIR),-L$(i)/lib) $(PBSDIR) $(OTHERLIBDIR)
FLIB += $(MKLLIB) $(METISLIB) $(CGNSLIB) $(HDF5LIB) $(SZIPLIB) $(PBSLIB) $(OTHERLIB)

FCBASE = $(REALLOCFLAG) $(TRACEFLAG) -g

ADDFLAGS = $(MKLFLAG) $(MOD) $(FINC) $(FCADD) $(CMPFLAGS)

COBJS := $(OBJDIR)/mycpu.o

OBJS := $(OBJDIR)/kind_types.o \
        $(OBJDIR)/main.o \
        $(OBJDIR)/order_mod.o \
        $(OBJDIR)/quadrature_mod.o \
        $(OBJDIR)/geovar_mod.o \
        $(OBJDIR)/flowvar_mod.o \
        $(OBJDIR)/ovar_mod.o \
        $(OBJDIR)/gmsh_mod.o \
        $(OBJDIR)/plot3d_mod.o \
        $(OBJDIR)/cgnstypes_mod.o \
        $(OBJDIR)/cgns_mod.o \
        $(OBJDIR)/io.o \
        $(OBJDIR)/metrics.o \
        $(OBJDIR)/poly.o \
        $(OBJDIR)/initialize.o \
        $(OBJDIR)/generic.o \
        $(OBJDIR)/bc_mod.o \
        $(OBJDIR)/bnd_profiles_mod.o \
        $(OBJDIR)/triangle.o \
        $(OBJDIR)/flux.o \
        $(OBJDIR)/filter_mod.o \
        $(OBJDIR)/time_mod.o \
        $(OBJDIR)/typesgen_mod.o \
        $(OBJDIR)/eqn_idx_mod.o \
        $(OBJDIR)/correction_mod.o \
        $(OBJDIR)/derivatives_mod.o \
        $(OBJDIR)/interpolation_mod.o \
        $(OBJDIR)/projection_mod.o \
        $(OBJDIR)/vandermonde_mod.o \
        $(OBJDIR)/metis.o \
        $(OBJDIR)/parallel.o \
        $(OBJDIR)/mappings.o \
        $(OBJDIR)/connectivity.o \
        $(OBJDIR)/limiters.o \
        $(OBJDIR)/mms_mod.o \
        $(OBJDIR)/channel_mod.o \
        $(OBJDIR)/restart_mod.o \
        $(OBJDIR)/postproc_mod.o \
        $(OBJDIR)/pbs_mod.o \
        $(OBJDIR)/cpu_info_mod.o \
	$(OBJDIR)/averaging_mod.o


HighlightName := | /bin/grep -hi --color '[^ ]\+\.\(f[0-9]*\|c\)\>'

.PHONY: opt
opt: OPTIM = $(OPTFLAG) $(ALIGN32)
opt: optimize

.PHONY: opt-host
opt-host: OPTIM = $(OPTFLAG) $(HOSTARCH) $(ALIGN32)
opt-host: optimize

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
westmere: OPTIM = $(OPTFLAG) $(WESOPT)
westmere: optimize

.PHONY: merope
merope: OPTIM = $(OPTFLAG) $(WESOPT)
merope: optimize

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
deplist: FCBASE += -gen-dep=deplist.txt
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
  FPE     := -fpe0 -fpe-all=0 -fp-stack-check
  MODFLAG := -module
  ARCHFLAG := -x
  WESARCH := SSE4.2
  SANARCH := AVX
  IVYARCH := CORE-AVX-I
  HASARCH := CORE-AVX2
  BROARCH := CORE-AVX2
  SKYARCH := CORE-AVX512
  FPPFLAG := -fpp
  IPFLAG  := -ip
  IPOFLAG := -ipo
  ALIGN32 := -align array32byte
  ALIGN64 := -align array64byte
  OPTFLAG := -O3 $(IPFLAG)
  DBGFLAG := -check all -debug all -ftrapuv $(FPE)
  TRACEFLAG := -traceback
  REALLOCFLAG := -assume realloc_lhs
  CMPFLAGS := -DSPECIAL_FOR_INTEL
 #CMPFLAGS := -DSPECIAL_FOR_INTEL -DDISABLE_QP -mkl=sequential
  AUTOTARGETARCH := $(ARCHFLAG)$(SANARCH) -ax$(SKYARCH),$(HASARCH),$(IVYARCH) $(ALIGN32)
  HOSTARCH := $(ARCHFLAG)Host
else ifeq "$(COMPFAMILY)" "gcc"
  # GCC flags
  MODFLAG := -J
  ARCHFLAG := -march=
  ######################################################
  #### USE THESE ARCH FLAGS FOR GCC 4.9.4 OR NEWER
  WESARCH := nehalem
  SANARCH := sandybridge
  IVYARCH := ivybridge
  HASARCH := haswell
  BROARCH := broadwell
  SKYARCH := skylake-avx512 # ONLY WORKS FOR 6.1.0 or newer
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
  CMPFLAGS := -DDISABLE_DTIO -DSPECIAL_FOR_GCC
  AUTOTARGETARCH := $(ARCHFLAG)$(SANARCH)
  AUTOTARGETARCH += -mtune=$(IVYARCH) -mtune=$(HASARCH) -mtune=$(BROARCH) -mtune=$(SKYARCH)
  HOSTARCH := $(ARCHFLAG)native
else ifeq "$(COMPFAMILY)" "pgi"
  # PGI flags
  MODFLAG := -module
  ARCHFLAG := -tp=
  WESARCH := nehalem
  SANARCH := sandybridge
  IVYARCH := sandybridge
  HASARCH := haswell
  BROARCH := haswell
  SKYARCH := haswell
  FPPFLAG := -Mpreprocess
  OPTFLAG := -O3
 #OPTFLAG := -O3 -fast -Mipa=fast,inline
  DBGFLAG := -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -Minform=inform -Ktrap=fp
  TRACEFLAG := -traceback
  REALLOCFLAG := -Mallocatable=03
  CMPFLAGS := -DSPECIAL_FOR_PGI -DDISABLE_QP -DDISABLE_DTIO -Mbackslash
  AUTOTARGETARCH := -tp=$(SANARCH),$(HASARCH)
  HOSTARCH :=
endif

WESOPT := $(ARCHFLAG)$(WESARCH) $(ALIGN32)
SANOPT := $(ARCHFLAG)$(SANARCH) $(ALIGN32)
IVYOPT := $(ARCHFLAG)$(IVYARCH) $(ALIGN32)
HASOPT := $(ARCHFLAG)$(HASARCH) $(ALIGN32)
BROOPT := $(ARCHFLAG)$(BROARCH) $(ALIGN32)
SKYOPT := $(ARCHFLAG)$(SKYARCH) $(ALIGN64)

include dep.mk
