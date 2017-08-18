
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
CFLAG = -O3 -g

HDF5LIB  := -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -ldl -lz
CGNSLIB  := -lcgns
METISLIB := -lmetis

ifdef SZIPDIR
  SZIPLIB := -lsz
endif

ifdef MKLFPP
  MKLFLAG := -mkl
  MKLINC  := -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
  MKLLIB  := -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
endif

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

TARGET ?= $(HOSTARCH)

FPP = $(FPPFLAG) $(METISFPP) $(PBSFPP) $(CGNSFPP) $(MPIFPP) $(MKLFPP) $(OTHERFPP)

INCDIR = $(sort $(CGNSDIR) $(METISDIR))
LIBDIR = $(sort $(HDF5DIR) $(SZIPDIR) $(CGNSDIR) $(METISDIR))

FINC = $(foreach i,$(INCDIR),-I$(i)/include) -I./ $(MKLINC) $(OTHERINCDIR)
FLIB = $(foreach i,$(LIBDIR),-L$(i)/lib) $(PBSDIR) $(OTHERLIBDIR)
FLIB += $(MKLLIB) $(METISLIB) $(CGNSLIB) $(HDF5LIB) $(SZIPLIB) $(PBSLIB) $(OTHERLIB)

FCBASE = $(TRACEFLAG) -g $(REALLOCFLAG)

ADDFLAGS = $(MKLFLAG) $(MOD) $(FINC) $(FCADD) $(CMPFLAGS)

COBJS := $(OBJDIR)/mycpu.o

OBJS := $(OBJDIR)/kind_types.o \
        $(OBJDIR)/Euler_tri.o \
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

.PHONY: O2-ip-opt
O2-ip-opt: OPTIM = -O2 $(OPTFLAG) $(TARGET)
O2-ip-opt: optimize

.PHONY: O3-ipo-opt
O3-ipo-opt: OPTIM = -O3 $(OPTFLAG) $(TARGET)
O3-ipo-opt: optimize

.PHONY: broadwell
broadwell: OPTIM = -O3 $(OPTFLAG) $(ARCHFLAG)$(BROARCH)
broadwell: optimize

.PHONY: haswell
haswell: OPTIM = -O3 $(OPTFLAG) $(ARCHFLAG)$(HASARCH)
haswell: optimize

.PHONY: ivybridge
ivybridge: OPTIM = -O3 $(OPTFLAG) $(ARCHFLAG)$(IVYARCH)
ivybridge: optimize

.PHONY: sandybridge
sandybridge: OPTIM = -O3 $(OPTFLAG) $(ARCHFLAG)$(SANARCH)
sandybridge: optimize

.PHONY: westmere
westmere: OPTIM = -O3 $(OPTFLAG) $(ARCHFLAG)$(WESARCH)
westmere: optimize

.PHONY: merope
merope: OPTIM = -O3 $(OPTFLAG) $(ARCHFLAG)$(WESARCH)
merope: optimize

.PHONY: optimize
optimize: FCFLAG = $(FCBASE) $(OPTIM) $(ADDFLAGS)
optimize: all

.PHONY: debug
debug: FPP += -DDEBUG_ON
debug: FCFLAG = $(FCBASE) $(DBGFLAG) $(ADDFLAGS)
debug: all


.PHONY: deplist
deplist: FCBASE += -gen-dep=deplist.txt
deplist: debug

.PHONY: time
time: Command-Beg := /usr/bin/time -f '%e'
time: Command-End  = 2>&1 | awk '{ printf "%-25s%8.2f seconds", "$<", $$1 }' $(HighlightName)
time: FCFLAG = $(FCBASE) -O3 $(OPTFLAG) $(TARGET) $(ADDFLAGS)
time: $(OBJS)
	$(FC) $(FPP) $(FCFLAG) -o $(EXECUTABLE) $(OBJS) $(FLIB)
	/bin/cp -f $(EXECUTABLE) $(EXEDIR)/$(EXECUTABLE)

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
##$(QUIET) $(Time-Beg) $(FC) $(FPP) $(FCFLAG) -c -o $@ $< $(Time-End) $(HighlightName)

.PHONY: link-to-haswell
link-to-haswell:
	if [ -d $(EXEDIR) ] ; then /bin/ln -sf recon_haswell $(EXEDIR)/recon_broadwell ; fi

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
  MODFLAG := -module
  ARCHFLAG := -x
  TUNEFLAG := -ax
  WESARCH := SSE4.2
  SANARCH := AVX
  IVYARCH := CORE-AVX-I
  HASARCH := CORE-AVX2
  BROARCH := CORE-AVX2
  FPPFLAG := -fpp
  OPTFLAG := -ip -align array32byte
  DBGFLAG := -check all -debug all -ftrapuv -fp-stack-check
  IPFLAG  := -ip
  IPOFLAG := -ipo
  TRACEFLAG := -traceback
  REALLOCFLAG := -assume realloc_lhs
  CMPFLAGS := -DSPECIAL_FOR_INTEL
  AUTOTARGETARCH := $(ARCHFLAG)$(SANARCH) -ax$(HASARCH),$(IVYARCH)
  HOSTARCH := $(ARCHFLAG)Host
else ifeq "$(COMPFAMILY)" "gcc"
  # GCC flags
  MODFLAG := -J
  ARCHFLAG := -march=
  TUNEFLAG := -mtune=
  WESARCH := westmere
  SANARCH := sandybridge
  IVYARCH := ivybridge
  HASARCH := haswell
  BROARCH := broadwell
  FPPFLAG := -cpp
  OPTFLAG := 
  DBGFLAG := -fcheck=all -Wall -Wextra -Wuninitialized
  TRACEFLAG := -fbacktrace
  REALLOCFLAG := -frealloc-lhs
  CMPFLAGS := -DDISABLE_DTIO -DSPECIAL_FOR_GCC
  EXTOBJS  := $(METISDIR)/include/metis.o
  AUTOTARGETARCH := $(ARCHFLAG)$(SANARCH)
  AUTOTARGETARCH += -mtune=$(IVYARCH) -mtune=$(HASARCH) -mtune=$(BROARCH)
  HOSTARCH := $(ARCHFLAG)native
else ifeq "$(COMPFAMILY)" "pgi"
  # PGI flags
  MODFLAG := -module
  ARCHFLAG := -tp=
  TUNEFLAG := -tp=
  WESARCH := nehalem
  SANARCH := sandybridge
  IVYARCH := sandybridge
  HASARCH := haswell
  BROARCH := haswell
  FPPFLAG := -Mpreprocess
  OPTFLAG :=
  DBGFLAG := -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -Minform=inform -Ktrap=fp
  TRACEFLAG := -traceback
  REALLOCFLAG := -Mallocatable=03
  CMPFLAGS := -DSPECIAL_FOR_PGI -DDISABLE_QP -DDISABLE_DTIO -Mbackslash
  EXTOBJS  := $(METISDIR)/include/metis.o
  AUTOTARGETARCH := -tp=$(SANARCH),$(HASARCH)
  HOSTARCH :=
endif

include dep.mk
