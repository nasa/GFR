
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

MOD = -module $(MODDIR)
FSFXFLG = .f90

CSFXFLG = .c
CFLAG = -O3 -g

HDF5LIB = -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -ldl

O2OPT = -O2 -ip
#O3OPT = -O3 -ipo
 O3OPT = -O3 -ip

TARGET = $(AUTOTARGETCPU) $(TARGETCPU)
ADDOPT = -debug inline-debug-info $(ALIGN)
#ADDOPT += -trace
#ADDOPT += -watch

EXEDIR := $(HOME)/bin

ifneq "$(MKLLIB)" ""
  MKLFLAG = -mkl
  MKLFPP = -DUSE_INTEL_MKL
  MKLINC = -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
endif

FPP = -fpp $(METISFPP) $(PBSFPP) $(CGNSFPP) $(MPIFPP) $(MKLFPP) $(OTHERFPP)

INCDIR = $(sort $(CGNSDIR) $(METISDIR))
LIBDIR = $(sort $(HDF5DIR) $(SZIPDIR) $(CGNSDIR) $(METISDIR) $(PBSDIR))

FINC = $(foreach i,$(INCDIR),-I$(i)/include) $(MKLINC) $(OTHERINCDIR)
FLIB = $(foreach i,$(LIBDIR),-L$(i)/lib) $(OTHERLIBDIR)
FLIB += $(MKLLIB) $(METISLIB) $(CGNSLIB) $(HDF5LIB) $(SZIPLIB) $(PBSLIB) $(OTHERLIB)

FCBASE = -traceback -g -assume realloc_lhs

ADDFLAGS = $(MKLFLAG) $(MOD) $(FINC) $(FCADD)

DEBUGFLAGS = -ftrapuv -check all -debug -fp-stack-check

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
        $(OBJDIR)/cpu_info_mod.o


HighlightName := | /bin/grep -hi --color '[^ ]\+\.\(f[0-9]*\|c\)\>'

.PHONY: O2-ip-opt
O2-ip-opt: OPTIM = $(O2OPT) $(TARGET) $(ADDOPT)
O2-ip-opt: optimize

.PHONY: O3-ipo-opt
O3-ipo-opt: OPTIM = $(O3OPT) $(TARGET) $(ADDOPT)
O3-ipo-opt: optimize

.PHONY: broadwell
broadwell: OPTIM = $(O3OPT) -xCORE-AVX2 $(ADDOPT)
broadwell: optimize

.PHONY: haswell
haswell: OPTIM = $(O3OPT) -xCORE-AVX2 $(ADDOPT)
haswell: optimize

.PHONY: ivybridge
ivybridge: OPTIM = $(O3OPT) -xCORE-AVX-I $(ADDOPT)
ivybridge: optimize

.PHONY: sandybridge
sandybridge: OPTIM = $(O3OPT) -xAVX $(ADDOPT)
sandybridge: optimize

.PHONY: westmere
westmere: OPTIM = $(O3OPT) -xSSE4.2 $(ADDOPT)
westmere: optimize

.PHONY: merope
merope: OPTIM = $(O3OPT) -xSSE4.2 $(ADDOPT)
merope: optimize

.PHONY: optimize
optimize: FCFLAG = $(FCBASE) $(OPTIM) $(ADDFLAGS)
optimize: all

.PHONY: debug
debug: FPP += -DDEBUG_ON
debug: FCFLAG = $(FCBASE) $(DEBUGFLAGS) $(ADDFLAGS)
debug: all


.PHONY: deplist
deplist: FCBASE += -gen-dep=deplist.txt
deplist: debug

.PHONY: time
time: Command-Beg := /usr/bin/time -f '%e'
time: Command-End  = 2>&1 | awk '{ printf "%-25s%8.2f seconds", "$<", $$1 }' $(HighlightName)
time: FCFLAG = $(FCBASE) $(O3OPT) $(TARGET) $(ADDOPT) $(ADDFLAGS)
time: $(OBJS)
	$(FC) $(FPP) $(FCFLAG) -o $(EXECUTABLE) $(OBJS) $(FLIB)
	/bin/cp -f $(EXECUTABLE) $(EXEDIR)/$(EXECUTABLE)

.PHONY: all
all: Command-Beg := echo
all: Command-End  = $(HighlightName); $(FC) $(FPP) $(FCFLAG) -c -o $@ $<
all: Command-Cend  = $(HighlightName); $(CC) $(CFLAG) -c -o $@ $<
all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS) $(COBJS)
	$(FC) $(FPP) $(FCFLAG) -o $(EXECUTABLE) $(OBJS) $(COBJS) $(FLIB)
	/bin/cp -f $(EXECUTABLE) $(EXEDIR)/$(EXECUTABLE)

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

.PHONY: clean
clean:
	/bin/rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod $(EXECUTABLE)

.PHONY: dump-variables
dump-variables:
	$(QUIET) echo "FC = $(FC)"
	$(QUIET) echo "FPP = $(FPP)"
	$(QUIET) echo "FCBASE = $(FCBASE)"
	$(QUIET) echo "O3OPT = $(O3OPT)"
	$(QUIET) echo "ADDOPT = $(ADDOPT)"
	$(QUIET) echo "ADDFLAGS = $(ADDFLAGS)"
	$(QUIET) echo "FCFLAG = $(FCFLAG)"
	$(QUIET) echo "OBJDIR = $(OBJDIR)"
	$(QUIET) echo "MODDIR = $(MODDIR)"
	$(QUIET) echo "EXECUTABLE = $(EXECUTABLE)"

include dep.mk
