
include system_variables.mk

OBJDIR := objects/$(BUILDDIR)
MODDIR := modules/$(BUILDDIR)
FUNDIR := Functions

vpath %.o $(OBJDIR)
vpath %.f90 ./
vpath %.mod $(MODDIR)

MOD := -module $(MODDIR)
FSFXFLG := .f90

HDF5LIB := -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz -ldl

O2OPT := -O2 -ip
O3OPT := -O3 -ipo

OPTIM := $(AUTOTARGETCPU) $(TARGETCPU) -debug inline-debug-info $(ALIGN)
#OPTIM += -trace
#OPTIM += -watch

EXEDIR := $(HOME)/bin

ifneq "$(MKLLIB)" ""
  MKLFLAG := -mkl
  MKLFPP := -DUSE_INTEL_MKL
  MKLINC := -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
endif

FPP := -fpp $(METISFPP) $(PBSFPP) -DUSE_MATRIX_PROCEDURES $(MKLFPP) $(OTHERFPP)

INCDIR := $(sort $(CGNSDIR) $(METISDIR))
LIBDIR := $(sort $(HDF5DIR) $(SZIPDIR) $(CGNSDIR) $(METISDIR) $(PBSDIR))

FINC := $(foreach i,$(INCDIR),-I$(i)/include) $(MKLINC) $(OTHERINCDIR)
FLIB := $(foreach i,$(LIBDIR),-L$(i)/lib) $(OTHERLIBDIR)
FLIB += $(MKLLIB) $(METISLIB) $(CGNSLIB) $(HDF5LIB) $(SZIPLIB) $(PBSLIB) $(OTHERLIB)

FCBASE := -ftrapuv -traceback -g

ADDFLAGS := $(MKLFLAG) $(MOD) $(FINC) $(FCADD)

DEBUGFLAGS := -check all -debug -fp-stack-check

OBJS := $(OBJDIR)/kind_types.o \
        $(OBJDIR)/Euler_tri.o \
        $(OBJDIR)/order_mod.o \
        $(OBJDIR)/quadrature_mod.o \
        $(OBJDIR)/geovar_mod.o \
        $(OBJDIR)/flowvar_mod.o \
        $(OBJDIR)/ovar_mod.o \
        $(OBJDIR)/gmsh_mod.o \
        $(OBJDIR)/plot3d_mod.o \
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
        $(OBJDIR)/restart_mod.o


O2-ip-opt: FCFLAG := $(FCBASE) $(O2OPT) $(OPTIM) $(ADDFLAGS)
O2-ip-opt: all


O3-ipo-opt: FCFLAG := $(FCBASE) $(O3OPT) $(OPTIM) $(ADDFLAGS)
O3-ipo-opt: all


debug: FPP += -DDEBUG_ON
debug: FCFLAG := $(FCBASE) $(DEBUGFLAGS) $(ADDFLAGS)
debug: all


deplist: FCBASE += -gen-dep=deplist.txt
deplist: debug


all: $(OBJS)
	$(FC) $(FPP) $(FCFLAG) -o $(EXECUTABLE) $(OBJS) $(FLIB)
	/bin/cp -f $(EXECUTABLE) $(EXEDIR)/$(EXECUTABLE)


%.o:%$(FSFXFLG)
	@echo $(FC) $(FPP) $(FCFLAG) -c -o $@ $< | /bin/grep -Phi --color '[^ ]+\.f[0-9]*';\
	$(FC) $(FPP) $(FCFLAG) -c -o $@ $<


$(OBJS):
	@echo $(FC) $(FPP) $(FCFLAG) -c -o $*.o $(*F)$(FSFXFLG) | /bin/grep -Phi --color '[^ ]+\.f[0-9]*';\
	$(FC) $(FPP) $(FCFLAG) -c -o $*.o $(*F)$(FSFXFLG)


.PHONY: clean
clean:
	/bin/rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod $(EXECUTABLE)


include dep.mk
