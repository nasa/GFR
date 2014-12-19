#
#
OBJDIR = objects/$(BUILDDIR)
MODDIR = modules/$(BUILDDIR)
FUNDIR = Functions

vpath %.o $(OBJDIR)
vpath %.f90 ./
vpath %.mod $(MODDIR)

FC = mpif90
MOD = -module $(MODDIR)
FSFXFLG = .f90

# MeTiS
METISLIB = -lmetis

# PBS Time
PBSTIME = -lpbs_time_left

# HDF5 Library
HDF5 = -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz -ldl -lm

# CGNS 3.2 Library
CGNS32 = -lcgns $(HDF5)
# CGNS 3.1 Library
CGNS31 = -lcgns -lpcgns $(HDF5)
# CGNS 2.5 Library
CGNS25 = -lcgns2

####################################
####################################
###   NAS SPECIFIC DEFINITIONS   ###
####################################
####################################

### HOMEDIR = /u/sspiegel

### CGNS31LIB = $(CGNS) $(HDF5) -lsz
### LOCALLIB = -L/nasa/szip/2.1/lib
### OTHERLIB = $(PBSTIME)

### OPTIM = -O3 -ipo -axAVX -xSSE4.1 #-ipo-jobs4

###  FPP = -fpp
###  FPP += -DPBS_ENV 
### #FPP += -DIO_FOR_NAS
### #FPP += -DCGNS_ADF

####################################
####################################
###   LAB SPECIFIC DEFINITIONS   ###
####################################
####################################

HOMEDIR = /home/sspiegel

LOCALLIB =
OTHERLIB =

 OPTIM = -O2 -ip -xSSE4.2 #-v #-watch
 IPO_OPT = -O3 -ipo #-ipo-jobs4
#OPTIM = -O2 -xSSE4.2 #-ipo-jobs4
#OPTIM = -O3 -ip -xSSE4.2
#OPTIM += -align array32byte

 FPP = -fpp
#FPP += -DPBS_ENV 
 FPP += -DIO_FOR_NAS
#FPP += -DCGNS_ADF
#FPP += -DUSE_INTEL_MKL
#FPP += -DCGNS_ADF2

##############################
##############################
###   COMMON DEFINITIONS   ###
##############################
##############################

 CGNSLIB = $(CGNS32)
#CGNSLIB = $(CGNS31)
#CGNSLIB = $(CGNS25)

EXEDIR = $(HOMEDIR)/bin
LIBDIR = $(HOMEDIR)/lib
INCDIR = $(HOMEDIR)/include

INCLUDE = -I$(INCDIR) -I/usr/local/include
LIBRARY = -L$(LIBDIR) -L/usr/local/lib $(LOCALLIB)
#LIBRARY += -L$(ADVISOR_XE_2013_DIR)/lib64 -ladvisor -ldl
#INCLUDE += -I$(ADVISOR_XE_2013_DIR)/include/intel64

FLIB = $(LIBRARY) $(METISLIB) $(CGNSLIB) $(OTHERLIB)

# Various compiler flag options

FCBASE = -ftrapuv -traceback -g #-warn stderrors
#FCBASE += -i8

TUNE  = -prof-gen $(OPTIM)
RETUNE = -prof-use $(OPTIM) $(IPOPT)
#
F08STD = -std08
F03STD = -std03
F95STD = -std95
F90STD = -std90
#
FCFLAG = $(FCBASE)
#FCFLAG += -multiple-processes=4
FCFLAG += $(FPP)
#FCFLAG += $(F08STD)
FCFLAG += $(MOD) $(INCLUDE)
#FCFLAG += -shared-intel
#FCFLAG += -parallel

#
#FCFLAG += -opt-report=2 -opt-report-file=medium.opt.log
#FCFLAG += -multiple-processes=6
#

# MKL options
#BLASLIB  = -lmkl_blas95_lp64
#MKLLIB += -lmkl_intel_lp64
#MKLLIB += -lmkl_sequential
#MKLLIB += -lmkl_intel_thread
#MKLLIB += -lmkl_core
#OPENMP = -openmp
#FLIB += -mkl $(BLASLIB) -Wl,--start-group $(MKLLIB) -Wl,--end-group  $(OPENMP) -lpthread
#FLIB += $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -openmp -lpthread -lm

#
OBJS = $(OBJDIR)/kind_types.o \
       $(OBJDIR)/Euler_tri.o \
       $(OBJDIR)/metvar_mod.o \
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
       $(OBJDIR)/startup.o \
       $(OBJDIR)/initialize.o \
       $(OBJDIR)/generic.o \
       $(OBJDIR)/bc_mod.o \
       $(OBJDIR)/triangle.o \
       $(OBJDIR)/flux.o \
       $(OBJDIR)/filter_mod.o \
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
       $(OBJDIR)/limiters.o
#
optimize: FCFLAG += $(OPTIM) $(ADDFLAGS)
optimize: all
#
ipo-opt: FCFLAG += $(OPTIM) $(IPO_OPT)
ipo-opt: all
#
debug: FCFLAG += -check all -debug -fp-stack-check -DDEBUG_ON
debug: all
#
optprofile: FCFLAG += -p -debug inline-debug-info
optprofile: optimize
#
debprofile: FCFLAG += -p -fno-inline-functions
debprofile: debug
#
profile: FCFLAG += -p -fno-inline-functions
profile: all
#
vtune: FCFLAG += -debug inline-debug-info
vtune: optimize
#
precise: FCFLAG += -fp-model precise -prec-sqrt -prec-div -assume protect_parens,minus0
precise: optimize
#
deplist: FCFLAG += -gen-dep=deplist.txt
deplist: debug

guide: FCFLAG += -parallel -O3 -ip -xSSE4.2 -guide=4 -guide-file=guide_files/$(*F).txt
guide: $(OBJS)
#
all: $(OBJS)
	$(FC) $(FCFLAG) -o $(EXECUTABLE) $(OBJS) $(FLIB)
	/bin/cp -f $(EXECUTABLE) $(EXEDIR)/$(EXECUTABLE)
#
%.o:%$(FSFXFLG)
	@echo $(FC) $(FCFLAG) -c -o $@ $< | /bin/grep -Phi --color '[^ ]+\.f[0-9]*';\
	$(FC) $(FCFLAG) -c -o $@ $<
#
$(OBJS):
	@echo $(FC) $(FCFLAG) -c -o $*.o $(*F)$(FSFXFLG) | /bin/grep -Phi --color '[^ ]+\.f[0-9]*';\
	$(FC) $(FCFLAG) -c -o $*.o $(*F)$(FSFXFLG)
#
basic: $(OBJS)
	@echo $(FC) $(FCFLAG) -o $(EXECUTABLE) $(OBJS) $(FLIB) | /bin/grep -Phi --color '[^ ]+\.f[0-9]*[^ ]*';\
	$(FC) $(FCFLAG) -o $(EXECUTABLE) $(OBJS) $(FLIB)
#
base: $(OBJS)
	$(FC) $(FCFLAG) -o $(EXECUTABLE) $(OBJS) $(FLIB)
#
retune: $(OBJS)
	$(FC) $(FCFLAG) $(RETUNE) -o $(EXECUTABLE) $(OBJS) $(FLIB)
#
tune: $(OBJS)
	$(FC) $(FCFLAG) $(TUNE) -o $(EXECUTABLE) $(OBJS) $(FLIB)
#
##prof: $(OBJS)
##	$(FC) $(FCFLAG) $(PROFILE) -o $(EXECUTABLE) $(OBJS) $(FLIB)
#
##optimize: $(OBJS)
##	$(FC) $(FCFLAG) $(OPTIMIZE) -o $(EXECUTABLE) $(OBJS) $(FLIB)
#
##precise: $(OBJS)
##	$(FC) $(FCPREC) -o $(EXECUTABLE) $(OBJS) $(FLIB)
#
f08: $(OBJS)
	$(FC) $(F08STD) $(MOD) -o $(EXECUTABLE) $(OBJS) $(FLIB)
#
f03: $(OBJS)
	$(FC) $(F03STD) $(MOD) -o $(EXECUTABLE) $(OBJS) $(FLIB)
#
f95: $(OBJS)
	$(FC) $(F95STD) $(MOD) -o $(EXECUTABLE) $(OBJS) $(FLIB)
#
f90: $(OBJS)
	$(FC) $(F90STD) $(MOD) -o $(EXECUTABLE) $(OBJS) $(FLIB)
#
test_poly: test_poly.f90 kind_types.f90 poly.f90
	$(FC) -check all -debug -g -traceback $(MOD) $(INCLUDE) -o test test_poly.f90
#
.PHONY: clean
clean:
	/bin/rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod $(EXECUTABLE)
#
include dep.mk
