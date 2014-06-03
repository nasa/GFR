#
#
EXECUTABLE = recon
OBJDIR = objects
MODDIR = modules
#EXEDIR = 
EXEDIR = /home/sspiegel/bin

vpath %.o $(OBJDIR)
vpath %.f90 ./
vpath %.mod $(MODDIR)

#MAKE = make
 MAKE = make -j 4

 FC = mpif90
#FC = h5pfc#-show
MOD = -module $(MODDIR)
FSFXFLG = .f90

# Include and library directories 
 INCDIR = -I/usr/local/include
#INCDIR = -I/u/sspiegel/include
 LOCALDIR = -L/usr/local/lib
#LOCALDIR = -L/u/sspiegel/lib -L/nasa/szip/2.1/lib

# MeTiS
METISLIB = -lmetis

# CGNS 3.1 Library
CGNS = -lcgns3
PCGNS = -lpcgns
HDF5 = -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz -ldl -lm
 CGNS31LIB = $(CGNS) $(PCGNS) $(HDF5)
#CGNS31LIB = $(CGNS) $(HDF5) -lsz

# CGNS 2.5 Library
CGNS25LIB = -lcgns2

# Linked Libraries
 FLIB = $(LOCALDIR) $(METISLIB) $(CGNS25LIB)
#FLIB = $(LOCALDIR) $(METISLIB) $(CGNS31LIB)
#FLIB += $(PBSTIME)

# Various compiler flag options
PRECISION = -fp-model precise -prec-sqrt -prec-div -assume protect_parens,minus0
 FCBASE = -vec-report0 -ftrapuv -traceback -g #-warn stderrors
#FCBASE = -ftrapuv -traceback -g -warn stderrors
 FPP = -fpp
#FPP += -DPBS_ENV 
 FPP += -DIO_FOR_NAS 
 FPP += -DCGNS_ADF 
#OPTIM = -O2 -ip -fpe0
 OPTWM = -O3 -ip -xSSE4.2
#OPTWM = -O3 -ipo -xSSE4.2
#OPTWM = -O3 -ipo -ipo-jobs4 -xSSE4.2
 OPTSB = -O3 -ipo -xAVX
 OPTIB = -O3 -ipo -xAVX
 OPTALL = -O3 -ipo -axAVX -xSSE4.1
 OPTIM = $(OPTWM)
DEBUG = -check all -debug -fp-stack-check

TUNE  = -prof-gen $(OPTIM)
RETUNE = -prof-use $(OPTIM) $(IPOPT)
PROFILE = -p -fno-inline-functions
#PROFILE = -p inline-debug-info $(OPTIM)
OPTIMIZE = -O3 $(IPOPT)
#
F08STD = -std08
F03STD = -std03
F95STD = -std95
F90STD = -std90
#
FCFLAG = $(FCBASE)
FCFLAG += $(FPP)
#FCFLAG += $(PRECISION)
#FCFLAG += $(DEBUG)
#FCFLAG += -gen-dep=deplist.txt
#FCFLAG += $(OPTIM)
 FCFLAG += $(F08STD)
#FCFLAG += -opt-report=2 -opt-report-file=medium.opt.log
#FCFLAG += -multiple-processes=6
#
FCFLAG += $(MOD) $(INCDIR)
#
#
OBJS = $(OBJDIR)/kind_types.o \
       $(OBJDIR)/Euler_tri.o \
       $(OBJDIR)/metvar_mod.o \
       $(OBJDIR)/ovar_mod.o \
       $(OBJDIR)/geovar_mod.o \
       $(OBJDIR)/flowvar_mod.o \
       $(OBJDIR)/io.o \
       $(OBJDIR)/metrics.o \
       $(OBJDIR)/poly.o \
       $(OBJDIR)/startup.o \
       $(OBJDIR)/initialize.o \
       $(OBJDIR)/generic.o \
       $(OBJDIR)/triangle.o \
       $(OBJDIR)/flux.o \
       $(OBJDIR)/parallel.o
#
optimize: FCFLAG += $(OPTIM)
optimize: all
#
debug: FCFLAG += $(DEBUG)
debug: all
#
optprofile: FCFLAG += $(OPTIM) $(PROFILE)
optprofile: all
#
debprofile: FCFLAG += $(DEBUG) $(PROFILE)
debprofile: all
#
profile: FCFLAG += $(PROFILE)
profile: all
#
vtune: FCFLAG += $(OPTIM) -debug inline-debug-info
vtune: all
#
precise: FCFLAG += $(OPTIM) $(PRECISION)
precise: all
#
all: $(OBJS)
	$(FC) $(FCFLAG) -o $(EXECUTABLE) $(OBJS) $(FLIB)
	/bin/cp -f $(EXECUTABLE) $(EXEDIR)/$(EXECUTABLE)
#
%.o:%$(FSFXFLG)
	@echo $(FC) $(FCFLAG) -c -o $*.o $(*F)$(FSFXFLG) | /bin/grep -Phi --color '[^ ]+\.f[0-9]*';\
	$(FC) $(FCFLAG) -c -o $*.o $(*F)$(FSFXFLG)
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
	$(FC) -check all -debug -g -traceback $(MOD) $(INCDIR) -o test test_poly.f90
#
.PHONY: clean
clean:
	/bin/rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod $(EXECUTABLE)
#
$(OBJDIR)/kind_types.o : kind_types.f90
$(OBJDIR)/Euler_tri.o : Euler_tri.f90 \
               $(OBJDIR)/kind_types.o \
               $(OBJDIR)/geovar_mod.o \
               $(OBJDIR)/ovar_mod.o \
               $(OBJDIR)/metvar_mod.o \
               $(OBJDIR)/flowvar_mod.o \
               $(OBJDIR)/io.o \
               $(OBJDIR)/metrics.o \
               $(OBJDIR)/startup.o \
               $(OBJDIR)/initialize.o \
               $(OBJDIR)/generic.o \
               $(OBJDIR)/flux.o \
               $(OBJDIR)/parallel.o
$(OBJDIR)/metvar_mod.o : metvar_mod.f90 $(OBJDIR)/kind_types.o
$(OBJDIR)/ovar_mod.o : ovar_mod.f90 $(OBJDIR)/kind_types.o
$(OBJDIR)/geovar_mod.o : geovar_mod.f90 $(OBJDIR)/kind_types.o
$(OBJDIR)/flowvar_mod.o : flowvar_mod.f90 $(OBJDIR)/kind_types.o
$(OBJDIR)/poly.o : poly.f90 $(OBJDIR)/kind_types.o
$(OBJDIR)/triangle.o : triangle.f90 \
             $(OBJDIR)/kind_types.o \
             $(OBJDIR)/startup.o \
             $(OBJDIR)/poly.o
$(OBJDIR)/parallel.o : parallel.f90 \
             $(OBJDIR)/kind_types.o \
             $(OBJDIR)/geovar_mod.o \
             $(OBJDIR)/flowvar_mod.o
$(OBJDIR)/metrics.o : metrics.f90 \
            $(OBJDIR)/kind_types.o \
            $(OBJDIR)/geovar_mod.o \
            $(OBJDIR)/ovar_mod.o \
            $(OBJDIR)/metvar_mod.o \
            $(OBJDIR)/poly.o \
            $(OBJDIR)/triangle.o \
            $(OBJDIR)/parallel.o \
            Functions/chop.f90
$(OBJDIR)/io.o : io.f90 \
       $(OBJDIR)/kind_types.o \
       $(OBJDIR)/geovar_mod.o \
       $(OBJDIR)/ovar_mod.o \
       $(OBJDIR)/metvar_mod.o \
       $(OBJDIR)/flowvar_mod.o \
       $(OBJDIR)/parallel.o \
       $(OBJDIR)/generic.o \
       Functions/vortex.f90
$(OBJDIR)/startup.o : startup.f90 \
            $(OBJDIR)/kind_types.o \
            $(OBJDIR)/geovar_mod.o \
            $(OBJDIR)/ovar_mod.o \
            $(OBJDIR)/metvar_mod.o \
            $(OBJDIR)/poly.o \
            Functions/chop.f90
$(OBJDIR)/initialize.o : initialize.f90 \
               $(OBJDIR)/kind_types.o \
               $(OBJDIR)/geovar_mod.o \
               $(OBJDIR)/ovar_mod.o \
               $(OBJDIR)/metvar_mod.o \
               $(OBJDIR)/flowvar_mod.o \
               $(OBJDIR)/generic.o \
               $(OBJDIR)/io.o \
               Functions/vortex.f90
$(OBJDIR)/generic.o : generic.f90 \
            $(OBJDIR)/kind_types.o \
            $(OBJDIR)/geovar_mod.o \
            $(OBJDIR)/ovar_mod.o \
            $(OBJDIR)/metvar_mod.o \
            $(OBJDIR)/flowvar_mod.o \
            $(OBJDIR)/parallel.o \
            Functions/vortex.f90
$(OBJDIR)/flux.o : flux.f90 \
         $(OBJDIR)/kind_types.o \
         $(OBJDIR)/geovar_mod.o \
         $(OBJDIR)/ovar_mod.o \
         $(OBJDIR)/metvar_mod.o \
         $(OBJDIR)/flowvar_mod.o \
         $(OBJDIR)/parallel.o \
         $(OBJDIR)/generic.o
