#
#
OBJDIR = objects
MODDIR = modules

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

OPTIM = -O3 -ipo -xSSE4.2 #-ipo-jobs4
#OPTIM = -O3 -ip -xSSE4.2

 FPP = -fpp
#FPP += -DPBS_ENV 
 FPP += -DIO_FOR_NAS
 FPP += -DCGNS_ADF

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

FLIB = $(LIBRARY) $(METISLIB) $(CGNSLIB) $(OTHERLIB)

# Various compiler flag options

FCBASE = -vec-report0 -ftrapuv -traceback -g #-warn stderrors

TUNE  = -prof-gen $(OPTIM)
RETUNE = -prof-use $(OPTIM) $(IPOPT)
#
F08STD = -std08
F03STD = -std03
F95STD = -std95
F90STD = -std90
#
FCFLAG = $(FCBASE) 
FCFLAG += $(FPP)
FCFLAG += $(F08STD)
FCFLAG += $(MOD) $(INCLUDE)
#FCFLAG += -shared-intel

#
#FCFLAG += -opt-report=2 -opt-report-file=medium.opt.log
#FCFLAG += -multiple-processes=6
#


#
OBJS = $(OBJDIR)/kind_types.o \
       $(OBJDIR)/Euler_tri.o \
       $(OBJDIR)/metvar_mod.o \
       $(OBJDIR)/ovar_mod.o \
       $(OBJDIR)/geovar_mod.o \
       $(OBJDIR)/flowvar_mod.o \
       $(OBJDIR)/gmsh_mod.o \
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
debug: FCFLAG += -check all -debug -fp-stack-check
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
	$(FC) -check all -debug -g -traceback $(MOD) $(INCLUDE) -o test test_poly.f90
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
$(OBJDIR)/gmsh_mod.o : gmsh_mod.f90 \
             $(OBJDIR)/kind_types.o \
             $(OBJDIR)/geovar_mod.o \
             $(OBJDIR)/ovar_mod.o
$(OBJDIR)/io.o : io.f90 \
       $(OBJDIR)/kind_types.o \
       $(OBJDIR)/geovar_mod.o \
       $(OBJDIR)/ovar_mod.o \
       $(OBJDIR)/metvar_mod.o \
       $(OBJDIR)/flowvar_mod.o \
       $(OBJDIR)/gmsh_mod.o \
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
