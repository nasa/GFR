# Lab

# Fortran MPI compiler
FC := mpiifort
CC := mpiicc

#COMPFAMILY := gcc
#COMPFAMILY := pgi
 COMPFAMILY := intel

# Library flags, NOTE: HDF5LIB is defined in slv.mk
SZIPLIB  := -lsz -lz
CGNSLIB  := -lcgns
METISLIB := -lmetis
PBSLIB   :=
#MKLLIB  := -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential

# Paths to directories
HDF5DIR  := $(HDF5_ROOT)
SZIPDIR  := $(SZIP_ROOT)
CGNSDIR  := $(CGNS_ROOT)
METISDIR := $(METIS_ROOT)
PBSDIR   :=

# Fortran Preprocessor Definitions
METISFPP := -DMETIS_5
PBSFPP   :=
CGNSFPP  := -DCGNS_3_3
MPIFPP   := -DUSE_MPI_F08

# Other system specific additions
OTHERLIB    :=# -ladvisor -ldl
OTHERLIBDIR :=# -L$(ADVISOR_XE_2016_DIR)/lib64
OTHERINCDIR :=# -I$(ADVISOR_XE_2016_DIR)/include/lib64
OTHERFPP    :=

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
  IPFLAG  := -ip
  IPOFLAG := -ipo
  TRACEFLAG := -traceback
  REALLOCFLAG := -assume realloc_lhs
  CHECKFLAG := -check all
  DEBUGFLAG := -debug -ftrapuv -fp-stack-check
  ADDOPT := -debug inline-debug-info -align array32byte
  EXTFLAGS :=
  OTHERINCDIR :=
  EXTOBJS  :=
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
  IPFLAG  := 
  IPOFLAG := 
  TRACEFLAG := -fbacktrace
  REALLOCFLAG := -frealloc-lhs
  CHECKFLAG := -fcheck=all
  DEBUGFLAG := -Wall -Wextra -Wuninitialized
  ADDOPT := 
  EXTFLAGS := -DDISABLE_DTIO -DSPECIAL_FOR_GCC
  OTHERINCDIR := -I./
  PBSLIB   :=
  EXTOBJS  :=
else ifeq "$(COMPFAMILY)" "pgi"
  # GCC flags
  MODFLAG := -module
  ARCHFLAG := -tp=
  TUNEFLAG := -tp=
  WESARCH := nahalem
  SANARCH := sandybridge
  IVYARCH := sandybridge
  HASARCH := haswell
  BROARCH := haswell
  FPPFLAG := -Mpreprocess
  IPFLAG  := 
  IPOFLAG := 
  TRACEFLAG := -traceback
  REALLOCFLAG := -Mallocatable=03
  CHECKFLAG := -Mbounds -Mchkfpstk -Mchkptr -Mchkstk
  DEBUGFLAG := -Minform=inform -Ktrap=fp
  ADDOPT := 
  EXTFLAGS := -DSPECIAL_FOR_PGI -DDISABLE_QP -DDISABLE_DTIO -Mbackslash
  OTHERINCDIR := -I./
  PBSLIB   := 
  EXTOBJS  := $(METIS_ROOT)/include/metis.o
endif


# Compilation flags for this architecture
TARGETCPU := $(ARCHFLAG)$(SANARCH)
AUTOTARGETCPU := $(TUNEFLAG)$(IVYARCH) $(TUNEFLAG)$(HASARCH) $(TUNEFLAG)$(BROARCH)
