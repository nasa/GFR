# Lab

# Fortran MPI compiler
FC := mpiifort

# Compilation flags for this specific architecture
TARGETCPU     := -xHost
AUTOTARGETCPU :=
ALIGN         := -align array16byte

# Library flags, NOTE: HDF5LIB is defined in slv.mk
SZIPLIB  := -lsz
CGNSLIB  := -lcgns
METISLIB := -lmetis5
PBSLIB   :=
#MKLLIB  := -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential

# Paths to directories
HDF5DIR  := $(HOME)
SZIPDIR  := $(HOME)
CGNSDIR  := $(HOME)
METISDIR := $(HOME)
PBSDIR   :=

# Fortran Preprocessor Definitions
METISFPP := -DMETIS_5
PBSFPP   :=

# Other system specific additions
OTHERLIB    :=# -ladvisor -ldl
OTHERLIBDIR :=# -L$(ADVISOR_XE_2016_DIR)/lib64
OTHERINCDIR :=# -I$(ADVISOR_XE_2016_DIR)/include/lib64
OTHERFPP    :=#
