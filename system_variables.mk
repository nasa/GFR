# Lab

# Fortran MPI compiler
FC := mpiifort

# Compilation flags for this specific architecture
TARGETCPU     := -xSSE4.2
AUTOTARGETCPU :=
ALIGN         := -align array32byte

# Library flags, NOTE: HDF5LIB is defined in slv.mk
SZIPLIB  := -lsz -lz
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
OTHERFPP    :=# -DUSE_TYPEBOUND_FACE_PACK
