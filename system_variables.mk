
###########################################################
###  Compiler/System specific definitions to build GFR  ###
###########################################################

# MPI compiler wrappers for Fortran and C
FC := mpif90
CC := mpicc

# Compiler family/brand
#COMPFAMILY := gcc
#COMPFAMILY := pgi
 COMPFAMILY := intel

# Directory to put Executables
EXEDIR := $(CURDIR)/bin

# Paths to install directories of required libraries
# not including lib or include sub-directories, e.g.,
#     HDF5DIR := /home/username/hdf5
#         where
#       $(HDF5DIR)/lib      => Contains all HDF5 shared and/or static libraries
#       $(HDF5DIR)/include  => Contains all HDF5 headers and Fortran modules
#
# NOTE: I use the environment variables $HDF5_ROOT, $CGNS_ROOT, $METIS_ROOT, 
#       and $SZIP_ROOT (if HDF5 is built with szip support) to store the paths
#       to these installed libraries.
HDF5DIR  := $(HDF5_ROOT)
SZIPDIR  := $(SZIP_ROOT)
CGNSDIR  := $(CGNS_ROOT)
METISDIR := $(METIS_ROOT)

# Fortran Preprocessor Definitions

# Select which version of the METIS library to use
#     METISFPP := -DMETIS_5 -DMETIS_32_BIT => METIS 5.1.0 (built 32-bit)
#     METISFPP := -DMETIS_5 -DMETIS_64_BIT => METIS 5.1.0 (built 64-bit)
#     METISFPP := -DMETIS_4                => METIS 4.0.3
METISFPP := -DMETIS_5 -DMETIS_32_BIT

# Enable/Disable code accessing the PBS environment
#     PBSFPP := -DPBS_ENV  => Enable  code for PBS environment
#     PBSFPP :=            => Disable code for PBS environment
PBSFPP :=

# Select which version of the CGNS library to use
#     CGNSFPP := -DCGNS_3_3  => CGNS 3.3.x
#     CGNSFPP :=             => CGNS 3.2.x
CGNSFPP := -DCGNS_3_3

# Select the MPI Fortran module to access
#     MPIFPP := -DUSE_MPI_F08  => USE MPI_F08
#     MPIFPP :=                => USE MPI
MPIFPP := -DUSE_MPI_F08
ifeq "$(COMPFAMILY)" "pgi"
MPIFPP :=
endif

# Select whether to link to the Intel Math Libraries
#     MKLFPP := -DUSE_INTEL_MKL  => Do     link to Intel MKL
#     MKLFPP :=                  => Do not link to Intel MKL
MKLFPP :=

# Other system specific additions
OTHERLIB    :=# -ladvisor -ldl
OTHERLIBDIR :=# -L$(ADVISOR_XE_2016_DIR)/lib64
OTHERINCDIR :=# -I$(ADVISOR_XE_2016_DIR)/include/lib64
OTHERFPP    :=

