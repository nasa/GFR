
###########################################################
###  Compiler/System specific definitions to build GFR  ###
###########################################################

# MPI compiler wrappers for Fortran and C
FC := $(GFR_FC)
CC := $(GFR_CC)

# Compiler family/brand
#    Working values: intel, oneapi
#    Future  values: aocc, cray, gnu, nvidia
COMPFAMILY := $(GFR_CompFamily)

# Batch Scheduler
#    Working values: pbs, "", not-defined (last two are identical)
#    Future  values: slurm
 BATCH_SCHEDULER := $(GFR_Scheduler)

# Directory to put Executables
EXEDIR := $(CURDIR)/bin

# Paths to install directories of required libraries
# not including lib or include sub-directories, e.g.,
#     HDF5DIR := /home/username/hdf5
#         where
#       $(HDF5DIR)/lib      => Contains all HDF5 shared and/or static libraries
#       $(HDF5DIR)/include  => Contains all HDF5 headers and Fortran modules
#
# NOTE: I use the environment variables $HDF5_ROOT, $CGNS_ROOT,
#       $METIS_ROOT, $PARMETIS_ROOT, and $MKLROOT to store the
#       paths to these installed libraries.
HDF5DIR     := $(HDF5_ROOT)
CGNSDIR     := $(CGNS_ROOT)
METISDIR    := $(METIS_ROOT)
PARMETISDIR := $(PARMETIS_ROOT)
MKLDIR      := $(MKLROOT)
#TECDIR     := /home/sspiegel/Downloads/tecio/tecio/install

###############################################################
# Fortran preprocessor definitions to enable/disable features #
###############################################################

# Toggle the use of OpenMP for thread parallelism
#     ENABLEOMP := 1         => Enable  OpenMP
#     ENABLEOMP :=           => Disable OpenMP
# NOTE: OpenMP currently only exists in certain non-solver parts of the code,
#       but should hopefully be added to the solver in the near future.
 ENABLEOMP := 
#ENABLEOMP := 1

# Set the integer type used by the METIS library
#     METISFPP :=                => METIS 5.1.0 (built 32-bit)
#     METISFPP := -DMETIS_64_BIT => METIS 5.1.0 (built 64-bit)
METISFPP := 

# Select whether to link to the Intel BLAS Math Library
#     ENABLE_MKL_BLAS := 1     => Do     link to Intel MKL BLAS library
#     ENABLE_MKL_BLAS :=       => Do not link to Intel MKL BLAS library
ENABLE_MKL_BLAS :=  

# Select whether to link to the Intel LAPACK Math Library
#     ENABLE_MKL_LAPACK := 1     => Do     link to Intel MKL LAPACK library
#     ENABLE_MKL_LAPACK :=       => Do not link to Intel MKL LAPACK library
ENABLE_MKL_LAPACK := 1

# Other system specific additions
OTHERLIB    :=# -ladvisor -ldl
OTHERLIBDIR :=# -L$(ADVISOR_XE_2016_DIR)/lib64
OTHERINCDIR :=# -I$(ADVISOR_XE_2016_DIR)/include/lib64
OTHERFPP    :=

##############################################################################
# Fortran preprocessor definitions that are automatically set based on above #
#  settings or definitions. Manually setting these is not recommended these  #
#   is NOT recommended unless you are certain you know what you are doing.   #
##############################################################################

# Enable ParMETIS functionality if the PARMETISDIR variable is defined
ifdef PARMETISDIR
  PARMETISFPP := -DENABLE_PARMETIS
endif

# Enable/Disable code accessing the PBS environment
#     PBSFPP := -DPBS_ENV  => Enable  code for PBS environment
#     PBSFPP :=            => Disable code for PBS environment
ifeq "$(BATCH_SCHEDULER)" "pbs"
  HPCFPP := -DPBS_ENV
endif

# Select the MPI Fortran module to access
#     MPIFPP := -DUSE_MPI_F08  => USE MPI_F08
#     MPIFPP :=                => USE MPI
MPIFPP := -DUSE_MPI_F08
ifneq "$(FC)" "ftn"
ifneq "$(COMPFAMILY)" "intel"
MPIFPP :=
endif
endif
