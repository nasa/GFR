#ifdef USE_MPI_F08

#define _MPI_MODULE_         mpi_f08
#define _MPI_DATA_TYPE_      type(MPI_DATATYPE)
#define _MPI_FILE_TYPE_      type(MPI_FILE)
#define _MPI_INFO_TYPE_      type(MPI_INFO)
#define _MPI_COMM_TYPE_      type(MPI_COMM)
#define _MPI_REQUEST_TYPE_   type(MPI_REQUEST)
#define _MPI_STATUS_TYPE_    type(MPI_STATUS)
#define _MPI_STATUS_SIZE_    1
#define _MPI_COMM_WORLD_VAL_ MPI_COMM_WORLD%MPI_VAL
#define _MPI_COMM_SELF_VAL_  MPI_COMM_SELF%MPI_VAL

#else

#define _MPI_MODULE_         mpi
#define _MPI_DATA_TYPE_      integer(INT_MPI)
#define _MPI_FILE_TYPE_      integer(INT_MPI)
#define _MPI_INFO_TYPE_      integer(INT_MPI)
#define _MPI_COMM_TYPE_      integer(INT_MPI)
#define _MPI_REQUEST_TYPE_   integer(INT_MPI)
#define _MPI_STATUS_TYPE_    integer(INT_MPI)
#define _MPI_STATUS_SIZE_    MPI_STATUS_SIZE
#define _MPI_COMM_WORLD_VAL_ MPI_COMM_WORLD
#define _MPI_COMM_SELF_VAL_  MPI_COMM_SELF

#endif
