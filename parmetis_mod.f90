module parmetis_mod
  !
#include <mpi_defs.h>
  !
  use module_kind_types, only : INT_MPI
  !
#ifdef USE_MPI_F08
  use _MPI_MODULE_, only : MPI_COMM
  use _MPI_MODULE_, only : MPI_DATATYPE
#endif
#ifdef METIS_64_BIT
  use _MPI_MODULE_, only : pm_mpi_inttyp => MPI_LONG_LONG_INT
#else
  use _MPI_MODULE_, only : pm_mpi_inttyp => MPI_INT
#endif
  !
  use metis_mod, only : idx_t, real_t
  !
  implicit none
  !
  private
  !
  public :: ParMETIS_PartKway
  public :: ParMETIS_PartGeomKway
  public :: ParMETIS_PartGeom
  public :: ParMETIS_RefineKway
  public :: ParMETIS_AdaptiveRepart
  public :: ParMETIS_Mesh2Dual
  public :: ParMETIS_PartMeshKway
  public :: ParMETIS_NodeND
  public :: ParMETIS_V32_NodeND
  public :: ParMETIS_SerialNodeND
  !
  interface
    !===========================================================================
    !===========================================================================
    function ParMETIS_V3_PartKway_C_Interface(vtxdist, xadj, adjncy, vwgt, &
                                              adjwgt, wgtflag, numflag, ncon, &
                                              nparts, tpwgts, ubvec, options, &
                                              edgecut, part, comm) &
                                         bind(c,name="ParMETIS_V3_PartKway")
      use, intrinsic :: iso_c_binding, only : c_int
      import :: idx_t, real_t
      integer(c_int) :: ParMETIS_V3_PartKway_C_Interface
      !
      integer(idx_t) :: vtxdist(*)
      integer(idx_t) :: xadj(*)
      integer(idx_t) :: adjncy(*)
      integer(idx_t) :: vwgt(*)
      integer(idx_t) :: adjwgt(*)
      integer(idx_t) :: wgtflag
      integer(idx_t) :: numflag
      integer(idx_t) :: ncon
      integer(idx_t) :: nparts
      real(real_t)   :: tpwgts(*)
      real(real_t)   :: ubvec(*)
      integer(idx_t) :: options(*)
      integer(idx_t) :: edgecut
      integer(idx_t) :: part(*)
      integer(c_int) :: comm
      !
    end function ParMETIS_V3_PartKway_C_Interface
    !===========================================================================
    !===========================================================================
    function ParMETIS_V3_PartGeomKway_C_Interface(vtxdist, xadj, adjncy, vwgt, &
                                                  adjwgt, wgtflag, numflag, &
                                                  ndims, xyz, ncon, nparts, &
                                                  tpwgts, ubvec, options, &
                                                  edgecut, part, comm) &
                                         bind(c,name="ParMETIS_V3_PartGeomKway")
      use, intrinsic :: iso_c_binding, only : c_int
      import :: idx_t, real_t
      integer(c_int) :: ParMETIS_V3_PartGeomKway_C_Interface
      !
      integer(idx_t) :: vtxdist(*)
      integer(idx_t) :: xadj(*)
      integer(idx_t) :: adjncy(*)
      integer(idx_t) :: vwgt(*)
      integer(idx_t) :: adjwgt(*)
      integer(idx_t) :: wgtflag
      integer(idx_t) :: numflag
      integer(idx_t) :: ndims
      real(real_t)   :: xyz(*)
      integer(idx_t) :: ncon
      integer(idx_t) :: nparts
      real(real_t)   :: tpwgts(*)
      real(real_t)   :: ubvec(*)
      integer(idx_t) :: options(*)
      integer(idx_t) :: edgecut
      integer(idx_t) :: part(*)
      integer(c_int) :: comm
      !
    end function ParMETIS_V3_PartGeomKway_C_Interface
    !===========================================================================
    !===========================================================================
    function ParMETIS_V3_PartGeom_C_Interface(vtxdist, ndims, xyz, part, comm) &
                                         bind(c,name="ParMETIS_V3_PartGeom")
      use, intrinsic :: iso_c_binding, only : c_int
      import :: idx_t, real_t
      integer(c_int) :: ParMETIS_V3_PartGeom_C_Interface
      !
      integer(idx_t) :: vtxdist(*)
      integer(idx_t) :: ndims
      real(real_t)   :: xyz(*)
      integer(idx_t) :: part(*)
      integer(c_int) :: comm
      !
    end function ParMETIS_V3_PartGeom_C_Interface
    !===========================================================================
    !===========================================================================
    function ParMETIS_V3_RefineKway_C_Interface(vtxdist, xadj, adjncy, vwgt, &
                                                adjwgt, wgtflag, numflag, &
                                                ncon, nparts, tpwgts, ubvec, &
                                                options, edgecut, part, comm) &
                                           bind(c,name="ParMETIS_V3_RefineKway")
      use, intrinsic :: iso_c_binding, only : c_int
      import :: idx_t, real_t
      integer(c_int) :: ParMETIS_V3_RefineKway_C_Interface
      !
      integer(idx_t) :: vtxdist(*)
      integer(idx_t) :: xadj(*)
      integer(idx_t) :: adjncy(*)
      integer(idx_t) :: vwgt(*)
      integer(idx_t) :: adjwgt(*)
      integer(idx_t) :: wgtflag
      integer(idx_t) :: numflag
      integer(idx_t) :: ncon
      integer(idx_t) :: nparts
      real(real_t)   :: tpwgts(*)
      real(real_t)   :: ubvec(*)
      integer(idx_t) :: options(*)
      integer(idx_t) :: edgecut
      integer(idx_t) :: part(*)
      integer(c_int) :: comm
      !
    end function ParMETIS_V3_RefineKway_C_Interface
    !===========================================================================
    !===========================================================================
    function ParMETIS_V3_AdaptiveRepart_C_Interface(vtxdist, xadj, adjncy, &
                                                    vwgt, vsize, adjwgt, &
                                                    wgtflag, numflag, ncon, &
                                                    nparts, tpwgts, ubvec, &
                                                    ipc2redist, options, &
                                                    edgecut, part, comm) &
                                       bind(c,name="ParMETIS_V3_AdaptiveRepart")
      use, intrinsic :: iso_c_binding, only : c_int
      import :: idx_t,real_t
      integer(c_int) :: ParMETIS_V3_AdaptiveRepart_C_Interface
      !
      integer(idx_t) :: vtxdist(*)
      integer(idx_t) :: xadj(*)
      integer(idx_t) :: adjncy(*)
      integer(idx_t) :: vwgt(*)
      integer(idx_t) :: vsize(*)
      integer(idx_t) :: adjwgt(*)
      integer(idx_t) :: wgtflag
      integer(idx_t) :: numflag
      integer(idx_t) :: ncon
      integer(idx_t) :: nparts
      real(real_t)   :: tpwgts(*)
      real(real_t)   :: ubvec(*)
      real(real_t)   :: ipc2redist
      integer(idx_t) :: options(*)
      integer(idx_t) :: edgecut
      integer(idx_t) :: part(*)
      integer(c_int) :: comm
      !
    end function ParMETIS_V3_AdaptiveRepart_C_Interface
    !===========================================================================
    !===========================================================================
    function ParMETIS_V3_Mesh2Dual_C_Interface(elmdist, eptr, eind, numflag, &
                                               ncommonnodes, xadj, adjncy, &
                                               comm) &
                                          bind(c,name="ParMETIS_V3_Mesh2Dual")
      use, intrinsic :: iso_c_binding, only : c_int,c_ptr
      import :: idx_t
      integer(c_int) :: ParMETIS_V3_Mesh2Dual_C_Interface
      !
      integer(idx_t) :: elmdist(*)
      integer(idx_t) :: eptr(*)
      integer(idx_t) :: eind(*)
      integer(idx_t) :: numflag
      integer(idx_t) :: ncommonnodes
      type(c_ptr)    :: xadj
      type(c_ptr)    :: adjncy
      integer(c_int) :: comm
      !
    end function ParMETIS_V3_Mesh2Dual_C_Interface
    !===========================================================================
    !===========================================================================
    function ParMETIS_V3_PartMeshKway_C_Interface(elmdist, eptr, eind, elmwgt, &
                                                  wgtflag, numflag, ncon, &
                                                  ncommonnodes, nparts, &
                                                  tpwgts, ubvec, options, &
                                                  edgecut, part, comm) &
                                         bind(c,name="ParMETIS_V3_PartMeshKway")
      use, intrinsic :: iso_c_binding, only : c_int
      import :: idx_t,real_t
      integer(c_int) :: ParMETIS_V3_PartMeshKway_C_Interface
      !
      integer(idx_t) :: elmdist(*)
      integer(idx_t) :: eptr(*)
      integer(idx_t) :: eind(*)
      integer(idx_t) :: elmwgt(*)
      integer(idx_t) :: wgtflag
      integer(idx_t) :: numflag
      integer(idx_t) :: ncon
      integer(idx_t) :: ncommonnodes
      integer(idx_t) :: nparts
      real(real_t)   :: tpwgts(*)
      real(real_t)   :: ubvec(*)
      integer(idx_t) :: options(*)
      integer(idx_t) :: edgecut
      integer(idx_t) :: part(*)
      integer(c_int) :: comm
      !
    end function ParMETIS_V3_PartMeshKway_C_Interface
    !===========================================================================
    !===========================================================================
    function ParMETIS_V3_NodeND_C_Interface(vtxdist, xadj, adjncy, numflag, &
                                            options, order, sizes, comm) &
                                       bind(c,name="ParMETIS_V3_NodeND")
      use, intrinsic :: iso_c_binding, only : c_int
      import :: idx_t
      integer(c_int) :: ParMETIS_V3_NodeND_C_Interface
      !
      integer(idx_t) :: vtxdist(*)
      integer(idx_t) :: xadj(*)
      integer(idx_t) :: adjncy(*)
      integer(idx_t) :: numflag
      integer(idx_t) :: options(*)
      integer(idx_t) :: order(*)
      integer(idx_t) :: sizes(*)
      integer(c_int) :: comm
      !
    end function ParMETIS_V3_NodeND_C_Interface
    !===========================================================================
    !===========================================================================
    function ParMETIS_V32_NodeND_C_Interface(vtxdist, xadj, adjncy, vwgt, &
                                             numflag, mtype, rtype, p_nseps, &
                                             s_nseps, ubfrac, seed, dbglvl, &
                                             order, sizes, comm) &
                                        bind(c,name="ParMETIS_V32_NodeND")
      use, intrinsic :: iso_c_binding, only : c_int
      import :: idx_t,real_t
      integer(c_int) :: ParMETIS_V32_NodeND_C_Interface
      !
      integer(idx_t) :: vtxdist(*)
      integer(idx_t) :: xadj(*)
      integer(idx_t) :: adjncy(*)
      integer(idx_t) :: vwgt(*)
      integer(idx_t) :: numflag
      integer(idx_t) :: mtype
      integer(idx_t) :: rtype
      integer(idx_t) :: p_nseps
      integer(idx_t) :: s_nseps
      real(real_t)   :: ubfrac
      integer(idx_t) :: seed
      integer(idx_t) :: dbglvl
      integer(idx_t) :: order(*)
      integer(idx_t) :: sizes(*)
      integer(c_int) :: comm
      !
    end function ParMETIS_V32_NodeND_C_Interface
    !===========================================================================
    !===========================================================================
    function ParMETIS_SerialNodeND_C_Interface(vtxdist, xadj, adjncy, numflag, &
                                               options, order, sizes, comm) &
                                          bind(c,name="ParMETIS_SerialNodeND")
      use, intrinsic :: iso_c_binding, only : c_int
      import :: idx_t
      integer(c_int) :: ParMETIS_SerialNodeND_C_Interface
      !
      integer(idx_t) :: vtxdist(*)
      integer(idx_t) :: xadj(*)
      integer(idx_t) :: adjncy(*)
      integer(idx_t) :: numflag
      integer(idx_t) :: options(*)
      integer(idx_t) :: order(*)
      integer(idx_t) :: sizes(*)
      integer(c_int) :: comm
      !
    end function ParMETIS_SerialNodeND_C_Interface
    !===========================================================================
    !===========================================================================
  end interface
  !
  private :: ParMETIS_V3_PartKway_C_Interface
  private :: ParMETIS_V3_PartGeomKway_C_Interface
  private :: ParMETIS_V3_PartGeom_C_Interface
  private :: ParMETIS_V3_RefineKway_C_Interface
  private :: ParMETIS_V3_AdaptiveRepart_C_Interface
  private :: ParMETIS_V3_Mesh2Dual_C_Interface
  private :: ParMETIS_V3_PartMeshKway_C_Interface
  private :: ParMETIS_V3_NodeND_C_Interface
  private :: ParMETIS_V32_NodeND_C_Interface
  private :: ParMETIS_SerialNodeND_C_Interface
  !
  !
  ! pmoptype_et
  !
  enum, bind(c)
    enumerator :: PARMETIS_OP_KMETIS
    enumerator :: PARMETIS_OP_GKMETIS
    enumerator :: PARMETIS_OP_GMETIS
    enumerator :: PARMETIS_OP_RMETIS
    enumerator :: PARMETIS_OP_AMETIS
    enumerator :: PARMETIS_OP_OMETIS
    enumerator :: PARMETIS_OP_M2DUAL
    enumerator :: PARMETIS_OP_MKMETIS
  end enum
  !
  integer(idx_t), parameter :: PARMETIS_MTYPE_LOCAL  = 1_idx_t
  integer(idx_t), parameter :: PARMETIS_MTYPE_GLOBAL = 2_idx_t
  !
  integer(idx_t), parameter :: PARMETIS_SRTYPE_GREEDY = 1_idx_t
  integer(idx_t), parameter :: PARMETIS_SRTYPE_2PHASE = 2_idx_t
  !
  integer(idx_t), parameter :: PARMETIS_PSR_COUPLED   = 1_idx_t
  integer(idx_t), parameter :: PARMETIS_PSR_UNCOUPLED = 2_idx_t
  !
  ! Debugging ParMETIS is achieved by setting
  ! PARMETIS_DEBUG_FLAG = 1
  ! PARMETIS_DEBUG_VAL = [any sum of the following DBGLVL parameters]
  !
  integer(idx_t), parameter :: PARMETIS_DEBUG_FLAG = 0_idx_t
 !integer(idx_t), parameter :: PARMETIS_DEBUG_FLAG = 1_idx_t
  !
  integer(idx_t), parameter :: PARMETIS_DBGLVL_TIME = 1_idx_t
  integer(idx_t), parameter :: PARMETIS_DBGLVL_INFO = 2_idx_t
  integer(idx_t), parameter :: PARMETIS_DBGLVL_PROGRESS = 4_idx_t
  integer(idx_t), parameter :: PARMETIS_DBGLVL_REFINEINFO = 8_idx_t
  integer(idx_t), parameter :: PARMETIS_DBGLVL_MATCHINFO = 16_idx_t
  integer(idx_t), parameter :: PARMETIS_DBGLVL_RMOVEINFO = 32_idx_t
  integer(idx_t), parameter :: PARMETIS_DBGLVL_REMAP = 64_idx_t
  integer(idx_t), parameter :: PARMETIS_DBGLVL_TWOHOP = 128_idx_t
  integer(idx_t), parameter :: PARMETIS_DBGLVL_DROPEDGES = 256_idx_t
  integer(idx_t), parameter :: PARMETIS_DBGLVL_FAST = 512_idx_t
  integer(idx_t), parameter :: PARMETIS_DBGLVL_ONDISK = 1024_idx_t
  !
  integer(idx_t), parameter :: PARMETIS_DEBUG_VAL = 0_idx_t
  !
  ! CREATE_TPWGTS: Logical constant determining whether or not tpwgts and ubvec
  !                are explicitly allocated before sent as an argument to the
  !                various parmetis functions. Leaving them deallocated should
  !                be equivalent to passing a NULL argument which supposedly
  !                tells parmetis to create equal partitions.
  !
  logical, parameter :: create_tpwgts = .false.
 !logical, parameter :: create_tpwgts = .true.
  !
  ! Module Scalars
  !
  integer(int_mpi), save :: my_pm_rank
  integer(int_mpi), save :: num_pm_ranks
  integer(int_mpi), save :: mpierr
  !
  ! Module Allocatable Arrays
  !
  integer(idx_t), allocatable :: options(:)
  !
contains
!
!###############################################################################
!
subroutine ParMETIS_PartKway(vtxdist,xadj,adjncy,vwgt,adjwgt, &
                             wgtflag,numflag,ncon,nparts, &
                             tpwgts,ubvec,edgecut,part,comm)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int
  !
  !.. Formal Arguments ..
  integer(idx_t),    intent(in) :: vtxdist(:)
  integer(idx_t),    intent(in) :: xadj(:)
  integer(idx_t),    intent(in) :: adjncy(:)
  integer(idx_t),    intent(in) :: vwgt(:)
  integer(idx_t),    intent(in) :: adjwgt(:)
  integer(idx_t),    intent(in) :: wgtflag
  integer(idx_t),    intent(in) :: numflag
  integer(idx_t),    intent(in) :: ncon
  integer(idx_t),    intent(in) :: nparts
  real(real_t),      intent(in) :: tpwgts(:)
  real(real_t),      intent(in) :: ubvec(:)
  integer(idx_t), intent(inout) :: edgecut
  integer(idx_t), intent(inout) :: part(:)
  _MPI_COMM_TYPE_,   intent(in) :: comm
  !
  !.. Local Scalars ..
  integer(c_int) :: ierr
  !
continue
#ifdef ENABLE_PARMETIS
  !
  options = [PARMETIS_DEBUG_FLAG,PARMETIS_DEBUG_VAL,0_idx_t] ! F03 AUTO-REALLOC
  !
  ierr = ParMETIS_V3_PartKway_C_Interface(vtxdist, &
                                          xadj,adjncy,vwgt,adjwgt, &
                                          wgtflag,numflag, &
                                          ncon,nparts,tpwgts,ubvec, &
                                          options,edgecut,part, &
                                          _MPI_TYPE_VAL_(comm))
  !
  if (allocated(options)) deallocate ( options )
  !
#endif
end subroutine ParMETIS_PartKway
!
!###############################################################################
!
subroutine ParMETIS_PartGeomKway(vtxdist,xadj,adjncy,vwgt,adjwgt, &
                                 wgtflag,numflag,ndims,xyz,ncon,nparts, &
                                 tpwgts,ubvec,edgecut,part,comm)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int
  !
  !.. Formal Arguments ..
  integer(idx_t),    intent(in) :: vtxdist(:)
  integer(idx_t),    intent(in) :: xadj(:)
  integer(idx_t),    intent(in) :: adjncy(:)
  integer(idx_t),    intent(in) :: vwgt(:)
  integer(idx_t),    intent(in) :: adjwgt(:)
  integer(idx_t),    intent(in) :: wgtflag
  integer(idx_t),    intent(in) :: numflag
  integer(idx_t),    intent(in) :: ndims
  real(real_t),      intent(in) :: xyz(:)
  integer(idx_t),    intent(in) :: ncon
  integer(idx_t),    intent(in) :: nparts
  real(real_t),      intent(in) :: tpwgts(:)
  real(real_t),      intent(in) :: ubvec(:)
  integer(idx_t), intent(inout) :: edgecut
  integer(idx_t), intent(inout) :: part(:)
  _MPI_COMM_TYPE_,   intent(in) :: comm
  !
  !.. Local Scalars ..
  integer(c_int) :: ierr
  !
continue
#ifdef ENABLE_PARMETIS
  !
  options = [PARMETIS_DEBUG_FLAG,PARMETIS_DEBUG_VAL,0_idx_t] ! F03 AUTO-REALLOC
  !
  ierr = ParMETIS_V3_PartGeomKway_C_Interface(vtxdist, &
                                              xadj,adjncy,vwgt,adjwgt, &
                                              wgtflag,numflag,ndims,xyz, &
                                              ncon,nparts,tpwgts,ubvec, &
                                              options,edgecut,part, &
                                              _MPI_TYPE_VAL_(comm))
  !
  if (allocated(options)) deallocate ( options )
  !
#endif
end subroutine ParMETIS_PartGeomKway
!
!###############################################################################
!
subroutine create_vtxdist(my_ncell,vtxdist,comm) 
  !
  !.. Use Statements .
  use _MPI_MODULE_, only : mpi_comm_size
  use _MPI_MODULE_, only : mpi_comm_rank
  use _MPI_MODULE_, only : mpi_allgather
  !
  !.. Formal Arguments ..
  integer(idx_t),               intent(in) :: my_ncell
  integer(idx_t), allocatable, intent(out) :: vtxdist(:)
  _MPI_COMM_TYPE_,              intent(in) :: comm
  !
  !.. Local Scalars ..
  integer :: n
  !
  !.. Local Allocatable Arrays ..
  integer(idx_t), allocatable :: cell_counts(:)
  !
continue
  !
  ! Get the number of parmetis ranks and the rank of this MPI process
  ! for this MPI communicator.
  !
  call mpi_comm_size(comm,num_pm_ranks,mpierr)
  !
  call mpi_comm_rank(comm,my_pm_rank,mpierr)
  !
  ! Allocate the temporary cell_counts array
  !
  allocate ( cell_counts(1:num_pm_ranks) , source=0_idx_t )
  !
  ! Use mpi_allgather to combine the individual cell counts
  !
  call mpi_allgather(my_ncell,1_int_mpi,pm_mpi_inttyp, &
                     cell_counts,1_int_mpi,pm_mpi_inttyp, &
                     comm,mpierr)
  !
  ! Allocate the vtxdist array
  !
  allocate ( vtxdist(1:num_pm_ranks+1) , source=0_idx_t )
  !
  ! Create the vtxdist array so it contains the global range of cells 
  ! local to each MPI process involved in the parallel partitioning:
  !   MPI process i contains the cells vtxdist(i+1) to vtxdist(i+2)-1.
  !
  vtxdist(1) = 1_idx_t
  do n = 1,num_pm_ranks
    vtxdist(n+1) = vtxdist(n) + cell_counts(n)
  end do
  !
  if (allocated(cell_counts)) deallocate ( cell_counts )
  !
end subroutine create_vtxdist
!
!###############################################################################
!
subroutine recombine_global_epart(vtxdist,part,epart,comm) 
  !
  !.. Use Statements .
  use _MPI_MODULE_, only : MPI_INTEGER
  use _MPI_MODULE_, only : MPI_IN_PLACE
  use _MPI_MODULE_, only : mpi_allgather
  !
  !.. Formal Arguments ..
  integer(idx_t),    intent(in) :: vtxdist(:)
  integer(idx_t),    intent(in) :: part(:)
  integer(idx_t), intent(inout) :: epart(:)
  _MPI_COMM_TYPE_,   intent(in) :: comm
  !
  !.. Local Scalars ..
  integer :: n,n1,n2
  integer(int_mpi) :: this_ncell,max_ncell
  !
  !.. Local Allocatable Arrays ..
  integer(idx_t), allocatable :: padded_epart(:)
  !
continue
  !
  ! Get the maximum number of cells on any given MPI process involved
  ! with the parallel partitioning.
  !
  max_ncell = 0_int_mpi
  do n = 2,size(vtxdist)
    this_ncell = int( vtxdist(n)-vtxdist(n-1) , kind=kind(max_ncell) )
    max_ncell = max( max_ncell , this_ncell )
  end do
  !
  ! Allocate a padded version of the global epart array where every MPI
  ! process has the same number of cells equal to the max number of cells.
  !
  allocate ( padded_epart(1:num_pm_ranks*max_ncell) , source=-99_idx_t )
  !
  ! Get the index for the first cell of this MPI process within the padded
  ! epart array, and the index for the last cell without padding.
  !
  n1 = my_pm_rank*max_ncell + 1
  n2 = n1-1 + size(part)
  !
  ! Copy the local piece of the global partitioning into the padded version
  ! of the global epart array.
  !
  padded_epart(n1:n2) = part(:)
  !
  ! Use mpi_allgather to combine the individual padded_epart arrays
  ! to create a global version.
  !
  call mpi_allgather(MPI_IN_PLACE,0_int_mpi,MPI_INTEGER, &
                     padded_epart,max_ncell,pm_mpi_inttyp, &
                     comm,mpierr)
  !
  ! Count the number of non-negative values in padded_epart. This should be
  ! equal to the total number of global cells for this simulation.
  !
  n = count( padded_epart(:) >= 0 )
  !
  ! If the number of non-negative values is correct, copy only the non-negative
  ! values into the epart array; otherwise, assign all cells the invalid
  ! partition number of -99.
  !
  if (size(epart) == n) then
    epart(:) = pack( padded_epart(:) , mask=padded_epart(:)>=0 )
  else
    epart(:) = -99_idx_t
  end if
  !
  if (allocated(padded_epart)) deallocate ( padded_epart )
  !
end subroutine recombine_global_epart
!
!###############################################################################
!
subroutine create_partition_weights(nparts,ncon,tpwgts,ubvec) 
  !
  !.. Formal Arguments ..
  integer(idx_t),             intent(in) :: nparts
  integer(idx_t),             intent(in) :: ncon
  real(real_t), allocatable, intent(out) :: tpwgts(:)
  real(real_t), allocatable, intent(out) :: ubvec(:)
  !
continue
  !
  if (create_tpwgts) then
    !
    allocate ( tpwgts(1:ncon*nparts) )
    !
    tpwgts(:) = 1.0_real_t / real(nparts,kind=real_t)
    !
    allocate ( ubvec(1:ncon) )
    !
    ! The serial METIS manual (v 5.1.0; March 30, 2013; section 5.8 on page 26)
    ! recommends the values
    !         1.001   if ncon==1
    !         1.01    otherwise
   !ubvec(:) = merge( 1.001_real_t , 1.01_real_t , ncon==1 )
    !
    ! The ParMETIS manual (v 4.0; March 30, 2013; section 5.1 on page 16)
    ! recommends quote:
    !   "A value of 1.05 for each of the ncon weights is recommended."
    ubvec(:) = 1.05_real_t
    !
  end if
  !
end subroutine create_partition_weights
!
!###############################################################################
!
pure &
function get_weight_flag(vwgt,adjwgt) result(return_value)
  !
  !.. Formal Arguments ..
  integer(idx_t), allocatable, intent(in) :: vwgt(:)
  integer(idx_t), allocatable, intent(in) :: adjwgt(:)
  !
  !.. Function Result ..
  integer(idx_t) :: return_value
  !
continue
  !
  ! The argument wgtflag for the parmetis functions tells parmetis what type of
  ! weighting to use to compute the partitions. The possible values are
  !
  !   wgtflag = 0   : no weighting
  !                   (vwgt is NOT allocated  ; adjwgt is NOT allocated)
  !   wgtflag = 1   : only use the edge weights
  !                   (vwgt is NOT allocated  ; adjwgt IS     allocated)
  !   wgtflag = 2   : only use the vertex weights
  !                   (vwgt IS     allocated  ; adjwgt is NOT allocated)
  !   wgtflag = 3   : use both the edge and vertex weights
  !                   (vwgt IS     allocated  ; adjwgt IS     allocated)
  !
  ! NOTE: The terms vertex and edge refer to components of the METIS dual graph
  !       which respectively correspond to grid cells and grid faces.
  !
  return_value = 0_idx_t
  if (allocated(adjwgt)) return_value = return_value + 1_idx_t
  if (allocated(vwgt))   return_value = return_value + 2_idx_t
  !
end function get_weight_flag
!
!###############################################################################
!
pure &
function get_number_of_vertex_weights(xadj,vwgt) result(return_value)
  !
  !.. Formal Arguments ..
  integer(idx_t),              intent(in) :: xadj(:)
  integer(idx_t), allocatable, intent(in) :: vwgt(:)
  !
  !.. Function Result ..
  integer(idx_t) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,n
  !
continue
  !
  n = size(xadj)-1 ! number of vertices in the dual graph (i.e., grid cells)
  !
  if (allocated(vwgt)) then
    ! the size of vwgt should be an integer multiple of the number of vertices
    ! in the dual graph, but make sure it is >= n
    i = max( size(vwgt) , n )
  else
    i = n ! give i an initial value
  end if
  !
  ! Make sure the returned value for ncon is >= 1
  !
  return_value = int( max( 1 , i/n ) , kind=kind(return_value) )
  !
end function get_number_of_vertex_weights
!
!###############################################################################
!
subroutine ParMETIS_PartGeom(vtxdist,ndims,xyz,part,comm)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int
  !
  !.. Formal Arguments ..
  integer(idx_t),    intent(in) :: vtxdist(:)
  integer(idx_t),    intent(in) :: ndims
  real(real_t),      intent(in) :: xyz(:)
  integer(idx_t), intent(inout) :: part(:)
  _MPI_COMM_TYPE_,   intent(in) :: comm
  !
  !.. Local Scalars ..
  integer(c_int) :: ierr
  !
continue
#ifdef ENABLE_PARMETIS
  !
  ierr = ParMETIS_V3_PartGeom_C_Interface(vtxdist,ndims,xyz,part, &
                                          _MPI_TYPE_VAL_(comm))
  !
#endif
end subroutine ParMETIS_PartGeom
!
!###############################################################################
!
subroutine ParMETIS_RefineKway(vtxdist,xadj,adjncy,vwgt,adjwgt, &
                               wgtflag,numflag,ncon,nparts,tpwgts, &
                               ubvec,edgecut,part,comm)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int
  !
  !.. Formal Arguments ..
  integer(idx_t),    intent(in) :: vtxdist(:)
  integer(idx_t),    intent(in) :: xadj(:)
  integer(idx_t),    intent(in) :: adjncy(:)
  integer(idx_t),    intent(in) :: vwgt(:)
  integer(idx_t),    intent(in) :: adjwgt(:)
  integer(idx_t),    intent(in) :: wgtflag
  integer(idx_t),    intent(in) :: numflag
  integer(idx_t),    intent(in) :: ncon
  integer(idx_t),    intent(in) :: nparts
  real(real_t),      intent(in) :: tpwgts(:)
  real(real_t),      intent(in) :: ubvec(:)
  integer(idx_t), intent(inout) :: edgecut
  integer(idx_t), intent(inout) :: part(:)
  _MPI_COMM_TYPE_,   intent(in) :: comm
  !
  !.. Local Scalars ..
  integer(c_int) :: ierr
  !
continue
#ifdef ENABLE_PARMETIS
  !
  options = [PARMETIS_DEBUG_FLAG,PARMETIS_DEBUG_VAL, & ! F2003 AUTO-REALLOC
             0_idx_t,0_idx_t]
  !
  ierr = ParMETIS_V3_RefineKway_C_Interface(vtxdist,xadj,adjncy,vwgt,adjwgt, &
                               wgtflag,numflag,ncon,nparts,tpwgts, &
                               ubvec,options,edgecut,part, &
                               _MPI_TYPE_VAL_(comm))
  !
  if (allocated(options)) deallocate ( options )
  !
#endif
end subroutine ParMETIS_RefineKway
!
!###############################################################################
!
subroutine ParMETIS_AdaptiveRepart(vtxdist,xadj,adjncy,vwgt,vsize,adjwgt, &
                                   wgtflag,numflag,ncon,nparts,tpwgts,ubvec, &
                                   ipc2redist,edgecut,part,comm)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int
  !
  !.. Formal Arguments ..
  integer(idx_t),    intent(in) :: vtxdist(:)
  integer(idx_t),    intent(in) :: xadj(:)
  integer(idx_t),    intent(in) :: adjncy(:)
  integer(idx_t),    intent(in) :: vwgt(:)
  integer(idx_t),    intent(in) :: vsize(:)
  integer(idx_t),    intent(in) :: adjwgt(:)
  integer(idx_t),    intent(in) :: wgtflag
  integer(idx_t),    intent(in) :: numflag
  integer(idx_t),    intent(in) :: ncon
  integer(idx_t),    intent(in) :: nparts
  real(real_t),      intent(in) :: tpwgts(:)
  real(real_t),      intent(in) :: ubvec(:)
  real(real_t),      intent(in) :: ipc2redist
  integer(idx_t), intent(inout) :: edgecut
  integer(idx_t), intent(inout) :: part(:)
  _MPI_COMM_TYPE_,   intent(in) :: comm
  !
  !.. Local Scalars ..
  integer(c_int) :: ierr
  !
continue
#ifdef ENABLE_PARMETIS
  !
  options = [PARMETIS_DEBUG_FLAG,PARMETIS_DEBUG_VAL, & ! F2003 AUTO-REALLOC
             0_idx_t,0_idx_t]
  !
  ierr = ParMETIS_V3_AdaptiveRepart_C_Interface(vtxdist,xadj,adjncy,vwgt, &
                                                vsize,adjwgt,wgtflag,numflag, &
                                                ncon,nparts,tpwgts,ubvec, &
                                                ipc2redist,options, &
                                                edgecut,part, &
                                                _MPI_TYPE_VAL_(comm))
  !
  if (allocated(options)) deallocate ( options )
  !
#endif
end subroutine ParMETIS_AdaptiveRepart
!
!###############################################################################
!
subroutine ParMETIS_Mesh2Dual(elmdist,eptr,eind,numflag,ncommonnodes, &
                              xadj,adjncy,comm)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int,c_ptr
  use, intrinsic :: iso_c_binding, only : c_f_pointer
  !
  use metis_mod, only : METIS_Free
  !
  !.. Formal Arguments ..
  integer(idx_t),                 intent(in) :: elmdist(:)
  integer(idx_t),                 intent(in) :: eptr(:)
  integer(idx_t),                 intent(in) :: eind(:)
  integer(idx_t),                 intent(in) :: numflag
  integer(idx_t),                 intent(in) :: ncommonnodes
  integer(idx_t), allocatable, intent(inout) :: xadj(:)
  integer(idx_t), allocatable, intent(inout) :: adjncy(:)
  _MPI_COMM_TYPE_,                intent(in) :: comm
  !
  !.. Local Scalars ..
  integer :: n,ne
  integer(c_int) :: ierr
  type(c_ptr) :: r_xadj
  type(c_ptr) :: r_adjncy
  !
  !.. Local Pointers ..
  integer(idx_t), pointer :: ptr(:) => null()
  !
continue
  !
  if (allocated(xadj)) deallocate ( xadj )
  if (allocated(adjncy)) deallocate ( adjncy )
  !
  ne = int(elmdist(size(elmdist)),kind=kind(ne))
  n = ne+1
  !
#ifdef ENABLE_PARMETIS
  !
  ierr = ParMETIS_V3_Mesh2Dual_C_Interface(elmdist,eptr,eind,numflag, &
                                           ncommonnodes,r_xadj,r_adjncy, &
                                           _MPI_TYPE_VAL_(comm))
  !
  call c_f_pointer( cptr=r_xadj , fptr=ptr , shape=[n] )
  !
  allocate ( xadj(1:n) )
  xadj(1:n) = ptr(1:n)
  ptr => null()
  !
  n = xadj(ne+1)
  call c_f_pointer( cptr=r_adjncy , fptr=ptr , shape=[n] )
  !
  allocate ( adjncy(1:n) )
  adjncy(1:n) = ptr(1:n)
  ptr => null()
  !
  ierr = METIS_Free(r_xadj)
  ierr = METIS_Free(r_adjncy)
  !
#else
  !
  allocate ( xadj(1:n) , source=0_idx_t )
  allocate ( adjncy(1:1) , source=0_idx_t )
  !
#endif
end subroutine ParMETIS_Mesh2Dual
!
!###############################################################################
!
subroutine ParMETIS_PartMeshKway(elmdist,eptr,eind,elmwgt, &
                                 wgtflag,numflag,ncon,ncommonnodes,nparts, &
                                 tpwgts,ubvec,edgecut,part,comm)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int
  !
  !.. Formal Arguments ..
  integer(idx_t),    intent(in) :: elmdist(:)
  integer(idx_t),    intent(in) :: eptr(:)
  integer(idx_t),    intent(in) :: eind(:)
  integer(idx_t),    intent(in) :: elmwgt(:)
  integer(idx_t),    intent(in) :: wgtflag
  integer(idx_t),    intent(in) :: numflag
  integer(idx_t),    intent(in) :: ncon
  integer(idx_t),    intent(in) :: ncommonnodes
  integer(idx_t),    intent(in) :: nparts
  real(real_t),      intent(in) :: tpwgts(:)
  real(real_t),      intent(in) :: ubvec(:)
  integer(idx_t), intent(inout) :: edgecut
  integer(idx_t), intent(inout) :: part(:)
  _MPI_COMM_TYPE_,   intent(in) :: comm
  !
  !.. Local Scalars ..
  integer(c_int) :: ierr
  !
continue
#ifdef ENABLE_PARMETIS
  !
  options = [PARMETIS_DEBUG_FLAG,PARMETIS_DEBUG_VAL,0_idx_t] ! F03 AUTO-REALLOC
  !
  ierr = ParMETIS_V3_PartMeshKway_C_Interface(elmdist,eptr,eind,elmwgt, &
                                              wgtflag,numflag,ncon, &
                                              ncommonnodes,nparts, &
                                              tpwgts,ubvec,options, &
                                              edgecut,part, &
                                              _MPI_TYPE_VAL_(comm))
  !
  if (allocated(options)) deallocate ( options )
  !
#endif
end subroutine ParMETIS_PartMeshKway
!
!###############################################################################
!
subroutine ParMETIS_NodeND(vtxdist,xadj,adjncy,numflag,order,sizes,comm)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int
  !
  !.. Formal Arguments ..
  integer(idx_t),    intent(in) :: vtxdist(:)
  integer(idx_t),    intent(in) :: xadj(:)
  integer(idx_t),    intent(in) :: adjncy(:)
  integer(idx_t),    intent(in) :: numflag
  integer(idx_t), intent(inout) :: order(:)
  integer(idx_t), intent(inout) :: sizes(:)
  _MPI_COMM_TYPE_,   intent(in) :: comm
  !
  !.. Local Scalars ..
  integer(c_int) :: ierr
  !
continue
#ifdef ENABLE_PARMETIS
  !
  options = [PARMETIS_DEBUG_FLAG,PARMETIS_DEBUG_VAL,0_idx_t] ! F03 AUTO-REALLOC
  !
  ierr = ParMETIS_V3_NodeND_C_Interface(vtxdist,xadj,adjncy,numflag, &
                                        options,order,sizes, &
                                        _MPI_TYPE_VAL_(comm))
  !
  if (allocated(options)) deallocate ( options )
  !
#endif
end subroutine ParMETIS_NodeND
!
!###############################################################################
!
subroutine ParMETIS_V32_NodeND(vtxdist,xadj,adjncy,vwgt,numflag,mtype,rtype, &
                               p_nseps,s_nseps,ubfrac,order,sizes,comm)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int
  !
  !.. Formal Arguments ..
  integer(idx_t),    intent(in) :: vtxdist(:)
  integer(idx_t),    intent(in) :: xadj(:)
  integer(idx_t),    intent(in) :: adjncy(:)
  integer(idx_t),    intent(in) :: vwgt(:)
  integer(idx_t),    intent(in) :: numflag
  integer(idx_t),    intent(in) :: mtype
  integer(idx_t),    intent(in) :: rtype
  integer(idx_t),    intent(in) :: p_nseps
  integer(idx_t),    intent(in) :: s_nseps
  real(real_t),      intent(in) :: ubfrac
  integer(idx_t), intent(inout) :: order(:)
  integer(idx_t), intent(inout) :: sizes(:)
  _MPI_COMM_TYPE_,   intent(in) :: comm
  !
  !.. Local Scalars ..
  integer(c_int) :: ierr
  integer(idx_t) :: seed,dbglvl
  !
continue
#ifdef ENABLE_PARMETIS
  !
  seed = 0_idx_t
  dbglvl = merge(PARMETIS_DEBUG_VAL,0_idx_t,PARMETIS_DEBUG_FLAG > 0) 
  !
  ierr = ParMETIS_V32_NodeND_C_Interface(vtxdist,xadj,adjncy,vwgt,numflag, &
                                         mtype,rtype,p_nseps,s_nseps, &
                                         ubfrac,seed,dbglvl,order,sizes, &
                                         _MPI_TYPE_VAL_(comm))
  !
#endif
end subroutine ParMETIS_V32_NodeND
!
!###############################################################################
!
subroutine ParMETIS_SerialNodeND(vtxdist,xadj,adjncy,numflag,order,sizes,comm)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int
  !
  !.. Formal Arguments ..
  integer(idx_t),    intent(in) :: vtxdist(:)
  integer(idx_t),    intent(in) :: xadj(:)
  integer(idx_t),    intent(in) :: adjncy(:)
  integer(idx_t),    intent(in) :: numflag
  integer(idx_t), intent(inout) :: order(:)
  integer(idx_t), intent(inout) :: sizes(:)
  _MPI_COMM_TYPE_,   intent(in) :: comm
  !
  !.. Local Scalars ..
  integer(c_int) :: ierr
  !
continue
#ifdef ENABLE_PARMETIS
  !
  options = [PARMETIS_DEBUG_FLAG,PARMETIS_DEBUG_VAL,0_idx_t] ! F03 AUTO-REALLOC
  !
  ierr = ParMETIS_SerialNodeND_C_Interface(vtxdist,xadj,adjncy,numflag, &
                                           options,order,sizes, &
                                           _MPI_TYPE_VAL_(comm))
  !
  if (allocated(options)) deallocate ( options )
  !
#endif
end subroutine ParMETIS_SerialNodeND
!
!###############################################################################
!
end module parmetis_mod
