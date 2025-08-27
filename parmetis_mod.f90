module parmetis_mod
  !
#include <mpi_defs.h>
  !
  use, intrinsic :: iso_c_binding, only : c_float,c_double
  use, intrinsic :: iso_c_binding, only : c_int32_t,c_int64_t
  !
#ifdef USE_MPI_F08
  use _MPI_MODULE_, only : MPI_COMM
#else
  use module_kind_types, only : INT_MPI
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
  integer, parameter :: PARMETIS_DBGLVL_TIME = 1
  integer, parameter :: PARMETIS_DBGLVL_INFO = 2
  integer, parameter :: PARMETIS_DBGLVL_PROGRESS = 4
  integer, parameter :: PARMETIS_DBGLVL_REFINEINFO = 8
  integer, parameter :: PARMETIS_DBGLVL_MATCHINFO = 16
  integer, parameter :: PARMETIS_DBGLVL_RMOVEINFO = 32
  integer, parameter :: PARMETIS_DBGLVL_REMAP = 64
  integer, parameter :: PARMETIS_DBGLVL_TWOHOP = 128
  integer, parameter :: PARMETIS_DBGLVL_DROPEDGES = 256
  integer, parameter :: PARMETIS_DBGLVL_FAST = 512
  integer, parameter :: PARMETIS_DBGLVL_ONDISK = 1024
  !
  integer(idx_t), parameter :: PARMETIS_DEBUG_VAL = 0_idx_t
  !
  !
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
  ierr = ParMETIS_V3_PartKway_C_Interface(vtxdist,xadj,adjncy,vwgt,adjwgt, &
                                          wgtflag,numflag,ncon,nparts, &
                                          tpwgts,ubvec,options,edgecut,part, &
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
  ierr = ParMETIS_V3_PartGeomKway_C_Interface(vtxdist,xadj,adjncy,vwgt,adjwgt, &
                                 wgtflag,numflag,ndims,xyz,ncon,nparts, &
                                 tpwgts,ubvec,options,edgecut,part, &
                                 _MPI_TYPE_VAL_(comm))
  !
  if (allocated(options)) deallocate ( options )
  !
#endif
end subroutine ParMETIS_PartGeomKway
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
