module metis5_mod
  !
  use, intrinsic :: iso_c_binding, only : c_float,c_double
  use, intrinsic :: iso_c_binding, only : c_int32_t,c_int64_t
  !
  implicit none
  !
  public
  !
  private :: c_float, c_double, c_int32_t, c_int64_t
  !
#ifdef METIS_64_BIT
  integer, parameter :: idx_t = c_int64_t
  integer, parameter :: real_t = c_double
#else
  integer, parameter :: idx_t = c_int32_t
  integer, parameter :: real_t = c_float
#endif
  !
  interface
    !===========================================================================
    !===========================================================================
    function METIS_PartGraphRecursive_C_Interface(nvtxs, ncon, xadj, adjncy, &
                                                  vwgt, vsize, adjwgt, nparts, &
                                                  tpwgts, ubvec, options, &
                                                  objval, part) &
                                         bind(c,name="METIS_PartGraphRecursive")
      use, intrinsic :: iso_c_binding, only : c_int
      import :: idx_t, real_t
      integer(c_int) :: METIS_PartGraphRecursive_C_Interface
      integer(idx_t) :: nvtxs
      integer(idx_t) :: ncon
      integer(idx_t) :: xadj(*)
      integer(idx_t) :: adjncy(*)
      integer(idx_t) :: vwgt(*)
      integer(idx_t) :: vsize(*)
      integer(idx_t) :: adjwgt(*)
      integer(idx_t) :: nparts
      real(real_t)   :: tpwgts(*)
      real(real_t)   :: ubvec(*)
      integer(idx_t) :: options(*)
      integer(idx_t) :: objval
      integer(idx_t) :: part(*)
    end function METIS_PartGraphRecursive_C_Interface
    !===========================================================================
    !===========================================================================
    function METIS_PartGraphKway_C_Interface(nvtxs, ncon, xadj, adjncy, vwgt, &
                                             vsize, adjwgt, nparts, tpwgts, &
                                             ubvec, options, objval, part) &
                                             bind(c,name="METIS_PartGraphKway")
      use, intrinsic :: iso_c_binding, only : c_int
      import :: idx_t, real_t
      integer(c_int) :: METIS_PartGraphKway_C_Interface
      integer(idx_t) :: nvtxs
      integer(idx_t) :: ncon
      integer(idx_t) :: xadj(*)
      integer(idx_t) :: adjncy(*)
      integer(idx_t) :: vwgt(*)
      integer(idx_t) :: vsize(*)
      integer(idx_t) :: adjwgt(*)
      integer(idx_t) :: nparts
      real(real_t)   :: tpwgts(*)
      real(real_t)   :: ubvec(*)
      integer(idx_t) :: options(*)
      integer(idx_t) :: objval
      integer(idx_t) :: part(*)
    end function METIS_PartGraphKway_C_Interface
    !===========================================================================
    !===========================================================================
    function METIS_PartMeshDual_C_Interface(ne, nn, eptr, eind, vwgt, vsize, &
                                            ncommon, nparts, tpwgts, options, &
                                            objval, epart, npart) &
                                            bind(c,name="METIS_PartMeshDual")
      use, intrinsic :: iso_c_binding, only : c_int
      import :: idx_t, real_t
      integer(c_int) :: METIS_PartMeshDual_C_Interface
      integer(idx_t) :: ne
      integer(idx_t) :: nn
      integer(idx_t) :: eptr(*)
      integer(idx_t) :: eind(*)
      integer(idx_t) :: vwgt(*)
      integer(idx_t) :: vsize(*)
      integer(idx_t) :: ncommon
      integer(idx_t) :: nparts
      real(real_t)   :: tpwgts(*)
      integer(idx_t) :: options(*)
      integer(idx_t) :: objval
      integer(idx_t) :: epart(*)
      integer(idx_t) :: npart(*)
    end function METIS_PartMeshDual_C_Interface
    !===========================================================================
    !===========================================================================
    function METIS_PartMeshNodal_C_Interface(ne, nn, eptr, eind, vwgt, vsize, &
                                             nparts, tpwgts, options, objval, &
                                             epart, npart) &
                                             bind(c,name="METIS_PartMeshNodal")
      use, intrinsic :: iso_c_binding, only : c_int
      import :: idx_t, real_t
      integer(c_int) :: METIS_PartMeshNodal_C_Interface
      integer(idx_t) :: ne
      integer(idx_t) :: nn
      integer(idx_t) :: eptr(*)
      integer(idx_t) :: eind(*)
      integer(idx_t) :: vwgt(*)
      integer(idx_t) :: vsize(*)
      integer(idx_t) :: nparts
      real(real_t)   :: tpwgts(*)
      integer(idx_t) :: options(*)
      integer(idx_t) :: objval
      integer(idx_t) :: epart(*)
      integer(idx_t) :: npart(*)
    end function METIS_PartMeshNodal_C_Interface
    !===========================================================================
    !===========================================================================
    function METIS_NodeND_C_Interface(nvtxs, xadj, adjncy, vwgt, &
                                      options, perm, iperm) &
                                      bind(c,name="METIS_NodeND")
      use, intrinsic :: iso_c_binding, only : c_int
      import :: idx_t
      integer(c_int) :: METIS_NodeND_C_Interface
      integer(idx_t) :: nvtxs
      integer(idx_t) :: xadj(*)
      integer(idx_t) :: adjncy(*)
      integer(idx_t) :: vwgt(*)
      integer(idx_t) :: options(*)
      integer(idx_t) :: perm(*)
      integer(idx_t) :: iperm(*)
    end function METIS_NodeND_C_Interface
    !===========================================================================
    !===========================================================================
    function METIS_MeshToDual_C_Interface(ne, nn, eptr, eind, ncommon, &
                                          numflag, r_xadj, r_adjncy) &
                                          bind(c,name="METIS_MeshToDual")
      use, intrinsic :: iso_c_binding, only : c_int,c_ptr
      import :: idx_t
      integer(c_int) :: METIS_MeshToDual_C_Interface
      integer(idx_t) :: ne
      integer(idx_t) :: nn
      integer(idx_t) :: eptr(*)
      integer(idx_t) :: eind(*)
      integer(idx_t) :: ncommon
      integer(idx_t) :: numflag
      type(c_ptr) :: r_xadj
      type(c_ptr) :: r_adjncy
    end function METIS_MeshToDual_C_Interface
    !===========================================================================
    !===========================================================================
    function METIS_MeshToNodal_C_Interface(ne, nn, eptr, eind, numflag, &
                                           r_xadj, r_adjncy) &
                                           bind(c,name="METIS_MeshToNodal")
      use, intrinsic :: iso_c_binding, only : c_int,c_ptr
      import :: idx_t
      integer(c_int) :: METIS_MeshToNodal_C_Interface
      integer(idx_t) :: ne
      integer(idx_t) :: nn
      integer(idx_t) :: eptr(*)
      integer(idx_t) :: eind(*)
      integer(idx_t) :: numflag
      type(c_ptr) :: r_xadj
      type(c_ptr) :: r_adjncy
    end function METIS_MeshToNodal_C_Interface
    !===========================================================================
    !===========================================================================
    function METIS_SetDefaultOptions_C_Interface(options) &
                                          bind(c,name="METIS_SetDefaultOptions")
      use, intrinsic :: iso_c_binding, only : c_int
      import :: idx_t
      integer(c_int) :: METIS_SetDefaultOptions_C_Interface
      integer(idx_t) :: options(*)
    end function METIS_SetDefaultOptions_C_Interface
    !===========================================================================
    !===========================================================================
    function METIS_Free(ptr) bind(c,name="METIS_Free")
      use, intrinsic :: iso_c_binding, only : c_int,c_ptr
      integer(c_int) :: METIS_Free
      type(c_ptr), value :: ptr
    end function METIS_Free
    !===========================================================================
    !===========================================================================
  end interface
  !
  private :: METIS_PartGraphRecursive_C_Interface
  private :: METIS_PartGraphKway_C_Interface
  private :: METIS_PartMeshDual_C_Interface
  private :: METIS_PartMeshNodal_C_Interface
  private :: METIS_NodeND_C_Interface
  private :: METIS_MeshToDual_C_Interface
  private :: METIS_MeshToNodal_C_Interface
  private :: METIS_SetDefaultOptions_C_Interface
  !
  !
  !
  integer, parameter :: METIS_NOPTIONS = 40
  !
  ! rstatus_et
  !
  enum, bind(c)
    enumerator :: METIS_OK = 1
    enumerator :: METIS_ERROR_INPUT = -2
    enumerator :: METIS_ERROR_MEMORY = -3
    enumerator :: METIS_ERROR = -4
  end enum
  !
  ! moptype_et
  !
  enum, bind(c)
    enumerator :: METIS_OP_PMETIS
    enumerator :: METIS_OP_KMETIS
    enumerator :: METIS_OP_OMETIS
  end enum
  !
  ! moptions_et
  !
  enum, bind(c)
    enumerator :: METIS_OPTION_PTYPE
    enumerator :: METIS_OPTION_OBJTYPE
    enumerator :: METIS_OPTION_CTYPE
    enumerator :: METIS_OPTION_IPTYPE
    enumerator :: METIS_OPTION_RTYPE
    enumerator :: METIS_OPTION_DBGLVL
    enumerator :: METIS_OPTION_NITER
    enumerator :: METIS_OPTION_NCUTS
    enumerator :: METIS_OPTION_SEED
    enumerator :: METIS_OPTION_NO2HOP
    enumerator :: METIS_OPTION_MINCONN
    enumerator :: METIS_OPTION_CONTIG
    enumerator :: METIS_OPTION_COMPRESS
    enumerator :: METIS_OPTION_CCORDER
    enumerator :: METIS_OPTION_PFACTOR
    enumerator :: METIS_OPTION_NSEPS
    enumerator :: METIS_OPTION_UFACTOR
    enumerator :: METIS_OPTION_NUMBERING
    enumerator :: METIS_OPTION_HELP
    enumerator :: METIS_OPTION_TPWGTS
    enumerator :: METIS_OPTION_NCOMMON
    enumerator :: METIS_OPTION_NOOUTPUT
    enumerator :: METIS_OPTION_BALANCE
    enumerator :: METIS_OPTION_GTYPE
    enumerator :: METIS_OPTION_UBVEC
  end enum
  !
  ! mptype_et
  !
  enum, bind(c)
    enumerator :: METIS_PTYPE_RB
    enumerator :: METIS_PTYPE_KWAY
  end enum
  !
  ! mgtype_et
  !
  enum, bind(c)
    enumerator :: METIS_GTYPE_DUAL
    enumerator :: METIS_GTYPE_NODAL
  end enum
  !
  ! mctype_et
  !
  enum, bind(c)
    enumerator :: METIS_CTYPE_RM
    enumerator :: METIS_CTYPE_SHEM
  end enum
  !
  ! miptype_et
  !
  enum, bind(c)
    enumerator :: METIS_IPTYPE_GROW
    enumerator :: METIS_IPTYPE_RANDOM
    enumerator :: METIS_IPTYPE_EDGE
    enumerator :: METIS_IPTYPE_NODE
    enumerator :: METIS_IPTYPE_METISRB
  end enum
  !
  ! mrtype_et
  !
  enum, bind(c)
    enumerator :: METIS_RTYPE_FM
    enumerator :: METIS_RTYPE_GREEDY
    enumerator :: METIS_RTYPE_SEP2SIDED
    enumerator :: METIS_RTYPE_SEP1SIDED
  end enum
  !
  ! mdbglvl_et
  !
  enum, bind(c)
    enumerator :: METIS_DBG_INFO = 1
    enumerator :: METIS_DBG_TIME = 2
    enumerator :: METIS_DBG_COARSEN = 4
    enumerator :: METIS_DBG_REFINE = 8
    enumerator :: METIS_DBG_IPART = 16
    enumerator :: METIS_DBG_MOVEINFO = 32
    enumerator :: METIS_DBG_SEPINFO = 64
    enumerator :: METIS_DBG_CONNINFO = 128
    enumerator :: METIS_DBG_CONTIGINFO = 256
    enumerator :: METIS_DBG_MEMORY = 2048
  end enum
  !
  ! mobjtype_et
  !
  enum, bind(c)
    enumerator :: METIS_OBJTYPE_CUT
    enumerator :: METIS_OBJTYPE_VOL
    enumerator :: METIS_OBJTYPE_NODE
  end enum
  !
contains
!
!###############################################################################
!
subroutine METIS_PartGraphRecursive(nvtxs,ncon,xadj,adjncy,vwgt,vsize,adjwgt, &
                                    nparts,tpwgts,ubvec,options,objval,part)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int
  !
  !.. Formal Arguments ..
  integer(idx_t),              intent(in) :: nvtxs
  integer(idx_t),              intent(in) :: ncon
  integer(idx_t),              intent(in) :: xadj(:)
  integer(idx_t),              intent(in) :: adjncy(:)
  integer(idx_t), allocatable, intent(in) :: vwgt(:)
  integer(idx_t), allocatable, intent(in) :: vsize(:)
  integer(idx_t), allocatable, intent(in) :: adjwgt(:)
  integer(idx_t),              intent(in) :: nparts
  real(real_t),   allocatable, intent(in) :: tpwgts(:)
  real(real_t),   allocatable, intent(in) :: ubvec(:)
  integer(idx_t),              intent(in) :: options(:)
  integer(idx_t),           intent(inout) :: objval
  integer(idx_t),           intent(inout) :: part(:)
  !
  !.. Local Scalars ..
  integer(c_int) :: ierr
  !
continue
  !
  ierr = METIS_PartGraphRecursive_C_Interface(nvtxs,ncon,xadj,adjncy,vwgt, &
                                              vsize,adjwgt,nparts,tpwgts, &
                                              ubvec,options,objval,part)
  !
end subroutine METIS_PartGraphRecursive
!
!###############################################################################
!
subroutine METIS_PartGraphKway(nvtxs,ncon,xadj,adjncy,vwgt,vsize,adjwgt, &
                               nparts,tpwgts,ubvec,options,objval,part)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int
  !
  !.. Formal Arguments ..
  integer(idx_t),              intent(in) :: nvtxs
  integer(idx_t),              intent(in) :: ncon
  integer(idx_t),              intent(in) :: xadj(:)
  integer(idx_t),              intent(in) :: adjncy(:)
  integer(idx_t), allocatable, intent(in) :: vwgt(:)
  integer(idx_t), allocatable, intent(in) :: vsize(:)
  integer(idx_t), allocatable, intent(in) :: adjwgt(:)
  integer(idx_t),              intent(in) :: nparts
  real(real_t),   allocatable, intent(in) :: tpwgts(:)
  real(real_t),   allocatable, intent(in) :: ubvec(:)
  integer(idx_t),              intent(in) :: options(:)
  integer(idx_t),           intent(inout) :: objval
  integer(idx_t),           intent(inout) :: part(:)
  !
  !.. Local Scalars ..
  integer(c_int) :: ierr
  !
continue
  !
  ierr = METIS_PartGraphKway_C_Interface(nvtxs,ncon,xadj,adjncy,vwgt, &
                                         vsize,adjwgt,nparts,tpwgts, &
                                         ubvec,options,objval,part)
  !
end subroutine METIS_PartGraphKway
!
!###############################################################################
!
subroutine METIS_PartMeshDual(ne,nn,eptr,eind,vwgt,vsize,ncommon,nparts, &
                              tpwgts,options,objval,epart,npart)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int
  !
  !.. Formal Arguments ..
  integer(idx_t),    intent(in) :: ne
  integer(idx_t),    intent(in) :: nn
  integer(idx_t),    intent(in) :: eptr(:)
  integer(idx_t),    intent(in) :: eind(:)
  integer(idx_t),    intent(in) :: vwgt(:)
  integer(idx_t),    intent(in) :: vsize(:)
  integer(idx_t),    intent(in) :: ncommon
  integer(idx_t),    intent(in) :: nparts
  real(real_t),      intent(in) :: tpwgts(:)
  integer(idx_t),    intent(in) :: options(:)
  integer(idx_t), intent(inout) :: objval
  integer(idx_t), intent(inout) :: epart(:)
  integer(idx_t), intent(inout) :: npart(:)
  !
  !.. Local Scalars ..
  integer(c_int) :: ierr
  !
continue
  !
  ierr = METIS_PartMeshDual_C_Interface(ne,nn,eptr,eind,vwgt,vsize, &
                                        ncommon,nparts,tpwgts,options, &
                                        objval,epart,npart)
  !
end subroutine METIS_PartMeshDual
!
!###############################################################################
!
subroutine METIS_PartMeshNodal(ne,nn,eptr,eind,vwgt,vsize,nparts, &
                               tpwgts,options,objval,epart,npart)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int
  !
  !.. Formal Arguments ..
  integer(idx_t),    intent(in) :: ne
  integer(idx_t),    intent(in) :: nn
  integer(idx_t),    intent(in) :: eptr(:)
  integer(idx_t),    intent(in) :: eind(:)
  integer(idx_t),    intent(in) :: vwgt(:)
  integer(idx_t),    intent(in) :: vsize(:)
  integer(idx_t),    intent(in) :: nparts
  real(real_t),      intent(in) :: tpwgts(:)
  integer(idx_t),    intent(in) :: options(:)
  integer(idx_t), intent(inout) :: objval
  integer(idx_t), intent(inout) :: epart(:)
  integer(idx_t), intent(inout) :: npart(:)
  !
  !.. Local Scalars ..
  integer(c_int) :: ierr
  !
continue
  !
  ierr = METIS_PartMeshNodal_C_Interface(ne,nn,eptr,eind,vwgt,vsize,nparts, &
                                         tpwgts,options,objval,epart,npart)
  !
end subroutine METIS_PartMeshNodal
!
!###############################################################################
!
subroutine METIS_NodeND(nvtxs,xadj,adjncy,vwgt,options,perm,iperm)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int
  !
  !.. Formal Arguments ..
  integer(idx_t),    intent(in) :: nvtxs
  integer(idx_t),    intent(in) :: xadj(:)
  integer(idx_t),    intent(in) :: adjncy(:)
  integer(idx_t),    intent(in) :: vwgt(:)
  integer(idx_t),    intent(in) :: options(:)
  integer(idx_t), intent(inout) :: perm(:)
  integer(idx_t), intent(inout) :: iperm(:)
  !
  !.. Local Scalars ..
  integer(c_int) :: ierr
  !
continue
  !
  ierr = METIS_NodeND_C_Interface(nvtxs,xadj,adjncy,vwgt,options,perm,iperm)
  !
end subroutine METIS_NodeND
!
!###############################################################################
!
subroutine METIS_MeshToDual(ne,nn,eptr,eind,ncommon,numflag,xadj,adjncy)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int,c_ptr
  use, intrinsic :: iso_c_binding, only : c_f_pointer
  !
  !.. Formal Arguments ..
  integer(idx_t),                 intent(in) :: ne
  integer(idx_t),                 intent(in) :: nn
  integer(idx_t),                 intent(in) :: eptr(:)
  integer(idx_t),                 intent(in) :: eind(:)
  integer(idx_t),                 intent(in) :: ncommon
  integer(idx_t),                 intent(in) :: numflag
  integer(idx_t), allocatable, intent(inout) :: xadj(:)
  integer(idx_t), allocatable, intent(inout) :: adjncy(:)
  !
  !.. Local Scalars ..
  integer :: n
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
  ierr = METIS_MeshToDual_C_Interface(ne,nn,eptr,eind,ncommon, &
                                      numflag,r_xadj,r_adjncy)
  !
  n = ne+1
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
end subroutine METIS_MeshToDual
!
!###############################################################################
!
subroutine METIS_MeshToNodal(ne,nn,eptr,eind,numflag,xadj,adjncy)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int,c_ptr
  use, intrinsic :: iso_c_binding, only : c_f_pointer
  !
  !.. Formal Arguments ..
  integer(idx_t),                 intent(in) :: ne
  integer(idx_t),                 intent(in) :: nn
  integer(idx_t),                 intent(in) :: eptr(:)
  integer(idx_t),                 intent(in) :: eind(:)
  integer(idx_t),                 intent(in) :: numflag
  integer(idx_t), allocatable, intent(inout) :: xadj(:)
  integer(idx_t), allocatable, intent(inout) :: adjncy(:)
  !
  !.. Local Scalars ..
  integer :: n
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
  ierr = METIS_MeshToNodal_C_Interface(ne,nn,eptr,eind,numflag,r_xadj,r_adjncy)
  !
  n = nn+1
  call c_f_pointer( cptr=r_xadj , fptr=ptr , shape=[n] )
  !
  allocate ( xadj(1:n) )
  xadj(1:n) = ptr(1:n)
  ptr => null()
  !
  n = xadj(nn+1)
  call c_f_pointer( cptr=r_adjncy , fptr=ptr , shape=[n] )
  !
  allocate ( adjncy(1:n) )
  adjncy(1:n) = ptr(1:n)
  ptr => null()
  !
  ierr = METIS_Free(r_xadj)
  ierr = METIS_Free(r_adjncy)
  !
end subroutine METIS_MeshToNodal
!
!###############################################################################
!
subroutine METIS_SetDefaultOptions(options)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int
  !
  !.. Formal Arguments ..
  integer(idx_t), intent(inout) :: options(:)
  !
  !.. Local Scalars ..
  integer(c_int) :: ierr
  !
continue
  !
  ierr = METIS_SetDefaultOptions_C_Interface(options)
  !
end subroutine METIS_SetDefaultOptions
!
!###############################################################################
!
end module metis5_mod
