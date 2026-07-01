!###############################################################################
!> This module provides support for Metis 5.
!###############################################################################
module metis_mod
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
  private :: create_partition_weights
  private :: check_METIS_error_value
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
      use, intrinsic :: iso_c_binding, only : c_int,c_ptr
      import :: idx_t, real_t
      integer(c_int) :: METIS_PartGraphRecursive_C_Interface
      integer(idx_t) :: nvtxs
      integer(idx_t) :: ncon
      integer(idx_t) :: xadj(*)
      integer(idx_t) :: adjncy(*)
      type(c_ptr), value :: vwgt
      type(c_ptr), value :: vsize
      type(c_ptr), value :: adjwgt
      integer(idx_t) :: nparts
      type(c_ptr), value :: tpwgts
      type(c_ptr), value :: ubvec
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
      use, intrinsic :: iso_c_binding, only : c_int,c_ptr
      import :: idx_t, real_t
      integer(c_int) :: METIS_PartGraphKway_C_Interface
      integer(idx_t) :: nvtxs
      integer(idx_t) :: ncon
      integer(idx_t) :: xadj(*)
      integer(idx_t) :: adjncy(*)
      type(c_ptr), value :: vwgt
      type(c_ptr), value :: vsize
      type(c_ptr), value :: adjwgt
      integer(idx_t) :: nparts
      type(c_ptr), value :: tpwgts
      type(c_ptr), value :: ubvec
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
  ! Debugging METIS is achieved by setting METIS_DEBUG_VAL > 0
  ! METIS_DEBUG_VAL = [any sum of the following DBGLVL parameters]
  !
  integer(idx_t), parameter :: METIS_DEBUG_VAL = int( 0 &
                                                     !+ METIS_DBG_INFO &
                                                     !+ METIS_DBG_TIME &
                                                     !+ METIS_DBG_COARSEN &
                                                     !+ METIS_DBG_REFINE &
                                                     !+ METIS_DBG_IPART &
                                                     !+ METIS_DBG_MOVEINFO &
                                                     !+ METIS_DBG_SEPINFO &
                                                     !+ METIS_DBG_CONNINFO &
                                                     !+ METIS_DBG_CONTIGINFO &
                                                     !+ METIS_DBG_MEMORY &
                                                      , kind=idx_t )
  !
  logical, parameter :: ACTIVATE_METIS_DEBUGGING = (METIS_DEBUG_VAL > 0_idx_t)
  !
  ! mobjtype_et
  !
  enum, bind(c)
    enumerator :: METIS_OBJTYPE_CUT
    enumerator :: METIS_OBJTYPE_VOL
    enumerator :: METIS_OBJTYPE_NODE
  end enum
  !
  ! Module Arrays
  !
  integer(idx_t), public, save :: options(0:METIS_NOPTIONS-1)
  !
contains
!
!###############################################################################
!
subroutine METIS_PartGraph(metis_option,xadj,adjncy, &
                           vwgt,adjwgt,nparts,objval,part)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int,c_loc,c_ptr,c_null_ptr
  use, intrinsic :: iso_fortran_env, only : iout => output_unit
  !
  !.. Formal Arguments ..
  integer,                             intent(in) :: metis_option
 !integer(idx_t),                      intent(in) :: nvtxs
 !integer(idx_t),                      intent(in) :: ncon
  integer(idx_t),                      intent(in) :: xadj(:)
  integer(idx_t),                      intent(in) :: adjncy(:)
  integer(idx_t), allocatable, target, intent(in) :: vwgt(:)
 !integer(idx_t), allocatable, target, intent(in) :: vsize(:)
  integer(idx_t), allocatable, target, intent(in) :: adjwgt(:)
  integer(idx_t),                      intent(in) :: nparts
 !real(real_t),   allocatable, target, intent(in) :: tpwgts(:)
 !real(real_t),   allocatable, target, intent(in) :: ubvec(:)
  integer(idx_t),                   intent(inout) :: objval
  integer(idx_t),                   intent(inout) :: part(:)
  !
  !.. Local Scalars ..
  integer(c_int) :: ierr
  integer(idx_t) :: ncon,nvtxs
  !
  !.. Local Allocatable Arrays ..
  integer(idx_t), allocatable, target :: vsize(:)
  real(real_t),   allocatable, target :: tpwgts(:)
  real(real_t),   allocatable, target :: ubvec(:)
  !
  !.. Local c_ptr Scalars ..
  type(c_ptr) :: c_vwgt
  type(c_ptr) :: c_vsize
  type(c_ptr) :: c_adjwgt
  type(c_ptr) :: c_tpwgts
  type(c_ptr) :: c_ubvec
  !
continue
  !
  ! Set up the options array for METIS
  !
  call set_METIS_options
  !
  nvtxs = size(xadj,kind=idx_t) - 1_idx_t
  !
  ncon = get_number_of_vertex_weights(xadj,vwgt)
  !
  c_tpwgts = c_null_ptr
  c_ubvec = c_null_ptr
  c_vsize = c_null_ptr
  c_vwgt = c_null_ptr
  c_adjwgt = c_null_ptr
  !
  if (any(metis_option == [4,5,6])) then
    call create_partition_weights(nparts,ncon,nvtxs,tpwgts,ubvec,vsize)
    if (allocated(vwgt)) c_vwgt = c_loc(vwgt)
    if (allocated(adjwgt)) c_adjwgt = c_loc(adjwgt)
    if (allocated(vsize)) c_vsize = c_loc(vsize)
    if (allocated(tpwgts)) c_tpwgts = c_loc(tpwgts)
    if (allocated(ubvec)) c_ubvec = c_loc(ubvec)
  end if
  !
  ! The following logic should default to metis_option=3 if the
  ! input metis_option isnt valid (a value of 1-6)
  !
  if (any(metis_option == [1,6])) then
    !
    ierr = METIS_PartGraphRecursive_C_Interface(nvtxs,ncon,xadj,adjncy, &
                                                c_vwgt,c_vsize,c_adjwgt, &
                                                nparts,c_tpwgts,c_ubvec, &
                                                options,objval,part)
    !
    ! Check the error code returned by the METIS function
    !
    call check_METIS_error_value(ierr,"METIS_PartGraphRecursive_C_Interface", &
                                 __LINE__,__FILE__)
    !
  else
    !
    if (any(metis_option == [2,5])) then
      options(METIS_OPTION_OBJTYPE) = METIS_OBJTYPE_CUT
    else
      options(METIS_OPTION_OBJTYPE) = METIS_OBJTYPE_VOL
    end if
    !
    ierr = METIS_PartGraphKway_C_Interface(nvtxs,ncon,xadj,adjncy,c_vwgt, &
                                           c_vsize,c_adjwgt,nparts,c_tpwgts, &
                                           c_ubvec,options,objval,part)
    !
    ! Check the error code returned by the METIS function
    !
    call check_METIS_error_value(ierr,"METIS_PartGraphKway_C_Interface", &
                                 __LINE__,__FILE__)
    !
  end if
  write (iout,'(a,i0)') &
    "Returned from METIS_PartGraphKway_C_Interface with error code = ",ierr
  flush (iout)
  !
  ! Reset the C address pointers to null
  !
  c_vwgt   = c_null_ptr
  c_adjwgt = c_null_ptr
  c_vsize = c_null_ptr
  c_tpwgts = c_null_ptr
  c_ubvec  = c_null_ptr
  !
  ! Deallocate the local arrays
  !
  if (allocated(vsize )) deallocate ( vsize )
  if (allocated(tpwgts)) deallocate ( tpwgts )
  if (allocated(ubvec )) deallocate ( ubvec )
  !
end subroutine METIS_PartGraph
!
!###############################################################################
!
subroutine METIS_PartGraphRecursive(nvtxs,ncon,xadj,adjncy,vwgt,adjwgt, &
                                    nparts,objval,part)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int,c_loc,c_ptr,c_null_ptr
  !
  !.. Formal Arguments ..
  integer(idx_t),                      intent(in) :: nvtxs
  integer(idx_t),                      intent(in) :: ncon
  integer(idx_t),                      intent(in) :: xadj(:)
  integer(idx_t),                      intent(in) :: adjncy(:)
  integer(idx_t), allocatable, target, intent(in) :: vwgt(:)
 !integer(idx_t), allocatable, target, intent(in) :: vsize(:)
  integer(idx_t), allocatable, target, intent(in) :: adjwgt(:)
  integer(idx_t),                      intent(in) :: nparts
 !real(real_t),   allocatable, target, intent(in) :: tpwgts(:)
 !real(real_t),   allocatable, target, intent(in) :: ubvec(:)
  integer(idx_t),                   intent(inout) :: objval
  integer(idx_t),                   intent(inout) :: part(:)
  !
  !.. Local Scalars ..
  integer(c_int) :: ierr
  !
  !.. Local Allocatable Arrays ..
  integer(idx_t), allocatable, target :: vsize(:)
  real(real_t),   allocatable, target :: tpwgts(:)
  real(real_t),   allocatable, target :: ubvec(:)
  !
  !.. Local c_ptr Scalars ..
  type(c_ptr) :: c_vwgt
  type(c_ptr) :: c_vsize
  type(c_ptr) :: c_adjwgt
  type(c_ptr) :: c_tpwgts
  type(c_ptr) :: c_ubvec
  !
continue
  !
  ! Set up the options array for METIS
  !
  call set_METIS_options
  !
  c_tpwgts = c_null_ptr
  c_ubvec = c_null_ptr
  c_vsize = c_null_ptr
  c_vwgt = c_null_ptr
  c_adjwgt = c_null_ptr
  !
  call create_partition_weights(nparts,ncon,nvtxs,tpwgts,ubvec,vsize)
  if (allocated(vwgt)) c_vwgt = c_loc(vwgt)
  if (allocated(adjwgt)) c_adjwgt = c_loc(adjwgt)
  if (allocated(vsize)) c_vsize = c_loc(vsize)
  if (allocated(tpwgts)) c_tpwgts = c_loc(tpwgts)
  if (allocated(ubvec)) c_ubvec = c_loc(ubvec)
  !
  ierr = METIS_PartGraphRecursive_C_Interface(nvtxs,ncon,xadj,adjncy,c_vwgt, &
                                              c_vsize,c_adjwgt,nparts,c_tpwgts, &
                                              c_ubvec,options,objval,part)
  !
  ! Check the error code returned by the METIS function
  !
  call check_METIS_error_value(ierr,"METIS_PartGraphRecursive_C_Interface", &
                               __LINE__,__FILE__)
  !
  ! Reset the C address pointers to null
  !
  c_vwgt   = c_null_ptr
  c_adjwgt = c_null_ptr
  c_vsize = c_null_ptr
  c_tpwgts = c_null_ptr
  c_ubvec  = c_null_ptr
  !
  ! Deallocate the local arrays
  !
  if (allocated(vsize )) deallocate ( vsize )
  if (allocated(tpwgts)) deallocate ( tpwgts )
  if (allocated(ubvec )) deallocate ( ubvec )
  !
end subroutine METIS_PartGraphRecursive
!
!###############################################################################
!
subroutine METIS_PartGraphKway(nvtxs,ncon,xadj,adjncy,vwgt,adjwgt, &
                               nparts,objval,part)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int,c_loc,c_ptr,c_null_ptr
  use, intrinsic :: iso_fortran_env, only : iout => output_unit
  !
  !.. Formal Arguments ..
  integer(idx_t),                      intent(in) :: nvtxs
  integer(idx_t),                      intent(in) :: ncon
  integer(idx_t),                      intent(in) :: xadj(:)
  integer(idx_t),                      intent(in) :: adjncy(:)
  integer(idx_t), allocatable, target, intent(in) :: vwgt(:)
 !integer(idx_t), allocatable, target, intent(in) :: vsize(:)
  integer(idx_t), allocatable, target, intent(in) :: adjwgt(:)
  integer(idx_t),                      intent(in) :: nparts
 !real(real_t),   allocatable, target, intent(in) :: tpwgts(:)
 !real(real_t),   allocatable, target, intent(in) :: ubvec(:)
  integer(idx_t),                   intent(inout) :: objval
  integer(idx_t),                   intent(inout) :: part(:)
  !
  !.. Local Scalars ..
  integer(c_int) :: ierr
  !
  !.. Local Allocatable Arrays ..
  integer(idx_t), allocatable, target :: vsize(:)
  real(real_t),   allocatable, target :: tpwgts(:)
  real(real_t),   allocatable, target :: ubvec(:)
  !
  !.. Local c_ptr Scalars ..
  type(c_ptr) :: c_vwgt
  type(c_ptr) :: c_vsize
  type(c_ptr) :: c_adjwgt
  type(c_ptr) :: c_tpwgts
  type(c_ptr) :: c_ubvec
  !
continue
  !
  ! Set up the options array for METIS
  !
  call set_METIS_options
  !
  c_tpwgts = c_null_ptr
  c_ubvec = c_null_ptr
  c_vsize = c_null_ptr
  c_vwgt = c_null_ptr
  c_adjwgt = c_null_ptr
  !
  call create_partition_weights(nparts,ncon,nvtxs,tpwgts,ubvec,vsize)
  if (allocated(vwgt)) c_vwgt = c_loc(vwgt)
  if (allocated(adjwgt)) c_adjwgt = c_loc(adjwgt)
  if (allocated(vsize)) c_vsize = c_loc(vsize)
  if (allocated(tpwgts)) c_tpwgts = c_loc(tpwgts)
  if (allocated(ubvec)) c_ubvec = c_loc(ubvec)
  !
  ierr = METIS_PartGraphKway_C_Interface(nvtxs,ncon,xadj,adjncy,c_vwgt, &
                                         c_vsize,c_adjwgt,nparts,c_tpwgts, &
                                         c_ubvec,options,objval,part)
  !
  ! Check the error code returned by the METIS function
  !
  call check_METIS_error_value(ierr,"METIS_PartGraphKway_C_Interface", &
                               __LINE__,__FILE__)
  write (iout,'(a,i0)') &
    "Returned from METIS_PartGraphKway_C_Interface with error code = ",ierr
  flush (iout)
  !
  ! Reset the C address pointers to null
  !
  c_vwgt   = c_null_ptr
  c_adjwgt = c_null_ptr
  c_vsize = c_null_ptr
  c_tpwgts = c_null_ptr
  c_ubvec  = c_null_ptr
  !
  ! Deallocate the local arrays
  !
  if (allocated(vsize )) deallocate ( vsize )
  if (allocated(tpwgts)) deallocate ( tpwgts )
  if (allocated(ubvec )) deallocate ( ubvec )
  !
end subroutine METIS_PartGraphKway
!
!###############################################################################
!
subroutine METIS_PartMeshDual(ne,nn,eptr,eind,vwgt,vsize,ncommon,nparts, &
                              tpwgts,objval,epart,npart)
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
  integer(idx_t), intent(inout) :: objval
  integer(idx_t), intent(inout) :: epart(:)
  integer(idx_t), intent(inout) :: npart(:)
  !
  !.. Local Scalars ..
  integer(c_int) :: ierr
  !
continue
  !
  ! Set up the options array for METIS
  !
  call set_METIS_options
  !
  ierr = METIS_PartMeshDual_C_Interface(ne,nn,eptr,eind,vwgt,vsize, &
                                        ncommon,nparts,tpwgts,options, &
                                        objval,epart,npart)
  !
  ! Check the error code returned by the METIS function
  !
  call check_METIS_error_value(ierr,"METIS_PartMeshDual_C_Interface", &
                               __LINE__,__FILE__)
  !
end subroutine METIS_PartMeshDual
!
!###############################################################################
!
subroutine METIS_PartMeshNodal(ne,nn,eptr,eind,vwgt,vsize,nparts, &
                               tpwgts,objval,epart,npart)
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
  integer(idx_t), intent(inout) :: objval
  integer(idx_t), intent(inout) :: epart(:)
  integer(idx_t), intent(inout) :: npart(:)
  !
  !.. Local Scalars ..
  integer(c_int) :: ierr
  !
continue
  !
  ! Set up the options array for METIS
  !
  call set_METIS_options
  !
  ierr = METIS_PartMeshNodal_C_Interface(ne,nn,eptr,eind,vwgt,vsize,nparts, &
                                         tpwgts,options,objval,epart,npart)
  !
  ! Check the error code returned by the METIS function
  !
  call check_METIS_error_value(ierr,"METIS_PartMeshNodal_C_Interface", &
                               __LINE__,__FILE__)
  !
end subroutine METIS_PartMeshNodal
!
!###############################################################################
!
subroutine METIS_NodeND(nvtxs,xadj,adjncy,vwgt,perm,iperm)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int
  !
  !.. Formal Arguments ..
  integer(idx_t),    intent(in) :: nvtxs
  integer(idx_t),    intent(in) :: xadj(:)
  integer(idx_t),    intent(in) :: adjncy(:)
  integer(idx_t),    intent(in) :: vwgt(:)
  integer(idx_t), intent(inout) :: perm(:)
  integer(idx_t), intent(inout) :: iperm(:)
  !
  !.. Local Scalars ..
  integer(c_int) :: ierr
  !
continue
  !
  ! Set up the options array for METIS
  !
  call set_METIS_options
  !
  ierr = METIS_NodeND_C_Interface(nvtxs,xadj,adjncy,vwgt,options,perm,iperm)
  !
  ! Check the error code returned by the METIS function
  !
  call check_METIS_error_value(ierr,"METIS_NodeND_C_Interface", &
                               __LINE__,__FILE__)
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
  ! Check the error code returned by the METIS function
  !
  call check_METIS_error_value(ierr,"METIS_MeshToDual_C_Interface", &
                               __LINE__,__FILE__)
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
  call check_METIS_error_value(ierr,"METIS_Free",__LINE__,__FILE__)
  ierr = METIS_Free(r_adjncy)
  call check_METIS_error_value(ierr,"METIS_Free",__LINE__,__FILE__)
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
  ! Check the error code returned by the METIS function
  !
  call check_METIS_error_value(ierr,"METIS_MeshToNodal_C_Interface", &
                               __LINE__,__FILE__)
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
  call check_METIS_error_value(ierr,"METIS_Free",__LINE__,__FILE__)
  ierr = METIS_Free(r_adjncy)
  call check_METIS_error_value(ierr,"METIS_Free",__LINE__,__FILE__)
  !
end subroutine METIS_MeshToNodal
!
!###############################################################################
!
subroutine set_METIS_options()
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int
  !
  !.. Local Scalars ..
  integer(c_int) :: ierr
  integer(c_int) :: default_dbg_opt
  !
continue
  !
  ! First set the METIS options to their default values
  !
  ierr = METIS_SetDefaultOptions_C_Interface(options)
  !
  ! Check the error code returned by the METIS function
  !
  call check_METIS_error_value(ierr,"METIS_SetDefaultOptions_C_Interface", &
                               __LINE__,__FILE__)
  !
  ! Customize some of the options
  !
  ! Set numbering to Fortran style (begins at 1)
  options(METIS_OPTION_NUMBERING) = 1
  ! Force contiguous partitions
  options(METIS_OPTION_CONTIG) = 1
  ! Explicitly minimize the maximum connectivity
  options(METIS_OPTION_MINCONN) = 1
  !
  ! Set the debugging level for METIS
  !
  default_dbg_opt = options(METIS_OPTION_DBGLVL)
  if (ACTIVATE_METIS_DEBUGGING) then
    if (default_dbg_opt > 0_idx_t) then
      options(METIS_OPTION_DBGLVL) = ior(default_dbg_opt,METIS_DEBUG_VAL)
    else
      options(METIS_OPTION_DBGLVL) = METIS_DEBUG_VAL
    end if
  end if
  !
end subroutine set_METIS_options
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
  !
  !
  ierr = METIS_SetDefaultOptions_C_Interface(options)
  !
  ! Check the error code returned by the METIS function
  !
  call check_METIS_error_value(ierr,"METIS_SetDefaultOptions_C_Interface", &
                               __LINE__,__FILE__)
  !
end subroutine METIS_SetDefaultOptions
!
!###############################################################################
!
subroutine create_partition_weights(nparts,ncon,nvtxs,tpwgts,ubvec,vsize)
  !
  !.. Formal Arguments ..
  integer(idx_t),             intent(in) :: nparts
  integer(idx_t),             intent(in) :: ncon
  integer(idx_t),             intent(in) :: nvtxs
  real(real_t),   allocatable, intent(out) :: tpwgts(:)
  real(real_t),   allocatable, intent(out) :: ubvec(:)
  integer(idx_t), allocatable, intent(out) :: vsize(:)
  !
continue
  !
  allocate ( tpwgts(1:ncon*nparts) )
  !
  tpwgts(:) = 1.0_real_t / real(nparts,kind=real_t)
  tpwgts(:) = tpwgts(:) / sum(tpwgts(:))
  !
  allocate ( ubvec(1:ncon) )
  !
  ! The serial METIS manual (v 5.1.0; March 30, 2013; section 5.8 on page 26)
  ! recommends the values
  !         1.001   if ncon==1
  !         1.01    otherwise
  ubvec(:) = merge( 1.001_real_t , 1.01_real_t , ncon==1 )
  !
  allocate ( vsize(1:nvtxs) )
  vsize(:) = 1_idx_t
  !
end subroutine create_partition_weights
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
subroutine check_METIS_error_value(ierr,pname,line_number,file_name)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding,   only : c_int
  use, intrinsic :: iso_fortran_env, only : iout => output_unit
  !
  !.. Formal Arguments ..
  integer(c_int),   intent(in) :: ierr
  character(len=*), intent(in) :: pname
  integer,          intent(in) :: line_number
  character(len=*), intent(in) :: file_name
  !
  !.. Local Scalars ..
  !
continue
  !
  if (ierr == int(METIS_OK,kind=c_int)) then
    !
    if (ACTIVATE_METIS_DEBUGGING) then
      write (iout,10) trim(adjustl(pname)), &
                      line_number, &
                      trim(adjustl(file_name))
    end if
    !
  else
    !
    write (iout,20) trim(adjustl(pname)), &
                    line_number, &
                    trim(adjustl(file_name))
    !
    if (ierr == int(METIS_ERROR_INPUT,kind=c_int)) then
      write (iout,21)
    else if (ierr == int(METIS_ERROR_MEMORY,kind=c_int)) then
      write (iout,22)
    else if (ierr == int(METIS_ERROR,kind=c_int)) then
      write (iout,23)
    else
      write (iout,24) ierr
    end if
    !
    write (iout,29)
    !
  end if
  !
  10 format (/,2x,"The call to METIS function ",a,/, &
               2x,"on line number ",i0," of file ",a," was successful!",/)
  !
  20 format (//,80("#"),/, &
                34("#"),"  ERROR!!!  ",34("#"),/, &
                80("#"),/, &
                2x,"THE CALL TO METIS FUNCTION ",a,/, &
                2x,"ON LINE NUMBER ",i0," OF FILE ",a)
  !
  21 format (2x,"RETURNED AN ERROR INDICATING SOMETHING WAS",/, &
             2x,"WRONG WITH THE INPUTS SUPPLIED TO THE FUNCTION!")
  22 format (2x,"RETURNED A MEMORY ERROR! GOOD LUCK FIGURING THAT OUT!")
  23 format (2x,"RETURNED A GENERAL ERROR! THANK YOU METIS ", &
                "FOR BEING SO HELPFUL",/, &
             2x,"WITH FIGURING OUT WHY YOU DO NOT WANT TO ", &
                "WORK WELL WITH OTHERS!")
  24 format (2x,"RETURNED AN UNKNOWN ERROR VALUE!",/, &
             2x,"        IERR = ",i0)
  !
  29 format (80("#"),//)
  !
end subroutine check_METIS_error_value
!
!###############################################################################
!
end module metis_mod
