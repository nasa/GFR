module octtree_mod
!
! Base module for performing an octtree based search through xyz data.
!
!#define PROFILE_ON
! -------------------------------------------------------------------------
! #ifdef compiler derectives:
!
! OCTID             Store a unique id in each octnode (useful for debugging).
! DEBUG_OCTDEALLOC  Print debug info while deallocating memory.
! DEBUG_OCTFILL     Print debug info while filling the octree.
! DEBUG_OCTFIND     Print debug info while finding a point from the octtree.
! DEBUG_OCTPRINT    Print octreee to file "octtree.print" during fill
!                   process.  Updates after every 10th point is added.
! -------------------------------------------------------------------------
!
  use module_kind_types, only : wp, alloc_error, error_message
  use module_kind_types, only : start_accumulation
  use module_kind_types, only : stop_accumulation
  use module_kind_types, only : reset_accumulated_time
  use module_kind_types, only : report_accumulated_time
  use module_kind_types, only : start_timer,stop_timer
  use module_kind_types, only : entering_procedure
  use module_kind_types, only : leaving_procedure
  use module_kind_types, only : mypnum,mpi_wtime,iout,zero
  use geovar, only : nr
  use octstat_mod
  use octroot_mod
  use octedat_mod
  use octnode_mod
  implicit none
  private
  public :: octtree_dealloc
  public :: octtree_init
  public :: octtree_addxyz
  public :: octtree_find_nearest_xyz
!
! -------------------------------------------------------------------------
! Abstract base octtree type.
  type, public :: octtree_type
    class(octedat_type) ,pointer :: octedat
    type(octroot_type)           :: octroot
    type(octnode_type)           :: octnode
! contains
!   final :: octtree_final
  end type octtree_type
!
contains
!
!##########################################################################
  subroutine octtree_final(octtree)
  !
  ! Finalize the octtree_type.
  !
  implicit none
  !
  !.. Formal Arguments ..
  type(octtree_type) ,intent(inout) :: octtree
  !
  !.. Local Scalars ..
  logical :: dealloc_octedat
  integer :: status
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octtree_final"
  !
continue
  !
!-DEBUG_BEG:
#ifdef DEBUG_OCTDEALLOC
  write(iout,*) 'DEBUG[octtree_final] Call octtree_dealloc.'
  flush(iout)
#endif
!-DEBUG_END:
  dealloc_octedat = .true.
  call octtree_dealloc(octtree,dealloc_octedat,status)
!-DEBUG_BEG:
#ifdef DEBUG_OCTDEALLOC
  write(iout,*) 'DEBUG[octtree_final] Leaving subroutine.'
  flush(iout)
#endif
!-DEBUG_END:
  return
  end subroutine octtree_final
!##########################################################################
  subroutine octtree_dealloc(octtree,dealloc_octedat,status)
  !
  ! Deallocate the octtree and all subnode data.
  !
  implicit none
  !
  !.. Formal Arguments ..
  class(octtree_type) ,intent(inout) :: octtree
  logical             ,intent(in)    :: dealloc_octedat
  integer             ,intent(out)   :: status
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octtree_dealloc"
  !
continue
!-DEBUG_BEG:
#ifdef DEBUG_OCTDEALLOC
  write(iout,*) 'DEBUG[octtree_dealloc] Deallocating octtree components.'
  flush(iout)
#endif
!-DEBUG_END:
  !
  status = OCT_STATUS_OK
  !
  ! Clear the octnode info.
  call octnode_dealloc(octtree%octnode,status)
  !
  ! Clear the octroot info.
  call octtree%octroot%dealloc(status)
  !
  ! Clear the octedat info.
  if ( dealloc_octedat ) call octtree%octedat%dealloc(status)
  !
!-DEBUG_BEG:
#ifdef DEBUG_OCTDEALLOC
  write(iout,*) 'DEBUG[octtree_dealloc] Leaving subroutine.'
  flush(iout)
#endif
!-DEBUG_END:
  return
  end subroutine octtree_dealloc
!##########################################################################
  subroutine octtree_init(octtree,bbox,status)
  !
  ! Initialize the empty octtree but do not polulate.
  ! Initialize octtree%octedat elsewhere, since we don't know anything
  ! about octedat structure here.
  !
  implicit none
  !
  !.. Formal Arguments ..
  type(octtree_type)         ,intent(inout) :: octtree
  real(wp)                   ,intent(in)    :: bbox(1:nr,1:2)
  integer                    ,intent(out)   :: status
  !
  !.. Local Scalars ..
  logical :: dealloc_octedat
  integer :: maxoctdata,maxoctnodes
  integer :: l,k
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octtree_init"
  !
continue
  !
  status = OCT_STATUS_OK  ! no error yet
  !
  ! Make sure this octtree is empty.
  dealloc_octedat = .false. ! don't erase octedat data
  call octtree_dealloc(octtree,dealloc_octedat,status)
  if ( status .ne. OCT_STATUS_OK ) return
  !
  ! Init the octroot.
  maxoctdata  = -1       ! use default from octnode_mod
  maxoctnodes = 2**nr
  call octroot_init(octtree%octroot,maxoctdata,maxoctnodes,status)
  if ( status .ne. OCT_STATUS_OK ) return
  !
  ! Init the initial octnode.
  call octnode_init(octtree%octnode,octtree%octroot, &
                    bbox,status) ! ,parent=null()
  if ( status .ne. OCT_STATUS_OK ) return
  !
! Push data into octtree.
! For a simple array src_xyzdata(3,npts):
!   do k = 1,size(src_xyzdata,dim=2)
!     call octnode_addxyz(octtree%octnode,octtree%octroot,k,status)
!     if ( status .ne. OCT_STATUS_OK ) return
!   end do
! However, we 
! call ed_data%octtree_fill(octtree,
  !
  return
  end subroutine octtree_init
!##########################################################################
  subroutine octtree_addxyz(octtree,ixyz_pt,xyz_pt,status)
  !
  ! Add an xyz point to this octtree.
  !
  implicit none
  !
  !.. Formal Arguments ..
  type(octtree_type)                    ,intent(inout) :: octtree
  integer              ,dimension(:)    ,intent(in)    :: ixyz_pt
  real(wp)             ,dimension(1:nr) ,intent(in)    :: xyz_pt
  integer                               ,intent(out)   :: status
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octtree_addxyz"
  !
continue
  !
  call octnode_addxyz(octtree%octnode,octtree%octroot,octtree%octedat, &
                          ixyz_pt,xyz_pt,status)
  return
  end subroutine octtree_addxyz
!##########################################################################
  subroutine octtree_find_containing_octnode(octtree, &
               xyz_pt, containing_octnode,status)
  !
  ! Within this octtree, find the smallest octnode or sub-octnode whose bbox
  ! contains xyz_pt.  This node will have no children, and may or may not
  ! have data.
  !
  implicit none
  !
  !.. Formal Arguments ..
  type(octtree_type)                    ,intent(in)  :: octtree
  real(wp) ,dimension(:)                ,intent(in)  :: xyz_pt
  type(octnode_type) ,pointer           ,intent(out) :: containing_octnode
  integer                               ,intent(out) :: status
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octtree_find_containing_octnode"
  !
continue
  !
  ! Set default containing octnode
  containing_octnode => null()
  !
  ! Search the octtree starting from the head octnode.
  call octnode_find_containing_octnode(octtree%octnode, &
         octtree%octroot,octtree%octedat, &
         xyz_pt,containing_octnode,status)
  !
  return
  end subroutine octtree_find_containing_octnode
!##########################################################################
  subroutine octtree_find_nearest_xyz(octtree,xyz_pt, ixyz_pt,status,tid,acc)
  !
  ! Find point in octtree that is nearest to desired xyz_pt.
  ! Return the index info used to retrieve the corresponding xyz point
  ! from routine octtree%octedat%evalxyz_pt(idx,xyz,status).
  !
  implicit none
  !
  !.. Formal Arguments ..
  type(octtree_type)               ,intent(in)  :: octtree
  real(wp)           ,dimension(:) ,intent(in)  :: xyz_pt
  integer            ,dimension(:) ,intent(out) :: ixyz_pt
  integer                          ,intent(out) :: status
  !
  integer, optional, intent(in) :: tid
  integer, optional, intent(in) :: acc
  !
  !.. Local Scalars ..
  type(octnode_type) ,pointer                   :: containing_octnode
  type(octnode_type) ,pointer                   :: best_octnode
  integer                                       :: best_octixyz
  real(wp)                                      :: best_dist2
  real(wp)                                      :: bbox_diag
  integer                                       :: ierr
!-DEBUG_BEG:
#ifdef DEBUG_OCTFIND
  integer l
#endif
!-DEBUG_END:
  !
  !.. Local Arrays ..
  integer            ,dimension(:) ,allocatable :: best_ixyz
  real(wp)           ,dimension(1:nr,1:2)       :: xyz_bbox
  real(wp)           ,dimension(1:nr)           :: bbox_size
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octtree_find_nearest_xyz"
  !
continue
  !
  if (present(acc)) then
    if (acc == 1) then
      call debug_timer(reset_accumulated_time,"Crap 1", &
                       acc=1,pid=tid)
      call debug_timer(reset_accumulated_time,"octnode_find_containing_node", &
                       acc=2,pid=tid)
      call debug_timer(reset_accumulated_time,"Crap 2", &
                       acc=3,pid=tid)
      call debug_timer(reset_accumulated_time,"octnode_find_nearest_xyz", &
                       acc=4,pid=tid)
      call debug_timer(reset_accumulated_time,"Crap 3", &
                       acc=5,pid=tid)
      return
    else if (acc == 2) then
      call debug_timer(report_accumulated_time,"Crap 1", &
                       acc=1,pid=tid)
      call debug_timer(report_accumulated_time,"octnode_find_containing_node", &
                       acc=2,pid=tid)
      call debug_timer(report_accumulated_time,"Crap 2", &
                       acc=3,pid=tid)
      call debug_timer(report_accumulated_time,"octnode_find_nearest_xyz", &
                       acc=4,pid=tid)
      call debug_timer(report_accumulated_time,"Crap 3", &
                       acc=5,pid=tid)
      return
    end if
  end if

  ! Set default values, in case no point is found.
  call debug_timer(start_accumulation,"Crap 1", &
                   acc=1,pid=tid)
  status = OCT_STATUS_NOT_FOUND
  containing_octnode => null()
  best_octnode => null()
  best_octixyz = 0
  best_dist2   = huge(best_dist2)
  ixyz_pt = 0
  !
  allocate( best_ixyz(octtree%octedat%maxindices) , stat=ierr , &
            errmsg=error_message )
  call alloc_error(pname,"best_ixyz",1,__LINE__,__FILE__,ierr,error_message)
  if ( ierr .ne. 0 ) then
    status = OCT_STATUS_MEMORY_ERROR
    return
  endif
  best_ixyz = 0
  !
  call debug_timer(stop_accumulation,"Crap 1", &
                   acc=1,pid=tid)
  call debug_timer(start_accumulation,"octnode_find_containing_octnode", &
                   acc=2,pid=tid)
  ! Find the smallest octnode containing xyz_pt.
  call octnode_find_containing_octnode(octtree%octnode, &
         octtree%octroot,octtree%octedat, xyz_pt,containing_octnode,status)
  call debug_timer(stop_accumulation,"octnode_find_containing_octnode", &
                   acc=2,pid=tid)
  !
  if ( status .eq. OCT_STATUS_OK ) then
    call debug_timer(start_accumulation,"Crap 2", &
                     acc=3,pid=tid)
    best_octnode => containing_octnode
    ! Now have the smallest octnode containing xyz_pt.
    ! Nearest point could be in the containing_octnode or an adjacent octnode
    ! (either a sibling octnode or a distance cousin from the octtree).
    !
    ! Know there is at least one point in containing_octnode or its siblings.
    ! Its distance is within twice the diagonal of containing_octnode bbox.
    ! Therefore, create a shpere (or cubic bounding box) of r=2*diagonal
    ! that is centered at xyz_pt, and search all octnodes that intersect it
    ! to see if there is a closer point.
    bbox_size = containing_octnode%bbox(:,2) - &
                containing_octnode%bbox(:,1)
    bbox_diag = sqrt(sum(bbox_size * bbox_size))
    xyz_bbox(:,1) = xyz_pt(:) - 2*bbox_diag
    xyz_bbox(:,2) = xyz_pt(:) + 2*bbox_diag
    !
    call debug_timer(stop_accumulation,"Crap 2", &
                     acc=3,pid=tid)
    ! Search octtree octnodes for point nearest to xyz_pt.
    call debug_timer(start_accumulation,"octnode_find_nearest_xyz", &
                     acc=4,pid=tid)
    call octnode_find_nearest_xyz(octtree%octnode,octtree%octroot, &
           octtree%octedat,xyz_pt,xyz_bbox, &
           best_octnode,best_octixyz,best_ixyz,best_dist2,status)
    if ( status .eq. OCT_STATUS_OK ) then
      ixyz_pt = best_ixyz
    end if
    call debug_timer(stop_accumulation,"octnode_find_nearest_xyz", &
                     acc=4,pid=tid)
  end if
  call debug_timer(start_accumulation,"Crap 3", &
                   acc=5,pid=tid)
  !
#ifdef DEBUG_OCTFIND
#ifdef OCTID
!-DEBUG_BEG:
  if ( associated(containing_octnode) ) then
    write(*,*) '  Containing octnode = ', containing_octnode%octid
    write(*,*) '    bbox_min         = ',(containing_octnode%bbox(l,1), l=1,nr)
    write(*,*) '    bbox_max         = ',(containing_octnode%bbox(l,2), l=1,nr)
  else
    write(*,*) '  Containing octnode = null'
  end if
  if ( associated(best_octnode) ) then
    write(*,*) '  Best       octnode = ',       best_octnode%octid
    write(*,*) '    bbox_min         = ',(      best_octnode%bbox(l,1), l=1,nr)
    write(*,*) '    bbox_max         = ',(      best_octnode%bbox(l,2), l=1,nr)
  else
    write(*,*) '  Best       octnode = null'
  end if
!-DEBUG_END:
#endif
#endif
  !
  if ( allocated(best_ixyz) ) then
    deallocate( best_ixyz , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"best_ixyz",2,__LINE__,__FILE__,ierr,error_message)
  end if
  call debug_timer(stop_accumulation,"Crap 3", &
                   acc=5,pid=tid)
  !
  return
  end subroutine octtree_find_nearest_xyz
!##########################################################################
  subroutine octtree_print(octtree,filename,status)
  !
  ! Print the octtree to file for debugging.
  !
  implicit none
  !
  !.. Formal Arguments ..
  type(octtree_type)                    ,intent(in)    :: octtree
  character (len=*)                     ,intent(in)    :: filename
  integer                               ,intent(out)   :: status
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octtree_print"
  !
continue
  !
  status = OCT_STATUS_OK
  !
  call octnode_print(octtree%octnode,octtree%octedat,status,filename)
  !
  return
  end subroutine octtree_print
!
!##########################################################################
!
#include "Functions/debug_timer.f90"
!
!##########################################################################
!
end module octtree_mod
