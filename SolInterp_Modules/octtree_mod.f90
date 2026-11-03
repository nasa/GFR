module octtree_mod
  !
  ! Base module for performing an octtree based search through xyz data.
  !
  !#define PROFILE_ON
  ! -------------------------------------------------------------------------
  ! #ifdef compiler derectives:
  !
  ! OCTID             Store a unique id in each octnode (useful for debugging).
  !                   OCTID IS NOW PERMANENT
  ! DEBUG_OCTDEALLOC  Print debug info while deallocating memory.
  ! DEBUG_OCTFILL     Print debug info while filling the octree.
  ! DEBUG_OCTFIND     Print debug info while finding a point from the octtree.
  ! DEBUG_OCTPRINT    Print octreee to file "octtree.print" during fill
  !                   process.  Updates after every 10th point is added.
  ! -------------------------------------------------------------------------
  !
  !.. Use Statements ..
  use module_kind_types, only : wp,lk,alloc_error,error_message
  use module_kind_types, only : true,fals
  use module_kind_types, only : start_accumulation
  use module_kind_types, only : stop_accumulation
  use module_kind_types, only : reset_accumulated_time
  use module_kind_types, only : report_accumulated_time
  use module_kind_types, only : start_timer,stop_timer
  use module_kind_types, only : entering_procedure
  use module_kind_types, only : leaving_procedure
  use module_kind_types, only : debug_timer_disabled
  use module_kind_types, only : mypnum,mpi_wtime,iout,zero
  use module_kind_types, only : gfr_timer
  use module_kind_types, only : glb_root
  use module_kind_types, only : uppercase ! needed for Functions/debug_timer.f90
  use octstat_mod
  use octroot_mod
  use octedat_mod
  use octnode_mod
  !
  implicit none
  !
  private
  !
  public :: octtree_dealloc
  public :: octtree_init
  public :: octtree_addxyz
  public :: octtree_find_nearest_xyz
  !
  !
  ! Abstract base octtree type.
  !
  type, public :: octtree_type
    class(octedat_type), pointer :: octedat
    type(octroot_type)           :: octroot
    type(octnode_type)           :: octnode
 !contains
 !  final :: octtree_final
  end type octtree_type
  !
contains
!
!###############################################################################
!
subroutine octtree_final(octtree)
  !
  ! Finalize the octtree_type.
  !
  !.. Formal Arguments ..
  type(octtree_type), intent(inout) :: octtree
  !
  !.. Local Scalars ..
  integer :: istat
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octtree_final"
  !
continue
  !
  !-DEBUG_BEG:
#ifdef DEBUG_OCTDEALLOC
  write (iout,*) 'DEBUG[octtree_final] Call octtree_dealloc.'
  flush (iout)
#endif
  !-DEBUG_END:
  !
  call octtree_dealloc(octtree,true,istat)
  !
  !-DEBUG_BEG:
#ifdef DEBUG_OCTDEALLOC
  write (iout,*) 'DEBUG[octtree_final] Leaving subroutine.'
  flush (iout)
#endif
  !-DEBUG_END:
  !
end subroutine octtree_final
!
!###############################################################################
!
subroutine octtree_dealloc(octtree,dealloc_octedat,istat)
  !
  ! Deallocate the octtree and all subnode data.
  !
  !.. Formal Arguments ..
  class(octtree_type), intent(inout) :: octtree
  logical(lk),            intent(in) :: dealloc_octedat
  integer,               intent(out) :: istat
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octtree_dealloc"
  !
continue
  !
  !-DEBUG_BEG:
#ifdef DEBUG_OCTDEALLOC
  write (iout,*) 'DEBUG[octtree_dealloc] Deallocating octtree components.'
  flush (iout)
#endif
  !-DEBUG_END:
  !
  istat = OCT_STATUS_OK
  !
  ! Clear the octnode info.
  !
  call octnode_dealloc(octtree%octnode,istat)
  !
  ! Clear the octroot info.
  !
  call octtree%octroot%dealloc(istat)
  !
  ! Clear the octedat info.
  !
  if (dealloc_octedat) call octtree%octedat%dealloc(istat)
  !
  !-DEBUG_BEG:
#ifdef DEBUG_OCTDEALLOC
  write (iout,*) 'DEBUG[octtree_dealloc] Leaving subroutine.'
  flush (iout)
#endif
  !-DEBUG_END:
  !
end subroutine octtree_dealloc
!
!###############################################################################
!
subroutine octtree_init(octtree,bbox,istat)
  !
  ! Initialize the empty octtree but do not polulate.
  ! Initialize octtree%octedat somewhere else, since we
  ! do not know anything about octedat structure here.
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr
  !
  !.. Formal Arguments ..
  type(octtree_type), intent(inout) :: octtree
  real(wp),              intent(in) :: bbox(1:nr,1:2)
  integer,              intent(out) :: istat
  !
  !.. Local Scalars ..
  integer :: l,k,maxoctdata,maxoctnodes
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octtree_init"
  !
continue
  !
  istat = OCT_STATUS_OK  ! no error yet
  !
  ! Make sure this octtree is empty, but do not erase octedat data
  !
  call octtree_dealloc(octtree,fals,istat)
  !
  if (istat /= OCT_STATUS_OK) return
  !
  ! Init the octroot.
  !
  maxoctdata  = -1       ! use default from octnode_mod
  maxoctnodes = 2**nr
  call octroot_init(octtree%octroot,maxoctdata,maxoctnodes,istat)
  !
  if (istat /= OCT_STATUS_OK) return
  !
  ! Init the initial octnode.
  !
  call octnode_init(octtree%octnode,octtree%octroot,bbox,istat) !,parent=null()
  !
  if (istat /= OCT_STATUS_OK) return
  !
  ! Push data into octtree.
  ! For a simple array src_xyzdata(3,npts):
  !   do k = 1,size(src_xyzdata,dim=2)
  !     call octnode_addxyz(octtree%octnode,octtree%octroot,k,istat)
  !     if (istat /= OCT_STATUS_OK) return
  !   end do
  ! However, we
  ! call ed_data%octtree_fill(octtree,
  !
end subroutine octtree_init
!
!###############################################################################
!
subroutine octtree_addxyz(octtree,ixyz_pt,xyz_pt,istat)
  !
  ! Add an xyz point to this octtree.
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr
  !
  !.. Formal Arguments ..
  type(octtree_type), intent(inout) :: octtree
  integer,               intent(in) :: ixyz_pt(:)
  real(wp),              intent(in) :: xyz_pt(1:nr)
  integer,              intent(out) :: istat
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octtree_addxyz"
  !
continue
  !
  call octnode_addxyz(octtree%octnode,octtree%octroot,octtree%octedat, &
                      ixyz_pt,xyz_pt,istat)
  !
end subroutine octtree_addxyz
!
!###############################################################################
!
subroutine octtree_find_containing_octnode(octtree,xyz_pt,containing_octnode, &
                                           istat)
  !
  ! Within this octtree, find the smallest octnode or sub-octnode whose bbox
  ! contains xyz_pt.  This node will have no children, and may or may not
  ! have data.
  !
  !.. Formal Arguments ..
  type(octtree_type),           intent(in) :: octtree
  real(wp),                     intent(in) :: xyz_pt(:)
  type(octnode_type), pointer, intent(out) :: containing_octnode
  integer,                     intent(out) :: istat
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octtree_find_containing_octnode"
  !
continue
  !
  ! Set default containing octnode
  !
  containing_octnode => null()
  !
  ! Search the octtree starting from the head octnode.
  !
  call octnode_find_containing_octnode(octtree%octnode, &
                                       octtree%octroot,octtree%octedat, &
                                       xyz_pt,containing_octnode,istat)
  !
end subroutine octtree_find_containing_octnode
!
!###############################################################################
!
subroutine octtree_find_nearest_xyz(octtree,xyz_pt, ixyz_pt,istat,tid,acc)
  !
  ! Find point in octtree that is nearest to desired xyz_pt.
  ! Return the index info used to retrieve the corresponding xyz point
  ! from routine octtree%octedat%evalxyz_pt(idx,xyz,istat).
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr
  !
  !.. Formal Arguments ..
  type(octtree_type), intent(in) :: octtree
  real(wp),           intent(in) :: xyz_pt(:)
  integer,           intent(out) :: ixyz_pt(:)
  integer,           intent(out) :: istat
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: tid
  integer, optional, intent(in) :: acc
  !
  !.. Local Scalars ..
  integer  :: best_octixyz
  real(wp) :: best_dist2
  real(wp) :: bbox_diag
  integer  :: ierr
  !
  !-DEBUG_BEG:
#ifdef DEBUG_OCTFIND
  integer l
#endif
  !-DEBUG_END:
  !
  !.. Local Arrays ..
  real(wp) :: xyz_bbox(1:nr,1:2)
  real(wp) :: bbox_size(1:nr)
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: best_ixyz(:)
  !
  !.. Local Pointers ..
  type(octnode_type), pointer :: containing_octnode
  type(octnode_type), pointer :: best_octnode
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octtree_find_nearest_xyz"
  !
continue
  !
#ifdef DEBUG_ON
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
  !
  ! Set default values, in case no point is found.
  !
  call debug_timer(start_accumulation,"Crap 1", &
                   acc=1,pid=tid)
#endif
  !
  istat = OCT_STATUS_NOT_FOUND
  containing_octnode => null()
  best_octnode => null()
  best_octixyz = 0
  best_dist2 = huge(best_dist2)
  ixyz_pt = 0
  !
  ! Allocate best_ixyz
  !
  allocate ( best_ixyz(octtree%octedat%maxindices) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"best_ixyz",1,__LINE__,__FILE__,ierr,error_message)
  !
#ifdef DEBUG_ON
  call debug_timer(stop_accumulation,"Crap 1", &
                   acc=1,pid=tid)
  !
  ! Find the smallest octnode containing xyz_pt.
  !
  call debug_timer(start_accumulation,"octnode_find_containing_octnode", &
                   acc=2,pid=tid)
#endif
  call octnode_find_containing_octnode(octtree%octnode,octtree%octroot, &
                                       octtree%octedat,xyz_pt, &
                                       containing_octnode,istat)
#ifdef DEBUG_ON
  call debug_timer(stop_accumulation,"octnode_find_containing_octnode", &
                   acc=2,pid=tid)
#endif
  !
  !
  !
  if (istat == OCT_STATUS_OK) then
    !
#ifdef DEBUG_ON
    call debug_timer(start_accumulation,"Crap 2", &
                     acc=3,pid=tid)
#endif
    !
    best_octnode => containing_octnode
    !
    ! Now have the smallest octnode containing xyz_pt.
    ! Nearest point could be in the containing_octnode or an adjacent octnode
    ! (either a sibling octnode or a distance cousin from the octtree).
    !
    ! Know there is at least one point in containing_octnode or its siblings.
    ! Its distance is within twice the diagonal of containing_octnode bbox.
    ! Therefore, create a shpere (or cubic bounding box) of r=2*diagonal
    ! that is centered at xyz_pt, and search all octnodes that intersect it
    ! to see if there is a closer point.
    !
    bbox_size = containing_octnode%bbox(:,2) - containing_octnode%bbox(:,1)
    bbox_diag = sqrt(sum(bbox_size * bbox_size))
    xyz_bbox(:,1) = xyz_pt(:) - 2*bbox_diag
    xyz_bbox(:,2) = xyz_pt(:) + 2*bbox_diag
    !
#ifdef DEBUG_ON
    call debug_timer(stop_accumulation,"Crap 2", &
                     acc=3,pid=tid)
    !
    ! Search octtree octnodes for point nearest to xyz_pt.
    !
    call debug_timer(start_accumulation,"octnode_find_nearest_xyz", &
                     acc=4,pid=tid)
#endif
    call octnode_find_nearest_xyz(octtree%octnode,octtree%octroot, &
           octtree%octedat,xyz_pt,xyz_bbox, &
           best_octnode,best_octixyz,best_ixyz,best_dist2,istat)
    if (istat == OCT_STATUS_OK) then
      ixyz_pt = best_ixyz
    end if
#ifdef DEBUG_ON
    call debug_timer(stop_accumulation,"octnode_find_nearest_xyz", &
                     acc=4,pid=tid)
#endif
    !
  end if
  !
#ifdef DEBUG_ON
  call debug_timer(start_accumulation,"Crap 3", &
                   acc=5,pid=tid)
#endif
  !
#ifdef DEBUG_OCTFIND
  !-DEBUG_BEG:
  if (associated(containing_octnode)) then
    write (iout,10) "Containing",containing_octnode%octid
    write (iout,11) (containing_octnode%bbox(l,1), l=1,nr)
    write (iout,12) (containing_octnode%bbox(l,2), l=1,nr)
  else
    write (iout,13) "Containing"
  end if
  !
  if (associated(best_octnode) then
    write (iout,10) "Best",best_octnode%octid
    write (iout,11) (best_octnode%bbox(l,1), l=1,nr)
    write (iout,12) (best_octnode%bbox(l,2), l=1,nr)
  else
    write (iout,13) "Best"
  end if
  !-DEBUG_END:
#endif
  !
  if (allocated(best_ixyz)) then
    deallocate ( best_ixyz , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"best_ixyz",2,__LINE__,__FILE__,ierr,error_message)
  end if
#ifdef DEBUG_ON
  call debug_timer(stop_accumulation,"Crap 3", &
                   acc=5,pid=tid)
#endif
  !
  ! Format Statements
  !
  10 format (2x,a10," octnode = ",i0)
  11 format (4x,"bbox_min",8x," = ",3(i0,:,1x))
  12 format (4x,"bbox_max",8x," = ",3(i0,:,1x))
  13 format (2x,a10," octnode = null")
  !
end subroutine octtree_find_nearest_xyz
!
!###############################################################################
!
subroutine octtree_print(octtree,filename,istat)
  !
  ! Print the octtree to file for debugging.
  !
  !.. Formal Arguments ..
  type(octtree_type), intent(in) :: octtree
  character(len=*),   intent(in) :: filename
  integer,           intent(out) :: istat
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octtree_print"
  !
continue
  !
  istat = OCT_STATUS_OK
  !
  call octnode_print(octtree%octnode,octtree%octedat,istat,filename)
  !
end subroutine octtree_print
!
!###############################################################################
!
#ifndef DEBUG_ON
#undef PROFILE_ON
#endif
!
#include "Functions/debug_timer.f90"
!
!###############################################################################
!
end module octtree_mod
