module octnode_mod
  !
  ! Base module for performing an octtree based search through xyz data.
  !
  ! We use octtree vs octree, since it better aligns with octroot and octnode.
  !
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
  !.. Use Statements ..
  use module_kind_types, only : wp,iout,alloc_error,error_message,half,io_error
  use octstat_mod
  use octroot_mod
  use octedat_mod
  !
  implicit none
  !
  private
  !
  public :: octnode_dealloc
  public :: octnode_init
  public :: octnode_addxyz
  public :: octnode_find_containing_octnode
  public :: octnode_find_nearest_xyz
  public :: octnode_print
  !
  ! -------------------------------------------------------------------------
  ! Octtree node type.  Stuff unique to each branch/leaf octnode.
  !   octid                   Unique identifier for each octnode.
  !   parent                  Parent oct node (null if head).
  !   child(1:MAXOCTNODES)    Child  oct nodes.
  !   bbox(1:nr,1:2)          The (min:max) bounding box of octnode.
  !   nxyzdata                Number of active elements in ixyzdata.
  !   ixyzdata(1:MAXINDICES,1:MAXOCTDATA)
  !                           Integer index to evaluate xyzdata(1:nr)
  !                           from whatever array/struct is in octedat.
  ! Note: child is a pointer because not all compilers support allocatable
  !       recursive components yet.
  type, public :: octnode_type
#ifdef OCTID
    integer                         :: octid = -1
#endif
    type(octnode_type), pointer     :: parent => null()
    type(octnode_type), pointer     :: child(:) => null()
    real(wp),           allocatable :: bbox(:,:)
    integer                         :: nxyzdata = 0
    integer,            allocatable :: ixyzdata(:,:)
 !contains
 !  final :: octnode_final
  end type octnode_type
  !
contains
!
!###############################################################################
!
subroutine octnode_final(this)
  !
  ! Finalize the octnode_type.
  !
  !.. Formal Arguments ..
  type(octnode_type) :: this
  !
  !.. Local Scalars ..
  integer :: istat
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octnode_final"
  !
continue
  !
  !-DEBUG_BEG:
#ifdef DEBUG_OCTDEALLOC
  write (iout,*) "DEBUG[octnode_final] Call octnode_dealloc."
  flush (iout)
#endif
  !-DEBUG_END:
  !
  call octnode_dealloc(this,istat)
  !
  !-DEBUG_BEG:
#ifdef DEBUG_OCTDEALLOC
  write (iout,*) "DEBUG[octnode_final] Leaving subroutine."
  flush (iout)
#endif
  !-DEBUG_END:
  !
end subroutine octnode_final
!
!###############################################################################
!
recursive &
subroutine octnode_dealloc(octnode,istat)
  !
  ! Deallocate the octnode and all child octnodes.
  !
  !.. Formal Arguments ..
  class(octnode_type) :: octnode
  integer             :: istat
  !
  !.. Local Scalars ..
  integer :: ierr,n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octnode_dealloc"
  !
continue
  !
  !-DEBUG_BEG:
#ifdef DEBUG_OCTDEALLOC
  write (iout,*) "DEBUG[octnode_dealloc] Deallocating octnode components."
  flush (iout)
#ifdef OCTID
  write (iout,*) "DEBUG[octnode_dealloc] octid = ",octnode%octid
  flush (iout)
#endif
#endif
  !-DEBUG_END:
  !
  istat = OCT_STATUS_OK
  !
  if (associated(octnode%parent)) nullify(octnode%parent)
  !
  if (associated(octnode%child)) then
    !
    do n = 1,size(octnode%child)
      call octnode_dealloc(octnode%child(n),ierr)
      if (ierr /= OCT_STATUS_OK) istat = ierr
    end do
    !
    deallocate ( octnode%child , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"octnode%child",2,__LINE__,__FILE__, &
                     ierr,error_message)
    !
  end if
  !
  if (allocated(octnode%bbox)) then
    deallocate ( octnode%bbox , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"octnode%bbox",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(octnode%ixyzdata)) then
    deallocate ( octnode%ixyzdata , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"octnode%ixyzdata",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  octnode%nxyzdata = 0
  !
  !-DEBUG_BEG:
#ifdef DEBUG_OCTDEALLOC
  write (iout,*) "DEBUG[octnode_dealloc] Leaving subroutine."
  flush (iout)
#endif
  !-DEBUG_END:
  !
end subroutine octnode_dealloc
!
!###############################################################################
!
subroutine octnode_init(octnode,octroot,bbox,istat,parent)
  !
  ! Initialize non-default components of an octnode.
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Formal Arguments ..
  type(octnode_type), intent(inout) :: octnode
  type(octroot_type), intent(inout) :: octroot
  real(wp),              intent(in) :: bbox(:,:)
  integer,              intent(out) :: istat
  !
  !.. Optional Arguments ..
  type(octnode_type), optional, target, intent(in) :: parent
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octnode_init"
  !
continue
  !
  istat = OCT_STATUS_OK
  !
  ! Add this octnode to the total count.
  !
  octroot%totoctnode = octroot%totoctnode + 1
  !
#ifdef OCTID
  ! If tracking octnode IDs then assign one.
  !
  octnode%octid = octroot%totoctnode
#endif
  !
  ! Associate the parent component to the parent arg or else null.
  !
  if (present(parent)) then
    octnode%parent => parent
  else
    nullify(octnode%parent)
  end if
  !
  ! Allocate space for and copy the octnode bounding box.
  !
  allocate ( octnode%bbox(1:nr,1:2) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"octnode%bbox",1,__LINE__,__FILE__,ierr,error_message)
  !
  if (ierr /= 0) then
    istat = OCT_STATUS_MEMORY_ERROR
  else
    octnode%bbox = bbox
  end if
  !
end subroutine octnode_init
!
!###############################################################################
!
recursive &
subroutine octnode_addxyz(octnode,octroot,octedat,ixyz_pt,xyz_pt,istat)
  !
  ! Add an xyz point to this octnode or one of the subnodes.
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Formal Arguments ..
  type(octnode_type), intent(inout) :: octnode
  type(octroot_type), intent(inout) :: octroot
  class(octedat_type),   intent(in) :: octedat
  integer,               intent(in) :: ixyz_pt(:)
  real(wp),              intent(in) :: xyz_pt(1:nr)
  integer,              intent(out) :: istat
  !
  !.. Local Scalars ..
  integer :: ierr,l,n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octnode_addxyz"
  !
continue
  !
  ! If xyz_pt not in this octnode, return with error.
  !
  istat = OCT_STATUS_NOT_FOUND
  do l = 1,nr
    if (xyz_pt(l) < octnode%bbox(l,1)) return
    if (xyz_pt(l) > octnode%bbox(l,2)) return
  end do
  !
  istat = OCT_STATUS_OK
  !
  !-DEBUG_BEG:
#ifdef DEBUG_OCTFILL
  write (iout,10) associated(octnode%child), allocated(octnode%ixyzdata), &
                  octnode%nxyzdata, octroot%maxoctdata
#ifdef OCTID
  write (iout,11) octnode%octid
#endif
  write (iout,12)
  flush (iout)
#endif
  !-DEBUG_END:
  !
  if (allocated(octnode%ixyzdata)) then
    !
    ! This octnode contains data.
    !
    if (octnode%nxyzdata == octroot%MAXOCTDATA) then
      !
      ! octnode data is full.  Create children.  Push exisiting data to them.
      ! Push of xyz_pt will occur in next if-block.
      !
      call octnode_makechildren(octnode,octroot,octedat,istat)
      if (istat /= OCT_STATUS_OK) return
      !
    end if
  end if
  !
  if (associated(octnode%child)) then
    !
    ! This octnode contains children.
    ! Send xyz_pt to each child octnode until we find the right one.
    !
    !-DEBUG_BEG:
#ifdef DEBUG_OCTFILL
    write (iout,13)
    flush (iout)
#endif
    !-DEBUG_END:
    !
    do n = 1,octroot%MAXOCTNODES
      !
      !-DEBUG_BEG:
#ifdef DEBUG_OCTFILL
#ifdef OCTID
      write (iout,14) octnode%octid, n, octnode%child(n)%octid
      flush (iout)
#endif
#endif
      !-DEBUG_END:
      !
      call octnode_addxyz(octnode%child(n),octroot,octedat, &
                          ixyz_pt,xyz_pt,istat)
      if (istat /= OCT_STATUS_NOT_FOUND) exit
      !
    end do
    !
  else
    !
    ! This octnode contains no children.  Check for data.
    !
    if (.not. allocated(octnode%ixyzdata)) then
      !
      ! This octnode contains no data yet.  Allocate space for it.
      !
      allocate ( octnode%ixyzdata(octedat%MAXINDICES,octroot%MAXOCTDATA) , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"octnode%ixyzdata",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      octnode%nxyzdata = 0
      !
    end if
    !
    ! This octnode contains data array with room for xyz_pt.
    ! Add xyz_pt and increment octnode data count.
    !
    if (octnode%nxyzdata < octroot%MAXOCTDATA) then
      !
      octnode%nxyzdata = octnode%nxyzdata + 1
      octnode%ixyzdata(:,octnode%nxyzdata) = ixyz_pt
      octroot%totoctdata = octroot%totoctdata + 1
      !
      !-DEBUG_BEG:
#ifdef DEBUG_OCTFILL
#ifdef OCTID
      write (iout,15) octnode%octid
      flush (iout)
#endif
#endif
      !-DEBUG_END:
      !
    else
      !
      ! Should never get here if logic above is correct for
      ! handling octnode data is full.
      !
      write (iout,16)
      stop
      !
    end if
    !
  end if
  !
  !-DEBUG_BEG:
#ifdef DEBUG_OCTFILL
  write (iout,17)
  flush (iout)
#endif
  !-DEBUG_END:
  !
  ! Format Statement
  !
  10 format ("DEBUG_BEG[octnode_addxyz]",/, &
             " associated(octnode%child)     = ",l1,/, &
             "  allocated(octnode%ixyzdata)  = ",l1,/, &
             "            octnode%nxyzdata   = ",i0,/, &
             "            octroot%maxoctdata = ",i0)
  11 format ("            octnode%octid      = ",i0)
  12 format ("DEBUG_END[octnode_addxyz]")
  13 format ("DEBUG[octnode_addxyz]  Send xyz to each child octnode.")
  14 format (1x,"DEBUG[octnode_addxyz]  Send xyz from parent/to child ", &
                i4,"(",i2,")/",i4)
  15 format (1x,"DEBUG[octnode_addxyz]  xyz added to octnode ",i4)
  16 format ("Octtree Error:  Overwrote data buffer")
  17 format ("DEBUG[octnode_addxyz]  Leaving octnode_addxyz.")
  !
end subroutine octnode_addxyz
!
!###############################################################################
!
subroutine octnode_makechildren(octnode,octroot,octedat,istat)
  !
  ! Create octnode children and push exisiting octnode data down to them.
  ! Assumes children do not already exist.
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Formal Arguments ..
  type(octnode_type), intent(inout) :: octnode
  type(octroot_type), intent(inout) :: octroot
  class(octedat_type),   intent(in) :: octedat
  integer,              intent(out) :: istat
  !
  !.. Local Scalars ..
  integer :: ierr,l,k,n,oct
  real(wp):: boxmid
  !
  !.. Local Arrays ..
  integer  :: bit(1:3)
  real(wp) :: bbox(1:nr,1:2)
  real(wp) :: subbbox(1:nr,1:2)
  real(wp) :: xyz_pt(1:nr)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: xyz_pts(:,:)
  integer,  allocatable :: ixyz_pt(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octnode_makechildren"
  !
continue
  !
  istat = OCT_STATUS_OK
  !
  !-DEBUG_BEG:
#ifdef DEBUG_OCTFILL
  write (iout,10)
#endif
  !-DEBUG_END:
  !
  ! Allocate space for the child nodes.
  !
  allocate ( octnode%child(1:octroot%MAXOCTNODES) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"octnode%child",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Get current bbox so we can subdivide.
  !
  bbox(:,:) = octnode%bbox(:,:)
  !
  !-DEBUG_BEG:
#ifdef DEBUG_OCTFILL
  write (iout,11)
#endif
  !-DEBUG_END:
  !
  ! Init the empty child nodes with parent and bbox info.
  !
  do n = 1,octroot%MAXOCTNODES
    bit(:) = 0
    oct = n
    if (oct > 4) then
      oct = oct - 4 ; bit(3) = 1
    end if
    if (oct > 2) then
      oct = oct - 2 ; bit(2) = 1
    end if
    if (oct > 1) then
      oct = oct - 1 ; bit(1) = 1
    end if
    do l = 1,nr
      boxmid = ( bbox(l,1) + bbox(l,2) ) * half
      if (bit(l) == 0) then
        subbbox(l,1) = bbox(l,1)
        subbbox(l,2) = boxmid
      else
        subbbox(l,1) = boxmid
        subbbox(l,2) = bbox(l,2)
      end if
    end do
    !
    ! Init this child node with parent and bbox info.
    !
    call octnode_init(octnode%child(n),octroot,subbbox,istat, &
                      parent=octnode)
    if (istat /= OCT_STATUS_OK) return
    !
  end do ! n
  !
  if (allocated(octnode%ixyzdata)) then
    !
    ! Push existing octnode data down to children.
    !
    ! Need to get xyz data for all points in ixyzdata.
    !
    allocate ( xyz_pts(nr,octnode%nxyzdata) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_pts",1,__LINE__,__FILE__,ierr,error_message)
    !
    call octedat%evalxyz_pts(octnode%nxyzdata,octnode%ixyzdata, &
                             xyz_pts,istat)
    !
    if (istat /= OCT_STATUS_OK) return
    !
    ! Now push each point down to child octnodes.
    !
    do k = 1,octnode%nxyzdata
      !
      ixyz_pt = octnode%ixyzdata(:,k)
      xyz_pt = xyz_pts(:,k)
      !
      do n = 1,octroot%MAXOCTNODES
        !
        call octnode_addxyz(octnode%child(n),octroot,octedat, &
                            ixyz_pt,xyz_pt,istat)
        !
        !-DEBUG_BEG:
#ifdef DEBUG_OCTFILL
        ! Print debug info when making children.
        if (k == 1) then
          if (n == 1) then
            write (iout,12) ixyz_pt
            write (iout,13) xyz_pt
          end if
          write (iout,14) "k,n,status",k,n,istat
          write (iout,14) "bboxmin",(octnode%child(n)%bbox(l,1), l=1,nr)
          write (iout,14) "bboxmax",(octnode%child(n)%bbox(l,2), l=1,nr)
          if (n == octroot%MAXOCTNODES) then
            write (iout,15)
          end if
        end if
#endif
        !-DEBUG_END:
        !
        if (istat == OCT_STATUS_OK) then
          ! Data pushed to child.  No need to try other children.
          exit  ! n-loop
        end if
        !
        if (istat /= OCT_STATUS_NOT_FOUND) then
          ! Error while pushing data to child.  Abort with status.
          return
        end if
        !
      end do  ! n-loop
      !
      !-DEBUG_BEG:
#ifdef DEBUG_OCTFILL
      if (istat /= OCT_STATUS_OK) then
        write (iout,16) octnode%nxyzdata,k,n,istat
      end if
#endif
      !-DEBUG_END:
      !
    end do ! k-loop
    !
    ! Deallocate temp index array for a single point.
    !
    deallocate ( xyz_pts , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_pts",2,__LINE__,__FILE__,ierr,error_message)
    !
    ! Deallocate existing octnode data.
    !
    deallocate ( octnode%ixyzdata , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"octnode%ixyzdatas",2,__LINE__,__FILE__, &
                     ierr,error_message)
    !
    ! Remove the existing octnode data from the global count.
    !
    octroot%totoctdata = octroot%totoctdata - octroot%MAXOCTDATA
    !
    ! Clear the octnode data count.
    !
    octnode%nxyzdata = 0
    !
  end if
  !
  ! Format Statements
  !
  10 format ("DEBUG[octnode_makechildren] Alloc child nodes.")
  11 format ("DEBUG[octnode_makechildren] Init  child nodes.")
  12 format ("DEBUG_BEG[octnode_makechildren]",/, &
             2x,"Push existing octnode data to children.",/, &
             2x,"ixyzdata = ",100(i0,:,1x))
  13 format (3x,"xyzdata = ",100(es14.6,:,1x))
  14 format (2x,a10," = ",3(i0,:,1x))
  15 format ("DEBUG_END[octnode_makechildren]")
  16 format ("DEBUG_BEG[octnode_makechildren]",/, &
             2x,"Should never get here!",/, &
             2x,"nxyzdata   = ",i0,/, &
             2x,"k,n,status = ",i0,", ",i0,", ",i0,/, &
             "DEBUG_END[octnode_makechildren]")
  !
end subroutine octnode_makechildren
!
!###############################################################################
!
recursive &
subroutine octnode_find_containing_octnode(octnode,octroot,octedat,xyz_pt, &
                                           containing_octnode,istat)
  !
  ! Within this octnode, find the smallest octnode or sub-octnode whose bbox
  ! contains xyz_pt.  This node will have no children, and may or may not
  ! have data.
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Formal Arguments ..
  type(octnode_type), target,   intent(in) :: octnode
  type(octroot_type),           intent(in) :: octroot
  class(octedat_type),          intent(in) :: octedat
  real(wp),                     intent(in) :: xyz_pt(:)
  type(octnode_type), pointer, intent(out) :: containing_octnode
  integer,                     intent(out) :: istat
  !
  !.. Local Scalars ..
  integer :: l,n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octnode_find_containing_octnode"
  !
continue
  !
  ! Set default containing octnode
  !
  containing_octnode => null()
  !
  ! If xyz_pt not within this octnode, return with error.
  !
  istat = OCT_STATUS_NOT_FOUND
  do l = 1,nr
    if (xyz_pt(l) < octnode%bbox(l,1)) return
    if (xyz_pt(l) > octnode%bbox(l,2)) return
  end do
  !
  if (associated(octnode%child)) then
    !
    ! This octnode has subnodes.  Search each one until a match is found.
    !
    do n = 1,octroot%MAXOCTNODES
      !
      call octnode_find_containing_octnode(octnode%child(n),octroot,octedat, &
                                           xyz_pt,containing_octnode,istat)
      if (istat == OCT_STATUS_OK) return
      !
    end do
    !
  else
    !
    ! This is the smallest octnode containing xyz_pt.
    !
    istat = OCT_STATUS_OK
    containing_octnode => octnode
    !
    !-DEBUG_BEG:
    !-DEBUG_CMT:  Increasing the container may be overly cautious
    !-DEBUG_CMT:  and increase cpu time.
#ifdef OCTBIGGERCONTAINER
    !
    ! If this octnode contains no data, set containing_octnode to parent.
    !
    if (.not. allocated(octnode%ixyzdata)) containing_octnode => octnode%parent
    !
#endif
    !-DEBUG_END:
    !
  end if
  !
end subroutine octnode_find_containing_octnode
!
!###############################################################################
!
recursive &
subroutine octnode_find_nearest_xyz(octnode,octroot,octedat,xyz_pt,xyz_bbox, &
                                    best_octnode,best_octixyz,best_ixyz, &
                                    best_dist2,istat)
  !
  ! Find point in octnode or subnodes that is nearest to desired xyz_pt.
  ! Only searches octnodes that intersect xyz_bbox that surrounds xyz_pt.
  !
  ! istat has these possible return values:
  !   OCT_STATUS_OK          Closer match to xyz_pt was found.
  !   OCT_STATUS_NOT_FOUND   Closer match to xyz_pt was not found.
  !   OCT_STATUS_DATA_ERROR  Error evaluating data within the octtree.
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Formal Arguments ..
  type(octnode_type), target,     intent(in) :: octnode
  type(octroot_type),             intent(in) :: octroot
  class(octedat_type),            intent(in) :: octedat
  real(wp),                       intent(in) :: xyz_pt(:)
  real(wp),                       intent(in) :: xyz_bbox(:,:)
  type(octnode_type), pointer, intent(inout) :: best_octnode
  integer,                     intent(inout) :: best_octixyz
  integer,                     intent(inout) :: best_ixyz(:)
  real(wp),                    intent(inout) :: best_dist2
  integer,                       intent(out) :: istat
  !
  !.. Local Scalars ..
  integer  :: ierr,l,k,n,octnode_status
  real(wp) :: dist_to_bbox,dist2
  !
  !.. Local Arrays ..
  real(wp) :: dxyz(1:nr)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: xyz_pts(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octnode_find_nearest_xyz"
  !
continue
  !
  ! Default status, unless a better match is found.
  !
  istat = OCT_STATUS_NOT_FOUND
  !
  ! If any min/max edge of xyz_bbox is right/left of max/min of octnode%bbox,
  ! then boxes to not overlap.
  !
  do l = 1,nr
    if (xyz_bbox(l,1) > octnode%bbox(l,2)) return
    if (xyz_bbox(l,2) < octnode%bbox(l,1)) return
  end do
  !
  !-COMMENT_BEG:
  !-COMMENT_CMT: This check may or may not improve search time.
  !-COMMENT_CMT: It is not necessary for the algorithm to work.
  !
  ! If the distance^2 from xyz_pt to the nearest edge of the bbox exceeds
  ! best_dist2, then any point in this octnode will exceed best_dist2.
  !
  dist_to_bbox = huge(dist_to_bbox)
  !
  do l = 1,nr
    dist_to_bbox =  min( dist_to_bbox, abs( xyz_pt(l) - octnode%bbox(l,1) ) , &
                                       abs( xyz_pt(l) - octnode%bbox(l,2) ) )
  end do
  !
  if (dist_to_bbox**2 > best_dist2) return
  !-COMMENT_END:
  !
  if (associated(octnode%child)) then
    ! This octnode has subnodes.  Test each one for a closer match to xyz_pt.
    ! If any subnode is better, indicate with status.
    !
    do n = 1,octroot%MAXOCTNODES
      !
      call octnode_find_nearest_xyz(octnode%child(n),octroot,octedat,xyz_pt, &
                                    xyz_bbox,best_octnode,best_octixyz, &
                                    best_ixyz,best_dist2,octnode_status)
      if (octnode_status == OCT_STATUS_OK) istat = OCT_STATUS_OK
      !
    end do
    !
  else if (allocated(octnode%ixyzdata)) then
    !
    ! This octnode has ixyzdata.  Test each point for a closer match to xyz_pt.
    ! If any is better, indicate with status and update best* variables.
    !
    ! Get xyz for all points in ixyzdata.
    !
    allocate ( xyz_pts(1:nr,1:octnode%nxyzdata) , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_pts",1,__LINE__,__FILE__,ierr,error_message)
    !
    call octedat%evalxyz_pts(octnode%nxyzdata,octnode%ixyzdata,xyz_pts,istat)
    !
    if (istat == OCT_STATUS_OK) then
      !
      do k = 1,octnode%nxyzdata
        !
        dxyz = xyz_pt(:) - xyz_pts(:,k)
        dist2 = dot_product(dxyz,dxyz)
        !
        if (dist2 < best_dist2) then
          best_octnode => octnode
          best_octixyz = k                     ! index within octnode data
          best_ixyz    = octnode%ixyzdata(:,k) ! index within octedat data
          best_dist2   = dist2                 ! distance squared
          istat = OCT_STATUS_OK
        end if
        !
      end do
      !
      deallocate ( xyz_pts , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"xyz_pts",2,__LINE__,__FILE__,ierr,error_message)
      !
    end if
    !
  end if
  !
end subroutine octnode_find_nearest_xyz
!
!###############################################################################
!
recursive &
subroutine octnode_print(octnode,octedat,istat,filename)
  !
  ! Print the octnode to file for debugging.
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Formal Arguments ..
  type(octnode_type),  intent(in) :: octnode
  class(octedat_type), intent(in) :: octedat
  integer,            intent(out) :: istat
  !
  !.. Optional Arguments ..
  character(len=*), optional, intent(in) :: filename
  !
  !.. Local Scalars ..
  integer :: ierr,l,k,n,io
  character(len=100) :: fname
  !
  !.. Local Arrays ..
  ! Note:  xyz_pts(1:nr,1:octnode%nxyzdata)
  real(wp), allocatable :: xyz_pts(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octnode_print"
  !
continue
  !
  istat = OCT_STATUS_OK
  !
  fname = "octnode_print.dat"
  if (present(filename)) then
    fname = trim(adjustl(filename))
  end if
  !
  open (newunit=io,file=trim(fname),form="formatted",status="replace", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  !
#ifdef OCTID
  ! Print octnode id.
  !
  write (io,10) octnode%octid
#endif
  !
  ! Print parent id.
  ! write (io,*) "parent = ", ???   ! how to handle root node?
  !
  ! Print bounding box.
  !
  if (allocated(octnode%bbox)) then
    write (io,11) "bboxmin",1,(octnode%bbox(l,1), l=1,nr)
    write (io,11) "bboxmax",2,(octnode%bbox(l,2), l=1,nr)
  else
    write (io,12) "bboxmin"
    write (io,12) "bboxmax"
  end if
  !
  ! Print data.
  !
  write (io,13) octnode%nxyzdata
  !
  if (allocated(octnode%ixyzdata) .and. octnode%nxyzdata > 0) then
    !
    allocate ( xyz_pts(nr,octnode%nxyzdata) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_pts",1,__LINE__,__FILE__,ierr,error_message)
    !
    call octedat%evalxyz_pts(octnode%nxyzdata,octnode%ixyzdata, &
                             xyz_pts,istat)
    !
    do k = 1,octnode%nxyzdata
      write (io,*) "           ", &
                   (octnode%ixyzdata(l,k), l=1,size(octnode%ixyzdata,dim=1)), &
                   (xyz_pts(l,k), l=1,nr)
    end do
    !
    deallocate ( xyz_pts , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_pts",2,__LINE__,__FILE__,ierr,error_message)
    !
  end if
  !
  ! Print children.
  !
  if (associated(octnode%child)) then
    write (io,14)
#ifdef OCTID
    do n = 1,size(octnode%child)
      write (io,15) octnode%child(n)%octid
    end do
#endif
    write (io,*)
    do n = 1,size(octnode%child)
      call octnode_print(octnode%child(n),octedat,istat)
    end do
  else
    write (io,16)
  end if
  !
  close (unit=io,iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Format Statements
  !
  10 format ("octid    = ",i0)
  11 format (a8," = ",i1,3(2x,i0,:))
  12 format (a8," = undefined")
  13 format ("nxyzdata = ",i0)
  14 format ("Contains children:")
  15 format (2x,i0)
  16 format ("Contains no children.",/)
  !
end subroutine octnode_print
!
!###############################################################################
!
end module octnode_mod
