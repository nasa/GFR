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
  use module_kind_types, only : wp, iout, alloc_error, error_message
  use geovar, only : nr
  use octstat_mod
  use octroot_mod
  use octedat_mod
  implicit none
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
    integer                                         :: octid    = -1
#endif
    type(octnode_type) ,pointer                     :: parent   => null()
    type(octnode_type) ,pointer     ,dimension(:)   :: child    => null()
    real(wp)           ,allocatable ,dimension(:,:) :: bbox
    integer                                         :: nxyzdata =  0
    integer            ,allocatable ,dimension(:,:) :: ixyzdata
! contains
!   final :: octnode_final
  end type octnode_type
!
contains
!
!##########################################################################
  subroutine octnode_final(this)
  !
  ! Finalize the octnode_type.
  !
  implicit none
  !
  !.. Formal Arguments ..
  type(octnode_type) :: this
  !
  !.. Local Scalars ..
  integer :: status
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octnode_final"
  !
continue
  !
!-DEBUG_BEG:
#ifdef DEBUG_OCTDEALLOC
  write(iout,*) 'DEBUG[octnode_final] Call octnode_dealloc.'
  flush(iout)
#endif
!-DEBUG_END:
  call octnode_dealloc(this,status)
!-DEBUG_BEG:
#ifdef DEBUG_OCTDEALLOC
  write(iout,*) 'DEBUG[octnode_final] Leaving subroutine.'
  flush(iout)
#endif
!-DEBUG_END:
  return
  end subroutine octnode_final
!##########################################################################
  recursive &
  subroutine octnode_dealloc(octnode,status)
  !
  ! Deallocate the octnode and all child octnodes.
  !
  implicit none
  !
  !.. Formal Arguments ..
  class(octnode_type) :: octnode
  integer             :: status
  !
  !.. Local Scalars ..
  integer :: ierr, n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octnode_dealloc"
  !
continue
!-DEBUG_BEG:
#ifdef DEBUG_OCTDEALLOC
  write(iout,*) 'DEBUG[octnode_dealloc] Deallocating octnode components.'
  flush(iout)
#ifdef OCTID
  write(iout,*) 'DEBUG[octnode_dealloc] octid = ',octnode%octid
  flush(iout)
#endif
#endif
!-DEBUG_END:
  !
  status = OCT_STATUS_OK
  if ( associated(octnode%parent) ) nullify(octnode%parent)
  if ( associated(octnode%child) ) then
    do n=1,size(octnode%child)
      call octnode_dealloc(octnode%child(n),ierr)
      if ( ierr .ne. OCT_STATUS_OK ) status = ierr
    end do
    deallocate( octnode%child , stat=ierr , errmsg=error_message )
    if ( ierr .ne. 0 ) status = OCT_STATUS_MEMORY_ERROR
    call alloc_error(pname,"octnode%child",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  if ( allocated(octnode%bbox) ) then
    deallocate( octnode%bbox , stat=ierr , errmsg=error_message )
    if ( ierr .ne. 0 ) status = OCT_STATUS_MEMORY_ERROR
    call alloc_error(pname,"octnode%bbox",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  if ( allocated(octnode%ixyzdata) ) then
    deallocate( octnode%ixyzdata , stat=ierr , errmsg=error_message )
    if ( ierr .ne. 0 ) status = OCT_STATUS_MEMORY_ERROR
    call alloc_error(pname,"octnode%ixyzdata",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  octnode%nxyzdata = 0
  !
!-DEBUG_BEG:
#ifdef DEBUG_OCTDEALLOC
  write(iout,*) 'DEBUG[octnode_dealloc] Leaving subroutine.'
  flush(iout)
#endif
!-DEBUG_END:
  return
  end subroutine octnode_dealloc
!##########################################################################
  subroutine octnode_init(octnode,octroot,bbox,status,parent)
  !
  ! Initialize non-default components of an octnode.
  !
  implicit none
  !
  !.. Formal Arguments ..
  type(octnode_type)                    ,intent(inout) :: octnode
  type(octroot_type)                    ,intent(inout) :: octroot
  real(wp)              ,dimension(:,:) ,intent(in)    :: bbox
  integer                               ,intent(out)   :: status
  type(octnode_type)  ,optional ,target ,intent(in)    :: parent
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octnode_init"
  !
continue
  !
  status = OCT_STATUS_OK
  !
  ! Add this octnode to the total count.
  octroot%totoctnode = octroot%totoctnode + 1
  !
#ifdef OCTID
  ! If tracking octnode id's then assign one.
  octnode%octid = octroot%totoctnode
#endif
  !
  ! Associate the parent component to the parent arg or else null.
  if ( present( parent ) ) then
    octnode%parent => parent
  else
    nullify(octnode%parent)
  end if
  !
  ! Allocate space for and copy the octnode bounding box.
  allocate( octnode%bbox(1:nr,1:2) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"octnode%bbox",1,__LINE__,__FILE__,ierr, &
                   error_message)
  if ( ierr .ne. 0 ) then
    status = OCT_STATUS_MEMORY_ERROR
  else
    octnode%bbox = bbox
  end if
  !
  return
  end subroutine octnode_init
!##########################################################################
  recursive &
  subroutine octnode_addxyz(octnode,octroot,octedat,ixyz_pt,xyz_pt,status)
  !
  ! Add an xyz point to this octnode or one of the subnodes.
  !
  implicit none
  !
  !.. Formal Arguments ..
  type(octnode_type)                    ,intent(inout) :: octnode
  type(octroot_type)                    ,intent(inout) :: octroot
  class(octedat_type)                   ,intent(in)    :: octedat
  integer              ,dimension(:)    ,intent(in)    :: ixyz_pt
  real(wp)             ,dimension(1:nr) ,intent(in)    :: xyz_pt
  integer                               ,intent(out)   :: status
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
  status = OCT_STATUS_NOT_FOUND
  do l = 1,nr
    if ( xyz_pt(l) .lt. octnode%bbox(l,1) ) return
    if ( xyz_pt(l) .gt. octnode%bbox(l,2) ) return
  end do
  !
  status = OCT_STATUS_OK
  !
!-DEBUG_BEG:
#ifdef DEBUG_OCTFILL
  write(iout,*) 'DEBUG_BEG[octnode_addxyz]'
  write(iout,*) ' associated(octnode%child)     = ',associated(octnode%child)
  write(iout,*) '  allocated(octnode%ixyzdata)  = ',allocated(octnode%ixyzdata)
  write(iout,*) '            octnode%nxyzdata   = ',          octnode%nxyzdata
  write(iout,*) '            octroot%maxoctdata = ',          octroot%maxoctdata
#ifdef OCTID
  write(iout,*) '            octnode%octid      = ',           octnode%octid
#endif
  write(iout,*) 'DEBUG_END[octnode_addxyz]'
  flush(iout)
#endif
!-DEBUG_END:
  !
  if ( allocated(octnode%ixyzdata) ) then
    ! This octnode contains data.
    if ( octnode%nxyzdata .eq. octroot%MAXOCTDATA ) then
      ! octnode data is full.  Create children.  Push exisiting data to them.
      ! Push of xyz_pt will occur in next if-block.
      call octnode_makechildren(octnode,octroot,octedat,status)
      if ( status .ne. OCT_STATUS_OK ) return
    end if
  end if
  !
  if ( associated(octnode%child) ) then
    ! This octnode contains children.
    ! Send xyz_pt to each child octnode until we find the right one.
!-DEBUG_BEG:
#ifdef DEBUG_OCTFILL
  write(iout,*) 'DEBUG[octnode_addxyz]  Send xyz to each child octnode.'
  flush(iout)
#endif
!-DEBUG_END:
    do n = 1,octroot%MAXOCTNODES
!-DEBUG_BEG:
#ifdef DEBUG_OCTFILL
#ifdef OCTID
  write(iout,'(x,a,i4,"(",i2,")/",i4)') &
    'DEBUG[octnode_addxyz]  Send xyz from parent/to child ', &
    octnode%octid, n, octnode%child(n)%octid
  flush(iout)
#endif
#endif
!-DEBUG_END:
      call octnode_addxyz(octnode%child(n),octroot,octedat, &
                          ixyz_pt,xyz_pt,status)
      if ( status .ne. OCT_STATUS_NOT_FOUND ) exit
    end do
  else
    ! This octnode contains no children.  Check for data.
    if ( .not. allocated(octnode%ixyzdata) ) then
      ! This octnode contains no data yet.  Allocate space for it.
      allocate( octnode%ixyzdata(octedat%MAXINDICES,octroot%MAXOCTDATA) , &
                stat=ierr , errmsg=error_message )
      call alloc_error(pname,"octnode%ixyzdata",1,__LINE__,__FILE__,ierr, &
                       error_message)
      if ( ierr .ne. 0 ) then
        status = OCT_STATUS_MEMORY_ERROR
        return
      end if
      octnode%nxyzdata = 0
    end if
    ! This octnode contains data array with room for xyz_pt.
    ! Add xyz_pt and increment octnode data count.
    if ( octnode%nxyzdata .lt. octroot%MAXOCTDATA ) then
      octnode%nxyzdata = octnode%nxyzdata + 1
      octnode%ixyzdata(:,octnode%nxyzdata) = ixyz_pt
      octroot%totoctdata = octroot%totoctdata + 1
!-DEBUG_BEG:
#ifdef DEBUG_OCTFILL
#ifdef OCTID
  write(iout,'(x,a,i4)') &
    'DEBUG[octnode_addxyz]  xyz added to octnode ', octnode%octid
  flush(iout)
#endif
#endif
!-DEBUG_END:
    else
      ! Should never get here if logic above is correct for
      ! handling octnode data is full.
      write(iout,*) 'Octtree Error:  Overwrote data buffer'
      stop
    end if
  end if
  !
!-DEBUG_BEG:
#ifdef DEBUG_OCTFILL
  write(iout,*) 'DEBUG[octnode_addxyz]  Leaving octnode_addxyz.'
  flush(iout)
#endif
!-DEBUG_END:
  return
  end subroutine octnode_addxyz
!##########################################################################
  subroutine octnode_makechildren(octnode,octroot,octedat,status)
  !
  ! Create octnode children and push exisiting octnode data down to them.
  ! Assumes children do not already exist.
  !
  implicit none
  !
  !.. Formal Arguments ..
  type(octnode_type)                    ,intent(inout) :: octnode
  type(octroot_type)                    ,intent(inout) :: octroot
  class(octedat_type)                   ,intent(in)    :: octedat
  integer                               ,intent(out)   :: status
  !
  !.. Local Parameters ..
  real(wp) ,parameter :: half = 0.5_wp
  !
  !.. Local Scalars ..
  integer ierr,l,k,n,oct
  real(wp):: boxmid
  !
  !.. Local Arrays ..
  integer ,dimension(1:3) :: bit
  real(wp),dimension(1:nr,1:2) :: bbox,subbbox
  real(wp),dimension(1:nr)             :: xyz_pt
  real(wp),dimension(:,:) ,allocatable :: xyz_pts
  integer ,dimension(:)   ,allocatable :: ixyz_pt
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octnode_makechildren"
  !
continue
  !
  status = OCT_STATUS_OK
  !
!-DEBUG_BEG:
#ifdef DEBUG_OCTFILL
  write(iout,*) 'DEBUG[octnode_makechildred] Alloc child nodes.'
#endif
!-DEBUG_END:
  !
  ! Allocate space for the child nodes.
  !
  allocate( octnode%child(1:octroot%MAXOCTNODES) , stat=ierr , &
            errmsg=error_message )
  call alloc_error(pname,"octnode%child",1,__LINE__,__FILE__,ierr, &
                   error_message)
  if ( ierr .ne. 0 ) then
    status = OCT_STATUS_MEMORY_ERROR
    return
  end if
  !
  ! Get current bbox so we can subdivide.
  !
  bbox(:,:) = octnode%bbox(:,:)
  !
!-DEBUG_BEG:
#ifdef DEBUG_OCTFILL
  write(iout,*) 'DEBUG[octnode_makechildred] Init  child nodes.'
#endif
!-DEBUG_END:
  !
  ! Init the empty child nodes with parent and bbox info.
  !
  do n = 1,octroot%MAXOCTNODES
    bit(:) = 0
    oct = n
    if ( oct .gt. 4 ) then
      oct = oct - 4 ; bit(3) = 1
    end if
    if ( oct .gt. 2 ) then
      oct = oct - 2 ; bit(2) = 1
    end if
    if ( oct .gt. 1 ) then
      oct = oct - 1 ; bit(1) = 1
    end if
    do l = 1,nr
      boxmid = ( bbox(l,1) + bbox(l,2) ) * half
      if ( bit(l) .eq. 0 ) then
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
    call octnode_init(octnode%child(n),octroot,subbbox,status, &
                      parent=octnode)
    if ( status .ne. OCT_STATUS_OK ) return
  end do ! n

  if ( allocated(octnode%ixyzdata) ) then
    ! Push existing octnode data down to children.
    !
    ! Need to get xyz data for all points in ixyzdata.
    allocate( xyz_pts(nr,octnode%nxyzdata) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_pts",1,__LINE__,__FILE__,ierr,error_message)
    if ( ierr .ne. 0 ) then
      status = OCT_STATUS_MEMORY_ERROR
      return
    endif
    call octedat%evalxyz_pts(octnode%nxyzdata,octnode%ixyzdata, &
                             xyz_pts,status)
    if ( status .ne. OCT_STATUS_OK ) return
    !
    ! Now push each point down to child octnodes.
    do k = 1,octnode%nxyzdata
      ixyz_pt = octnode%ixyzdata(:,k)
       xyz_pt = xyz_pts(:,k)
      do n = 1,octroot%MAXOCTNODES
        call octnode_addxyz(octnode%child(n),octroot,octedat, &
                            ixyz_pt,xyz_pt,status)
!
!-DEBUG_BEG:
#ifdef DEBUG_OCTFILL
  ! Print debug info when making children.
  if ( k .eq.  1 ) then
    if ( n .eq. 1 ) then
    write(*,*) 'DEBUG_BEG[octnode_makechildren]'
    write(*,*) '  Push existing octnode data to children.'
    write(*,*) '  ixyzdata = ',ixyz_pt
    write(*,*) '   xyzdata = ', xyz_pt
    end if
    write(*,*) '  k,n,status = ',k,n,status
    write(*,*) '  bboxmin    = ',(octnode%child(n)%bbox(l,1), l=1,nr)
    write(*,*) '  bboxmax    = ',(octnode%child(n)%bbox(l,2), l=1,nr)
    if ( n .eq. octroot%MAXOCTNODES ) then
    write(*,*) 'DEBUG_END[octnode_makechildren]'
    end if
  end if
#endif
!-DEBUG_END:
!
        if ( status .eq. OCT_STATUS_OK ) then
          ! Data pushed to child.  No need to try other children.
          exit  ! n-loop
        end if
        if ( status .ne. OCT_STATUS_NOT_FOUND ) then
          ! Error while pushing data to child.  Abort with status.
          return
        end if
      end do  ! n-loop
!
!-DEBUG_BEG:
#ifdef DEBUG_OCTFILL
  if ( status .ne. OCT_STATUS_OK ) then
    write(*,*) 'DEBUG_BEG[octnode_makechildren]'
    write(*,*) '  Should never get here!'
    write(*,*) '  nxyzdata   = ',octnode%nxyzdata
    write(*,*) '  k,n,status = ',k,n,status
    write(*,*) 'DEBUG_END[octnode_makechildren]'
  end if
#endif
!-DEBUG_END:
!
    end do ! k-loop
    ! Deallocate temp index array for a single point.
    deallocate( xyz_pts , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_pts",2,__LINE__,__FILE__,ierr,error_message)
    ! Deallocate existing octnode data.
    deallocate( octnode%ixyzdata , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"octnode%ixyzdatas",2,__LINE__,__FILE__,ierr, &
                     error_message)
    ! Remove the existing octnode data from the global count.
    octroot%totoctdata = octroot%totoctdata - octroot%MAXOCTDATA
    ! Clear the octnode data count.
    octnode%nxyzdata = 0
  end if
  !
  return
  end subroutine octnode_makechildren
!##########################################################################
  recursive &
  subroutine octnode_find_containing_octnode(octnode,octroot,octedat, &
               xyz_pt, containing_octnode,status)
  !
  ! Within this octnode, find the smallest octnode or sub-octnode whose bbox
  ! contains xyz_pt.  This node will have no children, and may or may not
  ! have data.
  !
  implicit none
  !
  !.. Formal Arguments ..
  type(octnode_type) ,target            ,intent(in)  :: octnode
  type(octroot_type)                    ,intent(in)  :: octroot
  class(octedat_type)                   ,intent(in)  :: octedat
  real(wp) ,dimension(:)                ,intent(in)  :: xyz_pt
  type(octnode_type) ,pointer           ,intent(out) :: containing_octnode
  integer                               ,intent(out) :: status
  !
  !.. Local Scalars ..
  integer l,n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octnode_find_containing_octnode"
  !
continue
  !
  ! Set default containing octnode
  containing_octnode => null()
  !
  ! If xyz_pt not within this octnode, return with error.
  status = OCT_STATUS_NOT_FOUND
  do l = 1,nr
    if ( xyz_pt(l) .lt. octnode%bbox(l,1) ) return
    if ( xyz_pt(l) .gt. octnode%bbox(l,2) ) return
  end do
!
  if ( associated(octnode%child) ) then
    ! This octnode has subnodes.  Search each one untill a match is found.
    do n = 1,octroot%MAXOCTNODES
      call octnode_find_containing_octnode(octnode%child(n), &
             octroot,octedat,xyz_pt,containing_octnode,status)
      if ( status .eq. OCT_STATUS_OK ) return
    end do
  else
    ! This is the smallest octnode containing xyz_pt.
    status = OCT_STATUS_OK
    containing_octnode => octnode
!
#ifdef OCTBIGGERCONTAINER
!-DEBUG_BEG:
!-DEBUG_CMT:  Increasing the container may be overly cautious 
!-DEBUG_CMT:  and increase cpu time.
!   If this octnode contains no data, set containing_octnode to parent.
    if ( .not. allocated(octnode%ixyzdata) ) &
      containing_octnode => octnode%parent
!-DEBUG_END:
#endif
!
  end if
  !
  return
  end subroutine octnode_find_containing_octnode
!##########################################################################
  recursive &
  subroutine octnode_find_nearest_xyz(octnode,octroot,octedat, &
               xyz_pt,xyz_bbox, &
               best_octnode,best_octixyz,best_ixyz,best_dist2,status)
  !
  ! Find point in octnode or subnodes that is nearest to desired xyz_pt.
  ! Only searches octnodes that intersect xyz_bbox that surrounds xyz_pt.
  !
  ! status has these possible return values:
  !   OCT_STATUS_OK          Closer match to xyz_pt was found.
  !   OCT_STATUS_NOT_FOUND   Closer match to xyz_pt was not found.
  !   OCT_STATUS_DATA_ERROR  Error evaluating data within the octtree.
  !
  implicit none
  !
  !.. Formal Arguments ..
  type(octnode_type)   ,target          ,intent(in)    :: octnode
  type(octroot_type)                    ,intent(in)    :: octroot
  class(octedat_type)                   ,intent(in)    :: octedat
  real(wp) ,dimension(:)                ,intent(in)    :: xyz_pt
  real(wp) ,dimension(:,:)              ,intent(in)    :: xyz_bbox
  type(octnode_type)   ,pointer         ,intent(inout) :: best_octnode
  integer                               ,intent(inout) :: best_octixyz
  integer              ,dimension(:)    ,intent(inout) :: best_ixyz
  real(wp)                              ,intent(inout) :: best_dist2
  integer                               ,intent(out)   :: status
  !
  !.. Local Scalars ..
  integer  :: ierr,l,k,n,octnode_status
  real(wp) :: dist_to_bbox,dist2
  !
  !.. Local Arrays ..
  real(wp) ,dimension(:,:)  ,allocatable :: xyz_pts
  real(wp) ,dimension(1:nr) :: dxyz
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octnode_find_nearest_xyz"
  !
continue
  !
  ! Default status, unless a better match is found.
  status = OCT_STATUS_NOT_FOUND
  !
  ! If any min/max edge of xyz_bbox is right/left of max/min of octnode%bbox,
  ! then boxes to not overlap.
  do l = 1,nr
    if ( xyz_bbox(l,1) .gt. octnode%bbox(l,2) ) return
    if ( xyz_bbox(l,2) .lt. octnode%bbox(l,1) ) return
  end do
  !
  !-COMMENT_BEG:
  !-COMMENT_CMT: This check may or may not improve search time.
  !-COMMENT_CMT: It is not necessary for the algorithm to work.
  !
  ! If the distance^2 from xyz_pt to the nearest edge of the bbox exceeds
  ! best_dist2, then any point in this octnode will exceed best_dist2.
  dist_to_bbox = huge(dist_to_bbox)
  do l = 1,nr
    dist_to_bbox =  min( dist_to_bbox, abs( xyz_pt(l) - octnode%bbox(l,1) ) &
                                     , abs( xyz_pt(l) - octnode%bbox(l,2) ) )
  end do
  if ( dist_to_bbox**2 .gt. best_dist2 ) return
  !-COMMENT_END:
  !
  if ( associated(octnode%child) ) then
    ! This octnode has subnodes.  Test each one for a closer match to xyz_pt.
    ! If any subnode is better, indicate with status.
    do n = 1,octroot%MAXOCTNODES
      call octnode_find_nearest_xyz(octnode%child(n),octroot,octedat, &
             xyz_pt,xyz_bbox, &
             best_octnode,best_octixyz,best_ixyz,best_dist2,octnode_status)
      if ( octnode_status .eq. OCT_STATUS_OK ) status = OCT_STATUS_OK
    end do

  else if ( allocated(octnode%ixyzdata) ) then
    ! This octnode has ixyzdata.  Test each point for a closer match to xyz_pt.
    ! If any is better, indicate with status and update best* variables.
    !
    ! Get xyz for all points in ixyzdata.
    allocate( xyz_pts(1:nr,1:octnode%nxyzdata) , stat=ierr , &
              errmsg=error_message )
    call alloc_error(pname,"xyz_pts",1,__LINE__,__FILE__,ierr,error_message)
    if ( ierr .ne. 0 ) then
      status = OCT_STATUS_MEMORY_ERROR
    else
      call octedat%evalxyz_pts(octnode%nxyzdata,octnode%ixyzdata, &
                               xyz_pts,status)
      if ( status .eq. OCT_STATUS_OK ) then
        do k = 1,octnode%nxyzdata
          dxyz = xyz_pt(:) - xyz_pts(:,k)
          dist2 = dot_product(dxyz,dxyz)
          if ( dist2 .lt. best_dist2 ) then
            best_octnode => octnode
            best_octixyz = k                     ! index within octnode data
            best_ixyz    = octnode%ixyzdata(:,k) ! index within octedat data
            best_dist2   = dist2                 ! distance squared
            status = OCT_STATUS_OK
          end if
        end do
        deallocate( xyz_pts , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"xyz_pts",2,__LINE__,__FILE__,ierr, &
                         error_message)
      end if
    end if
  end if
  !
  return
  end subroutine octnode_find_nearest_xyz
!##########################################################################
  recursive &
  subroutine octnode_print(octnode,octedat,status,filename)
  !
  ! Print the octnode to file for debugging.
  !
  implicit none
  !
  !.. Formal Arguments ..
  type(octnode_type)                    ,intent(in)    :: octnode
  class(octedat_type)                   ,intent(in)    :: octedat
  integer                               ,intent(out)   :: status
  character (len=*)           ,optional ,intent(in)    :: filename
  !
  !.. Local Scalars ..
  integer :: ierr,l,k,n
  !
  !.. Local Arrays ..
  ! Note:  xyz_pts(1:nr,1:octnode%nxyzdata)
  real(wp) ,allocatable ,dimension(:,:) :: xyz_pts
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "octnode_print"
  !
continue
  !
  status = OCT_STATUS_OK
  !
  if ( present(filename) ) then
    open(unit=20,file=trim(filename),form='formatted',status='unknown')
  end if
  !
#ifdef OCTID
  ! Print octnode id.
  write(20,*) 'octid    = ',octnode%octid
#endif
  !
  ! Print parent id.
  ! write(20,*) 'parent = ', ???   ! how to handle root node?
  !
  ! Print bounding box.
  if ( allocated(octnode%bbox) ) then
    write(20,*) 'bboxmin  = ',1,(octnode%bbox(l,1), l=1,nr)
    write(20,*) 'bboxmax  = ',2,(octnode%bbox(l,2), l=1,nr)
  else
    write(20,*) 'bboxmin  = undefined'
    write(20,*) 'bboxmax  = undefined'
  end if
  !
  ! Print data.
  write(20,*) 'nxyzdata = ',octnode%nxyzdata
  if ( allocated(octnode%ixyzdata) .and. octnode%nxyzdata .gt. 0 ) then
    allocate( xyz_pts(nr,octnode%nxyzdata) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_pts",1,__LINE__,__FILE__,ierr,error_message)
    call octedat%evalxyz_pts(octnode%nxyzdata,octnode%ixyzdata, &
                             xyz_pts,status)
    do k = 1,octnode%nxyzdata
      write(20,*) '           ', &
        (octnode%ixyzdata(l,k), l=1,size(octnode%ixyzdata,dim=1)), &
        ( xyz_pts(l,k), l=1,nr)
    end do
    deallocate( xyz_pts , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_pts",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Print children.
  if ( associated(octnode%child) ) then
    write(20,*) 'Contains children:'
#ifdef OCTID
    do n = 1,size(octnode%child)
      write(20,*) '  ',octnode%child(n)%octid
    end do
#endif
    write(20,*)
    do n = 1,size(octnode%child)
      call octnode_print(octnode%child(n),octedat,status)
    end do
  else
    write(20,*) 'Contains no children.'
    write(20,*)
  end if
  !
  if ( present(filename) ) then
    close(unit=20)
  end if
  return
  end subroutine octnode_print
!##########################################################################
!
end module octnode_mod
