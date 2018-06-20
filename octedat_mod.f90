module octedat_mod
!
! Base module for defining an abstract base type for storing
! encapsulated data for use with the octtree search.
!
  use module_kind_types, only : wp
  use octstat_mod                   ! will be needed by future derived types

  implicit none

  ! Define an abstract base type for encapsulated data which will be used to
  ! hide details of the actual xyz array/structure from the octtree routines.
  type ,abstract :: octedat_type
    ! Define two strings to identify what xyz data is contained.
    ! This is mainly for convenience/debugging.
    !   varname = name of variable being stored
    !   vardesc = description of variable (ie, a sentence)
    character(len=32)  :: varname = ""
    character(len=256) :: vardesc = ""
    !
    ! Pointer to xyz array/structure to be provided by the derived type.
    ! Cannot define here.  Don't know what the type/rank will be yet.
    ! Might require several different types (cells,faces,etc).
!   real(wp) ,dimension(:,:,...) ,pointer :: xyzdata    ! array
!   type(face_xyz)               ,pointer :: xyzdata    ! structure
    !
    ! Specify the maximun number of indices that will be used to
    ! access a given xyzdata element.  This value can be changed by
    ! the derived type.
    integer :: MAXINDICES = 1
    !
contains
    ! Functions to be provided by the derived type.  See below.
    procedure (octedat_dealloc_interface)    ,deferred :: dealloc
!   procedure (octedat_evalbbox_interface)   ,deferred :: evalbbox
    procedure (octedat_evalxyzpt_interface)  ,deferred :: evalxyz_pt
    procedure (octedat_evalxyzpts_interface) ,deferred :: evalxyz_pts
!   generic :: evalxyz => evalxyz_pt, evalxyz_pts
  end type octedat_type

  abstract interface
    ! Deallocate any memory inside the octedat structure.
    subroutine octedat_dealloc_interface(this,status)
      import :: octedat_type
      class(octedat_type)                 ,intent(inout) :: this
      integer                             ,intent(out)   :: status
    end subroutine octedat_dealloc_interface
  end interface
  
! abstract interface
!   ! Evaluate the bounding box, bbox(1:nr,1:2), with min/max values of
!   ! the octedat xyzdata array/struct.
!   ! Status indicates OCT_STATUS_OK or OCT_STATUS_DATA_ERROR.
!   subroutine octedat_evalbbox_interface(this,bbox,status)
!     import :: octedat_type
!     import :: wp
!     class(octedat_type)                 ,intent(in)    :: this
!     real(wp)            ,dimension(:,:) ,intent(out)   :: bbox
!     integer                             ,intent(out)   :: status
!   end subroutine octedat_evalbbox_interface
! end interface

  abstract interface
    ! Evaluate xyz(1:nr) at a single point given by idx(:).
    ! Elements of idx(:) serve as indices to octedat xyzdata array/struct.
    ! Status indicates OCT_STATUS_OK or OCT_STATUS_DATA_ERROR.
    subroutine octedat_evalxyzpt_interface(this,idx,xyz,status)
      import :: octedat_type
      import :: wp
      class(octedat_type)                 ,intent(in)    :: this
      integer             ,dimension(:)   ,intent(in)    :: idx
      real(wp)            ,dimension(:)   ,intent(out)   :: xyz
      integer                             ,intent(out)   :: status
    end subroutine octedat_evalxyzpt_interface
  end interface

  abstract interface
    ! Evaluate xyz(1:nr,1:npts) at a multiple point given by idx(:,1:npts).
    ! Elements of idx(:,ipt) serve as indices to octedat xyzdata array/struct.
    ! Status indicates OCT_STATUS_OK or OCT_STATUS_DATA_ERROR.
    subroutine octedat_evalxyzpts_interface(this,npts,idx,xyz,status)
      import :: octedat_type
      import :: wp
      class(octedat_type)                 ,intent(in)    :: this
      integer                             ,intent(in)    :: npts
      integer             ,dimension(:,:) ,intent(in)    :: idx
      real(wp)            ,dimension(:,:) ,intent(out)   :: xyz
      integer                             ,intent(out)   :: status
    end subroutine octedat_evalxyzpts_interface
  end interface

!##########################################################################

end module octedat_mod
