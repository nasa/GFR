module octroot_mod
  !
  ! Module for octree root info.  Stuff common to the entire octtree.
  !
  !.. Use Statements ..
  use octstat_mod
  !
  implicit none
  !
  private
  !
  public :: octroot_init
  !
  ! Octtree root type.  Stuff common to the entire octtree.
  ! Set some default values.  Each octtree can then alter its own instance.
  !
  type, public :: octroot_type
    integer :: totoctdata  =  0  ! global number of data points in tree
    integer :: totoctnode  =  0  ! global number of octnodes in tree
    integer :: MAXOCTDATA  = 16  ! max data pts in each octnode
    integer :: MAXOCTNODES =  8  ! max subnodes in each octnode
                                 ! Typically this is 2**nr, but can't set
                                 ! here because nr is not a parameter.
                                 ! Set during init instead, if desired.
  contains
    procedure :: init    => octroot_init
    procedure :: dealloc => octroot_dealloc
    final :: octroot_final, octroot_rank1_final
  end type octroot_type
  !
contains
!
!###############################################################################
!
subroutine octroot_final(this)
  !
  ! Finalize the octroot_type.
  !
  !.. Passed Argument from Invoking Object ..
  type(octroot_type) :: this
  !
  !.. Local Scalars ..
  integer :: istat
  !
continue
  !
  call octroot_dealloc(this,istat)
  !
end subroutine octroot_final
!
!###############################################################################
!
subroutine octroot_rank1_final(this)
  !
  ! Finalize an array of octroot_type.
  !
  !.. Passed Argument from Invoking Object ..
  type(octroot_type) :: this(:)
  !
  !.. Local Scalars ..
  integer :: n,istat
  !
continue
  !
  do n = 1,size(this)
    call octroot_dealloc(this(n),istat)
  end do
  !
end subroutine octroot_rank1_final
!
!###############################################################################
!
subroutine octroot_dealloc(octroot,istat)
  !
  ! Deallocare the octroot data.
  !
  !.. Passed Argument from Invoking Object ..
  class(octroot_type), intent(inout) :: octroot
  !
  !.. Formal Arguments ..
  integer, intent(out) :: istat
  !
continue
  !
  ! Nothing much to dealloc, just reset counters.
  octroot%totoctdata = 0
  octroot%totoctnode = 0
  !
  istat = OCT_STATUS_OK
  !
end subroutine octroot_dealloc
!
!###############################################################################
!
subroutine octroot_init(octroot,maxoctdata,maxoctnodes,istat)
  !
  ! Init octroot values.
  !
  !.. Passed Argument from Invoking Object ..
  class(octroot_type), intent(inout) :: octroot
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: maxoctdata
  integer,  intent(in) :: maxoctnodes
  integer, intent(out) :: istat
  !
continue
  !
  ! Negative input values will retain current (likely default) setting.
  if (maxoctdata  >= 0) octroot%MAXOCTDATA  = maxoctdata
  if (maxoctnodes >= 0) octroot%MAXOCTNODES = maxoctnodes
  !
  ! Reset counters.
  octroot%totoctdata = 0
  octroot%totoctnode = 0
  !
  istat = OCT_STATUS_OK
  !
end subroutine octroot_init
!
!###############################################################################
!
end module octroot_mod
