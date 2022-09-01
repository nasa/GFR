module element_type_mod
  !
  use module_kind_types, only : wp,Geom_Edge
  !
  implicit none
  !
  private
  !
  public :: assign_solpts_edge_ptr
  !
  ! element_t    : derived type that stores the location of the solution points
  !                and their respective quadrature weights within the standard
  !                region of a given cell geometry.
  !
  type, public :: element_t
    real(wp), allocatable :: pts(:,:)
    real(wp), allocatable :: wts(:)
    real(wp), allocatable :: face_pts(:,:)
    real(wp), allocatable :: face_wts(:)
  end type element_t
  !
  ! std_elem : Array giving the location of the solution points and their
  !            respective quadrature weights for all combinations of cell
  !            geometry and cell order that are possible in the current
  !            simulation
  !
  !  For a given elemental geometry, this_geom,
  !                   and its order, this_order :
  !
  !  std_elem(this_geom,this_order)%pts(l,k) : The location of the solution
  !                                            point 'k' in the direction 'l' of
  !                                            the local coordinate system
  !                                            within the standard region of the
  !                                            current element geometry.
  !
  !  std_elem(this_geom,this_order)%wts(k) :
  !                       The quadrature weight of the solution point 'k' within
  !                       the standard region of the current element geometry.
  !
  type(element_t), public, save, target, allocatable :: std_elem(:,:)
  !
  ! solpts_edge : pointer to the 'n_order' 1D edge solution points
  !
  real(wp), public, protected, pointer :: solpts_edge(:) => null()
  !
  ! face_elem : Array giving the location of the flux points and their
  !             respective quadrature weights for all combinations of cell
  !             geometry and cell order that are possible in the current
  !             simulation
  !
  !  For a given face geometry, this_geom,
  !              and its order, this_order :
  !
  !  face_elem(this_geom,this_order)%pts(l,k) :
  !                       The location of the flux point 'k' in the direction
  !                       'l' of the local coordinate system within the
  !                       standard region of the current face geometry.
  !
  !  face_elem(this_geom,this_order)%wts(k) :
  !                       The quadrature weight of the flux point 'k' within
  !                       the standard region of the current face geometry.
  !
  type(element_t), public, save, target, allocatable :: face_elem(:,:)
  !
contains
!
!###############################################################################
!
subroutine assign_solpts_edge_ptr(this_order)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: this_order
  !
continue
  !
  if (allocated(std_elem(Geom_Edge,this_order)%pts)) then
    solpts_edge => std_elem(Geom_Edge,this_order)%pts(1,:)
  else
    solpts_edge => null()
  end if
  !
end subroutine assign_solpts_edge_ptr
!
!###############################################################################
!
end module element_type_mod
