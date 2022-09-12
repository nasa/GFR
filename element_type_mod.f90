module element_type_mod
  !
  use module_kind_types, only : wp,true,fals,lk,intseq
  use module_kind_types, only : Geom_Min,Geom_Max
  use module_kind_types, only : Geom_Edge,Geom_Tria,Geom_Quad
  use module_kind_types, only : Geom_Tetr,Geom_Pyra,Geom_Pris,Geom_Hexa
  use module_kind_types, only : entering_procedure,leaving_procedure,debug_timer
  !
  implicit none
  !
  private
  !
  public :: find_valid_geometries
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
  ! GEOM_IS_USED : logical array that identifies whether or not a given
  !                element type could be used in the current simulation
  !
  logical(lk), public, protected, save :: geom_is_used(Geom_Min:Geom_Max)
  !
contains
!
!###############################################################################
!
subroutine find_valid_geometries()
  !
  !.. Use Statements ..
  use geovar,    only : cell_geom,cell_order
  use order_mod, only : n_min_geom,n_max_geom,q_order
  use order_mod, only : geom_solpts,geom_flxpts
  use order_mod, only : maxpts,maxSP,maxFP,maxQP,maxGP
  use order_mod, only : o_order,maxOP
  use order_mod, only : e_order,maxEP
  !
  !.. Local Scalars ..
  integer :: nc,ne,no,this_geom
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "find_valid_geometries"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Initialize geom_is_used to false
  !
  geom_is_used = fals
  !
  ! Loop through the geometries and mark which ones are explicitly used
  !
  do this_geom = Geom_Min,Geom_Max
    geom_is_used(this_geom) = any(cell_geom == this_geom)
  end do
  !
  ! Now go back and mark any lower dimension geometries needed by
  ! those geometries already marked (e.g., Geom_Edge needed by Geom_Quad)
  !
  ! Mark Geom_Quad as used if any 3D elements with quad faces are used
  !
  if (any( geom_is_used([Geom_Pyra,Geom_Pris,Geom_Hexa]) )) then
    geom_is_used(Geom_Quad) = true
  end if
  !
  ! Mark Geom_Tria as used if any 3D elements with triangle faces are used
  !
  if (any( geom_is_used([Geom_Tetr,Geom_Pyra,Geom_Pris]) )) then
    geom_is_used(Geom_Tria) = true
  end if
  !
  ! Mark Geom_Edge as used if any 2D elements are used
  ! NOTE: This should always be true for 2D or 3D
  !
  if (any( geom_is_used([Geom_Tria,Geom_Quad]) )) then
    geom_is_used(Geom_Edge) = true
  end if
  !
  ! Find the minimum and maximum geometry types
  !
  n_min_geom = minval( intseq(Geom_Min,Geom_Max) , mask=geom_is_used )
  n_max_geom = maxval( intseq(Geom_Min,Geom_Max) , mask=geom_is_used )
  !
  ! While we are here, find the maximum number of points in any cell
  !
  maxSP = 1
  maxFP = 1
  maxQP = 1
  maxGP = 1
  maxOP = 1
  maxEP = 1
  !
  do nc = 1,size(cell_geom)
    !
    maxSP = max( maxSP , geom_solpts(cell_geom(nc),cell_order(nc)) )
    maxQP = max( maxQP , geom_solpts(cell_geom(nc),q_order) )
    !
    maxFP = max( maxFP , geom_flxpts(cell_geom(nc),cell_order(nc)) )
    maxGP = max( maxGP , geom_flxpts(cell_geom(nc),q_order) )
    !
  end do
  !
  if (o_order > 0) then
    do nc = 1,size(cell_geom)
      no = max(cell_order(nc),o_order)
      maxOP = max( maxOP , geom_solpts(cell_geom(nc),no) )
    end do
  else
    maxOP = maxSP
  end if
  !
  if (e_order > 0) then
    do nc = 1,size(cell_geom)
      ne = max(cell_order(nc),e_order)
      maxEP = max( maxEP , geom_solpts(cell_geom(nc),ne) )
    end do
  else
    maxEP = maxSP
  end if
  !
  maxpts = max( maxSP , maxQP , maxFP , maxGP , maxOP , maxEP )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine find_valid_geometries
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
