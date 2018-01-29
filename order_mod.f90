module order_mod
  !
  !.. Use Statements ..
  use module_kind_types, only : Geom_Node,Geom_Edge,Geom_Tria,Geom_Quad
  use module_kind_types, only : Geom_Tetr,Geom_Pyra,Geom_Pris,Geom_Hexa
  use module_kind_types, only : Geom_Min,Geom_Max,Order_Min,Order_Max,lk
  !
  implicit none
  !
  private
  !
  ! #########################
  ! ### Public Procedures ###
  ! #########################
  !
  public :: initialize_geom_solpts
  public :: find_geom_max_interp_order
  !
  ! #########################
  ! ### Module Parameters ###
  ! #########################
  !
  ! ######################
  ! ### Public Scalars ###
  ! ######################
  !
  integer, public, save :: n_order = 3
  integer, public, save :: p_order = 3
  integer, public, save :: q_order = 3
  integer, public, save :: o_order = -3
  integer, public, save :: e_order = -3
  !
  integer, public, save :: n_min_geom
  integer, public, save :: n_max_geom
  integer, public, save :: n_min_order
  integer, public, save :: n_max_order
  !
  integer, public, save :: maxpts
  integer, public, save :: maxSP
  integer, public, save :: maxQP
  integer, public, save :: maxFP
  integer, public, save :: maxGP
  integer, public, save :: maxOP
  integer, public, save :: maxEP
  !
  ! Number of solution points for each cell geometry type
  !
  integer, public, protected, save :: n_sp1de
  integer, public, protected, save :: n_sp2dt
  integer, public, protected, save :: n_sp2dq
  !
  ! Number of face points for each face geometry type
  !
  integer, public, protected, save :: n_fp1de
  integer, public, protected, save :: n_fp2dt
  integer, public, protected, save :: n_fp2dq
  !
  ! Number of quadrature / over-integration
  ! solution points for each cell geometry type
  !
  integer, public, protected, save :: n_qp1de
  integer, public, protected, save :: n_qp2dt
  integer, public, protected, save :: n_qp2dq
  !
  ! Number of solution points for all geometries of orders 0 to max_order
  !
  integer, public, protected, save :: &
                              geom_solpts(Geom_Min:Geom_Max,Order_Min:Order_Max)
  !
  ! Number of face/flux points for all geometries of orders 0 to max_order
  !
  integer, public, protected, save :: &
                              geom_flxpts(Geom_Min:Geom_Max,Order_Min:Order_Max)
  !
  ! Number of edge/flux points for all geometries of orders 0 to max_order
  !
  integer, public, protected, save :: &
                              geom_edgpts(Geom_Min:Geom_Max,Order_Min:Order_Max)
  !
  ! geom_mio "geom_max_interp_order" gives the maximum interpolation order
  ! required for a given combination of cell geometry and cell order
  !
  integer, public, protected, save :: &
                              geom_mio(Geom_Min:Geom_Max,Order_Min:Order_Max)
  !
contains
!
!###############################################################################
!
subroutine initialize_geom_solpts(use_qdtr_face_pts)
  !
  !.. Formal Arguments ..
  logical(lk), intent(in) :: use_qdtr_face_pts
  !
  !.. Local Scalars ..
  integer :: i,j
  !
continue
  !
  do j = lbound(geom_solpts,dim=2),ubound(geom_solpts,dim=2)
    do i = lbound(geom_solpts,dim=1),ubound(geom_solpts,dim=1)
      geom_solpts(i,j) = Cell_Solpts(i,j)
      geom_flxpts(i,j) = Cell_Flxpts(i,j)
      geom_edgpts(i,j) = Cell_Edgpts(i,j)
    end do
  end do
  !
  n_sp1de = edge_solpts(n_order)
  n_sp2dt = tria_solpts(n_order)
  n_sp2dq = quad_solpts(n_order)
  !
  n_qp1de = edge_solpts(q_order)
  n_qp2dt = tria_solpts(q_order)
  n_qp2dq = quad_solpts(q_order)
  !
  if (use_qdtr_face_pts) then
    n_fp1de = n_qp1de
    n_fp2dt = n_qp2dt
    n_fp2dq = n_qp2dq
  else
    n_fp1de = n_sp1de
    n_fp2dt = n_sp2dt
    n_fp2dq = n_sp2dq
  end if
  !
  n_min_order = min( n_order , p_order , q_order )
  n_max_order = max( n_order , p_order , q_order )
  !
 !if (o_order > 0) then
 !  n_min_order = min( n_min_order , o_order )
 !  n_max_order = max( n_max_order , o_order )
 !end if
 !!
 !if (e_order > 0) then
 !  n_min_order = min( n_min_order , e_order )
 !  n_max_order = max( n_max_order , e_order )
 !end if
  !
end subroutine initialize_geom_solpts
!
!###############################################################################
!
subroutine find_geom_max_interp_order(cell_geom,cell_order,cells_on_left, &
                                      cells_on_right,order_of_faces)
  !
  !.. Formal Arguments ..
  integer, dimension(:), intent(in) :: cell_geom
  integer, dimension(:), intent(in) :: cell_order
  integer, dimension(:), intent(in) :: cells_on_left
  integer, dimension(:), intent(in) :: cells_on_right
  integer, dimension(:), intent(in) :: order_of_faces
  !
  !.. Local Scalars ..
  integer :: nf,face_order
  integer :: lgeom,lorder,left_cell
  integer :: rgeom,rorder,right_cell
  !
continue
  !
  geom_mio = -1
  !
  do nf = 1,size(order_of_faces)
    !
    left_cell = cells_on_left(nf)
    right_cell = cells_on_right(nf)
    !
    lgeom = cell_geom(left_cell)
    rgeom = cell_geom(right_cell)
    !
    lorder = cell_order(left_cell)
    rorder = cell_order(right_cell)
    !
    face_order = order_of_faces(nf)
    !
    geom_mio(lgeom,lorder) = max( geom_mio(lgeom,lorder) , face_order )
    geom_mio(rgeom,rorder) = max( geom_mio(rgeom,rorder) , face_order )
    !
  end do
  !
end subroutine find_geom_max_interp_order
!
!###############################################################################
!
include "Functions/cell_solpts.f90"
!
!###############################################################################
!
end module order_mod
