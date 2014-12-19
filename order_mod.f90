module order_mod
  !
  !.. Use Statements ..
  use module_kind_types, only : Geom_Node,Geom_Edge,Geom_Tria,Geom_Quad
  use module_kind_types, only : Geom_Tetr,Geom_Pyra,Geom_Pris,Geom_Hexa
  use module_kind_types, only : Geom_Min,Geom_Max,Order_Min,Order_Max
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
  integer, public, save :: all_orders = -3
  !
  integer, public, save :: n_min_geom
  integer, public, save :: n_max_geom
  integer, public, save :: n_min_order
  integer, public, save :: n_max_order
  !
  integer, public, save :: maxpts
  !
  ! Number of solution points for each cell geometry type
  !
  integer, public, protected, save :: n_sp1de
  integer, public, protected, save :: n_sp2dt
  integer, public, protected, save :: n_sp2dq
  integer, public, protected, save :: n_sptet
  integer, public, protected, save :: n_sppyr
  integer, public, protected, save :: n_sppri
  integer, public, protected, save :: n_sphex
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
  integer, public, protected, save :: n_qptet
  integer, public, protected, save :: n_qppyr
  integer, public, protected, save :: n_qppri
  integer, public, protected, save :: n_qphex
  !
  ! Number of solution points for all geometries of orders 0 to max_order
  !
  integer, public, protected, save :: geom_solpts(Geom_Min:Geom_Max,Order_Min:Order_Max)
  !
contains
!
!###############################################################################
!
subroutine initialize_geom_solpts(use_qdtr_face_pts)
  !
  !.. Formal Arguments ..
  logical, intent(in) :: use_qdtr_face_pts
  !
  !.. Local Scalars ..
  integer :: i,j
  !
continue
  !
  do j = lbound(geom_solpts,dim=2),ubound(geom_solpts,dim=2)
    do i = lbound(geom_solpts,dim=1),ubound(geom_solpts,dim=1)
      geom_solpts(i,j) = Cell_Solpts(i,j)
    end do
  end do
  !
  n_sp1de = edge_solpts(n_order)
  n_sp2dt = tria_solpts(n_order)
  n_sp2dq = quad_solpts(n_order)
  n_sptet = tetr_solpts(n_order)
  n_sppyr = pyra_solpts(n_order)
  n_sppri = pris_solpts(n_order)
  n_sphex = hexa_solpts(n_order)
  !
  n_qp1de = edge_solpts(q_order)
  n_qp2dt = tria_solpts(q_order)
  n_qp2dq = quad_solpts(q_order)
  n_qptet = tetr_solpts(q_order)
  n_qppyr = pyra_solpts(q_order)
  n_qppri = pris_solpts(q_order)
  n_qphex = hexa_solpts(q_order)
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
  if (o_order > 0) then
    n_min_order = min( n_min_order , o_order )
    n_max_order = max( n_max_order , o_order )
  end if
  !
end subroutine initialize_geom_solpts
!
!###############################################################################
!
include "Functions/cell_solpts.f90"
!
!###############################################################################
!
end module order_mod
