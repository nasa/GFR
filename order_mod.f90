module order_mod
  !
  !.. Use Statements ..
  use module_kind_types, only : true,fals,lk,intseq
  use module_kind_types, only : Geom_Node,Geom_Edge,Geom_Tria,Geom_Quad
  use module_kind_types, only : Geom_Tetr,Geom_Pyra,Geom_Pris,Geom_Hexa
  use module_kind_types, only : Geom_Min,Geom_Max,Order_Min,Order_Max
  use module_kind_types, only : debug_timer,entering_procedure,leaving_procedure
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
  public :: find_valid_geometries
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
  integer, public, save :: t_order = -1
  !
  integer, public, save :: n_min_geom
  integer, public, save :: n_max_geom
  integer, public, save :: n_min_order
  integer, public, save :: n_max_order
  !
  integer, public, save :: maxpts
  integer, public, save :: maxSP ! maximum number of cell solution points
  integer, public, save :: maxQP ! maximum number of cell quadrature points
  integer, public, save :: maxFP ! maximum number of flux points
  integer, public, save :: maxGP ! maximum number of flux quadrature points
  integer, public, save :: maxOP ! maximum number of cell output points
  integer, public, save :: maxEP ! maximum number of cell error points
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
  !$ACC DECLARE CREATE(geom_solpts)
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
  ! GEOM_IS_USED : logical array that identifies whether or not a given
  !                element type could be used in the current simulation
  !
  logical(lk), public, protected, save :: geom_is_used(Geom_Min:Geom_Max)
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
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_geom_solpts"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
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
  n_min_order = min( n_order , p_order , q_order , t_order )
  n_max_order = max( n_order , p_order , q_order , t_order )
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
 !!$ACC ENTER DATA COPYIN(geom_solpts)
  !$ACC UPDATE DEVICE(geom_solpts)
  !
  call debug_timer(leaving_procedure,pname)
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
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "find_geom_max_interp_order"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
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
  call debug_timer(leaving_procedure,pname)
  !
end subroutine find_geom_max_interp_order
!
!###############################################################################
!
subroutine find_valid_geometries()
  !
  !.. Use Statements ..
  use geovar, only : cell_geom,cell_order
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
include "Functions/cell_solpts.f90"
!
!###############################################################################
!
end module order_mod
