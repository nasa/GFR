module interpsolvar_mod

! Module to define variables used in interpolation of
! (grid1,sol1) -> (grid2,sol2).

  use module_kind_types, only : wp,lk
  use geovar,            only : geom_conn, conflict, fp_t
  use geovar,            only : cell_t, face_t, grid_t
  use geovar,            only : map_t
  use flowvar,           only : face_xyz,face_var,faceusp_dt,sol_nodes
  implicit none
  private
  public :: faceusp_dt,sol_nodes

! Grid2 copies of variables from geovar.
  integer, public, save :: g2_n_solpts
  integer, public, save :: g2_n_flxpts
  logical(lk), public, save :: g2_using_quadrature
  integer, public, save :: g2_nr
  integer, public, save :: g2_nnode
  integer, public, save :: g2_ncell
  integer, public, save :: g2_nface
  integer, public, save :: g2_nedge
  integer, public, save :: g2_nfbnd
  integer, public, save :: g2_n_totcel
  integer, public, save :: g2_n_totnod
  integer, public, save :: g2_n_totpts
  integer, public, save, allocatable :: g2_cell_geom(:)
  integer, public, save, allocatable :: g2_cell_order(:)
  real(wp), public, save, allocatable :: g2_xyz_nodes(:,:)
  real(wp), public, save, allocatable :: g2_xyz(:,:)
  integer, public, save :: g2_n_global_node
  integer, public, save :: g2_n_global_totnod
  integer, public, save :: g2_n_global_cell
  integer, public, save :: g2_n_global_totcel
  integer, public, save :: g2_n_global_edge
  integer, public, save :: g2_n_global_fbnd
  integer, public, save :: g2_n_global_solpts
  integer, public, save :: g2_n_global_totpts
  real(wp), public, save, allocatable :: g2_global_xyz_nodes(:,:)
  real(wp), public, save, allocatable :: g2_global_xyz(:,:)
  integer,  public, save, allocatable :: g2_global_nofc_ptr(:)
  integer,  public, save, allocatable :: g2_global_nofc(:)
  integer,  public, save, allocatable :: g2_global_pinc_ptr(:)
  integer,  public, save, allocatable :: g2_global_cell_geom(:)
  integer,  public, save, allocatable :: g2_global_cell_order(:)
  integer, public, save, allocatable :: g2_bface(:,:)
  type(conflict), public, save, allocatable :: g2_bc_conflict(:)
  integer, public, save, allocatable :: g2_nodes_of_cell_ptr(:)
  integer, public, save, allocatable :: g2_nodes_of_cell(:)
  type(geom_conn), public, save, allocatable :: g2_pst_geom(:,:)
  type(fp_t), public, save, allocatable :: g2_fp
  type(face_t), public, save, allocatable :: g2_face(:)
  type(cell_t), public, save, allocatable :: g2_cell(:)
  type(grid_t), public, save, allocatable :: g2_grid

! Grid2 copies of variables from flowvar.
  real(wp), public, save, allocatable :: g2_sol_nodes(:,:)
  real(wp), public, save, allocatable :: g2_usp(:,:)
  real(wp), public, save, allocatable :: g2_dusp(:,:,:)
  type(face_xyz), public, save, allocatable :: g2_facexyz_dt(:)
  type(face_var), public, save, allocatable :: g2_faceusp_dt(:)

  real(wp), public, save, allocatable :: g2_facexyz_array(:,:)
  real(wp), public, save, allocatable :: g2_faceusp_array(:,:,:)

! Grid2 copies of variables from ovar.
  character(len=32), public, save, allocatable :: g2_bc_names(:)

! Grid2 copies of variables from parallel.
  type(map_t), public, save, allocatable :: g2_cell_map(:)

!
!##########################################################################
end module interpsolvar_mod
