module interpsolvar_mod
  !
  ! Module to define variables used in interpolation of
  ! (grid1,sol1) -> (grid2,sol2).
  !
  !.. Use Statements ..
  use module_kind_types, only : wp,lk
  use geovar,            only : geom_conn, conflict, fp_t
  use geovar,            only : cell_t, face_t, grid_t
  use geovar,            only : map_t
  use geovar,            only : bc_grp_t
  use flowvar,           only : sol_nodes
 !use flowvar,           only : faceusp_array
 !use flowvar,           only : faceusp_dt
 !use flowvar,           only : face_xyz,face_var
  !
  use matrix_types_mod,  only : matrix_gpu_t
  use quadrature_mod,    only : element_type
  use vandermonde_mod,   only : vandermonde_matrix
  use derivatives_mod,   only : derivative_matrices
  use interpolation_mod, only : interpolation_matrix
  use interpolation_mod, only : outerpolation_matrix
  use projection_mod,    only : projection_matrix
  use filter_mod,        only : filter_t
  !
  implicit none
  !
  private
  !
  public :: sol_nodes
 !public :: faceusp_dt
 !public :: faceusp_array
  !
  ! Grid2 copies of variables from order_mod.
  integer, public, save :: g2_n_order
  integer, public, save :: g2_p_order
  integer, public, save :: g2_q_order
  integer, public, save :: g2_o_order
  integer, public, save :: g2_e_order
  integer, public, save :: g2_t_order
  !
  ! Grid2 copies of variables from geovar.
  integer,     public, save :: g2_n_solpts
  integer,     public, save :: g2_n_flxpts
  logical(lk), public, save :: g2_using_quadrature
  integer,     public, save :: g2_nr
  integer,     public, save :: g2_nnode
  integer,     public, save :: g2_ncell
  integer,     public, save :: g2_nface
  integer,     public, save :: g2_nedge
  integer,     public, save :: g2_nfbnd
  integer,     public, save :: g2_n_totcel
  integer,     public, save :: g2_n_totnod
  integer,     public, save :: g2_n_totpts
  !
  integer,  public, save, allocatable :: g2_cell_geom(:)
  integer,  public, save, allocatable :: g2_cell_order(:)
  real(wp), public, save, allocatable :: g2_xyz_nodes(:,:)
  real(wp), public, save, allocatable :: g2_xyz(:,:)
  !
  integer, public, save :: g2_n_global_node
  integer, public, save :: g2_n_global_totnod
  integer, public, save :: g2_n_global_cell
  integer, public, save :: g2_n_global_totcel
  integer, public, save :: g2_n_global_edge
  integer, public, save :: g2_n_global_fbnd
  integer, public, save :: g2_n_global_solpts
  integer, public, save :: g2_n_global_totpts
  !
  real(wp), public, save, allocatable :: g2_global_xyz_nodes(:,:)
  real(wp), public, save, allocatable :: g2_global_xyz(:,:)
  integer,  public, save, allocatable :: g2_global_nofc_ptr(:)
  integer,  public, save, allocatable :: g2_global_nofc(:)
  integer,  public, save, allocatable :: g2_global_pinc_ptr(:)
  integer,  public, save, allocatable :: g2_global_cell_geom(:)
  integer,  public, save, allocatable :: g2_global_cell_order(:)
  integer,  public, save, allocatable :: g2_bface(:,:)
  !
  type(conflict), public, save, allocatable :: g2_bc_conflict(:)
  !
  integer, public, save, allocatable :: g2_nodes_of_cell_ptr(:)
  integer, public, save, allocatable :: g2_nodes_of_cell(:)
  !
  type(geom_conn), public, save, allocatable :: g2_pst_geom(:,:)
  !
  type(fp_t), public, save, allocatable :: g2_fp
  !
  type(face_t), public, save, allocatable :: g2_face(:)
  type(cell_t), public, save, allocatable :: g2_cell(:)
  !
  type(grid_t), public, save, allocatable :: g2_grid
  !
  type(bc_grp_t), public, save, allocatable :: g2_bc_grps(:)
  !
  ! Grid2 copies of variables from flowvar.
  real(wp), public, save, allocatable :: g2_sol_nodes(:,:)
  real(wp), public, save, allocatable :: g2_usp(:,:)
  real(wp), public, save, allocatable :: g2_gpsp(:,:)
  real(wp), public, save, allocatable :: g2_dusp(:,:,:)
  !
 !type(face_xyz), public, save, allocatable :: g2_facexyz_dt(:)
 !type(face_var), public, save, allocatable :: g2_faceusp_dt(:)
  !
  real(wp), public, save, allocatable :: g2_facexyz_array(:,:)
  real(wp), public, save, allocatable :: g2_faceusp_array(:,:,:)
  !
  ! Grid2 copies of variables from ovar.
  character(len=32), public, save, allocatable :: g2_bc_names(:)
  !
  ! Grid2 copies of variables from parallel.
  type(map_t), public, save, allocatable :: g2_cell_map(:)
  !
  ! Grid2 copies of variables from quadrature_mod
  type(element_type), public, save, allocatable :: g2_std_elem(:,:)
  type(element_type), public, save, allocatable :: g2_face_elem(:,:)
  !
  ! Grid2 copies of variables from vandermonde_mod
  type(vandermonde_matrix), public, save, allocatable :: g2_vand(:,:)
  !
  ! Grid2 copies of variables from derivatives_mod
  type(derivative_matrices), public, save, allocatable :: g2_deriv(:,:)
  !
  ! Grid2 copies of variables from correction_mod
  type(matrix_gpu_t), public, save, allocatable :: g2_correct(:,:)
  !
  ! Grid2 copies of variables from interpolation_mod
  type(interpolation_matrix), public, save, allocatable :: g2_interp(:,:)
  type(outerpolation_matrix), public, save, allocatable :: g2_outerp(:,:)
  !
  ! Grid2 copies of variables from projection_mod
  type(projection_matrix), public, save, allocatable :: g2_projct(:,:)
  !
  ! Grid2 copies of variables from filter_mod
  type(filter_t), public, save, allocatable :: g2_filter(:,:)
  !
  ! Grid2 copies of variables from metrics_mod
  real(wp), public, save, allocatable :: g2_metrics_array(:,:,:)
  real(wp), public, save, allocatable :: g2_jac(:)
  real(wp), public, save, allocatable :: g2_wtdjac(:)
  real(wp), public, save, allocatable :: g2_facenrm(:,:)
  real(wp), public, save, allocatable :: g2_facejac(:,:)
  !
end module interpsolvar_mod
!
