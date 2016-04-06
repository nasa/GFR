module postprocessor_mod
  !
  !.. Global Use Statements ..
  use module_kind_types
  use eqn_idx, only : nq
  !
  implicit none
  !
  private
  !
  !
  ! ##################################
  ! ###  Module Public Procedures  ###
  ! ##################################
  !
  public :: get_output_connectivity
  !
  !
  ! ################################
  ! ###  Module Local Variables  ###
  ! ################################
  !
  integer, save :: lcell
  integer, save :: lnode
  integer, save :: lfbnd
  integer, save :: lface
  integer, save :: ledge
  !
  type :: pnt_map_t
    integer, allocatable :: sp(:)
    integer, allocatable :: fp(:)
    integer, allocatable :: ep(:)
    integer, allocatable :: np(:)
  end type pnt_map_t
  !
  type :: output_face_pts_t
    integer              :: face_idx
    integer, allocatable :: out_idx(:)
  end type output_face_pts_t
  !
  type :: output_edge_pts_t
    integer              :: edge_idx
    integer, allocatable :: out_idx(:)
  end type output_edge_pts_t
  !
  type :: output_node_pts_t
    integer :: node_idx
    integer :: out_idx
  end type output_node_pts_t
  !
  integer,                 allocatable :: output_sol_pts(:)
  !
  type(output_face_pts_t), allocatable :: output_face_pts(:)
  !
  type(output_edge_pts_t), allocatable :: output_edge_pts(:)
  !
  type(output_node_pts_t), allocatable :: output_node_pts(:)
  !
  type :: plot_func_t
    integer :: func
  end type plot_func_t
  !
  type(plot_func_t), public, parameter :: &
              plot_static_density     = plot_func_t(1), &
              plot_velocity_x         = plot_func_t(2), &
              plot_velocity_y         = plot_func_t(3), &
              plot_velocity_z         = plot_func_t(4), &
              plot_static_pressure    = plot_func_t(5), &
              plot_static_temperature = plot_func_t(6), &
              plot_momentum_x         = plot_func_t(7), &
              plot_momentum_y         = plot_func_t(8), &
              plot_momentum_z         = plot_func_t(9), &
              plot_entropy            = plot_func_t(10), &
              plot_mach_number        = plot_func_t(11), &
              plot_vorticity_x        = plot_func_t(12), &
              plot_vorticity_y        = plot_func_t(13), &
              plot_vorticity_z        = plot_func_t(14)
  !
  public :: get_plot_functions
  !
  !
  FINISHED SUBROUTINE MAP_TO_OUTPUT_POINTS AND THE RELATED FUNCTIONS FOR
  HEXAHEDRA CELLS. NOW NEED TO GO BACK AND FIGURE OUT HOW TO APPLY THIS
  FOR THE CONNECTIVITY AND FOR CREATING THE SOLUTION ARRAY SENT TO CGNS.
  !
  WHEN CREATING THE SOLUTION ARRAY, NEED TO FIGURE OUT WHAT TO DO ABOUT
  POINTS FOR WHICH AN ADJACENT PARTITION IS RESPONSIBLE. IN THESE CASES
  THERE IS A POINT IN THE SOLUTION ARRAY THAT NEEDS TO BE SKIPPED FOR A
  PARTITION THAT SHARES THIS POINT BUT IS NOT RESPONSIBLE FOR ITS SOLUTION.
  !
  !
  " NEED FIGURE OUT CONNECTIVITY OF THE GLOBAL GRID AND WHICH CELLS ARE"
  " RESPONSIBLE FOR WHICH OUTPUT POINTS."
  "  - NEED TO FIGURE THIS OUT FOR BOTH SERIAL AND PARALLEL CASES."
  "  - USE CELL(:)%NODE(:)%CELL_IDX, CELL(:)%EDGE(:)%CELL_IDX(:), AND"
  "    CELL(:)%FACE(:)%CELL_IDX(:) TO MAP THE LOCATIONS OF THE NODE, EDGE,"
  "    AND FACE VALUES WITHIN NODEUSP, NODEDUSP, EDGEUSP, EDGEDUSP, FACEUSP,"
  "    AND FACEDUSP TO THE CORRECT LOCATIONS WITHIN THE ARRAY THATS CREATED"
  "    TO COLLECT THE SOLUTIONS FOR SENDING TO CGNS."
  "  - TO REMOVE THE DUPLICATE POINT INDICES ON PARTITION BOUNDARIES IN THE"
  "    CGNS FILE, A POSSIBLE PATH COULD BE TO USE THE NODEUSP, EDGEUSP, AND"
  "    FACEUSP ARRAYS TO STORE A REAL VALUE OF A GIVEN POINT INDEX ON THE"
  "    PROCESSOR THAT IS RESPONSIBLE FOR WRITING THE SOLUTION AT THIS POINT,"
  "    AND ALL OTHER PROCESSORS STORE A NEGATIVE VALUE (OR SOMETHING THAT CAN"
  "    BE IDENTIFIED AS A POINT INDEX FROM THE WRONG PROCESSOR) AND USE THE"
  "    EXISTING COMMUNICATION/EXCHANGE CODE TO SYNCRONIZE THE CORRECT POINT"
  "    INDEX ON ALL PROCESSORS. THERE COULD BE SOME LOGIC THAT ONLY ADDS"
  "    POSITIVE VALUES FROM THE RECEIVE BUFFER BACK INTO NODEUSP, EDGEUSP,"
  "    OR FACEUSP. AFTER THAT, ITS JUST A MATTER OF CONVERTING THE REAL"
  "    VALUES TO INTEGERS AND USE THEM TO OVERWRITE THE PROPER POINT INDICES"
  "    IN THE CONNECTIVITY ARRAY THAT SENT INTO CGNS. JUST MAKE SURE NOT"
  "    TO INCLUDE ANY AVERAGING COEFFICIENTS."
  !
  " NEED TO FIGURE OUT HOW TO OVERRIDE THE INTERIOR SOLUTIONS FOR THE"
  " SOLUTIONS ON BOUNDARY FACES, EDGES, AND NODES."
  "   - NEED TO ADJUST THE LOCAL CELLS_WITH_NODE AND CELLS_WITH_EDGE ARRAYS"
  "     ON EACH PARTITION SO THAT IF A NODE OR AN EDGE ARE ON A BOUNDARY"
  "     FACES THE ONLY CELLS REMAINING WITHIN THESE ARRAYS ARE THE GHOST"
  "     CELLS/FACES THAT CONTRIBUTE TO THE SOLUTION AT THE NODE/EDGE."
  "       - NEED TO TAKE INTO ACCOUNT LOCAL BC PRIORITIES AS WELL."
  "   - NEED TO CREATE SOMETHING TO INTERPOLATE THE BOUNDARY FACE SOLUTION"
  "     FOR THE GHOST CELLS/FACES TO THE NODES/EDGES FOR THE BOUNDARY FACES."
  "     INTERPOLATING THE INTERIOR SOLUTION TO THOSE NODE/EDGE POINTS WOULD"
  "     NOT BE CORRECT, FOR EXAMPLE INTERPOLATING THE HOST CELL SOLUTION"
  "     TO A BOUNDARY FACE THAT HAS A NO-SLIP WALL BOUNDARY CONDITION."
  "   - THIS WILL PROBABLY NEED TO BE DONE BEFORE EXCHANGING SINCE THE"
  "     SOLUTION FROM AN ADJACENT PARTITION MIGHT OVERRIDE THE LOCAL"
  "     SOLUTION OR VICE-VERSA."
  "   + IF THIS IS DONE BEFORE THE EXCHANGE, WE NEED TO MAKE SURE THE"
  "     AVERAGING CORRECTION IS BASED ON THE NUMBER OF BOUNDARY FACES"
  "     THAT CONTAIN THE NODE/EDGE AND NOT THE NUMBER OF INTERNAL GRID"
  "     CELLS."
  !
  " NEED TO FIGURE OUT HOW TO CONSTRUCT THE GLOBAL EDGE AND NODE"
  " CONNECTIVITY SO THAT WE CAN CREATE THE MPI DATA TYPES FOR EXCHANGING"
  " EDGE AND NODE SOLUTIONS."
  "   - NEED TO CREATE THE EDGE DERIVED TYPE ARRAY"
  "       + MAKE SURE THAT THE INDICES IN EDGE(:)%CPU(:)%EDGPTS STARTS AT 0"
  "       - MIGHT NEED TO FIX THE ORDERING OF THE EDGE POINT INDICES THAT"
  "         ARE STORED IN EDGE(:)%CPU(:)%EDGPTS. CURRENTLY IT JUST ASSUMES"
  "         THE ORDERING IS [0:NUMBER OF EDGE POINTS-1]."
  "   + NEED TO CREATE THE NODE DERIVED TYPE ARRAY"
  "   + NEED TO CREATE THE EDGE_MAP ARRAY"
  "   + NEED TO CREATE THE EXCH_NEUSP ARRAY"
  "   - "
  "   - "
  !
  " NEED TO FIGURE OUT HOW TO COLLECT THE INTERIOR, FACE, EDGE, AND NODE"
  " SOLUTIONS TOGETHER AS ONE SOLUTION FOR A CELL BEFORE WRITING IT TO"
  " THE CGNS SOLUTION FILE."
  "   - USE THE SAME SUBCONNECTIVITY FROM BEFORE INCLUDING THE ADD POINTS"
  "     VALUE DEPENDING ON THE CONTINUOUS SOLUTION INPUT PARAMETER. WE JUST"
  "     NEED A SIMPLE INTEGER ARRAY THAT MAPS THE VARIOUS FACE, EDGE, AND"
  "     NODE POINTS TO THE CORRECT LOCATIONS WITHIN THE ARRAY THAT WILL"
  "     STORE THE SOLUTION FOR THE CELL."
  !
  " NEED TO DEALLOCATE THE EDGEXYZ AND FACEXYZ (EXCEPT WHEN USING MMS)"
  " ARRAYS AFTER WRITING THE FIRST GRID CGNS FILE SINCE THEY WILL NO"
  " LONGER BE NEEDED."
  !
  " THE LOCAL AVERAGING COEFFICIENTS SHOULD ACCOUNT ONLY FOR LOCAL CELLS."
  " THE BND_NODES AND BND_EDGES ARRAYS WILL CORRECT THE LOCAL AVERAGING"
  " COEFFICIENTS TO THE CORRECT GLOBAL AVERAGING COEFFICIENT."
  !
  !
  !
contains
!
!###############################################################################
!
subroutine get_output_connectivity(bface,cell_geom, &
                                   nodes_of_cell,nodes_of_cell_ptr)
  !
  !.. Use Statements ..
  use connectivity_mod, only : find_cells_surrounding_each_node
  use connectivity_mod, only : find_cells_surrounding_each_cell
  !
  !.. Formal Arguments ..
  integer, dimension(:,:), intent(in) :: bface
  integer, dimension(:),   intent(in) :: cell_geom
  integer, dimension(:),   intent(in) :: nodes_of_cell
  integer, dimension(:),   intent(in) :: nodes_of_cell_ptr
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:)   :: cells_surr_cell
  integer, allocatable, dimension(:)   :: cells_surr_cell_ptr
  integer, allocatable, dimension(:)   :: edges_with_node
  integer, allocatable, dimension(:)   :: edges_with_node_ptr
  integer, allocatable, dimension(:,:) :: nodes_on_edge
  !
  integer, allocatable, dimension(:)   :: cells_with_node
  integer, allocatable, dimension(:)   :: cells_with_node_ptr
  integer, allocatable, dimension(:)   :: cells_with_edge
  integer, allocatable, dimension(:)   :: cells_with_edge_ptr
  integer, allocatable, dimension(:,:) :: cells_on_face
  !
continue
  !
  ! Set the module size parameters
  !
  lnode = maxval(nodes_of_cell)
  lfbnd = size(bface,dim=2)
  lcell = size(cell_geom) - lfbnd
  !
  ! Get the cells surrounding each node
  !
  call find_cells_surrounding_each_node(bface, &
                                        nodes_of_cell,nodes_of_cell_ptr, &
                                        cells_with_node,cells_with_node_ptr)
  !
  ! Get the cells surrounding each cell
  !
  call find_cells_surrounding_each_cell(bface,cell_geom, &
                                        nodes_of_cell,nodes_of_cell_ptr, &
                                        cells_with_node,cells_with_node_ptr, &
                                        cells_surr_cell,cells_surr_cell_ptr)
  !
  ! Get the nodes defining each edge as well as the edges surrounding each node
  !
  call find_nodes_on_each_edge(cell_geom,nodes_of_cell,nodes_of_cell_ptr, &
                               cells_with_node,cells_with_node_ptr, &
                               edges_with_node,edges_with_node_ptr, &
                               nodes_on_edge)
  !
  ! Get the cells that surround each edge
  !
  call find_cells_surrounding_each_edge(nodes_on_edge, &
                                        nodes_of_cell,nodes_of_cell_ptr, &
                                        cells_with_node,cells_with_node_ptr, &
                                        cells_with_edge,cells_with_edge_ptr)
  !
  ! Get the cells on each face
  !
  call find_cells_on_each_face(bface,cells_surr_cell,cells_surr_cell_ptr, &
                               cells_on_face)
  !
  ! Get the unit cell subconnectivities for use when writing solution file
  !
  call get_unit_cell_subconnectivity
  !
end subroutine get_output_connectivity
!
!###############################################################################
!
subroutine find_cells_on_each_face(bface,cells_surr_cell,cells_surr_cell_ptr, &
                                   cells_on_face)
  !
  ! Procedure that creates the face arrays. Look in
  ! the module geovar for definitions of each array.
  !
  !.. Formal Arguments ..
  integer,              dimension(:,:),    intent(in) :: bface
  integer, allocatable, dimension(:),   intent(inout) :: cells_surr_cell
  integer, allocatable, dimension(:),   intent(inout) :: cells_surr_cell_ptr
  integer, allocatable, dimension(:,:), intent(inout) :: cells_on_face
  !
  !.. Local Scalars ..
  integer :: c1,c2,n1,n2,nf,nc,ierr
  integer :: face_counter
  integer :: this_cell,adj_cell
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "find_cells_on_each_face"
  !
continue
  !
  ! ###########
  ! ###########
  ! ###  1  ###
  ! ###########
  ! ###########
  !
  ! Find the total number of faces
  ! and then allocate the face array
  !
  lface = lfbnd
  do nc = 1,lcell
    !
    n1 = cells_surr_cell_ptr(nc)+1
    n2 = cells_surr_cell_ptr(nc+1)
    !
    lface = lface + count( cells_surr_cell(n1:n2)>nc .and. &
                           cells_surr_cell(n1:n2)<=lcell )
    !
  end do
  !
  allocate ( cells_on_face(1:2,1:lface) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_on_face",1,ierr,error_message)
  !
  ! ###########
  ! ###########
  ! ###  2  ###
  ! ###########
  ! ###########
  !
  ! Get the cells on each boundary face
  !
  do nf = 1,lfbnd
    !
    cells_on_face(1,nf) = bface( 2,nf) ! Host  cell on boundary face
    cells_on_face(2,nf) = bface(10,nf) ! Ghost cell on boundary face
    !
  end do
  !
  ! ###########
  ! ###########
  ! ###  3  ###
  ! ###########
  ! ###########
  !
  ! Get the cells for the interior faces
  !
  face_counter = lfbnd
  !
  do this_cell = 1,lcell
    !
    c1 = cells_surr_cell_ptr(this_cell)+1
    c2 = cells_surr_cell_ptr(this_cell+1)
    !
    do nc = c1,c2
      !
      adj_cell = cells_surr_cell(nc)
      !
      if (adj_cell>this_cell .and. adj_cell<=lcell) then
        !
        face_counter = face_counter + 1
        !
        ! Store the cells on this interior face
        !
        cells_on_face(1,face_counter) = this_cell
        cells_on_face(2,face_counter) = adj_cell
        !
      end if
      !
    end do
    !
  end do
  !
  ! We are finished with cells_surr_cell
  ! and cells_surr_cell_ptr so deallocate
  !
  deallocate ( cells_surr_cell , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_surr_cell",2,ierr,error_message)
  !
  deallocate ( cells_surr_cell_ptr , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_surr_cell_ptr",2,ierr,error_message)
  !
end subroutine find_cells_on_each_face
!
!###############################################################################
!
subroutine collect_solution_for_output(var_sp,var_fp,var_ep,var_np)
  !
  !.. Use Statements ..
  use geovar,  only : cell
  use flowvar, only : usp,dusp
  use flowvar, only : faceusp,facedusp
  use flowvar, only : nodeusp,nodedusp
  use flowvar, only : edgeusp,edgedusp
  !
  !.. Formal Arguments ..
  integer,                 intent(in) :: plot_function
  real(wp), dimension(:), intent(out) :: var
  !
  !.. Local Scalars ..
  integer :: k,l,n,nc,ne,nf,nn,np1,np2
  integer :: face_idx,edge_idx,node_idx
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "collect_solution_for_output"
  !
continue
  !
  cell_loop: do nc = 1,ncell
    !
    np1 = cell(nc)%beg_sp
    np2 = cell(nc)%end_sp
    !
    do n = np1,np2
      !
      l = output_sol_pts(n)
      !
      var(l) = compute_output_variable(plot_function,usp(:,n),dusp(:,:,n))
      !
    end do
    !
  end do cell_loop
  !
  !
  !
  face_loop: do nf = 1,size(output_face_pts)
    !
    face_idx = output_face_pts(nf)%face_idx
    !
    do k = 1,size(output_face_pts(nf)%out_idx)
      !
      l = output_face_pts(nf)%out_idx(k)
      !
      var(l) = compute_output_variable(plot_function, &
                                       faceusp(face_idx)%v(:,k,1), &
                                       facedusp(face_idx)%v(:,:,k,1))
      !
    end do
    !
  end do face_loop
  !
  !
  !
  edge_loop: do ne = 1,size(output_edge_pts)
    !
    edge_idx = output_edge_pts(ne)%edge_idx
    !
    do k = 1,size(output_edge_pts(ne)%out_idx)
      !
      l = output_edge_pts(ne)%out_idx(k)
      !
      var(l) = compute_output_variable(plot_function, &
                                       edgeusp(edge_idx)%v(:,k), &
                                       edgedusp(edge_idx)%v(:,:,k))
      !
    end do
    !
  end do edge_loop
  !
  !
  !
  node_loop: do nn = 1,size(output_node_pts)
    !
    node_idx = output_node_pts(nn)%node_idx
    !
    l = output_node_pts(nn)%out_idx
    !
    var(l) = compute_output_variable(plot_function, &
                                     nodeusp(node_idx)%v(:), &
                                     nodedusp(node_idx)%v(:,:))
    !
  end do node_loop
  !
end subroutine collect_solution_for_output
!
!###############################################################################
!
subroutine get_face_solutions
  !
  !.. Use Statements ..
  use geovar,    only : nfbnd,nface
  use flowvar,   only : faceusp,facedusp
  !
  use flux_mod,  only : interp_cv_to_flxpts
  use flux_mod,  only : interp_dusp_to_flxpts
  use module_bc, only : get_faceusp_bc
  use module_bc, only : get_facedusp_bc
  !
  !.. Local Scalars ..
  integer :: nf
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_face_solutions"
  !
continue
  !
  call faceusp%zero
  call facedusp%zero
  !
  call interp_cv_to_flxpts
  !
  call get_faceusp_bc
  !
  do nf = 1,nfbnd
    faceusp(nf)%v(:,:,1) = faceusp(nf)%v(:,:,2)
  end do
  !
  do nf = nfbnd+1,nface
    faceusp(nf)%v(:,:,1) = half * (faceusp(nf)%v(:,:,1) + faceusp(nf)%v(:,:,2))
  end do
  !
  !
  !
  call interp_dusp_to_flxpts
  !
  call get_facedusp_bc
  !
  do nf = 1,nfbnd
    facedusp(nf)%v(:,:,:,1) = facedusp(nf)%v(:,:,:,2)
  end do
  !
  do nf = nfbnd+1,nface
    facedusp(nf)%v(:,:,:,1) = half * (facedusp(nf)%v(:,:,:,1) + &
                                      facedusp(nf)%v(:,:,:,2))
  end do
  !
end subroutine get_face_solutions
!
!###############################################################################
!
subroutine get_edge_and_node_solutions
  !
  !.. Use Statements ..
  use order_mode, only : maxSP
  use geovar,     only : nr,ncell,cell
  use flowvar,    only : usp,dusp
  use flowvar,    only : edgeusp,edgedusp
  use flowvar,    only : nodeusp,nodedusp
  use interpolation_mod, only : interp
  !
  use parallel_mod, only : exch_edge_and_node_solutions
  !
  !.. Local Scalars ..
  integer :: nc,n1,n2,np,k,l,m,n,ne,nn,kf,lf
  integer :: this_geom,this_order,edge_order
  real(wp) :: ave
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nq)              :: upt
  real(wp), dimension(1:nr,1:nq)         :: dupt
  real(wp), dimension(1:maxSP,1:nq)      :: tusp
  real(wp), dimension(1:maxSP,1:nr,1:nq) :: tdusp
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_edge_and_node_solutions"
  !
continue
  !
  ! Zero out the nodeusp, nodedusp, edgeusp, and edgedusp arrays
  !
  call nodeusp%zero
  call nodedusp%zero
  !
  call edgeusp%zero
  call edgedusp%zero
  !
  ! Loop over all the cells on the current processor and compute the
  ! contribution to the edge and node solutions from each cell
  !
  cell_loop: do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
    n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
    np = n2-n1+1         ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! Create transpose of the usp array to reduce cache misses
    !
    do k = 1,np
      do m = 1,nq
        tusp(k,m) = usp(m,k+n1-1)
      end do
    end do
    !
    ! Create transpose of the dusp array to reduce cache misses
    !
    do k = 1,np
      do m = 1,nq
        do l = 1,nr
          tusp(k,l,m) = dusp(l,m,k+n1-1)
        end do
      end do
    end do
    !
    ! Interpolate the solution and gradients to the edge points
    ! of the cell and add to the edgeusp and edgedusp arrays
    !
    cell_edges: do n = 1,size(cell(nc)%edge)
      !
      ne = cell(nc)%edge(n)%idx
      ave = cell(nc)%edge(n)%ave_coef
      edge_order = cell(nc)%edge(n)%order
      !
      edge_pts: do kf = 1,size(cell(nc)%edge(n)%cell_idx)
        !
        lf = cell(nc)%edge(n)%cell_idx(kf)
        !
        do m = 1,nq
          do l = 1,nr
            dupt(l,m) = interp(this_geom,this_order)% &
                                  toEdge(edge_order)% &
                        dot(lf,tdusp(1:np,l,m))
          end do
        end do
        !
        do m = 1,nq
          upt(m) = interp(this_geom,this_order)% &
                             toEdge(edge_order)% &
                   dot(lf,tusp(1:np,m))
        end do
        !
        edgeusp(ne)%v(1:nq,kf) = edgeusp(ne)%v(1:nq,kf) + ave*upt(1:nq)
        !
        edgedusp(ne)%v(1:nr,1:nq,kf) = edgedusp(ne)%v(1:nr,1:nq,kf) + &
                                       ave*dupt(1:nr,1:nq)
        !
      end do edge_pts
      !
    end do cell_edges
    !
    ! Interpolate the solution and gradients to the nodes
    ! of the cell and add to the nodeusp and nodedusp arrays
    !
    cell_nodes: do n = 1,size(cell(nc)%node)
      !
      nn  = cell(nc)%node(n)%idx
      ave = cell(nc)%node(n)%ave_coef
      lf  = n
      !
      do m = 1,nq
        do l = 1,nr
          dupt(l,m) = interp(this_geom,this_order)% &
                                            toNode% &
                      dot(lf,tdusp(1:np,l,m))
        end do
      end do
      !
      do m = 1,nq
        upt(m) = interp(this_geom,this_order)% &
                                       toNode% &
                 dot(lf,tusp(1:np,m))
      end do
      !
      nodeusp(nn)%v(1:nq) = nodeusp(nn)%v(1:nq) + ave*upt(1:nq)
      !
      nodedusp(nn)%v(1:nr,1:nq) = nodedusp(nn)%v(1:nr,1:nq) + &
                                  ave*dupt(1:nr,1:nq)
      !
    end do cell_nodes
    !
  end do cell_loop
  !
  !
  !
  if (ncpu > 1) then
    call exch_edge_and_node_solutions(nodeusp,nodedusp,edgeusp,edgedusp)
  end if
  !
end subroutine get_edge_and_node_solutions
!
!###############################################################################
!
subroutine collect_cell_subconnectivity
  !
  !.. Use Statements ..
  use order_mod, only : gmin => n_min_geom
  use order_mod, only : gmax => n_max_geom
  use order_mod, only : omin => n_min_order
  use order_mod, only : omax => n_max_order
  !
  !.. Local Scalars ..
  integer :: nop,nsp,nfp,nep,nnp
  integer :: mop,msp,mfp,mep,mnp
  !
  !.. Local Arrays ..
  integer, dimension(1:ncpu) :: cpu_mxpts
  !
  !.. Local Allocatable Arrays ..
  type(pnt_map_t), allocatable :: map(:,:)
  !
  integer, allocatable :: cell_points_ptr(:)
  integer, allocatable :: cell_points(:)
  integer, allocatable :: cell_pts(:)
  integer, allocatable :: face_pt(:)
  integer, allocatable :: edge_pt(:)
  integer, allocatable :: node_pt(:)
  integer, allocatable :: remaining_pts(:)
  !
  logical, allocatable :: my_faces(:)
  logical, allocatable :: my_edges(:)
  logical, allocatable :: my_nodes(:)
  logical, allocatable :: node_already_done(:)
  logical, allocatable :: edge_already_done(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "collect_cell_subconnectivity"
  !
continue
  !
  ! Initialize the maximum number of output points to zero for all processors
  !
  cpu_mxpts = 0
  !
  allocate ( map(gmin:gmax,omin:omax) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"map",1,ierr,error_message)
  !
  !##########################################################################
  !##########################################################################
  !###                                                                    ###
  !###  1. First create the integer pointer array for the indices of all  ###
  !###     the points within each cell on this partition.                 ###
  !###                                                                    ###
  !##########################################################################
  !##########################################################################
  !
  allocate ( cell_points_ptr(1:ncell+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_points_ptr",1,ierr,error_message)
  !
  ! Loop through each cell and count the number of output points
  ! for each cell on this partition
  !
  mfp = 0
  mep = 0
  mnp = 0
  do nc = 1,ncell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Create the mapping arrays if they havent already been created for
    ! this combination of cell geometry and solution order
    !
    if (.not. allocated(map(this_geom,this_order)%sp)) then
      call map_to_output_points(this_geom,this_order,map(this_geom,this_order))
    end if
    !
    ! Total number of solution points for this cell
    !
    nsp = size( map(this_geom,this_order)%sp )
    msp = max(msp,nsp)
    !
    ! Total number of face points for this cell
    !
    nfp = size( map(this_geom,this_order)%fp )
    mfp = max(mfp,nfp)
    !
    ! Total number of edge points for this cell
    !
    nep = size( map(this_geom,this_order)%ep )
    mep = max(mep,nep)
    !
    ! Total number of node points for this cell
    !
    nnp = size( map(this_geom,this_order)%np )
    mnp = max(mnp,nnp)
    !
    ! Total number of output points for this cell
    !
    nop = nsp + nfp + nep + nnp
    mop = max(mop,nop)
    !
    ! Add the total number of output points for this cell to
    ! the cell_points_ptr array
    !
    cell_points_ptr(nc+1) = cell_points_ptr(nc) + nop
    !
  end do
  !
  ! Save the total number of number of output points for this processor
  !
  cpu_mxpts(mypnum+1) = last(cell_points_ptr)
  !
  ! Collect the total number of output points from all processors
  !
  if (ncpu > 1) then
    call mpi_allgather(MPI_IN_PLACE,0_int_mpi,MPI_DATATYPE_NULL, &
                       cpu_mxpts,1_int_mpi,mpi_inttyp, &
                       mpi_comm_world,mpierr)
  end if
  !
  ! Create the point offset for this processor
  !
  my_pts_offset = sum( cpu_mxpts(1:mypnum) )
  !
  !#########################################################################
  !#########################################################################
  !###                                                                   ###
  !###  2. Perform the necessary allocations and initializations before  ###
  !###     creating the initial point indices for each cell on this      ###
  !###     partition.                                                    ###
  !###                                                                   ###
  !#########################################################################
  !#########################################################################
  !
  ! Allocate the cell_points array which will contain the indices for the output
  ! points within each cell. This will initially contain the indices for each
  ! cell with all colocated node, edge, and face points on partition boundaries
  ! having a unique index from each cell containing the node, edge, or face.
  ! After exchanging the unique indices between adjacent partitions, the indices
  ! for these colocated points will be reduced to a single index that is shared
  ! by all cells that contain the point.
  !
  allocate ( cell_points(1:last(cell_points_ptr)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_points",1,ierr,error_message)
  !
  ! Allocate the arrays to mark nodes and edges that have already been completed
  !
  nn = size(nodeusp)
  allocate ( node_already_done(1:nn) , source=fals , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"node_already_done",1,ierr,error_message)
  !
  ne = size(edgeusp)
  allocate ( edge_already_done(1:ne) , source=fals , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"edge_already_done",1,ierr,error_message)
  !
  ! Allocate the cell_pts array to temporarily store the point indices of
  ! the current cell thats being operated on
  !
  allocate ( cell_pts(1:mop) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_pts",1,ierr,error_message)
  !
  allocate ( face_pt(1:mfp) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"face_pt",1,ierr,error_message)
  !
  allocate ( edge_pt(1:mep) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"edge_pt",1,ierr,error_message)
  !
  allocate ( node_pt(1:mnp) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"node_pt",1,ierr,error_message)
  !
  ! Zero the nodeusp, nodedusp, edgeusp, edgedusp, and faceusp array
  ! before we start storing point indices in these arrays
  !
  call nodeusp%zero
  call nodedusp%zero
  call edgeusp%zero
  call edgedusp%zero
  call faceusp%zero
  !
  !###########################################################################
  !###########################################################################
  !###                                                                     ###
  !###  3. Loop through the cells on this partition and find the indices   ###
  !###     for the output points of each cell and use the node and edge    ###
  !###     mask arrays that were just allocated to make sure that we only  ###
  !###     assign a point index for the first node or edge point found     ###
  !###     for the nodes and edges with multiple colocated points and      ###
  !###     that are only found on this partition.                          ###
  !###                                                                     ###
  !###########################################################################
  !###########################################################################
  !
  pack_cell_pts: do nc = 1,ncell
    !
    ! First initialize cell_pts to 0
    !
    cell_pts = 0
    face_pt  = 0
    edge_pt  = 0
    node_pt  = 0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    nop1 = cell_points_ptr(nc)+1
    nop2 = cell_points_ptr(nc+1)
    nop  = nop2-nop1+1
    !
    ! Total number of solution points for this cell
    !
    nsp = size( map(this_geom,this_order)%sp )
    !
    ! Find the total number of face points for this cell and create the
    ! map from the face point lf to the correct location in cell_pts
    !
    nfp = size( map(this_geom,this_order)%fp )
    !
    face_pt(1:nfp) = nsp + intseq(1,nfp)
    !
    ! Find the total number of edge points for this cell and create the
    ! map from the edge point lf to the correct location in cell_pts
    !
    nep = size( map(this_geom,this_order)%ep )
    !
    edge_pt(1:nep) = nsp + nfp + intseq(1,nep)
    !
    ! Find the total number of node points for this cell and create the
    ! map from the node point lf to the correct location in cell_pts
    !
    nnp = size( map(this_geom,this_order)%np )
    !
    node_pt(1:nnp) = nsp + nfp + nep + intseq(1,nnp)
    !
    ! Create the initial indices for the solution points in this cell
    !
    cell_pts(1:nsp) = nop1 + map(this_geom,this_order)%sp
    !
    ! Create the initial indices for the face points of this cell
    !
    cell_pts(face_pt) = nop1 + map(this_geom,this_order)%fp
    !
    ! Create the initial indices for the edge points of this cell
    !
    cell_pts(edge_pt) = nop1 + map(this_geom,this_order)%ep
    !
    ! Create the initial indices for the node points of this cell
    !
    cell_pts(node_pt) = nop1 + map(this_geom,this_order)%np
    !
    ! Adjust cell_pts for the offset of this processor
    !
    cell_pts(1:nop) = cell_pts(1:nop) + my_pts_offset
    !
    ! Pack up the face point indices of this cell into faceusp
    !
    pack_cell_face_pts: do n = 1,size(cell(nc)%face)
      !
      nf = cell(nc)%face(n)%idx
      sf = cell(nc)%face(n)%side
      !
      do kf = 1,size(cell(nc)%face(n)%cell_idx)
        !
        lf = cell(nc)%face(n)%cell_idx(kf)
        !
        faceusp(nf)%v(1,kf,sf) = real( cell_pts( face_pt(lf) ) , kind=wp )
        !
      end do
      !
    end do pack_cell_face_pts
    !
    ! Pack up the edge point indices of this cell into edgeusp
    !
    pack_cell_edge_pts: do n = 1,size(cell(nc)%edge)
      !
      ne = cell(nc)%edge(n)%idx
      !
      if (edge_already_done(ne)) cycle pack_cell_edge_pts
      !
      do kf = 1,size(cell(nc)%edge(n)%cell_idx)
        !
        lf = cell(nc)%edge(n)%cell_idx(kf)
        !
        edgeusp(ne)%v(1,kf) = real( cell_pts( edge_pt(lf) ) , kind=wp )
        !
      end do
      !
      edge_already_done(ne) = true
      !
    end do pack_cell_edge_pts
    !
    ! Pack up the node point indices of this cell into nodeusp
    !
    pack_cell_node_pts: do n = 1,size(cell(nc)%node)
      !
      nn = cell(nc)%node(n)%idx
      !
      if (node_already_done(nn)) cycle pack_cell_node_pts
      !
      lf = n
      !
      nodeusp(nn)%v(1) = real( cell_pts( node_pt(lf) ) , kind=wp )
      !
      node_already_done(nn) = true
      !
    end do pack_cell_node_pts
    !
    ! Save the output points for this cell, in correct order, and
    ! before any corrections for duplicates
    !
    cell_points(nop1:nop2) = cell_pts(1:nop)
    !
  end do pack_cell_pts
  !
  if (allocated(node_already_done)) then
    deallocate ( node_already_done , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"node_already_done",2,ierr,error_message)
  end if
  !
  if (allocated(edge_already_done)) then
    deallocate ( edge_already_done , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edge_already_done",2,ierr,error_message)
  end if
  !
  !###########################################################################
  !###########################################################################
  !###                                                                     ###
  !###  4. If using multiple partitions/processors, exchange the point     ###
  !###     indices for any nodes, edges, and faces, that are on partition  ###
  !###     boundaries.                                                     ###
  !###                                                                     ###
  !###########################################################################
  !###########################################################################
  !
  if (ncpu > 1) then
    call exch_connectivity(nodeusp,nodedusp,edgeusp,edgedusp,faceusp)
  end if
  !
  !##########################################################################
  !##########################################################################
  !###                                                                    ###
  !###  5. Loop through the cells on this partition and unpack the point  ###
  !###     indices from the nodeusp, edgeusp, and faceusp arrays.         ###
  !###                                                                    ###
  !##########################################################################
  !##########################################################################
  !
  do nc = 1,ncell
    !
    cell_pts = 0
    face_pt  = 0
    edge_pt  = 0
    node_pt  = 0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    nop1 = cell_points_ptr(nc)+1
    nop2 = cell_points_ptr(nc+1)
    nop  = nop2-nop1+1
    !
    ! Total number of solution points for this cell
    !
    nsp = size( map(this_geom,this_order)%sp )
    !
    ! Find the total number of face points for this cell and create the
    ! map from the face point lf to the correct location in cell_pts
    !
    nfp = size( map(this_geom,this_order)%fp )
    !
    face_pt(1:nfp) = nsp + intseq(1,nfp)
    !
    ! Find the total number of edge points for this cell and create the
    ! map from the edge point lf to the correct location in cell_pts
    !
    nep = size( map(this_geom,this_order)%ep )
    !
    edge_pt(1:nep) = nsp + nfp + intseq(1,nep)
    !
    ! Find the total number of node points for this cell and create the
    ! map from the node point lf to the correct location in cell_pts
    !
    nnp = size( map(this_geom,this_order)%np )
    !
    node_pt(1:nnp) = nsp + nfp + nep + intseq(1,nnp)
    !
    ! Copy the output points for this cell, in correct order, and
    ! before any corrections for duplicates
    !
    cell_pts(1:nop) = cell_points(nop1:nop2)
    !
    ! Unpack the face point indices of this cell from faceusp
    !
    unpack_cell_face_pts: do n = 1,size(cell(nc)%face)
      !
      nf = cell(nc)%face(n)%idx
      !
      nfp = 
      !
      do kf = 1,size(cell(nc)%face(n)%cell_idx)
        !
        lf = cell(nc)%face(n)%cell_idx(kf)
        !
        cell_pts( face_pt(lf) ) = nint( faceusp(nf)%v(1,kf,1) , &
                                        kind=kind(cell_pts) )
        !
      end do
      !
    end do unpack_cell_face_pts
    !
    ! Unpack the edge point indices of this cell from edgeusp
    !
    unpack_cell_edge_pts: do n = 1,size(cell(nc)%edge)
      !
      ne = cell(nc)%edge(n)%idx
      !
      do kf = 1,size(cell(nc)%edge(n)%cell_idx)
        !
        lf = cell(nc)%edge(n)%cell_idx(kf)
        !
        cell_pts( edge_pt(lf) ) = nint( edgeusp(ne)%v(1,kf) , &
                                        kind=kind(cell_pts) )
        !
      end do
      !
    end do unpack_cell_edge_pts
    !
    ! Unpack the node point indices of this cell from nodeusp
    !
    unpack_cell_node_pts: do n = 1,size(cell(nc)%node)
      !
      nn = cell(nc)%node(n)%idx
      !
      lf = n
      !
      cell_pts( node_pt(lf) ) = nint( nodeusp(nn)%v(1) , kind=kind(cell_pts) )
      !
    end do unpack_cell_node_pts
    !
    ! Copy the unpacked point indices into cell_points
    !
    cell_points(nop1:nop2) = cell_pts(1:nop)
    !
  end do
  !
  !#############################################################################
  !#############################################################################
  !###                                                                       ###
  !###  6.                                                                   ###
  !###                                                                       ###
  !###                                                                       ###
  !###                                                                       ###
  !###                                                                       ###
  !#############################################################################
  !#############################################################################
  !
  nf = size(faceusp)
  allocate ( my_faces(1:nf) , source=fals , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"my_faces",1,ierr,error_message
  !
  ne = size(edgeusp)
  allocate ( my_edges(1:ne) , source=fals , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"my_edges",1,ierr,error_message
  !
  nn = size(nodeusp)
  allocate ( my_nodes(1:nn) , source=fals , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"my_nodes",1,ierr,error_message
  !
  bad_indices = 0
  !
  ! Find the face points for which this processor is responsible
  !
  find_my_faces: do nf = 1,size(faceusp)
    !
    n = size(faceusp(nf)%v,dim=2)
    m = count(faceusp(nf)%v(1,:,1) > zero)
    !
    if (n == m) then
      !
      my_faces(nf) = true
      !
    else if (m > 0) then
      !
      if (bad_indices == 0) write (iout,1)
      !
      bad_indices = bad_indices + 1
      write (iout,2) nf,m,n
      !
    end if
    !
  end do find_my_faces
  !
  ! Find the edge points for which this processor is responsible
  !
  find_my_edges: do ne = 1,size(edgeusp)
    !
    n = size(edgeusp(ne)%v,dim=2)
    m = count(edgeusp(ne)%v(1,:) > zero)
    !
    if (n == m) then
      !
      my_edges(ne) = true
      !
    else if (m > 0) then
      !
      if (bad_indices == 0) write (iout,1)
      !
      bad_indices = bad_indices + 1
      write (iout,3) ne,m,n
      !
    end if
    !
  end do find_my_edges
  !
  ! Find the node points for which this processor is responsible
  !
  find_my_nodes: do nn = 1,size(nodeusp)
    !
    if (nodeusp(nn)%v(1) > zero) then
      my_nodes(nn) = true
    end if
    !
  end do find_my_nodes
  !
  !############################################################################
  !############################################################################
  !###                                                                      ###
  !###  7. Since all the colocated points with multiple unique indices for  ###
  !###     the same point were reduced to a single index, there is now a    ###
  !###     large number of point indices which are no longer used,          ###
  !###     creating 'holes' in the indexing. Now we will go through and     ###
  !###     remove these holes in the indexing so the point indices are a    ###
  !###     continuous interval.                                             ###
  !###                                                                      ###
  !############################################################################
  !############################################################################
  !
  ! Find the maximum point index
  !
  mxpt = maxval(cell_points)
  !
  ! Collectively find the maximum point index from all processors
  !
  if (ncpu > 1) then
    call mpi_allreduce(MPI_IN_PLACE,mxpt,1_int_mpi,mpi_inttyp, &
                       MPI_MAX,mpi_comm_world,mpierr)
  end if
  !
  ! Allocate the array to mark the points on each processor that still remain
  !
  allocate ( remaining_pts(1:mxpt) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"remaining_pts",1,ierr,error_message)
  !
  ! Mark all the remaining points on this processor
  ! NOTE: We need to use the absolute value of cell_points because some
  !       of the indices are still negative to mark points for which an
  !       adjacent processor is responsible for writing the solution at
  !       the point.
  !
  remaining_pts( abs(cell_points(:)) ) = 1
  !
  ! Collectively reduce the remaining points on all processors
  !
  if (ncpu > 1) then
    call mpi_allreduce(MPI_IN_PLACE,remaining_pts,int(mxpt,kind=int_mpi), &
                       mpi_inttyp,MPI_MAX,mpi_comm_world,mpierr)
  end if
  !
  ! Convert remaining_pts into an array that maps a point index in cell_points
  ! into the equivalent point index if all non-used point indices were removed
  !
 !do n = size(remaining_pts),1,-1
 !  if (remaining_pts(n) /= 0) then
 !    remaining_pts(n) = sum( remaining_pts(1:n-1) )
 !  end if
 !end do
  mxpt = 0
  do n = 1,size(remaining_pts)
    if (remaining_pts(n) /= 0) then
      mxpt = mxpt + 1
      remaining_pts(n) = mxpt
    end if
  end do
  !
  ! Now redo the point indices in cell_points to account
  ! for only the remaining point indices
  ! NOTE: We need to make sure that we keep the negative sign for indices
  !       within cell_points that mark points for which an adjacent partition
  !       is responsible for writing the solution at the point.
  !
  do n = 1,size(cell_points)
    cell_points(n) = sign( remaining_pts( abs(cell_points(n)) ) , &
                           cell_points(n) )
  end do
  !
  if (allocated(remaining_pts)) then
    deallocate ( remaining_pts , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"remaining_pts",2,ierr,error_message)
  end if
  !
  !###########################################################################
  !###########################################################################
  !###                                                                     ###
  !###  8. Create the arrays that map all the points that the current      ###
  !###     processor is responsible for writing to the CGNS solution file  ###
  !###     to the correct location within the array that collects the      ###
  !###     plot function values at each of these points and is then sent   ###
  !###     into the CGNS function for writing the solution.                ###
  !###                                                                     ###
  !###########################################################################
  !###########################################################################
  !
  allocate ( output_sol_pts(1:n_solpts) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"output_sol_pts",1,ierr,error_message)
  !
  nf = count( my_faces )
  allocate ( output_face_pts(1:nf) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"output_face_pts",1,ierr,error_message)
  !
  ne = count( my_edges )
  allocate ( output_edge_pts(1:ne) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"output_edge_pts",1,ierr,error_message)
  !
  nn = count( my_nodes )
  allocate ( output_node_pts(1:nn) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"output_node_pts",1,ierr,error_message)
  !
  !
  do nc = 1,ncell
    !
    cell_pts = 0
    face_pt  = 0
    edge_pt  = 0
    node_pt  = 0
    !
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    nop1 = cell_points_ptr(nc)+1
    nop2 = cell_points_ptr(nc+1)
    nop  = nop2-nop1+1
    !
    ! Total number of solution points for this cell
    !
    nsp = size( map(this_geom,this_order)%sp )
    !
    ! Find the total number of face points for this cell and create the
    ! map from the face point lf to the correct location in cell_pts
    !
    nfp = size( map(this_geom,this_order)%fp )
    !
    face_pt(1:nfp) = nsp + intseq(1,nfp)
    !
    ! Find the total number of edge points for this cell and create the
    ! map from the edge point lf to the correct location in cell_pts
    !
    nep = size( map(this_geom,this_order)%ep )
    !
    edge_pt(1:nep) = nsp + nfp + intseq(1,nep)
    !
    ! Find the total number of node points for this cell and create the
    ! map from the node point lf to the correct location in cell_pts
    !
    nnp = size( map(this_geom,this_order)%np )
    !
    node_pt(1:nnp) = nsp + nfp + nep + intseq(1,nnp)
    !
    ! Make a copy of the output points for this cell
    !
    cell_pts(1:nop) = cell_points(nop1:nop2)
    !
    ! Get the output point indices for the solution points of this cell
    !
    output_sol_pts(n1:n2) = cell_pts(1:nsp)
    !
    ! Get the output point indices for the face points of this cell    
    !
    ifc = 0
    create_output_face_pts: do n = 1,size(cell(nc)%face)
      !
      nf = cell(nc)%face(n)%idx
      !
      ! my_faces(nf) will be true if the current processor is responsible
      ! for the solution at the points on this face AND this face hasnt
      ! already been added to output_face_pts
      !
      if (.not. my_faces(nf)) cycle create_output_face_pts
      !
      ! Increment the face counter by one and then set
      ! my_faces(nf) to false so that this face is not
      ! found a second time by the other cell on this face
      !
      ifc = ifc + 1
      my_faces(nf) = fals
      !
      ! Save this face index to output_face_pts
      !
      output_face_pts(ifc)%face_idx = nf
      !
      ! Get the number of face points for this face
      !
      nfp = size(cell(nc)%face(n)%cell_idx)
      !
      ! Allocate the out_idx component of output_face_pts for this face
      !
      allocate ( output_face_pts(ifc)%out_idx(1:nfp) , source=0 , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) "face",ifc
      call alloc_error(pname,array_name,1,ierr,error_message)
      !
      do kf = 1,nfp
        !
        lf = cell(nc)%face(n)%cell_idx(kf)
        !
        output_face_pts(ifc)%out_idx(kf) = cell_pts( face_pt(lf) )
        !
      end do
      !
    end do create_output_face_pts
    !
    ! Get the output point indices for the edge points of this cell    
    !
    iec = 0
    create_output_edge_pts: do n = 1,size(cell(nc)%edge)
      !
      ne = cell(nc)%edge(n)%idx
      !
      ! my_edges(ne) will be true if the current processor is responsible
      ! for the solution at the points on this edge AND this edge hasnt
      ! already been added to output_edge_pts
      !
      if (.not. my_edges(ne)) cycle create_output_edge_pts
      !
      ! Increment the edge counter by one and then set
      ! my_edges(ne) to false so that this edge is not
      ! found a second time by other cells with this edge
      !
      iec = iec + 1
      my_edges(ne) = fals
      !
      ! Save this edge index to output_edge_pts
      !
      output_edge_pts(iec)%edge_idx = ne
      !
      ! Get the number of edge points for this edge
      !
      nep = size(cell(nc)%edge(n)%cell_idx)
      !
      ! Allocate the out_idx component of output_edge_pts for this edge
      !
      allocate ( output_edge_pts(iec)%out_idx(1:nep) , source=0 , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) "edge",iec
      call alloc_error(pname,array_name,1,ierr,error_message)
      !
      do kf = 1,nep
        !
        lf = cell(nc)%edge(n)%cell_idx(kf)
        !
        output_edge_pts(iec)%out_idx(kf) = cell_pts( edge_pt(lf) )
        !
      end do
      !
    end do create_output_edge_pts
    !
    ! Get the output point indices for the nodes of this cell    
    !
    inc = 0
    create_output_node_pts: do n = 1,size(cell(nc)%node)
      !
      nn = cell(nc)%node(n)%idx
      !
      ! my_nodes(nn) will be true if the current processor is responsible
      ! for the solution at this node AND this node hasnt already been
      ! added to output_node_pts
      !
      if (.not. my_nodes(nn)) cycle create_output_node_pts
      !
      ! Increment the node counter by one and then set
      ! my_nodes(nn) to false so that this node is not
      ! found a second time by other cells with this node
      !
      inc = inc + 1
      my_nodes(nn) = fals
      !
      ! Save this node index to output_node_pts
      !
      output_node_pts(inc)%node_idx = nn
      !
      ! Save the output index for this node to output_node_pts
      !
      lf = n
      !
      output_node_pts(inc)%out_idx = cell_pts( node_pt(lf) )
      !
      !
    end do create_output_node_pts
    !
  end do

  !
  !#######################################################
  !#######################################################
  !###                                                 ###
  !###  9. Deallocate the local arrays before leaving  ###
  !###                                                 ###
  !#######################################################
  !#######################################################
  !
  if (allocated(cell_points)) then
    deallocate ( cell_points , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cell_points",2,ierr,error_message)
  end if
  !
  if (allocated(cell_points_ptr)) then
    deallocate ( cell_points_ptr , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cell_points_ptr",2,ierr,error_message)
  end if
  !
  if (allocated(my_faces)) then
    deallocate ( my_faces , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"my_faces",2,ierr,error_message)
  end if
  !
  if (allocated(my_edges)) then
    deallocate ( my_edges , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"my_edges",2,ierr,error_message)
  end if
  !
  if (allocated(my_nodes)) then
    deallocate ( my_nodes , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"my_nodes",2,ierr,error_message)
  end if
  !
  if (allocated(map)) then
    deallocate ( map , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"map",2,ierr,error_message)
  end if
  !
  if (allocated(cell_pts)) then
    deallocate ( cell_pts , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cell_pts",2,ierr,error_message)
  end if
  !
  if (allocated(face_pt)) then
    deallocate ( face_pt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"face_pt",2,ierr,error_message)
  end if
  !
  if (allocated(edge_pt)) then
    deallocate ( edge_pt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edge_pt",2,ierr,error_message)
  end if
  !
  if (allocated(node_pt)) then
    deallocate ( node_pt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"node_pt",2,ierr,error_message)
  end if
  !
end subroutine collect_cell_subconnectivity
!
!###############################################################################
!
pure function compute_output_variable(plot_function,upt,dupt)
  !
  !.. Formal Arguments ..
  type(plot_func_t),        intent(in) :: plot_function
  real(wp), dimension(:),   intent(in) :: upt
  real(wp), dimension(:,:), intent(in) :: dupt
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
continue
  !
  if (plot_function == plot_static_density) then
    !
    return_value = upt(nec)
    !
  else if (plot_function == plot_velocity_x) then
    !
    return_value = upt(nmx) / upt(nec)
    !
  else if (plot_function == plot_velocity_y) then
    !
    return_value = upt(nmy) / upt(nec)
    !
  else if (plot_function == plot_velocity_z) then
    !
    return_value = upt(nmz) / upt(nec)
    !
  else if (plot_function == plot_momentum_x) then
    !
    return_value = upt(nmx)
    !
  else if (plot_function == plot_momentum_y) then
    !
    return_value = upt(nmy)
    !
  else if (plot_function == plot_momentum_z) then
    !
    return_value = upt(nmz)
    !
  else if (plot_function == plot_static_pressure) then
    !
    return_value = pressure_cv( upt )
    !
  else if (plot_function == plot_static_temperature) then
    !
    return_value = temperature_cv( upt )
    !
  else if (plot_function == plot_mach_number) then
    !
    return_value = mach_number_cv( upt )
    !
  else if (plot_function == plot_entropy) then
    !
    return_value = entropy_cv( upt )
    !
  else if (plot_function == plot_vorticity_x) then
    !
    return_value = vorticity_dusp( upt , dupt , 1 )
    !
  else if (plot_function == plot_vorticity_y) then
    !
    return_value = vorticity_dusp( upt , dupt , 2 )
    !
  else if (plot_function == plot_vorticity_z) then
    !
    return_value = vorticity_dusp( upt , dupt , 3 )
    !
  end if
  !
end function compute_output_variable
!
!###############################################################################
!
subroutine map_to_output_points(this_geom,this_order,map)
  !
  !.. Formal Arguments ..
  integer,            intent(in) :: this_geom
  integer,            intent(in) :: this_order
  type(pnt_map_t), intent(inout) :: map
  !
  !.. Local Scalars ..
  integer :: ne,nf,nep,nfp,nnp,nsp,ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "map_to_output_points"
  !
continue
  !
  !
  !###################################
  !###################################
  !#####          FACES          #####
  !###################################
  !###################################
  !
  ! Number of solution points for this combination
  ! of cell geometry and solution order
  !
  nsp = Cell_Solpts(this_geom,this_order)
  !
  ! Number of faces for this cell
  !
  nfp = 0
  do nf = 1,geom_faces(this_geom)
    nfp = nfp + Cell_Solpts( cell_face_geom(nf,this_geom) , this_order )
  end do
  !
  ! Number of edges for this cell
  !
  nep = 0
  if (any(this_geom == Geom_3D)) then
    do ne = 1,geom_edges(this_geom)
      nep = nep + Cell_Solpts(Geom_Edge,this_order)
    end do
  end if
  !
  ! Number of nodes for this cell
  !
  nnp = geom_nodes(this_geom)
  !
  ! Allocate the different components of the map derived type
  !
  allocate ( map%sp(1:nsp) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"map%sp",1,ierr,error_message)
  !
  allocate ( map%fp(1:nfp) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"map%sp",1,ierr,error_message)
  !
  allocate ( map%ep(1:nep) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"map%sp",1,ierr,error_message)
  !
  allocate ( map%np(1:nnp) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"map%sp",1,ierr,error_message)
  !
  ! Compute the maps from the local indices to the output indices
  !
  select case (this_geom)
      !
    case (Geom_Tria)
      !
      map%sp = get_tria_solpt_indices(this_order)
      map%fp = get_tria_face_indices(this_order)
      map%np = get_tria_node_indices(this_order)
      !
    case (Geom_Quad)
      !
      map%sp = get_quad_solpt_indices(this_order)
      map%fp = get_quad_face_indices(this_order)
      map%np = get_quad_node_indices(this_order)
      !
    case (Geom_Tetr)
      !
      map%sp = get_tetr_solpt_indices(this_order)
      map%fp = get_tetr_face_indices(this_order)
      map%ep = get_tetr_edge_indices(this_order)
      map%np = get_tetr_node_indices(this_order)
      !
    case (Geom_Pyra)
      !
      map%sp = get_pyra_solpt_indices(this_order)
      map%fp = get_pyra_face_indices(this_order)
      map%ep = get_pyra_edge_indices(this_order)
      map%np = get_pyra_node_indices(this_order)
      !
    case (Geom_Pris)
      !
      map%sp = get_pris_solpt_indices(this_order)
      map%fp = get_pris_face_indices(this_order)
      map%ep = get_pris_edge_indices(this_order)
      map%np = get_pris_node_indices(this_order)
      !
    case (Geom_Hexa)
      !
      map%sp = get_hexa_solpt_indices(this_order)
      map%fp = get_hexa_face_indices(this_order)
      map%ep = get_hexa_edge_indices(this_order)
      map%np = get_hexa_node_indices(this_order)
      !
    case default
      !
      write (error_message,1)
      call stop_recon(abort,pname,error_message)
      !
  end select
  !
  ! Finally, adjust the mapping indices so that they start at 0 instead of 1
  !
  map%sp = map%sp - 1
  map%fp = map%fp - 1
  map%ep = map%ep - 1
  map%np = map%np - 1
  !
  ! Format Statements
  !
  1 format ("ERROR! Encountered invalid geometry when trying to compute ", &
            "the arrays to map local point indices to the corresponding ", &
            "output indices.")
  !
end subroutine map_to_output_points
!
!###############################################################################
!
pure function get_hexa_solpt_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, dimension(1:(order+1)**3) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,n,np1d,nt1d
  !
continue
  !
  return_value = 0
  !
  np1d = order + 1
  nt1d = np1d + 2
  !
  l = 0
  n = nt1d*nt1d
  !
  do k = 1,np1d
    n = n + nt1d ! skip all points along the j-min face
    do j = 1,np1d
      n = n + 1 ! skip the point for the i-min face
      do i = 1,np1d
        l = l + 1
        n = n + 1
        return_value(l) = n
      end do
      n = n + 1 ! skip the point for the i-max face
    end do
    n = n + nt1d ! skip all points along the j-max face
  end do
  !
end function get_hexa_solpt_indices
!
!###############################################################################
!
pure function get_hexa_node_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, dimension(1:8) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,ip,jp,kp,np1d,nt1d,nt2d
  !
continue
  !
  return_value = 0
  !
  np1d = order + 1
  nt1d = np1d + 2
  nt2d = nt1d**2
  !
  l = 0
  do k = 1,2
    kp = (k-1)*(nt1d-1)
    do j = 1,2
      jp = (j-1)*(nt1d-1)
      do i = 1,2
        ip = (i-1)*(nt1d-1) - (2*i-3)*jp
        l = l + 1
        return_value(l) = kp*nt2d + jp*nt1d + ip + 1
      end do
    end do
  end do
 !return_value(1) = 1
 !return_value(2) = nt1d
 !return_value(3) = nt1d*nt1d
 !return_value(4) = nt1d*nt1d - nt1d + 1
 !return_value(5) = nt1d*nt1d*(nt1d-1) + 1
 !return_value(6) = nt1d*nt1d*(nt1d-1) + nt1d
 !return_value(7) = nt1d*nt1d*(nt1d-1) + nt1d*nt1d
 !return_value(8) = nt1d*nt1d*(nt1d-1) + nt1d*nt1d - nt1d + 1

  !
end function get_hexa_node_indices
!
!###############################################################################
!
pure function get_hexa_face_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, dimension(1:6*(order+1)*(order+1)) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,ip,iexp,istride
  integer :: j,jp,jexp,jstride
  integer :: k,kp,kexp,kstride
  integer :: l,nf,np1d,nt1d
  !
  !.. Local Parameters ..
  integer, parameter :: dx = 0 ! The i/x/xi   direction
  integer, parameter :: dy = 1 ! The j/y/eta  direction
  integer, parameter :: dz = 2 ! The k/z/zeta direction
  !
  ! fexp(1,:) :  gives the stride exponent for the direction in which the local
  !              face indices are consecutive
  !
  ! fexp(2,:) :  gives the stride exponent for the direction in which the local
  !              face indices are not consecutive
  !
  integer, parameter :: fexp(2,6) = reshape( [dy,dx, dx,dy, &
                                              dy,dz, dz,dy, &
                                              dy,dz, dz,dx], shape(fexp) )
  !
continue
  !
  return_value = 0
  !
  np1d = order + 1
  nt1d = np1d + 2
  !
  l = 0
  do nf = 1,6
    !
    k = merge( 1 , nt1d , any(nf == [Xmin3D_Face,Ymin3D_Face,ZMin3D_Face]) )
    !
    iexp = fexp(1,nf)
    jexp = fexp(2,nf)
    kexp = 3 - iexp - jexp
    !
    istride = nt1d**iexp
    jstride = nt1d**jexp
    kstride = nt1d**kexp
    !
    kp = (k-1)*kstride
    do j = 2,np1d+1
      jp = (j-1)*jstride
      do i = 2,np1d+1
        ip = (i-1)*istride
        l = l + 1
        return_value(l) = kp + jp + ip + 1
      end do
    end do
    !
  end do
  !
end function get_hexa_face_indices
!
!###############################################################################
!
pure function get_hexa_edge_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result ..
  integer, dimension(1:12*(order+1)) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,l,ne,idir,np1d,nt1d
  integer :: ic,ip,iexp,istride
  integer :: jc,jp,jexp,jstride
  integer :: kc,kp,kexp,kstride
  !
  !.. Local Parameters ..
  integer, parameter :: dx = 0 ! The i/x/xi   direction
  integer, parameter :: dy = 1 ! The j/y/eta  direction
  integer, parameter :: dz = 2 ! The k/z/zeta direction
  !
  integer, parameter :: fwd = 0 ! The edge-pt indices are in increasing order
  integer, parameter :: rev = 1 ! The edge-pt indices are in decreasing order
  !
  integer, parameter :: dmin = 0 ! The direction is a min face
  integer, parameter :: dmax = 1 ! The direction is a max face
  !
  ! eexp(1,:) : gives the stride exponent for the direction of
  !             lowest  dimension that is constant for the edge
  !
  ! eexp(2,:) : gives the stride exponent for the direction of
  !             highest dimension that is constant for the edge
  !
  integer, parameter :: eexp(2,12) = reshape( [dy,dz, dx,dz, dy,dz, dx,dz,  &
                                               dx,dy, dx,dy, dx,dy, dx,dy,  &
                                               dy,dz, dx,dz, dy,dz, dx,dz], &
                                               shape(eexp) )
  !
  ! ecoef(1,:) : gives the coefficient that puts the edge-pt indices in
  !              increasing or decreasing order for the non-constant
  !              direction of the edge
  !
  ! ecoef(2,:) : gives the coefficient that adjusts the edge-pt
  !              indices based on whether the edge is on a
  !              min or max face for the direction eexp(1,:)
  !
  ! ecoef(3,:) : gives the coefficient that adjusts the edge-pt
  !              indices based on whether the edge is on a
  !              min or max face for the direction eexp(2,:)
  !
  integer, parameter :: ecoef(3,12) = reshape( [fwd,dmin,dmin, fwd,dmax,dmin, &
                                                rev,dmax,dmin, rev,dmin,dmin, &
                                                fwd,dmin,dmin, fwd,dmax,dmin, &
                                                fwd,dmax,dmax, fwd,dmin,dmax, &
                                                fwd,dmin,dmax, fwd,dmax,dmax, &
                                                rev,dmax,dmax, rev,dmin,dmax], &
                                                shape(ecoef) )
  !
continue
  !
  return_value = 0
  !
  np1d = order + 1
  nt1d = np1d + 2
  !
  l = 0
  do ne = 1,12
    !
    ic = ecoef(1,ne)
    jc = ecoef(2,ne)
    kc = ecoef(3,ne)
    !
    jexp = eexp(1,ne)
    kexp = eexp(2,ne)
    iexp = 3 - kexp - jexp
    !
    istride = nt1d**iexp
    jstride = nt1d**jexp
    kstride = nt1d**kexp
    !
    ip = ic * (nt1d-1) * istride
    jp = jc * (nt1d-1) * jstride
    kp = kc * (nt1d-1) * kstride
    !
    do i = 2,np1d+1
      l = l + 1
      idir = (2*ic-1) * (1-i) * istride
      return_value(l) = kp + jp + ip + idir + 1
    end do
    !
  end do
  !
end function get_hexa_edge_indices
!
!###############################################################################
!
pure function get_tria_solpt_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Tria,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_tria_solpt_indices
!
!###############################################################################
!
pure function get_tria_face_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Tria,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_tria_face_indices
!
!###############################################################################
!
pure function get_tria_node_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Tria,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_tria_node_indices
!
!###############################################################################
!
pure function get_quad_solpt_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Quad,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_quad_solpt_indices
!
!###############################################################################
!
pure function get_quad_face_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Quad,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_quad_face_indices
!
!###############################################################################
!
pure function get_quad_node_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Quad,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_quad_node_indices
!
!###############################################################################
!
pure function get_tetr_solpt_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Tetr,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_tetr_solpt_indices
!
!###############################################################################
!
pure function get_tetr_face_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Tetr,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_tetr_face_indices
!
!###############################################################################
!
pure function get_tetr_edge_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Tetr,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_tetr_edge_indices
!
!###############################################################################
!
pure function get_tetr_node_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Tetr,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_tetr_node_indices
!
!###############################################################################
!
pure function get_pyra_solpt_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Pyra,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_pyra_solpt_indices
!
!###############################################################################
!
pure function get_pyra_face_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Pyra,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_pyra_face_indices
!
!###############################################################################
!
pure function get_pyra_edge_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Pyra,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_pyra_edge_indices
!
!###############################################################################
!
pure function get_pyra_node_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Pyra,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_pyra_node_indices
!
!###############################################################################
!
pure function get_pris_solpt_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Pris,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_pris_solpt_indices
!
!###############################################################################
!
pure function get_pris_face_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Pris,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_pris_face_indices
!
!###############################################################################
!
pure function get_pris_edge_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Pris,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_pris_edge_indices
!
!###############################################################################
!
pure function get_pris_node_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Pris,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_pris_node_indices
!
!###############################################################################
!
include "Functions/last.f90"
!
!###############################################################################
!
end module postprocessor_mod
