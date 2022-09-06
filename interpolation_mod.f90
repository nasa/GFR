module interpolation_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use matrix_types_mod, only : matrix
  use matrix_types_mod, only : create_matrix_sparsity
  use matrix_types_mod, only : output_sparsity
  use eqn_idx, only : nec,nmx,nmy,nmz,nee
  !
  implicit none
  !
  private
  !
  !
  !.. Public Procedures ..
  !
  public :: init_interpolation_matrices
  public :: init_outerpolation_matrices
  public :: interpolation_memory_usage
  public :: interp_cell_sol_to_xyz
  public :: move_interp_to_device
  public :: move_outerp_to_device
  public :: outerp_mv_mult
  public :: init_ref_vectors
  !
  !
  !
  ! interpolation_matrix : derived type that stores the matrices for
  !                        interpolating a polynomial quantity within a
  !                        cell to a face point or cell point for a given
  !                        combination of cell geometry and cell order
  !
  type, public :: interpolation_matrix
    type(matrix), allocatable :: toFace(:) ! interpolate to a face point
    type(matrix), allocatable :: toCell(:) ! interpolate to a cell point
    type(matrix), allocatable :: toEdge(:) ! interpolate to a edge point
    type(matrix), allocatable :: toNode    ! interpolate to a node
  end type interpolation_matrix
  !
  ! interp : interpolation matrix for all combinations of cell geometry
  !          and cell order that are possible in the current simulation
  !
  !     For cell 'n' :
  !
  !         The interpolation matrix (giving the contribution of each solution
  !         point within cell 'n' to the solution point 'k' that is currently
  !         of interest) is given by the expression
  !
  !         interp(this_geom,this_order)%toFace(to_order)%mat(:,k)
  !         interp(this_geom,this_order)%toCell(to_order)%mat(:,k)
  !
  !     Where :
  !
  !        this_geom = cell(n)%geom
  !        this_order = cell(n)%order
  !        to_order = order of face or cell to which we are interpolating
  !        k = if using toFace, the face solution point we are interpolating to
  !            if using toCell, the cell solution point we are interpolating to
  !
  type(interpolation_matrix), public, save, allocatable :: interp(:,:)
  !
  ! outerpolation_matrix : derived type that stores the "outerpolation" matrices
  !                        (short for "interpolation for output") which
  !                        interpolate the cell solution to specific quadrature
  !                        points.
  !
  type, public :: outerpolation_matrix
    ! output : Interpolation matrix from solution points to quadrature points
    !          used to output the solution for post-processing
    type(matrix), allocatable :: output
    ! init : Interpolation matrix from solution points to quadrature points
    !        used to compute initial conditions
    type(matrix), allocatable :: init
    ! error : Interpolation matrix from solution points to quadrature points
    !         used for error computations
    type(matrix), allocatable :: error
    ! err_wts : Array containing the quadrature weights
    !           for the error computations
    real(wp), allocatable :: err_wts(:)
  end type outerpolation_matrix
  !
  ! outerp : "outerpolation" matrix for all combinations of cell geometry
  !          and cell order that are possible in the current simulation.
  !          These matrices interpolate the solution from the cell solution
  !          points to a set of nodal points only meant for output and to
  !          visualize the solution.
  !
  type(outerpolation_matrix), public, save, allocatable :: outerp(:,:)
  !
  type :: spanning_vectors_t
    logical(lk) :: found = fals
    integer :: idx1
    integer :: idx2
    integer :: idx3
    real(wp), allocatable :: comp_a(:)
    real(wp), allocatable :: comp_b(:)
    real(wp), allocatable :: comp_c(:)
    real(wp), allocatable :: comp_r(:)
  end type spanning_vectors_t
  !
  type(spanning_vectors_t), save, allocatable :: ref_vectors(:,:)
  !
  interface outerp_mv_mult
    module procedure outerp_mv_mult_with_rv, &
                     outerp_mv_mult_in_place
  end interface outerp_mv_mult
  !
contains
!
!###############################################################################
!
subroutine init_ref_vectors()
  !
  !.. Use Statements ..
  use order_mod,      only : n_min_geom,n_min_order
  use order_mod,      only : n_max_geom,n_max_order
  use quadrature_mod, only : geom_is_used
  !
  !.. Local Scalars ..
  integer :: this_geom,this_order,ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_ref_vectors"
  !
continue
  !
  if (allocated(ref_vectors)) then
    deallocate ( ref_vectors , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ref_vectors",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( ref_vectors(n_min_geom:n_max_geom,n_min_order:n_max_order) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"ref_vectors",1,__LINE__,__FILE__,ierr,error_message)
  !
  order_loop: &
  do this_order = n_min_order,n_max_order
    geom_loop: &
    do this_geom = n_min_geom,n_max_geom
      !
      if (.not. geom_is_used(this_geom)) cycle geom_loop
      if (all(this_geom /= Geom_Valid)) cycle geom_loop
      if (this_geom == Geom_Node) cycle geom_loop
      !
      if (.not. ref_vectors(this_geom,this_order)%found) then
        ref_vectors(this_geom,this_order) = find_ref_vectors(this_geom, &
                                                             this_order)
      end if
      !
    end do geom_loop
  end do order_loop
  !
end subroutine init_ref_vectors
!
!###############################################################################
!
subroutine original_interp_cell_sol_to_xyz(nc,xyz_pt,return_value)
  !
  !.. Use Statements ..
  use order_mod,        only : maxSP
  use order_mod,        only : n_min_geom,n_min_order
  use order_mod,        only : n_max_geom,n_max_order
  use geovar,           only : nr,cell,xyz
  use flowvar,          only : usp
  use element_type_mod, only : std_elem,solpts_edge
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nc
  real(wp),    intent(in) :: xyz_pt(:)
  real(wp), intent(inout) :: return_value(:)
  !
  !.. Local Scalars ..
  integer  :: l,m,n,n1,n2,np,ierr
  integer  :: idx0,idx1,idx2,idx3
  integer  :: this_geom,this_order
  real(wp) :: errmag
  !
  !.. Local Arrays ..
  real(wp) :: xyz_cell(1:nr,1:maxSP)
  real(wp) :: phys_a(1:nr),comp_a(1:nr)
  real(wp) :: phys_b(1:nr),comp_b(1:nr)
  real(wp) :: phys_c(1:nr),comp_c(1:nr)
  real(wp) :: phys_r(1:nr),comp_r(1:nr)
  real(wp) :: del_abc(1:nr)
  real(wp) :: interp_mat(1:maxSP)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "original_interp_cell_sol_to_xyz"
  !
  integer,  parameter :: max_iterations = 50
 !real(wp), parameter :: error_goal = eps3
  real(wp), parameter :: error_goal = eps6
  !
continue
  !
  ! Information for the current cell
  !
  n1 = cell(nc)%beg_sp
  n2 = cell(nc)%end_sp
  np = n2-n1+1
  !
  this_geom = cell(nc)%geom
  this_order = cell(nc)%order
  !
  if (.not. allocated(ref_vectors)) then
    allocate ( ref_vectors(n_min_geom:n_max_geom,n_min_order:n_max_order) , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ref_vectors",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Get the indices for the reference points for this combination of
  ! geometry/solution-order if it hasnt been done yet
  !
  if (.not. ref_vectors(this_geom,this_order)%found) then
    ref_vectors(this_geom,this_order) = find_ref_vectors(this_geom,this_order)
  end if
  !
  ! Get the indices of the points to use in this cell to create the
  ! reference vectors that will used to find the computational space
  ! coordinates corresponding to the physical space coordinates xyz_pt
  !
  idx0 = 1
  idx1 = ref_vectors(this_geom,this_order)%idx1
  idx2 = ref_vectors(this_geom,this_order)%idx2
  idx3 = ref_vectors(this_geom,this_order)%idx3
  !
  ! Create a local copy of the physical space coordinates
  ! of the solution points for this cell
  !
  xyz_cell(1:nr,1:maxSP) = xyz(1:nr,n1:n2)
  !
  ! Get the physical space vector for the first point 'a'
  !
  phys_a(:) = xyz_cell(:,idx1) - xyz_cell(:,idx0)
  !
  ! Get the computational space vector for the first point 'a'
  !
  comp_a(:) = std_elem(this_geom,this_order)%pts(:,idx1) &
            - std_elem(this_geom,this_order)%pts(:,idx0)
  !
  ! Get the physical space vector for the second point 'b'
  !
  phys_b(:) = xyz_cell(:,idx2) - xyz_cell(:,idx0)
  !
  ! Get the computational space vector for the second point 'b'
  !
  comp_b(:) = std_elem(this_geom,this_order)%pts(:,idx2) &
            - std_elem(this_geom,this_order)%pts(:,idx0)
  !
  if (nr == 3) then
    !
    ! If this is 3D, get the physical space vector for the third points 'c'
    !
    phys_c(:) = xyz_cell(:,idx3) - xyz_cell(:,idx0)
    !
    ! If this is 3D, get the computational space vector for the third point 'c'
    !
    comp_c(:) = std_elem(this_geom,this_order)%pts(:,idx3) &
              - std_elem(this_geom,this_order)%pts(:,idx0)
    !
  end if
  !
  ! Compute the physical space vector from the "origin" location
  ! to the location to which we want to interpolate
  !
  phys_r(:) = xyz_pt(:) - xyz_cell(:,idx0)
  !
  ! Initialize the location of xyz_pt in computational space to the
  ! "origin" location used to create the other computational space vectors
  !
  comp_r(:) = std_elem(this_geom,this_order)%pts(:,idx0)
  !
  !
  iter_loop: do n = 1,max_iterations
    !
    del_abc(:) = interpolate_vector_solve(phys_r,phys_a,phys_b,phys_c)
    !
    if (nr == 3) then
      !
      comp_r(1) = comp_r(1) + del_abc(1)*comp_a(1) &
                            + del_abc(2)*comp_b(1) &
                            + del_abc(3)*comp_c(1)
      !
      comp_r(2) = comp_r(2) + del_abc(1)*comp_a(2) &
                            + del_abc(2)*comp_b(2) &
                            + del_abc(3)*comp_c(2)
      !
      comp_r(3) = comp_r(3) + del_abc(1)*comp_a(3) &
                            + del_abc(2)*comp_b(3) &
                            + del_abc(3)*comp_c(3)
      !
    else
      !
      comp_r(1) = comp_r(1) + del_abc(1)*comp_a(1) + del_abc(2)*comp_b(1)
      comp_r(2) = comp_r(2) + del_abc(1)*comp_a(2) + del_abc(2)*comp_b(2)
      !
    end if
    !
    ! Create a new interpolation matrix for the
    ! updated computational space coordinates in comp_r
    !
    select case (this_geom)
      !
      case (Geom_Edge)
        !
        interp_mat(1:np) = InterpToPoint_Edge(solpts_edge,comp_r)
        !
      case (Geom_Tria)
        !
        interp_mat(1:np) = &
                   InterpToPoint_Tria(std_elem(this_geom,this_order)%pts,comp_r)
        !
      case (Geom_Quad)
        !
        interp_mat(1:np) = InterpToPoint_Quad(solpts_edge,comp_r)
        !
      case (Geom_Tetr)
        !
        interp_mat(1:np) = &
                   InterpToPoint_Tetr(std_elem(this_geom,this_order)%pts,comp_r)
        !
      case (Geom_Pyra)
        !
        interp_mat(1:np) = &
                   InterpToPoint_Pyra(std_elem(this_geom,this_order)%pts,comp_r)
        !
      case (Geom_Pris)
        !
        interp_mat(1:np) = &
                   InterpToPoint_Pris(std_elem(this_geom,this_order)%pts,comp_r)
        !
      case (Geom_Hexa)
        !
        interp_mat(1:np) = InterpToPoint_Hexa(solpts_edge,comp_r)
        !
    end select
    !
    ! Find the physical space coordinates corresponding to the
    ! updated computational space coordinates in comp_r using
    ! this new interpolation matrix
    !
    do l = 1,nr
      phys_r(l) = xyz_pt(l) - dot( interp_mat(1:np) , xyz_cell(l,1:np) )
    end do
    !
    ! Check to see if we have converged
    !
    errmag = norm2( phys_r(1:nr) )
    !
    if (errmag < error_goal) then
      exit iter_loop
    end if
    !
  end do iter_loop
  !
  ! Use the last interpolation matrix that was created to interpolate
  ! the solution in the current cell to the input point
  !
  do m = 1,size(usp,dim=2)
    return_value(m) = dot( interp_mat(1:np) , usp(n1:n2,m) )
  end do
  !
end subroutine original_interp_cell_sol_to_xyz
!
!###############################################################################
!
subroutine interp_cell_sol_to_xyz(nc,xyz_pt,return_value)
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nc
  real(wp),    intent(in) :: xyz_pt(:)
  real(wp), intent(inout) :: return_value(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_cell_sol_to_xyz_3d"
  !
continue
  !
  if (nr == 3) then
    call interp_cell_sol_to_xyz_3d(nc,xyz_pt,return_value)
  else if (nr == 2) then
    call interp_cell_sol_to_xyz_2d(nc,xyz_pt,return_value)
  else
    call stop_gfr(abort,pname,__LINE__,__FILE__,"nr /= [23]")
  end if
  !
end subroutine interp_cell_sol_to_xyz
!
!###############################################################################
!
subroutine interp_cell_sol_to_xyz_3d(nc,xyz_pt,return_value)
  !
  !.. Use Statements ..
  use order_mod,        only : maxSP
  use order_mod,        only : n_min_geom,n_min_order
  use order_mod,        only : n_max_geom,n_max_order
  use geovar,           only : cell,xyz
  use flowvar,          only : usp
  use element_type_mod, only : std_elem,solpts_edge
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nc
  real(wp),    intent(in) :: xyz_pt(:)
  real(wp), intent(inout) :: return_value(:)
  !
  !.. Local Scalars ..
  integer  :: l,m,n,niter,n0,n1,n2,np,ierr
  integer  :: idx0,idx1,idx2,idx3
  integer  :: this_geom,this_order
  real(wp) :: errmag,a
  real(wp) :: x0,x1,x2,x3
  real(wp) :: y0,y1,y2,y3
  real(wp) :: z0,z1,z2,z3
  real(wp) :: valx,valy,valz
  real(wp) :: valsp1,valsp2,valsp3,valsp4,valsp5,valsp6
  real(wp) :: valec,valmx,valmy,valmz,valee
  real(wp) :: ax,ay,az,bx,by,bz,cx,cy,cz,rx,ry,rz
  real(wp) :: m11,m12,m13,m21,m22,m23,m31,m32,m33
  real(wp) :: det
  !
  !.. Local Arrays ..
 !real(wp) :: xyz_cell(1:3,1:maxSP)
  real(wp) :: phys_a(1:3),comp_a(1:3)
  real(wp) :: phys_b(1:3),comp_b(1:3)
  real(wp) :: phys_c(1:3),comp_c(1:3)
  real(wp) :: phys_r(1:3),comp_r(1:3)
  real(wp) :: del_abc(1:3)
  real(wp) :: interp_mat(1:maxSP)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_cell_sol_to_xyz_3d"
  !
  integer,  parameter :: max_iterations = 50
 !real(wp), parameter :: error_goal = eps3
  real(wp), parameter :: error_goal = eps6
  !
continue
  !
  ! Information for the current cell
  !
  n0 = cell(nc)%beg_sp - 1
  np = cell(nc)%end_sp - n0
  !
 !n1 = cell(nc)%beg_sp
 !n2 = cell(nc)%end_sp
 !np = n2-n1+1
  !
  this_geom = cell(nc)%geom
  this_order = cell(nc)%order
  !
  ! Get the indices of the points to use in this cell to create the
  ! reference vectors that will used to find the computational space
  ! coordinates corresponding to the physical space coordinates xyz_pt
  !
  idx0 = 1
  idx1 = ref_vectors(this_geom,this_order)%idx1
  idx2 = ref_vectors(this_geom,this_order)%idx2
  idx3 = ref_vectors(this_geom,this_order)%idx3
  !
  ! Create a local copy of the physical space coordinates
  ! of the solution points for this cell
  !
 !xyz_cell(1:nr,1:maxSP) = xyz(1:nr,n1:n2)
  x0 = xyz(1,idx0) ; y0 = xyz(2,idx0) ; z0 = xyz(3,idx0)
  x1 = xyz(1,idx1) ; y1 = xyz(2,idx1) ; z1 = xyz(3,idx1)
  x2 = xyz(1,idx2) ; y2 = xyz(2,idx2) ; z2 = xyz(3,idx2)
  x3 = xyz(1,idx3) ; y3 = xyz(2,idx3) ; z3 = xyz(3,idx3)
  !
  ! Get the physical space vector for the first point 'a'
  !
 !phys_a(:) = xyz_cell(:,idx1) - xyz_cell(:,idx0)
  phys_a(1) = x1 - x0
  phys_a(2) = y1 - y0
  phys_a(3) = z1 - z0
  !
  ! Get the computational space vector for the first point 'a'
  !
 !comp_a(:) = std_elem(this_geom,this_order)%pts(:,idx1) &
 !          - std_elem(this_geom,this_order)%pts(:,idx0)
  comp_a(1) = ref_vectors(this_geom,this_order)%comp_a(1)
  comp_a(2) = ref_vectors(this_geom,this_order)%comp_a(2)
  comp_a(3) = ref_vectors(this_geom,this_order)%comp_a(3)
  !
  ! Get the physical space vector for the second point 'b'
  !
 !phys_b(:) = xyz_cell(:,idx2) - xyz_cell(:,idx0)
  phys_b(1) = x2 - x0
  phys_b(2) = y2 - y0
  phys_b(3) = z2 - z0
  !
  ! Get the computational space vector for the second point 'b'
  !
 !comp_b(:) = std_elem(this_geom,this_order)%pts(:,idx2) &
 !          - std_elem(this_geom,this_order)%pts(:,idx0)
  comp_b(1) = ref_vectors(this_geom,this_order)%comp_b(1)
  comp_b(2) = ref_vectors(this_geom,this_order)%comp_b(2)
  comp_b(3) = ref_vectors(this_geom,this_order)%comp_b(3)
  !
  ! If this is 3D, get the physical space vector for the third points 'c'
  !
 !phys_c(:) = xyz_cell(:,idx3) - xyz_cell(:,idx0)
  phys_c(1) = x3 - x0
  phys_c(2) = y3 - y0
  phys_c(3) = z3 - z0
  !
  ! If this is 3D, get the computational space vector for the third point 'c'
  !
 !comp_c(:) = std_elem(this_geom,this_order)%pts(:,idx3) &
 !          - std_elem(this_geom,this_order)%pts(:,idx0)
  comp_c(1) = ref_vectors(this_geom,this_order)%comp_c(1)
  comp_c(2) = ref_vectors(this_geom,this_order)%comp_c(2)
  comp_c(3) = ref_vectors(this_geom,this_order)%comp_c(3)
  !
  ! Compute the physical space vector from the "origin" location
  ! to the location to which we want to interpolate
  !
 !phys_r(:) = xyz_pt(:) - xyz_cell(:,idx0)
  phys_r(1) = xyz_pt(1) - x0
  phys_r(2) = xyz_pt(2) - y0
  phys_r(3) = xyz_pt(3) - z0
  !
  ! Initialize the location of xyz_pt in computational space to the
  ! "origin" location used to create the other computational space vectors
  !
 !comp_r(:) = std_elem(this_geom,this_order)%pts(:,idx0)
  comp_r(1) = ref_vectors(this_geom,this_order)%comp_r(1)
  comp_r(2) = ref_vectors(this_geom,this_order)%comp_r(2)
  comp_r(3) = ref_vectors(this_geom,this_order)%comp_r(3)
  !
  del_abc(1) = zero
  del_abc(2) = zero
  del_abc(3) = zero
  !
  iter_loop: do niter = 1,max_iterations
    !
   !del_abc(:) = interpolate_vector_solve(phys_r,phys_a,phys_b,phys_c)
    !
    ax = phys_a(1) ; ay = phys_a(2) ; az = phys_a(3)
    bx = phys_b(1) ; by = phys_b(2) ; bz = phys_b(3)
    cx = phys_c(1) ; cy = phys_c(2) ; cz = phys_c(3)
    rx = phys_r(1) ; ry = phys_r(2) ; rz = phys_r(3)
    !
    m11 = by*cz - bz*cy
    m12 = bz*cx - bx*cz
    m13 = bx*cy - by*cx
    !
    m21 = az*cy - ay*cz
    m22 = ax*cz - az*cx
    m23 = ay*cx - ax*cy
    !
    m31 = ay*bz - az*by
    m32 = az*bx - ax*bz
    m33 = ax*by - ay*bx
    !
    det = ax*m11 + ay*m12 + az*m13
    !
    if (abs(det) > eps6) then
      !
      del_abc(1) = rx*m11 + ry*m12 + rz*m13
      del_abc(2) = rx*m21 + ry*m22 + rz*m23
      del_abc(3) = rx*m31 + ry*m32 + rz*m33
      !
      del_abc(1) = del_abc(1) / det
      del_abc(2) = del_abc(2) / det
      del_abc(3) = del_abc(3) / det
      !
    end if
    !
    comp_r(1) = comp_r(1) + del_abc(1)*comp_a(1) &
                          + del_abc(2)*comp_b(1) &
                          + del_abc(3)*comp_c(1)
    !
    comp_r(2) = comp_r(2) + del_abc(1)*comp_a(2) &
                          + del_abc(2)*comp_b(2) &
                          + del_abc(3)*comp_c(2)
    !
    comp_r(3) = comp_r(3) + del_abc(1)*comp_a(3) &
                          + del_abc(2)*comp_b(3) &
                          + del_abc(3)*comp_c(3)
    !
    ! Create a new interpolation matrix for the
    ! updated computational space coordinates in comp_r
    !
    select case (this_geom)
      !
      case (Geom_Edge)
        !
        interp_mat(1:np) = InterpToPoint_Edge(solpts_edge,comp_r)
        !
      case (Geom_Tria)
        !
        interp_mat(1:np) = &
                   InterpToPoint_Tria(std_elem(this_geom,this_order)%pts,comp_r)
        !
      case (Geom_Quad)
        !
        interp_mat(1:np) = InterpToPoint_Quad(solpts_edge,comp_r)
        !
      case (Geom_Tetr)
        !
        interp_mat(1:np) = &
                   InterpToPoint_Tetr(std_elem(this_geom,this_order)%pts,comp_r)
        !
      case (Geom_Pyra)
        !
        interp_mat(1:np) = &
                   InterpToPoint_Pyra(std_elem(this_geom,this_order)%pts,comp_r)
        !
      case (Geom_Pris)
        !
        interp_mat(1:np) = &
                   InterpToPoint_Pris(std_elem(this_geom,this_order)%pts,comp_r)
        !
      case (Geom_Hexa)
        !
        interp_mat(1:np) = InterpToPoint_Hexa(solpts_edge,comp_r)
        !
    end select
    !
    ! Find the physical space coordinates corresponding to the
    ! updated computational space coordinates in comp_r using
    ! this new interpolation matrix
    !
   !do l = 1,nr
   !  phys_r(l) = xyz_pt(l) - dot( interp_mat(1:np) , xyz_cell(l,1:np) )
   !end do
    valx = zero
    valy = zero
    valz = zero
    do n = 1,np
      a = interp_mat(n)
      valx = valx + a * xyz(1,n0+n)
      valy = valy + a * xyz(2,n0+n)
      valz = valz + a * xyz(3,n0+n)
    end do
    phys_r(1) = xyz_pt(1) - valx
    phys_r(2) = xyz_pt(2) - valy
    phys_r(3) = xyz_pt(3) - valz
    !
    ! Check to see if we have converged
    !
   !errmag = norm2( phys_r(1:nr) )
    errmag = sqrt( phys_r(1)*phys_r(1) + &
                   phys_r(2)*phys_r(2) + &
                   phys_r(3)*phys_r(3) )
    !
    if (errmag < error_goal) then
      exit iter_loop
    end if
    !
  end do iter_loop
  !
  ! Use the last interpolation matrix that was created to interpolate
  ! the solution in the current cell to the input point
  !
#define TEMPORARY_FOR_UTRC_CASE
#ifdef TEMPORARY_FOR_UTRC_CASE
  valsp1 = zero
  valsp2 = zero
  valsp3 = zero
  valsp4 = zero
  valsp5 = zero
  valsp6 = zero
  valec  = zero
  valmx  = zero
  valmy  = zero
  valmz  = zero
  valee  = zero
  do n = 1,np
    a = interp_mat(n)
    valsp1 = valsp1 + a * usp(n0+n,1)
    valsp2 = valsp2 + a * usp(n0+n,2)
    valsp3 = valsp3 + a * usp(n0+n,3)
    valsp4 = valsp4 + a * usp(n0+n,4)
    valsp5 = valsp5 + a * usp(n0+n,5)
    valsp6 = valsp6 + a * usp(n0+n,6)
    valec  = valec  + a * usp(n0+n,nec)
    valmx  = valmx  + a * usp(n0+n,nmx)
    valmy  = valmy  + a * usp(n0+n,nmy)
    valmz  = valmz  + a * usp(n0+n,nmz)
    valee  = valee  + a * usp(n0+n,nee)
  end do
  return_value(1)   = valsp1
  return_value(2)   = valsp2
  return_value(3)   = valsp3
  return_value(4)   = valsp4
  return_value(5)   = valsp5
  return_value(6)   = valsp6
  return_value(nec) = valec
  return_value(nmx) = valmx
  return_value(nmy) = valmy
  return_value(nmz) = valmz
  return_value(nee) = valee
#else
 !do m = 1,size(usp,dim=2)
 !  return_value(m) = dot( interp_mat(1:np) , usp(n1:n2,m) )
 !end do
  valec  = zero
  valmx  = zero
  valmy  = zero
  valmz  = zero
  valee  = zero
  do n = 1,np
    a = interp_mat(n)
    valec  = valec  + a * usp(n0+n,nec)
    valmx  = valmx  + a * usp(n0+n,nmx)
    valmy  = valmy  + a * usp(n0+n,nmy)
    valmz  = valmz  + a * usp(n0+n,nmz)
    valee  = valee  + a * usp(n0+n,nee)
  end do
  return_value(nec) = valec
  return_value(nmx) = valmx
  return_value(nmy) = valmy
  return_value(nmz) = valmz
  return_value(nee) = valee
#endif
  !
end subroutine interp_cell_sol_to_xyz_3d
!
!###############################################################################
!
subroutine interp_cell_sol_to_xyz_2d(nc,xyz_pt,return_value)
  !
  !.. Use Statements ..
  use order_mod,        only : maxSP
  use order_mod,        only : n_min_geom,n_min_order
  use order_mod,        only : n_max_geom,n_max_order
  use geovar,           only : cell,xyz
  use flowvar,          only : usp
  use element_type_mod, only : std_elem,solpts_edge
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nc
  real(wp),    intent(in) :: xyz_pt(:)
  real(wp), intent(inout) :: return_value(:)
  !
  !.. Local Scalars ..
  integer  :: l,m,n,niter,n0,n1,n2,np,ierr
  integer  :: idx0,idx1,idx2
  integer  :: this_geom,this_order
  real(wp) :: errmag,a
  real(wp) :: x0,x1,x2
  real(wp) :: y0,y1,y2
  real(wp) :: valx,valy
  real(wp) :: valsp1,valsp2,valsp3,valsp4,valsp5,valsp6
  real(wp) :: valec,valmx,valmy,valee
  real(wp) :: ax,ay,bx,by,cx,cy,rx,ry
  real(wp) :: m11,m12,m21,m22
  real(wp) :: det
  !
  !.. Local Arrays ..
 !real(wp) :: xyz_cell(1:3,1:maxSP)
  real(wp) :: phys_a(1:2),comp_a(1:2)
  real(wp) :: phys_b(1:2),comp_b(1:2)
  real(wp) :: phys_r(1:2),comp_r(1:2)
  real(wp) :: del_abc(1:2)
  real(wp) :: interp_mat(1:maxSP)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_cell_sol_to_xyz_2d"
  !
  integer,  parameter :: max_iterations = 50
 !real(wp), parameter :: error_goal = eps3
  real(wp), parameter :: error_goal = eps6
  !
continue
  !
  ! Information for the current cell
  !
  n0 = cell(nc)%beg_sp - 1
  np = cell(nc)%end_sp - n0
  !
 !n1 = cell(nc)%beg_sp
 !n2 = cell(nc)%end_sp
 !np = n2-n1+1
  !
  this_geom = cell(nc)%geom
  this_order = cell(nc)%order
  !
  if (.not. allocated(ref_vectors)) then
    allocate ( ref_vectors(n_min_geom:n_max_geom,n_min_order:n_max_order) , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ref_vectors",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Get the indices for the reference points for this combination of
  ! geometry/solution-order if it hasnt been done yet
  !
  if (.not. ref_vectors(this_geom,this_order)%found) then
    ref_vectors(this_geom,this_order) = find_ref_vectors(this_geom,this_order)
  end if
  !
  ! Get the indices of the points to use in this cell to create the
  ! reference vectors that will used to find the computational space
  ! coordinates corresponding to the physical space coordinates xyz_pt
  !
  idx0 = 1
  idx1 = ref_vectors(this_geom,this_order)%idx1
  idx2 = ref_vectors(this_geom,this_order)%idx2
  !
  ! Create a local copy of the physical space coordinates
  ! of the solution points for this cell
  !
 !xyz_cell(1:nr,1:maxSP) = xyz(1:nr,n1:n2)
  x0 = xyz(1,idx0) ; y0 = xyz(2,idx0)
  x1 = xyz(1,idx1) ; y1 = xyz(2,idx1)
  x2 = xyz(1,idx2) ; y2 = xyz(2,idx2)
  !
  ! Get the physical space vector for the first point 'a'
  !
 !phys_a(:) = xyz_cell(:,idx1) - xyz_cell(:,idx0)
  phys_a(1) = x1 - x0
  phys_a(2) = y1 - y0
  !
  ! Get the computational space vector for the first point 'a'
  !
 !comp_a(:) = std_elem(this_geom,this_order)%pts(:,idx1) &
 !          - std_elem(this_geom,this_order)%pts(:,idx0)
  comp_a(1) = ref_vectors(this_geom,this_order)%comp_a(1)
  comp_a(2) = ref_vectors(this_geom,this_order)%comp_a(2)
  !
  ! Get the physical space vector for the second point 'b'
  !
 !phys_b(:) = xyz_cell(:,idx2) - xyz_cell(:,idx0)
  phys_b(1) = x2 - x0
  phys_b(2) = y2 - y0
  !
  ! Get the computational space vector for the second point 'b'
  !
 !comp_b(:) = std_elem(this_geom,this_order)%pts(:,idx2) &
 !          - std_elem(this_geom,this_order)%pts(:,idx0)
  comp_b(1) = ref_vectors(this_geom,this_order)%comp_b(1)
  comp_b(2) = ref_vectors(this_geom,this_order)%comp_b(2)
  !
  ! Compute the physical space vector from the "origin" location
  ! to the location to which we want to interpolate
  !
 !phys_r(:) = xyz_pt(:) - xyz_cell(:,idx0)
  phys_r(1) = xyz_pt(1) - x0
  phys_r(2) = xyz_pt(2) - y0
  !
  ! Initialize the location of xyz_pt in computational space to the
  ! "origin" location used to create the other computational space vectors
  !
 !comp_r(:) = std_elem(this_geom,this_order)%pts(:,idx0)
  comp_r(1) = ref_vectors(this_geom,this_order)%comp_r(1)
  comp_r(2) = ref_vectors(this_geom,this_order)%comp_r(2)
  !
  del_abc(1) = zero
  del_abc(2) = zero
  !
  iter_loop: do niter = 1,max_iterations
    !
   !del_abc(:) = interpolate_vector_solve(phys_r,phys_a,phys_b,phys_c)
    !
    ax = phys_a(1) ; ay = phys_a(2)
    bx = phys_b(1) ; by = phys_b(2)
    rx = phys_r(1) ; ry = phys_r(2)
    !
    m11 = -by
    m12 =  bx
    !
    m21 =  ay
    m22 = -ax
    !
    det = ax*m11 + ay*m12
    !
    if (abs(det) > eps6) then
      !
      del_abc(1) = rx*m11 + ry*m12
      del_abc(2) = rx*m21 + ry*m22
      !
      del_abc(1) = del_abc(1) / det
      del_abc(2) = del_abc(2) / det
      !
    end if
    !
    comp_r(1) = comp_r(1) + del_abc(1)*comp_a(1) + del_abc(2)*comp_b(1)
    comp_r(2) = comp_r(2) + del_abc(1)*comp_a(2) + del_abc(2)*comp_b(2)
    !
    ! Create a new interpolation matrix for the
    ! updated computational space coordinates in comp_r
    !
    select case (this_geom)
      !
      case (Geom_Edge)
        !
        interp_mat(1:np) = InterpToPoint_Edge(solpts_edge,comp_r)
        !
      case (Geom_Tria)
        !
        interp_mat(1:np) = &
                   InterpToPoint_Tria(std_elem(this_geom,this_order)%pts,comp_r)
        !
      case (Geom_Quad)
        !
        interp_mat(1:np) = InterpToPoint_Quad(solpts_edge,comp_r)
        !
      case (Geom_Tetr)
        !
        interp_mat(1:np) = &
                   InterpToPoint_Tetr(std_elem(this_geom,this_order)%pts,comp_r)
        !
      case (Geom_Pyra)
        !
        interp_mat(1:np) = &
                   InterpToPoint_Pyra(std_elem(this_geom,this_order)%pts,comp_r)
        !
      case (Geom_Pris)
        !
        interp_mat(1:np) = &
                   InterpToPoint_Pris(std_elem(this_geom,this_order)%pts,comp_r)
        !
      case (Geom_Hexa)
        !
        interp_mat(1:np) = InterpToPoint_Hexa(solpts_edge,comp_r)
        !
    end select
    !
    ! Find the physical space coordinates corresponding to the
    ! updated computational space coordinates in comp_r using
    ! this new interpolation matrix
    !
   !do l = 1,nr
   !  phys_r(l) = xyz_pt(l) - dot( interp_mat(1:np) , xyz_cell(l,1:np) )
   !end do
    valx = zero
    valy = zero
    do n = 1,np
      a = interp_mat(n)
      valx = valx + a * xyz(1,n0+n)
      valy = valy + a * xyz(2,n0+n)
    end do
    phys_r(1) = xyz_pt(1) - valx
    phys_r(2) = xyz_pt(2) - valy
    !
    ! Check to see if we have converged
    !
   !errmag = norm2( phys_r(1:nr) )
    errmag = sqrt( phys_r(1)*phys_r(1) + phys_r(2)*phys_r(2) )
    !
    if (errmag < error_goal) then
      exit iter_loop
    end if
    !
  end do iter_loop
  !
  ! Use the last interpolation matrix that was created to interpolate
  ! the solution in the current cell to the input point
  !
#ifdef TEMPORARY_FOR_UTRC_CASE
  valsp1 = zero
  valsp2 = zero
  valsp3 = zero
  valsp4 = zero
  valsp5 = zero
  valsp6 = zero
  valec  = zero
  valmx  = zero
  valmy  = zero
  valee  = zero
  do n = 1,np
    a = interp_mat(n)
    valsp1 = valsp1 + a * usp(n0+n,1)
    valsp2 = valsp2 + a * usp(n0+n,2)
    valsp3 = valsp3 + a * usp(n0+n,3)
    valsp4 = valsp4 + a * usp(n0+n,4)
    valsp5 = valsp5 + a * usp(n0+n,5)
    valsp6 = valsp6 + a * usp(n0+n,6)
    valec  = valec  + a * usp(n0+n,nec)
    valmx  = valmx  + a * usp(n0+n,nmx)
    valmy  = valmy  + a * usp(n0+n,nmy)
    valee  = valee  + a * usp(n0+n,nee)
  end do
  return_value(1)   = valsp1
  return_value(2)   = valsp2
  return_value(3)   = valsp3
  return_value(4)   = valsp4
  return_value(5)   = valsp5
  return_value(6)   = valsp6
  return_value(nec) = valec
  return_value(nmx) = valmx
  return_value(nmy) = valmy
  return_value(nee) = valee
#else
 !do m = 1,size(usp,dim=2)
 !  return_value(m) = dot( interp_mat(1:np) , usp(n1:n2,m) )
 !end do
  valec  = zero
  valmx  = zero
  valmy  = zero
  valee  = zero
  do n = 1,np
    a = interp_mat(n)
    valec  = valec  + a * usp(n0+n,nec)
    valmx  = valmx  + a * usp(n0+n,nmx)
    valmy  = valmy  + a * usp(n0+n,nmy)
    valee  = valee  + a * usp(n0+n,nee)
  end do
  return_value(nec) = valec
  return_value(nmx) = valmx
  return_value(nmy) = valmy
  return_value(nee) = valee
#endif
  !
end subroutine interp_cell_sol_to_xyz_2d
!
!###############################################################################
!
pure function find_ref_vectors(this_geom,this_order) result(return_value)
  !
  !.. Use Statements ..
  use order_mod,        only : geom_solpts
  use element_type_mod, only : std_elem
  !
  !.. Formal Arguments ..
  integer, intent(in) :: this_geom
  integer, intent(in) :: this_order
  !
  !.. Function Result ..
  type(spanning_vectors_t) :: return_value
  !
  !.. Local Scalars ..
  integer :: np_edge,np_tria,np_quad,ndim
  integer :: idx0,idx1,idx2,idx3
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: comp_a(:)
  real(wp), allocatable :: comp_b(:)
  real(wp), allocatable :: comp_c(:)
  real(wp), allocatable :: comp_r(:)
  !
  !.. Local Parameters ..
 !logical(lk), parameter :: use_short_vectors = fals
  logical(lk), parameter :: use_short_vectors = true
  !
continue
  !
  idx0 = 1
  !
  np_edge = geom_solpts(Geom_Edge,this_order)
  np_tria = geom_solpts(Geom_Tria,this_order)
  np_quad = geom_solpts(Geom_Quad,this_order)
  !
  if (use_short_vectors) then
    !
    select case (this_geom)
      !
      case (Geom_Tria)
        !
        idx1 = 2
        idx2 = np_edge+1
        idx3 = 0
        ndim = 2
        !
      case (Geom_Quad)
        !
        idx1 = 2
        idx2 = np_edge+1
        idx3 = 0
        ndim = 2
        !
      case (Geom_Tetr)
        !
        idx1 = 2
        idx2 = np_edge+1
        idx3 = np_tria+1
        ndim = 3
        !
      case (Geom_Pyra)
        !
        idx1 = 2
        idx2 = np_edge+1
        idx3 = np_quad+1
        ndim = 3
        !
      case (Geom_Pris)
        !
        idx1 = 2
        idx2 = np_edge+1
        idx3 = np_tria+1
        ndim = 3
        !
      case (Geom_Hexa)
        !
        idx1 = 2
        idx2 = np_edge+1
        idx3 = np_quad+1
        ndim = 3
        !
      case default
        !
        idx1 = 2
        idx2 = 0
        idx3 = 0
        ndim = 1
        !
    end select
    !
  else
    !
    select case (this_geom)
      !
      case (Geom_Tria)
        !
        idx1 = np_edge
        idx2 = np_tria
        idx3 = 0
        ndim = 2
        !
      case (Geom_Quad)
        !
        idx1 = np_edge
        idx2 = np_edge*(np_edge-1) + 1
        idx3 = 0
        ndim = 2
        !
      case (Geom_Tetr)
        !
        idx1 = np_edge
        idx2 = np_tria
        idx3 = geom_solpts(Geom_Tetr,this_order)
        ndim = 3
        !
      case (Geom_Pyra)
        !
        idx1 = np_edge
        idx2 = np_edge*(np_edge-1) + 1
        idx3 = geom_solpts(Geom_Pyra,this_order)
        ndim = 3
        !
      case (Geom_Pris)
        !
        idx1 = np_edge
        idx2 = np_tria
        idx3 = np_tria*(np_edge-1) + 1
        ndim = 3
        !
      case (Geom_Hexa)
        !
        idx1 = np_edge
        idx2 = np_edge*(np_edge-1) + 1
        idx3 = np_quad*(np_edge-1) + 1
        ndim = 3
        !
      case default
        !
        idx1 = np_edge
        idx2 = 0
        idx3 = 0
        ndim = 1
        !
    end select
    !
  end if
  !
  comp_a = std_elem(this_geom,this_order)%pts(:,idx1) & ! F2003 AUTO-REALLOC
         - std_elem(this_geom,this_order)%pts(:,idx0)
  !
  if (idx2 > 0) then
    comp_b = std_elem(this_geom,this_order)%pts(:,idx2) & ! F2003 AUTO-REALLOC
           - std_elem(this_geom,this_order)%pts(:,idx0)
  end if
  !
  if (idx3 > 0) then
    comp_c = std_elem(this_geom,this_order)%pts(:,idx3) & ! F2003 AUTO-REALLOC
           - std_elem(this_geom,this_order)%pts(:,idx0)
  end if
  !
  comp_r = std_elem(this_geom,this_order)%pts(:,idx0)
  !
  return_value%idx1 = idx1
  return_value%idx2 = idx2
  return_value%idx3 = idx3
  !
  if (allocated(comp_a)) return_value%comp_a = comp_a ! F2003 AUTO-REALLOC
  if (allocated(comp_b)) return_value%comp_b = comp_b ! F2003 AUTO-REALLOC
  if (allocated(comp_c)) return_value%comp_c = comp_c ! F2003 AUTO-REALLOC
  if (allocated(comp_r)) return_value%comp_r = comp_r ! F2003 AUTO-REALLOC
  !
  return_value%found = true
  !
end function find_ref_vectors
!
!###############################################################################
!
pure function interpolate_vector_solve(r,a,b,c) result(return_value)
  !
  ! This function finds the coefficients that make a linear combination of
  ! some reference vectors equivalent to the vector r.
  !
  ! For 3D: three reference vectors a, b, and c, the solution for these
  !         coefficients becomes a linear equation of the form
  !
  !             coeff_a * ax + coeff_b * bx + coeff_c * cx = rx
  !             coeff_a * ay + coeff_b * by + coeff_c * cy = ry
  !             coeff_a * az + coeff_b * bz + coeff_c * cz = rz
  !
  !             or in matrix form:
  !
  !             | ax  bx  cx |   | coeff_a |   | rx |
  !             | ay  by  cy | X | coeff_b | = | ry |
  !             | az  bz  cz |   | coeff_c |   | rz |
  !
  ! For 2D: two reference vectors a and b, the solution for these
  !         coefficients becomes a linear equation of the form:
  !
  !            coeff_a * ax + coeff_b * bx = rx
  !            coeff_a * ay + coeff_b * by = ry
  !
  !            or in matrix form:
  !
  !            | ax  bx |   | coeff_a |   | rx |
  !            |        | X |         | = |    |
  !            | ay  by |   | coeff_b |   | ry |
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: r(:)
  real(wp), intent(in) :: a(:)
  real(wp), intent(in) :: b(:)
  real(wp), intent(in) :: c(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(r))
  !
  !.. Local Scalars ..
  integer  :: nr
  real(wp) :: ax,ay,az,bx,by,bz,cx,cy,cz,rx,ry,rz
  real(wp) :: m11,m12,m13,m21,m22,m23,m31,m32,m33
  real(wp) :: det
  !
continue
  !
  nr = size(r)
  !
  return_value(:) = zero
  !
  if (nr == 3) then
    !
    ax = a(1) ; ay = a(2) ; az = a(3)
    bx = b(1) ; by = b(2) ; bz = b(3)
    cx = c(1) ; cy = c(2) ; cz = c(3)
    rx = r(1) ; ry = r(2) ; rz = r(3)
    !
    m11 = by*cz - bz*cy
    m12 = bz*cx - bx*cz
    m13 = bx*cy - by*cx
    !
    m21 = az*cy - ay*cz
    m22 = ax*cz - az*cx
    m23 = ay*cx - ax*cy
    !
    m31 = ay*bz - az*by
    m32 = az*bx - ax*bz
    m33 = ax*by - ay*bx
    !
    det = ax*m11 + ay*m12 + az*m13
    !
    if (abs(det) > eps6) then
      !
      return_value(1) = rx*m11 + ry*m12 + rz*m13
      return_value(2) = rx*m21 + ry*m22 + rz*m23
      return_value(3) = rx*m31 + ry*m32 + rz*m33
      !
      return_value(1) = return_value(1) / det
      return_value(2) = return_value(2) / det
      return_value(3) = return_value(3) / det
      !
    end if
    !
  else
    !
    ax = a(1) ; ay = a(2)
    bx = b(1) ; by = b(2)
    rx = r(1) ; ry = r(2)
    !
    m11 = -by
    m12 =  bx
    !
    m21 =  ay
    m22 = -ax
    !
    det = ax*m11 + ay*m12
    !
    if (abs(det) > eps6) then
      !
      return_value(1) = rx*m11 + ry*m12
      return_value(2) = rx*m21 + ry*m22
      !
      return_value(1) = return_value(1) / det
      return_value(2) = return_value(2) / det
      !
    end if
    !
  end if
  !
end function interpolate_vector_solve
!
!###############################################################################
!
pure function InterpToPoint_Edge(xi_sp,xi_pt) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xi_sp(:)
  real(wp), intent(in) :: xi_pt(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(xi_sp))
  !
  !.. Local Scalars ..
  integer  :: i
  real(mp) :: val
  !
  !.. Local Arrays ..
  real(mp) :: xisp(1:size(xi_sp))
  real(mp) :: xipt(1:size(xi_pt))
  !
continue
  !
  xisp(:) = real( xi_sp(:) , kind=kind(xi_sp) )
  xipt(:) = real( xi_pt(:) , kind=kind(xi_pt) )
  !
  do i = 1,size(xisp)
    val = eval_LagrangePoly( [i] , xipt , xisp )
    return_value(i) = real(val,kind=kind(return_value))
  end do
  !
end function InterpToPoint_Edge
!
!###############################################################################
!
pure function InterpToPoint_Quad(xi_sp,xi_pt) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xi_sp(:)
  real(wp), intent(in) :: xi_pt(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(xi_sp)**2)
  !
  !.. Local Scalars ..
  integer  :: i,j,n,nsp
  real(mp) :: val
  !
  !.. Local Arrays ..
  real(mp) :: xisp(1:size(xi_sp))
  real(mp) :: xipt(1:size(xi_pt))
  !
continue
  !
  nsp = size(xi_sp)
  !
  xisp(:) = real( xi_sp(:) , kind=kind(xi_sp) )
  xipt(:) = real( xi_pt(:) , kind=kind(xi_pt) )
  !
  ! Get the coefficients for interpolating a polynomial from the solution
  ! points to the input point. These coefficients are simply the value of
  ! the 2D Lagrange polynomial for a specified solution point evaluated at
  ! the coordinates of the input point.
  !
  ! Get the value of the Lagrange polynomial for the
  ! input point at each of the interior solution points
  !
  n = 0
  do j = 1,nsp
    do i = 1,nsp
      n = n + 1
      val = eval_LagrangePoly( [i,j] , xipt , xisp )
      return_value(n) = real(val,kind=kind(return_value))
    end do
  end do
  !
end function InterpToPoint_Quad
!
!###############################################################################
!
pure function InterpToPoint_Tria(rs_sp,rs_pt) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: rs_sp(:,:)
  real(wp), intent(in) :: rs_pt(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(rs_sp,dim=2))
  !
  !.. Local Scalars ..
  integer :: np,nfp
  !
continue
  !
 !np = size(rs_sp,dim=2)
 !nfp = solpts2order(Geom_Tria,np) + 1
  !
  return_value(:) = zero
  !
  return_value(1) = one
  !
end function InterpToPoint_Tria
!
!###############################################################################
!
pure function InterpToPoint_Tetr(rs_sp,rs_pt) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: rs_sp(:,:)
  real(wp), intent(in) :: rs_pt(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(rs_sp,dim=2))
  !
  !.. Local Scalars ..
  integer :: i,n,nsp3d,nfp3d,nep3d
  !
continue
  !
  nsp3d = size(rs_sp,dim=2)
  !
  return_value(:) = zero
  !
  return_value(1) = one
  !
  ! Find the order for this tetrahedra
  !
  n = -1
  do i = 0,256
    if (Tetr_Solpts(i) == nsp3d) then
      n = i
      exit
    end if
  end do
  !
  ! If n is still -1, return with the interpolation matrix still being zero
  !
  if (n == -1) return
  !
  nfp3d = Cell_Flxpts(Geom_Tetr,n)
  nep3d = n + 1
  !
end function InterpToPoint_Tetr
!
!###############################################################################
!
pure function InterpToPoint_Pyra(rs_sp,rs_pt) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: rs_sp(:,:)
  real(wp), intent(in) :: rs_pt(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(rs_sp,dim=2))
  !
  !.. Local Scalars ..
  integer :: i,n,nfp3d,nsp3d,ierr
  !
continue
  !
  nsp3d = size(rs_sp,dim=2)
  !
  return_value(:) = zero
  !
  return_value(1) = one
  !
  ! Find the order for this tetrahedra
  !
  n = -1
  do i = 0,256
    if (Pyra_Solpts(i) == nsp3d) then
      n = i
      exit
    end if
  end do
  !
  nfp3d = Cell_Flxpts(Geom_Pyra,n)
  !
  ! If n is still -1, return with the interpolation matrix still being zero
  !
  if (n == -1) return
  !
end function InterpToPoint_Pyra
!
!###############################################################################
!
pure function InterpToPoint_Pris(rs_sp,rs_pt) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: rs_sp(:,:)
  real(wp), intent(in) :: rs_pt(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(rs_sp,dim=2))
  !
  !.. Local Scalars ..
  integer :: i,n,nfp3d,nsp3d,ierr
  !
continue
  !
  nsp3d = size(rs_sp,dim=2)
  !
  return_value(:) = zero
  !
  return_value(1) = one
  !
  ! Find the order for this tetrahedra
  !
  n = -1
  do i = 0,256
    if (Pris_Solpts(i) == nsp3d) then
      n = i
      exit
    end if
  end do
  !
  nfp3d = Cell_Flxpts(Geom_Pris,n)
  !
  ! If n is still -1, return with the interpolation matrix still being zero
  !
  if (n == -1) return
  !
end function InterpToPoint_Pris
!
!###############################################################################
!
pure function InterpToPoint_Hexa(xi_sp,xi_pt) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xi_sp(:)
  real(wp), intent(in) :: xi_pt(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(xi_sp)**3)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,n,nsp
  real(mp) :: val
  !
  !.. Local Arrays ..
  real(mp) :: xisp(1:size(xi_sp))
  real(mp) :: xipt(1:size(xi_pt))
  !
continue
  !
  nsp = size(xi_sp)
  !
  xisp(:) = real( xi_sp(:) , kind=kind(xi_sp) )
  xipt(:) = real( xi_pt(:) , kind=kind(xi_pt) )
  !
  ! Get the coefficients for interpolating a polynomial from the solution
  ! points to the input point. These coefficients are simply the value of
  ! the 3D Lagrange polynomial for a specified solution point evaluated at
  ! the coordinates of the input point.
  !
  ! Get the value of the Lagrange polynomial for the
  ! current point at each of the interior solution points
  !
  n = 0
  do k = 1,nsp
    do j = 1,nsp
      do i = 1,nsp
        n = n + 1
        val = eval_LagrangePoly( [i,j,k] , xipt , xisp )
        return_value(n) = real(val,kind=kind(return_value))
      end do
    end do
  end do
  !
end function InterpToPoint_Hexa
!
!###############################################################################
!
pure subroutine outerp_mv_mult_with_rv(outerp_mat,vector,return_value)
  !
  !.. Formal Arguments ..
  type(matrix), intent(in) :: outerp_mat
  real(wp),     intent(in) :: vector(:)
  real(wp),  intent(inout) :: return_value(:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n
  real(wp) :: a,val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "outerp_mv_mult_with_rv"
  !
continue
  !
  do n = 1,size(return_value)
    !
    val = zero
    !
    do i = outerp_mat%col_ptr(n)+1, &
           outerp_mat%col_ptr(n+1)
      !
      k = outerp_mat%col_idx(i)
      a = outerp_mat%col_val(i)
      !
      val = val + a * vector(k)
      !
    end do
    !
    return_value(n) = val
    !
  end do
  !
end subroutine outerp_mv_mult_with_rv
!
!###############################################################################
!
pure subroutine outerp_mv_mult_in_place(outerp_mat,vector)
  !
  !.. Formal Arguments ..
  type(matrix), intent(in) :: outerp_mat
  real(wp),  intent(inout) :: vector(:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n
  real(wp) :: a,val
  !
  !.. Local Arrays ..
  real(wp) :: tmp(1:size(outerp_mat%mat,dim=1))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "outerp_mv_mult_in_place"
  !
continue
  !
  do n = 1,size(outerp_mat%mat,dim=1)
    !
    val = zero
    !
    do i = outerp_mat%col_ptr(n)+1, &
           outerp_mat%col_ptr(n+1)
      !
      k = outerp_mat%col_idx(i)
      a = outerp_mat%col_val(i)
      !
      val = val + a * vector(k)
      !
    end do
    !
    tmp(n) = val
    !
  end do
  !
  vector(:) = tmp(:)
  !
end subroutine outerp_mv_mult_in_place
!
!###############################################################################
!
subroutine init_interpolation_matrices()
  !
  !.. Use Statements ..
  use ovar,             only : loc_quadrature_pts
  use order_mod,        only : geom_solpts,geom_flxpts
  use order_mod,        only : n_min_geom,n_max_geom
  use order_mod,        only : n_min_order,n_max_order
  use order_mod,        only : geom_mio,geom_edgpts
  use quadrature_mod,   only : geom_is_used
  use element_type_mod, only : std_elem,face_elem
  use polynomial_mod,   only : legendre_gauss_quadrature
  use polynomial_mod,   only : legendre_gauss_lobatto_quadrature
  !
  !.. Local Scalars ..
  integer :: gmin,gmax,omin,omax,ierr
  integer :: this_omin,this_omax
  integer :: nsolpts,this_order,this_geom
  integer :: nvolpts,cell_order
  integer :: nflxpts,face_order
  integer :: nedgpts,edge_order
  integer :: nnodpts,node_order
  character(len=200) :: array_name
  !
  !.. Local Allocatable Arrays ..
  real(mp), allocatable :: edgpts(:)
  real(mp), allocatable :: flxpts(:)
  real(mp), allocatable :: pts(:)
  real(mp), allocatable :: solpts(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_interpolation_matrices"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  gmin = n_min_geom
  gmax = n_max_geom
  omin = n_min_order
  omax = n_max_order
  !
  ! Allocate the interp array
  !
  allocate ( interp(gmin:gmax,omin:omax) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"interp",1,__LINE__,__FILE__,ierr,error_message)
  !
  geom_loop: do this_geom = gmin,gmax
    !
    if (.not. geom_is_used(this_geom)) cycle geom_loop
    if (all(this_geom /= Geom_Valid)) cycle geom_loop
    if (this_geom == Geom_Node) cycle geom_loop
    !
    order_loop: do this_order = omin,omax
      !
      ! Get the number of solution points for this cell
      !
      nsolpts = geom_solpts(this_geom,this_order)
      !
      ! Get the minimum and maximum interpolation order for
      ! this combination of cell geometry and cell order
      !
      this_omin = this_order
      this_omax = omax
      !
      ! Allocate the interp%toCell component of interp
      !
      allocate ( interp(this_geom,this_order)%toCell(this_omin:this_omax) , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"toCell"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Loop through all the interpolation orders for this combo
      ! and compute the corresponding interpolation matrices
      !
      cell_order_loop: do cell_order = this_omin,this_omax
        !
        ! Get the number of volume points for this cell order
        !
        nvolpts = geom_solpts(this_geom,cell_order)
        !
        ! Allocate the interpolation matrix for this cell order
        !
        allocate ( interp(this_geom,this_order)% &
                             toCell(cell_order)%mat(1:nvolpts,1:nsolpts) , &
                   source=zero , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"toCell",cell_order
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message)
        !
        pts = init_array(cell_order+1,real(zero,kind=kind(pts))) ! F2003 AUTO-REALLOC
        !
        if (loc_quadrature_pts == Legendre_Gauss) then
          !
          call legendre_gauss_quadrature(pts)
          !
        else if (loc_quadrature_pts == Legendre_Gauss_Lobatto) then
          !
          call legendre_gauss_lobatto_quadrature(pts)
          !
        else if (loc_quadrature_pts == Equi_Distant_Points) then
          !
          pts(:) = real(intseq(-cell_order,cell_order,2),kind=kind(pts)) / &
                   real(max(cell_order,1),kind=kind(pts))
          !
        else
          !
          write (error_message,3) "quadrature", Legendre_Gauss, &
                                  Legendre_Gauss_Lobatto, Equi_Distant_Points, &
                                  "loc_quadrature_pts",loc_quadrature_pts
          call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
          !
        end if
        !
        ! Create the matrix to interpolate to quadrature points
        !
        interp(this_geom,this_order)%toCell(cell_order)%mat(:,:) = &
                     InterpToCell(this_geom,this_order,nsolpts,nvolpts,pts)
        !
        ! Deallocate the pts array if needed
        !
        if (allocated(pts)) then
          deallocate ( pts , stat=ierr , errmsg=error_message )
          call alloc_error(pname,"pts",2,__LINE__,__FILE__,ierr, &
                           error_message,skip_alloc_pause)
        end if
        !
      end do cell_order_loop
      !
      ! Get the minimum and maximum interpolation order for
      ! this combination of cell geometry and cell order
      !
      this_omin = this_order
      this_omax = geom_mio(this_geom,this_order)
      !
      ! this_omax should not be smaller than this_omin. If for some reason
      ! it is, cycle to the next order for this geometry. The array geom_mio
      ! was initialized to -1 so this should also account for geometries that
      ! were never found while creating the array.
      !
      if (this_omax < this_omin) cycle order_loop
      !
      ! Allocate the interp%toFace component of interp
      !
      allocate ( interp(this_geom,this_order)%toFace(this_omin:this_omax) , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"toFace"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Loop through all the interpolation orders for this combo
      ! and compute the corresponding interpolation matrices
      !
      face_order_loop: do face_order = this_omin,this_omax
        !
        ! Get the number of face/flux points for this face order
        !
        nflxpts = geom_flxpts(this_geom,face_order)
        !
        ! Allocate the interpolation matrix for this face order
        !
        allocate ( interp(this_geom,this_order)% &
                             toFace(face_order)%mat(1:nsolpts,1:nflxpts) , &
                   source=zero , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"toFace",face_order
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message)
        !
        ! Compute the interpolation matrices depending on geometry type
        !
        select case (this_geom)
          case (Geom_Edge)
            !
            edgpts = real(std_elem(this_geom,this_order)%pts(1,:), &
                          kind=kind(edgpts))  ! F2003 AUTO-REALLOC
            !
            interp(this_geom,this_order)%toFace(face_order)%mat = &
                    InterpToFace_Edge( edgpts )
            !
          case (Geom_Tria)
            !
            solpts = real(std_elem(this_geom,this_order)%pts, &
                          kind=kind(solpts)) ! F2003 AUTO-REALLOC
            !
            interp(this_geom,this_order)%toFace(face_order)%mat = &
                    InterpToFace_Tria( solpts )
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,10) "tria",this_order,face_order
            call output_sparsity( interp(this_geom,this_order)% &
                                            toFace(face_order) , array_name )
#endif
            !
          case (Geom_Quad)
            !
            edgpts = real(std_elem(Geom_Edge,this_order)%pts(1,:), &
                          kind=kind(edgpts))  ! F2003 AUTO-REALLOC
            flxpts = real(face_elem(Geom_Edge,face_order)%pts(1,:), &
                          kind=kind(flxpts))  ! F2003 AUTO-REALLOC
            !
            interp(this_geom,this_order)%toFace(face_order)%mat = &
                    InterpToFace_Quad( edgpts , flxpts )
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,10) "quad",this_order,face_order
            call output_sparsity( interp(this_geom,this_order)% &
                                            toFace(face_order) , array_name )
#endif
            !
          case (Geom_Tetr)
            !
            solpts = real(std_elem(this_geom,this_order)%pts, &
                          kind=kind(solpts)) ! F2003 AUTO-REALLOC
            !
            interp(this_geom,this_order)%toFace(face_order)%mat = &
                    InterpToFace_Tetr( solpts )
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,10) "tetr",this_order,face_order
            call output_sparsity( interp(this_geom,this_order)% &
                                            toFace(face_order) , array_name )
#endif
            !
          case (Geom_Pyra)
            !
            solpts = real(std_elem(this_geom,this_order)%pts, &
                          kind=kind(solpts)) ! F2003 AUTO-REALLOC
            !
            interp(this_geom,this_order)%toFace(face_order)%mat = &
                    InterpToFace_Pyra( solpts )
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,10) "pyra",this_order,face_order
            call output_sparsity( interp(this_geom,this_order)% &
                                            toFace(face_order) , array_name )
#endif
            !
          case (Geom_Pris)
            !
            solpts = real(std_elem(this_geom,this_order)%pts, &
                          kind=kind(solpts)) ! F2003 AUTO-REALLOC
            !
            interp(this_geom,this_order)%toFace(face_order)%mat = &
                    InterpToFace_Pris( solpts )
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,10) "pris",this_order,face_order
            call output_sparsity( interp(this_geom,this_order)% &
                                            toFace(face_order) , array_name )
#endif
            !
          case (Geom_Hexa)
            !
            edgpts = real(std_elem(Geom_Edge,this_order)%pts(1,:), &
                          kind=kind(edgpts))  ! F2003 AUTO-REALLOC
            flxpts = real(face_elem(Geom_Edge,face_order)%pts(1,:), &
                          kind=kind(flxpts))  ! F2003 AUTO-REALLOC
            !
            interp(this_geom,this_order)%toFace(face_order)%mat = &
                    InterpToFace_Hexa( edgpts , flxpts )
#ifdef SPECIAL_OUTPUT_SPARSITY
            write (array_name,10) "hexa",this_order,face_order
            call output_sparsity( interp(this_geom,this_order)% &
                                            toFace(face_order) , array_name )
#endif
            !
          case default
            !
            write (error_message,2)
            call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
            !
        end select
        !
      end do face_order_loop
      !
      ! Allocate the interp%toEdge component of interp
      !
      allocate ( interp(this_geom,this_order)%toEdge(this_omin:this_omax) , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"toEdge"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Loop through all the interpolation orders for this combo
      ! and compute the corresponding interpolation matrices
      !
      edge_order_loop: do edge_order = this_omin,this_omax
        !
        ! Get the number of edge/flux points for this edge order
        !
        nedgpts = geom_edgpts(this_geom,edge_order)
        !
        ! If nedgpts is zero, skip to the next edge order. This will be
        ! the case for 1D and 2D geometries since the faces are the edges
        !
        if (nedgpts == 0) cycle edge_order_loop
        !
        ! Allocate the interpolation matrix for this edge order
        !
        allocate ( interp(this_geom,this_order)% &
                             toEdge(edge_order)%mat(1:nsolpts,1:nedgpts) , &
                   source=zero , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"toEdge",edge_order
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message)
        !
        ! Compute the interpolation matrices depending on geometry type
        !
        select case (this_geom)
          case (Geom_Tetr)
            !
            solpts = real(std_elem(this_geom,this_order)%pts, &
                          kind=kind(solpts)) ! F2003 AUTO-REALLOC
            !
            interp(this_geom,this_order)%toEdge(edge_order)%mat = &
                    InterpToEdge_Tetr( solpts )
            !
          case (Geom_Pyra)
            !
            solpts = real(std_elem(this_geom,this_order)%pts, &
                          kind=kind(solpts)) ! F2003 AUTO-REALLOC
            !
            interp(this_geom,this_order)%toEdge(edge_order)%mat = &
                    InterpToEdge_Pyra( solpts )
            !
          case (Geom_Pris)
            !
            solpts = real(std_elem(this_geom,this_order)%pts, &
                          kind=kind(solpts)) ! F2003 AUTO-REALLOC
            !
            interp(this_geom,this_order)%toEdge(edge_order)%mat = &
                    InterpToEdge_Pris( solpts )
            !
          case (Geom_Hexa)
            !
            edgpts = real(std_elem(Geom_Edge,this_order)%pts(1,:), &
                          kind=kind(edgpts))  ! F2003 AUTO-REALLOC
            flxpts = real(face_elem(Geom_Edge,edge_order)%pts(1,:), &
                          kind=kind(flxpts))  ! F2003 AUTO-REALLOC
            !
            interp(this_geom,this_order)%toEdge(edge_order)%mat = &
                    InterpToEdge_Hexa( edgpts , flxpts )
            !
          case default
            !
            write (error_message,2)
            call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
            !
        end select
        !
      end do edge_order_loop
      !
      ! There isnt really an order for nodes. The following loop exists
      ! simply so it can be skipped if this is an edge geometry
      !
      node_loop: do node_order = 1,1
        !
        ! Skip edge geometries since the nodes are the faces
        !
        if (this_geom == Geom_Edge) exit node_loop
        !
        ! Allocate the interp%toNode component of interp
        !
        allocate ( interp(this_geom,this_order)%toNode , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"toNode"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message)
        !
        ! Get the number of nodes for this cell
        !
        nnodpts = geom_nodes(this_geom)
        !
        ! Allocate the interpolation matrix for the nodes
        !
        allocate ( interp(this_geom,this_order)% &
                                         toNode%mat(1:nsolpts,1:nnodpts) , &
                   source=zero , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"toNode%mat"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message)
        !
        ! Compute the interpolation matrices depending on geometry type
        !
        select case (this_geom)
          case (Geom_Tria)
            !
            solpts = real(std_elem(this_geom,this_order)%pts, &
                          kind=kind(solpts)) ! F2003 AUTO-REALLOC
            !
            interp(this_geom,this_order)%toNode%mat = InterpToNode_Tria(solpts)
            !
          case (Geom_Quad)
            !
            edgpts = real(std_elem(Geom_Edge,this_order)%pts(1,:), &
                          kind=kind(edgpts))  ! F2003 AUTO-REALLOC
            !
            interp(this_geom,this_order)%toNode%mat = InterpToNode_Quad(edgpts)
            !
          case (Geom_Tetr)
            !
            solpts = real(std_elem(this_geom,this_order)%pts, &
                          kind=kind(solpts)) ! F2003 AUTO-REALLOC
            !
            interp(this_geom,this_order)%toNode%mat = InterpToNode_Tetr(solpts)
            !
          case (Geom_Pyra)
            !
            solpts = real(std_elem(this_geom,this_order)%pts, &
                          kind=kind(solpts)) ! F2003 AUTO-REALLOC
            !
            interp(this_geom,this_order)%toNode%mat = InterpToNode_Pyra(solpts)
            !
          case (Geom_Pris)
            !
            solpts = real(std_elem(this_geom,this_order)%pts, &
                          kind=kind(solpts)) ! F2003 AUTO-REALLOC
            !
            interp(this_geom,this_order)%toNode%mat = InterpToNode_Pris(solpts)
            !
          case (Geom_Hexa)
            !
            edgpts = real(std_elem(Geom_Edge,this_order)%pts(1,:), &
                          kind=kind(edgpts))  ! F2003 AUTO-REALLOC
            !
            interp(this_geom,this_order)%toNode%mat = InterpToNode_Hexa(edgpts)
            !
          case default
            !
            write (error_message,2)
            call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
            !
        end select
        !
      end do node_loop
      !
    end do order_loop
    !
  end do geom_loop
  !
  do this_geom = gmin,gmax
    do this_order = omin,omax
      if (allocated(interp(this_geom,this_order)%toCell)) then
       !call interp(this_geom,this_order)%toCell%check_for_sparsity
        call create_matrix_sparsity( interp(this_geom,this_order)%toCell )
      end if
      if (allocated(interp(this_geom,this_order)%toFace)) then
       !call interp(this_geom,this_order)%toFace%check_for_sparsity
        call create_matrix_sparsity( interp(this_geom,this_order)%toFace )
      end if
      if (allocated(interp(this_geom,this_order)%toEdge)) then
       !call interp(this_geom,this_order)%toEdge%check_for_sparsity
        call create_matrix_sparsity( interp(this_geom,this_order)%toEdge )
      end if
      if (allocated(interp(this_geom,this_order)%toNode)) then
       !call interp(this_geom,this_order)%toNode%check_for_sparsity
        call create_matrix_sparsity( interp(this_geom,this_order)%toNode )
      end if
    end do
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("interp(",a,",",i0,")%",a,:,"(",i0,")%mat")
  2 format (" A grid cell of an unknown geometry type was found!")
  3 format (2x,"Invalid option given for the location of the ",a," points!",/, &
            4x,"Valid options are :",4x,i0," => Legendre Gauss Points",/, &
            27x,i0," => Legendre Gauss Lobatto Points",/, &
            27x,i0," => Equidistant Points",/, &
            9x,"Option given : ",a," = ",i0)
  !
  10 format (a4,".interp.p",i0,".f",i0,".tec.dat")
  !
end subroutine init_interpolation_matrices
!
!###############################################################################
!
subroutine init_outerpolation_matrices()
  !
  !.. Use Statements ..
  use ovar,            only : loc_solution_pts,loc_output_pts
  use ovar,            only : loc_error_pts,loc_init_pts
  use order_mod,       only : n_order,o_order,e_order
  use order_mod,       only : geom_solpts,geom_flxpts
  use order_mod,       only : n_min_geom,n_max_geom
  use order_mod,       only : n_min_order,n_max_order
  use quadrature_mod,  only : geom_is_used
  use polynomial_mod,  only : legendre_gauss_quadrature
  use polynomial_mod,  only : legendre_gauss_lobatto_quadrature
  !
  !.. Local Scalars ..
  integer :: gmin,gmax,omin,omax,ierr
  integer :: this_geom,this_order,face_order
  integer :: this_omin,this_omax
  integer :: npts,opts,epts
  integer :: output_order,error_order
  logical(lk) :: init_matrix_is_needed
  logical(lk) :: output_matrix_is_needed
  logical(lk) :: error_matrix_is_needed
  character(len=200) :: array_name
  !
  !.. Local Allocatable Arrays ..
  real(mp), allocatable :: pts(:)
  real(mp), allocatable :: wts(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_outerpolation_matrices"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  gmin = n_min_geom
  gmax = n_max_geom
  omin = n_min_order
  omax = n_max_order
  !
  ! Allocate the outerp array
  !
  allocate ( outerp(gmin:gmax,omin:omax) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"outerp",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Check if we need to create the interpolation matrices within outerp
  !
  init_matrix_is_needed = fals
  output_matrix_is_needed = fals
  error_matrix_is_needed = fals
  !
  ! Check if we need the matrix to compute the initial conditions
  !
  if (loc_solution_pts /= loc_init_pts) then
    init_matrix_is_needed = true
  end if
  !
  ! Check if we need the matrix to output the solution for post-processing
  !
  if (loc_solution_pts /= loc_output_pts .or. o_order > n_order) then
    output_matrix_is_needed = true
  end if
  !
  ! Check if we need the matrix for computing the error
  !
  if (loc_solution_pts /= loc_error_pts .or. e_order > n_order) then
    error_matrix_is_needed = true
  end if
  !
  ! If none of these matrices are needed, just return from this subroutine
  !
  if (.not.init_matrix_is_needed .and. &
      .not.output_matrix_is_needed .and. &
      .not.error_matrix_is_needed) then
    call debug_timer(leaving_procedure,pname)
    return
  end if
  !
  ! Loop through the geometry/order combinations and create any needed matrices
  !
  geom_loop: do this_geom = gmin,gmax
    !
    if (.not. geom_is_used(this_geom)) cycle geom_loop
    if (all(this_geom /= Geom_Valid)) cycle geom_loop
    if (this_geom == Geom_Node) cycle geom_loop
    !
    order_loop: do this_order = omin,omax
      !
      ! Get the number of solution points for this cell
      !
      npts = geom_solpts(this_geom,this_order)
      !
      ! Create the init component of outerp if needed
      !
      if (init_matrix_is_needed) then
        !
        ! Allocate the init component of outerp
        !
        allocate ( outerp(this_geom,this_order)%init , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"init"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
        !
        ! Allocate the mat component of outerp%init
        !
        allocate ( outerp(this_geom,this_order)%init%mat(1:npts,1:npts) , &
                   source=zero , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"init%mat"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
        !
        ! Get the 1D location of the initialization points
        !
        ! F2003 AUTO-REALLOC
        pts = init_array(this_order+1,real(zero,kind=kind(pts)))
        !
        if (loc_init_pts == Legendre_Gauss) then
          !
          call legendre_gauss_quadrature(pts)
          !
        else if (loc_init_pts == Legendre_Gauss_Lobatto) then
          !
          call legendre_gauss_lobatto_quadrature(pts)
          !
        else
          !
          write (error_message,2) "initialization", &
                                  Legendre_Gauss, Legendre_Gauss_Lobatto, &
                                  "loc_init_pts",loc_init_pts
          call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
          !
        end if
        !
        ! Create the initialization matrix
        !
        outerp(this_geom,this_order)%init%mat(:,:) = &
                        InterpToCell(this_geom,this_order,npts,npts,pts)
        !
      end if
      !
      ! Deallocate the pts array if needed
      !
      if (allocated(pts)) then
        deallocate ( pts , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"pts",2,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
      end if
      !
      ! Create the output component of outerp if needed
      !
      if (output_matrix_is_needed) then
        !
        ! Get the order for output
        !
        output_order = max(this_order,o_order)
        !
        ! Get the number of solution points for this cell
        !
        opts = geom_solpts(this_geom,output_order)
        !
        ! Allocate the output component of outerp
        !
        allocate ( outerp(this_geom,this_order)%output , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"output"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
        !
        ! Allocate the mat component of outerp%output
        !
        allocate ( outerp(this_geom,this_order)%output%mat(1:opts,1:npts) , &
                   source=zero , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"output%mat"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
        !
        ! Get the 1D location of the output points
        !
        ! F2003 AUTO-REALLOC
        pts = init_array(output_order+1,real(zero,kind=kind(pts)))
        !
        if (loc_output_pts == Legendre_Gauss) then
          !
          call legendre_gauss_quadrature(pts)
          !
        else if (loc_output_pts == Legendre_Gauss_Lobatto) then
          !
          call legendre_gauss_lobatto_quadrature(pts)
          !
        else if (loc_output_pts == Equi_Distant_Points) then
          !
          pts(:) = real(intseq(-output_order,output_order,2),kind=kind(pts)) / &
                   real(max(output_order,1),kind=kind(pts))
          !
        else
          !
          write (error_message,3) "output", Legendre_Gauss, &
                                  Legendre_Gauss_Lobatto, Equi_Distant_Points, &
                                  "loc_output_pts",loc_output_pts
          call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
          !
        end if
        !
        ! Create the output matrix
        !
        outerp(this_geom,this_order)%output%mat(:,:) = &
                          InterpToCell(this_geom,this_order,npts,opts,pts)
        !
      end if
      !
      ! Deallocate the pts array if needed
      !
      if (allocated(pts)) then
        deallocate ( pts , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"pts",2,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
      end if
      !
      ! Create the error component of outerp if needed
      !
      if (error_matrix_is_needed) then
        !
        ! Get the order for computing error
        !
        error_order = max(this_order,e_order)
        !
        ! Get the number of solution points for this cell
        !
        epts = geom_solpts(this_geom,error_order)
        !
        ! Allocate the error component of outerp
        !
        allocate ( outerp(this_geom,this_order)%error , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"error"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
        !
        ! Allocate the mat component of outerp%error
        !
        allocate ( outerp(this_geom,this_order)%error%mat(1:epts,1:npts) , &
                   source=zero , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"error%mat"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
        !
        ! Get the 1D location of the error points
        ! and their corresponding quadrature weights
        !
        ! F2003 AUTO-REALLOC
        pts = init_array(error_order+1,real(zero,kind=kind(pts)))
        ! F2003 AUTO-REALLOC
        wts = init_array(error_order+1,real(zero,kind=kind(wts)))
        !
        if (loc_error_pts == Legendre_Gauss) then
          !
          call legendre_gauss_quadrature(pts,wts)
          !
        else if (loc_error_pts == Legendre_Gauss_Lobatto) then
          !
          call legendre_gauss_lobatto_quadrature(pts,wts)
          !
        else
          !
          write (error_message,2) "error", &
                                  Legendre_Gauss, Legendre_Gauss_Lobatto, &
                                  "loc_error_pts",loc_error_pts
          call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
          !
        end if
        !
        ! Create the error matrix
        !
        outerp(this_geom,this_order)%error%mat(:,:) = &
                         InterpToCell(this_geom,this_order,npts,epts,pts)
        !
        ! Allocate the error_wts component of outerp
        !
        allocate ( outerp(this_geom,this_order)%err_wts(1:epts) , &
                   source=zero , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"err_wts"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
        !
        ! Compute the error weights
        !
        outerp(this_geom,this_order)%err_wts(:) = &
                                               InterpWeights(this_geom,epts,wts)
        !
      end if
      !
      ! Deallocate the pts array if needed
      !
      if (allocated(pts)) then
        deallocate ( pts , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"pts",2,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
      end if
      !
      if (allocated(wts)) then
        deallocate ( wts , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"wts",2,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
      end if
      !
    end do order_loop
    !
  end do geom_loop
  !
  do this_geom = gmin,gmax
    do this_order = omin,omax
      if (allocated(outerp(this_geom,this_order)%init)) then
        call create_matrix_sparsity( outerp(this_geom,this_order)%init )
      end if
      if (allocated(outerp(this_geom,this_order)%output)) then
        call create_matrix_sparsity( outerp(this_geom,this_order)%output )
      end if
      if (allocated(outerp(this_geom,this_order)%error)) then
        call create_matrix_sparsity( outerp(this_geom,this_order)%error )
      end if
    end do
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("outerp(",a,",",i0,")%",a)
  2 format (2x,"Invalid option given for the location of the ",a," points!",/, &
            4x,"Valid options are :",4x,i0," => Legendre Gauss Points",/, &
            27x,i0," => Legendre Gauss Lobatto Points",/, &
            9x,"Option given : ",a," = ",i0)
  3 format (2x,"Invalid option given for the location of the ",a," points!",/, &
            4x,"Valid options are :",4x,i0," => Legendre Gauss Points",/, &
            27x,i0," => Legendre Gauss Lobatto Points",/, &
            27x,i0," => Equidistant Points",/, &
            9x,"Option given : ",a," = ",i0)
  !
end subroutine init_outerpolation_matrices
!
!###############################################################################
!
function InterpToCell(this_geom,this_order,opts,npts,pts) result(return_value)
  !
  ! THIS_GEOM  : geometry of the current cell
  ! THIS_ORDER : order of the solution polynomial for the current cell
  ! OPTS       : number of "old" solution points we will be interpolating FROM
  ! NPTS       : number of "new" solution points we will be interpolating TO
  ! PTS        : array of the 1D points we will be interpolating to
  !
  !.. Use Statements ..
  use element_type_mod, only : std_elem
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: this_geom
  integer,  intent(in) :: this_order
  integer,  intent(in) :: opts
  integer,  intent(in) :: npts
  real(mp), intent(in) :: pts(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:npts,1:opts)
  !
  !.. Local Allocatable Arrays ..
  real(mp), allocatable :: edgpts(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "InterptoCell"
  !
continue
  !
  ! Compute the outerpolation matrices depending on geometry type
  !
  select case (this_geom)
    case (Geom_Edge)
      !
      ! F2003 AUTO-REALLOC
      edgpts = real(std_elem(Geom_Edge,this_order)%pts(1,:), &
                    kind=kind(edgpts))
      !
      return_value(:,:) = InterpToCell_Edge( edgpts , pts )
      !
    case (Geom_Tria)
      !
      ! Not sure about triangles here so just use the identity matrix
      !
      return_value(:,:) = identity_matrix( return_value )
      !
    case (Geom_Quad)
      !
      ! F2003 AUTO-REALLOC
      edgpts = real(std_elem(Geom_Edge,this_order)%pts(1,:), &
                    kind=kind(edgpts))
      !
      return_value(:,:) = InterpToCell_Quad( edgpts , pts )
      !
    case (Geom_Tetr)
      !
      ! Not sure about tetrahedra here so just use the identity matrix
      !
      return_value(:,:) = identity_matrix( return_value )
      !
    case (Geom_Pyra)
      !
      ! Not sure about pyramids here so just use the identity matrix
      !
      return_value(:,:) = identity_matrix( return_value )
      !
    case (Geom_Pris)
      !
      ! Not sure about prisms here so just use the identity matrix
      !
      return_value(:,:) = identity_matrix( return_value )
      !
    case (Geom_Hexa)
      !
      ! F2003 AUTO-REALLOC
      edgpts = real(std_elem(Geom_Edge,this_order)%pts(1,:), &
                    kind=kind(edgpts))
      !
      return_value(:,:) = InterpToCell_Hexa( edgpts , pts )
      !
    case default
      !
      write (error_message,1)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      !
  end select
  !
  ! Format Statements
  !
  1 format (" A grid cell of an unknown geometry type was found!")
  !
end function InterpToCell
!
!###############################################################################
!
function InterpWeights(this_geom,npts,wts) result(return_value)
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: this_geom
  integer,  intent(in) :: npts
  real(mp), intent(in) :: wts(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:npts)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "InterpWeights"
  !
continue
  !
  ! Compute the outerpolation matrices depending on geometry type
  !
  select case (this_geom)
    case (Geom_Edge)
      return_value(:) = wts
    case (Geom_Tria)
      return_value(:) = zero
    case (Geom_Quad)
      return_value(:) = InterpWeights_Quad( wts )
    case (Geom_Tetr)
      return_value(:) = zero
    case (Geom_Pyra)
      return_value(:) = zero
    case (Geom_Pris)
      return_value(:) = zero
    case (Geom_Hexa)
      return_value(:) = InterpWeights_Hexa( wts )
    case default
      !
      write (error_message,1)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      !
  end select
  !
  ! Format Statements
  !
  1 format (" A grid cell of an unknown geometry type was found!")
  !
end function InterpWeights
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function General_Interp_Matrix(this_geom,this_order,xyz) result(return_value)
  !
  !.. Use Statements ..
  use vandermonde_mod, only : vand
  use vandermonde_mod, only : General_Basis_Matrix
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: this_geom
  integer,  intent(in) :: this_order
  real(wp), intent(in) :: xyz(:,:)
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: bmat(:,:)
  real(wp), allocatable :: vmat(:,:)
  real(wp), allocatable :: imat(:,:)
  !
continue
  !
  ! Get the basis matrix of the basis functions evaluated at points xyz
  !
  bmat = General_Basis_Matrix(this_geom,this_order,xyz) ! F2003 AUTO-REALLOC
  !
  ! Extract the nodal-to-modal Vandermonde matrix
  !
  vmat = vand(this_geom,this_order)%nodal2modal%mat ! F2003 AUTO-REALLOC
  !
  ! Compute the interpolation matrix
  !
  imat = matmul(transpose(vmat),transpose(bmat)) ! F2003 AUTO-REALLOC
  !
  ! Use the chop function to copy imat into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( real(imat,kind=qp) ) ! F2003 AUTO-REALLOC
  !
end function General_Interp_Matrix
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function generic_interpolation_matrix(basis,newbasis,wmat) &
                                    result(return_value)
  !
  !.. Formal Arguments ..
  real(mp), intent(in) :: basis(:,:)
  real(mp), intent(in) :: newbasis(:,:)
  real(mp), intent(in) :: wmat(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(wmat,1),1:size(newbasis,1))
  !
  !.. Local Arrays ..
  real(mp) :: interp_mat(1:size(wmat,1),1:size(newbasis,1))
  real(mp) :: basisT(1:size(basis,2),1:size(basis,1))
  real(mp) :: pmat(1:size(basis,2),1:size(basis,2))
  real(mp) :: tmat(1:size(basis,2),1:size(wmat,2))
  !
continue
  !
  basisT = matmul( transpose(basis) , wmat )
  basisT = chop( basisT , basisT )
  !
  pmat = matmul( basisT , basis )
  pmat = chop( pmat , pmat )
  !
#ifdef SPECIAL_FOR_INTEL
  pmat = invert_matrix_lapack95( pmat )
#else
  pmat = invert_matrix( pmat )
#endif
  pmat = chop( pmat , pmat )
  !
  tmat = matmul( pmat , basisT )
  tmat = chop( tmat , tmat )
  !
  interp_mat = transpose( matmul( newbasis , tmat ) )
  interp_mat = chop( interp_mat , interp_mat )
  !
  ! Use the chop function to copy interp_mat into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value(:,:) = chop( interp_mat )
  !
end function generic_interpolation_matrix
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function InterpToFace_Edge(xi) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  real(mp), intent(in) :: xi(:)
  !
  !.. Function Return Value ..
  real(wp) :: return_value(1:size(xi),1:2)
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Arrays ..
  real(mp) :: interp_mat(1:size(xi),1:2)
  !
  !.. Local Parameters ..
  real(mp), parameter :: one = 1.0_mp
  !
continue
  !
  do i = 1,size(xi)
    interp_mat(i,1) = eval_LagrangePoly([i],[-one],xi)
    interp_mat(i,2) = eval_LagrangePoly([i],[one],xi)
  end do
  !
  ! Use the chop function to copy interp_mat into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( interp_mat )
  !
end function InterpToFace_Edge
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function InterpToCell_Edge(xi,xi_interp) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  real(mp), intent(in) :: xi(:)
  real(mp), intent(in) :: xi_interp(:)
  !
  !.. Function Return Value ..
  real(wp) :: return_value(1:size(xi),1:size(xi))
  !
  !.. Local Scalars ..
  integer :: i,j,npts
  !
  !.. Local Arrays ..
  real(mp) :: interp_mat(1:size(xi),1:size(xi))
  !
continue
  !
  npts = size(xi)
  !
  ! Evaluate the Lagrange polynomial at each interpolation point
  !
  do j = 1,npts
    do i = 1,npts
      interp_mat(i,j) = eval_LagrangePoly([j],[xi_interp(i)],xi)
    end do
  end do
  !
  ! Use the chop function to copy interp_mat into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( interp_mat )
  !
end function InterpToCell_Edge
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function InterpToNode_Quad(xi_sp) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  real(mp), intent(in) :: xi_sp(:)
  !
  !.. Function Return Value ..
  real(wp) :: return_value(1:size(xi_sp)**2,1:4)
  !
  !.. Local Scalars ..
  integer :: i,j,n,np,nfp,nsp
  !
  !.. Local Arrays ..
  real(mp) :: xy(1:2)
  real(mp) :: nodes(1:2,1:4)
  real(mp) :: interp_mat(1:size(xi_sp)**2,1:4)
  !
continue
  !
  nsp = size(xi_sp)
  !
  ! Fill in the nodes for the reference element
  !
  nodes(1:2,1) = real([-one,-one],kind=kind(nodes))
  nodes(1:2,2) = real([ one,-one],kind=kind(nodes))
  nodes(1:2,3) = real([ one, one],kind=kind(nodes))
  nodes(1:2,4) = real([-one, one],kind=kind(nodes))
  !
  ! For each cell node, get the coefficients for interpolating
  ! a polynomial from the solution points to the node.
  ! These coefficients are simply the value of the 2D Lagrange
  ! polynomial for a specified solution point evaluated at
  ! the coordinates of the cell node.
  !
  ! Loops over the nodes of the cell
  !
  do np = 1,4
    !
    ! Loop to get the value of the Lagrange polynomial for the
    ! current node at each of the interior solution points
    !
    xy(1:2) = nodes(1:2,np)
    !
    n = 0
    do j = 1,nsp
      do i = 1,nsp
        n = n + 1
        interp_mat(n,np) = eval_LagrangePoly([i,j],xy,xi_sp)
      end do
    end do
    !
  end do
  !
  ! Use the chop function to copy interp_mat into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( interp_mat )
  !
end function InterpToNode_Quad
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function InterpToFace_Quad(xi_sp,xi_fp) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  real(mp), intent(in) :: xi_sp(:)
  real(mp), intent(in) :: xi_fp(:)
  !
  !.. Function Return Value ..
  real(wp) :: return_value(1:size(xi_sp)**2,1:4*size(xi_fp))
  !
  !.. Local Scalars ..
  integer :: i,j,n,nf,nfp,nsp
  !
  !.. Local Arrays ..
  real(mp) :: xy(1:2)
  real(mp) :: interp_mat(1:size(xi_sp)**2,1:4*size(xi_fp))
  !
  !.. Local Parameters ..
  real(mp), parameter :: one = 1.0_mp
  !
continue
  !
  nfp = size(xi_fp)
  nsp = size(xi_sp)
  !
  ! For each flux point, get the coefficients for interpolating
  ! a polynomial from the solution points to the flux point.
  ! These coefficients are simply the value of the 2D Lagrange
  ! polynomial for a specified solution point evaluated at
  ! the coordinates of the flux point.
  !
  ! Face 1 - bottom/south edge of quad
  !
  xy(2) = -one
  do n = 1,nfp
    !
    nf = n
    !
    xy(1) = xi_fp(n)
    !
    do j = 1,nsp
      do i = 1,nsp
        interp_mat(nsp*(j-1)+i,nf) = eval_LagrangePoly([i,j],xy,xi_sp)
      end do
    end do
    !
  end do
  !
  ! Face 2 - right/east edge of quad
  !
  xy(1) = one
  do n = 1,nfp
    !
    nf = n + nfp
    !
    xy(2) = xi_fp(n)
    !
    do j = 1,nsp
      do i = 1,nsp
        interp_mat(nsp*(j-1)+i,nf) = eval_LagrangePoly([i,j],xy,xi_sp)
      end do
    end do
    !
  end do
  !
  ! Face 3 - top/north edge of quad
  !
  xy(2) = one
  do n = 1,nfp
    !
    nf = n + 2*nfp
    !
    xy(1) = xi_fp(nfp+1-n)
    !
    do j = 1,nsp
      do i = 1,nsp
        interp_mat(nsp*(j-1)+i,nf) = eval_LagrangePoly([i,j],xy,xi_sp)
      end do
    end do
    !
  end do
  !
  ! Face 4 - left/west edge of quad
  !
  xy(1) = -one
  do n = 1,nfp
    !
    nf = n + 3*nfp
    !
    xy(2) = xi_fp(nfp+1-n)
    !
    do j = 1,nsp
      do i = 1,nsp
        interp_mat(nsp*(j-1)+i,nf) = eval_LagrangePoly([i,j],xy,xi_sp)
      end do
    end do
    !
  end do
  !
  ! Use the chop function to copy interp_mat into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( interp_mat )
  !
end function InterpToFace_Quad
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function InterpToCell_Quad(xi_old,xi_new) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  real(mp), intent(in) :: xi_old(:)
  real(mp), intent(in) :: xi_new(:)
  !
  !.. Function Return Value ..
  real(wp) :: return_value(1:size(xi_new)**2,1:size(xi_old)**2)
  !
  !.. Local Scalars ..
  integer :: o,oi,oj,opts
  integer :: n,ni,nj,npts
  !
  !.. Local Arrays ..
  real(mp) :: xy(1:2)
  real(mp) :: interp_mat(1:size(xi_new)**2,1:size(xi_old)**2)
  !
continue
  !
  opts = size(xi_old)
  npts = size(xi_new)
  !
  ! Evaluate the Lagrange polynomial at each interpolation point
  !
  do oj = 1,opts
    do oi = 1,opts
      !
      o = (oj-1)*opts + oi
      !
      do nj = 1,npts
        !
        xy(2) = xi_new(nj)
        !
        do ni = 1,npts
          !
          xy(1) = xi_new(ni)
          !
          n = (nj-1)*npts + ni
          !
          interp_mat(n,o) = eval_LagrangePoly([oi,oj],xy,xi_old)
          !
        end do
      end do
      !
    end do
  end do
  !
  ! Use the chop function to copy interp_mat into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( interp_mat )
  !
end function InterpToCell_Quad
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function InterpWeights_Quad(wts) result(return_value)
  !
  !.. Formal Arguments ..
  real(mp), intent(in) :: wts(:)
  !
  !.. Function Return Value ..
  real(wp) :: return_value(1:size(wts)**2)
  !
  !.. Local Scalars ..
  integer :: m,mi,mj,npts
  !
  !.. Local Arrays ..
  real(mp) :: weights(1:size(wts)**2)
  !
continue
  !
  npts = size(wts)
  !
  ! Evaluate the tensor product of the 1D weights
  !
  do mj = 1,npts
    do mi = 1,npts
      !
      m = (mj-1)*npts + mi
      !
      weights(m) = wts(mi)*wts(mj)
      !
    end do
  end do
  !
  ! Use the chop function to copy weights into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( weights )
  !
end function InterpWeights_Quad
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function InterpToNode_Tria(rs) result(return_value)
  !
  !.. Formal Arguments ..
  real(mp), intent(in) :: rs(:,:)
  !
  !.. Function Return Argument ..
  real(wp) :: return_value(1:size(rs,2),1:3)
  !
  !.. Local Scalars ..
  integer :: np,nfp
  !
continue
  !
  np = size(rs,2)
  nfp = solpts2order(Geom_Tria,np) + 1
  !
  return_value(:,:) = zero
  !
  ! Node 1
  !
  return_value(1,1) = one
  !
  ! Node 2
  !
  return_value(nfp,2) = one
  !
  ! Node 3
  !
  return_value(np,3) = one
  !
end function InterpToNode_Tria
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function InterpToFace_Tria(rs) result(return_value)
  !
  !.. Formal Arguments ..
  real(mp), intent(in) :: rs(:,:)
  !
  !.. Function Return Argument ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,np,nfp
  !
  !.. Local Arrays ..
  integer, allocatable :: Fmask(:,:)
  !
continue
  !
  np = size(rs,2)
  nfp = solpts2order(Geom_Tria,np) + 1
  !
  return_value = init_array(np,3*nfp,zero) ! F2003 AUTO-REALLOC
  !
  Fmask = init_array(nfp,3,0) ! F2003 AUTO-REALLOC
  !
  ! Face 1 - bottom edge of triangle
  !
  do i = 1,nfp
    Fmask(i,1) = i
  end do
  !
  ! Face 2 - diagonal edge of triangle
  !
  j = 0
  do i = 1,nfp
    j = j + nfp + 1 - i
    Fmask(i,2) = j
  end do
  !
  ! Face 3 - vertical edge of triangle
  !
  j = np + 1
  do i = 1,nfp
    j = j - i
    Fmask(i,3) = j
  end do
  !
  do j = 1,3
    do i = 1,nfp
      return_value(Fmask(i,j),nfp*(j-1)+i) = one
    end do
  end do
  !
end function InterpToFace_Tria
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function InterpToNode_Tetr(rst) result(return_value)
  !
  !.. Formal Arguments ..
  real(mp), intent(in) :: rst(:,:)
  !
  !.. Function Return Argument ..
  real(wp) :: return_value(1:size(rst,2),1:4)
  !
  !.. Local Scalars ..
  integer :: i,n,nsp3d,nfp3d,nep3d
  !
continue
  !
  nsp3d = size(rst,2)
  !
  return_value(:,:) = zero
  !
  ! Find the order for this tetrahedra
  !
  n = -1
  do i = 0,256
    if (Tetr_Solpts(i) == nsp3d) then
      n = i
      exit
    end if
  end do
  !
  ! If n is still -1, return with the interpolation matrix still being zero
  !
  if (n == -1) return
  !
  nfp3d = Cell_Flxpts(Geom_Tetr,n)
  nep3d = n + 1
  !
  ! Node 1
  !
  return_value(1,1) = one
  !
  ! Node 2
  !
  return_value(nep3d,2) = one
  !
  ! Node 3
  !
  return_value(nfp3d,3) = one
  !
  ! Node 4
  !
  return_value(nsp3d,4) = one
  !
end function InterpToNode_Tetr
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function InterpToEdge_Tetr(rst) result(return_value)
  !
  !.. Formal Arguments ..
  real(mp), intent(in) :: rst(:,:)
  !
  !.. Function Return Argument ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,n,nep3d,nsp3d,ierr
  !
continue
  !
  nsp3d = size(rst,2)
  !
  ! Find the order for this tetrahedra
  !
  n = -1
  do i = 0,256
    if (Tetr_Solpts(i) == nsp3d) then
      n = i
      exit
    end if
  end do
  !
  nep3d = geom_edges(Geom_Tetr) * (n + 1)
  !
  ! Allocate the return array containing the
  ! working precision interpolation matrix
  !
  allocate ( return_value(1:nsp3d,1:nep3d) , source=zero , stat=ierr )
  !
  ! If n is still -1, return with the interpolation matrix still being zero
  !
  if (n == -1) return
  !
end function InterpToEdge_Tetr
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function InterpToFace_Tetr(rst) result(return_value)
  !
  !.. Formal Arguments ..
  real(mp), intent(in) :: rst(:,:)
  !
  !.. Function Return Argument ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,n,nfp3d,nsp3d,ierr
  !
continue
  !
  nsp3d = size(rst,2)
  !
  ! Find the order for this tetrahedra
  !
  n = -1
  do i = 0,256
    if (Tetr_Solpts(i) == nsp3d) then
      n = i
      exit
    end if
  end do
  !
  nfp3d = Cell_Flxpts(Geom_Tetr,n)
  !
  ! Allocate the return array containing the
  ! working precision interpolation matrix
  !
  allocate ( return_value(1:nsp3d,1:nfp3d) , source=zero , stat=ierr )
  !
  ! If n is still -1, return with the interpolation matrix still being zero
  !
  if (n == -1) return
  !
end function InterpToFace_Tetr
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function InterpToNode_Pyra(rst) result(return_value)
  !
  !.. Formal Arguments ..
  real(mp), intent(in) :: rst(:,:)
  !
  !.. Function Return Argument ..
  real(wp) :: return_value(1:size(rst,2),1:5)
  !
  !.. Local Scalars ..
  integer :: i,n,nfp3d,nsp3d,ierr
  !
continue
  !
  nsp3d = size(rst,2)
  !
  return_value(:,:) = zero
  !
  ! Find the order for this tetrahedra
  !
  n = -1
  do i = 0,256
    if (Pyra_Solpts(i) == nsp3d) then
      n = i
      exit
    end if
  end do
  !
  nfp3d = Cell_Flxpts(Geom_Pyra,n)
  !
  ! If n is still -1, return with the interpolation matrix still being zero
  !
  if (n == -1) return
  !
end function InterpToNode_Pyra
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function InterpToEdge_Pyra(rst) result(return_value)
  !
  !.. Formal Arguments ..
  real(mp), intent(in) :: rst(:,:)
  !
  !.. Function Return Argument ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,n,nep3d,nsp3d,ierr
  !
continue
  !
  nsp3d = size(rst,2)
  !
  ! Find the order for this tetrahedra
  !
  n = -1
  do i = 0,256
    if (Pyra_Solpts(i) == nsp3d) then
      n = i
      exit
    end if
  end do
  !
  nep3d = geom_edges(Geom_Pyra) * (n + 1)
  !
  ! Allocate the return array containing the
  ! working precision interpolation matrix
  !
  allocate ( return_value(1:nsp3d,1:nep3d) , source=zero , stat=ierr )
  !
  ! If n is still -1, return with the interpolation matrix still being zero
  !
  if (n == -1) return
  !
end function InterpToEdge_Pyra
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function InterpToFace_Pyra(rst) result(return_value)
  !
  !.. Formal Arguments ..
  real(mp), intent(in) :: rst(:,:)
  !
  !.. Function Return Argument ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,n,nfp3d,nsp3d,ierr
  !
continue
  !
  nsp3d = size(rst,2)
  !
  ! Find the order for this tetrahedra
  !
  n = -1
  do i = 0,256
    if (Pyra_Solpts(i) == nsp3d) then
      n = i
      exit
    end if
  end do
  !
  nfp3d = Cell_Flxpts(Geom_Pyra,n)
  !
  ! Allocate the return array containing the
  ! working precision interpolation matrix
  !
  allocate ( return_value(1:nsp3d,1:nfp3d) , source=zero , stat=ierr )
  !
  ! If n is still -1, return with the interpolation matrix still being zero
  !
  if (n == -1) return
  !
end function InterpToFace_Pyra
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function InterpToNode_Pris(rst) result(return_value)
  !
  !.. Formal Arguments ..
  real(mp), intent(in) :: rst(:,:)
  !
  !.. Function Return Argument ..
  real(wp) :: return_value(1:size(rst,2),1:6)
  !
  !.. Local Scalars ..
  integer :: i,n,nfp3d,nsp3d,ierr
  !
continue
  !
  nsp3d = size(rst,2)
  !
  return_value(:,:) = zero
  !
  ! Find the order for this tetrahedra
  !
  n = -1
  do i = 0,256
    if (Pris_Solpts(i) == nsp3d) then
      n = i
      exit
    end if
  end do
  !
  nfp3d = Cell_Flxpts(Geom_Pris,n)
  !
  ! If n is still -1, return with the interpolation matrix still being zero
  !
  if (n == -1) return
  !
end function InterpToNode_Pris
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function InterpToEdge_Pris(rst) result(return_value)
  !
  !.. Formal Arguments ..
  real(mp), intent(in) :: rst(:,:)
  !
  !.. Function Return Argument ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,n,nep3d,nsp3d,ierr
  !
continue
  !
  nsp3d = size(rst,2)
  !
  ! Find the order for this tetrahedra
  !
  n = -1
  do i = 0,256
    if (Pris_Solpts(i) == nsp3d) then
      n = i
      exit
    end if
  end do
  !
  nep3d = geom_edges(Geom_Pris) * (n + 1)
  !
  ! Allocate the return array containing the
  ! working precision interpolation matrix
  !
  allocate ( return_value(1:nsp3d,1:nep3d) , source=zero , stat=ierr )
  !
  ! If n is still -1, return with the interpolation matrix still being zero
  !
  if (n == -1) return
  !
end function InterpToEdge_Pris
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function InterpToFace_Pris(rst) result(return_value)
  !
  !.. Formal Arguments ..
  real(mp), intent(in) :: rst(:,:)
  !
  !.. Function Return Argument ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,n,nfp3d,nsp3d,ierr
  !
continue
  !
  nsp3d = size(rst,2)
  !
  ! Find the order for this tetrahedra
  !
  n = -1
  do i = 0,256
    if (Pris_Solpts(i) == nsp3d) then
      n = i
      exit
    end if
  end do
  !
  nfp3d = Cell_Flxpts(Geom_Pris,n)
  !
  ! Allocate the return array containing the
  ! working precision interpolation matrix
  !
  allocate ( return_value(1:nsp3d,1:nfp3d) , source=zero , stat=ierr )
  !
  ! If n is still -1, return with the interpolation matrix still being zero
  !
  if (n == -1) return
  !
end function InterpToFace_Pris
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function InterpToNode_Hexa(xi_sp) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  real(mp), intent(in) :: xi_sp(:)
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,k,n,np,ierr
  integer :: nsp1d,nsp3d
  !
  !.. Local Arrays ..
  real(mp) :: xyz(1:3)
  real(mp) :: nodes(1:3,1:8)
  !
  !.. Local Allocatable Arrays ..
  real(mp), allocatable :: interp_mat(:,:)
  !
continue
  !
  nsp1d = size(xi_sp)
  nsp3d = nsp1d**3
  !
  ! Fill in the nodes for the reference element
  !
  nodes(1:3,1) = real([-one,-one,-one],kind=kind(nodes))
  nodes(1:3,2) = real([ one,-one,-one],kind=kind(nodes))
  nodes(1:3,3) = real([ one, one,-one],kind=kind(nodes))
  nodes(1:3,4) = real([-one, one,-one],kind=kind(nodes))
  nodes(1:3,5) = real([-one,-one, one],kind=kind(nodes))
  nodes(1:3,6) = real([ one,-one, one],kind=kind(nodes))
  nodes(1:3,7) = real([ one, one, one],kind=kind(nodes))
  nodes(1:3,8) = real([-one, one, one],kind=kind(nodes))
  !
  ! Create the interpolation matrix using quad precision
  !
  ! F2003 AUTO-REALLOC
  interp_mat = init_array(nsp3d,8,real(zero,kind=kind(interp_mat)))
  !
  ! For each cell node, get the coefficients for interpolating
  ! a polynomial from the solution points to the node.
  ! These coefficients are simply the value of the 3D Lagrange
  ! polynomial for a specified solution point evaluated at
  ! the coordinates of the cell node.
  !
  ! Loops over the nodes of the cell
  !
  do np = 1,8
    !
    ! Loop to get the value of the Lagrange polynomial for the
    ! current node at each of the interior solution points
    !
    xyz(1:3) = nodes(1:3,np)
    !
    n = 0
    do k = 1,nsp1d
      do j = 1,nsp1d
        do i = 1,nsp1d
          n = n + 1
          interp_mat(n,np) = eval_LagrangePoly([i,j,k],xyz,xi_sp)
        end do
      end do
    end do
    !
  end do
  !
  ! Use the chop function to copy interp_mat into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( interp_mat ) ! F2003 AUTO-REALLOC
  !
end function InterpToNode_Hexa
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function InterpToEdge_Hexa(xi_sp,xi_ep) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  real(mp), intent(in) :: xi_sp(:)
  real(mp), intent(in) :: xi_ep(:)
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,k,ifa,jfa,kfa,n,ne,ierr
  integer :: nep1d,nep3d
  integer :: nsp1d,nsp3d
  integer :: this_edge,n1,n2,l,lbeg,lend,lint,lp
  !
  !.. Local Arrays ..
  real(mp) :: xyz(1:3)
  !
  !.. Local Allocatable Arrays ..
  real(mp), allocatable :: interp_mat(:,:)
  !
  !.. Local Parameters ..
  real(mp), parameter :: one = 1.0_mp
  integer, parameter :: unit_cell(1:3,1:8) = reshape( [ 0,0,0, 1,0,0, &
                                                        1,1,0, 0,1,0, &
                                                        0,0,1, 1,0,1, &
                                                        1,1,1, 0,1,1], &
                                                      shape(unit_cell) )
  logical(lk), parameter :: use_condensed_loop = true
 !logical(lk), parameter :: use_condensed_loop = fals
  !
continue
  !
  nep1d = size(xi_ep)
  nsp1d = size(xi_sp)
  !
  nep3d = 12*nep1d
  nsp3d = nsp1d**3
  !
  ! F2003 AUTO-REALLOC
  interp_mat = init_array(nsp3d,nep3d,real(zero,kind=kind(interp_mat)))
  !
  if (use_condensed_loop) then
    !
    ne = 0
    !
    do this_edge = 1,geom_edges(Geom_Hexa)
      !
      n1 = edge_nodes(1,this_edge,Geom_Hexa) + 1
      n2 = edge_nodes(2,this_edge,Geom_Hexa) + 1
      !
      xyz(1) = merge( -one , one , unit_cell(1,n1) == 0 )
      xyz(2) = merge( -one , one , unit_cell(2,n1) == 0 )
      xyz(3) = merge( -one , one , unit_cell(3,n1) == 0 )
      !
      if (unit_cell(1,n1) /= unit_cell(1,n2)) then
        l = 1
      else if (unit_cell(2,n1) /= unit_cell(2,n2)) then
        l = 2
      else if (unit_cell(3,n1) /= unit_cell(3,n2)) then
        l = 3
      end if
      !
      lbeg = max( 1 , nep1d*unit_cell(l,n1) )
      lend = max( 1 , nep1d*unit_cell(l,n2) )
      lint = merge( 1 , -1 , unit_cell(l,n1) < unit_cell(l,n2) )
      !
      do lp = lbeg,lend,lint
        !
        ne = ne + 1 ! Proceed to next point on this edge
        !
        xyz(l) = xi_ep(lp)
        !
        n = 0
        do k = 1,nsp1d
          do j = 1,nsp1d
            do i = 1,nsp1d
              n = n + 1
              interp_mat(n,ne) = eval_LagrangePoly([i,j,k],xyz,xi_sp)
            end do
          end do
        end do
        !
      end do
      !
    end do
    !
  else
    !
    ! For each edge point, get the coefficients for interpolating
    ! a polynomial from the solution points to the edge point.
    ! These coefficients are simply the value of the 3D Lagrange
    ! polynomial for a specified solution point evaluated at
    ! the coordinates of the edge point.
    !
    ! #####################################################
    ! #####################################################
    ! #####################################################
    ! #####  Edge  1 - x/xi  min and z/zeta min edge  #####
    ! #####################################################
    ! #####################################################
    ! #####################################################
    !
    ne = 0*nep1d+1 ! First point on this edge
    !
    xyz(1) = -one ! x/xi   coordinate of current edge point
    xyz(3) = -one ! z/zeta coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do jfa = 1,nep1d
      !
      xyz(2) = xi_ep(jfa) ! y/eta coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat(n,ne) = eval_LagrangePoly([i,j,k],xyz,xi_sp)
          end do
        end do
      end do
      !
      ne = ne + 1 ! Proceed to next point on this face
      !
    end do
    !
    ! #####################################################
    ! #####################################################
    ! #####################################################
    ! #####  Edge  2 - y/eta max and z/zeta min edge  #####
    ! #####################################################
    ! #####################################################
    ! #####################################################
    !
    ne = 1*nep1d+1 ! First point on this edge
    !
    xyz(2) = +one ! y/eta  coordinate of current edge point
    xyz(3) = -one ! z/zeta coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do ifa = 1,nep1d
      !
      xyz(1) = xi_ep(ifa) ! x/xi  coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat(n,ne) = eval_LagrangePoly([i,j,k],xyz,xi_sp)
          end do
        end do
      end do
      !
      ne = ne + 1 ! Proceed to next point on this face
      !
    end do
    !
    ! #####################################################
    ! #####################################################
    ! #####################################################
    ! #####  Edge  3 - x/xi  max and z/zeta min edge  #####
    ! #####################################################
    ! #####################################################
    ! #####################################################
    !
    ne = 2*nep1d+1 ! First point on this edge
    !
    xyz(1) = +one ! x/xi   coordinate of current edge point
    xyz(3) = -one ! z/zeta coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do jfa = 1,nep1d
      !
      xyz(2) = xi_ep(jfa) ! y/eta coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat(n,ne) = eval_LagrangePoly([i,j,k],xyz,xi_sp)
          end do
        end do
      end do
      !
      ne = ne + 1 ! Proceed to next point on this face
      !
    end do
    !
    ! #####################################################
    ! #####################################################
    ! #####################################################
    ! #####  Edge  4 - y/eta min and z/zeta min edge  #####
    ! #####################################################
    ! #####################################################
    ! #####################################################
    !
    ne = 3*nep1d+1 ! First point on this edge
    !
    xyz(2) = -one ! y/eta  coordinate of current edge point
    xyz(3) = -one ! z/zeta coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do ifa = 1,nep1d
      !
      xyz(1) = xi_ep(ifa) ! x/xi  coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat(n,ne) = eval_LagrangePoly([i,j,k],xyz,xi_sp)
          end do
        end do
      end do
      !
      ne = ne + 1 ! Proceed to next point on this face
      !
    end do
    !
    ! #####################################################
    ! #####################################################
    ! #####################################################
    ! #####  Edge  5 - y/eta min and z/zeta max edge  #####
    ! #####################################################
    ! #####################################################
    ! #####################################################
    !
    ne = 4*nep1d+1 ! First point on this edge
    !
    xyz(2) = -one ! y/eta  coordinate of current edge point
    xyz(3) = +one ! z/zeta coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do ifa = 1,nep1d
      !
      xyz(1) = xi_ep(ifa) ! x/xi  coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat(n,ne) = eval_LagrangePoly([i,j,k],xyz,xi_sp)
          end do
        end do
      end do
      !
      ne = ne + 1 ! Proceed to next point on this face
      !
    end do
    !
    ! #####################################################
    ! #####################################################
    ! #####################################################
    ! #####  Edge  6 - x/xi  max and z/zeta max edge  #####
    ! #####################################################
    ! #####################################################
    ! #####################################################
    !
    ne = 5*nep1d+1 ! First point on this edge
    !
    xyz(1) = +one ! x/xi   coordinate of current edge point
    xyz(3) = +one ! z/zeta coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do jfa = 1,nep1d
      !
      xyz(2) = xi_ep(jfa) ! y/eta coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat(n,ne) = eval_LagrangePoly([i,j,k],xyz,xi_sp)
          end do
        end do
      end do
      !
      ne = ne + 1 ! Proceed to next point on this face
      !
    end do
    !
    ! #####################################################
    ! #####################################################
    ! #####################################################
    ! #####  Edge  7 - y/eta max and z/zeta max edge  #####
    ! #####################################################
    ! #####################################################
    ! #####################################################
    !
    ne = 6*nep1d+1 ! First point on this edge
    !
    xyz(2) = +one ! y/eta  coordinate of current edge point
    xyz(3) = +one ! z/zeta coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do ifa = 1,nep1d
      !
      xyz(1) = xi_ep(ifa) ! x/xi  coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat(n,ne) = eval_LagrangePoly([i,j,k],xyz,xi_sp)
          end do
        end do
      end do
      !
      ne = ne + 1 ! Proceed to next point on this face
      !
    end do
    !
    ! #####################################################
    ! #####################################################
    ! #####################################################
    ! #####  Edge  8 - x/xi  min and z/zeta max edge  #####
    ! #####################################################
    ! #####################################################
    ! #####################################################
    !
    ne = 7*nep1d+1 ! First point on this edge
    !
    xyz(1) = -one ! x/xi   coordinate of current edge point
    xyz(3) = +one ! z/zeta coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do jfa = 1,nep1d
      !
      xyz(2) = xi_ep(jfa) ! y/eta coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat(n,ne) = eval_LagrangePoly([i,j,k],xyz,xi_sp)
          end do
        end do
      end do
      !
      ne = ne + 1 ! Proceed to next point on this face
      !
    end do
    !
    ! #####################################################
    ! #####################################################
    ! #####################################################
    ! #####  Edge  9 - x/xi  max and y/eta  max edge  #####
    ! #####################################################
    ! #####################################################
    ! #####################################################
    !
    ne = 8*nep1d+1 ! First point on this edge
    !
    xyz(1) = +one ! x/xi   coordinate of current edge point
    xyz(2) = +one ! y/eta  coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do kfa = 1,nep1d
      !
      xyz(3) = xi_ep(kfa) ! z/zeta coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat(n,ne) = eval_LagrangePoly([i,j,k],xyz,xi_sp)
          end do
        end do
      end do
      !
      ne = ne + 1 ! Proceed to next point on this face
      !
    end do
    !
    ! #####################################################
    ! #####################################################
    ! #####################################################
    ! #####  Edge 10 - x/xi  max and y/eta  min edge  #####
    ! #####################################################
    ! #####################################################
    ! #####################################################
    !
    ne = 9*nep1d+1 ! First point on this edge
    !
    xyz(1) = +one ! x/xi   coordinate of current edge point
    xyz(2) = -one ! y/eta  coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do kfa = 1,nep1d
      !
      xyz(3) = xi_ep(kfa) ! z/zeta coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat(n,ne) = eval_LagrangePoly([i,j,k],xyz,xi_sp)
          end do
        end do
      end do
      !
      ne = ne + 1 ! Proceed to next point on this face
      !
    end do
    !
    ! #####################################################
    ! #####################################################
    ! #####################################################
    ! #####  Edge 11 - x/xi  min and y/eta  min edge  #####
    ! #####################################################
    ! #####################################################
    ! #####################################################
    !
    ne = 10*nep1d+1 ! First point on this edge
    !
    xyz(1) = -one ! x/xi   coordinate of current edge point
    xyz(2) = -one ! y/eta  coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do kfa = 1,nep1d
      !
      xyz(3) = xi_ep(kfa) ! z/zeta coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat(n,ne) = eval_LagrangePoly([i,j,k],xyz,xi_sp)
          end do
        end do
      end do
      !
      ne = ne + 1 ! Proceed to next point on this face
      !
    end do
    !
    ! #####################################################
    ! #####################################################
    ! #####################################################
    ! #####  Edge 12 - x/xi  min and y/eta  max edge  #####
    ! #####################################################
    ! #####################################################
    ! #####################################################
    !
    ne = 11*nep1d+1 ! First point on this edge
    !
    xyz(1) = -one ! x/xi   coordinate of current edge point
    xyz(2) = +one ! y/eta  coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do kfa = 1,nep1d
      !
      xyz(3) = xi_ep(kfa) ! z/zeta coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat(n,ne) = eval_LagrangePoly([i,j,k],xyz,xi_sp)
          end do
        end do
      end do
      !
      ne = ne + 1 ! Proceed to next point on this face
      !
    end do
    !
  end if
  !
  ! Use the chop function to copy interp_mat into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( interp_mat ) ! F2003 AUTO-REALLOC
  !
end function InterpToEdge_Hexa
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function InterpToFace_Hexa(xi_sp,xi_fp) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  real(mp), intent(in) :: xi_sp(:)
  real(mp), intent(in) :: xi_fp(:)
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,k,ifa,jfa,kfa,n,nf,ierr
  integer :: nfp1d,nfp2d,nfp3d
  integer :: nsp1d,nsp2d,nsp3d
  !
  !.. Local Arrays ..
  real(mp) :: xyz(1:3)
  !
  !.. Local Allocatable Arrays ..
  real(mp), allocatable :: interp_mat(:,:)
  !
  !.. Local Parameters ..
  real(mp), parameter :: one = 1.0_mp
  !
continue
  !
  ! #########################################################################
  ! #########################################################################
  ! #########################################################################
  ! #########################################################################
  ! #####  NOTE: THE ORDER OF THE FLUX POINTS IS SUCH THAT THE LOWEST   #####
  ! #####        COORDINATE DIRECTION IS THE INNER MOST LOOP, i.e.,     #####
  ! #####          FOR AN X-FACE, INNER LOOP OVER Y, OUTER LOOP OVER Z  #####
  ! #####          FOR A  Y-FACE, INNER LOOP OVER X, OUTER LOOP OVER Z  #####
  ! #####          FOR A  Z-FACE, INNER LOOP OVER X, OUTER LOOP OVER Y  #####
  ! #########################################################################
  ! #########################################################################
  ! #########################################################################
  ! #########################################################################
  !
  nfp1d = size(xi_fp)
  nsp1d = size(xi_sp)
  !
  nfp2d = nfp1d*nfp1d
  nsp2d = nsp1d*nsp1d
  !
  nfp3d = 6*nfp2d
  nsp3d = nsp1d**3
  !
  ! F2003 AUTO-REALLOC
  interp_mat = init_array(nsp3d,nfp3d,real(zero,kind=kind(interp_mat)))
  !
  ! For each flux point, get the coefficients for interpolating
  ! a polynomial from the solution points to the flux point.
  ! These coefficients are simply the value of the 3D Lagrange
  ! polynomial for a specified solution point evaluated at
  ! the coordinates of the flux point.
  !
  ! ############################################################
  ! ############################################################
  ! ############################################################
  ! #####  Face 1 - z/zeta min face of quad (bottom/down)  #####
  ! ############################################################
  ! ############################################################
  ! ############################################################
  !
  nf = 0*nfp2d+1 ! First point on this face
  !
  ! Loops over points on this face
  !
  xyz(3) = -one ! z/zeta coordinate of current face/flux point
  !
  do jfa = 1,nfp1d
    !
    xyz(2) = xi_fp(jfa) ! y/eta coordinate of current face/flux point
    !
    do ifa = 1,nfp1d
      !
      xyz(1) = xi_fp(ifa) ! x/xi coordinate of current face/flux point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! face/flux point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1 ! Equivalent to: n = (k-1)*nsp2d + (j-1)*nsp1d + i
            interp_mat(n,nf) = eval_LagrangePoly([i,j,k],xyz,xi_sp)
          end do
        end do
      end do
      !
      nf = nf + 1 ! Proceed to next point on this face
      !
    end do
  end do
  !
  ! #######################################################
  ! #######################################################
  ! #######################################################
  ! #####  Face 2 - z/zeta max face of quad (top/up)  #####
  ! #######################################################
  ! #######################################################
  ! #######################################################
  !
  nf = 1*nfp2d+1 ! First point on this face
  !
  ! Loops over points on this face
  !
  xyz(3) = one ! z/zeta coordinate of current face/flux point
  !
  do jfa = 1,nfp1d
    !
    xyz(2) = xi_fp(jfa) ! y/eta coordinate of current face/flux point
    !
    do ifa = 1,nfp1d
      !
      xyz(1) = xi_fp(ifa) ! x/xi coordinate of current face/flux point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! face/flux point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1 ! Equivalent to: n = (k-1)*nsp2d + (j-1)*nsp1d + i
            interp_mat(n,nf) = eval_LagrangePoly([i,j,k],xyz,xi_sp)
          end do
        end do
      end do
      !
      nf = nf + 1 ! Proceed to next point on this face
      !
    end do
  end do
  !
  ! #########################################################
  ! #########################################################
  ! #########################################################
  ! #####  Face 3 - x/xi max face of quad (right/east)  #####
  ! #########################################################
  ! #########################################################
  ! #########################################################
  !
  nf = 2*nfp2d+1 ! First point on this face
  !
  ! Loops over points on this face
  !
  xyz(1) = one ! x/xi coordinate of current face/flux point
  !
  do kfa = 1,nfp1d
    !
    xyz(3) = xi_fp(kfa) ! z/zeta coordinate of current face/flux point
    !
    do jfa = 1,nfp1d
      !
      xyz(2) = xi_fp(jfa) ! y/eta coordinate of current face/flux point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! face/flux point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1 ! Equivalent to: n = (k-1)*nsp2d + (j-1)*nsp1d + i
            interp_mat(n,nf) = eval_LagrangePoly([i,j,k],xyz,xi_sp)
          end do
        end do
      end do
      !
      nf = nf + 1 ! Proceed to next point on this face
      !
    end do
  end do
  !
  ! ########################################################
  ! ########################################################
  ! ########################################################
  ! #####  Face 4 - x/xi min face of quad (left/west)  #####
  ! ########################################################
  ! ########################################################
  ! ########################################################
  !
  nf = 3*nfp2d+1 ! First point on this face
  !
  ! Loops over points on this face
  !
  xyz(1) = -one ! x/xi coordinate of current face/flux point
  !
  do kfa = 1,nfp1d
    !
    xyz(3) = xi_fp(kfa) ! z/zeta coordinate of current face/flux point
    !
    do jfa = 1,nfp1d
      !
      xyz(2) = xi_fp(jfa) ! y/eta coordinate of current face/flux point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! face/flux point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1 ! Equivalent to: n = (k-1)*nsp2d + (j-1)*nsp1d + i
            interp_mat(n,nf) = eval_LagrangePoly([i,j,k],xyz,xi_sp)
          end do
        end do
      end do
      !
      nf = nf + 1 ! Proceed to next point on this face
      !
    end do
  end do
  !
  ! ###########################################################
  ! ###########################################################
  ! ###########################################################
  ! #####  Face 5 - y/eta min face of quad (front/south)  #####
  ! ###########################################################
  ! ###########################################################
  ! ###########################################################
  !
  nf = 4*nfp2d+1 ! First point on this face
  !
  ! Loops over points on this face
  !
  xyz(2) = -one ! y/eta coordinate of current face/flux point
  !
  do kfa = 1,nfp1d
    !
    xyz(3) = xi_fp(kfa) ! z/zeta coordinate of current face/flux point
    !
    do ifa = 1,nfp1d
      !
      xyz(1) = xi_fp(ifa) ! x/xi coordinate of current face/flux point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! face/flux point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1 ! Equivalent to: n = (k-1)*nsp2d + (j-1)*nsp1d + i
            interp_mat(n,nf) = eval_LagrangePoly([i,j,k],xyz,xi_sp)
          end do
        end do
      end do
      !
      nf = nf + 1 ! Proceed to next point on this face
      !
    end do
  end do
  !
  ! ##########################################################
  ! ##########################################################
  ! ##########################################################
  ! #####  Face 6 - y/eta max face of quad (back/north)  #####
  ! ##########################################################
  ! ##########################################################
  ! ##########################################################
  !
  nf = 5*nfp2d+1 ! First point on this face
  !
  ! Loops over points on this face
  !
  xyz(2) = one ! y/eta coordinate of current face/flux point
  !
  do kfa = 1,nfp1d
    !
    xyz(3) = xi_fp(kfa) ! z/zeta coordinate of current face/flux point
    !
    do ifa = 1,nfp1d
      !
      xyz(1) = xi_fp(ifa) ! x/xi coordinate of current face/flux point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! face/flux point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1 ! Equivalent to: n = (k-1)*nsp2d + (j-1)*nsp1d + i
            interp_mat(n,nf) = eval_LagrangePoly([i,j,k],xyz,xi_sp)
          end do
        end do
      end do
      !
      nf = nf + 1 ! Proceed to next point on this face
      !
    end do
  end do
  !
  ! Use the chop function to copy interp_mat into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( interp_mat ) ! F2003 AUTO-REALLOC
  !
end function InterpToFace_Hexa
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function InterpToCell_Hexa(xi_old,xi_new) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  real(mp), intent(in) :: xi_old(:)
  real(mp), intent(in) :: xi_new(:)
  !
  !.. Function Return Value ..
  real(wp) :: return_value(1:size(xi_new)**3,1:size(xi_old)**3)
  !
  !.. Local Scalars ..
  integer :: o,oi,oj,ok,opts
  integer :: n,ni,nj,nk,npts
  !
  !.. Local Arrays ..
  real(mp) :: xyz(1:3)
  real(mp) :: interp_mat(1:size(xi_new)**3,1:size(xi_old)**3)
  !
continue
  !
  opts = size(xi_old)
  npts = size(xi_new)
  !
  ! Evaluate the Lagrange polynomial at each interpolation point
  !
  do ok = 1,opts
    do oj = 1,opts
      do oi = 1,opts
        !
        o = (ok-1)*opts*opts + (oj-1)*opts + oi
        !
        do nk = 1,npts
          !
          xyz(3) = xi_new(nk)
          !
          do nj = 1,npts
            !
            xyz(2) = xi_new(nj)
            !
            do ni = 1,npts
              !
              xyz(1) = xi_new(ni)
              !
              n = (nk-1)*npts*npts + (nj-1)*npts + ni
              !
              interp_mat(n,o) = eval_LagrangePoly([oi,oj,ok],xyz,xi_old)
              !
            end do
          end do
        end do
        !
      end do
    end do
  end do
  !
  ! Use the chop function to copy interp_mat into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( interp_mat )
  !
end function InterpToCell_Hexa
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function InterpWeights_Hexa(wts) result(return_value)
  !
  !.. Formal Arguments ..
  real(mp), intent(in) :: wts(:)
  !
  !.. Function Return Value ..
  real(wp) :: return_value(1:size(wts)**3)
  !
  !.. Local Scalars ..
  integer :: m,mi,mj,mk,npts
  !
  !.. Local Arrays ..
  real(mp) :: weights(1:size(wts)**3)
  !
continue
  !
  npts = size(wts)
  !
  ! Evaluate the tensor product of the 1D weights
  !
  do mk = 1,npts
    do mj = 1,npts
      do mi = 1,npts
        !
        m = (mk-1)*npts*npts + (mj-1)*npts + mi
        !
        weights(m) = wts(mi)*wts(mj)*wts(mk)
        !
      end do
    end do
  end do
  !
  ! Use the chop function to copy weights into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( weights )
  !
end function InterpWeights_Hexa
!
!###############################################################################
!
subroutine interpolation_memory_usage(iunit)
  !
  !.. Use Statements ..
  use module_memory
#ifndef DISABLE_DTIO
  use matrix_types_mod, only : matrix_memory_usage
#endif
  !
  !.. Formal Arguments ..
  integer, optional, intent(in) :: iunit
  !
  !.. Local Scalars ..
  integer :: iou,i,j,k
  character(len=36) :: array_name
  type(memory) :: array,dt
  type(memory) :: total
  !
#ifndef DISABLE_DTIO
continue
  !
  iou = iout
  if (present(iunit)) iou = iunit
  !
  call total%reset
  !
  write (iou,1)
  !
  ! Check the memory of the arrays within this module
  !
  ! Memory usage for the interp array of type interpolation_matrix
  !
  if (allocated(interp)) then
    !
    call array%set( storage_size(interp,kind=inttype) , &
                            size(interp,kind=inttype) )
    array_name = "interp"
    write (iou,2) array_name,array%get_units()
    !
    call dt%set(array)
    !
    do j = lbound(interp,dim=2),ubound(interp,dim=2)
      do i = lbound(interp,dim=1),ubound(interp,dim=1)
        !
        if (allocated(interp(i,j)%toFace)) then
          !
          call array%set( storage_size(interp(i,j)%toFace,kind=inttype) , &
                                  size(interp(i,j)%toFace,kind=inttype) )
          !
          write (array_name,4) "interp",Geom_Name(i),j,"toFace"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          do k = lbound(interp(i,j)%toFace,dim=1), &
                 ubound(interp(i,j)%toFace,dim=1)
            !
            array = matrix_memory_usage( interp(i,j)%toFace(k) )
            !
            write (array_name,4) "interp",Geom_Name(i),j,"toFace",k
            write (iou,2) array_name,array%get_units()
            !
            call dt%add(array)
            !
          end do
          !
        end if
        !
        if (allocated(interp(i,j)%toCell)) then
          !
          call array%set( storage_size(interp(i,j)%toCell,kind=inttype) , &
                                  size(interp(i,j)%toCell,kind=inttype) )
          !
          write (array_name,4) "interp",Geom_Name(i),j,"toCell"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          do k = lbound(interp(i,j)%toCell,dim=1), &
                 ubound(interp(i,j)%toCell,dim=1)
            !
            array = matrix_memory_usage( interp(i,j)%toCell(k) )
            !
            write (array_name,4) "interp",Geom_Name(i),j,"toCell",k
            write (iou,2) array_name,array%get_units()
            !
            call dt%add(array)
            !
          end do
          !
        end if
        !
        if (allocated(interp(i,j)%toEdge)) then
          !
          call array%set( storage_size(interp(i,j)%toEdge,kind=inttype) , &
                                  size(interp(i,j)%toEdge,kind=inttype) )
          !
          write (array_name,4) "interp",Geom_Name(i),j,"toEdge"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          do k = lbound(interp(i,j)%toEdge,dim=1), &
                 ubound(interp(i,j)%toEdge,dim=1)
            !
            array = matrix_memory_usage( interp(i,j)%toEdge(k) )
            !
            write (array_name,4) "interp",Geom_Name(i),j,"toEdge",k
            write (iou,2) array_name,array%get_units()
            !
            call dt%add(array)
            !
          end do
          !
        end if
        !
        if (allocated(interp(i,j)%toNode)) then
          !
          call array%set( storage_size(interp(i,j)%toNode,kind=inttype) , &
                                                    int(1,kind=inttype) )
          !
          write (array_name,4) "interp",Geom_Name(i),j,"toNode"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          array = matrix_memory_usage( interp(i,j)%toNode )
          !
          write (array_name,4) "interp",Geom_Name(i),j,"toNode%mat"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
        end if
        !
      end do
    end do
    !
    call total%add(dt)
    write (array_name,5) "interp"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "interp"
    !
  end if
  !
  ! Memory usage for the outerp array of type outerpolation_matrix
  !
  if (allocated(outerp)) then
    !
    call array%set( storage_size(outerp,kind=inttype) , &
                            size(outerp,kind=inttype) )
    array_name = "outerp"
    write (iou,2) array_name,array%get_units()
    !
    call dt%set(array)
    !
    do j = lbound(outerp,dim=2),ubound(outerp,dim=2)
      do i = lbound(outerp,dim=1),ubound(outerp,dim=1)
        !
        if (allocated(outerp(i,j)%init)) then
          !
          call array%set( storage_size(outerp(i,j)%init,kind=inttype) )
          !
          write (array_name,4) "outerp",Geom_Name(i),j,"init"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          array = matrix_memory_usage( outerp(i,j)%init )
          !
          write (array_name,4) "outerp",Geom_Name(i),j,"init%mat"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
        end if
        !
        if (allocated(outerp(i,j)%output)) then
          !
          call array%set( storage_size(outerp(i,j)%output,kind=inttype) )
          !
          write (array_name,4) "outerp",Geom_Name(i),j,"output"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          array = matrix_memory_usage( outerp(i,j)%output )
          !
          write (array_name,4) "outerp",Geom_Name(i),j,"output%mat"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
        end if
        !
        if (allocated(outerp(i,j)%error)) then
          !
          call array%set( storage_size(outerp(i,j)%error,kind=inttype) )
          !
          write (array_name,4) "outerp",Geom_Name(i),j,"error"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          array = matrix_memory_usage( outerp(i,j)%error )
          !
          write (array_name,4) "outerp",Geom_Name(i),j,"error%mat"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
        end if
        !
        if (allocated(outerp(i,j)%err_wts)) then
          !
          call array%set( storage_size(outerp(i,j)%err_wts,kind=inttype) , &
                                  size(outerp(i,j)%err_wts,kind=inttype) )
          !
          write (array_name,4) "outerp",Geom_Name(i),j,"err_wts"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
        end if
        !
      end do
    end do
    !
    call total%add(dt)
    write (array_name,5) "outerp"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "outerp"
    !
  end if
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module interpolation_mod")
  2 format (" Array: ",a,dt)
  3 format (" Total memory = ",f8.2,1x,a,:,a)
  4 format (a,"(",a,",",i0,")%",a,:,"(",i0,")%mat")
  5 format ("'",a,"'")
  6 format (" Array: ",a," memory = ",i0,1x,a,a)
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  !
#endif
end subroutine interpolation_memory_usage
!
!###############################################################################
!
subroutine move_interp_to_device()
  !
  !.. Local Scalars ..
  integer :: i,i1,i2,j,j1,j2,k,k1,k2,m1,m2,n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "move_interp_to_device"
  !
continue
  !
  if (allocated(interp)) then
    !$ACC ENTER DATA COPYIN(interp)
    i1 = lbound(interp,dim=1)
    i2 = ubound(interp,dim=1)
    j1 = lbound(interp,dim=2)
    j2 = ubound(interp,dim=2)
    !$ACC ENTER DATA COPYIN(interp(i1:i2,j1:j2))
    do j = j1,j2
      do i = i1,i2
        !
        !---------------------------------------------------------------
        !
        if (allocated(interp(i,j)%toFace)) then
          !$ACC ENTER DATA COPYIN(interp(i,j)%toFace)
          k1 = lbound(interp(i,j)%toFace,dim=1)
          k2 = ubound(interp(i,j)%toFace,dim=1)
          !$ACC ENTER DATA COPYIN(interp(i,j)%toFace(k1:k2))
          do k = k1,k2
            !
            if (allocated(interp(i,j)%toFace(k)%mat)) then
              !$ACC ENTER DATA COPYIN(interp(i,j)%toFace(k)%mat)
              m1 = lbound(interp(i,j)%toFace(k)%mat,dim=1)
              m2 = ubound(interp(i,j)%toFace(k)%mat,dim=1)
              n1 = lbound(interp(i,j)%toFace(k)%mat,dim=2)
              n2 = ubound(interp(i,j)%toFace(k)%mat,dim=2)
              !$ACC ENTER DATA COPYIN(interp(i,j)%toFace(k)%mat(m1:m2,n1:n2))
            end if
            !
            if (allocated(interp(i,j)%toFace(k)%row_ptr)) then
              !$ACC ENTER DATA COPYIN(interp(i,j)%toFace(k)%row_ptr)
              m1 = lbound(interp(i,j)%toFace(k)%row_ptr,dim=1)
              m2 = ubound(interp(i,j)%toFace(k)%row_ptr,dim=1)
              !$ACC ENTER DATA COPYIN(interp(i,j)%toFace(k)%row_ptr(m1:m2))
            end if
            !
            if (allocated(interp(i,j)%toFace(k)%row_idx)) then
              !$ACC ENTER DATA COPYIN(interp(i,j)%toFace(k)%row_idx)
              m1 = lbound(interp(i,j)%toFace(k)%row_idx,dim=1)
              m2 = ubound(interp(i,j)%toFace(k)%row_idx,dim=1)
              !$ACC ENTER DATA COPYIN(interp(i,j)%toFace(k)%row_idx(m1:m2))
            end if
            !
            if (allocated(interp(i,j)%toFace(k)%row_val)) then
              !$ACC ENTER DATA COPYIN(interp(i,j)%toFace(k)%row_val)
              m1 = lbound(interp(i,j)%toFace(k)%row_val,dim=1)
              m2 = ubound(interp(i,j)%toFace(k)%row_val,dim=1)
              !$ACC ENTER DATA COPYIN(interp(i,j)%toFace(k)%row_val(m1:m2))
            end if
            !
            if (allocated(interp(i,j)%toFace(k)%col_ptr)) then
              !$ACC ENTER DATA COPYIN(interp(i,j)%toFace(k)%col_ptr)
              m1 = lbound(interp(i,j)%toFace(k)%col_ptr,dim=1)
              m2 = ubound(interp(i,j)%toFace(k)%col_ptr,dim=1)
              !$ACC ENTER DATA COPYIN(interp(i,j)%toFace(k)%col_ptr(m1:m2))
            end if
            !
            if (allocated(interp(i,j)%toFace(k)%col_idx)) then
              !$ACC ENTER DATA COPYIN(interp(i,j)%toFace(k)%col_idx)
              m1 = lbound(interp(i,j)%toFace(k)%col_idx,dim=1)
              m2 = ubound(interp(i,j)%toFace(k)%col_idx,dim=1)
              !$ACC ENTER DATA COPYIN(interp(i,j)%toFace(k)%col_idx(m1:m2))
            end if
            !
            if (allocated(interp(i,j)%toFace(k)%col_val)) then
              !$ACC ENTER DATA COPYIN(interp(i,j)%toFace(k)%col_val)
              m1 = lbound(interp(i,j)%toFace(k)%col_val,dim=1)
              m2 = ubound(interp(i,j)%toFace(k)%col_val,dim=1)
              !$ACC ENTER DATA COPYIN(interp(i,j)%toFace(k)%col_val(m1:m2))
            end if
            !
          end do
          !
        end if
        !
        !---------------------------------------------------------------
        !
        if (allocated(interp(i,j)%toCell)) then
          !$ACC ENTER DATA COPYIN(interp(i,j)%toCell)
          k1 = lbound(interp(i,j)%toCell,dim=1)
          k2 = ubound(interp(i,j)%toCell,dim=1)
          !$ACC ENTER DATA COPYIN(interp(i,j)%toCell(k1:k2))
          do k = k1,k2
            !
            if (allocated(interp(i,j)%toCell(k)%mat)) then
              !$ACC ENTER DATA COPYIN(interp(i,j)%toCell(k)%mat)
              m1 = lbound(interp(i,j)%toCell(k)%mat,dim=1)
              m2 = ubound(interp(i,j)%toCell(k)%mat,dim=1)
              n1 = lbound(interp(i,j)%toCell(k)%mat,dim=2)
              n2 = ubound(interp(i,j)%toCell(k)%mat,dim=2)
              !$ACC ENTER DATA COPYIN(interp(i,j)%toCell(k)%mat(m1:m2,n1:n2))
            end if
            !
            if (allocated(interp(i,j)%toCell(k)%row_ptr)) then
              !$ACC ENTER DATA COPYIN(interp(i,j)%toCell(k)%row_ptr)
              m1 = lbound(interp(i,j)%toCell(k)%row_ptr,dim=1)
              m2 = ubound(interp(i,j)%toCell(k)%row_ptr,dim=1)
              !$ACC ENTER DATA COPYIN(interp(i,j)%toCell(k)%row_ptr(m1:m2))
            end if
            !
            if (allocated(interp(i,j)%toCell(k)%row_idx)) then
              !$ACC ENTER DATA COPYIN(interp(i,j)%toCell(k)%row_idx)
              m1 = lbound(interp(i,j)%toCell(k)%row_idx,dim=1)
              m2 = ubound(interp(i,j)%toCell(k)%row_idx,dim=1)
              !$ACC ENTER DATA COPYIN(interp(i,j)%toCell(k)%row_idx(m1:m2))
            end if
            !
            if (allocated(interp(i,j)%toCell(k)%row_val)) then
              !$ACC ENTER DATA COPYIN(interp(i,j)%toCell(k)%row_val)
              m1 = lbound(interp(i,j)%toCell(k)%row_val,dim=1)
              m2 = ubound(interp(i,j)%toCell(k)%row_val,dim=1)
              !$ACC ENTER DATA COPYIN(interp(i,j)%toCell(k)%row_val(m1:m2))
            end if
            !
            if (allocated(interp(i,j)%toCell(k)%col_ptr)) then
              !$ACC ENTER DATA COPYIN(interp(i,j)%toCell(k)%col_ptr)
              m1 = lbound(interp(i,j)%toCell(k)%col_ptr,dim=1)
              m2 = ubound(interp(i,j)%toCell(k)%col_ptr,dim=1)
              !$ACC ENTER DATA COPYIN(interp(i,j)%toCell(k)%col_ptr(m1:m2))
            end if
            !
            if (allocated(interp(i,j)%toCell(k)%col_idx)) then
              !$ACC ENTER DATA COPYIN(interp(i,j)%toCell(k)%col_idx)
              m1 = lbound(interp(i,j)%toCell(k)%col_idx,dim=1)
              m2 = ubound(interp(i,j)%toCell(k)%col_idx,dim=1)
              !$ACC ENTER DATA COPYIN(interp(i,j)%toCell(k)%col_idx(m1:m2))
            end if
            !
            if (allocated(interp(i,j)%toCell(k)%col_val)) then
              !$ACC ENTER DATA COPYIN(interp(i,j)%toCell(k)%col_val)
              m1 = lbound(interp(i,j)%toCell(k)%col_val,dim=1)
              m2 = ubound(interp(i,j)%toCell(k)%col_val,dim=1)
              !$ACC ENTER DATA COPYIN(interp(i,j)%toCell(k)%col_val(m1:m2))
            end if
            !
          end do
          !
        end if
        !
        !---------------------------------------------------------------
        !
        if (allocated(interp(i,j)%toEdge)) then
          !$ACC ENTER DATA COPYIN(interp(i,j)%toEdge)
          k1 = lbound(interp(i,j)%toEdge,dim=1)
          k2 = ubound(interp(i,j)%toEdge,dim=1)
          !$ACC ENTER DATA COPYIN(interp(i,j)%toEdge(k1:k2))
          do k = k1,k2
            !
            if (allocated(interp(i,j)%toEdge(k)%mat)) then
              !$ACC ENTER DATA COPYIN(interp(i,j)%toEdge(k)%mat)
              m1 = lbound(interp(i,j)%toEdge(k)%mat,dim=1)
              m2 = ubound(interp(i,j)%toEdge(k)%mat,dim=1)
              n1 = lbound(interp(i,j)%toEdge(k)%mat,dim=2)
              n2 = ubound(interp(i,j)%toEdge(k)%mat,dim=2)
              !$ACC ENTER DATA COPYIN(interp(i,j)%toEdge(k)%mat(m1:m2,n1:n2))
            end if
            !
            if (allocated(interp(i,j)%toEdge(k)%row_ptr)) then
              !$ACC ENTER DATA COPYIN(interp(i,j)%toEdge(k)%row_ptr)
              m1 = lbound(interp(i,j)%toEdge(k)%row_ptr,dim=1)
              m2 = ubound(interp(i,j)%toEdge(k)%row_ptr,dim=1)
              !$ACC ENTER DATA COPYIN(interp(i,j)%toEdge(k)%row_ptr(m1:m2))
            end if
            !
            if (allocated(interp(i,j)%toEdge(k)%row_idx)) then
              !$ACC ENTER DATA COPYIN(interp(i,j)%toEdge(k)%row_idx)
              m1 = lbound(interp(i,j)%toEdge(k)%row_idx,dim=1)
              m2 = ubound(interp(i,j)%toEdge(k)%row_idx,dim=1)
              !$ACC ENTER DATA COPYIN(interp(i,j)%toEdge(k)%row_idx(m1:m2))
            end if
            !
            if (allocated(interp(i,j)%toEdge(k)%row_val)) then
              !$ACC ENTER DATA COPYIN(interp(i,j)%toEdge(k)%row_val)
              m1 = lbound(interp(i,j)%toEdge(k)%row_val,dim=1)
              m2 = ubound(interp(i,j)%toEdge(k)%row_val,dim=1)
              !$ACC ENTER DATA COPYIN(interp(i,j)%toEdge(k)%row_val(m1:m2))
            end if
            !
            if (allocated(interp(i,j)%toEdge(k)%col_ptr)) then
              !$ACC ENTER DATA COPYIN(interp(i,j)%toEdge(k)%col_ptr)
              m1 = lbound(interp(i,j)%toEdge(k)%col_ptr,dim=1)
              m2 = ubound(interp(i,j)%toEdge(k)%col_ptr,dim=1)
              !$ACC ENTER DATA COPYIN(interp(i,j)%toEdge(k)%col_ptr(m1:m2))
            end if
            !
            if (allocated(interp(i,j)%toEdge(k)%col_idx)) then
              !$ACC ENTER DATA COPYIN(interp(i,j)%toEdge(k)%col_idx)
              m1 = lbound(interp(i,j)%toEdge(k)%col_idx,dim=1)
              m2 = ubound(interp(i,j)%toEdge(k)%col_idx,dim=1)
              !$ACC ENTER DATA COPYIN(interp(i,j)%toEdge(k)%col_idx(m1:m2))
            end if
            !
            if (allocated(interp(i,j)%toEdge(k)%col_val)) then
              !$ACC ENTER DATA COPYIN(interp(i,j)%toEdge(k)%col_val)
              m1 = lbound(interp(i,j)%toEdge(k)%col_val,dim=1)
              m2 = ubound(interp(i,j)%toEdge(k)%col_val,dim=1)
              !$ACC ENTER DATA COPYIN(interp(i,j)%toEdge(k)%col_val(m1:m2))
            end if
            !
          end do
          !
        end if
        !
        !---------------------------------------------------------------
        !
        if (allocated(interp(i,j)%toNode)) then
          !$ACC ENTER DATA COPYIN(interp(i,j)%toNode)
          if (allocated(interp(i,j)%toNode%mat)) then
            !$ACC ENTER DATA COPYIN(interp(i,j)%toNode%mat)
            m1 = lbound(interp(i,j)%toNode%mat,dim=1)
            m2 = ubound(interp(i,j)%toNode%mat,dim=1)
            n1 = lbound(interp(i,j)%toNode%mat,dim=2)
            n2 = ubound(interp(i,j)%toNode%mat,dim=2)
            !$ACC ENTER DATA COPYIN(interp(i,j)%toNode%mat(m1:m2,n1:n2))
          end if
          !
          if (allocated(interp(i,j)%toNode%row_ptr)) then
            !$ACC ENTER DATA COPYIN(interp(i,j)%toNode%row_ptr)
            m1 = lbound(interp(i,j)%toNode%row_ptr,dim=1)
            m2 = ubound(interp(i,j)%toNode%row_ptr,dim=1)
            !$ACC ENTER DATA COPYIN(interp(i,j)%toNode%row_ptr(m1:m2))
          end if
          !
          if (allocated(interp(i,j)%toNode%row_idx)) then
            !$ACC ENTER DATA COPYIN(interp(i,j)%toNode%row_idx)
            m1 = lbound(interp(i,j)%toNode%row_idx,dim=1)
            m2 = ubound(interp(i,j)%toNode%row_idx,dim=1)
            !$ACC ENTER DATA COPYIN(interp(i,j)%toNode%row_idx(m1:m2))
          end if
          !
          if (allocated(interp(i,j)%toNode%row_val)) then
            !$ACC ENTER DATA COPYIN(interp(i,j)%toNode%row_val)
            m1 = lbound(interp(i,j)%toNode%row_val,dim=1)
            m2 = ubound(interp(i,j)%toNode%row_val,dim=1)
            !$ACC ENTER DATA COPYIN(interp(i,j)%toNode%row_val(m1:m2))
          end if
          !
          if (allocated(interp(i,j)%toNode%col_ptr)) then
            !$ACC ENTER DATA COPYIN(interp(i,j)%toNode%col_ptr)
            m1 = lbound(interp(i,j)%toNode%col_ptr,dim=1)
            m2 = ubound(interp(i,j)%toNode%col_ptr,dim=1)
            !$ACC ENTER DATA COPYIN(interp(i,j)%toNode%col_ptr(m1:m2))
          end if
          !
          if (allocated(interp(i,j)%toNode%col_idx)) then
            !$ACC ENTER DATA COPYIN(interp(i,j)%toNode%col_idx)
            m1 = lbound(interp(i,j)%toNode%col_idx,dim=1)
            m2 = ubound(interp(i,j)%toNode%col_idx,dim=1)
            !$ACC ENTER DATA COPYIN(interp(i,j)%toNode%col_idx(m1:m2))
          end if
          !
          if (allocated(interp(i,j)%toNode%col_val)) then
            !$ACC ENTER DATA COPYIN(interp(i,j)%toNode%col_val)
            m1 = lbound(interp(i,j)%toNode%col_val,dim=1)
            m2 = ubound(interp(i,j)%toNode%col_val,dim=1)
            !$ACC ENTER DATA COPYIN(interp(i,j)%toNode%col_val(m1:m2))
          end if
        end if
        !
        !---------------------------------------------------------------
        !
      end do
    end do
  end if
  !
end subroutine move_interp_to_device
!
!###############################################################################
!
subroutine move_outerp_to_device()
  !
  !.. Local Scalars ..
  integer :: i,i1,i2,j,j1,j2,m1,m2,n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "move_outerp_to_device"
  !
continue
  !
  if (allocated(outerp)) then
    !$ACC ENTER DATA COPYIN(outerp)
    i1 = lbound(outerp,dim=1)
    i2 = ubound(outerp,dim=1)
    j1 = lbound(outerp,dim=2)
    j2 = ubound(outerp,dim=2)
    !$ACC ENTER DATA COPYIN(outerp(i1:i2,j1:j2))
    do j = j1,j2
      do i = i1,i2
        !
        !---------------------------------------------------------------
        !
        if (allocated(outerp(i,j)%output)) then
          !$ACC ENTER DATA COPYIN(outerp(i,j)%output)
          if (allocated(outerp(i,j)%output%mat)) then
            !$ACC ENTER DATA COPYIN(outerp(i,j)%output%mat)
            m1 = lbound(outerp(i,j)%output%mat,dim=1)
            m2 = ubound(outerp(i,j)%output%mat,dim=1)
            n1 = lbound(outerp(i,j)%output%mat,dim=2)
            n2 = ubound(outerp(i,j)%output%mat,dim=2)
            !$ACC ENTER DATA COPYIN(outerp(i,j)%output%mat(m1:m2,n1:n2))
          end if
          !
          if (allocated(outerp(i,j)%output%row_ptr)) then
            !$ACC ENTER DATA COPYIN(outerp(i,j)%output%row_ptr)
            m1 = lbound(outerp(i,j)%output%row_ptr,dim=1)
            m2 = ubound(outerp(i,j)%output%row_ptr,dim=1)
            !$ACC ENTER DATA COPYIN(outerp(i,j)%output%row_ptr(m1:m2))
          end if
          !
          if (allocated(outerp(i,j)%output%row_idx)) then
            !$ACC ENTER DATA COPYIN(outerp(i,j)%output%row_idx)
            m1 = lbound(outerp(i,j)%output%row_idx,dim=1)
            m2 = ubound(outerp(i,j)%output%row_idx,dim=1)
            !$ACC ENTER DATA COPYIN(outerp(i,j)%output%row_idx(m1:m2))
          end if
          !
          if (allocated(outerp(i,j)%output%row_val)) then
            !$ACC ENTER DATA COPYIN(outerp(i,j)%output%row_val)
            m1 = lbound(outerp(i,j)%output%row_val,dim=1)
            m2 = ubound(outerp(i,j)%output%row_val,dim=1)
            !$ACC ENTER DATA COPYIN(outerp(i,j)%output%row_val(m1:m2))
          end if
          !
          if (allocated(outerp(i,j)%output%col_ptr)) then
            !$ACC ENTER DATA COPYIN(outerp(i,j)%output%col_ptr)
            m1 = lbound(outerp(i,j)%output%col_ptr,dim=1)
            m2 = ubound(outerp(i,j)%output%col_ptr,dim=1)
            !$ACC ENTER DATA COPYIN(outerp(i,j)%output%col_ptr(m1:m2))
          end if
          !
          if (allocated(outerp(i,j)%output%col_idx)) then
            !$ACC ENTER DATA COPYIN(outerp(i,j)%output%col_idx)
            m1 = lbound(outerp(i,j)%output%col_idx,dim=1)
            m2 = ubound(outerp(i,j)%output%col_idx,dim=1)
            !$ACC ENTER DATA COPYIN(outerp(i,j)%output%col_idx(m1:m2))
          end if
          !
          if (allocated(outerp(i,j)%output%col_val)) then
            !$ACC ENTER DATA COPYIN(outerp(i,j)%output%col_val)
            m1 = lbound(outerp(i,j)%output%col_val,dim=1)
            m2 = ubound(outerp(i,j)%output%col_val,dim=1)
            !$ACC ENTER DATA COPYIN(outerp(i,j)%output%col_val(m1:m2))
          end if
        end if
        !
        !---------------------------------------------------------------
        !
        if (allocated(outerp(i,j)%init)) then
          !$ACC ENTER DATA COPYIN(outerp(i,j)%init)
          if (allocated(outerp(i,j)%init%mat)) then
            !$ACC ENTER DATA COPYIN(outerp(i,j)%init%mat)
            m1 = lbound(outerp(i,j)%init%mat,dim=1)
            m2 = ubound(outerp(i,j)%init%mat,dim=1)
            n1 = lbound(outerp(i,j)%init%mat,dim=2)
            n2 = ubound(outerp(i,j)%init%mat,dim=2)
            !$ACC ENTER DATA COPYIN(outerp(i,j)%init%mat(m1:m2,n1:n2))
          end if
          !
          if (allocated(outerp(i,j)%init%row_ptr)) then
            !$ACC ENTER DATA COPYIN(outerp(i,j)%init%row_ptr)
            m1 = lbound(outerp(i,j)%init%row_ptr,dim=1)
            m2 = ubound(outerp(i,j)%init%row_ptr,dim=1)
            !$ACC ENTER DATA COPYIN(outerp(i,j)%init%row_ptr(m1:m2))
          end if
          !
          if (allocated(outerp(i,j)%init%row_idx)) then
            !$ACC ENTER DATA COPYIN(outerp(i,j)%init%row_idx)
            m1 = lbound(outerp(i,j)%init%row_idx,dim=1)
            m2 = ubound(outerp(i,j)%init%row_idx,dim=1)
            !$ACC ENTER DATA COPYIN(outerp(i,j)%init%row_idx(m1:m2))
          end if
          !
          if (allocated(outerp(i,j)%init%row_val)) then
            !$ACC ENTER DATA COPYIN(outerp(i,j)%init%row_val)
            m1 = lbound(outerp(i,j)%init%row_val,dim=1)
            m2 = ubound(outerp(i,j)%init%row_val,dim=1)
            !$ACC ENTER DATA COPYIN(outerp(i,j)%init%row_val(m1:m2))
          end if
          !
          if (allocated(outerp(i,j)%init%col_ptr)) then
            !$ACC ENTER DATA COPYIN(outerp(i,j)%init%col_ptr)
            m1 = lbound(outerp(i,j)%init%col_ptr,dim=1)
            m2 = ubound(outerp(i,j)%init%col_ptr,dim=1)
            !$ACC ENTER DATA COPYIN(outerp(i,j)%init%col_ptr(m1:m2))
          end if
          !
          if (allocated(outerp(i,j)%init%col_idx)) then
            !$ACC ENTER DATA COPYIN(outerp(i,j)%init%col_idx)
            m1 = lbound(outerp(i,j)%init%col_idx,dim=1)
            m2 = ubound(outerp(i,j)%init%col_idx,dim=1)
            !$ACC ENTER DATA COPYIN(outerp(i,j)%init%col_idx(m1:m2))
          end if
          !
          if (allocated(outerp(i,j)%init%col_val)) then
            !$ACC ENTER DATA COPYIN(outerp(i,j)%init%col_val)
            m1 = lbound(outerp(i,j)%init%col_val,dim=1)
            m2 = ubound(outerp(i,j)%init%col_val,dim=1)
            !$ACC ENTER DATA COPYIN(outerp(i,j)%init%col_val(m1:m2))
          end if
        end if
        !
        !---------------------------------------------------------------
        !
        if (allocated(outerp(i,j)%error)) then
          !$ACC ENTER DATA COPYIN(outerp(i,j)%error)
          if (allocated(outerp(i,j)%error%mat)) then
            !$ACC ENTER DATA COPYIN(outerp(i,j)%error%mat)
            m1 = lbound(outerp(i,j)%error%mat,dim=1)
            m2 = ubound(outerp(i,j)%error%mat,dim=1)
            n1 = lbound(outerp(i,j)%error%mat,dim=2)
            n2 = ubound(outerp(i,j)%error%mat,dim=2)
            !$ACC ENTER DATA COPYIN(outerp(i,j)%error%mat(m1:m2,n1:n2))
          end if
          !
          if (allocated(outerp(i,j)%error%row_ptr)) then
            !$ACC ENTER DATA COPYIN(outerp(i,j)%error%row_ptr)
            m1 = lbound(outerp(i,j)%error%row_ptr,dim=1)
            m2 = ubound(outerp(i,j)%error%row_ptr,dim=1)
            !$ACC ENTER DATA COPYIN(outerp(i,j)%error%row_ptr(m1:m2))
          end if
          !
          if (allocated(outerp(i,j)%error%row_idx)) then
            !$ACC ENTER DATA COPYIN(outerp(i,j)%error%row_idx)
            m1 = lbound(outerp(i,j)%error%row_idx,dim=1)
            m2 = ubound(outerp(i,j)%error%row_idx,dim=1)
            !$ACC ENTER DATA COPYIN(outerp(i,j)%error%row_idx(m1:m2))
          end if
          !
          if (allocated(outerp(i,j)%error%row_val)) then
            !$ACC ENTER DATA COPYIN(outerp(i,j)%error%row_val)
            m1 = lbound(outerp(i,j)%error%row_val,dim=1)
            m2 = ubound(outerp(i,j)%error%row_val,dim=1)
            !$ACC ENTER DATA COPYIN(outerp(i,j)%error%row_val(m1:m2))
          end if
          !
          if (allocated(outerp(i,j)%error%col_ptr)) then
            !$ACC ENTER DATA COPYIN(outerp(i,j)%error%col_ptr)
            m1 = lbound(outerp(i,j)%error%col_ptr,dim=1)
            m2 = ubound(outerp(i,j)%error%col_ptr,dim=1)
            !$ACC ENTER DATA COPYIN(outerp(i,j)%error%col_ptr(m1:m2))
          end if
          !
          if (allocated(outerp(i,j)%error%col_idx)) then
            !$ACC ENTER DATA COPYIN(outerp(i,j)%error%col_idx)
            m1 = lbound(outerp(i,j)%error%col_idx,dim=1)
            m2 = ubound(outerp(i,j)%error%col_idx,dim=1)
            !$ACC ENTER DATA COPYIN(outerp(i,j)%error%col_idx(m1:m2))
          end if
          !
          if (allocated(outerp(i,j)%error%col_val)) then
            !$ACC ENTER DATA COPYIN(outerp(i,j)%error%col_val)
            m1 = lbound(outerp(i,j)%error%col_val,dim=1)
            m2 = ubound(outerp(i,j)%error%col_val,dim=1)
            !$ACC ENTER DATA COPYIN(outerp(i,j)%error%col_val(m1:m2))
          end if
        end if
        !
        !---------------------------------------------------------------
        !
      end do
    end do
  end if
  !
end subroutine move_outerp_to_device
!
!###############################################################################
!
include "Functions/solpts2order.f90"
include "Functions/cell_solpts.f90"
!
!###############################################################################
!
end module interpolation_mod
!
!###############################################################################
!
