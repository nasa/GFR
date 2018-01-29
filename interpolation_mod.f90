module interpolation_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use generic_types_mod, only : matrix
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
  type(interpolation_matrix), public, save, target, allocatable :: interp(:,:)
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
  type(outerpolation_matrix), public, save, target, allocatable :: outerp(:,:)
  !
  type :: spanning_vectors_t
    logical(lk) :: found = fals
    integer :: idx1
    integer :: idx2
    integer :: idx3
  end type spanning_vectors_t
  !
  type(spanning_vectors_t), save, allocatable :: ref_vectors(:,:)
  !
contains
!
!###############################################################################
!
subroutine interpolate_restart_solution(nc,xyz_pt,return_value)
  !
  !.. Use Statements ..
  use order_mod,      only : maxSP
  use order_mod,      only : n_min_geom,n_min_order
  use order_mod,      only : n_max_geom,n_max_order
  use geovar,         only : nr,cell,xyz
  use flowvar,        only : usp
  use quadrature_mod, only : std_elem,solpts_edge
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
  !.. Local Pointers ..
  real(wp), pointer, contiguous :: solpts(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interpolate_restart_solution"
  !
  integer, parameter :: max_iterations = 10
  !
continue
  !
  ! Initialize the local pointers to disassociated
  !
  solpts => null()
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
        solpts => std_elem(this_geom,this_order)%pts
        !
        interp_mat(1:np) = InterpToPoint_Tria(solpts,comp_r)
        !
        if (associated(solpts)) solpts => null()
        !
      case (Geom_Quad)
        !
        interp_mat(1:np) = InterpToPoint_Quad(solpts_edge,comp_r)
        !
      case (Geom_Tetr)
        !
        solpts => std_elem(this_geom,this_order)%pts
        !
        interp_mat(1:np) = InterpToPoint_Tetr(solpts,comp_r)
        !
        if (associated(solpts)) solpts => null()
        !
      case (Geom_Pyra)
        !
        solpts => std_elem(this_geom,this_order)%pts
        !
        interp_mat(1:np) = InterpToPoint_Pyra(solpts,comp_r)
        !
        if (associated(solpts)) solpts => null()
        !
      case (Geom_Pris)
        !
        solpts => std_elem(this_geom,this_order)%pts
        !
        interp_mat(1:np) = InterpToPoint_Pris(solpts,comp_r)
        !
        if (associated(solpts)) solpts => null()
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
    if (errmag < eps3) then
      exit iter_loop
    end if
    !
  end do iter_loop
  !
  ! Use the last interpolation matrix that was created to interpolate
  ! the solution in the current cell to the input point
  !
  do m = 1,size(usp,dim=1)
    return_value(m) = dot( interp_mat(1:np) , usp(m,n1:n2) )
  end do
  !
end subroutine interpolate_restart_solution
!
!###############################################################################
!
pure function find_ref_vectors(this_geom,this_order) result(return_value)
  !
  !.. Use Statements ..
  use order_mod, only : geom_solpts
  !
  !.. Formal Arguments ..
  integer, intent(in) :: this_geom
  integer, intent(in) :: this_order
  !
  !.. Function Result ..
  type(spanning_vectors_t) :: return_value
  !
  !.. Local Scalars ..
  integer :: np_edge,np_tria,np_quad
  !
  !.. Local Parameters ..
 !logical(lk), parameter :: use_short_vectors = fals
  logical(lk), parameter :: use_short_vectors = true
  !
continue
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
        return_value%idx1 = 2
        return_value%idx2 = np_edge+1
        return_value%idx3 = 0
        !
      case (Geom_Quad)
        !
        return_value%idx1 = 2
        return_value%idx2 = np_edge+1
        return_value%idx3 = 0
        !
      case (Geom_Tetr)
        !
        return_value%idx1 = 2
        return_value%idx2 = np_edge+1
        return_value%idx3 = np_tria+1
        !
      case (Geom_Pyra)
        !
        return_value%idx1 = 2
        return_value%idx2 = np_edge+1
        return_value%idx3 = np_quad+1
        !
      case (Geom_Pris)
        !
        return_value%idx1 = 2
        return_value%idx2 = np_edge+1
        return_value%idx3 = np_tria+1
        !
      case (Geom_Hexa)
        !
        return_value%idx1 = 2
        return_value%idx2 = np_edge+1
        return_value%idx3 = np_quad+1
        !
      case default
        !
        return_value%idx1 = 2
        return_value%idx2 = 0
        return_value%idx3 = 0
        !
    end select
    !
  else
    !
    select case (this_geom)
      !
      case (Geom_Tria)
        !
        return_value%idx1 = np_edge
        return_value%idx2 = np_tria
        return_value%idx3 = 0
        !
      case (Geom_Quad)
        !
        return_value%idx1 = np_edge
        return_value%idx2 = np_edge*(np_edge-1) + 1
        return_value%idx3 = 0
        !
      case (Geom_Tetr)
        !
        return_value%idx1 = np_edge
        return_value%idx2 = np_tria
        return_value%idx3 = geom_solpts(Geom_Tetr,this_order)
        !
      case (Geom_Pyra)
        !
        return_value%idx1 = np_edge
        return_value%idx2 = np_edge*(np_edge-1) + 1
        return_value%idx3 = geom_solpts(Geom_Pyra,this_order)
        !
      case (Geom_Pris)
        !
        return_value%idx1 = np_edge
        return_value%idx2 = np_tria
        return_value%idx3 = np_tria*(np_edge-1) + 1
        !
      case (Geom_Hexa)
        !
        return_value%idx1 = np_edge
        return_value%idx2 = np_edge*(np_edge-1) + 1
        return_value%idx3 = np_quad*(np_edge-1) + 1
        !
      case default
        !
        return_value%idx1 = np_edge
        return_value%idx2 = 0
        return_value%idx3 = 0
        !
    end select
    !
  end if
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
  integer :: i
  !
continue
  !
  do i = 1,size(xi_sp)
    return_value(i) = eval_LagrangePoly( [i] , xi_pt , xi_sp )
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
  integer :: i,j,n,nsp
  !
continue
  !
  nsp = size(xi_sp)
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
      return_value(n) = eval_LagrangePoly( [i,j] , xi_pt , xi_sp )
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
 !nfp = np2n( np ) + 1
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
  integer :: i,j,k,n,nsp
  !
continue
  !
  nsp = size(xi_sp)
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
        return_value(n) = eval_LagrangePoly( [i,j,k] , xi_pt , xi_sp )
      end do
    end do
  end do
  !
end function InterpToPoint_Hexa
!
!###############################################################################
!
subroutine init_interpolation_matrices()
  !
  !.. Use Statements ..
  use order_mod,      only : geom_solpts,geom_flxpts
  use order_mod,      only : n_min_geom,n_max_geom
  use order_mod,      only : n_min_order,n_max_order
  use order_mod,      only : geom_mio,geom_edgpts
  use quadrature_mod, only : geom_is_used
  use quadrature_mod, only : std_elem,face_elem
  !
  !.. Local Scalars ..
  integer :: gmin,gmax,omin,omax,ierr
  integer :: this_omin,this_omax
  integer :: nsolpts,this_order,this_geom
  integer :: nflxpts,face_order
  integer :: nedgpts,edge_order
  integer :: nnodpts,node_order
  character(len=200) :: array_name
  !
  !.. Local Pointers ..
  type(interpolation_matrix), pointer :: this_interp
  type(matrix), pointer :: this_face
  type(matrix), pointer :: this_edge
  type(matrix), pointer :: this_node
 !real(wp), pointer, contiguous :: solpts(:,:)
  real(wp), pointer :: solpts(:,:)
 !real(wp), pointer, contiguous :: flxpts(:,:)
 !real(wp), pointer :: flxpts(:,:)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: edgpts(:)
  real(wp), allocatable :: flxpts(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_interpolation_matrices"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Initialize the local pointers to disassociated
  !
  this_interp => null()
  this_face   => null()
  this_edge   => null()
  this_node   => null()
  solpts      => null()
 !flxpts      => null()
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
      ! Assign local pointer as alias to simplify the code
      !
      this_interp => interp(this_geom,this_order)
      !
      ! Get the number of solution points for this cell
      !
      nsolpts = geom_solpts(this_geom,this_order)
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
      allocate ( this_interp%toFace(this_omin:this_omax) , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"toFace"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Loop through all the interpolation orders for this combo
      ! and compute the corresponding interpolation matrices
      !
      face_order_loop: do face_order = this_omin,this_omax
        !
        ! Assign local pointers to these matrices
        ! as aliases to simplify the code
        !
        this_face => this_interp%toFace(face_order)
        !
        ! Get the number of face/flux points for this face order
        !
        nflxpts = geom_flxpts(this_geom,face_order)
        !
        ! Allocate the interpolation matrix for this face order
        !
        allocate ( this_face%mat(1:nsolpts,1:nflxpts) , &
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
            edgpts = std_elem(this_geom,this_order)%pts(1,:) ! F03 AUTO-REALLOC
           !solpts => std_elem(this_geom,this_order)%pts
            !
            this_face%mat = InterpToFace_Edge( edgpts )
           !this_face%mat = InterpToFace_Edge( solpts(1,:) )
            !
          case (Geom_Tria)
            !
            solpts => std_elem(this_geom,this_order)%pts
            !
            this_face%mat = InterpToFace_Tria( solpts )
            !
          case (Geom_Quad)
            !
            edgpts = std_elem(Geom_Edge,this_order)%pts(1,:) ! F03 AUTO-REALLOC
           !solpts => std_elem(Geom_Edge,this_order)%pts
            flxpts = face_elem(Geom_Edge,face_order)%pts(1,:) ! F03 AUTO-REALLOC
           !flxpts => face_elem(Geom_Edge,face_order)%pts
            !
            this_face%mat = InterpToFace_Quad( edgpts , flxpts )
           !this_face%mat = InterpToFace_Quad( solpts(1,:) , flxpts(1,:) )
            !
          case (Geom_Tetr)
            !
            solpts => std_elem(this_geom,this_order)%pts
            !
            this_face%mat = InterpToFace_Tetr( solpts )
            !
          case (Geom_Pyra)
            !
            solpts => std_elem(this_geom,this_order)%pts
            !
            this_face%mat = InterpToFace_Pyra( solpts )
            !
          case (Geom_Pris)
            !
            solpts => std_elem(this_geom,this_order)%pts
            !
            this_face%mat = InterpToFace_Pris( solpts )
            !
          case (Geom_Hexa)
            !
            edgpts = std_elem(Geom_Edge,this_order)%pts(1,:) ! F03 AUTO-REALLOC
           !solpts => std_elem(Geom_Edge,this_order)%pts
            flxpts = face_elem(Geom_Edge,face_order)%pts(1,:) ! F03 AUTO-REALLOC
           !flxpts => face_elem(Geom_Edge,face_order)%pts
            !
            this_face%mat = InterpToFace_Hexa( edgpts , flxpts )
           !this_face%mat = InterpToFace_Hexa( solpts(1,:) , flxpts(1,:) )
            !
          case default
            !
            write (error_message,2)
            call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
            !
        end select
        !
        ! Remove all associations of the local pointers before continuing
        !
        if (associated(this_face)) this_face => null()
        if (associated(solpts   )) solpts    => null()
       !if (associated(flxpts   )) flxpts    => null()
        !
      end do face_order_loop
      !
      ! Allocate the interp%toEdge component of interp
      !
      allocate ( this_interp%toEdge(this_omin:this_omax) , &
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
        ! Assign local pointers to these matrices
        ! as aliases to simplify the code
        !
        this_edge => this_interp%toEdge(edge_order)
        !
        ! Allocate the interpolation matrix for this edge order
        !
        allocate ( this_edge%mat(1:nsolpts,1:nedgpts) , &
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
            solpts => std_elem(this_geom,this_order)%pts
            !
            this_edge%mat = InterpToEdge_Tetr( solpts )
            !
          case (Geom_Pyra)
            !
            solpts => std_elem(this_geom,this_order)%pts
            !
            this_edge%mat = InterpToEdge_Pyra( solpts )
            !
          case (Geom_Pris)
            !
            solpts => std_elem(this_geom,this_order)%pts
            !
            this_edge%mat = InterpToEdge_Pris( solpts )
            !
          case (Geom_Hexa)
            !
            edgpts = std_elem(Geom_Edge,this_order)%pts(1,:) ! F03 AUTO-REALLOC
           !solpts => std_elem(Geom_Edge,this_order)%pts
            flxpts = face_elem(Geom_Edge,edge_order)%pts(1,:) ! F03 AUTO-REALLOC
           !flxpts => face_elem(Geom_Edge,edge_order)%pts
            !
            this_edge%mat = InterpToEdge_Hexa( edgpts , flxpts )
           !this_edge%mat = InterpToEdge_Hexa( solpts(1,:) , flxpts(1,:) )
            !
          case default
            !
            write (error_message,2)
            call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
            !
        end select
        !
        ! Remove all associations of the local pointers before continuing
        !
        if (associated(this_edge)) this_edge => null()
        if (associated(solpts   )) solpts    => null()
       !if (associated(flxpts   )) flxpts    => null()
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
        allocate ( this_interp%toNode , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"toNode"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message)
        !
        ! Assign local pointers to these matrices
        ! as aliases to simplify the code
        !
        this_node => this_interp%toNode
        !
        ! Get the number of nodes for this cell
        !
        nnodpts = geom_nodes(this_geom)
        !
        ! Allocate the interpolation matrix for the nodes
        !
        allocate ( this_node%mat(1:nsolpts,1:nnodpts) , &
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
            solpts => std_elem(this_geom,this_order)%pts
            !
            this_node%mat = InterpToNode_Tria( solpts )
            !
          case (Geom_Quad)
            !
            edgpts = std_elem(Geom_Edge,this_order)%pts(1,:) ! F03 AUTO-REALLOC
           !solpts => std_elem(Geom_Edge,this_order)%pts
            !
            this_node%mat = InterpToNode_Quad( edgpts )
           !this_node%mat = InterpToNode_Quad( solpts(1,:) )
            !
          case (Geom_Tetr)
            !
            solpts => std_elem(this_geom,this_order)%pts
            !
            this_node%mat = InterpToNode_Tetr( solpts )
            !
          case (Geom_Pyra)
            !
            solpts => std_elem(this_geom,this_order)%pts
            !
            this_node%mat = InterpToNode_Pyra( solpts )
            !
          case (Geom_Pris)
            !
            solpts => std_elem(this_geom,this_order)%pts
            !
            this_node%mat = InterpToNode_Pris( solpts )
            !
          case (Geom_Hexa)
            !
            edgpts = std_elem(Geom_Edge,this_order)%pts(1,:) ! F03 AUTO-REALLOC
           !solpts => std_elem(Geom_Edge,this_order)%pts
            !
            this_node%mat = InterpToNode_Hexa( edgpts )
           !this_node%mat = InterpToNode_Hexa( solpts(1,:) )
            !
          case default
            !
            write (error_message,2)
            call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
            !
        end select
        !
        ! Remove all associations of the local pointers before continuing
        !
        if (associated(solpts   )) solpts    => null()
       !if (associated(flxpts   )) flxpts    => null()
        if (associated(this_node)) this_node => null()
        !
      end do node_loop
      !
      ! Remove all associations of the local pointers before continuing
      !
      if (associated(this_interp)) this_interp => null()
      !
    end do order_loop
    !
  end do geom_loop
  !
  do this_geom = gmin,gmax
    do this_order = omin,omax
      if (allocated(interp(this_geom,this_order)%toFace)) then
        call interp(this_geom,this_order)%toFace%check_for_sparsity
      end if
      if (allocated(interp(this_geom,this_order)%toEdge)) then
        call interp(this_geom,this_order)%toEdge%check_for_sparsity
      end if
      if (allocated(interp(this_geom,this_order)%toNode)) then
        call interp(this_geom,this_order)%toNode%check_for_sparsity
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
  !.. Local Pointers ..
  type(outerpolation_matrix), pointer :: this_outerp
 !real(wp), pointer, contiguous :: imat(:,:)
  real(wp), pointer :: imat(:,:)
 !real(wp), pointer, contiguous :: omat(:,:)
  real(wp), pointer :: omat(:,:)
 !real(wp), pointer, contiguous :: emat(:,:)
  real(wp), pointer :: emat(:,:)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: pts(:)
  real(wp), allocatable :: wts(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_outerpolation_matrices"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Initialize the local pointers to disassociated
  !
  this_outerp => null()
  imat        => null()
  omat        => null()
  emat        => null()
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
      ! Assign local pointer as alias to simplify the code
      !
      this_outerp => outerp(this_geom,this_order)
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
        allocate ( this_outerp%init , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"init"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
        !
        ! Allocate the mat component of outerp%init
        !
        allocate ( this_outerp%init%mat(1:npts,1:npts) , source=zero , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"init%mat"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
        !
        ! Assign local pointersto these matrices
        ! as aliases to simplify the code
        !
        imat => this_outerp%init%mat
        !
        ! Get the 1D location of the initialization points
        !
        allocate ( pts(1:this_order+1) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"pts",1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
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
        imat(:,:) = InterpToCell(this_geom,this_order,npts,npts,pts)
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
        allocate ( this_outerp%output , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"output"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
        !
        ! Allocate the mat component of outerp%output
        !
        allocate ( this_outerp%output%mat(1:opts,1:npts) , source=zero , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"output%mat"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
        !
        ! Assign local pointersto these matrices
        ! as aliases to simplify the code
        !
        omat => this_outerp%output%mat
        !
        ! Get the 1D location of the output points
        !
        allocate ( pts(1:output_order+1) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"pts",1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
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
          pts(:) = real(intseq(-output_order,output_order,2),kind=wp) / &
                   real(max(output_order,1),kind=wp)
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
        omat(:,:) = InterpToCell(this_geom,this_order,npts,opts,pts)
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
        allocate ( this_outerp%error , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"error"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
        !
        ! Allocate the mat component of outerp%error
        !
        allocate ( this_outerp%error%mat(1:epts,1:npts) , source=zero , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"error%mat"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
        !
        ! Assign local pointersto these matrices
        ! as aliases to simplify the code
        !
        emat => this_outerp%error%mat
        !
        ! Get the 1D location of the error points
        ! and their corresponding quadrature weights
        !
        allocate ( pts(1:error_order+1) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"pts",1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
        !
        allocate ( wts(1:error_order+1) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"wts",1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
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
        emat(:,:) = InterpToCell(this_geom,this_order,npts,epts,pts)
        !
        ! Allocate the error_wts component of outerp
        !
        allocate ( this_outerp%err_wts(1:epts) , source=zero , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"err_wts"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
        !
        ! Compute the error weights
        !
        this_outerp%err_wts(:) = InterpWeights(this_geom,epts,wts)
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
      ! Remove all associations of the local pointers before continuing
      !
      if (associated(this_outerp)) this_outerp => null()
      if (associated(imat       )) imat        => null()
      if (associated(omat       )) omat        => null()
      if (associated(emat       )) emat        => null()
      !
    end do order_loop
    !
  end do geom_loop
  !
  do this_geom = gmin,gmax
    do this_order = omin,omax
      if (allocated(outerp(this_geom,this_order)%init)) then
        call outerp(this_geom,this_order)%init%check_for_sparsity
      end if
      if (allocated(outerp(this_geom,this_order)%output)) then
        call outerp(this_geom,this_order)%output%check_for_sparsity
      end if
      if (allocated(outerp(this_geom,this_order)%error)) then
        call outerp(this_geom,this_order)%error%check_for_sparsity
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
  use quadrature_mod, only : std_elem
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: this_geom
  integer,  intent(in) :: this_order
  integer,  intent(in) :: opts
  integer,  intent(in) :: npts
  real(wp), intent(in) :: pts(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:npts,1:opts)
  !
  !.. Local Pointers ..
 !real(wp), pointer, contiguous :: solpts(:,:)
  real(wp), pointer :: solpts(:,:)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: edgpts(:)
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
      edgpts = std_elem(Geom_Edge,this_order)%pts(1,:) ! F2003 AUTO-REALLOC
     !solpts => std_elem(Geom_Edge,this_order)%pts
      !
      return_value(:,:) = InterpToCell_Edge( edgpts , pts )
     !return_value(:,:) = InterpToCell_Edge( solpts(1,:) , pts )
      !
    case (Geom_Tria)
      !
      ! Not sure about triangles here so just use the identity matrix
      !
      return_value(:,:) = identity_matrix( return_value )
      !
    case (Geom_Quad)
      !
      edgpts = std_elem(Geom_Edge,this_order)%pts(1,:) ! F2003 AUTO-REALLOC
     !solpts => std_elem(Geom_Edge,this_order)%pts
      !
      return_value(:,:) = InterpToCell_Quad( edgpts , pts )
     !return_value(:,:) = InterpToCell_Quad( solpts(1,:) , pts )
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
      edgpts = std_elem(Geom_Edge,this_order)%pts(1,:) ! F2003 AUTO-REALLOC
     !solpts => std_elem(Geom_Edge,this_order)%pts
      !
      return_value(:,:) = InterpToCell_Hexa( edgpts , pts )
     !return_value(:,:) = InterpToCell_Hexa( solpts(1,:) , pts )
      !
    case default
      !
      write (error_message,1)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      !
  end select
  !
  if (associated(solpts)) solpts => null()
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
  real(wp), intent(in) :: wts(:)
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
pure function generic_interpolation_matrix(basis,newbasis,wmat) &
                                    result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: basis
  real(wp), dimension(:,:), intent(in) :: newbasis
  real(wp), dimension(:,:), intent(in) :: wmat
  !
  !.. Function Result ..
  real(wp), dimension(1:size(wmat,1),1:size(newbasis,1)) :: return_value
  !
  !.. Local Arrays ..
  real(qp), dimension(1:size(wmat,1),1:size(newbasis,1)) :: interp_mat_qp
  real(qp), dimension(1:size(basis,2),1:size(basis,1)) :: basisT_qp
  real(qp), dimension(1:size(basis,2),1:size(basis,2)) :: pmat_qp
  real(qp), dimension(1:size(basis,2),1:size(wmat,2)) :: tmat_qp
  real(qp), dimension(1:size(basis,1),1:size(basis,2)) :: basis_qp
  real(qp), dimension(1:size(newbasis,1),1:size(newbasis,2)) :: newbasis_qp
  real(qp), dimension(1:size(wmat,1),1:size(wmat,2)) :: wmat_qp
  !
continue
  !
  basis_qp = real( basis , kind=qp )
  newbasis_qp = real( newbasis , kind=qp )
  wmat_qp = real( wmat , kind=qp )
  !
  basisT_qp = matmul( transpose(basis_qp) , wmat_qp )
  basisT_qp = chop( basisT_qp , basisT_qp )
  !
  pmat_qp = matmul( basisT_qp , basis_qp )
  pmat_qp = chop( pmat_qp , pmat_qp )
  !
  pmat_qp = invert_matrix( pmat_qp )
  pmat_qp = chop( pmat_qp , pmat_qp )
  !
  tmat_qp = matmul( pmat_qp , basisT_qp )
  tmat_qp = chop( tmat_qp , tmat_qp )
  !
  interp_mat_qp = transpose( matmul( newbasis_qp , tmat_qp ) )
  interp_mat_qp = chop( interp_mat_qp , interp_mat_qp )
  !
  ! Use the chop function to copy interp_mat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( interp_mat_qp )
  !
end function generic_interpolation_matrix
!
!###############################################################################
!
pure function InterpToFace_Edge(xi) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(wp), dimension(1:size(xi),1:2) :: return_value
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Arrays ..
  real(qp), dimension(1:size(xi)) :: xi_qp
  real(qp), dimension(1:size(xi),1:2) :: interp_mat_qp
  !
continue
  !
  ! Create a quad precision copy of the array giving the solution points
  !
  xi_qp = real( xi , kind=qp )
  !
  do i = 1,size(xi)
    interp_mat_qp(i,1) = eval_LagrangePoly([i],[-qone],xi_qp)
    interp_mat_qp(i,2) = eval_LagrangePoly([i],[qone],xi_qp)
  end do
  !
  ! Use the chop function to copy interp_mat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( interp_mat_qp )
  !
end function InterpToFace_Edge
!
!###############################################################################
!
pure function InterpToCell_Edge(xi,xi_interp) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xi
  real(wp), dimension(:), intent(in) :: xi_interp
  !
  !.. Function Return Value ..
  real(wp), dimension(1:size(xi),1:size(xi)) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,npts
  !
  !.. Local Arrays ..
  real(qp), dimension(1:size(xi)) :: xi_qp
  real(qp), dimension(1:size(xi)) :: pts_qp
  real(qp), dimension(1:size(xi),1:size(xi)) :: interp_mat_qp
  !
continue
  !
  npts = size(xi)
  !
  ! Create a quad precision copy of the array giving the solution points
  !
  xi_qp = real( xi , kind=qp )
  !
  ! Create a quad precision copy of the array giving the interpolation points
  !
  pts_qp = real( xi_interp , kind=qp )
  !
  ! Evaluate the Lagrange polynomial at each interpolation point
  !
  do j = 1,npts
    do i = 1,npts
      interp_mat_qp(i,j) = eval_LagrangePoly([j],[pts_qp(i)],xi_qp)
    end do
  end do
  !
  ! Use the chop function to copy interp_mat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( interp_mat_qp )
  !
end function InterpToCell_Edge
!
!###############################################################################
!
pure function InterpToNode_Quad(xi_sp) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xi_sp
  !
  !.. Function Return Value ..
  real(wp), dimension(1:size(xi_sp)**2,1:4) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,n,np,nfp,nsp
  !
  !.. Local Arrays ..
  real(qp), dimension(1:2) :: xy_qp
  real(qp), dimension(1:2,1:4) :: nodes_qp
  real(qp), dimension(1:size(xi_sp)) :: solpts_qp
  real(qp), dimension(1:size(xi_sp)**2,1:4) :: interp_mat_qp
  !
continue
  !
  nsp = size(xi_sp)
  !
  ! Fill in the nodes for the reference element
  !
  nodes_qp(1:2,1) = [-qone,-qone]
  nodes_qp(1:2,2) = [ qone,-qone]
  nodes_qp(1:2,3) = [ qone, qone]
  nodes_qp(1:2,4) = [-qone, qone]
  !
  ! Create quad precision versions of the arrays containing
  ! the locations of the solution points
  !
  solpts_qp = real( xi_sp , kind=qp )
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
    xy_qp(1:2) = nodes_qp(1:2,np)
    !
    n = 0
    do j = 1,nsp
      do i = 1,nsp
        n = n + 1
        interp_mat_qp(n,np) = eval_LagrangePoly([i,j],xy_qp,solpts_qp)
      end do
    end do
    !
  end do
  !
  ! Use the chop function to copy interp_mat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( interp_mat_qp )
  !
end function InterpToNode_Quad
!
!###############################################################################
!
pure function InterpToFace_Quad(xi_sp,xi_fp) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xi_sp
  real(wp), dimension(:), intent(in) :: xi_fp
  !
  !.. Function Return Value ..
  real(wp), dimension(1:size(xi_sp)**2,1:4*size(xi_fp)) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,n,nf,nfp,nsp
  !
  !.. Local Arrays ..
  real(qp), dimension(1:2) :: xy_qp
  real(qp), dimension(1:size(xi_sp)) :: solpts_qp
  real(qp), dimension(1:size(xi_fp)) :: flxpts_qp
  real(qp), dimension(1:size(xi_sp)**2,1:4*size(xi_fp)) :: interp_mat_qp
  !
continue
  !
  nfp = size(xi_fp)
  nsp = size(xi_sp)
  !
  ! Create quad precision versions of the arrays containing
  ! the locations of the solution and flux points
  !
  solpts_qp = real( xi_sp , kind=qp )
  flxpts_qp = real( xi_fp , kind=qp )
  !
  ! For each flux point, get the coefficients for interpolating
  ! a polynomial from the solution points to the flux point.
  ! These coefficients are simply the value of the 2D Lagrange
  ! polynomial for a specified solution point evaluated at
  ! the coordinates of the flux point.
  !
  ! Face 1 - bottom/south edge of quad
  !
  xy_qp(2) = -qone
  do n = 1,nfp
    !
    nf = n
    !
    xy_qp(1) = flxpts_qp(n)
    !
    do j = 1,nsp
      do i = 1,nsp
        interp_mat_qp(nsp*(j-1)+i,nf) = eval_LagrangePoly([i,j],xy_qp,solpts_qp)
      end do
    end do
    !
  end do
  !
  ! Face 2 - right/east edge of quad
  !
  xy_qp(1) = qone
  do n = 1,nfp
    !
    nf = n + nfp
    !
    xy_qp(2) = flxpts_qp(n)
    !
    do j = 1,nsp
      do i = 1,nsp
        interp_mat_qp(nsp*(j-1)+i,nf) = eval_LagrangePoly([i,j],xy_qp,solpts_qp)
      end do
    end do
    !
  end do
  !
  ! Face 3 - top/north edge of quad
  !
  xy_qp(2) = qone
  do n = 1,nfp
    !
    nf = n + 2*nfp
    !
    xy_qp(1) = flxpts_qp(nfp+1-n)
    !
    do j = 1,nsp
      do i = 1,nsp
        interp_mat_qp(nsp*(j-1)+i,nf) = eval_LagrangePoly([i,j],xy_qp,solpts_qp)
      end do
    end do
    !
  end do
  !
  ! Face 4 - left/west edge of quad
  !
  xy_qp(1) = -qone
  do n = 1,nfp
    !
    nf = n + 3*nfp
    !
    xy_qp(2) = flxpts_qp(nfp+1-n)
    !
    do j = 1,nsp
      do i = 1,nsp
        interp_mat_qp(nsp*(j-1)+i,nf) = eval_LagrangePoly([i,j],xy_qp,solpts_qp)
      end do
    end do
    !
  end do
  !
  ! Use the chop function to copy interp_mat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( interp_mat_qp )
  !
end function InterpToFace_Quad
!
!###############################################################################
!
pure function InterpToCell_Quad(xi_old,xi_new) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xi_old
  real(wp), dimension(:), intent(in) :: xi_new
  !
  !.. Function Return Value ..
  real(wp), dimension(1:size(xi_new)**2,1:size(xi_old)**2) :: return_value
  !
  !.. Local Scalars ..
  integer :: o,oi,oj,opts
  integer :: n,ni,nj,npts
  !
  !.. Local Arrays ..
  real(qp), dimension(1:2) :: xy_qp
  real(qp), dimension(1:size(xi_old)) :: xi_old_qp
  real(qp), dimension(1:size(xi_new)) :: xi_new_qp
  real(qp), dimension(1:size(xi_new)**2,1:size(xi_old)**2) :: interp_mat_qp
  !
continue
  !
  opts = size(xi_old)
  npts = size(xi_new)
  !
  ! Create a quad precision copy of the array giving the solution points
  !
  xi_old_qp = real( xi_old , kind=qp )
  !
  ! Create a quad precision copy of the array giving the interpolation points
  !
  xi_new_qp = real( xi_new , kind=qp )
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
        xy_qp(2) = xi_new_qp(nj)
        !
        do ni = 1,npts
          !
          xy_qp(1) = xi_new_qp(ni)
          !
          n = (nj-1)*npts + ni
          !
          interp_mat_qp(n,o) = eval_LagrangePoly([oi,oj],xy_qp,xi_old_qp)
          !
        end do
      end do
      !
    end do
  end do
  !
  ! Use the chop function to copy interp_mat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( interp_mat_qp )
  !
end function InterpToCell_Quad
!
!###############################################################################
!
pure function InterpWeights_Quad(wts) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: wts
  !
  !.. Function Return Value ..
  real(wp), dimension(1:size(wts)**2) :: return_value
  !
  !.. Local Scalars ..
  integer :: m,mi,mj,npts
  !
  !.. Local Arrays ..
  real(qp), dimension(1:size(wts)) :: wts_qp
  real(qp), dimension(1:size(wts)**2) :: weights_qp
  !
continue
  !
  npts = size(wts)
  !
  ! Create a quad precision copy of the array giving the 1D quadrature weights
  !
  wts_qp = real( wts , kind=qp )
  !
  ! Evaluate the tensor product of the 1D weights
  !
  do mj = 1,npts
    do mi = 1,npts
      !
      m = (mj-1)*npts + mi
      !
      weights_qp(m) = wts_qp(mi)*wts_qp(mj)
      !
    end do
  end do
  !
  ! Use the chop function to copy weights_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( weights_qp )
  !
end function InterpWeights_Quad
!
!###############################################################################
!
pure function InterpToNode_Tria(rs) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Function Return Argument ..
  real(wp), dimension(1:size(rs,2),1:3) :: return_value
  !
  !.. Local Scalars ..
  integer :: np,nfp
  !
continue
  !
  np = size(rs,2)
  nfp = np2n( np ) + 1
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
pure function InterpToFace_Tria(rs) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Function Return Argument ..
  real(wp), dimension(1:size(rs,2),1:3*(np2n(size(rs,2))+1)) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,np,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:np2n(size(rs,2))+1,1:3) :: Fmask
  !
continue
  !
  np = size(rs,2)
  nfp = np2n( np ) + 1
  !
  return_value(:,:) = zero
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
pure function InterpToNode_Tetr(rs) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Function Return Argument ..
  real(wp), dimension(1:size(rs,2),1:4) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,n,nsp3d,nfp3d,nep3d
  !
continue
  !
  nsp3d = size(rs,2)
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
pure function InterpToEdge_Tetr(rs) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Function Return Argument ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,n,nep3d,nsp3d,ierr
  !
continue
  !
  nsp3d = size(rs,2)
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
pure function InterpToFace_Tetr(rs) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Function Return Argument ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,n,nfp3d,nsp3d,ierr
  !
continue
  !
  nsp3d = size(rs,2)
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
pure function InterpToNode_Pyra(rs) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Function Return Argument ..
  real(wp), dimension(1:size(rs,2),1:5) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,n,nfp3d,nsp3d,ierr
  !
continue
  !
  nsp3d = size(rs,2)
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
pure function InterpToEdge_Pyra(rs) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Function Return Argument ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,n,nep3d,nsp3d,ierr
  !
continue
  !
  nsp3d = size(rs,2)
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
pure function InterpToFace_Pyra(rs) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Function Return Argument ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,n,nfp3d,nsp3d,ierr
  !
continue
  !
  nsp3d = size(rs,2)
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
pure function InterpToNode_Pris(rs) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Function Return Argument ..
  real(wp), dimension(1:size(rs,2),1:6) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,n,nfp3d,nsp3d,ierr
  !
continue
  !
  nsp3d = size(rs,2)
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
pure function InterpToEdge_Pris(rs) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Function Return Argument ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,n,nep3d,nsp3d,ierr
  !
continue
  !
  nsp3d = size(rs,2)
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
pure function InterpToFace_Pris(rs) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Function Return Argument ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,n,nfp3d,nsp3d,ierr
  !
continue
  !
  nsp3d = size(rs,2)
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
pure function InterpToNode_Hexa(xi_sp) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xi_sp
  !
  !.. Function Return Value ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,k,n,np,ierr
  integer :: nsp1d,nsp3d
  !
  !.. Local Arrays ..
  real(qp), dimension(1:3) :: xyz_qp
  real(qp), dimension(1:3,1:8) :: nodes_qp
  real(qp), dimension(1:size(xi_sp)) :: solpts_qp
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: interp_mat_qp(:,:)
  !
continue
  !
  nsp1d = size(xi_sp)
  nsp3d = nsp1d**3
  !
  ! Fill in the nodes for the reference element
  !
  nodes_qp(1:3,1) = [-qone,-qone,-qone]
  nodes_qp(1:3,2) = [ qone,-qone,-qone]
  nodes_qp(1:3,3) = [ qone, qone,-qone]
  nodes_qp(1:3,4) = [-qone, qone,-qone]
  nodes_qp(1:3,5) = [-qone,-qone, qone]
  nodes_qp(1:3,6) = [ qone,-qone, qone]
  nodes_qp(1:3,7) = [ qone, qone, qone]
  nodes_qp(1:3,8) = [-qone, qone, qone]
  !
  ! Allocate the return array containing the
  ! working precision interpolation matrix
  !
  allocate ( return_value(1:nsp3d,1:8) , source=zero , stat=ierr )
  !
  ! Create the interpolation matrix using quad precision
  !
  allocate ( interp_mat_qp(1:nsp3d,1:8) , source=qzero , stat=ierr )
  !
  ! Create quad precision versions of the arrays containing
  ! the locations of the solution points
  !
  solpts_qp = real( xi_sp , kind=qp )
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
    xyz_qp(1:3) = nodes_qp(1:3,np)
    !
    n = 0
    do k = 1,nsp1d
      do j = 1,nsp1d
        do i = 1,nsp1d
          n = n + 1
          interp_mat_qp(n,np) = eval_LagrangePoly([i,j,k],xyz_qp,solpts_qp)
        end do
      end do
    end do
    !
  end do
  !
  ! Use the chop function to copy interp_mat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value(:,:) = chop( interp_mat_qp(:,:) )
  !
end function InterpToNode_Hexa
!
!###############################################################################
!
pure function InterpToEdge_Hexa(xi_sp,xi_ep) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xi_sp
  real(wp), dimension(:), intent(in) :: xi_ep
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
  real(qp), dimension(1:3) :: xyz_qp
  real(qp), dimension(1:size(xi_sp)) :: solpts_qp
  real(qp), dimension(1:size(xi_ep)) :: edgpts_qp
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: interp_mat_qp(:,:)
  !
  !.. Local Parameters ..
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
  ! Allocate the return array containing the
  ! working precision interpolation matrix
  !
  allocate ( return_value(1:nsp3d,1:nep3d) , source=zero , stat=ierr )
  !
  ! Create the interpolation matrix using quad precision
  !
  allocate ( interp_mat_qp(1:nsp3d,1:nep3d) , source=qzero , stat=ierr )
  !
  ! Create quad precision versions of the arrays containing
  ! the locations of the solution and edge points
  !
  solpts_qp = real( xi_sp , kind=qp )
  edgpts_qp = real( xi_ep , kind=qp )
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
      xyz_qp(1) = merge( -qone , qone , unit_cell(1,n1) == 0 )
      xyz_qp(2) = merge( -qone , qone , unit_cell(2,n1) == 0 )
      xyz_qp(3) = merge( -qone , qone , unit_cell(3,n1) == 0 )
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
        xyz_qp(l) = edgpts_qp(lp)
        !
        n = 0
        do k = 1,nsp1d
          do j = 1,nsp1d
            do i = 1,nsp1d
              n = n + 1
              interp_mat_qp(n,ne) = eval_LagrangePoly([i,j,k],xyz_qp,solpts_qp)
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
    xyz_qp(1) = -qone ! x/xi   coordinate of current edge point
    xyz_qp(3) = -qone ! z/zeta coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do jfa = 1,nep1d
      !
      xyz_qp(2) = edgpts_qp(jfa) ! y/eta coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat_qp(n,ne) = eval_LagrangePoly([i,j,k],xyz_qp,solpts_qp)
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
    xyz_qp(2) = +qone ! y/eta  coordinate of current edge point
    xyz_qp(3) = -qone ! z/zeta coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do ifa = 1,nep1d
      !
      xyz_qp(1) = edgpts_qp(ifa) ! x/xi  coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat_qp(n,ne) = eval_LagrangePoly([i,j,k],xyz_qp,solpts_qp)
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
    xyz_qp(1) = +qone ! x/xi   coordinate of current edge point
    xyz_qp(3) = -qone ! z/zeta coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do jfa = 1,nep1d
      !
      xyz_qp(2) = edgpts_qp(jfa) ! y/eta coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat_qp(n,ne) = eval_LagrangePoly([i,j,k],xyz_qp,solpts_qp)
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
    xyz_qp(2) = -qone ! y/eta  coordinate of current edge point
    xyz_qp(3) = -qone ! z/zeta coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do ifa = 1,nep1d
      !
      xyz_qp(1) = edgpts_qp(ifa) ! x/xi  coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat_qp(n,ne) = eval_LagrangePoly([i,j,k],xyz_qp,solpts_qp)
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
    xyz_qp(2) = -qone ! y/eta  coordinate of current edge point
    xyz_qp(3) = +qone ! z/zeta coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do ifa = 1,nep1d
      !
      xyz_qp(1) = edgpts_qp(ifa) ! x/xi  coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat_qp(n,ne) = eval_LagrangePoly([i,j,k],xyz_qp,solpts_qp)
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
    xyz_qp(1) = +qone ! x/xi   coordinate of current edge point
    xyz_qp(3) = +qone ! z/zeta coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do jfa = 1,nep1d
      !
      xyz_qp(2) = edgpts_qp(jfa) ! y/eta coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat_qp(n,ne) = eval_LagrangePoly([i,j,k],xyz_qp,solpts_qp)
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
    xyz_qp(2) = +qone ! y/eta  coordinate of current edge point
    xyz_qp(3) = +qone ! z/zeta coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do ifa = 1,nep1d
      !
      xyz_qp(1) = edgpts_qp(ifa) ! x/xi  coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat_qp(n,ne) = eval_LagrangePoly([i,j,k],xyz_qp,solpts_qp)
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
    xyz_qp(1) = -qone ! x/xi   coordinate of current edge point
    xyz_qp(3) = +qone ! z/zeta coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do jfa = 1,nep1d
      !
      xyz_qp(2) = edgpts_qp(jfa) ! y/eta coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat_qp(n,ne) = eval_LagrangePoly([i,j,k],xyz_qp,solpts_qp)
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
    xyz_qp(1) = +qone ! x/xi   coordinate of current edge point
    xyz_qp(2) = +qone ! y/eta  coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do kfa = 1,nep1d
      !
      xyz_qp(3) = edgpts_qp(kfa) ! z/zeta coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat_qp(n,ne) = eval_LagrangePoly([i,j,k],xyz_qp,solpts_qp)
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
    xyz_qp(1) = +qone ! x/xi   coordinate of current edge point
    xyz_qp(2) = -qone ! y/eta  coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do kfa = 1,nep1d
      !
      xyz_qp(3) = edgpts_qp(kfa) ! z/zeta coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat_qp(n,ne) = eval_LagrangePoly([i,j,k],xyz_qp,solpts_qp)
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
    xyz_qp(1) = -qone ! x/xi   coordinate of current edge point
    xyz_qp(2) = -qone ! y/eta  coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do kfa = 1,nep1d
      !
      xyz_qp(3) = edgpts_qp(kfa) ! z/zeta coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat_qp(n,ne) = eval_LagrangePoly([i,j,k],xyz_qp,solpts_qp)
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
    xyz_qp(1) = -qone ! x/xi   coordinate of current edge point
    xyz_qp(2) = +qone ! y/eta  coordinate of current edge point
    !
    ! Loops over points on this edge
    !
    do kfa = 1,nep1d
      !
      xyz_qp(3) = edgpts_qp(kfa) ! z/zeta coordinate of current edge point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! edge point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1
            interp_mat_qp(n,ne) = eval_LagrangePoly([i,j,k],xyz_qp,solpts_qp)
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
  ! Use the chop function to copy interp_mat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( interp_mat_qp )
  !
end function InterpToEdge_Hexa
!
!###############################################################################
!
pure function InterpToFace_Hexa(xi_sp,xi_fp) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xi_sp
  real(wp), dimension(:), intent(in) :: xi_fp
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
  real(qp), dimension(1:3) :: xyz_qp
  real(qp), dimension(1:size(xi_sp)) :: solpts_qp
  real(qp), dimension(1:size(xi_fp)) :: flxpts_qp
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: interp_mat_qp(:,:)
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
  ! Allocate the return array containing the
  ! working precision interpolation matrix
  !
  allocate ( return_value(1:nsp3d,1:nfp3d) , source=zero , stat=ierr )
  !
  ! Create the interpolation matrix using quad precision
  !
  allocate ( interp_mat_qp(1:nsp3d,1:nfp3d) , source=qzero , stat=ierr )
  !
  ! Create quad precision versions of the arrays containing
  ! the locations of the solution and flux points
  !
  solpts_qp = real( xi_sp , kind=qp )
  flxpts_qp = real( xi_fp , kind=qp )
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
  xyz_qp(3) = -qone ! z/zeta coordinate of current face/flux point
  !
  do jfa = 1,nfp1d
    !
    xyz_qp(2) = flxpts_qp(jfa) ! y/eta coordinate of current face/flux point
    !
    do ifa = 1,nfp1d
      !
      xyz_qp(1) = flxpts_qp(ifa) ! x/xi coordinate of current face/flux point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! face/flux point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1 ! Equivalent to: n = (k-1)*nsp2d + (j-1)*nsp1d + i
            interp_mat_qp(n,nf) = eval_LagrangePoly([i,j,k],xyz_qp,solpts_qp)
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
  xyz_qp(3) = qone ! z/zeta coordinate of current face/flux point
  !
  do jfa = 1,nfp1d
    !
    xyz_qp(2) = flxpts_qp(jfa) ! y/eta coordinate of current face/flux point
    !
    do ifa = 1,nfp1d
      !
      xyz_qp(1) = flxpts_qp(ifa) ! x/xi coordinate of current face/flux point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! face/flux point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1 ! Equivalent to: n = (k-1)*nsp2d + (j-1)*nsp1d + i
            interp_mat_qp(n,nf) = eval_LagrangePoly([i,j,k],xyz_qp,solpts_qp)
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
  xyz_qp(1) = qone ! x/xi coordinate of current face/flux point
  !
  do kfa = 1,nfp1d
    !
    xyz_qp(3) = flxpts_qp(kfa) ! z/zeta coordinate of current face/flux point
    !
    do jfa = 1,nfp1d
      !
      xyz_qp(2) = flxpts_qp(jfa) ! y/eta coordinate of current face/flux point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! face/flux point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1 ! Equivalent to: n = (k-1)*nsp2d + (j-1)*nsp1d + i
            interp_mat_qp(n,nf) = eval_LagrangePoly([i,j,k],xyz_qp,solpts_qp)
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
  xyz_qp(1) = -qone ! x/xi coordinate of current face/flux point
  !
  do kfa = 1,nfp1d
    !
    xyz_qp(3) = flxpts_qp(kfa) ! z/zeta coordinate of current face/flux point
    !
    do jfa = 1,nfp1d
      !
      xyz_qp(2) = flxpts_qp(jfa) ! y/eta coordinate of current face/flux point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! face/flux point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1 ! Equivalent to: n = (k-1)*nsp2d + (j-1)*nsp1d + i
            interp_mat_qp(n,nf) = eval_LagrangePoly([i,j,k],xyz_qp,solpts_qp)
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
  xyz_qp(2) = -qone ! y/eta coordinate of current face/flux point
  !
  do kfa = 1,nfp1d
    !
    xyz_qp(3) = flxpts_qp(kfa) ! z/zeta coordinate of current face/flux point
    !
    do ifa = 1,nfp1d
      !
      xyz_qp(1) = flxpts_qp(ifa) ! x/xi coordinate of current face/flux point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! face/flux point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1 ! Equivalent to: n = (k-1)*nsp2d + (j-1)*nsp1d + i
            interp_mat_qp(n,nf) = eval_LagrangePoly([i,j,k],xyz_qp,solpts_qp)
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
  xyz_qp(2) = qone ! y/eta coordinate of current face/flux point
  !
  do kfa = 1,nfp1d
    !
    xyz_qp(3) = flxpts_qp(kfa) ! z/zeta coordinate of current face/flux point
    !
    do ifa = 1,nfp1d
      !
      xyz_qp(1) = flxpts_qp(ifa) ! x/xi coordinate of current face/flux point
      !
      ! Loop to get the value of the Lagrange polynomial for the current
      ! face/flux point at each of the interior solution points
      !
      n = 0
      do k = 1,nsp1d
        do j = 1,nsp1d
          do i = 1,nsp1d
            n = n + 1 ! Equivalent to: n = (k-1)*nsp2d + (j-1)*nsp1d + i
            interp_mat_qp(n,nf) = eval_LagrangePoly([i,j,k],xyz_qp,solpts_qp)
          end do
        end do
      end do
      !
      nf = nf + 1 ! Proceed to next point on this face
      !
    end do
  end do
  !
  ! Use the chop function to copy interp_mat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( interp_mat_qp )
  !
end function InterpToFace_Hexa
!
!###############################################################################
!
pure function InterpToCell_Hexa(xi_old,xi_new) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xi_old
  real(wp), dimension(:), intent(in) :: xi_new
  !
  !.. Function Return Value ..
  real(wp), dimension(1:size(xi_new)**3,1:size(xi_old)**3) :: return_value
  !
  !.. Local Scalars ..
  integer :: o,oi,oj,ok,opts
  integer :: n,ni,nj,nk,npts
  !
  !.. Local Arrays ..
  real(qp), dimension(1:3) :: xyz_qp
  real(qp), dimension(1:size(xi_old)) :: xi_old_qp
  real(qp), dimension(1:size(xi_new)) :: xi_new_qp
  real(qp), dimension(1:size(xi_new)**3,1:size(xi_old)**3) :: interp_mat_qp
  !
continue
  !
  opts = size(xi_old)
  npts = size(xi_new)
  !
  ! Create a quad precision copy of the array giving the solution points
  !
  xi_old_qp = real( xi_old , kind=qp )
  !
  ! Create a quad precision copy of the array giving the interpolation points
  !
  xi_new_qp = real( xi_new , kind=qp )
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
          xyz_qp(3) = xi_new_qp(nk)
          !
          do nj = 1,npts
            !
            xyz_qp(2) = xi_new_qp(nj)
            !
            do ni = 1,npts
              !
              xyz_qp(1) = xi_new_qp(ni)
              !
              n = (nk-1)*npts*npts + (nj-1)*npts + ni
              !
              interp_mat_qp(n,o) = eval_LagrangePoly([oi,oj,ok], &
                                                     xyz_qp,xi_old_qp)
              !
            end do
          end do
        end do
        !
      end do
    end do
  end do
  !
  ! Use the chop function to copy interp_mat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( interp_mat_qp )
  !
end function InterpToCell_Hexa
!
!###############################################################################
!
pure function InterpWeights_Hexa(wts) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: wts
  !
  !.. Function Return Value ..
  real(wp), dimension(1:size(wts)**3) :: return_value
  !
  !.. Local Scalars ..
  integer :: m,mi,mj,mk,npts
  !
  !.. Local Arrays ..
  real(qp), dimension(1:size(wts)) :: wts_qp
  real(qp), dimension(1:size(wts)**3) :: weights_qp
  !
continue
  !
  npts = size(wts)
  !
  ! Create a quad precision copy of the array giving the 1D quadrature weights
  !
  wts_qp = real( wts , kind=qp )
  !
  ! Evaluate the tensor product of the 1D weights
  !
  do mk = 1,npts
    do mj = 1,npts
      do mi = 1,npts
        !
        m = (mk-1)*npts*npts + (mj-1)*npts + mi
        !
        weights_qp(m) = wts_qp(mi)*wts_qp(mj)*wts_qp(mk)
        !
      end do
    end do
  end do
  !
  ! Use the chop function to copy weights_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( weights_qp )
  !
end function InterpWeights_Hexa
!
!###############################################################################
!
subroutine interpolation_memory_usage(iunit)
  !
  !.. Use Statements ..
  use module_memory
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
            array = interp(i,j)%toFace(k)%memory_usage()
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
            array = interp(i,j)%toCell(k)%memory_usage()
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
            array = interp(i,j)%toEdge(k)%memory_usage()
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
          array = interp(i,j)%toNode%memory_usage()
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
          array = outerp(i,j)%init%memory_usage()
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
          array = outerp(i,j)%output%memory_usage()
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
          array = outerp(i,j)%error%memory_usage()
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
include "Functions/np2n.f90"
include "Functions/cell_solpts.f90"
!
!###############################################################################
!
end module interpolation_mod
!
!###############################################################################
!
