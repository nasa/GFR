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
    ! gauss : Interpolate to gauss quadrature points
    type(matrix), allocatable :: gauss
    ! lobatto : Interpolate to lobatto quadrature points
    type(matrix), allocatable :: lobatto
    ! equidist : Interpolate to equidistant points
    type(matrix), allocatable :: equidist
    ! output : Pointer to quadrature points used to
    !          output the solution for post-processing
    type(matrix), pointer :: output => null()
    ! init : Pointer to quadrature points used
    !        to compute initial conditions
    type(matrix), pointer :: init => null()
    ! error : Pointer to quadrature points used for error computations
    type(matrix), pointer :: error => null()
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
contains
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
  real(wp), pointer, contiguous :: solpts(:,:)
  real(wp), pointer, contiguous :: flxpts(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_interpolation_matrices"
  !
continue
  !
  ! Initialize the local pointers to disassociated
  !
  this_interp => null()
  this_face   => null()
  this_edge   => null()
  this_node   => null()
  solpts      => null()
  flxpts      => null()
  !
  gmin = n_min_geom
  gmax = n_max_geom
  omin = n_min_order
  omax = n_max_order
  !
  ! Allocate the interp array
  !
  allocate ( interp(gmin:gmax,omin:omax) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"interp",1,ierr,error_message)
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
      call alloc_error(pname,array_name,1,ierr,error_message)
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
        call alloc_error(pname,array_name,1,ierr,error_message)
        !
        ! Compute the interpolation matrices depending on geometry type
        !
        select case (this_geom)
          case (Geom_Edge)
            !
            solpts => std_elem(this_geom,this_order)%pts
            !
            this_face%mat = InterpToFace_Edge( solpts(1,:) )
            !
          case (Geom_Tria)
            !
            solpts => std_elem(this_geom,this_order)%pts
            !
            this_face%mat = InterpToFace_Tria( solpts )
            !
          case (Geom_Quad)
            !
            solpts => std_elem(Geom_Edge,this_order)%pts
            flxpts => face_elem(Geom_Edge,face_order)%pts
            !
            this_face%mat = InterpToFace_Quad( solpts(1,:) , flxpts(1,:) )
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
            solpts => std_elem(Geom_Edge,this_order)%pts
            flxpts => face_elem(Geom_Edge,face_order)%pts
            !
            this_face%mat = InterpToFace_Hexa( solpts(1,:) , flxpts(1,:) )
            !
          case default
            !
            write (error_message,2)
            call stop_recon(stop_mpi,pname,error_message)
            !
        end select
        !
        ! Remove all associations of the local pointers before continuing
        !
        if (associated(this_face)) this_face => null()
        if (associated(solpts   )) solpts    => null()
        if (associated(flxpts   )) flxpts    => null()
        !
      end do face_order_loop
      !
      ! Allocate the interp%toEdge component of interp
      !
      allocate ( this_interp%toEdge(this_omin:this_omax) , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"toEdge"
      call alloc_error(pname,array_name,1,ierr,error_message)
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
        call alloc_error(pname,array_name,1,ierr,error_message)
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
            solpts => std_elem(Geom_Edge,this_order)%pts
            flxpts => face_elem(Geom_Edge,edge_order)%pts
            !
            this_edge%mat = InterpToEdge_Hexa( solpts(1,:) , flxpts(1,:) )
            !
          case default
            !
            write (error_message,2)
            call stop_recon(stop_mpi,pname,error_message)
            !
        end select
        !
        ! Remove all associations of the local pointers before continuing
        !
        if (associated(this_edge)) this_edge => null()
        if (associated(solpts   )) solpts    => null()
        if (associated(flxpts   )) flxpts    => null()
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
        call alloc_error(pname,array_name,1,ierr,error_message)
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
        call alloc_error(pname,array_name,1,ierr,error_message)
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
            solpts => std_elem(Geom_Edge,this_order)%pts
            !
            this_node%mat = InterpToNode_Quad( solpts(1,:) )
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
            solpts => std_elem(Geom_Edge,this_order)%pts
            !
            this_node%mat = InterpToNode_Hexa( solpts(1,:) )
            !
          case default
            !
            write (error_message,2)
            call stop_recon(stop_mpi,pname,error_message)
            !
        end select
        !
        ! Remove all associations of the local pointers before continuing
        !
        if (associated(solpts   )) solpts    => null()
        if (associated(flxpts   )) flxpts    => null()
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
  use ovar,            only : iop_solpts,loc_output_pts
  use ovar,            only : loc_error_pts,loc_init_pts
  use order_mod,       only : geom_solpts,geom_flxpts
  use order_mod,       only : n_min_geom,n_max_geom
  use order_mod,       only : n_min_order,n_max_order
  use order_mod,       only : geom_mio
  use quadrature_mod,  only : geom_is_used
  use quadrature_mod,  only : std_elem,face_elem
  use polynomial_mod,  only : legendre_gauss_quadrature
  use polynomial_mod,  only : legendre_gauss_lobatto_quadrature
  !
  !.. Local Scalars ..
  integer :: gmin,gmax,omin,omax,ierr
  integer :: this_geom,this_order,face_order
  integer :: this_omin,this_omax
  integer :: npts
  logical :: gauss_matrix_is_needed
  logical :: lobatto_matrix_is_needed
  logical :: equidist_matrix_is_needed
  logical :: error_weights_are_needed
  character(len=200) :: array_name
  !
  !.. Local Pointers ..
  type(outerpolation_matrix), pointer :: this_outerp
  real(wp), pointer, contiguous :: gmat(:,:)
  real(wp), pointer, contiguous :: lmat(:,:)
  real(wp), pointer, contiguous :: emat(:,:)
  real(wp), pointer, contiguous :: solpts(:,:)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: gauss_pts(:)
  real(wp), allocatable :: gauss_wts(:)
  real(wp), allocatable :: lobatto_pts(:)
  real(wp), allocatable :: lobatto_wts(:)
  real(wp), allocatable :: equidist_pts(:)
  real(wp), allocatable :: error_wts(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_outerpolation_matrices"
  !
continue
  !
  ! Initialize the local pointers to disassociated
  !
  this_outerp => null()
  gmat        => null()
  lmat        => null()
  emat        => null()
  solpts      => null()
  !
  gmin = n_min_geom
  gmax = n_max_geom
  omin = n_min_order
  omax = n_max_order
  !
  ! Allocate the outerp array
  !
  allocate ( outerp(gmin:gmax,omin:omax) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"outerp",1,ierr,error_message)
  !
  ! Check if we need to create the interpolation matrices within outerp
  !
  gauss_matrix_is_needed = fals
  lobatto_matrix_is_needed = fals
  error_weights_are_needed = fals
  !
  ! Check if we need either matrix for computing the error
  !
  if (iop_solpts /= loc_error_pts) then
    error_weights_are_needed = true
    if (loc_error_pts == Legendre_Gauss) gauss_matrix_is_needed = true
    if (loc_error_pts == Legendre_Gauss_Lobatto) lobatto_matrix_is_needed = true
  end if
  !
  ! Check if we need either matrix to compute the initial conditions
  !
  if (iop_solpts /= loc_init_pts) then
    if (loc_init_pts == Legendre_Gauss) gauss_matrix_is_needed = true
    if (loc_init_pts == Legendre_Gauss_Lobatto) lobatto_matrix_is_needed = true
  end if
  !
  ! Check if we need either matrix to output the solution for post-processing
  !
  if (iop_solpts /= loc_output_pts) then
    if (loc_output_pts == Legendre_Gauss) then
      gauss_matrix_is_needed = true
    end if
    if (loc_output_pts == Legendre_Gauss_Lobatto) then
      lobatto_matrix_is_needed = true
    end if
    if (loc_output_pts == Equi_Distant_Points) then
      equidist_matrix_is_needed = true
    end if
  end if
  !
  ! If neither matrix is needed, just return from this subroutine
  !
  if (.not.gauss_matrix_is_needed .and. .not.lobatto_matrix_is_needed) return
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
      ! Allocate the outerp%gauss component of outerp if needed
      !
      if (gauss_matrix_is_needed) then
        !
        allocate ( this_outerp%gauss , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"gauss"
        call alloc_error(pname,array_name,1,ierr,error_message,skip_alloc_pause)
        !
        allocate ( this_outerp%gauss%mat(1:npts,1:npts) , &
                   source=zero , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"gauss%mat"
        call alloc_error(pname,array_name,1,ierr,error_message)
        !
        ! Assign local pointers to these matrices
        ! as aliases to simplify the code
        !
        gmat => this_outerp%gauss%mat
        !
        ! Pre-compute the 1D Legendre Gauss points and their
        ! corresponding quadrature weights for this order
        !
        allocate ( gauss_pts(1:this_order+1) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"gauss_pts",1,ierr,error_message, &
                         skip_alloc_pause)
        !
        allocate ( gauss_wts(1:this_order+1) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"gauss_wts",1,ierr,error_message, &
                         skip_alloc_pause)
        !
        call legendre_gauss_quadrature(gauss_pts,gauss_wts)
        !
      end if
      !
      ! Allocate the outerp%lobatto component of outerp if needed
      !
      if (lobatto_matrix_is_needed) then
        !
        allocate ( this_outerp%lobatto , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"lobatto"
        call alloc_error(pname,array_name,1,ierr,error_message,skip_alloc_pause)
        !
        allocate ( this_outerp%lobatto%mat(1:npts,1:npts) , &
                   source=zero , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"lobatto%mat"
        call alloc_error(pname,array_name,1,ierr,error_message)
        !
        ! Assign local pointers to these matrices
        ! as aliases to simplify the code
        !
        lmat => this_outerp%lobatto%mat
        !
        ! Pre-compute the 1D Legendre Gauss Lobatto points and their
        ! corresponding quadrature weights for this order
        !
        allocate ( lobatto_pts(1:this_order+1) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"lobatto_pts",1,ierr,error_message, &
                         skip_alloc_pause)
        !
        allocate ( lobatto_wts(1:this_order+1) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"lobatto_wts",1,ierr,error_message, &
                         skip_alloc_pause)
        !
        call legendre_gauss_lobatto_quadrature(lobatto_pts,lobatto_wts)
        !
      end if
      !
      ! Allocate the outerp%equidist component of outerp if needed
      !
      if (equidist_matrix_is_needed) then
        !
        allocate ( this_outerp%equidist , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"equidist"
        call alloc_error(pname,array_name,1,ierr,error_message,skip_alloc_pause)
        !
        allocate ( this_outerp%equidist%mat(1:npts,1:npts) , &
                   source=zero , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"equidist%mat"
        call alloc_error(pname,array_name,1,ierr,error_message)
        !
        ! Assign local pointers to these matrices
        ! as aliases to simplify the code
        !
        emat => this_outerp%equidist%mat
        !
        ! Pre-compute the 1D equidistant points for this order
        ! NOTE: Equidistant points should only be used as output points
        !       where the quadrature weights arent needed, so we dont
        !       need to compute the weights for the equidistant points.
        !       
        allocate ( equidist_pts(1:this_order+1) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"equidist_pts",1,ierr,error_message, &
                         skip_alloc_pause)
        !
        equidist_pts(:) = real(intseq(-this_order,this_order,2),kind=wp) / &
                          real(max(this_order,1),kind=wp)
        !
      end if
      !
      ! Allocate the err_wts component of outerp along with the weights
      ! of the 1D error points if the error weights are needed
      !
      if (error_weights_are_needed) then
        !
        allocate ( this_outerp%err_wts(1:npts) , source=zero , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"err_wts"
        call alloc_error(pname,array_name,1,ierr,error_message)
        !
        allocate ( error_wts(1:this_order+1) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"error_wts",1,ierr,error_message, &
                         skip_alloc_pause)
        !
        if (loc_error_pts == Legendre_Gauss) then
          error_wts = gauss_wts
        else if (loc_error_pts == Legendre_Gauss_Lobatto) then
          error_wts = lobatto_wts
        end if
        !
      end if
      !
      ! Compute the outerpolation matrices depending on geometry type
      !
      select case (this_geom)
        case (Geom_Edge)
          !
          solpts => std_elem(Geom_Edge,this_order)%pts
          !
          if (gauss_matrix_is_needed) then
            gmat = InterpToCell_Edge( solpts(1,:) , gauss_pts )
          end if
          !
          if (lobatto_matrix_is_needed) then
            lmat = InterpToCell_Edge( solpts(1,:) , lobatto_pts )
          end if
          !
          if (equidist_matrix_is_needed) then
            emat = InterpToCell_Edge( solpts(1,:) , equidist_pts )
          end if
          !
          if (error_weights_are_needed) then
            this_outerp%err_wts = error_wts
          end if
          !
        case (Geom_Tria)
          !
          ! Not sure about triangles here so just use the identity matrix
          !
          if (gauss_matrix_is_needed) gmat = identity_matrix( gmat )
          if (lobatto_matrix_is_needed) lmat = identity_matrix( lmat )
          if (equidist_matrix_is_needed) emat = identity_matrix( emat )
          !
        case (Geom_Quad)
          !
          solpts => std_elem(Geom_Edge,this_order)%pts
          !
          if (gauss_matrix_is_needed) then
            gmat = InterpToCell_Quad( solpts(1,:) , gauss_pts )
          end if
          !
          if (lobatto_matrix_is_needed) then
            lmat = InterpToCell_Quad( solpts(1,:) , lobatto_pts )
          end if
          !
          if (equidist_matrix_is_needed) then
            emat = InterpToCell_Quad( solpts(1,:) , equidist_pts )
          end if
          !
          if (error_weights_are_needed) then
            this_outerp%err_wts = InterpWeights_Quad( error_wts )
          end if
          !
        case (Geom_Tetr)
          !
          ! Not sure about tetrahedra here so just use the identity matrix
          !
          if (gauss_matrix_is_needed) gmat = identity_matrix( gmat )
          if (lobatto_matrix_is_needed) lmat = identity_matrix( lmat )
          if (equidist_matrix_is_needed) emat = identity_matrix( emat )
          !
        case (Geom_Pyra)
          !
          ! Not sure about pyramids here so just use the identity matrix
          !
          if (gauss_matrix_is_needed) gmat = identity_matrix( gmat )
          if (lobatto_matrix_is_needed) lmat = identity_matrix( lmat )
          if (equidist_matrix_is_needed) emat = identity_matrix( emat )
          !
        case (Geom_Pris)
          !
          ! Not sure about prisms here so just use the identity matrix
          !
          if (gauss_matrix_is_needed) gmat = identity_matrix( gmat )
          if (lobatto_matrix_is_needed) lmat = identity_matrix( lmat )
          if (equidist_matrix_is_needed) emat = identity_matrix( emat )
          !
        case (Geom_Hexa)
          !
          solpts => std_elem(Geom_Edge,this_order)%pts
          !
          if (gauss_matrix_is_needed) then
            gmat = InterpToCell_Hexa( solpts(1,:) , gauss_pts )
          end if
          !
          if (lobatto_matrix_is_needed) then
            lmat = InterpToCell_Hexa( solpts(1,:) , lobatto_pts )
          end if
          !
          if (equidist_matrix_is_needed) then
            emat = InterpToCell_Hexa( solpts(1,:) , equidist_pts )
          end if
          !
          if (error_weights_are_needed) then
            this_outerp%err_wts = InterpWeights_Hexa( error_wts )
          end if
          !
        case default
          !
          write (error_message,2)
          call stop_recon(stop_mpi,pname,error_message)
          !
      end select
      !
      ! Now assign the pointers within outerp as needed
      !
      ! If needed, assign the pointer to the matrix for computing the error
      !
      if (iop_solpts /= loc_error_pts) then
        if (loc_error_pts == Legendre_Gauss) then
          this_outerp%error => this_outerp%gauss
        end if
        if (loc_error_pts == Legendre_Gauss_Lobatto) then
          this_outerp%error => this_outerp%lobatto
        end if
      end if
      !
      ! If needed, assign the pointer to the matrix used
      ! to compute the initial conditions
      !
      if (iop_solpts /= loc_init_pts) then
        if (loc_init_pts == Legendre_Gauss) then
          this_outerp%init => this_outerp%gauss
        end if
        if (loc_init_pts == Legendre_Gauss_Lobatto) then
          this_outerp%init => this_outerp%lobatto
        end if
      end if
      !
      ! If needed, assign the pointer to the matrix used
      ! to output the solution for post-processing
      !
      if (iop_solpts /= loc_output_pts) then
        if (loc_output_pts == Legendre_Gauss) then
          this_outerp%output => this_outerp%gauss
        end if
        if (loc_output_pts == Legendre_Gauss_Lobatto) then
          this_outerp%output => this_outerp%lobatto
        end if
        if (loc_output_pts == Equi_Distant_Points) then
          this_outerp%output => this_outerp%equidist
        end if
      end if
      !
      ! Deallocate any of the local arrays that were allocated
      !
      if (allocated(gauss_pts)) then
        deallocate ( gauss_pts , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"gauss_pts",2,ierr,error_message, &
                         skip_alloc_pause)
      end if
      if (allocated(gauss_wts)) then
        deallocate ( gauss_wts , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"gauss_wts",2,ierr,error_message, &
                         skip_alloc_pause)
      end if
      if (allocated(lobatto_pts)) then
        deallocate ( lobatto_pts , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"lobatto_pts",2,ierr,error_message, &
                         skip_alloc_pause)
      end if
      if (allocated(lobatto_wts)) then
        deallocate ( lobatto_wts , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"lobatto_wts",2,ierr,error_message, &
                         skip_alloc_pause)
      end if
      if (allocated(equidist_pts)) then
        deallocate ( equidist_pts , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"equidist_pts",2,ierr,error_message, &
                         skip_alloc_pause)
      end if
      if (allocated(error_wts)) then
        deallocate ( error_wts , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"error_wts",2,ierr,error_message, &
                         skip_alloc_pause)
      end if
      !
      ! Remove all associations of the local pointers before continuing
      !
      if (associated(this_outerp)) this_outerp => null()
      if (associated(gmat       )) gmat        => null()
      if (associated(lmat       )) lmat        => null()
      if (associated(emat       )) emat        => null()
      if (associated(solpts     )) solpts      => null()
      !
    end do order_loop
    !
  end do geom_loop
  !
  do this_geom = gmin,gmax
    do this_order = omin,omax
      if (allocated(outerp(this_geom,this_order)%gauss)) then
        call outerp(this_geom,this_order)%gauss%check_for_sparsity
      end if
      if (allocated(outerp(this_geom,this_order)%lobatto)) then
        call outerp(this_geom,this_order)%lobatto%check_for_sparsity
      end if
      if (allocated(outerp(this_geom,this_order)%equidist)) then
        call outerp(this_geom,this_order)%equidist%check_for_sparsity
      end if
    end do
  end do
  !
  1 format ("outerp(",a,",",i0,")%",a)
  2 format (" A grid cell of an unknown geometry type was found!")
  !
end subroutine init_outerpolation_matrices
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
pure function InterpToCell_Quad(xi,xi_interp) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xi
  real(wp), dimension(:), intent(in) :: xi_interp
  !
  !.. Function Return Value ..
  real(wp), dimension(1:size(xi)**2,1:size(xi)**2) :: return_value
  !
  !.. Local Scalars ..
  integer :: m,mi,mj,n,ni,nj,npts
  !
  !.. Local Arrays ..
  real(qp), dimension(1:2) :: xy_qp
  real(qp), dimension(1:size(xi)) :: xi_qp
  real(qp), dimension(1:size(xi)) :: pts_qp
  real(qp), dimension(1:size(xi)**2,1:size(xi)**2) :: interp_mat_qp
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
  do mj = 1,npts
    do mi = 1,npts
      !
      m = (mj-1)*npts + mi
      !
      do nj = 1,npts
        !
        xy_qp(2) = pts_qp(nj)
        !
        do ni = 1,npts
          !
          xy_qp(1) = pts_qp(ni)
          !
          n = (nj-1)*npts + ni
          !
          interp_mat_qp(n,m) = eval_LagrangePoly([mi,mj],xy_qp,xi_qp)
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
  return_value = chop( interp_mat_qp )
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
  logical, parameter :: use_condensed_loop = true
 !logical, parameter :: use_condensed_loop = fals
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
pure function InterpToCell_Hexa(xi,xi_interp) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xi
  real(wp), dimension(:), intent(in) :: xi_interp
  !
  !.. Function Return Value ..
  real(wp), dimension(1:size(xi)**3,1:size(xi)**3) :: return_value
  !
  !.. Local Scalars ..
  integer :: m,mi,mj,mk,n,ni,nj,nk,npts
  !
  !.. Local Arrays ..
  real(qp), dimension(1:3) :: xy_qp
  real(qp), dimension(1:size(xi)) :: xi_qp
  real(qp), dimension(1:size(xi)) :: pts_qp
  real(qp), dimension(1:size(xi)**3,1:size(xi)**3) :: interp_mat_qp
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
  do mk = 1,npts
    do mj = 1,npts
      do mi = 1,npts
        !
        m = (mk-1)*npts*npts + (mj-1)*npts + mi
        !
        do nk = 1,npts
          !
          xy_qp(3) = pts_qp(nk)
          !
          do nj = 1,npts
            !
            xy_qp(2) = pts_qp(nj)
            !
            do ni = 1,npts
              !
              xy_qp(1) = pts_qp(ni)
              !
              n = (nk-1)*npts*npts + (nj-1)*npts + ni
              !
              interp_mat_qp(n,m) = eval_LagrangePoly([mi,mj,mk],xy_qp,xi_qp)
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
        if (allocated(outerp(i,j)%gauss)) then
          !
          call array%set( storage_size(outerp(i,j)%gauss,kind=inttype) )
          !
          write (array_name,4) "outerp",Geom_Name(i),j,"gauss"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          array = outerp(i,j)%gauss%memory_usage()
          !
          write (array_name,4) "outerp",Geom_Name(i),j,"gauss%mat"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
        end if
        !
        if (allocated(outerp(i,j)%lobatto)) then
          !
          call array%set( storage_size(outerp(i,j)%lobatto,kind=inttype) )
          !
          write (array_name,4) "outerp",Geom_Name(i),j,"lobatto"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          array = outerp(i,j)%lobatto%memory_usage()
          !
          write (array_name,4) "outerp",Geom_Name(i),j,"lobatto%mat"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
        end if
        !
        if (allocated(outerp(i,j)%equidist)) then
          !
          call array%set( storage_size(outerp(i,j)%equidist,kind=inttype) )
          !
          write (array_name,4) "outerp",Geom_Name(i),j,"equidist"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          array = outerp(i,j)%equidist%memory_usage()
          !
          write (array_name,4) "outerp",Geom_Name(i),j,"equidist%mat"
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
