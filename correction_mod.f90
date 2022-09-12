module correction_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use matrix_types_mod, only : matrix
  use matrix_types_mod, only : create_matrix_sparsity
  use matrix_types_mod, only : output_sparsity
  use matrix_types_mod, only : output_matrix
  !
  implicit none
  !
  private
  !
  !
  !.. Public Procedures ..
  !
  public :: init_correction_matrices
  public :: correction_memory_usage
  public :: move_correct_to_device
  !
  !
  !
  ! correct : correction matrix for all combinations of cell geometry and
  !           and cell order that are possible in the current simulation
  !
  !     For cell 'n' :
  !
  !         The correction matrix (giving the contribution to each solution
  !         point within cell 'n' from the interface flux point 'k' that is
  !         currently of interest) is given by the expression
  !
  !         correct(this_geom,this_order)%mat(:,k)
  !
  !     Where :
  !
  !        this_geom = cell(n)%geom
  !        this_order = cell(n)%order
  !        k = the flx-pt on a face of cell 'n' that is currently of interest
  !
  type(matrix), public, save, allocatable :: correct(:,:)
  !
  interface face_interpolation_matrix
    module procedure face_interpolation_matrix, &
                     face_interpolation_matrix_from_qdtr
  end interface face_interpolation_matrix
  !
contains
!
!###############################################################################
!
subroutine init_correction_matrices()
  !
  !.. Use Statements ..
  use order_mod,        only : n_min_geom,n_max_geom
  use order_mod,        only : n_min_order,n_max_order
  use element_type_mod, only : geom_is_used
  !
  !.. Local Scalars ..
  integer :: gmin,gmax,omin,omax
  integer :: this_geom,this_order,ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_correction_matrices"
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
  ! Allocate the correct array
  !
  allocate ( correct(gmin:gmax,omin:omax) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"correct",1,__LINE__,__FILE__,ierr,error_message)
  !
  geom_loop: do this_geom = gmin,gmax
    !
    if (.not. geom_is_used(this_geom)) cycle geom_loop
    if (all(this_geom /= Geom_Valid)) cycle geom_loop
    if (this_geom == Geom_Node) cycle geom_loop
    !
    order_loop: do this_order = omin,omax
      !
      ! Compute the correction matrix for this geometry/order combo
      !
      correct(this_geom,this_order)%mat = & ! F2003 AUTO-REALLOC
                             General_Correction_Matrix( this_geom , this_order )
      !
      error_message = "Generalized Correction"
      call output_matrix(correct(this_geom,this_order), &
                         error_message,this_geom,this_order, &
                         output_transpose=.true.)
      !
#ifdef SPECIAL_OUTPUT_SPARSITY
      write (error_message,1) Geom_Short_Name(this_geom),this_order
      call output_sparsity( correct(this_geom,this_order) , error_message )
#endif
      !
    end do order_loop
    !
  end do geom_loop
  !
  ! Check for sparsity within each of the correction matrices
  !
  call create_matrix_sparsity( correct )
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (a4,".correct.p",i0,".tec.dat")
  !
end subroutine init_correction_matrices
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function CorrectionMatrix_Edge(xi,eval_deriv) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : correction_function
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xi
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: eval_deriv
  !
  !.. Function Result ..
  real(wp), dimension(1:size(xi)) :: return_value
  !
  !.. Local Scalars ..
  integer  :: n,k,j,npts,ndeg
  logical(lk) :: compute_derivative
  real(qp) :: cn,cnm1,cnm2
  !
  !.. Local Arrays ..
  real(qp), dimension(1:size(xi)) :: x
  real(qp), dimension(1:size(xi)) :: cmat
  real(qp), dimension(0:size(xi)) :: coef_c
  real(qp), dimension(0:size(xi)) :: coef_dc
  !
  real(qp), dimension( 0:size(xi), 1:size(xi)) :: xexp
  real(qp), dimension(-1:size(xi),-1:size(xi)) :: coef_Le
  real(qp), dimension( 0:size(xi), 0:size(xi)) :: coef_RR
  !
continue
  !
  ! #########################################################################
  !
  !       THE LOOP INDICES IN THIS FUNCTION ARE USED FOR THE FOLLOWING
  !       SPECIFIC TASKS
  !
  !    n : index used for looping through polynomial degrees / exponents
  !    k : index used for looping over the solution points
  !    j : index used for looping through the coefficients of a polynomial
  !
  ! #########################################################################
  !
  npts = size(xi) ! number of solution points
  ndeg = size(xi) ! polynomial degree of the correction function
  !
  compute_derivative = true
  if (present(eval_deriv)) then
    compute_derivative = eval_deriv
  end if
  !
  ! Initialize the function result
  !
  return_value(:) = zero
  !
  ! Initialize the local arrays
  !
  x(:) = real( xi(:) , kind=qp )
  !
  cmat(:)      = qzero
  coef_c(:)    = qzero
  coef_dc(:)   = qzero
  coef_Le(:,:) = qzero
  coef_RR(:,:) = qzero
  !
  ! Initialize specific coefficients of Legendre polynomial
  !
  coef_Le(0,0) = qone
  coef_Le(1,1) = qone
  !
  ! Get the value of each nodal coordinate raised to all possible powers
  !   using xk = x(k), xexp(0:neg,k) contains the values of
  !        [ xk^0, xk^1, xk^2, ..., xk^(ndeg-1), xk^ndeg]
  !
  xexp(0,:) = qone
  do k = 1,npts
    do n = 1,ndeg
      xexp(n,k) = xexp(n-1,k) * x(k)
    end do
  end do
  !
  ! Get the coefficients for every Legendre polynomial up to degree ndeg
  !
  do n = 2,ndeg
    !
    cnm1 = real( 2*n-1 , kind=qp ) / real( n , kind=qp )
    cnm2 = real(   n-1 , kind=qp ) / real( n , kind=qp )
    !
    ! Loop to get the coefficients for Legendre polynomial of degree n
    !
    do j = 0,n
      coef_Le(j,n) = cnm1*coef_Le(j-1,n-1) - cnm2*coef_Le(j,n-2)
    end do
    !
  end do
  !
  ! Get the coefficients for every right Radau polynomial up to degree ndeg
  !
  do n = 1,ndeg
    !
    cn = qhalf * real( (-1)**n , kind=qp )
    !
    ! Loop to get the coefficients for right Radau polynomial of degree n
    !
    do j = 0,n
      coef_RR(j,n) = cn * (coef_Le(j,n) - coef_Le(j,n-1))
    end do
    !
  end do
  !
  ! Get the coefficients for the correction function
  !
  if (correction_function == Correction_g2) then
    !
    ! Use the g2 correction function
    !
    cn   = real( ndeg-1 , kind=qp ) / real( 2*ndeg-1 , kind=qp )
    cnm1 = real( ndeg   , kind=qp ) / real( 2*ndeg-1 , kind=qp )
    !
    do j = 0,ndeg
      coef_c(j) = cn*coef_RR(j,ndeg) + cnm1*coef_RR(j,ndeg-1)
    end do
    !
  else if (correction_function == Correction_gGauss) then
    !
    ! Use the gGauss correction function
    !
    cn   = real( ndeg-1 , kind=qp ) / real( 2*ndeg-1 , kind=qp )
    cnm1 = real( ndeg   , kind=qp ) / real( 2*ndeg-1 , kind=qp )
    !
    do j = 0,ndeg
      coef_c(j) = cn*coef_RR(j,ndeg) + cnm1*coef_RR(j,ndeg-1)
    end do
    !
  else
    !
    ! Use the right Radau polynomial for the correction function
    !
    do j = 0,ndeg
      coef_c(j) = coef_RR(j,ndeg)
    end do
    !
  end if
  !
  if (compute_derivative) then
    !
    ! Compute the coefficients for the derivative of the correction function
    !
    do j = 1,ndeg
      coef_dc(j-1) = real( j , kind=qp ) * coef_c(j)
    end do
    !
    ! Finally, get the value for the derivative of the
    ! correction function at each 1D solution point
    !
    do k = 1,npts
      cmat(k) = -sum( xexp(0:ndeg-1,k) * coef_dc(0:ndeg-1) )
    end do
    !
  else
    !
    ! Finally get the value of the correction function
    ! at each 1D solution point
    !
    do k = 1,npts
      cmat(k) = sum( xexp(0:ndeg,k) * coef_c(0:ndeg) )
    end do
    !
  end if
  !
  ! Use the chop function to copy cmat into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( cmat )
  !
end function CorrectionMatrix_Edge
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function General_Correction_Matrix(this_geom,this_order) result(return_value)
  !
  !.. Use Statements ..
  use order_mod,        only : geom_solpts,geom_flxpts
  use element_type_mod, only : element_t
  use element_type_mod, only : std_elem,face_elem
  use quadrature_mod,   only : get_edge_quadrature
  use quadrature_mod,   only : get_tria_quadrature
  use quadrature_mod,   only : get_quad_quadrature
  use quadrature_mod,   only : get_ref_elem_node_idx_for_face
  use quadrature_mod,   only : map_face_pts
  use vandermonde_mod,  only : vand
  use vandermonde_mod,  only : General_Basis_Matrix
  !
  !.. Formal Arguments ..
  integer, intent(in) :: this_geom
  integer, intent(in) :: this_order
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: n1,n2,nsp,nfp,ndim
  integer :: this_face,face_geom
  !
  logical :: needsEdge,needsTria,needsQuad
  !
  type(element_t) :: edge_qdtr
  type(element_t) :: tria_qdtr
  type(element_t) :: quad_qdtr
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: ip(:)
  real(wp), allocatable :: rs(:,:)
  real(wp), allocatable :: node_xyz(:,:)
  real(wp), allocatable :: face_pts(:,:)
  real(wp), allocatable :: bmat(:,:)
  real(wp), allocatable :: smat(:,:)
  real(wp), allocatable :: bwimat(:,:)
  real(wp), allocatable :: cmat(:,:)
  !
  real(wp), allocatable :: imat(:,:)
  real(wp), allocatable :: wimat(:,:)
  real(wp), allocatable :: wimatEdge(:,:)
  real(wp), allocatable :: wimatTria(:,:)
  real(wp), allocatable :: wimatQuad(:,:)
  !
  type(matrix) :: tmp
  integer :: i,j
  !
continue
  !
  nsp = geom_solpts(this_geom,this_order)
  nfp = geom_flxpts(this_geom,this_order)
  !
  ! Special logic if this is for an edge element
  !
  if (this_geom == Geom_Edge) then
    !
    return_value = init_array(nsp,nfp,zero) ! F2003 AUTO-REALLOC
    !
    edge_qdtr = std_elem(Geom_Edge,this_order)
    !
    return_value(:,1) = CorrectionMatrix_Edge( edge_qdtr%pts(1,:) )
    return_value(:,2) = return_value(nsp:1:-1,1)
    !
    return
    !
  end if
 !!
 !needsEdge = any( this_geom == [Geom_Tria,Geom_Quad] )
 !needsTria = any( this_geom == [Geom_Tetr,Geom_Pyra,Geom_Pris] )
 !needsQuad = any( this_geom == [Geom_Pyra,Geom_Pris,Geom_Hexa] )
  !
  needsEdge = any( cell_face_geom(:,this_geom) == Geom_Edge )
  needsTria = any( cell_face_geom(:,this_geom) == Geom_Tria )
  needsQuad = any( cell_face_geom(:,this_geom) == Geom_Quad )
  !
  ! Compute the wimatEdge matrix if this element has edge faces.
  !
  if (needsEdge) then
    !
    ! Get the coordinates and weights for the edge-face
    ! quadrature/integration points
    !
    edge_qdtr = get_edge_quadrature(this_order,Legendre_Gauss)
    !
    ! Compute the interpolation matrix from the face flux points
    ! to the face quadrature/integration points
    !
    imat = face_interpolation_matrix(Geom_Edge, & ! F2003 AUTO-REALLOC
                                     face_elem(Geom_Edge,this_order), &
                                     edge_qdtr)
    !
    ! Compute diag(wts)*imat
    !
    ! F2003 AUTO-REALLOC
    wimatEdge = spread(edge_qdtr%wts,dim=1,ncopies=size(imat,dim=1))*imat
    !
  end if
  !
  ! Compute the wimatTria matrix if this element has triangle faces.
  !
  if (needsTria) then
    !
    ! Get the coordinates and weights for the triangle-face
    ! quadrature/integration points
    !
    tria_qdtr = get_tria_quadrature(this_order,WilliamsShunn_TriPoints)
    !
    ! Compute the interpolation matrix from the face flux points
    ! to the face quadrature/integration points
    !
    imat = face_interpolation_matrix(Geom_Tria, & ! F2003 AUTO-REALLOC
                                     face_elem(Geom_Tria,this_order), &
                                     tria_qdtr)
    !
    ! Compute diag(wts)*imat
    !
    ! F2003 AUTO-REALLOC
    wimatTria = spread(tria_qdtr%wts,dim=1,ncopies=size(imat,dim=1))*imat
    !
  end if
  !
  ! Compute the wimatQuad matrix if this element has quad faces.
  !
  if (needsQuad) then
    !
    ! Get the coordinates and weights for the quad-face
    ! quadrature/integration points
    !
    quad_qdtr = get_quad_quadrature(this_order,Legendre_Gauss)
    !
    ! Compute the interpolation matrix from the face flux points
    ! to the face quadrature/integration points
    !
    imat = face_interpolation_matrix(Geom_Quad, & ! F2003 AUTO-REALLOC
                                     face_elem(Geom_Quad,this_order), &
                                     quad_qdtr)
    !
    ! Compute diag(wts)*imat
    !
    ! F2003 AUTO-REALLOC
    wimatQuad = spread(quad_qdtr%wts,dim=1,ncopies=size(imat,dim=1))*imat
    !
  end if
  !
  ! Now we loop over the faces to compute the sections of bwimat
  !
  bwimat = init_array(nfp,nsp,zero) ! F2003 AUTO-REALLOC
  !
  n2 = 0
  ndim = geom_dimen(this_geom)
  !
  do this_face = 1,geom_faces(this_geom)
    !
    ! Get the geometry of the current face
    !
    face_geom = cell_face_geom(this_face,this_geom)
    !
    ! Load up rs and wimat depending what geometry the current face is
    !
    if (face_geom == Geom_Edge) then
      rs = edge_qdtr%pts ! F2003 AUTO-REALLOC
      wimat = wimatEdge ! F2003 AUTO-REALLOC
    else if (face_geom == Geom_Tria) then
      rs = tria_qdtr%pts ! F2003 AUTO-REALLOC
      wimat = wimatTria ! F2003 AUTO-REALLOC
    else if (face_geom == Geom_Quad) then
      rs = quad_qdtr%pts ! F2003 AUTO-REALLOC
      wimat = wimatQuad ! F2003 AUTO-REALLOC
    end if
    !
    ! Get the extents of the face points for the current face within the
    ! range of all face points across all faces of the current cell
    !
    n1 = n2 + 1
    n2 = n2 + size(rs,dim=2)
    !
    ! Get the node indices of the current face for the reference element
    !
    ip = get_ref_elem_node_idx_for_face(this_geom,this_face) ! F2003 AUTO-REALLOC
    !
    ! Get the coordinates of the nodes for the current face
    !
    node_xyz = ref_elem_nodes_xyz(1:ndim,ip(:),this_geom) ! F2003 AUTO-REALLOC
    !
    ! Map the face integration points onto the
    ! current face of the reference element
    !
    face_pts = map_face_pts(node_xyz,rs) ! F2003 AUTO-REALLOC
 !  !
 !  write (error_message,1) "face_pts",this_face
 !  tmp%mat = transpose(face_pts) ! F2003 AUTO-REALLOC
 !  call output_matrix(tmp,error_message,this_geom,this_order)
    !
    ! Compute a basis matrix containing the values of the cell orthonormal
    ! basis functions evaluated at the face integration points
    !
    bmat = General_Basis_Matrix(this_geom,this_order,face_pts) ! F2003 AUTO-REALLOC
 !  !
 !  write (error_message,1) "bmat",this_face
 !  tmp%mat = transpose(bmat) ! F2003 AUTO-REALLOC
 !  call output_matrix(tmp,error_message,this_geom,this_order)
    !
    smat = matmul( transpose(bmat) , transpose(wimat) ) ! F2003 AUTO-REALLOC
 !  !
 !  write (error_message,1) "smat",this_face
 !  tmp%mat = transpose(smat) ! F2003 AUTO-REALLOC
 !  call output_matrix(tmp,error_message,this_geom,this_order)
    !
    bwimat(n1:n2,1:nsp) = transpose(smat(:,:))
    !
  end do
 !1 format (a," (Face #",i1,")")
 !!
 !error_message = "bwimat"
 !tmp%mat = bwimat ! F2003 AUTO-REALLOC
 !call output_matrix(tmp,error_message,this_geom,this_order)
 !!
 !error_message = "vand"
 !call output_matrix(vand(this_geom,this_order)%modal2nodal, &
 !                   error_message,this_geom,this_order)
  !
  ! F2003 AUTO-REALLOC
  cmat = matmul( bwimat , transpose(vand(this_geom,this_order)%modal2nodal%mat) )
 !!
 !error_message = "cmat"
 !tmp%mat = cmat ! F2003 AUTO-REALLOC
 !call output_matrix(tmp,error_message,this_geom,this_order)
  !
  return_value = chop( real( transpose(cmat) , kind=qp ) ) ! F2003 AUTO-REALLOC
  !
end function General_Correction_Matrix
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function face_interpolation_matrix(face_geom,from_order,from_loc, &
                                   to_order,to_loc) result(return_value)
  !
  ! This function is needed because we cant use General_Interp_Matrix here in
  ! the situation where the face flux points use different quadrature points
  ! from the interior solution points of the cell. The matrices stored in
  ! vand(cell_geom,cell_order) are based on the input loc_solution_pts and the
  ! Vandermonde matrix needed here needs to be based of of loc_flux_pts.
  !
  !.. Use Statements ..
  use order_mod,        only : geom_solpts,geom_flxpts
  use element_type_mod, only : element_t
  use quadrature_mod,   only : get_edge_quadrature
  use quadrature_mod,   only : get_tria_quadrature
  use quadrature_mod,   only : get_quad_quadrature
  use vandermonde_mod,  only : vandermonde_inverse
  use vandermonde_mod,  only : General_Basis_Matrix
  !
  !.. Formal Arguments ..
  integer, intent(in) :: face_geom
  integer, intent(in) :: from_order
  integer, intent(in) :: from_loc
  integer, intent(in) :: to_order
  integer, intent(in) :: to_loc
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: npf,npt
  !
  type(element_t) :: from_qdtr
  type(element_t) :: to_qdtr
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: vand(:,:)
  real(wp), allocatable :: vinv(:,:)
  real(wp), allocatable :: bmat(:,:)
  !
continue
  !
  npf = geom_solpts(face_geom,from_order)
  npt = geom_solpts(face_geom,to_order)
  !
  return_value = init_array(npf,npt,zero) ! F2003 AUTO-REALLOC
  !
  select case (face_geom)
    !
    case (Geom_Edge)
      !
      from_qdtr = get_edge_quadrature(from_order,from_loc)
      to_qdtr = get_edge_quadrature(to_order,to_loc)
      !
    case (Geom_Tria)
      !
      from_qdtr = get_tria_quadrature(from_order,from_loc)
      to_qdtr = get_tria_quadrature(to_order,to_loc)
      !
    case (Geom_Quad)
      !
      from_qdtr = get_quad_quadrature(from_order,from_loc)
      to_qdtr = get_quad_quadrature(to_order,to_loc)
      !
    case default
      !
      return
      !
  end select
  !
  ! F2003 AUTO-REALLOC
  vand = General_Basis_Matrix(face_geom,from_order,nodal_points=from_qdtr%pts)
  vinv = Vandermonde_Inverse(vand) ! F2003 AUTO-REALLOC
  !
  ! F2003 AUTO-REALLOC
  bmat = General_Basis_Matrix(face_geom,from_order,to_qdtr%pts)
  !
  return_value(1:npf,1:npt) = matmul( transpose(vinv) , transpose(bmat) )
  !
end function face_interpolation_matrix
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function face_interpolation_matrix_from_qdtr(face_geom,from_qdtr,to_qdtr) &
                                      result(return_value)
  !
  ! This function is needed because we cant use General_Interp_Matrix here in
  ! the situation where the face flux points use different quadrature points
  ! from the interior solution points of the cell. The matrices stored in
  ! vand(cell_geom,cell_order) are based on the input loc_solution_pts and the
  ! Vandermonde matrix needed here needs to be based of of loc_flux_pts.
  !
  !.. Use Statements ..
  use element_type_mod, only : element_t
  use vandermonde_mod,  only : vandermonde_inverse
  use vandermonde_mod,  only : General_Basis_Matrix
  !
  !.. Formal Arguments ..
  integer,         intent(in) :: face_geom
  type(element_t), intent(in) :: from_qdtr
  type(element_t), intent(in) :: to_qdtr
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: npf,npt,from_order
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: vand(:,:)
  real(wp), allocatable :: vinv(:,:)
  real(wp), allocatable :: bmat(:,:)
  !
continue
  !
  npf = size(from_qdtr%wts)
  npt = size(to_qdtr%wts)
  !
  from_order = solpts2order(face_geom,npf)
  !
  return_value = init_array(npf,npt,zero) ! F2003 AUTO-REALLOC
  !
  ! F2003 AUTO-REALLOC
  vand = General_Basis_Matrix(face_geom,from_order,xyz=from_qdtr%pts)
  vinv = Vandermonde_Inverse(vand) ! F2003 AUTO-REALLOC
  !
  ! F2003 AUTO-REALLOC
  bmat = General_Basis_Matrix(face_geom,from_order,to_qdtr%pts)
  !
  return_value(1:npf,1:npt) = matmul( transpose(vinv) , transpose(bmat) )
  !
end function face_interpolation_matrix_from_qdtr
!
!###############################################################################
!
subroutine correction_memory_usage(iunit)
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
  integer :: iou,i,j
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
  ! Memory usage for the correct array of type matrix
  !
  if (allocated(correct)) then
    !
    call dt%set( storage_size(correct,kind=inttype) , &
                         size(correct,kind=inttype) )
    !
    do j = lbound(correct,dim=2),ubound(correct,dim=2)
      do i = lbound(correct,dim=1),ubound(correct,dim=1)
        !
        array = matrix_memory_usage( correct(i,j) )
        !
        write (array_name,4) "correct",Geom_Name(i),j
        write (iou,2) array_name,array%get_units()
        !
        call dt%add(array)
        !
      end do
    end do
    !
    call total%add(dt)
    write (array_name,5) "correct"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "correct"
    !
  end if
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module correction_mod")
  2 format (" Array: ",a,dt)
  3 format (" Total memory = ",f8.2,1x,a)
  4 format (a,"(",a,",",i0,")%mat")
  5 format ("'",a,"'")
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  !
#endif
end subroutine correction_memory_usage
!
!###############################################################################
!
subroutine move_correct_to_device()
  !
  !.. Local Scalars ..
  integer :: i,i1,i2,j,j1,j2,m1,m2,n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "move_correct_to_device"
  !
continue
  !
  if (allocated(correct)) then
    !$ACC ENTER DATA COPYIN(correct)
    i1 = lbound(correct,dim=1)
    i2 = ubound(correct,dim=1)
    j1 = lbound(correct,dim=2)
    j2 = ubound(correct,dim=2)
    !$ACC ENTER DATA COPYIN(correct(i1:i2,j1:j2))
    do j = j1,j2
      do i = i1,i2
        !
        if (allocated(correct(i,j)%mat)) then
          !$ACC ENTER DATA COPYIN(correct(i,j)%mat)
          m1 = lbound(correct(i,j)%mat,dim=1)
          m2 = ubound(correct(i,j)%mat,dim=1)
          n1 = lbound(correct(i,j)%mat,dim=2)
          n2 = ubound(correct(i,j)%mat,dim=2)
          !$ACC ENTER DATA COPYIN(correct(i,j)%mat(m1:m2,n1:n2))
        end if
        !
        if (allocated(correct(i,j)%row_ptr)) then
          !$ACC ENTER DATA COPYIN(correct(i,j)%row_ptr)
          m1 = lbound(correct(i,j)%row_ptr,dim=1)
          m2 = ubound(correct(i,j)%row_ptr,dim=1)
          !$ACC ENTER DATA COPYIN(correct(i,j)%row_ptr(m1:m2))
        end if
        !
        if (allocated(correct(i,j)%row_idx)) then
          !$ACC ENTER DATA COPYIN(correct(i,j)%row_idx)
          m1 = lbound(correct(i,j)%row_idx,dim=1)
          m2 = ubound(correct(i,j)%row_idx,dim=1)
          !$ACC ENTER DATA COPYIN(correct(i,j)%row_idx(m1:m2))
        end if
        !
        if (allocated(correct(i,j)%row_val)) then
          !$ACC ENTER DATA COPYIN(correct(i,j)%row_val)
          m1 = lbound(correct(i,j)%row_val,dim=1)
          m2 = ubound(correct(i,j)%row_val,dim=1)
          !$ACC ENTER DATA COPYIN(correct(i,j)%row_val(m1:m2))
        end if
        !
        if (allocated(correct(i,j)%col_ptr)) then
          !$ACC ENTER DATA COPYIN(correct(i,j)%col_ptr)
          m1 = lbound(correct(i,j)%col_ptr,dim=1)
          m2 = ubound(correct(i,j)%col_ptr,dim=1)
          !$ACC ENTER DATA COPYIN(correct(i,j)%col_ptr(m1:m2))
        end if
        !
        if (allocated(correct(i,j)%col_idx)) then
          !$ACC ENTER DATA COPYIN(correct(i,j)%col_idx)
          m1 = lbound(correct(i,j)%col_idx,dim=1)
          m2 = ubound(correct(i,j)%col_idx,dim=1)
          !$ACC ENTER DATA COPYIN(correct(i,j)%col_idx(m1:m2))
        end if
        !
        if (allocated(correct(i,j)%col_val)) then
          !$ACC ENTER DATA COPYIN(correct(i,j)%col_val)
          m1 = lbound(correct(i,j)%col_val,dim=1)
          m2 = ubound(correct(i,j)%col_val,dim=1)
          !$ACC ENTER DATA COPYIN(correct(i,j)%col_val(m1:m2))
        end if
        !
      end do
    end do
  end if
  !
end subroutine move_correct_to_device
!
!###############################################################################
!
include "Functions/solpts2order.f90"
!
!###############################################################################
!
end module correction_mod
!
!###############################################################################
!
