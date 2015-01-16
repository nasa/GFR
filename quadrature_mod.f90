module quadrature_mod
  !
  use module_kind_types
  !
  implicit none
  !
  private
  !
  public :: find_valid_geometries
  public :: init_geom_quadrature_rules
  public :: init_face_quadrature_rules
  !
  ! element_type : derived type that stores the location of the solution points
  !                and their respective quadrature weights within the standard
  !                region of a given cell geometry.
  !
  type, public :: element_type
    real(wp), allocatable :: pts(:,:)
    real(wp), allocatable :: wts(:)
  end type element_type
  !
  ! std_elem : Array giving the location of the solution points and their
  !            respective quadrature weights for all combinations of cell
  !            geometry and cell order that are possible in the current
  !            simulation
  !
  !  For a given elemental geometry, this_geom,
  !                   and its order, this_order :
  !
  !  std_elem(this_geom,this_order)%pts(l,k) : The location of the solution
  !                                            point 'k' in the direction 'l' of
  !                                            the local coordinate system
  !                                            within the standard region of the
  !                                            current element geometry.
  !
  !  std_elem(this_geom,this_order)%wts(k) : The quadrature weight of the solution
  !                                          point 'k' within the standard region
  !                                          of the current element geometry.
  !
  type(element_type), public, save, target, allocatable :: std_elem(:,:)
  !
  ! face_elem : Array giving the location of the flux points and their
  !             respective quadrature weights for all combinations of cell
  !             geometry and cell order that are possible in the current
  !             simulation
  !
  !  For a given face geometry, this_geom,
  !              and its order, this_order :
  !
  !  face_elem(this_geom,this_order)%pts(l,k) : The location of the flux
  !                                             point 'k' in the direction 'l' of
  !                                             the local coordinate system
  !                                             within the standard region of the
  !                                             current face geometry.
  !
  !  face_elem(this_geom,this_order)%wts(k) : The quadrature weight of the flux
  !                                           point 'k' within the standard region
  !                                           of the current face geometry.
  !
  type(element_type), public, save, target, allocatable :: face_elem(:,:)
  !
  !
  logical, public, protected, save :: geom_is_used(Geom_Min:Geom_Max)
  !
  interface get_triangle_weights_at_solpts
    module procedure get_triangle_weights_at_solpts_DP, &
                     get_triangle_weights_at_solpts_QP
  end interface get_triangle_weights_at_solpts
  !
contains
!
!###############################################################################
!
subroutine find_valid_geometries()
  !
  !.. Use Statements ..
  use geovar,    only : cell_geom,cell_order
  use order_mod, only : n_min_geom,n_max_geom
  use order_mod, only : geom_solpts,maxpts
  !
  !.. Local Scalars ..
  integer :: nc,this_geom
  !
continue
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
  n_min_geom = minval( [Geom_Min:Geom_Max] , mask=geom_is_used )
  n_max_geom = maxval( [Geom_Min:Geom_Max] , mask=geom_is_used )
  !
  ! While we are here, find the maximum number of points in any cell
  !
  maxpts = 1
  do nc = 1,size(cell_geom)
    maxpts = max( maxpts , geom_solpts(cell_geom(nc),cell_order(nc)) )
  end do
  !
end subroutine find_valid_geometries
!
!###############################################################################
!
subroutine init_geom_quadrature_rules()
  !
  !.. Use Statements ..
  use order_mod, only : geom_solpts
  use order_mod, only : n_min_geom,n_max_geom
  use order_mod, only : n_min_order,n_max_order
  use ovar,      only : iop_solpts
  use ovar,      only : iop_tripts
  !
  !.. Local Scalars ..
  integer :: n,npts,ndim,ierr
  integer :: this_geom,this_order
  character(len=100) :: array_name
  character(len=200) :: error_message
  !
  !.. Local Pointers ..
  type(element_type), pointer :: this_elem
  type(element_type), pointer :: edge_elem
  type(element_type), pointer :: tria_elem
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_geom_quadrature_rules"
  !
continue
  !
  ! Initialize this_elem to disassociated
  !
  this_elem => null()
  edge_elem => null()
  tria_elem => null()
  !
  ! Allocate the std_elem array
  !
  allocate ( std_elem(Geom_Min:Geom_Max,n_min_order:n_max_order) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"std_elem",1,ierr,error_message)
  !
  ! Loop through all the possible orders and get the quadrature rules
  ! for the most basic geometries that are used (Geom_Edge and Geom_Tria).
  ! We need to do this before all other geometry types because the
  ! remaining geometries are dependent on these quadrature rules. 
  ! NOTE: The pts array within the type element_type is only needed for these
  !       two geometries (maybe tetrahedra as well) so they will be the only
  !       ones in which this array will be allocated and created.
  !
  this_geom = Geom_Edge
  !
  basic_geom_loop: do n = 1,2
    !
    basic_order_loop: do this_order = n_min_order,n_max_order
      !
      ! Assign the pointer 'this_elem' to std_elem(this_geom,this_order)
      ! to simplify the remaining code in this loop and make it more legible
      !
      this_elem => std_elem(this_geom,this_order)
      !
      ! Spatial dimension of the current geometry
      ! NOTE: Make sure that ndim is at minimum 1, primarily for Geom_Node
      !
      ndim = max( 1 , geom_dimen(this_geom) )
      !
      ! Number of solution points for the current geometry/order combination
      !
      npts = geom_solpts(this_geom,this_order)
      !
      ! Allocate the wts component of the current geometry/order combination
      !
      allocate ( this_elem%wts(1:npts) , stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"wts"
      call alloc_error(pname,array_name,1,ierr,error_message)
      !
      ! Allocate the pts component of the current geometry/order combination
      !
      allocate ( this_elem%pts(1:ndim,1:npts) , stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"pts"
      call alloc_error(pname,array_name,1,ierr,error_message)
      !
      ! Fill in the solution points and quadrature weights
      ! for an edge with the current order
      !
      if (this_geom == Geom_Edge) then
        this_elem = get_edge_quadrature(this_order,iop_solpts)
      else if (this_geom == Geom_Tria) then
        this_elem = get_tria_quadrature(this_order,iop_tripts)
      end if
      !
      ! Disassociate this_elem before continuing to the next loop combo
      !
      if (associated(this_elem)) this_elem => null()
      !
    end do basic_order_loop
    !
    ! Change this_geom to Geom_Tria so the second pass through
    ! basic_order_loop works on triangle geometry
    !
    this_geom = Geom_Tria
    !
    if (.not. geom_is_used(this_geom)) exit basic_geom_loop
    !
  end do basic_geom_loop
  !
  ! Loop through all the possible combinations of cell geometry and order
  ! and get the quadrature rules for each one
  !
  order_loop: do this_order = n_min_order,n_max_order
    !
    ! Assign pointers to the edge and triangle elements of the std_elem
    ! array in order to simplify the remaining code in this loop and
    ! make it more legible.
    !
    edge_elem => std_elem(Geom_Edge,this_order)
    tria_elem => std_elem(Geom_Tria,this_order)
    !
    geom_loop: do this_geom = n_min_geom,n_max_geom
      !
      ! Cycle to the next geometry if the current one is:
      !   1. not used
      !   2. an edge or triangle since these have already been done
      !   3. not valid
      !
      if (.not. geom_is_used(this_geom)) cycle geom_loop
      if (any(this_geom == [Geom_Edge,Geom_Tria])) cycle geom_loop
      if (all(this_geom /= Geom_Valid)) cycle geom_loop
      !
      !
      ! Assign the pointer 'this_elem' to std_elem(this_geom,this_order)
      ! to simplify the remaining code in this loop and make it more legible
      !
      this_elem => std_elem(this_geom,this_order)
      !
      ! Spatial dimension of the current geometry
      ! NOTE: Make sure that ndim is at minimum 1, primarily for Geom_Node
      !
      ndim = max( 1 , geom_dimen(this_geom) )
      !
      ! Number of solution points for the current geometry/order combination
      !
      npts = geom_solpts(this_geom,this_order)
      !
      ! Allocate the wts component of the current geometry/order combination
      !
      allocate ( this_elem%wts(1:npts) , stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"wts"
      call alloc_error(pname,array_name,1,ierr,error_message)
      !
      ! We dont want to allocate the pts component for this geometry since
      ! it will never be used and will just be a waste of memory.
      !
     !allocate ( this_elem%pts(1:ndim,1:npts) , stat=ierr , errmsg=error_message )
     !write (array_name,1) Geom_Name(this_geom),this_order,"pts"
     !call alloc_error(pname,array_name,1,ierr,error_message)
      !
      ! Fill in the quadrature weights for the current geometry/order combination
      !
      select case (this_geom)
        case (Geom_Node)
          !
          this_elem%wts = one
          !
        case (Geom_Quad)
          !
          this_elem = get_quad_quadrature(edge_elem)
          !
        case (Geom_Tetr)
          !
          write (error_message,3) "Geom_Tetr"
          call stop_recon(stop_mpi,pname,error_message)
          this_elem = get_tetr_quadrature(this_order)
          !
        case (Geom_Pyra)
          !
          write (error_message,3) "Geom_Pyra"
          call stop_recon(stop_mpi,pname,error_message)
          this_elem = get_pyra_quadrature(this_order)
          !
        case (Geom_Pris)
          !
          write (error_message,3) "Geom_Pris"
          call stop_recon(stop_mpi,pname,error_message)
          this_elem = get_pris_quadrature(this_order)
          !
        case (Geom_Hexa)
          !
          this_elem = get_hexa_quadrature(edge_elem)
          !
        case default
          !
          write (error_message,2)
          call stop_recon(stop_mpi,pname,error_message)
          !
      end select
      !
      ! Disassociate the pointer 'this_elem' before
      ! continuing to the next geometry
      !
      if (associated(this_elem)) this_elem => null()
      !
    end do geom_loop
    !
    ! Disassociate the pointers 'edge_elem' and 'tria_elem' before
    ! continuing to the next order
    !
    if (associated(edge_elem)) edge_elem => null()
    if (associated(tria_elem)) tria_elem => null()
    !
  end do order_loop
  !
  ! Format Statements
  !
  1 format ("std_elem(",a,",",i0,")%",a)
  2 format (" A grid cell of an unknown geometry type was found!")
  3 format (" A grid cell of geometry type '",a,"' was encountered while", &
            " trying to compute the quadrature rules for each geometry type.", &
            " This geometry type is currently invalid because it has not", &
            " yet been fully implemented.")
  !
end subroutine init_geom_quadrature_rules
!
!###############################################################################
!
subroutine init_face_quadrature_rules()
  !
  !.. Use Statements ..
  use order_mod, only : geom_solpts
  use order_mod, only : n_min_order,n_max_order
  use ovar,      only : iop_flxpts
  use ovar,      only : iop_tripts
  !
  !.. Local Scalars ..
  integer :: n,npts,ndim,ierr
  integer :: gmin,gmax,this_order
  logical :: there_are_quad_faces
  logical :: there_are_tria_faces
  character(len=100) :: array_name
  character(len=200) :: error_message
  !
  !.. Local Pointers ..
  type(element_type), pointer :: this_elem
  type(element_type), pointer :: edge_elem
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_face_quadrature_rules"
  !
continue
  !
  ! Initialize this_elem to disassociated
  !
  this_elem => null()
  edge_elem => null()
  !
  ! Determine if there are any of the geometries used for this simulation
  ! are 3D to see if we need to find the quadrature rules for 2D faces.
  !
  gmin = Geom_Edge
  gmax = Geom_Edge
  !
  there_are_tria_faces = fals
  there_are_quad_faces = fals
  !
  if (any(geom_is_used(Geom_3D))) then
    gmax = maxval( Geom_2D , mask=geom_is_used(Geom_2D) )
    there_are_tria_faces = geom_is_used(Geom_Tria)
    there_are_quad_faces = geom_is_used(Geom_Quad)
  end if
  !
  ! Allocate the face_elem array
  !
  allocate ( face_elem(gmin:gmax,n_min_order:n_max_order) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"face_elem",1,ierr,error_message)
  !
  ! Loop through all the possible orders and get the edge quadrature rules.
  !
  edge_loop: do this_order = n_min_order,n_max_order
    !
    ! Assign the pointer 'this_elem' to face_elem(Geom_Edge,this_order)
    ! to simplify the remaining code in this loop and make it more legible
    !
    this_elem => face_elem(Geom_Edge,this_order)
    !
    ! Spatial dimension of the current geometry
    ! NOTE: Make sure that ndim is at minimum 1, primarily for Geom_Node
    !
    ndim = max( 1 , geom_dimen(Geom_Edge) )
    !
    ! Number of solution points for the current geometry/order combination
    !
    npts = geom_solpts(Geom_Edge,this_order)
    !
    ! Allocate the wts component of the current geometry/order combination
    !
    allocate ( this_elem%wts(1:npts) , stat=ierr , errmsg=error_message )
    write (array_name,1) Geom_Name(Geom_Edge),this_order,"wts"
    call alloc_error(pname,array_name,1,ierr,error_message)
    !
    ! Allocate the pts component of the current geometry/order combination
    !
    allocate ( this_elem%pts(1:ndim,1:npts) , stat=ierr , errmsg=error_message )
    write (array_name,1) Geom_Name(Geom_Edge),this_order,"pts"
    call alloc_error(pname,array_name,1,ierr,error_message)
    !
    ! Fill in the solution points and quadrature weights
    ! for an edge with the current order
    !
    this_elem = get_edge_quadrature(this_order,iop_flxpts)
    !
    ! Disassociate this_elem before continuing to the next loop combo
    !
    if (associated(this_elem)) this_elem => null()
    !
  end do edge_loop
  !
  ! Loop through all the possible orders and get the triangle quadrature rules
  ! if there are any triangle faces exist for this simulation
  !
  if (there_are_tria_faces) then
    !
    tria_loop: do this_order = n_min_order,n_max_order
      !
      ! Assign the pointer 'this_elem' to face_elem(Geom_Tria,this_order)
      ! to simplify the remaining code in this loop and make it more legible
      !
      this_elem => face_elem(Geom_Tria,this_order)
      !
      ! Spatial dimension of the current geometry
      ! NOTE: Make sure that ndim is at minimum 1, primarily for Geom_Node
      !
      ndim = max( 1 , geom_dimen(Geom_Tria) )
      !
      ! Number of solution points for the current geometry/order combination
      !
      npts = geom_solpts(Geom_Tria,this_order)
      !
      ! Allocate the wts component of the current geometry/order combination
      !
      allocate ( this_elem%wts(1:npts) , stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(Geom_Tria),this_order,"wts"
      call alloc_error(pname,array_name,1,ierr,error_message)
      !
      ! Allocate the pts component of the current geometry/order combination
      !
      allocate ( this_elem%pts(1:ndim,1:npts) , stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(Geom_Tria),this_order,"pts"
      call alloc_error(pname,array_name,1,ierr,error_message)
      !
      ! Fill in the solution points and quadrature weights
      ! for a triangle face with the current order
      !
      this_elem = get_tria_quadrature(this_order,iop_tripts)
      !
      ! Disassociate this_elem before continuing to the next loop combo
      !
      if (associated(this_elem)) this_elem => null()
      !
    end do tria_loop
    !
  end if
  !
  ! Loop through all the possible orders and get the quad quadrature rules
  ! if there are any quad faces exist for this simulation
  !
  if (there_are_quad_faces) then
    !
    quad_loop: do this_order = n_min_order,n_max_order
      !
      ! Assign a pointer to the edge element of the face_elem array in order
      ! to simplify the remaining code in this loop and make it more legible.
      !
      edge_elem => face_elem(Geom_Edge,this_order)
      !
      ! Assign the pointer 'this_elem' to face_elem(Geom_Quad,this_order)
      ! to simplify the remaining code in this loop and make it more legible
      !
      this_elem => face_elem(Geom_Quad,this_order)
      !
      ! Spatial dimension of the current geometry
      ! NOTE: Make sure that ndim is at minimum 1, primarily for Geom_Node
      !
      ndim = max( 1 , geom_dimen(Geom_Quad) )
      !
      ! Number of solution points for the current geometry/order combination
      !
      npts = geom_solpts(Geom_Quad,this_order)
      !
      ! Allocate the wts component of the current geometry/order combination
      !
      allocate ( this_elem%wts(1:npts) , stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(Geom_Quad),this_order,"wts"
      call alloc_error(pname,array_name,1,ierr,error_message)
     !!
     !! Allocate the pts component of the current geometry/order combination
     !!
     !allocate ( this_elem%pts(1:ndim,1:npts) , stat=ierr , errmsg=error_message )
     !write (array_name,1) Geom_Name(Geom_Quad),this_order,"pts"
     !call alloc_error(pname,array_name,1,ierr,error_message)
      !
      ! Fill in the solution points and quadrature weights
      ! for a quad face with the current order
      !
      this_elem = get_quad_quadrature(edge_elem)
      !
      ! Disassociate this_elem before continuing to the next loop combo
      !
      if (associated(this_elem)) this_elem => null()
      if (associated(edge_elem)) edge_elem => null()
      !
    end do quad_loop
    !
  end if
  !
  ! Format Statements
  !
  1 format ("face_elem(",a,",",i0,")%",a)
  !
end subroutine init_face_quadrature_rules
!
!###############################################################################
!
pure function get_edge_quadrature(n,loc) result(return_value)
  !
  use order_mod,      only : geom_solpts
  use polynomial_mod, only : nodes_legendre_gauss
  use polynomial_mod, only : nodes_legendre_gauss_lobatto
  use polynomial_mod, only : weights_legendre_gauss
  use polynomial_mod, only : weights_legendre_gauss_lobatto
  !
  integer, intent(in) :: n
  integer, intent(in) :: loc
  !
  type(element_type) :: return_value
  !
  integer :: ndim,npts,ierr
  !
  real(qp), allocatable :: pts(:)
  real(qp), allocatable :: wts(:)
  real(qp), allocatable :: val_at_pts(:)
  !
continue
  !
  ndim = geom_dimen(Geom_Edge)
  npts = geom_solpts(Geom_Edge,n)
  !
  allocate ( pts(1:npts) , source=qzero , stat=ierr )
  allocate ( wts(1:npts) , source=qzero , stat=ierr )
  allocate ( val_at_pts(1:npts) , source=qzero , stat=ierr )
  !
  select case (loc)
    !
    case (Legendre_Gauss)
      !
      ! Edge quadrature using Legendre Gauss nodes
      !
      call nodes_legendre_gauss(npts,pts,val_at_pts)
      call weights_legendre_gauss(npts,pts,val_at_pts,wts)
      !
    case (Legendre_Gauss_Lobatto)
      !
      ! Edge quadrature using Legendre Gauss Lobatto nodes
      !
      call nodes_legendre_gauss_lobatto(n,pts,val_at_pts)
      call weights_legendre_gauss_lobatto(n,val_at_pts,wts)
      !
    case default
      !
      ! Default to Legendre Gauss points
      !
      call nodes_legendre_gauss(npts,pts,val_at_pts)
      call weights_legendre_gauss(npts,pts,val_at_pts,wts)
      !
  end select
  !
  deallocate ( val_at_pts , stat=ierr )
  !
  allocate ( return_value%pts(1:ndim,1:npts) , source=zero , stat=ierr )
  allocate ( return_value%wts(1:npts) , source=zero , stat=ierr )
  !
  return_value%pts(1,:) = chop( pts )
  return_value%wts = chop( wts )
  !
  deallocate ( pts , stat=ierr )
  deallocate ( wts , stat=ierr )
  !
end function get_edge_quadrature
!
!###############################################################################
!
pure function get_tria_quadrature(n,loc) result(return_value)
  !
  use order_mod,      only : geom_solpts
  use triangle_mod,   only : TriNodes2D_DG
  use triangle_mod,   only : TriNodes2D_HT
  !
  integer, intent(in) :: n
  integer, intent(in) :: loc
  !
  type(element_type) :: return_value
  !
  integer :: ndim,npts,ierr
  !
  real(qp), allocatable :: pts(:,:)
  real(qp), allocatable :: wts(:)
  !
continue
  !
  ndim = geom_dimen(Geom_Tria)
  npts = geom_solpts(Geom_Tria,n)
  !
  allocate ( pts(1:ndim,1:npts) , source=qzero , stat=ierr )
  allocate ( wts(1:npts) , source=qzero , stat=ierr )
  !
  select case (loc)
    !
    case (DG_TriPoints)
      !
      ! Alpha optimized triangle quadrature from Hesthaven and Warburton
      !
      pts = TriNodes2D_DG(n,qzero)
      !
    case (HT_TriPoints)
      !
      ! 1D Lobatto quadrature extrapolated to a barycentric triangle
      ! coordinate system
      !
      pts = TriNodes2D_HT(n,qzero)
      !
    case default
      !
      ! Default to the barycentric Lobatto quadrature
      !
      pts = TriNodes2D_HT(n,qzero)
      !
  end select
  !
  ! The method for computing the triangle quadrature weights is
  ! independent of the method for getting the point locations
  !
  wts = get_triangle_weights_at_solpts(pts)
  !
  allocate ( return_value%pts(1:ndim,1:npts) , source=zero , stat=ierr )
  allocate ( return_value%wts(1:npts) , source=zero , stat=ierr )
  !
  return_value%pts = chop( pts )
  return_value%wts = chop( wts )
  !
  deallocate ( pts , stat=ierr )
  deallocate ( wts , stat=ierr )
  !
end function get_tria_quadrature
!
!###############################################################################
!
pure function get_quad_quadrature(edge_elem) result(return_value)
  !
  use order_mod, only : geom_solpts
  !
  type(element_type), intent(in) :: edge_elem
  !
  type(element_type) :: return_value
  !
  integer :: i,j,l,n,ndim,npts,ierr
  !
 !real(qp), allocatable :: pts(:,:)
  real(qp), allocatable :: wts(:)
  !
continue
  !
  n = size(edge_elem%wts) - 1
  !
 !ndim = geom_dimen(Geom_Quad)
  npts = geom_solpts(Geom_Quad,n)
 !!
 !! Compute the tensor product of the 1D solution points
 !!
 !allocate ( pts(1:ndim,1:npts) , source=qzero , stat=ierr )
 !!
 !l = 0
 !do j = 1,size(edge_elem%pts,dim=2)
 !  do i = 1,size(edge_elem%pts,dim=2)
 !    l = l + 1
 !    pts(1,l) = real( edge_elem%pts(1,i) , kind=qp )
 !    pts(2,l) = real( edge_elem%pts(1,j) , kind=qp )
 !  end do
 !end do
 !!
 !allocate ( return_value%pts(1:ndim,1:npts) , source=zero , stat=ierr )
 !!
 !return_value%pts = chop( pts )
 !!
 !deallocate ( pts , stat=ierr )
  !
  ! Compute the tensor product of the 1D weights
  !
  allocate ( wts(1:npts) , source=qzero , stat=ierr )
  !
  l = 0
  do j = 1,size(edge_elem%wts)
    do i = 1,size(edge_elem%wts)
      l = l + 1
      wts(l) = real( edge_elem%wts(i) , kind=qp ) *  &
               real( edge_elem%wts(j) , kind=qp )
    end do
  end do
  !
  allocate ( return_value%wts(1:npts) , source=zero , stat=ierr )
  !
  return_value%wts = chop( wts )
  !
  deallocate ( wts , stat=ierr )
  !
end function get_quad_quadrature
!
!###############################################################################
!
pure function get_tetr_quadrature(n) result(return_value)
  !
  use order_mod, only : geom_solpts
  !
  integer, intent(in) :: n
  !
  type(element_type) :: return_value
  !
  integer :: npts,ierr
  !
continue
  !
  npts = geom_solpts(Geom_Tetr,n)
  !
  allocate ( return_value%wts(1:npts) , source=zero , stat=ierr )
  !
  return_value%wts = one / real(npts,kind=wp)
  !
end function get_tetr_quadrature
!
!###############################################################################
!
pure function get_pyra_quadrature(n) result(return_value)
  !
  use order_mod, only : geom_solpts
  !
  integer, intent(in) :: n
  !
  type(element_type) :: return_value
  !
  integer :: npts,ierr
  !
continue
  !
  npts = geom_solpts(Geom_Pyra,n)
  !
  allocate ( return_value%wts(1:npts) , source=zero , stat=ierr )
  !
  return_value%wts = one / real(npts,kind=wp)
  !
end function get_pyra_quadrature
!
!###############################################################################
!
pure function get_pris_quadrature(n) result(return_value)
  !
  use order_mod, only : geom_solpts
  !
  integer, intent(in) :: n
  !
  type(element_type) :: return_value
  !
  integer :: npts,ierr
  !
continue
  !
  npts = geom_solpts(Geom_Pris,n)
  !
  allocate ( return_value%wts(1:npts) , source=zero , stat=ierr )
  !
  return_value%wts = one / real(npts,kind=wp)
  !
end function get_pris_quadrature
!
!###############################################################################
!
pure function get_hexa_quadrature(edge_elem) result(return_value)
  !
  use order_mod, only : geom_solpts
  !
  type(element_type), intent(in) :: edge_elem
  !
  type(element_type) :: return_value
  !
  integer :: i,j,k,l,n,ndim,npts,ierr
  !
 !real(qp), allocatable :: pts(:,:)
  real(qp), allocatable :: wts(:)
  !
continue
  !
  n = size(edge_elem%wts) - 1
  !
 !ndim = geom_dimen(Geom_Hexa)
  npts = geom_solpts(Geom_Hexa,n)
 !!
 !! Compute the tensor product of the 1D solution points
 !!
 !allocate ( pts(1:ndim,1:npts) , source=qzero , stat=ierr )
 !!
 !l = 0
 !do k = 1,size(edge_elem%pts,dim=2)
 !  do j = 1,size(edge_elem%pts,dim=2)
 !    do i = 1,size(edge_elem%pts,dim=2)
 !      l = l + 1
 !      pts(1,l) = real( edge_elem%pts(1,i) , kind=qp )
 !      pts(2,l) = real( edge_elem%pts(1,j) , kind=qp )
 !      pts(3,l) = real( edge_elem%pts(1,k) , kind=qp )
 !    end do
 !  end do
 !end do
 !!
 !allocate ( return_value%pts(1:ndim,1:npts) , source=zero , stat=ierr )
 !!
 !return_value%pts = chop( pts )
 !!
 !deallocate ( pts , stat=ierr )
  !
  ! Compute the tensor product of the 1D weights
  !
  allocate ( wts(1:npts) , source=qzero , stat=ierr )
  !
  l = 0
  do k = 1,size(edge_elem%wts)
    do j = 1,size(edge_elem%wts)
      do i = 1,size(edge_elem%wts)
        l = l + 1
        wts(l) = real( edge_elem%wts(i) , kind=qp ) * &
                 real( edge_elem%wts(j) , kind=qp ) * &
                 real( edge_elem%wts(k) , kind=qp )
      end do
    end do
  end do
  !
  allocate ( return_value%wts(1:npts) , source=zero , stat=ierr )
  !
  return_value%wts = chop( wts )
  !
  deallocate ( wts , stat=ierr )
  !
end function get_hexa_quadrature
!
!###############################################################################
!
pure function get_triangle_weights_at_solpts_DP(rs) result(return_value)
  !
  ! NOTE : This seems to fail for large n_order (>15)
  !
  use polynomial_mod, only : nodes_legendre_gauss
  use polynomial_mod, only : weights_legendre_gauss
  !
  real(dp), dimension(:,:), intent(in) :: rs
  !
  real(dp), dimension(1:size(rs,2)) :: return_value
  !
  integer  :: i,j,k,n,m,np
  integer  :: nspts,nfpts
  real(dp) :: g,h,x,y
  !
  integer, parameter :: nquad = 20
  !
  real(dp), dimension(1:nquad)                   :: xi
  real(dp), dimension(1:nquad)                   :: dxi
  real(dp), dimension(1:nquad)                   :: weights
  real(dp), dimension(1:size(rs,2),1:size(rs,2)) :: V
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: half = 0.5_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  ! On entry, r and s need to be defined on the interval [0,1]
  !
  nspts = size(rs,2)
  nfpts = np2nfp(nspts)
  !
  ! Compute the coefficients of the Vandermonde matrix for
  ! the Lagrange polynomials using monomial basis functions
  !
  do n = 1,nspts
    k = 0
    do i = 1,nfpts
      do j = 1,i
        k = k + 1
        V(k,n) = rs(1,n)**(i-j) * rs(2,n)**(j-1)
      end do
    end do
  end do
  !
  ! Now get the inverse of this matrix for use later
  !
  V = invert_matrix(V)
  !
  ! Get the Gauss quadrature points on the interval [-1,1] and
  ! the corresponding weights necessary to perform the integration
  !
  call nodes_legendre_gauss(nquad,xi,dxi)
  call weights_legendre_gauss(nquad,xi,dxi,weights)
  !
  ! Reset the quadrature points and weights to the interval [0,1]
  !
  xi(:) = half*(one+xi(:))
  weights(:) = half*weights(:)
  !
  ! Initialize the triangle weights to zero
  !
  return_value(:) = zero
  !
  ! Loop over each solution point performing a double integral on the
  ! corresponding Lagrange polynomial over the domain of the triangle
  ! i.e., ipt    : current solution point
  !       L(ipt) : Lagrange polynomial for ipt
  !
  !   /\ 1    /\ 1-xi
  !   |       |
  !   | d(xi) | d(eta) [L(ipt)]
  !   |       |
  !  \/ 0    \/ 0
  !
  do np = 1,nspts
    !
    ! Outer integral loop
    !
    outer_integral: do n = 1,nquad
      !
      ! Inner integral loop
      !
      h = zero
      inner_integral: do m = 1,nquad
        !
        x = xi(n)
        y = (one-xi(n))*xi(m)
        !
        ! Evaluate the polynomial at the point x,y
        !
        k = 0
        g = zero
        do i = 1,nfpts
          do j = 1,i
            k = k + 1
            g = g + V(np,k) * x**(i-j) * y**(j-1)
          end do
        end do
        !
        ! Add the weighted polynomial value to the inner integral result
        !
        h = h + (one-xi(n)) * weights(m) * g
        !
      end do inner_integral
      !
      ! Add the weighted inner integral result to the outer integral result
      !
      return_value(np) = return_value(np) + weights(n) * h
      !
    end do outer_integral
    !
  end do
  !
end function get_triangle_weights_at_solpts_DP
!
!###############################################################################
!
pure function get_triangle_weights_at_solpts_QP(rs) result(return_value)
  !
  ! NOTE : This seems to fail for large n_order (>15)
  !
  use polynomial_mod, only : nodes_legendre_gauss
  use polynomial_mod, only : weights_legendre_gauss
  !
  real(qp), dimension(:,:), intent(in) :: rs
  !
  real(qp), dimension(1:size(rs,2)) :: return_value
  !
  integer  :: i,j,k,n,m,np
  integer  :: nspts,nfpts
  real(qp) :: g,h,x,y
  !
  integer, parameter :: nquad = 20
  !
  real(qp), dimension(1:nquad)                   :: xi
  real(qp), dimension(1:nquad)                   :: dxi
  real(qp), dimension(1:nquad)                   :: weights
  real(qp), dimension(1:size(rs,2),1:size(rs,2)) :: V
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: half = 0.5_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  ! On entry, r and s need to be defined on the interval [0,1]
  !
  nspts = size(rs,2)
  nfpts = np2nfp(nspts)
  !
  ! Compute the coefficients of the Vandermonde matrix for
  ! the Lagrange polynomials using monomial basis functions
  !
  do n = 1,nspts
    k = 0
    do i = 1,nfpts
      do j = 1,i
        k = k + 1
        V(k,n) = rs(1,n)**(i-j) * rs(2,n)**(j-1)
      end do
    end do
  end do
  !
  ! Now get the inverse of this matrix for use later
  !
  V = invert_matrix(V)
  !
  ! Get the Gauss quadrature points on the interval [-1,1] and
  ! the corresponding weights necessary to perform the integration
  !
  call nodes_legendre_gauss(nquad,xi,dxi)
  call weights_legendre_gauss(nquad,xi,dxi,weights)
  !
  ! Reset the quadrature points and weights to the interval [0,1]
  !
  xi(:) = half*(one+xi(:))
  weights(:) = half*weights(:)
  !
  ! Initialize the triangle weights to zero
  !
  return_value(:) = zero
  !
  ! Loop over each solution point performing a double integral on the
  ! corresponding Lagrange polynomial over the domain of the triangle
  ! i.e., ipt    : current solution point
  !       L(ipt) : Lagrange polynomial for ipt
  !
  !   /\ 1    /\ 1-xi
  !   |       |
  !   | d(xi) | d(eta) [L(ipt)]
  !   |       |
  !  \/ 0    \/ 0
  !
  do np = 1,nspts
    !
    ! Outer integral loop
    !
    outer_integral: do n = 1,nquad
      !
      ! Inner integral loop
      !
      h = zero
      inner_integral: do m = 1,nquad
        !
        x = xi(n)
        y = (one-xi(n))*xi(m)
        !
        ! Evaluate the polynomial at the point x,y
        !
        k = 0
        g = zero
        do i = 1,nfpts
          do j = 1,i
            k = k + 1
            g = g + V(np,k) * x**(i-j) * y**(j-1)
          end do
        end do
        !
        ! Add the weighted polynomial value to the inner integral result
        !
        h = h + (one-xi(n)) * weights(m) * g
        !
      end do inner_integral
      !
      ! Add the weighted inner integral result to the outer integral result
      !
      return_value(np) = return_value(np) + weights(n) * h
      !
    end do outer_integral
    !
  end do
  !
end function get_triangle_weights_at_solpts_QP
!
!###############################################################################
!
include "Functions/np2nfp.f90"
!
!###############################################################################
!
end module quadrature_mod
