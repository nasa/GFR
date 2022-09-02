module quadrature_mod
  !
!!!#define PROFILE_ON
!!!#define DISABLE_PURE
  !
  use module_kind_types, global_debug_timer =>      debug_timer, &
                    pure_global_debug_timer => pure_debug_timer
  !
  use element_type_mod, only : element_t
  !
  implicit none
  !
  private
  !
  public :: find_valid_geometries
  public :: init_geom_quadrature_rules
  public :: init_face_quadrature_rules
  public :: get_edge_quadrature
  public :: get_tria_quadrature
  public :: get_quad_quadrature
  public :: get_ref_elem_node_idx_for_face
  public :: map_face_pts
  public :: quadrature_memory_usage
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
subroutine find_valid_geometries()
  !
  !.. Use Statements ..
  use geovar,    only : cell_geom,cell_order
  use order_mod, only : n_min_geom,n_max_geom,q_order
  use order_mod, only : geom_solpts,geom_flxpts
  use order_mod, only : maxpts,maxSP,maxFP,maxQP,maxGP
  use order_mod, only : o_order,maxOP
  use order_mod, only : e_order,maxEP
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
subroutine init_geom_quadrature_rules()
  !
  !.. Use Statements ..
  use order_mod, only : n_order,geom_solpts
  use order_mod, only : n_min_geom,n_max_geom
  use order_mod, only : n_min_order,n_max_order
  use ovar,      only : loc_solution_pts
  use ovar,      only : loc_triangle_pts
  use ovar,      only : loc_tetrahedra_pts
  !
  use element_type_mod, only : std_elem
  use element_type_mod, only : assign_solpts_edge_ptr
  !
  !.. Local Scalars ..
  integer :: n,ierr
  integer :: this_geom,this_order
  character(len=300) :: timer_message
  !
  !.. Local Pointers ..
  type(element_t), pointer :: this_elem
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_geom_quadrature_rules"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
 !call dump_solpts
  !
  ! Initialize this_elem to disassociated
  !
  this_elem => null()
  !
  ! Allocate the std_elem array
  !
  if (.not. allocated(std_elem)) then
    allocate ( std_elem(Geom_Min:Geom_Max,n_min_order:n_max_order) , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"std_elem",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Loop through all the possible orders and get the quadrature rules
  ! for the most basic geometries that are used (Geom_Edge and Geom_Tria).
  ! We need to do this before all other geometry types because the
  ! remaining geometries are dependent on these quadrature rules.
  ! NOTE: The pts array within the type element_t is only needed for these
  !       two geometries (maybe tetrahedra as well) so they will be the only
  !       ones in which this array will be allocated and created.
  !
  call check_for_bad_number(__LINE__,__FILE__)
  !
  this_geom = Geom_Edge
  !
  basic_geom_loop: do n = 1,2
    !
    basic_order_loop: do this_order = n_min_order,n_max_order
      write (timer_message,1) trim(adjustl(Geom_Name(this_geom))),this_order
      call debug_timer(start_timer,timer_message)
      !
      ! Assign the pointer 'this_elem' to std_elem(this_geom,this_order)
      ! to simplify the remaining code in this loop and make it more legible
      !
      this_elem => std_elem(this_geom,this_order)
      !
      ! If wts is already allocated, then the rest of this loop was done before
      !
      if (allocated(this_elem%wts)) then
        if (associated(this_elem)) this_elem => null()
        cycle basic_order_loop
      end if
      !
      ! Fill in the solution points and quadrature weights
      ! for an edge with the current order
      !
      if (this_geom == Geom_Edge) then
        this_elem = get_edge_quadrature(this_order,loc_solution_pts)
      else if (this_geom == Geom_Tria) then
        this_elem = get_tria_quadrature(this_order,loc_triangle_pts)
      end if
      call output_quadrature(this_elem,this_geom,this_order,"Cell")
      !
      ! Disassociate this_elem before continuing to the next loop combo
      !
      if (associated(this_elem)) this_elem => null()
      call debug_timer(stop_timer,timer_message)
      call check_for_bad_number(__LINE__,__FILE__)
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
      write (timer_message,1) trim(adjustl(Geom_Name(this_geom))),this_order
      call debug_timer(start_timer,timer_message)
      !
      ! Assign the pointer 'this_elem' to std_elem(this_geom,this_order)
      ! to simplify the remaining code in this loop and make it more legible
      !
      this_elem => std_elem(this_geom,this_order)
      !
      ! If wts is already allocated, then the rest of this loop was done before
      !
      if (allocated(this_elem%wts)) then
        if (associated(this_elem)) this_elem => null()
        cycle geom_loop
      end if
      !
      ! Fill in the quadrature points and weights for the
      ! current geometry/order combination
      !
      select case (this_geom)
        case (Geom_Node)
          !
          this_elem%pts = init_array(1,1,zero) ! F2003 AUTO-REALLOC
          this_elem%wts = init_array(1,one)    ! F2003 AUTO-REALLOC
          !
        case (Geom_Quad)
          !
          this_elem = get_quad_quadrature(this_order,loc_solution_pts)
          !
        case (Geom_Tetr)
          !
          this_elem = get_tetr_quadrature(this_order,loc_tetrahedra_pts)
          !
        case (Geom_Pyra)
          !
          this_elem = get_pyra_quadrature(this_order,loc_solution_pts)
          !
        case (Geom_Pris)
          !
          this_elem = get_pris_quadrature(this_order,loc_triangle_pts, &
                                                     loc_solution_pts)
          !
        case (Geom_Hexa)
          !
          this_elem = get_hexa_quadrature(this_order,loc_solution_pts)
          !
        case default
          !
          write (error_message,2)
          call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
          !
      end select
      !
      call output_quadrature(this_elem,this_geom,this_order,"Cell")
      !
      ! Disassociate the pointer 'this_elem' before
      ! continuing to the next geometry
      !
      if (associated(this_elem)) this_elem => null()
      call debug_timer(stop_timer,timer_message)
      call check_for_bad_number(__LINE__,__FILE__)
      !
    end do geom_loop
  end do order_loop
  !
  ! Before we leave, assign the pointer solpts_edge to its target
  !
  call assign_solpts_edge_ptr(n_order)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("Geometry = ",a,"  ;  Order = ",i0)
  2 format (" A grid cell of an unknown geometry type was found!")
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
  use ovar,      only : loc_flux_pts
  use ovar,      only : loc_triangle_pts
  !
  use element_type_mod, only : face_elem
  !
  !.. Local Scalars ..
  integer :: n,ierr
  integer :: gmin,gmax,this_order
  logical(lk) :: there_are_quad_faces
  logical(lk) :: there_are_tria_faces
  character(len=300) :: timer_message
  !
  !.. Local Pointers ..
  type(element_t), pointer :: this_elem
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_face_quadrature_rules"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Initialize this_elem to disassociated
  !
  this_elem => null()
  !
  ! Determine if there are any of the geometries used for this simulation
  ! are 3D to see if we need to find the quadrature rules for 2D faces.
  !
  gmin = Geom_Edge
  gmax = maxval( Geom_2D )
  !
  there_are_tria_faces = fals
  there_are_quad_faces = fals
  !
  if (any(geom_is_used(Geom_3D))) then
   !gmax = maxval( Geom_2D , mask=geom_is_used(Geom_2D) )
    there_are_tria_faces = geom_is_used(Geom_Tria)
    there_are_quad_faces = geom_is_used(Geom_Quad)
  end if
  !
  ! Allocate the face_elem array if not already available
  !
  if (.not. allocated(face_elem)) then
    allocate ( face_elem(gmin:gmax,n_min_order:n_max_order) , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"face_elem",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Loop through all the possible orders and get the edge quadrature rules.
  !
  call check_for_bad_number(__LINE__,__FILE__)
  !
  edge_loop: do this_order = n_min_order,n_max_order
    write (timer_message,1) trim(adjustl(Geom_Name(Geom_Edge))),this_order
    call debug_timer(start_timer,timer_message)
    !
    ! Assign the pointer 'this_elem' to face_elem(Geom_Edge,this_order)
    ! to simplify the remaining code in this loop and make it more legible
    !
    this_elem => face_elem(Geom_Edge,this_order)
    !
    ! If wts is already allocated, then the rest of this loop was done before
    !
    if (allocated(this_elem%wts)) then
      if (associated(this_elem)) this_elem => null()
      cycle edge_loop
    end if
    !
    ! Fill in the solution points and quadrature weights
    ! for an edge with the current order
    !
    this_elem = get_edge_quadrature(this_order,loc_flux_pts)
    !
    call output_quadrature(this_elem,Geom_Edge,this_order,"Face")
    !
    ! Disassociate this_elem before continuing to the next loop combo
    !
    if (associated(this_elem)) this_elem => null()
    call debug_timer(stop_timer,timer_message)
    call check_for_bad_number(__LINE__,__FILE__)
    !
  end do edge_loop
  !
  ! Loop through all the possible orders and get the triangle quadrature rules
  ! if there are any triangle faces exist for this simulation
  !
  if (there_are_tria_faces) then
    !
    tria_loop: do this_order = n_min_order,n_max_order
      write (timer_message,1) trim(adjustl(Geom_Name(Geom_Tria))),this_order
      call debug_timer(start_timer,timer_message)
      !
      ! Assign the pointer 'this_elem' to face_elem(Geom_Tria,this_order)
      ! to simplify the remaining code in this loop and make it more legible
      !
      this_elem => face_elem(Geom_Tria,this_order)
      !
      ! If wts is already allocated, then the rest of this loop was done before
      !
      if (allocated(this_elem%wts)) then
        if (associated(this_elem)) this_elem => null()
        cycle tria_loop
      end if
      !
      ! Fill in the solution points and quadrature weights
      ! for a triangle face with the current order
      !
      this_elem = get_tria_quadrature(this_order,loc_triangle_pts)
      !
      call output_quadrature(this_elem,Geom_Tria,this_order,"Face")
      !
      ! Disassociate this_elem before continuing to the next loop combo
      !
      if (associated(this_elem)) this_elem => null()
      call debug_timer(stop_timer,timer_message)
      call check_for_bad_number(__LINE__,__FILE__)
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
      write (timer_message,1) trim(adjustl(Geom_Name(Geom_Quad))),this_order
      call debug_timer(start_timer,timer_message)
      !
      ! Assign the pointer 'this_elem' to face_elem(Geom_Quad,this_order)
      ! to simplify the remaining code in this loop and make it more legible
      !
      this_elem => face_elem(Geom_Quad,this_order)
      !
      ! If wts is already allocated, then the rest of this loop was done before
      !
      if (allocated(this_elem%wts)) then
        if (associated(this_elem)) this_elem => null()
        cycle quad_loop
      end if
      !
      ! Fill in the solution points and quadrature weights
      ! for a quad face with the current order
      !
      this_elem = get_quad_quadrature(this_order,loc_flux_pts)
      !
      call output_quadrature(this_elem,Geom_Quad,this_order,"Face")
      !
      ! Disassociate this_elem before continuing to the next loop combo
      !
      if (associated(this_elem)) this_elem => null()
      call debug_timer(stop_timer,timer_message)
      call check_for_bad_number(__LINE__,__FILE__)
      !
    end do quad_loop
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("Geometry = ",a,"  ;  Order = ",i0)
  !
end subroutine init_face_quadrature_rules
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_generic_quadrature_wts(this_geom,this_order,pts) &
                             result(return_value)
  !
  !.. Use Statements ..
 !use vandermonde_mod, only : Compute_Basis_Matrix
  use vandermonde_mod, only : General_Basis_Matrix
  use vandermonde_mod, only : Vandermonde_Inverse
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: this_geom
  integer,  intent(in) :: this_order
  real(wp), intent(in) :: pts(:,:)
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: vmat(:,:)
  real(wp), allocatable :: vmatI(:,:)
  real(wp), allocatable :: mmatI(:,:)
  !
continue
  !
  ! Compute the Vandermonde matrix for the solution points
  !
 !vmat = Compute_Basis_Matrix(this_geom,this_order,pts) ! F2003 AUTO-REALLOC
  vmat = General_Basis_Matrix(this_geom,this_order,pts) ! F2003 AUTO-REALLOC
  !
  ! Compute the inverse of the Vandermonde matrix
  !
  vmatI = Vandermonde_Inverse( vmat ) ! F2003 AUTO-REALLOC
  !
  ! Compute the inverse mass matrix
  !
  mmatI = matmul( transpose(vmatI) , vmatI ) ! F2003 AUTO-REALLOC
  !
  ! The quadrature weights are just the sum of each row of mmatI
  ! NOTE: mmatI is a square, symmetric matrix, so it doesnt matter whether
  !       we sum the rows or columns
  !
  return_value = sum( mmatI , dim=1 ) ! F2003 AUTO-REALLOC
  !
end function get_generic_quadrature_wts
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_edge_quadrature(n,location) result(return_value)
  !
  use order_mod,      only : geom_solpts
  use polynomial_mod, only : legendre_gauss_quadrature
  use polynomial_mod, only : legendre_gauss_lobatto_quadrature
  !
  integer, intent(in) :: n
  integer, intent(in) :: location
  !
  type(element_t) :: return_value
  !
  integer :: ndim,npts,ierr
  !
  real(mp), allocatable :: pts(:)
  real(mp), allocatable :: wts(:)
  !
continue
  !
  ndim = geom_dimen(Geom_Edge)
  npts = geom_solpts(Geom_Edge,n)
  !
  pts = init_array(npts,real(zero,kind=kind(pts))) ! F2003 AUTO-REALLOC
  wts = init_array(npts,real(zero,kind=kind(wts))) ! F2003 AUTO-REALLOC
  !
  select case (location)
    !
    case (Legendre_Gauss)
      !
      ! Edge quadrature using Legendre Gauss nodes
      !
      call legendre_gauss_quadrature(pts,wts)
#ifdef DISABLE_PURE
      call check_for_bad_number(__LINE__,__FILE__)
#endif
      !
    case (Legendre_Gauss_Lobatto)
      !
      ! Edge quadrature using Legendre Gauss Lobatto nodes
      !
      call legendre_gauss_lobatto_quadrature(pts,wts)
#ifdef DISABLE_PURE
      call check_for_bad_number(__LINE__,__FILE__)
#endif
      !
    case default
      !
      ! Default to Legendre Gauss points
      !
      call legendre_gauss_quadrature(pts,wts)
#ifdef DISABLE_PURE
      call check_for_bad_number(__LINE__,__FILE__)
#endif
      !
  end select
  !
  return_value%pts = reshape( chop(pts) , [ndim,npts] ) ! F2003 AUTO-REALLOC
  return_value%wts = chop(wts) ! F2003 AUTO-REALLOC
  !
  ! NOTE: The following call to get_generic_quadrature_wts
  !       will give the same set of quadrature weights
  ! F2003 AUTO-REALLOC
 !return_value%wts = get_generic_quadrature_wts(Geom_Edge,n,return_value%pts)
  !
end function get_edge_quadrature
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_tria_quadrature(p,location) result(return_value)
  !
  !.. Use Statements ..
  use triangle_mod, only : TriNodes2D
  !
  !.. Formal Arguments ..
  integer, intent(in) :: p
  integer, intent(in) :: location
  !
  !.. Function Result ..
  type(element_t) :: return_value
  !
  !.. Local Allocatable Arrays ..
  real(mp), allocatable :: pts(:,:)
  !
continue
  !
  ! Get the solution points
  !
  pts = TriNodes2D(p,location) ! F2003 AUTO-REALLOC
  !
  ! Save the solution points into return_value
  !
  return_value%pts = real(pts,kind=kind(return_value%pts)) ! F2003 AUTO-REALLOC
  !
  ! Compute the quadrature weights using the inverse mass matrix
  !
  ! F2003 AUTO-REALLOC
  return_value%wts = get_generic_quadrature_wts(Geom_Tria,p,return_value%pts)
  !
end function get_tria_quadrature
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_quad_quadrature(n,location) result(return_value)
  !
  use order_mod, only : geom_solpts
  !
  integer, intent(in) :: n
  integer, intent(in) :: location
  !
  type(element_t) :: return_value
  !
  integer :: i,j,l,ndim,npts,ierr
  !
  real(mp), allocatable :: pts(:,:)
  real(mp), allocatable :: wts(:)
  !
  type(element_t) :: edge_elem
  !
continue
  !
  edge_elem = get_edge_quadrature(n,location)
  !
  ndim = geom_dimen(Geom_Quad)
  npts = geom_solpts(Geom_Quad,n)
  !
  ! Compute the tensor product of the 1D solution points
  !
  pts = init_array(ndim,npts,real(zero,kind=kind(pts))) ! F2003 AUTO-REALLOC
  !
  l = 0
  do j = 1,size(edge_elem%pts,dim=2)
    do i = 1,size(edge_elem%pts,dim=2)
      l = l + 1
      pts(1,l) = real( edge_elem%pts(1,i) , kind=kind(pts) )
      pts(2,l) = real( edge_elem%pts(1,j) , kind=kind(pts) )
    end do
  end do
  !
  return_value%pts = chop( pts ) ! F2003 AUTO-REALLOC
  !
  ! Compute the tensor product of the 1D weights
  !
  wts = init_array(npts,real(zero,kind=kind(wts))) ! F2003 AUTO-REALLOC
  !
  l = 0
  do j = 1,size(edge_elem%wts)
    do i = 1,size(edge_elem%wts)
      l = l + 1
      wts(l) = real( edge_elem%wts(i) , kind=kind(wts) ) *  &
               real( edge_elem%wts(j) , kind=kind(wts) )
    end do
  end do
  !
  return_value%wts = chop( wts ) ! F2003 AUTO-REALLOC
  !
  ! NOTE: The following call to get_generic_quadrature_wts
  !       will give the same set of quadrature weights
  ! F2003 AUTO-REALLOC
 !return_value%wts = get_generic_quadrature_wts(Geom_Quad,n,return_value%pts)
  !
end function get_quad_quadrature
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_tetr_quadrature(p,location) result(return_value)
  !
  !.. Use Statements ..
  use tetrahedra_mod, only : TetNodes3D
  !
  !.. Formal Arguments ..
  integer, intent(in) :: p
  integer, intent(in) :: location
  !
  !.. Function Result ..
  type(element_t) :: return_value
  !
  !.. Local Allocatable Arrays ..
  real(mp), allocatable :: pts(:,:)
  real(mp), allocatable :: wts(:)
  !
continue
  !
  ! Get the solution points
  !
  pts = TetNodes3D(p,location) ! F2003 AUTO-REALLOC
 !pts = TetNodes3D(p,location,use_Gauss) ! F2003 AUTO-REALLOC
  !
  ! Save the solution points into return_value
  !
  return_value%pts = real(pts,kind=kind(return_value%pts)) ! F2003 AUTO-REALLOC
  !
  ! Compute the quadrature weights using the inverse mass matrix
  !
  ! F2003 AUTO-REALLOC
  return_value%wts = get_generic_quadrature_wts(Geom_Tetr,p,return_value%pts)
  !
end function get_tetr_quadrature
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_pyra_quadrature(p,location) result(return_value)
  !
  !.. Use Statements ..
  use pyramid_mod, only : PyrNodes3D
  !
  !.. Formal Arguments ..
  integer, intent(in) :: p
  integer, intent(in) :: location
  !
  !.. Function Result ..
  type(element_t) :: return_value
  !
  !.. Local Scalars ..
  integer :: node_loc
  !
  !.. Local Allocatable Arrays ..
  real(mp), allocatable :: pts(:,:)
  !
continue
  !
  node_loc = Legendre_Gauss
  if (location == Legendre_Gauss_Lobatto) then
    node_loc = Legendre_Gauss_Lobatto
  end if
  !
  ! Get the solution points
  !
  pts = PyrNodes3D(p,node_loc) ! F2003 AUTO-REALLOC
  !
  ! Save the solution points into return_value
  !
  return_value%pts = real(pts,kind=kind(return_value%pts)) ! F2003 AUTO-REALLOC
  !
  ! Compute the quadrature weights using the inverse mass matrix
  !
  ! F2003 AUTO-REALLOC
  return_value%wts = get_generic_quadrature_wts(Geom_Pyra,p,return_value%pts)
  !
end function get_pyra_quadrature
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_pris_quadrature(p,tria_loc,edge_loc) result(return_value)
  !
  !.. Use Statements ..
  use prism_mod, only : PriNodes3D
  !
  !.. Formal Arguments ..
  integer, intent(in) :: p
  integer, intent(in) :: tria_loc
  integer, intent(in) :: edge_loc
  !
  !.. Function Result ..
  type(element_t) :: return_value
  !
  !.. Local Allocatable Arrays ..
  real(mp), allocatable :: pts(:,:)
  !
continue
  !
  ! Get the solution points
  !
  pts = PriNodes3D(p,tria_loc,edge_loc) ! F2003 AUTO-REALLOC
 !pts = PriNodes3D(p,tria_loc,edge_loc,tria_edge_loc) ! F2003 AUTO-REALLOC
  !
  ! Save the solution points into return_value
  !
  return_value%pts = real(pts,kind=kind(return_value%pts)) ! F2003 AUTO-REALLOC
  !
  ! Compute the quadrature weights using the inverse mass matrix
  !
  ! F2003 AUTO-REALLOC
  return_value%wts = get_generic_quadrature_wts(Geom_Pris,p,return_value%pts)
  !
end function get_pris_quadrature
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_hexa_quadrature(n,location) result(return_value)
  !
  use order_mod, only : geom_solpts
  !
  integer, intent(in) :: n
  integer, intent(in) :: location
  !
  type(element_t) :: return_value
  !
  integer :: i,j,k,l,ndim,npts,ierr
  !
  real(mp), allocatable :: pts(:,:)
  real(mp), allocatable :: wts(:)
  !
  type(element_t) :: edge_elem
  !
continue
  !
  edge_elem = get_edge_quadrature(n,location)
  !
  ndim = geom_dimen(Geom_Hexa)
  npts = geom_solpts(Geom_Hexa,n)
  !
  ! Compute the tensor product of the 1D solution points
  !
  pts = init_array(ndim,npts,real(zero,kind=kind(pts))) ! F2003 AUTO-REALLOC
  !
  l = 0
  do k = 1,size(edge_elem%pts,dim=2)
    do j = 1,size(edge_elem%pts,dim=2)
      do i = 1,size(edge_elem%pts,dim=2)
        l = l + 1
        pts(1,l) = real( edge_elem%pts(1,i) , kind=kind(pts) )
        pts(2,l) = real( edge_elem%pts(1,j) , kind=kind(pts) )
        pts(3,l) = real( edge_elem%pts(1,k) , kind=kind(pts) )
      end do
    end do
  end do
  !
  return_value%pts = chop( pts ) ! F2003 AUTO-REALLOC
  !
  ! Compute the tensor product of the 1D weights
  !
  wts = init_array(npts,real(zero,kind=kind(wts))) ! F2003 AUTO-REALLOC
  !
  l = 0
  do k = 1,size(edge_elem%wts)
    do j = 1,size(edge_elem%wts)
      do i = 1,size(edge_elem%wts)
        l = l + 1
        wts(l) = real( edge_elem%wts(i) , kind=kind(wts) ) * &
                 real( edge_elem%wts(j) , kind=kind(wts) ) * &
                 real( edge_elem%wts(k) , kind=kind(wts) )
      end do
    end do
  end do
  !
  return_value%wts = chop( wts ) ! F2003 AUTO-REALLOC
  !
  ! NOTE: The following call to get_generic_quadrature_wts
  !       will give the same set of quadrature weights
  ! F2003 AUTO-REALLOC
 !return_value%wts = get_generic_quadrature_wts(Geom_Hexa,n,return_value%pts)
  !
end function get_hexa_quadrature
!
!###############################################################################
!
subroutine get_std_elem_face_pts()
  !
  !.. Use Statements ..
  use element_type_mod, only : std_elem
  !
  !.. Local Scalars ..
  integer :: this_geom,this_order
  !
continue
  !
  do this_geom = lbound(std_elem,dim=1),ubound(std_elem,dim=1)
    do this_order = lbound(std_elem,dim=2),ubound(std_elem,dim=2)
      !
      if (allocated(std_elem(this_geom,this_order)%pts)) then
        call get_face_pts(std_elem(this_geom,this_order),this_geom,this_order)
      end if
      !
    end do
  end do
  !
end subroutine get_std_elem_face_pts
!
!###############################################################################
!
pure &
function get_ref_elem_node_idx_for_face(this_geom,this_face) &
                                 result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: this_geom
  integer, intent(in) :: this_face
  !
  !.. Function Result ..
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: nfnods
  !
continue
  !
  nfnods = num_face_nodes(this_face,this_geom)
  !
  return_value = side_nodes(1:nfnods,this_face,this_geom)+1 ! F2003 AUTO-REALLOC
  !
  ! We need to switch the node indices for some special cases
  !
  if (this_geom == Geom_Tetr) then
    !
    if (this_face == 3) then
      return_value = return_value([2,1,3]) ! F2003 AUTO-REALLOC
    else if (this_face == 4) then
      return_value = return_value([1,3,2]) ! F2003 AUTO-REALLOC
    end if
    !
  else if (this_geom == Geom_Pyra) then
    !
    if (any(this_face == [3,4])) then
      return_value = return_value([2,1,3]) ! F2003 AUTO-REALLOC
    else if (this_face == 5) then
      return_value = return_value([1,4,3,2]) ! F2003 AUTO-REALLOC
    end if
  else if (this_geom == Geom_Pris) then
    !
    if (this_face == 1) then
      return_value = return_value([1,3,2]) ! F2003 AUTO-REALLOC
    else if (this_face == 5) then
      return_value = return_value([1,4,3,2]) ! F2003 AUTO-REALLOC
    end if
    !
  else if (this_geom == Geom_Hexa) then
    !
    if (any(this_face == [1,4,6])) then
      return_value = return_value([1,4,3,2]) ! F2003 AUTO-REALLOC
    end if
    !
  end if
  !
end function get_ref_elem_node_idx_for_face
!
!###############################################################################
!
pure &
subroutine get_face_pts(this_elem,this_geom,this_order)
  !
  !.. Use Statements ..
  use order_mod,        only : geom_flxpts
  use element_type_mod, only : face_elem
  !
  !.. Formal Arguments ..
  integer, intent(in) :: this_geom
  integer, intent(in) :: this_order
  type(element_t), intent(inout) :: this_elem
  !
  !.. Local Scalars ..
  integer :: if1,if2,ndim,ntfp,nfp,nfnods
  integer :: this_face,face_geom
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: ip(:)
  real(wp), allocatable :: rs(:,:)
  real(wp), allocatable :: node_xyz(:,:)
  real(wp), allocatable :: face_pts(:,:)
  !
continue
  !
  if (this_geom == Geom_Edge) then
    this_elem%face_pts = reshape( [-one,one] , [1,2] ) ! F2003 AUTO-REALLOC
    this_elem%face_wts = [one,one] ! F2003 AUTO-REALLOC
    return
  end if
  !
  ndim = geom_dimen(this_geom)
  ntfp = geom_flxpts(this_geom,this_order)
  !
  this_elem%face_pts = init_array(ndim,ntfp,zero) ! F2003 AUTO-REALLOC
  this_elem%face_wts = init_array(ntfp,zero) ! F2003 AUTO-REALLOC
  !
  if2 = 0
  !
  do this_face = 1,geom_faces(this_geom)
    !
    face_geom = cell_face_geom(this_face,this_geom)
    !
    rs = face_elem(face_geom,this_order)%pts ! F2003 AUTO-REALLOC
    !
    ip = get_ref_elem_node_idx_for_face(this_geom,this_face) ! F03 AUTO-REALLOC
    !
    node_xyz = ref_elem_nodes_xyz(1:ndim,ip(:),this_geom) ! F2003 AUTO-REALLOC
    !
    face_pts = map_face_pts(node_xyz,rs) ! F2003 AUTO-REALLOC
    !
    if1 = if2 + 1
    if2 = if2 + size(face_pts,dim=2)
    !
    this_elem%face_pts(1:ndim,if1:if2) = face_pts(:,:)
    this_elem%face_wts(if1:if2) = face_elem(face_geom,this_order)%wts(:)
    !
  end do
  !
end subroutine get_face_pts
!
!###############################################################################
!
pure &
function map_face_pts(v,rs) result(return_value)
  !
  ! This function maps face points defined on a
  ! reference element to arbitrary face points
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: v(:,:)
  real(wp), intent(in) :: rs(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value( 1:size(v,dim=1) , 1:size(rs,dim=2) )
  !
continue
  !
  if (size(v,dim=2) == 2) then
    return_value = map_edge_pts(v,rs(1,:))
  else if (size(v,dim=2) == 3) then
    return_value = map_tria_pts(v,rs)
  else if (size(v,dim=2) == 4) then
    return_value = map_quad_pts(v,rs)
  else
    return_value = zero
  end if
  !
end function map_face_pts
!
!###############################################################################
!
pure &
function map_edge_pts(v,xi) result(return_value)
  !
  ! This function maps the points xi defined on the reference edge [-1,+1]
  ! to the arbitrary edge defined by [v(:,1),v(:,2)]
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: v(:,:)
  real(wp), intent(in) :: xi(:)
  !
  !.. Function Result ..
  real(wp) :: return_value( 1:size(v,dim=1) , 1:size(xi) )
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: x,xp1,xm1
  !
continue
  !
  do n = 1,size(xi)
    !
    x = xi(n)
    !
    xp1 = x + one
    xm1 = x - one
    !
    return_value(:,n) = half*( v(:,2)*xp1 - &
                               v(:,1)*xm1 )
    !
  end do
  !
end function map_edge_pts
!
!###############################################################################
!
pure &
function map_tria_pts(v,rs) result(return_value)
  !
  ! This function maps the points rs defined on the reference triangle
  ! to the arbitrary triangle defined by v
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: v(:,:)
  real(wp), intent(in) :: rs(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value( 1:size(v,dim=1) , 1:size(rs,dim=2) )
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: r,s,rp1,sp1,rps
  !
continue
  !
  do n = 1,size(rs,dim=2)
    !
    r = rs(1,n)
    s = rs(2,n)
    !
    rp1 = r + one
    sp1 = s + one
    rps = r + s
    !
    return_value(:,n) = half*( v(:,3)*sp1 + &
                               v(:,2)*rp1 - &
                               v(:,1)*rps )
    !
  end do
  !
end function map_tria_pts
!
!###############################################################################
!
pure &
function map_quad_pts(v,rs) result(return_value)
  !
  ! This function maps the points rs defined on the reference quad
  ! to the arbitrary quad defined by v
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: v(:,:)
  real(wp), intent(in) :: rs(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value( 1:size(v,dim=1) , 1:size(rs,dim=2) )
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: r,s,opr,omr,ops,oms
  !
continue
  !
  do n = 1,size(rs,dim=2)
    !
    r = rs(1,n)
    s = rs(2,n)
    !
    opr = one + r
    omr = one - r
    !
    ops = one + s
    oms = one - s
    !
    return_value(:,n) = one4*( v(:,1)*omr*oms + &
                               v(:,2)*opr*oms + &
                               v(:,3)*opr*ops + &
                               v(:,4)*omr*ops )
    !
  end do
  !
end function map_quad_pts
!
!###############################################################################
!
subroutine output_quadrature(this_elem,this_geom,this_order,this_type)
  !
  !.. Use Statements ..
  use ovar, only : mout => matrix_output_unit
  !
  !.. Formal Arguments ..
  class(element_t), intent(in) :: this_elem
  integer,          intent(in) :: this_geom
  integer,          intent(in) :: this_order
  character(len=*), intent(in) :: this_type
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "output_quadrature"
  !
continue
  !
  if (mout == iout) return
  !
  if (exclusive_element >= 0) then
    if (this_geom /= exclusive_element) then
      return
    end if
  end if
  !
  write (mout,1) this_order, &
                 trim(adjustl(Geom_Full_Name(this_geom))), &
                 trim(adjustl(this_type))
  !
  write (mout,2) "Solution"
  call output_pts_and_wts(this_geom,this_elem%pts,this_elem%wts)
  !
  if (allocated(this_elem%face_pts) .and. allocated(this_elem%face_wts)) then
    write (iout,2) "  Flux"
    call output_pts_and_wts(this_geom,this_elem%face_pts,this_elem%face_wts)
  else
    write (iout,3) merge("TRUE ","FALSE",allocated(this_elem%face_pts)), &
                   merge("TRUE ","FALSE",allocated(this_elem%face_wts))
  end if
  !
  1 format (/, &
            3(80("#"),/),/, &
            3x,"for a P",i0,1x,a,1x,a," element",/,/, &
            80("="))
  2 format (/, &
            2x,19("*"),/, &
            4x,a," Points",/, &
            2x,19("*"),/, &
            3x,"Reference Coordinates : Quadrature Weight",/, &
            80("-"))
  3 format (/, &
            2x,59("*"),/, &
            12x,"allocated( this_elem%face_pts ) = ",a,/, &
            12x,"allocated( this_elem%face_wts ) = ",a,/, &
            2x,59("*"))
  !
end subroutine output_quadrature
!
!###############################################################################
!
subroutine output_pts_and_wts(this_geom,pts,wts,mout)
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: this_geom
  real(wp), intent(in) :: pts(:,:)
  real(wp), intent(in) :: wts(:)
  integer,  intent(in) :: mout
  !
  !.. Local Scalars ..
  integer :: i,j,n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "output_pts_and_wts"
  !
continue
  !
  n1 = size(pts,dim=1)
  n2 = size(pts,dim=2)
  !
  if (geom_dimen(this_geom) == 1) then
    do j = 1,n2
      write (mout,1) (pts(i,j),i=1,n1), wts(j)
    end do
    write (mout,4) sum(wts(:))
  else if (geom_dimen(this_geom) == 2) then
    do j = 1,n2
      write (mout,2) (pts(i,j),i=1,n1), wts(j)
    end do
    write (mout,5) sum(wts(:))
  else if (geom_dimen(this_geom) == 3) then
    do j = 1,n2
      write (mout,3) (pts(i,j),i=1,n1), wts(j)
    end do
    write (mout,6) sum(wts(:))
  end if
  !
  1 format (1x,f10.6,"  :  ",f10.6)
  2 format (2(1x,f10.6),"  :  ",f10.6)
  3 format (3(1x,f10.6),"  :  ",f10.6)
  4 format (16x,f10.6)
  5 format (27x,f10.6)
  6 format (38x,f10.6)
  !
end subroutine output_pts_and_wts
!
!###############################################################################
!
subroutine quadrature_memory_usage(iunit)
  !
  !.. Use Statements ..
  use module_memory
  use element_type_mod, only : std_elem
  use element_type_mod, only : face_elem
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
  if (allocated(std_elem)) then
    !
    call dt%set( storage_size(std_elem,kind=inttype) , &
                         size(std_elem,kind=inttype) )
    !
    do j = lbound(std_elem,dim=2),ubound(std_elem,dim=2)
      do i = lbound(std_elem,dim=1),ubound(std_elem,dim=1)
        !
        if (allocated(std_elem(i,j)%pts)) then
          !
          call array%set( storage_size(std_elem(i,j)%pts,kind=inttype) , &
                                  size(std_elem(i,j)%pts,kind=inttype) )
          !
          write (array_name,4) "std_elem",Geom_Name(i),j,"pts"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
        end if
        !
        if (allocated(std_elem(i,j)%wts)) then
          !
          call array%set( storage_size(std_elem(i,j)%wts,kind=inttype) , &
                                  size(std_elem(i,j)%wts,kind=inttype) )
          !
          write (array_name,4) "std_elem",Geom_Name(i),j,"wts"
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
    write (array_name,5) "std_elem"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "std_elem"
    !
  end if
  !
  if (allocated(face_elem)) then
    !
    call dt%set( storage_size(face_elem,kind=inttype) , &
                         size(face_elem,kind=inttype) )
    !
    do j = lbound(face_elem,dim=2),ubound(face_elem,dim=2)
      do i = lbound(face_elem,dim=1),ubound(face_elem,dim=1)
        !
        if (allocated(face_elem(i,j)%pts)) then
          !
          call array%set( storage_size(face_elem(i,j)%pts,kind=inttype) , &
                                  size(face_elem(i,j)%pts,kind=inttype) )
          !
          write (array_name,4) "face_elem",Geom_Name(i),j,"pts"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
        end if
        !
        if (allocated(face_elem(i,j)%wts)) then
          !
          call array%set( storage_size(face_elem(i,j)%wts,kind=inttype) , &
                                  size(face_elem(i,j)%wts,kind=inttype) )
          !
          write (array_name,4) "face_elem",Geom_Name(i),j,"wts"
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
    write (array_name,5) "face_elem"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "face_elem"
    !
  end if
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module quadrature_mod")
  2 format (" Array: ",a,dt)
  3 format (" Total memory = ",f8.2,1x,a)
  4 format (a,"(",a,",",i0,")%",a)
  5 format ("'",a,"'")
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  !
#endif
end subroutine quadrature_memory_usage
!
!###############################################################################
!
subroutine dump_solpts()
  !
  !.. Use Statements ..
  use order_mod, only : n_order,geom_solpts
  use order_mod, only : n_min_geom,n_max_geom
  use order_mod, only : n_min_order,n_max_order
  !
  !.. Local Scalars ..
  integer :: ierr,location,io
  integer :: this_geom,this_order
  !
  !.. Local Allocatable Arrays ..
  type(element_t), target, allocatable, dimension(:,:) :: edge_elem
  type(element_t), target, allocatable, dimension(:,:) :: quad_elem
  type(element_t), target, allocatable, dimension(:,:) :: hexa_elem
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "dump_solpts"
  character(len=*), parameter :: filename = "geom_pts.dat"
  integer, parameter :: min_order = 1
  integer, parameter :: max_order = 8
  !
continue
  !
  open (newunit=io,file=filename,status="replace",action="write", &
        position="rewind",form="formatted",iostat=ierr,iomsg=error_message)
  call io_error(pname,filename,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (io,1)
  !
  allocate ( edge_elem(min_order:max_order,1:3) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"edge_elem",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( quad_elem(min_order:max_order,1:3) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"quad_elem",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( hexa_elem(min_order:max_order,1:3) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"hexa_elem",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Loop through all the possible orders and get the quadrature rules
  ! for the most basic geometries that are used (Geom_Edge and Geom_Tria).
  ! We need to do this before all other geometry types because the
  ! remaining geometries are dependent on these quadrature rules.
  ! NOTE: The pts array within the type element_t is only needed for these
  !       two geometries (maybe tetrahedra as well) so they will be the only
  !       ones in which this array will be allocated and created.
  !
  this_geom = Geom_Edge
  !
  do location = 1,3
    !
    basic_order_loop: do this_order = min_order,max_order
      !
      ! Fill in the solution points and quadrature weights
      ! for an edge with the current order
      !
      edge_elem(this_order,location) = get_edge_points(this_order,location)
      !
    end do basic_order_loop
    !
  end do
  !
  ! Change this_geom to Geom_Tria so the second pass through
  ! basic_order_loop works on triangle geometry
  !
  this_geom = Geom_Quad
  !
  ! Loop through all the possible combinations of cell geometry and order
  ! and get the quadrature rules for each one
  !
  quad_loc_loop: do location = 1,3
    !
    quad_order_loop: do this_order = min_order,max_order
      !
      ! Fill in the quadrature weights for the
      ! current geometry/order combination
      !
      quad_elem(this_order,location) = &
                               get_quad_points( edge_elem(this_order,location) )
      !
    end do quad_order_loop
    !
  end do quad_loc_loop
  !
  ! Change this_geom to Geom_Tria so the second pass through
  ! basic_order_loop works on triangle geometry
  !
  this_geom = Geom_Hexa
  !
  ! Loop through all the possible combinations of cell geometry and order
  ! and get the quadrature rules for each one
  !
  hexa_loc_loop: do location = 1,3
    !
    hexa_order_loop: do this_order = min_order,max_order
      !
      ! Fill in the quadrature weights for the
      ! current geometry/order combination
      !
      hexa_elem(this_order,location) = &
                               get_hexa_points( edge_elem(this_order,location) )
      !
    end do hexa_order_loop
    !
  end do hexa_loc_loop
  !
  this_geom = Geom_Quad
  do this_order = min_order,max_order
    call write_points(quad_elem(this_order,3),this_geom,this_order,3)
    call write_points(quad_elem(this_order,2),this_geom,this_order,2)
    call write_points(quad_elem(this_order,1),this_geom,this_order,1)
  end do
  !
  this_geom = Geom_Hexa
  do this_order = min_order,max_order
    call write_points(hexa_elem(this_order,3),this_geom,this_order,3)
    call write_points(hexa_elem(this_order,2),this_geom,this_order,2)
    call write_points(hexa_elem(this_order,1),this_geom,this_order,1)
  end do
  !
  close (io,iostat=ierr,iomsg=error_message)
  call io_error(pname,filename,2,__LINE__,__FILE__,ierr,error_message)
  !
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,"dumping solution points")
  !
  ! Format Statements
  !
  1 format ('VARIABLES = "X","Y","Z"')
  !
  contains
  !
  !#############################################################################
  !
#ifndef DISABLE_PURE
  pure &
#endif
  function get_edge_points(n,nloc) result(return_value)
    !
    !.. Use Statements ..
    use order_mod,      only : geom_solpts
    use polynomial_mod, only : legendre_gauss_quadrature
    use polynomial_mod, only : legendre_gauss_lobatto_quadrature
    !
    !.. Formal Arguments ..
    integer, intent(in) :: n
    integer, intent(in) :: nloc
    !
    !.. Function Result ..
    type(element_t) :: return_value
    !
    !.. Local Scalars ..
    integer :: nd,np,ierr
    !
    !.. Local Allocatable Arrays ..
    real(mp), allocatable :: pts(:)
    real(mp), allocatable :: wts(:)
    !
  continue
    !
    nd = geom_dimen(Geom_Edge)
    np = geom_solpts(Geom_Edge,n)
    !
    pts = init_array(np,real(zero,kind=kind(pts))) ! F2003 AUTO-REALLOC
    wts = init_array(np,real(zero,kind=kind(wts))) ! F2003 AUTO-REALLOC
    !
    select case (nloc)
      !
      case (Legendre_Gauss)
        !
        ! Edge quadrature using Legendre Gauss nodes
        !
        call legendre_gauss_quadrature(pts,wts)
        !
      case (Legendre_Gauss_Lobatto)
        !
        ! Edge quadrature using Legendre Gauss Lobatto nodes
        !
        call legendre_gauss_lobatto_quadrature(pts,wts)
        !
      case (3)
        !
        ! Edge quadrature using [-1,Legendre Gauss nodes,+1]
        !
        call legendre_gauss_quadrature(pts,wts)
        !
    end select
    !
    if (nloc == 3) then
      !
      ! F2003 AUTO-REALLOC
      return_value%pts = reshape( [-one,chop(pts),+one] , [nd,np+2] )
      return_value%wts = [zero,chop(wts),zero] ! F2003 AUTO-REALLOC
      !
    else
      !
      return_value%pts = reshape( chop(pts) , [nd,np] ) ! F2003 AUTO-REALLOC
      return_value%wts = chop(wts) ! F2003 AUTO-REALLOC
      !
    end if
    !
  end function get_edge_points
  !
  !#############################################################################
  !
#ifndef DISABLE_PURE
  pure &
#endif
  function get_quad_points(edge) result(return_value)
    !
    !.. Use Statements ..
    use order_mod, only : geom_solpts
    !
    !.. Formal Arguments ..
    type(element_t), intent(in) :: edge
    !
    !.. Function Result ..
    type(element_t) :: return_value
    !
    !.. Local Scalars ..
    integer :: i,j,l,nd,np,ierr
    !
    !.. Local Allocatable Arrays ..
    real(mp), allocatable :: pts(:,:)
    !
  continue
    !
    nd = geom_dimen(Geom_Quad)
    np = size(edge%pts,dim=2)**nd
    !
    ! Compute the tensor product of the 1D solution points
    !
    pts = init_array(nd,np,real(zero,kind=kind(pts))) ! F2003 AUTO-REALLOC
    !
    l = 0
    do j = 1,size(edge%pts,dim=2)
      do i = 1,size(edge%pts,dim=2)
        l = l + 1
        pts(1,l) = real( edge%pts(1,i) , kind=kind(pts) )
        pts(2,l) = real( edge%pts(1,j) , kind=kind(pts) )
      end do
    end do
    !
    return_value%pts = chop( pts ) ! F2003 AUTO-REALLOC
    !
    return_value%wts = init_array(np,zero) ! F2003 AUTO-REALLOC
    l = 0
    do j = 1,size(edge%wts)
      do i = 1,size(edge%wts)
        l = l + 1
        return_value%wts(l) = real( l , kind=kind(return_value%wts) )
      end do
    end do
    !
  end function get_quad_points
  !
  !#############################################################################
  !
#ifndef DISABLE_PURE
  pure &
#endif
  function get_hexa_points(edge) result(return_value)
    !
    !.. Use Statements ..
    use order_mod, only : geom_solpts
    !
    !.. Formal Arguments ..
    type(element_t), intent(in) :: edge
    !
    !.. Function Result ..
    type(element_t) :: return_value
    !
    !.. Local Scalars ..
    integer :: i,j,k,l,nd,np,ierr
    !
    !.. Local Allocatable Arrays ..
    real(mp), allocatable :: pts(:,:)
    !
  continue
    !
    nd = geom_dimen(Geom_Hexa)
    np = size(edge%pts,dim=2)**nd
    !
    ! Compute the tensor product of the 1D solution points
    !
    pts = init_array(nd,np,real(zero,kind=kind(pts))) ! F2003 AUTO-REALLOC
    !
    l = 0
    do k = 1,size(edge%pts,dim=2)
      do j = 1,size(edge%pts,dim=2)
        do i = 1,size(edge%pts,dim=2)
          l = l + 1
          pts(1,l) = real( edge%pts(1,i) , kind=kind(pts) )
          pts(2,l) = real( edge%pts(1,j) , kind=kind(pts) )
          pts(3,l) = real( edge%pts(1,k) , kind=kind(pts) )
        end do
      end do
    end do
    !
    return_value%pts = chop( pts ) ! F2003 AUTO-REALLOC
    !
    return_value%wts = init_array(np,zero) ! F2003 AUTO-REALLOC
    l = 0
    do k = 1,size(edge%wts)
      do j = 1,size(edge%wts)
        do i = 1,size(edge%wts)
          l = l + 1
          return_value%wts(l) = real( l , kind=kind(return_value%wts) )
        end do
      end do
    end do
    !
  end function get_hexa_points
  !
  !#############################################################################
  !
  subroutine write_points(elem,geom,order,nloc)
    !
    !.. Formal Arguments ..
    integer,         intent(in) :: geom
    integer,         intent(in) :: order
    integer,         intent(in) :: nloc
    type(element_t), intent(in) :: elem
    !
    !.. Local Scalars ..
    integer :: i,j,k,l,m,n,ierr
    integer :: nc,nc1d,nd,np,np1d,nnpc
    integer :: ncells,totcel,totpts
    !
    character(len=20) :: zonetype
    character(len=40) :: zonename
    !
    !.. Local Allocatable Arrays ..
    integer,  allocatable, dimension(:,:,:) :: ipelem
    real(wp), allocatable, dimension(:,:,:) :: xyz
    real(wp), allocatable, dimension(:,:)   :: xyz_offset
    !
    !.. Local Parameters ..
    character(len=*), parameter :: pname = "write_points"
    character(len=*), parameter :: cloc(1:3) = ["Gauss     ", &
                                                "Lobatto   ", &
                                                "Gauss+Ends"]
    real(wp), parameter :: ds = one/ten
    !
  continue
    !
    nd = size(elem%pts,dim=1)
    np = size(elem%pts,dim=2)
    !
    if (geom == Geom_Quad) then
      !
      np1d = nint( sqrt( real(np,kind=wp) ) )
      if (np1d*np1d /= np) then
        write (error_message,100) "Quad",np,np1d
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
      nc1d = np1d-1
      nc = nc1d*nc1d
      !
      zonetype = "FEQUADRILATERAL"
      ncells = 4
      nnpc = 4
      !
    else if (geom == Geom_Hexa) then
      !
      np1d = nint( (real(np,kind=wp)**one3) )
      if (np1d*np1d*np1d /= np) then
        write (error_message,100) "Hexa",np,np1d
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
      nc1d = np1d-1
      nc = nc1d*nc1d*nc1d
      !
      zonetype = "FEBRICK"
      ncells = 8
      nnpc = 8
      !
    else
      !
      write (error_message,101)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      !
    end if
    !
    totcel = nc*ncells
    totpts = np*ncells
    !
    !
    !
    ipelem     = init_array(nnpc,nc,ncells,0) ! F2003 AUTO-REALLOC
    xyz        = init_array(3,np,ncells,zero) ! F2003 AUTO-REALLOC
    xyz_offset = init_array(3,ncells,zero)    ! F2003 AUTO-REALLOC
    !
    m = 0
    do n = 1,ncells
      !
      l = 0
      do k = 1,merge(nc1d,1,geom == Geom_Hexa)
        do j = 1,nc1d
          do i = 1,nc1d
            !
            l = l + 1
            ipelem(1,l,n) = m + np1d*np1d*(k-1) + np1d*(j-1) + i
            ipelem(2,l,n) = m + np1d*np1d*(k-1) + np1d*(j-1) + i + 1
            ipelem(3,l,n) = m + np1d*np1d*(k-1) + np1d*(j  ) + i + 1
            ipelem(4,l,n) = m + np1d*np1d*(k-1) + np1d*(j  ) + i
            if (geom == Geom_Hexa) then
              ipelem(5,l,n) = m + np1d*np1d*(k  ) + np1d*(j-1) + i
              ipelem(6,l,n) = m + np1d*np1d*(k  ) + np1d*(j-1) + i + 1
              ipelem(7,l,n) = m + np1d*np1d*(k  ) + np1d*(j  ) + i + 1
              ipelem(8,l,n) = m + np1d*np1d*(k  ) + np1d*(j  ) + i
            end if
            !
          end do
        end do
      end do
      m = m + np
    end do
    !
    n = 0
    do k = 1,merge(2,1,geom==Geom_Hexa)
      do j = 1,2
        do i = 1,2
          n = n + 1
          xyz_offset(1,n) = real( 2*i - 3 , kind=wp ) * (one + ds)
          xyz_offset(2,n) = real( 2*j - 3 , kind=wp ) * (one + ds)
          xyz_offset(3,n) = real( 2*k - 3 , kind=wp ) * (one + ds)
        end do
      end do
    end do
    !
    do n = 1,ncells
      do l = 1,np
        xyz(1:nd,l,n) = elem%pts(1:nd,l) + xyz_offset(1:nd,n)
      end do
    end do
    !
    !
    write (zonename,1) order,trim(Geom_Name(geom)),trim(cloc(nloc))
    !
    write (io,2) trim(adjustl(zonename)),totpts,totcel, &
                 trim(adjustl(zonetype))
    !
    do i = 1,3
      write (io,3) ((xyz(i,l,n), l=1,np), n=1,ncells)
    end do
    !
    do n = 1,ncells
      do l = 1,nc
        write (io,4) (ipelem(i,l,n), i=1,nnpc)
      end do
    end do
    !
    !
    if (allocated(ipelem)) then
      deallocate ( ipelem , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"ipelem",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    if (allocated(xyz)) then
      deallocate ( xyz , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"xyz",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    if (allocated(xyz_offset)) then
      deallocate ( xyz_offset , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"xyz_offset",2,__LINE__,__FILE__,ierr, &
                       error_message)
    end if
    !
    ! Format Statements
    !
    1 format ("P",i0," - ",a," at ",a)
    2 format ('ZONE T="',a,'", N=',i0,', E=',i0, &
              ', DATAPACKING=BLOCK, ZONETYPE=',a)
    3 format (5es16.8)
    4 format (8(1x,i6))
    !
    100 format ("  Error finding np1d for ",a,"!",/, &
                "       np   = ",i0,/,"       np1d = ",i0)
    101 format ("  ERROR! Invalid element type!")
    !
  end subroutine write_points
  !
  !#############################################################################
  !
end subroutine dump_solpts
!
!###############################################################################
!
#include "Functions/check_for_bad_number.f90"
#include "Functions/debug_timer.f90"
!
!###############################################################################
!
end module quadrature_mod
!
!###############################################################################
!
