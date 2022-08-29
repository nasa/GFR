module quadrature_mod
  !
!!!#define PROFILE_ON
!!!#define DISABLE_PURE
  !
  use module_kind_types, global_debug_timer =>      debug_timer, &
                    pure_global_debug_timer => pure_debug_timer
  !
  implicit none
  !
  private
  !
  public :: find_valid_geometries
  public :: init_geom_quadrature_rules
  public :: init_face_quadrature_rules
  public :: quadrature_memory_usage
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
  !  std_elem(this_geom,this_order)%wts(k) :
  !                       The quadrature weight of the solution point 'k' within
  !                       the standard region of the current element geometry.
  !
  type(element_type), public, save, target, allocatable :: std_elem(:,:)
  !
  ! solpts_edge : pointer to the 'n_order' 1D edge solution points
  !
  real(wp), public, protected, pointer :: solpts_edge(:) => null()
  !
  ! face_elem : Array giving the location of the flux points and their
  !             respective quadrature weights for all combinations of cell
  !             geometry and cell order that are possible in the current
  !             simulation
  !
  !  For a given face geometry, this_geom,
  !              and its order, this_order :
  !
  !  face_elem(this_geom,this_order)%pts(l,k) :
  !                       The location of the flux point 'k' in the direction
  !                       'l' of the local coordinate system within the
  !                       standard region of the current face geometry.
  !
  !  face_elem(this_geom,this_order)%wts(k) :
  !                       The quadrature weight of the flux point 'k' within
  !                       the standard region of the current face geometry.
  !
  type(element_type), public, save, target, allocatable :: face_elem(:,:)
  !
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
  character(len=:), parameter :: pname = "find_valid_geometries"
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
  !.. Local Scalars ..
  integer :: n,npts,ndim,ierr
  integer :: this_geom,this_order
  character(len=100) :: array_name
  character(len=300) :: timer_message
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
  call debug_timer(entering_procedure,pname)
  !
 !call dump_solpts
  !
  ! Initialize this_elem to disassociated
  !
  this_elem => null()
  edge_elem => null()
  tria_elem => null()
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
  ! NOTE: The pts array within the type element_type is only needed for these
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
      write (timer_message,100) trim(adjustl(Geom_Name(this_geom))),this_order
      100 format ("Geometry = ",a,"  ;  Order = ",i0)
      call debug_timer(start_timer,timer_message)
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
      ! Allocate the wts component of the current geometry/order combination.
      ! If already allocated, then the rest of this loop was done before
      !
      if (allocated(this_elem%wts)) then
        if (associated(this_elem)) this_elem => null()
        cycle basic_order_loop
      end if
      !
      allocate ( this_elem%wts(1:npts) , stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"wts"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Allocate the pts component of the current geometry/order combination
      !
      allocate ( this_elem%pts(1:ndim,1:npts) , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"pts"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
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
      write (timer_message,100) trim(adjustl(Geom_Name(this_geom))),this_order
      call debug_timer(start_timer,timer_message)
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
      ! Allocate the wts component of the current geometry/order combination.
      ! If already allocated, then the rest of this loop was done before
      !
      if (allocated(this_elem%wts)) then
        if (associated(this_elem)) this_elem => null()
        cycle geom_loop
      end if
      !
      allocate ( this_elem%wts(1:npts) , stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"wts"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Allocate the pts component of the current geometry/order combination
      !
      allocate ( this_elem%pts(1:ndim,1:npts) , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"pts"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Fill in the quadrature points and weights for the
      ! current geometry/order combination
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
          this_elem = get_hexa_quadrature(edge_elem)
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
    !
    ! Disassociate the pointers 'edge_elem' and 'tria_elem' before
    ! continuing to the next order
    !
    if (associated(edge_elem)) edge_elem => null()
    if (associated(tria_elem)) tria_elem => null()
    !
  end do order_loop
  !
  ! Before we leave, assign the pointer solpts_edge to its target
  !
  solpts_edge => std_elem(Geom_Edge,n_order)%pts(1,:)
  !
  call debug_timer(leaving_procedure,pname)
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
  use ovar,      only : loc_flux_pts
  use ovar,      only : loc_triangle_pts
  !
  !.. Local Scalars ..
  integer :: n,npts,ndim,ierr
  integer :: gmin,gmax,this_order
  logical(lk) :: there_are_quad_faces
  logical(lk) :: there_are_tria_faces
  character(len=100) :: array_name
  character(len=300) :: timer_message
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
  call debug_timer(entering_procedure,pname)
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
    write (timer_message,100) trim(adjustl(Geom_Name(Geom_Edge))),this_order
    100 format ("Geometry = ",a,"  ;  Order = ",i0)
    call debug_timer(start_timer,timer_message)
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
    ! Allocate the wts component of the current geometry/order combination.
    ! If already allocated, then the rest of this loop was done before
    !
    if (allocated(this_elem%wts)) then
      if (associated(this_elem)) this_elem => null()
      cycle edge_loop
    end if
    !
    allocate ( this_elem%wts(1:npts) , stat=ierr , errmsg=error_message )
    write (array_name,1) Geom_Name(Geom_Edge),this_order,"wts"
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Allocate the pts component of the current geometry/order combination
    !
    allocate ( this_elem%pts(1:ndim,1:npts) , stat=ierr , errmsg=error_message )
    write (array_name,1) Geom_Name(Geom_Edge),this_order,"pts"
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
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
      write (timer_message,100) trim(adjustl(Geom_Name(Geom_Tria))),this_order
      call debug_timer(start_timer,timer_message)
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
      ! Allocate the wts component of the current geometry/order combination.
      ! If already allocated, then the rest of this loop was done before
      !
      if (allocated(this_elem%wts)) then
        if (associated(this_elem)) this_elem => null()
        cycle tria_loop
      end if
      !
      allocate ( this_elem%wts(1:npts) , stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(Geom_Tria),this_order,"wts"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Allocate the pts component of the current geometry/order combination
      !
      allocate ( this_elem%pts(1:ndim,1:npts) , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(Geom_Tria),this_order,"pts"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
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
      write (timer_message,100) trim(adjustl(Geom_Name(Geom_Quad))),this_order
      call debug_timer(start_timer,timer_message)
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
      ! Allocate the wts component of the current geometry/order combination.
      ! If already allocated, then the rest of this loop was done before
      !
      if (allocated(this_elem%wts)) then
        if (associated(this_elem)) this_elem => null()
        cycle quad_loop
      end if
      !
      allocate ( this_elem%wts(1:npts) , stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(Geom_Quad),this_order,"wts"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
     !!
     !! Allocate the pts component of the current geometry/order combination
     !!
     !allocate ( this_elem%pts(1:ndim,1:npts) , &
     !           stat=ierr , errmsg=error_message )
     !write (array_name,1) Geom_Name(Geom_Quad),this_order,"pts"
     !call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Fill in the solution points and quadrature weights
      ! for a quad face with the current order
      !
      this_elem = get_quad_quadrature(edge_elem)
      !
      call output_quadrature(this_elem,Geom_Quad,this_order,"Face")
      !
      ! Disassociate this_elem before continuing to the next loop combo
      !
      if (associated(this_elem)) this_elem => null()
      if (associated(edge_elem)) edge_elem => null()
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
  1 format ("face_elem(",a,",",i0,")%",a)
  !
end subroutine init_face_quadrature_rules
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_quadrature_wts(this_geom,this_order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: this_geom
  integer, intent(in) :: this_order
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: i,ierr
  !
continue
  !
  ! NEED TO FIGURE OUT WHAT DO IF QUADRATURE POINTS ARE NOT THE SAME AS THOSE
  ! USED TO COMPUTE THE INFORMATION STORED IN THE STD_ELEM ARRAY.
  !
  allocate ( return_value(1:this_order+1) , stat=ierr )
  !
  return_value = one/real(this_order+1,kind=wp)
  !
end function get_quadrature_wts
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
  type(element_type) :: return_value
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
  use triangle_mod, only : get_triangle_weights_at_solpts
  use triangle_mod, only : old_get_triangle_weights_at_solpts
  !
  !.. Formal Arguments ..
  integer, intent(in) :: p
  integer, intent(in) :: location
  !
  !.. Function Result ..
  type(element_type) :: return_value
  !
  !.. Local Allocatable Arrays ..
  real(mp), allocatable :: pts(:,:)
  real(mp), allocatable :: wts(:)
  !
continue
  !
  ! Get the solution points
  !
  pts = TriNodes2D(p,location) ! F2003 AUTO-REALLOC
  !
  ! The method for computing the triangle quadrature weights is
  ! independent of the method for getting the point locations
  !
  wts = get_triangle_weights_at_solpts(pts)     ! F2003 AUTO-REALLOC
 !wts = old_get_triangle_weights_at_solpts(pts) ! F2003 AUTO-REALLOC
  !
 !return_value%pts = chop( pts ) ! F2003 AUTO-REALLOC
 !return_value%wts = chop( wts ) ! F2003 AUTO-REALLOC
  !
  return_value%pts = real(pts,kind=kind(return_value%pts)) ! F2003 AUTO-REALLOC
  return_value%wts = real(wts,kind=kind(return_value%wts)) ! F2003 AUTO-REALLOC
  !
end function get_tria_quadrature
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_quad_quadrature(edge_elem) result(return_value)
  !
  use order_mod, only : geom_solpts
  !
  type(element_type), intent(in) :: edge_elem
  !
  type(element_type) :: return_value
  !
  integer :: i,j,l,n,ndim,npts,ierr
  !
  real(mp), allocatable :: pts(:,:)
  real(mp), allocatable :: wts(:)
  !
continue
  !
  n = size(edge_elem%wts) - 1
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
  use tetrahedra_mod, only : get_tetrahedron_weights_at_solpts
  !
  !.. Formal Arguments ..
  integer, intent(in) :: p
  integer, intent(in) :: location
  !
  !.. Function Result ..
  type(element_type) :: return_value
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
  ! The method for computing the tetrahedra quadrature weights is
  ! independent of the method for getting the point locations
  !
  wts = get_tetrahedron_weights_at_solpts(pts)     ! F2003 AUTO-REALLOC
  !
 !return_value%pts = chop( pts ) ! F2003 AUTO-REALLOC
 !return_value%wts = chop( wts ) ! F2003 AUTO-REALLOC
  !
  return_value%pts = real(pts,kind=kind(return_value%pts)) ! F2003 AUTO-REALLOC
  return_value%wts = real(wts,kind=kind(return_value%wts)) ! F2003 AUTO-REALLOC
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
  use pyramid_mod, only : get_pyramid_weights_at_solpts
  !
  !.. Formal Arguments ..
  integer, intent(in) :: p
  integer, intent(in) :: location
  !
  !.. Function Result ..
  type(element_type) :: return_value
  !
  !.. Local Scalars ..
  integer :: node_loc
  !
  !.. Local Allocatable Arrays ..
  real(mp), allocatable :: pts(:,:)
  real(mp), allocatable :: wts(:)
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
  ! The method for computing the pyramid quadrature weights is
  ! independent of the method for getting the point locations
  !
  wts = get_pyramid_weights_at_solpts(pts)     ! F2003 AUTO-REALLOC
  !
 !return_value%pts = chop( pts ) ! F2003 AUTO-REALLOC
 !return_value%wts = chop( wts ) ! F2003 AUTO-REALLOC
  !
  return_value%pts = real(pts,kind=kind(return_value%pts)) ! F2003 AUTO-REALLOC
  return_value%wts = real(wts,kind=kind(return_value%wts)) ! F2003 AUTO-REALLOC
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
  use prism_mod, only : get_prism_weights_at_solpts
  !
  !.. Formal Arguments ..
  integer, intent(in) :: p
  integer, intent(in) :: tria_loc
  integer, intent(in) :: edge_loc
  !
  !.. Function Result ..
  type(element_type) :: return_value
  !
  !.. Local Allocatable Arrays ..
  real(mp), allocatable :: pts(:,:)
  real(mp), allocatable :: wts(:)
  !
continue
  !
  ! Get the solution points
  !
  pts = PriNodes3D(p,tria_loc,edge_loc) ! F2003 AUTO-REALLOC
 !pts = PriNodes3D(p,tria_loc,edge_loc,tria_edge_loc) ! F2003 AUTO-REALLOC
  !
  ! The method for computing the prism quadrature weights is
  ! independent of the method for getting the point locations
  !
  wts = get_prism_weights_at_solpts(pts)     ! F2003 AUTO-REALLOC
  !
 !return_value%pts = chop( pts ) ! F2003 AUTO-REALLOC
 !return_value%wts = chop( wts ) ! F2003 AUTO-REALLOC
  !
  return_value%pts = real(pts,kind=kind(return_value%pts)) ! F2003 AUTO-REALLOC
  return_value%wts = real(wts,kind=kind(return_value%wts)) ! F2003 AUTO-REALLOC
  !
end function get_pris_quadrature
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_hexa_quadrature(edge_elem) result(return_value)
  !
  use order_mod, only : geom_solpts
  !
  type(element_type), intent(in) :: edge_elem
  !
  type(element_type) :: return_value
  !
  integer :: i,j,k,l,n,ndim,npts,ierr
  !
  real(mp), allocatable :: pts(:,:)
  real(mp), allocatable :: wts(:)
  !
continue
  !
  n = size(edge_elem%wts) - 1
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
end function get_hexa_quadrature
!
!###############################################################################
!
subroutine output_quadrature(this_elem,this_geom,this_order,this_type)
  !
  !.. Use Statements ..
  use ovar, only : mout => matrix_output_unit
  !
  !.. Formal Arguments ..
  class(element_type), intent(in) :: this_elem
  integer,             intent(in) :: this_geom
  integer,             intent(in) :: this_order
  character(len=*),    intent(in) :: this_type
  !
  !.. Local Scalars ..
  integer :: i,j,n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "output_quadrature"
  !
continue
  !
  if (mout == iout) return
  !
  write (mout,1) this_order, trim(adjustl(Geom_Full_Name(this_geom))), &
                 trim(adjustl(this_type))
  !
  n1 = size(this_elem%pts,dim=1)
  n2 = size(this_elem%pts,dim=2)
  !
  if (geom_dimen(this_geom) == 1) then
    do j = 1,n2
      write (mout,2) (this_elem%pts(i,j),i=1,n1), this_elem%wts(j)
    end do
  else if (geom_dimen(this_geom) == 2) then
    do j = 1,n2
      write (mout,3) (this_elem%pts(i,j),i=1,n1), this_elem%wts(j)
    end do
  else if (geom_dimen(this_geom) == 3) then
    do j = 1,n2
      write (mout,4) (this_elem%pts(i,j),i=1,n1), this_elem%wts(j)
    end do
  end if
  !
  1 format (/, &
            80("#"),/, &
            3x,"Solution Point Coordinates : Quadrature Weight",/, &
            3x,"for a P",i0,1x,a,1x,a," element",/, &
            80("-"),/)
  2 format (1x,es14.7,"  :  ",es14.7)
  3 format (2(1x,es14.7),"  :  ",es14.7)
  4 format (3(1x,es14.7),"  :  ",es14.7)

  !
end subroutine output_quadrature
!
!###############################################################################
!
subroutine quadrature_memory_usage(iunit)
  !
  !.. Use Statements ..
  use module_memory
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
  integer :: npts,ndim,ierr,location,io
  integer :: this_geom,this_order
  character(len=100) :: array_name
  !
  !.. Local Pointers ..
  type(element_type), pointer :: this_elem
  type(element_type), pointer :: ref_edge
  !
  !.. Local Allocatable Arrays ..
  type(element_type), target, allocatable, dimension(:,:) :: edge_elem
  type(element_type), target, allocatable, dimension(:,:) :: quad_elem
  type(element_type), target, allocatable, dimension(:,:) :: hexa_elem
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
  write (io,10)
  !
  ! Initialize this_elem to disassociated
  !
  this_elem => null()
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
  ! NOTE: The pts array within the type element_type is only needed for these
  !       two geometries (maybe tetrahedra as well) so they will be the only
  !       ones in which this array will be allocated and created.
  !
  this_geom = Geom_Edge
  !
  do location = 1,3
    !
    basic_order_loop: do this_order = min_order,max_order
      !
      ! Assign the pointer 'this_elem' to std_elem(this_geom,this_order)
      ! to simplify the remaining code in this loop and make it more legible
      !
      this_elem => edge_elem(this_order,location)
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
      if (location == 3) npts = npts + 2
      !
      ! Allocate the wts component of the current geometry/order combination
      !
      allocate ( this_elem%wts(1:npts) , stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"wts"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Allocate the pts component of the current geometry/order combination
      !
      allocate ( this_elem%pts(1:ndim,1:npts) , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"pts"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Fill in the solution points and quadrature weights
      ! for an edge with the current order
      !
      this_elem = get_edge_points(this_order,location)
      !
      ! Disassociate this_elem before continuing to the next loop combo
      !
      if (associated(this_elem)) this_elem => null()
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
      ! Assign pointers to the edge and triangle elements of the std_elem
      ! array in order to simplify the remaining code in this loop and
      ! make it more legible.
      !
      ref_edge => edge_elem(this_order,location)
      !
      ! Assign the pointer 'this_elem' to std_elem(this_geom,this_order)
      ! to simplify the remaining code in this loop and make it more legible
      !
      this_elem => quad_elem(this_order,location)
      !
      ! Spatial dimension of the current geometry
      ! NOTE: Make sure that ndim is at minimum 1, primarily for Geom_Node
      !
      ndim = max( 1 , geom_dimen(this_geom) )
      !
      ! Number of solution points for the current geometry/order combination
      !
      npts = size(ref_edge%wts)**ndim
      !
      ! Allocate the wts component of the current geometry/order combination
      !
      allocate ( this_elem%wts(1:npts) , stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"wts"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Allocate the pts component of the current geometry/order combination
      !
      allocate ( this_elem%pts(1:ndim,1:npts) , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"pts"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Fill in the quadrature weights for the
      ! current geometry/order combination
      !
      this_elem = get_quad_points(ref_edge)
      !
      ! Disassociate the pointer 'this_elem' before
      ! continuing to the next geometry
      !
      if (associated(this_elem)) this_elem => null()
      !
      ! Disassociate the pointers 'edge_elem' and 'tria_elem' before
      ! continuing to the next order
      !
      if (associated(ref_edge)) ref_edge => null()
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
      ! Assign pointers to the edge and triangle elements of the std_elem
      ! array in order to simplify the remaining code in this loop and
      ! make it more legible.
      !
      ref_edge => edge_elem(this_order,location)
      !
      ! Assign the pointer 'this_elem' to std_elem(this_geom,this_order)
      ! to simplify the remaining code in this loop and make it more legible
      !
      this_elem => hexa_elem(this_order,location)
      !
      ! Spatial dimension of the current geometry
      ! NOTE: Make sure that ndim is at minimum 1, primarily for Geom_Node
      !
      ndim = max( 1 , geom_dimen(this_geom) )
      !
      ! Number of solution points for the current geometry/order combination
      !
      npts = size(ref_edge%wts)**ndim
      !
      ! Allocate the wts component of the current geometry/order combination
      !
      allocate ( this_elem%wts(1:npts) , stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"wts"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Allocate the pts component of the current geometry/order combination
      !
      allocate ( this_elem%pts(1:ndim,1:npts) , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"pts"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Fill in the quadrature weights for the
      ! current geometry/order combination
      !
      this_elem = get_hexa_points(ref_edge)
      !
      ! Disassociate the pointer 'this_elem' before
      ! continuing to the next geometry
      !
      if (associated(this_elem)) this_elem => null()
      !
      ! Disassociate the pointers 'edge_elem' and 'tria_elem' before
      ! continuing to the next order
      !
      if (associated(ref_edge)) ref_edge => null()
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
  1 format ("std_elem(",a,",",i0,")%",a)
  2 format (" A grid cell of an unknown geometry type was found!")
  3 format (" A grid cell of geometry type '",a,"' was encountered while", &
            " trying to compute the quadrature rules for each geometry type.", &
            " This geometry type is currently invalid because it has not", &
            " yet been fully implemented.")
  !
  10 format ('VARIABLES = "X","Y","Z"')
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
    type(element_type) :: return_value
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
    type(element_type), intent(in) :: edge
    !
    !.. Function Result ..
    type(element_type) :: return_value
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
    type(element_type), intent(in) :: edge
    !
    !.. Function Result ..
    type(element_type) :: return_value
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
    integer,            intent(in) :: geom
    integer,            intent(in) :: order
    integer,            intent(in) :: nloc
    type(element_type), intent(in) :: elem
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
    allocate ( ipelem(1:nnpc,1:nc,1:ncells) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ipelem",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( xyz(1:3,1:np,1:ncells) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( xyz_offset(1:3,1:ncells) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_offset",1,__LINE__,__FILE__,ierr,error_message)
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
