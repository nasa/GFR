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
  !        this_geom = cell_geom(n)
  !        this_order = cell_order(n)
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
    ! Interpolate to gauss quadrature points
    type(matrix), allocatable :: gauss
    ! Interpolate to lobatto quadrature points
    type(matrix), allocatable :: lobatto
    ! Pointer to quadrature points used for error computations
    real(wp), pointer, contiguous :: error(:,:) => null()
    ! Pointer to quadrature points used to output the solution for post-processing
    real(wp), pointer, contiguous :: output(:,:) => null()
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
  use geovar,         only : cell_geom,cell_order
  use order_mod,      only : geom_solpts,geom_flxpts
  use order_mod,      only : n_min_geom,n_max_geom
  use order_mod,      only : n_min_order,n_max_order
  use order_mod,      only : geom_mio
  use quadrature_mod, only : geom_is_used
  use quadrature_mod, only : std_elem,face_elem
  !
  !.. Local Scalars ..
  integer :: gmin,gmax,omin,omax,ierr
  integer :: this_geom,this_order,face_order
  integer :: this_omin,this_omax
  integer :: nsolpts,nflxpts
  character(len=200) :: array_name
  character(len=200) :: error_message
  !
  !.. Local Pointers ..
  type(interpolation_matrix), pointer :: this_interp
  type(matrix), pointer :: this_face
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
        ! Assign local pointers to these matrices as aliases to simplify the code
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
        if (associated(this_face  )) this_face   => null()
        if (associated(solpts     )) solpts      => null()
        if (associated(flxpts     )) flxpts      => null()
        !
      end do face_order_loop
      !
      ! Remove all associations of the local pointers before continuing
      !
      if (associated(this_interp)) this_interp => null()
      !
    end do order_loop
    !
  end do geom_loop
  !
  1 format ("interp(",a,",",i0,")%",a,:,"(",i0,")%mat")
  2 format (" A grid cell of an unknown geometry type was found!")
  !
end subroutine init_interpolation_matrices
!
!###############################################################################
!
subroutine init_outerpolation_matrices()
end subroutine init_outerpolation_matrices
!
!###############################################################################
!
pure function generic_interpolation_matrix(basis,newbasis,wmat) result(return_value)
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
  integer  :: i,j,n,nf,nfp,nsp
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
pure function InterpToFace_Tria(rs) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Function Return Argument ..
  real(wp), dimension(1:size(rs,2),1:3*np2nfp(size(rs,2))) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,np,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:np2nfp(size(rs,2)),1:3) :: Fmask
  !
continue
  !
  np = size(rs,2)
  nfp = np2nfp( np )
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
  do ifa = 1,nfp1d
    !
    xyz_qp(1) = flxpts_qp(ifa) ! x/xi coordinate of current face/flux point
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
  do jfa = 1,nfp1d
    !
    xyz_qp(2) = flxpts_qp(jfa) ! y/eta coordinate of current face/flux point
    !
    do kfa = 1,nfp1d
      !
      xyz_qp(3) = flxpts_qp(kfa) ! z/zeta coordinate of current face/flux point
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
  do ifa = 1,nfp1d
    !
    xyz_qp(1) = flxpts_qp(ifa) ! x/xi coordinate of current face/flux point
    !
    do kfa = 1,nfp1d
      !
      xyz_qp(3) = flxpts_qp(kfa) ! z/zeta coordinate of current face/flux point
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
include "Functions/np2nfp.f90"
include "Functions/cell_solpts.f90"
!
!###############################################################################
!
end module interpolation_mod
!
!###############################################################################
!
