module projection_mod
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
  public :: init_projection_matrices
  public :: project_to_porder
  public :: projection_memory_usage
  !
  !
  !
  ! projection_matrix : derived type that stores the matrices for
  !                     projecting a polynomial quantity within a
  !                     cell to a face point or cell point for a given
  !                     combination of cell geometry and cell order
  !
  type, public :: projection_matrix
    ! limiter : projection operator used to 'limit' a polynomial function
    !           (defined by the nodal values of the polynomial given at the
    !           solution points of a cell) down to a lower degree polynomial
    !           than the cell order
    type(matrix), allocatable :: limiter(:)
    ! CellOI : projection operator used to project a polynomial function,
    !             defined by the nodal values of the polynomial given at a
    !             set of quadrature points in a cell (or face),
    !             down to a polynomial function of the same degree as the
    !             solution polynomial of the cell (or face),
    !             defined by the nodal values of the polynomial given at the
    !             solution points of the cell (or flux points of the face).
    type(matrix), allocatable :: CellOI(:)
    ! FaceOI : projection operator used to project a polynomial function,
    !             defined by the nodal values of the polynomial given at a
    !             set of quadrature points in a cell (or face),
    !             down to a polynomial function of the same degree as the
    !             solution polynomial of the cell (or face),
    !             defined by the nodal values of the polynomial given at the
    !             solution points of the cell (or flux points of the face).
    type(matrix), allocatable :: FaceOI(:)
  end type projection_matrix
  !
  ! projct : projection matrix for all combinations of cell geometry and
  !          and cell order that are possible in the current simulation
  !
  !     For cell 'n' (or face 'nf') :
  !
  !         The projection matrix, giving the contribution of each quadrature
  !         point within cell 'n' (or face 'nf') to the solution point (or
  !         flux point) 'k' that is currently of interest, is given by the
  !         expression
  !
  !
  !     OVER-INTEGRATING :
  !
  !         IF: over-integrating a polynomial function in a cell
  !
  !             projct(this_geom,this_order)%CellOI(quadrature_order)%mat(:,k)
  !
  !               this_geom = cell(n)%geom
  !               this_order = cell(n)%order
  !               quadrature_order = q_order
  !
  !         ELSE IF: over-integrating a polynomial function on a face
  !
  !             projct(this_geom,this_order)%FaceOI(quadrature_order)%mat(:,k)
  !
  !               this_geom = face(nf)%geom
  !               this_order = face(nf)%order
  !               quadrature_order = q_order
  !
  !     PROJECTION LIMITING :
  !
  !         projct(this_geom,this_order)%limiter(lower_order)%mat(:,k)
  !
  !         IF: limiting a polynomial function in a cell
  !
  !           this_geom = cell(n)%geom
  !           this_order = cell(n)%order
  !           lower_order = reduced order desired
  !
  !         ELSE IF: limiting a polynomial function on a face
  !
  !           this_geom = face(nf)%geom
  !           this_order = face(nf)%order
  !           lower_order = reduced order desired
  !
  type(projection_matrix), public, save, target, allocatable :: projct(:,:)
  !
  integer, save :: current_elem_geom
  integer, save :: current_elem_order
  integer, save :: current_proj_order
  integer, save :: iu = 21
  !
contains
!
!###############################################################################
!
subroutine init_projection_matrices()
  !
  !.. Use Statements ..
  use geovar,          only : nr
  use order_mod,       only : geom_solpts
  use order_mod,       only : n_min_geom,n_max_geom
  use order_mod,       only : n_min_order,n_max_order
  use quadrature_mod,  only : geom_is_used
  use quadrature_mod,  only : std_elem,face_elem
  use vandermonde_mod, only : vand
  use vandermonde_mod, only : Compute_Basis_Matrix
  !
  !.. Local Scalars ..
  integer :: gmin,gmax,omin,omax,ierr
  integer :: this_geom,elem_order
  integer :: cq_order,fq_order,lmtr_order
  integer :: nsp,nfp,ncqp,nfqp,nlp
  integer :: qdtr_omin,qdtr_omax
  integer :: lmtr_omin,lmtr_omax
  !
  character(len=200) :: array_name
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: bmat(:,:)
  real(wp), allocatable :: vmat(:,:)
  !
  !.. Local Pointers ..
  type(matrix),            pointer :: proj
  type(projection_matrix), pointer :: this_projct
 !real(wp), contiguous,    pointer :: pmat(:,:)
  real(wp),                pointer :: pmat(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_projection_matrices"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  pmat        => null()
  proj        => null()
  this_projct => null()
  !
  gmin = n_min_geom
  gmax = n_max_geom
  omin = n_min_order
  omax = n_max_order
  !
  ! Allocate the projct array
  !
  allocate ( projct(gmin:gmax,omin:omax) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"projct",1,__LINE__,__FILE__,ierr,error_message)
  !
  geom_loop: do this_geom = gmin,gmax
    !
    current_elem_geom = this_geom
    !
    if (.not. geom_is_used(this_geom)) cycle geom_loop
    if (all(this_geom /= Geom_Valid)) cycle geom_loop
    if (this_geom == Geom_Node) cycle geom_loop
    !
    elem_order_loop: do elem_order = omin,omax
      !
      current_elem_order = elem_order
      !
      ! Assign local pointer as alias to simplify the code
      !
      this_projct => projct(this_geom,elem_order)
      !
      ! Get the number of solution points for this cell
      !
      nsp = geom_solpts(this_geom,elem_order)
      !
      !###################################################
      !###################################################
      !#####  OVER-INTEGRATION PROJECTION OPERATORS  #####
      !###################################################
      !###################################################
      !
      !===========================================================
      !  OVER-INTEGRATION (cellOI) PROJECTION OPERATOR FOR CELLS
      !===========================================================
      !
      ! Get the minimum and maximum over-integration projection
      ! order for this combination of cell geometry and cell order
      !
      qdtr_omin = elem_order+1
      qdtr_omax = omax
      !
      ! Allocate the projct%CellOI component of projct
      !
      allocate ( this_projct%CellOI(qdtr_omin:qdtr_omax) , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),elem_order,"CellOI"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Loop through all the projection orders for this combo
      ! and compute the corresponding projection matrices
      !
      cell_qdtr_order_loop: do cq_order = qdtr_omin,qdtr_omax
        !
        current_proj_order = cq_order
        !
        ! Get the number of quadrature points for this quadrature order
        !
        ncqp = geom_solpts(this_geom,cq_order)
        !
        ! Compute the basis matrix for the quadrature points in the
        ! solution polynomial space
        ! NOTE : Dont use F2003 auto-reallocation. This makes sure the shapes
        !        of the array and the function result are the same
        !
        allocate ( bmat(1:ncqp,1:nsp) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"bmat",1,__LINE__,__FILE__,ierr,error_message)
        !
        bmat(:,:) = Compute_Basis_Matrix(this_geom,cq_order, &
                                         basis_order=elem_order)
        !
        ! Assign local pointer as alias to simplify the code
        !
        proj => this_projct%CellOI(cq_order)
        !
        ! Allocate the projection matrix for this quadrature order
        ! NOTE : proj is actually the transpose of the projection matrix
        ! NOTE : Dont use F2003 auto-reallocation. This makes sure the shapes
        !        of the array and the function result are the same
        !
        allocate ( proj%mat(1:ncqp,1:nsp) , source=zero , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),elem_order,"CellOI",cq_order
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message)
        !
        ! Assign local pointer as alias to simplify the code
        !
        pmat => proj%mat
        !
        ! Compute the projection matrix
        !
        pmat(:,:) = Compute_Projection_Matrix( &
                        basis=bmat, &
                        weights=std_elem(this_geom,cq_order)%wts, &
                        vand=vand(this_geom,elem_order)%modal2nodal%mat)
        !
        if (associated(pmat)) pmat => null()
        if (associated(proj)) proj => null()
        !
        ! Deallocate bmat for the next cell quadrature order
        !
        if (allocated(bmat)) then
          deallocate ( bmat , stat=ierr , errmsg=error_message )
          call alloc_error(pname,"bmat",2,__LINE__,__FILE__,ierr,error_message)
        end if
        !
      end do cell_qdtr_order_loop
      !
      !===========================================================
      !  OVER-INTEGRATION (faceOI) PROJECTION OPERATOR FOR FACES
      !===========================================================
      !
      if (geom_dimen(this_geom) < nr) then
        !
        ! This is also a face geometry in this simulation, so we need to
        ! create face over-integration projection matrices since the locations
        ! of the flux points and the solution points might not be the same
        !
        ! Copy nsp over to nfp, just to prevent confusion here since we are
        ! dealing with a face element
        !
        nfp = nsp
        !
        ! Allocate the projct%FaceOI component of projct
        !
        allocate ( this_projct%FaceOI(qdtr_omin:qdtr_omax) , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),elem_order,"FaceOI"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message)
        !
        ! Loop through all the projection orders for this combo
        ! and compute the corresponding projection matrices
        !
        face_qdtr_order_loop: do fq_order = qdtr_omin,qdtr_omax
          !
          current_proj_order = fq_order
          !
          ! Get the number of quadrature points for this quadrature order
          !
          nfqp = geom_solpts(this_geom,fq_order)
          !
          ! Compute the Vandermonde matrix for the flux points
          ! NOTE : Dont use F2003 auto-reallocation. This makes sure the shapes
          !        of the array and the function result are the same
          !
          allocate ( vmat(1:nfp,1:nfp) , source=zero , &
                     stat=ierr , errmsg=error_message )
          call alloc_error(pname,"vmat",1,__LINE__,__FILE__,ierr,error_message)
          !
          vmat(:,:) = Compute_Basis_Matrix( &
                          this_geom,elem_order, &
                          nodal_points=face_elem(this_geom,elem_order)%pts)
          !
          ! Compute the basis matrix for the quadrature points in the
          ! solution polynomial space
          ! NOTE : Dont use F2003 auto-reallocation. This makes sure the shapes
          !        of the array and the function result are the same
          !
          allocate ( bmat(1:nfqp,1:nfp) , source=zero , &
                     stat=ierr , errmsg=error_message )
          call alloc_error(pname,"bmat",1,__LINE__,__FILE__,ierr,error_message)
          !
          bmat(:,:) = Compute_Basis_Matrix( &
                          this_geom,fq_order,basis_order=elem_order, &
                          nodal_points=face_elem(this_geom,fq_order)%pts)
          !
          ! Assign local pointer as alias to simplify the code
          !
          proj => this_projct%FaceOI(fq_order)
          !
          ! Allocate the projection matrix for this quadrature order
          ! NOTE : proj is actually the transpose of the projection matrix
          ! NOTE : Dont use F2003 auto-reallocation. This makes sure the shapes
          !        of the array and the function result are the same
          !
          allocate ( proj%mat(1:nfqp,1:nfp) , source=zero , &
                     stat=ierr , errmsg=error_message )
          write (array_name,1) Geom_Name(this_geom),elem_order,"FaceOI",fq_order
          call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                           error_message)
          !
          ! Assign local pointer as alias to simplify the code
          !
          pmat => proj%mat
          !
          ! Compute the projection matrix
          !
          pmat(:,:) = Compute_Projection_Matrix( &
                          basis=bmat, &
                          weights=face_elem(this_geom,fq_order)%wts, &
                          vand=vmat)
          !
          if (associated(pmat)) pmat => null()
          if (associated(proj)) proj => null()
          !
          ! Deallocate bmat and vmat for the next face quadrature order
          !
          if (allocated(bmat)) then
            deallocate ( bmat , stat=ierr , errmsg=error_message )
            call alloc_error(pname,"bmat",2,__LINE__,__FILE__,ierr, &
                             error_message)
          end if
          !
          if (allocated(vmat)) then
            deallocate ( vmat , stat=ierr , errmsg=error_message )
            call alloc_error(pname,"vmat",2,__LINE__,__FILE__,ierr, &
                             error_message)
          end if
          !
        end do face_qdtr_order_loop
        !
      end if
      !
      !#####################################################
      !#####################################################
      !#####  LIMITING (limiter) PROJECTION OPERATORS  #####
      !#####################################################
      !#####################################################
      !
      ! Get the minimum and maximum over-integration projection
      ! order for this combination of cell geometry and cell order
      !
      lmtr_omin = 0
      lmtr_omax = elem_order-1
      !
      ! Allocate the projct%limiter component of projct
      !
      allocate ( this_projct%limiter(lmtr_omin:lmtr_omax) , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),elem_order,"limiter"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Loop through all the projection orders for this combo
      ! and compute the corresponding projection matrices
      !
      cell_order_loop: do lmtr_order = lmtr_omin,lmtr_omax
        !
        current_proj_order = lmtr_order
        !
        ! Get the number of modes for the reduce limiter order
        !
        nlp = geom_solpts(this_geom,lmtr_order)
        !
        ! Compute the basis matrix for the solution points in the
        ! reduced order polynomial space
        ! NOTE : Dont use F2003 auto-reallocation. This makes sure the shapes
        !        of the array and the function result are the same
        !
        allocate ( bmat(1:nsp,1:nlp) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"bmat",1,__LINE__,__FILE__,ierr,error_message)
        !
        bmat(:,:) = Compute_Basis_Matrix(this_geom,elem_order, &
                                         basis_order=lmtr_order)
        !
        ! Assign local pointer as alias to simplify the code
        !
        proj => this_projct%limiter(lmtr_order)
        !
        ! Allocate the projection matrix for this cell order
        ! NOTE : proj is actually the transpose of the projection matrix
        ! NOTE : Dont use F2003 auto-reallocation. This makes sure the shapes
        !        of the array and the function result are the same
        !
        allocate ( proj%mat(1:nsp,1:nsp) , source=zero , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),elem_order, &
                             "limiter",lmtr_order
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message)
        !
        ! Assign local pointer as alias to simplify the code
        !
        pmat => proj%mat
        !
        ! Compute the projection matrices depending on geometry type
        !
        pmat(:,:) = Compute_Projection_Matrix( &
                        basis=bmat, &
                        weights=std_elem(this_geom,elem_order)%wts)
        !
        if (associated(pmat)) pmat => null()
        if (associated(proj)) proj => null()
        !
        ! Deallocate bmat for the next limiter order
        !
        if (allocated(bmat)) then
          deallocate ( bmat , stat=ierr , errmsg=error_message )
          call alloc_error(pname,"bmat",2,__LINE__,__FILE__,ierr, &
                           error_message)
        end if
        !
      end do cell_order_loop
      !
      if (associated(pmat       )) pmat        => null()
      if (associated(proj       )) proj        => null()
      if (associated(this_projct)) this_projct => null()
      !
    end do elem_order_loop
    !
  end do geom_loop
  !
  do this_geom = gmin,gmax
    do elem_order = omin,omax
      if (allocated(projct(this_geom,elem_order)%CellOI)) then
        call projct(this_geom,elem_order)%CellOI%check_for_sparsity
      end if
      if (allocated(projct(this_geom,elem_order)%FaceOI)) then
        call projct(this_geom,elem_order)%FaceOI%check_for_sparsity
      end if
      if (allocated(projct(this_geom,elem_order)%limiter)) then
        call projct(this_geom,elem_order)%limiter%check_for_sparsity
      end if
    end do
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("projct(",a,",",i0,")%",a,:,"(",i0,")%mat")
  !
end subroutine init_projection_matrices
!
!###############################################################################
!
subroutine project_to_porder(cell,array,order)
  !
  !... calculates d(utd)/dt
  !
  !.. Use Statements ..
  use order_mod,  only : maxpts,p_order
  use geovar,     only : cell_t
  !
  !.. Formal Arguments ..
  type(cell_t), dimension(:),    intent(in) :: cell
  real(wp),   dimension(:,:), intent(inout) :: array
  integer,          optional,    intent(in) :: order
  !
  !.. Local Scalars ..
  integer :: k,m,n,nc,n1,n2,np,nvar
  integer :: this_geom,this_order,proj_order
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxpts,1:maxpts) :: pmat
  real(wp), dimension(1:maxpts,1:size(array,dim=1)) :: tmp
  !
continue
  !
  proj_order = p_order
  if (present(order)) proj_order = order
  !
  nvar = size(array,dim=1)
  !
  do nc = 1,size(cell)
    !
    n1 = cell(nc)%beg_sp ! beginning index for solpts in cell nc
    n2 = cell(nc)%end_sp ! ending    index for solpts in cell nc
    np = n2-n1+1 ! number of solution points in this cell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    if (proj_order < this_order) then
      !
      ! Store the transpose of usp for this cell to reduce cache misses
      !
      do k = 1,np
        do m = 1,nvar
          tmp(k,m) = array(m,n1-1+k)
        end do
      end do
      !
      ! Project the solution down to proj_order
      ! at each solution point of this cell
      !
      do k = 1,np
        !
        n = k+n1-1 ! total index of solution point within array
        !
        array(1:nvar,n) = projct(this_geom,this_order)% &
                                   limiter(proj_order)% &
                          dot_mat( k , tmp(1:np,1:nvar) )
       !do m = 1,nvar
       !  array(m,n) = projct(this_geom,this_order)% &
       !                        limiter(proj_order)% &
       !               dot( k , tmp(1:np,m) )
       !end do
        !
      end do
      !
    end if
    !
  end do
  !
end subroutine project_to_porder
!
!###############################################################################
!
!!!#define ENABLE_PROJECTION_OUTPUT
#ifndef ENABLE_PROJECTION_OUTPUT
pure &
#endif
function Compute_Projection_Matrix(basis,weights,vand) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: basis(:,:)
  real(wp), intent(in) :: weights(:)
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: vand(:,:)
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j
  !
  !.. Local Arrays ..
  real(qp) :: bmat(1:size(basis,dim=1), &
                   1:size(basis,dim=2))
  real(qp) :: bTxW(1:size(basis,dim=2), &
                   1:size(basis,dim=1))
  real(qp) :: bTxWxb(1:size(basis,dim=2), &
                     1:size(basis,dim=2))
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: pmat(:,:)
  real(qp), allocatable :: vmat(:,:)
  !
continue
  !
  bmat(:,:) = real( basis , kind=kind(bmat) )
  !
#ifdef ENABLE_PROJECTION_OUTPUT
  write (iu,1) "current_elem_geom  = ",current_elem_geom
  write (iu,1) "current_elem_order = ",current_elem_order
  write (iu,1) "current_proj_order = ",current_proj_order
  !
  write (iu,2) "bmat", &
               ((bmat(i,j),i=1,size(bmat,dim=1)),j=1,size(bmat,dim=2))
#endif
  !
  ! bTxW = transpose(basis) x diag(weights)
  !
  do j = 1,size(bmat,dim=2)
    do i = 1,size(bmat,dim=1)
      bTxW(j,i) = bmat(i,j) * real(weights(i),kind=kind(bTxW))
    end do
  end do
  !
  ! Make sure all values in bTxW that are less than
  ! epsilon for working precision are set to zero
  !
  bTxW = chop( bTxW , bTxW )
#ifdef ENABLE_PROJECTION_OUTPUT
  write (iu,2) "bTxW", &
               ((bTxW(i,j),i=1,size(bTxW,dim=1)),j=1,size(bTxW,dim=2))
#endif
  !
  ! bTxWxb = transpose(basis) x diag(weights) x basis
  !
  bTxWxb = matmul( bTxW , bmat )
  !
  ! Make sure all values in bTxWxb that are less than
  ! epsilon for working precision are set to zero
  !
  bTxWxb = chop( bTxWxb , bTxWxb )
#ifdef ENABLE_PROJECTION_OUTPUT
  write (iu,2) "bTxWxb", &
               ((bTxWxb(i,j),i=1,size(bTxWxb,dim=1)),j=1,size(bTxWxb,dim=2))
#endif
  !
  ! Get the inverse of bTxWxb
  ! bTxWxb = inverse( transpose(basis) x diag(weights) x basis )
  !
  bTxWxb = invert_matrix( bTxWxb )
  !
  ! Make sure all values in bTxWxb that are less than
  ! epsilon for working precision are set to zero
  !
  bTxWxb = chop( bTxWxb , bTxWxb )
#ifdef ENABLE_PROJECTION_OUTPUT
  write (iu,2) "inverse(bTxWxb)", &
               ((bTxWxb(i,j),i=1,size(bTxWxb,dim=1)),j=1,size(bTxWxb,dim=2))
#endif
  !
  ! Compute the final projection matrix
  !
  if (present(vand)) then
    !
    ! pmat = vand x bTxWxb x bTxW
    !
    vmat = real( vand , kind=kind(vmat) ) ! USING F2003 AUTO-REALLOCATION
#ifdef ENABLE_PROJECTION_OUTPUT
  write (iu,2) "vand", &
               ((vand(i,j),i=1,size(vand,dim=1)),j=1,size(vand,dim=2))
#endif
    !
    pmat = matmul( vmat , matmul( bTxWxb , bTxW ) ) ! F2003 AUTO-REALLOCATION
    !
  else
    !
    ! pmat = basis x bTxWxb x bTxW
    !
    pmat = matmul( bmat , matmul( bTxWxb , bTxW ) ) ! F2003 AUTO-REALLOCATION
    !
  end if
  !
  ! Make sure all values in pmat that are less than
  ! epsilon for working precision are set to zero
  !
  pmat(:,:) = chop( pmat , pmat )
#ifdef ENABLE_PROJECTION_OUTPUT
  write (iu,2) "pmat", &
               ((pmat(i,j),i=1,size(pmat,dim=1)),j=1,size(pmat,dim=2))
#endif
  !
  ! Return the projection matrix in working precsion
  ! NOTE : The transpose of the projection matrix is returned because the
  !        second dimension (or rows) of the projection matrix is used for
  !        dot products with solution vectors.  This reduces the number of
  !        cache misses whenever the projection matrix is used.
  !
  return_value = chop( transpose( pmat ) ) ! USING F2003 AUTO-REALLOCATION
#ifdef ENABLE_PROJECTION_OUTPUT
  write (iu,2) "return_value",&
               ((return_value(i,j),i=1,size(return_value,dim=1)), &
                                   j=1,size(return_value,dim=2))
  !
  iu = iu + 1
  1 format (a,i0)
  2 format (/,a,/,100(6(1x,g12.4,:),/))
 !2 format (/,a,/,*(6(1x,g12.4,:),/))
#endif
  !
end function Compute_Projection_Matrix
!
!###############################################################################
!
subroutine projection_memory_usage(iunit)
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
  if (allocated(projct)) then
    !
    call dt%set( storage_size(projct,kind=inttype) , &
                         size(projct,kind=inttype) )
    !
    do j = lbound(projct,dim=2),ubound(projct,dim=2)
      do i = lbound(projct,dim=1),ubound(projct,dim=1)
        !
        if (allocated(projct(i,j)%CellOI)) then
          !
          call array%set( storage_size(projct(i,j)%CellOI,kind=inttype) , &
                                  size(projct(i,j)%CellOI,kind=inttype) )
          !
          write (array_name,4) Geom_Name(i),j,"CellOI"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          do k = lbound(projct(i,j)%CellOI,dim=1), &
                 ubound(projct(i,j)%CellOI,dim=1)
            !
            array = projct(i,j)%CellOI(k)%memory_usage()
            !
            write (array_name,4) Geom_Name(i),j,"CellOI",k
            write (iou,2) array_name,array%get_units()
            !
            call dt%add(array)
            !
          end do
          !
        end if
        !
        if (allocated(projct(i,j)%FaceOI)) then
          !
          call array%set( storage_size(projct(i,j)%FaceOI,kind=inttype) , &
                                  size(projct(i,j)%FaceOI,kind=inttype) )
          !
          write (array_name,4) Geom_Name(i),j,"FaceOI"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          do k = lbound(projct(i,j)%FaceOI,dim=1), &
                 ubound(projct(i,j)%FaceOI,dim=1)
            !
            array = projct(i,j)%FaceOI(k)%memory_usage()
            !
            write (array_name,4) Geom_Name(i),j,"FaceOI",k
            write (iou,2) array_name,array%get_units()
            !
            call dt%add(array)
            !
          end do
          !
        end if
        !
        if (allocated(projct(i,j)%limiter)) then
          !
          call array%set( storage_size(projct(i,j)%limiter,kind=inttype) , &
                                  size(projct(i,j)%limiter,kind=inttype) )
          !
          write (array_name,4) Geom_Name(i),j,"limiter"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          do k = lbound(projct(i,j)%limiter,dim=1), &
                 ubound(projct(i,j)%limiter,dim=1)
            !
            array = projct(i,j)%limiter(k)%memory_usage()
            !
            write (array_name,4) Geom_Name(i),j,"limiter",k
            write (iou,2) array_name,array%get_units()
            !
            call dt%add(array)
            !
          end do
          !
        end if
        !
      end do
    end do
    !
    call total%add(dt)
    write (array_name,5) "projct"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "projct"
    !
  end if
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module projection_mod")
  2 format (" Array: ",a,dt)
  3 format (" Total memory = ",f8.2,1x,a)
  4 format ("projct(",a,",",i0,")%",a,:,"(",i0,")%mat")
  5 format ("'",a,"'")
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  !
#endif
end subroutine projection_memory_usage
!
!###############################################################################
!
end module projection_mod
!
!###############################################################################
!
