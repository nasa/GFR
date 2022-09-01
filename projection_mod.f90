module projection_mod
  !
!!!#define PROFILE_ON
  !
  !.. Use Statements ..
  use module_kind_types, global_debug_timer =>      debug_timer, &
                    pure_global_debug_timer => pure_debug_timer
  !
  use matrix_types_mod, only : matrix => matrix_gpu_t
  use matrix_types_mod, only : create_matrix_sparsity
  !
  use eqn_idx, only : ncs,nec,nmx,nmy,nmz,nmb,nme,nee,nq
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
  public :: move_projct_to_device
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
  type(projection_matrix), public, save, allocatable :: projct(:,:)
  !
  integer, save :: current_elem_geom
  integer, save :: current_elem_order
  integer, save :: current_proj_order
  integer, save :: iu = 21
  !
  !
  ! #########################################
  ! #####   MODULE GENERIC INTERFACES   #####
  ! #########################################
  !
  interface project_to_porder
    module procedure project_to_porder_vector, &
                     project_to_porder_array
  end interface project_to_porder
  !
contains
!
!###############################################################################
!
subroutine init_projection_matrices()
  !
  !.. Use Statements ..
  use geovar,          only : nr
  use order_mod,       only : geom_solpts,p_order
  use order_mod,       only : n_min_geom,n_max_geom
  use order_mod,       only : n_min_order,n_max_order
  use quadrature_mod,  only : geom_is_used
  use quadrature_mod,  only : std_elem,face_elem
  use vandermonde_mod, only : vand
  use vandermonde_mod, only : Compute_Basis_Matrix
  use ovar,            only : Limiter_Option
  use ovar,            only : Limiter_Filter_Option
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
  character(len=300) :: timer_message
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: bmat(:,:)
  real(wp), allocatable :: vmat(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_projection_matrices"
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
      allocate ( projct(this_geom,elem_order)%CellOI(qdtr_omin:qdtr_omax) , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),elem_order,"CellOI"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Loop through all the projection orders for this combo
      ! and compute the corresponding projection matrices
      !
      cell_qdtr_order_loop: do cq_order = qdtr_omin,qdtr_omax
        write (timer_message,100) trim(adjustl(Geom_Name(this_geom))), &
                                  elem_order,cq_order
        100 format ("CellOI Geometry = ",a,"  ;  Order = ",i0, &
                    "  ;  cq_order = ",i0)
        call debug_timer(start_timer,timer_message)
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
        ! Allocate the projection matrix for this quadrature order
        ! NOTE : proj is actually the transpose of the projection matrix
        ! NOTE : Dont use F2003 auto-reallocation. This makes sure the shapes
        !        of the array and the function result are the same
        !
        allocate ( projct(this_geom,elem_order)% &
                               CellOI(cq_order)%mat(1:ncqp,1:nsp) , &
                   source=zero , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),elem_order,"CellOI",cq_order
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message)
        !
        ! Compute the projection matrix
        !
        projct(this_geom,elem_order)%CellOI(cq_order)%mat(:,:) = &
                  Compute_Projection_Matrix( &
                        basis=bmat, &
                        weights=std_elem(this_geom,cq_order)%wts, &
                        vand=vand(this_geom,elem_order)%modal2nodal%mat)
        !
        ! Deallocate bmat for the next cell quadrature order
        !
        if (allocated(bmat)) then
          deallocate ( bmat , stat=ierr , errmsg=error_message )
          call alloc_error(pname,"bmat",2,__LINE__,__FILE__,ierr,error_message)
        end if
        call debug_timer(stop_timer,timer_message)
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
        allocate ( projct(this_geom,elem_order)%FaceOI(qdtr_omin:qdtr_omax) , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),elem_order,"FaceOI"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message)
        !
        ! Loop through all the projection orders for this combo
        ! and compute the corresponding projection matrices
        !
        face_qdtr_order_loop: do fq_order = qdtr_omin,qdtr_omax
          write (timer_message,101) trim(adjustl(Geom_Name(this_geom))), &
                                    elem_order,fq_order
          101 format ("FaceOI Geometry = ",a,"  ;  Order = ",i0, &
                      "  ;  fq_order = ",i0)
          call debug_timer(start_timer,timer_message)
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
          ! Allocate the projection matrix for this quadrature order
          ! NOTE : proj is actually the transpose of the projection matrix
          ! NOTE : Dont use F2003 auto-reallocation. This makes sure the shapes
          !        of the array and the function result are the same
          !
          allocate ( projct(this_geom,elem_order)% &
                                 FaceOI(fq_order)%mat(1:nfqp,1:nfp) , &
                     source=zero , stat=ierr , errmsg=error_message )
          write (array_name,1) Geom_Name(this_geom),elem_order,"FaceOI",fq_order
          call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                           error_message)
          !
          ! Compute the projection matrix
          !
          projct(this_geom,elem_order)%FaceOI(fq_order)%mat(:,:) = &
                    Compute_Projection_Matrix( &
                          basis=bmat, &
                          weights=face_elem(this_geom,fq_order)%wts, &
                          vand=vmat)
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
          call debug_timer(stop_timer,timer_message)
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
      if (Limiter_Option == 0) then
        lmtr_omin = elem_order-1
      else if (Limiter_Option >= 1) then
        if (Limiter_Filter_Option /= 0) then
          lmtr_omin = elem_order-1
        end if
      end if
      lmtr_omin = max( 0 , min(lmtr_omin,p_order) )
      !
      lmtr_omax = elem_order-1
      !
      ! Allocate the projct%limiter component of projct
      !
      allocate ( projct(this_geom,elem_order)%limiter(lmtr_omin:lmtr_omax) , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),elem_order,"limiter"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Loop through all the projection orders for this combo
      ! and compute the corresponding projection matrices
      !
      cell_order_loop: do lmtr_order = lmtr_omin,lmtr_omax
        write (timer_message,102) trim(adjustl(Geom_Name(this_geom))), &
                                  elem_order,lmtr_order
        102 format ("Projection Limiter Basis Geometry = ",a, &
                    "  ;  Order = ",i0, &
                    "  ;  lmtr_order = ",i0)
        call debug_timer(start_timer,timer_message)
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
        call debug_timer(stop_timer,timer_message)
        write (timer_message,103) trim(adjustl(Geom_Name(this_geom))), &
                                  elem_order,lmtr_order
        103 format ("Projection Limiter Matrix Geometry = ",a, &
                    "  ;  Order = ",i0, &
                    "  ;  lmtr_order = ",i0)
        call debug_timer(start_timer,timer_message)
        !
        ! Allocate the projection matrix for this cell order
        ! NOTE : proj is actually the transpose of the projection matrix
        ! NOTE : Dont use F2003 auto-reallocation. This makes sure the shapes
        !        of the array and the function result are the same
        !
        allocate ( projct(this_geom,elem_order)% &
                            limiter(lmtr_order)%mat(1:nsp,1:nsp) , &
                   source=zero , stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),elem_order, &
                             "limiter",lmtr_order
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message)
        !
        ! Compute the projection matrices depending on geometry type
        !
        projct(this_geom,elem_order)%limiter(lmtr_order)%mat(:,:) = &
                  Compute_Projection_Matrix( &
                        basis=bmat, &
                        weights=std_elem(this_geom,elem_order)%wts)
        call debug_timer(stop_timer,timer_message)
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
    end do elem_order_loop
    !
  end do geom_loop
  !
  do this_geom = gmin,gmax
    do elem_order = omin,omax
      if (allocated(projct(this_geom,elem_order)%CellOI)) then
        write (timer_message,200) trim(adjustl(Geom_Name(this_geom))),elem_order
        200 format ("CellOI sparsity check:  Geometry = ",a,"  ;  Order = ",i0)
        call debug_timer(start_timer,timer_message)
        call create_matrix_sparsity( projct(this_geom,elem_order)%CellOI )
        call debug_timer(stop_timer,timer_message)
      end if
      if (allocated(projct(this_geom,elem_order)%FaceOI)) then
        write (timer_message,201) trim(adjustl(Geom_Name(this_geom))),elem_order
        201 format ("FaceOI sparsity check:  Geometry = ",a,"  ;  Order = ",i0)
        call debug_timer(start_timer,timer_message)
        call create_matrix_sparsity( projct(this_geom,elem_order)%FaceOI )
        call debug_timer(stop_timer,timer_message)
      end if
      if (allocated(projct(this_geom,elem_order)%limiter)) then
        write (timer_message,202) trim(adjustl(Geom_Name(this_geom))),elem_order
        202 format ("Limiter sparsity check:  Geometry = ",a,"  ;  Order = ",i0)
        call debug_timer(start_timer,timer_message)
        call create_matrix_sparsity( projct(this_geom,elem_order)%limiter )
        call debug_timer(stop_timer,timer_message)
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
subroutine project_to_porder_vector(vector,proj_order)
  !
  !.. Use Statements ..
  use order_mod,  only : maxSP,p_order
  use geovar,     only : ncell,cell
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: vector(:)
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: proj_order
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order,order
  real(wp) :: a,val
  !
  !.. Local Arrays ..
  real(wp) :: old_vector(1:maxSP)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "project_to_porder_vector"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  order = p_order
  if (present(proj_order)) then
    order = proj_order
  end if
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(vector,cell,projct) &
  !$ACC&              PRIVATE(old_vector)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    if (order < this_order) then
      !
      ! Save the current, pre-projection section of vector
      ! for this cell into the temporary array old_vector.
      !
      !$ACC LOOP VECTOR
      do n = 1,np
        old_vector(n) = vector(n0+n)
      end do
      !
      ! Project the vector down to order at each solution point of this cell
      !
      !$ACC LOOP VECTOR
      do n = 1,np
        !
        val = zero
        !
        do i = projct(this_geom,this_order)%limiter(order)%row_ptr(n)+1, &
               projct(this_geom,this_order)%limiter(order)%row_ptr(n+1)
          !
          k = projct(this_geom,this_order)%limiter(order)%row_idx(i)
          a = projct(this_geom,this_order)%limiter(order)%row_val(i)
          !
          val = val + a * old_vector(k)
          !
        end do
        !
        vector(n0+n) = val
        !
      end do
      !
    end if
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine project_to_porder_vector
!
!###############################################################################
!
subroutine project_to_porder_array(array,proj_order,update_gpsp)
  !
  !.. Use Statements ..
  use order_mod, only : p_order
  use geovar,    only : n_solpts
  use ovar,      only : force_nd_subroutines
  use flowvar,   only : gpsp
  !
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: array(:,:)
  !
  !.. Optional Arguments ..
  integer,     optional, intent(in) :: proj_order
  logical(lk), optional, intent(in) :: update_gpsp
  !
  !.. Local Scalars ..
  integer     :: nr,nvar,order
  logical(lk) :: gpsp_update_requested
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "project_to_porder_array"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  gpsp_update_requested = fals
  if (present(update_gpsp)) then
    gpsp_update_requested = update_gpsp
  end if
  !
  nr = nme-nmb+1
  nvar = size(array,dim=2)
  !
  order = p_order
  if (present(proj_order)) then
    order = proj_order
  end if
  !
  if (force_nd_subroutines) then
    !
    call project_to_porder_nd(1,nvar,array,order)
    !
  else
    !
    if (nr == 3) then
      call project_to_porder_3d(array,order)
    else if (nr == 2) then
      call project_to_porder_2d(array,order)
    else
      call project_to_porder_1d(array,order)
    end if
    !
    if (ncs > 0) then
      call project_to_porder_nd(1,ncs,array,order)
    end if
    !
  end if
  !
  if (gpsp_update_requested) then
    call update_thermo_properties(1,n_solpts,array,gpsp, &
                                  loc="cell",caller_name=pname)
    call update_transport_properties(1,n_solpts,array,gpsp)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine project_to_porder_array
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine project_to_porder_3d(array,order)
  !
  !.. Use Statements ..
  use order_mod, only : maxSP
  use geovar,    only : ncell,cell
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: order
  real(wp), intent(inout) :: array(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,var_ec,var_mx,var_my,var_mz,var_ee
  !
  !.. Local Arrays ..
  real(wp) :: old_array(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "project_to_porder_3d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(array,cell,projct) &
  !$ACC&              PRIVATE(old_array)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    if (order < this_order) then
      !
      ! Save the current, pre-projection section of array
      ! for this cell into the temporary array old_array.
      !
      !$ACC LOOP VECTOR
      do n = 1,np
        old_array(n,nec) = array(n0+n,nec)
        old_array(n,nmx) = array(n0+n,nmx)
        old_array(n,nmy) = array(n0+n,nmy)
        old_array(n,nmz) = array(n0+n,nmz)
        old_array(n,nee) = array(n0+n,nee)
      end do
      !
      ! Project the variables to order at each solution point of this cell
      !
      !$ACC LOOP VECTOR
      do n = 1,np
        !
        var_ec = zero
        var_mx = zero
        var_my = zero
        var_mz = zero
        var_ee = zero
        !
        do i = projct(this_geom,this_order)%limiter(order)%row_ptr(n)+1, &
               projct(this_geom,this_order)%limiter(order)%row_ptr(n+1)
          !
          k = projct(this_geom,this_order)%limiter(order)%row_idx(i)
          a = projct(this_geom,this_order)%limiter(order)%row_val(i)
          !
          var_ec = var_ec + a * old_array(k,nec)
          var_mx = var_mx + a * old_array(k,nmx)
          var_my = var_my + a * old_array(k,nmy)
          var_mz = var_mz + a * old_array(k,nmz)
          var_ee = var_ee + a * old_array(k,nee)
          !
        end do
        !
        array(n0+n,nec) = var_ec
        array(n0+n,nmx) = var_mx
        array(n0+n,nmy) = var_my
        array(n0+n,nmz) = var_mz
        array(n0+n,nee) = var_ee
        !
      end do
      !
    end if
    !
  end do
  !
end subroutine project_to_porder_3d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine project_to_porder_2d(array,order)
  !
  !.. Use Statements ..
  use order_mod, only : maxSP
  use geovar,    only : ncell,cell
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: order
  real(wp), intent(inout) :: array(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,var_ec,var_mx,var_my,var_ee
  !
  !.. Local Arrays ..
  real(wp) :: old_array(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "project_to_porder_2d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(array,cell,projct) &
  !$ACC&              PRIVATE(old_array)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    if (order < this_order) then
      !
      ! Save the current, pre-projection section of array
      ! for this cell into the temporary array old_array.
      !
      !$ACC LOOP VECTOR
      do n = 1,np
        old_array(n,nec) = array(n0+n,nec)
        old_array(n,nmx) = array(n0+n,nmx)
        old_array(n,nmy) = array(n0+n,nmy)
        old_array(n,nee) = array(n0+n,nee)
      end do
      !
      ! Project the variables to order at each solution point of this cell
      !
      !$ACC LOOP VECTOR
      do n = 1,np
        !
        var_ec = zero
        var_mx = zero
        var_my = zero
        var_ee = zero
        !
        do i = projct(this_geom,this_order)%limiter(order)%row_ptr(n)+1, &
               projct(this_geom,this_order)%limiter(order)%row_ptr(n+1)
          !
          k = projct(this_geom,this_order)%limiter(order)%row_idx(i)
          a = projct(this_geom,this_order)%limiter(order)%row_val(i)
          !
          var_ec = var_ec + a * old_array(k,nec)
          var_mx = var_mx + a * old_array(k,nmx)
          var_my = var_my + a * old_array(k,nmy)
          var_ee = var_ee + a * old_array(k,nee)
          !
        end do
        !
        array(n0+n,nec) = var_ec
        array(n0+n,nmx) = var_mx
        array(n0+n,nmy) = var_my
        array(n0+n,nee) = var_ee
        !
      end do
      !
    end if
    !
  end do
  !
end subroutine project_to_porder_2d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine project_to_porder_1d(array,order)
  !
  !.. Use Statements ..
  use order_mod, only : maxSP
  use geovar,    only : ncell,cell
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: order
  real(wp), intent(inout) :: array(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,var_ec,var_mx,var_ee
  !
  !.. Local Arrays ..
  real(wp) :: old_array(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "project_to_porder_1d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(array,cell,projct) &
  !$ACC&              PRIVATE(old_array)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    if (order < this_order) then
      !
      ! Save the current, pre-projection section of array
      ! for this cell into the temporary array old_array.
      !
      !$ACC LOOP VECTOR
      do n = 1,np
        old_array(n,nec) = array(n0+n,nec)
        old_array(n,nmx) = array(n0+n,nmx)
        old_array(n,nee) = array(n0+n,nee)
      end do
      !
      ! Project the variables to order at each solution point of this cell
      !
      !$ACC LOOP VECTOR
      do n = 1,np
        !
        var_ec = zero
        var_mx = zero
        var_ee = zero
        !
        do i = projct(this_geom,this_order)%limiter(order)%row_ptr(n)+1, &
               projct(this_geom,this_order)%limiter(order)%row_ptr(n+1)
          !
          k = projct(this_geom,this_order)%limiter(order)%row_idx(i)
          a = projct(this_geom,this_order)%limiter(order)%row_val(i)
          !
          var_ec = var_ec + a * old_array(k,nec)
          var_mx = var_mx + a * old_array(k,nmx)
          var_ee = var_ee + a * old_array(k,nee)
          !
        end do
        !
        array(n0+n,nec) = var_ec
        array(n0+n,nmx) = var_mx
        array(n0+n,nee) = var_ee
        !
      end do
      !
    end if
    !
  end do
  !
end subroutine project_to_porder_1d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine project_to_porder_nd(mbeg,mend,array,order)
  !
  !.. Use Statements ..
  use order_mod, only : maxSP
  use geovar,    only : ncell,cell
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: mbeg,mend
  integer,     intent(in) :: order
  real(wp), intent(inout) :: array(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,m,n,nc,np,n0,nvar
  integer  :: this_geom,this_order
  real(wp) :: a,val
  !
  !.. Local Arrays ..
  real(wp) :: old_array(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "project_to_porder_nd"
  !
continue
  !
  nvar = size(array,dim=2)
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(array,cell,projct) &
  !$ACC&              PRIVATE(old_array)
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    if (order < this_order) then
      !
      ! Save the current, pre-projection section of array
      ! for this cell into the temporary array old_array.
      !
      !$ACC LOOP VECTOR COLLAPSE(2)
      do m = mbeg,mend
        do n = 1,np
          old_array(n,m) = array(n0+n,m)
        end do
      end do
      !
      ! Project the variables to order at each solution point of this cell
      !
      !$ACC LOOP VECTOR COLLAPSE(2)
      do m = mbeg,mend
        do n = 1,np
          !
          val = zero
          !
          do i = projct(this_geom,this_order)%limiter(order)%row_ptr(n)+1, &
                 projct(this_geom,this_order)%limiter(order)%row_ptr(n+1)
            !
            k = projct(this_geom,this_order)%limiter(order)%row_idx(i)
            a = projct(this_geom,this_order)%limiter(order)%row_val(i)
            !
            val = val + a * old_array(k,m)
            !
          end do
          !
          array(n0+n,m) = val
          !
        end do
      end do
      !
    end if
    !
  end do
  !
end subroutine project_to_porder_nd
!
!###############################################################################
!
!!!#define ENABLE_PROJECTION_OUTPUT
#ifndef ENABLE_PROJECTION_OUTPUT
pure &
#endif
function Compute_Projection_Matrix(basis,weights,vand) result(return_value)
  !
#ifdef SPECIAL_FOR_INTEL
  !.. Use Statements ..
  use blas95, only : gemm
#endif
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
  real(mp) :: bmat(1:size(basis,dim=1), &
                   1:size(basis,dim=2))
  real(mp) :: bTxW(1:size(basis,dim=2), &
                   1:size(basis,dim=1))
  real(mp) :: bTxWxb(1:size(basis,dim=2), &
                     1:size(basis,dim=2))
  !
  !.. Local Allocatable Arrays ..
  real(mp), allocatable :: pmat(:,:)
  real(mp), allocatable :: vmat(:,:)
  real(wp), allocatable :: tmp(:,:)
  real(wp), allocatable :: tmp_mat(:,:)
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
#ifdef SPECIAL_FOR_INTEL
  tmp_mat = init_array(size(bTxWxb,dim=1), & ! F2003 AUTO-REALLOC
                       size(bTxWxB,dim=2),zero)
  call gemm( real(bTxW,kind=wp) , real(bmat,kind=wp) , tmp_mat )
  bTxWxb(:,:) = real(tmp_mat,kind=kind(bTxWxb))
 !call gemm( bTxW , bmat , bTxWxb )
#else
  bTxWxb = matmul( bTxW , bmat )
#endif
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
#ifdef SPECIAL_FOR_INTEL
  bTxWxb = invert_matrix_lapack95( bTxWxb )
#else
  bTxWxb = invert_matrix( bTxWxb )
#endif
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
#ifdef SPECIAL_FOR_INTEL
    tmp = init_array(size(bTxWxb,dim=1), & ! F2003 AUTO-REALLOC
                     size(bTxW,dim=2),zero)
    tmp_mat = init_array(size(vmat,dim=1), & ! F2003 AUTO-REALLOC
                         size(bTxW,dim=2),zero)
    call gemm( real(bTxWxb,kind=wp) , real(bTxW,kind=wp) , tmp )
    call gemm( real(vmat,kind=wp) , tmp , tmp_mat )
    pmat = real(tmp_mat,kind=kind(pmat)) ! F2003 AUTO-REALLOC
   !call gemm( bTxWxb , bTxW , tmp )
   !call gemm( vmat , tmp , pmat )
#else
    pmat = matmul( vmat , matmul( bTxWxb , bTxW ) ) ! F2003 AUTO-REALLOCATION
#endif
    !
  else
    !
    ! pmat = basis x bTxWxb x bTxW
    !
#ifdef SPECIAL_FOR_INTEL
    tmp = init_array(size(bTxWxb,dim=1), & ! F2003 AUTO-REALLOC
                     size(bTxW,dim=2),zero)
    tmp_mat = init_array(size(bmat,dim=1), & ! F2003 AUTO-REALLOC
                         size(bTxW,dim=2),zero)
    call gemm( real(bTxWxb,kind=wp) , real(bTxW,kind=wp) , tmp )
    call gemm( real(bmat,kind=wp) , tmp , tmp_mat )
    pmat = real(tmp_mat,kind=kind(pmat)) ! F2003 AUTO-REALLOC
   !call gemm( bTxWxb , bTxW , tmp )
   !call gemm( bmat , tmp , pmat )
#else
    pmat = matmul( bmat , matmul( bTxWxb , bTxW ) ) ! F2003 AUTO-REALLOCATION
#endif
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
            array = matrix_memory_usage( projct(i,j)%CellOI(k) )
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
            array = matrix_memory_usage( projct(i,j)%FaceOI(k) )
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
            array = matrix_memory_usage( projct(i,j)%limiter(k) )
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
subroutine move_projct_to_device()
  !
  !.. Local Scalars ..
  integer :: i,i1,i2,j,j1,j2,k,k1,k2,m1,m2,n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "move_projct_to_device"
  !
continue
  !
  if (allocated(projct)) then
    !$ACC ENTER DATA COPYIN(projct)
    i1 = lbound(projct,dim=1)
    i2 = ubound(projct,dim=1)
    j1 = lbound(projct,dim=2)
    j2 = ubound(projct,dim=2)
    !$ACC ENTER DATA COPYIN(projct(i1:i2,j1:j2))
    do j = j1,j2
      do i = i1,i2
        !
        !---------------------------------------------------------------
        !
        if (allocated(projct(i,j)%limiter)) then
          !$ACC ENTER DATA COPYIN(projct(i,j)%limiter)
          k1 = lbound(projct(i,j)%limiter,dim=1)
          k2 = ubound(projct(i,j)%limiter,dim=1)
          !$ACC ENTER DATA COPYIN(projct(i,j)%limiter(k1:k2))
          do k = k1,k2
            !
            if (allocated(projct(i,j)%limiter(k)%mat)) then
              !$ACC ENTER DATA COPYIN(projct(i,j)%limiter(k)%mat)
              m1 = lbound(projct(i,j)%limiter(k)%mat,dim=1)
              m2 = ubound(projct(i,j)%limiter(k)%mat,dim=1)
              n1 = lbound(projct(i,j)%limiter(k)%mat,dim=2)
              n2 = ubound(projct(i,j)%limiter(k)%mat,dim=2)
              !$ACC ENTER DATA COPYIN(projct(i,j)%limiter(k)%mat(m1:m2,n1:n2))
            end if
            !
            if (allocated(projct(i,j)%limiter(k)%row_ptr)) then
              !$ACC ENTER DATA COPYIN(projct(i,j)%limiter(k)%row_ptr)
              m1 = lbound(projct(i,j)%limiter(k)%row_ptr,dim=1)
              m2 = ubound(projct(i,j)%limiter(k)%row_ptr,dim=1)
              !$ACC ENTER DATA COPYIN(projct(i,j)%limiter(k)%row_ptr(m1:m2))
            end if
            !
            if (allocated(projct(i,j)%limiter(k)%row_idx)) then
              !$ACC ENTER DATA COPYIN(projct(i,j)%limiter(k)%row_idx)
              m1 = lbound(projct(i,j)%limiter(k)%row_idx,dim=1)
              m2 = ubound(projct(i,j)%limiter(k)%row_idx,dim=1)
              !$ACC ENTER DATA COPYIN(projct(i,j)%limiter(k)%row_idx(m1:m2))
            end if
            !
            if (allocated(projct(i,j)%limiter(k)%row_val)) then
              !$ACC ENTER DATA COPYIN(projct(i,j)%limiter(k)%row_val)
              m1 = lbound(projct(i,j)%limiter(k)%row_val,dim=1)
              m2 = ubound(projct(i,j)%limiter(k)%row_val,dim=1)
              !$ACC ENTER DATA COPYIN(projct(i,j)%limiter(k)%row_val(m1:m2))
            end if
            !
            if (allocated(projct(i,j)%limiter(k)%col_ptr)) then
              !$ACC ENTER DATA COPYIN(projct(i,j)%limiter(k)%col_ptr)
              m1 = lbound(projct(i,j)%limiter(k)%col_ptr,dim=1)
              m2 = ubound(projct(i,j)%limiter(k)%col_ptr,dim=1)
              !$ACC ENTER DATA COPYIN(projct(i,j)%limiter(k)%col_ptr(m1:m2))
            end if
            !
            if (allocated(projct(i,j)%limiter(k)%col_idx)) then
              !$ACC ENTER DATA COPYIN(projct(i,j)%limiter(k)%col_idx)
              m1 = lbound(projct(i,j)%limiter(k)%col_idx,dim=1)
              m2 = ubound(projct(i,j)%limiter(k)%col_idx,dim=1)
              !$ACC ENTER DATA COPYIN(projct(i,j)%limiter(k)%col_idx(m1:m2))
            end if
            !
            if (allocated(projct(i,j)%limiter(k)%col_val)) then
              !$ACC ENTER DATA COPYIN(projct(i,j)%limiter(k)%col_val)
              m1 = lbound(projct(i,j)%limiter(k)%col_val,dim=1)
              m2 = ubound(projct(i,j)%limiter(k)%col_val,dim=1)
              !$ACC ENTER DATA COPYIN(projct(i,j)%limiter(k)%col_val(m1:m2))
            end if
            !
          end do
          !
        end if
        !
        !---------------------------------------------------------------
        !
        if (allocated(projct(i,j)%CellOI)) then
          !$ACC ENTER DATA COPYIN(projct(i,j)%CellOI)
          k1 = lbound(projct(i,j)%CellOI,dim=1)
          k2 = ubound(projct(i,j)%CellOI,dim=1)
          !$ACC ENTER DATA COPYIN(projct(i,j)%CellOI(k1:k2))
          do k = k1,k2
            !
            if (allocated(projct(i,j)%CellOI(k)%mat)) then
              !$ACC ENTER DATA COPYIN(projct(i,j)%CellOI(k)%mat)
              m1 = lbound(projct(i,j)%CellOI(k)%mat,dim=1)
              m2 = ubound(projct(i,j)%CellOI(k)%mat,dim=1)
              n1 = lbound(projct(i,j)%CellOI(k)%mat,dim=2)
              n2 = ubound(projct(i,j)%CellOI(k)%mat,dim=2)
              !$ACC ENTER DATA COPYIN(projct(i,j)%CellOI(k)%mat(m1:m2,n1:n2))
            end if
            !
            if (allocated(projct(i,j)%CellOI(k)%row_ptr)) then
              !$ACC ENTER DATA COPYIN(projct(i,j)%CellOI(k)%row_ptr)
              m1 = lbound(projct(i,j)%CellOI(k)%row_ptr,dim=1)
              m2 = ubound(projct(i,j)%CellOI(k)%row_ptr,dim=1)
              !$ACC ENTER DATA COPYIN(projct(i,j)%CellOI(k)%row_ptr(m1:m2))
            end if
            !
            if (allocated(projct(i,j)%CellOI(k)%row_idx)) then
              !$ACC ENTER DATA COPYIN(projct(i,j)%CellOI(k)%row_idx)
              m1 = lbound(projct(i,j)%CellOI(k)%row_idx,dim=1)
              m2 = ubound(projct(i,j)%CellOI(k)%row_idx,dim=1)
              !$ACC ENTER DATA COPYIN(projct(i,j)%CellOI(k)%row_idx(m1:m2))
            end if
            !
            if (allocated(projct(i,j)%CellOI(k)%row_val)) then
              !$ACC ENTER DATA COPYIN(projct(i,j)%CellOI(k)%row_val)
              m1 = lbound(projct(i,j)%CellOI(k)%row_val,dim=1)
              m2 = ubound(projct(i,j)%CellOI(k)%row_val,dim=1)
              !$ACC ENTER DATA COPYIN(projct(i,j)%CellOI(k)%row_val(m1:m2))
            end if
            !
            if (allocated(projct(i,j)%CellOI(k)%col_ptr)) then
              !$ACC ENTER DATA COPYIN(projct(i,j)%CellOI(k)%col_ptr)
              m1 = lbound(projct(i,j)%CellOI(k)%col_ptr,dim=1)
              m2 = ubound(projct(i,j)%CellOI(k)%col_ptr,dim=1)
              !$ACC ENTER DATA COPYIN(projct(i,j)%CellOI(k)%col_ptr(m1:m2))
            end if
            !
            if (allocated(projct(i,j)%CellOI(k)%col_idx)) then
              !$ACC ENTER DATA COPYIN(projct(i,j)%CellOI(k)%col_idx)
              m1 = lbound(projct(i,j)%CellOI(k)%col_idx,dim=1)
              m2 = ubound(projct(i,j)%CellOI(k)%col_idx,dim=1)
              !$ACC ENTER DATA COPYIN(projct(i,j)%CellOI(k)%col_idx(m1:m2))
            end if
            !
            if (allocated(projct(i,j)%CellOI(k)%col_val)) then
              !$ACC ENTER DATA COPYIN(projct(i,j)%CellOI(k)%col_val)
              m1 = lbound(projct(i,j)%CellOI(k)%col_val,dim=1)
              m2 = ubound(projct(i,j)%CellOI(k)%col_val,dim=1)
              !$ACC ENTER DATA COPYIN(projct(i,j)%CellOI(k)%col_val(m1:m2))
            end if
            !
          end do
          !
        end if
        !
        !---------------------------------------------------------------
        !
        if (allocated(projct(i,j)%FaceOI)) then
          !$ACC ENTER DATA COPYIN(projct(i,j)%FaceOI)
          k1 = lbound(projct(i,j)%FaceOI,dim=1)
          k2 = ubound(projct(i,j)%FaceOI,dim=1)
          !$ACC ENTER DATA COPYIN(projct(i,j)%FaceOI(k1:k2))
          do k = k1,k2
            !
            if (allocated(projct(i,j)%FaceOI(k)%mat)) then
              !$ACC ENTER DATA COPYIN(projct(i,j)%FaceOI(k)%mat)
              m1 = lbound(projct(i,j)%FaceOI(k)%mat,dim=1)
              m2 = ubound(projct(i,j)%FaceOI(k)%mat,dim=1)
              n1 = lbound(projct(i,j)%FaceOI(k)%mat,dim=2)
              n2 = ubound(projct(i,j)%FaceOI(k)%mat,dim=2)
              !$ACC ENTER DATA COPYIN(projct(i,j)%FaceOI(k)%mat(m1:m2,n1:n2))
            end if
            !
            if (allocated(projct(i,j)%FaceOI(k)%row_ptr)) then
              !$ACC ENTER DATA COPYIN(projct(i,j)%FaceOI(k)%row_ptr)
              m1 = lbound(projct(i,j)%FaceOI(k)%row_ptr,dim=1)
              m2 = ubound(projct(i,j)%FaceOI(k)%row_ptr,dim=1)
              !$ACC ENTER DATA COPYIN(projct(i,j)%FaceOI(k)%row_ptr(m1:m2))
            end if
            !
            if (allocated(projct(i,j)%FaceOI(k)%row_idx)) then
              !$ACC ENTER DATA COPYIN(projct(i,j)%FaceOI(k)%row_idx)
              m1 = lbound(projct(i,j)%FaceOI(k)%row_idx,dim=1)
              m2 = ubound(projct(i,j)%FaceOI(k)%row_idx,dim=1)
              !$ACC ENTER DATA COPYIN(projct(i,j)%FaceOI(k)%row_idx(m1:m2))
            end if
            !
            if (allocated(projct(i,j)%FaceOI(k)%row_val)) then
              !$ACC ENTER DATA COPYIN(projct(i,j)%FaceOI(k)%row_val)
              m1 = lbound(projct(i,j)%FaceOI(k)%row_val,dim=1)
              m2 = ubound(projct(i,j)%FaceOI(k)%row_val,dim=1)
              !$ACC ENTER DATA COPYIN(projct(i,j)%FaceOI(k)%row_val(m1:m2))
            end if
            !
            if (allocated(projct(i,j)%FaceOI(k)%col_ptr)) then
              !$ACC ENTER DATA COPYIN(projct(i,j)%FaceOI(k)%col_ptr)
              m1 = lbound(projct(i,j)%FaceOI(k)%col_ptr,dim=1)
              m2 = ubound(projct(i,j)%FaceOI(k)%col_ptr,dim=1)
              !$ACC ENTER DATA COPYIN(projct(i,j)%FaceOI(k)%col_ptr(m1:m2))
            end if
            !
            if (allocated(projct(i,j)%FaceOI(k)%col_idx)) then
              !$ACC ENTER DATA COPYIN(projct(i,j)%FaceOI(k)%col_idx)
              m1 = lbound(projct(i,j)%FaceOI(k)%col_idx,dim=1)
              m2 = ubound(projct(i,j)%FaceOI(k)%col_idx,dim=1)
              !$ACC ENTER DATA COPYIN(projct(i,j)%FaceOI(k)%col_idx(m1:m2))
            end if
            !
            if (allocated(projct(i,j)%FaceOI(k)%col_val)) then
              !$ACC ENTER DATA COPYIN(projct(i,j)%FaceOI(k)%col_val)
              m1 = lbound(projct(i,j)%FaceOI(k)%col_val,dim=1)
              m2 = ubound(projct(i,j)%FaceOI(k)%col_val,dim=1)
              !$ACC ENTER DATA COPYIN(projct(i,j)%FaceOI(k)%col_val(m1:m2))
            end if
            !
          end do
          !
        end if
        !
        !---------------------------------------------------------------
        !
      end do
    end do
  end if
  !
end subroutine move_projct_to_device
!
!###############################################################################
!
#include "Functions/debug_timer.f90"
!
!###############################################################################
!
end module projection_mod
!
!###############################################################################
!
