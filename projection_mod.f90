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
    type(matrix), allocatable :: toFace(:) ! project to a face point
    type(matrix), allocatable :: toCell(:) ! project to a cell point
  end type projection_matrix
  !
  ! projct : projection matrix for all combinations of cell geometry and 
  !          and cell order that are possible in the current simulation
  !
  !     For cell 'n' :
  !
  !         The projection matrix (giving the contribution of each solution
  !         point within cell 'n' to the solution point 'k' that is currently
  !         of interest) is given by the expression
  !
  !         projct(this_geom,from_order)%toFace(to_order)%mat(:,k)
  !         projct(this_geom,from_order)%toCell(to_order)%mat(:,k)
  !
  !     Where :
  !
  !        this_geom = cell(n)%geom
  !
  !   IF: over-integrating by projecting the residual (i.e. p_order < n_order)
  !        from_order = cell(n)%order
  !        to_order   = p_order
  !   ELSE IF: over-integrating by using quadrature (i.e. q_order > n_order)
  !        from_order = q_order
  !        to_order   = cell(n)%order
  !   END IF
  !        k = if using toFace, the face solution point we are projecting to
  !            if using toCell, the cell solution point we are projecting to
  !
  type(projection_matrix), public, save, target, allocatable :: projct(:,:)
  !
contains
!
!###############################################################################
!
subroutine init_projection_matrices()
  !
  !.. Use Statements ..
  use geovar,         only : cell
  use quadrature_mod, only : geom_is_used
  use quadrature_mod, only : std_elem
  use vandermonde_mod, only : vand
  !
  !.. Local Scalars ..
  integer :: ierr
  !
continue
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
  integer :: k,n,l,n1,n2,nsp,nvar
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
  do n = 1,size(cell)
    !
    n1 = cell(n)%beg_sp
    n2 = cell(n)%end_sp
    nsp = n2-n1+1
    !
    this_geom = cell(n)%geom
    this_order = cell(n)%order
    !
    if (proj_order < this_order) then
      !
      do k = 1,nsp
        do l = 1,nvar
          tmp(k,l) = array(l,n1-1+k)
        end do
      end do
      !
      pmat(1:nsp,1:nsp) = projct(this_geom,this_order)% &
                                    toCell(proj_order)% &
                          mat(1:nsp,1:nsp)
      !
      do k = 1,nsp
        do l = 1,nvar
          array(l,n1-1+k) = dot( pmat(1:nsp,k) , tmp(1:nsp,l) )
        end do
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
pure function generic_projection_matrix(basis,wmat) result(return_value)
  !
  real(qp), dimension(:,:), intent(in) :: basis
  real(qp), dimension(:,:), intent(in) :: wmat
  !
  real(qp), dimension(1:size(wmat,1),1:size(wmat,2)) :: return_value
  !
  real(qp), dimension(1:size(basis,2),1:size(basis,1)) :: basisT
  real(qp), dimension(1:size(basis,2),1:size(basis,2)) :: pmat
  real(qp), dimension(1:size(wmat,1),1:size(wmat,2)) :: qmat
  !
continue
  !
  basisT(:,:) = qzero
  pmat(:,:) = qzero
  qmat(:,:) = qzero
  !
  basisT = matmul( transpose(basis) , wmat )
  basisT = chop( basisT , basisT )
  !
  pmat = matmul( basisT , basis )
  pmat = chop( pmat , pmat )
  !
  pmat = invert_matrix( pmat )
  pmat = chop( pmat , pmat )
  !
  qmat = matmul( basis , matmul( pmat , basisT ) )
  qmat = chop( qmat , qmat )
  !
  return_value = transpose( qmat )
  return_value = chop( return_value , return_value )
  !
end function generic_projection_matrix
!
!###############################################################################
!
!!subroutine interp2quadrature(vand_1d_qp,vand_1d_inv_qp, &
!!                             vand_quad_qp,vand_quad_inv_qp)
!!  !
!!  !.. Use Statements ..
!!  use order_mod, only : n_order,p_order,o_order,n_sp1de,n_sp2dt,n_sp2dq
!!  use order_mod, only : q_order,n_qp1de,n_qp2dt,n_qp2dq
!!  !
!!  use geovar, only : using_quadrature
!!  !
!!  use ovar,   only : iop_qdtpts
!!  !
!!  use metrics_mod, only : dqr_LPmat
!!  use metrics_mod, only : dqs_LPmat
!!  !
!!  use qpvar, only : qdtpts_1d_qp
!!  use qpvar, only : gauss_1d_qp,lobatto_1d_qp
!!  use qpvar, only : gauss_wts_2d_qp,lobatto_wts_2d_qp
!!  use qpvar, only : dqr_LPmat_qp
!!  use qpvar, only : dqs_LPmat_qp
!!  !
!!  use polynomial_mod, only : nodes_legendre_gauss
!!  use polynomial_mod, only : weights_legendre_gauss
!!  use polynomial_mod, only : nodes_legendre_gauss_lobatto
!!  use polynomial_mod, only : weights_legendre_gauss_lobatto
!!  !
!!  use quadrature_mod, only : std_elem
!!  !
!!  !.. Formal Arguments ..
!!  real(qp), dimension(:,:), intent(in) :: vand_1d_qp
!!  real(qp), dimension(:,:), intent(in) :: vand_1d_inv_qp
!!  real(qp), dimension(:,:), intent(in) :: vand_quad_qp
!!  real(qp), dimension(:,:), intent(in) :: vand_quad_inv_qp
!!  !
!!  !.. Local Scalars ..
!!  integer :: i,j,k,l,n,ierr
!!  integer :: li,lj,ki,kj,nsp,nqp,npp,nop
!!  integer :: limit_order
!!  character(len=32) :: array_name
!!  !
!!  !.. Local Arrays ..
!!  real(qp) :: solpts_1d_qp(1:n_sp1de)
!!  real(qp) :: solwts_1d_qp(1:n_sp1de)
!!  real(qp), dimension(1:n_qp1de) :: vp_qdtr_qp
!!  real(qp), dimension(1:abs(o_order)+1) :: vp_out_qp
!!  real(qp), dimension(1:abs(o_order)+1) :: outpts_1d_qp
!!  real(qp), dimension(1:n_qp1de,1:n_qp1de) :: vand1D_qdtr_qp
!!  real(qp), dimension(1:n_qp1de,1:n_qp1de) :: vand1D_qdtr_inv_qp
!!  real(qp), dimension(1:n_qp1de,1:n_qp1de) :: grad_vand1D_qdtr_qp
!!  real(qp), dimension(1:n_qp1de,1:n_sp1de) :: vand1D_f2g_qp
!!  real(qp), dimension(1:n_sp2dq,1:n_sp2dq) :: wt_mat_qp
!!  real(qp), dimension(1:n_sp2dq,1:n_sp2dq) :: vand2D_lim_qp
!!  real(qp), dimension(1:n_qp2dq,1:n_sp2dq,1:2) :: vand2D_s2q_qp
!!  real(qp), dimension(1:n_sp2dq,1:(p_order+1)**2,1:2) :: vand2D_p2s_qp
!!  real(qp), dimension(1:(abs(o_order)+1)**2,1:n_sp2dq,1:2) :: vand2D_s2o_qp
!!  real(qp) :: vand2D_p2o_qp(1:(abs(o_order)+1)**2,1:(p_order+1)**2,1:2)
!!  !
!!  !.. Local Allocatable Arrays From Module qpvar ..
!!  real(qp), allocatable :: qdtwts_1d_qp(:)
!!  real(qp), allocatable :: invM_mat_qp(:,:,:)
!!  real(qp), allocatable :: dq1D_mat_qp(:,:)
!!  real(qp), allocatable :: qdtwts_2d_qp(:,:)
!!  real(qp), allocatable :: vand_qdtr_qp(:,:,:)
!!  real(qp), allocatable :: vand_qdtr_inv_qp(:,:,:)
!!  real(qp), allocatable :: interp_s2o_qp(:,:,:,:)
!!  real(qp), allocatable :: dqr_mat_qp(:,:,:)
!!  real(qp), allocatable :: dqs_mat_qp(:,:,:)
!!  real(qp), allocatable :: dqr_CRmat_qp(:,:,:)
!!  real(qp), allocatable :: dqs_CRmat_qp(:,:,:)
!!  real(qp), allocatable :: projct_limiter_qp(:,:,:,:)
!!  real(qp), allocatable :: projct_g2f_qp(:,:)
!!  real(qp), allocatable :: projct_q2s_qp(:,:,:)
!!  real(qp), allocatable :: projct_s2p_qp(:,:,:,:)
!!  real(qp), allocatable :: interp_f2g_qp(:,:)
!!  real(qp), allocatable :: interp_s2q_qp(:,:,:)
!!  real(qp), allocatable :: interp_s2g_qp(:,:,:,:)
!!  !
!!  !
!!  !.. Local Allocatable Arrays From Module geovar ..
!!  real(wp), allocatable :: qdtwts_1d(:)
!!  real(wp), allocatable :: invM_mat(:,:,:)
!!  real(wp), allocatable :: qdtwts_2d(:,:)
!!  real(wp), allocatable :: vand_qdtr(:,:,:)
!!  real(wp), allocatable :: vand_qdtr_inv(:,:,:)
!!  real(wp), allocatable :: interp_s2o(:,:,:,:)
!!  real(wp), allocatable :: projct_limiter(:,:,:,:)
!!  real(wp), allocatable :: projct_g2f(:,:)
!!  real(wp), allocatable :: projct_q2s(:,:,:)
!!  real(wp), allocatable :: projct_s2p(:,:,:,:)
!!  real(wp), allocatable :: interp_f2g(:,:)
!!  real(wp), allocatable :: interp_s2q(:,:,:)
!!  real(wp), allocatable :: interp_s2g(:,:,:,:)
!!  !
!!  !.. Local Allocatable Arrays From Module metrics_mod ..
!!  real(wp), allocatable :: dq1D_mat(:,:)
!!  real(wp), allocatable :: dqr_mat(:,:,:)
!!  real(wp), allocatable :: dqs_mat(:,:,:)
!!  real(wp), allocatable :: dqr_CRmat(:,:,:)
!!  real(wp), allocatable :: dqs_CRmat(:,:,:)
!!  !
!!  !.. Local Pointers ..
!!  real(wp), pointer, contiguous :: weights(:)
!!  !
!!  !.. Local Parameters ..
!!  logical, parameter :: array_debugging = fals
!! !logical, parameter :: array_debugging = true
!!  character(len=*), parameter :: pname = "interp2quadrature"
!!  !
!!continue
!!  !
!!  weights => null()
!!  !
!!  vp_qdtr_qp(:) = qzero
!!  vp_out_qp(:) = qzero
!!  outpts_1d_qp(:) = qzero
!!  vand1D_qdtr_qp(:,:) = qzero
!!  vand1D_qdtr_inv_qp(:,:) = qzero
!!  grad_vand1D_qdtr_qp(:,:) = qzero
!!  vand1D_f2g_qp(:,:) = qzero
!!  vand2D_s2q_qp(:,:,:) = qzero
!!  vand2D_p2s_qp(:,:,:) = qzero
!!  vand2D_s2o_qp(:,:,:) = qzero
!!  vand2D_p2o_qp(:,:,:) = qzero
!!  wt_mat_qp(:,:) = qzero
!!  !
!!  npp = (p_order+1)**2
!!  nop = (abs(o_order)+1)**2
!!  !
!!  solpts_1d_qp = real(solpts_1d,kind=qp)
!!  solwts_1d_qp = real(solwts_1d,kind=qp)
!!  !
!!  weights => std_elem(Geom_Quad,n_order)%wts
!!  !
!!  if (using_quadrature .and. mypnum == 0) write (iout,1)
!!  !
!!  ! ###########################################################
!!  ! ###########################################################
!!  ! ##                                                       ##
!!  ! ##  Allocate all the arrays needed for quadrature in 2D  ##
!!  ! ##                                                       ##
!!  ! ###########################################################
!!  ! ###########################################################
!!  !
!!  ! #####################################################################
!!  ! #####################################################################
!!  ! ##     USE QUAD PRECISION ARRAYS TO CREATE ALL OF THE ARRAYS       ##
!!  ! ##   PREVIOUSLY ALLOCATED. WE WILL THEN PERFORM A CHOP TO ZERO     ##
!!  ! ##  EVERYTHING WITH AN ABSOLUTE VALUE SMALLER THAN EPSILON(0_WP).  ##
!!  ! #####################################################################
!!  ! #####################################################################
!!  !
!!  ! Allocate the quad precision arrays
!!  !
!!  allocate ( qdtpts_1d_qp(1:n_qp1de) , source=qzero , stat=ierr )
!!  call alloc_error(pname,"qdtpts_1d_qp",1,ierr)
!!  allocate ( qdtwts_1d_qp(1:n_qp1de) , source=qzero , stat=ierr )
!!  call alloc_error(pname,"qdtwts_1d_qp",1,ierr)
!!  allocate ( dq1D_mat_qp(1:n_qp1de,1:n_qp1de), source=qzero , stat=ierr )
!!  call alloc_error(pname,"dq1D_mat_qp",1,ierr)
!!  allocate ( interp_f2g_qp(1:n_sp1de,1:n_qp1de) , source=qzero , stat=ierr )
!!  call alloc_error(pname,"interp_f2g_qp",1,ierr)
!!  allocate ( projct_g2f_qp(1:n_qp1de,1:n_sp1de) , source=qzero , stat=ierr )
!!  call alloc_error(pname,"projct_g2f_qp",1,ierr)
!!  !
!!  allocate ( qdtwts_2d_qp(1:n_qp2dq,1:2) , source=qzero , stat=ierr )
!!  call alloc_error(pname,"qdtwts_2d_qp",1,ierr)
!!  allocate ( vand_qdtr_qp(1:n_qp2dq,1:n_qp2dq,1:2), source=qzero , stat=ierr )
!!  call alloc_error(pname,"vand_qdtr_qp",1,ierr)
!!  allocate ( vand_qdtr_inv_qp(1:n_qp2dq,1:n_qp2dq,1:2), source=qzero , &
!!             stat=ierr )
!!  call alloc_error(pname,"vand_qdtr_inv_qp",1,ierr)
!!  allocate ( interp_s2q_qp(1:n_sp2dq,1:n_qp2dq,1:2) , source=qzero , &
!!             stat=ierr )
!!  call alloc_error(pname,"interp_s2q_qp",1,ierr)
!!  allocate ( interp_s2o_qp(1:n_sp2dq,1:nop,1:2,1:3) , source=qzero , &
!!             stat=ierr )
!!  call alloc_error(pname,"interp_s2o_qp",1,ierr)
!!  allocate ( projct_q2s_qp(1:n_qp2dq,1:n_sp2dq,1:2) , source=qzero , &
!!             stat=ierr )
!!  call alloc_error(pname,"projct_q2s_qp",1,ierr)
!!  allocate ( projct_s2p_qp(1:n_sp2dq,1:n_sp2dq,1:2,1:5) , source=qzero , &
!!             stat=ierr )
!!  call alloc_error(pname,"projct_s2p_qp",1,ierr)
!!  allocate ( projct_limiter_qp(1:n_sp2dq,1:n_sp2dq,1:2,1:n_order) , &
!!             source=qzero , stat=ierr )
!!  call alloc_error(pname,"projct_limiter_qp",1,ierr)
!!  allocate ( dqr_mat_qp(1:n_qp2dq,1:n_qp2dq,1:2), source=qzero , stat=ierr )
!!  call alloc_error(pname,"dqr_mat_qp",1,ierr)
!!  allocate ( dqs_mat_qp(1:n_qp2dq,1:n_qp2dq,1:2), source=qzero , stat=ierr )
!!  call alloc_error(pname,"dqs_mat_qp",1,ierr)
!!  allocate ( dqr_LPmat_qp(1:n_qp2dq,1:n_sp2dq,1:2), source=qzero , stat=ierr )
!!  call alloc_error(pname,"dqr_LPmat_qp",1,ierr)
!!  allocate ( dqs_LPmat_qp(1:n_qp2dq,1:n_sp2dq,1:2), source=qzero , stat=ierr )
!!  call alloc_error(pname,"dqs_LPmat_qp",1,ierr)
!!  allocate ( dqr_CRmat_qp(1:n_sp2dq,1:n_qp2dq,1:2), source=qzero , stat=ierr )
!!  call alloc_error(pname,"dqr_CRmat_qp",1,ierr)
!!  allocate ( dqs_CRmat_qp(1:n_sp2dq,1:n_qp2dq,1:2), source=qzero , stat=ierr )
!!  call alloc_error(pname,"dqs_CRmat_qp",1,ierr)
!!  allocate ( invM_mat_qp(1:n_sp2dq,1:n_sp2dq,1:2) , source=qzero , stat=ierr )
!!  call alloc_error(pname,"invM_mat_qp",1,ierr)
!!  allocate ( interp_s2g_qp(1:n_sp2dq,1:4*n_qp1de,1:2) , source=qzero , &
!!             stat=ierr )
!!  call alloc_error(pname,"interp_s2g_qp",1,ierr)
!!  !
!!  ! ####################################
!!  ! ####################################
!!  ! ##                                ##
!!  ! ##  Section for quadrature in 1D  ##
!!  ! ##                                ##
!!  ! ####################################
!!  ! ####################################
!!  !
!!  ! 1D quadrature nodes and weights
!!  !
!!  if (iop_qdtpts == Legendre_Gauss_Lobatto) then
!!    call nodes_legendre_gauss_lobatto(q_order,qdtpts_1d_qp,vp_qdtr_qp)
!!    call weights_legendre_gauss_lobatto(q_order,vp_qdtr_qp,qdtwts_1d_qp)
!!  else
!!    call nodes_legendre_gauss(n_qp1de,qdtpts_1d_qp,vp_qdtr_qp)
!!    call weights_legendre_gauss(n_qp1de,qdtpts_1d_qp,vp_qdtr_qp,qdtwts_1d_qp)
!!  end if
!!  qdtpts_1d_qp = chop( qdtpts_1d_qp , qdtpts_1d_qp )
!!  qdtwts_1d_qp = chop( qdtwts_1d_qp , qdtwts_1d_qp )
!!  !
!!  ! 1D output nodes
!!  !
!!  call nodes_legendre_gauss_lobatto(abs(o_order),outpts_1d_qp,vp_out_qp)
!!  outpts_1d_qp = chop( outpts_1d_qp , outpts_1d_qp )
!!  !
!!  ! 1D Derivative matrix for the quadrature points
!!  !
!!  vand1D_qdtr_qp(:,:) = Vandermonde1D(qdtpts_1d_qp)
!!  vand1D_qdtr_qp = chop( vand1D_qdtr_qp , vand1D_qdtr_qp )
!!  !
!!  vand1D_qdtr_inv_qp(:,:) = invert_matrix( vand1D_qdtr_qp )
!!  vand1D_qdtr_inv_qp = chop( vand1D_qdtr_inv_qp , vand1D_qdtr_inv_qp )
!!  !
!!  grad_vand1D_qdtr_qp(:,:) = GradVandermonde1D(qdtpts_1d_qp)
!!  grad_vand1D_qdtr_qp = chop( grad_vand1D_qdtr_qp , grad_vand1D_qdtr_qp )
!!  !
!!  dq1D_mat_qp(:,:) = matmul( grad_vand1D_qdtr_qp , vand1D_qdtr_inv_qp )
!!  dq1D_mat_qp = chop( dq1D_mat_qp , dq1D_mat_qp )
!!  !
!!  ! 1D Vandermonde matrix giving the value for each of the
!!  ! #(n_sp1de) modes at the #(n_qp1de) quadrature points
!!  !
!!  vand1d_f2g_qp(:,:) = Vandermonde1D_N(qdtpts_1d_qp,n_sp1de)
!!  vand1d_f2g_qp = chop( vand1d_f2g_qp , vand1d_f2g_qp )
!!  !
!!  ! Matrix for interpolating 1D nodal values at face
!!  ! flux/solution points to face quadrature points
!!  !
!!  interp_f2g_qp(:,:) = transpose( matmul( vand1d_f2g_qp , vand_1d_inv_qp ) )
!!  interp_f2g_qp = chop( interp_f2g_qp , interp_f2g_qp )
!!  !
!!  ! Matrix for projecting nodal values at face quadrature
!!  ! points to the face flux/solution points
!!  !
!!  projct_g2f_qp(:,:) = matmul( vand1D_f2g_qp , transpose(vand_1D_qp) )
!!  projct_g2f_qp = chop( projct_g2f_qp , projct_g2f_qp )
!!  do i = 1,n_sp1de
!!    projct_g2f_qp(:,i) = projct_g2f_qp(:,i) * qdtwts_1d_qp(:)
!!  end do
!!  projct_g2f_qp = chop( projct_g2f_qp , projct_g2f_qp )
!!  !
!!  !
!!  ! ################################################
!!  ! ################################################
!!  ! ##                                            ##
!!  ! ##  Section for quadrature on triangle cells  ##
!!  ! ##                                            ##
!!  ! ################################################
!!  ! ################################################
!!  !
!!  nsp = n_sp2dt
!!  nqp = n_qp2dt
!!  !
!!  ! #####################################################
!!  ! #####################################################
!!  ! ##                                                 ##
!!  ! ##  Section for quadrature on quadrilateral cells  ##
!!  ! ##                                                 ##
!!  ! #####################################################
!!  ! #####################################################
!!  !
!!  nsp = n_sp2dq
!!  nqp = n_qp2dq
!!  !
!!  ! 2D quadrature weights for the quadrature points
!!  !
!!  k = 0
!!  do i = 1,n_qp1de
!!    do j = 1,n_qp1de
!!      k = k + 1
!!      qdtwts_2d_qp(k,2) = qdtwts_1d_qp(i)*qdtwts_1d_qp(j)
!!    end do
!!  end do
!!  qdtwts_2d_qp = chop( qdtwts_2d_qp , qdtwts_2d_qp )
!!  !
!!  do i = 1,n_qp2dq
!!    wt_mat_qp(i,i) = qdtwts_2d_qp(i,2)
!!  end do
!!  wt_mat_qp = chop( wt_mat_qp , wt_mat_qp )
!!  !
!!  ! Full Vandermonde matrix for the quadrature nodes and modes
!!  !
!!  vand_qdtr_qp(1:nqp,1:nqp,2) = Vandermonde2D_quad(qdtpts_1d_qp)
!!  vand_qdtr_qp = chop( vand_qdtr_qp , vand_qdtr_qp )
!!  !
!!  vand_qdtr_inv_qp(1:nqp,1:nqp,2) = invert_matrix(vand_qdtr_qp(1:nqp,1:nqp,2))
!!  vand_qdtr_inv_qp = chop( vand_qdtr_inv_qp , vand_qdtr_inv_qp )
!!  !
!!  ! ###########################################
!!  !  2D matrix for project solution to p_order
!!  ! ###########################################
!!  !
!!  ! Project at solution points
!!  !
!!  vand2D_p2s_qp(1:nsp,1:npp,2) = Vandermonde2D_quad_N(solpts_1d_qp,p_order)
!!  vand2D_p2s_qp = chop( vand2D_p2s_qp , vand2D_p2s_qp )
!!  !
!!  projct_s2p_qp(1:nsp,1:nsp,2,1) = matmul( vand2D_p2s_qp(1:nsp,1:npp,2) , &
!!                                 transpose(vand2D_p2s_qp(1:nsp,1:npp,2)) )
!!  projct_s2p_qp = chop( projct_s2p_qp , projct_s2p_qp )
!!  do i = 1,nsp
!!    projct_s2p_qp(i,1:nsp,2,1) = projct_s2p_qp(i,1:nsp,2,1) * &
!!                                 real(weights(i),kind=qp)
!!  end do
!!  projct_s2p_qp = chop( projct_s2p_qp , projct_s2p_qp )
!!  !
!!  ! Project at Gauss points
!!  !
!!  vand2D_p2s_qp(1:nsp,1:npp,2) = Vandermonde2D_quad_N(gauss_1d_qp,p_order)
!!  vand2D_p2s_qp = chop( vand2D_p2s_qp , vand2D_p2s_qp )
!!  !
!!  projct_s2p_qp(1:nsp,1:nsp,2,2) = &
!!    matmul( vand2D_p2s_qp(1:nsp,1:npp,2) , &
!!    transpose(vand2D_p2s_qp(1:nsp,1:npp,2)) )
!!  projct_s2p_qp = chop( projct_s2p_qp , projct_s2p_qp )
!!  do i = 1,nsp
!!    projct_s2p_qp(i,1:nsp,2,2) = projct_s2p_qp(i,1:nsp,2,2) * &
!!                                 gauss_wts_2d_qp(i)
!!  end do
!!  projct_s2p_qp = chop( projct_s2p_qp , projct_s2p_qp )
!!  !
!!  ! Project at Lobatto points
!!  !
!!  vand2D_p2s_qp(1:nsp,1:npp,2) = Vandermonde2D_quad_N(lobatto_1d_qp,p_order)
!!  vand2D_p2s_qp = chop( vand2D_p2s_qp , vand2D_p2s_qp )
!!  !
!!  projct_s2p_qp(1:nsp,1:nsp,2,3) = &
!!    matmul( vand2D_p2s_qp(1:nsp,1:npp,2) , &
!!            transpose(vand2D_p2s_qp(1:nsp,1:npp,2)) )
!!  projct_s2p_qp = chop( projct_s2p_qp , projct_s2p_qp )
!!  do i = 1,nsp
!!    projct_s2p_qp(i,1:nsp,2,3) = projct_s2p_qp(i,1:nsp,2,3) * &
!!                                 lobatto_wts_2d_qp(i)
!!  end do
!!  projct_s2p_qp = chop( projct_s2p_qp , projct_s2p_qp )
!!  !
!!  ! Index 4 used to be radau projection, but Im not sure what index 5 was.
!!  ! I just copied over the solution and Gauss projections, respectively.
!!  !
!!  projct_s2p_qp(1:nsp,1:nsp,2,4) = projct_s2p_qp(1:nsp,1:nsp,2,1)
!!  projct_s2p_qp(1:nsp,1:nsp,2,5) = projct_s2p_qp(1:nsp,1:nsp,2,2)
!!  !
!!  ! ###########################################
!!  ! Create the projection limiter matrices
!!  ! ###########################################
!!  !
!!  do n = 1,n_order
!!    !
!!    limit_order = n_order - n
!!    !
!!    ! ------------------
!!    ! For quadrilaterals
!!    ! ------------------
!!    !
!!    nsp = n_sp2dq
!!    npp = (limit_order+1)**2
!!    !
!!    ! Get the 2D Vandermonde matrix for the current projection limiter order
!!    !
!!    vand2D_lim_qp(1:nsp,1:npp) = Vandermonde2D_quad_N(solpts_1d_qp, &
!!                                                      limit_order)
!!    vand2D_lim_qp = chop( vand2D_lim_qp , vand2D_lim_qp )
!!    !
!!    ! Compute the matrix multiplication of Vand*transpose(Vand)
!!    !
!!    projct_limiter_qp(1:nsp,1:nsp,2,n) = &
!!      matmul( vand2D_lim_qp(1:nsp,1:npp) , &
!!              transpose(vand2D_lim_qp(1:nsp,1:npp)) )
!!    projct_limiter_qp = chop( projct_limiter_qp , projct_limiter_qp )
!!    !
!!    ! Finally add in the quadrature weights for the solution points to
!!    ! complete the projection matrices for the current limiter order
!!    !
!!    do i = 1,nsp
!!      projct_limiter_qp(i,1:nsp,2,n) = projct_limiter_qp(i,1:nsp,2,n) * &
!!                                       real(weights(i),kind=qp)
!!    end do
!!    projct_limiter_qp = chop( projct_limiter_qp , projct_limiter_qp )
!!    !
!!  end do
!!  !
!!  ! ###########################################
!!  ! 2D matrix for interpolating nodal values at SPs to QPs
!!  ! ###########################################
!!  !
!!  vand2D_s2q_qp(1:nqp,1:nsp,2) = Vandermonde2D_quad_N(qdtpts_1d_qp,n_order)
!!  vand2D_s2q_qp = chop( vand2D_s2q_qp , vand2D_s2q_qp )
!!  !
!!  interp_s2q_qp(1:nsp,1:nqp,2) = &
!!        transpose( matmul( vand2D_s2q_qp(1:nqp,1:nsp,2) , vand_quad_inv_qp ) )
!!  interp_s2q_qp = chop( interp_s2q_qp , interp_s2q_qp )
!!  !
!!  ! ###########################################
!!  ! Matrix for interpolating cell nodal values at SPs to face QPs
!!  ! ###########################################
!!  !
!!  interp_s2g_qp(1:nsp,1:4*n_qp1de,2) = InterpToFace2D_quad(qdtpts_1d_qp, &
!!                                                           solpts_1d_qp)
!!  interp_s2g_qp = chop( interp_s2g_qp , interp_s2g_qp )
!!  !
!!  ! ###########################################
!!  ! Matrix for interpolating cell nodal values at SPs to output SPs
!!  ! ###########################################
!!  !
!!  vand2D_p2o_qp(1:nop,1:npp,2) = Vandermonde2D_quad_N(outpts_1d_qp,p_order)
!!  vand2D_p2o_qp = chop( vand2D_p2o_qp , vand2D_p2o_qp )
!! !call output2d(true,real(vand2D_p2o_qp(1:nop,1:npp,2),kind=wp), &
!!                2,2,2,'vand2d_p2o')
!! !call output2d(true,real(vand2D_p2s_qp(1:nsp,1:npp,2),kind=wp), &
!!                2,2,2,'vand2d_p2s')
!!  !
!!  interp_s2o_qp(1:nsp,1:nop,2,1) = &
!!                       interpolation_matrix(vand2d_p2s_qp(1:nsp,1:npp,2), &
!!                                            vand2d_p2o_qp(1:nop,1:npp,2), &
!!                                            wt_mat_qp(1:nsp,1:nsp))
!!  !
!!  vand2D_s2o_qp(1:nop,1:nsp,2) = Vandermonde2D_quad_N(outpts_1d_qp,n_order)
!!  vand2D_s2o_qp = chop( vand2D_s2o_qp , vand2D_s2o_qp )
!!  !
!!  interp_s2o_qp(1:nsp,1:nop,2,2) = &
!!        transpose( matmul( vand2D_s2o_qp(1:nop,1:nsp,2) , vand_quad_inv_qp ) )
!!  interp_s2o_qp = chop( interp_s2o_qp , interp_s2o_qp )
!!  !
!!  interp_s2o_qp(1:nsp,1:nop,2,3) = &
!!   matmul( vand2D_p2s_qp(1:nsp,1:npp,2) , &
!!           transpose(vand2d_p2o_qp(1:nop,1:npp,2)) )
!!  !transpose(matmul( vand2D_p2o_qp(1:nop,1:npp,2) , &
!!                     transpose(vand2d_p2s_qp(1:nsp,1:npp,2)) ))
!!  interp_s2o_qp = chop( interp_s2o_qp , interp_s2o_qp )
!!  do i = 1,nop
!!    interp_s2o_qp(1:nsp,i,2,3) = interp_s2o_qp(1:nsp,i,2,3) * &
!!                                 real(weights(1:nsp),kind=qp)
!!  end do
!!  interp_s2o_qp = chop( interp_s2o_qp , interp_s2o_qp )
!!  !
!!  ! ###########################################
!!  ! 2D matrix to project nodal values at QPs to SPs
!!  ! ###########################################
!!  !
!!  projct_q2s_qp(1:nqp,1:nsp,2) = &
!!        matmul( vand2D_s2q_qp(1:nqp,1:nsp,2) , transpose(vand_quad_qp) )
!!  projct_q2s_qp = chop( projct_q2s_qp , projct_q2s_qp )
!!  do i = 1,nsp
!!    projct_q2s_qp(1:nqp,i,2) = projct_q2s_qp(1:nqp,i,2) * &
!!                               qdtwts_2d_qp(1:nqp,2)
!!  end do
!!  projct_q2s_qp = chop( projct_q2s_qp , projct_q2s_qp )
!!  !
!!  ! ###########################################
!!  ! 2D Derivative matrix for the quadrature points
!!  ! ###########################################
!!  !
!!  call Dmatrix2D_quad(dq1d_mat_qp,dqr_mat_qp(1:nqp,1:nqp,2), &
!!                                  dqs_mat_qp(1:nqp,1:nqp,2))
!!  dqr_mat_qp = chop( dqr_mat_qp , dqr_mat_qp )
!!  dqs_mat_qp = chop( dqs_mat_qp , dqs_mat_qp )
!!  !
!!  ! Pre-multiply the derivative matrices by the projection
!!  ! matrix for when using the derivative of the Lagrange polynomial
!!  !
!!  dqr_LPmat_qp(1:nqp,1:nsp,2) = &
!!      matmul( transpose(dqr_mat_qp(1:nqp,1:nqp,2)) , &
!!              projct_q2s_qp(1:nqp,1:nsp,2) )
!!  dqs_LPmat_qp(1:nqp,1:nsp,2) = &
!!      matmul( transpose(dqs_mat_qp(1:nqp,1:nqp,2)) , &
!!              projct_q2s_qp(1:nqp,1:nsp,2) )
!!  dqr_LPmat_qp = chop( dqr_LPmat_qp , dqr_LPmat_qp )
!!  dqs_LPmat_qp = chop( dqs_LPmat_qp , dqs_LPmat_qp )
!!  !
!!  ! Pre-multiply the derivative matrices by the
!!  ! interpolation matrix for when using the chain-rule
!!  !
!!  dqr_CRmat_qp(1:nsp,1:nqp,2) = &
!!      matmul( interp_s2q_qp(1:nsp,1:nqp,2) , &
!!              transpose(dqr_mat_qp(1:nqp,1:nqp,2)) )
!!  dqs_CRmat_qp(1:nsp,1:nqp,2) = &
!!      matmul( interp_s2q_qp(1:nsp,1:nqp,2) , &
!!              transpose(dqs_mat_qp(1:nqp,1:nqp,2)) )
!!  dqr_CRmat_qp = chop( dqr_CRmat_qp , dqr_CRmat_qp )
!!  dqs_CRmat_qp = chop( dqs_CRmat_qp , dqs_CRmat_qp )
!!  !
!!  ! Inverse of the mass matrix
!!  !
!!  invM_mat_qp(1:nsp,1:nsp,2) = matmul( vand_quad_qp , transpose(vand_quad_qp))
!!  invM_mat_qp = chop( invM_mat_qp , invM_mat_qp )
!!  !
!!  ! ##############################################
!!  ! ##############################################
!!  ! ##   CHOP THE QUAD PRECISION ARRAYS AND     ##
!!  ! ##  STORE IN THE WORKING PRECISION ARRAYS.  ##
!!  ! ##############################################
!!  ! ##############################################
!!  !
!!  ! Allocate the working precision arrays
!!  !
!!  allocate ( qdtpts_1d(1:n_qp1de) , source=zero , stat=ierr )
!!  call alloc_error(pname,"qdtpts_1d",1,ierr)
!!  allocate ( qdtwts_1d(1:n_qp1de) , source=zero , stat=ierr )
!!  call alloc_error(pname,"qdtwts_1d",1,ierr)
!!  allocate ( dq1D_mat(1:n_qp1de,1:n_qp1de), source=zero , stat=ierr )
!!  call alloc_error(pname,"dq1D_mat",1,ierr)
!!  allocate ( interp_f2g(1:n_sp1de,1:n_qp1de) , source=zero , stat=ierr )
!!  call alloc_error(pname,"interp_f2g",1,ierr)
!!  allocate ( projct_g2f(1:n_qp1de,1:n_sp1de) , source=zero , stat=ierr )
!!  call alloc_error(pname,"projct_g2f",1,ierr)
!!  !
!!  allocate ( qdtwts_2d(1:n_qp2dq,1:2) , source=zero , stat=ierr )
!!  call alloc_error(pname,"qdtwts_2d",1,ierr)
!!  allocate ( vand_qdtr(1:n_qp2dq,1:n_qp2dq,1:2), source=zero , stat=ierr )
!!  call alloc_error(pname,"vand_qdtr",1,ierr)
!!  allocate ( vand_qdtr_inv(1:n_qp2dq,1:n_qp2dq,1:2), source=zero , stat=ierr )
!!  call alloc_error(pname,"vand_qdtr_inv",1,ierr)
!!  allocate ( interp_s2q(1:n_sp2dq,1:n_qp2dq,1:2) , source=zero , stat=ierr )
!!  call alloc_error(pname,"interp_s2q",1,ierr)
!!  allocate ( interp_s2o(1:n_sp2dq,1:nop,1:2,1:3) , source=zero , stat=ierr )
!!  call alloc_error(pname,"interp_s2o",1,ierr)
!!  allocate ( projct_q2s(1:n_qp2dq,1:n_sp2dq,1:2) , source=zero , stat=ierr )
!!  call alloc_error(pname,"projct_q2s",1,ierr)
!!  allocate ( projct_s2p(1:n_sp2dq,1:n_sp2dq,1:2,1:5) , source=zero,stat=ierr )
!!  call alloc_error(pname,"projct_s2p",1,ierr)
!!  allocate ( dqr_mat(1:n_qp2dq,1:n_qp2dq,1:2), source=zero , stat=ierr )
!!  call alloc_error(pname,"dqr_mat",1,ierr)
!!  allocate ( dqs_mat(1:n_qp2dq,1:n_qp2dq,1:2), source=zero , stat=ierr )
!!  call alloc_error(pname,"dqs_mat",1,ierr)
!!  allocate ( dqr_LPmat(1:n_qp2dq,1:n_sp2dq,1:2), source=zero , stat=ierr )
!!  call alloc_error(pname,"dqr_LPmat",1,ierr)
!!  allocate ( dqs_LPmat(1:n_qp2dq,1:n_sp2dq,1:2), source=zero , stat=ierr )
!!  call alloc_error(pname,"dqs_LPmat",1,ierr)
!!  allocate ( dqr_CRmat(1:n_sp2dq,1:n_qp2dq,1:2), source=zero , stat=ierr )
!!  call alloc_error(pname,"dqr_CRmat",1,ierr)
!!  allocate ( dqs_CRmat(1:n_sp2dq,1:n_qp2dq,1:2), source=zero , stat=ierr )
!!  call alloc_error(pname,"dqs_CRmat",1,ierr)
!!  allocate ( invM_mat(1:n_sp2dq,1:n_sp2dq,1:2) , source=zero , stat=ierr )
!!  call alloc_error(pname,"invM_mat",1,ierr)
!!  allocate ( interp_s2g(1:n_sp2dq,1:4*n_qp1de,1:2) , source=zero , stat=ierr )
!!  call alloc_error(pname,"interp_s2g",1,ierr)
!!  allocate ( projct_limiter(1:n_sp2dq,1:n_sp2dq,1:2,1:n_order) , &
!!             source=zero , stat=ierr )
!!  call alloc_error(pname,"projct_limiter",1,ierr)
!!  !
!!  qdtpts_1d      = chop( qdtpts_1d_qp )
!!  qdtwts_1d      = chop( qdtwts_1d_qp )
!!  dq1D_mat       = chop( dq1D_mat_qp )
!!  interp_f2g     = chop( interp_f2g_qp )
!!  projct_g2f     = chop( projct_g2f_qp )
!!  qdtwts_2d      = chop( qdtwts_2d_qp )
!!  vand_qdtr      = chop( vand_qdtr_qp )
!!  vand_qdtr_inv  = chop( vand_qdtr_inv_qp )
!!  interp_s2q     = chop( interp_s2q_qp )
!!  interp_s2o     = chop( interp_s2o_qp )
!!  interp_s2g     = chop( interp_s2g_qp )
!!  projct_q2s     = chop( projct_q2s_qp )
!!  projct_s2p     = chop( projct_s2p_qp )
!!  projct_limiter = chop( projct_limiter_qp )
!!  dqr_mat        = chop( dqr_mat_qp )
!!  dqs_mat        = chop( dqs_mat_qp )
!!  dqr_LPmat      = chop( dqr_LPmat_qp )
!!  dqs_LPmat      = chop( dqs_LPmat_qp )
!!  dqr_CRmat      = chop( dqr_CRmat_qp )
!!  dqs_CRmat      = chop( dqs_CRmat_qp )
!!  invM_mat       = chop( invM_mat_qp )
!!  !
!!  call output2d(array_debugging,interp_f2g(:,:), &
!!                2,2,3,'interp_f2g')
!!  call output2d(array_debugging,projct_g2f(:,:), &
!!                2,2,3,'projct_g2f')
!!  call output2d(array_debugging,interp_s2q(1:nsp,1:nqp,2), &
!!                2,2,2,'interp_s2q')
!!  call output2d(array_debugging,interp_s2o(1:nsp,1:nqp,2,1), &
!!                2,2,2,'interp_s2o')
!! !call output2d(true,interp_s2o(1:nsp,1:nop,2,1), &
!!                2,2,2,'interp_s2o1')
!! !call output2d(true,interp_s2o(1:nsp,1:nop,2,2), &
!!                2,2,2,'interp_s2o2')
!! !call output2d(true,interp_s2o(1:nsp,1:nop,2,3), &
!!                2,2,2,'interp_s2o3')
!!  call output2d(array_debugging,interp_s2g(1:nsp,1:4*n_qp1de,2), &
!!                2,1,2,'interp_s2g')
!!  call output2d(array_debugging,projct_q2s(1:nqp,1:nsp,2), &
!!                2,2,2,'projct_q2s')
!!  call output2d(array_debugging,projct_s2p(1:nsp,1:nsp,2,1), &
!!                2,2,2,'projct_s2p')
!! !call output2d(true,projct_s2p(1:nsp,1:nsp,2,1), &
!!                2,2,2,'projct_s2p1')
!! !call output2d(true,projct_s2p(1:nsp,1:nsp,2,5), &
!!                2,2,2,'projct_s2p5')
!!  do n = n_order,1,-1
!!    write (array_name,'("projct_limiter",i0)') n
!!   !call output2d(true,           projct_limiter(1:nsp,1:nsp,2,n),wp), &
!!                  2,2,2,array_name)
!!    call output2d(array_debugging,projct_limiter(1:nsp,1:nsp,2,n),wp), &
!!                  2,2,2,array_name)
!!  end do
!!  call output2d(array_debugging,dqr_mat(1:nqp,1:nqp,2),2,2,2,'dqr_mat')
!!  call output2d(array_debugging,dqs_mat(1:nqp,1:nqp,2),2,2,2,'dqs_mat')
!!  call output2d(array_debugging,dqr_LPmat(1:nqp,1:nsp,2),2,2,2,'dqr_LPmat')
!!  call output2d(array_debugging,dqs_LPmat(1:nqp,1:nsp,2),2,2,2,'dqs_LPmat')
!!  call output2d(array_debugging,dqr_CRmat(1:nsp,1:nqp,2),2,2,2,'dqr_CRmat')
!!  call output2d(array_debugging,dqs_CRmat(1:nsp,1:nqp,2),2,2,2,'dqs_CRmat')
!!  call output2d(array_debugging,invM_mat(1:nsp,1:nsp,2),2,2,2,'invM_mat')
!!  !
!!  ! Deallocate the local allocatable arrays before we leave
!!  !
!!  deallocate ( dq1D_mat_qp , stat=ierr )
!!  call alloc_error(pname,"dq1D_mat_qp",2,ierr)
!!  deallocate ( qdtwts_2d_qp , stat=ierr )
!!  call alloc_error(pname,"qdtwts_2d_qp",2,ierr)
!!  deallocate ( vand_qdtr_qp , stat=ierr )
!!  call alloc_error(pname,"vand_qdtr_qp",2,ierr)
!!  deallocate ( vand_qdtr_inv_qp , stat=ierr )
!!  call alloc_error(pname,"vand_qdtr_inv_qp",2,ierr)
!!  deallocate ( interp_s2o_qp , stat=ierr )
!!  call alloc_error(pname,"interp_s2o_qp",2,ierr)
!!  deallocate ( dqr_mat_qp , stat=ierr )
!!  call alloc_error(pname,"dqr_mat_qp",2,ierr)
!!  deallocate ( dqs_mat_qp , stat=ierr )
!!  call alloc_error(pname,"dqs_mat_qp",2,ierr)
!!  deallocate ( dqr_CRmat_qp , stat=ierr )
!!  call alloc_error(pname,"dqr_CRmat_qp",2,ierr)
!!  deallocate ( dqs_CRmat_qp , stat=ierr )
!!  call alloc_error(pname,"dqs_CRmat_qp",2,ierr)
!!  deallocate ( invM_mat_qp , stat=ierr )
!!  call alloc_error(pname,"invM_mat_qp",2,ierr)
!!  deallocate ( dq1D_mat , stat=ierr )
!!  call alloc_error(pname,"dq1D_mat",2,ierr)
!!  deallocate ( qdtwts_2d , stat=ierr )
!!  call alloc_error(pname,"qdtwts_2d",2,ierr)
!!  deallocate ( vand_qdtr , stat=ierr )
!!  call alloc_error(pname,"vand_qdtr",2,ierr)
!!  deallocate ( vand_qdtr_inv , stat=ierr )
!!  call alloc_error(pname,"vand_qdtr_inv",2,ierr)
!!  deallocate ( interp_s2o , stat=ierr )
!!  call alloc_error(pname,"interp_s2o",2,ierr)
!!  deallocate ( dqr_mat , stat=ierr )
!!  call alloc_error(pname,"dqr_mat",2,ierr)
!!  deallocate ( dqs_mat , stat=ierr )
!!  call alloc_error(pname,"dqs_mat",2,ierr)
!!  deallocate ( dqr_CRmat , stat=ierr )
!!  call alloc_error(pname,"dqr_CRmat",2,ierr)
!!  deallocate ( dqs_CRmat , stat=ierr )
!!  call alloc_error(pname,"dqs_CRmat",2,ierr)
!!  deallocate ( invM_mat , stat=ierr )
!!  call alloc_error(pname,"invM_mat",2,ierr)
!!  deallocate ( qdtwts_1d_qp , stat=ierr )
!!  call alloc_error(pname,"qdtwts_1d_qp",2,ierr)
!!  deallocate ( qdtwts_1d , stat=ierr )
!!  call alloc_error(pname,"qdtwts_1d",2,ierr)
!!  deallocate ( interp_f2g_qp , stat=ierr )
!!  call alloc_error(pname,"interp_f2g_qp",2,ierr)
!!  deallocate ( projct_g2f_qp , stat=ierr )
!!  call alloc_error(pname,"projct_g2f_qp",2,ierr)
!!  deallocate ( interp_s2q_qp , stat=ierr )
!!  call alloc_error(pname,"interp_s2q_qp",2,ierr)
!!  deallocate ( projct_q2s_qp , stat=ierr )
!!  call alloc_error(pname,"projct_q2s_qp",2,ierr)
!!  deallocate ( projct_s2p_qp , stat=ierr )
!!  call alloc_error(pname,"projct_s2p_qp",2,ierr)
!!  deallocate ( interp_s2g_qp , stat=ierr )
!!  call alloc_error(pname,"interp_s2g_qp",2,ierr)
!!  deallocate ( projct_limiter_qp , stat=ierr )
!!  call alloc_error(pname,"projct_limiter_qp",2,ierr)
!!  deallocate ( interp_f2g , stat=ierr )
!!  call alloc_error(pname,"interp_f2g",2,ierr)
!!  deallocate ( projct_g2f , stat=ierr )
!!  call alloc_error(pname,"projct_g2f",2,ierr)
!!  deallocate ( interp_s2q , stat=ierr )
!!  call alloc_error(pname,"interp_s2q",2,ierr)
!!  deallocate ( projct_q2s , stat=ierr )
!!  call alloc_error(pname,"projct_q2s",2,ierr)
!!  deallocate ( projct_s2p , stat=ierr )
!!  call alloc_error(pname,"projct_s2p",2,ierr)
!!  deallocate ( interp_s2g , stat=ierr )
!!  call alloc_error(pname,"interp_s2g",2,ierr)
!!  deallocate ( projct_limiter , stat=ierr )
!!  call alloc_error(pname,"projct_limiter",2,ierr)
!!  !
!!  if (associated(weights)) weights => null()
!!  !
!!  ! Format Statements
!!  !
!!  1 format (/," ************************************************",/, &
!!              "     USING OVER-INTEGRATION FOR DE-ALIASING",/, &
!!              " ************************************************",/)
!!  !
!!end subroutine interp2quadrature
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
        if (allocated(projct(i,j)%toFace)) then
          !
          call array%set( storage_size(projct(i,j)%toFace,kind=inttype) , &
                                  size(projct(i,j)%toFace,kind=inttype) )
          !
          write (array_name,4) Geom_Name(i),j,"toFace"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          do k = lbound(projct(i,j)%toFace,dim=1), &
                 ubound(projct(i,j)%toFace,dim=1)
            !
            array = projct(i,j)%toFace(k)%memory_usage()
            !
            write (array_name,4) Geom_Name(i),j,"toFace",k
            write (iou,2) array_name,array%get_units()
            !
            call dt%add(array)
            !
          end do
          !
        end if
        !
        if (allocated(projct(i,j)%toCell)) then
          !
          call array%set( storage_size(projct(i,j)%toCell,kind=inttype) , &
                                  size(projct(i,j)%toCell,kind=inttype) )
          !
          write (array_name,4) Geom_Name(i),j,"toCell"
          write (iou,2) array_name,array%get_units()
          !
          call dt%add(array)
          !
          do k = lbound(projct(i,j)%toCell,dim=1), &
                 ubound(projct(i,j)%toCell,dim=1)
            !
            array = projct(i,j)%toCell(k)%memory_usage()
            !
            write (array_name,4) Geom_Name(i),j,"toCell",k
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
end subroutine projection_memory_usage
!
!###############################################################################
!
end module projection_mod
!
!###############################################################################
!
