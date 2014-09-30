module startup_mod
  !
  ! Module that contains various functions and routines
  ! used to set up the computational problem
  !
  use module_kind_types
  use polynomial_mod, only : JacobiP => normalized_jacobi_poly
  !
  implicit none
  !
  private
  !
  ! Public Procedures
  !
  public :: get_quadrature_data
  public :: old_get_quadrature_data
  public :: Vandermonde1D
  !
  ! Local logical that controls whether array/debugging information is output
  !
  logical, parameter :: array_debugging = fals
 !logical, parameter :: array_debugging = true
  !
  ! Generic Interfaces for the different Functions
  !
  interface CellVolume1D
    module procedure CellVolume1D_DP, CellVolume1D_QP
  end interface CellVolume1D
  !
  interface CorrectionMatrix1D
    module procedure CorrectionMatrix1D_DP, CorrectionMatrix1D_QP
  end interface CorrectionMatrix1D
  !
  interface CorrectionMatrix2D_quad
    module procedure CorrectionMatrix2D_quad_DP, CorrectionMatrix2D_quad_QP
  end interface CorrectionMatrix2D_quad
  !
  interface CorrectionMatrix2D_quad_dg
    module procedure CorrectionMatrix2D_quad_dg_DP, CorrectionMatrix2D_quad_dg_QP
  end interface CorrectionMatrix2D_quad_dg
  !
  interface CorrectionMatrix2D_tria
    module procedure CorrectionMatrix2D_tria_DP, CorrectionMatrix2D_tria_QP
  end interface CorrectionMatrix2D_tria
  !
  interface d2lagrange
    module procedure d2lagrange_DP, d2lagrange_QP
  end interface d2lagrange
  !
  interface dlagrange
    module procedure dlagrange_DP, dlagrange_QP
  end interface dlagrange
  !
  interface Dmatrix1D
    module procedure Dmatrix1D_DP, Dmatrix1D_QP
  end interface Dmatrix1D
  !
  interface Dmatrix2D_quad
    module procedure Dmatrix2D_quad_DP, Dmatrix2D_quad_QP
  end interface Dmatrix2D_quad
  !
  interface Dmatrix2D_quad_weak
    module procedure Dmatrix2D_quad_weak_DP, Dmatrix2D_quad_weak_QP
  end interface Dmatrix2D_quad_weak
  !
  interface Dmatrix2D_tria
    module procedure Dmatrix2D_tria_DP, Dmatrix2D_tria_QP
  end interface Dmatrix2D_tria
  !
  interface Dmatrix2D_tria_weak
    module procedure Dmatrix2D_tria_weak_DP, Dmatrix2D_tria_weak_QP
  end interface Dmatrix2D_tria_weak
  !
  interface eval_LagrangePoly1D
    module procedure eval_LagrangePoly1D_DP, eval_LagrangePoly1D_QP
  end interface eval_LagrangePoly1D
  !
  interface eval_LagrangePoly2D
    module procedure eval_LagrangePoly2D_DP, eval_LagrangePoly2D_QP
  end interface eval_LagrangePoly2D
  !
  interface Filter2D_quad
    module procedure Filter2D_quad_DP, Filter2D_quad_QP
  end interface Filter2D_quad
  !
  interface Filter2D_tria
    module procedure Filter2D_tria_DP, Filter2D_tria_QP
  end interface Filter2D_tria
  !
  interface GradJacobiP
    module procedure GradJacobiP_DP, GradJacobiP_QP
  end interface GradJacobiP
  !
  interface GradVandermonde1D
    module procedure GradVandermonde1D_DP, GradVandermonde1D_QP
  end interface GradVandermonde1D
  !
  interface GradVandermonde2D
    module procedure GradVandermonde2D_DP, GradVandermonde2D_QP
  end interface GradVandermonde2D
  !
  interface LagrangePoly1D
    module procedure LagrangePoly1D_DP, LagrangePoly1D_QP
  end interface LagrangePoly1D
  !
  interface LagrangePoly2D_quad
    module procedure LagrangePoly2D_quad_DP, LagrangePoly2D_quad_QP
  end interface LagrangePoly2D_quad
  !
  interface LagrangePoly2D_tria
    module procedure LagrangePoly2D_tria_DP, LagrangePoly2D_tria_QP
  end interface LagrangePoly2D_tria
  !
  interface rstoab
    module procedure rstoab_DP, rstoab_QP
  end interface rstoab
  !
  interface Vandermonde1D
    module procedure Vandermonde1D_DP, Vandermonde1D_QP
  end interface Vandermonde1D
  !
  interface Vandermonde1D_N
    module procedure Vandermonde1D_N_DP, Vandermonde1D_N_QP
  end interface Vandermonde1D_N
  !
  interface Vandermonde2D_quad
    module procedure Vandermonde2D_quad_DP, Vandermonde2D_quad_QP
  end interface Vandermonde2D_quad
  !
  interface Vandermonde2D_quad_N
    module procedure Vandermonde2D_quad_N_DP, Vandermonde2D_quad_N_QP
  end interface Vandermonde2D_quad_N
  !
  interface Vandermonde2D_tria
    module procedure Vandermonde2D_tria_DP, Vandermonde2D_tria_QP
  end interface Vandermonde2D_tria
  !
  ! Karniadakis, 1993; Nodes, Modes, and Flow Codes; Physics Today.
  ! "Electronic component speeds and densities have improved by a factor of more
  !  than 10^5 in the last half-century. This development is unrivaled in other
  !  fields of human endeavor; if automobiles had undergone similar improvements,
  !  today a Cadillac would sell for less than a penny, or it would be capable
  !  of a peak speed in excess of 1% of the speed of light, or one gallon of gas
  !  would suffice for about ten trips to the Moon."
  !
contains
!
!###############################################################################
!
subroutine get_quadrature_data()
  !
  !.. Use Statements ..
  use geovar, only : nr,n_order,n_sp1de,n_sp2dt,n_sp2dq
  use geovar, only : solpts_1d,solwts_1d
  use geovar, only : rs,flxpts_1d,interp_s2f,filter
  use geovar, only : vand_1d,vand_1d_inv,vand_inv
  use geovar, only : vand_quad,vand_quad_inv
  use geovar, only : vand_tria,vand_tria_inv
  use metvar, only : cormat,d1_mat,dr_mat,ds_mat,dr_hat,ds_hat,cormat_dg
  !
  use qpvar,  only : d1_mat_qp,dr_mat_qp,ds_mat_qp,dr_hat_qp,ds_hat_qp
  use qpvar,  only : vand_1d_qp,vand_1d_inv_qp
  use qpvar,  only : vand_quad_qp,vand_quad_inv_qp
  use qpvar,  only : vand_tria_qp,vand_tria_inv_qp
  use qpvar,  only : filter_qp,cormat_qp,cormat_dg_qp
  use qpvar,  only : rs_qp,interp_s2f_qp
  use qpvar,  only : solpts_1d_qp,solwts_1d_qp,flxpts_1d_qp
  use qpvar,  only : deallocate_all_qp_arrays
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,n,nf,ierr
  integer :: li,lj,ki,kj
  !
  integer :: nrand
  integer, allocatable, dimension(:) :: rseed
  !
  !.. Local Quad Precision Arrays ..
  real(qp), dimension(1:n_sp1de,1:n_sp1de)      :: vand_temp_qp
  real(qp), dimension(1:nr,1:n_sp2dt)           :: adjd_rs_qp
  real(qp), dimension(1:n_sp1de,1:n_sp1de)      :: grad_vand_1d_qp
  real(qp), dimension(1:n_sp2dt,1:n_sp2dt,1:nr) :: grad_vand_tria_qp
  !
  !.. Local Working Precision Arrays ..
  real(wp), dimension(1:1,1:n_sp2dq)       :: sol,mode
  real(wp), dimension(1:n_sp1de,1:n_sp1de) :: vand_temp
  real(wp), dimension(1:n_sp1de,1:n_sp1de) :: grad_vand_1d
  !
continue
  !
  sol(:,:) = zero
  mode(:,:) = zero
  adjd_rs_qp(:,:) = qzero
  grad_vand_1d_qp(:,:) = qzero
  grad_vand_tria_qp(:,:,:) = qzero
  !
  grad_vand_1d(:,:) = zero
  !
  ! #####################################################################
  ! #####################################################################
  ! ##     USE QUAD PRECISION ARRAYS TO CREATE ALL OF THE ARRAYS       ##
  ! ##   PREVIOUSLY ALLOCATED. WE WILL THEN PERFORM A CHOP TO ZERO     ##
  ! ##  EVERYTHING WITH AN ABSOLUTE VALUE SMALLER THAN EPSILON(0_WP).  ##
  ! #####################################################################
  ! #####################################################################
  !
  ! Allocate the quad precision arrays
  !
  allocate ( vand_1d_qp(1:n_sp1de,1:n_sp1de) , source=qzero , stat=ierr )
  call alloc_error("get_quadrature_data","vand_1d_qp",1,ierr)
  allocate ( vand_tria_qp(1:n_sp2dt,1:n_sp2dt) , source=qzero , stat=ierr )
  call alloc_error("get_quadrature_data","vand_tria_qp",1,ierr)
  allocate ( vand_quad_qp(1:n_sp2dq,1:n_sp2dq) , source=qzero , stat=ierr )
  call alloc_error("get_quadrature_data","vand_quad_qp",1,ierr)
  allocate ( vand_1d_inv_qp(1:n_sp1de,1:n_sp1de) , source=qzero , stat=ierr )
  call alloc_error("get_quadrature_data","vand_1d_inv_qp",1,ierr)
  allocate ( vand_tria_inv_qp(1:n_sp2dt,1:n_sp2dt) , source=qzero , stat=ierr )
  call alloc_error("get_quadrature_data","vand_tria_inv_qp",1,ierr)
  allocate ( vand_quad_inv_qp(1:n_sp2dq,1:n_sp2dq) , source=qzero , stat=ierr )
  call alloc_error("get_quadrature_data","vand_quad_inv_qp",1,ierr)
  !
  allocate ( d1_mat_qp(1:n_sp1de,1:n_sp1de) , source=qzero , stat=ierr )
  call alloc_error("get_quadrature_data","d1_mat_qp",1,ierr)
  allocate ( dr_mat_qp(1:n_sp2dq,1:n_sp2dq,1:2) , source=qzero , stat=ierr )
  call alloc_error("get_quadrature_data","dr_mat_qp",1,ierr)
  allocate ( ds_mat_qp(1:n_sp2dq,1:n_sp2dq,1:2) , source=qzero , stat=ierr )
  call alloc_error("get_quadrature_data","ds_mat_qp",1,ierr)
  allocate ( dr_hat_qp(1:n_sp2dq,1:n_sp2dq,1:2) , source=qzero , stat=ierr )
  call alloc_error("get_quadrature_data","dr_hat_qp",1,ierr)
  allocate ( ds_hat_qp(1:n_sp2dq,1:n_sp2dq,1:2) , source=qzero , stat=ierr )
  call alloc_error("get_quadrature_data","ds_hat_qp",1,ierr)
  !
  allocate ( interp_s2f_qp(1:n_sp2dq,1:4*n_sp1de,1:2) , source=qzero , stat=ierr )
  call alloc_error("get_quadrature_data","interp_s2f_qp",1,ierr)
  !
  allocate ( cormat_qp(1:n_sp2dq,1:4*n_sp1de,1:2) , source=qzero , stat=ierr )
  call alloc_error("get_quadrature_data","cormat_qp",1,ierr)
  allocate ( cormat_dg_qp(1:n_sp2dq,1:4*n_sp1de,1:2) , source=qzero , stat=ierr )
  call alloc_error("get_quadrature_data","cormat_dg_qp",1,ierr)
  !
  allocate ( filter_qp(1:n_sp2dq,1:n_sp2dq,1:2) , source=qzero , stat=ierr )
  call alloc_error("get_quadrature_data","filter_qp",1,ierr)
  !
  ! Adjust r and s to be on the interval [-1,1]. The routines
  ! Vandermonde2d_tria, GradVandermonde2D, and CorrectionMatrix2D_tria
  ! need to be adjusted to work for r and s on the interval [0,1].
  !
  adjd_rs_qp(:,:) = qtwo*rs_qp(:,:) - qone
  adjd_rs_qp = chop_qp( adjd_rs_qp )
  !
  ! 1D Vandermonde matrix and inverse
  !
  vand_1d_qp = Vandermonde1d(solpts_1d_qp)
  vand_1d_qp = chop_qp( vand_1d_qp )
  !
  vand_1d_inv_qp = invert_matrix(vand_1d_qp)
  vand_1d_inv_qp = chop_qp( vand_1d_inv_qp )
  !
  ! 2D triangle Vandermonde matrix and inverse
  !
  vand_tria_qp = Vandermonde2d_tria(adjd_rs_qp)
  vand_tria_qp = chop_qp( vand_tria_qp )
  !
  vand_tria_inv_qp = invert_matrix(vand_tria_qp)
  vand_tria_inv_qp = chop_qp( vand_tria_inv_qp )
  !
  ! 2D quadrilateral Vandermonde matrix and inverse
  !
  vand_quad_qp = Vandermonde2d_quad(solpts_1d_qp)
  vand_quad_qp = chop_qp( vand_quad_qp )
  !
  vand_quad_inv_qp = invert_matrix(vand_quad_qp)
  vand_quad_inv_qp = chop_qp( vand_quad_inv_qp )
  !
  ! Get the derivative of the 1D Vandermonde matrix
  !
  grad_vand_1d_qp = GradVandermonde1D(solpts_1d_qp)
  grad_vand_1d_qp = chop_qp( grad_vand_1d_qp )
  !
  ! Get the derivative of the Vandermonde matrix and store in return_value
  !
  grad_vand_tria_qp = GradVandermonde2D(adjd_rs_qp)
  grad_vand_tria_qp = chop_qp( grad_vand_tria_qp )
  !
  ! Get the 1D derivative matrix
  !
  d1_mat_qp = matmul( grad_vand_1d_qp , vand_1d_inv_qp )
  d1_mat_qp = chop_qp( d1_mat_qp )
  !
  ! Allocation of arrays from the metvar module
  !
  cormat_qp(1:n_sp2dt,1:3*n_sp1de,1) = CorrectionMatrix2D_tria(adjd_rs_qp)
  cormat_qp(1:n_sp2dq,1:4*n_sp1de,2) = CorrectionMatrix2D_quad(solpts_1d_qp)
  cormat_qp = chop_qp( cormat_qp )
  !
  ! Get the derivative matrices for the strong form of the equations
  !
  call Dmatrix2d_tria(grad_vand_tria_qp,vand_tria_inv_qp, &
                      dr_mat_qp(1:n_sp2dt,1:n_sp2dt,1), &
                      ds_mat_qp(1:n_sp2dt,1:n_sp2dt,1))
  !
  call Dmatrix2d_quad(d1_mat_qp,dr_mat_qp(:,:,2),ds_mat_qp(:,:,2))
  dr_mat_qp = chop_qp( dr_mat_qp )
  ds_mat_qp = chop_qp( ds_mat_qp )
  !
  ! Get the derivative matrices for the weak form of the equations
  !
  call Dmatrix2D_tria_weak(adjd_rs_qp,vand_tria_qp, &
                           dr_hat_qp(1:n_sp2dt,1:n_sp2dt,1), &
                           ds_hat_qp(1:n_sp2dt,1:n_sp2dt,1))
  !
  call Dmatrix2d_quad_weak(solwts_1d_qp,d1_mat_qp, &
                           dr_hat_qp(:,:,2),ds_hat_qp(:,:,2))
  dr_hat_qp = chop_qp( dr_hat_qp )
  ds_hat_qp = chop_qp( ds_hat_qp )
  !
  ! Get the Lagrange polynomial to project the variables/fluxes from the
  ! solution points in the cell interior to the flux points on the cell face
  !
  interp_s2f_qp(1:n_sp2dt,1:3*n_sp1de,1) = LagrangePoly2D_tria(rs_qp)
  interp_s2f_qp(1:n_sp2dq,1:4*n_sp1de,2) = LagrangePoly2D_quad(flxpts_1d_qp,solpts_1d_qp)
  interp_s2f_qp = chop_qp( interp_s2f_qp )
  !
  ! Get the correction matrix for DG methods
  !
  cormat_dg_qp(1:n_sp2dq,1:4*n_sp1de,2) = &
              CorrectionMatrix2D_quad_dg(solpts_1d_qp,solwts_1d_qp)
  cormat_dg_qp = chop_qp( cormat_dg_qp )
  !
  ! Get the filter matrices
  !
  filter_qp(1:n_sp2dt,1:n_sp2dt,1) = Filter2D_tria(vand_tria_qp,vand_tria_inv_qp)
  filter_qp(1:n_sp2dq,1:n_sp2dq,2) = Filter2D_quad(vand_quad_qp,vand_quad_inv_qp)
  filter_qp = chop_qp( filter_qp )
  !
  ! Allocate the working precision arrays
  !
  allocate ( vand_1d(1:n_sp1de,1:n_sp1de) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","vand_1d",1,ierr)
  allocate ( vand_tria(1:n_sp2dt,1:n_sp2dt) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","vand_tria",1,ierr)
  allocate ( vand_quad(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","vand_quad",1,ierr)
  allocate ( vand_1d_inv(1:n_sp1de,1:n_sp1de) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","vand_1d_inv",1,ierr)
  allocate ( vand_tria_inv(1:n_sp2dt,1:n_sp2dt) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","vand_tria_inv",1,ierr)
  allocate ( vand_quad_inv(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","vand_quad_inv",1,ierr)
  allocate ( vand_inv(1:n_sp2dq,1:n_sp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","vand_inv",1,ierr)
  !
  allocate ( d1_mat(1:n_sp1de,1:n_sp1de) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","d1_mat",1,ierr)
  allocate ( dr_mat(1:n_sp2dq,1:n_sp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","dr_mat",1,ierr)
  allocate ( ds_mat(1:n_sp2dq,1:n_sp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","ds_mat",1,ierr)
  allocate ( dr_hat(1:n_sp2dq,1:n_sp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","dr_hat",1,ierr)
  allocate ( ds_hat(1:n_sp2dq,1:n_sp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","ds_hat",1,ierr)
  !
  allocate ( interp_s2f(1:n_sp2dq,1:4*n_sp1de,1:2) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","interp_s2f",1,ierr)
  !
  allocate ( cormat(1:n_sp2dq,1:4*n_sp1de,1:2) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","cormat",1,ierr)
  allocate ( cormat_dg(1:n_sp2dq,1:4*n_sp1de,1:2) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","cormat_dg",1,ierr)
  !
  allocate ( filter(1:n_sp2dq,1:n_sp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","filter",1,ierr)
  !
  ! ##############################################
  ! ##############################################
  ! ##   CHOP THE QUAD PRECISION ARRAYS AND     ##
  ! ##  STORE IN THE WORKING PRECISION ARRAYS.  ##
  ! ##############################################
  ! ##############################################
  !
  vand_1d        = chop( vand_1d_qp )
  vand_1d_inv    = chop( vand_1d_inv_qp )
  vand_tria      = chop( vand_tria_qp )
  vand_tria_inv  = chop( vand_tria_inv_qp )
  vand_quad      = chop( vand_quad_qp )
  vand_quad_inv  = chop( vand_quad_inv_qp )
  grad_vand_1d   = chop( grad_vand_1d_qp )
  d1_mat         = chop( d1_mat_qp )
  cormat         = chop( cormat_qp )
  dr_mat         = chop( dr_mat_qp )
  ds_mat         = chop( ds_mat_qp )
  dr_hat         = chop( dr_hat_qp )
  ds_hat         = chop( ds_hat_qp )
  interp_s2f     = chop( interp_s2f_qp )
  cormat_dg      = chop( cormat_dg_qp )
  filter         = chop( filter_qp )
  !
  vand_inv(1:n_sp2dt,1:n_sp2dt,1) = vand_tria_inv
  vand_inv(1:n_sp2dq,1:n_sp2dq,2) = vand_quad_inv
  !
  ! Output some of the arrays if debugging
  !
 !call output2d(true,vand_quad(:,:),2,2,0,'vand_quad')
 !call output2d(true,vand_quad(:,:),2,2,2,'vand_quad2')
  call output2d(array_debugging,vand_quad(:,:),2,2,0,'vand_quad')
 !call output2d(true,vand_quad_inv(:,:),2,2,0,'vand_quad_inv')
 !call output2d(true,vand_quad_inv(:,:),2,2,2,'vand_quad_inv2')
  call output2d(array_debugging,vand_quad_inv(:,:),2,2,0,'vand_quad_inv')
  !
  call output2d(array_debugging,interp_s2f(1:n_sp2dq,1:4*n_sp1de,2),2,1,2,'interp_s2f')
  !
  call output2d(array_debugging,vand_1d(:,:),2,2,0,'vand_1d')
  call output2d(array_debugging,vand_1d_inv(:,:),2,2,0,'vand_1d_inv')
  call output2d(array_debugging,grad_vand_1d(:,:),2,2,0,'grad_vand_1d')
  call output2d(array_debugging,d1_mat(:,:),2,2,0,'d1_mat')
  call output2d(array_debugging,cormat(1:n_sp2dq,1:4*n_sp1de,2),2,1,2,'cormat')
  call output2d(array_debugging,cormat_dg(1:n_sp2dq,1:4*n_sp1de,2),2,1,2,'cormat_dg')
  call output2d(array_debugging,filter(1:n_sp2dq,1:n_sp2dq,2),2,2,2,'filter')
  !
  ! Output the solution and flux points to a tecplot file for visualization
  !
 !call write_quadrature
  !
  call gauss_lobatto_matrices
  !
  call interp2quadrature
  !
  ! #################################################################
  ! #################################################################
  ! ##    WE NO LONGER NEED ANY OF THE QUAD PRECISION ARRAYS SO    ##
  ! ##  CALL THE SUBROUTINE TO DEALLOCATE ALL THAT ARE ALLOCATED.  ##
  ! #################################################################
  ! #################################################################
  !
  call deallocate_all_qp_arrays()
  !
 !write (101,*) "Quad Vandermonde Matrix"
 !1 format (a,i0,a,i0,a)
 !2 format (*(es12.4))
 !do lj = 0,n_order
 !  do li = 0,n_order
 !    l = n_sp1de*lj + li + 1
 !    write (101,*)
 !    write (101,1) "Basis function for term (x^",li,")*(y^",lj,")"
 !    do kj = n_sp1de,1,-1
 !      write (101,2) (vand_quad(n_sp1de*(kj-1)+ki,l),ki=1,n_sp1de)
 !    end do
 !  end do
 !end do
 !!
 !!
 !write (100,*) "Inverse of Quad Vandermonde Matrix"
 !do lj = 0,n_order
 !  do li = 0,n_order
 !    l = n_sp1de*lj + li + 1
 !    write (100,*)
 !    write (100,1) "Basis function for term (x^",li,")*(y^",lj,")"
 !    do kj = n_sp1de,1,-1
 !      write (100,2) (vand_quad_inv(n_sp1de*(kj-1)+ki,l),ki=1,n_sp1de)
 !    end do
 !  end do
 !end do
 !!
 !!
 !write (101,*) "Correction Function"
 !write (102,*) "DG Correction Function"
 !do nf = 1,4*n_sp1de
 !  !
 !  if (nf <= 1*n_sp1de) then
 !    write (101,5) nf
 !    write (102,5) nf
 !  else if (nf <= 2*n_sp1de) then
 !    write (101,6) nf - n_sp1de
 !    write (102,6) nf - n_sp1de
 !  else if (nf <= 3*n_sp1de) then
 !    write (101,7) n_sp1de+1-(nf-2*n_sp1de)
 !    write (102,7) n_sp1de+1-(nf-2*n_sp1de)
 !  else
 !    write (101,8) n_sp1de+1-(nf-3*n_sp1de)
 !    write (102,8) n_sp1de+1-(nf-3*n_sp1de)
 !  end if
 !  do j = n_sp1de,1,-1
 !    write (101,4) (cormat((j-1)*n_sp1de+i,nf,2),i=1,n_sp1de)
 !    write (102,4) (cormat_dg((j-1)*n_sp1de+i,nf,2),i=1,n_sp1de)
 !  end do
 !  !
 !end do
 !4 format (*(es12.4))
 !5 format (/," Bottom/South edge for xi = ",i0)
 !6 format (/," Right/East edge for eta = ",i0)
 !7 format (/," Top/North edge for xi = ",i0)
 !8 format (/," Left/West edge for eta = ",i0)
 !do j = 1,3*n_sp1de
 !  write (201,1) (cormat(i,j,1),i=1,n_sp2dt)
 !end do
 !write (102,*) "Dr matrix"
 !write (103,*) "Ds matrix"
 !do j = 1,n_sp2dq
 !  write (102,1) (dr_mat(i,j,2),i=1,n_sp2dq)
 !  write (103,1) (ds_mat(i,j,2),i=1,n_sp2dq)
 !end do
 !write (106,*) "Dr matrix"
 !write (107,*) "Ds matrix"
 !do l = 1,n_sp1de
 !  do k = 1,n_sp1de
 !    write (106,2) k,l
 !    write (107,2) k,l
 !    do j = 1,n_sp1de
 !      write (106,3) (dr_mat( (n_sp1de-j)*n_sp1de+i , (l-1)*n_sp1de+k , 2 ),i=1,n_sp1de)
 !      write (107,3) (ds_mat( (n_sp1de-j)*n_sp1de+i , (l-1)*n_sp1de+k , 2 ),i=1,n_sp1de)
 !    end do
 !  end do
 !end do
 !2 format (/,"i = ",i2,", j = ",i2)
 !3 format (*(es12.4))
 !!
 !write (110,*) "interp_s2f"
 !do l = 1,4
 !  do k = 1,n_sp1de
 !    write (110,4) l,k
 !    do j = 1,n_sp1de
 !      write (110,3) (interp_s2f( (n_sp1de-j)*n_sp1de+i , (l-1)*n_sp1de+k , 2 ),i=1,n_sp1de)
 !    end do
 !  end do
 !end do
 !4 format (/,"Face = ",i2,", i = ",i2)
 !!
 !write (108,*) "1D Dr"
 !write (109,*) "1D Ds"
 !do j = 1,n_sp1de
 !  write (108,3) (d1_mat(i,n_sp1de+1-j),i=1,n_sp1de)
 !  write (109,3) (d1_mat(n_sp1de+1-j,i),i=1,n_sp1de)
 !end do
 !write (*,*)
 !write (*,*)
 !do i = 1,n_sp1de
 !  write (*,3) (dlagrange(k,solpts_1d(i),solpts_1d),k=1,n_sp1de)
 !end do
 !write (*,*)
 !write (*,*)
 !!
 !do j = 1,n_sp2dt
 !  write (202,1) (dr_mat(i,j,1),i=1,n_sp2dt)
 !  write (203,1) (ds_mat(i,j,1),i=1,n_sp2dt)
 !end do
 !write (105,*) "interp_s2f matrix"
 !do j = 1,4*n_sp1de
 !  write (*,1) (interp_s2f(i,j,2),i=1,n_sp2dq)
 !end do
 !!
 !1 format (<n_sp2dq>es12.4)
 !!
 !call stop_recon(stop_mpi,"get_quadrature_data")
 !!
end subroutine get_quadrature_data
!
!###############################################################################
!
subroutine interp2quadrature()
  !
  !.. Use Statements ..
  use geovar, only : n_order,p_order,o_order,n_sp1de,n_sp2dt,n_sp2dq
  use geovar, only : q_order,n_qp1de,n_qp2dt,n_qp2dq
  use geovar, only : qdtpts_1d,qdtwts_1d,qdtwts_2d
  use geovar, only : invM_mat,interp_s2o
  use geovar, only : interp_f2g,projct_g2f,vand_qdtr
  use geovar, only : interp_s2q,projct_q2s,vand_qdtr_inv
  use geovar, only : interp_s2g,using_quadrature,projct_s2p
  use ovar,   only : iop_qdtpts
  use metvar, only : dr_mat,ds_mat,dq1D_mat
  use metvar, only : dqr_mat,dqr_LPmat,dqr_CRmat
  use metvar, only : dqs_mat,dqs_LPmat,dqs_CRmat
  !
  use qpvar, only : qdtpts_1d_qp,qdtwts_1d_qp,qdtwts_2d_qp
  use qpvar, only : solpts_1d_qp,solwts_1d_qp,invM_mat_qp
  use qpvar, only : gauss_1d_qp,lobatto_1d_qp,radau_1d_qp
  use qpvar, only : gauss_wts_2d_qp,lobatto_wts_2d_qp,radau_wts_2d_qp
  use qpvar, only : vand_1d_qp,vand_1d_inv_qp,vand_quad_qp,vand_quad_inv_qp
  use qpvar, only : interp_f2g_qp,projct_g2f_qp,vand_qdtr_qp
  use qpvar, only : interp_s2q_qp,projct_q2s_qp,vand_qdtr_inv_qp
  use qpvar, only : interp_s2g_qp,projct_s2p_qp,weights2d_qp
  use qpvar, only : interp_s2o_qp
  use qpvar, only : dr_mat_qp,ds_mat_qp,dq1D_mat_qp
  use qpvar, only : dqr_mat_qp,dqr_LPmat_qp,dqr_CRmat_qp
  use qpvar, only : dqs_mat_qp,dqs_LPmat_qp,dqs_CRmat_qp
  !
  use polynomial_mod, only : nodes_legendre_gauss, &
                             weights_legendre_gauss
  use polynomial_mod, only : nodes_legendre_gauss_lobatto, &
                             weights_legendre_gauss_lobatto
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,ierr
  integer :: li,lj,ki,kj,nsp,nqp,npp,nop
  !
  !.. Local Arrays ..
  real(qp), dimension(1:n_qp1de) :: vp_qdtr_qp
  real(qp), dimension(1:abs(o_order)+1) :: vp_out_qp
  real(qp), dimension(1:abs(o_order)+1) :: outpts_1d_qp
  real(qp), dimension(1:n_qp1de,1:n_qp1de) :: vand1D_qdtr_qp
  real(qp), dimension(1:n_qp1de,1:n_qp1de) :: vand1D_qdtr_inv_qp
  real(qp), dimension(1:n_qp1de,1:n_qp1de) :: grad_vand1D_qdtr_qp
  real(qp), dimension(1:n_qp1de,1:n_sp1de) :: vand1D_f2g_qp
  real(qp), dimension(1:n_sp2dq,1:n_sp2dq) :: wt_mat_qp
  real(qp), dimension(1:n_qp2dq,1:n_sp2dq,1:2) :: vand2D_s2q_qp
  real(qp), dimension(1:n_sp2dq,1:(p_order+1)**2,1:2) :: vand2D_p2s_qp
  real(qp), dimension(1:(abs(o_order)+1)**2,1:n_sp2dq,1:2) :: vand2D_s2o_qp
  real(qp), dimension(1:(abs(o_order)+1)**2,1:(p_order+1)**2,1:2) :: vand2D_p2o_qp
  !
continue
  !
  vp_qdtr_qp(:) = qzero
  vp_out_qp(:) = qzero
  outpts_1d_qp(:) = qzero
  vand1D_qdtr_qp(:,:) = qzero
  vand1D_qdtr_inv_qp(:,:) = qzero
  grad_vand1D_qdtr_qp(:,:) = qzero
  vand1D_f2g_qp(:,:) = qzero
  vand2D_s2q_qp(:,:,:) = qzero
  vand2D_p2s_qp(:,:,:) = qzero
  vand2D_s2o_qp(:,:,:) = qzero
  vand2D_p2o_qp(:,:,:) = qzero
  wt_mat_qp(:,:) = qzero
  !
  npp = (p_order+1)**2
  nop = (abs(o_order)+1)**2
  !
  if (using_quadrature .and. mypnum == 0) write (iout,1)
  !
  ! ###########################################################
  ! ###########################################################
  ! ##                                                       ##
  ! ##  Allocate all the arrays needed for quadrature in 2D  ##
  ! ##                                                       ##
  ! ###########################################################
  ! ###########################################################
  !
  ! #####################################################################
  ! #####################################################################
  ! ##     USE QUAD PRECISION ARRAYS TO CREATE ALL OF THE ARRAYS       ##
  ! ##   PREVIOUSLY ALLOCATED. WE WILL THEN PERFORM A CHOP TO ZERO     ##
  ! ##  EVERYTHING WITH AN ABSOLUTE VALUE SMALLER THAN EPSILON(0_WP).  ##
  ! #####################################################################
  ! #####################################################################
  !
  ! Allocate the quad precision arrays
  !
  allocate ( qdtpts_1d_qp(1:n_qp1de) , source=qzero , stat=ierr )
  call alloc_error("interp2quadrature","qdtpts_1d_qp",1,ierr)
  allocate ( qdtwts_1d_qp(1:n_qp1de) , source=qzero , stat=ierr )
  call alloc_error("interp2quadrature","qdtwts_1d_qp",1,ierr)
  allocate ( dq1D_mat_qp(1:n_qp1de,1:n_qp1de), source=qzero , stat=ierr )
  call alloc_error("interp2quadrature","dq1D_mat_qp",1,ierr)
  allocate ( interp_f2g_qp(1:n_sp1de,1:n_qp1de) , source=qzero , stat=ierr )
  call alloc_error("interp2quadrature","interp_f2g_qp",1,ierr)
  allocate ( projct_g2f_qp(1:n_qp1de,1:n_sp1de) , source=qzero , stat=ierr )
  call alloc_error("interp2quadrature","projct_g2f_qp",1,ierr)
  !
  allocate ( qdtwts_2d_qp(1:n_qp2dq,1:2) , source=qzero , stat=ierr )
  call alloc_error("interp2quadrature","qdtwts_2d_qp",1,ierr)
  allocate ( vand_qdtr_qp(1:n_qp2dq,1:n_qp2dq,1:2), source=qzero , stat=ierr )
  call alloc_error("interp2quadrature","vand_qdtr_qp",1,ierr)
  allocate ( vand_qdtr_inv_qp(1:n_qp2dq,1:n_qp2dq,1:2), source=qzero , stat=ierr )
  call alloc_error("interp2quadrature","vand_qdtr_inv_qp",1,ierr)
  allocate ( interp_s2q_qp(1:n_sp2dq,1:n_qp2dq,1:2) , source=qzero , stat=ierr )
  call alloc_error("interp2quadrature","interp_s2q_qp",1,ierr)
  allocate ( interp_s2o_qp(1:n_sp2dq,1:nop,1:2,1:3) , source=qzero , stat=ierr )
  call alloc_error("interp2quadrature","interp_s2o_qp",1,ierr)
  allocate ( projct_q2s_qp(1:n_qp2dq,1:n_sp2dq,1:2) , source=qzero , stat=ierr )
  call alloc_error("interp2quadrature","projct_q2s_qp",1,ierr)
  allocate ( projct_s2p_qp(1:n_sp2dq,1:n_sp2dq,1:2,1:5) , source=qzero , stat=ierr )
  call alloc_error("interp2quadrature","projct_s2p_qp",1,ierr)
  allocate ( dqr_mat_qp(1:n_qp2dq,1:n_qp2dq,1:2), source=qzero , stat=ierr )
  call alloc_error("interp2quadrature","dqr_mat_qp",1,ierr)
  allocate ( dqs_mat_qp(1:n_qp2dq,1:n_qp2dq,1:2), source=qzero , stat=ierr )
  call alloc_error("interp2quadrature","dqs_mat_qp",1,ierr)
  allocate ( dqr_LPmat_qp(1:n_qp2dq,1:n_sp2dq,1:2), source=qzero , stat=ierr )
  call alloc_error("interp2quadrature","dqr_LPmat_qp",1,ierr)
  allocate ( dqs_LPmat_qp(1:n_qp2dq,1:n_sp2dq,1:2), source=qzero , stat=ierr )
  call alloc_error("interp2quadrature","dqs_LPmat_qp",1,ierr)
  allocate ( dqr_CRmat_qp(1:n_sp2dq,1:n_qp2dq,1:2), source=qzero , stat=ierr )
  call alloc_error("interp2quadrature","dqr_CRmat_qp",1,ierr)
  allocate ( dqs_CRmat_qp(1:n_sp2dq,1:n_qp2dq,1:2), source=qzero , stat=ierr )
  call alloc_error("interp2quadrature","dqs_CRmat_qp",1,ierr)
  allocate ( invM_mat_qp(1:n_sp2dq,1:n_sp2dq,1:2) , source=qzero , stat=ierr )
  call alloc_error("interp2quadrature","invM_mat_qp",1,ierr)
  allocate ( interp_s2g_qp(1:n_sp2dq,1:4*n_qp1de,1:2) , source=qzero , stat=ierr )
  call alloc_error("interp2quadrature","interp_s2g_qp",1,ierr)
  !
  ! ####################################
  ! ####################################
  ! ##                                ##
  ! ##  Section for quadrature in 1D  ##
  ! ##                                ##
  ! ####################################
  ! ####################################
  !
  ! 1D quadrature nodes and weights
  !
  if (iop_qdtpts == Legendre_Gauss_Lobatto) then
    call nodes_legendre_gauss_lobatto(q_order,qdtpts_1d_qp,vp_qdtr_qp)
    call weights_legendre_gauss_lobatto(q_order,vp_qdtr_qp,qdtwts_1d_qp)
  else
    call nodes_legendre_gauss(n_qp1de,qdtpts_1d_qp,vp_qdtr_qp)
    call weights_legendre_gauss(n_qp1de,qdtpts_1d_qp,vp_qdtr_qp,qdtwts_1d_qp)
  end if
  qdtpts_1d_qp = chop_qp( qdtpts_1d_qp )
  qdtwts_1d_qp = chop_qp( qdtwts_1d_qp )
  !
  ! 1D output nodes
  !
  call nodes_legendre_gauss_lobatto(abs(o_order),outpts_1d_qp,vp_out_qp)
  outpts_1d_qp = chop_qp( outpts_1d_qp )
  !
  ! 1D Derivative matrix for the quadrature points
  !
  vand1D_qdtr_qp(:,:) = Vandermonde1D(qdtpts_1d_qp)
  vand1D_qdtr_qp = chop_qp( vand1D_qdtr_qp )
  !
  vand1D_qdtr_inv_qp(:,:) = invert_matrix( vand1D_qdtr_qp )
  vand1D_qdtr_inv_qp = chop_qp( vand1D_qdtr_inv_qp )
  !
  grad_vand1D_qdtr_qp(:,:) = GradVandermonde1D(qdtpts_1d_qp)
  grad_vand1D_qdtr_qp = chop_qp( grad_vand1D_qdtr_qp )
  !
  dq1D_mat_qp(:,:) = matmul( grad_vand1D_qdtr_qp , vand1D_qdtr_inv_qp )
  dq1D_mat_qp = chop_qp( dq1D_mat_qp )
  !
  ! 1D Vandermonde matrix giving the value for each of the
  ! #(n_sp1de) modes at the #(n_qp1de) quadrature points
  !
  vand1d_f2g_qp(:,:) = Vandermonde1D_N(qdtpts_1d_qp,n_sp1de)
  vand1d_f2g_qp = chop_qp( vand1d_f2g_qp )
  !
  ! Matrix for interpolating 1D nodal values at face
  ! flux/solution points to face quadrature points
  !
  interp_f2g_qp(:,:) = transpose( matmul( vand1d_f2g_qp , vand_1d_inv_qp ) )
  interp_f2g_qp = chop_qp( interp_f2g_qp )
  !
  ! Matrix for projecting nodal values at face quadrature
  ! points to the face flux/solution points
  !
  projct_g2f_qp(:,:) = matmul( vand1D_f2g_qp , transpose(vand_1D_qp) )
  projct_g2f_qp = chop_qp( projct_g2f_qp )
  do i = 1,n_sp1de
    projct_g2f_qp(:,i) = projct_g2f_qp(:,i) * qdtwts_1d_qp(:)
  end do
  projct_g2f_qp = chop_qp( projct_g2f_qp )
  !
  !
  ! ################################################
  ! ################################################
  ! ##                                            ##
  ! ##  Section for quadrature on triangle cells  ##
  ! ##                                            ##
  ! ################################################
  ! ################################################
  !
  nsp = n_sp2dt
  nqp = n_qp2dt
  !
  ! #####################################################
  ! #####################################################
  ! ##                                                 ##
  ! ##  Section for quadrature on quadrilateral cells  ##
  ! ##                                                 ##
  ! #####################################################
  ! #####################################################
  !
  nsp = n_sp2dq
  nqp = n_qp2dq
  !
  ! 2D quadrature weights for the quadrature points
  !
  k = 0
  do i = 1,n_qp1de
    do j = 1,n_qp1de
      k = k + 1
      qdtwts_2d_qp(k,2) = qdtwts_1d_qp(i)*qdtwts_1d_qp(j)
    end do
  end do
  qdtwts_2d_qp = chop_qp( qdtwts_2d_qp )
  !
  do i = 1,n_qp2dq
    wt_mat_qp(i,i) = qdtwts_2d_qp(i,2)
  end do
  wt_mat_qp = chop_qp( wt_mat_qp )
  !
  ! Full Vandermonde matrix for the quadrature nodes and modes
  !
  vand_qdtr_qp(1:nqp,1:nqp,2) = Vandermonde2D_quad(qdtpts_1d_qp)
  vand_qdtr_qp = chop_qp( vand_qdtr_qp )
  !
  vand_qdtr_inv_qp(1:nqp,1:nqp,2) = invert_matrix(vand_qdtr_qp(1:nqp,1:nqp,2))
  vand_qdtr_inv_qp = chop_qp( vand_qdtr_inv_qp )
  !
  ! ###########################################
  !  2D matrix for project solution to p_order
  ! ###########################################
  !
  ! Project at solution points
  !
  vand2D_p2s_qp(1:nsp,1:npp,2) = Vandermonde2D_quad_N(solpts_1d_qp,p_order)
  vand2D_p2s_qp = chop_qp( vand2D_p2s_qp )
  !
  projct_s2p_qp(1:nsp,1:nsp,2,1) = &
    matmul( vand2D_p2s_qp(1:nsp,1:npp,2) , transpose(vand2D_p2s_qp(1:nsp,1:npp,2)) )
  projct_s2p_qp = chop_qp( projct_s2p_qp )
  do i = 1,nsp
    projct_s2p_qp(i,1:nsp,2,1) = projct_s2p_qp(i,1:nsp,2,1) * weights2d_qp(i,2)
  end do
  projct_s2p_qp = chop_qp( projct_s2p_qp )
  !
  ! Project at Gauss points
  !
  vand2D_p2s_qp(1:nsp,1:npp,2) = Vandermonde2D_quad_N(gauss_1d_qp,p_order)
  vand2D_p2s_qp = chop_qp( vand2D_p2s_qp )
  !
  projct_s2p_qp(1:nsp,1:nsp,2,2) = &
    matmul( vand2D_p2s_qp(1:nsp,1:npp,2) , transpose(vand2D_p2s_qp(1:nsp,1:npp,2)) )
  projct_s2p_qp = chop_qp( projct_s2p_qp )
  do i = 1,nsp
    projct_s2p_qp(i,1:nsp,2,2) = projct_s2p_qp(i,1:nsp,2,2) * gauss_wts_2d_qp(i)
  end do
  projct_s2p_qp = chop_qp( projct_s2p_qp )
  !
  ! Project at Lobatto points
  !
  vand2D_p2s_qp(1:nsp,1:npp,2) = Vandermonde2D_quad_N(lobatto_1d_qp,p_order)
  vand2D_p2s_qp = chop_qp( vand2D_p2s_qp )
  !
  projct_s2p_qp(1:nsp,1:nsp,2,3) = &
    matmul( vand2D_p2s_qp(1:nsp,1:npp,2) , transpose(vand2D_p2s_qp(1:nsp,1:npp,2)) )
  projct_s2p_qp = chop_qp( projct_s2p_qp )
  do i = 1,nsp
    projct_s2p_qp(i,1:nsp,2,3) = projct_s2p_qp(i,1:nsp,2,3) * lobatto_wts_2d_qp(i)
  end do
  projct_s2p_qp = chop_qp( projct_s2p_qp )
  !
  ! Project at Radau points
  !
  vand2D_p2s_qp(1:nsp,1:npp,2) = Vandermonde2D_quad_N(radau_1d_qp,p_order)
  vand2D_p2s_qp = chop_qp( vand2D_p2s_qp )
  !
  projct_s2p_qp(1:nsp,1:nsp,2,4) = &
    matmul( vand2D_p2s_qp(1:nsp,1:npp,2) , transpose(vand2D_p2s_qp(1:nsp,1:npp,2)) )
  projct_s2p_qp = chop_qp( projct_s2p_qp )
  do i = 1,nsp
    projct_s2p_qp(i,1:nsp,2,4) = projct_s2p_qp(i,1:nsp,2,4) * radau_wts_2d_qp(i)
  end do
  projct_s2p_qp = chop_qp( projct_s2p_qp )
  !
  projct_s2p_qp(1:nsp,1:nsp,2,5) = projection_matrix(vand2d_p2s_qp(1:nsp,1:npp,2), &
                                                     wt_mat_qp(1:nsp,1:nsp))
  !
  !
  ! ###########################################
  ! 2D matrix for interpolating nodal values at SPs to QPs
  ! ###########################################
  !
  vand2D_s2q_qp(1:nqp,1:nsp,2) = Vandermonde2D_quad_N(qdtpts_1d_qp,n_order)
  vand2D_s2q_qp = chop_qp( vand2D_s2q_qp )
  !
  interp_s2q_qp(1:nsp,1:nqp,2) = &
        transpose( matmul( vand2D_s2q_qp(1:nqp,1:nsp,2) , vand_quad_inv_qp ) )
  interp_s2q_qp = chop_qp( interp_s2q_qp )
  !
  ! ###########################################
  ! Matrix for interpolating cell nodal values at SPs to face QPs
  ! ###########################################
  !
  interp_s2g_qp(1:nsp,1:4*n_qp1de,2) = LagrangePoly2D_quad(qdtpts_1d_qp,solpts_1d_qp)
  interp_s2g_qp = chop_qp( interp_s2g_qp )
  !
  ! ###########################################
  ! Matrix for interpolating cell nodal values at SPs to output SPs
  ! ###########################################
  !
  vand2D_p2o_qp(1:nop,1:npp,2) = Vandermonde2D_quad_N(outpts_1d_qp,p_order)
  vand2D_p2o_qp = chop_qp( vand2D_p2o_qp )
 !call output2d(true,real(vand2D_p2o_qp(1:nop,1:npp,2),kind=wp),2,2,2,'vand2d_p2o')
 !call output2d(true,real(vand2D_p2s_qp(1:nsp,1:npp,2),kind=wp),2,2,2,'vand2d_p2s')
  !
  interp_s2o_qp(1:nsp,1:nop,2,1) = interpolation_matrix(vand2d_p2s_qp(1:nsp,1:npp,2), &
                                                        vand2d_p2o_qp(1:nop,1:npp,2), &
                                                        wt_mat_qp(1:nsp,1:nsp))
  !
  vand2D_s2o_qp(1:nop,1:nsp,2) = Vandermonde2D_quad_N(outpts_1d_qp,n_order)
  vand2D_s2o_qp = chop_qp( vand2D_s2o_qp )
  !
  interp_s2o_qp(1:nsp,1:nop,2,2) = &
        transpose( matmul( vand2D_s2o_qp(1:nop,1:nsp,2) , vand_quad_inv_qp ) )
  interp_s2o_qp = chop_qp( interp_s2o_qp )
  !
  interp_s2o_qp(1:nsp,1:nop,2,3) = &
   matmul( vand2D_p2s_qp(1:nsp,1:npp,2) , transpose(vand2d_p2o_qp(1:nop,1:npp,2)) )
  !transpose(matmul( vand2D_p2o_qp(1:nop,1:npp,2) , transpose(vand2d_p2s_qp(1:nsp,1:npp,2)) ))
  interp_s2o_qp = chop_qp( interp_s2o_qp )
  do i = 1,nop
    interp_s2o_qp(1:nsp,i,2,3) = interp_s2o_qp(1:nsp,i,2,3)*weights2d_qp(1:nsp,2)
  end do
  interp_s2o_qp = chop_qp( interp_s2o_qp )
  !
  ! ###########################################
  ! 2D matrix to project nodal values at QPs to SPs
  ! ###########################################
  !
  projct_q2s_qp(1:nqp,1:nsp,2) = &
        matmul( vand2D_s2q_qp(1:nqp,1:nsp,2) , transpose(vand_quad_qp) )
  projct_q2s_qp = chop_qp( projct_q2s_qp )
  do i = 1,nsp
    projct_q2s_qp(1:nqp,i,2) = projct_q2s_qp(1:nqp,i,2) * qdtwts_2d_qp(1:nqp,2)
  end do
  projct_q2s_qp = chop_qp( projct_q2s_qp )
  !
  ! ###########################################
  ! 2D Derivative matrix for the quadrature points
  ! ###########################################
  !
  call Dmatrix2D_quad(dq1d_mat_qp,dqr_mat_qp(1:nqp,1:nqp,2),dqs_mat_qp(1:nqp,1:nqp,2))
  dqr_mat_qp = chop_qp( dqr_mat_qp )
  dqs_mat_qp = chop_qp( dqs_mat_qp )
  !
  ! Pre-multiply the derivative matrices by the projection
  ! matrix for when using the derivative of the Lagrange polynomial
  !
  dqr_LPmat_qp(1:nqp,1:nsp,2) = &
      matmul( transpose(dqr_mat_qp(1:nqp,1:nqp,2)) , projct_q2s_qp(1:nqp,1:nsp,2) )
  dqs_LPmat_qp(1:nqp,1:nsp,2) = &
      matmul( transpose(dqs_mat_qp(1:nqp,1:nqp,2)) , projct_q2s_qp(1:nqp,1:nsp,2) )
  dqr_LPmat_qp = chop_qp( dqr_LPmat_qp )
  dqs_LPmat_qp = chop_qp( dqs_LPmat_qp )
  !
  ! Pre-multiply the derivative matrices by the
  ! interpolation matrix for when using the chain-rule
  !
  dqr_CRmat_qp(1:nsp,1:nqp,2) = &
      matmul( interp_s2q_qp(1:nsp,1:nqp,2) , transpose(dqr_mat_qp(1:nqp,1:nqp,2)) )
  dqs_CRmat_qp(1:nsp,1:nqp,2) = &
      matmul( interp_s2q_qp(1:nsp,1:nqp,2) , transpose(dqs_mat_qp(1:nqp,1:nqp,2)) )
  dqr_CRmat_qp = chop_qp( dqr_CRmat_qp )
  dqs_CRmat_qp = chop_qp( dqs_CRmat_qp )
  !
  ! Inverse of the mass matrix
  !
  invM_mat_qp(1:nsp,1:nsp,2) = matmul( vand_quad_qp , transpose(vand_quad_qp) )
  invM_mat_qp = chop_qp( invM_mat_qp )
  !
  ! ##############################################
  ! ##############################################
  ! ##   CHOP THE QUAD PRECISION ARRAYS AND     ##
  ! ##  STORE IN THE WORKING PRECISION ARRAYS.  ##
  ! ##############################################
  ! ##############################################
  !
  ! Allocate the working precision arrays
  !
  allocate ( qdtpts_1d(1:n_qp1de) , source=zero , stat=ierr )
  call alloc_error("interp2quadrature","qdtpts_1d",1,ierr)
  allocate ( qdtwts_1d(1:n_qp1de) , source=zero , stat=ierr )
  call alloc_error("interp2quadrature","qdtwts_1d",1,ierr)
  allocate ( dq1D_mat(1:n_qp1de,1:n_qp1de), source=zero , stat=ierr )
  call alloc_error("interp2quadrature","dq1D_mat",1,ierr)
  allocate ( interp_f2g(1:n_sp1de,1:n_qp1de) , source=zero , stat=ierr )
  call alloc_error("interp2quadrature","interp_f2g",1,ierr)
  allocate ( projct_g2f(1:n_qp1de,1:n_sp1de) , source=zero , stat=ierr )
  call alloc_error("interp2quadrature","projct_g2f",1,ierr)
  !
  allocate ( qdtwts_2d(1:n_qp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("interp2quadrature","qdtwts_2d",1,ierr)
  allocate ( vand_qdtr(1:n_qp2dq,1:n_qp2dq,1:2), source=zero , stat=ierr )
  call alloc_error("interp2quadrature","vand_qdtr",1,ierr)
  allocate ( vand_qdtr_inv(1:n_qp2dq,1:n_qp2dq,1:2), source=zero , stat=ierr )
  call alloc_error("interp2quadrature","vand_qdtr_inv",1,ierr)
  allocate ( interp_s2q(1:n_sp2dq,1:n_qp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("interp2quadrature","interp_s2q",1,ierr)
  allocate ( interp_s2o(1:n_sp2dq,1:nop,1:2,1:3) , source=zero , stat=ierr )
  call alloc_error("interp2quadrature","interp_s2o",1,ierr)
  allocate ( projct_q2s(1:n_qp2dq,1:n_sp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("interp2quadrature","projct_q2s",1,ierr)
  allocate ( projct_s2p(1:n_sp2dq,1:n_sp2dq,1:2,1:5) , source=zero , stat=ierr )
  call alloc_error("interp2quadrature","projct_s2p",1,ierr)
  allocate ( dqr_mat(1:n_qp2dq,1:n_qp2dq,1:2), source=zero , stat=ierr )
  call alloc_error("interp2quadrature","dqr_mat",1,ierr)
  allocate ( dqs_mat(1:n_qp2dq,1:n_qp2dq,1:2), source=zero , stat=ierr )
  call alloc_error("interp2quadrature","dqs_mat",1,ierr)
  allocate ( dqr_LPmat(1:n_qp2dq,1:n_sp2dq,1:2), source=zero , stat=ierr )
  call alloc_error("interp2quadrature","dqr_LPmat",1,ierr)
  allocate ( dqs_LPmat(1:n_qp2dq,1:n_sp2dq,1:2), source=zero , stat=ierr )
  call alloc_error("interp2quadrature","dqs_LPmat",1,ierr)
  allocate ( dqr_CRmat(1:n_sp2dq,1:n_qp2dq,1:2), source=zero , stat=ierr )
  call alloc_error("interp2quadrature","dqr_CRmat",1,ierr)
  allocate ( dqs_CRmat(1:n_sp2dq,1:n_qp2dq,1:2), source=zero , stat=ierr )
  call alloc_error("interp2quadrature","dqs_CRmat",1,ierr)
  allocate ( invM_mat(1:n_sp2dq,1:n_sp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("interp2quadrature","invM_mat",1,ierr)
  allocate ( interp_s2g(1:n_sp2dq,1:4*n_qp1de,1:2) , source=zero , stat=ierr )
  call alloc_error("interp2quadrature","interp_s2g",1,ierr)
  !
  qdtpts_1d     = chop( qdtpts_1d_qp )
  qdtwts_1d     = chop( qdtwts_1d_qp )
  dq1D_mat      = chop( dq1D_mat_qp )
  interp_f2g    = chop( interp_f2g_qp )
  projct_g2f    = chop( projct_g2f_qp )
  qdtwts_2d     = chop( qdtwts_2d_qp )
  vand_qdtr     = chop( vand_qdtr_qp )
  vand_qdtr_inv = chop( vand_qdtr_inv_qp )
  interp_s2q    = chop( interp_s2q_qp )
  interp_s2o    = chop( interp_s2o_qp )
  interp_s2g    = chop( interp_s2g_qp )
  projct_q2s    = chop( projct_q2s_qp )
  projct_s2p    = chop( projct_s2p_qp )
  dqr_mat       = chop( dqr_mat_qp )
  dqs_mat       = chop( dqs_mat_qp )
  dqr_LPmat     = chop( dqr_LPmat_qp )
  dqs_LPmat     = chop( dqs_LPmat_qp )
  dqr_CRmat     = chop( dqr_CRmat_qp )
  dqs_CRmat     = chop( dqs_CRmat_qp )
  invM_mat      = chop( invM_mat_qp )
  !
  call output2d(array_debugging,interp_f2g(:,:),2,2,3,'interp_f2g')
  call output2d(array_debugging,projct_g2f(:,:),2,2,3,'projct_g2f')
  call output2d(array_debugging,interp_s2q(1:nsp,1:nqp,2),2,2,2,'interp_s2q')
  call output2d(array_debugging,interp_s2o(1:nsp,1:nqp,2,1),2,2,2,'interp_s2o')
 !call output2d(true,interp_s2o(1:nsp,1:nop,2,1),2,2,2,'interp_s2o1')
 !call output2d(true,interp_s2o(1:nsp,1:nop,2,2),2,2,2,'interp_s2o2')
 !call output2d(true,interp_s2o(1:nsp,1:nop,2,3),2,2,2,'interp_s2o3')
  call output2d(array_debugging,interp_s2g(1:nsp,1:4*n_qp1de,2),2,1,2,'interp_s2g')
  call output2d(array_debugging,projct_q2s(1:nqp,1:nsp,2),2,2,2,'projct_q2s')
  call output2d(array_debugging,projct_s2p(1:nsp,1:nsp,2,1),2,2,2,'projct_s2p')
 !call output2d(true,projct_s2p(1:nsp,1:nsp,2,1),2,2,2,'projct_s2p1')
 !call output2d(true,projct_s2p(1:nsp,1:nsp,2,5),2,2,2,'projct_s2p5')
  call output2d(array_debugging,dr_mat(1:nsp,1:nsp,2),2,2,2,'dr_mat')
  call output2d(array_debugging,ds_mat(1:nsp,1:nsp,2),2,2,2,'ds_mat')
  call output2d(array_debugging,dqr_mat(1:nqp,1:nqp,2),2,2,2,'dqr_mat')
  call output2d(array_debugging,dqs_mat(1:nqp,1:nqp,2),2,2,2,'dqs_mat')
  call output2d(array_debugging,dqr_LPmat(1:nqp,1:nsp,2),2,2,2,'dqr_LPmat')
  call output2d(array_debugging,dqs_LPmat(1:nqp,1:nsp,2),2,2,2,'dqs_LPmat')
  call output2d(array_debugging,dqr_CRmat(1:nsp,1:nqp,2),2,2,2,'dqr_CRmat')
  call output2d(array_debugging,dqs_CRmat(1:nsp,1:nqp,2),2,2,2,'dqs_CRmat')
  call output2d(array_debugging,invM_mat(1:nsp,1:nsp,2),2,2,2,'invM_mat')
  !
 !write (121,5) "Solution point weights 1D"
 !write (121,4) (solwts_1d(i),i=1,n_sp1de)
 !write (121,5) "Quadrature weights 1D"
 !write (121,4) (qdtwts_1d(i),i=1,n_sp1de)
 !write (121,5) "Inverse solution point weights 1D"
 !write (121,4) (one/solwts_1d(i),i=1,n_sp1de)
 !write (121,5) "Inverse quadrature weights 1D"
 !write (121,4) (one/qdtwts_1d(i),i=1,n_sp1de)
 !4 format (/,*(f12.6))
 !5 format (/,a)
 !!
 !write (104,*) "Dqr matrix"
 !write (105,*) "Dqs matrix"
 !!
 !write (120,*) "weights"
 !do j = 1,n_sp1de
 !  write (120,3) (weights2d( (n_sp1de-j)*n_sp1de+i , 2 ),i=1,n_sp1de)
 !end do
 !write (120,*) "Inverse weights"
 !do j = 1,n_sp1de
 !  write (120,3) (one/weights2d( (n_sp1de-j)*n_sp1de+i , 2 ),i=1,n_sp1de)
 !end do
 !write (120,*) "Quadrature weights"
 !do j = 1,n_qp1de
 !  write (120,3) (qdtwts_2d( (n_qp1de-j)*n_qp1de+i ),i=1,n_qp1de)
 !end do
 !write (120,*) "Inverse mass matrix"
 !do l = 1,n_sp1de
 !  do k = 1,n_sp1de
 !    write (104,2) k,l
 !    write (105,2) k,l
 !    write (120,2) k,l
 !    do j = 1,n_qp1de
 !      write (104,3) (dqr_mat( (n_qp1de-j)*n_qp1de+i , (l-1)*n_sp1de+k , 2 ),i=1,n_qp1de)
 !      write (105,3) (dqs_mat( (n_qp1de-j)*n_qp1de+i , (l-1)*n_sp1de+k , 2 ),i=1,n_qp1de)
 !    end do
 !    do j = 1,n_sp1de
 !      write (120,3) (invM_mat( (n_sp1de-j)*n_sp1de+i , (l-1)*n_sp1de+k , 2 ),i=1,n_sp1de)
 !    end do
 !  end do
 !end do
 !2 format (/,"i = ",i2,", j = ",i2)
 !3 format (*(es12.4))
 !flush(120)
 !!
 !write (111,*) "Lagrange on Quadrature Points"
 !do l = 1,4
 !  do k = 1,n_sp1de
 !    write (111,4) l,k
 !    do j = 1,n_qp1de
 !      write (111,3) (interp_s2g( (n_qp1de-j)*n_qp1de+i , (l-1)*n_sp1de+k ),i=1,n_qp1de)
 !    end do
 !  end do
 !end do
 !4 format (/,"Face = ",i2,", i = ",i2)
 !write (102,*) "Quadrature Interpolation Vandermonde Matrix"
 !do lj = 0,n_order
 !  do li = 0,n_order
 !    l = n_sp1de*lj + li + 1
 !    write (102,*)
 !    write (102,11) "Basis function for term (x^",li,")*(y^",lj,")"
 !    do kj = n_qp1de,1,-1
 !      write (102,12) (vand_qdtr(n_qp1de*(kj-1)+ki,l),ki=1,n_qp1de)
 !    end do
 !  end do
 !end do
 !!
 !write (103,*) "Quadrature Interpolation Matrix"
 !do lj = 0,n_order
 !  do li = 0,n_order
 !    l = n_sp1de*lj + li + 1
 !    write (103,*)
 !    write (103,11) "Basis function for term (x^",li,")*(y^",lj,")"
 !    do kj = n_qp1de,1,-1
 !      write (103,12) (interp_s2q(n_qp1de*(kj-1)+ki,l),ki=1,n_qp1de)
 !    end do
 !  end do
 !end do
 !!
 !if (mypnum == 0) call write_lagrange_tecplot
 !!
 !10 format (100es14.6)
  11 format (a,i0,a,i0,a)
  12 format (*(es12.4))
  !
  1 format (/," ************************************************",/, &
              "     USING OVER-INTEGRATION FOR DE-ALIASING",/, &
              " ************************************************",/)
 !!
end subroutine interp2quadrature
!
!###############################################################################
!
subroutine gauss_lobatto_matrices()
  !
  !.. Use Statements ..
  use geovar, only : n_order,n_sp1de,n_sp2dt,n_sp2dq
  use geovar, only : exact2numer,numer2exact
  use geovar, only : radau2numer,numer2radau
  use geovar, only : gauss2numer,numer2gauss
  use geovar, only : lobatto2numer,numer2lobatto
  use ovar,   only : error_loc,itestcase
  !
  use qpvar, only : gauss_1d_qp,lobatto_1d_qp,radau_1d_qp
  use qpvar, only : vand_quad_qp,vand_quad_inv_qp
  use qpvar, only : exact2numer_qp,numer2exact_qp
  use qpvar, only : radau2numer_qp,numer2radau_qp
  use qpvar, only : gauss2numer_qp,numer2gauss_qp
  use qpvar, only : lobatto2numer_qp,numer2lobatto_qp
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable, dimension(:,:) :: gauss_vand_qp
  real(qp), allocatable, dimension(:,:) :: gauss_vand_inv_qp
  real(qp), allocatable, dimension(:,:) :: lobatto_vand_qp
  real(qp), allocatable, dimension(:,:) :: lobatto_vand_inv_qp
  real(qp), allocatable, dimension(:,:) :: radau_vand_qp
  real(qp), allocatable, dimension(:,:) :: radau_vand_inv_qp
  !
continue
  !
  ! Allocate the quad precision arrays
  !
  allocate ( radau_vand_qp(1:n_sp2dq,1:n_sp2dq) , source=qzero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","radau_vand_qp",1,ierr)
  allocate ( radau_vand_inv_qp(1:n_sp2dq,1:n_sp2dq) , source=qzero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","radau_vand_inv_qp",1,ierr)
  !
  allocate ( gauss_vand_qp(1:n_sp2dq,1:n_sp2dq) , source=qzero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","gauss_vand_qp",1,ierr)
  allocate ( gauss_vand_inv_qp(1:n_sp2dq,1:n_sp2dq) , source=qzero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","gauss_vand_inv_qp",1,ierr)
  !
  allocate ( lobatto_vand_qp(1:n_sp2dq,1:n_sp2dq) , source=qzero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","lobatto_vand_qp",1,ierr)
  allocate ( lobatto_vand_inv_qp(1:n_sp2dq,1:n_sp2dq) , source=qzero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","lobatto_vand_inv_qp",1,ierr)
  !
  allocate ( radau2numer_qp(1:n_sp2dq,1:n_sp2dq) , source=qzero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","radau2numer_qp",1,ierr)
  allocate ( numer2radau_qp(1:n_sp2dq,1:n_sp2dq) , source=qzero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","numer2radau_qp",1,ierr)
  !
  allocate ( gauss2numer_qp(1:n_sp2dq,1:n_sp2dq) , source=qzero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","gauss2numer_qp",1,ierr)
  allocate ( numer2gauss_qp(1:n_sp2dq,1:n_sp2dq) , source=qzero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","numer2gauss_qp",1,ierr)
  !
  allocate ( lobatto2numer_qp(1:n_sp2dq,1:n_sp2dq) , source=qzero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","lobatto2numer_qp",1,ierr)
  allocate ( numer2lobatto_qp(1:n_sp2dq,1:n_sp2dq) , source=qzero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","numer2lobatto_qp",1,ierr)
  !
  allocate ( exact2numer_qp(1:n_sp2dq,1:n_sp2dq) , source=qzero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","exact2numer_qp",1,ierr)
  allocate ( numer2exact_qp(1:n_sp2dq,1:n_sp2dq) , source=qzero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","numer2exact_qp",1,ierr)
  !
  ! Get the Vandermonde matrix and its inverse for the right Radau nodes
  !
  radau_vand_qp = Vandermonde2d_quad(radau_1d_qp)
  radau_vand_qp = chop_qp( radau_vand_qp )
  !
  radau_vand_inv_qp = invert_matrix(radau_vand_qp)
  radau_vand_inv_qp = chop_qp( radau_vand_inv_qp )
  !
  ! Get the Vandermonde matrix and its inverse for the Gauss nodes
  !
  gauss_vand_qp = Vandermonde2d_quad(gauss_1d_qp)
  gauss_vand_qp = chop_qp( gauss_vand_qp )
  !
  gauss_vand_inv_qp = invert_matrix(gauss_vand_qp)
  gauss_vand_inv_qp = chop_qp( gauss_vand_inv_qp )
  !
  ! Get the Vandermonde matrix and its inverse for the Lobatto nodes
  !
  lobatto_vand_qp = Vandermonde2d_quad(lobatto_1d_qp)
  lobatto_vand_qp = chop_qp( lobatto_vand_qp )
  !
  lobatto_vand_inv_qp = invert_matrix(lobatto_vand_qp)
  lobatto_vand_inv_qp = chop_qp( lobatto_vand_inv_qp )
  !
  ! Now get the conversion matrices through matrix multiplication
  !
 !if (any(itestcase == Advec_Problems)) then
    !
    radau2numer_qp = matmul( vand_quad_qp , radau_vand_inv_qp )
    radau2numer_qp = chop_qp( radau2numer_qp )
    !
    numer2radau_qp = matmul( radau_vand_qp , vand_quad_inv_qp )
    numer2radau_qp = chop_qp( numer2radau_qp )
    !
    gauss2numer_qp = matmul( vand_quad_qp , gauss_vand_inv_qp )
    gauss2numer_qp = chop_qp( gauss2numer_qp )
    !
    numer2gauss_qp = matmul( gauss_vand_qp , vand_quad_inv_qp )
    numer2gauss_qp = chop_qp( numer2gauss_qp )
    !
    lobatto2numer_qp = matmul( vand_quad_qp , lobatto_vand_inv_qp )
    lobatto2numer_qp = chop_qp( lobatto2numer_qp )
    !
    numer2lobatto_qp = matmul( lobatto_vand_qp , vand_quad_inv_qp )
    numer2lobatto_qp = chop_qp( numer2lobatto_qp )
    !
 !else
    !
    numer2exact_qp = merge( matmul(gauss_vand_qp,vand_quad_inv_qp), &
                            matmul(radau_vand_qp,vand_quad_inv_qp), &
                            error_loc > 0 )
    numer2exact_qp = chop_qp( numer2exact_qp )
    !
    exact2numer_qp = merge( matmul(vand_quad_qp,gauss_vand_inv_qp), &
                            matmul(vand_quad_qp,radau_vand_inv_qp), &
                            error_loc > 0 )
    exact2numer_qp = chop_qp( exact2numer_qp )
    !
 !end if
  !
  ! Deallocate the arrays that arent needed anymore
  !
  deallocate ( gauss_vand_qp , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","gauss_vand_qp",2,ierr)
  deallocate ( gauss_vand_inv_qp , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","gauss_vand_inv_qp",2,ierr)
  deallocate ( lobatto_vand_qp , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","lobatto_vand_qp",2,ierr)
  deallocate ( lobatto_vand_inv_qp , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","lobatto_vand_inv_qp",2,ierr)
  deallocate ( radau_vand_qp , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","radau_vand_qp",2,ierr)
  deallocate ( radau_vand_inv_qp , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","radau_vand_inv_qp",2,ierr)
  !
  ! Allocate the working precision arrays
  !
  allocate ( radau2numer(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","radau2numer",1,ierr)
  allocate ( numer2radau(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","numer2radau",1,ierr)
  !
  allocate ( gauss2numer(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","gauss2numer",1,ierr)
  allocate ( numer2gauss(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","numer2gauss",1,ierr)
  !
  allocate ( lobatto2numer(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","lobatto2numer",1,ierr)
  allocate ( numer2lobatto(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","numer2lobatto",1,ierr)
  !
  allocate ( exact2numer(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","exact2numer",1,ierr)
  allocate ( numer2exact(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","numer2exact",1,ierr)
  !
 !if (any(itestcase == Advec_Problems)) then
    radau2numer   = chop( radau2numer_qp )
    numer2radau   = chop( numer2radau_qp )
    gauss2numer   = chop( gauss2numer_qp )
    numer2gauss   = chop( numer2gauss_qp )
    lobatto2numer = chop( lobatto2numer_qp )
    numer2lobatto = chop( numer2lobatto_qp )
 !else
    numer2exact   = chop( numer2exact_qp )
    exact2numer   = chop( exact2numer_qp )
 !end if
  call output2d(array_debugging,radau2numer(:,:),2,2,2,'radau2numer')
  call output2d(array_debugging,numer2radau(:,:),2,2,2,'numer2radau')
  call output2d(array_debugging,gauss2numer(:,:),2,2,2,'gauss2numer')
  call output2d(array_debugging,numer2gauss(:,:),2,2,2,'numer2gauss')
  call output2d(array_debugging,lobatto2numer(:,:),2,2,2,'lobatto2numer')
  call output2d(array_debugging,numer2lobatto(:,:),2,2,2,'numer2lobatto')
  call output2d(array_debugging,numer2exact(:,:),2,2,2,'numer2exact')
  call output2d(array_debugging,exact2numer(:,:),2,2,2,'exact2numer')
  !
  1 format ('VARIABLES = "X" "Y"')
  2 format ('ZONE T="Reference Quad", I=2, J=2')
  3 format ('ZONE T="Gauss Quad", I=',i0,', J=',i0)
  4 format ('ZONE T="Lobatto Quad", I=',i0,', J=',i0)
  5 format ('ZONE T="Reference Tri", N=3, E=1, DATAPACKING=POINT, ZONETYPE=FETRIANGLE')
  6 format ('ZONE T="Gauss Tri"')
  7 format ('ZONE T="Lobatto Tri"')
  8 format (2es23.15)
  9 format (3i4)
  !
end subroutine gauss_lobatto_matrices
!
!###############################################################################
!
pure function projection_matrix(basis,wmat) result(return_value)
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
  basisT = chop_qp( basisT )
  !
  pmat = matmul( basisT , basis )
  pmat = chop_qp( pmat )
  !
  pmat = invert_matrix( pmat )
  pmat = chop_qp( pmat )
  !
  qmat = matmul( basis , matmul( pmat , basisT ) )
  qmat = chop_qp( qmat )
  !
  return_value = transpose( qmat )
  return_value = chop_qp( return_value )
  !
end function projection_matrix
!
!###############################################################################
!
pure function interpolation_matrix(basis,newbasis,wmat) result(return_value)
  !
  real(qp), dimension(:,:), intent(in) :: basis
  real(qp), dimension(:,:), intent(in) :: newbasis
  real(qp), dimension(:,:), intent(in) :: wmat
  !
  real(qp), dimension(1:size(wmat,1),1:size(newbasis,1)) :: return_value
  !
  real(qp), dimension(1:size(basis,2),1:size(basis,1)) :: basisT
  real(qp), dimension(1:size(basis,2),1:size(basis,2)) :: pmat
  real(qp), dimension(1:size(basis,2),1:size(wmat,2)) :: tmat
  !
continue
  !
  basisT(:,:) = qzero
  pmat(:,:) = qzero
  tmat(:,:) = qzero
  !
  basisT = matmul( transpose(basis) , wmat )
  basisT = chop_qp( basisT )
  !
  pmat = matmul( basisT , basis )
  pmat = chop_qp( pmat )
  !
  pmat = invert_matrix( pmat )
  pmat = chop_qp( pmat )
  !
  tmat = matmul( pmat , basisT )
  tmat = chop_qp( tmat )
  !
  return_value = transpose( matmul( newbasis , tmat ) )
  return_value = chop_qp( return_value )
  !
end function interpolation_matrix
!
!###############################################################################
!
pure function CellVolume1D_DP(nk,np,x,Dr) result(return_value)
  !
  integer, intent(in) :: nk
  integer, intent(in) :: np
  !
  real(dp), dimension(1:np,1:nk), intent(in) :: x
  real(dp), dimension(1:np,1:np), intent(in) :: Dr
  !
  real(dp), dimension(1:nk) :: return_value
  !
  integer :: i
  !
  integer, parameter :: lp = dp
  !
continue
  !
  do i = 1,nk
    return_value(i) = dot_product(Dr(1,:),x(:,i))
  end do
  !
end function CellVolume1D_DP
!
!###############################################################################
!
pure function CellVolume1D_QP(nk,np,x,Dr) result(return_value)
  !
  integer, intent(in) :: nk
  integer, intent(in) :: np
  !
  real(qp), dimension(1:np,1:nk), intent(in) :: x
  real(qp), dimension(1:np,1:np), intent(in) :: Dr
  !
  real(qp), dimension(1:nk) :: return_value
  !
  integer :: i
  !
  integer, parameter :: lp = qp
  !
continue
  !
  do i = 1,nk
    return_value(i) = dot_product(Dr(1,:),x(:,i))
  end do
  !
end function CellVolume1D_QP
!
!###############################################################################
!
pure function CorrectionMatrix1D_DP(xi) result(return_value)
  !
  !.. Formal Arguments ..
  real(dp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(dp), dimension(1:size(xi)) :: return_value
  !
  !.. Local Scalars ..
  integer   :: i,j,nfp
  real(dp) :: c1,c2
  !
  !.. Local Arrays ..
  real(dp), dimension(0:size(xi))              :: cdRR
  real(dp), dimension(-1:size(xi),-1:size(xi)) :: coefLe
  !
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: half = 0.5_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  nfp = size(xi)
  !
  return_value(:) = zero
  !
  ! Initialize the coefficients of Legendre polynomial up to degree nfp
  !
  coefLe(:,:) = zero
  coefLe(0,0) = one
  coefLe(1,1) = one
  !
  ! Get the coefficients for every Legendre polynomial up to degree nfp
  !
  do i = 2,nfp
    do j = 0,i
      !
      c1 = real(2*i-1,kind=lp)/real(i,kind=lp)
      c2 = real(  i-1,kind=lp)/real(i,kind=lp)
      !
      coefLe(i,j) = c1*coefLe(i-1,j-1) - c2*coefLe(i-2,j)
      !
    end do
  end do
  !
  ! Get the coefficients for the derivative of
  ! the right Radau polynomial of degree nfp
  !
  c1 = half * real((-1)**nfp,kind=lp)
  do i = 1,nfp
    cdRR(i-1) = c1 * real(i,kind=lp) * (coefLe(nfp,i) - coefLe(nfp-1,i))
  end do
  !
  ! Finally, get the value for the derivative of the
  ! right Radau polynomial at each 1D solution point
  !
  do i = 1,nfp
    return_value(i) = -sum( (xi(i)**(/0:nfp-1/)) * cdRR(0:nfp-1) )
  end do
  !
end function CorrectionMatrix1D_DP
!
!###############################################################################
!
pure function CorrectionMatrix1D_QP(xi) result(return_value)
  !
  !.. Formal Arguments ..
  real(qp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(qp), dimension(1:size(xi)) :: return_value
  !
  !.. Local Scalars ..
  integer   :: i,j,nfp
  real(qp) :: c1,c2
  !
  !.. Local Arrays ..
  real(qp), dimension(0:size(xi))              :: cdRR
  real(qp), dimension(-1:size(xi),-1:size(xi)) :: coefLe
  !
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: half = 0.5_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  nfp = size(xi)
  !
  return_value(:) = zero
  !
  ! Initialize the coefficients of Legendre polynomial up to degree nfp
  !
  coefLe(:,:) = zero
  coefLe(0,0) = one
  coefLe(1,1) = one
  !
  ! Get the coefficients for every Legendre polynomial up to degree nfp
  !
  do i = 2,nfp
    do j = 0,i
      !
      c1 = real(2*i-1,kind=lp)/real(i,kind=lp)
      c2 = real(  i-1,kind=lp)/real(i,kind=lp)
      !
      coefLe(i,j) = c1*coefLe(i-1,j-1) - c2*coefLe(i-2,j)
      !
    end do
  end do
  !
  ! Get the coefficients for the derivative of
  ! the right Radau polynomial of degree nfp
  !
  c1 = half * real((-1)**nfp,kind=lp)
  do i = 1,nfp
    cdRR(i-1) = c1 * real(i,kind=lp) * (coefLe(nfp,i) - coefLe(nfp-1,i))
  end do
  !
  ! Finally, get the value for the derivative of the
  ! right Radau polynomial at each 1D solution point
  !
  do i = 1,nfp
    return_value(i) = -sum( (xi(i)**(/0:nfp-1/)) * cdRR(0:nfp-1) )
  end do
  !
end function CorrectionMatrix1D_QP
!
!###############################################################################
!
pure function CorrectionMatrix2D_quad_dg_DP(pts1d,wts1d) result(return_value)
  !
  !.. Formal Arguments ..
  real(dp), dimension(:), intent(in) :: pts1d
  real(dp), dimension(:), intent(in) :: wts1d
  !
  !.. Function Return Value ..
  real(dp), dimension(1:size(pts1d)**2,1:4*size(pts1d)) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,k,nf,np
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  return_value(:,:) = zero
  !
  np = size(pts1d)
  !
  ! Face 1 - bottom/south edge of quad
  !
  nf = 1
  do i = 1,np
    do j = 1,np
      k = (j-1)*np + i
      return_value(k,nf) = eval_LagrangePoly1D(j,-one,pts1d) / wts1d(j)
    end do
    nf = nf + 1
  end do
  !
  ! Face 2 - right/east edge of quad
  !
  nf = np+1
  do j = 1,np
    do i = np,1,-1
      k = (j-1)*np + i
      return_value(k,nf) = eval_LagrangePoly1D(i,one,pts1d) / wts1d(i)
    end do
    nf = nf + 1
  end do
  !
  ! Face 3 - top/north edge of quad
  !
  nf = 2*np+1
  do i = np,1,-1
    do j = np,1,-1
      k = (j-1)*np + i
      return_value(k,nf) = eval_LagrangePoly1D(j,one,pts1d) / wts1d(j)
    end do
    nf = nf + 1
  end do
  !
  ! Face 4 - left/west edge of quad
  !
  nf = 3*np+1
  do j = np,1,-1
    do i = 1,np
      k = (j-1)*np + i
      return_value(k,nf) = eval_LagrangePoly1D(i,-one,pts1d) / wts1d(i)
    end do
    nf = nf + 1
  end do
  !
end function CorrectionMatrix2D_quad_dg_DP
!
!###############################################################################
!
pure function CorrectionMatrix2D_quad_dg_QP(pts1d,wts1d) result(return_value)
  !
  !.. Formal Arguments ..
  real(qp), dimension(:), intent(in) :: pts1d
  real(qp), dimension(:), intent(in) :: wts1d
  !
  !.. Function Return Value ..
  real(qp), dimension(1:size(pts1d)**2,1:4*size(pts1d)) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,k,nf,np
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  return_value(:,:) = zero
  !
  np = size(pts1d)
  !
  ! Face 1 - bottom/south edge of quad
  !
  nf = 1
  do i = 1,np
    do j = 1,np
      k = (j-1)*np + i
      return_value(k,nf) = eval_LagrangePoly1D(j,-one,pts1d) / wts1d(j)
    end do
    nf = nf + 1
  end do
  !
  ! Face 2 - right/east edge of quad
  !
  nf = np+1
  do j = 1,np
    do i = np,1,-1
      k = (j-1)*np + i
      return_value(k,nf) = eval_LagrangePoly1D(i,one,pts1d) / wts1d(i)
    end do
    nf = nf + 1
  end do
  !
  ! Face 3 - top/north edge of quad
  !
  nf = 2*np+1
  do i = np,1,-1
    do j = np,1,-1
      k = (j-1)*np + i
      return_value(k,nf) = eval_LagrangePoly1D(j,one,pts1d) / wts1d(j)
    end do
    nf = nf + 1
  end do
  !
  ! Face 4 - left/west edge of quad
  !
  nf = 3*np+1
  do j = np,1,-1
    do i = 1,np
      k = (j-1)*np + i
      return_value(k,nf) = eval_LagrangePoly1D(i,-one,pts1d) / wts1d(i)
    end do
    nf = nf + 1
  end do
  !
end function CorrectionMatrix2D_quad_dg_QP
!
!###############################################################################
!
pure function CorrectionMatrix2D_quad_DP(xi) result(return_value)
  !
  !.. Formal Arguments ..
  real(dp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(dp), dimension(1:size(xi)**2,1:4*size(xi)) :: return_value
  !
  !.. Local Scalars ..
  integer   :: i,j,nfp
  !
  !.. Local Arrays ..
  integer,   dimension(1:size(xi),1:size(xi)) :: Fmask
  real(dp), dimension(1:size(xi)) :: V1D
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  nfp = size(xi)
  !
  ! Initialize the correction matrix to zero
  !
  return_value(:,:) = zero
  !
  ! 1D Correction matrix
  !
  V1D = CorrectionMatrix1D(xi)
  !
  ! Face 1 - bottom/south edge of quad
  !
  do j = 1,nfp
    do i = 1,nfp
      Fmask(i,j) = j + (i-1)*nfp
    end do
  end do
  !
  do j = 1,nfp
    do i = 1,nfp
      return_value(Fmask(i,j),j) = V1D(i)
    end do
  end do
  !
  ! Face 2 - right/east edge of quad
  !
  do j = 1,nfp
    do i = 1,nfp
      Fmask(i,j) = (j-1)*nfp + nfp - i + 1
    end do
  end do
  !
  do j = 1,nfp
    do i = 1,nfp
      return_value(Fmask(i,j),j+nfp) = V1D(i)
    end do
  end do
  !
  ! Face 3 - top/north edge of quad
  !
  do j = 1,nfp
    do i = 1,nfp
      Fmask(i,j) = (nfp-i)*nfp + nfp - j + 1
    end do
  end do
  !
  do j = 1,nfp
    do i = 1,nfp
      return_value(Fmask(i,j),j+2*nfp) = V1D(i)
    end do
  end do
  !
  ! Face 4 - left/west edge of quad
  !
  do j = 1,nfp
    do i = 1,nfp
      Fmask(i,j) = (nfp-j)*nfp + i
    end do
  end do
  !
  do j = 1,nfp
    do i = 1,nfp
      return_value(Fmask(i,j),j+3*nfp) = V1D(i)
    end do
  end do
  !
end function CorrectionMatrix2D_quad_DP
!
!###############################################################################
!
pure function CorrectionMatrix2D_quad_QP(xi) result(return_value)
  !
  !.. Formal Arguments ..
  real(qp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(qp), dimension(1:size(xi)**2,1:4*size(xi)) :: return_value
  !
  !.. Local Scalars ..
  integer   :: i,j,nfp
  !
  !.. Local Arrays ..
  integer,   dimension(1:size(xi),1:size(xi)) :: Fmask
  real(qp), dimension(1:size(xi)) :: V1D
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  nfp = size(xi)
  !
  ! Initialize the correction matrix to zero
  !
  return_value(:,:) = zero
  !
  ! 1D Correction matrix
  !
  V1D = CorrectionMatrix1D(xi)
  !
  ! Face 1 - bottom/south edge of quad
  !
  do j = 1,nfp
    do i = 1,nfp
      Fmask(i,j) = j + (i-1)*nfp
    end do
  end do
  !
  do j = 1,nfp
    do i = 1,nfp
      return_value(Fmask(i,j),j) = V1D(i)
    end do
  end do
  !
  ! Face 2 - right/east edge of quad
  !
  do j = 1,nfp
    do i = 1,nfp
      Fmask(i,j) = (j-1)*nfp + nfp - i + 1
    end do
  end do
  !
  do j = 1,nfp
    do i = 1,nfp
      return_value(Fmask(i,j),j+nfp) = V1D(i)
    end do
  end do
  !
  ! Face 3 - top/north edge of quad
  !
  do j = 1,nfp
    do i = 1,nfp
      Fmask(i,j) = (nfp-i)*nfp + nfp - j + 1
    end do
  end do
  !
  do j = 1,nfp
    do i = 1,nfp
      return_value(Fmask(i,j),j+2*nfp) = V1D(i)
    end do
  end do
  !
  ! Face 4 - left/west edge of quad
  !
  do j = 1,nfp
    do i = 1,nfp
      Fmask(i,j) = (nfp-j)*nfp + i
    end do
  end do
  !
  do j = 1,nfp
    do i = 1,nfp
      return_value(Fmask(i,j),j+3*nfp) = V1D(i)
    end do
  end do
  !
end function CorrectionMatrix2D_quad_QP
!
!###############################################################################
!
pure function CorrectionMatrix2D_tria_DP(rs) result(return_value)
  !
  ! RS : r and s coordinates for the standard
  !      triangle element on the interval [-1,1]
  !
  real(dp), dimension(:,:), intent(in) :: rs
  !
  real(dp), dimension(1:size(rs,2),1:3*np2nfp(size(rs,2))) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,n,np,nfp
  !
  !.. Local Arrays ..
  integer,  dimension(1:np2nfp(size(rs,2)),1:3) :: Fmask
  real(dp), dimension(1:np2nfp(size(rs,2)),1:np2nfp(size(rs,2))) :: V1D
  real(dp), dimension(1:size(rs,2),1:size(rs,2)) :: Vtri
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: two  = 2.0_lp
  !
continue
  !
  return_value = zero
  !
  np  = size(rs,2)
  n   = np2n(np)
  nfp = np2nfp(np)
  !
  ! Compute the 2D Vandermonde matrix for the standard triangle element
  ! NOTE: r and s should be already be adjusted to the interval [-1,1] which is
  !       what the function Vandermonde2D_tria is expecting meaning we dont
  !       want to adjust them again.
  !
  Vtri(:,:) = Vandermonde2D_tria( rs(:,:) )
  !
  ! HTs method for getting Fmask
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
  ! DG Book's method for getting Fmask
  !
 !Fmask(:,1) = pack( (/(i,i=1,np)/) , abs(rs(2,:) + one    ) < eps12 )
 !Fmask(:,2) = pack( (/(i,i=1,np)/) , abs(rs(1,:) + rs(2,:)) < eps12 )
 !Fmask(:,3) = pack( (/(i,i=1,np)/) , abs(rs(1,:) + one    ) < eps12 )
  !
  ! Face 1
  !
  V1D = Vandermonde1D(rs(1,Fmask(:,1)))
  V1D = matmul( V1D , transpose(V1D) )
  return_value(Fmask(:,1),1:nfp) = invert_matrix(V1D)
  !
  ! Face 2
  !
  V1D = Vandermonde1D(rs(1,Fmask(:,2)))
  V1D = matmul( V1D , transpose(V1D) )
  return_value(Fmask(:,2),nfp+1:2*nfp) = invert_matrix(V1D)
  !
  ! Face 3
  !
  V1D = Vandermonde1D(rs(2,Fmask(:,3)))
  V1D = matmul( V1D , transpose(V1D) )
  return_value(Fmask(:,3),2*nfp+1:3*nfp) = invert_matrix(V1D)
  !
  ! Finally, compute the correction matrix
  !
  return_value = matmul( transpose(Vtri) , return_value )
  return_value = two*matmul( Vtri , return_value )
  !
end function CorrectionMatrix2D_tria_DP
!
!###############################################################################
!
pure function CorrectionMatrix2D_tria_QP(rs) result(return_value)
  !
  ! RS : r and s coordinates for the standard
  !      triangle element on the interval [-1,1]
  !
  real(qp), dimension(:,:), intent(in) :: rs
  !
  real(qp), dimension(1:size(rs,2),1:3*np2nfp(size(rs,2))) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,n,np,nfp
  !
  !.. Local Arrays ..
  integer,  dimension(1:np2nfp(size(rs,2)),1:3) :: Fmask
  real(qp), dimension(1:np2nfp(size(rs,2)),1:np2nfp(size(rs,2))) :: V1D
  real(qp), dimension(1:size(rs,2),1:size(rs,2)) :: Vtri
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: two  = 2.0_lp
  !
continue
  !
  return_value = zero
  !
  np  = size(rs,2)
  n   = np2n(np)
  nfp = np2nfp(np)
  !
  ! Compute the 2D Vandermonde matrix for the standard triangle element
  ! NOTE: r and s should be already be adjusted to the interval [-1,1] which is
  !       what the function Vandermonde2D_tria is expecting meaning we dont
  !       want to adjust them again.
  !
  Vtri(:,:) = Vandermonde2D_tria( rs(:,:) )
  !
  ! HTs method for getting Fmask
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
  ! DG Book's method for getting Fmask
  !
 !Fmask(:,1) = pack( (/(i,i=1,np)/) , abs(rs(2,:) + one    ) < eps12 )
 !Fmask(:,2) = pack( (/(i,i=1,np)/) , abs(rs(1,:) + rs(2,:)) < eps12 )
 !Fmask(:,3) = pack( (/(i,i=1,np)/) , abs(rs(1,:) + one    ) < eps12 )
  !
  ! Face 1
  !
  V1D = Vandermonde1D(rs(1,Fmask(:,1)))
  V1D = matmul( V1D , transpose(V1D) )
  return_value(Fmask(:,1),1:nfp) = invert_matrix(V1D)
  !
  ! Face 2
  !
  V1D = Vandermonde1D(rs(1,Fmask(:,2)))
  V1D = matmul( V1D , transpose(V1D) )
  return_value(Fmask(:,2),nfp+1:2*nfp) = invert_matrix(V1D)
  !
  ! Face 3
  !
  V1D = Vandermonde1D(rs(2,Fmask(:,3)))
  V1D = matmul( V1D , transpose(V1D) )
  return_value(Fmask(:,3),2*nfp+1:3*nfp) = invert_matrix(V1D)
  !
  ! Finally, compute the correction matrix
  !
  return_value = matmul( transpose(Vtri) , return_value )
  return_value = two*matmul( Vtri , return_value )
  !
end function CorrectionMatrix2D_tria_QP
!
!###############################################################################
!
pure function d2lagrange_DP(m,x,xi) result(return_value)
  !
  ! Get value of second derivative of
  ! m-th Lagrange polynomial at point x.
  !
  !                   nfp                nfp                nfp
  !  d2               ___                ___               _____
  ! --- [L_m (x)]  =  \         1        \         1        | |   x   - x(j)
  ! dx2               /__  -----------   /__  -----------   | |  -----------
  !                        x(m) - x(l)        x(m) - x(n)   j=1  x(m) - x(j)
  !                   l=1                n=1               j/=l
  !                  l/=m               n/=l               j/=n
  !                                     n/=m               j/=m
  !
  ! m  : index where the value at xi(m) is 1
  ! X  : x location at which to evaluate the second
  !      derivative of the m-th Lagrange polynomial 
  ! XI : array of all the interpolation points
  !
  !.. Formal Arguments ..
  integer,                 intent(in) :: m
  real(dp),               intent(in) :: x
  real(dp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(dp) :: return_value
  !
  !.. Local Scalars ..
  integer :: l,n,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  nfp = size(xi)
  !
  imask(:) = (/ (l,l=1,nfp) /)
  !
  return_value = zero
  !
  if (nfp > 2) then
    !
    do l = 1,nfp
      if (l == m) cycle
      do n = 1,nfp
        if ( (n == m) .or. (n == l) ) cycle
        return_value = return_value + &
                       product( x    -xi(1:nfp) , &
                                 (imask/=m .and. imask/=l .and. imask/=n) ) / &
                       product( xi(m)-xi(1:nfp) , &
                                 (imask/=m .and. imask/=l .and. imask/=n) ) / &
                       (xi(m) - xi(n))
      end do
      return_value = return_value / (xi(m) - xi(l))
    end do
    !
  end if
  !
end function d2lagrange_DP
!
!###############################################################################
!
pure function d2lagrange_QP(m,x,xi) result(return_value)
  !
  ! Get value of second derivative of
  ! m-th Lagrange polynomial at point x.
  !
  !                   nfp                nfp                nfp
  !  d2               ___                ___               _____
  ! --- [L_m (x)]  =  \         1        \         1        | |   x   - x(j)
  ! dx2               /__  -----------   /__  -----------   | |  -----------
  !                        x(m) - x(l)        x(m) - x(n)   j=1  x(m) - x(j)
  !                   l=1                n=1               j/=l
  !                  l/=m               n/=l               j/=n
  !                                     n/=m               j/=m
  !
  ! m  : index where the value at xi(m) is 1
  ! X  : x location at which to evaluate the second
  !      derivative of the m-th Lagrange polynomial 
  ! XI : array of all the interpolation points
  !
  !.. Formal Arguments ..
  integer,                 intent(in) :: m
  real(qp),               intent(in) :: x
  real(qp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(qp) :: return_value
  !
  !.. Local Scalars ..
  integer :: l,n,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  nfp = size(xi)
  !
  imask(:) = (/ (l,l=1,nfp) /)
  !
  return_value = zero
  !
  if (nfp > 2) then
    !
    do l = 1,nfp
      if (l == m) cycle
      do n = 1,nfp
        if ( (n == m) .or. (n == l) ) cycle
        return_value = return_value + &
                       product( x    -xi(1:nfp) , &
                                 (imask/=m .and. imask/=l .and. imask/=n) ) / &
                       product( xi(m)-xi(1:nfp) , &
                                 (imask/=m .and. imask/=l .and. imask/=n) ) / &
                       (xi(m) - xi(n))
      end do
      return_value = return_value / (xi(m) - xi(l))
    end do
    !
  end if
  !
end function d2lagrange_QP
!
!###############################################################################
!
pure function dlagrange_DP(m,x,xi) result(return_value)
  !
  ! Get the value of the derivative of the
  ! m-th Lagrange polynomial at the point x.
  !
  !                  nfp                nfp
  !   d              ___               _____
  !  -- [L_m (x)] =  \         1        | |   x   - x(j)
  !  dx              /__  -----------   | |  -----------
  !                       x(m) - x(l)   j=1  x(m) - x(j)
  !                  l=1               j/=l
  !                 l/=m               j/=m
  !
  ! m  : index where the value at xi(m) is 1
  ! X  : x location at which to evaluate the derivative
  !      of the m-th Lagrange polynomial 
  ! XI : array of all the interpolation points
  !
  !.. Formal Arguments ..
  integer,                 intent(in) :: m
  real(dp),               intent(in) :: x
  real(dp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(dp) :: return_value
  !
  !.. Local Scalars ..
  integer :: l,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  nfp = size(xi)
  !
  imask(:) = (/ (l,l=1,nfp) /)
  !
  return_value = zero
  !
  if (nfp > 1) then
    !
    do l = 1,nfp
      if (l == m) cycle
      return_value = return_value + &
                     product( x    -xi(1:nfp) , (imask/=m .and. imask/=l) ) / &
                     product( xi(m)-xi(1:nfp) , (imask/=m .and. imask/=l) ) / &
                     (xi(m) - xi(l))
    end do
    !
  end if
  !
end function dlagrange_DP
!
!###############################################################################
!
pure function dlagrange_QP(m,x,xi) result(return_value)
  !
  ! Get the value of the derivative of the
  ! m-th Lagrange polynomial at the point x.
  !
  !                  nfp                nfp
  !   d              ___               _____
  !  -- [L_m (x)] =  \         1        | |   x   - x(j)
  !  dx              /__  -----------   | |  -----------
  !                       x(m) - x(l)   j=1  x(m) - x(j)
  !                  l=1               j/=l
  !                 l/=m               j/=m
  !
  ! m  : index where the value at xi(m) is 1
  ! X  : x location at which to evaluate the derivative
  !      of the m-th Lagrange polynomial 
  ! XI : array of all the interpolation points
  !
  !.. Formal Arguments ..
  integer,                 intent(in) :: m
  real(qp),               intent(in) :: x
  real(qp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(qp) :: return_value
  !
  !.. Local Scalars ..
  integer :: l,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  nfp = size(xi)
  !
  imask(:) = (/ (l,l=1,nfp) /)
  !
  return_value = zero
  !
  if (nfp > 1) then
    !
    do l = 1,nfp
      if (l == m) cycle
      return_value = return_value + &
                     product( x    -xi(1:nfp) , (imask/=m .and. imask/=l) ) / &
                     product( xi(m)-xi(1:nfp) , (imask/=m .and. imask/=l) ) / &
                     (xi(m) - xi(l))
    end do
    !
  end if
  !
end function dlagrange_QP
!
!###############################################################################
!
pure function Dmatrix1D_DP(V,Vr) result(return_value)
  !
  real(dp), dimension(:,:), intent(in) :: V
  real(dp), dimension(:,:), intent(in) :: Vr
  !
  real(dp), dimension(1:size(V,1),1:size(V,2)) :: return_value
  !
  integer, parameter :: lp = dp
  !
continue
  !
  return_value = invert_matrix(V)
  return_value = matmul(Vr,return_value)
  !
end function Dmatrix1D_DP
!
!###############################################################################
!
pure function Dmatrix1D_QP(V,Vr) result(return_value)
  !
  real(qp), dimension(:,:), intent(in) :: V
  real(qp), dimension(:,:), intent(in) :: Vr
  !
  real(qp), dimension(1:size(V,1),1:size(V,2)) :: return_value
  !
  integer, parameter :: lp = qp
  !
continue
  !
  return_value = invert_matrix(V)
  return_value = matmul(Vr,return_value)
  !
end function Dmatrix1D_QP
!
!###############################################################################
!
subroutine Dmatrix2D_quad_DP(d1_mat,dr_mat,ds_mat)
  !
  !.. Formal Arguments ..
  real(dp), dimension(:,:),    intent(in) :: d1_mat
  real(dp), dimension(:,:), intent(inout) :: dr_mat
  real(dp), dimension(:,:), intent(inout) :: ds_mat
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,n
  integer :: lbeg,lend
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  n = size(d1_mat,1)
  !
  ! Initialize the quad parts of the derivative matrices to zero
  !
  dr_mat(:,:) = zero
  ds_mat(:,:) = zero
  !
  ! First the r/xi direction
  !
  k = 0
  do j = 1,n
    lbeg = (j-1)*n + 1
    lend = j*n
    do i = 1,n
      k = k + 1
      dr_mat(lbeg:lend,k) = d1_mat(i,:)
    end do
  end do
  !
  ! Next the s/eta direction
  !
  k = 0
  do j = 1,n
    do i = 1,n
      k = k + 1
      do l = 1,n
        ds_mat(i+(l-1)*n,k) = d1_mat(j,l)
      end do
    end do
  end do
  !
end subroutine Dmatrix2D_quad_DP
!
!###############################################################################
!
subroutine Dmatrix2D_quad_QP(d1_mat,dr_mat,ds_mat)
  !
  !.. Formal Arguments ..
  real(qp), dimension(:,:),    intent(in) :: d1_mat
  real(qp), dimension(:,:), intent(inout) :: dr_mat
  real(qp), dimension(:,:), intent(inout) :: ds_mat
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,n
  integer :: lbeg,lend
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  n = size(d1_mat,1)
  !
  ! Initialize the quad parts of the derivative matrices to zero
  !
  dr_mat(:,:) = zero
  ds_mat(:,:) = zero
  !
  ! First the r/xi direction
  !
  k = 0
  do j = 1,n
    lbeg = (j-1)*n + 1
    lend = j*n
    do i = 1,n
      k = k + 1
      dr_mat(lbeg:lend,k) = d1_mat(i,:)
    end do
  end do
  !
  ! Next the s/eta direction
  !
  k = 0
  do j = 1,n
    do i = 1,n
      k = k + 1
      do l = 1,n
        ds_mat(i+(l-1)*n,k) = d1_mat(j,l)
      end do
    end do
  end do
  !
end subroutine Dmatrix2D_quad_QP
!
!###############################################################################
!
subroutine Dmatrix2D_tria_DP(grad_vand,vand_inv,dr_mat,ds_mat)
  !
  !.. Formal Arguments ..
  real(dp), dimension(:,:),    intent(in) :: vand_inv
  real(dp), dimension(:,:,:),  intent(in) :: grad_vand
  real(dp), dimension(:,:), intent(inout) :: dr_mat
  real(dp), dimension(:,:), intent(inout) :: ds_mat
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: two  = 2.0_lp
  !
continue
  !
  ! Initialize the triangle parts of the derivative matrices to zero
  !
  dr_mat(:,:) = zero
  ds_mat(:,:) = zero
  !
  ! Compute the derivative matrix
  !
  dr_mat(:,:) = two * matmul( grad_vand(:,:,1) , vand_inv(:,:) )
  ds_mat(:,:) = two * matmul( grad_vand(:,:,2) , vand_inv(:,:) )
  !
  ! Transpose these matrices to get them in the right ordering
  !
  dr_mat(:,:) = transpose( dr_mat(:,:) )
  ds_mat(:,:) = transpose( ds_mat(:,:) )
  !
end subroutine Dmatrix2D_tria_DP
!
!###############################################################################
!
subroutine Dmatrix2D_tria_QP(grad_vand,vand_inv,dr_mat,ds_mat)
  !
  !.. Formal Arguments ..
  real(qp), dimension(:,:),    intent(in) :: vand_inv
  real(qp), dimension(:,:,:),  intent(in) :: grad_vand
  real(qp), dimension(:,:), intent(inout) :: dr_mat
  real(qp), dimension(:,:), intent(inout) :: ds_mat
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: two  = 2.0_lp
  !
continue
  !
  ! Initialize the triangle parts of the derivative matrices to zero
  !
  dr_mat(:,:) = zero
  ds_mat(:,:) = zero
  !
  ! Compute the derivative matrix
  !
  dr_mat(:,:) = two * matmul( grad_vand(:,:,1) , vand_inv(:,:) )
  ds_mat(:,:) = two * matmul( grad_vand(:,:,2) , vand_inv(:,:) )
  !
  ! Transpose these matrices to get them in the right ordering
  !
  dr_mat(:,:) = transpose( dr_mat(:,:) )
  ds_mat(:,:) = transpose( ds_mat(:,:) )
  !
end subroutine Dmatrix2D_tria_QP
!
!###############################################################################
!
pure subroutine Dmatrix2D_tria_weak_DP(rs,V,dr_hat,ds_hat)
  !
  real(dp), dimension(:,:),    intent(in) :: rs
  real(dp), dimension(:,:),    intent(in) :: V
  real(dp), dimension(:,:), intent(inout) :: dr_hat
  real(dp), dimension(:,:), intent(inout) :: ds_hat
  !
  !.. Local Arrays ..
  real(dp), dimension(1:size(V,1),1:size(V,2)) :: invVVT
  real(dp), dimension(1:size(V,1),1:size(V,2),1:size(rs,1)) :: dmat
  !
  integer, parameter :: lp = dp
  !
continue
  !
  ! Precompute the inverse of V*V_Transpose
  !
  invVVT = transpose(V)
  invVVT = matmul(V,invVVT)
  invVVT = invert_matrix(invVVT)
  !
  ! Get the derivative of the Vandermonde matrix and store in dmat
  !
  dmat = GradVandermonde2D(rs)
  !
  ! Compute the weak form of the derivative matrix
  !
  dr_hat(:,:) = transpose( dmat(:,:,1) )
  dr_hat(:,:) = matmul( V(:,:) , dr_hat(:,:) )
  dr_hat(:,:) = matmul( dr_hat(:,:) , invVVT )
  !
  ds_hat(:,:) = transpose( dmat(:,:,2) )
  ds_hat(:,:) = matmul( V(:,:) , ds_hat(:,:) )
  ds_hat(:,:) = matmul( ds_hat(:,:) , invVVT )
  !
end subroutine Dmatrix2D_tria_weak_DP
!
!###############################################################################
!
subroutine Dmatrix2D_tria_weak_QP(rs,V,dr_hat,ds_hat)
  !
  real(qp), dimension(:,:),    intent(in) :: rs
  real(qp), dimension(:,:),    intent(in) :: V
  real(qp), dimension(:,:), intent(inout) :: dr_hat
  real(qp), dimension(:,:), intent(inout) :: ds_hat
  !
  !.. Local Arrays ..
  real(qp), dimension(1:size(V,1),1:size(V,2)) :: invVVT
  real(qp), dimension(1:size(V,1),1:size(V,2),1:size(rs,1)) :: dmat
  !
  integer, parameter :: lp = qp
  !
continue
  !
  ! Precompute the inverse of V*V_Transpose
  !
  invVVT = transpose(V)
  invVVT = matmul(V,invVVT)
  invVVT = invert_matrix(invVVT)
  !
  ! Get the derivative of the Vandermonde matrix and store in dmat
  !
  dmat = GradVandermonde2D(rs)
  !
  ! Compute the weak form of the derivative matrix
  !
  dr_hat(:,:) = transpose( dmat(:,:,1) )
  dr_hat(:,:) = matmul( V(:,:) , dr_hat(:,:) )
  dr_hat(:,:) = matmul( dr_hat(:,:) , invVVT )
  !
  ds_hat(:,:) = transpose( dmat(:,:,2) )
  ds_hat(:,:) = matmul( V(:,:) , ds_hat(:,:) )
  ds_hat(:,:) = matmul( ds_hat(:,:) , invVVT )
  !
end subroutine Dmatrix2D_tria_weak_QP
!
!###############################################################################
!
subroutine Dmatrix2d_quad_weak_DP(wts1d,d1_mat,dr_hat,ds_hat)
  !
  ! Create the weak form of the derivative matrices
  !
  !.. Formal Arguments ..
  real(dp), dimension(:),      intent(in) :: wts1d
  real(dp), dimension(:,:),    intent(in) :: d1_mat
  real(dp), dimension(:,:), intent(inout) :: dr_hat
  real(dp), dimension(:,:), intent(inout) :: ds_hat
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,np,lbeg,lend
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  !
  !.. Local Arrays ..
  real(lp), dimension(1:size(wts1d),1:size(wts1d)) :: d1_hat
  !
continue
  !
  np = size(wts1d)
  !
  ! Initialize the weak form of the derivative matrices to zero
  !
  dr_hat(:,:) = zero
  ds_hat(:,:) = zero
  !
  ! Create the weak form of the 1D derivative matrix
  !
  do j = 1,np
    do i = 1,np
      d1_hat(i,j) = -d1_mat(j,i) * wts1d(j) / wts1d(i)
    end do
  end do
  !
  ! Create the weak form of the quadrilateral r/xi derivative matrix
  !
  k = 0
  do j = 1,np
    lbeg = (j-1)*np + 1
    lend = j*np
    do i = 1,np
      k = k + 1
      dr_hat(lbeg:lend,k) = d1_hat(i,:)
    end do
  end do
  !
  ! Create the weak form of the quadrilateral s/eta derivative matrix
  !
  k = 0
  do j = 1,np
    do i = 1,np
      k = k + 1
      do l = 1,np
        ds_hat(i+(l-1)*np,k) = d1_hat(j,l)
      end do
    end do
  end do
  !
end subroutine Dmatrix2d_quad_weak_DP
!
!###############################################################################
!
subroutine Dmatrix2d_quad_weak_QP(wts1d,d1_mat,dr_hat,ds_hat)
  !
  ! Create the weak form of the derivative matrices
  !
  !.. Formal Arguments ..
  real(qp), dimension(:),      intent(in) :: wts1d
  real(qp), dimension(:,:),    intent(in) :: d1_mat
  real(qp), dimension(:,:), intent(inout) :: dr_hat
  real(qp), dimension(:,:), intent(inout) :: ds_hat
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,np,lbeg,lend
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: zero = 0.0_lp
  !
  !.. Local Arrays ..
  real(lp), dimension(1:size(wts1d),1:size(wts1d)) :: d1_hat
  !
continue
  !
  np = size(wts1d)
  !
  ! Initialize the weak form of the derivative matrices to zero
  !
  dr_hat(:,:) = zero
  ds_hat(:,:) = zero
  !
  ! Create the weak form of the 1D derivative matrix
  !
  do j = 1,np
    do i = 1,np
      d1_hat(i,j) = -d1_mat(j,i) * wts1d(j) / wts1d(i)
    end do
  end do
  !
  ! Create the weak form of the quadrilateral r/xi derivative matrix
  !
  k = 0
  do j = 1,np
    lbeg = (j-1)*np + 1
    lend = j*np
    do i = 1,np
      k = k + 1
      dr_hat(lbeg:lend,k) = d1_hat(i,:)
    end do
  end do
  !
  ! Create the weak form of the quadrilateral s/eta derivative matrix
  !
  k = 0
  do j = 1,np
    do i = 1,np
      k = k + 1
      do l = 1,np
        ds_hat(i+(l-1)*np,k) = d1_hat(j,l)
      end do
    end do
  end do
  !
end subroutine Dmatrix2d_quad_weak_QP
!
!###############################################################################
!
pure function eval_LagrangePoly1D_DP(k,x,xi) result(return_value)
  !
  ! Get the value of the 1D Lagrange polynomial at point x
  !
  ! K  : The value of the Lagrange polynomial is 1 at xi(k)
  ! X  : Location at which to evaluate the k-th Lagrange polynomial
  ! XI : Array of all the interpolation points
  !
  ! This evaluates the Lagrange polynomial at the location x.
  ! The Lagrange polynomial corresponds to a function that
  ! is zero at all iterpolation points xi(i) for i=1,nfp
  ! except at xi(k) where the function is equal to 1.
  !
  !.. Formal Arguments ..
  integer,                 intent(in) :: k
  real(dp),               intent(in) :: x
  real(dp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(dp) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
  !
  integer, parameter :: lp = dp
  !
continue
  !
  nfp = size(xi)
  imask(:) = (/ (i,i=1,nfp) /)
  !
  return_value = product( x    -xi(1:nfp) , imask /= k ) / &
                 product( xi(k)-xi(1:nfp) , imask /= k )
  !
end function eval_LagrangePoly1D_DP
!
!###############################################################################
!
pure function eval_LagrangePoly1D_QP(k,x,xi) result(return_value)
  !
  ! Get the value of the 1D Lagrange polynomial at point x
  !
  ! K  : The value of the Lagrange polynomial is 1 at xi(k)
  ! X  : Location at which to evaluate the k-th Lagrange polynomial
  ! XI : Array of all the interpolation points
  !
  ! This evaluates the Lagrange polynomial at the location x.
  ! The Lagrange polynomial corresponds to a function that
  ! is zero at all iterpolation points xi(i) for i=1,nfp
  ! except at xi(k) where the function is equal to 1.
  !
  !.. Formal Arguments ..
  integer,                 intent(in) :: k
  real(qp),               intent(in) :: x
  real(qp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(qp) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
  !
  integer, parameter :: lp = qp
  !
continue
  !
  nfp = size(xi)
  imask(:) = (/ (i,i=1,nfp) /)
  !
  return_value = product( x    -xi(1:nfp) , imask /= k ) / &
                 product( xi(k)-xi(1:nfp) , imask /= k )
  !
end function eval_LagrangePoly1D_QP
!
!###############################################################################
!
pure function eval_LagrangePoly2D_DP(k,l,x,y,xi,eta) result(return_value)
  !
  ! Get value of the 2D Lagrange polynomial at point x,y
  !
  ! K   : The value of the Lagrange polynomial in coordinate XI  is 1 at xi(k)
  ! L   : The value of the Lagrange polynomial in coordinate ETA is 1 at eta(l)
  ! X   : Location at which to evaluate the k-th Lagrange polynomial in coordinate XI
  ! Y   : Location at which to evaluate the l-th Lagrange polynomial in coordinate ETA
  ! XI  : Array of all the interpolation points in XI  coordinate
  ! ETA : Array of all the interpolation points in ETA coordinate
  !
  ! This evaluates the 2D Lagrange polynomial at the location x,y.
  ! The 2D Lagrange polynomial corresponds to a function that
  ! is zero at all iterpolation points [xi(i),eta(j)] for [i=1,ifp;j=1,jfp]
  ! except at [xi(k),eta(l)] where the function is equal to 1.
  !
  !.. Formal Arguments ..
  integer,                 intent(in) :: k
  integer,                 intent(in) :: l
  real(dp),               intent(in) :: x
  real(dp),               intent(in) :: y
  real(dp), dimension(:), intent(in) :: xi
  real(dp), dimension(:), intent(in) :: eta
  !
  !.. Function Return Value ..
  real(dp) :: return_value
  !
  integer, parameter :: lp = dp
  !
continue
  !
  return_value = eval_LagrangePoly1D(k,x,xi) * eval_LagrangePoly1D(l,y,eta)
  !
end function eval_LagrangePoly2D_DP
!
!###############################################################################
!
pure function eval_LagrangePoly2D_QP(k,l,x,y,xi,eta) result(return_value)
  !
  ! Get value of the 2D Lagrange polynomial at point x,y
  !
  ! K   : The value of the Lagrange polynomial in coordinate XI  is 1 at xi(k)
  ! L   : The value of the Lagrange polynomial in coordinate ETA is 1 at eta(l)
  ! X   : Location at which to evaluate the k-th Lagrange polynomial in coordinate XI
  ! Y   : Location at which to evaluate the l-th Lagrange polynomial in coordinate ETA
  ! XI  : Array of all the interpolation points in XI  coordinate
  ! ETA : Array of all the interpolation points in ETA coordinate
  !
  ! This evaluates the 2D Lagrange polynomial at the location x,y.
  ! The 2D Lagrange polynomial corresponds to a function that
  ! is zero at all iterpolation points [xi(i),eta(j)] for [i=1,ifp;j=1,jfp]
  ! except at [xi(k),eta(l)] where the function is equal to 1.
  !
  !.. Formal Arguments ..
  integer,                 intent(in) :: k
  integer,                 intent(in) :: l
  real(qp),               intent(in) :: x
  real(qp),               intent(in) :: y
  real(qp), dimension(:), intent(in) :: xi
  real(qp), dimension(:), intent(in) :: eta
  !
  !.. Function Return Value ..
  real(qp) :: return_value
  !
  integer, parameter :: lp = qp
  !
continue
  !
  return_value = eval_LagrangePoly1D(k,x,xi) * eval_LagrangePoly1D(l,y,eta)
  !
end function eval_LagrangePoly2D_QP
!
!###############################################################################
!
pure function Filter2D_tria_DP(vand,vand_inv) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Use Statements ..
  use geovar, only : n_order,n_sp2dt
  use ovar,   only : Filter_Order, Filter_Cutoff,Filter_Alpha
  use ovar,   only : Filter_etac
  !
  !.. Formal Arguments ..
  real(dp), dimension(1:n_sp2dt,1:n_sp2dt), intent(in) :: vand
  real(dp), dimension(1:n_sp2dt,1:n_sp2dt), intent(in) :: vand_inv
  !
  !.. Function Return Value ..
  real(dp), dimension(1:n_sp2dt,1:n_sp2dt) :: return_value
  !
  !.. Local Scalars ..
  integer  :: i,j,k
  real(dp) :: rd,rn,falpha,fetac
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one = 1.0_lp
  !
continue
  !
  ! Covert Filter_Alpha to local precision
  !
  falpha = real( Filter_Alpha , kind=lp )
  fetac  = real( Filter_etac , kind=lp )
  !
  ! Initialize the local filter matrices to identity matrices
  !
  return_value(:,:) = zero
  do i = 1,n_sp2dt
    return_value(i,i) = one
  end do
  !
  ! Compute the filter diagonal for triangles
  !
  k = 1
  do i = 0,n_order
    do j = 0,n_order-i
      !
      if (i+j >= Filter_Cutoff) then
        !
        rn = real(i,kind=lp)/real(n_order,kind=lp) - fetac
        rd = one - fetac
        !
        return_value(k,k) = exp(-falpha * (rn/rd)**Filter_Order)
        !
      end if
      !
      k = k + 1
      !
    end do
  end do
  !
  ! Compute the final filter matrix for triangles
  !
  return_value = matmul( vand , return_value )
  return_value = matmul( return_value , vand_inv )
  !
end function Filter2D_tria_DP
!
!###############################################################################
!
pure function Filter2D_tria_QP(vand,vand_inv) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Use Statements ..
  use geovar, only : n_order,n_sp2dt
  use ovar,   only : Filter_Order, Filter_Cutoff,Filter_Alpha
  use ovar,   only : Filter_etac
  !
  !.. Formal Arguments ..
  real(qp), dimension(1:n_sp2dt,1:n_sp2dt), intent(in) :: vand
  real(qp), dimension(1:n_sp2dt,1:n_sp2dt), intent(in) :: vand_inv
  !
  !.. Function Return Value ..
  real(qp), dimension(1:n_sp2dt,1:n_sp2dt) :: return_value
  !
  !.. Local Scalars ..
  integer  :: i,j,k
  real(qp) :: rd,rn,falpha,fetac
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one = 1.0_lp
  !
continue
  !
  ! Covert Filter_Alpha to local precision
  !
  falpha = real( Filter_Alpha , kind=lp )
  fetac  = real( Filter_etac , kind=lp )
  !
  ! Initialize the local filter matrices to identity matrices
  !
  return_value(:,:) = zero
  do i = 1,n_sp2dt
    return_value(i,i) = one
  end do
  !
  ! Compute the filter diagonal for triangles
  !
  k = 1
  do i = 0,n_order
    do j = 0,n_order-i
      !
      if (i+j >= Filter_Cutoff) then
        !
        rn = real(i,kind=lp)/real(n_order,kind=lp) - fetac
        rd = one - fetac
        !
        return_value(k,k) = exp(-falpha * (rn/rd)**Filter_Order)
        !
      end if
      !
      k = k + 1
      !
    end do
  end do
  !
  ! Compute the final filter matrix for triangles
  !
  return_value = matmul( vand , return_value )
  return_value = matmul( return_value , vand_inv )
  !
end function Filter2D_tria_QP
!
!###############################################################################
!
pure function Filter2D_quad_DP(vand,vand_inv) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Use Statements ..
  use geovar, only : n_order,n_sp1de,n_sp2dq
  use ovar,   only : Filter_Option,Filter_Order, &
                     Filter_Cutoff,Filter_Alpha
  use ovar,   only : Filter_etac
  !
  !.. Formal Arguments ..
  real(dp), dimension(1:n_sp2dq,1:n_sp2dq), intent(in) :: vand
  real(dp), dimension(1:n_sp2dq,1:n_sp2dq), intent(in) :: vand_inv
  !
  !.. Function Return Value ..
  real(dp), dimension(1:n_sp2dq,1:n_sp2dq) :: return_value
  !
  !.. Local Scalars ..
  integer  :: i,j,k
  real(dp) :: rd,rn,falpha,fetac
  !
  !.. Local Arrays ..
  real(dp), dimension(1:n_sp1de) :: filter_1D
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one = 1.0_lp
  !
continue
  !
  ! Covert Filter_Alpha to local precision
  !
  falpha = real( Filter_Alpha , kind=lp )
  fetac  = real( Filter_etac , kind=lp )
  !
  ! Initialize the local filter matrices to identity matrices
  !
  filter_1D(:) = one
  !
  return_value(:,:) = zero
  do i = 1,n_sp2dq
    return_value(i,i) = one
  end do
  !
  ! Compute the filter diagonal for quadrilaterals
  !
  if (Filter_Option == 1) then
    !
    ! Option 1 - Tensor product of the 1D filter
    !           (Equivalent to using a 1D filter in each direction)
    !
    ! First get the 1D filter diagonal
    !
    do i = 0,n_order
      if (i >= Filter_Cutoff) then
        !
        rn = real(i,kind=lp)/real(n_order,kind=lp) - fetac
        rd = one - fetac
        !
        filter_1D(i+1) = exp(-falpha * (rn/rd)**Filter_Order)
        !
      end if
    end do
    !
    ! Now fill the filter diagonal for quads
    !
    do j = 0,n_order
      do i = 0,n_order
        !
        k = n_sp1de*j + i + 1
        !
        return_value(k,k) = filter_1D(i+1)*filter_1D(j+1)
        !
      end do
    end do
    !
  else if (Filter_Option == 2) then
    !
    ! Option 2 - The quad basis functions actually contains cross polynomial
    !            terms up to order 2*p (as opposed to just 1*p for triangles).
    !            This filter applies the cut-off frequency based on 2*p as
    !            opposed to basing it on just p in 1D like Option 1.
    !
    do i = 0,n_order
      do j = 0,n_order
        !
        if (i+j >= Filter_Cutoff) then
          !
          k = n_sp1de*j + i + 1
          !
          rn = real(i+j-Filter_Cutoff,kind=lp)
          rd = real(2*n_order-Filter_Cutoff,kind=lp)
          !
          return_value(k,k) = exp(-falpha * (rn/rd)**Filter_Order)
          !
        end if
        !
      end do
    end do
    !
  end if
  !
  ! Compute the final filter matrix for quadrilaterals
  !
  return_value = matmul( vand , return_value )
  return_value = matmul( return_value , vand_inv )
  !
end function Filter2D_quad_DP
!
!###############################################################################
!
pure function Filter2D_quad_QP(vand,vand_inv) result(return_value)
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Use Statements ..
  use geovar, only : n_order,n_sp1de,n_sp2dq
  use ovar,   only : Filter_Option,Filter_Order, &
                     Filter_Cutoff,Filter_Alpha
  use ovar,   only : Filter_etac
  !
  !.. Formal Arguments ..
  real(qp), dimension(1:n_sp2dq,1:n_sp2dq), intent(in) :: vand
  real(qp), dimension(1:n_sp2dq,1:n_sp2dq), intent(in) :: vand_inv
  !
  !.. Function Return Value ..
  real(qp), dimension(1:n_sp2dq,1:n_sp2dq) :: return_value
  !
  !.. Local Scalars ..
  integer  :: i,j,k
  real(qp) :: rd,rn,falpha,fetac
  !
  !.. Local Arrays ..
  real(qp), dimension(1:n_sp1de) :: filter_1D
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one = 1.0_lp
  !
continue
  !
  ! Covert Filter_Alpha to local precision
  !
  falpha = real( Filter_Alpha , kind=lp )
  fetac  = real( Filter_etac , kind=lp )
  !
  ! Initialize the local filter matrices to identity matrices
  !
  filter_1D(:) = one
  !
  return_value(:,:) = zero
  do i = 1,n_sp2dq
    return_value(i,i) = one
  end do
  !
  ! Compute the filter diagonal for quadrilaterals
  !
  if (Filter_Option == 1) then
    !
    ! Option 1 - Tensor product of the 1D filter
    !           (Equivalent to using a 1D filter in each direction)
    !
    ! First get the 1D filter diagonal
    !
    do i = 0,n_order
      if (i >= Filter_Cutoff) then
        !
        rn = real(i,kind=lp)/real(n_order,kind=lp) - fetac
        rd = one - fetac
        !
        filter_1D(i+1) = exp(-falpha * (rn/rd)**Filter_Order)
        !
      end if
    end do
    !
    ! Now fill the filter diagonal for quads
    !
    do j = 0,n_order
      do i = 0,n_order
        !
        k = n_sp1de*j + i + 1
        !
        return_value(k,k) = filter_1D(i+1)*filter_1D(j+1)
        !
      end do
    end do
    !
  else if (Filter_Option == 2) then
    !
    ! Option 2 - The quad basis functions actually contains cross polynomial
    !            terms up to order 2*p (as opposed to just 1*p for triangles).
    !            This filter applies the cut-off frequency based on 2*p as
    !            opposed to basing it on just p in 1D like Option 1.
    !
    do i = 0,n_order
      do j = 0,n_order
        !
        if (i+j >= Filter_Cutoff) then
          !
          k = n_sp1de*j + i + 1
          !
          rn = real(i+j-Filter_Cutoff,kind=lp)
          rd = real(2*n_order-Filter_Cutoff,kind=lp)
          !
          return_value(k,k) = exp(-falpha * (rn/rd)**Filter_Order)
          !
        end if
        !
      end do
    end do
    !
  end if
  !
  ! Compute the final filter matrix for quadrilaterals
  !
  return_value = matmul( vand , return_value )
  return_value = matmul( return_value , vand_inv )
  !
end function Filter2D_quad_QP
!
!###############################################################################
!
pure function GradJacobiP_DP(n,alpha,beta,x) result(return_value)
  !
  integer,   intent(in) :: n
  real(dp), intent(in) :: alpha
  real(dp), intent(in) :: beta
  real(dp), intent(in) :: x
  !
  real(dp) :: return_value
  !
  !.. Local Scalars ..
  real(dp) :: rn,a1,b1,ab1
  !
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !
  if (n == 0) then
    return_value = zero
  else
    !
    rn  = real(n,kind=lp)
    a1  = alpha + one
    b1  = beta  + one
    ab1 = alpha + b1
    !
    return_value = sqrt(rn*(rn+ab1)) * JacobiP(n-1,a1,b1,x)
    !
  end if
  !
end function GradJacobiP_DP
!
!###############################################################################
!
pure function GradJacobiP_QP(n,alpha,beta,x) result(return_value)
  !
  integer,   intent(in) :: n
  real(qp), intent(in) :: alpha
  real(qp), intent(in) :: beta
  real(qp), intent(in) :: x
  !
  real(qp) :: return_value
  !
  !.. Local Scalars ..
  real(qp) :: rn,a1,b1,ab1
  !
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !
  if (n == 0) then
    return_value = zero
  else
    !
    rn  = real(n,kind=lp)
    a1  = alpha + one
    b1  = beta  + one
    ab1 = alpha + b1
    !
    return_value = sqrt(rn*(rn+ab1)) * JacobiP(n-1,a1,b1,x)
    !
  end if
  !
end function GradJacobiP_QP
!
!###############################################################################
!
pure function GradVandermonde1D_DP(xi) result(return_value)
  !
  real(dp), dimension(:), intent(in) :: xi
  !
  real(dp), dimension(1:size(xi),1:size(xi)) :: return_value
  !
  integer   :: i,j
  real(dp) :: rj
  !
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  return_value(:,:) = zero
  !
  do j = 2,size(xi)
    do i = 1,size(xi)
      rj = real(j-1,kind=lp)
      return_value(i,j) = sqrt(rj*(rj+one)) * JacobiP(j-2,one,one,xi(i))
    end do
  end do
  !
end function GradVandermonde1D_DP
!
!###############################################################################
!
pure function GradVandermonde1D_QP(xi) result(return_value)
  !
  real(qp), dimension(:), intent(in) :: xi
  !
  real(qp), dimension(1:size(xi),1:size(xi)) :: return_value
  !
  integer   :: i,j
  real(qp) :: rj
  !
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  return_value(:,:) = zero
  !
  do j = 2,size(xi)
    do i = 1,size(xi)
      rj = real(j-1,kind=lp)
      return_value(i,j) = sqrt(rj*(rj+one)) * JacobiP(j-2,one,one,xi(i))
    end do
  end do
  !
end function GradVandermonde1D_QP
!
!###############################################################################
!
pure function GradVandermonde2D_DP(rs) result(return_value)
  !
  real(dp), dimension(:,:), intent(in) :: rs
  !
  real(dp), dimension(1:size(rs,2),1:size(rs,2),1:size(rs,1)) :: return_value
  !
  !.. Local Scalars ..
  integer   :: i,im1,j,k,m,n,np
  real(dp) :: ri,twoiph,alpha,alphb,beta
  real(dp) :: fa,gb,dfa,dgb,h1mb,dmdr,dmds
  !
  !.. Local Arrays ..
  real(dp), dimension(1:2) :: ab
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero  = 0.0_lp
  real(lp), parameter :: half  = 0.5_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  real(lp), parameter :: sqrt2 = sqrt(two)
  !
continue
  !
  np = size(rs,2)
  n  = np2n(np)
  !
  m = 1
  do i = 0,n
    !
    ri   = real(i,kind=lp)
    im1  = max(i-1,0)
    !
    twoiph = sqrt2 * two**i
    !
    alpha = zero
    alphb = real(2*i+1,kind=lp)
    beta  = zero
    !
    do j = 0,n-i
      !
      do k = 1,np
        !
        ! Get a,b from r,s
        !
        ab = rstoab(rs(1:2,k))
        !
        fa = JacobiP(i,alpha,beta,ab(1))
        gb = JacobiP(j,alphb,beta,ab(2))
        !
        dfa = GradJacobiP(i,alpha,beta,ab(1))
        dgb = GradJacobiP(j,alphb,beta,ab(2))
        !
        h1mb = half*(one-ab(2))
        !
        dmdr = dfa * gb * h1mb**im1
        dmds = half*dmdr*(one+ab(1)) + fa*(dgb*h1mb**i - half*ri*gb*h1mb**im1)
        !
        return_value(k,m,1) = twoiph * dmdr
        return_value(k,m,2) = twoiph * dmds
        !
      end do
      !
      m = m + 1
      !
    end do
  end do
  !
end function GradVandermonde2D_DP
!
!###############################################################################
!
pure function GradVandermonde2D_QP(rs) result(return_value)
  !
  real(qp), dimension(:,:), intent(in) :: rs
  !
  real(qp), dimension(1:size(rs,2),1:size(rs,2),1:size(rs,1)) :: return_value
  !
  !.. Local Scalars ..
  integer   :: i,im1,j,k,m,n,np
  real(qp) :: ri,twoiph,alpha,alphb,beta
  real(qp) :: fa,gb,dfa,dgb,h1mb,dmdr,dmds
  !
  !.. Local Arrays ..
  real(qp), dimension(1:2) :: ab
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: zero  = 0.0_lp
  real(lp), parameter :: half  = 0.5_lp
  real(lp), parameter :: one   = 1.0_lp
  real(lp), parameter :: two   = 2.0_lp
  real(lp), parameter :: sqrt2 = sqrt(two)
  !
continue
  !
  np = size(rs,2)
  n  = np2n(np)
  !
  m = 1
  do i = 0,n
    !
    ri   = real(i,kind=lp)
    im1  = max(i-1,0)
    !
    twoiph = sqrt2 * two**i
    !
    alpha = zero
    alphb = real(2*i+1,kind=lp)
    beta  = zero
    !
    do j = 0,n-i
      !
      do k = 1,np
        !
        ! Get a,b from r,s
        !
        ab = rstoab(rs(1:2,k))
        !
        fa = JacobiP(i,alpha,beta,ab(1))
        gb = JacobiP(j,alphb,beta,ab(2))
        !
        dfa = GradJacobiP(i,alpha,beta,ab(1))
        dgb = GradJacobiP(j,alphb,beta,ab(2))
        !
        h1mb = half*(one-ab(2))
        !
        dmdr = dfa * gb * h1mb**im1
        dmds = half*dmdr*(one+ab(1)) + fa*(dgb*h1mb**i - half*ri*gb*h1mb**im1)
        !
        return_value(k,m,1) = twoiph * dmdr
        return_value(k,m,2) = twoiph * dmds
        !
      end do
      !
      m = m + 1
      !
    end do
  end do
  !
end function GradVandermonde2D_QP
!
!###############################################################################
!
pure function LagrangePoly1D_DP(xi) result(return_value)
  !
  !.. Formal Arguments ..
  real(dp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(dp), dimension(1:size(xi)) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  nfp = size(xi)
  imask(:) = (/ (i,i=1,nfp) /)
  !
  do i = 1,nfp
    return_value(i) = product( -one -xi(1:nfp) , imask /= i ) / &
                      product( xi(i)-xi(1:nfp) , imask /= i )
  end do
  !
end function LagrangePoly1D_DP
!
!###############################################################################
!
pure function LagrangePoly1D_QP(xi) result(return_value)
  !
  !.. Formal Arguments ..
  real(qp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(qp), dimension(1:size(xi)) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  nfp = size(xi)
  imask(:) = (/ (i,i=1,nfp) /)
  !
  do i = 1,nfp
    return_value(i) = product( -one -xi(1:nfp) , imask /= i ) / &
                      product( xi(i)-xi(1:nfp) , imask /= i )
  end do
  !
end function LagrangePoly1D_QP
!
!###############################################################################
!
pure function LagrangePoly2D_quad_DP(xi_fp,xi_sp) result(return_value)
  !
  !.. Formal Arguments ..
  real(dp), dimension(:), intent(in) :: xi_fp
  real(dp), dimension(:), intent(in) :: xi_sp
  !
  !.. Function Return Value ..
  real(dp), dimension(1:size(xi_sp)**2,1:4*size(xi_fp)) :: return_value
  !
  !.. Local Scalars ..
  integer   :: i,j,n,nf,nfp,nsp
  real(dp) :: x,y
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  nfp = size(xi_fp)
  nsp = size(xi_sp)
  !
  ! For each flux point, get the coefficients for interpolating
  ! a polynomial from the solution points to the flux point.
  ! These coefficients are simply the value of the 2D Lagrange
  ! polynomial for a specified solution point evaluated at
  ! the coordinates of the flux point.
  !
  ! Face 1 - bottom/south edge of quad
  !
  do n = 1,nfp
    !
    nf = n
    !
    x = xi_fp(n)
    y = -one
    !
    do j = 1,nsp
      do i = 1,nsp
        return_value(nsp*(j-1)+i,nf) = eval_LagrangePoly2D(i,j,x,y,xi_sp,xi_sp)
      end do
    end do
    !
  end do
  !
  ! Face 2 - right/east edge of quad
  !
  do n = 1,nfp
    !
    nf = n + nfp
    !
    x = one
    y = xi_fp(n)
    !
    do j = 1,nsp
      do i = 1,nsp
        return_value(nsp*(j-1)+i,nf) = eval_LagrangePoly2D(i,j,x,y,xi_sp,xi_sp)
      end do
    end do
    !
  end do
  !
  ! Face 3 - top/north edge of quad
  !
  do n = 1,nfp
    !
    nf = n + 2*nfp
    !
    x = xi_fp(nfp+1-n)
    y = one
    !
    do j = 1,nsp
      do i = 1,nsp
        return_value(nsp*(j-1)+i,nf) = eval_LagrangePoly2D(i,j,x,y,xi_sp,xi_sp)
      end do
    end do
    !
  end do
  !
  ! Face 4 - left/west edge of quad
  !
  do n = 1,nfp
    !
    nf = n + 3*nfp
    !
    x = -one
    y = xi_fp(nfp+1-n)
    !
    do j = 1,nsp
      do i = 1,nsp
        return_value(nsp*(j-1)+i,nf) = eval_LagrangePoly2D(i,j,x,y,xi_sp,xi_sp)
      end do
    end do
    !
  end do
  !
end function LagrangePoly2D_quad_DP
!
!###############################################################################
!
pure function LagrangePoly2D_quad_QP(xi_fp,xi_sp) result(return_value)
  !
  !.. Formal Arguments ..
  real(qp), dimension(:), intent(in) :: xi_fp
  real(qp), dimension(:), intent(in) :: xi_sp
  !
  !.. Function Return Value ..
  real(qp), dimension(1:size(xi_sp)**2,1:4*size(xi_fp)) :: return_value
  !
  !.. Local Scalars ..
  integer   :: i,j,n,nf,nfp,nsp
  real(qp) :: x,y
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  nfp = size(xi_fp)
  nsp = size(xi_sp)
  !
  ! For each flux point, get the coefficients for interpolating
  ! a polynomial from the solution points to the flux point.
  ! These coefficients are simply the value of the 2D Lagrange
  ! polynomial for a specified solution point evaluated at
  ! the coordinates of the flux point.
  !
  ! Face 1 - bottom/south edge of quad
  !
  do n = 1,nfp
    !
    nf = n
    !
    x = xi_fp(n)
    y = -one
    !
    do j = 1,nsp
      do i = 1,nsp
        return_value(nsp*(j-1)+i,nf) = eval_LagrangePoly2D(i,j,x,y,xi_sp,xi_sp)
      end do
    end do
    !
  end do
  !
  ! Face 2 - right/east edge of quad
  !
  do n = 1,nfp
    !
    nf = n + nfp
    !
    x = one
    y = xi_fp(n)
    !
    do j = 1,nsp
      do i = 1,nsp
        return_value(nsp*(j-1)+i,nf) = eval_LagrangePoly2D(i,j,x,y,xi_sp,xi_sp)
      end do
    end do
    !
  end do
  !
  ! Face 3 - top/north edge of quad
  !
  do n = 1,nfp
    !
    nf = n + 2*nfp
    !
    x = xi_fp(nfp+1-n)
    y = one
    !
    do j = 1,nsp
      do i = 1,nsp
        return_value(nsp*(j-1)+i,nf) = eval_LagrangePoly2D(i,j,x,y,xi_sp,xi_sp)
      end do
    end do
    !
  end do
  !
  ! Face 4 - left/west edge of quad
  !
  do n = 1,nfp
    !
    nf = n + 3*nfp
    !
    x = -one
    y = xi_fp(nfp+1-n)
    !
    do j = 1,nsp
      do i = 1,nsp
        return_value(nsp*(j-1)+i,nf) = eval_LagrangePoly2D(i,j,x,y,xi_sp,xi_sp)
      end do
    end do
    !
  end do
  !
end function LagrangePoly2D_quad_QP
!
!###############################################################################
!
pure function LagrangePoly2D_tria_DP(rs) result(return_value)
  !
  !.. Formal Arguments ..
  real(dp), dimension(:,:), intent(in) :: rs
  !
  !.. Function Return Argument ..
  real(dp), dimension(1:size(rs,2),1:3*np2nfp(size(rs,2))) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,np,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:np2nfp(size(rs,2)),1:3) :: Fmask
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
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
end function LagrangePoly2D_tria_DP
!
!###############################################################################
!
pure function LagrangePoly2D_tria_QP(rs) result(return_value)
  !
  !.. Formal Arguments ..
  real(qp), dimension(:,:), intent(in) :: rs
  !
  !.. Function Return Argument ..
  real(qp), dimension(1:size(rs,2),1:3*np2nfp(size(rs,2))) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,np,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:np2nfp(size(rs,2)),1:3) :: Fmask
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
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
end function LagrangePoly2D_tria_QP
!
!###############################################################################
!
pure function rstoab_DP(rs) result(return_value)
  !
  real(dp), dimension(1:2), intent(in) :: rs
  !
  real(dp), dimension(1:2) :: return_value
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: two  = 2.0_lp
  !
continue
  !
  return_value(1) = two * (one+rs(1)) * reciprocal(one-rs(2)) - one
  return_value(2) = rs(2)
  !
end function rstoab_DP
!
!###############################################################################
!
pure function rstoab_QP(rs) result(return_value)
  !
  real(qp), dimension(1:2), intent(in) :: rs
  !
  real(qp), dimension(1:2) :: return_value
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: two  = 2.0_lp
  !
continue
  !
  return_value(1) = two * (one+rs(1)) * reciprocal(one-rs(2)) - one
  return_value(2) = rs(2)
  !
end function rstoab_QP
!
!###############################################################################
!
pure function Vandermonde1D_DP(xi) result(return_value)
  !
  real(dp), dimension(:), intent(in) :: xi
  !
  real(dp), dimension(1:size(xi),1:size(xi)) :: return_value
  !
  integer :: i,j
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  return_value(:,:) = zero
  !
  do j = 1,size(xi)
    do i = 1,size(xi)
      return_value(i,j) = JacobiP(j-1,zero,zero,xi(i))
    end do
  end do
  !
end function Vandermonde1D_DP
!
!###############################################################################
!
pure function Vandermonde1D_QP(xi) result(return_value)
  !
  real(qp), dimension(:), intent(in) :: xi
  !
  real(qp), dimension(1:size(xi),1:size(xi)) :: return_value
  !
  integer :: i,j
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  return_value(:,:) = zero
  !
  do j = 1,size(xi)
    do i = 1,size(xi)
      return_value(i,j) = JacobiP(j-1,zero,zero,xi(i))
    end do
  end do
  !
end function Vandermonde1D_QP
!
!###############################################################################
!
pure function Vandermonde1D_N_DP(xi,N) result(return_value)
  !
  integer,                 intent(in) :: N
  real(dp), dimension(:), intent(in) :: xi
  !
  real(dp), dimension(1:size(xi),1:N) :: return_value
  !
  integer :: i,j
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  return_value(:,:) = zero
  !
  do j = 1,N
    do i = 1,size(xi)
      return_value(i,j) = JacobiP(j-1,zero,zero,xi(i))
    end do
  end do
  !
end function Vandermonde1D_N_DP
!
!###############################################################################
!
pure function Vandermonde1D_N_QP(xi,N) result(return_value)
  !
  integer,                 intent(in) :: N
  real(qp), dimension(:), intent(in) :: xi
  !
  real(qp), dimension(1:size(xi),1:N) :: return_value
  !
  integer :: i,j
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  return_value(:,:) = zero
  !
  do j = 1,N
    do i = 1,size(xi)
      return_value(i,j) = JacobiP(j-1,zero,zero,xi(i))
    end do
  end do
  !
end function Vandermonde1D_N_QP
!
!###############################################################################
!
pure function Vandermonde2D_quad_DP(xi) result(return_value)
  !
  real(dp), dimension(:), intent(in) :: xi
  !
  real(dp), dimension(1:size(xi)**2,1:size(xi)**2) :: return_value
  !
  !.. Local Scalars ..
  integer   :: k,ki,kj,l,li,lj
  integer   :: n_order,n_sp1de,n_sp2dq
  real(dp) :: x,y,hx,hy
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  n_sp1de = size(xi)
  n_sp2dq = n_sp1de*n_sp1de
  n_order = n_sp1de-1
  !
  do kj = 1,n_sp1de
    do ki = 1,n_sp1de
      !
      k = n_sp1de*(kj-1) + ki
      !
      x = xi(ki)
      y = xi(kj)
      !
      do lj = 0,n_order
        do li = 0,n_order
          !
          l = n_sp1de*lj + li + 1
          !
          hx = JacobiP(li,zero,zero,x)
          hy = JacobiP(lj,zero,zero,y)
          !
          return_value(k,l) = hx*hy
          !
        end do
      end do
      !
    end do
  end do
  !
end function Vandermonde2D_quad_DP
!
!###############################################################################
!
pure function Vandermonde2D_quad_QP(xi) result(return_value)
  !
  real(qp), dimension(:), intent(in) :: xi
  !
  real(qp), dimension(1:size(xi)**2,1:size(xi)**2) :: return_value
  !
  !.. Local Scalars ..
  integer   :: k,ki,kj,l,li,lj
  integer   :: n_order,n_sp1de,n_sp2dq
  real(qp) :: x,y,hx,hy
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  n_sp1de = size(xi)
  n_sp2dq = n_sp1de*n_sp1de
  n_order = n_sp1de-1
  !
  do kj = 1,n_sp1de
    do ki = 1,n_sp1de
      !
      k = n_sp1de*(kj-1) + ki
      !
      x = xi(ki)
      y = xi(kj)
      !
      do lj = 0,n_order
        do li = 0,n_order
          !
          l = n_sp1de*lj + li + 1
          !
          hx = JacobiP(li,zero,zero,x)
          hy = JacobiP(lj,zero,zero,y)
          !
          return_value(k,l) = hx*hy
          !
        end do
      end do
      !
    end do
  end do
  !
end function Vandermonde2D_quad_QP
!
!###############################################################################
!
pure function Vandermonde2D_quad_N_DP(xi,N) result(return_value)
  !
  integer,                 intent(in) :: N
  real(dp), dimension(:), intent(in) :: xi
  !
  real(dp), dimension(1:size(xi)**2,1:(N+1)**2) :: return_value
  !
  !.. Local Scalars ..
  integer   :: k,ki,kj,l,li,lj
  integer   :: nsp,nso,nsq,nqp,nqo,nqq
  real(dp) :: x,y,hx,hy
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  nqp = size(xi)
  nqo = nqp-1
  nqq = nqp*nqp
  !
  nsp = N+1
  nso = N
  nsq = nsp*nsp
  !
  do kj = 1,nqp
    do ki = 1,nqp
      !
      ! k is the 1D index of the solution point at
      ! which each basis function is being evaluated
      !
      k = nqp*(kj-1) + ki
      !
      x = xi(ki)
      y = xi(kj)
      !
      do lj = 0,nso
        do li = 0,nso
          !
          ! l is the index of the current basis function being evaluated
          !
          l = nsp*lj + li + 1
          !
          hx = JacobiP(li,zero,zero,x)
          hy = JacobiP(lj,zero,zero,y)
          !
          return_value(k,l) = hx*hy
          !
        end do
      end do
      !
    end do
  end do
  !
end function Vandermonde2D_quad_N_DP
!
!###############################################################################
!
pure function Vandermonde2D_quad_N_QP(xi,N) result(return_value)
  !
  integer,                 intent(in) :: N
  real(qp), dimension(:), intent(in) :: xi
  !
  real(qp), dimension(1:size(xi)**2,1:(N+1)**2) :: return_value
  !
  !.. Local Scalars ..
  integer   :: k,ki,kj,l,li,lj
  integer   :: nsp,nso,nsq,nqp,nqo,nqq
  real(qp) :: x,y,hx,hy
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  nqp = size(xi)
  nqo = nqp-1
  nqq = nqp*nqp
  !
  nsp = N+1
  nso = N
  nsq = nsp*nsp
  !
  do kj = 1,nqp
    do ki = 1,nqp
      !
      ! k is the 1D index of the solution point at
      ! which each basis function is being evaluated
      !
      k = nqp*(kj-1) + ki
      !
      x = xi(ki)
      y = xi(kj)
      !
      do lj = 0,nso
        do li = 0,nso
          !
          ! l is the index of the current basis function being evaluated
          !
          l = nsp*lj + li + 1
          !
          hx = JacobiP(li,zero,zero,x)
          hy = JacobiP(lj,zero,zero,y)
          !
          return_value(k,l) = hx*hy
          !
        end do
      end do
      !
    end do
  end do
  !
end function Vandermonde2D_quad_N_QP
!
!###############################################################################
!
pure function Vandermonde2D_tria_DP(rs) result(return_value)
  !
  real(dp), dimension(:,:), intent(in) :: rs
  !
  real(dp), dimension(1:size(rs,2),1:size(rs,2)) :: return_value
  !
  !.. Local Scalars ..
  integer   :: i,j,k,n,np,m
  real(dp) :: h1,h2,alpha,beta
  !
  !.. Local Arrays ..
  real(dp), dimension(1:2) :: ab
  !
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: two  = 2.0_lp
  !
continue
  !
  np = size(rs,2)
  n  = np2n(np)
  !
  m = 1
  do i = 0,n
    do j = 0,n-i
      !
      ! Evaluate the basis function, psi(m), at the np solution points
      ! where m = j + i*(n+1) + 1 - i*(i-1)/2
      !
      do k = 1,np
        !
        ! Get a,b from r,s
        !
        ab = rstoab(rs(1:2,k))
        !
        ! Get the value of the Jacobi polynomial for a and b
        !
        alpha = zero
        beta  = zero
        h1 = JacobiP(i,alpha,beta,ab(1))
        !
        alpha = real(2*i+1,kind=lp)
        beta  = zero
        h2 = JacobiP(j,alpha,beta,ab(2))
        !
        ! Evaluate the psi(m) at node k
        !
        return_value(k,m) = sqrt(two) * h1 * h2 * (one-ab(2))**i
        !
      end do
      !
      m = m + 1
      !
    end do
  end do
  !
end function Vandermonde2D_tria_DP
!
!###############################################################################
!
pure function Vandermonde2D_tria_QP(rs) result(return_value)
  !
  real(qp), dimension(:,:), intent(in) :: rs
  !
  real(qp), dimension(1:size(rs,2),1:size(rs,2)) :: return_value
  !
  !.. Local Scalars ..
  integer   :: i,j,k,n,np,m
  real(qp) :: h1,h2,alpha,beta
  !
  !.. Local Arrays ..
  real(qp), dimension(1:2) :: ab
  !
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: two  = 2.0_lp
  !
continue
  !
  np = size(rs,2)
  n  = np2n(np)
  !
  m = 1
  do i = 0,n
    do j = 0,n-i
      !
      ! Evaluate the basis function, psi(m), at the np solution points
      ! where m = j + i*(n+1) + 1 - i*(i-1)/2
      !
      do k = 1,np
        !
        ! Get a,b from r,s
        !
        ab = rstoab(rs(1:2,k))
        !
        ! Get the value of the Jacobi polynomial for a and b
        !
        alpha = zero
        beta  = zero
        h1 = JacobiP(i,alpha,beta,ab(1))
        !
        alpha = real(2*i+1,kind=lp)
        beta  = zero
        h2 = JacobiP(j,alpha,beta,ab(2))
        !
        ! Evaluate the psi(m) at node k
        !
        return_value(k,m) = sqrt(two) * h1 * h2 * (one-ab(2))**i
        !
      end do
      !
      m = m + 1
      !
    end do
  end do
  !
end function Vandermonde2D_tria_QP
!
!###############################################################################
!
subroutine output2d(output_requested,array,dim1,dim2,ctype,arrayname)
  !
  !.. Use Statements ..
  use ifport
  !
  !.. Formal Arguments ..
  logical,                  intent(in) :: output_requested
  integer,                  intent(in) :: dim1,dim2,ctype
  real(wp), dimension(:,:), intent(in) :: array
  character(len=*),         intent(in) :: arrayname
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,np1,np2
  integer :: ni,nj,nk,nl,ierr
  integer :: exitstat,cmndstat
  character(len=1000) :: cmndmsg,command
  character(len=6+len_trim(adjustl(arrayname))) :: fname
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(array,1),1:size(array,2)) :: copy
  !
  !.. Local Constants ..
  integer, parameter :: face = 1
  integer, parameter :: cell = 2
  integer, parameter :: iu = 18421
  !
continue
  !
  exitstat = 0
  cmndstat = 0
  cmndmsg = " "
  !
  if (mypnum == 0 .and. (output_requested .or. array_debugging)) then
    !
    fname = "array."//trim(adjustl(arrayname))
    open(iu,file=fname,status="replace",action="write",form="formatted",iostat=ierr)
    call io_error('output2d',fname,1,ierr)
    !
    np1 = size(array,1)
    np2 = size(array,2)
    !
    copy(:,:) = array(:,:)
    !
    write (iu,1) trim(adjustl(arrayname))
    !
    ! Default values
    !
    ni = np1
    nj = 1
    nk = np2
    nl = 1
    !
    if (ctype == 1) then
      !
      ni = np2n(np1)
      nj = ni
      nk = np2n(np2)
      nl = nk
      !
    else if (ctype == 2) then
      !
      if (dim1 == face) then
        ni = np1/4
        nj = 4
      else if (dim1 == cell) then
        ni = nint( sqrt(real(np1,kind=wp)) )
        nj = ni
      end if
      !
      if (dim2 == face) then
        nk = np2/4
        nl = 4
      else if (dim2 == cell) then
        nk = nint( sqrt(real(np2,kind=wp)) )
        nl = nk
      end if
      !
      do l = 1,nl
        do k = 1,nk
          write (iu,2) k,l
          do j = 1,nj
            write (iu,3) ( copy( (nj-j)*ni+i , (l-1)*nk+k ) , i=1,ni )
          end do
        end do
      end do
      !
    else
      !
      do j = 1,np1
        write (iu,3) ( copy(j,i) , i=1,np2 )
      end do
      !
    end if
    !
    flush (iu)
    !
    close (iu,iostat=ierr)
    call io_error('output2d',fname,2,ierr)
    !
    command = "sed -i -e 's/-\<0\.000\>/   -  /g' -e 's/\<0\.000\>/  -  /g' "// &
              trim(adjustl(fname))
    write (iout,4) "Executing: ",trim(adjustl(command))
    !
    exitstat = system( trim(command) )
   !call execute_command_line(trim(command),true,exitstat,cmndstat,cmndmsg)
    !
    if (exitstat /= 0 .or. cmndstat /= 0) then
      write (iout,99) exitstat,cmndstat,trim(adjustl(cmndmsg))
    end if
    !
  end if
  !
  1 format (/,a,/)
  2 format (/,"i = ",i2,",j = ",i2)
  3 format (*(g12.4))
  4 format (a,a)
  11 format (a,i0,a,i0,a)
  99 format (/," Error attempting to use command line!",/, &
               " exit status    = ",i0,/, &
               " command status = ",i0,/, &
               " error message  = ",a,//)
  !
end subroutine output2d
!
!###############################################################################
!
subroutine write_lagrange_tecplot()
  !
  !.. Use Statements ..
  use geovar, only : n_sp1de,n_sp2dq,n_qp1de,n_qp2dq
  use geovar, only : solpts_1d,qdtpts_1d,interp_s2f,interp_s2g
  use metvar, only : dr_mat,ds_mat
  use metvar, only : dqr_mat,dqs_mat
  !
  !.. Local Scalars ..
  integer   :: i,j,k,ierr
  integer, parameter :: iotec = 4928
  !
continue
  !
  open (iotec,file="lagrange.tec.dat",status="replace", &
              action="write",form="formatted",iostat=ierr)
  call io_error("write_tecplot_partitions","lagrange.tec.dat",1,ierr)
  !
  write (iotec,1)
  ! 
  write (iotec,3) 1,n_sp1de,n_sp1de
  write (iotec,2) ((solpts_1d(i),i=1,n_sp1de),j=1,n_sp1de)
  write (iotec,2) ((solpts_1d(j),i=1,n_sp1de),j=1,n_sp1de)
  write (iotec,2) ((interp_s2f((j-1)*n_sp1de+i,1,2),i=1,n_sp1de),j=1,n_sp1de)
  !
  do k = 2,4*n_sp1de
    write (iotec,4) k,n_sp1de,n_sp1de
    write (iotec,2) ((interp_s2f((j-1)*n_sp1de+i,k,2),i=1,n_sp1de),j=1,n_sp1de)
  end do
  !
  do k = 1,n_sp2dq
    write (iotec,5) k,n_sp1de,n_sp1de
    write (iotec,2) ((dr_mat((j-1)*n_sp1de+i,k,2)/maxval(dqr_mat(:,k,2)),i=1,n_sp1de),j=1,n_sp1de)
  end do
  !
  do k = 1,n_sp2dq
    write (iotec,6) k,n_sp1de,n_sp1de
    write (iotec,2) ((ds_mat((j-1)*n_sp1de+i,k,2)/maxval(dqs_mat(:,k,2)),i=1,n_sp1de),j=1,n_sp1de)
  end do
  !
  write (iotec,7) 1,n_qp1de,n_qp1de
  write (iotec,2) ((qdtpts_1d(i),i=1,n_qp1de),j=1,n_qp1de)
  write (iotec,2) ((qdtpts_1d(j),i=1,n_qp1de),j=1,n_qp1de)
  write (iotec,2) ((interp_s2g((j-1)*n_qp1de+i,1,2),i=1,n_qp1de),j=1,n_qp1de)
  !
  do k = 2,4*n_sp1de
    write (iotec,8) k,n_qp1de,n_qp1de,4*n_sp1de+2*n_sp1de*n_sp1de+1
    write (iotec,2) ((interp_s2g((j-1)*n_qp1de+i,k,2),i=1,n_qp1de),j=1,n_qp1de)
  end do
  !
  do k = 1,n_sp2dq
    write (iotec,9) k,n_qp1de,n_qp1de,4*n_sp1de+2*n_sp1de*n_sp1de+1
    write (iotec,2) ((dqr_mat((j-1)*n_qp1de+i,k,2)/maxval(dqr_mat(:,k,2)),i=1,n_qp1de),j=1,n_qp1de)
  end do
  !
  do k = 1,n_sp2dq
    write (iotec,10) k,n_qp1de,n_qp1de,4*n_sp1de+2*n_sp1de*n_sp1de+1
    write (iotec,2) ((dqs_mat((j-1)*n_qp1de+i,k,2)/maxval(dqs_mat(:,k,2)),i=1,n_qp1de),j=1,n_qp1de)
  end do
  !
  close (iotec,iostat=ierr)
  call io_error("write_tecplot_partitions","lagrange.tec.dat",2,ierr)
  !
  !
  !
  1 format ('VARIABLES = "X", "Y", "Z"')
  2 format (5es18.10)
  3 format ('ZONE T="La-SP - ',i0,'", I=',i0,', J=',i0,', K=1, DATAPACKING=BLOCK')
  4 format ('ZONE T="La-SP - ',i0,'", I=',i0,', J=',i0, &
            ', K=1, DATAPACKING=BLOCK, VARSHARELIST=([1,2]=1)')
  5 format ('ZONE T="Dr-SP - ',i0,'", I=',i0,', J=',i0, &
            ', K=1, DATAPACKING=BLOCK, VARSHARELIST=([1,2]=1)')
  6 format ('ZONE T="Ds-SP - ',i0,'", I=',i0,', J=',i0, &
            ', K=1, DATAPACKING=BLOCK, VARSHARELIST=([1,2]=1)')
  7 format ('ZONE T="La-QP - ',i0,'", I=',i0,', J=',i0,', K=1, DATAPACKING=BLOCK')
  8 format ('ZONE T="La-QP - ',i0,'", I=',i0,', J=',i0, &
            ', K=1, DATAPACKING=BLOCK, VARSHARELIST=([1,2]=',i0,')')
  9 format ('ZONE T="Dr-QP - ',i0,'", I=',i0,', J=',i0, &
            ', K=1, DATAPACKING=BLOCK, VARSHARELIST=([1,2]=',i0,')')
 10 format ('ZONE T="Ds-QP - ',i0,'", I=',i0,', J=',i0, &
            ', K=1, DATAPACKING=BLOCK, VARSHARELIST=([1,2]=',i0,')')
  !
end subroutine write_lagrange_tecplot
!
!###############################################################################
!
subroutine write_quadrature()
  !
  !.. Use Statements ..
  use geovar, only : n_sp1de,solpts_1d,flxpts_1d
  !
  !.. Local Scalars ..
  integer :: i,j,ierr
  !
  !.. Local Constants ..
  integer, parameter :: iqd = 18
  !
continue
  !
  open(iqd,file="quad.tec.dat",status="replace",action="write",iostat=ierr)
  call io_error("write_quadrature","quad.tec.dat",1,ierr)
  !
  write (iqd,*) 'VARIABLES = "X" "Y"'
  write (iqd,*) 'ZONE T="Solution Points"'
  do j = 1,n_sp1de
    do i = 1,n_sp1de
      write (iqd,1) c0_25*( -one*(one-solpts_1d(j))*(one-solpts_1d(i)) + &
                             one*(one-solpts_1d(j))*(one+solpts_1d(i)) + &
                             one*(one+solpts_1d(j))*(one+solpts_1d(i)) - &
                             one*(one+solpts_1d(j))*(one-solpts_1d(i)) ), &
                    c0_25*( -one*(one-solpts_1d(j))*(one-solpts_1d(i)) - &
                             one*(one-solpts_1d(j))*(one+solpts_1d(i)) + &
                             one*(one+solpts_1d(j))*(one+solpts_1d(i)) + &
                             one*(one+solpts_1d(j))*(one-solpts_1d(i)) )
    end do
  end do
  write (iqd,*) 'ZONE T="Flux Points"'
  do i = 1,n_sp1de
    write (iqd,1) flxpts_1d(i),-one
  end do
  do i = 1,n_sp1de
    write (iqd,1) one,flxpts_1d(i)
  end do
  do i = 1,n_sp1de
    write (iqd,1) flxpts_1d(n_sp1de+1-i),one
  end do
  do i = 1,n_sp1de
    write (iqd,1) -one,flxpts_1d(n_sp1de+1-i)
  end do
  !
  close (iqd,iostat=ierr)
  call io_error("write_quadrature","quad.tec.dat",2,ierr)
  !
  1 format (25es12.4)
  !
end subroutine write_quadrature
!
!###############################################################################
!
include "Functions/np2n.f90"
include "Functions/np2nfp.f90"
include "Functions/chop.f90"
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
!  OLD (NON-QUAD PRECISION) VERSIONS OF:
!         get_quadrature_data
!         gauss_lobatto_matrices
!         interp2quadrature
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine old_get_quadrature_data()
  !
  !.. Use Statements ..
  use geovar, only : nr,n_sp1de,n_sp2dt,n_sp2dq,solpts_1d,solwts_1d
  use geovar, only : rs,flxpts_1d,n_order,interp_s2f,filter
  use geovar, only : vand_1d,vand_1d_inv
  use geovar, only : vand_quad,vand_quad_inv
  use geovar, only : vand_tria,vand_tria_inv
  use metvar, only : cormat,d1_mat,dr_mat,ds_mat,dr_hat,ds_hat,cormat_dg
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,n,nf,ierr
  integer :: li,lj,ki,kj
  real(wp) :: time0,time1
  !
  !.. Local Allocatable Arrays ..
  real(wp), dimension(1:nr,1:n_sp2dt) :: adjd_rs
  real(wp), dimension(1:n_sp1de,1:n_sp1de) :: grad_vand_1d
  real(wp), dimension(1:n_sp2dt,1:n_sp2dt,1:nr) :: grad_vand_tria
  !
continue
  !
  adjd_rs(:,:) = zero
  grad_vand_1d(:,:) = zero
  grad_vand_tria(:,:,:) = zero
  !
  ! Allocate the working precision arrays
  !
  allocate ( vand_1d(1:n_sp1de,1:n_sp1de) , source=zero , stat=ierr )
  call alloc_error("old_get_quadrature_data","vand_1d",1,ierr)
  allocate ( vand_tria(1:n_sp2dt,1:n_sp2dt) , source=zero , stat=ierr )
  call alloc_error("old_get_quadrature_data","vand_tria",1,ierr)
  allocate ( vand_quad(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("old_get_quadrature_data","vand_quad",1,ierr)
  allocate ( vand_1d_inv(1:n_sp1de,1:n_sp1de) , source=zero , stat=ierr )
  call alloc_error("old_get_quadrature_data","vand_1d_inv",1,ierr)
  allocate ( vand_tria_inv(1:n_sp2dt,1:n_sp2dt) , source=zero , stat=ierr )
  call alloc_error("old_get_quadrature_data","vand_tria_inv",1,ierr)
  allocate ( vand_quad_inv(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("old_get_quadrature_data","vand_quad_inv",1,ierr)
  !
  allocate ( d1_mat(1:n_sp1de,1:n_sp1de) , source=zero , stat=ierr )
  call alloc_error("old_get_quadrature_data","d1_mat",1,ierr)
  allocate ( dr_mat(1:n_sp2dq,1:n_sp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("old_get_quadrature_data","dr_mat",1,ierr)
  allocate ( ds_mat(1:n_sp2dq,1:n_sp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("old_get_quadrature_data","ds_mat",1,ierr)
  allocate ( dr_hat(1:n_sp2dq,1:n_sp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("old_get_quadrature_data","dr_hat",1,ierr)
  allocate ( ds_hat(1:n_sp2dq,1:n_sp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("old_get_quadrature_data","ds_hat",1,ierr)
  !
  allocate ( interp_s2f(1:n_sp2dq,1:4*n_sp1de,1:2) , source=zero , stat=ierr )
  call alloc_error("old_get_quadrature_data","interp_s2f",1,ierr)
  !
  allocate ( cormat(1:n_sp2dq,1:4*n_sp1de,1:2) , source=zero , stat=ierr )
  call alloc_error("old_get_quadrature_data","cormat",1,ierr)
  allocate ( cormat_dg(1:n_sp2dq,1:4*n_sp1de,1:2) , source=zero , stat=ierr )
  call alloc_error("old_get_quadrature_data","cormat_dg",1,ierr)
  !
  allocate ( filter(1:n_sp2dq,1:n_sp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("old_get_quadrature_data","filter",1,ierr)
  !
  ! Adjust r and s to be on the interval [-1,1]. The routines
  ! Vandermonde2d_tria, GradVandermonde2D, and CorrectionMatrix2D_tria
  ! need to be adjusted to work for r and s on the interval [0,1].
  !
  adjd_rs(:,:) = two*rs(:,:) - one
  !
  ! 1D Vandermonde matrix and inverse
  !
  vand_1d = Vandermonde1d(solpts_1d)
  vand_1d_inv = invert_matrix(vand_1d)
  !
  ! 2D triangle Vandermonde matrix and inverse
  !
  vand_tria = Vandermonde2d_tria(adjd_rs)
  vand_tria_inv = invert_matrix(vand_tria)
  !
  ! 2D quadrilateral Vandermonde matrix and inverse
  !
  vand_quad = Vandermonde2d_quad(solpts_1d)
  vand_quad_inv = invert_matrix(vand_quad)
  call output2d(array_debugging,vand_quad(:,:),2,2,0,'vand_quad')
  call output2d(array_debugging,vand_quad_inv(:,:),2,2,0,'vand_quad_inv')
  !
  ! Get the derivative of the 1D Vandermonde matrix
  !
  grad_vand_1d = GradVandermonde1D(solpts_1d)
  !
  ! Get the derivative of the Vandermonde matrix
  !
  grad_vand_tria = GradVandermonde2D(adjd_rs)
  !
  ! Get the 1D derivative matrix
  !
  d1_mat = matmul( grad_vand_1d , vand_1d_inv )
  !
  ! Get the derivative matrices for the strong form of the equations
  !
  call Dmatrix2d_tria(grad_vand_tria,vand_tria_inv, &
                      dr_mat(1:n_sp2dt,1:n_sp2dt,1), &
                      ds_mat(1:n_sp2dt,1:n_sp2dt,1))
  !
  call Dmatrix2d_quad(d1_mat,dr_mat(:,:,2),ds_mat(:,:,2))
  !
  ! Get the derivative matrices for the weak form of the equations
  !
  call Dmatrix2d_tria_weak(adjd_rs,vand_tria, &
                           dr_hat(1:n_sp2dt,1:n_sp2dt,1), &
                           ds_hat(1:n_sp2dt,1:n_sp2dt,1))
  !
  call Dmatrix2d_quad_weak(solwts_1d,d1_mat,dr_hat(:,:,2),ds_hat(:,:,2))
  !
  ! Get the Lagrange polynomial to project the variables/fluxes from the
  ! solution points in the cell interior to the flux points on the cell face
  !
  interp_s2f(1:n_sp2dt,1:3*n_sp1de,1) = LagrangePoly2D_tria(rs)
  interp_s2f(1:n_sp2dq,1:4*n_sp1de,2) = LagrangePoly2D_quad(flxpts_1d,solpts_1d)
  !
  call output2d(array_debugging,interp_s2f(1:n_sp2dq,1:4*n_sp1de,2),2,1,2,'interp_s2f')
  !
  ! Get the correction matrix 
  !
  cormat(1:n_sp2dt,1:3*n_sp1de,1) = CorrectionMatrix2D_tria(adjd_rs)
  cormat(1:n_sp2dq,1:4*n_sp1de,2) = CorrectionMatrix2D_quad(solpts_1d)
  !
  ! Get the correction matrix for DG methods
  !
  cormat_dg(1:n_sp2dq,1:4*n_sp1de,2) = &
           CorrectionMatrix2D_quad_dg(solpts_1d,solwts_1d)
  !
  ! Get the filter matrices
  !
  filter(1:n_sp2dt,1:n_sp2dt,1) = Filter2D_tria(vand_tria,vand_tria_inv)
  filter(1:n_sp2dq,1:n_sp2dq,2) = Filter2D_quad(vand_quad,vand_quad_inv)
  !
  ! Get the matrices to switch between the Gauss and Lobatto nodes
  !
  call old_gauss_lobatto_matrices
  !
  call old_interp2quadrature
  !
  call output2d(array_debugging,vand_1d(:,:),2,2,0,'vand_1d')
  call output2d(array_debugging,vand_1d_inv(:,:),2,2,0,'vand_1d_inv')
  call output2d(array_debugging,grad_vand_1d(:,:),2,2,0,'grad_vand_1d')
  call output2d(array_debugging,d1_mat(:,:),2,2,0,'d1_mat')
  call output2d(array_debugging,cormat(1:n_sp2dq,1:4*n_sp1de,2),2,1,2,'cormat')
  call output2d(array_debugging,cormat_dg(1:n_sp2dq,1:4*n_sp1de,2),2,1,2,'cormat_dg')
  call output2d(array_debugging,filter(1:n_sp2dq,1:n_sp2dq,2),2,2,2,'filter')
  !
 !!
 !! Output the solution and flux points to a tecplot file for visualization
 !!
 !call write_quadrature
 !!
 !write (101,*) "Quad Vandermonde Matrix"
 !1 format (a,i0,a,i0,a)
 !2 format (*(es12.4))
 !do lj = 0,n_order
 !  do li = 0,n_order
 !    l = n_sp1de*lj + li + 1
 !    write (101,*)
 !    write (101,1) "Basis function for term (x^",li,")*(y^",lj,")"
 !    do kj = n_sp1de,1,-1
 !      write (101,2) (vand_quad(n_sp1de*(kj-1)+ki,l),ki=1,n_sp1de)
 !    end do
 !  end do
 !end do
 !!
 !!
 !write (100,*) "Inverse of Quad Vandermonde Matrix"
 !do lj = 0,n_order
 !  do li = 0,n_order
 !    l = n_sp1de*lj + li + 1
 !    write (100,*)
 !    write (100,1) "Basis function for term (x^",li,")*(y^",lj,")"
 !    do kj = n_sp1de,1,-1
 !      write (100,2) (vand_quad_inv(n_sp1de*(kj-1)+ki,l),ki=1,n_sp1de)
 !    end do
 !  end do
 !end do
 !!
 !!
 !write (101,*) "Correction Function"
 !write (102,*) "DG Correction Function"
 !do nf = 1,4*n_sp1de
 !  !
 !  if (nf <= 1*n_sp1de) then
 !    write (101,5) nf
 !    write (102,5) nf
 !  else if (nf <= 2*n_sp1de) then
 !    write (101,6) nf - n_sp1de
 !    write (102,6) nf - n_sp1de
 !  else if (nf <= 3*n_sp1de) then
 !    write (101,7) n_sp1de+1-(nf-2*n_sp1de)
 !    write (102,7) n_sp1de+1-(nf-2*n_sp1de)
 !  else
 !    write (101,8) n_sp1de+1-(nf-3*n_sp1de)
 !    write (102,8) n_sp1de+1-(nf-3*n_sp1de)
 !  end if
 !  do j = n_sp1de,1,-1
 !    write (101,4) (cormat((j-1)*n_sp1de+i,nf,2),i=1,n_sp1de)
 !    write (102,4) (cormat_dg((j-1)*n_sp1de+i,nf,2),i=1,n_sp1de)
 !  end do
 !  !
 !end do
 !4 format (*(es12.4))
 !5 format (/," Bottom/South edge for xi = ",i0)
 !6 format (/," Right/East edge for eta = ",i0)
 !7 format (/," Top/North edge for xi = ",i0)
 !8 format (/," Left/West edge for eta = ",i0)
 !do j = 1,3*n_sp1de
 !  write (201,1) (cormat(i,j,1),i=1,n_sp2dt)
 !end do
 !write (102,*) "Dr matrix"
 !write (103,*) "Ds matrix"
 !do j = 1,n_sp2dq
 !  write (102,1) (dr_mat(i,j,2),i=1,n_sp2dq)
 !  write (103,1) (ds_mat(i,j,2),i=1,n_sp2dq)
 !end do
 !write (106,*) "Dr matrix"
 !write (107,*) "Ds matrix"
 !do l = 1,n_sp1de
 !  do k = 1,n_sp1de
 !    write (106,2) k,l
 !    write (107,2) k,l
 !    do j = 1,n_sp1de
 !      write (106,3) (dr_mat( (n_sp1de-j)*n_sp1de+i , (l-1)*n_sp1de+k , 2 ),i=1,n_sp1de)
 !      write (107,3) (ds_mat( (n_sp1de-j)*n_sp1de+i , (l-1)*n_sp1de+k , 2 ),i=1,n_sp1de)
 !    end do
 !  end do
 !end do
 !2 format (/,"i = ",i2,", j = ",i2)
 !3 format (*(es12.4))
 !!
 !write (110,*) "interp_s2f"
 !do l = 1,4
 !  do k = 1,n_sp1de
 !    write (110,4) l,k
 !    do j = 1,n_sp1de
 !      write (110,3) (interp_s2f( (n_sp1de-j)*n_sp1de+i , (l-1)*n_sp1de+k , 2 ),i=1,n_sp1de)
 !    end do
 !  end do
 !end do
 !4 format (/,"Face = ",i2,", i = ",i2)
 !!
 !write (108,*) "1D Dr"
 !write (109,*) "1D Ds"
 !do j = 1,n_sp1de
 !  write (108,3) (d1_mat(i,n_sp1de+1-j),i=1,n_sp1de)
 !  write (109,3) (d1_mat(n_sp1de+1-j,i),i=1,n_sp1de)
 !end do
 !write (*,*)
 !write (*,*)
 !do i = 1,n_sp1de
 !  write (*,3) (dlagrange(k,solpts_1d(i),solpts_1d),k=1,n_sp1de)
 !end do
 !write (*,*)
 !write (*,*)
 !!
 !do j = 1,n_sp2dt
 !  write (202,1) (dr_mat(i,j,1),i=1,n_sp2dt)
 !  write (203,1) (ds_mat(i,j,1),i=1,n_sp2dt)
 !end do
 !write (105,*) "interp_s2f matrix"
 !do j = 1,4*n_sp1de
 !  write (*,1) (interp_s2f(i,j,2),i=1,n_sp2dq)
 !end do
 !!
 !1 format (<n_sp2dq>es12.4)
 !!
 !call stop_recon(stop_mpi,"old_get_quadrature_data")
 !!
end subroutine old_get_quadrature_data
!
!###############################################################################
!
subroutine old_gauss_lobatto_matrices()
  !
  !.. Use Statements ..
  use geovar, only : n_order,n_sp1de,n_sp2dt,n_sp2dq
  use geovar, only : gauss_1d,lobatto_1d,radau_1d
  use geovar, only : exact2numer,numer2exact
  use geovar, only : radau2numer,numer2radau
  use geovar, only : gauss2numer,numer2gauss
  use geovar, only : lobatto2numer,numer2lobatto
  use geovar, only : vand_quad,vand_quad_inv
  use ovar,   only : error_loc,itestcase
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,ierr
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable, dimension(:,:) :: gauss_vand
  real(wp), allocatable, dimension(:,:) :: gauss_vand_inv
  real(wp), allocatable, dimension(:,:) :: lobatto_vand
  real(wp), allocatable, dimension(:,:) :: lobatto_vand_inv
  real(wp), allocatable, dimension(:,:) :: radau_vand
  real(wp), allocatable, dimension(:,:) :: radau_vand_inv
  !
  !.. Local Parameters ..
  real(wp), dimension(1:2,1:3) :: &
        tnodes = reshape([zero,zero,one,zero,zero,one],[2,3])
  real(wp), dimension(1:2,1:4) :: &
        qnodes = reshape([-one,-one,one,-one,one,one,-one,one],[2,4])
  !
continue
  !
  ! Get the Vandermonde matrix and its inverse for the right Radau nodes
  !
  allocate ( radau_vand(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("old_gauss_lobatto_matrices","radau_vand",1,ierr)
  allocate ( radau_vand_inv(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("old_gauss_lobatto_matrices","radau_vand_inv",1,ierr)
  !
  radau_vand = Vandermonde2d_quad(radau_1d)
  radau_vand_inv = invert_matrix(radau_vand)
  !
  ! Get the Vandermonde matrix and its inverse for the Gauss nodes
  !
  allocate ( gauss_vand(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("old_gauss_lobatto_matrices","gauss_vand",1,ierr)
  allocate ( gauss_vand_inv(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("old_gauss_lobatto_matrices","gauss_vand_inv",1,ierr)
  !
  gauss_vand = Vandermonde2d_quad(gauss_1d)
  gauss_vand_inv = invert_matrix(gauss_vand)
  !
  ! Get the Vandermonde matrix and its inverse for the Lobatto nodes
  !
  allocate ( lobatto_vand(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("old_gauss_lobatto_matrices","lobatto_vand",1,ierr)
  allocate ( lobatto_vand_inv(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("old_gauss_lobatto_matrices","lobatto_vand_inv",1,ierr)
  !
  lobatto_vand = Vandermonde2d_quad(lobatto_1d)
  lobatto_vand_inv = invert_matrix(lobatto_vand)
  !
  ! Now get the conversion matrices through matrix multiplication
  !
 !if (any(itestcase == Advec_Problems)) then
    !
    allocate ( radau2numer(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
    call alloc_error("old_gauss_lobatto_matrices","radau2numer",1,ierr)
    allocate ( numer2radau(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
    call alloc_error("old_gauss_lobatto_matrices","numer2radau",1,ierr)
    !
    radau2numer = matmul( vand_quad , radau_vand_inv )
    numer2radau = matmul( radau_vand , vand_quad_inv )
    !
    allocate ( gauss2numer(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
    call alloc_error("old_gauss_lobatto_matrices","gauss2numer",1,ierr)
    allocate ( numer2gauss(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
    call alloc_error("old_gauss_lobatto_matrices","numer2gauss",1,ierr)
    !
    gauss2numer = matmul( vand_quad , gauss_vand_inv )
    numer2gauss = matmul( gauss_vand , vand_quad_inv )
    !
    allocate ( lobatto2numer(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
    call alloc_error("old_gauss_lobatto_matrices","lobatto2numer",1,ierr)
    allocate ( numer2lobatto(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
    call alloc_error("old_gauss_lobatto_matrices","numer2lobatto",1,ierr)
    !
    lobatto2numer = matmul( vand_quad , lobatto_vand_inv )
    numer2lobatto = matmul( lobatto_vand , vand_quad_inv )
    !
 !else
    !
    allocate ( exact2numer(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
    call alloc_error("old_gauss_lobatto_matrices","exact2numer",1,ierr)
    allocate ( numer2exact(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
    call alloc_error("old_gauss_lobatto_matrices","numer2exact",1,ierr)
    !
    numer2exact = merge( matmul(gauss_vand,vand_quad_inv), &
                         matmul(radau_vand,vand_quad_inv), &
                         error_loc > 0 )
    !
    exact2numer = merge( matmul(vand_quad,gauss_vand_inv), &
                         matmul(vand_quad,radau_vand_inv), &
                         error_loc > 0 )
    !
 !end if
  !
  deallocate ( gauss_vand , stat=ierr )
  call alloc_error("old_gauss_lobatto_matrices","gauss_vand",2,ierr)
  deallocate ( gauss_vand_inv , stat=ierr )
  call alloc_error("old_gauss_lobatto_matrices","gauss_vand_inv",2,ierr)
  deallocate ( lobatto_vand , stat=ierr )
  call alloc_error("old_gauss_lobatto_matrices","lobatto_vand",2,ierr)
  deallocate ( lobatto_vand_inv , stat=ierr )
  call alloc_error("old_gauss_lobatto_matrices","lobatto_vand_inv",2,ierr)
  deallocate ( radau_vand , stat=ierr )
  call alloc_error("old_gauss_lobatto_matrices","radau_vand",2,ierr)
  deallocate ( radau_vand_inv , stat=ierr )
  call alloc_error("old_gauss_lobatto_matrices","radau_vand_inv",2,ierr)
  !
  call output2d(array_debugging,radau2numer(:,:),2,2,2,'radau2numer')
  call output2d(array_debugging,numer2radau(:,:),2,2,2,'numer2radau')
  call output2d(array_debugging,gauss2numer(:,:),2,2,2,'gauss2numer')
  call output2d(array_debugging,numer2gauss(:,:),2,2,2,'numer2gauss')
  call output2d(array_debugging,lobatto2numer(:,:),2,2,2,'lobatto2numer')
  call output2d(array_debugging,numer2lobatto(:,:),2,2,2,'numer2lobatto')
  call output2d(array_debugging,numer2exact(:,:),2,2,2,'numer2exact')
  call output2d(array_debugging,exact2numer(:,:),2,2,2,'exact2numer')
  !
  1 format ('VARIABLES = "X" "Y"')
  2 format ('ZONE T="Reference Quad", I=2, J=2')
  3 format ('ZONE T="Gauss Quad", I=',i0,', J=',i0)
  4 format ('ZONE T="Lobatto Quad", I=',i0,', J=',i0)
  5 format ('ZONE T="Reference Tri", N=3, E=1, DATAPACKING=POINT, ZONETYPE=FETRIANGLE')
  6 format ('ZONE T="Gauss Tri"')
  7 format ('ZONE T="Lobatto Tri"')
  8 format (2es23.15)
  9 format (3i4)
  !
end subroutine old_gauss_lobatto_matrices
!
!###############################################################################
!
subroutine old_interp2quadrature()
  !
  !.. Use Statements ..
  use geovar, only : n_order,n_sp1de,n_sp2dt,n_sp2dq
  use geovar, only : q_order,n_qp1de,n_qp2dt,n_qp2dq
  use geovar, only : qdtpts_1d,qdtwts_1d,qdtwts_2d
  use geovar, only : solpts_1d,solwts_1d,invM_mat,weights2d
  use geovar, only : interp_f2g,projct_g2f,vand_qdtr
  use geovar, only : interp_s2q,projct_q2s,vand_qdtr_inv
  use geovar, only : interp_s2g,using_quadrature
  use geovar, only : vand_1d,vand_1d_inv
  use geovar, only : vand_quad,vand_quad_inv
  use ovar,   only : iop_qdtpts
  use metvar, only : dr_mat,ds_mat,dq1D_mat
  use metvar, only : dqr_mat,dqr_LPmat,dqr_CRmat
  use metvar, only : dqs_mat,dqs_LPmat,dqs_CRmat
  use polynomial_mod, only : nodes_legendre_gauss, &
                             weights_legendre_gauss
  use polynomial_mod, only : nodes_legendre_gauss_lobatto, &
                             weights_legendre_gauss_lobatto
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,ierr
  integer :: li,lj,ki,kj,nsp,nqp
  !
  !.. Local Arrays ..
  real(wp), dimension(1:n_qp1de) :: vp
  real(wp), dimension(1:n_sp2dq,1:n_sp2dq) :: tmpmat
  real(wp), dimension(1:n_qp1de,1:n_qp1de) :: vand1D_qp
  real(wp), dimension(1:n_qp1de,1:n_qp1de) :: dv1D_qp
  real(wp), dimension(1:n_qp1de,1:n_sp1de) :: vand1D_f2g
  real(wp), dimension(1:n_qp2dq,1:n_sp2dq,1:2) :: vand2D_s2q
  !
continue
  !
  if (using_quadrature .and. mypnum == 0) write (iout,1)
  !
  ! ###########################################################
  ! ###########################################################
  ! ##                                                       ##
  ! ##  Allocate all the arrays needed for quadrature in 2D  ##
  ! ##                                                       ##
  ! ###########################################################
  ! ###########################################################
  !
  allocate ( qdtpts_1d(1:n_qp1de) , source=zero , stat=ierr )
  call alloc_error("old_interp2quadrature","qdtpts_1d",1,ierr)
  allocate ( qdtwts_1d(1:n_qp1de) , source=zero , stat=ierr )
  call alloc_error("old_interp2quadrature","qdtwts_1d",1,ierr)
  allocate ( dq1D_mat(1:n_qp1de,1:n_qp1de), source=zero , stat=ierr )
  call alloc_error("old_interp2quadrature","dq1D_mat",1,ierr)
  allocate ( interp_f2g(1:n_sp1de,1:n_qp1de) , source=zero , stat=ierr )
  call alloc_error("old_interp2quadrature","interp_f2g",1,ierr)
  allocate ( projct_g2f(1:n_qp1de,1:n_sp1de) , source=zero , stat=ierr )
  call alloc_error("old_interp2quadrature","projct_g2f",1,ierr)
  !
  ! ####################################
  ! ####################################
  ! ##                                ##
  ! ##  Section for quadrature in 1D  ##
  ! ##                                ##
  ! ####################################
  ! ####################################
  !
  ! 1D quadrature nodes and weights
  !
  if (iop_qdtpts == Legendre_Gauss_Lobatto) then
    call nodes_legendre_gauss_lobatto(q_order,qdtpts_1d,vp)
    call weights_legendre_gauss_lobatto(q_order,vp,qdtwts_1d)
  else
    call nodes_legendre_gauss(n_qp1de,qdtpts_1d,vp)
    call weights_legendre_gauss(n_qp1de,qdtpts_1d,vp,qdtwts_1d)
  end if
  !
  ! 1D Derivative matrix for the quadrature points
  !
  vand1d_qp(:,:) = Vandermonde1D(qdtpts_1d)
  dv1d_qp(:,:) = GradVandermonde1D(qdtpts_1d)
  dq1D_mat(:,:) = Dmatrix1D(vand1d_qp,dv1d_qp)
  !
  ! 1D Vandermonde matrix giving the value for each of the
  ! #(n_sp1de) modes at the #(n_qp1de) quadrature points
  !
  vand1d_f2g(:,:) = Vandermonde1D_N(qdtpts_1d,n_sp1de)
  !
  ! Matrix for interpolating 1D nodal values at face
  ! flux/solution points to face quadrature points
  !
  interp_f2g(:,:) = transpose( matmul( vand1d_f2g , vand_1d_inv ) )
  !
  call output2d(array_debugging,interp_f2g(:,:),2,2,3,'interp_f2g')
  !
  ! Matrix for projecting nodal values at face quadrature
  ! points to the face flux/solution points
  !
  projct_g2f(:,:) = matmul( vand1D_f2g , transpose(vand_1D) )
  do i = 1,n_sp1de
    projct_g2f(:,i) = projct_g2f(:,i) * qdtwts_1d(:)
  end do
  !
  call output2d(array_debugging,projct_g2f(:,:),2,2,3,'projct_g2f')
  !
  ! ###########################################################
  ! ###########################################################
  ! ##                                                       ##
  ! ##  Allocate all the arrays needed for quadrature in 2D  ##
  ! ##                                                       ##
  ! ###########################################################
  ! ###########################################################
  !
  allocate ( qdtwts_2d(1:n_qp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("old_interp2quadrature","qdtwts_2d",1,ierr)
  allocate ( vand_qdtr(1:n_qp2dq,1:n_qp2dq,1:2), source=zero , stat=ierr )
  call alloc_error("old_interp2quadrature","vand_qdtr",1,ierr)
  allocate ( vand_qdtr_inv(1:n_qp2dq,1:n_qp2dq,1:2), source=zero , stat=ierr )
  call alloc_error("old_interp2quadrature","vand_qdtr_inv",1,ierr)
  allocate ( interp_s2q(1:n_sp2dq,1:n_qp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("old_interp2quadrature","interp_s2q",1,ierr)
  allocate ( projct_q2s(1:n_qp2dq,1:n_sp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("old_interp2quadrature","projct_q2s",1,ierr)
  allocate ( dqr_mat(1:n_qp2dq,1:n_qp2dq,1:2), source=zero , stat=ierr )
  call alloc_error("old_interp2quadrature","dqr_mat",1,ierr)
  allocate ( dqs_mat(1:n_qp2dq,1:n_qp2dq,1:2), source=zero , stat=ierr )
  call alloc_error("old_interp2quadrature","dqs_mat",1,ierr)
  allocate ( dqr_LPmat(1:n_qp2dq,1:n_sp2dq,1:2), source=zero , stat=ierr )
  call alloc_error("old_interp2quadrature","dqr_LPmat",1,ierr)
  allocate ( dqs_LPmat(1:n_qp2dq,1:n_sp2dq,1:2), source=zero , stat=ierr )
  call alloc_error("old_interp2quadrature","dqs_LPmat",1,ierr)
  allocate ( dqr_CRmat(1:n_sp2dq,1:n_qp2dq,1:2), source=zero , stat=ierr )
  call alloc_error("old_interp2quadrature","dqr_CRmat",1,ierr)
  allocate ( dqs_CRmat(1:n_sp2dq,1:n_qp2dq,1:2), source=zero , stat=ierr )
  call alloc_error("old_interp2quadrature","dqs_CRmat",1,ierr)
  allocate ( invM_mat(1:n_sp2dq,1:n_sp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("old_interp2quadrature","invM_mat",1,ierr)
  allocate ( interp_s2g(1:n_sp2dq,1:4*n_qp1de,1:2) , source=zero , stat=ierr )
  call alloc_error("old_interp2quadrature","interp_s2g",1,ierr)
  !
  ! ################################################
  ! ################################################
  ! ##                                            ##
  ! ##  Section for quadrature on triangle cells  ##
  ! ##                                            ##
  ! ################################################
  ! ################################################
  !
  nsp = n_sp2dt
  nqp = n_qp2dt
  !
  ! #####################################################
  ! #####################################################
  ! ##                                                 ##
  ! ##  Section for quadrature on quadrilateral cells  ##
  ! ##                                                 ##
  ! #####################################################
  ! #####################################################
  !
  nsp = n_sp2dq
  nqp = n_qp2dq
  !
  ! 2D quadrature weights for the quadrature points
  !
  k = 0
  do i = 1,n_qp1de
    do j = 1,n_qp1de
      k = k + 1
      qdtwts_2d(k,2) = qdtwts_1d(i)*qdtwts_1d(j)
    end do
  end do
  !
  ! Full Vandermonde matrix for the quadrature nodes and modes
  !
  vand_qdtr(1:nqp,1:nqp,2) = Vandermonde2D_quad(qdtpts_1d)
  vand_qdtr_inv(1:nqp,1:nqp,2) = invert_matrix(vand_qdtr(1:nqp,1:nqp,2))
  !
  ! 2D matrix for interpolating nodal values at SPs to QPs
  !
  vand2D_s2q(1:nqp,1:nsp,2) = Vandermonde2D_quad_N(qdtpts_1d,n_order)
  interp_s2q(1:nsp,1:nqp,2) = &
        transpose( matmul( vand2D_s2q(1:nqp,1:nsp,2) , vand_quad_inv ) )
  !
  call output2d(array_debugging,interp_s2q(1:nsp,1:nqp,2),2,2,2,'interp_s2q')
  !
  ! Matrix for interpolating cell nodal values at SPs to face QPs
  !
  interp_s2g(1:nsp,1:4*n_qp1de,2) = LagrangePoly2D_quad(qdtpts_1d,solpts_1d)
  !
  call output2d(array_debugging,interp_s2g(1:nsp,1:4*n_qp1de,2),2,1,2,'interp_s2g')
  !
  ! 2D matrix to project nodal values at QPs to SPs
  !
  projct_q2s(1:nqp,1:nsp,2) = &
        matmul( vand2D_s2q(1:nqp,1:nsp,2) , transpose(vand_quad) )
  do i = 1,nsp
    projct_q2s(1:nqp,i,2) = projct_q2s(1:nqp,i,2) * qdtwts_2d(1:nqp,2)
  end do
  !
  call output2d(array_debugging,projct_q2s(1:nqp,1:nsp,2),2,2,2,'projct_q2s')
  !
  ! 2D Derivative matrix for the quadrature points
  !
  call output2d(array_debugging,dr_mat(1:nsp,1:nsp,2),2,2,2,'dr_mat')
  call output2d(array_debugging,ds_mat(1:nsp,1:nsp,2),2,2,2,'ds_mat')
  !
  call Dmatrix2D_quad(dq1d_mat,dqr_mat(1:nqp,1:nqp,2),dqs_mat(1:nqp,1:nqp,2))
  !
  call output2d(array_debugging,dqr_mat(1:nqp,1:nqp,2),2,2,2,'dqr_mat')
  call output2d(array_debugging,dqs_mat(1:nqp,1:nqp,2),2,2,2,'dqs_mat')
  !
  ! Pre-multiply the derivative matrices by the projection
  ! matrix for when using the derivative of the Lagrange polynomial
  !
  dqr_LPmat(1:nqp,1:nsp,2) = &
      matmul( transpose(dqr_mat(1:nqp,1:nqp,2)) , projct_q2s(1:nqp,1:nsp,2) )
  dqs_LPmat(1:nqp,1:nsp,2) = &
      matmul( transpose(dqs_mat(1:nqp,1:nqp,2)) , projct_q2s(1:nqp,1:nsp,2) )
  !
  call output2d(array_debugging,dqr_LPmat(1:nqp,1:nsp,2),2,2,2,'dqr_LPmat')
  call output2d(array_debugging,dqs_LPmat(1:nqp,1:nsp,2),2,2,2,'dqs_LPmat')
  !
  ! Pre-multiply the derivative matrices by the
  ! interpolation matrix for when using the chain-rule
  !
  dqr_CRmat(1:nsp,1:nqp,2) = &
      matmul( interp_s2q(1:nsp,1:nqp,2) , transpose(dqr_mat(1:nqp,1:nqp,2)) )
  dqs_CRmat(1:nsp,1:nqp,2) = &
      matmul( interp_s2q(1:nsp,1:nqp,2) , transpose(dqs_mat(1:nqp,1:nqp,2)) )
  !
  call output2d(array_debugging,dqr_CRmat(1:nsp,1:nqp,2),2,2,2,'dqr_CRmat')
  call output2d(array_debugging,dqs_CRmat(1:nsp,1:nqp,2),2,2,2,'dqs_CRmat')
  !
  ! Inverse of the mass matrix
  !
  invM_mat(1:nsp,1:nsp,2) = matmul( vand_quad , transpose(vand_quad) )
  !
  call output2d(array_debugging,invM_mat(1:nsp,1:nsp,2),2,2,2,'invM_mat')
  !
 !write (121,5) "Solution point weights 1D"
 !write (121,4) (solwts_1d(i),i=1,n_sp1de)
 !write (121,5) "Quadrature weights 1D"
 !write (121,4) (qdtwts_1d(i),i=1,n_sp1de)
 !write (121,5) "Inverse solution point weights 1D"
 !write (121,4) (one/solwts_1d(i),i=1,n_sp1de)
 !write (121,5) "Inverse quadrature weights 1D"
 !write (121,4) (one/qdtwts_1d(i),i=1,n_sp1de)
 !4 format (/,*(f12.6))
 !5 format (/,a)
 !!
 !write (104,*) "Dqr matrix"
 !write (105,*) "Dqs matrix"
 !!
 !write (120,*) "weights"
 !do j = 1,n_sp1de
 !  write (120,3) (weights2d( (n_sp1de-j)*n_sp1de+i , 2 ),i=1,n_sp1de)
 !end do
 !write (120,*) "Inverse weights"
 !do j = 1,n_sp1de
 !  write (120,3) (one/weights2d( (n_sp1de-j)*n_sp1de+i , 2 ),i=1,n_sp1de)
 !end do
 !write (120,*) "Quadrature weights"
 !do j = 1,n_qp1de
 !  write (120,3) (qdtwts_2d( (n_qp1de-j)*n_qp1de+i ),i=1,n_qp1de)
 !end do
 !write (120,*) "Inverse mass matrix"
 !do l = 1,n_sp1de
 !  do k = 1,n_sp1de
 !    write (104,2) k,l
 !    write (105,2) k,l
 !    write (120,2) k,l
 !    do j = 1,n_qp1de
 !      write (104,3) (dqr_mat( (n_qp1de-j)*n_qp1de+i , (l-1)*n_sp1de+k , 2 ),i=1,n_qp1de)
 !      write (105,3) (dqs_mat( (n_qp1de-j)*n_qp1de+i , (l-1)*n_sp1de+k , 2 ),i=1,n_qp1de)
 !    end do
 !    do j = 1,n_sp1de
 !      write (120,3) (invM_mat( (n_sp1de-j)*n_sp1de+i , (l-1)*n_sp1de+k , 2 ),i=1,n_sp1de)
 !    end do
 !  end do
 !end do
 !2 format (/,"i = ",i2,", j = ",i2)
 !3 format (*(es12.4))
 !flush(120)
 !!
 !write (111,*) "Lagrange on Quadrature Points"
 !do l = 1,4
 !  do k = 1,n_sp1de
 !    write (111,4) l,k
 !    do j = 1,n_qp1de
 !      write (111,3) (interp_s2g( (n_qp1de-j)*n_qp1de+i , (l-1)*n_sp1de+k ),i=1,n_qp1de)
 !    end do
 !  end do
 !end do
 !4 format (/,"Face = ",i2,", i = ",i2)
 !write (102,*) "Quadrature Interpolation Vandermonde Matrix"
 !do lj = 0,n_order
 !  do li = 0,n_order
 !    l = n_sp1de*lj + li + 1
 !    write (102,*)
 !    write (102,11) "Basis function for term (x^",li,")*(y^",lj,")"
 !    do kj = n_qp1de,1,-1
 !      write (102,12) (vand_qdtr(n_qp1de*(kj-1)+ki,l),ki=1,n_qp1de)
 !    end do
 !  end do
 !end do
 !!
 !write (103,*) "Quadrature Interpolation Matrix"
 !do lj = 0,n_order
 !  do li = 0,n_order
 !    l = n_sp1de*lj + li + 1
 !    write (103,*)
 !    write (103,11) "Basis function for term (x^",li,")*(y^",lj,")"
 !    do kj = n_qp1de,1,-1
 !      write (103,12) (interp_s2q(n_qp1de*(kj-1)+ki,l),ki=1,n_qp1de)
 !    end do
 !  end do
 !end do
 !!
 !if (mypnum == 0) call write_lagrange_tecplot
 !!
 !10 format (100es14.6)
  11 format (a,i0,a,i0,a)
  12 format (*(es12.4))
  !
  1 format (/," ************************************************",/, &
              "     USING OVER-INTEGRATION FOR DE-ALIASING",/, &
              " ************************************************",/)
 !!
end subroutine old_interp2quadrature
!
!###############################################################################
!
end module startup_mod
