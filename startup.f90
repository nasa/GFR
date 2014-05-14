module startup_mod
  !
  ! Module that contains various functions and routines
  ! used to set up the computational problem
  !
  use module_kind_types
  use polynomial_mod, only : JacobiP => normalized_jacobi_poly
  use geovar,         only : vand_1d,   vand_1d_inv, &
                             vand_tria, vand_tria_inv, &
                             vand_quad, vand_quad_inv
  !
  implicit none
  !
  private
  !
 !include "kind_defs.h"
 !include "constants.h"
  !
  public :: get_quadrature_data
  public :: Vandermonde1D
 !public :: Vandermonde2D_tria
 !public :: Vandermonde2D_quad
 !public :: GradVandermonde1D
 !public :: GradVandermonde2D
 !public :: Dmatrix1D
 !public :: Dmatrix2D_tria
 !public :: Dmatrix2D_weak
 !public :: Dmatrix2D_quad
  public :: invert_matrix
  public :: determinant
 !public :: CellVolume1D
 !public :: CorrectionMatrix1D
 !public :: CorrectionMatrix2D_tria
 !public :: CorrectionMatrix2D_quad
 !public :: rstoab
  !
  real(vbp), save, allocatable, dimension(:,:)   :: grad_vand_1d
  real(vbp), save, allocatable, dimension(:,:,:) :: grad_vand_tria
  !
  logical, parameter :: array_debugging = fals
 !logical, parameter :: array_debugging = true
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
  use geovar, only : nr,n_sp1de,n_sp2dt,n_sp2dq,solpts_1d
  use geovar, only : rs,flxpts_1d,n_order,interp_s2f
  use metvar, only : cormat,d1_mat,dr_mat,ds_mat,dr_hat,ds_hat,cormat_dg
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,n,nf,ierr
  integer :: li,lj,ki,kj
  real(vbp) :: time0,time1
  !
  !.. Local Allocatable Arrays ..
  real(vbp), allocatable, dimension(:,:) :: adjd_rs
  !
continue
  !
  ! Allocation of arrays local to the module
  !
  allocate ( d1_mat(1:n_sp1de,1:n_sp1de) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","d1_mat",1,ierr)
  allocate ( adjd_rs(1:size(rs,1),1:size(rs,2)) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","adjd_rs",1,ierr)
  !
  ! Adjust r and s to be on the interval [-1,1]. The routines
  ! Vandermonde2d_tria, GradVandermonde2D, and CorrectionMatrix2D_tria
  ! need to be adjusted to work for r and s on the interval [0,1].
  !
  adjd_rs(:,:) = two*rs(:,:) - one
  !
  ! Allocate and compute the various
  ! Vandermonde matrices and their inverses
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
  allocate ( grad_vand_1d(1:n_sp1de,1:n_sp1de) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","grad_vand_1d",1,ierr)
  allocate ( grad_vand_tria(1:n_sp2dt,1:n_sp2dt,1:nr) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","grad_vand_tria",1,ierr)
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
  call output2d(vand_quad(:,:),2,2,0,'vand_quad')
  call output2d(vand_quad_inv(:,:),2,2,0,'vand_quad_inv')
  !
  ! Get the derivative of the 1D Vandermonde matrix
  !
  grad_vand_1d = GradVandermonde1D(solpts_1d)
  !
  ! Get the derivative of the Vandermonde matrix and store in return_value
  !
  grad_vand_tria = GradVandermonde2D(adjd_rs)
  !
  ! Get the 1D derivative matrix
  !
  d1_mat = matmul( grad_vand_1d , vand_1d_inv )
  !
  ! Allocation of arrays from the metvar module
  !
  allocate ( dr_mat(1:n_sp2dq,1:n_sp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","dr_mat",1,ierr)
  allocate ( ds_mat(1:n_sp2dq,1:n_sp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","ds_mat",1,ierr)
  allocate ( dr_hat(1:n_sp2dq,1:n_sp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","dr_hat",1,ierr)
  allocate ( ds_hat(1:n_sp2dq,1:n_sp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","ds_hat",1,ierr)
  allocate ( cormat(1:n_sp2dq,1:4*n_sp1de,1:2) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","cormat",1,ierr)
  !
  cormat(1:n_sp2dt,1:3*n_sp1de,1) = CorrectionMatrix2D_tria(adjd_rs)
  cormat(1:n_sp2dq,1:4*n_sp1de,2) = CorrectionMatrix2D_quad(solpts_1d)
  !
  call Dmatrix2d_tria
  !
  call Dmatrix2d_quad(d1_mat,dr_mat(:,:,2),ds_mat(:,:,2))
  !
  call Dmatrix2d_weak
  !
  ! Get the Lagrange polynomial to project the variables/fluxes from the
  ! solution points in the cell interior to the flux points on the cell face
  !
  allocate ( interp_s2f(1:n_sp2dq,1:4*n_sp1de,1:2) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","interp_s2f",1,ierr)
  !
  interp_s2f(1:n_sp2dt,1:3*n_sp1de,1) = LagrangePoly2D_tria(rs)
  interp_s2f(1:n_sp2dq,1:4*n_sp1de,2) = LagrangePoly2D_quad(flxpts_1d,solpts_1d)
  !
  call output2d(interp_s2f(1:n_sp2dq,1:4*n_sp1de,2),2,1,2,'interp_s2f')
  !
  ! Output the solution and flux points to a tecplot file for visualization
  !
 !call write_quadrature
  !
  deallocate ( adjd_rs , stat=ierr )
  call alloc_error("get_quadrature_data","adjd_rs",2,ierr)
  !
  ! Get the correction matrix for DG methods
  !
  call CorrectionMatrix2D_dg
  !
  ! Get the filter matrices
  !
  call Filter2D
  !
  ! Get the matrices to switch between the Gauss and Lobatto nodes
  !
  call gauss_lobatto_matrices
  !
  call interp2quadrature
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
pure function Vandermonde1D(xi) result(return_value)
  !
  real(vbp), dimension(:), intent(in) :: xi
  !
  real(vbp), dimension(1:size(xi),1:size(xi)) :: return_value
  !
  integer :: i,j
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
end function Vandermonde1D
!
!###############################################################################
!
pure function Vandermonde2D_tria(rs) result(return_value)
  !
  real(vbp), dimension(:,:), intent(in) :: rs
  !
  real(vbp), dimension(1:size(rs,2),1:size(rs,2)) :: return_value
  !
  !.. Local Scalars ..
  integer   :: i,j,k,n,np,m
  real(vbp) :: h1,h2,alpha,beta
  real(vbp) :: denom,dennz,denrt
  !
  !.. Local Arrays ..
  real(vbp), dimension(1:2) :: ab
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
        alpha = real(2*i+1,vbp)
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
end function Vandermonde2D_tria
!
!###############################################################################
!
pure function Vandermonde2D_quad(xi) result(return_value)
  !
  real(vbp), dimension(:), intent(in) :: xi
  !
  real(vbp), dimension(1:size(xi)**2,1:size(xi)**2) :: return_value
  !
  !.. Local Scalars ..
  integer   :: k,ki,kj,l,li,lj
  integer   :: n_order,n_sp1de,n_sp2dq
  real(vbp) :: x,y,hx,hy
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
  ! This isn't set up correctly so just make
  ! it an identity matrix for the time being
  !
 !return_value(:,:) = zero
 !do i = 1,np**2
 !  return_value(i,i) = one
 !end do
  !
end function Vandermonde2D_quad
!
!###############################################################################
!
pure function GradVandermonde2D(rs) result(return_value)
  !
  real(vbp), dimension(:,:), intent(in) :: rs
  !
  real(vbp), dimension(1:size(rs,2),1:size(rs,2),1:size(rs,1)) :: return_value
  !
  !.. Local Scalars ..
  integer   :: i,im1,j,k,m,n,np
  real(vbp) :: ri,twoiph,alpha,alphb,beta
  real(vbp) :: fa,gb,dfa,dgb,h1mb,dmdr,dmds
  !
  !.. Local Arrays ..
  real(vbp), dimension(1:2) :: ab
  !
  !.. Constants ..
  real(vbp), parameter :: sqrt2 = sqrt(two)
  !
continue
  !
  np = size(rs,2)
  n  = np2n(np)
  !
  m = 1
  do i = 0,n
    !
    ri   = real(i,vbp)
    im1  = max(i-1,0)
    !
    twoiph = sqrt2 * two**i
    !
    alpha = zero
    alphb = real(2*i+1,vbp)
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
end function GradVandermonde2D
!
!###############################################################################
!
subroutine Dmatrix2D_tria()
  !
  !.. Use Statements ..
  use geovar, only : n_sp2dt
  use metvar, only : dr_mat,ds_mat
  !
  !.. Local Scalars ..
  integer :: i
  !
continue
  !
  ! Initialize the triangle parts of the derivative matrices to zero
  !
  dr_mat(:,:,1) = zero
  ds_mat(:,:,1) = zero
  !
  ! Compute the derivative matrix
  !
  dr_mat(1:n_sp2dt,1:n_sp2dt,1) = two*matmul( grad_vand_tria(:,:,1) , vand_tria_inv )
  ds_mat(1:n_sp2dt,1:n_sp2dt,1) = two*matmul( grad_vand_tria(:,:,2) , vand_tria_inv )
  !
  ! Transpose these matrices to get them in the right ordering
  !
  dr_mat(1:n_sp2dt,1:n_sp2dt,1) = transpose(dr_mat(1:n_sp2dt,1:n_sp2dt,1))
  ds_mat(1:n_sp2dt,1:n_sp2dt,1) = transpose(ds_mat(1:n_sp2dt,1:n_sp2dt,1))
  !
end subroutine Dmatrix2D_tria
!
!###############################################################################
!
subroutine Dmatrix2D_quad(d1_mat,dr_mat,ds_mat)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:),    intent(in) :: d1_mat
  real(wp), dimension(:,:), intent(inout) :: dr_mat
  real(wp), dimension(:,:), intent(inout) :: ds_mat
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,n
  integer :: lbeg,lend
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
end subroutine Dmatrix2D_quad
!
!###############################################################################
!
subroutine Dmatrix2d_weak()
  !
  ! Create the weak form of the derivative matrices
  !
  !.. Use Statements ..
  use geovar, only : nr,n_sp1de,solwts_1d
  use metvar, only : d1_mat,dr_hat,ds_hat
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,n,lbeg,lend
  !
  !.. Local Arrays ..
  real(vbp), dimension(1:n_sp1de,1:n_sp1de) :: d1_hat
  !
continue
  !
  ! Initialize the weak form of the derivative matrices to zero
  !
  dr_hat(:,:,:) = zero
  ds_hat(:,:,:) = zero
  !
  ! Create the weak form of the 1D derivative matrix
  !
  do n = 1,n_sp1de
    do j = 1,n_sp1de
      d1_hat(j,n) = -d1_mat(n,j)*solwts_1d(n)/solwts_1d(j)
    end do
  end do
  !
  ! Create the weak form of the quadrilateral r/xi derivative matrix
  !
  k = 0
  do j = 1,n_sp1de
    lbeg = (j-1)*n_sp1de + 1
    lend = j*n_sp1de
    do i = 1,n_sp1de
      k = k + 1
      dr_hat(lbeg:lend,k,2) = d1_hat(i,:)
    end do
  end do
  !
  ! Create the weak form of the quadrilateral s/eta derivative matrix
  !
  k = 0
  do j = 1,n_sp1de
    do i = 1,n_sp1de
      k = k + 1
      do l = 1,n_sp1de
        ds_hat(i+(l-1)*n_sp1de,k,2) = d1_hat(j,l)
      end do
    end do
  end do
 !!
 !! Create the weak form of the triangle derivative matrices
 !!
 !write (111,*) "1D Dr Weak"
 !write (112,*) "1D Ds Weak"
 !do j = 1,n_sp1de
 !  write (111,3) (d1_hat(i,n_sp1de+1-j),i=1,n_sp1de)
 !  write (112,3) (d1s_hat(i,n_sp1de+1-j),i=1,n_sp1de)
 !end do
 !!
 !write (115,*) "1D Derivative Matrix"
 !do j = 1,n_sp1de
 !  write (115,4) (d1_mat(i,n_sp1de+1-j),i=1,n_sp1de)
 !end do
 !write (115,*)
 !write (115,*)
 !do j = 1,n_sp1de
 !  write (115,4) (d1r_mat(i,n_sp1de+1-j),i=1,n_sp1de)
 !end do
 !write (113,*) "2D Dr Weak"
 !write (114,*) "2D Ds Weak"
 !do l = 1,n_sp1de
 !  do k = 1,n_sp1de
 !    write (113,2) k,l
 !    write (114,2) k,l
 !    do j = 1,n_sp1de
 !      write (113,3) (dr_hat( (l-1)*n_sp1de+k , (n_sp1de-j)*n_sp1de+i , 2 ),i=1,n_sp1de)
 !      write (114,3) (ds_hat( (l-1)*n_sp1de+k , (n_sp1de-j)*n_sp1de+i , 2 ),i=1,n_sp1de)
 !    end do
 !  end do
 !end do
 !!
 !2 format (/,"i = ",i2,", j = ",i2)
 !3 format (<n_sp1de>es12.4)
 !4 format (<n_sp1de>es22.14)
  !
end subroutine Dmatrix2d_weak
!
!###############################################################################
!
pure function Dmatrix2D_tria_weak(rs,V) result(return_value)
  !
  real(vbp), dimension(:,:), intent(in) :: rs
  real(vbp), dimension(:,:), intent(in) :: V
  !
  real(vbp), dimension(1:size(V,1),1:size(V,2),1:size(rs,1)) :: return_value
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Arrays ..
  real(vbp), dimension(1:size(V,1),1:size(V,2)) :: invVVT
  !
continue
  !
  ! Precompute the inverse of V*V_Transpose
  !
  invVVT = transpose(V)
  invVVT = matmul(V,invVVT)
  invVVT = invert_matrix(invVVT)
  !
  ! Get the derivative of the Vandermonde matrix and store in return_value
  !
  return_value = GradVandermonde2D(rs)
  !
  ! Compute the weak form of the derivative matrix
  !
  do i = 1,size(rs,1)
    return_value(:,:,i) = transpose(return_value(:,:,i))
    return_value(:,:,i) = matmul( V , return_value(:,:,i) )
    return_value(:,:,i) = matmul( return_value(:,:,i) , invVVT )
  end do
  !
end function Dmatrix2D_tria_weak
!
!###############################################################################
!
pure function GradJacobiP(n,alpha,beta,x) result(return_value)
  !
  integer,   intent(in) :: n
  real(vbp), intent(in) :: alpha
  real(vbp), intent(in) :: beta
  real(vbp), intent(in) :: x
  !
  real(vbp) :: return_value
  !
  !.. Local Scalars ..
  real(vbp) :: rn,a1,b1,ab1
  !
continue
  !
  if (n == 0) then
    return_value = zero
  else
    !
    rn  = real(n,vbp)
    a1  = alpha + one
    b1  = beta  + one
    ab1 = alpha + b1
    !
    return_value = sqrt(rn*(rn+ab1)) * JacobiP(n-1,a1,b1,x)
    !
  end if
  !
end function GradJacobiP
!
!###############################################################################
!
pure function GradVandermonde1D(xi) result(return_value)
  !
  real(vbp), dimension(:), intent(in) :: xi
  !
  real(vbp), dimension(1:size(xi),1:size(xi)) :: return_value
  !
  integer   :: i,j
  real(vbp) :: rj
  !
continue
  !
  return_value(:,:) = zero
  !
  do j = 2,size(xi)
    do i = 1,size(xi)
      rj = real(j-1,vbp)
      return_value(i,j) = sqrt(rj*(rj+one)) * JacobiP(j-2,one,one,xi(i))
    end do
  end do
  !
end function GradVandermonde1D
!
!###############################################################################
!
pure function Dmatrix1D(V,Vr) result(return_value)
  !
  real(vbp), dimension(:,:), intent(in) :: V
  real(vbp), dimension(:,:), intent(in) :: Vr
  !
  real(vbp), dimension(1:size(V,1),1:size(V,2)) :: return_value
  !
continue
  !
  return_value = invert_matrix(V)
  return_value = matmul(Vr,return_value)
  !
end function Dmatrix1D
!
!###############################################################################
!
pure function CellVolume1D(nk,np,x,Dr) result(return_value)
  !
  integer, intent(in) :: nk
  integer, intent(in) :: np
  !
  real(vbp), dimension(1:np,1:nk), intent(in) :: x
  real(vbp), dimension(1:np,1:np), intent(in) :: Dr
  !
  real(vbp), dimension(1:nk) :: return_value
  !
  integer :: i
  !
continue
  !
  do i = 1,nk
    return_value(i) = dot_product(Dr(1,:),x(:,i))
  end do
  !
end function CellVolume1D
!
!###############################################################################
!
!pure function CorrectionMatrix1D(V) result(return_value)
!  !
!  real(vbp), dimension(:,:), intent(in) :: V
!  !
!  real(vbp), dimension(size(V,1),2) :: return_value
!  !
!  real(vbp), dimension(size(V,1),2) :: Emat
!  !
!continue
!  !
!  Emat = zero
!  Emat(1,1) = one
!  Emat(size(V,1),2) = one
!  !
!  return_value = matmul(V,matmul(transpose(V),Emat))
!  !
!end function CorrectionMatrix1D
!
!###############################################################################
!
pure function CorrectionMatrix2D_tria(rs) result(return_value)
  !
  real(vbp), dimension(:,:), intent(in) :: rs
  !
  real(vbp), dimension(1:size(rs,2),1:3*np2nfp(size(rs,2))) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,n,np,nfp
  !
  !.. Local Arrays ..
  integer,   dimension(1:np2nfp(size(rs,2)),1:3) :: Fmask
  real(vbp), dimension(1:np2nfp(size(rs,2)),1:np2nfp(size(rs,2))) :: V1D
  !
continue
  !
  return_value = zero
  !
  np  = size(rs,2)
  n   = np2n(np)
  nfp = np2nfp(np)
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
  return_value = matmul( transpose(vand_tria) , return_value )
  return_value = two*matmul( vand_tria , return_value )
  !
end function CorrectionMatrix2D_tria
!
!###############################################################################
!
pure function CorrectionMatrix2D_quad(xi) result(return_value)
  !
  !.. Formal Arguments ..
  real(vbp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(vbp), dimension(1:size(xi)**2,1:4*size(xi)) :: return_value
  !
  !.. Local Scalars ..
  integer   :: i,j,nfp
  !
  !.. Local Arrays ..
  integer,   dimension(1:size(xi),1:size(xi)) :: Fmask
  real(vbp), dimension(1:size(xi)) :: V1D
 !real(vbp), dimension(1:size(xi),1:size(xi)) :: V1D
  !
continue
  !
  nfp = size(xi)
  !
  ! Initialize the correction matrix to zero
  !
  return_value(:,:) = zero
  !
  ! 1D Vandermonde matrix
  !
 !do i = 1,nfp
 !  V1D(i) = dot_product( vand_1D(1,1:nfp) , vand_1D(i,1:nfp) )
 !end do
 !V1D = matmul( vand_1D , transpose(vand_1d) )
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
 !    return_value(Fmask(i,j),j) = V1D(i,1)
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
 !    return_value(Fmask(i,j),j+nfp) = V1D(i,1)
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
 !    return_value(Fmask(i,j),j+2*nfp) = V1D(i,1)
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
 !    return_value(Fmask(i,j),j+3*nfp) = V1D(i,1)
      return_value(Fmask(i,j),j+3*nfp) = V1D(i)
    end do
  end do
  !
end function CorrectionMatrix2D_quad
!
!###############################################################################
!
subroutine CorrectionMatrix2D_dg()
  !
  !.. Use Statements ..
  use geovar, only : n_sp1de,n_sp2dq,solwts_1d,solpts_1d,interp_s2f
  use metvar, only : cormat_dg
  !
  !.. Local Scalars ..
  integer   :: i,j,k,nf,ierr
  !
continue
  !
  ! Allocate the correction matrix for the DG method
  !
  allocate ( cormat_dg(1:n_sp2dq,1:4*n_sp1de,1:2) , source=zero , stat=ierr )
  call alloc_error("CorrectionMatrix2D_dg","cormat_dg",1,ierr)
  !
  ! Face 1 - bottom/south edge of quad
  !
  nf = 1
  do i = 1,n_sp1de
    do j = 1,n_sp1de
      k = (j-1)*n_sp1de + i
      cormat_dg(k,nf,2) = eval_LagrangePoly1D(j,-one,solpts_1d)/solwts_1d(j)
    end do
    nf = nf + 1
  end do
  !
  ! Face 2 - right/east edge of quad
  !
  nf = n_sp1de+1
  do j = 1,n_sp1de
    do i = n_sp1de,1,-1
      k = (j-1)*n_sp1de + i
      cormat_dg(k,nf,2) = eval_LagrangePoly1D(i,one,solpts_1d)/solwts_1d(i)
    end do
    nf = nf + 1
  end do
  !
  ! Face 3 - top/north edge of quad
  !
  nf = 2*n_sp1de+1
  do i = n_sp1de,1,-1
    do j = n_sp1de,1,-1
      k = (j-1)*n_sp1de + i
      cormat_dg(k,nf,2) = eval_LagrangePoly1D(j,one,solpts_1d)/solwts_1d(j)
    end do
    nf = nf + 1
  end do
  !
  ! Face 4 - left/west edge of quad
  !
  nf = 3*n_sp1de+1
  do j = n_sp1de,1,-1
    do i = 1,n_sp1de
      k = (j-1)*n_sp1de + i
      cormat_dg(k,nf,2) = eval_LagrangePoly1D(i,-one,solpts_1d)/solwts_1d(i)
    end do
    nf = nf + 1
  end do
  !
 !do nf = 3*n_sp1de+1,4*n_sp1de
 !  !
 !  j = 
 !  do j = 1,n_sp1de
 !    do i = 1,n_sp1de
 !      k = (j-1)*n_sp1de + i
 !     !cormat_dg(k,nf,2) = interp_s2f(k,nf,2)/solwts_1d(i)
 !      cormat_dg(k,nf,2) = eval_LagrangePoly1D(i,-one,solpts_1d)/solwts_1d(i)
 !    end do
 !  end do
 !  !
 !end do
 !!
 !write (98,*) "DGSEM Correction Function"
 !write (99,*) "interp_s2f interpolation polynomial"
 !do nf = 1,4*n_sp1de
 !  !
 !  if (nf <= 1*n_sp1de) then
 !    write (98,2) nf
 !    write (99,2) nf
 !  else if (nf <= 2*n_sp1de) then
 !    write (98,3) nf - n_sp1de
 !    write (99,3) nf - n_sp1de
 !  else if (nf <= 3*n_sp1de) then
 !    write (98,4) n_sp1de+1-(nf-2*n_sp1de)
 !    write (99,4) n_sp1de+1-(nf-2*n_sp1de)
 !  else
 !    write (98,5) n_sp1de+1-(nf-3*n_sp1de)
 !    write (99,5) n_sp1de+1-(nf-3*n_sp1de)
 !  end if
 !  !
 !  do j = n_sp1de,1,-1
 !    write (98,1) (cormat_dg((j-1)*n_sp1de+i,nf,2),i=1,n_sp1de)
 !    write (99,1) (interp_s2f((j-1)*n_sp1de+i,nf,2),i=1,n_sp1de)
 !  end do
 !  !
 !end do
 !2 format (/," Bottom/South edge for xi = ",i0)
 !3 format (/," Right/East edge for eta = ",i0)
 !4 format (/," Top/North edge for xi = ",i0)
 !5 format (/," Left/West edge for eta = ",i0)
 !!
 !write (97,1) (solwts_1d(i),i=1,n_sp1de)
 !1 format (<n_sp1de>es12.4)
 !!
 !write (*,*)
 !do i = 1,n_sp1de
 !  write (*,*) eval_LagrangePoly1D(i,-one,solpts_1d), &
 !              eval_LagrangePoly1D(i,one,solpts_1d)
 !end do
 !write (*,*)
  !
end subroutine CorrectionMatrix2D_dg
!
!###############################################################################
!
pure function LagrangePoly1D(xi) result(return_value)
  !
  !.. Formal Arguments ..
  real(vbp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(vbp), dimension(1:size(xi)) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
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
end function LagrangePoly1D
!
!###############################################################################
!
pure function LagrangePoly2D_quad(xi_fp,xi_sp) result(return_value)
  !
  !.. Formal Arguments ..
  real(vbp), dimension(:), intent(in) :: xi_fp
  real(vbp), dimension(:), intent(in) :: xi_sp
  !
  !.. Function Return Value ..
  real(vbp), dimension(1:size(xi_sp)**2,1:4*size(xi_fp)) :: return_value
  !
  !.. Local Scalars ..
  integer   :: i,j,n,nf,nfp,nsp
  real(vbp) :: x,y
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
end function LagrangePoly2D_quad
!
!###############################################################################
!
pure function LagrangePoly2D_tria(rs) result(return_value)
  !
  !.. Formal Arguments ..
  real(vbp), dimension(:,:), intent(in) :: rs
  !
  !.. Function Return Argument ..
  real(vbp), dimension(1:size(rs,2),1:3*np2nfp(size(rs,2))) :: return_value
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
end function LagrangePoly2D_tria
!
!###############################################################################
!
pure function eval_LagrangePoly1D(k,x,xi) result(return_value)
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
  real(vbp),               intent(in) :: x
  real(vbp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(vbp) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
  !
continue
  !
  nfp = size(xi)
  imask(:) = (/ (i,i=1,nfp) /)
  !
  return_value = product( x    -xi(1:nfp) , imask /= k ) / &
                 product( xi(k)-xi(1:nfp) , imask /= k )
  !
end function eval_LagrangePoly1D
!
!###############################################################################
!
pure function eval_LagrangePoly2D(k,l,x,y,xi,eta) result(return_value)
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
  real(vbp),               intent(in) :: x
  real(vbp),               intent(in) :: y
  real(vbp), dimension(:), intent(in) :: xi
  real(vbp), dimension(:), intent(in) :: eta
  !
  !.. Function Return Value ..
  real(vbp) :: return_value
  !
continue
  !
  return_value = eval_LagrangePoly1D(k,x,xi) * eval_LagrangePoly1D(l,y,eta)
  !
end function eval_LagrangePoly2D
!
!###############################################################################
!
pure function dlagrange(m,x,xi) result(return_value)
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
  real(vbp),               intent(in) :: x
  real(vbp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(vbp) :: return_value
  !
  !.. Local Scalars ..
  integer :: l,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
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
end function dlagrange
!
!###############################################################################
!
pure function d2lagrange(m,x,xi) result(return_value)
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
  real(vbp),               intent(in) :: x
  real(vbp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(vbp) :: return_value
  !
  !.. Local Scalars ..
  integer :: l,n,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:size(xi)) :: imask
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
end function d2lagrange
!
!###############################################################################
!
pure function CorrectionMatrix1D(xi) result(return_value)
  !
  !.. Formal Arguments ..
  real(vbp), dimension(:), intent(in) :: xi
  !
  !.. Function Return Value ..
  real(vbp), dimension(1:size(xi)) :: return_value
  !
  !.. Local Scalars ..
  integer   :: i,j,nfp
  real(vbp) :: c1,c2
  !
  !.. Local Arrays ..
  real(vbp), dimension(0:size(xi))              :: cdRR
  real(vbp), dimension(-1:size(xi),-1:size(xi)) :: coefLe
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
      c1 = real(2*i-1,vbp)/real(i,vbp)
      c2 = real(  i-1,vbp)/real(i,vbp)
      !
      coefLe(i,j) = c1*coefLe(i-1,j-1) - c2*coefLe(i-2,j)
      !
    end do
  end do
  !
  ! Get the coefficients for the derivative of
  ! the right Radau polynomial of degree nfp
  !
  c1 = half * real((-1)**nfp,vbp)
  do i = 1,nfp
    cdRR(i-1) = c1 * real(i,vbp) * (coefLe(nfp,i) - coefLe(nfp-1,i))
  end do
  !
  ! Finally, get the value for the derivative of the
  ! right Radau polynomial at each 1D solution point
  !
  do i = 1,nfp
    return_value(i) = -sum( (xi(i)**(/0:nfp-1/)) * cdRR(0:nfp-1) )
  end do
  !
end function CorrectionMatrix1D
!
!###############################################################################
!
subroutine Filter2D()
  !
  ! Routine to compute the spectral filter matrices
  !
  !.. Use Statements ..
  use geovar, only : n_order,n_sp1de,n_sp2dt,n_sp2dq,filter
  use ovar,   only : Filter_Option,Filter_Order, &
                     Filter_Cutoff,Filter_Alpha
  !
  !.. Local Scalars ..
  integer   :: i,j,k,ierr
  real(vbp) :: rd,rn
  !
  !.. Local Arrays ..
  real(vbp), dimension(1:n_sp1de)           :: filter_1D
  real(vbp), dimension(1:n_sp2dt,1:n_sp2dt) :: filter_tria
  real(vbp), dimension(1:n_sp2dq,1:n_sp2dq) :: filter_quad
  !
continue
  !
  ! Allocate the filter array
  !
  allocate ( filter(1:n_sp2dq,1:n_sp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("Filter2D","filter",1,ierr)
  !
  ! Initialize the local filter matrices to identity matrices
  !
  filter_1D(:) = one
  !
  filter_tria(:,:) = zero
  do i = 1,n_sp2dt
    filter_tria(i,i) = one
  end do
  !
  filter_quad(:,:) = zero
  do i = 1,n_sp2dq
    filter_quad(i,i) = one
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
        rn = real(i+j-Filter_Cutoff,vbp)
        rd = real(n_order-Filter_Cutoff,vbp)
        !
        filter_tria(k,k) = exp(-Filter_Alpha * (rn/rd)**Filter_Order)
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
  filter_tria = matmul( vand_tria , filter_tria )
  filter_tria = matmul( filter_tria , vand_tria_inv )
  !
  ! Compute the filter diagonal for quadrilaterals
  !
  if (Filter_Option == 1) then
    !
    ! Option 1 - Equivalent to using a 1D filter in each direction
    !
    ! First get the 1D filter diagonal
    !
    do i = 0,n_order
      if (i >= Filter_Cutoff) then
        !
        rn = real(i-Filter_Cutoff,vbp)
        rd = real(n_order-Filter_Cutoff,vbp)
        !
        filter_1D(i+1) = exp(-Filter_Alpha * (rn/rd)**Filter_Order)
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
        filter_quad(k,k) = filter_1D(i+1)*filter_1D(j+1)
        !
      end do
    end do
    !
   !1 format (<n_sp1de>es22.12)
   !write (103,*) "Quad filter option 1"
   !do j = 0,n_order
   !  write (103,1) (filter_quad(n_sp1de*j+i+1,n_sp1de*j+i+1),i=0,n_order)
   !end do
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
          rn = real(i+j-Filter_Cutoff,vbp)
          rd = real(2*n_order-Filter_Cutoff,vbp)
          !
          filter_quad(k,k) = exp(-Filter_Alpha * (rn/rd)**Filter_Order)
          !
        end if
        !
      end do
    end do
   !write (104,*) "Quad filter option 2"
   !do j = 0,n_order
   !  write (104,1) (filter_quad(n_sp1de*j+i+1,n_sp1de*j+i+1),i=0,n_order)
   !end do
    !
  end if
 !stop
  !
  ! Compute the final filter matrix for quadrilaterals
  !
  filter_quad = matmul( vand_quad , filter_quad )
  filter_quad = matmul( filter_quad , vand_quad_inv )
  !
  ! Finally, store the triangle and quad filters in the common filter array
  !
  filter(1:n_sp2dt,1:n_sp2dt,1) = filter_tria
  filter(1:n_sp2dq,1:n_sp2dq,2) = filter_quad
  !
 !1 format (<n_sp2dq>es14.6)
 !do j = 1,n_sp2dt
 !  write (103,1) (filter(i,j,1),i=1,n_sp2dt)
 !end do
 !do j = 1,n_sp2dq
 !  write (104,1) (filter(i,j,2),i=1,n_sp2dq)
 !end do
  !
end subroutine Filter2D
!
!###############################################################################
!
!pure function FaceToCellConnectivity1D(nface) result(return_value)
!  !
!  integer, intent(in) :: nface
!  !
!  integer, dimension(2,nface) :: return_value
!  !
!  integer :: i
!  !
!continue
!  !
!  do i = 1,nface
!    return_value(1,i) = i-1
!    return_value(2,i) = i
!  end do
!  !
!end function FaceToCellConnectivity1D
!!
!!###############################################################################
!!
!pure function CellToCellConnectivity1D(ncell) result(return_value)
!  !
!  integer, intent(in) :: ncell
!  !
!  integer, dimension(2,ncell) :: return_value
!  !
!  integer :: i
!  !
!continue
!  !
!  do i = 1,ncell
!    return_value(1,i) = i-1
!    return_value(2,i) = i+1
!  end do
!  !
!end function CellToCellConnectivity1D
!!
!!###############################################################################
!!
!pure function FaceNormals1D(nface) result(return_value)
!  !
!  integer, intent(in) :: nface
!  !
!  real(vbp), dimension(nface) :: return_value
!  !
!  integer :: i
!  !
!continue
!  !
!  do i = 1,nface
!    return_value(i) = one
!  end do
!  !
!end function FaceNormals1D
!
!###############################################################################
!
pure function determinant(a) result(return_value)
  !
  real(vbp), dimension(:,:), intent(in) :: a
  !
  real(vbp) :: return_value
  !
  integer :: i,n
  !
  integer,   dimension(size(a,1))           :: indx
  real(vbp), dimension(size(a,1),size(a,1)) :: b
  !
continue
  !
  n = size(a,1)
  !
  b(:,:) = a(:,:)
  !
  call ludcmp(b,n,indx,return_value)
  !
  do i = 1,n
    return_value = return_value * b(i,i)
  end do
  !
end function determinant
!
!###############################################################################
!
pure function invert_matrix(a) result(return_value)
  !
  real(vbp), dimension(:,:), intent(in) :: a
  !
  real(vbp), dimension(size(a,1),size(a,1)) :: return_value
  !
  integer   :: np,i,ii,j,k,ll,imax
  real(vbp) :: summ
  !
  integer,   dimension(size(a,1))           :: indx
  real(vbp), dimension(size(a,1))           :: temp
  real(vbp), dimension(size(a,1))           :: vv
  real(vbp), dimension(size(a,1),size(a,1)) :: b
  !
continue
  !
  np = size(a,1)
  !
  b(:,:) = a(:,:)
  !
  return_value(:,:) = zero
  do i = 1,np
    return_value(i,i) = one
  end do
  !
  ! Perform the LU decomposition
  !
  vv = one/maxval(abs(b),dim=2)
  !
  do j = 1,np
    !
    imax = (j-1) + sum(maxloc(vv(j:np)*abs(b(j:np,j))))
    if (j /= imax) then
      temp(:) = b(imax,:)
      b(imax,:) = b(j,:)
      b(j,:) = temp(:)
      vv(imax) = vv(j)
    end if
    !
    indx(j) = imax
    !
    if (b(j,j) == zero) b(j,j) = rndoff
    !
    b(j+1:np,j) = b(j+1:np,j)/b(j,j)
    b(j+1:np,j+1:np) = b(j+1:np,j+1:np) &
                     - spread(b(j+1:np,j),dim=2,ncopies=np-j)* &
                       spread(b(j,j+1:np),dim=1,ncopies=np-j)
    !
  end do
  !
  ! Perform back substitution step
  !
  do k = 1,np
    !
    ii = 0
    do i = 1,np
      ll = indx(i)
      summ = return_value(ll,k)
      return_value(ll,k) = return_value(i,k)
      if (ii /= 0) then
        summ = summ - dot_product(b(i,ii:i-1),return_value(ii:i-1,k))
      else if (summ /= zero) then
        ii = i
      end if
      return_value(i,k) = summ
    end do
    !
    do i = np,1,-1
      summ = dot_product(b(i,i+1:np),return_value(i+1:np,k))
      return_value(i,k) = (return_value(i,k) - summ)/b(i,i)
    end do
    !
  end do
  !
end function invert_matrix
!
!###############################################################################
!
pure function invert_matrix90(a) result(return_value)
  !
  real(vbp), dimension(:,:), intent(in) :: a
  !
  real(vbp), dimension(size(a,1),size(a,1)) :: return_value
  !
  integer   :: np,i,k
  real(vbp) :: d
  !
  integer,   dimension(size(a,1))           :: indx
  real(vbp), dimension(size(a,1),size(a,1)) :: b
  !
continue
  !
  np = size(a,1)
  !
  b(:,:) = a(:,:)
  !
  return_value(:,:) = zero
  do i = 1,np
    return_value(i,i) = one
  end do
  !
  ! Perform LU decomposition
  !
  call ludcmp(b,np,indx,d)
  !
  ! Perform back substitution step
  !
  do k = 1,np
    call lubksb(b,np,indx,return_value(1:np,k))
  end do
  !
end function invert_matrix90
!
!###############################################################################
!
pure subroutine ludcmp(a,np,indx,d)
  !
  integer,                            intent(in) :: np
  integer,   dimension(1:np),        intent(out) :: indx
  real(vbp),                         intent(out) :: d
  real(vbp), dimension(1:np,1:np), intent(inout) :: a
  !
  integer :: i,j,k,imax
  !
  real(vbp), dimension(1:np) :: vv,temp
  !
continue
  !
  d = one
  vv = one/maxval(abs(a),dim=2)
  !
  do j = 1,np
    !
    imax = (j-1) + sum(maxloc(vv(j:np)*abs(a(j:np,j))))
    if (j /= imax) then
      temp(:) = a(imax,:)
      a(imax,:) = a(j,:)
      a(j,:) = temp(:)
      d = -d
      vv(imax) = vv(j)
    end if
    !
    indx(j) = imax
    !
    if (a(j,j) == zero) a(j,j) = rndoff
    !
    a(j+1:np,j) = a(j+1:np,j)/a(j,j)
    a(j+1:np,j+1:np) = a(j+1:np,j+1:np) - spread(a(j+1:np,j),dim=2,ncopies=np-j)* &
                                          spread(a(j,j+1:np),dim=1,ncopies=np-j)
    !
  end do
  !
end subroutine ludcmp
!
!###############################################################################
!
pure subroutine lubksb(a,np,indx,y)
  !
  integer,                         intent(in) :: np
  integer,   dimension(1:np),      intent(in) :: indx
  real(vbp), dimension(1:np,1:np), intent(in) :: a
  real(vbp), dimension(1:np),   intent(inout) :: y
  !
  integer :: i,ii,j,ll
  !
  real(vbp) :: summ
  !
continue
  !
  ii = 0
  do i = 1,np
    ll = indx(i)
    summ = y(ll)
    y(ll) = y(i)
    if (ii /= 0) then
      summ = summ - dot_product(a(i,ii:i-1),y(ii:i-1))
    else if (summ /= zero) then
      ii = i
    end if
    y(i) = summ
  end do
  !
  do i = np,1,-1
    y(i) = (y(i) - dot_product(a(i,i+1:np),y(i+1:np)))/a(i,i)
  end do
  !
end subroutine lubksb
!
!###############################################################################
!
pure function rstoab(rs) result(return_value)
  !
  real(vbp), dimension(1:2), intent(in) :: rs
  !
  real(vbp), dimension(1:2) :: return_value
  !
continue
  !
  return_value(1) = two * (one+rs(1)) * reciprocal(one-rs(2)) - one
  return_value(2) = rs(2)
  !
end function rstoab
!
!###############################################################################
!
subroutine gauss_lobatto_matrices()
  !
  !.. Use Statements ..
  use geovar, only : n_order,n_sp1de,n_sp2dt,n_sp2dq
  use geovar, only : gauss_1d,lobatto_1d,radau_1d
  use geovar, only : exact2numer,numer2exact
  use geovar, only : radau2numer,numer2radau
  use geovar, only : gauss2numer,numer2gauss
  use geovar, only : lobatto2numer,numer2lobatto
  use ovar,   only : error_loc,itestcase
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,ierr
  !
  !.. Local Allocatable Arrays ..
  real(vbp), allocatable, dimension(:,:) :: gauss_vand
  real(vbp), allocatable, dimension(:,:) :: gauss_vand_inv
  real(vbp), allocatable, dimension(:,:) :: lobatto_vand
  real(vbp), allocatable, dimension(:,:) :: lobatto_vand_inv
  real(vbp), allocatable, dimension(:,:) :: radau_vand
  real(vbp), allocatable, dimension(:,:) :: radau_vand_inv
  !
  !.. Local Parameters ..
  real(vbp), dimension(1:2,1:3) :: &
        tnodes = reshape([zero,zero,one,zero,zero,one],[2,3])
  real(vbp), dimension(1:2,1:4) :: &
        qnodes = reshape([-one,-one,one,-one,one,one,-one,one],[2,4])
  !
continue
  !
  ! Get the Vandermonde matrix and its inverse for the right Radau nodes
  !
  allocate ( radau_vand(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","radau_vand",1,ierr)
  allocate ( radau_vand_inv(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","radau_vand_inv",1,ierr)
  !
  radau_vand = Vandermonde2d_quad(radau_1d)
  radau_vand_inv = invert_matrix(radau_vand)
  !
  ! Get the Vandermonde matrix and its inverse for the Gauss nodes
  !
  allocate ( gauss_vand(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","gauss_vand",1,ierr)
  allocate ( gauss_vand_inv(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","gauss_vand_inv",1,ierr)
  !
  gauss_vand = Vandermonde2d_quad(gauss_1d)
  gauss_vand_inv = invert_matrix(gauss_vand)
  !
  ! Get the Vandermonde matrix and its inverse for the Lobatto nodes
  !
  allocate ( lobatto_vand(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","lobatto_vand",1,ierr)
  allocate ( lobatto_vand_inv(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","lobatto_vand_inv",1,ierr)
  !
  lobatto_vand = Vandermonde2d_quad(lobatto_1d)
  lobatto_vand_inv = invert_matrix(lobatto_vand)
  !
  ! Now get the conversion matrices through matrix multiplication
  !
 !if (any(itestcase == Advec_Problems)) then
    !
    allocate ( radau2numer(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
    call alloc_error("gauss_lobatto_matrices","radau2numer",1,ierr)
    allocate ( numer2radau(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
    call alloc_error("gauss_lobatto_matrices","numer2radau",1,ierr)
    !
    radau2numer = matmul( vand_quad , radau_vand_inv )
    numer2radau = matmul( radau_vand , vand_quad_inv )
    !
    allocate ( gauss2numer(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
    call alloc_error("gauss_lobatto_matrices","gauss2numer",1,ierr)
    allocate ( numer2gauss(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
    call alloc_error("gauss_lobatto_matrices","numer2gauss",1,ierr)
    !
    gauss2numer = matmul( vand_quad , gauss_vand_inv )
    numer2gauss = matmul( gauss_vand , vand_quad_inv )
    !
    allocate ( lobatto2numer(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
    call alloc_error("gauss_lobatto_matrices","lobatto2numer",1,ierr)
    allocate ( numer2lobatto(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
    call alloc_error("gauss_lobatto_matrices","numer2lobatto",1,ierr)
    !
    lobatto2numer = matmul( vand_quad , lobatto_vand_inv )
    numer2lobatto = matmul( lobatto_vand , vand_quad_inv )
    !
 !else
    !
    allocate ( exact2numer(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
    call alloc_error("gauss_lobatto_matrices","exact2numer",1,ierr)
    allocate ( numer2exact(1:n_sp2dq,1:n_sp2dq) , source=zero , stat=ierr )
    call alloc_error("gauss_lobatto_matrices","numer2exact",1,ierr)
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
  call alloc_error("gauss_lobatto_matrices","gauss_vand",2,ierr)
  deallocate ( gauss_vand_inv , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","gauss_vand_inv",2,ierr)
  deallocate ( lobatto_vand , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","lobatto_vand",2,ierr)
  deallocate ( lobatto_vand_inv , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","lobatto_vand_inv",2,ierr)
  deallocate ( radau_vand , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","radau_vand",2,ierr)
  deallocate ( radau_vand_inv , stat=ierr )
  call alloc_error("gauss_lobatto_matrices","radau_vand_inv",2,ierr)
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
pure function Vandermonde1D_N(xi,N) result(return_value)
  !
  integer,                 intent(in) :: N
  real(vbp), dimension(:), intent(in) :: xi
  !
  real(vbp), dimension(1:size(xi),1:N) :: return_value
  !
  integer :: i,j
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
end function Vandermonde1D_N
!
!###############################################################################
!
pure function Vandermonde2D_quad_N(xi,N) result(return_value)
  !
  integer,                 intent(in) :: N
  real(vbp), dimension(:), intent(in) :: xi
  !
  real(vbp), dimension(1:size(xi)**2,1:(N+1)**2) :: return_value
  !
  !.. Local Scalars ..
  integer   :: k,ki,kj,l,li,lj
  integer   :: nsp,nso,nsq,nqp,nqo,nqq
  real(vbp) :: x,y,hx,hy
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
end function Vandermonde2D_quad_N
!
!###############################################################################
!
subroutine interp2quadrature()
  !
  !.. Use Statements ..
  use geovar, only : n_order,n_sp1de,n_sp2dt,n_sp2dq
  use geovar, only : n_order_qp,n_qp1de,n_qp2dt,n_qp2dq
  use geovar, only : qdtpts_1d,qdtwts_1d,qdtwts_2d
  use geovar, only : solpts_1d,solwts_1d,invM_mat,weights2d
  use geovar, only : interp_f2g,projct_g2f,vand_qdtr
  use geovar, only : interp_s2q,projct_q2s,vand_qdtr_inv
  use geovar, only : interp_s2g,using_quadrature
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
  if (using_quadrature) write (iout,1)
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
    call nodes_legendre_gauss_lobatto(n_order_qp,qdtpts_1d,vp)
    call weights_legendre_gauss_lobatto(n_order_qp,vp,qdtwts_1d)
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
  call output2d(interp_f2g(:,:),2,2,3,'interp_f2g')
  !
  ! Matrix for projecting nodal values at face quadrature
  ! points to the face flux/solution points
  !
  projct_g2f(:,:) = matmul( vand1D_f2g , transpose(vand_1D) )
  do i = 1,n_sp1de
    projct_g2f(:,i) = projct_g2f(:,i) * qdtwts_1d(:)
  end do
  !
  call output2d(projct_g2f(:,:),2,2,3,'projct_g2f')
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
  call alloc_error("interp2quadrature","qdtwts_2d",1,ierr)
  allocate ( vand_qdtr(1:n_qp2dq,1:n_qp2dq,1:2), source=zero , stat=ierr )
  call alloc_error("interp2quadrature","vand_qdtr",1,ierr)
  allocate ( vand_qdtr_inv(1:n_qp2dq,1:n_qp2dq,1:2), source=zero , stat=ierr )
  call alloc_error("interp2quadrature","vand_qdtr_inv",1,ierr)
  allocate ( interp_s2q(1:n_sp2dq,1:n_qp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("interp2quadrature","interp_s2q",1,ierr)
  allocate ( projct_q2s(1:n_qp2dq,1:n_sp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("interp2quadrature","projct_q2s",1,ierr)
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
  call output2d(interp_s2q(1:nsp,1:nqp,2),2,2,2,'interp_s2q')
  !
  ! Matrix for interpolating cell nodal values at SPs to face QPs
  !
  interp_s2g(1:nsp,1:4*n_qp1de,2) = LagrangePoly2D_quad(qdtpts_1d,solpts_1d)
  !
  call output2d(interp_s2g(1:nsp,1:4*n_qp1de,2),2,1,2,'interp_s2g')
  !
  ! 2D matrix to project nodal values at QPs to SPs
  !
  projct_q2s(1:nqp,1:nsp,2) = &
        matmul( vand2D_s2q(1:nqp,1:nsp,2) , transpose(vand_quad) )
  do i = 1,nsp
    projct_q2s(1:nqp,i,2) = projct_q2s(1:nqp,i,2) * qdtwts_2d(1:nqp,2)
  end do
  !
  call output2d(projct_q2s(1:nqp,1:nsp,2),2,2,2,'projct_q2s')
  !
  ! 2D Derivative matrix for the quadrature points
  !
  call output2d(dr_mat(1:nsp,1:nsp,2),2,2,2,'dr_mat')
  call output2d(ds_mat(1:nsp,1:nsp,2),2,2,2,'ds_mat')
  !
  call Dmatrix2D_quad(dq1d_mat,dqr_mat(1:nqp,1:nqp,2),dqs_mat(1:nqp,1:nqp,2))
  !
  call output2d(dqr_mat(1:nqp,1:nqp,2),2,2,2,'dqr_mat')
  call output2d(dqs_mat(1:nqp,1:nqp,2),2,2,2,'dqs_mat')
  !
  ! Pre-multiply the derivative matrices by the projection
  ! matrix for when using the derivative of the Lagrange polynomial
  !
  dqr_LPmat(1:nqp,1:nsp,2) = &
      matmul( transpose(dqr_mat(1:nqp,1:nqp,2)) , projct_q2s(1:nqp,1:nsp,2) )
  dqs_LPmat(1:nqp,1:nsp,2) = &
      matmul( transpose(dqs_mat(1:nqp,1:nqp,2)) , projct_q2s(1:nqp,1:nsp,2) )
  !
  call output2d(dqr_LPmat(1:nqp,1:nsp,2),2,2,2,'dqr_LPmat')
  call output2d(dqs_LPmat(1:nqp,1:nsp,2),2,2,2,'dqs_LPmat')
  !
  ! Pre-multiply the derivative matrices by the
  ! interpolation matrix for when using the chain-rule
  !
  dqr_CRmat(1:nsp,1:nqp,2) = &
      matmul( interp_s2q(1:nsp,1:nqp,2) , transpose(dqr_mat(1:nqp,1:nqp,2)) )
  dqs_CRmat(1:nsp,1:nqp,2) = &
      matmul( interp_s2q(1:nsp,1:nqp,2) , transpose(dqs_mat(1:nqp,1:nqp,2)) )
  !
  call output2d(dqr_CRmat(1:nsp,1:nqp,2),2,2,2,'dqr_CRmat')
  call output2d(dqs_CRmat(1:nsp,1:nqp,2),2,2,2,'dqs_CRmat')
  !
  ! Inverse of the mass matrix
  !
  invM_mat(1:nsp,1:nsp,2) = matmul( vand_quad , transpose(vand_quad) )
  where (abs(invM_mat) <= eps12)  invM_mat = zero
  !
  call output2d(invM_mat(1:nsp,1:nsp,2),2,2,2,'invM_mat')
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
subroutine output2d(array,dim1,dim2,ctype,arrayname)
  !
  !.. Use Statements ..
  use ifport
  !
  !.. Formal Arguments ..
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
  if (mypnum == 0 .and. array_debugging) then
    !
    fname = "array."//trim(adjustl(arrayname))
    open(iu,file=fname,status="replace",action="write",form="formatted",iostat=ierr)
    call io_error('output2d',fname,1,ierr)
    !
    np1 = size(array,1)
    np2 = size(array,2)
    !
    copy(:,:) = array(:,:)
    where (abs(copy) <= eps12) copy = zero
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
    command = "sed -i 's/\<0\.000\>/  -  /g' "//trim(adjustl(fname))
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
include "Functions/deltafun.f90"
include "Functions/reciprocal.f90"
include "Functions/np2n.f90"
include "Functions/np2nfp.f90"
include "Functions/chop.f90"
!
!###############################################################################
!
end module startup_mod
