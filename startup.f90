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
  use geovar, only : nr,n_sp1de,n_sp2dt,n_sp2dq,solpts_1d,rs,flxpts_1d,n_order
  use metvar, only : cormat,d1_mat,dr_mat,ds_mat,dr_hat,ds_hat,lagrange
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
  call Dmatrix2d_quad
  !
  call Dmatrix2d_weak
  !
  ! Get the Lagrange polynomial to project the variables/fluxes from the
  ! solution points in the cell interior to the flux points on the cell face
  !
  allocate ( lagrange(1:n_sp2dq,1:4*n_sp1de,1:2) , source=zero , stat=ierr )
  call alloc_error("get_quadrature_data","lagrange",1,ierr)
  !
  lagrange(1:n_sp2dt,1:3*n_sp1de,1) = LagrangePoly2D_tria(rs)
  lagrange(1:n_sp2dq,1:4*n_sp1de,2) = LagrangePoly2D_quad(flxpts_1d,solpts_1d)
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
 !write (101,*) "Quad Vandermonde Matrix"
 !1 format (a,i0,a,i0,a)
 !2 format (<n_sp1de>es12.4)
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
 !write (102,*) "Inverse of Quad Vandermonde Matrix"
 !do lj = 0,n_order
 !  do li = 0,n_order
 !    l = n_sp1de*lj + li + 1
 !    write (102,*)
 !    write (102,1) "Basis function for term (x^",li,")*(y^",lj,")"
 !    do kj = n_sp1de,1,-1
 !      write (102,2) (vand_quad_inv(n_sp1de*(kj-1)+ki,l),ki=1,n_sp1de)
 !    end do
 !  end do
 !end do
 !stop
 !!
 !write (101,*) "DG Correction Function"
 !do nf = 1,4*n_sp1de
 !  !
 !  if (nf <= 1*n_sp1de) then
 !    write (101,5) nf
 !  else if (nf <= 2*n_sp1de) then
 !    write (101,6) nf - n_sp1de
 !  else if (nf <= 3*n_sp1de) then
 !    write (101,7) n_sp1de+1-(nf-2*n_sp1de)
 !  else
 !    write (101,8) n_sp1de+1-(nf-3*n_sp1de)
 !  end if
 !  do j = n_sp1de,1,-1
 !    write (101,4) (cormat((j-1)*n_sp1de+i,nf,2),i=1,n_sp1de)
 !  end do
 !  !
 !end do
 !4 format (<n_sp1de>es12.4)
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
 !3 format (<n_sp1de>es12.4)
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
 !write (105,*) "Lagrange matrix"
 !do j = 1,4*n_sp1de
 !  write (*,1) (lagrange(i,j,2),i=1,n_sp2dq)
 !end do
 !!
 !1 format (<n_sp2dq>es12.4)
  !
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
subroutine Dmatrix2D_quad()
  !
  !.. Use Statements ..
  use geovar, only : nr,n_sp1de,n_sp2dq
  use metvar, only : d1_mat,dr_mat,ds_mat
  !
  !.. Local Scalars ..
  integer :: i,j,k,l
  integer :: lbeg,lend
  !
continue
  !
  ! Initialize the quad parts of the derivative matrices to zero
  !
  dr_mat(:,:,2) = zero
  ds_mat(:,:,2) = zero
  !
  ! First the r/xi direction
  !
  k = 0
  do j = 1,n_sp1de
    lbeg = (j-1)*n_sp1de + 1
    lend = j*n_sp1de
    do i = 1,n_sp1de
      k = k + 1
      dr_mat(lbeg:lend,k,2) = d1_mat(i,:)
    end do
  end do
  !
  ! Next the s/eta direction
  !
  k = 0
  do j = 1,n_sp1de
    do i = 1,n_sp1de
      k = k + 1
      do l = 1,n_sp1de
        ds_mat(i+(l-1)*n_sp1de,k,2) = d1_mat(j,l)
      end do
    end do
  end do
  !
 !do k = 1,n_sp2dq
 !  do j = 1,n_sp1de
 !    write (104,'(10es12.4)') (dr_mat(n_sp1de*(j-1)+i,k,2),i=1,n_sp1de)
 !    write (105,'(10es12.4)') (ds_mat(n_sp1de*(j-1)+i,k,2),i=1,n_sp1de)
 !  end do
 !  write (104,*)
 !  write (105,*)
 !end do
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
  use geovar, only : n_sp1de,n_sp2dq,solwts_1d,solpts_1d
  use metvar, only : lagrange,cormat_dg
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
 !     !cormat_dg(k,nf,2) = lagrange(k,nf,2)/solwts_1d(i)
 !      cormat_dg(k,nf,2) = eval_LagrangePoly1D(i,-one,solpts_1d)/solwts_1d(i)
 !    end do
 !  end do
 !  !
 !end do
 !!
 !write (98,*) "DGSEM Correction Function"
 !write (99,*) "Lagrange interpolation polynomial"
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
 !    write (99,1) (lagrange((j-1)*n_sp1de+i,nf,2),i=1,n_sp1de)
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
include "Functions/deltafun.f90"
include "Functions/reciprocal.f90"
include "Functions/np2n.f90"
include "Functions/np2nfp.f90"
!
!###############################################################################
!
end module startup_mod
