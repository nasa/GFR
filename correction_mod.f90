module correction_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use matrix_types_mod, only : matrix
  use matrix_types_mod, only : create_matrix_sparsity
  use matrix_types_mod, only : output_sparsity
  use matrix_types_mod, only : output_matrix
  !
  implicit none
  !
  private
  !
  !
  !.. Public Procedures ..
  !
  public :: init_correction_matrices
  public :: correction_memory_usage
  public :: move_correct_to_device
  !
  !
  !
  ! correct : correction matrix for all combinations of cell geometry and
  !           and cell order that are possible in the current simulation
  !
  !     For cell 'n' :
  !
  !         The correction matrix (giving the contribution to each solution
  !         point within cell 'n' from the interface flux point 'k' that is
  !         currently of interest) is given by the expression
  !
  !         correct(this_geom,this_order)%mat(:,k)
  !
  !     Where :
  !
  !        this_geom = cell(n)%geom
  !        this_order = cell(n)%order
  !        k = the flx-pt on a face of cell 'n' that is currently of interest
  !
  type(matrix), public, save, allocatable :: correct(:,:)
  !
  interface face_interpolation_matrix
    module procedure face_interpolation_matrix, &
                     face_interpolation_matrix_from_qdtr
  end interface face_interpolation_matrix
  !
contains
!
!###############################################################################
!
subroutine init_correction_matrices()
  !
  !.. Use Statements ..
  use order_mod,        only : geom_solpts,geom_flxpts
  use order_mod,        only : n_min_geom,n_max_geom
  use order_mod,        only : n_min_order,n_max_order
  use order_mod,        only : n_order
  use quadrature_mod,   only : geom_is_used
  use element_type_mod, only : std_elem
  !
  !.. Local Scalars ..
  integer :: gmin,gmax,omin,omax
  integer :: this_geom,this_order
  integer :: nsolpts,nflxpts,ierr
  integer :: n,n1,n2,i
  character(len=200) :: array_name
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: edgpts(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_correction_matrices"
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
  ! Allocate the correct array
  !
  allocate ( correct(gmin:gmax,omin:omax) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"correct",1,__LINE__,__FILE__,ierr,error_message)
  !
  geom_loop: do this_geom = gmin,gmax
    !
    if (.not. geom_is_used(this_geom)) cycle geom_loop
    if (all(this_geom /= Geom_Valid)) cycle geom_loop
    if (this_geom == Geom_Node) cycle geom_loop
    !
    order_loop: do this_order = omin,omax
      !
      ! Get the number of solution points for this cell
      !
      nsolpts = geom_solpts(this_geom,this_order)
      !
      ! Get the number of face/flux points for this cell
      !
      nflxpts = geom_flxpts(this_geom,this_order)
      !
      ! Allocate the correction matrix
      !
     !allocate ( this_correct%mat(1:nsolpts,1:nflxpts) , &
      allocate ( correct(this_geom,this_order)%mat(1:nsolpts,1:nflxpts) , &
                 source=zero , stat=ierr , errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Compute the correction matrices depending on geometry type
      !
      select case (this_geom)
        case (Geom_Edge)
          !
          edgpts = std_elem(this_geom,this_order)%pts(1,:) ! F2003 AUTO-REALLOC
          !
          correct(this_geom,this_order)%mat(:,1) = &
                     CorrectionMatrix_Edge( edgpts )
          !
          correct(this_geom,this_order)%mat(:,2) = &
                     correct(this_geom,this_order)%mat(nsolpts:1:-1,1)
          !
          error_message = "Correction"
          call output_matrix(correct(this_geom,this_order), &
                             error_message,this_geom,this_order, &
                             output_transpose=.true.)
          !
        case (Geom_Tria)
          !
          correct(this_geom,this_order)%mat = &
                     CorrectionMatrix_Tria( std_elem(this_geom,this_order)%pts )
#ifdef SPECIAL_OUTPUT_SPARSITY
          write (array_name,10) "tria",this_order
          call output_sparsity( correct(this_geom,this_order) , array_name )
#endif
          !
        case (Geom_Quad)
          !
          edgpts = std_elem(Geom_Edge,this_order)%pts(1,:) ! F2003 AUTO-REALLOC
          !
          correct(this_geom,this_order)%mat = CorrectionMatrix_Quad( edgpts )
#ifdef SPECIAL_OUTPUT_SPARSITY
          write (array_name,10) "quad",this_order
          call output_sparsity( correct(this_geom,this_order) , array_name )
#endif
          !
          error_message = "Correction"
          call output_matrix(correct(this_geom,this_order), &
                             error_message,this_geom,this_order, &
                             output_transpose=.true.)
          !
        case (Geom_Tetr)
          !
          correct(this_geom,this_order)%mat = &
                     CorrectionMatrix_Tetr( std_elem(this_geom,this_order)%pts )
#ifdef SPECIAL_OUTPUT_SPARSITY
          write (array_name,10) "tetr",this_order
          call output_sparsity( correct(this_geom,this_order) , array_name )
#endif
          !
        case (Geom_Pyra)
          !
          correct(this_geom,this_order)%mat = &
                     CorrectionMatrix_Pyra( std_elem(this_geom,this_order)%pts )
#ifdef SPECIAL_OUTPUT_SPARSITY
          write (array_name,10) "pyra",this_order
          call output_sparsity( correct(this_geom,this_order) , array_name )
#endif
          !
        case (Geom_Pris)
          !
          correct(this_geom,this_order)%mat = &
                     CorrectionMatrix_Pris( std_elem(this_geom,this_order)%pts )
#ifdef SPECIAL_OUTPUT_SPARSITY
          write (array_name,10) "pris",this_order
          call output_sparsity( correct(this_geom,this_order) , array_name )
#endif
          !
        case (Geom_Hexa)
          !
          edgpts = std_elem(Geom_Edge,this_order)%pts(1,:) ! F2003 AUTO-REALLOC
          !
          correct(this_geom,this_order)%mat = CorrectionMatrix_Hexa( edgpts )
#ifdef SPECIAL_OUTPUT_SPARSITY
          write (array_name,10) "hexa",this_order
          call output_sparsity( correct(this_geom,this_order) , array_name )
#endif
          !
          error_message = "Correction"
          call output_matrix(correct(this_geom,this_order), &
                             error_message,this_geom,this_order, &
                             output_transpose=.true.)
          !
        case default
          !
          write (error_message,2)
          call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
          !
      end select
      !
    end do order_loop
    !
  end do geom_loop
  !
  ! Check for sparsity within each of the correction matrices
  !
 !call correct%check_for_sparsity
  call create_matrix_sparsity( correct )
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("correct(",a,",",i0,")%mat")
  2 format (" A grid cell of an unknown geometry type was found!")
  !
  10 format (a4,".correct.p",i0,".tec.dat")
  !
end subroutine init_correction_matrices
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function CorrectionMatrix_Edge(xi,eval_deriv) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : correction_function
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xi
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: eval_deriv
  !
  !.. Function Result ..
  real(wp), dimension(1:size(xi)) :: return_value
  !
  !.. Local Scalars ..
  integer  :: n,k,j,npts,ndeg
  logical(lk) :: compute_derivative
  real(qp) :: cn,cnm1,cnm2
  !
  !.. Local Arrays ..
  real(qp), dimension(1:size(xi)) :: x
  real(qp), dimension(1:size(xi)) :: cmat
  real(qp), dimension(0:size(xi)) :: coef_c
  real(qp), dimension(0:size(xi)) :: coef_dc
  !
  real(qp), dimension( 0:size(xi), 1:size(xi)) :: xexp
  real(qp), dimension(-1:size(xi),-1:size(xi)) :: coef_Le
  real(qp), dimension( 0:size(xi), 0:size(xi)) :: coef_RR
  !
continue
  !
  ! #########################################################################
  !
  !       THE LOOP INDICES IN THIS FUNCTION ARE USED FOR THE FOLLOWING
  !       SPECIFIC TASKS
  !
  !    n : index used for looping through polynomial degrees / exponents
  !    k : index used for looping over the solution points
  !    j : index used for looping through the coefficients of a polynomial
  !
  ! #########################################################################
  !
  npts = size(xi) ! number of solution points
  ndeg = size(xi) ! polynomial degree of the correction function
  !
  compute_derivative = true
  if (present(eval_deriv)) then
    compute_derivative = eval_deriv
  end if
  !
  ! Initialize the function result
  !
  return_value(:) = zero
  !
  ! Initialize the local arrays
  !
  x(:) = real( xi(:) , kind=qp )
  !
  cmat(:)      = qzero
  coef_c(:)    = qzero
  coef_dc(:)   = qzero
  coef_Le(:,:) = qzero
  coef_RR(:,:) = qzero
  !
  ! Initialize specific coefficients of Legendre polynomial
  !
  coef_Le(0,0) = qone
  coef_Le(1,1) = qone
  !
  ! Get the value of each nodal coordinate raised to all possible powers
  !   using xk = x(k), xexp(0:neg,k) contains the values of
  !        [ xk^0, xk^1, xk^2, ..., xk^(ndeg-1), xk^ndeg]
  !
  xexp(0,:) = qone
  do k = 1,npts
    do n = 1,ndeg
      xexp(n,k) = xexp(n-1,k) * x(k)
    end do
  end do
  !
  ! Get the coefficients for every Legendre polynomial up to degree ndeg
  !
  do n = 2,ndeg
    !
    cnm1 = real( 2*n-1 , kind=qp ) / real( n , kind=qp )
    cnm2 = real(   n-1 , kind=qp ) / real( n , kind=qp )
    !
    ! Loop to get the coefficients for Legendre polynomial of degree n
    !
    do j = 0,n
      coef_Le(j,n) = cnm1*coef_Le(j-1,n-1) - cnm2*coef_Le(j,n-2)
    end do
    !
  end do
  !
  ! Get the coefficients for every right Radau polynomial up to degree ndeg
  !
  do n = 1,ndeg
    !
    cn = qhalf * real( (-1)**n , kind=qp )
    !
    ! Loop to get the coefficients for right Radau polynomial of degree n
    !
    do j = 0,n
      coef_RR(j,n) = cn * (coef_Le(j,n) - coef_Le(j,n-1))
    end do
    !
  end do
  !
  ! Get the coefficients for the correction function
  !
  if (correction_function == Correction_g2) then
    !
    ! Use the g2 correction function
    !
    cn   = real( ndeg-1 , kind=qp ) / real( 2*ndeg-1 , kind=qp )
    cnm1 = real( ndeg   , kind=qp ) / real( 2*ndeg-1 , kind=qp )
    !
    do j = 0,ndeg
      coef_c(j) = cn*coef_RR(j,ndeg) + cnm1*coef_RR(j,ndeg-1)
    end do
    !
  else if (correction_function == Correction_gGauss) then
    !
    ! Use the gGauss correction function
    !
    cn   = real( ndeg-1 , kind=qp ) / real( 2*ndeg-1 , kind=qp )
    cnm1 = real( ndeg   , kind=qp ) / real( 2*ndeg-1 , kind=qp )
    !
    do j = 0,ndeg
      coef_c(j) = cn*coef_RR(j,ndeg) + cnm1*coef_RR(j,ndeg-1)
    end do
    !
  else
    !
    ! Use the right Radau polynomial for the correction function
    !
    do j = 0,ndeg
      coef_c(j) = coef_RR(j,ndeg)
    end do
    !
  end if
  !
  if (compute_derivative) then
    !
    ! Compute the coefficients for the derivative of the correction function
    !
    do j = 1,ndeg
      coef_dc(j-1) = real( j , kind=qp ) * coef_c(j)
    end do
    !
    ! Finally, get the value for the derivative of the
    ! correction function at each 1D solution point
    !
    do k = 1,npts
      cmat(k) = -sum( xexp(0:ndeg-1,k) * coef_dc(0:ndeg-1) )
    end do
    !
  else
    !
    ! Finally get the value of the correction function
    ! at each 1D solution point
    !
    do k = 1,npts
      cmat(k) = sum( xexp(0:ndeg,k) * coef_c(0:ndeg) )
    end do
    !
  end if
  !
  ! Use the chop function to copy cmat into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( cmat )
  !
end function CorrectionMatrix_Edge
!
!###############################################################################
!
function CorrectionMatrix_Edge_LP(pts1d,wts1d) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly1D
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: pts1d
  real(wp), dimension(:), intent(in) :: wts1d
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: n,np,ierr
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: cmat_qp(:)
  real(qp), allocatable :: pts1d_qp(:)
  real(qp), allocatable :: wts1d_qp(:)
  !
continue
  !
  np = size(pts1d)
  !
  ! Allocate the return array containing the working
  ! precision correction matrix
  !
  allocate ( return_value(1:np) , source=zero , stat=ierr )
  !
  ! Create the correction matrix using quad precision and make quad precision
  ! copies of the 1D solution points and their corresponding quadrature weights
  !
  allocate ( cmat_qp(1:np) , source=qzero , stat=ierr )
  !
  pts1d_qp = real(pts1d,kind=kind(pts1d_qp)) ! F2003 AUTO-REALLOC
  wts1d_qp = real(wts1d,kind=kind(wts1d_qp)) ! F2003 AUTO-REALLOC
  !
  !
  !
  do n = 1,np
    cmat_qp(n) = eval_LagrangePoly1D(n,-qone,pts1d_qp) / wts1d_qp(n)
  end do
  !
  ! Use the chop function to copy cmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( cmat_qp )
  !
end function CorrectionMatrix_Edge_LP
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function CorrectionMatrix_Quad_LP(pts1d,wts1d) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly1D
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: pts1d
  real(wp), dimension(:), intent(in) :: wts1d
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,k,nf,np,nsp,nfp,ierr
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: cmat_qp(:,:)
  real(qp), allocatable :: pts1d_qp(:)
  real(qp), allocatable :: wts1d_qp(:)
  !
continue
  !
  np = size(pts1d)
  nsp = np*np
  nfp = 4*np
  !
  ! Allocate the return array containing the working
  ! precision correction matrix
  !
  allocate ( return_value(1:nsp,1:nfp) , source=zero , stat=ierr )
  !
  ! Create the correction matrix using quad precision and make quad precision
  ! copies of the 1D solution points and their corresponding quadrature weights
  !
  allocate ( cmat_qp(1:nsp,1:nfp) , source=qzero , stat=ierr )
  !
  pts1d_qp = real(pts1d,kind=kind(pts1d_qp)) ! F2003 AUTO-REALLOC
  wts1d_qp = real(wts1d,kind=kind(wts1d_qp)) ! F2003 AUTO-REALLOC
  !
  ! Face 1 - bottom/south edge of quad
  !
  nf = 1
  do i = 1,np
    do j = 1,np
      k = (j-1)*np + i
      cmat_qp(k,nf) = eval_LagrangePoly1D(j,-qone,pts1d_qp) / wts1d_qp(j)
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
      cmat_qp(k,nf) = eval_LagrangePoly1D(i,qone,pts1d_qp) / wts1d_qp(i)
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
      cmat_qp(k,nf) = eval_LagrangePoly1D(j,qone,pts1d_qp) / wts1d_qp(j)
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
      cmat_qp(k,nf) = eval_LagrangePoly1D(i,-qone,pts1d_qp) / wts1d_qp(i)
    end do
    nf = nf + 1
  end do
  !
  ! Use the chop function to copy cmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( cmat_qp )
  !
end function CorrectionMatrix_Quad_LP
!
!###############################################################################
!
function CorrectionMatrix_Quad_2D_LP(pts1d,wts1d) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly1D
  use polynomial_mod, only : eval_LagrangePoly2D
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: pts1d
  real(wp), dimension(:), intent(in) :: wts1d
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,k,n,nf,np,nsp,nfp,ierr
  real(qp) :: x,y
  !
  !.. Local Allocatable Arrays ..
  real(qp), allocatable :: cmat_qp(:,:)
  real(qp), allocatable :: pts1d_qp(:)
  real(qp), allocatable :: wts1d_qp(:)
  !
continue
  !
  np = size(pts1d)
  nsp = np*np
  nfp = 4*np
  !
  ! Allocate the return array containing the working
  ! precision correction matrix
  !
  allocate ( return_value(1:nsp,1:nfp) , source=zero , stat=ierr )
  !
  ! Create the correction matrix using quad precision and make quad precision
  ! copies of the 1D solution points and their corresponding quadrature weights
  !
  allocate ( cmat_qp(1:nsp,1:nfp) , source=qzero , stat=ierr )
  !
  pts1d_qp = real(pts1d,kind=kind(pts1d_qp)) ! F2003 AUTO-REALLOC
  wts1d_qp = real(wts1d,kind=kind(wts1d_qp)) ! F2003 AUTO-REALLOC
  !
  ! Face 1 - bottom/south edge of quad
  !
  write (iout,*)
  write (iout,*)
  write (iout,*)
  1 format (2x,i1,2x,"cmat(",i2,",",i2,")", &
               " = eval_LagrangePoly2D(",i1,",",i1,",", &
               f4.1,",",f4.1,",pts1d,pts1d) / wts1d(",i1,") = ",f11.6)
  !
  nf = 1
  y = -qone
  do n = 1,np
    x = pts1d_qp(n)
    do j = 1,np
      do i = 1,np
        k = (j-1)*np + i
        cmat_qp(k,nf) = eval_LagrangePoly2D(i,j,x,y,pts1d_qp,pts1d_qp) / &
                                                                     wts1d_qp(j)
        write (iout,1) n,k,nf,i,j,x,y,j,cmat_qp(k,nf)
      end do
    end do
    write (iout,*)
    nf = nf + 1
  end do
  !
 !nf = 1
 !do i = 1,np
 !  do j = 1,np
 !    k = (j-1)*np + i
 !    cmat_qp(k,nf) = eval_LagrangePoly1D(j,-qone,pts1d_qp) / wts1d_qp(j)
 !    write (iout,1) k,nf,j,"-1.0",j,cmat_qp(k,nf)
 !  end do
 !  nf = nf + 1
 !end do
  !
  ! Face 2 - right/east edge of quad
  !
  nf = np+1
  x = qone
  do n = 1,np
    y = pts1d_qp(n)
    do j = 1,np
      do i = 1,np
        k = (j-1)*np + i
        cmat_qp(k,nf) = eval_LagrangePoly2D(i,j,x,y,pts1d_qp,pts1d_qp) / &
                                                                     wts1d_qp(i)
        write (iout,1) n,k,nf,i,j,x,y,i,cmat_qp(k,nf)
      end do
    end do
    write (iout,*)
    nf = nf + 1
  end do
  !
 !nf = np+1
 !do j = 1,np
 !  do i = np,1,-1
 !    k = (j-1)*np + i
 !    cmat_qp(k,nf) = eval_LagrangePoly1D(i,qone,pts1d_qp) / wts1d_qp(i)
 !    write (iout,1) k,nf,i," 1.0",i,cmat_qp(k,nf)
 !  end do
 !  nf = nf + 1
 !end do
  !
  ! Face 3 - top/north edge of quad
  !
  nf = 2*np+1
  y = qone
  do n = 1,np
    x = pts1d_qp(np+1-n)
    do j = 1,np
      do i = 1,np
        k = (j-1)*np + i
        cmat_qp(k,nf) = eval_LagrangePoly2D(i,j,x,y,pts1d_qp,pts1d_qp) / &
                                                                     wts1d_qp(j)
        write (iout,1) np+1-n,k,nf,i,j,x,y,j,cmat_qp(k,nf)
      end do
    end do
    write (iout,*)
    nf = nf + 1
  end do
  !
 !nf = 2*np+1
 !do i = np,1,-1
 !  do j = np,1,-1
 !    k = (j-1)*np + i
 !    cmat_qp(k,nf) = eval_LagrangePoly1D(j,qone,pts1d_qp) / wts1d_qp(j)
 !    write (iout,1) k,nf,j," 1.0",j,cmat_qp(k,nf)
 !  end do
 !  nf = nf + 1
 !end do
  !
  ! Face 4 - left/west edge of quad
  !
  nf = 3*np+1
  x = -qone
  do n = 1,np
    y = pts1d_qp(np+1-n)
    do j = 1,np
      do i = 1,np
        k = (j-1)*np + i
        cmat_qp(k,nf) = eval_LagrangePoly2D(i,j,x,y,pts1d_qp,pts1d_qp) / &
                                                                     wts1d_qp(i)
        write (iout,1) np+1-n,k,nf,i,j,x,y,i,cmat_qp(k,nf)
      end do
    end do
    write (iout,*)
    nf = nf + 1
  end do
  !
 !nf = 3*np+1
 !do j = np,1,-1
 !  do i = 1,np
 !    k = (j-1)*np + i
 !    cmat_qp(k,nf) = eval_LagrangePoly1D(i,-qone,pts1d_qp) / wts1d_qp(i)
 !    write (iout,1) k,nf,i,"-1.0",i,cmat_qp(k,nf)
 !  end do
 !  nf = nf + 1
 !end do
  !
  write (iout,*)
  write (iout,*)
  write (iout,*)
  !
  ! Use the chop function to copy cmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( cmat_qp )
  !
end function CorrectionMatrix_Quad_2D_LP
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function CorrectionMatrix_Quad_from_Edge(xi,eval_deriv) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xi
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: eval_deriv
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,k,np,nf,nsp,nfp,ierr
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: cmat1D(:)
  !
continue
  !
  ! ###########################################################################
  ! ###########################################################################
  ! ####  NOTE: No need to use chopping function in this procedure. The    ####
  ! ####        quad correction matrix is really just the 1D correction    ####
  ! ####        matrix applied to each face/flux point of the standard     ####
  ! ####        quad element. The 1D correction matrix is obtained by the  ####
  ! ####        function CorrectionMatrix_Edge which applies the chopping  ####
  ! ####        function to its result so there is no need to repeat this  ####
  ! ####        operation again.                                           ####
  ! ###########################################################################
  ! ###########################################################################
  !
  np = size(xi)
  nsp = np*np
  nfp = 4*np
  !
  ! Allocate the return array containing the working
  ! precision correction matrix
  !
  allocate ( return_value(1:nsp,1:nfp) , source=zero , stat=ierr )
  !
  ! Create the 1D correction matrix
  !
  allocate ( cmat1D(1:np) , source=zero , stat=ierr )
  !
  ! 1D Correction matrix
  !
  cmat1D = CorrectionMatrix_Edge(xi,eval_deriv)
  !
  ! Face 1 - y/eta min face of quad (bottom/south)
  !
  nf = 1
  do j = 1,np
    do i = 1,np
      k = j + (i-1)*np
      return_value(k,nf) = cmat1D(i)
    end do
    nf = nf + 1
  end do
  !
  ! Face 2 - x/xi max face of quad (right/east)
  !
  nf = np+1
  do j = 1,np
    do i = 1,np
      k = (j-1)*np + np - i + 1
      return_value(k,nf) = cmat1D(i)
    end do
    nf = nf + 1
  end do
  !
  ! Face 3 - y/eta max face of quad (top/north)
  !
  nf = 2*np+1
  do j = 1,np
    do i = 1,np
      k = (np-i)*np + np - j + 1
      return_value(k,nf) = cmat1D(i)
    end do
    nf = nf + 1
  end do
  !
  ! Face 4 - x/xi min face of quad (left/west)
  !
  nf = 3*np+1
  do j = 1,np
    do i = 1,np
      k = (np-j)*np + i
      return_value(k,nf) = cmat1D(i)
    end do
    nf = nf + 1
  end do
  !
end function CorrectionMatrix_Quad_from_Edge
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function General_Correction_Matrix(this_geom,this_order) result(return_value)
  !
  !.. Use Statements ..
  use order_mod,        only : geom_solpts,geom_flxpts
  use element_type_mod, only : element_type
  use element_type_mod, only : std_elem,face_elem
  use quadrature_mod,   only : get_edge_quadrature
  use quadrature_mod,   only : get_tria_quadrature
  use quadrature_mod,   only : get_quad_quadrature
  use quadrature_mod,   only : get_ref_elem_node_idx_for_face
  use quadrature_mod,   only : map_face_pts
  use vandermonde_mod,  only : vand
  use vandermonde_mod,  only : General_Basis_Matrix
  !
  !.. Formal Arguments ..
  integer, intent(in) :: this_geom
  integer, intent(in) :: this_order
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: n1,n2,nsp,nfp,ndim
  integer :: this_face,face_geom
  !
  logical :: needsEdge,needsTria,needsQuad
  !
  type(element_type) :: edge_qdtr
  type(element_type) :: tria_qdtr
  type(element_type) :: quad_qdtr
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: ip(:)
  real(wp), allocatable :: rs(:,:)
  real(wp), allocatable :: node_xyz(:,:)
  real(wp), allocatable :: face_pts(:,:)
  real(wp), allocatable :: bmat(:,:)
  real(wp), allocatable :: smat(:,:)
  real(wp), allocatable :: bwimat(:,:)
  real(wp), allocatable :: cmat(:,:)
  !
  real(wp), allocatable :: imat(:,:)
  real(wp), allocatable :: wimat(:,:)
  real(wp), allocatable :: wimatEdge(:,:)
  real(wp), allocatable :: wimatTria(:,:)
  real(wp), allocatable :: wimatQuad(:,:)
  !
  type(matrix) :: tmp
  integer :: i,j
  !
continue
  !
  nsp = geom_solpts(this_geom,this_order)
  nfp = geom_flxpts(this_geom,this_order)
  !
  ! Special logic if this is for an edge element
  !
  if (this_geom == Geom_Edge) then
    !
    return_value = init_array(nsp,nfp,zero) ! F2003 AUTO-REALLOC
    !
    edge_qdtr = std_elem(Geom_Edge,this_order)
    !
    return_value(:,1) = CorrectionMatrix_Edge( edge_qdtr%pts(1,:) )
    return_value(:,2) = return_value(nsp:1:-1,1)
    !
    return
    !
  end if
 !!
 !needsEdge = any( this_geom == [Geom_Tria,Geom_Quad] )
 !needsTria = any( this_geom == [Geom_Tetr,Geom_Pyra,Geom_Pris] )
 !needsQuad = any( this_geom == [Geom_Pyra,Geom_Pris,Geom_Hexa] )
  !
  needsEdge = any( cell_face_geom(:,this_geom) == Geom_Edge )
  needsTria = any( cell_face_geom(:,this_geom) == Geom_Tria )
  needsQuad = any( cell_face_geom(:,this_geom) == Geom_Quad )
  !
  ! Compute the wimatEdge matrix if this element has edge faces.
  !
  if (needsEdge) then
    !
    ! Get the coordinates and weights for the edge-face
    ! quadrature/integration points
    !
    edge_qdtr = get_edge_quadrature(this_order,Legendre_Gauss)
    !
    ! Compute the interpolation matrix from the face flux points
    ! to the face quadrature/integration points
    !
    imat = face_interpolation_matrix(Geom_Edge, & ! F2003 AUTO-REALLOC
                                     face_elem(Geom_Edge,this_order), &
                                     edge_qdtr)
    !
    ! Compute diag(wts)*imat
    !
    ! F2003 AUTO-REALLOC
    wimatEdge = spread(edge_qdtr%wts,dim=1,ncopies=size(imat,dim=1))*imat
    !
  end if
  !
  ! Compute the wimatTria matrix if this element has triangle faces.
  !
  if (needsTria) then
    !
    ! Get the coordinates and weights for the triangle-face
    ! quadrature/integration points
    !
    tria_qdtr = get_tria_quadrature(this_order,WilliamsShunn_TriPoints)
    !
    ! Compute the interpolation matrix from the face flux points
    ! to the face quadrature/integration points
    !
    imat = face_interpolation_matrix(Geom_Tria, & ! F2003 AUTO-REALLOC
                                     face_elem(Geom_Tria,this_order), &
                                     tria_qdtr)
    !
    ! Compute diag(wts)*imat
    !
    ! F2003 AUTO-REALLOC
    wimatTria = spread(tria_qdtr%wts,dim=1,ncopies=size(imat,dim=1))*imat
    !
  end if
  !
  ! Compute the wimatQuad matrix if this element has quad faces.
  !
  if (needsQuad) then
    !
    ! Get the coordinates and weights for the quad-face
    ! quadrature/integration points
    !
    quad_qdtr = get_quad_quadrature(this_order,Legendre_Gauss)
    !
    ! Compute the interpolation matrix from the face flux points
    ! to the face quadrature/integration points
    !
    imat = face_interpolation_matrix(Geom_Quad, & ! F2003 AUTO-REALLOC
                                     face_elem(Geom_Quad,this_order), &
                                     quad_qdtr)
    !
    ! Compute diag(wts)*imat
    !
    ! F2003 AUTO-REALLOC
    wimatQuad = spread(quad_qdtr%wts,dim=1,ncopies=size(imat,dim=1))*imat
    !
  end if
  !
  ! Now we loop over the faces to compute the sections of bwimat
  !
  bwimat = init_array(nfp,nsp,zero) ! F2003 AUTO-REALLOC
  !
  n2 = 0
  ndim = geom_dimen(this_geom)
  !
  do this_face = 1,geom_faces(this_geom)
    !
    ! Get the geometry of the current face
    !
    face_geom = cell_face_geom(this_face,this_geom)
    !
    ! Load up rs and wimat depending what geometry the current face is
    !
    if (face_geom == Geom_Edge) then
      rs = edge_qdtr%pts ! F2003 AUTO-REALLOC
      wimat = wimatEdge ! F2003 AUTO-REALLOC
    else if (face_geom == Geom_Tria) then
      rs = tria_qdtr%pts ! F2003 AUTO-REALLOC
      wimat = wimatTria ! F2003 AUTO-REALLOC
    else if (face_geom == Geom_Quad) then
      rs = quad_qdtr%pts ! F2003 AUTO-REALLOC
      wimat = wimatQuad ! F2003 AUTO-REALLOC
    end if
    !
    ! Get the extents of the face points for the current face within the
    ! range of all face points across all faces of the current cell
    !
    n1 = n2 + 1
    n2 = n2 + size(rs,dim=2)
    !
    ! Get the node indices of the current face for the reference element
    !
    ip = get_ref_elem_node_idx_for_face(this_geom,this_face) ! F2003 AUTO-REALLOC
    !
    ! Get the coordinates of the nodes for the current face
    !
    node_xyz = ref_elem_nodes_xyz(1:ndim,ip(:),this_geom) ! F2003 AUTO-REALLOC
    !
    ! Map the face integration points onto the
    ! current face of the reference element
    !
    face_pts = map_face_pts(node_xyz,rs) ! F2003 AUTO-REALLOC
 !  !
 !  write (error_message,1) "face_pts",this_face
 !  tmp%mat = transpose(face_pts) ! F2003 AUTO-REALLOC
 !  call output_matrix(tmp,error_message,this_geom,this_order)
    !
    ! Compute a basis matrix containing the values of the cell orthonormal
    ! basis functions evaluated at the face integration points
    !
    bmat = General_Basis_Matrix(this_geom,this_order,face_pts) ! F2003 AUTO-REALLOC
 !  !
 !  write (error_message,1) "bmat",this_face
 !  tmp%mat = transpose(bmat) ! F2003 AUTO-REALLOC
 !  call output_matrix(tmp,error_message,this_geom,this_order)
    !
    smat = matmul( transpose(bmat) , transpose(wimat) ) ! F2003 AUTO-REALLOC
 !  !
 !  write (error_message,1) "smat",this_face
 !  tmp%mat = transpose(smat) ! F2003 AUTO-REALLOC
 !  call output_matrix(tmp,error_message,this_geom,this_order)
    !
    bwimat(n1:n2,1:nsp) = transpose(smat(:,:))
    !
  end do
 !1 format (a," (Face #",i1,")")
 !!
 !error_message = "bwimat"
 !tmp%mat = bwimat ! F2003 AUTO-REALLOC
 !call output_matrix(tmp,error_message,this_geom,this_order)
 !!
 !error_message = "vand"
 !call output_matrix(vand(this_geom,this_order)%modal2nodal, &
 !                   error_message,this_geom,this_order)
  !
  ! F2003 AUTO-REALLOC
  cmat = matmul( bwimat , transpose(vand(this_geom,this_order)%modal2nodal%mat) )
 !!
 !error_message = "cmat"
 !tmp%mat = cmat ! F2003 AUTO-REALLOC
 !call output_matrix(tmp,error_message,this_geom,this_order)
  !
  return_value = chop( real( transpose(cmat) , kind=qp ) ) ! F2003 AUTO-REALLOC
  !
end function General_Correction_Matrix
!
!###############################################################################
!
function CorrectionMatrix_Quad(this_order) result(return_value)
  !
  !.. Use Statements ..
  use order_mod,        only : geom_solpts,geom_flxpts
  use element_type_mod, only : element_type
  use element_type_mod, only : face_elem
  use quadrature_mod,   only : get_edge_quadrature
  use quadrature_mod,   only : get_ref_elem_node_idx_for_face
  use quadrature_mod,   only : map_face_pts
  use vandermonde_mod,  only : vand
  use vandermonde_mod,  only : General_Basis_Matrix
  !
  !.. Formal Arguments ..
  integer, intent(in) :: this_order
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: n1,n2,nsp,nfp,ndim
  integer :: this_face,face_geom
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: ip(:)
  real(wp), allocatable :: imat_edge(:,:)
  real(wp), allocatable :: wimat(:,:)
  real(wp), allocatable :: rs(:,:)
  real(wp), allocatable :: node_xyz(:,:)
  real(wp), allocatable :: face_pts(:,:)
  real(wp), allocatable :: bmat(:,:)
  real(wp), allocatable :: smat(:,:)
  real(wp), allocatable :: bwimat(:,:)
  real(wp), allocatable :: cmat(:,:)
  !
  type(element_type) :: face_qdtr
  type(matrix) :: tmp
  integer :: i,j
  !
continue
  !
  nsp = geom_solpts(Geom_Quad,this_order)
  nfp = geom_flxpts(Geom_Quad,this_order)
  !
  ! Get the coordinates and weights for the face quadrature/integration points
  !
  face_qdtr = get_edge_quadrature(this_order,Legendre_Gauss)
  !
  ! Compute the interpolation matrix from the face flux points
  ! to the face quadrature/integration points
  !
  ! F2003 AUTO-REALLOC
  imat_edge = face_interpolation_matrix(Geom_Edge, &
                                        face_elem(Geom_Edge,this_order), &
                                        face_qdtr)
 !!
 !error_message = "edge-face interpolation (imat_edge)"
 !tmp%mat = imat_edge ! F2003 AUTO-REALLOC
 !call output_matrix(tmp,error_message,Geom_Edge,this_order)
 !!
 !! Compute the interpolation matrix from the face flux points
 !! to the face quadrature/integration points
 !!
 !! F2003 AUTO-REALLOC
 !imat_edge = face_interpolation_matrix(Geom_Edge,this_order,loc_flux_pts, &
 !                                                this_order,Legendre_Gauss)
 !!
 !error_message = "edge-face interpolation (imat_edge)"
 !tmp%mat = imat_edge ! F2003 AUTO-REALLOC
 !call output_matrix(tmp,error_message,Geom_Edge,this_order)
  !
  ! Compute diag(wts)*imat
  !
  ! F2003 AUTO-REALLOC
  wimat = spread(face_qdtr%wts,dim=1,ncopies=size(imat_edge,dim=1))*imat_edge
 !!
 !error_message = "weighted edge-face interpolation (wimat)"
 !tmp%mat = wimat ! F2003 AUTO-REALLOC
 !call output_matrix(tmp,error_message,Geom_Edge,this_order)
  !
  bwimat = init_array(nfp,nsp,zero) ! F2003 AUTO-REALLOC
  !
  n2 = 0
  ndim = geom_dimen(Geom_Quad)
  !
  do this_face = 1,geom_faces(Geom_Quad)
    !
    face_geom = cell_face_geom(this_face,Geom_Quad)
    !
    rs = face_qdtr%pts ! F2003 AUTO-REALLOC
    !
    n1 = n2 + 1
    n2 = n2 + size(rs,dim=2)
    !
    ip = get_ref_elem_node_idx_for_face(Geom_Quad,this_face) ! F2003 AUTO-REALLOC
    !
    node_xyz = ref_elem_nodes_xyz(1:ndim,ip(:),Geom_Quad) ! F2003 AUTO-REALLOC
    !
    face_pts = map_face_pts(node_xyz,rs) ! F2003 AUTO-REALLOC
 !  !
 !  write (error_message,1) "face_pts",this_face
 !  tmp%mat = transpose(face_pts) ! F2003 AUTO-REALLOC
 !  call output_matrix(tmp,error_message,Geom_Quad,this_order)
    !
    bmat = General_Basis_Matrix(Geom_Quad,this_order,face_pts) ! F2003 AUTO-REALLOC
 !  !
 !  write (error_message,1) "bmat",this_face
 !  tmp%mat = transpose(bmat) ! F2003 AUTO-REALLOC
 !  call output_matrix(tmp,error_message,Geom_Quad,this_order)
    !
    smat = matmul( transpose(bmat) , transpose(wimat) ) ! F2003 AUTO-REALLOC
 !  !
 !  write (error_message,1) "smat",this_face
 !  tmp%mat = transpose(smat) ! F2003 AUTO-REALLOC
 !  call output_matrix(tmp,error_message,Geom_Quad,this_order)
    !
    bwimat(n1:n2,1:nsp) = transpose(smat(:,:))
    !
  end do
 !1 format (a," (Face #",i1,")")
 !!
 !error_message = "bwimat"
 !tmp%mat = bwimat ! F2003 AUTO-REALLOC
 !call output_matrix(tmp,error_message,Geom_Quad,this_order)
 !!
 !error_message = "vand"
 !call output_matrix(vand(Geom_Quad,this_order)%modal2nodal, &
 !                   error_message,Geom_Quad,this_order)
  !
  ! F2003 AUTO-REALLOC
  cmat = matmul( bwimat , transpose(vand(Geom_Quad,this_order)%modal2nodal%mat) )
 !!
 !error_message = "cmat"
 !tmp%mat = cmat ! F2003 AUTO-REALLOC
 !call output_matrix(tmp,error_message,Geom_Quad,this_order)
  !
  return_value = chop( real( transpose(cmat) , kind=qp ) ) ! F2003 AUTO-REALLOC
  !
end function CorrectionMatrix_Quad
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function face_interpolation_matrix(face_geom,from_order,from_loc, &
                                   to_order,to_loc) result(return_value)
  !
  ! This function is needed because we cant use General_Interp_Matrix here in
  ! the situation where the face flux points use different quadrature points
  ! from the interior solution points of the cell. The matrices stored in
  ! vand(cell_geom,cell_order) are based on the input loc_solution_pts and the
  ! Vandermonde matrix needed here needs to be based of of loc_flux_pts.
  !
  !.. Use Statements ..
  use order_mod,        only : geom_solpts,geom_flxpts
  use element_type_mod, only : element_type
  use quadrature_mod,   only : get_edge_quadrature
  use quadrature_mod,   only : get_tria_quadrature
  use quadrature_mod,   only : get_quad_quadrature
  use vandermonde_mod,  only : Compute_Basis_Matrix
  use vandermonde_mod,  only : vandermonde_inverse
  use vandermonde_mod,  only : General_Basis_Matrix
  !
  !.. Formal Arguments ..
  integer, intent(in) :: face_geom
  integer, intent(in) :: from_order
  integer, intent(in) :: from_loc
  integer, intent(in) :: to_order
  integer, intent(in) :: to_loc
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: npf,npt
  !
  type(element_type) :: from_qdtr
  type(element_type) :: to_qdtr
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: vand(:,:)
  real(wp), allocatable :: vinv(:,:)
  real(wp), allocatable :: bmat(:,:)
  !
continue
  !
  npf = geom_solpts(face_geom,from_order)
  npt = geom_solpts(face_geom,to_order)
  !
  return_value = init_array(npf,npt,zero) ! F2003 AUTO-REALLOC
  !
  select case (face_geom)
    !
    case (Geom_Edge)
      !
      from_qdtr = get_edge_quadrature(from_order,from_loc)
      to_qdtr = get_edge_quadrature(to_order,to_loc)
      !
    case (Geom_Tria)
      !
      from_qdtr = get_tria_quadrature(from_order,from_loc)
      to_qdtr = get_tria_quadrature(to_order,to_loc)
      !
    case (Geom_Quad)
      !
      from_qdtr = get_quad_quadrature(from_order,from_loc)
      to_qdtr = get_quad_quadrature(to_order,to_loc)
      !
    case default
      !
      return
      !
  end select
  !
  ! F2003 AUTO-REALLOC
  vand = Compute_Basis_Matrix(face_geom,from_order,nodal_points=from_qdtr%pts)
  vinv = Vandermonde_Inverse(vand) ! F2003 AUTO-REALLOC
  !
  ! F2003 AUTO-REALLOC
  bmat = General_Basis_Matrix(face_geom,from_order,to_qdtr%pts)
  !
  return_value(1:npf,1:npt) = matmul( transpose(vinv) , transpose(bmat) )
  !
end function face_interpolation_matrix
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function face_interpolation_matrix_from_qdtr(face_geom,from_qdtr,to_qdtr) &
                                      result(return_value)
  !
  ! This function is needed because we cant use General_Interp_Matrix here in
  ! the situation where the face flux points use different quadrature points
  ! from the interior solution points of the cell. The matrices stored in
  ! vand(cell_geom,cell_order) are based on the input loc_solution_pts and the
  ! Vandermonde matrix needed here needs to be based of of loc_flux_pts.
  !
  !.. Use Statements ..
  use element_type_mod, only : element_type
  use vandermonde_mod,  only : vandermonde_inverse
  use vandermonde_mod,  only : General_Basis_Matrix
  !
  !.. Formal Arguments ..
  integer,            intent(in) :: face_geom
  type(element_type), intent(in) :: from_qdtr
  type(element_type), intent(in) :: to_qdtr
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: npf,npt,from_order
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: vand(:,:)
  real(wp), allocatable :: vinv(:,:)
  real(wp), allocatable :: bmat(:,:)
  !
continue
  !
  npf = size(from_qdtr%wts)
  npt = size(to_qdtr%wts)
  !
  from_order = solpts2order(face_geom,npf)
  !
  return_value = init_array(npf,npt,zero) ! F2003 AUTO-REALLOC
  !
  ! F2003 AUTO-REALLOC
  vand = General_Basis_Matrix(face_geom,from_order,xyz=from_qdtr%pts)
  vinv = Vandermonde_Inverse(vand) ! F2003 AUTO-REALLOC
  !
  ! F2003 AUTO-REALLOC
  bmat = General_Basis_Matrix(face_geom,from_order,to_qdtr%pts)
  !
  return_value(1:npf,1:npt) = matmul( transpose(vinv) , transpose(bmat) )
  !
end function face_interpolation_matrix_from_qdtr
!
!###############################################################################
!
function CorrectionMatrix_Tria(this_order) result(return_value)
  !
  ! RS : r and s coordinates for the standard
  !      triangle element on the interval [-1,1]
  !
  !.. Use Statements ..
  use order_mod,        only : geom_solpts,geom_flxpts
  use ovar,             only : loc_flux_pts
  use element_type_mod, only : element_type
  use element_type_mod, only : face_elem
  use quadrature_mod,   only : get_edge_quadrature
  use quadrature_mod,   only : get_ref_elem_node_idx_for_face
  use quadrature_mod,   only : map_face_pts
  use vandermonde_mod,  only : vand
  use vandermonde_mod,  only : General_Basis_Matrix
  !
  !.. Formal Arguments ..
  integer, intent(in) :: this_order
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: n1,n2,nsp,nfp,ndim
  integer :: this_face,face_geom
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: ip(:)
  real(wp), allocatable :: imat_edge(:,:)
  real(wp), allocatable :: wimat(:,:)
  real(wp), allocatable :: rs(:,:)
  real(wp), allocatable :: node_xyz(:,:)
  real(wp), allocatable :: face_pts(:,:)
  real(wp), allocatable :: bmat(:,:)
  real(wp), allocatable :: smat(:,:)
  real(wp), allocatable :: bwimat(:,:)
  real(wp), allocatable :: cmat(:,:)
  !
  type(element_type) :: face_qdtr
  type(matrix) :: tmp
  !
continue
  !
  nsp = geom_solpts(Geom_Tria,this_order)
  nfp = geom_flxpts(Geom_Tria,this_order)
  !
  ! Get the coordinates and weights for the face quadrature/integration points
  !
  face_qdtr = get_edge_quadrature(this_order,Legendre_Gauss)
  !
  ! Compute the interpolation matrix from the face flux points
  ! to the face quadrature/integration points
  !
  ! F2003 AUTO-REALLOC
  imat_edge = face_interpolation_matrix(Geom_Edge, &
                                        face_elem(Geom_Edge,this_order), &
                                        face_qdtr)
 !!
 !error_message = "edge-face interpolation (imat_edge)"
 !tmp%mat = imat_edge ! F2003 AUTO-REALLOC
 !call output_matrix(tmp,error_message,Geom_Edge,this_order)
 !!
 !! Compute the interpolation matrix from the face flux points
 !! to the face quadrature/integration points
 !!
 !! F2003 AUTO-REALLOC
 !imat_edge = face_interpolation_matrix(Geom_Edge,this_order,loc_flux_pts, &
 !                                                this_order,Legendre_Gauss)
 !!
 !error_message = "edge-face interpolation (imat_edge)"
 !tmp%mat = imat_edge ! F2003 AUTO-REALLOC
 !call output_matrix(tmp,error_message,Geom_Edge,this_order)
  !
  ! Compute diag(wts)*imat
  !
  ! F2003 AUTO-REALLOC
  wimat = spread(face_qdtr%wts,dim=1,ncopies=size(imat_edge,dim=1))*imat_edge
 !!
 !error_message = "weighted edge-face interpolation (wimat)"
 !tmp%mat = wimat ! F2003 AUTO-REALLOC
 !call output_matrix(tmp,error_message,Geom_Edge,this_order)
  !
  bwimat = init_array(nfp,nsp,zero) ! F2003 AUTO-REALLOC
  !
  n2 = 0
  ndim = geom_dimen(Geom_Tria)
  !
  do this_face = 1,geom_faces(Geom_Tria)
    !
    face_geom = cell_face_geom(this_face,Geom_Tria)
    !
    rs = face_qdtr%pts ! F2003 AUTO-REALLOC
    !
    n1 = n2 + 1
    n2 = n2 + size(rs,dim=2)
    !
    ip = get_ref_elem_node_idx_for_face(Geom_Tria,this_face) ! F2003 AUTO-REALLOC
    !
    node_xyz = ref_elem_nodes_xyz(1:ndim,ip(:),Geom_Tria) ! F2003 AUTO-REALLOC
    !
    face_pts = map_face_pts(node_xyz,rs) ! F2003 AUTO-REALLOC
 !  !
 !  write (error_message,1) "face_pts",this_face
 !  tmp%mat = transpose(face_pts) ! F2003 AUTO-REALLOC
 !  call output_matrix(tmp,error_message,Geom_Tria,this_order)
    !
    bmat = General_Basis_Matrix(Geom_Tria,this_order,face_pts) ! F2003 AUTO-REALLOC
 !  !
 !  write (error_message,1) "bmat",this_face
 !  tmp%mat = transpose(bmat) ! F2003 AUTO-REALLOC
 !  call output_matrix(tmp,error_message,Geom_Tria,this_order)
    !
    smat = matmul( transpose(bmat) , transpose(wimat) ) ! F2003 AUTO-REALLOC
 !  !
 !  write (error_message,1) "smat",this_face
 !  tmp%mat = transpose(smat) ! F2003 AUTO-REALLOC
 !  call output_matrix(tmp,error_message,Geom_Tria,this_order)
    !
    bwimat(n1:n2,1:nsp) = transpose(smat(:,:))
    !
  end do
 !1 format (a," (Face #",i1,")")
 !!
 !error_message = "bwimat"
 !tmp%mat = bwimat ! F2003 AUTO-REALLOC
 !call output_matrix(tmp,error_message,Geom_Tria,this_order)
 !!
 !error_message = "vand"
 !call output_matrix(vand(Geom_Tria,this_order)%modal2nodal, &
 !                   error_message,Geom_Tria,this_order)
  !
  ! F2003 AUTO-REALLOC
  cmat = matmul( bwimat , transpose(vand(Geom_Tria,this_order)%modal2nodal%mat) )
 !!
 !error_message = "cmat"
 !tmp%mat = cmat ! F2003 AUTO-REALLOC
 !call output_matrix(tmp,error_message,Geom_Tria,this_order)
  !
  return_value = chop( real( transpose(cmat) , kind=qp ) ) ! F2003 AUTO-REALLOC
  !
end function CorrectionMatrix_Tria
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function CorrectionMatrix_Tria_backup(rs,eval_deriv) result(return_value)
  !
  ! RS : r and s coordinates for the standard
  !      triangle element on the interval [-1,1]
  !
  !.. Use Statements ..
  use vandermonde_mod, only : vand
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: eval_deriv
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,n1,n2,np,nsp,nfp
  integer :: this_order,ierr
  logical(lk) :: compute_derivative
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: Fmask(:,:)
  real(qp), allocatable :: cmat_qp(:,:)
  real(qp), allocatable :: mmat1D_qp(:,:)
  real(qp), allocatable :: vedge_qp(:,:)
  real(qp), allocatable :: vtria_qp(:,:)
  !
continue
  !
  nsp = size(rs,dim=2)
  this_order = solpts2order(Geom_Tria,nsp)
  np = this_order+1
  nfp = 3*np
  !
  compute_derivative = true
  if (present(eval_deriv)) then
    compute_derivative = eval_deriv
  end if
  !
  ! Allocate the return array containing the working
  ! precision correction matrix
  !
  allocate ( return_value(1:nsp,1:nfp) , source=zero , stat=ierr )
  !
  ! Create the correction matrix using quad precision and allocate the
  ! inverse mass matrix and Fmask arrays.
  !
  allocate ( cmat_qp(1:nsp,1:nfp) , source=qzero , stat=ierr )
  allocate ( mmat1D_qp(1:np,1:np) , source=qzero , stat=ierr )
  allocate ( Fmask(1:np,1:3) , source=0 , stat=ierr )
  !
  ! Make quad precision copies of the Edge and Tria Vandermonde matrices.
  !
  allocate ( vedge_qp(1:np,1:np) , source=qzero , stat=ierr )
  allocate ( vtria_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  !
  vedge_qp = real( vand(Geom_Edge,this_order)%modal2nodal%mat , kind=qp )
  vtria_qp = real( vand(Geom_Tria,this_order)%modal2nodal%mat , kind=qp )
  !
  ! Face 1 - bottom edge of triangle
  !
  do i = 1,np
    Fmask(i,1) = i
  end do
  !
  ! Face 2 - diagonal edge of triangle
  !
  j = 0
  do i = 1,np
    j = j + np + 1 - i
    Fmask(i,2) = j
  end do
  !
  ! Face 3 - vertical edge of triangle
  !
  j = nsp + 1
  do i = 1,np
    j = j - i
    Fmask(i,3) = j
  end do
  !
  ! DG Book's method for getting Fmask
  !
 !Fmask(:,1) = pack( (/(i,i=1,nsp)/) , abs(rs(2,:) + one    ) < eps12 )
 !Fmask(:,2) = pack( (/(i,i=1,nsp)/) , abs(rs(1,:) + rs(2,:)) < eps12 )
 !Fmask(:,3) = pack( (/(i,i=1,nsp)/) , abs(rs(1,:) + one    ) < eps12 )
  !
  ! Compute the inverse mass matrix along an edge
  !
#ifdef SPECIAL_FOR_INTEL
  mmat1D_qp = invert_matrix_lapack95( matmul( vedge_qp , transpose(vedge_qp) ) )
#else
  mmat1D_qp = invert_matrix( matmul( vedge_qp , transpose(vedge_qp) ) )
#endif
  !
  ! Face 1
  !
  n1 = 1  ! index of first node on this face
  n2 = np ! index of last  node on this face
  !
  cmat_qp(Fmask(:,1),n1:n2) = mmat1D_qp
  !
  ! Face 2
  !
  n1 = np+1 ! index of first node on this face
  n2 = 2*np ! index of last  node on this face
  !
  cmat_qp(Fmask(:,2),n1:n2) = mmat1D_qp
  !
  ! Face 3
  !
  n1 = 2*np+1 ! index of first node on this face
  n2 = 3*np   ! index of last  node on this face
  !
  cmat_qp(Fmask(:,3),n1:n2) = mmat1D_qp
  !
  ! Finally, compute the correction matrix
  !
  cmat_qp = matmul( transpose(vtria_qp) , cmat_qp )
  cmat_qp = matmul( vtria_qp , cmat_qp )
  if (.not. tria_is_biunit) then
    cmat_qp = qtwo*cmat_qp
  end if
  !
  ! Use the chop function to copy cmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( cmat_qp )
  !
end function CorrectionMatrix_Tria_backup
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function CorrectionMatrix_Tetr(rs,eval_deriv) result(return_value)
  !
  ! RS : r and s coordinates for the standard
  !      triangle element on the interval [-1,1]
  !
  !.. Use Statements ..
  use vandermonde_mod, only : vand
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: eval_deriv
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,n1,n2,np,nsp,nfp
  integer :: this_order,ierr
  logical(lk) :: compute_derivative
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: Fmask(:,:)
  real(qp), allocatable :: cmat_qp(:,:)
  real(qp), allocatable :: mmat1D_qp(:,:)
  real(qp), allocatable :: vedge_qp(:,:)
  real(qp), allocatable :: vtria_qp(:,:)
  !
continue
  !
  nsp = size(rs,dim=2)
  this_order = solpts2order(Geom_Tetr,nsp)
  np = this_order+1
  nfp = 3*np
  !
  compute_derivative = true
  if (present(eval_deriv)) then
    compute_derivative = eval_deriv
  end if
  !
  ! Allocate the return array containing the working
  ! precision correction matrix
  !
  allocate ( return_value(1:nsp,1:nfp) , source=zero , stat=ierr )
  !
  ! Create the correction matrix using quad precision and allocate the
  ! inverse mass matrix and Fmask arrays.
  !
  allocate ( cmat_qp(1:nsp,1:nfp) , source=qzero , stat=ierr )
  allocate ( mmat1D_qp(1:np,1:np) , source=qzero , stat=ierr )
  allocate ( Fmask(1:np,1:3) , source=0 , stat=ierr )
  !
  ! Make quad precision copies of the Edge and Tria Vandermonde matrices.
  !
  allocate ( vedge_qp(1:np,1:np) , source=qzero , stat=ierr )
  allocate ( vtria_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  !
  vedge_qp = real( vand(Geom_Edge,this_order)%modal2nodal%mat , kind=qp )
  vtria_qp = real( vand(Geom_Tria,this_order)%modal2nodal%mat , kind=qp )
  !
  ! Face 1 - bottom edge of triangle
  !
  do i = 1,np
    Fmask(i,1) = i
  end do
  !
  ! Face 2 - diagonal edge of triangle
  !
  j = 0
  do i = 1,np
    j = j + np + 1 - i
    Fmask(i,2) = j
  end do
  !
  ! Face 3 - vertical edge of triangle
  !
  j = nsp + 1
  do i = 1,np
    j = j - i
    Fmask(i,3) = j
  end do
  !
  ! DG Book's method for getting Fmask
  !
 !Fmask(:,1) = pack( (/(i,i=1,nsp)/) , abs(rs(2,:) + one    ) < eps12 )
 !Fmask(:,2) = pack( (/(i,i=1,nsp)/) , abs(rs(1,:) + rs(2,:)) < eps12 )
 !Fmask(:,3) = pack( (/(i,i=1,nsp)/) , abs(rs(1,:) + one    ) < eps12 )
  !
  ! Compute the inverse mass matrix along an edge
  !
#ifdef SPECIAL_FOR_INTEL
  mmat1D_qp = invert_matrix_lapack95( matmul( vedge_qp , transpose(vedge_qp) ) )
#else
  mmat1D_qp = invert_matrix( matmul( vedge_qp , transpose(vedge_qp) ) )
#endif
  !
  ! Face 1
  !
  n1 = 1  ! index of first node on this face
  n2 = np ! index of last  node on this face
  !
  cmat_qp(Fmask(:,1),n1:n2) = mmat1D_qp
  !
  ! Face 2
  !
  n1 = np+1 ! index of first node on this face
  n2 = 2*np ! index of last  node on this face
  !
  cmat_qp(Fmask(:,2),n1:n2) = mmat1D_qp
  !
  ! Face 3
  !
  n1 = 2*np+1 ! index of first node on this face
  n2 = 3*np   ! index of last  node on this face
  !
  cmat_qp(Fmask(:,3),n1:n2) = mmat1D_qp
  !
  ! Finally, compute the correction matrix
  !
  cmat_qp = matmul( transpose(vtria_qp) , cmat_qp )
  cmat_qp = qtwo*matmul( vtria_qp , cmat_qp )
  !
  ! Use the chop function to copy cmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( cmat_qp )
  !
end function CorrectionMatrix_Tetr
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function CorrectionMatrix_Pyra(rs,eval_deriv) result(return_value)
  !
  ! RS : r and s coordinates for the standard
  !      triangle element on the interval [-1,1]
  !
  !.. Use Statements ..
  use vandermonde_mod, only : vand
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: eval_deriv
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,n1,n2,np,nsp,nfp
  integer :: this_order,ierr
  logical(lk) :: compute_derivative
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: Fmask(:,:)
  real(qp), allocatable :: cmat_qp(:,:)
  real(qp), allocatable :: mmat1D_qp(:,:)
  real(qp), allocatable :: vedge_qp(:,:)
  real(qp), allocatable :: vtria_qp(:,:)
  !
continue
  !
  nsp = size(rs,dim=2)
  this_order = solpts2order(Geom_Pyra,nsp)
  np = this_order+1
  nfp = 3*np
  !
  compute_derivative = true
  if (present(eval_deriv)) then
    compute_derivative = eval_deriv
  end if
  !
  ! Allocate the return array containing the working
  ! precision correction matrix
  !
  allocate ( return_value(1:nsp,1:nfp) , source=zero , stat=ierr )
  !
  ! Create the correction matrix using quad precision and allocate the
  ! inverse mass matrix and Fmask arrays.
  !
  allocate ( cmat_qp(1:nsp,1:nfp) , source=qzero , stat=ierr )
  allocate ( mmat1D_qp(1:np,1:np) , source=qzero , stat=ierr )
  allocate ( Fmask(1:np,1:3) , source=0 , stat=ierr )
  !
  ! Make quad precision copies of the Edge and Tria Vandermonde matrices.
  !
  allocate ( vedge_qp(1:np,1:np) , source=qzero , stat=ierr )
  allocate ( vtria_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  !
  vedge_qp = real( vand(Geom_Edge,this_order)%modal2nodal%mat , kind=qp )
  vtria_qp = real( vand(Geom_Tria,this_order)%modal2nodal%mat , kind=qp )
  !
  ! Face 1 - bottom edge of triangle
  !
  do i = 1,np
    Fmask(i,1) = i
  end do
  !
  ! Face 2 - diagonal edge of triangle
  !
  j = 0
  do i = 1,np
    j = j + np + 1 - i
    Fmask(i,2) = j
  end do
  !
  ! Face 3 - vertical edge of triangle
  !
  j = nsp + 1
  do i = 1,np
    j = j - i
    Fmask(i,3) = j
  end do
  !
  ! DG Book's method for getting Fmask
  !
 !Fmask(:,1) = pack( (/(i,i=1,nsp)/) , abs(rs(2,:) + one    ) < eps12 )
 !Fmask(:,2) = pack( (/(i,i=1,nsp)/) , abs(rs(1,:) + rs(2,:)) < eps12 )
 !Fmask(:,3) = pack( (/(i,i=1,nsp)/) , abs(rs(1,:) + one    ) < eps12 )
  !
  ! Compute the inverse mass matrix along an edge
  !
#ifdef SPECIAL_FOR_INTEL
  mmat1D_qp = invert_matrix_lapack95( matmul( vedge_qp , transpose(vedge_qp) ) )
#else
  mmat1D_qp = invert_matrix( matmul( vedge_qp , transpose(vedge_qp) ) )
#endif
  !
  ! Face 1
  !
  n1 = 1  ! index of first node on this face
  n2 = np ! index of last  node on this face
  !
  cmat_qp(Fmask(:,1),n1:n2) = mmat1D_qp
  !
  ! Face 2
  !
  n1 = np+1 ! index of first node on this face
  n2 = 2*np ! index of last  node on this face
  !
  cmat_qp(Fmask(:,2),n1:n2) = mmat1D_qp
  !
  ! Face 3
  !
  n1 = 2*np+1 ! index of first node on this face
  n2 = 3*np   ! index of last  node on this face
  !
  cmat_qp(Fmask(:,3),n1:n2) = mmat1D_qp
  !
  ! Finally, compute the correction matrix
  !
  cmat_qp = matmul( transpose(vtria_qp) , cmat_qp )
  cmat_qp = qtwo*matmul( vtria_qp , cmat_qp )
  !
  ! Use the chop function to copy cmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( cmat_qp )
  !
end function CorrectionMatrix_Pyra
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function CorrectionMatrix_Pris(rs,eval_deriv) result(return_value)
  !
  ! RS : r and s coordinates for the standard
  !      triangle element on the interval [-1,1]
  !
  !.. Use Statements ..
  use vandermonde_mod, only : vand
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rs
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: eval_deriv
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,n1,n2,np,nsp,nfp
  integer :: this_order,ierr
  logical(lk) :: compute_derivative
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: Fmask(:,:)
  real(qp), allocatable :: cmat_qp(:,:)
  real(qp), allocatable :: mmat1D_qp(:,:)
  real(qp), allocatable :: vedge_qp(:,:)
  real(qp), allocatable :: vtria_qp(:,:)
  !
continue
  !
  nsp = size(rs,dim=2)
  this_order = solpts2order(Geom_Pris,nsp)
  np = this_order+1
  nfp = 3*np
  !
  compute_derivative = true
  if (present(eval_deriv)) then
    compute_derivative = eval_deriv
  end if
  !
  ! Allocate the return array containing the working
  ! precision correction matrix
  !
  allocate ( return_value(1:nsp,1:nfp) , source=zero , stat=ierr )
  !
  ! Create the correction matrix using quad precision and allocate the
  ! inverse mass matrix and Fmask arrays.
  !
  allocate ( cmat_qp(1:nsp,1:nfp) , source=qzero , stat=ierr )
  allocate ( mmat1D_qp(1:np,1:np) , source=qzero , stat=ierr )
  allocate ( Fmask(1:np,1:3) , source=0 , stat=ierr )
  !
  ! Make quad precision copies of the Edge and Tria Vandermonde matrices.
  !
  allocate ( vedge_qp(1:np,1:np) , source=qzero , stat=ierr )
  allocate ( vtria_qp(1:nsp,1:nsp) , source=qzero , stat=ierr )
  !
  vedge_qp = real( vand(Geom_Edge,this_order)%modal2nodal%mat , kind=qp )
  vtria_qp = real( vand(Geom_Tria,this_order)%modal2nodal%mat , kind=qp )
  !
  ! Face 1 - bottom edge of triangle
  !
  do i = 1,np
    Fmask(i,1) = i
  end do
  !
  ! Face 2 - diagonal edge of triangle
  !
  j = 0
  do i = 1,np
    j = j + np + 1 - i
    Fmask(i,2) = j
  end do
  !
  ! Face 3 - vertical edge of triangle
  !
  j = nsp + 1
  do i = 1,np
    j = j - i
    Fmask(i,3) = j
  end do
  !
  ! DG Book's method for getting Fmask
  !
 !Fmask(:,1) = pack( (/(i,i=1,nsp)/) , abs(rs(2,:) + one    ) < eps12 )
 !Fmask(:,2) = pack( (/(i,i=1,nsp)/) , abs(rs(1,:) + rs(2,:)) < eps12 )
 !Fmask(:,3) = pack( (/(i,i=1,nsp)/) , abs(rs(1,:) + one    ) < eps12 )
  !
  ! Compute the inverse mass matrix along an edge
  !
#ifdef SPECIAL_FOR_INTEL
  mmat1D_qp = invert_matrix_lapack95( matmul( vedge_qp , transpose(vedge_qp) ) )
#else
  mmat1D_qp = invert_matrix( matmul( vedge_qp , transpose(vedge_qp) ) )
#endif
  !
  ! Face 1
  !
  n1 = 1  ! index of first node on this face
  n2 = np ! index of last  node on this face
  !
  cmat_qp(Fmask(:,1),n1:n2) = mmat1D_qp
  !
  ! Face 2
  !
  n1 = np+1 ! index of first node on this face
  n2 = 2*np ! index of last  node on this face
  !
  cmat_qp(Fmask(:,2),n1:n2) = mmat1D_qp
  !
  ! Face 3
  !
  n1 = 2*np+1 ! index of first node on this face
  n2 = 3*np   ! index of last  node on this face
  !
  cmat_qp(Fmask(:,3),n1:n2) = mmat1D_qp
  !
  ! Finally, compute the correction matrix
  !
  cmat_qp = matmul( transpose(vtria_qp) , cmat_qp )
  cmat_qp = qtwo*matmul( vtria_qp , cmat_qp )
  !
  ! Use the chop function to copy cmat_qp into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( cmat_qp )
  !
end function CorrectionMatrix_Pris
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function CorrectionMatrix_Hexa(xi,eval_deriv) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xi
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: eval_deriv
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: nf,nep,nfp,nbp,nsp
  integer :: i,j,k,l,ierr
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: cmat1D(:)
  !
continue
  !
  ! #########################################################################
  ! #########################################################################
  ! #########################################################################
  ! #########################################################################
  ! #####  NOTE: THE ORDER OF THE FLUX POINTS IS SUCH THAT THE LOWEST   #####
  ! #####        COORDINATE DIRECTION IS THE INNER MOST LOOP, i.e.,     #####
  ! #####          FOR AN X-FACE, INNER LOOP OVER Y, OUTER LOOP OVER Z  #####
  ! #####          FOR A  Y-FACE, INNER LOOP OVER X, OUTER LOOP OVER Z  #####
  ! #####          FOR A  Z-FACE, INNER LOOP OVER X, OUTER LOOP OVER Y  #####
  ! #########################################################################
  ! #########################################################################
  ! #########################################################################
  ! #########################################################################
  !
  ! ###########################################################################
  ! ###########################################################################
  ! ####  NOTE: No need to use chopping function in this procedure. The    ####
  ! ####        hexa correction matrix is really just the 1D correction    ####
  ! ####        matrix applied to each face/flux point of the standard     ####
  ! ####        hexa element. The 1D correction matrix is obtained by the  ####
  ! ####        function CorrectionMatrix_Edge which applies the chopping  ####
  ! ####        function to its result so there is no need to repeat this  ####
  ! ####        operation again.                                           ####
  ! ###########################################################################
  ! ###########################################################################
  !
  nep = size(xi)
  nfp = nep*nep
  nbp = 6*nfp
  nsp = nep*nep*nep
  !
  ! Allocate the return array containing the working
  ! precision correction matrix
  !
  allocate ( return_value(1:nsp,1:nbp) , source=zero , stat=ierr )
  !
  ! Create the 1D correction matrix
  !
  allocate ( cmat1D(1:nep) , source=zero , stat=ierr )
  !
  ! 1D Correction matrix
  !
  cmat1D = CorrectionMatrix_Edge(xi,eval_deriv)
  !
  ! ############################################################
  ! ############################################################
  ! #####  Face 1 - z/zeta min face of quad (bottom/down)  #####
  ! ############################################################
  ! ############################################################
  !
  nf = 0*nfp+1 ! First point on this face
  !
  ! Loops over points on this face
  !
  do j = 1,nep
    do i = 1,nep
      !
      ! Loop over interior solution points in z/zeta direction
      !
      do k = 1,nep
        l = i + (j-1)*nep + (k-1)*nfp
        return_value(l,nf) = cmat1D(k)
      end do
      !
      nf = nf + 1 ! Proceed to next point on this face
      !
    end do
  end do
  !
  ! #######################################################
  ! #######################################################
  ! #####  Face 2 - z/zeta max face of quad (top/up)  #####
  ! #######################################################
  ! #######################################################
  !
  nf = 1*nfp+1 ! First point on this face
  !
  ! Loops over points on this face
  !
  do j = 1,nep
    do i = 1,nep
      !
      ! Loop over interior solution points in z/zeta direction
      !
      do k = 1,nep
        l = i + (j-1)*nep + nfp*(nep-k)
        return_value(l,nf) = cmat1D(k)
      end do
      !
      nf = nf + 1 ! Proceed to next point on this face
      !
    end do
  end do
  !
  ! #########################################################
  ! #########################################################
  ! #####  Face 3 - x/xi max face of quad (right/east)  #####
  ! #########################################################
  ! #########################################################
  !
  nf = 2*nfp+1 ! First point on this face
  !
  ! Loops over points on this face
  !
  do k = 1,nep
    do j = 1,nep
      !
      ! Loop over interior solution points in x/xi direction
      !
      do i = 1,nep
        l = j*nep + (k-1)*nfp - i + 1
        return_value(l,nf) = cmat1D(i)
      end do
      !
      nf = nf + 1 ! Proceed to next point on this face
      !
    end do
  end do
  !
  ! ########################################################
  ! ########################################################
  ! #####  Face 4 - x/xi min face of quad (left/west)  #####
  ! ########################################################
  ! ########################################################
  !
  nf = 3*nfp+1 ! First point on this face
  !
  ! Loops over points on this face
  !
  do k = 1,nep
    do j = 1,nep
      !
      ! Loop over interior solution points in x/xi direction
      !
      do i = 1,nep
        l = (j-1)*nep + (k-1)*nfp + i
        return_value(l,nf) = cmat1D(i)
      end do
      !
      nf = nf + 1 ! Proceed to next point on this face
      !
    end do
  end do
  !
  ! ###########################################################
  ! ###########################################################
  ! #####  Face 5 - y/eta min face of quad (front/south)  #####
  ! ###########################################################
  ! ###########################################################
  !
  nf = 4*nfp+1 ! First point on this face
  !
  ! Loops over points on this face
  !
  do k = 1,nep
    do i = 1,nep
      !
      ! Loop over interior solution points in y/eta direction
      !
      do j = 1,nep
        l = (k-1)*nfp + i + (j-1)*nep
        return_value(l,nf) = cmat1D(j)
      end do
      !
      nf = nf + 1 ! Proceed to next point on this face
      !
    end do
  end do
  !
  ! ##########################################################
  ! ##########################################################
  ! #####  Face 6 - y/eta max face of quad (back/north)  #####
  ! ##########################################################
  ! ##########################################################
  !
  nf = 5*nfp+1 ! First point on this face
  !
  ! Loops over points on this face
  !
  do k = 1,nep
    do i = 1,nep
      !
      ! Loop over interior solution points in y/eta direction
      !
      do j = 1,nep
        l = (k-1)*nfp + i + (nep-j)*nep
        return_value(l,nf) = cmat1D(j)
      end do
      !
      nf = nf + 1 ! Proceed to next point on this face
      !
    end do
  end do
  !
end function CorrectionMatrix_Hexa
!
!###############################################################################
!
subroutine correction_memory_usage(iunit)
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
  ! Memory usage for the correct array of type matrix
  !
  if (allocated(correct)) then
    !
    call dt%set( storage_size(correct,kind=inttype) , &
                         size(correct,kind=inttype) )
    !
    do j = lbound(correct,dim=2),ubound(correct,dim=2)
      do i = lbound(correct,dim=1),ubound(correct,dim=1)
        !
        array = matrix_memory_usage( correct(i,j) )
        !
        write (array_name,4) "correct",Geom_Name(i),j
        write (iou,2) array_name,array%get_units()
        !
        call dt%add(array)
        !
      end do
    end do
    !
    call total%add(dt)
    write (array_name,5) "correct"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "correct"
    !
  end if
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module correction_mod")
  2 format (" Array: ",a,dt)
  3 format (" Total memory = ",f8.2,1x,a)
  4 format (a,"(",a,",",i0,")%mat")
  5 format ("'",a,"'")
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  !
#endif
end subroutine correction_memory_usage
!
!###############################################################################
!
subroutine move_correct_to_device()
  !
  !.. Local Scalars ..
  integer :: i,i1,i2,j,j1,j2,m1,m2,n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "move_correct_to_device"
  !
continue
  !
  if (allocated(correct)) then
    !$ACC ENTER DATA COPYIN(correct)
    i1 = lbound(correct,dim=1)
    i2 = ubound(correct,dim=1)
    j1 = lbound(correct,dim=2)
    j2 = ubound(correct,dim=2)
    !$ACC ENTER DATA COPYIN(correct(i1:i2,j1:j2))
    do j = j1,j2
      do i = i1,i2
        !
        if (allocated(correct(i,j)%mat)) then
          !$ACC ENTER DATA COPYIN(correct(i,j)%mat)
          m1 = lbound(correct(i,j)%mat,dim=1)
          m2 = ubound(correct(i,j)%mat,dim=1)
          n1 = lbound(correct(i,j)%mat,dim=2)
          n2 = ubound(correct(i,j)%mat,dim=2)
          !$ACC ENTER DATA COPYIN(correct(i,j)%mat(m1:m2,n1:n2))
        end if
        !
        if (allocated(correct(i,j)%row_ptr)) then
          !$ACC ENTER DATA COPYIN(correct(i,j)%row_ptr)
          m1 = lbound(correct(i,j)%row_ptr,dim=1)
          m2 = ubound(correct(i,j)%row_ptr,dim=1)
          !$ACC ENTER DATA COPYIN(correct(i,j)%row_ptr(m1:m2))
        end if
        !
        if (allocated(correct(i,j)%row_idx)) then
          !$ACC ENTER DATA COPYIN(correct(i,j)%row_idx)
          m1 = lbound(correct(i,j)%row_idx,dim=1)
          m2 = ubound(correct(i,j)%row_idx,dim=1)
          !$ACC ENTER DATA COPYIN(correct(i,j)%row_idx(m1:m2))
        end if
        !
        if (allocated(correct(i,j)%row_val)) then
          !$ACC ENTER DATA COPYIN(correct(i,j)%row_val)
          m1 = lbound(correct(i,j)%row_val,dim=1)
          m2 = ubound(correct(i,j)%row_val,dim=1)
          !$ACC ENTER DATA COPYIN(correct(i,j)%row_val(m1:m2))
        end if
        !
        if (allocated(correct(i,j)%col_ptr)) then
          !$ACC ENTER DATA COPYIN(correct(i,j)%col_ptr)
          m1 = lbound(correct(i,j)%col_ptr,dim=1)
          m2 = ubound(correct(i,j)%col_ptr,dim=1)
          !$ACC ENTER DATA COPYIN(correct(i,j)%col_ptr(m1:m2))
        end if
        !
        if (allocated(correct(i,j)%col_idx)) then
          !$ACC ENTER DATA COPYIN(correct(i,j)%col_idx)
          m1 = lbound(correct(i,j)%col_idx,dim=1)
          m2 = ubound(correct(i,j)%col_idx,dim=1)
          !$ACC ENTER DATA COPYIN(correct(i,j)%col_idx(m1:m2))
        end if
        !
        if (allocated(correct(i,j)%col_val)) then
          !$ACC ENTER DATA COPYIN(correct(i,j)%col_val)
          m1 = lbound(correct(i,j)%col_val,dim=1)
          m2 = ubound(correct(i,j)%col_val,dim=1)
          !$ACC ENTER DATA COPYIN(correct(i,j)%col_val(m1:m2))
        end if
        !
      end do
    end do
  end if
  !
end subroutine move_correct_to_device
!
!###############################################################################
!
include "Functions/solpts2order.f90"
!
!###############################################################################
!
end module correction_mod
!
!###############################################################################
!
