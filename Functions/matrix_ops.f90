pure function invert_matrix_SP(a) result(return_value)
  !
  integer, parameter :: lp = sp
  !
  real(lp), dimension(:,:), intent(in) :: a
  !
  real(lp), dimension(size(a,1),size(a,1)) :: return_value
  !
  integer  :: np,i,ii,j,k,ll,imax
  real(lp) :: summ
  !
  integer,  dimension(size(a,1))           :: indx
  real(lp), dimension(size(a,1))           :: temp
  real(lp), dimension(size(a,1))           :: vv
  real(lp), dimension(size(a,1),size(a,1)) :: b
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one = 1.0_lp
  real(lp), parameter :: rndoff = epsilon(zero)
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
        summ = summ - dot(b(i,ii:i-1),return_value(ii:i-1,k))
      else if (summ /= zero) then
        ii = i
      end if
      return_value(i,k) = summ
    end do
    !
    do i = np,1,-1
      summ = dot(b(i,i+1:np),return_value(i+1:np,k))
      return_value(i,k) = (return_value(i,k) - summ)/b(i,i)
    end do
    !
  end do
  !
end function invert_matrix_SP
!
!###############################################################################
!
pure function invert_matrix_DP(a) result(return_value)
  !
  integer, parameter :: lp = dp
  !
  real(lp), dimension(:,:), intent(in) :: a
  !
  real(lp), dimension(size(a,1),size(a,1)) :: return_value
  !
  integer  :: np,i,ii,j,k,ll,imax
  real(lp) :: summ
  !
  integer,  dimension(size(a,1))           :: indx
  real(lp), dimension(size(a,1))           :: temp
  real(lp), dimension(size(a,1))           :: vv
  real(lp), dimension(size(a,1),size(a,1)) :: b
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one = 1.0_lp
  real(lp), parameter :: rndoff = epsilon(zero)
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
        summ = summ - dot(b(i,ii:i-1),return_value(ii:i-1,k))
      else if (summ /= zero) then
        ii = i
      end if
      return_value(i,k) = summ
    end do
    !
    do i = np,1,-1
      summ = dot(b(i,i+1:np),return_value(i+1:np,k))
      return_value(i,k) = (return_value(i,k) - summ)/b(i,i)
    end do
    !
  end do
  !
end function invert_matrix_DP
!
!###############################################################################
!
pure function invert_matrix_QP(a) result(return_value)
  !
  integer, parameter :: lp = qp
  !
  real(lp), dimension(:,:), intent(in) :: a
  !
  real(lp), dimension(size(a,1),size(a,1)) :: return_value
  !
  integer  :: np,i,ii,j,k,ll,imax
  real(lp) :: summ
  !
  integer,  dimension(size(a,1))           :: indx
  real(lp), dimension(size(a,1))           :: temp
  real(lp), dimension(size(a,1))           :: vv
  real(lp), dimension(size(a,1),size(a,1)) :: b
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one = 1.0_lp
  real(lp), parameter :: rndoff = epsilon(zero)
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
end function invert_matrix_QP
!
!###############################################################################
!
pure function invert_matrix90_SP(a) result(return_value)
  !
  integer, parameter :: lp = sp
  !
  real(lp), dimension(:,:), intent(in) :: a
  !
  real(lp), dimension(size(a,1),size(a,1)) :: return_value
  !
  integer   :: np,i,k
  real(lp) :: d
  !
  integer,  dimension(size(a,1))           :: indx
  real(lp), dimension(size(a,1),size(a,1)) :: b
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
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
end function invert_matrix90_SP
!
!###############################################################################
!
pure function invert_matrix90_DP(a) result(return_value)
  !
  integer, parameter :: lp = dp
  !
  real(lp), dimension(:,:), intent(in) :: a
  !
  real(lp), dimension(size(a,1),size(a,1)) :: return_value
  !
  integer   :: np,i,k
  real(lp) :: d
  !
  integer,  dimension(size(a,1))           :: indx
  real(lp), dimension(size(a,1),size(a,1)) :: b
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
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
end function invert_matrix90_DP
!
!###############################################################################
!
pure function invert_matrix90_QP(a) result(return_value)
  !
  integer, parameter :: lp = qp
  !
  real(lp), dimension(:,:), intent(in) :: a
  !
  real(lp), dimension(size(a,1),size(a,1)) :: return_value
  !
  integer   :: np,i,k
  real(lp) :: d
  !
  integer,  dimension(size(a,1))           :: indx
  real(lp), dimension(size(a,1),size(a,1)) :: b
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
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
end function invert_matrix90_QP
!
!###############################################################################
!
pure subroutine ludcmp_SP(a,np,indx,d)
  !
  integer, parameter :: lp = sp
  !
  integer,                           intent(in) :: np
  integer,  dimension(1:np),        intent(out) :: indx
  real(lp),                         intent(out) :: d
  real(lp), dimension(1:np,1:np), intent(inout) :: a
  !
  integer :: j,imax
  !
  real(lp), dimension(1:np) :: vv,temp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one = 1.0_lp
  real(lp), parameter :: rndoff = epsilon(zero)
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
    a(j+1:np,j+1:np) = a(j+1:np,j+1:np) - &
                       spread(a(j+1:np,j),dim=2,ncopies=np-j)* &
                       spread(a(j,j+1:np),dim=1,ncopies=np-j)
    !
  end do
  !
end subroutine ludcmp_SP
!
!###############################################################################
!
pure subroutine ludcmp_DP(a,np,indx,d)
  !
  integer, parameter :: lp = dp
  !
  integer,                           intent(in) :: np
  integer,  dimension(1:np),        intent(out) :: indx
  real(lp),                         intent(out) :: d
  real(lp), dimension(1:np,1:np), intent(inout) :: a
  !
  integer :: j,imax
  !
  real(lp), dimension(1:np) :: vv,temp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one = 1.0_lp
  real(lp), parameter :: rndoff = epsilon(zero)
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
    a(j+1:np,j+1:np) = a(j+1:np,j+1:np) - &
                       spread(a(j+1:np,j),dim=2,ncopies=np-j)* &
                       spread(a(j,j+1:np),dim=1,ncopies=np-j)
    !
  end do
  !
end subroutine ludcmp_DP
!
!###############################################################################
!
pure subroutine ludcmp_QP(a,np,indx,d)
  !
  integer, parameter :: lp = qp
  !
  integer,                           intent(in) :: np
  integer,  dimension(1:np),        intent(out) :: indx
  real(lp),                         intent(out) :: d
  real(lp), dimension(1:np,1:np), intent(inout) :: a
  !
  integer :: j,imax
  !
  real(lp), dimension(1:np) :: vv,temp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one = 1.0_lp
  real(lp), parameter :: rndoff = epsilon(zero)
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
    a(j+1:np,j+1:np) = a(j+1:np,j+1:np) - &
                       spread(a(j+1:np,j),dim=2,ncopies=np-j)* &
                       spread(a(j,j+1:np),dim=1,ncopies=np-j)
    !
  end do
  !
end subroutine ludcmp_QP
!
!###############################################################################
!
pure subroutine lubksb_SP(a,np,indx,y)
  !
  integer, parameter :: lp = sp
  !
  integer,                        intent(in) :: np
  integer,  dimension(1:np),      intent(in) :: indx
  real(lp), dimension(1:np,1:np), intent(in) :: a
  real(lp), dimension(1:np),   intent(inout) :: y
  !
  integer :: i,ii,ll
  !
  real(lp) :: summ
  !
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  ii = 0
  do i = 1,np
    ll = indx(i)
    summ = y(ll)
    y(ll) = y(i)
    if (ii /= 0) then
      summ = summ - dot(a(i,ii:i-1),y(ii:i-1))
    else if (summ /= zero) then
      ii = i
    end if
    y(i) = summ
  end do
  !
  do i = np,1,-1
    y(i) = (y(i) - dot(a(i,i+1:np),y(i+1:np)))/a(i,i)
  end do
  !
end subroutine lubksb_SP
!
!###############################################################################
!
pure subroutine lubksb_DP(a,np,indx,y)
  !
  integer, parameter :: lp = dp
  !
  integer,                        intent(in) :: np
  integer,  dimension(1:np),      intent(in) :: indx
  real(lp), dimension(1:np,1:np), intent(in) :: a
  real(lp), dimension(1:np),   intent(inout) :: y
  !
  integer :: i,ii,ll
  !
  real(lp) :: summ
  !
  real(lp), parameter :: zero = 0.0_lp
  !
continue
  !
  ii = 0
  do i = 1,np
    ll = indx(i)
    summ = y(ll)
    y(ll) = y(i)
    if (ii /= 0) then
      summ = summ - dot(a(i,ii:i-1),y(ii:i-1))
    else if (summ /= zero) then
      ii = i
    end if
    y(i) = summ
  end do
  !
  do i = np,1,-1
    y(i) = (y(i) - dot(a(i,i+1:np),y(i+1:np)))/a(i,i)
  end do
  !
end subroutine lubksb_DP
!
!###############################################################################
!
pure subroutine lubksb_QP(a,np,indx,y)
  !
  integer, parameter :: lp = qp
  !
  integer,                        intent(in) :: np
  integer,  dimension(1:np),      intent(in) :: indx
  real(lp), dimension(1:np,1:np), intent(in) :: a
  real(lp), dimension(1:np),   intent(inout) :: y
  !
  integer :: i,ii,ll
  !
  real(lp) :: summ
  !
  real(lp), parameter :: zero = 0.0_lp
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
end subroutine lubksb_QP
!
!###############################################################################
!
pure function determinant_SP(a) result(return_value)
  !
  integer, parameter :: lp = sp
  !
  real(lp), dimension(:,:), intent(in) :: a
  !
  real(lp) :: return_value
  !
  integer :: i,n
  !
  integer,   dimension(size(a,1))           :: indx
  real(lp), dimension(size(a,1),size(a,1)) :: b
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
end function determinant_SP
!
!###############################################################################
!
pure function determinant_DP(a) result(return_value)
  !
  integer, parameter :: lp = dp
  !
  real(lp), dimension(:,:), intent(in) :: a
  !
  real(lp) :: return_value
  !
  integer :: i,n
  !
  integer,   dimension(size(a,1))           :: indx
  real(lp), dimension(size(a,1),size(a,1)) :: b
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
end function determinant_DP
!
!###############################################################################
!
pure function determinant_QP(a) result(return_value)
  !
  integer, parameter :: lp = qp
  !
  real(lp), dimension(:,:), intent(in) :: a
  !
  real(lp) :: return_value
  !
  integer :: i,n
  !
  integer,   dimension(size(a,1))           :: indx
  real(lp), dimension(size(a,1),size(a,1)) :: b
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
end function determinant_QP
!
!###############################################################################
!
pure function identity_matrix_dim_SP(a,n1,opt_n2) result(return_value)
  !
  integer, parameter :: lp = SP
  !
  real(lp), intent(in) :: a
  integer, intent(in) :: n1
  integer, optional, intent(in) :: opt_n2
  !
  real(lp), allocatable :: return_value(:,:)
  !
  integer :: i,j,n2,ierr
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one = 1.0_lp
  !
continue
  !
  n2 = n1
  if (present(opt_n2)) n2 = opt_n2
  !
  allocate ( return_value(1:n1,1:n2) , source=zero , stat=ierr )
  !
  do i = 1,min(n1,n2)
    return_value(i,i) = one
  end do
  !
end function identity_matrix_dim_SP
!
!###############################################################################
!
pure function identity_matrix_dim_DP(a,n1,opt_n2) result(return_value)
  !
  integer, parameter :: lp = DP
  !
  real(lp), intent(in) :: a
  integer, intent(in) :: n1
  integer, optional, intent(in) :: opt_n2
  !
  real(lp), allocatable :: return_value(:,:)
  !
  integer :: i,j,n2,ierr
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one = 1.0_lp
  !
continue
  !
  n2 = n1
  if (present(opt_n2)) n2 = opt_n2
  !
  allocate ( return_value(1:n1,1:n2) , source=zero , stat=ierr )
  !
  do i = 1,min(n1,n2)
    return_value(i,i) = one
  end do
  !
end function identity_matrix_dim_DP
!
!###############################################################################
!
pure function identity_matrix_dim_QP(a,n1,opt_n2) result(return_value)
  !
  integer, parameter :: lp = QP
  !
  real(lp), intent(in) :: a
  integer, intent(in) :: n1
  integer, optional, intent(in) :: opt_n2
  !
  real(lp), allocatable :: return_value(:,:)
  !
  integer :: i,j,n2,ierr
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one = 1.0_lp
  !
continue
  !
  n2 = n1
  if (present(opt_n2)) n2 = opt_n2
  !
  allocate ( return_value(1:n1,1:n2) , source=zero , stat=ierr )
  !
  do i = 1,min(n1,n2)
    return_value(i,i) = one
  end do
  !
end function identity_matrix_dim_QP
!
!###############################################################################
!
pure function identity_matrix_mat_SP(array) result(return_value)
  !
  integer, parameter :: lp = SP
  !
  real(lp), dimension(:,:), intent(in) :: array
  !
  real(lp), dimension(1:size(array,dim=1),1:size(array,dim=2)) :: return_value
  !
  integer :: i,n1,n2
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one = 1.0_lp
  !
continue
  !
  n1 = size(array,dim=1)
  n2 = size(array,dim=2)
  !
  return_value = zero
  !
  do i = 1,min(n1,n2)
    return_value(i,i) = one
  end do
  !
end function identity_matrix_mat_SP
!
!###############################################################################
!
pure function identity_matrix_mat_DP(array) result(return_value)
  !
  integer, parameter :: lp = DP
  !
  real(lp), dimension(:,:), intent(in) :: array
  !
  real(lp), dimension(1:size(array,dim=1),1:size(array,dim=2)) :: return_value
  !
  integer :: i,n1,n2
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one = 1.0_lp
  !
continue
  !
  n1 = size(array,dim=1)
  n2 = size(array,dim=2)
  !
  return_value = zero
  !
  do i = 1,min(n1,n2)
    return_value(i,i) = one
  end do
  !
end function identity_matrix_mat_DP
!
!###############################################################################
!
pure function identity_matrix_mat_QP(array) result(return_value)
  !
  integer, parameter :: lp = QP
  !
  real(lp), dimension(:,:), intent(in) :: array
  !
  real(lp), dimension(1:size(array,dim=1),1:size(array,dim=2)) :: return_value
  !
  integer :: i,n1,n2
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one = 1.0_lp
  !
continue
  !
  n1 = size(array,dim=1)
  n2 = size(array,dim=2)
  !
  return_value = zero
  !
  do i = 1,min(n1,n2)
    return_value(i,i) = one
  end do
  !
end function identity_matrix_mat_QP
!
!###############################################################################
!
