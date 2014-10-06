pure function invert_matrix_DP(a) result(return_value)
  !
  real(dp), dimension(:,:), intent(in) :: a
  !
  real(dp), dimension(size(a,1),size(a,1)) :: return_value
  !
  integer   :: np,i,ii,j,k,ll,imax
  real(dp) :: summ
  !
  integer,   dimension(size(a,1))           :: indx
  real(dp), dimension(size(a,1))           :: temp
  real(dp), dimension(size(a,1))           :: vv
  real(dp), dimension(size(a,1),size(a,1)) :: b
  !
  integer, parameter :: lp = dp
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
  real(qp), dimension(:,:), intent(in) :: a
  !
  real(qp), dimension(size(a,1),size(a,1)) :: return_value
  !
  integer   :: np,i,ii,j,k,ll,imax
  real(qp) :: summ
  !
  integer,   dimension(size(a,1))           :: indx
  real(qp), dimension(size(a,1))           :: temp
  real(qp), dimension(size(a,1))           :: vv
  real(qp), dimension(size(a,1),size(a,1)) :: b
  !
  integer, parameter :: lp = qp
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
pure function invert_matrix90_DP(a) result(return_value)
  !
  real(dp), dimension(:,:), intent(in) :: a
  !
  real(dp), dimension(size(a,1),size(a,1)) :: return_value
  !
  integer   :: np,i,k
  real(dp) :: d
  !
  integer,   dimension(size(a,1))           :: indx
  real(dp), dimension(size(a,1),size(a,1)) :: b
  !
  integer, parameter :: lp = dp
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
  real(qp), dimension(:,:), intent(in) :: a
  !
  real(qp), dimension(size(a,1),size(a,1)) :: return_value
  !
  integer   :: np,i,k
  real(qp) :: d
  !
  integer,   dimension(size(a,1))           :: indx
  real(qp), dimension(size(a,1),size(a,1)) :: b
  !
  integer, parameter :: lp = qp
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
pure subroutine ludcmp_DP(a,np,indx,d)
  !
  integer,                            intent(in) :: np
  integer,   dimension(1:np),        intent(out) :: indx
  real(dp),                         intent(out) :: d
  real(dp), dimension(1:np,1:np), intent(inout) :: a
  !
  integer :: j,imax
  !
  real(dp), dimension(1:np) :: vv,temp
  !
  integer, parameter :: lp = dp
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
    a(j+1:np,j+1:np) = a(j+1:np,j+1:np) - spread(a(j+1:np,j),dim=2,ncopies=np-j)* &
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
  integer,                            intent(in) :: np
  integer,   dimension(1:np),        intent(out) :: indx
  real(qp),                         intent(out) :: d
  real(qp), dimension(1:np,1:np), intent(inout) :: a
  !
  integer :: j,imax
  !
  real(qp), dimension(1:np) :: vv,temp
  !
  integer, parameter :: lp = qp
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
    a(j+1:np,j+1:np) = a(j+1:np,j+1:np) - spread(a(j+1:np,j),dim=2,ncopies=np-j)* &
                                          spread(a(j,j+1:np),dim=1,ncopies=np-j)
    !
  end do
  !
end subroutine ludcmp_QP
!
!###############################################################################
!
pure subroutine lubksb_DP(a,np,indx,y)
  !
  integer,                         intent(in) :: np
  integer,   dimension(1:np),      intent(in) :: indx
  real(dp), dimension(1:np,1:np), intent(in) :: a
  real(dp), dimension(1:np),   intent(inout) :: y
  !
  integer :: i,ii,ll
  !
  real(dp) :: summ
  !
  integer, parameter :: lp = dp
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
  integer,                         intent(in) :: np
  integer,   dimension(1:np),      intent(in) :: indx
  real(qp), dimension(1:np,1:np), intent(in) :: a
  real(qp), dimension(1:np),   intent(inout) :: y
  !
  integer :: i,ii,ll
  !
  real(qp) :: summ
  !
  integer, parameter :: lp = qp
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
pure function determinant_DP(a) result(return_value)
  !
  real(dp), dimension(:,:), intent(in) :: a
  !
  real(dp) :: return_value
  !
  integer :: i,n
  !
  integer,   dimension(size(a,1))           :: indx
  real(dp), dimension(size(a,1),size(a,1)) :: b
  !
  integer, parameter :: lp = dp
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
  real(qp), dimension(:,:), intent(in) :: a
  !
  real(qp) :: return_value
  !
  integer :: i,n
  !
  integer,   dimension(size(a,1))           :: indx
  real(qp), dimension(size(a,1),size(a,1)) :: b
  !
  integer, parameter :: lp = qp
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
