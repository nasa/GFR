pure function intseq_no_interval_i1(seq_beg,seq_end) result(return_value)
  !
  use iso_fortran_env, only : lp => int8
  !
  integer(lp), intent(in) :: seq_beg
  integer(lp), intent(in) :: seq_end
  !
  integer(lp) :: return_value(1:seq_end-seq_beg+1_lp)
  !
  integer(lp) :: n
  !
continue
  return_value = [(n,n=seq_beg,seq_end)]
end function intseq_no_interval_i1
!
!###############################################################################
!
pure function intseq_with_interval_i1(seq_beg,seq_end,seq_interval) &
                               result(return_value)
  !
  use iso_fortran_env, only : lp => int8
  !
  integer(lp), intent(in) :: seq_beg
  integer(lp), intent(in) :: seq_end
  integer(lp), intent(in) :: seq_interval
  !
  integer(lp), allocatable :: return_value(:)
  !
  integer(lp) :: loc_interval,n
  integer :: ierr
  !
continue
  loc_interval = 1_lp
  if (seq_interval /= 0_lp) loc_interval = seq_interval
  n = (seq_end-seq_beg)/loc_interval + 1_lp
  allocate ( return_value(1:n) , stat=ierr )
  return_value = [(n,n=seq_beg,seq_end,loc_interval)]
end function intseq_with_interval_i1
!
!###############################################################################
!
pure function intseq_rev_dir_i1(seq_end,rev_dir) result(return_value)
  !
  use iso_fortran_env, only : lp => int8
  !
  integer(lp), intent(in) :: seq_end
  logical(lk), intent(in) :: rev_dir
  !
  integer(lp) :: return_value(1:seq_end)
  !
  integer(lp) :: n
  !
continue
  if (rev_dir) then
    return_value = [(n,n=seq_end,1,-1)]
  else
    return_value = [(n,n=1,seq_end)]
  end if
end function intseq_rev_dir_i1
!
!###############################################################################
!
pure function intseq_no_interval_i2(seq_beg,seq_end) result(return_value)
  !
  use iso_fortran_env, only : lp => int16
  !
  integer(lp), intent(in) :: seq_beg
  integer(lp), intent(in) :: seq_end
  !
  integer(lp) :: return_value(1:seq_end-seq_beg+1_lp)
  !
  integer(lp) :: n
  !
continue
  return_value = [(n,n=seq_beg,seq_end)]
end function intseq_no_interval_i2
!
!###############################################################################
!
pure function intseq_with_interval_i2(seq_beg,seq_end,seq_interval) &
                               result(return_value)
  !
  use iso_fortran_env, only : lp => int16
  !
  integer(lp), intent(in) :: seq_beg
  integer(lp), intent(in) :: seq_end
  integer(lp), intent(in) :: seq_interval
  !
  integer(lp), allocatable :: return_value(:)
  !
  integer(lp) :: loc_interval,n
  integer :: ierr
  !
continue
  loc_interval = 1_lp
  if (seq_interval /= 0_lp) loc_interval = seq_interval
  n = (seq_end-seq_beg)/loc_interval + 1_lp
  allocate ( return_value(1:n) , stat=ierr )
  return_value = [(n,n=seq_beg,seq_end,loc_interval)]
end function intseq_with_interval_i2
!
!###############################################################################
!
pure function intseq_rev_dir_i2(seq_end,rev_dir) result(return_value)
  !
  use iso_fortran_env, only : lp => int16
  !
  integer(lp), intent(in) :: seq_end
  logical(lk), intent(in) :: rev_dir
  !
  integer(lp) :: return_value(1:seq_end)
  !
  integer(lp) :: n
  !
continue
  if (rev_dir) then
    return_value = [(n,n=seq_end,1,-1)]
  else
    return_value = [(n,n=1,seq_end)]
  end if
end function intseq_rev_dir_i2
!
!###############################################################################
!
pure function intseq_no_interval_i4(seq_beg,seq_end) result(return_value)
  !
  use iso_fortran_env, only : lp => int32
  !
  integer(lp), intent(in) :: seq_beg
  integer(lp), intent(in) :: seq_end
  !
  integer(lp) :: return_value(1:seq_end-seq_beg+1_lp)
  !
  integer(lp) :: n
  !
continue
  return_value = [(n,n=seq_beg,seq_end)]
end function intseq_no_interval_i4
!
!###############################################################################
!
pure function intseq_with_interval_i4(seq_beg,seq_end,seq_interval) &
                               result(return_value)
  !
  use iso_fortran_env, only : lp => int32
  !
  integer(lp), intent(in) :: seq_beg
  integer(lp), intent(in) :: seq_end
  integer(lp), intent(in) :: seq_interval
  !
  integer(lp), allocatable :: return_value(:)
  !
  integer(lp) :: loc_interval,n
  integer :: ierr
  !
continue
  loc_interval = 1_lp
  if (seq_interval /= 0_lp) loc_interval = seq_interval
  n = (seq_end-seq_beg)/loc_interval + 1_lp
  allocate ( return_value(1:n) , stat=ierr )
  return_value = [(n,n=seq_beg,seq_end,loc_interval)]
end function intseq_with_interval_i4
!
!###############################################################################
!
pure function intseq_rev_dir_i4(seq_end,rev_dir) result(return_value)
  !
  use iso_fortran_env, only : lp => int32
  !
  integer(lp), intent(in) :: seq_end
  logical(lk), intent(in) :: rev_dir
  !
  integer(lp) :: return_value(1:seq_end)
  !
  integer(lp) :: n
  !
continue
  if (rev_dir) then
    return_value = [(n,n=seq_end,1,-1)]
  else
    return_value = [(n,n=1,seq_end)]
  end if
end function intseq_rev_dir_i4
!
!###############################################################################
!
pure function intseq_no_interval_i8(seq_beg,seq_end) result(return_value)
  !
  use iso_fortran_env, only : lp => int64
  !
  integer(lp), intent(in) :: seq_beg
  integer(lp), intent(in) :: seq_end
  !
  integer(lp) :: return_value(1:seq_end-seq_beg+1_lp)
  !
  integer(lp) :: n
  !
continue
  return_value = [(n,n=seq_beg,seq_end)]
end function intseq_no_interval_i8
!
!###############################################################################
!
pure function intseq_with_interval_i8(seq_beg,seq_end,seq_interval) &
                               result(return_value)
  !
  use iso_fortran_env, only : lp => int64
  !
  integer(lp), intent(in) :: seq_beg
  integer(lp), intent(in) :: seq_end
  integer(lp), intent(in) :: seq_interval
  !
  integer(lp), allocatable :: return_value(:)
  !
  integer(lp) :: loc_interval,n
  integer :: ierr
  !
continue
  loc_interval = 1_lp
  if (seq_interval /= 0_lp) loc_interval = seq_interval
  n = (seq_end-seq_beg)/loc_interval + 1_lp
  allocate ( return_value(1:n) , stat=ierr )
  return_value = [(n,n=seq_beg,seq_end,loc_interval)]
end function intseq_with_interval_i8
!
!###############################################################################
!
pure function intseq_rev_dir_i8(seq_end,rev_dir) result(return_value)
  !
  use iso_fortran_env, only : lp => int64
  !
  integer(lp), intent(in) :: seq_end
  logical(lk), intent(in) :: rev_dir
  !
  integer(lp) :: return_value(1:seq_end)
  !
  integer(lp) :: n
  !
continue
  if (rev_dir) then
    return_value = [(n,n=seq_end,1,-1)]
  else
    return_value = [(n,n=1,seq_end)]
  end if
end function intseq_rev_dir_i8
!
!###############################################################################
!
