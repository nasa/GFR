pure function intseq_no_interval(seq_beg,seq_end) result(return_value)
  integer, intent(in) :: seq_beg
  integer, intent(in) :: seq_end
  integer :: return_value(1:seq_end-seq_beg+1)
  integer :: n
continue
  return_value = [(n,n=seq_beg,seq_end)]
end function intseq_no_interval
pure function intseq_with_interval(seq_beg,seq_end,seq_interval) &
                            result(return_value)
  integer, intent(in) :: seq_beg
  integer, intent(in) :: seq_end
  integer, intent(in) :: seq_interval
  integer, allocatable :: return_value(:)
  integer :: iint,ierr,n
continue
  iint = 1
  if (seq_interval /= 0) iint = seq_interval
  n = (seq_end-seq_beg)/iint + 1  
  allocate ( return_value(1:n) , stat=ierr )
  return_value = [(n,n=seq_beg,seq_end,iint)]
end function intseq_with_interval
