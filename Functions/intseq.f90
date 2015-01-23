pure function intseq(seq_beg,seq_end,seq_interval) result(return_value)
  integer, intent(in) :: seq_beg
  integer, intent(in) :: seq_end
  integer, optional, intent(in) :: seq_interval
  integer, allocatable :: return_value(:)
  integer :: iint,ierr,n
continue
  iint = 1
  if (present(seq_interval)) then
    if (seq_interval /= 0) iint = seq_interval
  end if
  n = (seq_end-seq_beg)/iint + 1  
  allocate ( return_value(1:n) , stat=ierr )
  return_value = [(n,n=seq_beg,seq_end,iint)]
end function intseq
