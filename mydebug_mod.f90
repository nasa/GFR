module mydebug_mod

! Module to help with debugging.
  use module_kind_types
  implicit none

contains

!##########################################################################
  subroutine mydebug(caller_name,line_num,file_nam,query)

! Debug the global solution variables.

  use eqn_idx, only :  nq
  use geovar,  only :  nnode, ncell, nface, nfbnd, nr
  use geovar,  only :  grid, xyz_nodes, xyz
  use flowvar, only :  facexyz, faceusp, sol_nodes, usp
  implicit none

  !.. Formal Arguments ..
  character(len=*),           intent(in) :: caller_name
  integer,                    intent(in) :: line_num
  character(len=*),           intent(in) :: file_nam
  character(len=*), optional, intent(in) :: query

  !.. Local ..
  integer :: n, nf
  real, dimension(1:nq) :: a_minval, a_maxval
  character(len=:) ,allocatable :: local_query

  ! Define a default list of variables to query.
  ! Override with optional list from arguments, if present.
  local_query = ''  ! to prevent gfortran warning
  local_query = '|nnode|ncell|nface|nfbnd|' // &
                '|facexyz|xyz_nodes|grid%xyz|xyz|' // &
                '|faceusp|sol_nodes|grid%sol|usp|'
  if ( present(query) ) local_query = query

  write(iout,*) 'DEBUG_BEG:[mydebug]'
  write(iout,*) '  Called from subroutine ',trim(caller_name), &
                ', line ',line_num,' in file ',trim(file_nam)

! Print info about nnode.
  if ( index(local_query,'|nnode|') .gt. 0 ) then
    write(iout,*) '  nnode     = ',nnode
  endif

! Print info about ncell.
  if ( index(local_query,'|ncell|') .gt. 0 ) then
    write(iout,*) '  ncell     = ',ncell
  endif

! Print info about ncell.
  if ( index(local_query,'|nface|') .gt. 0 ) then
    write(iout,*) '  nface     = ',nface
  endif

! Print info about ncell.
  if ( index(local_query,'|nfbnd|') .gt. 0 ) then
    write(iout,*) '  nfbnd     = ',nfbnd
  endif

! Print info about facexyz.
  if ( index(local_query,'|facexyz|') .gt. 0 ) then
  if ( allocated( facexyz ) ) then
    write(iout,*) '  facexyz   =   allocated to size:', size(facexyz)
    a_minval =  huge(1.0)
    a_maxval = -huge(1.0)
    do nf=1,size(facexyz)     ! size should equal nface
      if ( allocated( facexyz(nf)%v ) ) then
        do n=1,nr
          a_minval(n) = min(a_minval(n), minval(facexyz(nf)%v(n,:)))
          a_maxval(n) = max(a_maxval(n), minval(facexyz(nf)%v(n,:)))
        enddo
      endif
    enddo
    write(iout,*) '  n, min, max = '
    do n=1,nr
      write(iout,*) n,a_minval(n),a_maxval(n)
    enddo
  else
    write(iout,*) '  facexyz  = unallocated'
  endif
  endif

! Print info about xyz_nodes.
  if ( index(local_query,'|xyz_nodes|') .gt. 0 ) then
  if ( allocated( xyz_nodes ) ) then
    write(iout,*) '  xyz_nodes =   allocated to size:', &
                  size(xyz_nodes,dim=1), size(xyz_nodes,dim=2)
    write(iout,*) '  n, min, max = '
    do n=1,size(xyz_nodes,dim=1)
      write(iout,*) n,minval(xyz_nodes(n,:)),maxval(xyz_nodes(n,:))
    enddo
  else
    write(iout,*) '  xyz_nodes = unallocated'
  endif
  endif

! Print info about grid%xyz.
  if ( index(local_query,'|grid%xyz|') .gt. 0 ) then
  if ( allocated( grid ) ) then
    write(iout,*) '  grid      =   allocated'
    if ( allocated( grid%xyz ) ) then
      write(iout,*) '  grid%xyz  =   allocated to size:', &
                    size(grid%xyz,dim=1), size(grid%xyz,dim=2)
      write(iout,*) '  n, min, max = '
      do n=1,size(grid%xyz,dim=1)
        write(iout,*) n,minval(grid%xyz(n,:)),maxval(grid%xyz(n,:))
      enddo
    else
      write(iout,*) '  grid%xyz  = unallocated'
    endif
  else
    write(iout,*) '  grid      = unallocated'
  endif
  endif

! Print info about xyz.
  if ( index(local_query,'|xyz|') .gt. 0 ) then
  if ( allocated( xyz ) ) then
    write(iout,*) '  xyz       =   allocated to size ', &
                     size(xyz,dim=1), size(xyz,dim=2)
    write(iout,*) '  n, min, max = '
    do n=1,size(xyz,dim=1)
      write(iout,*) n,minval(xyz(n,:)),maxval(xyz(n,:))
    enddo
  else
    write(iout,*) '  xyz       = unallocated'
  endif
  endif



! Print info about faceusp.
  if ( index(local_query,'|faceusp|') .gt. 0 ) then
  if ( allocated( faceusp ) ) then
    write(iout,*) '  faceusp   =   allocated to size:', size(faceusp)
    a_minval =  huge(1.0)
    a_maxval = -huge(1.0)
    do nf=1,size(faceusp)     ! size should equal nface
      if ( allocated( faceusp(nf)%v ) ) then
        do n=1,nq
          a_minval(n) = min(a_minval(n), minval(faceusp(nf)%v(n,:,:)))
          a_maxval(n) = max(a_maxval(n), minval(faceusp(nf)%v(n,:,:)))
        enddo
      endif
    enddo
    write(iout,*) '  n, min, max = '
    do n=1,nq
      write(iout,*) n,a_minval(n),a_maxval(n)
    enddo
  else
    write(iout,*) '  faceusp   = unallocated'
  endif
  endif

! Print info about sol_nodes.
  if ( index(local_query,'|sol_nodes|') .gt. 0 ) then
  if ( allocated( sol_nodes ) ) then
    write(iout,*) '  sol_nodes =   allocated to size:', &
                  size(sol_nodes,dim=1), size(sol_nodes,dim=2)
    write(iout,*) '  n, min, max = '
    do n=1,nq
      write(iout,*) n,minval(sol_nodes(n,:)),maxval(sol_nodes(n,:))
    enddo
  else
    write(iout,*) '  sol_nodes = unallocated'
  endif
  endif

! Print info about grid%sol.
  if ( index(local_query,'|grid%sol|') .gt. 0 ) then
  if ( allocated( grid ) ) then
    write(iout,*) '  grid      =   allocated'
    if ( allocated( grid%sol ) ) then
      write(iout,*) '  grid%sol  =   allocated to size:', &
                    size(grid%sol,dim=1), size(grid%sol,dim=2)
      write(iout,*) '  n, min, max = '
      do n=1,nq
        write(iout,*) n,minval(grid%sol(n,:)),maxval(grid%sol(n,:))
      enddo
    else
      write(iout,*) '  grid%sol  = unallocated'
    endif
  else
    write(iout,*) '  grid      = unallocated'
  endif
  endif

! Print info about usp.
  if ( index(local_query,'|usp|') .gt. 0 ) then
  if ( allocated( usp ) ) then
    write(iout,*) '  usp       =   allocated to size ', &
                  size(usp,dim=1), size(usp,dim=2)
    write(iout,*) '  n, min, max = '
    do n=1,nq
      write(iout,*) n,minval(usp(n,:)),maxval(usp(n,:))
    enddo
  else
    write(iout,*) '  usp       = unallocated'
  endif
  endif


  write(iout,*) 'DEBUG_END:[mydebug]'
  flush(iout)

! Print same info for g2_* variables.
  call mydebug_g2(caller_name,line_num,file_nam,local_query)

  return
  end subroutine mydebug
!##########################################################################
  subroutine mydebug_g2(caller_name,line_num,file_nam,query)

! Debug the global solution variables for g2 grid.

  use eqn_idx, only :  nq
  use geovar,           only :  nr
  use interpsolvar_mod, only :  g2_nnode, g2_ncell, g2_nface, g2_nfbnd
  use interpsolvar_mod, only :  g2_grid, g2_xyz_nodes, g2_xyz
  use interpsolvar_mod, only :  g2_facexyz, g2_faceusp, g2_sol_nodes, g2_usp
  implicit none

  !.. Formal Arguments ..
  character(len=*),           intent(in) :: caller_name
  integer,                    intent(in) :: line_num
  character(len=*),           intent(in) :: file_nam
  character(len=*), optional, intent(in) :: query

  !.. Local ..
  integer :: n, nf
  real, dimension(1:nq) :: a_minval, a_maxval
  character(len=:) ,allocatable :: local_query

  ! Define a default list of variables to query.
  ! Override with optional list from arguments, if present.
  local_query = ''  ! to prevent gfortran warning
  local_query = '|nnode|ncell|nface|nfbnd|' // &
                '|facexyz|xyz_nodes|grid%xyz|xyz|' // &
                '|faceusp|sol_nodes|grid%sol|usp|'
  if ( present(query) ) local_query = query

  write(iout,*) 'DEBUG_BEG:[mydebug_g2]'
  write(iout,*) '  Called from subroutine ',trim(caller_name), &
                ', line ',line_num,' in file ',trim(file_nam)

! Print info about nnode.
  if ( index(local_query,'|nnode|') .gt. 0 ) then
    write(iout,*) '  g2_nnode     = ',g2_nnode
  endif

! Print info about ncell.
  if ( index(local_query,'|ncell|') .gt. 0 ) then
    write(iout,*) '  g2_ncell     = ',g2_ncell
  endif

! Print info about ncell.
  if ( index(local_query,'|nface|') .gt. 0 ) then
    write(iout,*) '  g2_nface     = ',g2_nface
  endif

! Print info about ncell.
  if ( index(local_query,'|nfbnd|') .gt. 0 ) then
    write(iout,*) '  g2_nfbnd     = ',g2_nfbnd
  endif

! Print info about facexyz.
  if ( index(local_query,'|facexyz|') .gt. 0 ) then
  if ( allocated( g2_facexyz ) ) then
    write(iout,*) '  g2_facexyz   =   allocated to size:', size(g2_facexyz)
    a_minval =  huge(1.0)
    a_maxval = -huge(1.0)
    do nf=1,size(g2_facexyz)     ! size should equal nface
      if ( allocated( g2_facexyz(nf)%v ) ) then
        do n=1,nr
          a_minval(n) = min(a_minval(n), minval(g2_facexyz(nf)%v(n,:)))
          a_maxval(n) = max(a_maxval(n), minval(g2_facexyz(nf)%v(n,:)))
        enddo
      endif
    enddo
    write(iout,*) '  n, min, max = '
    do n=1,nr
      write(iout,*) n,a_minval(n),a_maxval(n)
    enddo
  else
    write(iout,*) '  g2_facexyz  = unallocated'
  endif
  endif

! Print info about xyz_nodes.
  if ( index(local_query,'|xyz_nodes|') .gt. 0 ) then
  if ( allocated( g2_xyz_nodes ) ) then
    write(iout,*) '  g2_xyz_nodes =   allocated to size:', &
                  size(g2_xyz_nodes,dim=1), size(g2_xyz_nodes,dim=2)
    write(iout,*) '  n, min, max = '
    do n=1,size(g2_xyz_nodes,dim=1)
      write(iout,*) n,minval(g2_xyz_nodes(n,:)),maxval(g2_xyz_nodes(n,:))
    enddo
  else
    write(iout,*) '  g2_xyz_nodes = unallocated'
  endif
  endif

! Print info about grid%xyz.
  if ( index(local_query,'|grid%xyz|') .gt. 0 ) then
  if ( allocated( g2_grid ) ) then
    write(iout,*) '  g2_grid      =   allocated'
    if ( allocated( g2_grid%xyz ) ) then
      write(iout,*) '  g2_grid%xyz  =   allocated to size:', &
                    size(g2_grid%xyz,dim=1), size(g2_grid%xyz,dim=2)
      write(iout,*) '  n, min, max = '
      do n=1,size(g2_grid%xyz,dim=1)
        write(iout,*) n,minval(g2_grid%xyz(n,:)),maxval(g2_grid%xyz(n,:))
      enddo
    else
      write(iout,*) '  g2_grid%xyz  = unallocated'
    endif
  else
    write(iout,*) '  g2_grid      = unallocated'
  endif
  endif

! Print info about xyz.
  if ( index(local_query,'|xyz|') .gt. 0 ) then
  if ( allocated( g2_xyz ) ) then
    write(iout,*) '  g2_xyz       =   allocated to size ', &
                     size(g2_xyz,dim=1), size(g2_xyz,dim=2)
    write(iout,*) '  n, min, max = '
    do n=1,size(g2_xyz,dim=1)
      write(iout,*) n,minval(g2_xyz(n,:)),maxval(g2_xyz(n,:))
    enddo
  else
    write(iout,*) '  g2_xyz       = unallocated'
  endif
  endif



! Print info about faceusp.
  if ( index(local_query,'|faceusp|') .gt. 0 ) then
  if ( allocated( g2_faceusp ) ) then
    write(iout,*) '  g2_faceusp   =   allocated to size:', size(g2_faceusp)
    a_minval =  huge(1.0)
    a_maxval = -huge(1.0)
    do nf=1,size(g2_faceusp)     ! size should equal nface
      if ( allocated( g2_faceusp(nf)%v ) ) then
        do n=1,nq
          a_minval(n) = min(a_minval(n), minval(g2_faceusp(nf)%v(n,:,:)))
          a_maxval(n) = max(a_maxval(n), minval(g2_faceusp(nf)%v(n,:,:)))
        enddo
      endif
    enddo
    write(iout,*) '  n, min, max = '
    do n=1,nq
      write(iout,*) n,a_minval(n),a_maxval(n)
    enddo
  else
    write(iout,*) '  g2_faceusp   = unallocated'
  endif
  endif

! Print info about sol_nodes.
  if ( index(local_query,'|sol_nodes|') .gt. 0 ) then
  if ( allocated( g2_sol_nodes ) ) then
    write(iout,*) '  g2_sol_nodes =   allocated to size:', &
                  size(g2_sol_nodes,dim=1), size(g2_sol_nodes,dim=2)
    write(iout,*) '  n, min, max = '
    do n=1,nq
      write(iout,*) n,minval(g2_sol_nodes(n,:)),maxval(g2_sol_nodes(n,:))
    enddo
  else
    write(iout,*) '  g2_sol_nodes = unallocated'
  endif
  endif

! Print info about grid%sol.
  if ( index(local_query,'|grid%sol|') .gt. 0 ) then
  if ( allocated( g2_grid ) ) then
    write(iout,*) '  g2_grid      =   allocated'
    if ( allocated( g2_grid%sol ) ) then
      write(iout,*) '  g2_grid%sol  =   allocated to size:', &
                    size(g2_grid%sol,dim=1), size(g2_grid%sol,dim=2)
      write(iout,*) '  n, min, max = '
      do n=1,nq
        write(iout,*) n,minval(g2_grid%sol(n,:)),maxval(g2_grid%sol(n,:))
      enddo
    else
      write(iout,*) '  g2_grid%sol  = unallocated'
    endif
  else
    write(iout,*) '  g2_grid      = unallocated'
  endif
  endif

! Print info about usp.
  if ( index(local_query,'|usp|') .gt. 0 ) then
  if ( allocated( g2_usp ) ) then
    write(iout,*) '  g2_usp       =   allocated to size ', &
                  size(g2_usp,dim=1), size(g2_usp,dim=2)
    write(iout,*) '  n, min, max = '
    do n=1,nq
      write(iout,*) n,minval(g2_usp(n,:)),maxval(g2_usp(n,:))
    enddo
  else
    write(iout,*) '  g2_usp       = unallocated'
  endif
  endif


  write(iout,*) 'DEBUG_END:[mydebug_g2]'
  flush(iout)

  return
  end subroutine mydebug_g2
!##########################################################################
end module mydebug_mod
