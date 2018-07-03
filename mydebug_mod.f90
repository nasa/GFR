module mydebug_mod
  !
  ! Module to help with debugging.
  !
  use module_kind_types
  !
  implicit none
  !
  private
  !
  public :: mydebug
  public :: mydebug_g2
  !
contains
!
!###############################################################################
!
subroutine mydebug(caller_name,line_num,file_nam,query)
  !
  ! Debug the global solution variables.
  !
  use eqn_idx, only : nq
  use geovar,  only : nnode, ncell, nface, nfbnd, nr
  use geovar,  only : grid, xyz_nodes, xyz
  use flowvar, only : facexyz_dt, faceusp_dt, sol_nodes, usp
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: caller_name
  integer,          intent(in) :: line_num
  character(len=*), intent(in) :: file_nam
  !
  !.. Optional Arguments ..
  character(len=*), optional, intent(in) :: query
  !
  !.. Local Scalars ..
  integer :: n, nf
  !
  !.. Local Arrays ..
  real(wp) :: a_minval(1:nq)
  real(wp) :: a_maxval(1:nq)
  !
  !.. Local Allocatable Arrays ..
  character(len=:), allocatable :: local_query
  !
continue
  !
  if (mypnum == glb_root) then
    !
    ! Define a default list of variables to query.
    ! Override with optional list from arguments, if present.
    local_query = ''  ! to prevent gfortran warning
    local_query = '|nnode|ncell|nface|nfbnd|' // &
                  '|facexyz|xyz_nodes|grid%xyz|xyz|' // &
                  '|faceusp|sol_nodes|grid%sol|usp|'
    !
    if (present(query)) local_query = query
    !
    write (iout,*) 'DEBUG_BEG:[mydebug]'
    write (iout,6) trim(adjustl(caller_name)), line_num, trim(adjustl(file_nam))
    !
    ! Print info about nnode.
    !
    if (index(local_query,'|nnode|') > 0) then
      write (iout,5) "nnode",nnode
    end if
    !
    ! Print info about ncell.
    !
    if (index(local_query,'|ncell|') > 0) then
      write (iout,5) "ncell",ncell
    end if
    !
    ! Print info about ncell.
    !
    if (index(local_query,'|nface|') > 0) then
      write (iout,5) "nface",nface
    end if
    !
    ! Print info about ncell.
    !
    if (index(local_query,'|nfbnd|') > 0) then
      write (iout,5) "nfbnd",nfbnd
    end if
    !
    ! Print info about facexyz.
    !
    if (index(local_query,'|facexyz|') > 0) then
      if (allocated(facexyz_dt)) then
        write (iout,3) "facexyz_dt", size(facexyz_dt)
        a_minval =  huge(one)
        a_maxval = -huge(one)
        do nf = 1,size(facexyz_dt)     ! size should equal nface
          if (allocated(facexyz_dt(nf)%v)) then
            do n = 1,nr
              a_minval(n) = min(a_minval(n), minval(facexyz_dt(nf)%v(n,:)))
              a_maxval(n) = max(a_maxval(n), minval(facexyz_dt(nf)%v(n,:)))
            end do
          end if
        end do
        write (iout,4)
        do n = 1,nr
          write (iout,7) n,a_minval(n),a_maxval(n)
        end do
      else
        write (iout,2) "facexyz_dt"
      end if
    end if
    !
    ! Print info about xyz_nodes.
    !
    if (index(local_query,'|xyz_nodes|') > 0) then
      if (allocated(xyz_nodes)) then
        write (iout,3) "xyz_nodes", shape(xyz_nodes)
        write (iout,4)
        do n = 1,size(xyz_nodes,dim=1)
          write (iout,7) n,minval(xyz_nodes(n,:)),maxval(xyz_nodes(n,:))
        end do
      else
        write (iout,2) "xyz_nodes"
      end if
    end if
    !
    ! Print info about grid%xyz.
    !
    if (index(local_query,'|grid%xyz|') > 0) then
      if (allocated(grid)) then
        write (iout,1) "grid"
        if (allocated(grid%xyz)) then
          write (iout,3) "grid%xyz", shape(grid%xyz)
          write (iout,4)
          do n = 1,size(grid%xyz,dim=1)
            write (iout,7) n,minval(grid%xyz(n,:)),maxval(grid%xyz(n,:))
          end do
        else
          write (iout,2) "grid%xyz"
        end if
      else
        write (iout,2) "grid"
      end if
    end if
    !
    ! Print info about xyz.
    !
    if (index(local_query,'|xyz|') > 0) then
      if (allocated(xyz)) then
        write (iout,3) "xyz", shape(xyz)
        write (iout,4)
        do n = 1,size(xyz,dim=1)
          write (iout,7) n,minval(xyz(n,:)),maxval(xyz(n,:))
        end do
      else
        write (iout,2) "xyz"
      end if
    end if
    !
    ! Print info about faceusp.
    !
    if (index(local_query,'|faceusp|') > 0) then
      if (allocated(faceusp_dt)) then
        write (iout,3) "faceusp_dt", size(faceusp_dt)
        a_minval =  huge(one)
        a_maxval = -huge(one)
        do nf = 1,size(faceusp_dt)     ! size should equal nface
          if (allocated(faceusp_dt(nf)%v)) then
            do n = 1,nq
              a_minval(n) = min(a_minval(n), minval(faceusp_dt(nf)%v(n,:,:)))
              a_maxval(n) = max(a_maxval(n), minval(faceusp_dt(nf)%v(n,:,:)))
            end do
          end if
        end do
        write (iout,4)
        do n = 1,nq
          write (iout,7) n,a_minval(n),a_maxval(n)
        end do
      else
        write (iout,2) "faceusp_dt"
      end if
    end if
    !
    ! Print info about sol_nodes.
    !
    if (index(local_query,'|sol_nodes|') > 0) then
      if (allocated(sol_nodes)) then
        write (iout,3) "sol_nodes", shape(sol_nodes)
        write (iout,4)
        do n = 1,nq
          write (iout,7) n,minval(sol_nodes(n,:)),maxval(sol_nodes(n,:))
        end do
      else
        write (iout,2) "sol_nodes"
      end if
    end if
    !
    ! Print info about grid%sol.
    !
    if (index(local_query,'|grid%sol|') > 0) then
      if (allocated(grid)) then
        write (iout,1) "grid"
        if (allocated(grid%sol)) then
          write (iout,3) "grid%sol", shape(grid%sol)
          write (iout,4)
          do n = 1,nq
            write (iout,7) n,minval(grid%sol(n,:)),maxval(grid%sol(n,:))
          end do
        else
          write (iout,2) "grid%sol"
        end if
      else
        write (iout,2) "grid"
      end if
    end if
    !
    ! Print info about usp.
    !
    if (index(local_query,'|usp|') > 0) then
      if (allocated(usp)) then
        write (iout,3) "usp", shape(usp)
        write (iout,4)
        do n = 1,nq
          write (iout,7) n,minval(usp(n,:)),maxval(usp(n,:))
        end do
      else
        write (iout,2) "usp"
      end if
    end if
    !
    write (iout,*) 'DEBUG_END:[mydebug]'
    flush (iout)
    !
  end if
  !
  ! Print same info for g2_* variables.
  !
  call mydebug_g2(caller_name,line_num,file_nam,local_query)
  !
  ! Format Statements
  !
  1 format (2x,a,t14," =   allocated")
  2 format (2x,a,t14," = unallocated")
  3 format (2x,a,t14," =   allocated to size: ",10(i0,:,", "))
  4 format (2x,"n, min, max = ")
  5 format (2x,a,t12," = ",i0)
  6 format (2x,"Called from subroutine ",a,", line ",i0," in file ",a)
  7 format (2x,i1,2(2x,es14.6))
  !
end subroutine mydebug
!
!###############################################################################
!
subroutine mydebug_g2(caller_name,line_num,file_nam,query)
  !
  ! Debug the global solution variables for g2 grid.
  !
  use eqn_idx,          only : nq
  use geovar,           only : nr
  use interpsolvar_mod, only : g2_nnode, g2_ncell, g2_nface, g2_nfbnd
  use interpsolvar_mod, only : g2_grid, g2_xyz_nodes, g2_xyz
  use interpsolvar_mod, only : g2_facexyz_dt, g2_faceusp_dt
  use interpsolvar_mod, only : g2_sol_nodes, g2_usp
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: caller_name
  integer,          intent(in) :: line_num
  character(len=*), intent(in) :: file_nam
  !
  !.. Optional Arguments ..
  character(len=*), optional, intent(in) :: query
  !
  !.. Local Scalars ..
  integer :: n, nf
  !
  !.. Local Arrays ..
  real(wp) :: a_minval(1:nq)
  real(wp) :: a_maxval(1:nq)
  !
  !.. Local Allocatable Arrays ..
  character(len=:), allocatable :: local_query
  !
continue
  !
  if (mypnum == glb_root) then
    !
    ! Define a default list of variables to query.
    ! Override with optional list from arguments, if present.
    local_query = ''  ! to prevent gfortran warning
    local_query = '|nnode|ncell|nface|nfbnd|' // &
                  '|facexyz|xyz_nodes|grid%xyz|xyz|' // &
                  '|faceusp|sol_nodes|grid%sol|usp|'
    !
    if (present(query)) local_query = query
    !
    write (iout,*) 'DEBUG_BEG:[mydebug_g2]'
    write (iout,6) trim(adjustl(caller_name)), line_num, trim(adjustl(file_nam))
    !
    ! Print info about nnode.
    !
    if (index(local_query,'|nnode|') > 0) then
      write (iout,5) "g2_nnode",g2_nnode
    end if
    !
    ! Print info about ncell.
    !
    if (index(local_query,'|ncell|') > 0) then
      write (iout,5) "g2_ncell",g2_ncell
    end if
    !
    ! Print info about ncell.
    !
    if (index(local_query,'|nface|') > 0) then
      write (iout,5) "g2_nface",g2_nface
    end if
    !
    ! Print info about ncell.
    !
    if (index(local_query,'|nfbnd|') > 0) then
      write (iout,5) "g2_nfbnd",g2_nfbnd
    end if
    !
    ! Print info about facexyz.
    !
    if (index(local_query,'|facexyz|') > 0) then
      if (allocated(g2_facexyz_dt)) then
        write (iout,3) "g2_facexyz_dt", size(g2_facexyz_dt)
        a_minval =  huge(one)
        a_maxval = -huge(one)
        do nf = 1,size(g2_facexyz_dt)     ! size should equal nface
          if (allocated(g2_facexyz_dt(nf)%v)) then
            do n = 1,nr
              a_minval(n) = min(a_minval(n), minval(g2_facexyz_dt(nf)%v(n,:)))
              a_maxval(n) = max(a_maxval(n), minval(g2_facexyz_dt(nf)%v(n,:)))
            end do
          end if
        end do
        write (iout,4)
        do n = 1,nr
          write (iout,7) n,a_minval(n),a_maxval(n)
        end do
      else
        write (iout,2) "g2_facexyz_dt"
      end if
    end if
    !
    ! Print info about xyz_nodes.
    !
    if (index(local_query,'|xyz_nodes|') > 0) then
      if (allocated(g2_xyz_nodes)) then
        write (iout,3) "g2_xyz_nodes", shape(g2_xyz_nodes)
        write (iout,4)
        do n = 1,size(g2_xyz_nodes,dim=1)
          write (iout,7) n,minval(g2_xyz_nodes(n,:)),maxval(g2_xyz_nodes(n,:))
        end do
      else
        write (iout,2) "g2_xyz_nodes"
      end if
    end if
    !
    ! Print info about grid%xyz.
    !
    if (index(local_query,'|grid%xyz|') > 0) then
      if (allocated(g2_grid)) then
        write (iout,1) "g2_grid"
        if (allocated(g2_grid%xyz)) then
          write (iout,3) "g2_grid%xyz", shape(g2_grid%xyz)
          write (iout,4)
          do n = 1,size(g2_grid%xyz,dim=1)
            write (iout,7) n,minval(g2_grid%xyz(n,:)),maxval(g2_grid%xyz(n,:))
          end do
        else
          write (iout,2) "g2_grid%xyz"
        end if
      else
        write (iout,2) "g2_grid"
      end if
    end if
    !
    ! Print info about xyz.
    !
    if (index(local_query,'|xyz|') > 0) then
      if (allocated(g2_xyz)) then
        write (iout,3) "g2_xyz", shape(g2_xyz)
        write (iout,4)
        do n = 1,size(g2_xyz,dim=1)
          write (iout,7) n,minval(g2_xyz(n,:)),maxval(g2_xyz(n,:))
        end do
      else
        write (iout,2) "g2_xyz"
      end if
    end if
    !
    ! Print info about faceusp.
    !
    if (index(local_query,'|faceusp|') > 0) then
      if (allocated(g2_faceusp_dt)) then
        write (iout,3) "g2_faceusp_dt", size(g2_faceusp_dt)
        a_minval =  huge(one)
        a_maxval = -huge(one)
        do nf = 1,size(g2_faceusp_dt)     ! size should equal nface
          if (allocated(g2_faceusp_dt(nf)%v)) then
            do n = 1,nq
              a_minval(n) = min(a_minval(n), minval(g2_faceusp_dt(nf)%v(n,:,:)))
              a_maxval(n) = max(a_maxval(n), minval(g2_faceusp_dt(nf)%v(n,:,:)))
            end do
          end if
        end do
        write (iout,4)
        do n = 1,nq
          write (iout,7) n,a_minval(n),a_maxval(n)
        end do
      else
        write (iout,2) "g2_faceusp_dt"
      end if
    end if
    !
    ! Print info about sol_nodes.
    !
    if (index(local_query,'|sol_nodes|') > 0) then
      if (allocated(g2_sol_nodes)) then
        write (iout,3) "g2_sol_nodes", shape(g2_sol_nodes)
        write (iout,4)
        do n = 1,nq
          write (iout,7) n,minval(g2_sol_nodes(n,:)),maxval(g2_sol_nodes(n,:))
        end do
      else
        write (iout,2) "g2_sol_nodes"
      end if
    end if
    !
    ! Print info about grid%sol.
    !
    if (index(local_query,'|grid%sol|') > 0) then
      if (allocated(g2_grid)) then
        write (iout,1) "g2_grid"
        if (allocated(g2_grid%sol)) then
          write (iout,3) "g2_grid%sol", shape(g2_grid%sol)
          write (iout,4)
          do n = 1,nq
            write (iout,7) n,minval(g2_grid%sol(n,:)),maxval(g2_grid%sol(n,:))
          end do
        else
          write (iout,2) "g2_grid%sol"
        end if
      else
        write (iout,2) "g2_grid"
      end if
    end if
    !
    ! Print info about usp.
    !
    if (index(local_query,'|usp|') > 0) then
      if (allocated(g2_usp)) then
        write (iout,3) "g2_usp", shape(g2_usp)
        write (iout,4)
        do n = 1,nq
          write (iout,7) n,minval(g2_usp(n,:)),maxval(g2_usp(n,:))
        end do
      else
        write (iout,2) "g2_usp"
      end if
    end if
    !
    write (iout,*) 'DEBUG_END:[mydebug_g2]'
    flush (iout)
    !
  end if
  !
  ! Format Statements
  !
  1 format (2x,a,t17," =   allocated")
  2 format (2x,a,t17," = unallocated")
  3 format (2x,a,t17," =   allocated to size: ",10(i0,:,", "))
  4 format (2x,"n, min, max = ")
  5 format (2x,a,t15," = ",i0)
  6 format (2x,"Called from subroutine ",a,", line ",i0," in file ",a)
  7 format (2x,i1,2(2x,es14.6))
  !
end subroutine mydebug_g2
!
!###############################################################################
!
end module mydebug_mod
