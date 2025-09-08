!###############################################################################
!> Module for finding and limiting cells with bad thermo properties
!> (pressure, temperature, density) in either the volume, face, or node points.
!###############################################################################
module thermofix_mod
  !
  !.. Use Statements ..
  use module_kind_types, global_debug_timer =>      debug_timer, &
                    pure_global_debug_timer => pure_debug_timer
  !
  use eqn_idx_mod, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq
  use eqn_idx_mod, only : ncs,npres,ntemp,ng
  !
  use ovar, only : force_nd_subroutines
  !
  implicit none
  !
  private
  !
  !
  !.. Public Procedures ..
  !
  public  :: check_cell_thermo_properties  ! obsolete
  public  :: check_face_thermo_properties  ! obsolete
  public  :: check_thermo_properties
  !
  !
  !.. Private Procedures ..
  !
  private :: find_bad_thermo_cells
  private :: mend_bad_thermo_cells
  private :: reinterpolate_bad_cells
  private :: reinterpolate_bad_cells_3d
  private :: reinterpolate_bad_cells_2d
  private :: reinterpolate_bad_cells_1d
  private :: reinterpolate_bad_cells_nd
  private :: emergency_dynamic_filter_bad_cells
  private :: emergency_dynamic_filter_bad_cells_3d
  private :: emergency_dynamic_filter_bad_cells_2d
  private :: emergency_dynamic_filter_bad_cells_1d
  private :: emergency_dynamic_filter_bad_cells_nd
  private :: emergency_limit_bad_cells
  private :: emergency_limit_bad_cells_3d
  private :: emergency_limit_bad_cells_2d
  private :: emergency_limit_bad_cells_1d
  private :: emergency_limit_bad_cells_nd
  private :: emergency_limit2_bad_cells
  private :: emergency_scale_bad_cells
  private :: emergency_scale_bad_cells_3d
  private :: emergency_scale_bad_cells_2d
  private :: emergency_scale_bad_cells_1d
  private :: emergency_scale_bad_cells_nd
  !
contains
!
!###############################################################################
!> Check volume points for bad thermo properties
!> (pressure, temperature, density) and limit the cell.
!>
!> This routine is obsolete. Use check_face_properties instead.
!###############################################################################
subroutine check_cell_thermo_properties()
  !
  !.. Use Statements ..
  use geomvar_mod, only : ncell,cell
  use flowvar_mod, only : usp,gpsp
  use ovar,        only : thermofix,rhoref,tref,aref
  !
  use, intrinsic :: ieee_arithmetic, only : ieee_is_negative
  !
  !.. Local Scalars ..
  integer :: nc,n1,n2,n
  logical(lk) :: bad_cell_found
  real(wp) :: min_dens_limit,min_temp_limit,min_pres_limit
  !
  !.. Local Allocatable Arrays ..
  logical(lk), allocatable :: bad_cell(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_cell_thermo_properties"
  !
continue
  !
  if (ncs > 0) return
  !
  min_dens_limit = thermofix%min_density / rhoref
  min_temp_limit = thermofix%min_temperature / tref
  min_pres_limit = thermofix%min_pressure / (rhoref*aref*aref)
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Find any cells with solution points that have
  ! a negative density, temperature, or pressure
  !
  bad_cell_found = fals
  !
  bad_cell = init_array(ncell,fals) ! F2003 AUTO-REALLOC
  !$ACC DATA CREATE(bad_cell)
  !
  cell_loop: &
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    !
    do n = n1,n2
     !if ( ieee_is_negative(  usp(n,nec)   ) .or. &
     !     ieee_is_negative( gpsp(n,npres) ) .or. &
     !     ieee_is_negative( gpsp(n,ntemp) ) ) then
      if (  usp(n,nec)   < min_dens_limit .or. &
           gpsp(n,ntemp) < min_temp_limit .or. &
           gpsp(n,npres) < min_pres_limit ) then
        bad_cell(nc)   = true
        bad_cell_found = true
        cycle cell_loop
      end if
    end do
    !
  end do cell_loop
  !$ACC UPDATE DEVICE(bad_cell)
  !
  if (bad_cell_found) then
    !
    ! Warning message
    !
    write (iout,1) mypnum
    flush (iout)
    !
    ! Use the emergency procedures to limit any cells with
    ! a negative density, temperature, or pressure.
    ! Update solution, thermo and transport values at volume points
    ! and face points.
    !
    call mend_bad_thermo_cells(bad_cell,mend_faces=true,mend_iter=0)
  end if
  !
  !$ACC END DATA
  if (allocated(bad_cell)) deallocate ( bad_cell )
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/,2x,"WARNING!!! CPU # ",i6," FOUND A CELL WITH ", &
                 "A NEGATIVE DENSITY, TEMPERATURE, OR PRESSURE,",/, &
              2x,"WARNING!!!",14x,"AND IS ATTEMPTING TO USE THE ", &
                 "EMERGENCY LIMITER TO CONTINUE RUNNING!",/)
  !
end subroutine check_cell_thermo_properties
!
!###############################################################################
!> Check face points for bad thermo properties
!> (pressure, temperature, density), limit the cell, and update the face.
!>
!> This routine is obsolete. Use check_face_properties instead.
!###############################################################################
subroutine check_face_thermo_properties()
  !
  !.. Use Statements ..
  use geomvar_mod, only : ncell,nfbnd,nface,face
  use flowvar_mod, only : faceusp,facegpsp
  use ovar,        only : thermofix,rhoref,tref,aref
  !
  use, intrinsic :: ieee_arithmetic, only : ieee_is_negative
  !
  !.. Local Scalars ..
  integer :: nf,f1,f2,n,left_cell,right_cell
  logical(lk) :: bad_cell_found
  real(wp) :: min_dens_limit,min_temp_limit,min_pres_limit
  !
  !.. Local Allocatable Arrays ..
  logical(lk), allocatable :: bad_cell(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_face_thermo_properties"
  !
continue
  !
  if (ncs > 0) return
  !
  min_dens_limit = thermofix%min_density / rhoref
  min_temp_limit = thermofix%min_temperature / tref
  min_pres_limit = thermofix%min_pressure / (rhoref*aref*aref)
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Find any cells with solution points that have
  ! a negative density, temperature, or pressure
  !
  bad_cell_found = fals
  !
  bad_cell = init_array(ncell,fals) ! F2003 AUTO-REALLOC
  !$ACC DATA CREATE(bad_cell)
  !
  ! Search the boundary faces for any negative density, temperature, or pressure
  !
  bnd_face_loop: &
  do nf = 1,nfbnd
    !
    f1 = face(nf)%beg_fp
    f2 = face(nf)%end_fp
    !
    left_cell = face(nf)%left%cell
    !
    if (.not. bad_cell(left_cell)) then
      !
      do n = f1,f2
       !if ( ieee_is_negative(  faceusp(n,nec  ,1) ) .or. &
       !     ieee_is_negative( facegpsp(n,ntemp,1) ) .or. &
       !     ieee_is_negative( facegpsp(n,npres,1) ) ) then
        if (  faceusp(n,nec  ,1) < min_dens_limit .or. &
             facegpsp(n,ntemp,1) < min_temp_limit .or. &
             facegpsp(n,npres,1) < min_pres_limit ) then
          bad_cell(left_cell) = true
          bad_cell_found      = true
          cycle bnd_face_loop
        end if
      end do
      !
    end if
    !
  end do bnd_face_loop
  !
  ! Search the interior faces for any negative density, temperature, or pressure
  !
  int_face_loop: &
  do nf = nfbnd+1,nface
    !
    f1 = face(nf)%beg_fp
    f2 = face(nf)%end_fp
    !
    left_cell = face(nf)%left%cell
    right_cell = face(nf)%right%cell
    !
    ! Check the left side of this face
    !
    left_cell_check: &
    if (.not. bad_cell(left_cell)) then
      !
      do n = f1,f2
       !if ( ieee_is_negative(  faceusp(n,nec  ,1) ) .or. &
       !     ieee_is_negative( facegpsp(n,ntemp,1) ) .or. &
       !     ieee_is_negative( facegpsp(n,npres,1) ) ) then
        if (  faceusp(n,nec  ,1) < min_dens_limit .or. &
             facegpsp(n,ntemp,1) < min_temp_limit .or. &
             facegpsp(n,npres,1) < min_pres_limit ) then
          bad_cell(left_cell) = true
          bad_cell_found      = true
          exit left_cell_check
        end if
      end do
      !
    end if left_cell_check
    !
    ! Check the right side of this face
    !
    right_cell_check: &
    if (.not. bad_cell(right_cell)) then
      !
      do n = f1,f2
       !if ( ieee_is_negative(  faceusp(n,nec  ,2) ) .or. &
       !     ieee_is_negative( facegpsp(n,ntemp,2) ) .or. &
       !     ieee_is_negative( facegpsp(n,npres,2) ) ) then
        if (  faceusp(n,nec  ,2) < min_dens_limit .or. &
             facegpsp(n,ntemp,2) < min_temp_limit .or. &
             facegpsp(n,npres,2) < min_pres_limit ) then
          bad_cell(right_cell) = true
          bad_cell_found       = true
          exit right_cell_check
        end if
      end do
      !
    end if right_cell_check
    !
  end do int_face_loop
  !$ACC UPDATE DEVICE(bad_cell)
  !
  if (bad_cell_found) then
    !
    ! Warning message
    !
    write (iout,1) mypnum
    flush (iout)
    !
    ! Use the emergency procedures to limit any cells with
    ! a negative density, temperature, or pressure.
    ! Update solution, thermo and transport values at volume points
    ! and face points.
    !
    call mend_bad_thermo_cells(bad_cell,mend_faces=true,mend_iter=0)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/,2x,"WARNING!!! CPU # ",i6," FOUND A FACE WITH ", &
                 "A NEGATIVE DENSITY, TEMPERATURE, OR PRESSURE,",/, &
              2x,"WARNING!!!",14x,"AND IS ATTEMPTING TO USE THE ", &
                 "EMERGENCY LIMITER TO CONTINUE RUNNING!",/)
  !
end subroutine check_face_thermo_properties
!
!###############################################################################
!
!> Check cells, faces, nodes for bad thermo properties
!> (pressure, temperature, density), apply a correction, and update values.
!>
!> Updates the solution, thermo and transport values:
!> - at volume points in     usp,     gpsp.
!> - at face   points in faceusp, facegpsp, if requested.
!>
!> Assumes that the solution and thermo values for the volume and face are
!> current upon entry.
!
!###############################################################################
subroutine check_thermo_properties(check_cells,check_faces,check_nodes, &
                                   update_faces,verbose)
  !
  !.. Use Statements ..
  use geomvar_mod,  only : ncell
  use cell_ave_mod, only : allocate_cell_ave, deallocate_cell_ave, &
                           compute_cell_ave_sol
  use ovar,         only : thermofix
!-DEBUG_BEG:
  use geomvar_mod,  only : cell, xyz_nodes
  use flowvar_mod,  only : usp,gpsp
  use cell_ave_mod, only : cell_ave_sol
  use ovar,         only : itcur,current_rk_stage
  real(wp) :: xyz_ctr(3)
  integer  :: nc,n,n1,n2,nn,nr,num_cell_nodes
!-DEBUG_END:
  !
  !.. Arguments ..
  logical(lk), intent(in) :: check_cells  !< Check cell volumes.
  logical(lk), intent(in) :: check_faces  !< Check cell faces.
  logical(lk), intent(in) :: check_nodes  !< Check cell nodes.
  logical(lk), intent(in) :: update_faces !< Update face values.
  logical(lk), intent(in) :: verbose      !< Output warnings.
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_thermo_properties"
  !
  !.. Local Allocatable Arrays ..
  logical(lk), allocatable :: bad_cell(:)
  !
  !.. Local Scalars ..
  integer     :: iter_mend, max_iter_mend, num_bad_Cell
  logical(lk) :: bad_cell_found, bad_face_found, bad_node_found, any_bad_found
  logical(lk) :: need_cell_ave, mend_faces
  !
continue
  !
  ! Skip thermo check
  !
  ! Multi-species should work here (untested)
 !if ( ncs > 0 ) return  ! multi-species flow
  !
  if ( thermofix%method == TFM_OFF ) return
  if ( .not. (check_cells .or. check_faces .or. check_nodes) ) return
  !
  ! Cell ave does not change with repair, so compute once now if needed.
  !
  need_cell_ave = ( thermofix%method == TFM_SCALE )
!-DEBUG_BEG:
!-DEBUG_CMT: need cell ave to print below
  need_cell_ave = true
!-DEBUG_END:
  if ( need_cell_ave ) then
    call allocate_cell_ave()
    call  compute_cell_ave_sol()
  end if
  !
  ! Determine whether to mend face values too.
  !
  mend_faces = ( check_faces .or. update_faces )
  !
  ! Init bad_cell(:)=true so all cells are checked/repaired.
  !
  bad_cell = init_array(ncell,true) ! F2003 AUTO-REALLOC
  !
  loop_iter_mend : &
  do iter_mend = 1, thermofix%iterations
    !
    call find_bad_thermo_cells( check_cells, check_faces, &
                                check_nodes, bad_cell, &
                                bad_cell_found, bad_face_found, bad_node_found )
    !
    any_bad_found = ( bad_cell_found .or. bad_face_found .or. bad_node_found )
    if ( .not. any_bad_found ) then  ! all cells good
      exit loop_iter_mend
    end if
    !
    ! Write message
    !
    if ( verbose == true ) then
      num_bad_cell = count(bad_cell==true)
      if (iter_mend == 1) write (iout,10) itcur, current_rk_stage
      write (iout,1) iter_mend, mypnum, num_bad_cell
!-DEBUG_BEG:
!-DEBUG_CMT: Print more bad_cell info
      !
      ! If only 1 bad cell, write cell centroid xyz.
      !
!     if ( num_bad_cell == 1 .and. iter_mend == thermofix%iterations ) then
      if ( num_bad_cell == 1 ) then
        do nc = 1, ncell  ! cell loop
          if ( bad_cell(nc) ) then
            xyz_ctr(:) = zero
            nr = nme-nmb+1
            n1 = cell(nc)%beg_sp        ! beginning sol-pt in this cell
            n2 = cell(nc)%end_sp        ! ending    sol-pt in this cell
            num_cell_nodes = size(cell(nc)%node)
            do n = 1,num_cell_nodes
              nn = cell(nc)%node(n)%idx   ! index in nodeusp
              xyz_ctr(1:nr) = xyz_ctr(1:nr) + xyz_nodes(1:nr,nn)
            end do
            !
            xyz_ctr(:) = xyz_ctr(:) / real(num_cell_nodes,wp)
            write(iout,2) mypnum,nc, &
                          bad_cell_found,bad_face_found,bad_node_found
            write(iout,3) xyz_ctr(1:3)
            write(iout,4) minval(  usp(n1:n2,nec  ) ), &
                          minval( gpsp(n1:n2,npres) ), &
                          minval( gpsp(n1:n2,ntemp) )
            write(iout,5) cell_ave_sol(nec,nc)
2           format(10x,'part,cell,bad_cell,bad_face,bad_node =',i6,i8,3l2)
3           format(10x,'Cell x,y,z =',3es13.5)
4           format(10x,'Cell r,p,t =',3es13.5)
5           format(10x,'Ave  r     =',3es13.5)
            exit  ! exit cell loop
          end if
        end do
      end if
!-DEBUG_END:
      flush (iout)
    end if
    !
    ! Update solution, thermo, transport for volume points only in bad cells.
    ! Update face points too if requested.
    !
    call mend_bad_thermo_cells(bad_cell,mend_faces,iter_mend)
    !
  end do loop_iter_mend
  !
  ! Deallocate memory.
  !
  if (allocated(bad_cell)) deallocate ( bad_cell )
  if ( need_cell_ave ) then
    call deallocate_cell_ave()
  end if
  !
  ! Format Statements
  !
  1 format (6x,"MEND ITER # ",i0,": CPU # ",i6," IS REPAIRING ",i0," CELLS.")
  10 format (2x,"THERMOFIX WARNING!!! ITER # ",i0," (RK Stage ",i0,"):")
  !
end subroutine check_thermo_properties
!
!###############################################################################
!
!> Find cells with bad thermo properties in cells, faces, and/or nodes.
!>
!> - Checks density, pressure, and temperature.
!> - Assumes  the volume solution and thermo in     usp and     gpsp is current.
!> - Assumes  the face   solution and thermo in faceusp and facegpsp is current.
!> - Computes the node   solution and thermo as needed.
!
!###############################################################################
subroutine find_bad_thermo_cells(check_cells,check_faces,check_nodes,bad_cell, &
                                 bad_cell_found,bad_face_found,bad_node_found)
  !
  !.. Use Statements ..
  use geomvar_mod,     only : ncell,cell,face
  use flowvar_mod,     only : usp,gpsp,faceusp,facegpsp
  use ovar,            only : thermofix,rhoref,tref,aref
  use thermo_mod,      only : update_thermo_properties
  use order_mod,       only : maxSP
  use mlp_mod,         only : interp_cv_to_nodepts
  !
  !.. Arguments ..
  logical(lk), intent(in)    :: check_cells  !< Check cell volumes.
  logical(lk), intent(in)    :: check_faces  !< Check cell faces.
  logical(lk), intent(in)    :: check_nodes  !< Check cell nodes.
  !> On input,  marks the cells to be checked.
  !> On output, marks the bad cells that were found.
  logical(lk), intent(inout) :: bad_cell(:)
  logical(lk), intent(out)   :: bad_cell_found !< Found bad volume value.
  logical(lk), intent(out)   :: bad_face_found !< Found bad face   value.
  logical(lk), intent(out)   :: bad_node_found !< Found bad node   value.
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "find_bad_thermo_cells"
  integer, parameter :: maxNP = 8   ! Max node points per cell
  !
  !.. Local Arrays ..
  !> conserved variables at node points of one cell
  real(wp), dimension(1:maxNP,1:nq) :: thiscell_nodeusp
  !> gas property variables at node points of one cell
  real(wp), dimension(1:maxNP,1:ng) :: thiscell_nodegpsp
  !
  !.. Local Scalars ..
  integer :: this_geom,this_order
  integer :: num_cell_faces,num_cell_nodes
  integer :: nc,n0,n1,n2,np,n,m
  integer :: ncf,nf,sf,f1,f2
  real(wp) :: min_dens_limit,min_pres_limit,min_temp_limit
  !
continue
  !
  bad_cell_found = fals
  bad_face_found = fals
  bad_node_found = fals
  !
  if ( .not. (check_cells .or. check_faces .or. check_nodes) ) then
    ! Nothing to check. Mark cells good and exit.
    bad_cell(:) = fals
    return
  end if
  !
  min_dens_limit = thermofix%min_density / rhoref
  min_temp_limit = thermofix%min_temperature / tref
  min_pres_limit = thermofix%min_pressure / (rhoref*aref*aref)
  !
  loop_cells: &
  do nc = 1, ncell
    if ( bad_cell(nc) ) then
      !
      ! cell info
      !
      this_geom  = cell(nc)%geom  ! geometry of this cell
      this_order = cell(nc)%order ! order of this cell
      n1 = cell(nc)%beg_sp        ! beginning sol-pt in this cell
      n2 = cell(nc)%end_sp        ! ending    sol-pt in this cell
      n0 = cell(nc)%beg_sp - 1    ! offset index for sol-pts in this cell
      np = cell(nc)%end_sp - n0   ! number of sol-pts in this cell
      num_cell_faces = size(cell(nc)%face)
      num_cell_nodes = size(cell(nc)%node)
      !
      ! Check cell volume
      !
      if ( check_cells ) then
        loop_vol_points: &
        do n = n1,n2
          if ( gpsp(n,npres) <= min_pres_limit .or. &
               gpsp(n,ntemp) <= min_temp_limit .or. &
                usp(n,nec)   <= min_dens_limit ) then
            bad_cell_found = true
            cycle loop_cells
          end if
        end do loop_vol_points ! n
      end if ! check_cells
      !
      !
      ! Check cell faces
      !
      if ( check_faces ) then
        loop_faces: &
        do ncf = 1,num_cell_faces
          nf  = cell(nc)%face(ncf)%idx  ! index of the current face
          sf  = cell(nc)%face(ncf)%side ! side of face that cell nc is on
          f1  = face(nf)%beg_fp
          f2  = face(nf)%end_fp
          loop_face_points: &
          do n = f1,f2
            ! Test face points of this cell for negative values.
            if ( facegpsp(n,npres,sf) <= min_pres_limit .or. &
                 facegpsp(n,ntemp,sf) <= min_temp_limit .or. &
                  faceusp(n,nec  ,sf) <= min_dens_limit ) then
              bad_face_found = true
              cycle loop_cells
            end if
          end do loop_face_points ! n
        end do loop_faces ! ncf
      end if ! check_faces
      !
      !
      ! Check cell nodes
      !
      if ( check_nodes ) then
        !
        ! Update flow for node points of this cell.
        call interp_cv_to_nodepts( cell(nc), n0, usp, thiscell_nodeusp )
        !
        ! Update thermo for node points of this cell.
        call update_thermo_properties(1,num_cell_nodes, &
                                      thiscell_nodeusp,thiscell_nodegpsp, &
                                      loc='node',caller_name=pname)
        !
        ! Test node points of this cell for negative values.
        loop_node_points : &
        do n = 1,num_cell_nodes
          if ( thiscell_nodegpsp(n,npres) <= min_pres_limit .or. &
               thiscell_nodegpsp(n,ntemp) <= min_temp_limit .or. &
                thiscell_nodeusp(n,nec)   <= min_dens_limit ) then
            bad_node_found = true
            cycle loop_cells
          end if
        end do loop_node_points ! n
        !
      end if ! check_nodes
      !
      ! Passed all checks
      !
      bad_cell(nc) = fals
      !
    end if ! bad_cell
  end do loop_cells
  !
end subroutine find_bad_thermo_cells
!
!###############################################################################
!
!> Mend cells with bad thermo properties.
!>
!> Applies an emergency procedure to the volume solution.
!> Updates the solution, thermo, and transport properties  
!> - for volume points in     usp and     gpsp.  
!> - for face   points in faceusp and facegpsp, if requested.
!
!###############################################################################
subroutine mend_bad_thermo_cells(bad_cell,mend_faces,mend_iter)
  !
  !.. Use Statements ..
  use flowvar_mod, only : usp,gpsp,faceusp,facegpsp
  use ovar,        only : thermofix
  !
  !.. Arguments ..
  !> Marks the cells to be repaired.
  logical(lk), intent(inout) :: bad_cell(:)
  !
  !> Update the solution, thermo and transport properties at face points
  !> of bad cells.
  logical(lk), intent(in)    :: mend_faces
  !
  !> Current mend iteration
  integer, intent(in)        :: mend_iter
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "mend_bad_thermo_cells"
  !
continue
  !
  ! Apply correction to the solution at volume points of bad cells.
  !
  select case (thermofix%method)
  case (TFM_OFF)
    return
  case (TFM_LIMIT)
    call emergency_limit_bad_cells(bad_cell,usp)
  case (-TFM_LIMIT)
    call emergency_dynamic_filter_bad_cells(bad_cell,usp,mend_iter)
  case (TFM_SCALE)
    call emergency_scale_bad_cells(bad_cell,usp)
  case (TFM_LIMI2)
    call emergency_limit2_bad_cells(bad_cell,usp,gpsp)
  case default
    write (error_message,100) thermofix%method
    100 format( "Unexpected value, thermo_limiter_method = ", i0 )
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end select
  !
  ! Update thermo and transport properties at volume points of bad cells.
  !
  call bad_cell_update_properties_sp(bad_cell,usp,gpsp)
  !
  ! Update solution, thermo and transport properties at face points of
  ! bad cells.
  !
  if ( mend_faces ) then
    call reinterpolate_bad_cells(bad_cell,usp,faceusp,facegpsp)
  end if
  !
end subroutine mend_bad_thermo_cells
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
!###############################################################################
!> For cells marked as bad, reinterpolate the cell volume solution data to the
!> face points and update the face thermo and transport property values.
!###############################################################################
subroutine reinterpolate_bad_cells(bad_cell,usp,faceusp,facegpsp)
  !
  use geomvar_mod,   only : nfbnd,nface,face
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  !.. Optional Arguments ..
  logical(lk),    intent(in) :: bad_cell(:)
  real(wp),       intent(in) :: usp(:,:)
  real(wp),    intent(inout) :: faceusp(:,:,:)
  real(wp),    intent(inout) :: facegpsp(:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr,nf,sf,f1,f2
  !
  !.. Local Allocatable Arrays
  logical(lk), allocatable :: faces_to_update(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "reinterpolate_bad_cells"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = nme-nmb+1
  !
  faces_to_update = init_array(nface,2,fals) ! F2003 AUTO-REALLOC
  !$ACC DATA CREATE(faces_to_update)
  !
  if (force_nd_subroutines) then
    !
    call reinterpolate_bad_cells_nd(1,nq,bad_cell,usp, &
                                    faceusp,faces_to_update)
    !
  else
    !
    if (nr == 3) then
      call reinterpolate_bad_cells_3d(bad_cell,usp, &
                                      faceusp,faces_to_update)
    else if (nr == 2) then
      call reinterpolate_bad_cells_2d(bad_cell,usp, &
                                      faceusp,faces_to_update)
    else
      call reinterpolate_bad_cells_1d(bad_cell,usp, &
                                      faceusp,faces_to_update)
    end if
    !
    if (ncs > 0) then
      call reinterpolate_bad_cells_nd(1,ncs,bad_cell,usp, &
                                      faceusp,faces_to_update)
    end if
    !
  end if
  !$ACC UPDATE HOST(faces_to_update)
  !
  ! Update boundary faces, side-of-face=1 for interior solution.
  ! Do not update......... side-of-face=2 for boundary values.
  do nf = 1,nfbnd  ! loop boundary faces
    if (faces_to_update(nf,1)) then
      f1 = face(nf)%beg_fp
      f2 = face(nf)%end_fp
      call update_thermo_properties(f1,f2, &
                                    faceusp(:,:,1),facegpsp(:,:,1), &
                                    loc="face",caller_name=pname)
      call update_transport_properties(f1,f2, &
                                       faceusp(:,:,1),facegpsp(:,:,1))
    end if
  end do
  !
  ! Update interior faces, side-of-face=1:2 for left:right solutions.
  do sf = 1,2  ! loop side of face
    do nf = nfbnd+1,nface  ! loop interior faces
      if (faces_to_update(nf,sf)) then
        f1 = face(nf)%beg_fp
        f2 = face(nf)%end_fp
        call update_thermo_properties(f1,f2, &
                                      faceusp(:,:,sf),facegpsp(:,:,sf), &
                                      loc="face",caller_name=pname)
        call update_transport_properties(f1,f2, &
                                         faceusp(:,:,sf),facegpsp(:,:,sf))
      end if
    end do
  end do
  !
  !$ACC END DATA
  if (allocated(faces_to_update)) deallocate ( faces_to_update )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine reinterpolate_bad_cells
!
!###############################################################################
!> For cells marked as bad, reinterpolate the cell volume solution data to the
!> face points. Only updates solution, not thermo or transport.
!>
!> Optimized routine for 3d.
!###############################################################################
pure &
subroutine reinterpolate_bad_cells_3d(bad_cell,usp, &
                                      faceusp,faces_to_update)
  !
  use geomvar_mod,       only : ncell,cell,fp
  use interpolation_mod, only : interp
  use geomvar_mod,       only : face
  !
  !.. Optional Arguments ..
  logical(lk),    intent(in) :: bad_cell(:)
  real(wp),       intent(in) :: usp(:,:)
  real(wp),    intent(inout) :: faceusp(:,:,:)
  logical(lk), intent(inout) :: faces_to_update(:,:)
  !
  !.. Local Scalars ..
  integer :: nc,n0,ncf,f0
  integer :: nf,kf,lf,sf,i,j,nfp
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: fp_offset,face_rot
  integer :: num_cell_faces
  real(wp) :: a,face_ec,face_mx,face_my,face_mz,face_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "reinterpolate_bad_cells_3d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(faceusp) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(fp) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(usp) &
  !$ACC&                   PRESENT(interp) &
  !$ACC&                   PRESENT(bad_cell,faces_to_update)
  cell_loop: &
  do nc = 1,ncell
    !
    if (.not. bad_cell(nc)) cycle cell_loop
    !
    n0 = cell(nc)%beg_sp-1
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! Now interpolate the interior solution gradients to the flux points
    !
    ! Loop over the faces of this cell
    !
    num_cell_faces = size(cell(nc)%face)
    !
    do ncf = 1,num_cell_faces
      !
      nf = cell(nc)%face(ncf)%idx  ! index of the current face
      sf = cell(nc)%face(ncf)%side ! side of face that cell nc is on
      !
      ! Mark this side of this face as needing
      ! an update of the thermo properties
      !
      faces_to_update(nf,sf) = true
      !
      face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
      face_order = cell(nc)%face(ncf)%order ! order of current face
      !
      fp_offset = cell(nc)%face(ncf)%fp_offset
      face_rot = cell(nc)%face(ncf)%rotation
      !
      f0 = face(nf)%beg_fp - 1
      nfp = face(nf)%end_fp - f0
      !
      !$ACC LOOP VECTOR
      do kf = 1,nfp
        !
        ! kf: index of the flx-pt local to this face
        ! lf: index of flx-pt within interp%mat
        !
        lf = fp_offset + fp%geord(face_geom,face_order)%rot(face_rot)%pts(kf)
        !
        face_ec = zero
        face_mx = zero
        face_my = zero
        face_mz = zero
        face_ee = zero
        !
        do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
               interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
          !
          j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
          a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
          !
          face_ec = face_ec + a * usp(n0+j,nec)
          face_mx = face_mx + a * usp(n0+j,nmx)
          face_my = face_my + a * usp(n0+j,nmy)
          face_mz = face_mz + a * usp(n0+j,nmz)
          face_ee = face_ee + a * usp(n0+j,nee)
          !
        end do
        !
        faceusp(f0+kf,nec,sf) = face_ec
        faceusp(f0+kf,nmx,sf) = face_mx
        faceusp(f0+kf,nmy,sf) = face_my
        faceusp(f0+kf,nmz,sf) = face_mz
        faceusp(f0+kf,nee,sf) = face_ee
        !
      end do
      !
    end do
    !
  end do cell_loop
  !
end subroutine reinterpolate_bad_cells_3d
!
!###############################################################################
!> For cells marked as bad, reinterpolate the cell volume solution data to the
!> face points. Only updates solution, not thermo or transport.
!>
!> Optimized routine for 2d.
!###############################################################################
pure &
subroutine reinterpolate_bad_cells_2d(bad_cell,usp, &
                                      faceusp,faces_to_update)
  !
  use geomvar_mod,       only : ncell,cell,fp
  use interpolation_mod, only : interp
  use geomvar_mod,       only : face
  !
  !.. Optional Arguments ..
  logical(lk),    intent(in) :: bad_cell(:)
  real(wp),       intent(in) :: usp(:,:)
  real(wp),    intent(inout) :: faceusp(:,:,:)
  logical(lk), intent(inout) :: faces_to_update(:,:)
  !
  !.. Local Scalars ..
  integer :: nc,n0,ncf,f0
  integer :: nf,kf,lf,sf,i,j,nfp
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: fp_offset,face_rot
  integer :: num_cell_faces
  real(wp) :: a,face_ec,face_mx,face_my,face_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "reinterpolate_bad_cells_2d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(faceusp) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(fp) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(usp) &
  !$ACC&                   PRESENT(interp) &
  !$ACC&                   PRESENT(bad_cell,faces_to_update)
  cell_loop: &
  do nc = 1,ncell
    !
    if (.not. bad_cell(nc)) cycle cell_loop
    !
    n0 = cell(nc)%beg_sp-1
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! Now interpolate the interior solution gradients to the flux points
    !
    ! Loop over the faces of this cell
    !
    num_cell_faces = size(cell(nc)%face)
    !
    do ncf = 1,num_cell_faces
      !
      nf = cell(nc)%face(ncf)%idx  ! index of the current face
      sf = cell(nc)%face(ncf)%side ! side of face that cell nc is on
      !
      ! Mark this side of this face as needing
      ! an update of the thermo properties
      !
      faces_to_update(nf,sf) = true
      !
      face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
      face_order = cell(nc)%face(ncf)%order ! order of current face
      !
      fp_offset = cell(nc)%face(ncf)%fp_offset
      face_rot = cell(nc)%face(ncf)%rotation
      !
      f0 = face(nf)%beg_fp - 1
      nfp = face(nf)%end_fp - f0
      !
      !$ACC LOOP VECTOR
      do kf = 1,nfp
        !
        ! kf: index of the flx-pt local to this face
        ! lf: index of flx-pt within interp%mat
        !
        lf = fp_offset + fp%geord(face_geom,face_order)%rot(face_rot)%pts(kf)
        !
        face_ec = zero
        face_mx = zero
        face_my = zero
        face_ee = zero
        !
        do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
               interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
          !
          j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
          a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
          !
          face_ec = face_ec + a * usp(n0+j,nec)
          face_mx = face_mx + a * usp(n0+j,nmx)
          face_my = face_my + a * usp(n0+j,nmy)
          face_ee = face_ee + a * usp(n0+j,nee)
          !
        end do
        !
        faceusp(f0+kf,nec,sf) = face_ec
        faceusp(f0+kf,nmx,sf) = face_mx
        faceusp(f0+kf,nmy,sf) = face_my
        faceusp(f0+kf,nee,sf) = face_ee
        !
      end do
      !
    end do
    !
  end do cell_loop
  !
end subroutine reinterpolate_bad_cells_2d
!
!###############################################################################
!> For cells marked as bad, reinterpolate the cell volume solution data to the
!> face points. Only updates solution, not thermo or transport.
!>
!> Optimized routine for 1d.
!###############################################################################
pure &
subroutine reinterpolate_bad_cells_1d(bad_cell,usp, &
                                      faceusp,faces_to_update)
  !
  use geomvar_mod,       only : ncell,cell,fp
  use interpolation_mod, only : interp
  use geomvar_mod,       only : face
  !
  !.. Optional Arguments ..
  logical(lk),    intent(in) :: bad_cell(:)
  real(wp),       intent(in) :: usp(:,:)
  real(wp),    intent(inout) :: faceusp(:,:,:)
  logical(lk), intent(inout) :: faces_to_update(:,:)
  !
  !.. Local Scalars ..
  integer :: nc,n0,ncf,f0
  integer :: nf,kf,lf,sf,i,j,nfp
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: fp_offset,face_rot
  integer :: num_cell_faces
  real(wp) :: a,face_ec,face_mx,face_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "reinterpolate_bad_cells_1d"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(faceusp) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(fp) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(usp) &
  !$ACC&                   PRESENT(interp) &
  !$ACC&                   PRESENT(bad_cell,faces_to_update)
  cell_loop: &
  do nc = 1,ncell
    !
    if (.not. bad_cell(nc)) cycle cell_loop
    !
    n0 = cell(nc)%beg_sp-1
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! Now interpolate the interior solution gradients to the flux points
    !
    ! Loop over the faces of this cell
    !
    num_cell_faces = size(cell(nc)%face)
    !
    do ncf = 1,num_cell_faces
      !
      nf = cell(nc)%face(ncf)%idx  ! index of the current face
      sf = cell(nc)%face(ncf)%side ! side of face that cell nc is on
      !
      ! Mark this side of this face as needing
      ! an update of the thermo properties
      !
      faces_to_update(nf,sf) = true
      !
      face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
      face_order = cell(nc)%face(ncf)%order ! order of current face
      !
      fp_offset = cell(nc)%face(ncf)%fp_offset
      face_rot = cell(nc)%face(ncf)%rotation
      !
      f0 = face(nf)%beg_fp - 1
      nfp = face(nf)%end_fp - f0
      !
      !$ACC LOOP VECTOR
      do kf = 1,nfp
        !
        ! kf: index of the flx-pt local to this face
        ! lf: index of flx-pt within interp%mat
        !
        lf = fp_offset + fp%geord(face_geom,face_order)%rot(face_rot)%pts(kf)
        !
        face_ec = zero
        face_mx = zero
        face_ee = zero
        !
        do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
               interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
          !
          j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
          a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
          !
          face_ec = face_ec + a * usp(n0+j,nec)
          face_mx = face_mx + a * usp(n0+j,nmx)
          face_ee = face_ee + a * usp(n0+j,nee)
          !
        end do
        !
        faceusp(f0+kf,nec,sf) = face_ec
        faceusp(f0+kf,nmx,sf) = face_mx
        faceusp(f0+kf,nee,sf) = face_ee
        !
      end do
      !
    end do
    !
  end do cell_loop
  !
end subroutine reinterpolate_bad_cells_1d
!
!###############################################################################
!> For cells marked as bad, reinterpolate the cell volume solution data to the
!> face points. Only updates solution, not thermo or transport.
!>
!> General routine for nd.
!###############################################################################
pure &
subroutine reinterpolate_bad_cells_nd(mbeg,mend,bad_cell,usp, &
                                      faceusp,faces_to_update)
  !
  use geomvar_mod,       only : ncell,cell,fp
  use interpolation_mod, only : interp
  use geomvar_mod,       only : face
  !
  !.. Optional Arguments ..
  integer,        intent(in) :: mbeg,mend
  logical(lk),    intent(in) :: bad_cell(:)
  real(wp),       intent(in) :: usp(:,:)
  real(wp),    intent(inout) :: faceusp(:,:,:)
  logical(lk), intent(inout) :: faces_to_update(:,:)
  !
  !.. Local Scalars ..
  integer :: m,nc,n0,ncf,f0
  integer :: nf,kf,lf,sf,i,j,nfp
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: fp_offset,face_rot
  integer :: num_cell_faces
  real(wp) :: a,val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "reinterpolate_bad_cells_nd"
  !
continue
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(faceusp) &
  !$ACC&                   PRESENT(face) &
  !$ACC&                   PRESENT(fp) &
  !$ACC&                   PRESENT(cell) &
  !$ACC&                   PRESENT(usp) &
  !$ACC&                   PRESENT(interp) &
  !$ACC&                   PRESENT(bad_cell,faces_to_update)
  cell_loop: &
  do nc = 1,ncell
    !
    if (.not. bad_cell(nc)) cycle cell_loop
    !
    n0 = cell(nc)%beg_sp-1
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! Now interpolate the interior solution gradients to the flux points
    !
    ! Loop over the faces of this cell
    !
    num_cell_faces = size(cell(nc)%face)
    !
    do ncf = 1,num_cell_faces
      !
      nf = cell(nc)%face(ncf)%idx  ! index of the current face
      sf = cell(nc)%face(ncf)%side ! side of face that cell nc is on
      !
      ! Mark this side of this face as needing
      ! an update of the thermo properties
      !
      faces_to_update(nf,sf) = true
      !
      face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
      face_order = cell(nc)%face(ncf)%order ! order of current face
      !
      fp_offset = cell(nc)%face(ncf)%fp_offset
      face_rot = cell(nc)%face(ncf)%rotation
      !
      f0 = face(nf)%beg_fp - 1
      nfp = face(nf)%end_fp - f0
      !
      !$ACC LOOP VECTOR COLLAPSE(2)
      do kf = 1,nfp
        do m = mbeg,mend
          !
          ! kf: index of the flx-pt local to this face
          ! lf: index of flx-pt within interp%mat
          !
          lf = fp_offset + fp%geord(face_geom,face_order)%rot(face_rot)%pts(kf)
          !
          val = zero
          !
          do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1,&
                 interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
            !
            j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
            a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
            !
            val = val + a * usp(n0+j,m)
            !
          end do
          !
          faceusp(f0+kf,m,sf) = val
          !
        end do
      end do
      !
    end do
    !
  end do cell_loop
  !
end subroutine reinterpolate_bad_cells_nd
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
!###############################################################################
!> Apply the dynamically scaled sponge filter to bad cells and update
!> the volume solution in usp.
!###############################################################################
subroutine emergency_dynamic_filter_bad_cells(bad_cell,usp,mend_iter)
  !
  !.. Use Statements ..
  use ovar, only : thermofix
  !
  !.. Formal Arguments ..
  logical(lk), intent(in) :: bad_cell(:)
  real(wp), intent(inout) :: usp(:,:)
  integer,     intent(in) :: mend_iter
  !
  !.. Local Scalars ..
  integer  :: nr
  real(wp) :: df_coef
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "emergency_dynamic_filter_bad_cells"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = nme-nmb+1
  !
  ! Get the dynamic filter coefficient, which is a value between 0 and 1.
  ! A value of 0 will do no filtering, and a value of 1 will apply the
  ! full filter.
  !
  if (any([mend_iter,thermofix%iterations] == 0)) then
    !
    ! Some logic to set df_coef to 1 if mend_iter is 0.
    ! This is just in case someone ends up using the
    ! check_cell_thermo_properties or check_face_thermo_properties
    ! subroutines which are currently not used anywhere.
    !
    df_coef = one
    !
  else
    !
    ! We will use mend_iter/thermofix%iterations as the dynamic filter
    ! coefficient. For example, if thermofix%iterations is 10, the dynamic
    ! filter coefficient will be [0.1, 0.2, 0.3, ..., 1.0] as mend_iter
    ! progresses from [1, 2, 3, ..., 10].
    !
    df_coef = real(mend_iter,kind=wp)/real(thermofix%iterations,kind=wp)
    !
  end if
  !
  if (force_nd_subroutines) then
    !
    call emergency_dynamic_filter_bad_cells_nd(1,nq,bad_cell,usp,df_coef)
    !
  else
    !
    if (nr == 3) then
      call emergency_dynamic_filter_bad_cells_3d(bad_cell,usp,df_coef)
    else if (nr == 2) then
      call emergency_dynamic_filter_bad_cells_2d(bad_cell,usp,df_coef)
    else
      call emergency_dynamic_filter_bad_cells_1d(bad_cell,usp,df_coef)
    end if
    !
    if (ncs > 0) then
      call emergency_dynamic_filter_bad_cells_nd(1,ncs,bad_cell,usp,df_coef)
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine emergency_dynamic_filter_bad_cells
!
!###############################################################################
!> Apply the dynamically scaled sponge filter to bad cells and update
!> the volume solution in usp.
!>
!> Optimized routine for 3d.
!###############################################################################
pure &
subroutine emergency_dynamic_filter_bad_cells_3d(bad_cell,usp,df_coef)
  !
  !.. Use Statements ..
  use order_mod,   only : maxSP
  use geomvar_mod, only : ncell,cell
  use filter_mod,  only : filter
  !
  !.. Formal Arguments ..
  logical(lk), intent(in) :: bad_cell(:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp),    intent(in) :: df_coef
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,filt_ec,filt_mx,filt_my,filt_mz,filt_ee,s,oms
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "emergency_dynamic_filter_bad_cells_3d"
  !
continue
  !
  ! Save the dynamic filter coefficient to local variable s and precompute 1-s.
  ! NOTE: This is mostly done to keep this code as similar as possible to the
  !       dynamic filtering code in the limiters module.
  !
  s = df_coef
  oms = one-s
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,bad_cell,filter) &
  !$ACC&              PRIVATE(old_usp)
  filt_cell_loop: &
  do nc = 1,ncell
    !
    ! If the current cell is not troubled, skip to the next one
    !
    if (.not.bad_cell(nc)) then
      cycle filt_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      old_usp(n,nec) = usp(n0+n,nec)
      old_usp(n,nmx) = usp(n0+n,nmx)
      old_usp(n,nmy) = usp(n0+n,nmy)
      old_usp(n,nmz) = usp(n0+n,nmz)
      old_usp(n,nee) = usp(n0+n,nee)
    end do
    !
    ! Apply the filter emergency limiter by multplying the filter
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      filt_ec = oms*old_usp(n,nec)
      filt_mx = oms*old_usp(n,nmx)
      filt_my = oms*old_usp(n,nmy)
      filt_mz = oms*old_usp(n,nmz)
      filt_ee = oms*old_usp(n,nee)
      !
      do i = filter(this_geom,this_order)%sponge%col_ptr(n)+1, &
             filter(this_geom,this_order)%sponge%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%sponge%col_idx(i)
        a = filter(this_geom,this_order)%sponge%col_val(i)
        !
        ! Use the scaled marker to linearly scale the strength of the 
        ! filter. The scaling is applied to all elements of the base
        ! filter matrix, and a (1-scaling) term needs to be added to
        ! the diagonal elements of the matrix.
        ! NOTE: Fmat = Vmat.f.Vinv   
        !          where : Fmat = base filter matrix
        !                  f = diagonal matrix of the filter coefficients
        !       Scaled(Fmat) = SFmat = Vmat.[(1-s)*I + s*f).Vinv
        !          where : s = scaled marker (0=no filter, 1=full filter)
        !                  I = identity matrix
        !       SFmat = (1-s)*Vmat.I.Vinv + s*Vmat.f.Vinv
        !             = (1-s)*I + s*Fmat
        !       We can use the (1-s)*I part as the initialization of the
        !       filter sums because that value is independent of the actual
        !       filter matrix.
        !
        filt_ec = filt_ec + a * s * old_usp(k,nec)
        filt_mx = filt_mx + a * s * old_usp(k,nmx)
        filt_my = filt_my + a * s * old_usp(k,nmy)
        filt_mz = filt_mz + a * s * old_usp(k,nmz)
        filt_ee = filt_ee + a * s * old_usp(k,nee)
        !
      end do
      !
      usp(n0+n,nec) = filt_ec
      usp(n0+n,nmx) = filt_mx
      usp(n0+n,nmy) = filt_my
      usp(n0+n,nmz) = filt_mz
      usp(n0+n,nee) = filt_ee
      !
    end do
    !
  end do filt_cell_loop
  !
end subroutine emergency_dynamic_filter_bad_cells_3d
!
!###############################################################################
!> Apply the dynamically scaled sponge filter to bad cells and update
!> the volume solution in usp.
!>
!> Optimized routine for 2d.
!###############################################################################
pure &
subroutine emergency_dynamic_filter_bad_cells_2d(bad_cell,usp,df_coef)
  !
  !.. Use Statements ..
  use order_mod,   only : maxSP
  use geomvar_mod, only : ncell,cell
  use filter_mod,  only : filter
  !
  !.. Formal Arguments ..
  logical(lk), intent(in) :: bad_cell(:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp),    intent(in) :: df_coef
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,filt_ec,filt_mx,filt_my,filt_ee,s,oms
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "emergency_dynamic_filter_bad_cells_2d"
  !
continue
  !
  ! Save the dynamic filter coefficient to local variable s and precompute 1-s.
  ! NOTE: This is mostly done to keep this code as similar as possible to the
  !       dynamic filtering code in the limiters module.
  !
  s = df_coef
  oms = one-s
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,bad_cell,filter) &
  !$ACC&              PRIVATE(old_usp)
  filt_cell_loop: &
  do nc = 1,ncell
    !
    ! If the current cell is not troubled, skip to the next one
    !
    if (.not.bad_cell(nc)) then
      cycle filt_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      old_usp(n,nec) = usp(n0+n,nec)
      old_usp(n,nmx) = usp(n0+n,nmx)
      old_usp(n,nmy) = usp(n0+n,nmy)
      old_usp(n,nee) = usp(n0+n,nee)
    end do
    !
    ! Apply the filter emergency limiter by multplying the filter
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      filt_ec = oms*old_usp(n,nec)
      filt_mx = oms*old_usp(n,nmx)
      filt_my = oms*old_usp(n,nmy)
      filt_ee = oms*old_usp(n,nee)
      !
      do i = filter(this_geom,this_order)%sponge%col_ptr(n)+1, &
             filter(this_geom,this_order)%sponge%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%sponge%col_idx(i)
        a = filter(this_geom,this_order)%sponge%col_val(i)
        !
        ! Use the scaled marker to linearly scale the strength of the 
        ! filter. The scaling is applied to all elements of the base
        ! filter matrix, and a (1-scaling) term needs to be added to
        ! the diagonal elements of the matrix.
        ! NOTE: Fmat = Vmat.f.Vinv   
        !          where : Fmat = base filter matrix
        !                  f = diagonal matrix of the filter coefficients
        !       Scaled(Fmat) = SFmat = Vmat.[(1-s)*I + s*f).Vinv
        !          where : s = scaled marker (0=no filter, 1=full filter)
        !                  I = identity matrix
        !       SFmat = (1-s)*Vmat.I.Vinv + s*Vmat.f.Vinv
        !             = (1-s)*I + s*Fmat
        !       We can use the (1-s)*I part as the initialization of the
        !       filter sums because that value is independent of the actual
        !       filter matrix.
        !
        filt_ec = filt_ec + a * s * old_usp(k,nec)
        filt_mx = filt_mx + a * s * old_usp(k,nmx)
        filt_my = filt_my + a * s * old_usp(k,nmy)
        filt_ee = filt_ee + a * s * old_usp(k,nee)
        !
      end do
      !
      usp(n0+n,nec) = filt_ec
      usp(n0+n,nmx) = filt_mx
      usp(n0+n,nmy) = filt_my
      usp(n0+n,nee) = filt_ee
      !
    end do
    !
  end do filt_cell_loop
  !
end subroutine emergency_dynamic_filter_bad_cells_2d
!
!###############################################################################
!> Apply the dynamically scaled sponge filter to bad cells and update
!> the volume solution in usp.
!>
!> Optimized routine for 1d.
!###############################################################################
pure &
subroutine emergency_dynamic_filter_bad_cells_1d(bad_cell,usp,df_coef)
  !
  !.. Use Statements ..
  use order_mod,   only : maxSP
  use geomvar_mod, only : ncell,cell
  use filter_mod,  only : filter
  !
  !.. Formal Arguments ..
  logical(lk), intent(in) :: bad_cell(:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp),    intent(in) :: df_coef
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,filt_ec,filt_mx,filt_ee,s,oms
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "emergency_dynamic_filter_bad_cells_1d"
  !
continue
  !
  ! Save the dynamic filter coefficient to local variable s and precompute 1-s.
  ! NOTE: This is mostly done to keep this code as similar as possible to the
  !       dynamic filtering code in the limiters module.
  !
  s = df_coef
  oms = one-s
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,bad_cell,filter) &
  !$ACC&              PRIVATE(old_usp)
  filt_cell_loop: &
  do nc = 1,ncell
    !
    ! If the current cell is not troubled, skip to the next one
    !
    if (.not.bad_cell(nc)) then
      cycle filt_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      old_usp(n,nec) = usp(n0+n,nec)
      old_usp(n,nmx) = usp(n0+n,nmx)
      old_usp(n,nee) = usp(n0+n,nee)
    end do
    !
    ! Apply the filter emergency limiter by multplying the filter
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      filt_ec = oms*old_usp(n,nec)
      filt_mx = oms*old_usp(n,nmx)
      filt_ee = oms*old_usp(n,nee)
      !
      do i = filter(this_geom,this_order)%sponge%col_ptr(n)+1, &
             filter(this_geom,this_order)%sponge%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%sponge%col_idx(i)
        a = filter(this_geom,this_order)%sponge%col_val(i)
        !
        ! Use the scaled marker to linearly scale the strength of the 
        ! filter. The scaling is applied to all elements of the base
        ! filter matrix, and a (1-scaling) term needs to be added to
        ! the diagonal elements of the matrix.
        ! NOTE: Fmat = Vmat.f.Vinv   
        !          where : Fmat = base filter matrix
        !                  f = diagonal matrix of the filter coefficients
        !       Scaled(Fmat) = SFmat = Vmat.[(1-s)*I + s*f).Vinv
        !          where : s = scaled marker (0=no filter, 1=full filter)
        !                  I = identity matrix
        !       SFmat = (1-s)*Vmat.I.Vinv + s*Vmat.f.Vinv
        !             = (1-s)*I + s*Fmat
        !       We can use the (1-s)*I part as the initialization of the
        !       filter sums because that value is independent of the actual
        !       filter matrix.
        !
        filt_ec = filt_ec + a * s * old_usp(k,nec)
        filt_mx = filt_mx + a * s * old_usp(k,nmx)
        filt_ee = filt_ee + a * s * old_usp(k,nee)
        !
      end do
      !
      usp(n0+n,nec) = filt_ec
      usp(n0+n,nmx) = filt_mx
      usp(n0+n,nee) = filt_ee
      !
    end do
    !
  end do filt_cell_loop
  !
end subroutine emergency_dynamic_filter_bad_cells_1d
!
!###############################################################################
!> Apply the dynamically scaled sponge filter to bad cells and update
!> the volume solution in usp.
!>
!> General routine for nd.
!###############################################################################
pure &
subroutine emergency_dynamic_filter_bad_cells_nd(mbeg,mend,bad_cell,usp,df_coef)
  !
  !.. Use Statements ..
  use order_mod,   only : maxSP
  use geomvar_mod, only : ncell,cell
  use filter_mod,  only : filter
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: mbeg,mend
  logical(lk), intent(in) :: bad_cell(:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp),    intent(in) :: df_coef
  !
  !.. Local Scalars ..
  integer  :: i,k,m,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,val,s,oms
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "emergency_dynamic_filter_bad_cells_nd"
  !
continue
  !
  ! Save the dynamic filter coefficient to local variable s and precompute 1-s.
  ! NOTE: This is mostly done to keep this code as similar as possible to the
  !       dynamic filtering code in the limiters module.
  !
  s = df_coef
  oms = one-s
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,bad_cell,filter) &
  !$ACC&              PRIVATE(old_usp)
  filt_cell_loop: &
  do nc = 1,ncell
    !
    ! If the current cell is not troubled, skip to the next one
    !
    if (.not.bad_cell(nc)) then
      cycle filt_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do m = mbeg,mend
      do n = 1,np
        old_usp(n,m) = usp(n0+n,m)
      end do
    end do
    !
    ! Apply the filter emergency limiter by multplying the filter
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do m = mbeg,mend
      do n = 1,np
        !
        val = oms*old_usp(n,m)
        !
        do i = filter(this_geom,this_order)%sponge%col_ptr(n)+1, &
               filter(this_geom,this_order)%sponge%col_ptr(n+1)
          !
          k = filter(this_geom,this_order)%sponge%col_idx(i)
          a = filter(this_geom,this_order)%sponge%col_val(i)
          !
          ! Use the scaled marker to linearly scale the strength of the 
          ! filter. The scaling is applied to all elements of the base
          ! filter matrix, and a (1-scaling) term needs to be added to
          ! the diagonal elements of the matrix.
          ! NOTE: Fmat = Vmat.f.Vinv   
          !          where : Fmat = base filter matrix
          !                  f = diagonal matrix of the filter coefficients
          !       Scaled(Fmat) = SFmat = Vmat.[(1-s)*I + s*f).Vinv
          !          where : s = scaled marker (0=no filter, 1=full filter)
          !                  I = identity matrix
          !       SFmat = (1-s)*Vmat.I.Vinv + s*Vmat.f.Vinv
          !             = (1-s)*I + s*Fmat
          !       We can use the (1-s)*I part as the initialization of the
          !       filter sums because that value is independent of the actual
          !       filter matrix.
          !
          val = val + a * s * old_usp(k,m)
          !
        end do
        !
        usp(n0+n,m) = val
        !
      end do
    end do
    !
  end do filt_cell_loop
  !
end subroutine emergency_dynamic_filter_bad_cells_nd
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
!###############################################################################
!> Apply the sponge filter to bad cells and update the volume solution in usp.
!###############################################################################
subroutine emergency_limit_bad_cells(bad_cell,usp)
  !
  !.. Use Statements ..
  use geomvar_mod, only : ncell,cell
  !
  !.. Formal Arguments ..
  logical(lk), intent(in) :: bad_cell(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: nc,nr,n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "emergency_limit_bad_cells"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = nme-nmb+1
  !
  if (force_nd_subroutines) then
    !
    call emergency_limit_bad_cells_nd(1,nq,bad_cell,usp)
    !
  else
    !
    if (nr == 3) then
      call emergency_limit_bad_cells_3d(bad_cell,usp)
    else if (nr == 2) then
      call emergency_limit_bad_cells_2d(bad_cell,usp)
    else
      call emergency_limit_bad_cells_1d(bad_cell,usp)
    end if
    !
    if (ncs > 0) then
      call emergency_limit_bad_cells_nd(1,ncs,bad_cell,usp)
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine emergency_limit_bad_cells
!
!###############################################################################
!> Apply the sponge filter to bad cells and update the volume solution in usp.
!>
!> Optimized routine for 3d.
!###############################################################################
pure &
subroutine emergency_limit_bad_cells_3d(bad_cell,usp)
  !
  !.. Use Statements ..
  use order_mod,   only : maxSP
  use geomvar_mod, only : ncell,cell
  use filter_mod,  only : filter
  !
  !.. Formal Arguments ..
  logical(lk), intent(in) :: bad_cell(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,filt_ec,filt_mx,filt_my,filt_mz,filt_ee
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "emergency_limit_bad_cells_3d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,bad_cell,filter) &
  !$ACC&              PRIVATE(old_usp)
  filt_cell_loop: &
  do nc = 1,ncell
    !
    ! If the current cell is not troubled, skip to the next one
    !
    if (.not.bad_cell(nc)) then
      cycle filt_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      old_usp(n,nec) = usp(n0+n,nec)
      old_usp(n,nmx) = usp(n0+n,nmx)
      old_usp(n,nmy) = usp(n0+n,nmy)
      old_usp(n,nmz) = usp(n0+n,nmz)
      old_usp(n,nee) = usp(n0+n,nee)
    end do
    !
    ! Apply the filter emergency limiter by multplying the filter
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      filt_ec = zero
      filt_mx = zero
      filt_my = zero
      filt_mz = zero
      filt_ee = zero
      !
      do i = filter(this_geom,this_order)%sponge%col_ptr(n)+1, &
             filter(this_geom,this_order)%sponge%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%sponge%col_idx(i)
        a = filter(this_geom,this_order)%sponge%col_val(i)
        !
        filt_ec = filt_ec + a * old_usp(k,nec)
        filt_mx = filt_mx + a * old_usp(k,nmx)
        filt_my = filt_my + a * old_usp(k,nmy)
        filt_mz = filt_mz + a * old_usp(k,nmz)
        filt_ee = filt_ee + a * old_usp(k,nee)
        !
      end do
      !
      usp(n0+n,nec) = filt_ec
      usp(n0+n,nmx) = filt_mx
      usp(n0+n,nmy) = filt_my
      usp(n0+n,nmz) = filt_mz
      usp(n0+n,nee) = filt_ee
      !
    end do
    !
  end do filt_cell_loop
  !
end subroutine emergency_limit_bad_cells_3d
!
!###############################################################################
!> Apply the sponge filter to bad cells and update the volume solution in usp.
!>
!> Optimized routine for 2d.
!###############################################################################
pure &
subroutine emergency_limit_bad_cells_2d(bad_cell,usp)
  !
  !.. Use Statements ..
  use order_mod,   only : maxSP
  use geomvar_mod, only : ncell,cell
  use filter_mod,  only : filter
  !
  !.. Formal Arguments ..
  logical(lk), intent(in) :: bad_cell(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,filt_ec,filt_mx,filt_my,filt_ee
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "emergency_limit_bad_cells_2d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,bad_cell,filter) &
  !$ACC&              PRIVATE(old_usp)
  filt_cell_loop: &
  do nc = 1,ncell
    !
    ! If the current cell is not troubled, skip to the next one
    !
    if (.not.bad_cell(nc)) then
      cycle filt_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      old_usp(n,nec) = usp(n0+n,nec)
      old_usp(n,nmx) = usp(n0+n,nmx)
      old_usp(n,nmy) = usp(n0+n,nmy)
      old_usp(n,nee) = usp(n0+n,nee)
    end do
    !
    ! Apply the filter emergency limiter by multplying the filter
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      filt_ec = zero
      filt_mx = zero
      filt_my = zero
      filt_ee = zero
      !
      do i = filter(this_geom,this_order)%sponge%col_ptr(n)+1, &
             filter(this_geom,this_order)%sponge%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%sponge%col_idx(i)
        a = filter(this_geom,this_order)%sponge%col_val(i)
        !
        filt_ec = filt_ec + a * old_usp(k,nec)
        filt_mx = filt_mx + a * old_usp(k,nmx)
        filt_my = filt_my + a * old_usp(k,nmy)
        filt_ee = filt_ee + a * old_usp(k,nee)
        !
      end do
      !
      usp(n0+n,nec) = filt_ec
      usp(n0+n,nmx) = filt_mx
      usp(n0+n,nmy) = filt_my
      usp(n0+n,nee) = filt_ee
      !
    end do
    !
  end do filt_cell_loop
  !
end subroutine emergency_limit_bad_cells_2d
!
!###############################################################################
!> Apply the sponge filter to bad cells and update the volume solution in usp.
!>
!> Optimized routine for 1d.
!###############################################################################
pure &
subroutine emergency_limit_bad_cells_1d(bad_cell,usp)
  !
  !.. Use Statements ..
  use order_mod,   only : maxSP
  use geomvar_mod, only : ncell,cell
  use filter_mod,  only : filter
  !
  !.. Formal Arguments ..
  logical(lk), intent(in) :: bad_cell(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,filt_ec,filt_mx,filt_ee
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "emergency_limit_bad_cells_1d"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,bad_cell,filter) &
  !$ACC&              PRIVATE(old_usp)
  filt_cell_loop: &
  do nc = 1,ncell
    !
    ! If the current cell is not troubled, skip to the next one
    !
    if (.not.bad_cell(nc)) then
      cycle filt_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      old_usp(n,nec) = usp(n0+n,nec)
      old_usp(n,nmx) = usp(n0+n,nmx)
      old_usp(n,nee) = usp(n0+n,nee)
    end do
    !
    ! Apply the filter emergency limiter by multplying the filter
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR
    do n = 1,np
      !
      filt_ec = zero
      filt_mx = zero
      filt_ee = zero
      !
      do i = filter(this_geom,this_order)%sponge%col_ptr(n)+1, &
             filter(this_geom,this_order)%sponge%col_ptr(n+1)
        !
        k = filter(this_geom,this_order)%sponge%col_idx(i)
        a = filter(this_geom,this_order)%sponge%col_val(i)
        !
        filt_ec = filt_ec + a * old_usp(k,nec)
        filt_mx = filt_mx + a * old_usp(k,nmx)
        filt_ee = filt_ee + a * old_usp(k,nee)
        !
      end do
      !
      usp(n0+n,nec) = filt_ec
      usp(n0+n,nmx) = filt_mx
      usp(n0+n,nee) = filt_ee
      !
    end do
    !
  end do filt_cell_loop
  !
end subroutine emergency_limit_bad_cells_1d
!
!###############################################################################
!> Apply the sponge filter to bad cells and update the volume solution in usp.
!>
!> General routine for nd.
!###############################################################################
pure &
subroutine emergency_limit_bad_cells_nd(mbeg,mend,bad_cell,usp)
  !
  !.. Use Statements ..
  use order_mod,   only : maxSP
  use geomvar_mod, only : ncell,cell
  use filter_mod,  only : filter
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: mbeg,mend
  logical(lk), intent(in) :: bad_cell(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,m,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,val
  !
  !.. Local Arrays ..
  real(wp) :: old_usp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "emergency_limit_bad_cells_nd"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,usp,bad_cell,filter) &
  !$ACC&              PRIVATE(old_usp)
  filt_cell_loop: &
  do nc = 1,ncell
    !
    ! If the current cell is not troubled, skip to the next one
    !
    if (.not.bad_cell(nc)) then
      cycle filt_cell_loop
    end if
    !
    n0 = cell(nc)%beg_sp - 1  ! offset index for sol-pts in this cell
    np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry of this cell
    this_order = cell(nc)%order ! order of this cell
    !
    ! Store a copy of usp for this cell because the individual
    ! solutions at the solution points within the cell will be
    ! altered during the projection computation.
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do m = mbeg,mend
      do n = 1,np
        old_usp(n,m) = usp(n0+n,m)
      end do
    end do
    !
    ! Apply the filter emergency limiter by multplying the filter
    ! matrix with the solution in the current cell
    !
    !$ACC LOOP VECTOR COLLAPSE(2)
    do m = mbeg,mend
      do n = 1,np
        !
        val = zero
        !
        do i = filter(this_geom,this_order)%sponge%col_ptr(n)+1, &
               filter(this_geom,this_order)%sponge%col_ptr(n+1)
          !
          k = filter(this_geom,this_order)%sponge%col_idx(i)
          a = filter(this_geom,this_order)%sponge%col_val(i)
          !
          val = val + a * old_usp(k,m)
          !
        end do
        !
        usp(n0+n,m) = val
        !
      end do
    end do
    !
  end do filt_cell_loop
  !
end subroutine emergency_limit_bad_cells_nd
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
!###############################################################################
!> Scale down the higher modes of bad cells and
!> update the volume solution in usp.
!>
!> Keep mode 0, but scale modes 1+ by factor=[0:1]. Xu uses 0.5.
!>   q = qave + (q-qave)*factor = factor*usp + (1-factor)*qave
!> Scales the solution towards qave.
!>
!> Ref:
!>   Xu (2009) doi:10.1016/j.jcp.2008.11.025
!###############################################################################
subroutine emergency_scale_bad_cells(bad_cell,usp)
  !
  !.. Use Statements ..
  use geomvar_mod,  only : ncell,cell
  use thermo_mod,   only : update_thermo_properties
  !
  !.. Arguments ..
  !> Marks the cells to be repaired.
  logical(lk), intent(in)    :: bad_cell(:)
  !
  !> Conserved variables at solution points for all cells.
  real(wp),    intent(inout) :: usp(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "emergency_scale_bad_cells"
  !
  !.. Local Scalars ..
  integer     :: nc,nr,n1,n2
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = nme-nmb+1
  !
  if (force_nd_subroutines) then
    !
    call emergency_scale_bad_cells_nd(1,nq,bad_cell,usp)
    !
  else
    !
    if (nr == 3) then
      call emergency_scale_bad_cells_3d(bad_cell,usp)
    else if (nr == 2) then
      call emergency_scale_bad_cells_2d(bad_cell,usp)
    else
      call emergency_scale_bad_cells_1d(bad_cell,usp)
    end if
    !
    if (ncs > 0) then
      call emergency_scale_bad_cells_nd(1,ncs,bad_cell,usp)
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine emergency_scale_bad_cells
!
!###############################################################################
!> Scale down the higher modes of bad cells and
!> update the volume solution in usp.
!>
!> Optimized routine for 3d.
!###############################################################################
subroutine emergency_scale_bad_cells_3d(bad_cell,usp)
  !
  !.. Use Statements ..
  use geomvar_mod,  only : ncell,cell
  use thermo_mod,   only : update_thermo_properties
  use cell_ave_mod, only : cell_ave_sol
  use ovar,         only : thermofix
  !
  !.. Arguments ..
  !> Marks the cells to be repaired.
  logical(lk), intent(in)    :: bad_cell(:)
  !
  !> Conserved variables at solution points for all cells.
  real(wp),    intent(inout) :: usp(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "emergency_scale_bad_cells_3d"
  !
  !.. Local Scalars ..
  integer  :: n,nc,n0,np
  real(wp) :: f1  ! solution scaling factor
  real(wp) :: f2  ! f2=1-f1
  !
continue
  !
  ! Scale the solution.
  !
  f1 = thermofix%scale_factor
  f2 = one - f1
  !
  loop_cells: &
  do nc = 1,ncell
    if ( bad_cell(nc) ) then
      n0 = cell(nc)%beg_sp - 1    ! offset index for sol-pts in this cell
      np = cell(nc)%end_sp - n0   ! number of sol-pts in this cell
      do n = 1,np  ! loop solution points
        usp(n0+n,nec) = f1*usp(n0+n,nec) + f2*cell_ave_sol(nec,nc)
        usp(n0+n,nmx) = f1*usp(n0+n,nmx) + f2*cell_ave_sol(nmx,nc)
        usp(n0+n,nmy) = f1*usp(n0+n,nmy) + f2*cell_ave_sol(nmy,nc)
        usp(n0+n,nmz) = f1*usp(n0+n,nmz) + f2*cell_ave_sol(nmz,nc)
        usp(n0+n,nee) = f1*usp(n0+n,nee) + f2*cell_ave_sol(nee,nc)
      end do ! n
    end if
  end do loop_cells
  !
end subroutine emergency_scale_bad_cells_3d
!
!###############################################################################
!> Scale down the higher modes of bad cells and
!> update the volume solution in usp.
!>
!> Optimized routine for 2d.
!###############################################################################
subroutine emergency_scale_bad_cells_2d(bad_cell,usp)
  !
  !.. Use Statements ..
  use geomvar_mod,  only : ncell,cell
  use thermo_mod,   only : update_thermo_properties
  use cell_ave_mod, only : cell_ave_sol
  use ovar,         only : thermofix
  !
  !.. Arguments ..
  !> Marks the cells to be repaired.
  logical(lk), intent(in)    :: bad_cell(:)
  !
  !> Conserved variables at solution points for all cells.
  real(wp),    intent(inout) :: usp(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "emergency_scale_bad_cells_2d"
  !
  !.. Local Scalars ..
  integer  :: n,nc,n0,np
  real(wp) :: f1  ! solution scaling factor
  real(wp) :: f2  ! f2=1-f1
  !
continue
  !
  ! Scale the solution.
  !
  f1 = thermofix%scale_factor
  f2 = one - f1
  !
  loop_cells: &
  do nc = 1,ncell
    if ( bad_cell(nc) ) then
      n0 = cell(nc)%beg_sp - 1    ! offset index for sol-pts in this cell
      np = cell(nc)%end_sp - n0   ! number of sol-pts in this cell
      do n = 1,np  ! loop solution points
        usp(n0+n,nec) = f1*usp(n0+n,nec) + f2*cell_ave_sol(nec,nc)
        usp(n0+n,nmx) = f1*usp(n0+n,nmx) + f2*cell_ave_sol(nmx,nc)
        usp(n0+n,nmy) = f1*usp(n0+n,nmy) + f2*cell_ave_sol(nmy,nc)
        usp(n0+n,nee) = f1*usp(n0+n,nee) + f2*cell_ave_sol(nee,nc)
      end do ! n
    end if
  end do loop_cells
  !
end subroutine emergency_scale_bad_cells_2d
!
!###############################################################################
!> Scale down the higher modes of bad cells and
!> update the volume solution in usp.
!>
!> Optimized routine for 1d.
!###############################################################################
subroutine emergency_scale_bad_cells_1d(bad_cell,usp)
  !
  !.. Use Statements ..
  use geomvar_mod,  only : ncell,cell
  use thermo_mod,   only : update_thermo_properties
  use cell_ave_mod, only : cell_ave_sol
  use ovar,         only : thermofix
  !
  !.. Arguments ..
  !> Marks the cells to be repaired.
  logical(lk), intent(in)    :: bad_cell(:)
  !
  !> Conserved variables at solution points for all cells.
  real(wp),    intent(inout) :: usp(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "emergency_scale_bad_cells_1d"
  !
  !.. Local Scalars ..
  integer  :: n,nc,n0,np
  real(wp) :: f1  ! solution scaling factor
  real(wp) :: f2  ! f2=1-f1
  !
continue
  !
  ! Scale the solution.
  !
  f1 = thermofix%scale_factor
  f2 = one - f1
  !
  loop_cells: &
  do nc = 1,ncell
    if ( bad_cell(nc) ) then
      n0 = cell(nc)%beg_sp - 1    ! offset index for sol-pts in this cell
      np = cell(nc)%end_sp - n0   ! number of sol-pts in this cell
      do n = 1,np  ! loop solution points
        usp(n0+n,nec) = f1*usp(n0+n,nec) + f2*cell_ave_sol(nec,nc)
        usp(n0+n,nmx) = f1*usp(n0+n,nmx) + f2*cell_ave_sol(nmx,nc)
        usp(n0+n,nee) = f1*usp(n0+n,nee) + f2*cell_ave_sol(nee,nc)
      end do ! n
    end if
  end do loop_cells
  !
end subroutine emergency_scale_bad_cells_1d
!
!###############################################################################
!> Scale down the higher modes of bad cells and
!> update the volume solution in usp.
!>
!> General routine for nd.
!###############################################################################
subroutine emergency_scale_bad_cells_nd(mbeg,mend,bad_cell,usp)
  !
  !.. Use Statements ..
  use geomvar_mod,  only : ncell,cell
  use thermo_mod,   only : update_thermo_properties
  use cell_ave_mod, only : cell_ave_sol
  use ovar,         only : thermofix
  !
  !.. Arguments ..
  integer,     intent(in)    :: mbeg        !< Beginning solution variable.
  integer,     intent(in)    :: mend        !< Ending    solution variable.
  logical(lk), intent(in)    :: bad_cell(:) !< Marks the cells to be repaired.
  !
  !> Conserved variables at solution points for all cells.
  real(wp),    intent(inout) :: usp(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "emergency_scale_bad_cells_nd"
  !
  !.. Local Scalars ..
  integer  :: m,n,nc,n0,np
  real(wp) :: f1  ! solution scaling factor
  real(wp) :: f2  ! f2=1-f1
  !
continue
  !
  ! Scale the solution.
  !
  f1 = thermofix%scale_factor
  f2 = one - f1
  !
  loop_cells: &
  do nc = 1,ncell
    if ( bad_cell(nc) ) then
      n0 = cell(nc)%beg_sp - 1    ! offset index for sol-pts in this cell
      np = cell(nc)%end_sp - n0   ! number of sol-pts in this cell
      do m = mbeg,mend  ! loop solution variables
        do n = 1,np  ! loop solution points
          usp(n0+n,m) = f1*usp(n0+n,m) + f2*cell_ave_sol(m,nc)
        end do ! n
      end do ! m
    end if
  end do loop_cells
  !
end subroutine emergency_scale_bad_cells_nd
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
!###############################################################################
!> Positivity Preserving limiting procedure of Sheshadri (thesis stanford)
!>
!> Scale density to a threshold, update solution
!> Scale conserved variables based off of a pressure threshold
!> *scale about the cell average
!>
!> Ref:
!> stacks.stanford.edu/file/nk877ng0918/Abhishek_thesis_e_version2-augmented.pdf
!###############################################################################
subroutine emergency_limit2_bad_cells(bad_cell,usp,gpsp)
  !
  !.. Use Statements ..
  use geomvar_mod,  only : ncell,cell,face
  use cell_ave_mod, only : cell_ave_sol, compute_cell_ave_var
  use ovar,         only : thermofix, rgasref_nd, rhoref, aref
  use flowvar_mod,  only : faceusp,facegpsp
  !
  !.. Arguments ..
  logical(lk), intent(in)    :: bad_cell(:) !< Marks the cells to be repaired.
  real(wp),    intent(inout) :: gpsp(:,:)
  !
  !> Conserved variables at solution points for all cells.
  real(wp),    intent(inout) :: usp(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "emergency_limit2_bad_cells"
  !
  !.. Local Scalars ..
  integer  :: m,n,nc,n1,n2,ncf,nf,f1,f2,sf
  real(wp) :: rhomin, pmin ! minimum density/pressure
  real(wp) :: rhoave ! average density
  real(wp), allocatable, dimension(:) :: pave ! average pressure
  real(wp), dimension(ncs) :: ysp
  real(wp) :: theta,tau ! scaling param
  real(wp) :: asdf
  real(wp) :: rholimit, plimit 
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  allocate( pave(ncell), source=zero )
  !
  rholimit = thermofix%min_density / rhoref
  plimit = thermofix%min_pressure / (rhoref*aref*aref)
  !
  loop_cells: &
  do nc = 1,ncell
    if ( bad_cell(nc) ) then
      !
      n1 = cell(nc)%beg_sp ! beginning solution point index for this cell
      n2 = cell(nc)%end_sp ! ending    solution point index for this cell
      !
      ! Get the minimum value of density from the solution points
      !
      rhomin = minval( usp(n1:n2,nec) )
      !
      ! Get the minimum value of density from the face/flux points
      !
      do ncf = 1,size(cell(nc)%face)
        nf = cell(nc)%face(ncf)%idx  ! index of current cell on face
        sf = cell(nc)%face(ncf)%side ! side of face that nc is on
        f1 = face(nf)%beg_fp ! beginning face point index for this face
        f2 = face(nf)%end_fp ! ending    face point index for this face
        rhomin = min( rhomin , minval(faceusp(f1:f2,nec,sf)) )
      end do
      !
      ! Extract the cell average density
      !
      rhoave = max(rholimit,cell_ave_sol(nec,nc))
      !
      ! Find the theta_1 scaling parameter
      theta = min(one, (rhoave - rholimit) / (rhoave - rhomin))
      !
      ! Scale the conserved variables to accound for the new density
      do n = n1,n2  ! loop solution points
        !
        if (ncs > 0) then
          ! calculate mass fractions
          !do m = 1,ncs
          !  ysp(m) = usp(n,m)/max(eps6, usp(n,nec))
          !end do
          !asdf = sum(ysp)
          !do m = 1,ncs
          !  ysp(m) = ysp(m) / asdf
          !end do
          !
          !do m = 1,ncs  ! loop solution variables
          !  ! limit species density
          !  usp(n,m) = theta*(usp(n,m)-rhoave*ysp(m)) + rhoave*ysp(m)
          !  !
          !end do
        end if
        !
        if (ncs > 0) then
          ! realculate species density
          do m = 1,ncs
            usp(n,m) = usp(n,m)*(theta*(usp(n,nec)-rhoave) + rhoave)/usp(n,nec)
          end do
        end if
        !
        ! limit mixture density
        usp(n,nec) = theta*(usp(n,nec)-rhoave) + rhoave
      end do ! n
      !
    end if
  end do loop_cells
  !
  call debug_timer(leaving_procedure,pname)
  return
  !
  call bad_cell_update_properties_sp(bad_cell,usp,gpsp)
  !
  call reinterpolate_bad_cells(bad_cell,usp,faceusp,facegpsp)
  !
  call compute_cell_ave_var(nq+1, pave)
  !
  ! Now limit based on pressure (i.e. entropy limiting)
  !
  loop_cells2: &
  do nc = 1,ncell
    if ( bad_cell(nc) ) then
      !
      n1 = cell(nc)%beg_sp ! beginning solution point index for this cell
      n2 = cell(nc)%end_sp ! ending    solution point index for this cell
      !
      ! Get the minimum value of pressure from the solution points
      !
      pmin = minval( gpsp(n1:n2,npres) )
      !
      ! Get the minimum value of pressure from the face/flux points
      !
      do ncf = 1,size(cell(nc)%face)
        nf = cell(nc)%face(ncf)%idx  ! index of current cell on face
        sf = cell(nc)%face(ncf)%side ! side of face that nc is on
        f1 = face(nf)%beg_fp ! beginning face point index for this face
        f2 = face(nf)%end_fp ! ending    face point index for this face
        pmin = min( pmin , minval(facegpsp(f1:f2,npres,sf)) )
      end do
      !
      ! tau parameter
      tau = min(zero, pmin-plimit)
      !
      theta = tau / (tau - (pave(nc)-plimit)) 
      !
      ! scale solution
      do n = n1,n2  ! loop solution points
        do m = 1,nq  ! loop solution variables
          usp(n,m) = usp(n,m) + theta*(cell_ave_sol(m,nc) - usp(n,m))
        end do ! m
      end do ! n
      !
    end if
  end do loop_cells2
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine emergency_limit2_bad_cells
!
!###############################################################################
!> Update thermo and transport properties at volume points of bad cells.
!###############################################################################
subroutine bad_cell_update_properties_sp(bad_cell,usp,gpsp)
  !
  !.. Use Statements ..
  use geomvar_mod,   only : ncell,cell
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  !.. Formal Arguments ..
  logical(lk), intent(in) :: bad_cell(:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: gpsp(:,:)
  !
  !.. Local Scalars ..
  integer :: nc,n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "bad_cell_update_properties"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  !
  ! Update the thermodynamic properties for any cells that had to be limited
  ! Update the transport     properties for any cells that had to be limited
  !
  do nc = 1,ncell
    if (bad_cell(nc)) then
      n1 = cell(nc)%beg_sp
      n2 = cell(nc)%end_sp
      call update_thermo_properties(n1,n2,usp,gpsp, &
                                    loc="cell",caller_name=pname)
      call update_transport_properties(n1,n2,usp,gpsp)
    end if
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine bad_cell_update_properties_sp
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
#include "Functions/debug_timer.f90"
!
!###############################################################################
!
end module thermofix_mod
