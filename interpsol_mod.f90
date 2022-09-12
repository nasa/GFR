module interpsol_mod
  !
  ! Module to handle interpolation of (grid1,sol1) -> (grid2,sol2).
  !
!#define PROFILE_ON
  !
  !.. Use Statements ..
  use module_kind_types, global_debug_timer =>      debug_timer, &
                    pure_global_debug_timer => pure_debug_timer
  use interpsolvar_mod
  use eqn_idx, only : nec,nmx,nmy,nmz,nee,nq
  !
  implicit none
  !
  private
  !
  public :: mv_global_to_grid2
  public :: mv_grid2_to_global
  public :: interpsol
  public :: interpsol_debug
  !
contains
!
!###############################################################################
!
subroutine mv_global_to_grid2()
  !
  ! Move the global variables to grid2, leaving globals deallocated.
  !
  !.. Use Statements ..
  use geovar
  use flowvar
  use ovar,   only : bc_names
  !
  use element_type_mod,  only : std_elem,face_elem
  use vandermonde_mod,   only : vand
  use derivatives_mod,   only : deriv
  use correction_mod,    only : correct
  use interpolation_mod, only : interp,outerp
  use projection_mod,    only : projct
  use filter_mod,        only : filter
  !
  use metrics_mod, only : facenrm,facejac
  use metrics_mod, only : jac,wtdjac
  use metrics_mod, only : metrics_array
  !
  use order_mod, only : n_order,p_order,q_order,o_order,e_order,t_order
  !
continue
  !
  !-DEBUG_BEG:
  write (iout,*) "DEBUG[mv_global_to_grid2] Move global vars to grid2."
  flush (iout)
  !-DEBUG_END:
  !
  ! Variables from order_mod.
  g2_n_order = n_order
  g2_p_order = p_order
  g2_q_order = q_order
  g2_o_order = o_order
  g2_e_order = e_order
  g2_t_order = t_order
  !
  ! Variables from geovar.
  g2_n_solpts         = n_solpts
  g2_n_flxpts         = n_flxpts
  g2_using_quadrature = using_quadrature
  g2_nr               = nr
  g2_nnode            = nnode
  g2_ncell            = ncell
  g2_nface            = nface
  g2_nedge            = nedge
  g2_nfbnd            = nfbnd
  g2_n_totcel         = n_totcel
  g2_n_totnod         = n_totnod
  g2_n_totpts         = n_totpts
  !
  call move_alloc( from=cell_geom,  to=g2_cell_geom  )
  call move_alloc( from=cell_order, to=g2_cell_order )
  call move_alloc( from=xyz_nodes,  to=g2_xyz_nodes  )
  call move_alloc( from=xyz,        to=g2_xyz        )
  !
  g2_n_global_node    = n_global_node
  g2_n_global_totnod  = n_global_totnod
  g2_n_global_cell    = n_global_cell
  g2_n_global_totcel  = n_global_totcel
  g2_n_global_edge    = n_global_edge
  g2_n_global_fbnd    = n_global_fbnd
  g2_n_global_solpts  = n_global_solpts
  g2_n_global_totpts  = n_global_totpts
  !
  call move_alloc( from=global_xyz_nodes,  to=g2_global_xyz_nodes  )
  call move_alloc( from=global_xyz,        to=g2_global_xyz        )
  call move_alloc( from=global_nofc_ptr,   to=g2_global_nofc_ptr   )
  call move_alloc( from=global_nofc,       to=g2_global_nofc       )
  call move_alloc( from=global_pinc_ptr,   to=g2_global_pinc_ptr   )
  call move_alloc( from=global_cell_geom,  to=g2_global_cell_geom  )
  call move_alloc( from=global_cell_order, to=g2_global_cell_order )
  call move_alloc( from=bface,             to=g2_bface             )
  !
  call move_alloc( from=bc_conflict, to=g2_bc_conflict )
  !
  call move_alloc( from=nodes_of_cell_ptr, to=g2_nodes_of_cell_ptr )
  call move_alloc( from=nodes_of_cell,     to=g2_nodes_of_cell     )
  !
  call move_alloc( from=pst_geom, to=g2_pst_geom )
  !
  call move_alloc( from=fp, to=g2_fp )
  !
  call move_alloc( from=face, to=g2_face )
  call move_alloc( from=cell, to=g2_cell )
  !
  call move_alloc( from=grid, to=g2_grid )
  !
  call move_alloc( from=cell_map, to=g2_cell_map )
  !
  call move_alloc( from=bc_grps, to=g2_bc_grps )
  !
  ! Variables from flowvar.
  call move_alloc( from=usp,           to=g2_usp           )
  call move_alloc( from=gpsp,          to=g2_gpsp          )
  call move_alloc( from=dusp,          to=g2_dusp          )
 !call move_alloc( from=facexyz_dt,    to=g2_facexyz_dt    )
  call move_alloc( from=facexyz_array, to=g2_facexyz_array )
  !
  ! Variables from ovar.
  call move_alloc( from=bc_names, to=g2_bc_names )
  !
  ! Variables from quadrature_mod
  call move_alloc( from=std_elem,  to=g2_std_elem  )
  call move_alloc( from=face_elem, to=g2_face_elem )
  !
  ! Variables from vandermonde_mod
  call move_alloc( from=vand, to=g2_vand )
  !
  ! Variables from derivatives_mod
  call move_alloc( from=deriv, to=g2_deriv )
  !
  ! Variables from correction_mod
  call move_alloc( from=correct, to=g2_correct )
  !
  ! Variables from interpolation_mod
  call move_alloc( from=interp, to=g2_interp )
  call move_alloc( from=outerp, to=g2_outerp )
  !
  ! Variables from projection_mod
  call move_alloc( from=projct, to=g2_projct )
  !
  ! Variables from filter_mod
  call move_alloc( from=filter, to=g2_filter )
  !
  ! Variables from metrics_mod
  call move_alloc( from=metrics_array, to=g2_metrics_array )
  call move_alloc( from=jac,           to=g2_jac           )
  call move_alloc( from=wtdjac,        to=g2_wtdjac        )
  call move_alloc( from=facenrm,       to=g2_facenrm       )
  call move_alloc( from=facejac,       to=g2_facejac       )
  !
  !-DEBUG_BEG:
  write (iout,*) "DEBUG[mv_global_to_grid2] Leaving subroutine."
  flush (iout)
  !-DEBUG_END:
  !
  !
end subroutine mv_global_to_grid2
!
!###############################################################################
!
subroutine mv_grid2_to_global()
  !
  ! Move the grid2 variables to global, leaving grid2 deallocated.
  !
  !.. Use Statements ..
  use geovar
  use flowvar
  use ovar,    only : bc_names
  !
  use element_type_mod,  only : std_elem,face_elem
  use vandermonde_mod,   only : vand
  use derivatives_mod,   only : deriv
  use correction_mod,    only : correct
  use interpolation_mod, only : interp,outerp
  use projection_mod,    only : projct
  use filter_mod,        only : filter
  !
  use metrics_mod, only : facenrm,facejac
  use metrics_mod, only : jac,wtdjac
  use metrics_mod, only : metrics_array
  !
  use order_mod, only : n_order,p_order,q_order
  use order_mod, only : o_order,e_order,t_order
  use order_mod, only : initialize_geom_solpts
  use order_mod, only : find_valid_geometries
  !
continue
  !
  !-DEBUG_BEG:
  write (iout,*) "DEBUG[mv_grid2_to_global] Move grid2 vars to global."
  flush (iout)
  !-DEBUG_END:
  !
  n_order = g2_n_order
  p_order = g2_p_order
  q_order = g2_q_order
  o_order = g2_o_order
  e_order = g2_e_order
  t_order = g2_t_order
  !
  n_solpts         = g2_n_solpts
  n_flxpts         = g2_n_flxpts
  using_quadrature = g2_using_quadrature
  nr               = g2_nr
  nnode            = g2_nnode
  ncell            = g2_ncell
  nface            = g2_nface
  nedge            = g2_nedge
  nfbnd            = g2_nfbnd
  n_totcel         = g2_n_totcel
  n_totnod         = g2_n_totnod
  n_totpts         = g2_n_totpts
  !
  call move_alloc( from=g2_cell_geom,  to=cell_geom  )
  call move_alloc( from=g2_cell_order, to=cell_order )
  !
  call initialize_geom_solpts(fals)
  call find_valid_geometries()
  !
  call move_alloc( from=g2_xyz_nodes,  to=xyz_nodes  )
  call move_alloc( from=g2_xyz,        to=xyz        )
  !
  n_global_node    = g2_n_global_node
  n_global_totnod  = g2_n_global_totnod
  n_global_cell    = g2_n_global_cell
  n_global_totcel  = g2_n_global_totcel
  n_global_edge    = g2_n_global_edge
  n_global_fbnd    = g2_n_global_fbnd
  n_global_solpts  = g2_n_global_solpts
  n_global_totpts  = g2_n_global_totpts
  !
  call move_alloc( from=g2_global_xyz_nodes,  to=global_xyz_nodes  )
  call move_alloc( from=g2_global_xyz,        to=global_xyz        )
  call move_alloc( from=g2_global_nofc_ptr,   to=global_nofc_ptr   )
  call move_alloc( from=g2_global_nofc,       to=global_nofc       )
  call move_alloc( from=g2_global_pinc_ptr,   to=global_pinc_ptr   )
  call move_alloc( from=g2_global_cell_geom,  to=global_cell_geom  )
  call move_alloc( from=g2_global_cell_order, to=global_cell_order )
  call move_alloc( from=g2_bface,             to=bface             )
  !
  call move_alloc( from=g2_bc_conflict, to=bc_conflict )
  !
  call move_alloc( from=g2_nodes_of_cell_ptr, to=nodes_of_cell_ptr )
  call move_alloc( from=g2_nodes_of_cell,     to=nodes_of_cell     )
  !
  call move_alloc( from=g2_pst_geom, to=pst_geom )
  !
  call move_alloc( from=g2_fp, to=fp )
  !
  call move_alloc( from=g2_face, to=face )
  call move_alloc( from=g2_cell, to=cell )
  !
  call move_alloc( from=g2_grid, to=grid )
  !
  call move_alloc( from=g2_cell_map, to=cell_map )
  !
  call move_alloc( from=g2_bc_grps, to=bc_grps )
  !
  !
  !-DEBUG_BEG:
  write (iout,*) "DEBUG_BEG[mv_grid2_to_global]"
  write (iout,*) "  Before move_alloc(from=g2_usp,to=usp):"
  write (iout,2001) "    g2_usp(1,1:nq) = ",g2_usp(1,1:nq)
  write (iout,2001) "       usp(1,1:nq) = ",   usp(1,1:nq)
  write (iout,2002) "g2_usp", &
                    lbound(g2_usp,dim=1),ubound(g2_usp,dim=1), &
                    lbound(g2_usp,dim=2),ubound(g2_usp,dim=2)
  write (iout,2002) "usp", &
                    lbound(usp,dim=1),ubound(usp,dim=1), &
                    lbound(usp,dim=2),ubound(usp,dim=2)
  flush (iout)
  !-DEBUG_END:
  !
  ! Variables from flowvar.
  call move_alloc( from=g2_usp,           to=usp           )
  call move_alloc( from=g2_gpsp,          to=gpsp          )
  call move_alloc( from=g2_dusp,          to=dusp          )
 !call move_alloc( from=g2_facexyz_dt,    to=facexyz_dt    )
  call move_alloc( from=g2_facexyz_array, to=facexyz_array )
  !
  !-DEBUG_BEG:
  write (iout,*) "  After  move_alloc(from=g2_usp,to=usp):"
  !
  if (allocated(g2_usp)) then
    write (iout,2001) "    g2_usp(1,1:nq) = ",g2_usp(1,1:nq)
    write (iout,2002) "g2_usp", &
                      lbound(g2_usp,dim=1),ubound(g2_usp,dim=1), &
                      lbound(g2_usp,dim=2),ubound(g2_usp,dim=2)
  else
    write (iout,*)    "    g2_usp(1,1:nq) = unallocated"
  end if
  !
  if (allocated(usp)) then
    write (iout,2001) "       usp(1,1:nq) = ",   usp(1,1:nq)
    write (iout,2002) "usp", &
                      lbound(usp,dim=1),ubound(usp,dim=1), &
                      lbound(usp,dim=2),ubound(usp,dim=2)
  else
    write (iout,*)    "    usp(1,1:nq) = unallocated"
  end if
  !
  write (iout,*) "DEBUG_END[mv_grid2_to_global]"
  flush (iout)
  !-DEBUG_END:
  !
  ! Variables from ovar.
  call move_alloc( from=g2_bc_names, to=bc_names )
  !
  ! Variables from quadrature_mod
  call move_alloc( from=g2_std_elem,  to=std_elem  )
  call move_alloc( from=g2_face_elem, to=face_elem )
  !
  ! Variables from vandermonde_mod
  call move_alloc( from=g2_vand, to=vand )
  !
  ! Variables from derivatives_mod
  call move_alloc( from=g2_deriv, to=deriv )
  !
  ! Variables from correction_mod
  call move_alloc( from=g2_correct, to=correct )
  !
  ! Variables from interpolation_mod
  call move_alloc( from=g2_interp, to=interp )
  call move_alloc( from=g2_outerp, to=outerp )
  !
  ! Variables from projection_mod
  call move_alloc( from=g2_projct, to=projct )
  !
  ! Variables from filter_mod
  call move_alloc( from=g2_filter, to=filter )
  !
  ! Variables from metrics_mod
  call move_alloc( from=g2_metrics_array, to=metrics_array )
  call move_alloc( from=g2_jac,           to=jac           )
  call move_alloc( from=g2_wtdjac,        to=wtdjac        )
  call move_alloc( from=g2_facenrm,       to=facenrm       )
  call move_alloc( from=g2_facejac,       to=facejac       )
  !
  !-DEBUG_BEG:
  write (iout,*) "DEBUG[mv_grid2_to_global] Leaving subroutine."
  flush (iout)
  !-DEBUG_END:
  !
  ! Format Statements
  !
  2001 format (x,a,e20.6,:,/,(20x,e20.6))
  2002 format ("  Bounds of ",a," = (",i0,":",i0,",",i0,":",i0,")")
  !
end subroutine mv_grid2_to_global
!
!###############################################################################
!
subroutine partition_usp_to_g2()
  !
  !.. Use Statements ..
  use geovar,  only : cell
  use flowvar, only : usp
  !
  !.. Local Scalars ..
  integer :: g1,g2,l1,l2,mypart,loc_cell,glb_cell
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "partition_usp_to_g2"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  mypart = mypnum + 1
  !
  do loc_cell = 1,size(g2_cell_map(mypart)%loc_to_glb)
    !
    glb_cell = g2_cell_map(mypart)%loc_to_glb(loc_cell)
    !
    l1 = g2_cell(loc_cell)%beg_sp
    l2 = g2_cell(loc_cell)%end_sp
    !
    g1 = cell(glb_cell)%beg_sp
    g2 = cell(glb_cell)%end_sp
    !
    g2_usp(l1:l2,1:nq) = usp(g1:g2,1:nq)
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine partition_usp_to_g2
!
!###############################################################################
!
subroutine interpsol()
  !
  ! Interpolate solution onto grid2.
  !
  !.. Use Statements ..
  use module_kind_types, only : dot,lk,fals
  use geovar,            only : ncell,nface,nfbnd,nr
  use geovar,            only : cell,face,fp
  use flux_mod,          only : interp_cv_to_flxpts  ! DEBUG: for debug output
  use metrics_mod,       only : geofa,facenrm
  use flowvar,           only : usp
 !use flowvar,           only : facexyz_dt
  use flowvar,           only : facexyz_array
  use order_mod,         only : maxSP,geom_solpts
  use octedat_facexyz_mod
  use octtree_mod
  use octstat_mod
  use interpolation_mod, only : interp_cell_sol_to_xyz, interp
  use interpolation_mod, only : init_ref_vectors
  !-DEBUG_BEG:
 !use ifport,            only : sleep
  use mydebug_mod
  !-DEBUG_END:
#ifdef ENABLE_OPENMP
  use omp_lib, only : omp_get_thread_num
  use omp_lib, only : omp_get_num_threads
#endif
  !
  !.. Local Scalars ..
  type(octedat_facexyz_type), target :: octedat
  type(octtree_type)                 :: octtree
  !
  integer :: istat,acc
  integer :: i,j,ifp,k,l,nc,ncf,nf,nfp,sf,f0
  integer :: n0,n1,n2,np,knp,m,lf
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: fp_offset,face_rot
  integer :: donor_cell
  !
  real(wp) :: a,val,face_dot
  real(wp) :: valsp1,valsp2,valsp3,valsp4,valsp5,valsp6
  real(wp) :: valec,valmx,valmy,valmz,valee
  !
  logical(lk) :: dealloc_octedat
  !
  !-DEBUG_BEG:
  integer :: lcell_min, lcell_max, rcell_min, rcell_max
  integer :: size_cell, size_face
  integer :: size_cell2, size_face2
  logical :: debug_print
  integer :: hours,minutes,seconds,pts_completed
  real(wp) :: start_time,new_time,elapsed_time,percent_complete
  !-DEBUG_END:
  !
  !.. Local Arrays ..
  integer  :: ixyz_pt(1:3) ! = (cell,face,facept)
  real(wp) :: xyz_pt(1:nr)
  real(wp) :: g2_xyz_pt(1:nr)
  real(wp) :: face_to_xyz(1:nr)
  real(wp) :: bbox(1:nr,1:2)
  real(wp) :: tusp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interpsol"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
 !if (ncpu > 1) then
 !  call partition_usp_to_g2
 !  call debug_timer(leaving_procedure,pname)
 !  return
 !end if
  !
  !-DEBUG_BEG:
  ! Need faceusp for debug output.  Must compute it first.
  call interp_cv_to_flxpts()
  ! Print global array info for debugging.
  call mydebug(pname,__LINE__,__FILE__)
  !-DEBUG_END:
  !
  ! Size the bounding box to contain the dst grid volume points.
  !
  if (mypnum == glb_root) then
    write (iout,10) pname, "Size bounding box"  ! DEBUG
    flush (iout)
  end if
  !
  bbox(:,1) = g2_xyz(:,1)
  bbox(:,2) = g2_xyz(:,1)
  !
  do k = 2,size(g2_xyz,dim=2)
    do l = 1,nr
      bbox(l,1) = min( bbox(l,1), g2_xyz(l,k) )
      bbox(l,2) = max( bbox(l,2), g2_xyz(l,k) )
    end do
  end do
  !
  ! Increase the bounding box to also contain the donor grid face points.
  !
  if (mypnum == glb_root) then
    write (iout,10) pname, "Resize bounding box"  ! DEBUG
    flush (iout)
  end if
  !
  do nf = 1,nface
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
   !nfp = size(facexyz_dt(nf)%v,dim=2)
    !
    do ifp = 1,nfp
      do l = 1,nr
       !bbox(l,1) = min( bbox(l,1), facexyz_dt(nf)%v(l,ifp) )
       !bbox(l,2) = max( bbox(l,2), facexyz_dt(nf)%v(l,ifp) )
        bbox(l,1) = min( bbox(l,1), facexyz_array(l,f0+ifp) )
        bbox(l,2) = max( bbox(l,2), facexyz_array(l,f0+ifp) )
      end do
    end do
    !
  end do
  !
  ! Init the encapsulated data for the octtree.
  !
  if (mypnum == glb_root) then
    write (iout,10) pname, "Init octedat"  ! DEBUG
    flush (iout)
  end if
  !
 !octedat%varname = "facexyz_dt(1:nface)%v(1:nr,1:nfp)"
  octedat%varname = "facexyz_array(1:nr,1:n_flxpts)"
  octedat%vardesc = "Get facexyz given cell,face,facept = (nc,ncf,nfp)."
 !octedat%facexyz => facexyz
 !octedat%cell    => cell
  octedat%MAXINDICES = 3
  !
  ! Encode octedat data using 3 indices:
  !   idx(1) = nc   ! cell no.
  !   idx(2) = ncf  ! cell face no.
  !   idx(3) = ifp  ! cell face facept no.
  !
  ! Init the empty octtree.
  !
  if (mypnum == glb_root) then
    write (iout,10) pname, "Init empty octtree"  ! DEBUG
    flush (iout)
  end if
  !
  octtree%octedat => octedat
  call octtree_init(octtree,bbox,istat)
  !
  !-DEBUG_BEG:
  if (mypnum == glb_root) then
    write (iout,10) pname, "Filling octree for use in interpolation"
    flush (iout)
  end if
  !-DEBUG_END:
  !
  ! Fill the octree with the donor grid face points.
  !
  do nc = 1,ncell  ! cell loop
    !
    !-DEBUG_BEG:
    !-DEBUG_CMT: Print percent complete status every 10%.
    if (mypnum == glb_root) then
      if (nc > 1 .and. (ncell-1)*10 > 100) then
        if (mod( nc-1, ((ncell-1)*10)/100 ) == 0) then
          write (iout,201) int( 100.0*(nc-1.0)/(ncell-1.0) + 0.5 )
          flush (iout)
          201 format ("Octree fill ",i3," percent complete.")
        end if
      end if
    end if
    !-DEBUG_END:
    !
    do ncf = 1,size(cell(nc)%face)  ! cell face loop
      !
      nf = cell(nc)%face(ncf)%idx   ! global face index
      sf = cell(nc)%face(ncf)%side
      !
      ! Only left faces to avoid duplicate pts from right face.
      !
      if (sf == 1) then
        !
        face_geom  = cell(nc)%face(ncf)%geom
        face_order = cell(nc)%face(ncf)%order
        nfp = geom_solpts(face_geom,face_order)  ! current cell face pts
        !
        do ifp = 1,nfp  ! current cell face point loop
          !
          !-DEBUG_BEG:
#ifdef DEBUG_OCTFILL
          if (mypnum == glb_root) then
            write (iout,400) pname,nc,ncell,ncf,size(cell(nc)%face),ifp,nfp
            flush (iout)
            400 format (1x,"DEBUG[",a,"]  nc,ncf,ifp = ",3(i4,'/',i4,:,2x))
          end if
#endif
          !-DEBUG_BEG:
          !   call sleep(1)  ! Sleep # seconds.
          !-DEBUG_END:
          !
          ! Load ixyz_pt(:) with the index values the octtree will track.
          !
          ixyz_pt(1) = nc   ! cell no.
          ixyz_pt(2) = ncf  ! cell face no.
          ixyz_pt(3) = ifp  ! cell face facept no.
          f0 = face(nf)%beg_fp - 1
          !
          xyz_pt = facexyz_array(:,f0+ifp)
         !xyz_pt = facexyz_dt(nf)%v(:,ifp)
          !
          call octtree_addxyz(octtree,ixyz_pt,xyz_pt,istat)
          if (istat /= OCT_STATUS_OK) then
            if (mypnum == glb_root) then
              write (iout,*) "ERROR:  Failed to init octree for interpolation."
              write (iout,*) "          status = ",istat
              call stop_gfr(abort,pname,__LINE__,__FILE__)
            end if
            ! Have all the host roots call mpi_barrier here so that the
            ! MPI processes that arent the global root do not continue
            ! past this point before the global root has a chance to
            ! call mpi_abort within the above conditional.
            call mpi_barrier(host_roots_comm,mpierr)
          end if
        end do  ! cell face point loop
      end if  ! Left face
    end do  ! cell face loop
  end do  ! cell loop
  !
  ! Interpolate solution onto each volume point in the destination grid.
  !
  !-DEBUG_BEG:
  if (mypnum == glb_root) then
    write (iout,10) pname, "Octtree fill complete"
    write (iout,*)
    write (iout,10) pname, "Performing interpolation"
    write (iout,*) "DEBUG_BEG:[interpsol]"
    size_cell=0; if (allocated(cell)) size_cell=size(cell)
    size_face=0; if (allocated(face)) size_face=size(face)
    size_cell2=0; if (allocated(g2_cell)) size_cell2=size(g2_cell)
    size_face2=0; if (allocated(g2_face)) size_face2=size(g2_face)
    write (iout,*) "      ncell/nface/nfbnd = ",ncell," / ",nface," / ",nfbnd
    write (iout,*) "   g2_ncell/nface/nfbnd = ",g2_ncell," / ",g2_nface, &
                                                       " / ",g2_nfbnd
    write (iout,*) "        size(cell/face) = ",size_cell," / ",size_face
    write (iout,*) "     g2_size(cell/face) = ",size_cell2," / ",size_face2
    !
    lcell_min = huge(0)
    rcell_min = huge(0)
    lcell_max = 0
    rcell_max = 0
    do nc = 1,ncell
      do ncf = 1,size(cell(nc)%face)
        nf = cell(nc)%face(ncf)%idx
        lcell_min = min(lcell_min, face(nf)%left%cell)
        rcell_min = min(rcell_min, face(nf)%right%cell)
        lcell_max = max(lcell_max, face(nf)%left%cell)
        rcell_max = max(rcell_max, face(nf)%right%cell)
      end do
    end do
    write (iout,*) "       lcell_min/max    = ",lcell_min," / ",lcell_max
    write (iout,*) "       rcell_min/max    = ",rcell_min," / ",rcell_max
    !
    lcell_min = huge(0)
    rcell_min = huge(0)
    lcell_max = 0
    rcell_max = 0
    do nc = 1,g2_ncell
      do ncf = 1,size(g2_cell(nc)%face)
        nf = g2_cell(nc)%face(ncf)%idx
        lcell_min = min(lcell_min, g2_face(nf)%left%cell)
        rcell_min = min(rcell_min, g2_face(nf)%right%cell)
        lcell_max = max(lcell_max, g2_face(nf)%left%cell)
        rcell_max = max(rcell_max, g2_face(nf)%right%cell)
      end do
    end do
    write (iout,*) "    g2_lcell_min/max    = ",lcell_min," / ",lcell_max
    write (iout,*) "    g2_rcell_min/max    = ",rcell_min," / ",rcell_max
    write (iout,*) "    g2_n_solpts         = ",g2_n_solpts
    write (iout,*) " Looping over g2_n_solpts in destination grid."
    !
    flush (iout)
  end if
  !
  ! Initialize the ref_vectors array needed for interpolation
  !
  call init_ref_vectors()
  !
!#define QUICK_SKIP_INTERPOLATE
#ifdef QUICK_SKIP_INTERPOLATE
  do m = 1,nq
    do k = 1,g2_n_solpts
      g2_usp(k,m) = usp(1,m)
    end do
  end do
#else
#define ASSUME_3D
  !-DEBUG_END:
  ! Loop over volume points in dst grid.
 !open (unit=6,form="formatted",carriagecontrol="fortran",buffered="no")
  !$OMP PARALLEL &
  !$OMP&         DEFAULT(SHARED) &
  !$OMP&         PRIVATE(k,g2_xyz_pt,nc,ncf,ifp,nf,face_to_xyz,donor_cell, &
  !$OMP&                 n0,n1,n2,np,this_geom,this_order,knp,m,tusp,tid, &
  !$OMP&                 face_geom,face_order,lf,istat,ixyz_pt,xyz_pt,val, &
  !$OMP&                 start_time,new_time,elapsed_time,percent_complete, &
  !$OMP&                 hours,minutes,seconds,pts_completed,fp_offset, &
  !$OMP&                 face_rot,face_dot, &
  !$OMP&                 valsp1,valsp2,valsp3,valsp4,valsp5,valsp6, &
  !$OMP&                 valec,valmx,valmy,valmz,valee)
  !
  call debug_timer(reset_accumulated_time, &
                   "octtree_find_nearest_xyz",acc=1)
  call debug_timer(reset_accumulated_time, &
                   "LEFT cell interp_cell_sol_to_xyz",acc=2)
  call debug_timer(reset_accumulated_time, &
                   "RIGHT cell interp_cell_sol_to_xyz",acc=3)
  call debug_timer(reset_accumulated_time, &
                   "outside boundary face",acc=4)
  !
  start_time = mpi_wtime()
#ifdef ENABLE_OPENMP
  tid = omp_get_thread_num()
  nthreads = omp_get_num_threads()
#endif
  pts_completed = 0
  chunk_size = ceiling( real(g2_n_solpts,kind=wp) / real(nthreads,kind=wp) )
  !
  if (mypnum == glb_root) then
    if (tid == glb_root) then
      write (iout,thread_report) pname,nthreads
    end if
  end if
  !
  g2_xyz_pt(:) = zero
  ixyz_pt(:) = 0
  call octtree_find_nearest_xyz(octtree,g2_xyz_pt, ixyz_pt,istat, &
                                tid=tid,acc=1)
  !
  !$OMP DO SCHEDULE(RUNTIME)
 !!$OMP DO SCHEDULE(GUIDED,chunk_size/10)
 !!$OMP DO SCHEDULE(STATIC,chunk_size)
  do k = 1,g2_n_solpts
    !
    !-DEBUG_BEG:
    !-DEBUG_CMT: Print percent complete status every 10%.
    if (mypnum == glb_root) then
      if (tid == glb_root) then
        if (pts_completed > 0 .and. (g2_n_solpts-1)*10 > 100) then
         !if (mod( k-1, ((g2_n_solpts-1)*10)/100 ) == 0) then
         !  write (iout,301) int( 100.0*(k-1.0)/(g2_n_solpts-1.0) + 0.5 )
         !  flush (iout)
         !  301 format ("Interpolation ",i3," percent complete.")
         !end if
          if (k == g2_n_solpts) then
            write (iout,303)
            flush (iout)
            303 format (t1,"Interpolation  100.000 percent complete", &
                          " ( 0:00:00 remaining).")
          else if (mod( pts_completed-1, 1000 ) == 0) then
            new_time = mpi_wtime()
            elapsed_time = new_time - start_time
            percent_complete = real(pts_completed,kind=wp) / &
                               real(chunk_size,kind=wp)
            new_time = elapsed_time/percent_complete - elapsed_time
            hours = nint(new_time) / 3600
            minutes = (nint(new_time) - hours*3600) / 60
            seconds = nint(new_time) - hours*3600 - minutes*60
            write (iout,302) 100.0_wp*percent_complete,hours,minutes,seconds
            flush (iout)
            302 format ("Interpolation ",f8.3," percent complete ", &
                       "(",i0,":",i2.2,":",i2.2," remaining).")
          end if
        end if
      end if
    end if
    !-DEBUG_END:
    !
    ! Get volume point xyz in dst grid to interp sol onto.
    !
#ifdef ASSUME_3D
    g2_xyz_pt(1) = g2_xyz(1,k)
    g2_xyz_pt(2) = g2_xyz(2,k)
    g2_xyz_pt(3) = g2_xyz(3,k)
#else
    g2_xyz_pt = g2_xyz(:,k)
#endif
    !
    ! Find index of nearest face pt in src grid.
    !
#ifdef DEBUG_ON
    call debug_timer(start_accumulation,"octtree_find_nearest_xyz",acc=1)
#endif
    call octtree_find_nearest_xyz(octtree,g2_xyz_pt, ixyz_pt,istat, &
                                  tid=tid,acc=0)
#ifdef DEBUG_ON
    call debug_timer(stop_accumulation,"octtree_find_nearest_xyz",acc=1)
#endif
    !
    if (istat == OCT_STATUS_OK) then
      !
      nc  = ixyz_pt(1)  ! cell no.
      ncf = ixyz_pt(2)  ! cell face no.
      ifp = ixyz_pt(3)  ! cell face face_pt no.
      nf  = cell(nc)%face(ncf)%idx   ! global face index
      f0  = face(nf)%beg_fp - 1
      !
#ifdef ASSUME_3D
     !xyz_pt(1) = facexyz_array(1,f0+ifp)
     !xyz_pt(2) = facexyz_array(2,f0+ifp)
     !xyz_pt(3) = facexyz_array(3,f0+ifp)
     !!
     !face_to_xyz(1) = g2_xyz_pt(1) - xyz_pt(1)   ! vector from face to dst pt
     !face_to_xyz(2) = g2_xyz_pt(2) - xyz_pt(2)   ! vector from face to dst pt
     !face_to_xyz(3) = g2_xyz_pt(3) - xyz_pt(3)   ! vector from face to dst pt
     !!
     !face_dot = facenrm(1,f0+ifp)*face_to_xyz(1) + &
     !           facenrm(2,f0+ifp)*face_to_xyz(2) + &
     !           facenrm(3,f0+ifp)*face_to_xyz(3)
      ! 
      ! Dot product of face normal with face_to_xyz vector
      !
      face_dot = facenrm(1,f0+ifp)*(g2_xyz_pt(1) - facexyz_array(1,f0+ifp)) + &
                 facenrm(2,f0+ifp)*(g2_xyz_pt(2) - facexyz_array(2,f0+ifp)) + &
                 facenrm(3,f0+ifp)*(g2_xyz_pt(3) - facexyz_array(3,f0+ifp))
#else
      xyz_pt(:) = facexyz_array(:,f0+ifp)
     !xyz_pt(:) = facexyz_dt(nf)%v(:,ifp)
      face_to_xyz(:) = g2_xyz_pt(:) - xyz_pt(:)   ! vector from face to dst pt
      face_dot = dot( facenrm(:,f0+ifp), face_to_xyz(:) )
#endif
      !
      ! Dot face normal with face_to_xyz vector.
      ! Face normal always points toward right cell.
      !
     !if (dot( geofa(nf)%v(:,ifp), face_to_xyz(:) ) >= 0) then
      if (face_dot >= 0) then
        !
        ! xyz_pt is inside this (left) donor cell.
        ! Interpolate donor sol to this xyz_pt.
        !
#ifdef DEBUG_ON
        call debug_timer(start_accumulation, &
                         "LEFT cell interp_cell_sol_to_xyz",acc=2)
#endif
        donor_cell = face(nf)%left%cell
        !
        call interp_cell_sol_to_xyz(donor_cell,g2_xyz_pt, g2_usp(k,1:nq))
#ifdef DEBUG_ON
        call debug_timer(stop_accumulation, &
                         "LEFT cell interp_cell_sol_to_xyz",acc=2)
#endif
        !
      else if (nf > nfbnd) then
        !
        ! xyz_pt is inside adjacent (right) cell.
        ! Interpolate donor sol to this xyz_pt.
        !
#ifdef DEBUG_ON
        call debug_timer(start_accumulation, &
                         "RIGHT cell interp_cell_sol_to_xyz",acc=3)
#endif
        donor_cell = face(nf)%right%cell
        !
        call interp_cell_sol_to_xyz(donor_cell,g2_xyz_pt, g2_usp(k,1:nq))
#ifdef DEBUG_ON
        call debug_timer(stop_accumulation, &
                         "RIGHT cell interp_cell_sol_to_xyz",acc=3)
#endif
        !
      else
        !
#ifdef DEBUG_ON
        call debug_timer(start_accumulation,"outside boundary face",acc=4)
#endif
        ! xyz_pt is outside boundary face of cell.
        ! Interpolate donor volume sol to closest donor face pt.
        ! Similar to subroutine interp_cv_to_flxpts.
        !
        n0 = cell(nc)%beg_sp - 1
        np = cell(nc)%end_sp - n0 ! number of sol-pts in this cell
        !
       !n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
       !n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
       !np = n2-n1+1         ! number of sol-pts in this cell
        !
        this_geom = cell(nc)%geom   ! geometry type of this cell
        this_order = cell(nc)%order ! solution order of this cell
        !
        ! Interpolate volume solution to face point location.
        ! Assign that face point solution to g2 dst point.
        !
        face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
        face_order = cell(nc)%face(ncf)%order ! order of current face
        !
        fp_offset = cell(nc)%face(ncf)%fp_offset
        face_rot = cell(nc)%face(ncf)%rotation
        !
        ! ifp: index of the flx-pt local to this face
        ! lf: index of flx-pt within interp%mat
       !lf = cell(nc)%face(ncf)%get_fp_idx(ifp)
        !
        lf = fp_offset + fp%geord(face_geom,face_order)%rot(face_rot)%pts(ifp)
        !
#define TEMPORARY_FOR_UTRC_CASE
#ifdef TEMPORARY_FOR_UTRC_CASE
        valsp1 = zero
        valsp2 = zero
        valsp3 = zero
        valsp4 = zero
        valsp5 = zero
        valsp6 = zero
        valec  = zero
        valmx  = zero
        valmy  = zero
        valmz  = zero
        valee  = zero
        !
        do i = interp(this_geom,this_order)% &
                         toFace(face_order)%row_ptr(lf)+1, &
               interp(this_geom,this_order)% &
                         toFace(face_order)%row_ptr(lf+1)
          !
          j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
          a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
          !
          valsp1 = valsp1 + a * usp(n0+j,1)
          valsp2 = valsp2 + a * usp(n0+j,2)
          valsp3 = valsp3 + a * usp(n0+j,3)
          valsp4 = valsp4 + a * usp(n0+j,4)
          valsp5 = valsp5 + a * usp(n0+j,5)
          valsp6 = valsp6 + a * usp(n0+j,6)
          valec  = valec  + a * usp(n0+j,nec)
          valmx  = valmx  + a * usp(n0+j,nmx)
          valmy  = valmy  + a * usp(n0+j,nmy)
          valmz  = valmz  + a * usp(n0+j,nmz)
          valee  = valee  + a * usp(n0+j,nee)
          !
        end do
        !
        g2_usp(k,1)   = valsp1
        g2_usp(k,2)   = valsp2
        g2_usp(k,3)   = valsp3
        g2_usp(k,4)   = valsp4
        g2_usp(k,5)   = valsp5
        g2_usp(k,6)   = valsp6
        g2_usp(k,nec) = valec
        g2_usp(k,nmx) = valmx
        g2_usp(k,nmy) = valmy
        g2_usp(k,nmz) = valmz
        g2_usp(k,nee) = valee
#else
#ifdef ASSUME_3D
        valec  = zero
        valmx  = zero
        valmy  = zero
        valmz  = zero
        valee  = zero
        !
        do i = interp(this_geom,this_order)% &
                         toFace(face_order)%row_ptr(lf)+1, &
               interp(this_geom,this_order)% &
                         toFace(face_order)%row_ptr(lf+1)
          !
          j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
          a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
          !
          valec  = valec  + a * usp(n0+j,nec)
          valmx  = valmx  + a * usp(n0+j,nmx)
          valmy  = valmy  + a * usp(n0+j,nmy)
          valmz  = valmz  + a * usp(n0+j,nmz)
          valee  = valee  + a * usp(n0+j,nee)
          !
        end do
        !
        g2_usp(k,nec) = valec
        g2_usp(k,nmx) = valmx
        g2_usp(k,nmy) = valmy
        g2_usp(k,nmz) = valmz
        g2_usp(k,nee) = valee
#else
        do m = 1,nq
          !
          val = zero
          !
          do i = interp(this_geom,this_order)% &
                           toFace(face_order)%row_ptr(lf)+1, &
                 interp(this_geom,this_order)% &
                           toFace(face_order)%row_ptr(lf+1)
            !
            j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
            a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
            !
            val = val + a * usp(n0+j,m)
            !
          end do
          !
          g2_usp(k,m) = val
          !
        end do
#endif
#endif
        !
       !g2_usp(1:nq,k) = interp(this_geom,this_order)% &
       !                        toFace(face_order)% &
       !                        dot_mat( lf , tusp(1:np,1:nq) )
       !do m = 1,nq
       !  g2_usp(m,k) = interp(this_geom,this_order)% &
       !                       toFace(face_order)% &
       !                       dot( lf , tusp(1:np,m) )
       !end do
#ifdef DEBUG_ON
        call debug_timer(stop_accumulation,"outside boundary face",acc=4)
#endif
        !
      end if
      !
    else
      !
      write (error_message,1000)
      call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
      !
    end if
    !
    pts_completed = pts_completed + 1
    !
  end do  ! k-loop
  !$OMP END DO
  !
  call octtree_find_nearest_xyz(octtree,g2_xyz_pt, ixyz_pt,istat, &
                                tid=tid,acc=2)
  call debug_timer(report_accumulated_time, &
                   "octtree_find_nearest_xyz",acc=1)
  call debug_timer(report_accumulated_time, &
                   "LEFT cell interp_cell_sol_to_xyz",acc=2)
  call debug_timer(report_accumulated_time, &
                   "RIGHT cell interp_cell_sol_to_xyz",acc=3)
  call debug_timer(report_accumulated_time, &
                   "outside boundary face",acc=4)
  write (iout,1001) mypnum,tid
  1001 format ("CPU # ",i0," , Thread # ",i2," finished interpolation!")
  flush (iout)
  !$OMP END PARALLEL
  !
  1000 format ("ERROR:  Octree failed to locate nearest", &
               " point for interpolation.")
  !
#endif
  !
  !-DEBUG_BEG:
  if (mypnum == glb_root) then
    write (iout,*) "DEBUG_END:[interpsol]"
    flush (iout)
  end if
  !-DEBUG_END:
  !
  ! Cleanup
  !
  !-DEBUG_BEG:
  if (mypnum == glb_root) then
    write (iout,10) pname, "Deallocate octtree"
    flush (iout)
  end if
  !-DEBUG_END:
  !
  dealloc_octedat = fals         ! Do this manually instead.
  call octtree_dealloc(octtree,dealloc_octedat,istat)
  !
  !-DEBUG_BEG:
  if (mypnum == glb_root) then
    write (iout,10) pname, "Deallocate octtree now complete"
    flush (iout)
  end if
  !-DEBUG_END:
  !
  nullify(octtree%octedat)          ! Manually clear octtree%octedat
 !nullify(octedat%facexyz)
 !nullify(octedat%cell)
  octedat%varname = ""
  octedat%vardesc = ""
  !
  !-DEBUG_BEG:
  if (mypnum == glb_root) then
    write (iout,10) pname, "Leaving subroutine interpsol"
    flush (iout)
  end if
  !-DEBUG_END:
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  10 format (1x,"DEBUG[",a,"]",a,".")
  !
end subroutine interpsol
!
!###############################################################################
!
subroutine interpsol_debug()
  !
  ! Interpolate solution onto a specific xyz location.
  ! This routine is similar to interpsol and is used to debug interpolation
  ! at a specified xyz location (hard-coded below).
  !
  !.. Use Statements ..
  use module_kind_types, only : dot,lk,fals
  use geovar,            only : nnode,ncell,nface,nfbnd,nr
  use geovar,            only : cell,face,fp
  use geovar,            only : xyz
  use flux_mod,          only : interp_cv_to_flxpts  ! DEBUG: for debug output
  use metrics_mod,       only : geofa,facenrm
  use flowvar,           only : usp
 !use flowvar,           only : facexyz_dt
  use flowvar,           only : facexyz_array
  use flowvar,           only : faceusp_array
  use order_mod,         only : maxSP,geom_solpts
  use octedat_facexyz_mod
  use octtree_mod
  use octstat_mod
  use interpolation_mod, only : interp_cell_sol_to_xyz, interp
  !-DEBUG_BEG:
  use mydebug_mod
  !-DEBUG_END:
  !
  !.. Local Scalars ..
  type(octedat_facexyz_type), target :: octedat
  type(octtree_type)                 :: octtree
  !
  integer :: istat,io,ierr
  integer :: i,j,ifp,k,l,nc,ncf,nf,nfp,sf,f0
  integer :: n1,n2,np,knp,m,lf
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: fp_offset,face_rot
  integer :: donor_cell
  integer :: donor_lr_flag
  !
  logical(lk) :: dealloc_octedat
  !
  real(wp) :: a,val
  !
  !-DEBUG_BEG:
  integer :: lcell_min, lcell_max, rcell_min, rcell_max
  integer :: size_cell, size_face
  integer :: size_cell2, size_face2
  !
  logical(lk) :: debug_print
  !-DEBUG_END:
  !
  !.. Local Arrays ..
  integer  :: ixyz_pt(1:3) ! = (cell,face,facept)
  real(wp) :: xyz_pt(1:nr)
  real(wp) :: g2_xyz_pt(1:nr)
  real(wp) :: face_to_xyz(1:nr)
  real(wp) :: g2_usp_pt(1:nq)
  real(wp) :: bbox(1:nr,1:2)
  real(wp) :: tusp(1:maxSP,1:nq)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interpsol_debug"
  character(len=*), parameter :: fname = "interpsol_debug.tec"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !-DEBUG_BEG:
  ! Need faceusp for debug output.  Must compute it first.
  call interp_cv_to_flxpts()
  ! Print global array info for debugging.
  call mydebug(pname,__LINE__,__FILE__)
  !-DEBUG_END:
  !
  ! Size the bounding box to contain the dst grid volume points.
  write (iout,10) pname, "Size bounding box"  ! DEBUG
  flush (iout)
  bbox(:,1) = g2_xyz(:,1)
  bbox(:,2) = g2_xyz(:,1)
  do k = 2,size(g2_xyz,dim=2)
    do l = 1,nr
      bbox(l,1) = min( bbox(l,1), g2_xyz(l,k) )
      bbox(l,2) = max( bbox(l,2), g2_xyz(l,k) )
    end do
  end do
  !
  ! Increase the bounding box to also contain the donor grid face points.
  write (iout,10) pname, "Resize bounding box"  ! DEBUG
  flush (iout)
  do nf = 1,nface
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
   !nfp = size(facexyz_dt(nf)%v,dim=2)
    do ifp = 1,nfp
      do l = 1,nr
       !bbox(l,1) = min( bbox(l,1), facexyz_dt(nf)%v(l,ifp) )
       !bbox(l,2) = max( bbox(l,2), facexyz_dt(nf)%v(l,ifp) )
        bbox(l,1) = min( bbox(l,1), facexyz_array(l,f0+ifp) )
        bbox(l,2) = max( bbox(l,2), facexyz_array(l,f0+ifp) )
      end do
    end do
  end do
  !
  ! Init the encapsulated data for the octtree.
  write (iout,10) pname, "Init octedat"  ! DEBUG
  flush (iout)
 !octedat%varname = "facexyz_dt(1:nface)%v(1:nr,1:nfp)"
  octedat%varname = "facexyz_array(1:nr,1:n_flxpts)"
  octedat%vardesc = "Get facexyz given cell,face,facept = (nc,ncf,nfp)."
 !octedat%facexyz => facexyz
 !octedat%cell    => cell
  octedat%MAXINDICES = 3
  ! Encode octedat data using 3 indices:
  !   idx(1) = nc   ! cell no.
  !   idx(2) = ncf  ! cell face no.
  !   idx(3) = ifp  ! cell face facept no.
  !
  ! Init the empty octtree.
  write (iout,10) pname, "Init empty octtree"  ! DEBUG
  flush (iout)
  octtree%octedat => octedat
  call octtree_init(octtree,bbox,istat)
  !
  !-DEBUG_BEG:
  write (iout,10) pname, "Filling octree for use in interpolation"
  flush (iout)
  !-DEBUG_END:
  !
  ! Fill the octree with the donor grid face points.
  do nc = 1,ncell  ! cell loop
    !-DEBUG_BEG:
    !-DEBUG_CMT: Print percent complete status every 10%.
    if (nc > 1 .and. (ncell-1)*10 > 100) then
      if (mod( nc-1, ((ncell-1)*10)/100 ) == 0) then
        write (iout,201) int( 100.0*(nc-1.0)/(ncell-1.0) + 0.5 )
        flush (iout)
        201 format ("Octree fill ",i3," percent complete.")
      end if
    end if
    !-DEBUG_END:
    do ncf = 1,size(cell(nc)%face)  ! cell face loop
      nf = cell(nc)%face(ncf)%idx   ! global face index
      sf = cell(nc)%face(ncf)%side
      ! Only left faces to avoid duplicate pts from right face.
      if (sf == 1) then
        face_geom  = cell(nc)%face(ncf)%geom
        face_order = cell(nc)%face(ncf)%order
        nfp = geom_solpts(face_geom,face_order)  ! current cell face pts
        do ifp = 1,nfp  ! current cell face point loop
          !-DEBUG_BEG:
#ifdef DEBUG_OCTFILL
          write (iout,400) pname,nc,ncell, ncf,size(cell(nc)%face),ifp,nfp
          400 format (1x,"DEBUG[",a,"]  nc,ncf,ifp = ",3(i4,'/',i4,:,2x))
#endif
          ! Load ixyz_pt(:) with the index values the octtree will track.
          ixyz_pt(1) = nc   ! cell no.
          ixyz_pt(2) = ncf  ! cell face no.
          ixyz_pt(3) = ifp  ! cell face facept no.
          f0 = face(nf)%beg_fp - 1
          xyz_pt = facexyz_array(:,f0+ifp)
         !xyz_pt = facexyz_dt(nf)%v(:,ifp)
          call octtree_addxyz(octtree,ixyz_pt,xyz_pt,istat)
          if (istat /= OCT_STATUS_OK) then
            if (mypnum == glb_root) then
              write (iout,*) "ERROR:  Failed to init octree for interpolation."
              write (iout,*) "          status = ",istat
              call stop_gfr(abort,pname,__LINE__,__FILE__)
            end if
            ! Have all the host roots call mpi_barrier here so that the
            ! MPI processes that arent the global root do not continue
            ! past this point before the global root has a chance to
            ! call mpi_abort within the above conditional.
            call mpi_barrier(host_roots_comm,mpierr)
          end if
        end do  ! cell face point loop
      end if  ! Left face
    end do  ! cell face loop
  end do  ! cell loop
  !
  !
  ! For a specified xyz location,
  !   Determine donor cell.
  !   Print xyz/sol for volume points in donor cell.
  !   Interpolate donor solution onto specified xyz location.
  !
  !
  !-DEBUG_BEG:
  write (iout,10) pname, "Octtree fill complete"
  !
  !-DEBUG_CMT: Write additional cell/face info for both grids.
  write (iout,*)
  write (iout,10) pname, "Performing interpolation"
  write (iout,*) "DEBUG_BEG:[interpsol_debug]"
  size_cell=0; if (allocated(cell)) size_cell=size(cell)
  size_face=0; if (allocated(face)) size_face=size(face)
  size_cell2=0; if (allocated(g2_cell)) size_cell2=size(g2_cell)
  size_face2=0; if (allocated(g2_face)) size_face2=size(g2_face)
  write (iout,*) "      nnode             = ",nnode
  write (iout,*) "   g2_nnode             = ",g2_nnode
  write (iout,*) "      ncell/nface/nfbnd = ",ncell," / ",nface," / ",nfbnd
  write (iout,*) "   g2_ncell/nface/nfbnd = ",g2_ncell," / ",g2_nface, &
                                                     " / ",g2_nfbnd
  write (iout,*) "        size(cell/face) = ",size_cell," / ",size_face
  write (iout,*) "     g2_size(cell/face) = ",size_cell2," / ",size_face2
  !
  lcell_min = huge(0)
  rcell_min = huge(0)
  lcell_max = 0
  rcell_max = 0
  do nc = 1,ncell
    do ncf = 1,size(cell(nc)%face)
      nf = cell(nc)%face(ncf)%idx
      lcell_min = min(lcell_min, face(nf)%left%cell)
      rcell_min = min(rcell_min, face(nf)%right%cell)
      lcell_max = max(lcell_max, face(nf)%left%cell)
      rcell_max = max(rcell_max, face(nf)%right%cell)
    end do
  end do
  write (iout,*) "       lcell_min/max    = ",lcell_min," / ",lcell_max
  write (iout,*) "       rcell_min/max    = ",rcell_min," / ",rcell_max
  !
  lcell_min = huge(0)
  rcell_min = huge(0)
  lcell_max = 0
  rcell_max = 0
  do nc = 1,g2_ncell
    do ncf = 1,size(g2_cell(nc)%face)
      nf = g2_cell(nc)%face(ncf)%idx
      lcell_min = min(lcell_min, g2_face(nf)%left%cell)
      rcell_min = min(rcell_min, g2_face(nf)%right%cell)
      lcell_max = max(lcell_max, g2_face(nf)%left%cell)
      rcell_max = max(rcell_max, g2_face(nf)%right%cell)
    end do
  end do
  write (iout,*) "    g2_lcell_min/max    = ",lcell_min," / ",lcell_max
  write (iout,*) "    g2_rcell_min/max    = ",rcell_min," / ",rcell_max
  write (iout,*) "    g2_n_solpts         = ",g2_n_solpts
  write (iout,*) " Looping over g2_n_solpts in destination grid."
  !
  flush (iout)
  !-DEBUG_END:
  !
  !
  ! Loop over test points.
  do k = 1,1
    !
    ! Specify volume point xyz to interp sol onto.
    !
    if (nr == 3) then
      g2_xyz_pt(1) = 4.25_wp
      g2_xyz_pt(2) = 19.1107340438399191_wp
      g2_xyz_pt(3) = 0.567418334509571576_wp
    else
      g2_xyz_pt = zero
    end if
    !
    ! Find index of nearest face pt in src grid.
    call octtree_find_nearest_xyz(octtree,g2_xyz_pt, ixyz_pt,istat)
    !
    !-DEBUG_BEG:
    ! debug_print determines which k points debug info is printed for.
    debug_print = true
    if (debug_print) then
      write (iout,*) "  Destination point:"
      write (iout,21) "xyz = ",g2_xyz_pt(:)
      write (iout,*) "  octtree search status = ",istat
      flush (iout)
    end if
    !-DEBUG_END:
    !
    if (istat == OCT_STATUS_OK) then
      ! Decode cell info from ixyz_pt.
      nc  = ixyz_pt(1)  ! cell no.
      ncf = ixyz_pt(2)  ! cell face no.
      ifp = ixyz_pt(3)  ! cell face face_pt no.
      nf  = cell(nc)%face(ncf)%idx   ! global face index
      f0  = face(nf)%beg_fp - 1
      xyz_pt(:) = facexyz_array(:,f0+ifp)
     !xyz_pt(:) = facexyz_dt(nf)%v(:,ifp)
      face_to_xyz(:) = g2_xyz_pt(:) - xyz_pt(:)   ! vector from face to dst pt
      !
      !-DEBUG_BEG:
      !-DEBUG_CMT: Print info about nearest donor cell/face.
      if (debug_print) then
        write (iout,20) "nc         = ",nc
        write (iout,20) "ncf        = ",ncf
        write (iout,20) "ifp        = ",ifp
        write (iout,20) "nf         = ",nf
        write (iout,20) "r-cell     = ",face(nf)%right%cell
        write (iout,20) "l-cell     = ",face(nf)%left%cell
        write (iout,21) "g2_xyz_pt  = ",g2_xyz_pt(:)
        write (iout,21) "facexyz_pt = ",facexyz_array(:,f0+ifp)
       !write (iout,21) "facexyz_pt = ",facexyz_dt(nf)%v(:,ifp)
       !write (iout,21) "normal     = ",geofa(nf)%v(:,ifp)
        write (iout,21) "normal     = ",facenrm(:,f0+ifp)
        write (iout,21) "face2dst   = ",face_to_xyz(:)
       !write (iout,21) "dot(n,f2d) = ",dot(geofa(nf)%v(:,ifp), face_to_xyz(:))
        write (iout,21) "dot(n,f2d) = ",dot(facenrm(:,f0+ifp), face_to_xyz(:))
        ! Output xyz and usp for closest face point on donor grid.
        write (iout,20) "Donor grid/solution info."
       !write (iout,21) "facexyz_dt(:)          = ",facexyz_dt(nf)%v(1:nr,ifp)
        write (iout,21) "facexyz_array(:)       = ",facexyz_array(1:nr,f0+ifp)
       !write (iout,21) "faceusp_dt(:,left)     = ",faceusp_dt(nf)%v(1:nq,ifp,1)
        write (iout,21) "faceusp_array(:,left)  = ",faceusp_array(f0+ifp,1:nq,1)
       !write (iout,21) "faceusp_dt(:,right)    = ",faceusp_dt(nf)%v(1:nq,ifp,2)
        write (iout,21) "faceusp_array(:,right) = ",faceusp_array(f0+ifp,1:nq,2)
        flush (iout)
      end if
      !-DEBUG_END:
      !
      ! Dot face normal with face_to_xyz vector.
      ! Face normal always points toward right cell.
     !if (dot( geofa(nf)%v(:,ifp), face_to_xyz(:) ) <= 0) then
      if (dot( facenrm(:,f0+ifp), face_to_xyz(:) ) <= 0) then
        ! xyz_pt is inside this (left) donor cell.
        ! Interpolate donor sol to this xyz_pt.
        donor_lr_flag = 1
        !-DEBUG_BEG:
        if (debug_print) then
          write (iout,20) "Using left cell as donor."
          flush (iout)
        end if
        !-DEBUG_END:
        donor_cell = face(nf)%left%cell
        call interp_cell_sol_to_xyz(donor_cell,g2_xyz_pt,g2_usp_pt)
      else if (nf > nfbnd) then
        ! xyz_pt is inside adjacent (right) cell.
        ! Interpolate donor sol to this xyz_pt.
        donor_lr_flag = 2
        donor_cell = face(nf)%right%cell
        call interp_cell_sol_to_xyz(donor_cell,g2_xyz_pt,g2_usp_pt)
        !-DEBUG_BEG:
        if (debug_print) then
          write (iout,20) "Using right cell as donor."
          flush (iout)
        end if
        !-DEBUG_END:
      else
        ! xyz_pt is outside boundary face of cell.
        ! Interpolate donor volume sol to closest donor face pt.
        ! Similar to subroutine interp_cv_to_flxpts.
        !
        donor_lr_flag = 3
        !-DEBUG_BEG:
        if (debug_print) then
          write (iout,20) "Using closest donor face pt."
          flush (iout)
        end if
        !-DEBUG_END:
        donor_cell = nc
        n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
        n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
        np = n2-n1+1         ! number of sol-pts in this cell
        !
        this_geom = cell(nc)%geom   ! geometry type of this cell
        this_order = cell(nc)%order ! solution order of this cell
        !
        ! Create transpose of usp array to reduce cache misses.
        !
        do knp = 1,np
          do m = 1,nq
            tusp(knp,m) = usp(knp+n1-1,m)
          end do
        end do
        !
        ! Interpolate volume solution to face point location.
        ! Assign that face point solution to g2 dst point.
        !
        face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
        face_order = cell(nc)%face(ncf)%order ! order of current face
        !
        fp_offset = cell(nc)%face(ncf)%fp_offset
        face_rot = cell(nc)%face(ncf)%rotation
        !
        ! ifp: index of the flx-pt local to this face
        ! lf: index of flx-pt within interp%mat
       !lf = cell(nc)%face(ncf)%get_fp_idx(ifp)
        !
        lf = fp_offset + fp%geord(face_geom,face_order)%rot(face_rot)%pts(ifp)
        !
        do m = 1,nq
          !
          val = zero
          !
          do i = interp(this_geom,this_order)% &
                           toFace(face_order)%row_ptr(lf)+1, &
                 interp(this_geom,this_order)% &
                           toFace(face_order)%row_ptr(lf+1)
            !
            j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
            a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
            !
            val = val + a * tusp(j,m)
            !
          end do
          !
          g2_usp_pt(m) = val
          !
        end do
        !
       !g2_usp_pt(1:nq) = interp(this_geom,this_order)% &
       !                         toFace(face_order)% &
       !                         dot_mat( lf , tusp(1:np,1:nq) )
       !do m = 1,nq
       !  g2_usp_pt(m) = interp(this_geom,this_order)% &
       !                        toFace(face_order)% &
       !                        dot( lf , tusp(1:np,m) )
       !end do
        !
      end if
      !
      ! Write xyz/usp for all volume points in donor cell.
      open (newunit=io,file=fname,status="replace",action="write", &
                       iostat=ierr,iomsg=error_message)
      call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
      !
      nc = donor_cell
      n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
      n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
      np = n2-n1+1         ! number of sol-pts in this cell
      knp = int((real(np))**one3+half)
      write (io,90)
      write (io,91) g2_xyz_pt(:)
      write (io,92)
      write (io,93) knp,knp,knp
      do knp = 1,np
        write (io,94) xyz(1:nr,knp+n1-1), usp(knp+n1-1,1:nq)
      end do
      !
      close (io,iostat=ierr,iomsg=error_message)
      call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
      !
    else
      ! istat /= OCT_STATUS_OK
      if (mypnum == glb_root) then
        write (error_message,80)
        call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
      end if
      ! Have all the host roots call mpi_barrier here so that the
      ! MPI processes that arent the global root do not continue
      ! past this point before the global root has a chance to
      ! call mpi_abort within the above conditional.
      call mpi_barrier(host_roots_comm,mpierr)
    end if
    !
    !-DEBUG_BEG:
    !-DEBUG_CMT: Output xyz and usp for point on dest grid.
    if (debug_print) then
      write (iout,20) "Interpolated grid/solution info."
      write (iout,21) "   xyz(:)          = ",g2_xyz_pt(:)
      write (iout,21) "   usp(:)          = ",g2_usp_pt(:)
      write (iout,*)
      flush (iout)
    end if
    !-DEBUG_END:
    !
  end do  ! k-loop
  !
  !-DEBUG_BEG:
  write (iout,*) "DEBUG_END:[interpsol_debug]"
  flush (iout)
  !-DEBUG_END:
  !
  ! Cleanup
  !
  !-DEBUG_BEG:
  write (iout,10) pname, "Deallocate octtree"
  !-DEBUG_END:
  dealloc_octedat = fals         ! Do this manually instead.
  call octtree_dealloc(octtree,dealloc_octedat,istat)
  !-DEBUG_BEG:
  write (iout,10) pname, "Deallocate octtree now complete"
  !-DEBUG_END:
  nullify(octtree%octedat)          ! Manually clear octtree%octedat
 !nullify(octedat%facexyz)
 !nullify(octedat%cell)
  octedat%varname = ""
  octedat%vardesc = ""
  !
  !-DEBUG_BEG:
  write (iout,10) pname, "Leaving subroutine interpsol_debug"
  !-DEBUG_END:
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  10 format (1x,"DEBUG[",a,"]",a,".")
  !
  20 format (5x,a,:,i0)
  21 format (5x,a,:,10(e20.6,:,/,26x))
  !
  80 format ("ERROR: OCTREE FAILED TO LOCATE NEAREST POINT FOR INTERPOLATION.")
  !
  90 format ("# Solution info for the donor cell containing point:")
  91 format ("#   xyz = ",3e20.10)
  92 format ('Variables="x","y","z","usp1","usp2","usp3","usp4","usp5"')
  93 format ('Zone T="Donor Cell", I=',i3," J=",i3," K=",i3,' F="POINT"')
  94 format (2x,99(e20.10))
  !
end subroutine interpsol_debug
!
!###############################################################################
!
subroutine deallocate_solver_arrays()
  !
  ! Deallocate the solver arrays that are not needed, in order to have
  ! enough memory for interpolation.
  !
  !.. Use Statements ..
  use geovar,      only : face,grid
  use flowvar,     only : faceflx
 !use flowvar,     only : faceusp_dt,facedusp_dt,facexyz_dt
  use flowvar,     only : faceusp_array,facedusp_array,facexyz_array
  use flowvar,     only : residual,dtsp,uoldsp,rk_sum,rssprk
  use flowvar,     only : nodeusp_dt,nodedusp_dt
  use flowvar,     only : edgeusp_dt,edgedusp_dt,edgexyz_dt
  use flowvar,     only : global_usp,global_var
  use metrics_mod, only : metrics_dt,geofa
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "deallocate_solver_arrays"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (allocated(grid)) then
    deallocate ( grid , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(faceusp_array)) then
    deallocate ( faceusp_array , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"faceusp_array",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(facedusp_array)) then
    deallocate ( facedusp_array , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"facedusp_array",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
 !!
 !if (allocated(faceusp_dt)) then
 !  deallocate ( faceusp_dt , stat=ierr , errmsg=error_message )
 !  call alloc_error(pname,"faceusp_dt",2,__LINE__,__FILE__,ierr,error_message)
 !end if
 !!
 !if (allocated(facedusp_dt)) then
 !  deallocate ( facedusp_dt , stat=ierr , errmsg=error_message )
 !  call alloc_error(pname,"facedusp_dt",2,__LINE__,__FILE__,ierr,error_message)
 !end if
  !
  if (allocated(faceflx)) then
    deallocate ( faceflx , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"faceflx",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(facexyz_array)) then
    deallocate ( facexyz_array , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"facexyz_array",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
 !!
 !if (allocated(facexyz_dt)) then
 !  deallocate ( facexyz_dt , stat=ierr , errmsg=error_message )
 !  call alloc_error(pname,"facexyz_dt",2,__LINE__,__FILE__,ierr,error_message)
 !end if
  !
  if (allocated(residual)) then
    deallocate ( residual , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"residual",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(dtsp)) then
    deallocate ( dtsp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"dtsp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(uoldsp)) then
    deallocate ( uoldsp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"uoldsp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(rk_sum)) then
    deallocate ( rk_sum , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"rk_sum",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(rssprk)) then
    deallocate ( rssprk , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"rssprk",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(nodeusp_dt)) then
    deallocate ( nodeusp_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nodeusp_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(nodedusp_dt)) then
    deallocate ( nodedusp_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nodedusp_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(edgeusp_dt)) then
    deallocate ( edgeusp_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edgeusp_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(edgedusp_dt)) then
    deallocate ( edgedusp_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edgedusp_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(edgexyz_dt)) then
    deallocate ( edgexyz_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edgexyz_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(global_usp)) then
    deallocate ( global_usp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"global_usp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(global_var)) then
    deallocate ( global_var , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"global_var",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(metrics_dt)) then
    deallocate ( metrics_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"metrics_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(geofa)) then
    deallocate ( geofa , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"geofa",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine deallocate_solver_arrays
!
!###############################################################################
!
#include "Functions/debug_timer.f90"
!
!###############################################################################
!
end module interpsol_mod
