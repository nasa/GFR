module interpsol_mod

! Module to handle interpolation of (grid1,sol1) -> (grid2,sol2).

  use module_kind_types
! use module_kind_types, only : wp,lk,true,fals,error_message, &
!                               Geom_Edge,Geom_Tria,Geom_Quad,Geom_Tetr, &
!                               Geom_Pyra,Geom_Pris,Geom_Hexa, &
!                               dot
  use interpsolvar_mod
  implicit none

  private
  public mv_global_to_grid2
  public mv_grid2_to_global
  public interpsol
  public interpsol_debug

contains

!##########################################################################
!
subroutine mv_global_to_grid2()

! Move the global variables to grid2, leaving globals deallocated.

  use geovar
  use flowvar
  use ovar,              only : bc_names
  implicit none

!-DEBUG_BEG:
  write(iout,*) 'DEBUG[mv_global_to_grid2] Move global vars to grid2.'
  flush(iout)
!-DEBUG_END:

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
  g2_n_global_node    = n_global_node
  g2_n_global_totnod  = n_global_totnod
  g2_n_global_cell    = n_global_cell
  g2_n_global_totcel  = n_global_totcel
  g2_n_global_edge    = n_global_edge
  g2_n_global_fbnd    = n_global_fbnd
  g2_n_global_solpts  = n_global_solpts
  g2_n_global_totpts  = n_global_totpts
  call move_alloc( from=cell_geom         , to=g2_cell_geom         )
  call move_alloc( from=cell_order        , to=g2_cell_order        )
  call move_alloc( from=xyz_nodes         , to=g2_xyz_nodes         )
  call move_alloc( from=xyz               , to=g2_xyz               )
  call move_alloc( from=global_xyz_nodes  , to=g2_global_xyz_nodes  )
  call move_alloc( from=global_xyz        , to=g2_global_xyz        )
  call move_alloc( from=global_nofc_ptr   , to=g2_global_nofc_ptr   )
  call move_alloc( from=global_nofc       , to=g2_global_nofc       )
  call move_alloc( from=global_pinc_ptr   , to=g2_global_pinc_ptr   )
  call move_alloc( from=global_cell_geom  , to=g2_global_cell_geom  )
  call move_alloc( from=global_cell_order , to=g2_global_cell_order )
  call move_alloc( from=bface             , to=g2_bface             )
  call move_alloc( from=bc_conflict       , to=g2_bc_conflict       )
  call move_alloc( from=nodes_of_cell_ptr , to=g2_nodes_of_cell_ptr )
  call move_alloc( from=nodes_of_cell     , to=g2_nodes_of_cell     )
  call move_alloc( from=nodes_of_cell     , to=g2_nodes_of_cell     )
  call move_alloc( from=pst_geom          , to=g2_pst_geom          )
  call move_alloc( from=fp                , to=g2_fp                )
  call move_alloc( from=face              , to=g2_face              )
  call move_alloc( from=cell              , to=g2_cell              )
  call move_alloc( from=grid              , to=g2_grid              )
  call move_alloc( from=cell_map          , to=g2_cell_map          )
  
! Variables from flowvar.
  call move_alloc( from=usp               , to=g2_usp               )
  call move_alloc( from=dusp              , to=g2_dusp              )
  call move_alloc( from=facexyz           , to=g2_facexyz           )
  
! Variables from ovar.
  call move_alloc( from=bc_names          , to=g2_bc_names          )
  
!-DEBUG_BEG:
  write(iout,*) 'DEBUG[mv_global_to_grid2] Leaving subroutine.'
  flush(iout)
!-DEBUG_END:

  return
end subroutine mv_global_to_grid2
!
!##########################################################################
!
subroutine mv_grid2_to_global()

! Move the grid2 variables to global, leaving grid2 deallocated.

  use geovar
  use flowvar
  use ovar,              only : bc_names
  implicit none

!-DEBUG_BEG:
  write(iout,*) 'DEBUG[mv_grid2_to_global] Move grid2 vars to global.'
  flush(iout)
!-DEBUG_END:

! Variables from geovar.
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
  n_global_node    = g2_n_global_node
  n_global_totnod  = g2_n_global_totnod
  n_global_cell    = g2_n_global_cell
  n_global_totcel  = g2_n_global_totcel
  n_global_edge    = g2_n_global_edge
  n_global_fbnd    = g2_n_global_fbnd
  n_global_solpts  = g2_n_global_solpts
  n_global_totpts  = g2_n_global_totpts
  call move_alloc( from=g2_cell_geom         , to=cell_geom         )
  call move_alloc( from=g2_cell_order        , to=cell_order        )
  call move_alloc( from=g2_xyz_nodes         , to=xyz_nodes         )
  call move_alloc( from=g2_xyz               , to=xyz               )
  call move_alloc( from=g2_global_xyz_nodes  , to=global_xyz_nodes  )
  call move_alloc( from=g2_global_xyz        , to=global_xyz        )
  call move_alloc( from=g2_global_nofc_ptr   , to=global_nofc_ptr   )
  call move_alloc( from=g2_global_nofc       , to=global_nofc       )
  call move_alloc( from=g2_global_pinc_ptr   , to=global_pinc_ptr   )
  call move_alloc( from=g2_global_cell_geom  , to=global_cell_geom  )
  call move_alloc( from=g2_global_cell_order , to=global_cell_order )
  call move_alloc( from=g2_bface             , to=bface             )
  call move_alloc( from=g2_bc_conflict       , to=bc_conflict       )
  call move_alloc( from=g2_nodes_of_cell_ptr , to=nodes_of_cell_ptr )
  call move_alloc( from=g2_nodes_of_cell     , to=nodes_of_cell     )
  call move_alloc( from=g2_nodes_of_cell     , to=nodes_of_cell     )
  call move_alloc( from=g2_pst_geom          , to=pst_geom          )
  call move_alloc( from=g2_fp                , to=fp                )
  call move_alloc( from=g2_face              , to=face              )
  call move_alloc( from=g2_cell              , to=cell              )
  call move_alloc( from=g2_grid              , to=grid              )
  call move_alloc( from=g2_cell_map          , to=cell_map          )

!-DEBUG_BEG:
  write(iout,*) 'DEBUG_BEG[mv_grid2_to_global]'
  write(iout,*) '  Before move_alloc(from=g2_usp,to=usp):'
  write(iout,2001) '    g2_usp(:,1) = ',g2_usp(:,1)
  write(iout,2001) '       usp(:,1) = ',   usp(:,1)
  2001 format(x,a,e20.6,:,/,(20x,e20.6))
  flush(iout)
!-DEBUG_END:

! Variables from flowvar.
  call move_alloc( from=g2_usp               , to=usp               )
  call move_alloc( from=g2_dusp              , to=dusp              )
  call move_alloc( from=g2_facexyz           , to=facexyz           )
  
!-DEBUG_BEG:
  write(iout,*) '  After  move_alloc(from=g2_usp,to=usp):'
  write(iout,*)    '    g2_usp(:,1) = unallocated'
  write(iout,2001) '       usp(:,1) = ',   usp(:,1)
  write(iout,*) 'DEBUG_END[mv_grid2_to_global]'
  flush(iout)
!-DEBUG_END:

! Variables from ovar.
  call move_alloc( from=g2_bc_names          , to=bc_names          )
  
!-DEBUG_BEG:
  write(iout,*) 'DEBUG[mv_global_to_grid2] Leaving subroutine.'
  flush(iout)
!-DEBUG_END:

  return
end subroutine mv_grid2_to_global
!
!##########################################################################
!
subroutine interpsol()
  
! Interpolate solution onto grid2.

  use module_kind_types, only : dot
  use eqn_idx, only : nq
  use geovar,            only : ncell,nface,nfbnd,nr
  use geovar,            only : cell,face
  use flux_mod,          only : interp_cv_to_flxpts  ! DEBUG: for debug output
  use metrics_mod,       only : geofa
  use flowvar,           only : facexyz, usp
  use order_mod,         only : maxSP,geom_solpts
  use octedat_facexyz_mod
  use octtree_mod
  use interpolation_mod, only : interp_cell_sol_to_xyz, interp
!-DEBUG_BEG:
! use ifport           , only : sleep
  use mydebug_mod
!-DEBUG_END:
  implicit none

  !.. Local Scalars ..
  type(octedat_facexyz_type) ,target :: octedat
  type(octtree_type)                 :: octtree
  integer :: status
  integer :: ifp,k,l,nc,ncf,nf,nfp,sf
  integer :: n1,n2,np,knp,m,lf
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: donor_cell
  logical :: dealloc_octedat
!-DEBUG_BEG:
  integer :: lcell_min, lcell_max, rcell_min, rcell_max
  integer :: size_cell, size_face
  integer :: size_cell2, size_face2
  logical :: debug_print
!-DEBUG_END:

  !.. Local Arrays ..
  integer  ,dimension(3)    :: ixyz_pt ! = (cell,face,facept)
  real(wp) ,dimension(nr)   :: xyz_pt, g2_xyz_pt,face_to_xyz
  real(wp) ,dimension(nr,2) :: bbox
  real(wp), dimension(1:maxSP,1:nq) :: tusp
  !
 
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interpsol"

continue

!-DEBUG_BEG:
  ! Need faceusp for debug output.  Must compute it first.
  call interp_cv_to_flxpts()
  ! Print global array info for debugging.
  call mydebug(pname,__LINE__,__FILE__)
!-DEBUG_END:

! Size the bounding box to contain the dst grid volume points.
  write(iout,*) 'DEBUG[interpsol] Size bounding box.'  ! DEBUG
  flush(iout)
  bbox(:,1) = g2_xyz(:,1)
  bbox(:,2) = g2_xyz(:,1)
  do k = 2,size(g2_xyz,dim=2)
    do l = 1,nr
      bbox(l,1) = min( bbox(l,1), g2_xyz(l,k) )
      bbox(l,2) = max( bbox(l,2), g2_xyz(l,k) )
    end do
  end do

! Increase the bounding box to also contain the donor grid face points.
  write(iout,*) 'DEBUG[interpsol] Resize bounding box.'  ! DEBUG
  flush(iout)
  do nf = 1,nface
    do ifp = 1,size(facexyz(nf)%v,dim=2)
      do l = 1,nr
        bbox(l,1) = min( bbox(l,1), facexyz(nf)%v(l,ifp) )
        bbox(l,2) = max( bbox(l,2), facexyz(nf)%v(l,ifp) )
      end do
    end do
  end do

! Init the encapsulated data for the octtree.
  write(iout,*) 'DEBUG[interpsol] Init octedat.'  ! DEBUG
  flush(iout)
  octedat%varname = "facexyz(1:nface)%v(1:nr,1:nfp)"
  octedat%vardesc = "Get facexyz given cell,face,facept = (nc,ncf,nfp)."
! octedat%facexyz => facexyz
! octedat%cell    => cell
  octedat%MAXINDICES = 3
  ! Encode octedat data using 3 indices:
  !   idx(1) = nc   ! cell no.
  !   idx(2) = ncf  ! cell face no.
  !   idx(3) = ifp  ! cell face facept no.

! Init the empty octtree.
  write(iout,*) 'DEBUG[interpsol] Init empty octtree.'  ! DEBUG
  flush(iout)
  octtree%octedat => octedat
  call octtree_init(octtree,bbox,status)

!-DEBUG_BEG:
  write(iout,*) 'DEBUG[interpsol] Filling octree for use in interpolation.'
  flush(iout)
!-DEBUG_END:

! Fill the octree with the donor grid face points.
  do nc = 1,ncell  ! cell loop
!-DEBUG_BEG:
!-DEBUG_CMT: Print percent complete status every 10%.
    if ( nc .gt. 1 .and. (ncell-1)*10 > 100 ) then
      if ( mod( nc-1, ((ncell-1)*10)/100 ) == 0 ) then
        write(iout,201) int( 100.0*(nc-1.0)/(ncell-1.0) + 0.5 )
        flush(iout)
        201 format('Octree fill ',i3,' percent complete.')
      end if
    end if
!-DEBUG_END:
    do ncf = 1,size(cell(nc)%face)  ! cell face loop
      nf = cell(nc)%face(ncf)%idx   ! global face index
      sf = cell(nc)%face(ncf)%side
      ! Only left faces to avoid duplicate pts from right face.
      if ( sf .eq. 1 ) then
        face_geom  = cell(nc)%face(ncf)%geom
        face_order = cell(nc)%face(ncf)%order
        nfp = geom_solpts(face_geom,face_order)  ! current cell face pts
        do ifp = 1,nfp  ! current cell face point loop
!-DEBUG_BEG:
#ifdef DEBUG_OCTFILL
    write(iout,'(x,3a, i4,"/",i4,"  ", i4,"/",i4,"  ", i4,"/",i4)') &
      'DEBUG[',pname,']  nc,ncf,ifp = ', &
      nc,ncell, ncf,size(cell(nc)%face) ,ifp,nfp
#endif
!-DEBUG_BEG:
!   call sleep(1)  ! Sleep # seconds.
!-DEBUG_END:
          ! Load ixyz_pt(:) with the index values the octtree will track.
          ixyz_pt(1) = nc   ! cell no.
          ixyz_pt(2) = ncf  ! cell face no.
          ixyz_pt(3) = ifp  ! cell face facept no.
          xyz_pt = facexyz(nf)%v(:,ifp)
          call octtree_addxyz(octtree,ixyz_pt,xyz_pt,status)
          if ( status .ne. OCT_STATUS_OK ) then
            if ( mypnum == glb_root) then
              write(iout,*) 'ERROR:  Failed to init octree for interpolation.'
              write(iout,*) '          status = ',status
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
  write(iout,*) 'DEBUG[interpsol] Octtree fill complete.'
  write(iout,*) 
  write(iout,*) 'DEBUG[interpsol] Performing interpolation.'
  write(iout,*) 'DEBUG_BEG:[interpsol]'
  size_cell=0; if ( allocated(cell) ) size_cell=size(cell)
  size_face=0; if ( allocated(face) ) size_face=size(face)
  size_cell2=0; if ( allocated(g2_cell) ) size_cell2=size(g2_cell)
  size_face2=0; if ( allocated(g2_face) ) size_face2=size(g2_face)
  write(iout,*) '      ncell/nface/nfbnd = ',ncell,' / ',nface,' / ',nfbnd
  write(iout,*) '   g2_ncell/nface/nfbnd = ',g2_ncell,' / ',g2_nface, &
                                                     ' / ',g2_nfbnd
  write(iout,*) '        size(cell/face) = ',size_cell ,' / ',size_face
  write(iout,*) '     g2_size(cell/face) = ',size_cell2,' / ',size_face2
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
    enddo
  enddo
  write(iout,*) '       lcell_min/max    = ',lcell_min,' / ',lcell_max
  write(iout,*) '       rcell_min/max    = ',rcell_min,' / ',rcell_max
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
    enddo
  enddo
  write(iout,*) '    g2_lcell_min/max    = ',lcell_min,' / ',lcell_max
  write(iout,*) '    g2_rcell_min/max    = ',rcell_min,' / ',rcell_max
  write(iout,*) '    g2_n_solpts         = ',g2_n_solpts
  write(iout,*) ' Looping over g2_n_solpts in destination grid.'
  !
  flush(iout)
!-DEBUG_END:
! Loop over volume points in dst grid.
  do k = 1,g2_n_solpts
!-DEBUG_BEG:
!-DEBUG_CMT: Print percent complete status every 10%.
    if ( k .gt. 1 .and. (g2_n_solpts-1)*10 > 100 ) then
      if ( mod( k-1, ((g2_n_solpts-1)*10)/100 ) == 0 ) then
        write(iout,301) int( 100.0*(k-1.0)/(g2_n_solpts-1.0) + 0.5 )
        flush(iout)
        301 format('Interpolation ',i3,' percent complete.')
      end if
    end if
!-DEBUG_END:
    ! Get volume point xyz in dst grid to interp sol onto.
    g2_xyz_pt = g2_xyz(:,k)
    ! Find index of nearest face pt in src grid.
    call octtree_find_nearest_xyz(octtree,g2_xyz_pt, ixyz_pt,status)
!-DEBUG_BEG:
! debug_print determines which k points debug info is printed for.
  debug_print = .false.
! if ( k <= 1+9 ) debug_print = .true.
! if ( any( [1,10,100,500,1000,10000,100000,300000] == k ) ) &
!    debug_print = .true.
! if ( k >= g2_n_solpts-9 ) debug_print = .true.
  !
  if ( debug_print ) then
    write(iout,*) '  Destination point = k = ',k,' / ',g2_n_solpts
    write(iout,*) '  octtree search status = ',status
    flush(iout)
  endif
!-DEBUG_END:
    if ( status .eq. OCT_STATUS_OK ) then
      nc  = ixyz_pt(1)  ! cell no.
      ncf = ixyz_pt(2)  ! cell face no.
      ifp = ixyz_pt(3)  ! cell face face_pt no.
      nf  = cell(nc)%face(ncf)%idx   ! global face index
      xyz_pt(:) = facexyz(nf)%v(:,ifp)
      face_to_xyz(:) = g2_xyz_pt(:) - xyz_pt(:)   ! vector from face to dst pt
!-DEBUG_BEG:
if ( debug_print ) then
  write(iout,*) '    nc         = ',nc
  write(iout,*) '    ncf        = ',ncf
  write(iout,*) '    ifp        = ',ifp
  write(iout,*) '    nf         = ',nf
  write(iout,*) '    r-cell     = ',face(nf)%right%cell
  write(iout,*) '    l-cell     = ',face(nf)%left%cell
  write(iout,*) '    g2_xyz_pt  = ',g2_xyz_pt(:)
  write(iout,*) '    facexyz_pt = ',facexyz(nf)%v(:,ifp)
  write(iout,*) '    normal     = ',geofa(nf)%v(:,ifp)
  write(iout,*) '    face2dst   = ',face_to_xyz(:)
  write(iout,*) '    dot(n,f2d) = ',dot( geofa(nf)%v(:,ifp), face_to_xyz(:) )
  if ( allocated(g2_usp) ) then
  write(iout,*) '    g2_usp     = allocated'
  else
  write(iout,*) '    g2_usp     = NOT allocated'
  endif
  flush(iout)
  write(iout,*) '    size(g2_usp,dim=1) = ',size(g2_usp,dim=1)
  write(iout,*) '    size(g2_usp,dim=2) = ',size(g2_usp,dim=2)
  ! Output xyz and usp for closest face point on donor grid.
  2001 format(x,a,e20.6,:,/,(26x,e20.6))
  write(iout,*)    '    Donor grid/solution info.'
  write(iout,2001) '    facexyz(:)         = ',facexyz(nf)%v(1:nr,ifp)
  write(iout,2001) '    faceusp(:,left)    = ',faceusp(nf)%v(1:nq,ifp,1)
  write(iout,2001) '    faceusp(:,right)   = ',faceusp(nf)%v(1:nq,ifp,2)
  flush(iout)
endif
!-DEBUG_END:
      ! Dot face normal with face_to_xyz vector.
      ! Face normal always points toward right cell.
      if ( dot( geofa(nf)%v(:,ifp), face_to_xyz(:) ) .le. 0 ) then
        ! xyz_pt is inside this (left) donor cell.
        ! Interpolate donor sol to this xyz_pt.
!-DEBUG_BEG:
if ( debug_print ) then
  write(iout,*) '    Using left cell as donor.'
  flush(iout)
endif
!-DEBUG_END:
        donor_cell = face(nf)%left%cell
        call interp_cell_sol_to_xyz &
               (donor_cell,g2_xyz_pt, g2_usp(:,k))
      else if ( nf > nfbnd ) then
        ! xyz_pt is inside adjacent (right) cell.
        ! Interpolate donor sol to this xyz_pt.
        donor_cell = face(nf)%right%cell
        call interp_cell_sol_to_xyz &
               (donor_cell,g2_xyz_pt, g2_usp(:,k))
!-DEBUG_BEG:
if ( debug_print ) then
  write(iout,*) '    Using right cell as donor.'
  flush(iout)
endif
!-DEBUG_END:
      else
        ! xyz_pt is outside boundary face of cell.
        ! Interpolate donor volume sol to closest donor face pt.
        ! Similar to subroutine interp_cv_to_flxpts.
        !
!-DEBUG_BEG:
if ( debug_print ) then
  write(iout,*) '    Using closest donor face pt.'
  flush(iout)
endif
!-DEBUG_END:
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
            tusp(knp,m) = usp(m,knp+n1-1)
          end do
        end do
        !
        ! Interpolate volume solution to face point location.
        ! Assign that face point solution to g2 dst point.
        !
        face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
        face_order = cell(nc)%face(ncf)%order ! order of current face
        !
        ! ifp: index of the flx-pt local to this face
        ! lf: index of flx-pt within interp%mat
        lf = cell(nc)%face(ncf)%get_fp_idx(ifp)
        !
        g2_usp(1:nq,k) = interp(this_geom,this_order)% &
                                toFace(face_order)% &
                                dot_mat( lf , tusp(1:np,1:nq) )
       !do m = 1,nq
       !  g2_usp(m,k) = interp(this_geom,this_order)% &
       !                       toFace(face_order)% &
       !                       dot( lf , tusp(1:np,m) )
       !end do
        !
      end if
    else
      if ( mypnum == glb_root) then
        write(iout,*) 'ERROR:  Octree failed to locate nearest point ', &
                      'for interpolation.'
        call stop_gfr(abort,pname,__LINE__,__FILE__)
      end if
      ! Have all the host roots call mpi_barrier here so that the
      ! MPI processes that arent the global root do not continue
      ! past this point before the global root has a chance to
      ! call mpi_abort within the above conditional.
      call mpi_barrier(host_roots_comm,mpierr)
    end if
!-DEBUG_BEG:
  ! Output xyz and usp for point on dest grid.
if ( debug_print ) then
  write(iout,*)    '    Destination grid/solution info.'
  write(iout,2001) '    g2_xyz(:)          = ',g2_xyz_pt(:)
  write(iout,2001) '    g2_usp(:)          = ',g2_usp(:,k)
  write(iout,*)
  flush(iout)
  if ( .false. .and. k .eq. 10 ) then
    error_message = '[DEBUG] Only print first 10 points and stop.'
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  endif
endif
!-DEBUG_END:
  end do  ! k-loop
  
!-DEBUG_BEG:
  write(iout,*) 'DEBUG_END:[interpsol]'
  flush(iout)
!-DEBUG_END:
!
! Cleanup
!
!-DEBUG_BEG:
  write(iout,*) 'DEBUG[interpsol] Deallocate octtree.'
!-DEBUG_END:
  dealloc_octedat = .false.         ! Do this manually instead.
  call octtree_dealloc(octtree,dealloc_octedat,status)
!-DEBUG_BEG:
  write(iout,*) 'DEBUG[interpsol] Deallocate octtree now complete.'
!-DEBUG_END:
  nullify(octtree%octedat)          ! Manually clear octtree%octedat
! nullify(octedat%facexyz)
! nullify(octedat%cell)
  octedat%varname = ""
  octedat%vardesc = ""

!-DEBUG_BEG:
  write(iout,*) 'DEBUG[interpsol] Leaving subroutine interpsol.'
!-DEBUG_END:
  return
end subroutine interpsol
!
!##########################################################################
!
subroutine interpsol_debug()
  
! Interpolate solution onto a specific xyz location.
! This routine is similar to interpsol and is used to debug interpolation
! at a specified xyz location (hard-coded below).

  use module_kind_types, only : dot
  use eqn_idx, only : nq
  use geovar,            only : nnode,ncell,nface,nfbnd,nr
  use geovar,            only : cell,face
  use geovar,            only : xyz
  use flux_mod,          only : interp_cv_to_flxpts  ! DEBUG: for debug output
  use metrics_mod,       only : geofa
  use flowvar,           only : facexyz, usp
  use order_mod,         only : maxSP,geom_solpts
  use octedat_facexyz_mod
  use octtree_mod
  use interpolation_mod, only : interp_cell_sol_to_xyz, interp
!-DEBUG_BEG:
  use mydebug_mod
!-DEBUG_END:
  implicit none

  !.. Local Scalars ..
  type(octedat_facexyz_type) ,target :: octedat
  type(octtree_type)                 :: octtree
  integer :: status
  integer :: iotec
  integer :: ifp,k,l,nc,ncf,nf,nfp,sf
  integer :: n1,n2,np,knp,m,lf
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: donor_cell
  integer :: donor_lr_flag
  logical :: dealloc_octedat
!-DEBUG_BEG:
  integer :: lcell_min, lcell_max, rcell_min, rcell_max
  integer :: size_cell, size_face
  integer :: size_cell2, size_face2
  logical :: debug_print
!-DEBUG_END:

  !.. Local Arrays ..
  integer  ,dimension(3)    :: ixyz_pt ! = (cell,face,facept)
  real(wp) ,dimension(nr)   :: xyz_pt, g2_xyz_pt,face_to_xyz
  real(wp) ,dimension(nq)   ::         g2_usp_pt
  real(wp) ,dimension(nr,2) :: bbox
  real(wp), dimension(1:maxSP,1:nq) :: tusp
  !
 
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interpsol_debug"

continue

!-DEBUG_BEG:
  ! Need faceusp for debug output.  Must compute it first.
  call interp_cv_to_flxpts()
  ! Print global array info for debugging.
  call mydebug(pname,__LINE__,__FILE__)
!-DEBUG_END:

! Size the bounding box to contain the dst grid volume points.
  write(iout,*) 'DEBUG[interpsol_debug] Size bounding box.'  ! DEBUG
  flush(iout)
  bbox(:,1) = g2_xyz(:,1)
  bbox(:,2) = g2_xyz(:,1)
  do k = 2,size(g2_xyz,dim=2)
    do l = 1,nr
      bbox(l,1) = min( bbox(l,1), g2_xyz(l,k) )
      bbox(l,2) = max( bbox(l,2), g2_xyz(l,k) )
    end do
  end do

! Increase the bounding box to also contain the donor grid face points.
  write(iout,*) 'DEBUG[interpsol_debug] Resize bounding box.'  ! DEBUG
  flush(iout)
  do nf = 1,nface
    do ifp = 1,size(facexyz(nf)%v,dim=2)
      do l = 1,nr
        bbox(l,1) = min( bbox(l,1), facexyz(nf)%v(l,ifp) )
        bbox(l,2) = max( bbox(l,2), facexyz(nf)%v(l,ifp) )
      end do
    end do
  end do

! Init the encapsulated data for the octtree.
  write(iout,*) 'DEBUG[interpsol_debug] Init octedat.'  ! DEBUG
  flush(iout)
  octedat%varname = "facexyz(1:nface)%v(1:nr,1:nfp)"
  octedat%vardesc = "Get facexyz given cell,face,facept = (nc,ncf,nfp)."
! octedat%facexyz => facexyz
! octedat%cell    => cell
  octedat%MAXINDICES = 3
  ! Encode octedat data using 3 indices:
  !   idx(1) = nc   ! cell no.
  !   idx(2) = ncf  ! cell face no.
  !   idx(3) = ifp  ! cell face facept no.

! Init the empty octtree.
  write(iout,*) 'DEBUG[interpsol_debug] Init empty octtree.'  ! DEBUG
  flush(iout)
  octtree%octedat => octedat
  call octtree_init(octtree,bbox,status)

!-DEBUG_BEG:
  write(iout,*) 'DEBUG[interpsol_debug] Filling octree for use in interpolation.'
  flush(iout)
!-DEBUG_END:

! Fill the octree with the donor grid face points.
  do nc = 1,ncell  ! cell loop
!-DEBUG_BEG:
!-DEBUG_CMT: Print percent complete status every 10%.
    if ( nc .gt. 1 .and. (ncell-1)*10 > 100 ) then
      if ( mod( nc-1, ((ncell-1)*10)/100 ) == 0 ) then
        write(iout,201) int( 100.0*(nc-1.0)/(ncell-1.0) + 0.5 )
        flush(iout)
        201 format('Octree fill ',i3,' percent complete.')
      end if
    end if
!-DEBUG_END:
    do ncf = 1,size(cell(nc)%face)  ! cell face loop
      nf = cell(nc)%face(ncf)%idx   ! global face index
      sf = cell(nc)%face(ncf)%side
      ! Only left faces to avoid duplicate pts from right face.
      if ( sf .eq. 1 ) then
        face_geom  = cell(nc)%face(ncf)%geom
        face_order = cell(nc)%face(ncf)%order
        nfp = geom_solpts(face_geom,face_order)  ! current cell face pts
        do ifp = 1,nfp  ! current cell face point loop
!-DEBUG_BEG:
#ifdef DEBUG_OCTFILL
    write(iout,'(x,3a, i4,"/",i4,"  ", i4,"/",i4,"  ", i4,"/",i4)') &
      'DEBUG[',pname,']  nc,ncf,ifp = ', &
      nc,ncell, ncf,size(cell(nc)%face) ,ifp,nfp
#endif
          ! Load ixyz_pt(:) with the index values the octtree will track.
          ixyz_pt(1) = nc   ! cell no.
          ixyz_pt(2) = ncf  ! cell face no.
          ixyz_pt(3) = ifp  ! cell face facept no.
          xyz_pt = facexyz(nf)%v(:,ifp)
          call octtree_addxyz(octtree,ixyz_pt,xyz_pt,status)
          if ( status .ne. OCT_STATUS_OK ) then
            if ( mypnum == glb_root) then
              write(iout,*) 'ERROR:  Failed to init octree for interpolation.'
              write(iout,*) '          status = ',status
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
! For a specified xyz location,
!   Determine donor cell.
!   Print xyz/sol for volume points in donor cell.
!   Interpolate donor solution onto specified xyz location.
!
!
!-DEBUG_BEG:
  write(iout,*) 'DEBUG[interpsol_debug] Octtree fill complete.'
!
!-DEBUG_CMT: Write additional cell/face info for both grids.
  write(iout,*) 
  write(iout,*) 'DEBUG[interpsol_debug] Performing interpolation.'
  write(iout,*) 'DEBUG_BEG:[interpsol_debug]'
  size_cell=0; if ( allocated(cell) ) size_cell=size(cell)
  size_face=0; if ( allocated(face) ) size_face=size(face)
  size_cell2=0; if ( allocated(g2_cell) ) size_cell2=size(g2_cell)
  size_face2=0; if ( allocated(g2_face) ) size_face2=size(g2_face)
  write(iout,*) '      nnode             = ',nnode
  write(iout,*) '   g2_nnode             = ',g2_nnode
  write(iout,*) '      ncell/nface/nfbnd = ',ncell,' / ',nface,' / ',nfbnd
  write(iout,*) '   g2_ncell/nface/nfbnd = ',g2_ncell,' / ',g2_nface, &
                                                     ' / ',g2_nfbnd
  write(iout,*) '        size(cell/face) = ',size_cell ,' / ',size_face
  write(iout,*) '     g2_size(cell/face) = ',size_cell2,' / ',size_face2
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
    enddo
  enddo
  write(iout,*) '       lcell_min/max    = ',lcell_min,' / ',lcell_max
  write(iout,*) '       rcell_min/max    = ',rcell_min,' / ',rcell_max
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
    enddo
  enddo
  write(iout,*) '    g2_lcell_min/max    = ',lcell_min,' / ',lcell_max
  write(iout,*) '    g2_rcell_min/max    = ',rcell_min,' / ',rcell_max
  write(iout,*) '    g2_n_solpts         = ',g2_n_solpts
  write(iout,*) ' Looping over g2_n_solpts in destination grid.'
  !
  flush(iout)
!-DEBUG_END:
!
!
! Loop over test points.
  do k = 1,1
    ! Specify volume point xyz to interp sol onto.
    if ( nr .eq. 3 ) then
      g2_xyz_pt(1) = 4.25
      g2_xyz_pt(2) = 19.1107340438399191
      g2_xyz_pt(3) = 0.567418334509571576
    else
      g2_xyz_pt = 0.0
    endif
    !
    ! Find index of nearest face pt in src grid.
    call octtree_find_nearest_xyz(octtree,g2_xyz_pt, ixyz_pt,status)
!
!-DEBUG_BEG:
! debug_print determines which k points debug info is printed for.
  debug_print = .true.
  if ( debug_print ) then
    write(iout,*) '  Destination point:'
    write(iout,*) '    xyz = ',g2_xyz_pt(:)
    write(iout,*) '  octtree search status = ',status
    flush(iout)
  endif
!-DEBUG_END:
!
    if ( status .eq. OCT_STATUS_OK ) then
      ! Decode cell info from ixyz_pt.
      nc  = ixyz_pt(1)  ! cell no.
      ncf = ixyz_pt(2)  ! cell face no.
      ifp = ixyz_pt(3)  ! cell face face_pt no.
      nf  = cell(nc)%face(ncf)%idx   ! global face index
      xyz_pt(:) = facexyz(nf)%v(:,ifp)
      face_to_xyz(:) = g2_xyz_pt(:) - xyz_pt(:)   ! vector from face to dst pt
!
!-DEBUG_BEG:
!-DEBUG_CMT: Print info about nearest donor cell/face.
if ( debug_print ) then
  write(iout,*) '    nc         = ',nc
  write(iout,*) '    ncf        = ',ncf
  write(iout,*) '    ifp        = ',ifp
  write(iout,*) '    nf         = ',nf
  write(iout,*) '    r-cell     = ',face(nf)%right%cell
  write(iout,*) '    l-cell     = ',face(nf)%left%cell
  write(iout,*) '    g2_xyz_pt  = ',g2_xyz_pt(:)
  write(iout,*) '    facexyz_pt = ',facexyz(nf)%v(:,ifp)
  write(iout,*) '    normal     = ',geofa(nf)%v(:,ifp)
  write(iout,*) '    face2dst   = ',face_to_xyz(:)
  write(iout,*) '    dot(n,f2d) = ',dot( geofa(nf)%v(:,ifp), face_to_xyz(:) )
  ! Output xyz and usp for closest face point on donor grid.
  2001 format(x,a,e20.6,:,/,(26x,e20.6))
  write(iout,*)    '    Donor grid/solution info.'
  write(iout,2001) '    facexyz(:)         = ',facexyz(nf)%v(1:nr,ifp)
  write(iout,2001) '    faceusp(:,left)    = ',faceusp(nf)%v(1:nq,ifp,1)
  write(iout,2001) '    faceusp(:,right)   = ',faceusp(nf)%v(1:nq,ifp,2)
  flush(iout)
endif
!-DEBUG_END:
!
      ! Dot face normal with face_to_xyz vector.
      ! Face normal always points toward right cell.
      if ( dot( geofa(nf)%v(:,ifp), face_to_xyz(:) ) .le. 0 ) then
        ! xyz_pt is inside this (left) donor cell.
        ! Interpolate donor sol to this xyz_pt.
        donor_lr_flag = 1
!-DEBUG_BEG:
if ( debug_print ) then
  write(iout,*) '    Using left cell as donor.'
  flush(iout)
endif
!-DEBUG_END:
        donor_cell = face(nf)%left%cell
        call interp_cell_sol_to_xyz &
               (donor_cell,g2_xyz_pt, g2_usp_pt  )
      else if ( nf > nfbnd ) then
        ! xyz_pt is inside adjacent (right) cell.
        ! Interpolate donor sol to this xyz_pt.
        donor_lr_flag = 2
        donor_cell = face(nf)%right%cell
        call interp_cell_sol_to_xyz &
               (donor_cell,g2_xyz_pt, g2_usp_pt  )
!-DEBUG_BEG:
if ( debug_print ) then
  write(iout,*) '    Using right cell as donor.'
  flush(iout)
endif
!-DEBUG_END:
      else
        ! xyz_pt is outside boundary face of cell.
        ! Interpolate donor volume sol to closest donor face pt.
        ! Similar to subroutine interp_cv_to_flxpts.
        !
        donor_lr_flag = 3
!-DEBUG_BEG:
if ( debug_print ) then
  write(iout,*) '    Using closest donor face pt.'
  flush(iout)
endif
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
            tusp(knp,m) = usp(m,knp+n1-1)
          end do
        end do
        !
        ! Interpolate volume solution to face point location.
        ! Assign that face point solution to g2 dst point.
        !
        face_geom = cell(nc)%face(ncf)%geom   ! geom  of current face
        face_order = cell(nc)%face(ncf)%order ! order of current face
        !
        ! ifp: index of the flx-pt local to this face
        ! lf: index of flx-pt within interp%mat
        lf = cell(nc)%face(ncf)%get_fp_idx(ifp)
        !
        g2_usp_pt(1:nq) = interp(this_geom,this_order)% &
                                 toFace(face_order)% &
                                 dot_mat( lf , tusp(1:np,1:nq) )
       !do m = 1,nq
       !  g2_usp_pt(m) = interp(this_geom,this_order)% &
       !                        toFace(face_order)% &
       !                        dot( lf , tusp(1:np,m) )
       !end do
        !
      end if
      !
      ! Write xyz/usp for all volume points in donor cell.
      open(newunit=iotec, file='interpsol_debug.tec', &
           form='formatted', status='unknown')
        nc = donor_cell
        n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
        n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
        np = n2-n1+1         ! number of sol-pts in this cell
        knp = int((real(np))**(1.0D+00/3.0D+00)+0.5)
        write(iotec,2601)
        write(iotec,2602) g2_xyz_pt(:)
        write(iotec,2611)
        write(iotec,2612) knp,knp,knp
 2601   format('# Solution info for the donor cell containing point:')
 2602   format('#   xyz = ',3e20.10)
 2611   format('Variables="x","y","z","usp1","usp2","usp3","usp4","usp5"')
 2612   format('Zone T="Donor Cell", I=',i3,' J=',i3,' K=',i3,' F="POINT"')
        do knp = 1,np
          write(iotec,2621) xyz(1:nr,knp+n1-1), &
                            usp(1:nq,knp+n1-1)
 2621     format(2x,99(e20.10))
        end do
      close(iotec)
      !
    else
      ! status .ne. OCT_STATUS_OK
      if ( mypnum == glb_root) then
        write(iout,*) 'ERROR:  Octree failed to locate nearest point ', &
                      'for interpolation.'
        call stop_gfr(abort,pname,__LINE__,__FILE__)
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
if ( debug_print ) then
  write(iout,*)    '    Interpolated grid/solution info.'
  write(iout,2001) '       xyz(:)          = ',g2_xyz_pt(:)
  write(iout,2001) '       usp(:)          = ',g2_usp_pt(:)
  write(iout,*)
  flush(iout)
endif
!-DEBUG_END:
!
  end do  ! k-loop
  
!-DEBUG_BEG:
  write(iout,*) 'DEBUG_END:[interpsol_debug]'
  flush(iout)
!-DEBUG_END:
!
! Cleanup
!
!-DEBUG_BEG:
  write(iout,*) 'DEBUG[interpsol_debug] Deallocate octtree.'
!-DEBUG_END:
  dealloc_octedat = .false.         ! Do this manually instead.
  call octtree_dealloc(octtree,dealloc_octedat,status)
!-DEBUG_BEG:
  write(iout,*) 'DEBUG[interpsol_debug] Deallocate octtree now complete.'
!-DEBUG_END:
  nullify(octtree%octedat)          ! Manually clear octtree%octedat
! nullify(octedat%facexyz)
! nullify(octedat%cell)
  octedat%varname = ""
  octedat%vardesc = ""

!-DEBUG_BEG:
  write(iout,*) 'DEBUG[interpsol_debug] Leaving subroutine interpsol_debug.'
!-DEBUG_END:
  return
end subroutine interpsol_debug
!###############################################################################
!
subroutine deallocate_solver_arrays()
  !
  ! Deallocate the solver arrays that are not needed, in order to have
  ! enough memory for interpolation.
  !
  !.. Use Statements ..
  use geovar,      only : face,grid
  use flowvar,     only : faceusp,facedusp,faceflx,facexyz
  use flowvar,     only : residual,dtsp,uoldsp,usp_sum,rssprk
  use flowvar,     only : nodeusp,nodedusp,edgeusp,edgedusp,edgexyz
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
  if (allocated(faceusp)) then
    deallocate ( faceusp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"faceusp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(facedusp)) then
    deallocate ( facedusp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"facedusp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(faceflx)) then
    deallocate ( faceflx , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"faceflx",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(facexyz)) then
    deallocate ( facexyz , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"facexyz",2,__LINE__,__FILE__,ierr,error_message)
  end if
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
  if (allocated(usp_sum)) then
    deallocate ( usp_sum , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"usp_sum",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(rssprk)) then
    deallocate ( rssprk , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"rssprk",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(nodeusp)) then
    deallocate ( nodeusp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nodeusp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(nodedusp)) then
    deallocate ( nodedusp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nodedusp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(edgeusp)) then
    deallocate ( edgeusp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edgeusp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(edgedusp)) then
    deallocate ( edgedusp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edgedusp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(edgexyz)) then
    deallocate ( edgexyz , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edgexyz",2,__LINE__,__FILE__,ierr,error_message)
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
!##########################################################################
end module interpsol_mod
