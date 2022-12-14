module periodic_mod
  !
  !.. Use Statements ..
  use module_kind_types
  !
  implicit none
  !
  private
  !
  ! ##################################
  ! ###  Module Public Procedures  ###
  ! ##################################
  !
  public :: match_periodic_faces
  public :: check_periodic_list_indices
  !
  ! ##############################
  ! ###  Module Derived Types  ###
  ! ##############################
  !
  type :: bnd_face_t
    integer :: bc_idx
    integer :: geom
    integer, allocatable :: nodes(:)
  end type bnd_face_t
  !
  type :: periodic_face_t
    integer :: nf
    integer :: bc_idx
    integer :: rotation
  end type periodic_face_t
  !
  type(periodic_face_t), save, allocatable :: periodic_faces(:)
  !
contains
!
!###############################################################################
!
subroutine match_periodic_faces()
  !
  !.. Use Statements ..
  use geovar, only : nr,nfbnd,nnode,ncell
  use geovar, only : bface,xyz_nodes
  use ovar,   only : bc_input
  !
  use module_cgns_types, only : CGLEN
  !
  !.. Local Scalars ..
  integer :: l,n,nf,nc,npbc,npair,npf,npn,nfnods,this_pair,ierr
  integer :: i,ip,iidx,ipf,igeom,match_iidx,num_iidx_cells
  integer :: j,jp,jidx,jpf,jgeom,match_jidx,num_jidx_cells
  !
  logical(lk) :: error_found
  !
  character(len=CGLEN) :: match_jnam,jname,iname,match_inam
  !
  !.. Local Arrays ..
  real(wp) :: pbc_rotation_center(1:nr)
  real(wp) :: pbc_rotation_angle(1:nr)
  real(wp) :: pbc_translation(1:nr)
  real(wp) :: xyz_ip(1:nr)
  !
  integer :: ipn(1:4)
  integer :: jpn(1:4)
  integer :: mpn(1:4)
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: pbc_idxs(:)
  integer, allocatable :: pairs(:,:)
  integer, allocatable :: paired(:)
  integer, allocatable :: ifaces(:)
  integer, allocatable :: jfaces(:)
  integer, allocatable :: inodes(:)
  integer, allocatable :: jnodes(:)
  integer, allocatable :: node_match(:)
  integer, allocatable :: face_match(:)
  integer, allocatable :: marker(:)
  !
  real(wp), allocatable :: ds(:)
  !
  logical(lk), allocatable :: msk(:)
  logical(lk), allocatable :: matched(:)
  !
  character(len=CGLEN), allocatable :: pbc_names(:)
  !
  type(bnd_face_t), allocatable :: bnd_face(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "match_periodic_faces"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Return if bc_input isnt allocated
  !
  if (.not. allocated(bc_input)) then
    call debug_timer(leaving_procedure,pname)
    return
  end if
  !
  ! Return if no periodic boundaries were specified in bc_input
  !
  if (count(bc_input(:)%set_bc_value == bc_periodic) == 0) then
    call debug_timer(leaving_procedure,pname)
    return
  end if
  !
  ! Collect the boundary face element connectivity into bnd_face
  ! using nodes_of_cell_ptr and nodes_of_cell. This is needed because
  ! only the host_roots have created the grid derived type, and this
  ! will quickly create the subset of grid that we need here, and all
  ! mpi processes will be able to do this successfully. If we were to
  ! use grid to do this, only the host roots would be able to do this
  ! and then we would have to mpi_bcast the updates to bface to all
  ! non-root processes.
  !
  call get_boundary_face_nodes(bnd_face)
  !
  ! Create an array that gives the indices of only the elements of
  ! bc_input that have been set to a periodic boundary condition
  !
  pbc_idxs = pack( intseq(1,size(bc_input)) , & ! F2003 AUTO-REALLOC
                   bc_input(:)%set_bc_value == bc_periodic )
  !
  npbc = size(pbc_idxs) ! number of periodic boundary conditions
  npair = npbc/2
  !
  ! Make sure that the number of periodic boundaries is even
  ! because they should be specified in matching pairs
  !
  if (mod(npbc,2) /= 0) then
    if (mypnum == glb_root) then
      write (iout,100)
      write (iout,101) npbc
      write (iout,100)
    end if
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  end if
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! Loop through the periodic boundary conditions and check the values
  ! specified for finding the periodic match.
  !
  do i = 1,npbc
    !
    iidx = pbc_idxs(i)
    !
    match_jidx = bc_input(iidx)%periodic%matching_bc_idx
    match_jnam = trim(adjustl(bc_input(iidx)%periodic%matching_bc_name))
    !
    ! Check the values specified for finding the periodic match, and
    ! complete the matching information if needed
    !
    if (match_jidx /= 0 .and. len_trim(match_jnam) > 0) then
      !
      ! periodic%matching_bc_idx and periodic%matching_bc_name have both
      ! been specified. Make sure this information is consistent. If it
      ! isnt treat matching_bc_idx as the good information
      !
      jname = trim(adjustl(bc_input(match_jidx)%name))
      !
      if (match_jnam /= jname) then
        !
        if (mypnum == glb_root) write (iout,102)
        !
        bc_input(iidx)%periodic%matching_bc_name = bc_input(match_jidx)%name
        !
      end if
      !
    else if (match_jidx /= 0 .and. len_trim(match_jnam) == 0) then
      !
      ! Only periodic%matching_bc_idx has been specified, so fill in
      ! the information for periodic%matching_bc_name
      !
      bc_input(iidx)%periodic%matching_bc_name = bc_input(match_jidx)%name
      !
    else if (match_jidx == 0 .and. len_trim(match_jnam) > 0) then
      !
      ! Only periodic%matching_bc_name has been specified, so fill in
      ! the information for periodic%matching_bc_idx
      !
      find_matching_jidx: &
      do j = 1,npbc
        !
        if (i == j) cycle find_matching_jidx
        !
        jidx = pbc_idxs(j)
        !
        jname = trim(adjustl(bc_input(jidx)%name))
        !
        if (jname == match_jnam) then
          bc_input(iidx)%periodic%matching_bc_idx = jidx
          exit find_matching_jidx
        end if
        !
      end do find_matching_jidx
      !
      if (mypnum == glb_root) then
        write (iout,100)
        write (iout,103)
        write (iout,100)
      end if
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
      !
    else
      !
      ! Neither periodic%matching_bc_name or periodic%matching_bc_idx has
      ! been specified. This can only be the case if there are only two
      ! periodic boundary conditions. If that is the case, set the information
      ! to the other periodic boundary. If there are more than two periodic
      ! boundaries, exit with an error.
      !
      if (npbc == 2) then
        !
        if (i == 1) then
          jidx = pbc_idxs(2)
        else
          jidx = pbc_idxs(1)
        end if
        !
        bc_input(iidx)%periodic%matching_bc_name = bc_input(jidx)%name
        bc_input(iidx)%periodic%matching_bc_idx  = jidx
        !
      else
        !
        if (mypnum == glb_root) then
          write (iout,100)
          write (iout,104)
          write (iout,100)
        end if
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
        !
      end if
      !
    end if
    !
  end do
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! At this point, periodic%matching_bc_name and periodic%matching_bc_idx
  ! should both be consistent for each periodic boundary condition, so we now
  ! need to go through and find matches.
  !
  ! Create an array of the boundary condition names
  !
  pbc_names = bc_input(pbc_idxs(:))%name ! F2003 AUTO-REALLOC
  !
  pairs = init_array(2,npair,0)   ! F2003 AUTO-REALLOC
  paired = init_array(npbc,0)     ! F2003 AUTO-REALLOC
  matched = init_array(npbc,fals) ! F2003 AUTO-REALLOC
  !
  n = 1
  outer_pair_loop: &
  do i = 1,npbc
    !
    if (paired(i) > 0) cycle outer_pair_loop
    !
    iidx = pbc_idxs(i)
    !
    iname = trim(adjustl(pbc_names(i)))
    !
    match_jidx = bc_input(iidx)%periodic%matching_bc_idx
    match_jnam = trim(adjustl(bc_input(iidx)%periodic%matching_bc_name))
    !
    ! Use match_jidx to find the match
    !
    jname = trim(adjustl(bc_input(match_jidx)%name))
    !
    match_iidx = bc_input(match_jidx)%periodic%matching_bc_idx
    match_inam = trim(adjustl(bc_input(match_jidx)%periodic%matching_bc_name))
    !
    if (iname == match_inam .and. iidx == match_iidx) then
      !
      do l = 1,npbc
        if (pbc_idxs(l) == match_jidx) then
          j = l
          exit
        end if
      end do
      !
      pairs(1:2,n) = [i,j]
      paired([i,j]) = [j,i]
      matched([i,j]) = true
      !
      n = n + 1
      !
      cycle outer_pair_loop
      !
    else
      !
      if (mypnum == glb_root) write (iout,105)
      !
    end if
    !
  end do outer_pair_loop
  !
  ! Make sure the expected number of matched pairs were found
  !
  if (n-1 /= npair) then
    if (mypnum == glb_root) then
      write (iout,100)
      write (iout,106) n-1,npair
      write (iout,100)
    end if
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  end if
  !
  ! Make sure all periodic boundary groups found a matching pair
  !
  if (.not. all(paired > 0)) then
    if (mypnum == glb_root) then
      write (iout,100)
      write (iout,107)
      do n = 1,npbc
        iidx = pbc_idxs(n)
        jidx = pbc_idxs(paired(n))
        if (jidx > 0) then
          write (iout,108) n,trim(adjustl(bc_input(iidx)%name)), &
                             trim(adjustl(bc_input(jidx)%name))
        else
          write (iout,109) n,trim(adjustl(bc_input(iidx)%name))
        end if
      end do
      write (iout,100)
    end if
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  end if
  !
  ! Make sure that the number of boundary faces match between two matching
  ! periodic boundary groups
  !
  error_found = fals
  do n = 1,npair
    iidx = pbc_idxs( pairs(1,n) )
    jidx = pbc_idxs( pairs(2,n) )
    num_iidx_cells = count(bnd_face(:)%bc_idx == iidx)
    num_jidx_cells = count(bnd_face(:)%bc_idx == jidx)
    if (num_iidx_cells /= num_jidx_cells) then
      if (.not. error_found) then
        if (mypnum == glb_root) then
          write (iout,100)
          write (iout,110)
        end if
        error_found = true
      end if
      if (mypnum == glb_root) then
        write (iout,111) n,trim(adjustl(bc_input(iidx)%name)), &
                           num_iidx_cells, &
                           trim(adjustl(bc_input(jidx)%name)), &
                           num_jidx_cells
      end if
    end if
  end do
  !
  if (error_found) then
    if (mypnum == glb_root) write (iout,100)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  end if
  error_found = fals
  !
  ! We no longer need pbc_names or paired so deallocate
  !
  if (allocated(pbc_names)) then
    deallocate ( pbc_names , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"pbc_names",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(paired)) then
    deallocate ( paired , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"paired",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! Go through each periodic pair and find the matching boundary faces an
  ! save the necessary matching information into bface
  !
  match_pair_faces: &
  do this_pair = 1,npair
    !
    iidx = pbc_idxs( pairs(1,this_pair) )
    jidx = pbc_idxs( pairs(2,this_pair) )
    !
    pbc_rotation_center(1:nr) = bc_input(iidx)%periodic%RotationCenter(1:nr)
    pbc_rotation_angle(1:nr) = bc_input(iidx)%periodic%RotationAngle(1:nr)
    pbc_translation(1:nr) = bc_input(iidx)%periodic%Translation(1:nr)
    !
    ! F2003 AUTO-REALLOC
    ifaces = pack( intseq(1,nfbnd) , bnd_face(:)%bc_idx==iidx )
    ! F2003 AUTO-REALLOC
    jfaces = pack( intseq(1,nfbnd) , bnd_face(:)%bc_idx==jidx )
    !
    npf = size(ifaces) ! number of periodic faces
    !
    ! Find all the grid nodes for the i-faces
    !
    msk = init_array(nnode,fals) ! F2003 AUTO-REALLOC
    !
    do nf = 1,npf
      msk( bnd_face( ifaces(nf) )%nodes(:) ) = true
    end do
    inodes = pack( intseq(1,nnode) , msk ) ! F2003 AUTO-REALLOC
    !
    ! Find all the grid nodes for the j-faces
    !
    msk(:) = fals
    do nf = 1,npf
      msk( bnd_face( jfaces(nf) )%nodes(:) ) = true
    end do
    jnodes = pack( intseq(1,nnode) , msk ) ! F2003 AUTO-REALLOC
    !
    npn = size(inodes) ! number of periodic nodes
    !
    ! Loop through the i-nodes and find each matching j-node
    !
    node_match = init_array(nnode,0) ! F2003 AUTO-REALLOC
    !
    matched = init_array(npn,fals) ! F2003 AUTO-REALLOC
    ds = init_array(npn,zero) ! F2003 AUTO-REALLOC
    !
    inode_loop: &
    do i = 1,npn
      !
      ip = inodes(i)
      xyz_ip(1:nr) = periodic_transformation(xyz_nodes(1:nr,ip), &
                                             pbc_rotation_center, &
                                             pbc_rotation_angle, &
                                             pbc_translation)
      !
      ds(:) = huge(zero) ! reset all to a large number so we can find the min
      jnode_loop: &
      do j = 1,npn
        !
        ! Skip to the next j-node if this one has already been matched
        if (matched(j)) cycle jnode_loop
        !
        ds(j) = norm2( xyz_ip(1:nr) - xyz_nodes(1:nr,jnodes(j)) )
        !
      end do jnode_loop
      !
      j = minloc( ds , dim=1 )
      matched(j) = true
      !
      ! Mark jp negative if ds(j) is > eps6
      !
      jp = jnodes(j)
      if (ds(j) > eps6) then
        jp = -jp
      end if
      ! 
      node_match(ip) = jp
      !
    end do inode_loop
    !
    ! Now used the matched nodes to find the matching faces
    !
    face_match = init_array(nfbnd,0) ! F2003 AUTO-REALLOC
    !
    matched = init_array(npf,fals) ! F2003 AUTO-REALLOC
    marker = init_array(nnode,0) ! F2003 AUTO-REALLOC
    !
    iface_loop: &
    do i = 1,npf
      !
      ipf = ifaces(i) ! current i-face index from 1-nfbnd
      !
      igeom = bnd_face(ipf)%geom ! geometry of i-face
      nfnods = size(bnd_face(ipf)%nodes(:)) ! number of nodes for i-face
      !
      ! Get the nodes for the current i-face
      ipn(1:nfnods) = bnd_face(ipf)%nodes(:) ! nodes defining i-face
      !
      ! Use node_match to covert the i-nodes to their matching j-nodes
      mpn(1:nfnods) = node_match( ipn(1:nfnods) )
      !
      ! Mark all the matching j-nodes with the value 1.
      ! All other nodes are set to 0.
      marker( mpn(1:nfnods) ) = 1
      !
      ! Loop over all the j-faces until we find the j-face whose nodes
      ! put into marker sum to nfnods
      !
      jface_loop: &
      do j = 1,npf
        !
        ! Skip to next j-face if this one has already been matched
        if (matched(j)) cycle jface_loop
        !
        jpf = jfaces(j) ! current j-face index from 1-nfbnd
        !
        jgeom = bnd_face(jpf)%geom ! geometry of j-face
        !
        ! Skip to next j-face if the j-face geometry is not the same as i-face
        if (igeom /= jgeom) cycle jface_loop
        !
        ! Get the nodes for the current j-face
        jpn(1:nfnods) = bnd_face(jpf)%nodes(:) ! nodes defining j-face
        !
        ! Check the sum of marker using only the j-face nodes. If the sum is
        ! equal to nfnods, this j-face is the match for the current i-face.
        if (sum( marker( jpn(1:nfnods) ) ) == nfnods) then
          face_match(ipf) = jpf
          matched(j) = true
          ! NOTE: We need to exit the j-face loop so we can reset all of
          !       elements of marker to 0. We want to do that after the
          !       j-face loop just in case a match wasnt found.
          exit jface_loop
        end if
        !
      end do jface_loop
      !
      marker( mpn(1:nfnods) ) = 0
      !
    end do iface_loop
    !
    ! Now that we have the matching periodic boundary faces, save this
    ! information to bface
    !
    do i = 1,npf
      !
      ipf = ifaces(i)
      jpf = face_match(ipf)
      !
      bface(5,ipf) = jpf
      bface(6,ipf) = bface(2,jpf)
      bface(8,ipf) = bface(3,jpf)
      bface(9,ipf) = bface(4,jpf)
      !
      bface(5,jpf) = ipf
      bface(6,jpf) = bface(2,ipf)
      bface(8,jpf) = bface(3,ipf)
      bface(9,jpf) = bface(4,ipf)
      !
    end do
    !
  end do match_pair_faces
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! Kill everything if an error was found
  !
  if (error_found) then
    if (mypnum == glb_root) write (iout,100)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  end if
  error_found = fals
  !
  ! Deallocate any local arrays
  !
  if (allocated(ds)) then
    deallocate ( ds , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ds",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(pbc_idxs)) then
    deallocate ( pbc_idxs , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"pbc_idxs",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(pairs)) then
    deallocate ( pairs , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"pairs",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(ifaces)) then
    deallocate ( ifaces , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ifaces",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(jfaces)) then
    deallocate ( jfaces , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"jfaces",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(inodes)) then
    deallocate ( inodes , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"inodes",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(jnodes)) then
    deallocate ( jnodes , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"jnodes",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(node_match)) then
    deallocate ( node_match , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"node_match",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(face_match)) then
    deallocate ( face_match , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"face_match",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(marker)) then
    deallocate ( marker , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"marker",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(msk)) then
    deallocate ( msk , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"msk",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(matched)) then
    deallocate ( matched , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"matched",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(bnd_face)) then
    deallocate ( bnd_face , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"bnd_face",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! Format Statements
  !
  100 format (/,2x,70("*"),/,2x,31("*")," ERROR! ",31("*"),/,2x,70("*"),/)
  101 format (2x,"Number of periodic boundary groups (",i0,") is odd!")
  102 format ("match_jnam /= jname, so resetting matching_bc_name")
  103 format ("could not find matching periodic name")
  104 format ("neither match_jidx or match_jnam were specified and npbc > 2")
  105 format ("the specified matching periodic boundary doesnt have the",/, &
              "reverse matching information specified!")
  106 format (2x,"The number of matched pairs does not equal ", &
                 "the expected number of pairs!",/, &
              4x," matched pairs = ",i0,/, &
              4x,"expected pairs = ",i0)
  107 format (2x,"Not all of the periodic BC groups were able ", &
                 "to find a matching periodic BC group!",//, &
              4x,"BC GROUP #",t18,"GROUP NAME",t53,"MATCHING GROUP NAME")
  108 format (4x,i5,t17,"'",a,"'",t52,"'",a,"'")
  109 format (4x,i5,t17,"'",a,"'",t52,"no match found!")
  110 format (2x,"The number of cells between two matching periodic ", &
                 "boundary groups are not the same!",//, &
              4x,"PERIODIC PAIR #",t25,"GROUP-NAME",t62,"NCELLS")
  111 format (7x,i5,t24,"'",a,"'",t60,i8,/, &
                    t24,"'",a,"'",t60,i8)
  !
end subroutine match_periodic_faces
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function periodic_transformation(xyz,center,angle,translate) &
                          result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xyz(:)
  real(wp), intent(in) :: center(:)
  real(wp), intent(in) :: angle(:)
  real(wp), intent(in) :: translate(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(xyz))
  !
  !.. Local Scalars ..
  integer :: l
  !
continue
  !
  ! NOT SURE WHAT TO DO YET WITH ROTATION SO TRANSLATION IS ONLY THING
  ! USED AT THIS POINT
  !
  do l = 1,size(xyz)
    return_value(l) = xyz(l) + translate(l)
  end do
  !
end function periodic_transformation
!
!###############################################################################
!
subroutine get_boundary_face_nodes(bnd_face)
  !
  !.. Use Statements ..
  use geovar, only : nfbnd,bface
  use geovar, only : nodes_of_cell_ptr,nodes_of_cell
  !
  !.. Formal Arguments ..
  type(bnd_face_t), allocatable, intent(inout) :: bnd_face(:)
  !
  !.. Local Scalars ..
  integer :: nf,n1,nfnods,ierr
  integer :: host_cell,host_geom,host_side
  !
  !.. Local Arrays ..
  integer :: ip(1:4)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_boundary_face_nodes"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (allocated(bnd_face)) then
    deallocate ( bnd_face , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"bnd_face",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( bnd_face(1:nfbnd) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bnd_face",1,__LINE__,__FILE__,ierr,error_message)
  !
  do nf = 1,nfbnd
    !
    host_cell = bface(2,nf)
    host_geom = bface(3,nf)
    host_side = bface(4,nf)
    nfnods = num_face_nodes(host_side,host_geom)
    !
    n1 = nodes_of_cell_ptr(host_cell)+1
    !
    ip(1:nfnods) = side_nodes(1:nfnods,host_side,host_geom)
    !
    bnd_face(nf)%bc_idx = bface(6,nf)
    !
    bnd_face(nf)%geom = cell_face_geom(host_side,host_geom)
    !
    bnd_face(nf)%nodes = nodes_of_cell( n1 + ip(1:nfnods) ) ! F2003 AUTO-REALLOC
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_boundary_face_nodes
!
!###############################################################################
!
subroutine check_periodic_list_indices()
  !
  !.. Use Statements ..
  use geovar,            only : nr,nfbnd
  use geovar,            only : bface,cell,face,xyz,fp
  use order_mod,         only : geom_solpts,maxSP,maxFP
  !
  use interpolation_mod, only : interp
  !
  !.. Local Scalars ..
  integer :: nf,nfp,l1,l2,nlp,r1,r2,nrp,kf,lf,rf,l,i,j,k,pf
  integer :: face_geom,face_order
  integer :: left_cell,left_geom,left_side,left_order
  integer :: right_cell,right_geom,right_side,right_order
  real(wp) :: a,val
  !
  !.. Local Arrays ..
  real(wp) :: txyz(1:maxSP,1:nr)
  real(wp) :: left_facexyz(1:nr,1:maxFP)
  real(wp) :: right_facexyz(1:nr,1:maxFP)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_periodic_list_indices"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  !
  find_periodic_faces: &
  do nf = 1,nfbnd
    !
    if (bface(1,nf) /= bc_periodic) cycle find_periodic_faces
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    nfp = geom_solpts(face_geom,face_order)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  CREATE FACE XYZ FOR LEFT CELL ON CURRENT PERIODIC FACE
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    ! Interpolate the left cell solution points to the periodic boundary face
    !
    left_cell = bface(2,nf)
    left_geom = bface(3,nf)
    left_side = bface(4,nf)
    left_order = cell(left_cell)%order
    !
    l1 = cell(left_cell)%beg_sp
    l2 = cell(left_cell)%end_sp
    nlp = l2-l1+1
    !
    ! Create transpose of left cell xyz array to reduce cache misses
    !
    txyz(:,:) = zero
    do k = 1,nlp
      do l = 1,nr
        txyz(k,l) = xyz(l,k+l1-1)
      end do
    end do
    !
    ! Initialize the face point coordinates for the left cell on the face
    !
    left_facexyz(:,:) = zero
    !
    ! Interpolate the left cell solution point coordinates to the face points
    !
    left_interpolate: &
    do kf = 1,nfp
      !
      lf = face(nf)%left%fp_offset + fp%geord(face_geom,face_order)% &
                                        rot(face(nf)%left%rotation)% &
                                     pts(kf)
      !
      do l = 1,nr
        val = zero
        do i = interp(left_geom,left_order)%toFace(face_order)%row_ptr(lf)+1, &
               interp(left_geom,left_order)%toFace(face_order)%row_ptr(lf+1)
          !
          j = interp(left_geom,left_order)%toFace(face_order)%row_idx(i)
          a = interp(left_geom,left_order)%toFace(face_order)%row_val(i)
          !
          val = val + a * txyz(j,l)
          !
        end do
        left_facexyz(l,kf) = val
      end do
      !
    end do left_interpolate
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  CREATE FACE XYZ FOR RIGHT CELL ON CURRENT PERIODIC FACE
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    ! Interpolate the right cell solution points to the periodic boundary face
    !
    pf = bface(5,nf)
    right_cell = bface(2,pf)
    right_geom = bface(3,pf)
    right_side = bface(4,pf)
    right_order = cell(right_cell)%order
    !
    r1 = cell(right_cell)%beg_sp
    r2 = cell(right_cell)%end_sp
    nrp = r2-r1+1
    !
    ! Create transpose of left cell xyz array to reduce cache misses
    !
    txyz(:,:) = zero
    do k = 1,nrp
      do l = 1,nr
        txyz(k,l) = xyz(l,k+r1-1)
      end do
    end do
    !
    ! Initialize the face point coordinates for the right cell on the face
    !
    right_facexyz(:,:) = zero
    !
    ! Interpolate the right cell solution point coordinates to the face points
    !
    right_interpolate: &
    do kf = 1,nfp
      ! NOTE: The right cell on face nf is the left cell on face pf
      rf = face(pf)%left%fp_offset + fp%geord(face_geom,face_order)% &
                                        rot(face(pf)%left%rotation)% &
                                     pts(kf)
      !
      do l = 1,nr
        val = zero
        do i = interp(right_geom,right_order)%toFace(face_order)%row_ptr(rf)+1,&
               interp(right_geom,right_order)%toFace(face_order)%row_ptr(rf+1)
          !
          j = interp(right_geom,right_order)%toFace(face_order)%row_idx(i)
          a = interp(right_geom,right_order)%toFace(face_order)%row_val(i)
          !
          val = val + a * txyz(j,l)
          !
        end do
        right_facexyz(l,kf) = val
      end do
      !
    end do right_interpolate
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    !  CREATE FACE XYZ FOR RIGHT CELL ON CURRENT PERIODIC FACE
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
  end do find_periodic_faces
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine check_periodic_list_indices
!
!###############################################################################
!
end module periodic_mod
