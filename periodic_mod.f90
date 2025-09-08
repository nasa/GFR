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
  public :: get_pbc_matrix
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
  use geomvar_mod, only : nr,nfbnd,nnode,ncell
  use geomvar_mod, only : bface,xyz_nodes
  use geomvar_mod, only : nodes_of_cell,nodes_of_cell_ptr
  use ovar,        only : bc_input
  use ovar,        only : grid_scaling_factor
  !
  use cgns_types_mod, only : CGLEN
  !
  use connectivity_mod, only : init_tri_face_rotation_values
  use connectivity_mod, only : GetFluxPointRotation_Edge
  use connectivity_mod, only : GetFluxPointRotation_Tria
  use connectivity_mod, only : GetFluxPointRotation_Quad
  use connectivity_mod, only : GetFaceNodes_for_QuadFace
  !
  !.. Local Scalars ..
  integer :: l,n,nf,nc,npbc,npair,npf,npn,nfnods,this_pair,ierr
  integer :: i,ip,iidx,ipf,igeom,match_iidx,num_iidx_cells
  integer :: j,jp,jidx,jpf,jgeom,match_jidx,num_jidx_cells
  integer :: lg,ls,rg,rs,rotation_i2j,rotation_j2i
  integer :: i1,i2,j1,j2,ihost,jhost,np
  !
  character(len=1) :: ips,jps
  !
  logical(lk) :: error_found
  !
  character(len=CGLEN) :: match_jnam,jname,iname,match_inam
  !
  !.. Local Arrays ..
  real(wp) :: xyz_ip(1:nr),xyz_jp(1:nr),xyz_ipt(1:nr)
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
  integer, allocatable :: node_match_i2j(:)
  integer, allocatable :: node_match_j2i(:)
  integer, allocatable :: face_match(:)
  integer, allocatable :: marker(:)
  integer, allocatable :: ipf_nodes(:)
  integer, allocatable :: jpf_nodes(:)
  integer, allocatable :: left_nodes(:)
  integer, allocatable :: right_nodes(:)
  integer, allocatable :: ihost_nodes(:)
  integer, allocatable :: jhost_nodes(:)
  integer, allocatable :: lcell_nodes(:)
  integer, allocatable :: rcell_nodes(:)
  !
  real(wp), allocatable :: ds(:)
  real(wp), allocatable :: pbc_matrix(:,:)
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
  ! If any of the periodic boundary faces are triangles, we will need to
  ! find the face rotations for these faces
  !
  if (any(bnd_face(:)%geom == Geom_Tria)) then
    call init_tri_face_rotation_values()
  end if
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
    match_jnam = uppercase( match_jnam )
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
      jname = uppercase( jname )
      !
      if (match_jnam /= jname) then
        !
        if (mypnum == glb_root) then
          write (iout,102) iidx,match_jidx, &
                  trim(adjustl(bc_input(iidx)%name)), &
                  trim(adjustl(jname)), &
                  trim(adjustl(match_jnam)), &
                  trim(adjustl(bc_input(match_jidx)%periodic%matching_bc_name))
        end if
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
    iname = uppercase( iname )
    !
    match_jidx = bc_input(iidx)%periodic%matching_bc_idx
    match_jnam = trim(adjustl(bc_input(iidx)%periodic%matching_bc_name))
    match_jnam = uppercase( match_jnam )
    !
    ! Use match_jidx to find the match
    !
    jname = trim(adjustl(bc_input(match_jidx)%name))
    jname = uppercase( jname )
    !
    match_iidx = bc_input(match_jidx)%periodic%matching_bc_idx
    match_inam = trim(adjustl(bc_input(match_jidx)%periodic%matching_bc_name))
    match_inam = uppercase( match_inam )
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
  ! Go through each periodic pair and find the matching boundary faces and
  ! save the necessary matching information into bface
  !
  match_pair_faces: &
  do this_pair = 1,npair
    !
    iidx = pbc_idxs( pairs(1,this_pair) )
    jidx = pbc_idxs( pairs(2,this_pair) )
    !
    ! Get the transformation matrix for the periodic BC
    !
    pbc_matrix = get_pbc_matrix(bc_input(iidx)%periodic) ! F2003 AUTO-REALLOC
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
    node_match_i2j = init_array(nnode,0) ! F2003 AUTO-REALLOC
    node_match_j2i = init_array(nnode,0) ! F2003 AUTO-REALLOC
    !
    matched = init_array(npn,fals) ! F2003 AUTO-REALLOC
    ds = init_array(npn,zero) ! F2003 AUTO-REALLOC
    !
    inode_loop: &
    do i = 1,npn
      !
      ds(:) = huge(zero) ! reset all to a large number so we can find the min
      !
      ip = inodes(i)
      !
      ! Check to see if the i-node is also in the list of j-nodes.
      ! This would happen for a pole boundary condition (rotated periodic
      ! boundary conditions that meet at the rotation axis).
      !
      if (any( ip == pack(jnodes(:),mask=(.not.matched(:))) )) then
        !
        ! This i-node is also in the list of j-nodes,
        ! so this is most likely a pole boundary condition.
        ! Find the index of ip/jp within the jnodes array.
        !
        j = findloc(jnodes,ip,dim=1)
        !
        ! Manually set the distance for this index of jnode to zero
        ! to deal with the logic after the end of this conditional.
        !
        ds(j) = zero
        !
      else
        !
        xyz_ip(1:nr) = periodic_transformation(xyz_nodes(1:nr,ip),pbc_matrix)
        !
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
        !
      end if
      !
      matched(j) = true
      !
      ! Mark jp negative if ds(j) is > eps6
      !
      jp = jnodes(j)
      if (ds(j) > eps3) then
        jp = -jp
      end if
      ! 
      node_match_i2j(ip) = jp
      node_match_j2i(abs(jp)) = sign(ip,jp)
      !
    end do inode_loop
    !
    ! If any of the node matches have a negative value, write out all
    ! the information about the node matches and abort.
    ! 
    if (any(node_match_i2j(:) < 0) .or. any(node_match_j2i(:) < 0)) then
      ! 
      if (mypnum == glb_root) then
        matching_node_check: &
        do i = 1,npn
          ip = inodes(i)
          jp = abs(node_match_i2j(ip))
          jps = merge("+","-",node_match_i2j(ip)>0)
          ips = merge("+","-",node_match_j2i(jp)>0)
          xyz_ip(1:nr) = xyz_nodes(1:nr,ip)
          xyz_ipt(1:nr) = periodic_transformation(xyz_ip(1:nr),pbc_matrix)
          xyz_ip(1:nr) = xyz_ip(1:nr)/grid_scaling_factor
          xyz_ipt(1:nr) = xyz_ipt(1:nr)/grid_scaling_factor
          xyz_jp(1:nr) = xyz_nodes(1:nr,jp)/grid_scaling_factor
          if (nr == 3) then
            write (240+this_pair,8243) ips,ip,xyz_ip(1),xyz_ip(2),xyz_ip(3), &
                                          xyz_ipt(1),xyz_ipt(2),xyz_ipt(3), &
                                       jps,jp,xyz_jp(1),xyz_jp(2),xyz_jp(3), &
                                       norm2(xyz_ipt(1:nr)-xyz_jp(1:nr))
          else if (nr == 2) then
            write (240+this_pair,8242) ips,ip,xyz_ip(1),xyz_ip(2), &
                                          xyz_ipt(1),xyz_ipt(2), &
                                       jps,jp,xyz_jp(1),xyz_jp(2), &
                                       norm2(xyz_ipt(1:nr)-xyz_jp(1:nr))
          end if
        end do matching_node_check
      end if
      !
      write (error_message,112)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      !
      8242 format (2x,"ip = ",a1,i8," : [",es16.8,1(",",es16.8), &
                                  "] -> [",es16.8,1(",",es16.8),2x, &
                    "; jp = ",a1,i8," : [",es16.8,1(",",es16.8),4x, &
                    "=> ds = ",es20.12)
      8243 format (2x,"ip = ",a1,i8," : [",es16.8,2(",",es16.8), &
                                  "] -> [",es16.8,2(",",es16.8),2x, &
                    "; jp = ",a1,i8," : [",es16.8,2(",",es16.8),4x, &
                    "=> ds = ",es20.12)
      !
    end if
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
      ! Use node_match_i2j to covert the i-nodes to their matching j-nodes
      mpn(1:nfnods) = node_match_i2j( ipn(1:nfnods) )
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
      ! If this is a triangle face, we need to find the rotation between
      ! the two faces
      !
      rotation_i2j = 0
      rotation_j2i = 0
      !
      if (bnd_face(ipf)%geom == Geom_Edge) then
        !
        ipf_nodes = bnd_face(ipf)%nodes(:) ! F2003 AUTO-REALLOC
        jpf_nodes = bnd_face(jpf)%nodes(:) ! F2003 AUTO-REALLOC
        !
        ! First get the rotation from ipf (left) to jpf (right)
        !
        left_nodes = ipf_nodes(:) ! F2003 AUTO-REALLOC
        right_nodes = node_match_j2i( jpf_nodes(:) ) ! F2003 AUTO-REALLOC
        !
        rotation_i2j = GetFluxPointRotation_Edge(left_nodes,right_nodes)
        !
        ! Next get the rotation from jpf (left) to ipf (right)
        !
        left_nodes = jpf_nodes(:) ! F2003 AUTO-REALLOC
        right_nodes = node_match_i2j( ipf_nodes(:) ) ! F2003 AUTO-REALLOC
        !
        rotation_j2i = GetFluxPointRotation_Edge(left_nodes,right_nodes)
        !
      else if (bnd_face(ipf)%geom == Geom_Tria) then
        !
        ipf_nodes = bnd_face(ipf)%nodes(:) ! F2003 AUTO-REALLOC
        jpf_nodes = bnd_face(jpf)%nodes(:) ! F2003 AUTO-REALLOC
        !
        ! First get the rotation from ipf (left) to jpf (right)
        !
        left_nodes = ipf_nodes(:) ! F2003 AUTO-REALLOC
        right_nodes = node_match_j2i( jpf_nodes(:) ) ! F2003 AUTO-REALLOC
        !
        lg = bface(3,ipf) ! left geom
        ls = bface(4,ipf) ! side/face of left cell
        !
        rg = bface(8,ipf) ! right geom
        rs = bface(9,ipf) ! side/face of right cell
        !
        rotation_i2j = GetFluxPointRotation_Tria(left_nodes,lg,ls, &
                                                 right_nodes,rg,rs)
        !
        ! Next get the rotation from jpf (left) to ipf (right)
        !
        left_nodes = jpf_nodes(:) ! F2003 AUTO-REALLOC
        right_nodes = node_match_i2j( ipf_nodes(:) ) ! F2003 AUTO-REALLOC
        !
        lg = bface(3,jpf) ! left geom
        ls = bface(4,jpf) ! side/face of left cell
        !
        rg = bface(8,jpf) ! right geom
        rs = bface(9,jpf) ! side/face of right cell
        !
        rotation_j2i = GetFluxPointRotation_Tria(left_nodes,lg,ls, &
                                                 right_nodes,rg,rs)
        !
      else if (bnd_face(ipf)%geom == Geom_Quad) then
        !
        np = geom_nodes(bnd_face(ipf)%geom)
        !
        ihost = bface(2,ipf)
        jhost = bface(2,jpf)
        !
        i1 = nodes_of_cell_ptr(ihost)+1
        i2 = nodes_of_cell_ptr(ihost+1)
        !
        j1 = nodes_of_cell_ptr(jhost)+1
        j2 = nodes_of_cell_ptr(jhost+1)
        !
        ihost_nodes = nodes_of_cell(i1:i2) ! F2003 AUTO-REALLOC
        jhost_nodes = nodes_of_cell(j1:j2) ! F2003 AUTO-REALLOC
        !
        left_nodes = init_array(np,0) ! F2003 AUTO-REALLOC
        right_nodes = init_array(np,0) ! F2003 AUTO-REALLOC
        !
        ! First get the rotation from ipf (left) to jpf (right)
        !
        lg = bface(3,ipf) ! left geom
        ls = bface(4,ipf) ! side/face of left cell
        !
        rg = bface(8,ipf) ! right geom
        rs = bface(9,ipf) ! side/face of right cell
        !
        lcell_nodes = ihost_nodes(:) ! F2003 AUTO-REALLOC
        rcell_nodes = jhost_nodes(:) ! F2003 AUTO-REALLOC
        !
        ! Convert the j-face nodes to their matching i-face nodes
        !
        do n = 1,np
          l = side_nodes(n,rs,rg) + 1
          rcell_nodes(l) = node_match_j2i( jhost_nodes(l) )
        end do
       !write (mypnum+6100,931) (ihost_nodes(n), n=1,size(ihost_nodes)), &
       !                        (jhost_nodes(n), n=1,size(jhost_nodes))
       !write (mypnum+5100,931) (lcell_nodes(n), n=1,size(lcell_nodes)), &
       !                        (rcell_nodes(n), n=1,size(rcell_nodes)), &
       !                        ls,rs
       !931 format (2x,8i4,4x,8i4,4x,2i3)
        !
        call GetFaceNodes_for_QuadFace(lcell_nodes,ls, &
                                       rcell_nodes,rs, &
                                       left_nodes(1:np),right_nodes(1:np))
       !write (mypnum+5300,932) (left_nodes(n), n=1,np), &
       !                        (right_nodes(n), n=1,np)
       !932 format (2x,4i4,4x,4i4)
        !
        rotation_i2j = GetFluxPointRotation_Quad(left_nodes(1:np), &
                                                 right_nodes(1:np))
        !
        ! Next get the rotation from jpf (left) to ipf (right)
        !
        lg = bface(3,jpf) ! left geom
        ls = bface(4,jpf) ! side/face of left cell
        !
        rg = bface(8,jpf) ! right geom
        rs = bface(9,jpf) ! side/face of right cell
        !
        lcell_nodes = jhost_nodes(:) ! F2003 AUTO-REALLOC
        rcell_nodes = ihost_nodes(:) ! F2003 AUTO-REALLOC
        !
        ! Convert the j-face nodes to their matching i-face nodes
        !
        do n = 1,np
          l = side_nodes(n,rs,rg) + 1
          rcell_nodes(l) = node_match_i2j( ihost_nodes(l) )
        end do
       !write (mypnum+5200,931) (lcell_nodes(n), n=1,size(lcell_nodes)), &
       !                        (rcell_nodes(n), n=1,size(rcell_nodes)), &
       !                        ls,rs
        !
        call GetFaceNodes_for_QuadFace(lcell_nodes,ls, &
                                       rcell_nodes,rs, &
                                       left_nodes(1:np),right_nodes(1:np))
       !write (mypnum+5400,932) (left_nodes(n), n=1,np), &
       !                        (right_nodes(n), n=1,np)
        !
        rotation_j2i = GetFluxPointRotation_Quad(left_nodes(1:np), &
                                                 right_nodes(1:np))
        !
      end if
      !
      ! Its unsure if the first dimension of bface will be permanently
      ! increased to 11, but if it is, bface(11,:) is meant to store the
      ! bc_index value for the face, so save this information to bface.
      !
      ! rotation_i2j/j2i will be a value from 1-6 for a triangle face or
      ! 1-8 for a quad face. We are doubling up the info stored in
      ! bface(11,:) for periodic faces by adding 10000*rotation to the
      ! bc_idx value. This shouldnt be a problem unless there are more
      ! than 10000 bc_input given in the input file, which seems
      ! extremely unlikely.
      !
      if (size(bface,dim=1) >= 11) then
        bface(11,ipf) = 10000*rotation_i2j + iidx
        bface(11,jpf) = 10000*rotation_j2i + jidx
      end if
      !
    end do
    !
  end do match_pair_faces
  !
  ! Save some information to periodic_faces for later
  !
  npf = count( bface(1,:) == bc_periodic )
  !
  allocate ( periodic_faces(1:npf) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"periodic_faces",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  i = 0
  do nf = 1,nfbnd
    if (bface(1,nf) == bc_periodic) then
      i = i + 1
      periodic_faces(i)%nf = nf
      periodic_faces(i)%bc_idx = bnd_face(nf)%bc_idx
     !periodic_faces(i)%rotation = -99
      periodic_faces(i)%rotation = bface(11,nf)/10000
    end if
  end do
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
  if (allocated(node_match_i2j)) then
    deallocate ( node_match_i2j , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"node_match_i2j",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(node_match_j2i)) then
    deallocate ( node_match_j2i , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"node_match_j2i",2,__LINE__,__FILE__, &
                     ierr,error_message)
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
  102 format ("match_jnam /= jname, so resetting matching_bc_name",//, &
              8x,"ibc = ",i0,/, &
              8x,"jbc = bc_input(ibc)%periodic%matching_bc_idx = ",i0,//, &
              8x,"iname = bc_input(ibc)%name = ",a,/, &
              8x,"jname = bc_input(jbc)%name = ",a,//, &
              8x,"match_jnam = bc_input(ibc)%periodic%matching_bc_name = ",a,/,&
              8x,"match_inam = bc_input(jbc)%periodic%matching_bc_name = ",a,//)
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
  112 format (" A NEGATIVE NODE INDEX WAS FOUND IN EITHER NODE_MATCH_I2J OR", &
              " NODE_MATCH_J2I. THIS MEANS ONE OR MORE NODES ON A PAIR OF", &
              " PERIODIC BOUNDARY FACES WERE UNABLE TO FIND A MATCHING NODE", &
              " FROM THE OPPOSITE FACE. THIS IS MOST LIKELY DUE TO AN ERROR", &
              " IN THE VALUES SPECIFIED IN THE INPUT FILE. FOR EXAMPLE, THE", &
              " TRANSLATION DISTANCE OR ROTATION ANGLE MIGHT SWITCHED TO THE", &
              " OPPOSITE DIRECTION, OR BE SPECIFIED IN THE WRONG COORDINATE", &
              " DIRECTION.")
  !
end subroutine match_periodic_faces
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function periodic_transformation(xyz,pbc_matrix) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xyz(:)
  real(wp), intent(in) :: pbc_matrix(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(xyz))
  !
  !.. Local Scalars ..
  integer :: n,l
  !
continue
  !
  n = size(xyz)
  !
  do l = 1,n
    return_value(l) = dot_product(xyz(1:n),pbc_matrix(1:n,l)) + &
                      pbc_matrix(n+1,l)
  end do
  !
end function periodic_transformation
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_pbc_matrix(periodic,ignore_translation) result(return_value)
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr
  use ovar,        only : periodic_bc_t
  !
  !.. Formal Arguments ..
  type(periodic_bc_t), intent(in) :: periodic
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: ignore_translation
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  real(wp) :: rtn_angle,unx,uny,unz,xc,yc,zc
  real(wp) :: rxx,rxy,rxz,ryx,ryy,ryz,rzx,rzy,rzz
  real(wp) :: c,cx,cy,cz,s,sx,sy,sz
  !
  logical(lk) :: rotation_only
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: rtn_angles(:)
  real(wp), allocatable :: rtn_axis(:)
  real(wp), allocatable :: rmat(:,:)
  !
continue
  !
  rtn_angles = periodic%RotationAngle(1:3) ! F2003 AUTO-REALLOC
  rtn_axis = periodic%RotationAxis(1:3) ! F2003 AUTO-REALLOC
  rtn_angle = periodic%RotationAxis(4)
  !
  rmat = identity_matrix(one,4) ! F2003 AUTO-REALLOC
  !
  xc = periodic%RotationCenter(1)
  yc = periodic%RotationCenter(2)
  zc = periodic%RotationCenter(3)
  !
  if (any(abs(rtn_angles(:)) > eps6)) then
    !
    ! Use the angles of rotation about each coordinate axis to form
    ! the rotation matrix
    ! NOTE: This creates a single rotation matrix using three rotations
    !       about the Cartesian axes, with the three rotation angles not
    !       necessarily being equal. The rotations are performed in the
    !       order of x-axis, then y-axis, and finally z-axis.
    !
    cx = cos( rtn_angles(1)*degrees2radians )
    sx = sin( rtn_angles(1)*degrees2radians )
    !
    cy = cos( rtn_angles(2)*degrees2radians )
    sy = sin( rtn_angles(2)*degrees2radians )
    !
    cz = cos( rtn_angles(3)*degrees2radians )
    sz = sin( rtn_angles(3)*degrees2radians )
    !
    rxx =    cy*cz
    rxy = sx*sy*cz - cx*sz
    rxz = cx*sy*cz + sx*sz
    !
    ryx =    cy*sz
    ryy = sx*sy*sz + cx*cz
    ryz = cx*sy*sz - sx*cz
    !
    rzx =   -sy
    rzy = sx*cy
    rzz = cx*cy
    !
  else if (any(abs(rtn_axis(:)) > eps6) .and. abs(rtn_angle) > eps6) then
    !
    ! Use rotation axis and angle of rotation about the axis to form
    ! the rotation matrix
    ! NOTE: These equations came from a Word document from FUN3D titled
    !       Transformation Matrices. It references the 2009 AIAA paper
    !       "Recent Enhancements to the FUN3D Flow Solver for Moving-Mesh
    !       Applications" by Biedron and Thomas, AIAA 2009-1360.
    !
    c = cos( rtn_angle*degrees2radians )
    s = sin( rtn_angle*degrees2radians )
    !
    ! Compute the components of the normal vector used as the rotation axis
    !
    unx = rtn_axis(1)/norm2(rtn_axis)
    uny = rtn_axis(2)/norm2(rtn_axis)
    unz = rtn_axis(3)/norm2(rtn_axis)
    !
    rxx = (one-c)*unx*unx + c
    ryy = (one-c)*uny*uny + c
    rzz = (one-c)*unz*unz + c
    !
    rxy = (one-c)*unx*uny - unz*s
    rxz = (one-c)*unx*unz + uny*s
    ryz = (one-c)*uny*unz - unx*s
    !
    ryx = (one-c)*unx*uny + unz*s
    rzx = (one-c)*unx*unz - uny*s
    rzy = (one-c)*uny*unz + unx*s
    !
  else
    !
    rxx = one  ; rxy = zero ; rxz = zero ; xc = zero
    ryx = zero ; ryy = one  ; ryz = zero ; yc = zero
    rzx = zero ; rzy = zero ; rzz = one  ; zc = zero
    !
  end if
  !
  rmat(1,1) = rxx
  rmat(2,1) = rxy
  rmat(3,1) = rxz
  rmat(4,1) = xc - rxx*xc - rxy*yc - rxz*zc + periodic%Translation(1)
  !
  rmat(1,2) = ryx
  rmat(2,2) = ryy
  rmat(3,2) = ryz
  rmat(4,2) = yc - ryx*xc - ryy*yc - ryz*zc + periodic%Translation(2)
  !
  rmat(1,3) = rzx
  rmat(2,3) = rzy
  rmat(3,3) = rzz
  rmat(4,3) = zc - rzx*xc - rzy*yc - rzz*zc + periodic%Translation(3)
  !
  ! Copy the relevant parts of rmat to return_value
  !
  rotation_only = fals
  if (present(ignore_translation)) then
    rotation_only = ignore_translation
  end if
  !
  if (rotation_only) then
    return_value = rmat(1:nr,1:nr) ! F2003 AUTO-REALLOC
  else
    return_value = init_array(nr+1,nr,zero) ! F2003 AUTO-REALLOC
    return_value(1:nr,1:nr) = rmat(1:nr,1:nr)
    return_value(nr+1,1:nr) = rmat(4,1:nr)
  end if
  !
end function get_pbc_matrix
!
!###############################################################################
!
subroutine get_boundary_face_nodes(bnd_face)
  !
  !.. Use Statements ..
  use geomvar_mod, only : nfbnd,bface
  use geomvar_mod, only : nodes_of_cell_ptr,nodes_of_cell
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
  use geomvar_mod, only : bc_grps
  use geomvar_mod, only : bface,face,fp
  use geomvar_mod, only : cell_map
  use ovar,        only : bc_in
  !
  !.. Local Scalars ..
  integer :: nf,nfp,i,k,pf,ipf
  integer :: irot,p1,p2,n
  integer :: face_geom,face_order
  integer :: my_cpu,adj_cpu
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_periodic_list_indices"
  !
continue
  !
#ifndef PBS_ENV
#ifdef DEBUG_ON
  n = 0
  do nf = 1,size(bface,dim=2)
    if (bface(1,nf) == bc_cpu_bnd) then
      if (bface(11,nf) > 0) then
        if (bc_in(mod(bface(11,nf),10000))%bc_type == bc_periodic) then
          n = n + 1
        end if
      end if
    end if
  end do
  write (iout,981) mypnum,"MPI",n
  write (iout,981) mypnum,"local",count(bface(1,:)==bc_periodic)
  981 format ("CPU ",i6," : # of ",a5," periodic faces = ",i0)
#endif
#endif
  return
  if (ncpu > 1) return
  if (.not. allocated(periodic_faces)) return
  !
  call debug_timer(entering_procedure,pname)
  !
  do ipf = 1,size(periodic_faces)
    write (mypnum+200,801) ipf, &
                           periodic_faces(ipf)%nf, &
                           periodic_faces(ipf)%bc_idx, &
                           periodic_faces(ipf)%rotation
  end do
  if (ncpu > 1) then
  do nf = 1,size(bface,dim=2)
    write (mypnum+300,801) nf,(bface(k,nf), k=1,size(bface,dim=1))
    my_cpu = mypnum+1
    adj_cpu = bface(7,nf)+1
    write (mypnum+500,801) nf,bface(1,nf), &
                              cell_map(my_cpu)%loc_to_glb(bface(2,nf)), &
                              bface(3,nf), &
                              bface(4,nf), &
                              bface(5,nf), &
                              cell_map(adj_cpu)%loc_to_glb(bface(6,nf)), &
                              bface(7,nf), &
                              bface(8,nf), &
                              bface(9,nf), &
                              (bface(k,nf), k=10,size(bface,dim=1))
  end do
  end if
  801 format (i4,4x,20i5)
  !
  !
  redo_periodic_list: &
  do n = 1,size(bc_grps)
    !
    if (bc_grps(n)%bc_type == bc_periodic .and. &
        allocated(bc_grps(n)%periodic_list)) then
      !
      i = 0
      do ipf = 1,size(periodic_faces)
        !
        nf = periodic_faces(ipf)%nf
        irot = periodic_faces(ipf)%rotation
       !write (iout,'(10i5)') ipf,nf,irot,face(nf)%right%rotation
        !
        face_geom = face(nf)%geom
        face_order = face(nf)%order
        !
        face(nf)%right%rotation = irot
        !
        pf = bface(5,nf)
        !
        p1 = face(pf)%beg_fp
        p2 = face(pf)%end_fp
        nfp = p2-p1+1
        !
        do k = 1,nfp
          i = i + 1
          bc_grps(n)%periodic_list(i) = p1 + fp%geord(face_geom,face_order)% &
                                                                  rot(irot)% &
                                             pts(k)
        end do
        write (mypnum+4400,801) nf,irot,bface(11,nf)/10000
        !
      end do
      !
      exit redo_periodic_list
      !
    end if
    !
  end do redo_periodic_list
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine check_periodic_list_indices
!
!###############################################################################
!
end module periodic_mod
