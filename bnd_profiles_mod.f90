module bnd_profiles_mod
  !
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,ntk,ntl,nq
  !
  implicit none
  !
  private
  !
  !
  !
  ! ##################################
  ! ###  Module Public Procedures  ###
  ! ##################################
  !
  public :: init_bnd_profile_structures
  public :: compute_bnd_profiles
  public :: bnd_profiles_memory_usage
  !
  !
  !
  ! ###########################
  ! ###  Module Parameters  ###
  ! ###########################
  !
  character(len=*), parameter :: fname = "bndprof.tec.dat"
  !
  !
  !
  ! ##############################
  ! ###  Module Derived Types  ###
  ! ##############################
  !
  !
  ! PROFILE_QUANTITY: Derived type containing the name of a profile
  !                   quantity and the corresponding data
  !
  type, public :: profile_quantity
    character(len=100) :: qty_name
    real(wp), allocatable :: qty_data(:)
  end type profile_quantity
  !
  ! PROFILE: Derived type of a single profile containing multiple quantities
  !
  type, public :: profile
    character(len=100) :: prof_name
    integer :: group_num
    integer, allocatable :: faces(:)
    integer, allocatable :: connectivity(:,:)
    type(profile_quantity), allocatable :: qty(:)
  end type profile
  !
  !
  !
  ! ###################################
  ! ###  Module Allocatable Arrays  ###
  ! ###################################
  !
  ! BND_PROF : Solution profiles along each boundary
  ! SOL_PROF : Solution profiles within the interior
  !
  type(profile), public, save, allocatable, dimension(:) :: bnd_prof
  type(profile), public, save, allocatable, dimension(:) :: sol_prof
  !
  !
  !
  ! ################################
  ! ###  Module Local Variables  ###
  ! ################################
  !
  !
contains
!
!###############################################################################
!
subroutine compute_bnd_profiles
  !
  !.. Use Statements ..
  use order_mod,         only : maxSP,maxFP,geom_solpts
  use quadrature_mod,    only : std_elem
  use geovar,            only : nr,bface,face,cell,xyz_nodes
  use geovar,            only : nodes_of_cell,nodes_of_cell_ptr
  use ovar,              only : machref,rhoref,aref,muref,tref
  use ovar,              only : cpref_nd,Pr
  use ovar,              only : output_bnd_profiles
  use metrics_mod,       only : geofa
  use flowvar,           only : usp,dusp
  use interpolation_mod, only : interp
  use interpolation_mod, only : outerp
  !
  !.. Local Scalars ..
  integer :: i,k,l,m,n,nf,np,n1,n2,kf,kl
  integer :: host_cell,host_geom,host_order
  integer :: ndata,face_geom,face_order,nfp
  integer :: ioprf,ierr
  real(wp) :: temp,visc,tauwall
  real(wp) :: cfmag,rex,utau,yplus,qdot
  real(wp) :: dynamic_pres,blasius_cf
  real(wp) :: uinf,eta_coef,vy_non_coef
  real(wp) :: xloc,yloc,vx,vy,nuref
  real(wp) :: eta,vx_non,vy_non
  real(wp) :: cfdrag,cfdrag_total,plate_length
  real(wp) :: cfdrag_ideal,jacfp,rexmax
  logical(lk)  :: use_coefs
  integer  :: nn1,host_side
  !
  character(len=20) :: zonetype
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: unrm
  real(wp), dimension(1:nr) :: xyz_fp
  real(wp), dimension(1:nq) :: cvfp
  real(wp), dimension(1:nq) :: pvfp
  real(wp), dimension(1:nr,1:nq) :: dcvdx
  real(wp), dimension(1:nr,1:nq) :: dpvdx
  real(wp), dimension(1:maxSP,1:nq) :: tusp
  real(wp), dimension(1:maxSP,1:nr,1:nq) :: tdusp
  real(wp), dimension(1:maxSP,1:nr,1:nq) :: tdvsp
  real(wp), dimension(1:nr,1:nq,1:maxSP) :: dvsp
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_bnd_profiles"
  !
continue
  !
  ! If output_bnd_profiles is not set then return without doing anything
  !
  if (.not. output_bnd_profiles) return
  !
  call debug_timer(entering_procedure,pname)
  !
  dynamic_pres = half * machref * machref
  !
  cfdrag_total = zero
  rexmax = -huge(zero)
  !
  ! #########
  ! #########
  ! ##  1  ##
  ! #########
  ! #########
  !
  ! Compute the various profile quantities
  !
  write (iout,1)
  1 format(//,10x,"INTEGRATED SKIN FRICTION DRAG",/,80("-"))
  do n = 1,size(bnd_prof)
    !
    ndata = 0 ! initialize the index within qty_data for this boundary group
    !
    cfdrag = zero ! initialize the integrated skin friction drag
    !
    do i = 1,size(bnd_prof(n)%faces)
      !
      nf = bnd_prof(n)%faces(i) ! boundary face of current profile
      !
      face_geom = face(nf)%geom   ! geometry of boundary face
      face_order = face(nf)%order ! order    of boundary face
      !
      nfp = geom_solpts(face_geom,face_order) ! number of flx-pts
      !
      host_cell = face(nf)%left%cell ! host cell on current boundary face
      !
      host_geom = cell(host_cell)%geom   ! host cell geometry type
      host_order = cell(host_cell)%order ! host cell solution order
      !
      host_side = bface(4,nf)
      !
      n1 = cell(host_cell)%beg_sp ! beginning index for sol-pts in host cell
      n2 = cell(host_cell)%end_sp ! ending    index for sol-pts in host cell
      np = n2-n1+1                ! number of sol-pts in host cell
      !
      ! Create transposes of the usp and dusp arrays to reduce cache misses
      !
      nn1 = nodes_of_cell_ptr(host_cell)+1
      !
      jacfp = norm2( half * (xyz_nodes(:,nodes_of_cell(nn1+ &
                                       side_nodes(1,host_side,host_geom))) - &
                             xyz_nodes(:,nodes_of_cell(nn1+ &
                                       side_nodes(2,host_side,host_geom))) &
                             ) &
                   )
      !
     !dvsp(1:nr,1:nq,1:np) = grad_cv_to_grad_pv( usp(1:nq,n1:n2) , &
     !                                           dusp(1,1:nq,n1:n2), &
     !                                           dusp(2,1:nq,n1:n2))
      !
      tusp(1:np,1:nq) = transpose( usp(1:nq,n1:n2) )
      do k = 1,np
        tdusp(k,1:nr,1:nq) = dusp(1:nr,1:nq,n1-1+k)
     !  tdvsp(k,1:nr,1:nq) = dvsp(1:nr,1:nq,k)
      end do
      !
      ! Compute the profile quantities at each data point for this face
      !
      do k = 1,nfp
        !
        ! Index of flx-pt local to the current host cell
        kl = face(nf)%left%get_fp_idx(face_geom,face_order,k)
        ! Index within qty_data for current flx-pt
        kf = ndata + k
        !
        unrm(1:nr) = geofa(nf)%v(1:nr,k)
        !
        ! Interpolate the conservative variables from the interior sol-pts
        ! to the current flux point on face nf
        !
        do m = 1,nq
          cvfp(m) = interp(host_geom,host_order)% &
                              toFace(face_order)% &
                    dot( kl , tusp(1:np,m) )
          do l = 1,nr
            dcvdx(l,m) = interp(host_geom,host_order)% &
                                   toFace(face_order)% &
                         dot( kl , tdusp(1:np,l,m) )
          ! dpvdx(l,m) = interp(host_geom,host_order)% &
          !                        toFace(face_order)% &
          !              dot( kl , tdvsp(1:np,l,m) )
          end do
        end do
        !
        ! Compute the primitive variables and their
        ! gradients at the current flux point
        !
        pvfp(1:nq) = usp2v_sp( cvfp(1:nq) )
        dpvdx(1:nr,1:nq) = grad_cv_to_grad_pv_sp( cvfp(1:nq) , &
                                                  dcvdx(1:nr,1:nq) )
        !
        temp = temperature_pv_sp(pvfp)
        visc = viscosity_pv_sp(pvfp)
        !
        pvfp(nec) = pvfp(nec) * rhoref
        pvfp(nmb:nme) = pvfp(nmb:nme) * aref
        pvfp(nee) = pvfp(nee) * rhoref * aref * aref
        !
        do m = 1,nq
          bnd_prof(n)%qty(nr+m)%qty_data(kf) = pvfp(m)
        end do
        !
        bnd_prof(n)%qty(nr+nq+1)%qty_data(kf) = temp * tref
        bnd_prof(n)%qty(nr+nq+2)%qty_data(kf) = visc / muref
        !
        ! Compute the other flow variables for the profile
        !
        do m = 1,nr
          xyz_fp(m) = bnd_prof(n)%qty(m)%qty_data(kf)
        end do
        !
        if (any(bface(1,nf) == bc_walls)) then
          !
         !tauwall = visc * &
         !          normal_streamwise_velocity_gradient_pv(pvfp,dpvdx,unrm)
          tauwall = visc * dpvdx(2,nmx)
         !tauwall = muref * dpvdx(2,nmx)
          tauwall = abs(tauwall)
          cfmag = tauwall / dynamic_pres
          rex = rhoref * machref * aref * xyz_fp(1) / muref
         !rex = rhoref * machref * aref * xyz_fp(1) / (visc*aref*rhoref)
          utau = sqrt( tauwall / cvfp(nec) )
          yplus = cvfp(1) * utau * xyz_fp(2) / visc
          blasius_cf = 0.664_wp / sqrt(max(rex,small))
         !qdot = -dot(dpvdx(1:nr,nee),unrm(1:nr)) * cpref_nd * visc / Pr
          !
          bnd_prof(n)%qty(nr+nq+3)%qty_data(kf) = cfmag
          bnd_prof(n)%qty(nr+nq+4)%qty_data(kf) = rex
          bnd_prof(n)%qty(nr+nq+5)%qty_data(kf) = yplus
          bnd_prof(n)%qty(nr+nq+6)%qty_data(kf) = tauwall
          bnd_prof(n)%qty(nr+nq+7)%qty_data(kf) = utau
          bnd_prof(n)%qty(nr+nq+7)%qty_data(kf) = cfmag/blasius_cf
         !bnd_prof(n)%qty(nr+nq+8)%qty_data(kf) = qdot
          bnd_prof(n)%qty(nr+nq+8)%qty_data(kf) = blasius_cf
          !
          if (rex > rexmax) then
            rexmax = rex
            plate_length = xyz_fp(1)
          end if
          if (xyz_fp(1) > 1.0e-4_wp) then
            cfdrag = cfdrag + cfmag*jacfp*std_elem(face_geom,face_order)%wts(k)
          end if
          !
        else
          !
          bnd_prof(n)%qty(nr+nq+3)%qty_data(kf) = zero
          bnd_prof(n)%qty(nr+nq+4)%qty_data(kf) = zero
          bnd_prof(n)%qty(nr+nq+5)%qty_data(kf) = zero
          bnd_prof(n)%qty(nr+nq+6)%qty_data(kf) = zero
          bnd_prof(n)%qty(nr+nq+7)%qty_data(kf) = zero
          bnd_prof(n)%qty(nr+nq+8)%qty_data(kf) = zero
          !
        end if
        !
      end do
      !
      ! Update ndata for the next face in this boundary group
      !
      ndata = ndata + nfp
      !
    end do
    !
    if (cfdrag /= zero) then
      cfdrag = cfdrag!/ plate_length
      cfdrag_total = cfdrag_total + cfdrag
      write (iout,2) bnd_prof(n)%group_num, &
                     trim(adjustl(bnd_prof(n)%prof_name)), &
                     cfdrag
     !write (iout,'(a,es23.15)') "  plate_length = ",plate_length
    end if
    2 format ("Boundary group #",i0," - (",a," B.C.) => cfdrag = ",es23.15)
    !
  end do
  write (iout,3) cfdrag_total
  3 format (80("-"),/,20x,"Total cfdrag = ",es23.15)
  write (iout,4) 1.328_wp / sqrt(rexmax)
  4 format (20x,"Ideal cfdrag = ",es23.15,/,80("-"))
  !
  ! Compute the subsonic outflow profile
  !
  do n = 1,size(bnd_prof)
    !
    if (bc_string_to_integer(bnd_prof(n)%prof_name) == bc_sub_outflow) then
      !
      uinf = maxval( bnd_prof(n)%qty(nr+nmx)%qty_data ) * aref
      !
      do kf = 1,size(bnd_prof(n)%qty(1)%qty_data)
        !
        use_coefs = true
       !use_coefs = fals
        !
        nuref = muref / rhoref
        !
        eta_coef = merge(half,one,use_coefs)
        vy_non_coef = merge(two,one,use_coefs)
        !
        xloc = bnd_prof(n)%qty(1)%qty_data(kf)
        yloc = bnd_prof(n)%qty(2)%qty_data(kf)
        !
        vx = bnd_prof(n)%qty(nr+nmx)%qty_data(kf) * aref
        vy = bnd_prof(n)%qty(nr+nmy)%qty_data(kf) * aref
        rex = uinf * xloc / nuref
        !
        eta = yloc * sqrt( max(eta_coef*uinf/(nuref*xloc),zero) )
        vx_non = vx / uinf
        vy_non = vy * sqrt( vy_non_coef * max(rex,zero) ) / uinf
        !
        bnd_prof(n)%qty(nr+nq+3)%qty_data(kf) = eta
        bnd_prof(n)%qty(nr+nq+4)%qty_data(kf) = vx_non
        bnd_prof(n)%qty(nr+nq+5)%qty_data(kf) = vy_non
        bnd_prof(n)%qty(nr+nq+6)%qty_data(kf) = rex
        bnd_prof(n)%qty(nr+nq+7)%qty_data(kf) = vx
        bnd_prof(n)%qty(nr+nq+8)%qty_data(kf) = vy
        !
      end do
      !
    end if
    !
  end do
  !
  ! #########
  ! #########
  ! ##  2  ##
  ! #########
  ! #########
  !
  ! Go back through the face quantities and interpolate these values
  ! to the output solution points
  !
  do n = 1,size(bnd_prof)
    !
    ndata = 0 ! initialize the index within qty_data for this boundary group
    !
    do i = 1,size(bnd_prof(n)%faces)
      !
      nf = bnd_prof(n)%faces(i) ! boundary face of current profile
      !
      face_geom = face(nf)%geom   ! geometry of boundary face
      face_order = face(nf)%order ! order    of boundary face
      !
      nfp = geom_solpts(face_geom,face_order) ! number of flx-pts
      !
      if (allocated(outerp(face_geom,face_order)%output)) then
        !
        n1 = ndata + 1
        n2 = ndata + nfp
        !
        do m = 1,size(bnd_prof(n)%qty)
          bnd_prof(n)%qty(m)%qty_data(n1:n2) = outerp(face_geom,face_order)% &
                                                                     output% &
                                 mv_mult( bnd_prof(n)%qty(m)%qty_data(n1:n2) )
        end do
        !
      end if
      !
      ! Update ndata for the next face in this boundary group
      !
      ndata = ndata + nfp
      !
    end do
    !
  end do
  !
  ! #########
  ! #########
  ! ##  3  ##
  ! #########
  ! #########
  !
  ! Output the profiles in Tecplot format
  !
  open (newunit=ioprf,file=fname,status="replace",action="write", &
                      form="formatted",position="rewind", &
                      iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (ioprf,100) (trim(bnd_prof(1)%qty(m)%qty_name), &
                          m=1,size(bnd_prof(1)%qty) )
  !
  do n = 1,size(bnd_prof)
    !
    n1 = size(bnd_prof(n)%connectivity,dim=1)
    n2 = size(bnd_prof(n)%connectivity,dim=2)
    np = size(bnd_prof(n)%qty(1)%qty_data)
    !
    zonetype = merge("FEQUADRILATERAL","FELINESEG      ",nr==3)
    !
    nf = bnd_prof(n)%faces(1)
   !if (any(bface(1,nf) == bc_walls)) then
      write (ioprf,101) trim(bnd_prof(n)%prof_name), &
                        np,n2,trim(zonetype)
   !else
   !  write (ioprf,101) trim(bnd_prof(n)%prof_name), &
   !                    np,n2,trim(zonetype),nr+nq+3,nr+nq+8
   !end if
    !
    do m = 1,size(bnd_prof(n)%qty)
      write (ioprf,102) (bnd_prof(n)%qty(m)%qty_data(i), i=1,np)
    end do
    !
    do m = 1,n2
      write (ioprf,103) (bnd_prof(n)%connectivity(i,m), i=1,n1)
    end do
    !
  end do
  !
  close (ioprf,iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  100 format ('VARIABLES =',100(' "',a,'"',:))
 !100 format ('VARIABLES =',*(' "',a,'"',:))
  101 format ('ZONE T="',a,'", N=',i0,', E=',i0, &
              ', DATAPACKING=BLOCK, ZONETYPE=',a,:, &
              ', PASSIVEVARLIST=[',i0,'-',i0,']')
  102 format (5es18.10)
  103 format (4(1x,i0))
  !
end subroutine compute_bnd_profiles
!
!###############################################################################
!
subroutine init_bnd_profile_structures
  !
  !.. Use Statements ..
  use order_mod,         only : maxSP,geom_solpts
  use geovar,            only : nr,nfbnd,nbfai,bface
  use geovar,            only : xyz,face,cell,pst_geom
  use ovar,              only : output_bnd_profiles
  use interpolation_mod, only : interp
  use interpolation_mod, only : outerp
  !
  !.. Local Scalars ..
  integer :: i,k,kl,kf,l,n,nf,np,n1,n2
  integer :: nqty,ndata,nfbnd_nc,ierr
  integer :: min_grp,max_grp,nbcgrps
  integer :: group_id,num_group_faces
  integer :: face_geom,face_order,nfp
  integer :: host_cell,host_geom,host_order
  integer :: flow_qtys,wall_qtys,n_subcells
  integer :: iogrd
  character(len=200) :: array_name
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxSP,1:nr) :: txyz
  !
  !.. Local Allocatable Arrays ..
  logical(lk), allocatable, dimension(:) :: msk
  integer, allocatable, dimension(:) :: grp_ids
  integer, allocatable, dimension(:) :: bface_map
  integer, allocatable, dimension(:,:) :: bface_nc
  character(len=100), allocatable, dimension(:) :: flow_var
  !
  !.. Local Parameters ..
  character(len=*), parameter :: char_xyz(1:3) = ["X","Y","Z"]
  character(len=*), parameter :: pname = "init_bnd_profile_structures"
 !logical(lk), parameter :: test_profile_structures = true
  logical(lk), parameter :: test_profile_structures = fals
  !
continue
  !
  ! If output_bnd_profiles is not set then return without initializing anything
  !
  if (.not. output_bnd_profiles) return
  !
  ! For now, only allow boundary profiles to be created if using a single cpu
  !
  if (ncpu /= 1) then
    output_bnd_profiles = fals
    return
  end if
  !
  call debug_timer(entering_procedure,pname)
  !
  ! #########
  ! #########
  ! ##  1  ##
  ! #########
  ! #########
  !
  ! Get the number of flow quantities and their names
  !
  flow_qtys = nr + nq + 2
  wall_qtys = flow_qtys + 6
  !
  allocate ( flow_var(1:wall_qtys) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"flow_var",1,__LINE__,__FILE__,ierr,error_message)
  !
  do n = 1,nr
    flow_var(n) = "Coordinate"//char_xyz(n)
  end do
  flow_var(nr+1) = "Density"
  do n = 1,nr
    flow_var(nr+n+1) = "Velocity"//char_xyz(n)
  end do
  flow_var(nr+nq+0) = "Pressure"
  flow_var(nr+nq+1) = "Temperature"
  flow_var(nr+nq+2) = "Viscosity"
  flow_var(flow_qtys+1) = "SkinFrictionMagnitude"
  flow_var(flow_qtys+2) = "ReynoldsX"
  flow_var(flow_qtys+3) = "YPlus"
  flow_var(flow_qtys+4) = "WallStress"
  flow_var(flow_qtys+5) = "VelocityFriction"
 !flow_var(flow_qtys+6) = "WallHeatFlux"
  flow_var(flow_qtys+6) = "BlasiusCf"
  !
  ! #########
  ! #########
  ! ##  2  ##
  ! #########
  ! #########
  !
  ! Create the non-communication boundary face msk array
  !
  allocate ( msk(1:nfbnd) , source=fals , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",1,__LINE__,__FILE__,ierr,error_message)
  !
  msk(:) = (bface(1,:)/=bc_periodic .and. bface(1,:)/=bc_cpu_bnd)
  !
  ! Count the number of non-communication boundary faces
  !
  nfbnd_nc = count( msk )
  !
  ! Create the non-communication boundary face mapping array. We need this
  ! to give us the original index for the boundary face (i.e., from 1:nfbnd)
  ! when working on the bface_nc array since the face index within bface_nc
  ! is not necessarily the same as its corresponding index in bface.
  !
  allocate ( bface_map(1:nfbnd_nc) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface_map",1,__LINE__,__FILE__,ierr,error_message)
  !
  bface_map = pack( intseq(1,nfbnd) , mask = msk )
  !
  ! Deallocate msk because it is no longer needed
  !
  deallocate ( msk , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Create the bface array for only non-communication type boundaries
  !
  allocate ( bface_nc(1:nbfai,1:nfbnd_nc) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface_nc",1,__LINE__,__FILE__,ierr,error_message)
  !
  bface_nc = bface(:,bface_map)
  !
  ! Find the number of different boundary condition groups on this processor.
  !
  ! Get minimum and maximum group IDs in bface(5,:)
  ! for the non-communication boundaries.
  !
  min_grp = minval( bface_nc(5,:) )
  max_grp = maxval( bface_nc(5,:) )
  !
  ! Loop through the possible group IDs to see how many BC groups are
  ! present on this processor. At the same time, store any group IDs
  ! that are present.
  ! NOTE: The grp_ids array enters this loop deallocated and the first
  !       call to reallocate will allocate it for the first time. Each
  !       subsequent call to reallocate will increase its size as needed.
  !
  nbcgrps = 0
  do n = min_grp,max_grp
    if (any(bface(5,:) == n)) then
      nbcgrps = nbcgrps + 1
      call reallocate( grp_ids , nbcgrps )
      grp_ids(nbcgrps) = n
    end if
  end do
  !
  ! #########
  ! #########
  ! ##  3  ##
  ! #########
  ! #########
  !
  ! Allocate the top level of the bnd_prof array now that we
  ! know how many BC groups are on this processor
  !
  allocate ( bnd_prof(1:nbcgrps) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bnd_prof",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Loop back through the boundary groups and get the indices
  ! of the boundary faces belonging to each group ID.
  !
  do n = 1,size(grp_ids)
    !
    ! Store this group ID in the bnd_prof array
    !
    group_id = grp_ids(n)
    bnd_prof(n)%group_num = group_id
    !
    ! Count the number of boundary faces with this group ID.
    !
    num_group_faces = count(bface_nc(5,:) == group_id)
    !
    ! Allocate the faces part of the bnd_prof array for this group ID
    !
    allocate ( bnd_prof(n)%faces(1:num_group_faces) , source=0 , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) n,"faces"
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Store the indices for the boundary faces that have this group ID
    ! NOTE: We need to use bface_map as the packing array to map the face
    !       index within bface_nc back to the index for bface.
    !
    bnd_prof(n)%faces = pack( bface_map , mask = (bface_nc(5,:) == group_id) )
    !
    ! Get the name of the BC for this group ID
    ! using the first face in the group.
    !
    nf = bnd_prof(n)%faces(1)
    bnd_prof(n)%prof_name = trim( bc_integer_to_string(bface(1,nf)) )
    !
    ! Allocate the number of flow quantities
    !
    nqty = wall_qtys
   !nqty = merge( wall_qtys , flow_qtys , any(bface(1,nf) == bc_walls) )
    !
    allocate ( bnd_prof(n)%qty(1:nqty) , stat=ierr , errmsg=error_message )
    write (array_name,1) n,"qty"
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
    !
    ! Get the total number of data points for this boundary group by summing
    ! the number of face points for each face that belongs to this group.
    !
    ndata = 0
    do i = 1,size(bnd_prof(n)%faces)
      !
      nf = bnd_prof(n)%faces(i)
      !
      face_geom = face(nf)%geom   ! geometry of the current boundary face
      face_order = face(nf)%order ! order    of the current boundary face
      !
      nfp = geom_solpts(face_geom,face_order) ! number of flx-pts
      !
      ! Add the number of face points for this face to the total
      ! number of data points for this boundary group.
      !
      ndata = ndata + nfp
      !
    end do
    !
    ! Loop over the number of flow quantities and
    ! allocate the data locations for each quantity.
    !
    do i = 1,size(bnd_prof(n)%qty)
      !
      ! Store the name of the flow quantity
      !
      bnd_prof(n)%qty(i)%qty_name = flow_var(i)
      !
      ! Allocate the arrays containing the actual data for each flow quantity
      !
      allocate ( bnd_prof(n)%qty(i)%qty_data(1:ndata) , source=zero , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) n,"qty",i
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
    end do
    !
  end do
  !
  ! Deallocate flow_var, grp_ids, bface_nc, and
  ! bface_map because they no longer needed
  !
  deallocate ( grp_ids , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"grp_ids",2,__LINE__,__FILE__,ierr, &
                   error_message,skip_alloc_pause)
  deallocate ( bface_nc , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface_nc",2,__LINE__,__FILE__,ierr,error_message)
  deallocate ( bface_map , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface_map",2,__LINE__,__FILE__,ierr,error_message)
  deallocate ( flow_var , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"flow_var",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! #########
  ! #########
  ! ##  4  ##
  ! #########
  ! #########
  !
  ! Fill the coordinate quantity sections of boundary profile
  !
  do n = 1,size(bnd_prof)
    !
    ndata = 0 ! initialize the index within qty_data for this boundary group
    !
    do i = 1,size(bnd_prof(n)%faces)
      !
      nf = bnd_prof(n)%faces(i) ! boundary face of current profile
      !
      face_geom = face(nf)%geom   ! geometry of boundary face
      face_order = face(nf)%order ! order    of boundary face
      !
      nfp = geom_solpts(face_geom,face_order) ! number of flx-pts
      !
      host_cell = face(nf)%left%cell ! host cell on current boundary face
      !
      host_geom = cell(host_cell)%geom   ! host cell geometry type
      host_order = cell(host_cell)%order ! host cell solution order
      !
      n1 = cell(host_cell)%beg_sp ! beginning index for sol-pts in host cell
      n2 = cell(host_cell)%end_sp ! ending    index for sol-pts in host cell
      np = n2-n1+1                ! number of sol-pts in host cell
      !
      ! Create transpose of the xyz array to reduce cache misses
      !
      txyz(1:np,1:nr) = transpose( xyz(1:nr,n1:n2) )
      !
      do k = 1,nfp
        !
        ! Index of flx-pt local to the current host cell
        kl = face(nf)%left%get_fp_idx(face_geom,face_order,k)
        ! Index within qty_data for current flx-pt
        kf = ndata + k
        !
        ! Interpolate the coordinates of the interior sol-pts to the
        ! current flux point on face nf
        !
        do l = 1,nr
          bnd_prof(n)%qty(l)%qty_data(kf) = interp(host_geom,host_order)% &
                                                      toFace(face_order)% &
                                            dot( kl , txyz(1:np,l) )
        end do
        !
      end do
      !
      ! Update ndata for the next face in this boundary group
      !
      ndata = ndata + nfp
      !
    end do
    !
  end do
  !
  ! #########
  ! #########
  ! ##  5  ##
  ! #########
  ! #########
  !
  ! Now we need to create the face connectivities for each boundary group
  !
  do n = 1,size(bnd_prof)
    !
    ! First count the total number of subcells for this boundary group
    !
    n_subcells = 0
    !
    do i = 1,size(bnd_prof(n)%faces)
      !
      nf = bnd_prof(n)%faces(i) ! boundary face of current profile
      !
      face_geom = face(nf)%geom   ! geometry of boundary face
      face_order = face(nf)%order ! order    of boundary face
      !
      ! Get the number of subcells for this face geometry at this order and
      ! add it to the total number of subcells for this boundary group.
      ! The number of subcells is basically:
      !     face_order**(spatial dimension of face_geom)
      ! NOTE: Make sure that the number of subcells is at least 1,
      !       particularly for the case of face_order=0
      !
      n_subcells = n_subcells + max( 1 , face_order**geom_dimen(face_geom) )
      !
    end do
    !
    ! Allocate the connectivity array component of bnd_prof
    !
    ! Get the size of the first dimension of the connectivity array.
    ! It will be 2 for 2D or 4 for 3D.
    !
    i = merge(4,2,nr==3)
    !
    allocate ( bnd_prof(n)%connectivity(1:i,1:n_subcells) , source=0 , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) n,"connectivity"
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Now create the subcell connectivities for each face of this boundary group
    !
    np = 1 ! initialize the index within qty_data for this boundary group
    l  = 0 ! initialize the subcell index for this boundary group
    !
    do i = 1,size(bnd_prof(n)%faces)
      !
      nf = bnd_prof(n)%faces(i) ! boundary face of current profile
      !
      face_geom = face(nf)%geom   ! geometry of boundary face
      face_order = face(nf)%order ! order    of boundary face
      !
      nfp = geom_solpts(face_geom,face_order) ! number of flx-pts
      !
      ! Number of subcells for the current face
      !
      n_subcells = max( 1 , face_order**geom_dimen(face_geom) )
      !
      ! Loop through the number of subcells and add the subconnectivity
      ! to the first solution point of the current grid cell to get the
      ! connectivity for the current subcell.
      !
      do k = 1,n_subcells
        !
        l = l + 1
        !
        bnd_prof(n)%connectivity(:,l) = np + &
                             pst_geom(face_geom,face_order)%connectivity(:,k)
        !
      end do
      !
      np = np + nfp
      !
    end do
    !
  end do
  !
  ! #########
  ! #########
  ! ##  6  ##
  ! #########
  ! #########
  !
  ! Output the interpolated coordinates to Tecplot for debugging
  !
  if (test_profile_structures) then
    !
    ! First go back through the face quantities and interpolate these values
    ! to the output solution points
    !
    do n = 1,size(bnd_prof)
      !
      ndata = 0 ! initialize the index within qty_data for this boundary group
      !
      do i = 1,size(bnd_prof(n)%faces)
        !
        nf = bnd_prof(n)%faces(i) ! boundary face of current profile
        !
        face_geom = face(nf)%geom   ! geometry of boundary face
        face_order = face(nf)%order ! order    of boundary face
        !
        nfp = geom_solpts(face_geom,face_order) ! number of flx-pts
        !
        if (allocated(outerp(face_geom,face_order)%output)) then
          !
          n1 = ndata + 1
          n2 = ndata + nfp
          !
          do l = 1,nr
            bnd_prof(n)%qty(l)%qty_data(n1:n2) = outerp(face_geom,face_order)% &
                                                                       output% &
                                 mv_mult( bnd_prof(n)%qty(l)%qty_data(n1:n2) )
          end do
          !
        end if
        !
        ! Update ndata for the next face in this boundary group
        !
        ndata = ndata + nfp
        !
      end do
      !
    end do
    !
    ! Now write the profile coordinates to the Tecplot file
    !
    open (newunit=iogrd,file=fname,status="replace",action="write", &
                        form="formatted",position="rewind", &
                        iostat=ierr,iomsg=error_message)
    call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
    !
    write (iogrd,100) (trim(bnd_prof(1)%qty(l)%qty_name), l=1,nr)
    do n = 1,size(bnd_prof)
      !
      n1 = size(bnd_prof(n)%connectivity,dim=1)
      n2 = size(bnd_prof(n)%connectivity,dim=2)
      np = size(bnd_prof(n)%qty(1)%qty_data)
      !
      if (nr == 2) then
        write (iogrd,101) trim(bnd_prof(n)%prof_name),np,n2,"FELINESEG"
      else if (nr == 3) then
        write (iogrd,101) trim(bnd_prof(n)%prof_name),np,n2,"FEQUADRILATERAL"
      end if
      !
      do l = 1,nr
        write (iogrd,102) (bnd_prof(n)%qty(l)%qty_data(i), i=1,np)
      end do
      !
      do l = 1,n2
        write (iogrd,103) (bnd_prof(n)%connectivity(i,l),i=1,n1)
      end do
      !
    end do
    !
    close (iogrd,iostat=ierr,iomsg=error_message)
    call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
    !
    write (error_message,104)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("bnd_prof(",i0,")%",a,:,"(",i0,")%qty_data")
  2 format (2x,i0,5x,"'",a,"'")
  3 format ("BC Group# = ",i0,": ",a)
  100 format ('VARIABLES =',100(' "',a,'"',:))
  101 format ('ZONE T="',a,'", N=',i0,', E=',i0, &
              ', DATAPACKING=BLOCK, ZONETYPE=',a)
  102 format (5es18.10)
  103 format (4(1x,i0))
  104 format ("Temporary stop after writing bndprof.tec.dat")
  !
end subroutine init_bnd_profile_structures
!
!###############################################################################
!
subroutine bnd_profiles_memory_usage(iunit)
  !
  !.. Use Statements ..
  use module_memory
  !
  !.. Formal Arguments ..
  integer, optional, intent(in) :: iunit
  !
  !.. Local Scalars ..
  integer :: iou,i,j
  character(len=36) :: array_name
  type(memory) :: array,dt
  type(memory) :: total
  !
#ifndef DISABLE_DTIO
continue
  !
  iou = iout
  if (present(iunit)) iou = iunit
  !
  call total%reset
  !
  write (iou,1)
  !
  ! Check the memory of the arrays within this module
  !
  if (allocated(bnd_prof)) then
    !
    call dt%set( storage_size(bnd_prof,kind=inttype) , &
                         size(bnd_prof,kind=inttype) )
    !
    do i = lbound(bnd_prof,dim=1),ubound(bnd_prof,dim=1)
      !
      call dt%add( storage_size(bnd_prof(i)%prof_name,kind=inttype) )
      call dt%add( storage_size(bnd_prof(i)%group_num,kind=inttype) )
      !
      call array%set( storage_size(bnd_prof(i)%faces,kind=inttype) , &
                              size(bnd_prof(i)%faces,kind=inttype) )
      !
      write (array_name,4) "bnd_prof",i,"faces"
      write (iou,2) array_name,array%get_units()
      !
      call dt%add(array)
      !
      call array%set( storage_size(bnd_prof(i)%connectivity,kind=inttype) , &
                              size(bnd_prof(i)%connectivity,kind=inttype) )
      !
      write (array_name,4) "bnd_prof",i,"connectivity"
      write (iou,2) array_name,array%get_units()
      !
      call dt%add(array)
      !
      call array%set( storage_size(bnd_prof(i)%qty,kind=inttype) , &
                              size(bnd_prof(i)%qty,kind=inttype) )
      !
      call dt%add(array)
      !
      do j = lbound(bnd_prof(i)%qty,dim=1),ubound(bnd_prof(i)%qty,dim=1)
        !
        call dt%add( storage_size(bnd_prof(i)%qty(j)%qty_name,kind=inttype) )
        !
        call array%set( &
                   storage_size(bnd_prof(i)%qty(j)%qty_data,kind=inttype) , &
                           size(bnd_prof(i)%qty(j)%qty_data,kind=inttype) )
        !
        write (array_name,4) "bnd_prof",i,"qty",j
        write (iou,2) array_name,array%get_units()
        !
        call dt%add(array)
        !
      end do
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "bnd_prof"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "bnd_prof"
    !
  end if
  !
  ! Check the memory of the arrays within this module
  !
  if (allocated(sol_prof)) then
    !
    call dt%set( storage_size(sol_prof,kind=inttype) , &
                         size(sol_prof,kind=inttype) )
    !
    do i = lbound(sol_prof,dim=1),ubound(sol_prof,dim=1)
      !
      call dt%add( storage_size(sol_prof(i)%prof_name,kind=inttype) )
      call dt%add( storage_size(sol_prof(i)%group_num,kind=inttype) )
      !
      call array%set( storage_size(sol_prof(i)%faces,kind=inttype) , &
                              size(sol_prof(i)%faces,kind=inttype) )
      !
      write (array_name,4) "sol_prof",i,"faces"
      write (iou,2) array_name,array%get_units()
      !
      call dt%add(array)
      !
      call array%set( storage_size(sol_prof(i)%connectivity,kind=inttype) , &
                              size(sol_prof(i)%connectivity,kind=inttype) )
      !
      write (array_name,4) "sol_prof",i,"connectivity"
      write (iou,2) array_name,array%get_units()
      !
      call dt%add(array)
      !
      call array%set( storage_size(sol_prof(i)%qty,kind=inttype) , &
                              size(sol_prof(i)%qty,kind=inttype) )
      !
      call dt%add(array)
      !
      do j = lbound(sol_prof(i)%qty,dim=1),ubound(sol_prof(i)%qty,dim=1)
        !
        call dt%add( storage_size(sol_prof(i)%qty(j)%qty_name,kind=inttype) )
        !
        call array%set( &
                   storage_size(sol_prof(i)%qty(j)%qty_data,kind=inttype) , &
                           size(sol_prof(i)%qty(j)%qty_data,kind=inttype) )
        !
        write (array_name,4) "sol_prof",i,"qty",j
        write (iou,2) array_name,array%get_units()
        !
        call dt%add(array)
        !
      end do
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "sol_prof"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "sol_prof"
    !
  end if
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module bnd_profiles_mod")
  2 format (" Array: ",a,dt)
  3 format (" Total memory = ",f8.2,1x,a)
  4 format (a,"(",i0,")%",a,:,"(",i0,")%qty_data")
  5 format ("'",a,"'")
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  !
#endif
end subroutine bnd_profiles_memory_usage
!
!###############################################################################
!
include "Functions/usp2v.f90"
include "Functions/grad_cv_to_grad_pv.f90"
include "Functions/temperature.f90"
include "Functions/viscosity.f90"
include "Functions/normal_streamwise_velocity_gradient.f90"
include "Functions/bc_string_to_integer.f90"
!
!###############################################################################
!
end module bnd_profiles_mod
