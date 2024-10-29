module bc_flow_stats_mod
  !
!#define PROFILE_ON
  !
  !.. Use Statements ..
  use module_kind_types, global_debug_timer =>      debug_timer, &
                    pure_global_debug_timer => pure_debug_timer
  use eqn_idx_mod, only : nec,nmx,nmy,nmz,npres,ntemp
  !
  implicit none
  !
  private
  !
  public :: compute_bc_stats
  !
  !.. Local Module Parameters ..
  !
  real(wp), parameter :: noz_height = 2.68_wp*0.0254_wp
  real(wp), parameter :: inj_radius = 0.25_wp*0.0254_wp
  real(wp), parameter :: lip_height = 0.05_wp*0.0254_wp + noz_height
  !
  real(wp), parameter :: noz_area = noz_height**2
  real(wp), parameter :: inj_area = PI*(inj_radius**2)
  !
  integer, parameter :: nstats = 11
  !
  !.. Local Module Variables ..
  !
  integer, save :: bc_stats_are_not_initialized = true
  !
  integer, save :: n_bc_stats = 0
  !
  !.. Local Module Derived Types ..
  !
  type :: ave_t
    real(wp) :: r  = zero
    real(wp) :: p0 = zero
    real(wp) :: p  = zero
    real(wp) :: t0 = zero
    real(wp) :: t  = zero
  end type ave_t
  !
  type :: bc_stat_t
    integer  :: n
    real(wp) :: area
    real(wp) :: mdot
    type(ave_t) :: ave
    type(ave_t) :: mfr_ave
    integer, allocatable :: my_faces(:)
  end type bc_stat_t
  !
  !.. Local Module Arrays ..
  !
  type(bc_stat_t), allocatable, save :: bc_stats(:)
  !
  real(wp), allocatable :: all_stats(:,:)
  !
contains
!
!###############################################################################
!
subroutine compute_bc_stats()
  !
  !.. Use Statements ..
  use ovar,             only : gam,rhoref,aref,tref
  use ovar,             only : moving_ref_frame
  use geovar,           only : face
  use flowvar,          only : faceusp,facegpsp,facexyz
  use metrics_mod,      only : facenrm
  use element_type_mod, only : std_elem
  !
  !.. Formal Arguments ..
  !
  !.. Optional Arguments ..
  !
  !.. Local Scalars ..
  integer  :: ibc,i0,n,nf,sf,f0,nfp,k
  integer  :: face_geom,face_order
  real(wp) :: rsgn,wt,rnum,mdot
  real(wp) :: x,rnx,x_r,vx,vx_r
  real(wp) :: y,rny,y_r,vy,vy_r
  real(wp) :: z,rnz,z_r,vz,vz_r
  real(wp) :: r,p,t,q2,a2,m2,p0,t0,total_cond
  real(wp) :: gm1,gam_over_gm1,my_mdot
  !
  !.. Local Derived Type Scalars ..
  type(ave_t) :: my_ave
  type(ave_t) :: my_mfr_ave
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_bc_stats"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (bc_stats_are_not_initialized) then
    call init_bc_stats()
    bc_stats_are_not_initialized = fals
  end if
  !
  if (.not. allocated(bc_stats)) then
    return
    call debug_timer(leaving_procedure,pname)
  end if
  !
  gm1 = gam - one
  gam_over_gm1 = gam/gm1
  !
  all_stats(:,:) = zero
  !
  call moving_ref_frame%get_rotation(x_r,y_r,z_r)
  !
  ! Compute the mass flow rate of the injector outflow surface
  !
  loop_over_bcs: &
  do ibc = 1,size(bc_stats)
    !
    my_mdot = zero
    !
    loop_over_my_bc_faces: &
    do n = 1,size(bc_stats(ibc)%my_faces)
      !
      nf = abs(bc_stats(ibc)%my_faces(n))
      sf = merge( 1 , 2 , bc_stats(ibc)%my_faces(n) > 0 )
      !
      rsgn = merge( one , -one , sf == 1 )
      !
      face_geom = face(nf)%geom
      face_order = face(nf)%order
      !
      f0  = face(nf)%beg_fp - 1
      nfp = face(nf)%end_fp - f0
      !
      loop_over_face_fp: &
      do k = 1,nfp
        !
        x = facexyz(1,f0+k)
        y = facexyz(2,f0+k)
        z = facexyz(3,f0+k)
        !
        rnx = facenrm(1,f0+k)
        rny = facenrm(2,f0+k)
        rnz = facenrm(3,f0+k)
        !
        r = faceusp(f0+k,nec,sf)
        p = facegpsp(f0+k,npres,sf)
        t = facegpsp(f0+k,ntemp,sf)
        !
        vx_r = y_r*z - z_r*y
        vy_r = z_r*x - x_r*z
        vz_r = x_r*y - y_r*x
        !
        vx = faceusp(f0+k,nmx,sf)/r - vx_r
        vy = faceusp(f0+k,nmy,sf)/r - vy_r
        vz = faceusp(f0+k,nmz,sf)/r - vz_r
        !
        q2 = vx*vx + vy*vy + vz*vz
        a2 = max( gam*p/r , rndoff )
        !
        m2 = q2/a2
        !
        total_cond = one + half*gm1*m2
        !
        mdot = rsgn*r*(vx*rnx + vy*rny + vz*rnz)
        p0 = p*(total_cond**gam_over_gm1)
        t0 = t*total_cond
        !
        wt = std_elem(face_geom,face_order)%wts(k)
        !
        my_mdot = my_mdot + wt*mdot
        !
        my_ave%r  = my_ave%p0 + r
        my_ave%p0 = my_ave%p0 + p0
        my_ave%p  = my_ave%p  + p
        my_ave%t0 = my_ave%t0 + t0
        my_ave%t  = my_ave%t  + t
        !
        my_mfr_ave%r  = my_mfr_ave%p0 + mdot*r
        my_mfr_ave%p0 = my_mfr_ave%p0 + mdot*p0
        my_mfr_ave%p  = my_mfr_ave%p  + mdot*p
        my_mfr_ave%t0 = my_mfr_ave%t0 + mdot*t0
        my_mfr_ave%t  = my_mfr_ave%t  + mdot*t
        !
      end do loop_over_face_fp
      !
    end do loop_over_my_bc_faces
    !
    i0 = (ibc-1)*nstats
    !
    all_stats(i0+ 1,mypnum+1) = my_mdot
    all_stats(i0+ 2,mypnum+1) = my_ave%r
    all_stats(i0+ 3,mypnum+1) = my_ave%p0
    all_stats(i0+ 4,mypnum+1) = my_ave%p
    all_stats(i0+ 5,mypnum+1) = my_ave%t0
    all_stats(i0+ 6,mypnum+1) = my_ave%t
    all_stats(i0+ 7,mypnum+1) = my_mfr_ave%r
    all_stats(i0+ 8,mypnum+1) = my_mfr_ave%p0
    all_stats(i0+ 9,mypnum+1) = my_mfr_ave%p
    all_stats(i0+10,mypnum+1) = my_mfr_ave%t0
    all_stats(i0+11,mypnum+1) = my_mfr_ave%t
    !
  end do loop_over_bcs
  !
  if (ncpu > 1) then
    call mpi_allgather(MPI_IN_PLACE,0_int_mpi,MPI_DATATYPE_NULL, &
                       all_stats,size(all_stats,dim=1,kind=INT_MPI), &
                       mpi_flttyp,MPI_COMM_WORLD,mpierr)
  end if
  !
  do ibc = 1,size(bc_stats)
    !
    i0 = (ibc-1)*nstats
    !
    rnum = real(bc_stats(ibc)%n,kind=wp)
    !
    ! Compute the total mass flow rate by summing
    ! all the individual mass flow rates
    !
    mdot = sum(all_stats(i0+1,:)) 
    !
    ! Finish computing regular average and then dimensionalize
    !
    my_ave%r  = (sum(all_stats(i0+2,:)) / rnum) * rhoref
    my_ave%p0 = (sum(all_stats(i0+3,:)) / rnum) * rhoref*aref*aref
    my_ave%p  = (sum(all_stats(i0+4,:)) / rnum) * rhoref*aref*aref
    my_ave%t0 = (sum(all_stats(i0+5,:)) / rnum) * tref
    my_ave%t  = (sum(all_stats(i0+6,:)) / rnum) * tref
    !
    ! Finish computing the mass-flow-rated average and then dimensionalize
    ! NOTE: Need to make sure that the mdot average is completed before
    !       dimensionalizing either mdot or the average.
    !
    my_mfr_ave%r  = (sum(all_stats(i0+ 7,:)) / mdot) * rhoref
    my_mfr_ave%p0 = (sum(all_stats(i0+ 8,:)) / mdot) * rhoref*aref*aref
    my_mfr_ave%p  = (sum(all_stats(i0+ 9,:)) / mdot) * rhoref*aref*aref
    my_mfr_ave%t0 = (sum(all_stats(i0+10,:)) / mdot) * tref
    my_mfr_ave%t  = (sum(all_stats(i0+11,:)) / mdot) * tref
    !
    ! Dimensionalize the mass flow rate and save it to bc_stats
    !
    bc_stats(ibc)%mdot = mdot * rhoref*aref / my_ave%r
    !
    ! Save the averaged quantities to bc_stats
    !
    bc_stats(ibc)%ave = my_ave
    bc_stats(ibc)%mfr_ave = my_mfr_ave

  end do
  !
 !blowing_ratio = (inj_mdot/inj_area) / (noz_mdot/noz_area)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_bc_stats
!
!###############################################################################
!
subroutine init_bc_stats()
  !
  !.. Use Statements ..
  use geovar,           only : nfbnd,bface,cell,face,bc_grps
  use metrics_mod,      only : facenrm
  use element_type_mod, only : std_elem
  !
  !.. Local Scalars ..
  integer  :: n,m,m0,nc,ncf,nf,nfp,k,f0,sf,ierr
  integer  :: grp_idx,bc_type,bc_idx
  integer  :: n_inflow,n_outflow,n_extra
  integer  :: face_geom,face_order
  real(wp) :: sum_inflow_area,exp_inflow_area,area
  real(wp) :: sum_outflow_area,exp_outflow_area,my_area
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: msk(:)
  integer, allocatable :: my_faces(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_bc_stats"
  !
  logical(lk), parameter :: compute_surface_areas = fals
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  msk = init_array(nfbnd,0) ! F2003 AUTO-REALLOC
  !
  ! Loop through the BC groups and count the number of BCs to compute stats for
  !
  n_inflow = 0
  n_outflow = 0
  n_extra = 0
  !
  do grp_idx = 1,size(bc_grps)
    !
    bc_type = bc_grps(grp_idx)%bc_type
    bc_idx = bc_grps(grp_idx)%bc_idx
    !
    ! Skip any communication boundary faces
    if (any(bc_type == bc_comm)) cycle
    !
    if (any(bc_type == bc_outflows)) then
      ! One or more outflow BCs exist
      ! NOTE: mixed inflow/outflow are assumed outflow
      n_outflow = 1
    else if (any(bc_type == bc_inflows)) then
      ! One or more inflow BCs exist
      n_inflow = 1
    end if
    !
  end do
  !
  ! Add up the number of bc_stats groups
  !
  m = n_inflow + n_outflow + n_extra
  !
  if (m > 0) then
    allocate ( bc_stats(1:m) , stat=ierr )
  else
    call debug_timer(leaving_procedure,pname)
    return
  end if
  m = 0
  !
  ! Find the boundary faces in the inflow boundary groups
  !
  if (n_inflow > 0) then
    !
    m = m + 1
    msk(:) = 0
    !
    do grp_idx = 1,size(bc_grps)
      !
      bc_type = bc_grps(grp_idx)%bc_type
      bc_idx = bc_grps(grp_idx)%bc_idx
      !
      ! Make sure mixed inflow/outflow are assumed outflow
      if (any(bc_type == bc_outflows)) cycle
      ! Only continue if this is an inflow
      if (all(bc_type /= bc_inflows)) cycle
      !
      do nf = 1,nfbnd
        if (bface(1,nf) == bc_type) then
          if (bface(5,nf) == grp_idx) then
         !if (bface(6,nf) == bc_idx) then
            !
            nc = bface(2,nf)
            ncf = bface(4,nf)
            sf = cell(nc)%face(ncf)%side
            !
            msk(nf) = merge(1,-1,sf==1)
            !
          end if
        end if
      end do
      !
    end do
    !
    ! Get the list of boundary faces with an inflow BC
    !
    my_faces = pack( intseq(1,nfbnd) , msk /= 0 ) ! F2003 AUTO-REALLOC
    !
    ! Multiply the list of faces by msk which will make the face index
    ! either positive or negative depending on which side of the face
    ! the host cell is on.
    ! NOTE: I think msk should be only positive since the host cell should
    !       always be the left cell on the face
    !
    my_faces(:) = my_faces(:) * msk(my_faces(:))
    !
    ! Copy the local my_faces array to bc_stats(m) using move_alloc
    !
    call move_alloc( from=my_faces , to=bc_stats(m)%my_faces )
    !
  end if
  !
  ! Find the boundary faces in the outflow boundary groups
  !
  if (n_outflow > 0) then
    !
    m = m + 1
    msk(:) = 0
    !
    do grp_idx = 1,size(bc_grps)
      !
      bc_type = bc_grps(grp_idx)%bc_type
      bc_idx = bc_grps(grp_idx)%bc_idx
      !
      ! Only continue if this is an outflow
      if (all(bc_type /= bc_outflows)) cycle
      !
      do nf = 1,nfbnd
        if (bface(1,nf) == bc_type) then
          if (bface(5,nf) == grp_idx) then
         !if (bface(6,nf) == bc_idx) then
            !
            nc = bface(2,nf)
            ncf = bface(4,nf)
            sf = cell(nc)%face(ncf)%side
            !
            msk(nf) = merge(1,-1,sf==1)
            !
          end if
        end if
      end do
      !
    end do
    !
    ! Get the list of boundary faces with an outflow BC
    !
    my_faces = pack( intseq(1,nfbnd) , msk /= 0 ) ! F2003 AUTO-REALLOC
    !
    ! Multiply the list of faces by msk which will make the face index
    ! either positive or negative depending on which side of the face
    ! the host cell is on.
    ! NOTE: I think msk should be only positive since the host cell should
    !       always be the left cell on the face
    !
    my_faces(:) = my_faces(:) * msk(my_faces(:))
    !
    ! Copy the local my_faces array to bc_stats(m) using move_alloc
    !
    call move_alloc( from=my_faces , to=bc_stats(m)%my_faces )
    !
  end if
  !
 !if (n_extra > 0) then
 !  !
 !  m = m + 1
 !  msk(:) = 0
 !  !
 !  do grp_idx = 1,size(bc_grps)
 !    !
 !    bc_type = bc_grps(grp_idx)%bc_type
 !    bc_idx = bc_grps(grp_idx)%bc_idx
 !    !
 !    ! Only continue if this is an outflow
 !    if (all(bc_type /= bc_outflows)) cycle
 !    !
 !    do nf = 1,nfbnd
 !      if (bface(1,nf) == bc_type) then
 !        if (bface(5,nf) == grp_idx) then
 !       !if (bface(6,nf) == bc_idx) then
 !          !
 !          nc = bface(2,nf)
 !          ncf = bface(4,nf)
 !          sf = cell(nc)%face(ncf)%side
 !          !
 !          msk(nf) = merge(1,-1,sf==1)
 !          !
 !        end if
 !      end if
 !    end do
 !    !
 !  end do
 !  !
 !  my_faces = pack( intseq(1,nfbnd) , msk /= 0 ) ! F2003 AUTO-REALLOC
 !  !
 !  my_faces(:) = my_faces(:) * msk(my_faces(:))
 !  !
 !  call move_alloc( from=my_faces , to=bc_stats(m)%my_faces )
 !  !
 !end do
  !
  ! Allocate the array for collecting the stats from all MPI processes
  ! using mpi_allgather
  !
  n = nstats*size(bc_stats)
  !
  all_stats = init_array(n,ncpu,zero) ! F2003 AUTO-REALLOC
  !
  do m = 1,size(bc_stats)
    !
    ! Compute the area of the inflow faces on this MPI process
    !
    my_area = zero
    !
    do n = 1,size(bc_stats(m)%my_faces)
      !
      nf = abs(bc_stats(m)%my_faces(n))
      !
      face_geom = face(nf)%geom
      face_order = face(nf)%order
      !
      f0  = face(nf)%beg_fp - 1
      nfp = face(nf)%end_fp - f0
      !
      do k = 1,nfp
        my_area = my_area + norm2(facenrm(:,f0+k)) * &
                            std_elem(face_geom,face_order)%wts(k)
      end do
      !
    end do
    !
    ! Collect the individual surface area and the
    ! total number of boundary faces from all MPI processes
    !
    m0 = (m-1)*nstats
    !
    all_stats(m0+1,mypnum+1) = my_area
    all_stats(m0+2,mypnum+1) = real(size(bc_stats(m)%my_faces),kind=wp)
    !
  end do
  !
  if (ncpu > 1) then
    call mpi_allgather(MPI_IN_PLACE,0_int_mpi,MPI_DATATYPE_NULL, &
                       all_stats,size(all_stats,dim=1,kind=INT_MPI), &
                       mpi_flttyp,MPI_COMM_WORLD,mpierr)
  end if
  !
  ! Sum the total number of boundary faces and surface area for each bc_stat
  !
  n = 0
  area = zero
  do m = 1,size(bc_stats)
    !
    m0 = (m-1)*nstats
    !
    do k = 1,ncpu
      area = area + all_stats(m0+1,k)
      n = n + nint(all_stats(m0+2,k))
    end do
    !
    bc_stats(m)%n = n
    bc_stats(m)%area = area
    !
  end do
  !
  !
  !
  if (compute_surface_areas) then
    !
    ! Get the total inflow and outflow surface areas
    !
    sum_inflow_area = bc_stats(1)%area
    sum_outflow_area = bc_stats(2)%area
    !
    if (mypnum == glb_root) then
      !
      sum_inflow_area = sum_inflow_area / (0.0254_wp**2)
      sum_outflow_area = sum_outflow_area / (0.0254_wp**2)
      exp_inflow_area = PI*one4*half
      exp_outflow_area = sum_outflow_area / (0.0254_wp**2)
      !
      write (iout,1) sum_inflow_area,exp_inflow_area, &
                   ten2*abs(sum_inflow_area-exp_inflow_area)/exp_inflow_area, &
                   sum_outflow_area,exp_outflow_area, &
                   ten2*abs(sum_outflow_area-exp_outflow_area)/exp_outflow_area
      !
      flush (iout)
      !
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/, &
            80("="),/, &
            2x,"Surface",4x,"Computed Area",3x,"Expected Area",3x,"% Diff",/, &
            80("-"),/, &
            2x," Inflow : ",es14.6,2x,es14.6,2x,f9.2,/, &
            2x,"Outflow : ",es14.6,2x,es14.6,2x,f9.2,/, &
            80("="),/)
  !
end subroutine init_bc_stats
!
!###############################################################################
!
#include "Functions/debug_timer.f90"
!
!###############################################################################
!
end module bc_flow_stats_mod
