module metrics_mod
  !
  !.. Global Use Statements ..
  use module_kind_types
  !
  implicit none
  !
  private
  !
  !
  ! ##################################
  ! ###  Module Public Procedures  ###
  ! ##################################
  !
  public :: get_metrics
  public :: metrics_memory_usage
  !
  type, public :: face_nrm
    real(wp), allocatable :: v(:,:)
  end type face_nrm
  !
  ! geofa : normal of the face at each of the solution points on the face
  !
  type(face_nrm), public, save, allocatable :: geofa(:)
  !
  logical(lk), parameter :: use_dJ = fals
  !
  !     For solution point 'n' of cell 'nc',
  !         metrics_dt(nc)%met(n,l,i) gives the partial derivative
  !                             d(xi_l)
  !                             -------
  !                             d(x_i)
  !     where ( x_1, x_2, x_3) = (x,y,z) -> the physical coordinates
  !     where (xi_1,xi_2,xi_3) = (xi,eta,zeta) or (r,s,t) -> the computational
  !                           coordinates within the standard elemental region
  !           NOTE: the first index l gives the index of the numerator and
  !                 the second index i gives the index of the denominator
  !           NOTE: For 2D z and zeta (or t) do not exist
  !
  type, public :: metrics_t
    ! met : first derivative metrics terms at the solution points
    !       SEE ABOVE DESCRIPTION
    real(wp), allocatable :: met(:,:,:)
    ! jac : jacobian of the solution points
    real(wp), allocatable :: jac(:)
    ! dJ : gradient of the jacobian at the solution points
    real(wp), allocatable :: dJ(:,:)
  contains
    ! vol : type-bound procedure to compute volume of a grid cell by
    !       integerating the jacobian at the solution points using the
    !       appropriate quadrature rule
    procedure, public  :: vol => compute_cell_volume
    procedure, private :: create_cell_metrics
  end type metrics_t
  !
  type(metrics_t), public, save, allocatable :: metrics_dt(:)
  !
  !
  !
  ! DEFINITIONS FOR PREVIOUS VERSION OF METRICS ARRAY
  !     For solution point 'n', metrics(l,i,n,1) gives the partial derivative
  !                             d (x_l)
  !                             -------
  !                             d(xi_i)
  !                    whereas, metrics(l,i,n,2) gives the partial derivative
  !                             d(xi_l)
  !                             ------- * jacobian(n)
  !                             d(x_i)
  !     where ( x_1, x_2, x_3) = (x,y,z) -> the physical coordinates
  !     where (xi_1,xi_2,xi_3) = (xi,eta,zeta) or (r,s,t) -> the computational
  !                           coordinates within the standard elemental region
  !           NOTE: the first index l gives the index of the numerator and
  !                 the second index i gives the index of the denominator
  !           NOTE: For 2D z and zeta (or t) do not exist
  !           NOTE: The jacobian(n) term is implicitly included in the
  !                 metrics(:,:,:,2) part of the array because multiplication
  !                 by the jacobian is desired where this metric term is
  !                 predominantly used.
  !
  ! this is essentially metrics(:,:,n) = metrics(:,:,n,2) / jacobian(n)
  !                            d(xi_l)
  ! or simply metrics(l,i,n) = -------
  !                            d(x_i)
  !
  ! spmat : solution points matrix
  ! dximat : dx matrix
  ! detmat : dy matrix
  ! spcemat : curved edge CE (at bottom), solution points
  ! spceximat : CE metrics d_xi
  ! spceetmat : CE metrics d_eta
  ! spcexi2mat : CE metrics d_xi2
  ! spceet2mat : CE metrics d_eta2
  ! spcexietmat : CE metrics d_xieta
  !
  ! lagrange : gives the Lagrange polynomial coefficients to project
  !            variables/fluxes from the solution points in the cell
  !            interior to the flux points on the cell faces
  !
 !real(wp), public, save, allocatable :: spmat(:,:)       !(3,ksp2d)
 !real(wp), public, save, allocatable :: spcemat(:,:)     !(ksp1d+1,ksp2d)
 !real(wp), public, save, allocatable :: spceximat(:,:)   !(ksp1d+1,ksp2d)
 !real(wp), public, save, allocatable :: spceetmat(:,:)   !(ksp1d+1,ksp2d)
 !real(wp), public, save, allocatable :: spcexi2mat(:,:)  !(ksp1d+1,ksp2d)
 !real(wp), public, save, allocatable :: spceet2mat(:,:)  !(ksp1d+1,ksp2d)
 !real(wp), public, save, allocatable :: spcexietmat(:,:) !(ksp1d+1,ksp2d)
  !
  ! wsp  : quadrature weights
  ! xiSP : xi coordinates of the edge solution points
  !        !! xisp and xisp1d seemed to be duplicates
  ! xiCE :
  !
  ! Hesthaven and Warburton equivalencies
  ! wsp = weights
  !
 !real(wp), public, save, allocatable, dimension(:) :: wsp  !(ksp2d)
 !real(wp), public, save, allocatable, dimension(:) :: xiSP !(ksp1d)
 !real(wp), public, save, allocatable, dimension(:) :: xiCE !(ksp1d)
  !
  ! XYSP    : x,y coordinates of the cell solution points
  ! XYSPxi  : d_(x,y)/d_xi  of the cell solution points, i.e. d_x/d_xi, d_y/d_xi
  ! XYSPet  : d_(x,y)/d_eta of the cell solution points
  ! XYeSP   : x,y coordinates of the edge solution points !! not used
  ! XYeSPxi : d_(x,y)/d_xi of the edge solution points
  !
  ! Hesthaven and Warburton equivalencies
  ! XYSP(1,:,:)   = x
  ! XYSP(2,:,:)   = y
  ! XYSPxi(1,:,:) = xr = matmul(Dr,x) -> rx =  ys/J
  ! XYSPxi(2,:,:) = yr = matmul(Dr,y) -> ry = -xs/J
  ! XYSPet(1,:,:) = xs = matmul(Ds,x) -> sx = -yr/J
  ! XYSPet(2,:,:) = ys = matmul(Ds,y) -> sy =  xr/J
  !
 !real(wp), public, save, allocatable :: XYSP(:,:,:)    !(2,ksp2d,-nfbnd:ncell)
 !real(wp), public, save, allocatable :: XYSPxi(:,:,:)  !(2,ksp2d,-nfbnd:ncell)
 !real(wp), public, save, allocatable :: XYSPet(:,:,:)  !(2,ksp2d,-nfbnd:ncell)
 !real(wp), public, save, allocatable :: XYeSP(:,:,:)   !(2,ksp1d,1:nface)
 !real(wp), public, save, allocatable :: XYeSPxi(:,:,:) !(2,ksp1d,1:nface)
  !
  ! DxyDrs : d_(x,y)/d_(r,s) of the cell solution points
  !
 !real(wp), public, save, allocatable, dimension(:,:,:) :: DxyDrs
  !
  ! rjaSP   : jacobian of the solution points
  ! rjaSPxi : d_jacobian/d_xi of the solution points
  ! rjaSPet : d_jacobian/d_eta of the solution points
  !
  ! Hesthaven and Warburthon equivalencies
  ! rjaSP = J = -xs*yr + xr*ys
  !
 !real(wp), public, save, allocatable :: rjaSP(:,:)   !(ksp2d,-nfbnd:ncell)
 !real(wp), public, save, allocatable :: rjaSPxi(:,:) !(ksp2d,-nfbnd:ncell)
 !real(wp), public, save, allocatable :: rjaSPet(:,:) !(ksp2d,-nfbnd:ncell)
  !
  !... for boundary cells, if id_be(j) = 0, then the cell has a curved edge
  !
  ! be_cSPind(1,k,j) k = 1,ksp: the solution point index of cell be_ic(j)
  ! be_cSPind(2,k,j) k = 1,ksp: the solution point index of exterior cell
  ! be_cSPind(3,k,j) k = 1,ksp: the solution point index of cell be_pc(j)
  ! e_cSPind(1,k,j) k = 1,ksp: the solution point index for the left cell
  ! e_cSPind(2,k,j) k = 1,ksp: the solution point index for the right cell
 !integer, public, save, allocatable :: e_cSPind(:,:,:)  !(2,ksp1d,nface)
 !integer, public, save, allocatable :: be_cSPind(:,:,:) !(3,ksp1d,nfbnd)
  !
contains
!
!###############################################################################
!
subroutine get_metrics()
  !
  !.. Use Statements ..
  use geovar, only : nr,n_totpts,n_totcel
  use geovar, only : ncell,cell,xyz
  use geovar, only : nnode,nfbnd,nface,n_solpts
  !
  use derivatives_mod, only : deriv
  !
  use quadrature_mod, only : std_elem
  !
  use ovar, only : local_volume,total_volume
  use ovar, only : VortexGrid_random_perturb
  !
  use order_mod, only : maxSP
  !
  !.. Local Scalars ..
  integer  :: i,k,l,ierr,j1,j2,ibeg,iend
  integer  :: n,nc,n1,n2,np
  integer  :: this_geom,this_order
  integer  :: ref_dir,phys_dir
  integer  :: min_jsp,min_jcl,min_jac_cpu,min_vol_cpu
  integer  :: max_jsp,max_jcl,max_jac_cpu,max_vol_cpu
  integer  :: min_ncell_cpu,max_ncell_cpu
  integer  :: min_nnode_cpu,max_nnode_cpu
  integer  :: min_nfbnd_cpu,max_nfbnd_cpu
  integer  :: min_nface_cpu,max_nface_cpu
  integer  :: min_nspts_cpu,max_nspts_cpu
  integer  :: min_ncv,max_ncv,min_ncv_cpu,max_ncv_cpu
  real(r8) :: vol,loc_vol,tot_vol
  real(wp) :: ave_ncell,ave_nnode,ave_vol,ave_cell_vol
  real(wp) :: ave_nfbnd,ave_nface,ave_nspts
  real(wp) :: jac_min,jac_max,vol_min,vol_max
  !
  !.. Local Arrays ..
  integer  :: iom(1:4)
  integer  :: part_info(1:5)
  integer  :: all_part_info(1:5,1:ncpu)
  real(wp) :: txyz(1:maxSP,1:nr)
  real(wp) :: outmet(1:nr,1:nr)
  real(wp) :: jac_info(1:11)
  real(wp) :: all_jac_info(1:11,1:ncpu)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_metrics"
  character(len=*), parameter :: fname(4) = ["metrics.fwd.dat", &
                                             "metrics.rev.dat", &
                                             "uniform.fwd.dat", &
                                             "uniform.rev.dat"]
  !
 !logical(lk), parameter :: output_metrics = true
  logical(lk), parameter :: output_metrics = fals
  !
  character(len=*), parameter :: cxyz(3,2) = reshape( ["X","Y","Z", &
                                                       "R","S","T"], &
                                                      shape(cxyz) )
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  jac_info(:) = zero
  all_jac_info(:,:) = zero
  !
  if (VortexGrid_random_perturb) then
    ibeg = 1
    iend = 2
  else
    ibeg = 3
    iend = 4
  end if
  if (output_metrics) then
    do i = ibeg,iend
      open(newunit=iom(i),file=fname(i),status="replace",action="write", &
                          form="formatted",position="rewind", &
                          iostat=ierr,iomsg=error_message)
      call io_error(pname,fname(i),1,__LINE__,__FILE__,ierr,error_message)
    end do
  end if
  !
  ! Allocate the required metrics related arrays
  !
  allocate ( metrics_dt(1:ncell) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"metrics_dt",1,__LINE__,__FILE__,ierr,error_message)
  !
  local_volume = zero
  total_volume = zero
  !
  min_jsp = 1
  max_jsp = 1
  min_jcl = 1
  max_jcl = 1
  jac_min = huge(zero)
  jac_max = -huge(zero)
  !
  min_ncv = 1
  max_ncv = 1
  vol_min = huge(zero)
  vol_max = -huge(zero)
  !
  ! Loop over the cells to get the metrics at each solution point
  !
  metrics_loop: do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    np = n2-n1+1
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    if (output_metrics) then
      do i = ibeg,iend
        write (iom(i),21) nc
      end do
    end if
    !
    ! Create transpose of the coordinates for this cell to reduce cache misses
    !
    do k = n1,n2
      do l = 1,nr
        txyz(k-n1+1,l) = xyz(l,k)
      end do
    end do
    !
    ! Compute the metrics and jacobian for this cell
    !
    call metrics_dt(nc)%create_cell_metrics(txyz(1:np,1:nr), &
                                            this_geom,this_order,nc)
    !
    ! Sum the product of the Jacobian and quadrature weight for each point
    ! in this cell to compute the volume of the current cell
    !
    vol = metrics_dt(nc)%vol(this_geom,this_order)
    !
    ! Add the volume of the current cell to the total volume
    ! of the current local partition
    !
    local_volume = local_volume + vol
    !
    if (vol < vol_min) then
      vol_min = vol
      min_ncv = nc
    end if
    if (vol > vol_max) then
      vol_max = vol
      max_ncv = nc
    end if
    !
    if (minval(metrics_dt(nc)%jac) < jac_min) then
      min_jsp = minloc(metrics_dt(nc)%jac,dim=1)
      jac_min = metrics_dt(nc)%jac(min_jsp)
      min_jcl = nc
    end if
    if (maxval(metrics_dt(nc)%jac) > jac_max) then
      max_jsp = maxloc(metrics_dt(nc)%jac,dim=1)
      jac_max = metrics_dt(nc)%jac(max_jsp)
      max_jcl = nc
    end if
    !
    ! Output the metrics if needed
    !
    if (output_metrics) then
      do n = 1,np
        do i = ibeg,iend
          j1 = mod(i+1,2) + 1
          j2 = mod(i  ,2) + 1
          if (j1 == 1) then
            outmet = invert_matrix( metrics_dt(nc)%met(n,1:nr,1:nr) )
          else
            outmet = metrics_dt(nc)%met(n,1:nr,1:nr) * &
                     metrics_dt(nc)%jac(n)
          end if
          where (abs(outmet) < eps14) outmet = zero
          write (iom(i),22) n+n1-1
          write (iom(i),23) ((cxyz(k,j1),k=1,nr), &
                              cxyz(l,j2), &
                             (outmet(k,l), k=1,nr), l=1,nr)
        end do
      end do
    end if
    !
  end do metrics_loop
  !
  jac_info(1) = real(min_jsp,kind=wp)
  jac_info(2) = real(min_jcl,kind=wp)
  jac_info(3) = real(max_jsp,kind=wp)
  jac_info(4) = real(max_jcl,kind=wp)
  jac_info(5) = jac_min
  jac_info(6) = jac_max
  jac_info(7) = local_volume
  jac_info(8) = real(min_ncv,kind=wp)
  jac_info(9) = real(max_ncv,kind=wp)
  jac_info(10) = vol_min
  jac_info(11) = vol_max
  !
  part_info(1) = ncell
  part_info(2) = nnode
  part_info(3) = nfbnd
  part_info(4) = nface
  part_info(5) = n_solpts
  !
  if (ncpu > 1) then
    !
    call mpi_allgather(    jac_info,size(jac_info,kind=int_mpi),mpi_flttyp, &
                       all_jac_info,size(jac_info,kind=int_mpi),mpi_flttyp, &
                       MPI_COMM_WORLD,mpierr)
    !
    call mpi_allgather(    part_info,size(part_info,kind=int_mpi),mpi_inttyp, &
                       all_part_info,size(part_info,kind=int_mpi),mpi_inttyp, &
                       MPI_COMM_WORLD,mpierr)
    !
  else
    !
    all_jac_info(:,1) = jac_info(:)
    all_part_info(:,1) = part_info(:)
    !
  end if
  !
  min_jac_cpu = minloc( all_jac_info(5,:) , dim=1 )
  max_jac_cpu = maxloc( all_jac_info(6,:) , dim=1 )
  min_vol_cpu = minloc( all_jac_info(7,:) , dim=1 )
  max_vol_cpu = maxloc( all_jac_info(7,:) , dim=1 )
  min_ncv_cpu = minloc( all_jac_info(10,:) , dim=1 )
  max_ncv_cpu = maxloc( all_jac_info(11,:) , dim=1 )
  !
  min_jsp = nint( all_jac_info(1,min_jac_cpu) )
  min_jcl = nint( all_jac_info(2,min_jac_cpu) )
  max_jsp = nint( all_jac_info(3,max_jac_cpu) )
  max_jcl = nint( all_jac_info(4,max_jac_cpu) )
  !
  min_ncv = nint( all_jac_info(8,min_ncv_cpu) )
  max_ncv = nint( all_jac_info(9,max_ncv_cpu) )
  !
  min_ncell_cpu = minloc( all_part_info(1,:) , dim=1 )
  max_ncell_cpu = maxloc( all_part_info(1,:) , dim=1 )
  min_nnode_cpu = minloc( all_part_info(2,:) , dim=1 )
  max_nnode_cpu = maxloc( all_part_info(2,:) , dim=1 )
  min_nfbnd_cpu = minloc( all_part_info(3,:) , dim=1 )
  max_nfbnd_cpu = maxloc( all_part_info(3,:) , dim=1 )
  min_nface_cpu = minloc( all_part_info(4,:) , dim=1 )
  max_nface_cpu = maxloc( all_part_info(4,:) , dim=1 )
  min_nspts_cpu = minloc( all_part_info(5,:) , dim=1 )
  max_nspts_cpu = maxloc( all_part_info(5,:) , dim=1 )
  !
  ave_ncell = real( sum( all_part_info(1,:) ) , kind=kind(ave_ncell) ) / &
              real( size(all_part_info,dim=2) , kind=kind(ave_ncell) )
  ave_nnode = real( sum( all_part_info(2,:) ) , kind=kind(ave_nnode) ) / &
              real( size(all_part_info,dim=2) , kind=kind(ave_nnode) )
  ave_nfbnd = real( sum( all_part_info(3,:) ) , kind=kind(ave_nfbnd) ) / &
              real( size(all_part_info,dim=2) , kind=kind(ave_nfbnd) )
  ave_nface = real( sum( all_part_info(4,:) ) , kind=kind(ave_nface) ) / &
              real( size(all_part_info,dim=2) , kind=kind(ave_nface) )
  ave_nspts = real( sum( all_part_info(5,:) ) , kind=kind(ave_nspts) ) / &
              real( size(all_part_info,dim=2) , kind=kind(ave_nspts) )
  ave_vol   = sum( all_jac_info(7,:) ) / &
              real( size(all_jac_info(7,:)) , kind=kind(ave_vol) )
  !
  ave_cell_vol = sum( all_jac_info(7,:) ) / &
                 real( sum( all_part_info(1,:) ) , kind=kind(ave_cell_vol) )
  !
  if (mypnum == 0) then
    !
    write (iout,1) all_jac_info( 5,min_jac_cpu),min_jsp,min_jcl,min_jac_cpu-1, &
                   all_jac_info( 6,max_jac_cpu),max_jsp,max_jcl,max_jac_cpu-1, &
                   all_jac_info(10,min_ncv_cpu),min_ncv,min_ncv_cpu-1, &
                   all_jac_info(11,max_ncv_cpu),max_ncv,max_ncv_cpu-1, &
                   ave_cell_vol
    !
    write (iout,2) all_part_info(1,min_ncell_cpu),min_ncell_cpu, &
                   all_part_info(1,max_ncell_cpu),max_ncell_cpu, &
                   nint(ave_ncell), &
                   all_part_info(2,min_nnode_cpu),min_nnode_cpu, &
                   all_part_info(2,max_nnode_cpu),max_nnode_cpu, &
                   nint(ave_nnode), &
                   all_part_info(3,min_nfbnd_cpu),min_nfbnd_cpu, &
                   all_part_info(3,max_nfbnd_cpu),max_nfbnd_cpu, &
                   nint(ave_nfbnd), &
                   all_part_info(4,min_nface_cpu),min_nface_cpu, &
                   all_part_info(4,max_nface_cpu),max_nface_cpu, &
                   nint(ave_nface), &
                   all_part_info(5,min_nspts_cpu),min_nspts_cpu, &
                   all_part_info(5,max_nspts_cpu),max_nspts_cpu, &
                   nint(ave_nspts), &
                   all_jac_info(7,min_vol_cpu),min_vol_cpu, &
                   all_jac_info(7,max_vol_cpu),max_vol_cpu, &
                   ave_vol
    !
  end if
  !
  1 format (/, &
     " Min Interior  Jacobian =",es16.8, &
                 " @ SP ",i0,t52,"of Cell ",i0,t69,"on CPU # ",i0,/, &
     " Max Interior  Jacobian =",es16.8, &
                 " @ SP ",i0,t52,"of Cell ",i0,t69,"on CPU # ",i0,//, &
     " Min Interior Cell Volume =",es16.8," @ Cell ",i0,t62,"on CPU # ",i0,/, &
     " Max Interior Cell Volume =",es16.8," @ Cell ",i0,t62,"on CPU # ",i0,/, &
     " Ave Interior Cell Volume =",es16.8)
  !
  2 format (/, &
          " Min ncell = ",i0,t25,"on CPU # ",i0,/, &
          " Max ncell = ",i0,t25,"on CPU # ",i0,/, &
          " Ave ncell = ",i0,//, &
          " Min nnode = ",i0,t25,"on CPU # ",i0,/, &
          " Max nnode = ",i0,t25,"on CPU # ",i0,/, &
          " Ave nnode = ",i0,//, &
          " Min nfbnd = ",i0,t25,"on CPU # ",i0,/, &
          " Max nfbnd = ",i0,t25,"on CPU # ",i0,/, &
          " Ave nfbnd = ",i0,//, &
          " Min nface = ",i0,t25,"on CPU # ",i0,/, &
          " Max nface = ",i0,t25,"on CPU # ",i0,/, &
          " Ave nface = ",i0,//, &
          " Min n_solpts = ",i0,t28,"on CPU # ",i0,/, &
          " Max n_solpts = ",i0,t28,"on CPU # ",i0,/, &
          " Ave n_solpts = ",i0,//, &
          " Min partition volume =",es16.8," on CPU # ",i0,/, &
          " Max partition volume =",es16.8," on CPU # ",i0,/, &
          " Ave partition volume =",es16.8,/)
  !
  ! Get the face normal at the solution points on each face
  !
  call compute_face_normals
  !
  ! Reduce the volumes of each partition to the root
  ! processor to determine the volume of the global grid.
  !
  tot_vol = real(local_volume,kind=r8)
  loc_vol = real(local_volume,kind=r8)
  !
  if (ncpu > 1) then
    call mpi_allreduce(loc_vol,tot_vol,1_int_mpi,mpi_flttyp, &
                       MPI_SUM,MPI_COMM_WORLD,mpierr)
  end if
  !
  total_volume = real(tot_vol,kind=wp)
  local_volume = real(loc_vol,kind=wp)
  !
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  ! Close the files associated with output_metrics
  !
  if (output_metrics) then
    do i = ibeg,iend
      close (iom(i),iostat=ierr,iomsg=error_message)
      call io_error(pname,fname(i),2,__LINE__,__FILE__,ierr,error_message)
    end do
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  21 format (/,60("-"),/," Cell ",i0,/,60("-"))
  22 format (5x,"Solution Point ",i4)
  23 format (3(2x,"d(",a,",",a,",",a,")/d",a,"  :  ",3es14.6,:,/))
  !
end subroutine get_metrics
!
!###############################################################################
!
subroutine compute_face_normals
  !
  !.. Use Statements ..
  use order_mod,         only : maxSP,maxFP
  use order_mod,         only : geom_solpts
  use geovar,            only : ncell,cell
  use geovar,            only : nface,face
  use geovar,            only : nr,xyz
  use interpolation_mod, only : interp
  use derivatives_mod,   only : deriv
  !
  !.. Local Scalars ..
  integer :: n,nc,nf,ncf,n1,n2,np
  integer :: kf,lf,l,k,nfp,ndim
  integer :: left_geom,left_order
  integer :: face_geom,face_order
  integer :: inrm,ierr,npcount,nfpts
  character(len=200) :: array_name
  !
  !.. Local Arrays ..
  real(wp) :: dvec(1:3,1:3)
  real(wp) :: txyz(1:maxSP,1:nr)
  real(wp) :: xyz_face(1:maxFP,1:nr)
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: ipface(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_face_normals"
  character(len=*), parameter :: fname = "face_normals.tec.dat"
 !logical(lk), parameter :: output_face_normals = true
  logical(lk), parameter :: output_face_normals = fals
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Allocate the outer part of the geofa array
  !
  allocate ( geofa(1:nface) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"geofa",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Allocate the data components of the geofa array
  !
  do nf = 1,nface
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    nfpts = geom_solpts(face_geom,face_order)
    !
    allocate ( geofa(nf)%v(1:nr,1:nfpts) , source=zero , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) "geofa",nf
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
    !
  end do
  !
  if (output_face_normals) then
    !
    npcount = 0
    nfpts = 0
    do nf = 1,nface
      nfpts = nfpts + size(geofa(nf)%v,dim=2)
    end do
    !
    open (newunit=inrm,file=fname,status="replace",action="write", &
                       form="formatted",position="rewind",iostat=ierr, &
                       iomsg=error_message)
    call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
    !
    if (nr == 2) then
      !
      write (inrm,20) nfpts,nface
      !
      allocate ( ipface(1:2,1:nface) , source=0 , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"ipface",1,__LINE__,__FILE__,ierr,error_message)
      !
    else if (nr == 3) then
      !
      write (inrm,21) nfpts,nface
      !
      allocate ( ipface(1:4,1:nface) , source=0 , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"ipface",1,__LINE__,__FILE__,ierr,error_message)
      !
    end if
    !
  end if
  !
  cell_loop: do nc = 1,ncell
    !
    left_geom  = cell(nc)%geom
    left_order = cell(nc)%order
    !
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    np = n2-n1+1
    !
    do k = 1,np
      do l = 1,nr
        txyz(k,l) = xyz(l,k+n1-1)
      end do
    end do
    !
    ! Loop through the faces of the current cell
    ! to get the normals at each face point
    !
    cell_face_loop: do ncf = 1,size(cell(nc)%face)
      !
      if (cell(nc)%face(ncf)%side == 2) cycle cell_face_loop
      !
      nf = cell(nc)%face(ncf)%idx
      !
      face_geom  = cell(nc)%face(ncf)%geom
      face_order = cell(nc)%face(ncf)%order
      !
      nfp = geom_solpts(face_geom,face_order)
      !
      ! Interpolate the coordinates of the solution points to get the
      ! coordinates of the flux points on the current face.
      !
      do l = 1,nr
        !
        do kf = 1,nfp
          !
          lf = cell(nc)%face(ncf)%get_fp_idx(kf)
          !
          xyz_face(kf,l) = interp(left_geom,left_order)% &
                                     toFace(face_order)% &
                           dot( lf , txyz(1:np,l) )
          !
        end do
        !
      end do
      !
      ! If the current face is an edge, initialize the k components of the
      ! two vectors the dvec array to one and everything else to zero
      !
      if (face_geom == Geom_Edge) then
        dvec = zero
        dvec(3,:) = one
      end if
      !
      ! Now compute the derivative of the flux point coordinates
      !
      ndim = size(deriv(face_geom,face_order)%d)
      !
      do k = 1,nfp
        !
        do n = 1,ndim
          !
          ! If the current cell face is a x-min, y-max, or z-min face, we need
          ! to reverse the order of the derivatives so that the cross_product
          ! computed below results in a outward facing normal for the left cell
          ! on the face.
          !
          kf = n
          if (face_geom /= Geom_Edge) then
            if (any(ncf == [Xmin3D_Face,Ymax3D_Face,Zmin3D_Face])) then
              kf = ndim+1 - n
            end if
          end if
          !
          do l = 1,nr
            dvec(l,kf) = deriv(face_geom,face_order)% &
                                                d(n)% &
                         dot( k , xyz_face(1:nfp,l) )
          end do
          !
        end do
        !
        ! Finally, the length/area normal to the current face at the current
        ! flux point is just the cross product of the two vectors in dvec
        !
        dvec(:,3) = cross_product( dvec(:,1) , dvec(:,2) )
        !
        ! Zero/chop all the normals less than 10*roundoff
        !
        geofa(nf)%v(1:nr,k) = chop( real(dvec(1:nr,3),kind=qp)/qten ) * ten
        !
        if (output_face_normals) then
          !
          npcount = npcount + 1
          !
          if (nr == 2) then
            !
            if (k == 1) then
              ipface(1,nf) = npcount
            else if (k == nfp) then
              ipface(2,nf) = npcount
            end if
            !
            write (inrm,22) (xyz_face(k,l),l=1,nr), &
                            (geofa(nf)%v(l,k),l=1,nr)!,npcount
            !
          else if (nr == 3) then
            !
            if (k == 1) then
              ipface(1,nf) = npcount
            else if (k == face_order+1) then
              ipface(2,nf) = npcount
            else if (k == nfp) then
              ipface(3,nf) = npcount
            else if (face_geom == Geom_Quad) then
              if (k == (face_order+1)*face_order) then
                ipface(4,nf) = npcount
              end if
            end if
            !
            if (face_geom == Geom_Tria) ipface(4,nf) = ipface(3,nf)
            !
            write (inrm,23) (xyz_face(k,l),l=1,nr), &
                            (geofa(nf)%v(l,k),l=1,nr)!,npcount
            !
          end if
          !
        end if
        !
      end do
      !
    end do cell_face_loop
    !
  end do cell_loop
  !
  if (output_face_normals) then
    !
    do nf = 1,nface
      write (inrm,24) (ipface(l,nf),l=1,size(ipface,dim=1))
    end do
    !
    close (inrm,iostat=ierr,iomsg=error_message)
    call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
    !
    if (allocated(ipface)) then
      deallocate ( ipface , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"ipface",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1  format (a,"(",i0,")%v")
  20 format ('VARIABLES = "X", "Y", "n_x", "n_y"',/, &
             'ZONE T="Face Normals", N=',i0,', E=',i0, &
             ', DATAPACKING=POINT, ZONETYPE=FELINESEG')
  21 format ('VARIABLES = "X", "Y", "Z", "n_x", "n_y", "n_z"',/, &
             'ZONE T="Face Normals", N=',i0,', E=',i0, &
             ', DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL')
  22 format (2es23.15,1x,2es23.15,:,1x,i0)
  23 format (3es14.6,1x,3es14.6,:,1x,i0)
  24 format (4(1x,i0))
  !
end subroutine compute_face_normals
!
!###############################################################################
!
pure function compute_cell_volume(this,geom,order) result(return_value)
  !
  !.. Use Statements ..
  use quadrature_mod, only : std_elem
  !
  !.. Passed Argument from Invoking Object ..
  class(metrics_t), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer, intent(in) :: geom
  integer, intent(in) :: order
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
continue
  !
  ! Sum the product of the Jacobian and quadrature weight for each point
  ! in this cell to compute the volume of the current cell
  !
  return_value = dot( this%jac , std_elem(geom,order)%wts )
  !
end function compute_cell_volume
!
!###############################################################################
!
subroutine create_cell_metrics(this,txyz,geom,order,cell_idx)
  !
  !.. Use Statements ..
  use derivatives_mod, only : deriv
  !
  !.. Passed Argument from Invoking Object ..
  class(metrics_t), intent(inout) :: this
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: txyz(:,:)
  integer,  intent(in) :: geom
  integer,  intent(in) :: order
  integer,  intent(in) :: cell_idx
  !
  !.. Local Scalars ..
  integer :: n,np,ierr
  character(len=50) :: array_name
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_cell_metrics"
  !
continue
  !
  np = size(txyz,dim=1)
  !
  ! Compute the metrics for all points in this cell
  !
  ! USING F2003 AUTO-REALLOCATION
  this%met = deriv(geom,order)%grad( txyz(:,:) )
 !do ref_dir = 1,nr
 !  do phys_dir = 1,nr
 !    metrics(phys_dir,ref_dir,n1:n2) = deriv(this_geom,this_order)% &
 !                                                       d(ref_dir)% &
 !                                      vm_mult( txyz(1:np,phys_dir) )
 !  end do
 !end do
  !
  ! Compute the Jacobian of the metrics tensor for each point in this cell
  !
  allocate ( this%jac(1:np) , source=zero , stat=ierr , errmsg=error_message )
  write (array_name,1) cell_idx
  call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
  !
  do n = 1,np
    this%jac(n) = compute_jacobian( this%met(n,:,:) )
  end do
  !
  ! Compute the inverse metrics for each point in this cell
  !
  do n = 1,np
    this%met(n,:,:) = reverse_metrics( this%met(n,:,:) ) / &
                      this%jac(n)
  end do
  !
  ! Compute the first derivative of Jacobian for all
  ! points in this cell if the array is allocated
  !
  if (use_dJ) then
    ! USING F2003 AUTO-REALLOCATION
    this%dJ = deriv(geom,order)%grad( this%jac(:) )
   !do ref_dir = 1,nr
   !  dJ(ref_dir,n1:n2) = deriv(this_geom,this_order)% &
   !                                       d(ref_dir)% &
   !                      mv_mult( jacobian(n1:n2) )
   !end do
  end if
  !
  ! Format Statements
  !
  1 format ("metrics_dt(",i0,")%jac")
  !
end subroutine create_cell_metrics
!
!###############################################################################
!
subroutine metrics_memory_usage(iunit)
  !
  !.. Use Statements ..
  use module_memory
  !
  !.. Formal Arguments ..
  integer, optional, intent(in) :: iunit
  !
  !.. Local Scalars ..
  integer :: i,iou
  character(len=36) :: array_name
  type(memory) :: array,dt
  type(memory) :: metmet,metjac,metdj
  type(memory) :: total
  !
#ifndef DISABLE_DTIO
continue
  !
  iou = iout
  if (present(iunit)) iou = iunit
  !
  call total%reset
  call metmet%reset
  call metjac%reset
  call metdj%reset
  !
  write (iou,1)
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(geofa)) then
    !
    call dt%set( storage_size(geofa,kind=inttype) , &
                         size(geofa,kind=inttype) )
    !
    write (array_name,5) "geofa(:)"
    write (iou,2) array_name,dt%get_units()
    !
    call array%reset
    !
    do i = lbound(geofa,dim=1),ubound(geofa,dim=1)
      !
      call array%add( storage_size(geofa(i)%v,kind=inttype) , &
                              size(geofa(i)%v,kind=inttype) )
      !
    end do
    !
    write (array_name,5) "geofa(:)%v"
    write (iou,6) array_name,array%get_units()
    !
    call dt%add(array)
    !
    call total%add(dt)
    write (array_name,5) "geofa"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "geofa"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(metrics_dt)) then
    !
    call dt%set( storage_size(metrics_dt,kind=inttype) , &
                         size(metrics_dt,kind=inttype) )
    !
    write (array_name,5) "metrics_dt(:)"
    write (iou,2) array_name,dt%get_units()
    !
    call array%reset
    !
    do i = lbound(metrics_dt,dim=1),ubound(metrics_dt,dim=1)
      !
      if (allocated(metrics_dt(i)%met)) then
        call array%set( storage_size(metrics_dt(i)%met,kind=inttype) , &
                                size(metrics_dt(i)%met,kind=inttype) )
        call metmet%add(array)
        call dt%add(array)
      end if
      !
      if (allocated(metrics_dt(i)%jac)) then
        call array%set( storage_size(metrics_dt(i)%jac,kind=inttype) , &
                                size(metrics_dt(i)%jac,kind=inttype) )
        call metjac%add(array)
        call dt%add(array)
      end if
      !
      if (allocated(metrics_dt(i)%dJ)) then
        call array%set( storage_size(metrics_dt(i)%dJ,kind=inttype) , &
                                size(metrics_dt(i)%dJ,kind=inttype) )
        call metdJ%add(array)
        call dt%add(array)
      end if
      !
    end do
    !
    write (array_name,5) "metrics_dt(:)%met"
    write (iou,6) array_name,metmet%get_units()
    !
    write (array_name,5) "metrics_dt(:)%jac"
    write (iou,6) array_name,metjac%get_units()
    !
    write (array_name,5) "metrics_dt(:)%dJ"
    write (iou,6) array_name,metdJ%get_units()
    !
    call total%add(dt)
    write (array_name,5) "metrics_dt"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "metrics_dt"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module metrics_mod")
  2 format (" Array: ",a,dt)
  3 format (" Total memory = ",f8.2,1x,a)
  5 format ("'",a,"'")
  6 format (" Array Component: ",a,dt)
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  !
#endif
end subroutine metrics_memory_usage
!
!###############################################################################
!
end module metrics_mod
